use std::borrow::Cow;
use std::error::Error as StdError;
use std::fmt;
use std::path::{Path, PathBuf};
use std::{io, result};

use thiserror::Error;

pub type Result<T, E = ParseError> = result::Result<T, E>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
pub enum ParseErrorKind {
    #[error("I/O error")]
    Io,
    #[error("input was not valid UTF-8")]
    Utf8,
    #[error("token stream is invalid")]
    InvalidToken,
    #[error("loop tag/value mismatch")]
    LoopMismatch,
    #[error("numeric value missing")]
    MissingValue,
    #[error("invalid numeric literal")]
    InvalidNumber,
    #[error("_atom_site loop missing")]
    MissingAtomLoop,
    #[error("unsupported feature")]
    Unsupported,
}

#[derive(Debug)]
pub struct ParseError {
    kind: ParseErrorKind,
    message: Cow<'static, str>,
    path: Option<PathBuf>,
    source: Option<Box<dyn StdError + Send + Sync>>,
}

impl ParseError {
    pub fn new(kind: ParseErrorKind, message: impl Into<Cow<'static, str>>) -> Self {
        Self {
            kind,
            message: message.into(),
            path: None,
            source: None,
        }
    }

    pub fn kind(&self) -> ParseErrorKind {
        self.kind
    }

    pub fn message(&self) -> &str {
        &self.message
    }

    pub fn path(&self) -> Option<&Path> {
        self.path.as_deref()
    }

    pub fn with_path(mut self, path: impl AsRef<Path>) -> Self {
        self.path = Some(path.as_ref().to_path_buf());
        self
    }

    pub fn with_source<E>(mut self, source: E) -> Self
    where
        E: StdError + Send + Sync + 'static,
    {
        self.source = Some(Box::new(source));
        self
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.path {
            Some(path) => write!(f, "{}: {} ({})", self.kind, self.message, path.display()),
            None => write!(f, "{}: {}", self.kind, self.message),
        }
    }
}

impl StdError for ParseError {
    fn source(&self) -> Option<&(dyn StdError + 'static)> {
        self.source
            .as_ref()
            .map(|boxed| boxed.as_ref() as &(dyn StdError + 'static))
    }
}

impl From<io::Error> for ParseError {
    fn from(err: io::Error) -> Self {
        ParseError::new(ParseErrorKind::Io, err.to_string()).with_source(err)
    }
}

impl From<std::string::FromUtf8Error> for ParseError {
    fn from(err: std::string::FromUtf8Error) -> Self {
        ParseError::new(ParseErrorKind::Utf8, err.to_string()).with_source(err)
    }
}
