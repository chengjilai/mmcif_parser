CARGO ?= cargo
PYTHON ?= python3

.PHONY: fmt test docs docs-rust docs-python clean-docs

fmt:
	$(CARGO) fmt --all

test:
	$(CARGO) test --workspace

# Wrapper around the documentation script so CI/users can run `make docs`.
docs:
	scripts/generate_docs.sh

# Convenience targets for invoking only a subset of the docs pipeline.
docs-rust:
	$(CARGO) doc --workspace --no-deps

# The Python docs rely on `pdoc` and an importable mmcif-parser package.
docs-python:
	PYTHON_BIN=$(PYTHON) scripts/generate_docs.sh

clean-docs:
	rm -rf docs/api/python
