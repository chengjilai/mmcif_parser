use std::collections::BTreeSet;
use std::ffi::OsString;
use std::path::Path;
use std::process::Command;

use mmcif_core::parser::parse_file;
use mmcif_core::AtomSite;
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct BiopythonStats {
    atom_count: usize,
    residue_count: usize,
    chain_ids: Vec<String>,
    residue_names_sample: Vec<String>,
    first_atom_name: String,
    first_residue: String,
}

#[test]
fn compare_against_biopython_reference() {
    let crate_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let cif_path = crate_dir.join("../../11AS.cif.gz");
    assert!(
        cif_path.exists(),
        "11AS.cif.gz missing at {}",
        cif_path.display()
    );

    let stats = run_biopython_summary(crate_dir, &cif_path);
    let structure = parse_file(&cif_path).expect("mmcif-core parses 11AS");

    assert_eq!(structure.atom_count(), stats.atom_count, "atom count diverged");
    assert_eq!(
        structure.metadata.entry_id.as_deref(),
        Some("11AS"),
        "entry id mismatch"
    );

    let our_chain_ids = collect_chain_ids(structure.atoms());
    assert_eq!(our_chain_ids, stats.chain_ids, "chain set mismatch");

    let our_residue_count = collect_residue_keys(structure.atoms()).len();
    assert_eq!(our_residue_count, stats.residue_count, "residue count mismatch");

    let our_samples = collect_residue_names(structure.atoms());
    assert_eq!(our_samples, stats.residue_names_sample, "residue sample mismatch");

    let our_first_atom = structure
        .atoms()
        .first()
        .and_then(atom_name)
        .expect("first atom name");
    assert_eq!(our_first_atom, stats.first_atom_name, "first atom mismatch");

    let our_first_residue = structure
        .atoms()
        .first()
        .and_then(|atom| atom.label_comp_id.as_deref())
        .map(|s| s.trim().to_ascii_uppercase())
        .expect("first residue");
    assert_eq!(our_first_residue, stats.first_residue, "first residue mismatch");

    let cell = structure
        .metadata
        .unit_cell
        .as_ref()
        .expect("unit cell parameters present");
    assert_close(cell.a(), 52.900, 1e-3, "a length");
    assert_close(cell.b(), 126.200, 1e-3, "b length");
    assert_close(cell.c(), 52.780, 1e-3, "c length");
    assert_close(cell.alpha(), 90.0, 1e-4, "alpha angle");
    assert_close(cell.beta(), 105.34, 1e-3, "beta angle");
    assert_close(cell.gamma(), 90.0, 1e-4, "gamma angle");
}

fn collect_chain_ids(atoms: &[AtomSite]) -> Vec<String> {
    let mut chains = BTreeSet::new();
    for atom in atoms {
        if let Some(id) = atom
            .auth_asym_id
            .as_deref()
            .or(atom.label_asym_id.as_deref())
        {
            let trimmed = id.trim();
            if !trimmed.is_empty() {
                chains.insert(trimmed.to_string());
            }
        }
    }
    chains.into_iter().collect()
}

#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct ResidueKey {
    chain_id: Option<String>,
    seq_id: Option<i32>,
    comp_id: Option<String>,
}

fn collect_residue_keys(atoms: &[AtomSite]) -> BTreeSet<ResidueKey> {
    let mut residues = BTreeSet::new();
    for atom in atoms {
        let chain_id = atom
            .auth_asym_id
            .as_deref()
            .or(atom.label_asym_id.as_deref())
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty());
        let seq_id = atom.auth_seq_id.or(atom.label_seq_id);
        let comp_id = atom
            .label_comp_id
            .as_deref()
            .map(|s| s.trim().to_ascii_uppercase());
        residues.insert(ResidueKey {
            chain_id,
            seq_id,
            comp_id,
        });
    }
    residues
}

fn collect_residue_names(atoms: &[AtomSite]) -> Vec<String> {
    let mut names = BTreeSet::new();
    for atom in atoms {
        if let Some(comp) = atom.label_comp_id.as_deref() {
            let trimmed = comp.trim();
            if !trimmed.is_empty() {
                names.insert(trimmed.to_ascii_uppercase());
            }
        }
    }
    names.into_iter().take(10).collect()
}

fn atom_name(atom: &AtomSite) -> Option<String> {
    atom
        .label_atom_id
        .as_deref()
        .or(atom.auth_atom_id.as_deref())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
}

fn run_biopython_summary(crate_dir: &Path, cif_path: &Path) -> BiopythonStats {
    let python = pick_python_interpreter(crate_dir);
    let script_path = crate_dir.join("tests/biopython_summary.py");
    let output = Command::new(&python)
        .arg(&script_path)
        .arg(cif_path)
        .output()
        .unwrap_or_else(|err| panic!("failed to spawn {python:?}: {err}"));

    if !output.status.success() {
        panic!(
            "Biopython helper failed (status {status:?}).\nstdout:\n{stdout}\nstderr:\n{stderr}\nEnsure Biopython is installed (pip install biopython) or set BIOPYTHON_PYTHON to an interpreter that can import Bio.PDB.",
            status = output.status,
            stdout = String::from_utf8_lossy(&output.stdout),
            stderr = String::from_utf8_lossy(&output.stderr),
        );
    }

    serde_json::from_slice(&output.stdout).expect("valid Biopython JSON output")
}

fn pick_python_interpreter(crate_dir: &Path) -> OsString {
    if let Some(explicit) = std::env::var_os("BIOPYTHON_PYTHON") {
        return explicit;
    }
    let repo_root = crate_dir.join("../../");
    let venv_python = repo_root.join(".venv/bin/python");
    if venv_python.exists() {
        return venv_python.into_os_string();
    }
    OsString::from("python3")
}

fn assert_close(actual: f64, expected: f64, tolerance: f64, label: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{} mismatch: got {}, expected {} (tol {})",
        label,
        actual,
        expected,
        tolerance
    );
}
