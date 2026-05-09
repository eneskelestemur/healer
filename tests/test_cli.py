"""
Tests for CLI utilities and the run_enumeration pipeline.

Covers:
- load_input: SMILES string, CSV, SDF, invalid input
- parse_rules: rule string parsing
- run_enumeration: all three healer types (molecule, site, fragment)
  with both n_jobs=1 (sequential) and n_jobs=2 (parallel/loky)
"""
import pytest
import textwrap
import pandas as pd
from pathlib import Path

# conftest.py imports rdkit_monkey_patch first — no need to repeat here.

from healer.cli import load_input, parse_rules, run_enumeration
from healer.domain.bb_repository import BBRepository


# ---------------------------------------------------------------------------
# load_input
# ---------------------------------------------------------------------------

def test_load_input_direct_smiles():
    result = load_input("c1ccccc1")
    assert result == ["c1ccccc1"]


def test_load_input_direct_smiles_complex():
    """Multi-atom SMILES without a file path is returned as-is."""
    smi = "CC(=O)Oc1ccccc1C(=O)O"
    result = load_input(smi)
    assert result == [smi]


def test_load_input_csv_default_column(tmp_path):
    csv_file = tmp_path / "mols.csv"
    csv_file.write_text("smiles\nc1ccccc1\nCC(=O)O\n")
    result = load_input(str(csv_file))
    assert result == ["c1ccccc1", "CC(=O)O"]


def test_load_input_csv_custom_column(tmp_path):
    csv_file = tmp_path / "mols.csv"
    csv_file.write_text("mol\nc1ccccc1\nCC(=O)O\n")
    result = load_input(str(csv_file), column="mol")
    assert result == ["c1ccccc1", "CC(=O)O"]


def test_load_input_sdf(tmp_path):
    """SDF with two valid molecules is read correctly."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    sdf_file = tmp_path / "mols.sdf"
    writer = Chem.SDWriter(str(sdf_file))
    for smi in ["c1ccccc1", "CC(=O)O"]:
        mol = Chem.MolFromSmiles(smi)
        AllChem.Compute2DCoords(mol)
        writer.write(mol)
    writer.close()

    result = load_input(str(sdf_file))
    assert len(result) == 2


def test_load_input_invalid_raises():
    with pytest.raises(ValueError, match="Invalid SMILES or file not found"):
        load_input("NOT_A_SMILES_AND_NOT_A_FILE")


# ---------------------------------------------------------------------------
# parse_rules
# ---------------------------------------------------------------------------

def test_parse_rules_basic():
    result = parse_rules("MW:0:500,HBD:0:5")
    assert result == {"MW": (0, 500), "HBD": (0, 5)}


def test_parse_rules_single():
    result = parse_rules("TPSA:0:140")
    assert result == {"TPSA": (0, 140)}


def test_parse_rules_ignores_malformed_entry():
    """Entries without exactly 3 colon-separated parts are silently skipped."""
    result = parse_rules("MW:0:500,BAD,HBD:0:5")
    assert result == {"MW": (0, 500), "HBD": (0, 5)}


# ---------------------------------------------------------------------------
# Shared helpers for run_enumeration tests
# ---------------------------------------------------------------------------

_BASE_INIT = {
    "bb_source": "test",
    "reaction_tags": "all",
    "sim_threshold": 0.0,
    "max_bbs_per_frag": 10,
    "shuffle_bb_order": False,
    "verbose": 0,
}

_BASE_ENUMERATE = {
    "max_evals_per_comp": None,
    "max_products_per_comp": None,
    "max_total_products": 3,
    "n_jobs": 1,
}

_BASE_RESULTS = {
    "calc_similarity": False,
    "calc_properties": False,
}


def _molecule_init(repo: BBRepository) -> dict:
    return {**_BASE_INIT, "bb_repository": repo}


def _site_init(repo: BBRepository) -> dict:
    d = {**_BASE_INIT, "bb_repository": repo}
    del d["sim_threshold"]
    del d["max_bbs_per_frag"]
    d.update({"rules": {}, "struct_rules": []})
    return d


def _fragment_init(repo: BBRepository) -> dict:
    return {**_BASE_INIT, "bb_repository": repo}


# ---------------------------------------------------------------------------
# run_enumeration — molecule mode
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_jobs", [1, 2])
def test_run_enumeration_molecule(test_bb_repository: BBRepository, tmp_path: Path, n_jobs: int):
    """molecule mode: aspirin; produces a CSV with at least the query molecule row."""
    out = tmp_path / f"out_mol_n{n_jobs}.csv"
    run_enumeration(
        healer_type="molecule",
        smiles_list=["CC(=O)Oc1ccccc1C(=O)O"],  # aspirin
        init_kwargs=_molecule_init(test_bb_repository),
        query_kwargs={
            "n_compositions": 3,
            "retro_tree_depth": 1,
            "min_frag_size": 3,
            "randomize_compositions": False,
            "random_seed": -1,
        },
        enumerate_kwargs={**_BASE_ENUMERATE, "n_jobs": n_jobs},
        results_kwargs=_BASE_RESULTS,
        output_path=str(out),
        verbose=2,  # suppress outer tqdm
    )
    assert out.exists(), "Output CSV was not created"
    df = pd.read_csv(out)
    assert len(df) >= 1


# ---------------------------------------------------------------------------
# run_enumeration — site mode
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_jobs", [1, 2])
def test_run_enumeration_site(test_bb_repository: BBRepository, tmp_path: Path, n_jobs: int):
    """site mode: aspirin with permissive rules; produces a CSV."""
    out = tmp_path / f"out_site_n{n_jobs}.csv"
    run_enumeration(
        healer_type="site",
        smiles_list=["CC(=O)Oc1ccccc1C(=O)O"],  # aspirin
        init_kwargs=_site_init(test_bb_repository),
        query_kwargs={"reactive_sites": None},
        enumerate_kwargs={**_BASE_ENUMERATE, "n_jobs": n_jobs},
        results_kwargs=_BASE_RESULTS,
        output_path=str(out),
        verbose=2,
    )
    assert out.exists(), "Output CSV was not created"
    df = pd.read_csv(out)
    assert len(df) >= 1


# ---------------------------------------------------------------------------
# run_enumeration — fragment mode
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_jobs", [1, 2])
def test_run_enumeration_fragment(test_bb_repository: BBRepository, tmp_path: Path, n_jobs: int):
    """fragment mode: two-fragment SMILES; produces a CSV."""
    out = tmp_path / f"out_frag_n{n_jobs}.csv"
    run_enumeration(
        healer_type="fragment",
        # Two disconnected fragments: a primary amine and a carboxylic acid.
        # FragmentHEALER pairs BBs against each fragment and tries reactions.
        smiles_list=["NCC(=O)O.c1cccnc1"],
        init_kwargs=_fragment_init(test_bb_repository),
        query_kwargs={},  # FragmentHEALER.set_query_mol takes no extra kwargs
        enumerate_kwargs={**_BASE_ENUMERATE, "n_jobs": n_jobs},
        results_kwargs=_BASE_RESULTS,
        output_path=str(out),
        verbose=2,
    )
    assert out.exists(), "Output CSV was not created"
    df = pd.read_csv(out)
    assert len(df) >= 1


# ---------------------------------------------------------------------------
# Parallel path produces the same number of molecules as sequential
# ---------------------------------------------------------------------------

def test_molecule_parallel_matches_sequential(test_bb_repository: BBRepository, tmp_path: Path):
    """n_jobs=2 enumerates the same number of products as n_jobs=1."""
    query_kwargs = {
        "n_compositions": 3,
        "retro_tree_depth": 1,
        "min_frag_size": 3,
        "randomize_compositions": False,
        "random_seed": -1,
    }
    enum_base = {**_BASE_ENUMERATE, "max_total_products": 10}

    out1 = tmp_path / "seq.csv"
    out2 = tmp_path / "par.csv"

    run_enumeration(
        healer_type="molecule",
        smiles_list=["CC(=O)Oc1ccccc1C(=O)O"],
        init_kwargs=_molecule_init(test_bb_repository),
        query_kwargs=query_kwargs,
        enumerate_kwargs={**enum_base, "n_jobs": 1},
        results_kwargs=_BASE_RESULTS,
        output_path=str(out1),
        verbose=2,
    )
    run_enumeration(
        healer_type="molecule",
        smiles_list=["CC(=O)Oc1ccccc1C(=O)O"],
        init_kwargs=_molecule_init(test_bb_repository),
        query_kwargs=query_kwargs,
        enumerate_kwargs={**enum_base, "n_jobs": 2},
        results_kwargs=_BASE_RESULTS,
        output_path=str(out2),
        verbose=2,
    )

    n_seq = len(pd.read_csv(out1))
    n_par = len(pd.read_csv(out2))
    assert n_seq == n_par, (
        f"Sequential produced {n_seq} rows, parallel produced {n_par} rows"
    )


# ---------------------------------------------------------------------------
# Parallelization reduces wall time for large candidate sets
# ---------------------------------------------------------------------------

def test_parallelization_speedup(test_bb_repository: BBRepository, tmp_path: Path):
    """
    n_jobs=2 is faster than n_jobs=1 when the candidate set is large enough
    to amortise loky worker overhead.

    Uses max_bbs_per_frag=-1 (all 100 test BBs) to produce a large candidate
    set (~100 seeds × 100 BBs × n_reactions per composition), then runs two
    compositions uncapped.  A warm-up parallel call initialises the loky
    worker pool before the timed runs, so process-spawn cost is excluded.
    """
    import time

    # Override max_bbs_per_frag to use all 100 BBs — the only test that does this.
    init = {**_molecule_init(test_bb_repository), "max_bbs_per_frag": -1}
    query = {
        "n_compositions": 2,
        "retro_tree_depth": 1,
        "min_frag_size": 3,
        "randomize_compositions": False,
        "random_seed": -1,
    }
    enum_no_cap = {"max_evals_per_comp": None, "max_products_per_comp": None, "max_total_products": None}

    def timed(n_jobs: int, out: Path) -> float:
        t0 = time.perf_counter()
        run_enumeration(
            healer_type="molecule",
            smiles_list=["CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"],  # penicillin
            init_kwargs=init,
            query_kwargs=query,
            enumerate_kwargs={**enum_no_cap, "n_jobs": n_jobs},
            results_kwargs=_BASE_RESULTS,
            output_path=str(out),
            verbose=2,
        )
        return time.perf_counter() - t0

    # Warm up loky worker pool so startup cost is not charged to the timed runs.
    timed(2, tmp_path / "warmup.csv")

    t_seq = timed(1, tmp_path / "seq.csv")
    t_par = timed(2, tmp_path / "par.csv")

    assert t_par < t_seq, (
        f"Parallel (n_jobs=2, {t_par:.1f}s) should be faster than "
        f"sequential (n_jobs=1, {t_seq:.1f}s) with 100 BBs and no product cap."
    )

