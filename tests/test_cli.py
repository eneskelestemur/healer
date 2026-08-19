"""Tests for the command line interface."""

from pathlib import Path

import pandas as pd
import pytest

# conftest.py imports rdkit_monkey_patch first — no need to repeat here.
from healer.cli import load_input, parse_rules, run_enumeration
from healer.domain.bb_repository import BBRepository
from tests.conftest import ASPIRIN_SMILES


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


_BASE_INIT = {
    "bb_source": "test",
    "reaction_tags": "all",
    "sim_threshold": 0.0,
    "max_bbs_per_frag": 10,
    "shuffle_bb_order": False,
    "show_progress": False,
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


@pytest.mark.parametrize("n_jobs", [1, 2])
def test_run_enumeration_molecule(
    test_bb_repository: BBRepository, tmp_path: Path, n_jobs: int
):
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
        show_progress=False,
    )
    assert out.exists(), "Output CSV was not created"
    df = pd.read_csv(out)
    assert len(df) >= 1


@pytest.mark.parametrize("n_jobs", [1, 2])
def test_run_enumeration_site(
    test_bb_repository: BBRepository, tmp_path: Path, n_jobs: int
):
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
        show_progress=False,
    )
    assert out.exists(), "Output CSV was not created"
    df = pd.read_csv(out)
    assert len(df) >= 1


@pytest.mark.parametrize("n_jobs", [1, 2])
def test_run_enumeration_fragment(
    test_bb_repository: BBRepository, tmp_path: Path, n_jobs: int
):
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
        show_progress=False,
    )
    assert out.exists(), "Output CSV was not created"
    df = pd.read_csv(out)
    assert len(df) >= 1


def test_molecule_parallel_matches_sequential(
    test_bb_repository: BBRepository, tmp_path: Path
):
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
        show_progress=False,
    )
    run_enumeration(
        healer_type="molecule",
        smiles_list=["CC(=O)Oc1ccccc1C(=O)O"],
        init_kwargs=_molecule_init(test_bb_repository),
        query_kwargs=query_kwargs,
        enumerate_kwargs={**enum_base, "n_jobs": 2},
        results_kwargs=_BASE_RESULTS,
        output_path=str(out2),
        show_progress=False,
    )

    n_seq = len(pd.read_csv(out1))
    n_par = len(pd.read_csv(out2))
    assert n_seq == n_par, (
        f"Sequential produced {n_seq} rows, parallel produced {n_par} rows"
    )


@pytest.mark.slow
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
    enum_no_cap = {
        "max_evals_per_comp": None,
        "max_products_per_comp": None,
        "max_total_products": None,
    }

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
            show_progress=False,
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


class TestConfigFiles:
    def test_values_are_read_from_json(self, tmp_path):
        from healer.cli import load_config

        path = tmp_path / "run.json"
        path.write_text('{"bb_source": "test", "max_total": 5}')
        assert load_config(str(path)) == {"bb_source": "test", "max_total": 5}

    def test_command_line_values_win_over_the_config(self):
        import argparse

        from healer.cli import merge_args_with_config

        args = argparse.Namespace(bb_source="US_stock", max_total=None)
        merged = merge_args_with_config(args, {"bb_source": "test", "max_total": 5})

        assert merged.bb_source == "US_stock"
        assert merged.max_total == 5

    def test_unknown_config_keys_are_added(self):
        import argparse

        from healer.cli import merge_args_with_config

        merged = merge_args_with_config(argparse.Namespace(), {"extra": 1})
        assert merged.extra == 1


class TestKwargBuilders:
    def _args(self, **overrides):
        import argparse

        base = {
            "bb_source": "test",
            "reactions": "all",
            "shuffle": False,
            "quiet": True,
            "sim_threshold": 0.5,
            "max_bbs_per_frag": -1,
            "n_compositions": 10,
            "randomize": False,
            "seed": -1,
            "retro_depth": 1,
            "min_frag_size": 3,
            "reactive_sites": None,
            "rules": None,
            "struct_rules": None,
            "max_evals": None,
            "max_products": None,
            "max_total": None,
            "n_jobs": 1,
            "similarity": False,
            "properties": False,
        }
        base.update(overrides)
        return argparse.Namespace(**base)

    def test_all_reactions_stay_a_string(self):
        from healer.cli import get_init_kwargs

        assert get_init_kwargs(self._args(), "molecule")["reaction_tags"] == "all"

    def test_a_tag_list_is_split(self):
        from healer.cli import get_init_kwargs

        kwargs = get_init_kwargs(
            self._args(reactions="amide coupling,N-arylation"), "molecule"
        )
        assert kwargs["reaction_tags"] == ["amide coupling", "N-arylation"]

    def test_quiet_disables_progress(self):
        from healer.cli import get_init_kwargs

        assert (
            get_init_kwargs(self._args(quiet=True), "molecule")["show_progress"]
            is False
        )
        assert (
            get_init_kwargs(self._args(quiet=False), "molecule")["show_progress"]
            is None
        )

    def test_site_mode_carries_the_rules(self):
        from healer.cli import get_init_kwargs

        kwargs = get_init_kwargs(self._args(rules="MW:0:300"), "site")
        assert kwargs["rules"] == {"MW": (0, 300)}
        assert "sim_threshold" not in kwargs

    def test_query_kwargs_differ_by_mode(self):
        from healer.cli import get_query_kwargs

        assert "n_compositions" in get_query_kwargs(self._args(), "molecule")
        assert get_query_kwargs(self._args(reactive_sites=[1]), "site") == {
            "reactive_sites": [1]
        }
        assert get_query_kwargs(self._args(), "fragment") == {}

    def test_limits_and_output_options_are_forwarded(self):
        from healer.cli import get_enumerate_kwargs, get_results_kwargs

        args = self._args(max_total=5, n_jobs=2, similarity=True)
        assert get_enumerate_kwargs(args)["max_total_products"] == 5
        assert get_enumerate_kwargs(args)["n_jobs"] == 2
        assert get_results_kwargs(args)["calc_similarity"] is True


class TestParser:
    def test_each_subcommand_is_available(self):
        from healer.cli import build_parser

        parser = build_parser()
        for command in ("molecule", "site", "fragment", "view"):
            assert parser.parse_args([command, "CCO"]).command == command

    def test_defaults_match_the_documentation(self):
        from healer.cli import build_parser

        args = build_parser().parse_args(["molecule", "CCO"])
        assert args.bb_source == "US_stock"
        assert args.reactions == "all"
        assert args.n_jobs == 1
        assert args.quiet is False

    def test_reactive_sites_are_parsed_as_json(self):
        from healer.cli import build_parser

        args = build_parser().parse_args(["site", "CCO", "--reactive-sites", "[1,2]"])
        assert args.reactive_sites == [1, 2]

    def test_a_missing_command_is_rejected(self):
        from healer.cli import build_parser

        with pytest.raises(SystemExit):
            build_parser().parse_args([])


class TestViewCommand:
    def test_an_svg_is_written_when_asked(self, tmp_path):
        import argparse

        from healer.cli import cmd_view

        out = tmp_path / "mol.svg"
        cmd_view(argparse.Namespace(smiles="CCO", output=str(out)))
        assert out.exists() and out.read_bytes().startswith(b"<?xml")

    def test_an_invalid_smiles_is_reported(self, caplog):
        import argparse
        import logging

        from healer.cli import cmd_view

        with caplog.at_level(logging.ERROR, logger="healer.cli"):
            cmd_view(argparse.Namespace(smiles="not-a-molecule", output=None))
        assert "Invalid SMILES" in caplog.text


class TestRunSummary:
    def test_invalid_molecules_are_skipped_and_counted(
        self, test_bb_repository, tmp_path, caplog
    ):
        import logging

        out = tmp_path / "out.csv"
        with caplog.at_level(logging.WARNING, logger="healer.cli"):
            run_enumeration(
                healer_type="molecule",
                smiles_list=["not-a-molecule", ASPIRIN_SMILES],
                init_kwargs=_molecule_init(test_bb_repository),
                query_kwargs={
                    "n_compositions": 1,
                    "retro_tree_depth": 1,
                    "min_frag_size": 3,
                    "randomize_compositions": False,
                    "random_seed": -1,
                },
                enumerate_kwargs={**_BASE_ENUMERATE, "max_total_products": 2},
                results_kwargs=_BASE_RESULTS,
                output_path=str(out),
                show_progress=False,
            )
        assert "Skipping invalid SMILES" in caplog.text
        assert out.exists()

    def test_a_failing_molecule_does_not_stop_the_run(
        self, test_bb_repository, tmp_path, caplog
    ):
        import logging

        out = tmp_path / "out.csv"
        with caplog.at_level(logging.ERROR, logger="healer.cli"):
            run_enumeration(
                healer_type="fragment",
                smiles_list=["CCO", "c1ccccc1N.CC(=O)O"],  # first has too few fragments
                init_kwargs=_fragment_init(test_bb_repository),
                query_kwargs={},
                enumerate_kwargs={**_BASE_ENUMERATE, "max_total_products": 2},
                results_kwargs=_BASE_RESULTS,
                output_path=str(out),
                show_progress=False,
            )
        assert "Error processing" in caplog.text
        assert out.exists()
