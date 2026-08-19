"""Tests for building block preprocessing."""

import json
import zipfile

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from healer.scripts import preprocess_bb_source as prep


def write_sdf(path, smiles_list, props=None):
    writer = Chem.SDWriter(str(path))
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        AllChem.Compute2DCoords(mol)
        mol.SetProp("_Name", f"BB{i}")
        for key, value in (props or {}).items():
            mol.SetProp(key, value)
        writer.write(mol)
    writer.close()
    return path


class TestAnnotations:
    def test_a_carboxylic_acid_is_annotated(self):
        annotated = prep.add_rxn_annotations(Chem.MolFromSmiles("CC(=O)O"))
        parsed = json.loads(annotated.GetProp("rxn_annotations"))

        assert parsed, "an acid should match at least one template"
        assert all(isinstance(v, list) for v in parsed.values())

    def test_positions_are_reactant_indices(self):
        annotated = prep.add_rxn_annotations(Chem.MolFromSmiles("NCc1ccccc1"))
        parsed = json.loads(annotated.GetProp("rxn_annotations"))
        for positions in parsed.values():
            assert all(p in (0, 1) for p in positions)

    def test_an_unreactive_molecule_gets_an_empty_mapping(self):
        annotated = prep.add_rxn_annotations(Chem.MolFromSmiles("CC"))
        assert json.loads(annotated.GetProp("rxn_annotations")) == {}

    def test_only_matching_reactions_are_recorded(self):
        """
        Empty matches must be left out. BBRepository indexes on the keys, so
        storing them would put every block under every reaction.
        """
        annotated = prep.add_rxn_annotations(Chem.MolFromSmiles("CC(=O)O"))
        parsed = json.loads(annotated.GetProp("rxn_annotations"))

        assert all(positions for positions in parsed.values())
        assert len(parsed) < len(prep.REACTIONS)

    def test_reaction_filtering_stays_selective(self, tmp_path, all_reactions):
        """
        A single reaction must not match every block. It would if empty matches
        were stored, since BBRepository indexes on the annotation keys.
        """
        from healer.domain.bb_repository import BBRepository

        source = write_sdf(tmp_path / "in.sdf", ["CC(=O)O", "NCc1ccccc1", "c1ccccc1Br"])
        prep.main(str(source), output_dir=str(tmp_path), verbose=False)

        repo = BBRepository(source_path=str(tmp_path / "in_processed.sdf"))
        repo.load(show_progress=False)

        per_reaction = [len(repo.get_bbs_for_reactions([rxn])) for rxn in all_reactions]
        assert any(0 < n < len(repo) for n in per_reaction)

    def test_annotations_match_what_the_template_reports(self, amide_reaction):
        acid = Chem.MolFromSmiles("CC(=O)O")
        parsed = json.loads(prep.add_rxn_annotations(acid).GetProp("rxn_annotations"))
        if amide_reaction.name in parsed:
            assert parsed[amide_reaction.name] == amide_reaction.get_reactant_index(
                acid
            )


class TestSaltStripping:
    def test_the_largest_fragment_is_kept(self):
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O.[Na+].[Cl-]")
        largest = prep.remove_smaller_fragments(mol)
        assert Chem.MolToSmiles(largest) == "CC(=O)Oc1ccccc1C(=O)O"

    def test_a_single_component_is_returned_unchanged(self):
        mol = Chem.MolFromSmiles("CCO")
        assert Chem.MolToSmiles(prep.remove_smaller_fragments(mol)) == "CCO"


class TestSdfIteration:
    def test_records_are_split_on_the_delimiter(self, tmp_path):
        path = write_sdf(tmp_path / "in.sdf", ["CCO", "CCC", "c1ccccc1"])
        records = list(prep._iter_sdf_records(str(path)))
        assert len(records) == 3
        assert all(r.rstrip().endswith("$$$$") for r in records)

    def test_chunks_group_the_stream(self):
        chunks = list(prep._ichunk(range(7), 3))
        assert [len(c) for c in chunks] == [3, 3, 1]


class TestProcessMol:
    def test_a_record_round_trips_with_annotations(self, tmp_path):
        path = write_sdf(tmp_path / "in.sdf", ["CC(=O)O"])
        block = next(prep._iter_sdf_records(str(path)))

        processed = prep._process_mol(block)
        assert processed is not None
        assert "rxn_annotations" in processed
        assert processed.rstrip().endswith("$$$$")

    def test_an_unparsable_record_is_dropped(self):
        assert prep._process_mol("not a molblock\n$$$$\n") is None

    def test_original_properties_are_preserved(self, tmp_path):
        path = write_sdf(
            tmp_path / "in.sdf", ["CC(=O)O"], props={"URL": "https://x.test"}
        )
        block = next(prep._iter_sdf_records(str(path)))
        assert "https://x.test" in prep._process_mol(block)


class TestZipHandling:
    def test_a_plain_sdf_is_passed_through(self, tmp_path):
        path = write_sdf(tmp_path / "in.sdf", ["CCO"])
        assert prep.extract_zip_if_needed(str(path), verbose=False) == str(path)

    def test_an_sdf_is_extracted_from_a_zip(self, tmp_path):
        sdf = write_sdf(tmp_path / "in.sdf", ["CCO"])
        archive = tmp_path / "bundle.zip"
        with zipfile.ZipFile(archive, "w") as zf:
            zf.write(sdf, arcname="in.sdf")
        sdf.unlink()

        extracted = prep.extract_zip_if_needed(str(archive), verbose=False)
        assert extracted.endswith(".sdf")
        assert Chem.SDMolSupplier(extracted)[0] is not None

    def test_a_zip_without_an_sdf_is_an_error(self, tmp_path):
        archive = tmp_path / "empty.zip"
        with zipfile.ZipFile(archive, "w") as zf:
            zf.writestr("readme.txt", "no molecules here")
        with pytest.raises(FileNotFoundError):
            prep.extract_zip_if_needed(str(archive), verbose=False)


class TestEndToEnd:
    def test_processing_writes_an_annotated_file(self, tmp_path):
        source = write_sdf(tmp_path / "in.sdf", ["CC(=O)O", "NCc1ccccc1"])
        prep.main(str(source), output_dir=str(tmp_path), verbose=False)

        out = tmp_path / "in_processed.sdf"
        assert out.exists()

        mols = [m for m in Chem.SDMolSupplier(str(out)) if m is not None]
        assert len(mols) == 2
        assert all(m.HasProp("rxn_annotations") for m in mols)

    def test_blocks_matching_no_reaction_are_dropped(self, tmp_path):
        """They cannot take part in any enumeration, so they are not written."""
        source = write_sdf(tmp_path / "in.sdf", ["CC(=O)O", "CC"])
        prep.main(str(source), output_dir=str(tmp_path), verbose=False)

        mols = [m for m in Chem.SDMolSupplier(str(tmp_path / "in_processed.sdf")) if m]
        assert len(mols) == 1
        assert Chem.MolToSmiles(mols[0]) == "CC(=O)O"

    def test_the_output_loads_as_a_repository(self, tmp_path):
        """Preprocessed files must be consumable by BBRepository."""
        from healer.domain.bb_repository import BBRepository

        source = write_sdf(tmp_path / "in.sdf", ["CC(=O)O", "NCc1ccccc1"])
        prep.main(str(source), output_dir=str(tmp_path), verbose=False)

        repo = BBRepository(source_path=str(tmp_path / "in_processed.sdf"))
        repo.load(show_progress=False)
        assert len(repo) == 2
        assert repo._reaction_bb_indices

    def test_salts_are_stripped_during_processing(self, tmp_path):
        source = write_sdf(tmp_path / "in.sdf", ["CC(=O)O.[Na+].[Cl-]"])
        prep.main(str(source), output_dir=str(tmp_path), verbose=False)

        mol = next(iter(Chem.SDMolSupplier(str(tmp_path / "in_processed.sdf"))))
        assert len(Chem.GetMolFrags(mol)) == 1

    def test_parallel_workers_give_the_same_output(self, tmp_path):
        smiles = ["CC(=O)O", "NCc1ccccc1", "CCO", "c1ccccc1Br"]
        outputs = []
        for workers, name in ((1, "seq"), (2, "par")):
            folder = tmp_path / name
            folder.mkdir()
            source = write_sdf(folder / "in.sdf", smiles)
            prep.main(
                str(source), output_dir=str(folder), verbose=False, n_workers=workers
            )
            mols = [
                m for m in Chem.SDMolSupplier(str(folder / "in_processed.sdf")) if m
            ]
            outputs.append(sorted(m.GetProp("rxn_annotations") for m in mols))

        assert outputs[0] == outputs[1]
