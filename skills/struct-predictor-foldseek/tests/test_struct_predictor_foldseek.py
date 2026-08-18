"""Tests for the struct-predictor-foldseek skill."""
from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Path setup — make the skill importable without package install
# ---------------------------------------------------------------------------
SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

import struct_predictor_foldseek as spf
from struct_predictor_foldseek import (
    DEMO_STRUCTURE,
    DEFAULT_DB_CACHE,
    DEFAULT_MAX_HITS,
    DEFAULT_MIN_TMSCORE,
    DEMO_NAME,
    DISCLAIMER,
    _FOLDSEEK_COLS,
    _DB_ALIASES,
    _build_cmd,
    _build_parser,
    _filter_and_rank_hits,
    _foldseek_version,
    _get_db_path,
    _parse_hits,
    _resolve_databases,
    _tmscore_band,
    run_foldseek_search,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_hit(
    target: str = "1L2Y_A",
    tmscore: float = 0.97,
    rmsd: float = 0.42,
    pident: float = 1.0,
    evalue: float = 1.2e-15,
    alnlen: int = 20,
) -> dict:
    """Return a minimal hit dict as if parsed from Foldseek TSV."""
    return {
        "query": "Trpcage_model_0.cif",
        "target": target,
        "pident": pident,
        "alnlen": alnlen,
        "mismatch": 0,
        "gapopen": 0,
        "qstart": 1,
        "qend": 20,
        "tstart": 1,
        "tend": 20,
        "evalue": evalue,
        "bits": 100.0,
        "alntmscore": tmscore,
        "lddt": 0.85,
    }


def _write_tsv(tmp_path: Path, hits: list[dict], name: str = "results.tsv") -> Path:
    """Write hits to a TSV file formatted like Foldseek easy-search output."""
    path = tmp_path / name
    with open(path, "w", encoding="utf-8") as fh:
        for h in hits:
            row = "\t".join(str(h.get(c, "")) for c in _FOLDSEEK_COLS)
            fh.write(row + "\n")
    return path


# ---------------------------------------------------------------------------
# TestConstants
# ---------------------------------------------------------------------------


class TestConstants:
    def test_demo_structure_exists(self):
        assert DEMO_STRUCTURE.exists(), f"Demo CIF not found: {DEMO_STRUCTURE}"

    def test_demo_name(self):
        assert DEMO_NAME == "Trpcage"

    def test_disclaimer_not_empty(self):
        assert len(DISCLAIMER) > 20

    def test_foldseek_cols_has_tmscore(self):
        assert "tmscore" in _FOLDSEEK_COLS

    def test_foldseek_cols_has_rmsd(self):
        assert "rmsd" in _FOLDSEEK_COLS

    def test_db_aliases_contains_pdb(self):
        assert "pdb" in _DB_ALIASES

    def test_db_aliases_contains_afdb(self):
        assert "afdb" in _DB_ALIASES

    def test_demo_cif_has_atom_site(self):
        """The bundled demo CIF must contain ATOM records for Foldseek."""
        text = DEMO_STRUCTURE.read_text()
        assert "ATOM" in text
        assert "_atom_site" in text


# ---------------------------------------------------------------------------
# TestResolveDatabases
# ---------------------------------------------------------------------------


class TestResolveDatabases:
    def test_pdb_alias(self):
        result = _resolve_databases("pdb")
        assert result == [("pdb", "pdb")]

    def test_afdb_alias(self):
        result = _resolve_databases("afdb")
        assert result == [("afdb", "afdb50")]

    def test_multiple_databases(self):
        result = _resolve_databases("pdb,afdb")
        assert len(result) == 2
        aliases = [r[0] for r in result]
        assert "pdb" in aliases
        assert "afdb" in aliases

    def test_unknown_alias_raises(self):
        with pytest.raises(ValueError, match="Unknown database alias"):
            _resolve_databases("notadb")

    def test_whitespace_trimmed(self):
        result = _resolve_databases(" pdb , afdb ")
        assert len(result) == 2

    def test_empty_after_split_ignored(self):
        result = _resolve_databases("pdb,,afdb")
        assert len(result) == 2


# ---------------------------------------------------------------------------
# TestGetDbPath
# ---------------------------------------------------------------------------


class TestGetDbPath:
    def test_returns_expected_path(self, tmp_path):
        p = _get_db_path("pdb", "pdb", tmp_path)
        assert p == tmp_path / "pdb" / "pdb"


# ---------------------------------------------------------------------------
# TestParseHits
# ---------------------------------------------------------------------------


class TestParseHits:
    def test_empty_file_returns_empty_list(self, tmp_path):
        tsv = tmp_path / "empty.tsv"
        tsv.write_text("")
        assert _parse_hits(tsv) == []

    def test_nonexistent_file_returns_empty_list(self, tmp_path):
        tsv = tmp_path / "nosuchfile.tsv"
        assert _parse_hits(tsv) == []

    def test_single_hit_parsed(self, tmp_path):
        hit = _make_hit("1L2Y_A", tmscore=0.97, rmsd=0.42, pident=1.0)
        tsv = _write_tsv(tmp_path, [hit])
        hits = _parse_hits(tsv)
        assert len(hits) == 1
        assert hits[0]["target"] == "1L2Y_A"
        assert abs(hits[0]["tmscore"] - 0.97) < 1e-6
        assert abs(hits[0]["rmsd"] - 0.42) < 1e-6

    def test_multiple_hits_parsed(self, tmp_path):
        hits_in = [_make_hit(f"PDB_{i}", tmscore=0.5 + i * 0.1) for i in range(5)]
        tsv = _write_tsv(tmp_path, hits_in)
        hits_out = _parse_hits(tsv)
        assert len(hits_out) == 5

    def test_invalid_float_defaults_to_zero(self, tmp_path):
        tsv = tmp_path / "bad.tsv"
        # Write a row with "N/A" in the alntmscore column
        row_vals = ["query", "target", "0.5", "20", "0", "0", "1", "20", "1", "20", "1e-5", "100", "N/A", "0.85"]
        tsv.write_text("\t".join(row_vals) + "\n")
        hits = _parse_hits(tsv)
        assert hits[0]["alntmscore"] == 0.0

    def test_numeric_columns_are_floats_or_ints(self, tmp_path):
        hit = _make_hit()
        tsv = _write_tsv(tmp_path, [hit])
        parsed = _parse_hits(tsv)[0]
        assert isinstance(parsed["alntmscore"], float)
        assert isinstance(parsed["lddt"], float)
        assert isinstance(parsed["alnlen"], int)


# ---------------------------------------------------------------------------
# TestFilterAndRankHits
# ---------------------------------------------------------------------------


class TestFilterAndRankHits:
    def test_filters_below_threshold(self):
        hits = [_make_hit(f"H{i}", tmscore=i * 0.1) for i in range(11)]
        result = _filter_and_rank_hits(hits, min_tmscore=0.5, max_hits=100)
        for h in result:
            assert h["alntmscore"] >= 0.5

    def test_sorted_by_tmscore_descending(self):
        hits = [_make_hit(f"H{i}", tmscore=i * 0.1) for i in range(11)]
        result = _filter_and_rank_hits(hits, min_tmscore=0.0, max_hits=100)
        scores = [h["alntmscore"] for h in result]
        assert scores == sorted(scores, reverse=True)

    def test_max_hits_respected(self):
        hits = [_make_hit(f"H{i}", tmscore=0.9) for i in range(100)]
        result = _filter_and_rank_hits(hits, min_tmscore=0.0, max_hits=10)
        assert len(result) == 10

    def test_empty_input_returns_empty(self):
        assert _filter_and_rank_hits([], min_tmscore=0.3, max_hits=50) == []

    def test_all_filtered_returns_empty(self):
        hits = [_make_hit(f"H{i}", tmscore=0.1) for i in range(5)]
        result = _filter_and_rank_hits(hits, min_tmscore=0.5, max_hits=50)
        assert result == []


# ---------------------------------------------------------------------------
# TestTmscoreBand
# ---------------------------------------------------------------------------


class TestTmscoreBand:
    @pytest.mark.parametrize("score,expected_fragment", [
        (0.95, "Near-identical"),
        (0.80, "Very similar"),
        (0.60, "same fold family"),
        (0.40, "Possibly related"),
        (0.20, "Likely unrelated"),
    ])
    def test_band_labels(self, score, expected_fragment):
        assert expected_fragment.lower() in _tmscore_band(score).lower()

    def test_boundary_0_9(self):
        assert "Near-identical" in _tmscore_band(0.9)

    def test_boundary_0_7(self):
        assert "Very similar" in _tmscore_band(0.7)

    def test_boundary_0_5(self):
        assert "fold family" in _tmscore_band(0.5)

    def test_boundary_0_3(self):
        assert "Possibly" in _tmscore_band(0.3)


# ---------------------------------------------------------------------------
# TestBuildCmd
# ---------------------------------------------------------------------------


class TestBuildCmd:
    def test_demo_mode(self, tmp_path):
        cmd = _build_cmd(DEMO_STRUCTURE, tmp_path, "pdb", 0.3, 50, demo=True)
        assert "--demo" in cmd
        assert "--input" not in cmd

    def test_input_mode(self, tmp_path):
        query = tmp_path / "structure.cif"
        cmd = _build_cmd(query, tmp_path, "pdb", 0.3, 50, demo=False)
        assert "--input" in cmd
        assert str(query) in cmd

    def test_contains_databases(self, tmp_path):
        cmd = _build_cmd(DEMO_STRUCTURE, tmp_path, "pdb,afdb", 0.3, 50, demo=True)
        assert "pdb,afdb" in cmd

    def test_contains_output(self, tmp_path):
        cmd = _build_cmd(DEMO_STRUCTURE, tmp_path, "pdb", 0.3, 50, demo=True)
        assert str(tmp_path) in cmd

    def test_contains_min_tmscore(self, tmp_path):
        cmd = _build_cmd(DEMO_STRUCTURE, tmp_path, "pdb", 0.45, 50, demo=True)
        assert "0.45" in cmd


# ---------------------------------------------------------------------------
# TestCLIParser
# ---------------------------------------------------------------------------


class TestCLIParser:
    def test_demo_flag(self):
        parser = _build_parser()
        args = parser.parse_args(["--demo", "--output", "/tmp/out"])
        assert args.demo is True

    def test_input_flag(self):
        parser = _build_parser()
        args = parser.parse_args(["--input", "structure.cif", "--output", "/tmp/out"])
        assert args.input == "structure.cif"
        assert args.demo is False

    def test_databases_default(self):
        parser = _build_parser()
        args = parser.parse_args(["--demo", "--output", "/tmp/out"])
        assert args.databases == "pdb"

    def test_databases_custom(self):
        parser = _build_parser()
        args = parser.parse_args(["--demo", "--output", "/tmp/out", "--databases", "pdb,afdb"])
        assert args.databases == "pdb,afdb"

    def test_min_tmscore_default(self):
        parser = _build_parser()
        args = parser.parse_args(["--demo", "--output", "/tmp/out"])
        assert abs(args.min_tmscore - DEFAULT_MIN_TMSCORE) < 1e-9

    def test_min_tmscore_custom(self):
        parser = _build_parser()
        args = parser.parse_args(["--demo", "--output", "/tmp/out", "--min-tmscore", "0.5"])
        assert abs(args.min_tmscore - 0.5) < 1e-9

    def test_max_hits_custom(self):
        parser = _build_parser()
        args = parser.parse_args(["--demo", "--output", "/tmp/out", "--max-hits", "20"])
        assert args.max_hits == 20

    def test_output_required(self):
        parser = _build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["--demo"])


# ---------------------------------------------------------------------------
# TestRunFoldseekSearchValidation
# ---------------------------------------------------------------------------


class TestRunFoldseekSearchValidation:
    def test_no_input_no_demo_raises(self, tmp_path):
        with pytest.raises(ValueError, match="--input or --demo"):
            run_foldseek_search(input_path=None, output_dir=tmp_path, demo=False)

    def test_missing_query_file_raises(self, tmp_path):
        missing = tmp_path / "nosuchfile.cif"
        with patch("struct_predictor_foldseek.shutil.which", return_value="/usr/bin/foldseek"):
            with patch("struct_predictor_foldseek._foldseek_version", return_value="9.0"):
                with pytest.raises(FileNotFoundError):
                    run_foldseek_search(
                        input_path=missing, output_dir=tmp_path / "out", demo=False
                    )

    def test_unknown_database_raises(self, tmp_path):
        with patch("struct_predictor_foldseek.shutil.which", return_value="/usr/bin/foldseek"):
            with patch("struct_predictor_foldseek._foldseek_version", return_value="9.0"):
                with pytest.raises(ValueError, match="Unknown database alias"):
                    run_foldseek_search(
                        input_path=None,
                        output_dir=tmp_path / "out",
                        databases="notadb",
                        demo=True,
                    )

    def test_foldseek_not_on_path_raises(self, tmp_path):
        with patch("struct_predictor_foldseek.shutil.which", return_value=None):
            with pytest.raises(RuntimeError, match="foldseek not found"):
                run_foldseek_search(
                    input_path=None, output_dir=tmp_path / "out", demo=True
                )


# ---------------------------------------------------------------------------
# TestRunFoldseekSearchWithMockFoldseek
# ---------------------------------------------------------------------------


class TestRunFoldseekSearchWithMockFoldseek:
    """End-to-end tests mocking subprocess and database paths."""

    def _make_fake_db(self, db_cache: Path, alias: str, foldseek_name: str) -> Path:
        db_path = db_cache / alias / foldseek_name
        db_path.mkdir(parents=True, exist_ok=True)
        return db_path

    def _fake_run(self, hits: list[dict], results_tsv_arg_index: int = 3):
        """Return a side_effect function that writes fake TSV output."""
        def _side_effect(cmd, **kwargs):
            # cmd is [exe, "easy-search", query, db, results_tsv, tmp_dir, ...]
            # results_tsv is the 5th positional arg (index 4)
            results_tsv = Path(cmd[4])
            results_tsv.parent.mkdir(parents=True, exist_ok=True)
            with open(results_tsv, "w", encoding="utf-8") as fh:
                for h in hits:
                    row = "\t".join(str(h.get(c, "")) for c in _FOLDSEEK_COLS)
                    fh.write(row + "\n")
            mock_proc = MagicMock()
            mock_proc.returncode = 0
            mock_proc.stdout = ""
            mock_proc.stderr = ""
            return mock_proc
        return _side_effect

    def test_demo_end_to_end_no_db(self, tmp_path):
        """Demo with a missing database should warn but not crash; result should show 0 hits."""
        output_dir = tmp_path / "out"
        db_cache = tmp_path / "dbs"

        with patch("struct_predictor_foldseek.shutil.which", return_value="/usr/bin/foldseek"):
            with patch("struct_predictor_foldseek._foldseek_version", return_value="9.0"):
                result = run_foldseek_search(
                    input_path=None,
                    output_dir=output_dir,
                    databases="pdb",
                    demo=True,
                    db_cache=db_cache,
                )

        assert (output_dir / "report.md").exists()
        assert (output_dir / "result.json").exists()
        assert result["n_hits_raw"] == 0
        assert result["demo"] is True

    def test_demo_end_to_end_with_hits(self, tmp_path):
        """Demo with a present database returns the mocked hits."""
        output_dir = tmp_path / "out"
        db_cache = tmp_path / "dbs"
        self._make_fake_db(db_cache, "pdb", "pdb")

        fake_hits = [
            _make_hit("1L2Y_A", tmscore=0.97, rmsd=0.42),
            _make_hit("2JOF_A", tmscore=0.72, rmsd=1.10),
            _make_hit("3NIR_B", tmscore=0.25, rmsd=3.80),  # below default threshold
        ]

        with patch("struct_predictor_foldseek.shutil.which", return_value="/usr/bin/foldseek"):
            with patch("struct_predictor_foldseek._foldseek_version", return_value="9.0"):
                with patch(
                    "struct_predictor_foldseek.subprocess.run",
                    side_effect=self._fake_run(fake_hits),
                ):
                    result = run_foldseek_search(
                        input_path=None,
                        output_dir=output_dir,
                        databases="pdb",
                        min_tmscore=DEFAULT_MIN_TMSCORE,
                        demo=True,
                        db_cache=db_cache,
                    )

        assert result["n_hits_raw"] == 3
        assert result["n_hits_filtered"] == 2  # one below 0.3
        assert result["top_hit"]["target"] == "1L2Y_A"
        assert abs(result["top_hit"]["tmscore"] - 0.97) < 1e-6
        assert result["demo"] is True

    def test_report_md_contains_disclaimer(self, tmp_path):
        output_dir = tmp_path / "out"
        db_cache = tmp_path / "dbs"
        self._make_fake_db(db_cache, "pdb", "pdb")

        with patch("struct_predictor_foldseek.shutil.which", return_value="/usr/bin/foldseek"):
            with patch("struct_predictor_foldseek._foldseek_version", return_value="9.0"):
                with patch(
                    "struct_predictor_foldseek.subprocess.run",
                    side_effect=self._fake_run([_make_hit()]),
                ):
                    run_foldseek_search(
                        input_path=None,
                        output_dir=output_dir,
                        demo=True,
                        db_cache=db_cache,
                    )

        report = (output_dir / "report.md").read_text()
        assert DISCLAIMER in report

    def test_result_json_fields(self, tmp_path):
        output_dir = tmp_path / "out"
        db_cache = tmp_path / "dbs"
        self._make_fake_db(db_cache, "pdb", "pdb")

        with patch("struct_predictor_foldseek.shutil.which", return_value="/usr/bin/foldseek"):
            with patch("struct_predictor_foldseek._foldseek_version", return_value="9.0"):
                with patch(
                    "struct_predictor_foldseek.subprocess.run",
                    side_effect=self._fake_run([_make_hit()]),
                ):
                    run_foldseek_search(
                        input_path=None,
                        output_dir=output_dir,
                        demo=True,
                        db_cache=db_cache,
                    )

        result = json.loads((output_dir / "result.json").read_text())
        assert "query" in result
        assert "databases" in result
        assert "n_hits_raw" in result
        assert "n_hits_filtered" in result
        assert "top_hit" in result
        assert "hits" in result
        assert "demo" in result
        assert "disclaimer" in result
        assert "chat_summary_lines" in result

    def test_tables_written(self, tmp_path):
        output_dir = tmp_path / "out"
        db_cache = tmp_path / "dbs"
        self._make_fake_db(db_cache, "pdb", "pdb")

        with patch("struct_predictor_foldseek.shutil.which", return_value="/usr/bin/foldseek"):
            with patch("struct_predictor_foldseek._foldseek_version", return_value="9.0"):
                with patch(
                    "struct_predictor_foldseek.subprocess.run",
                    side_effect=self._fake_run([_make_hit()]),
                ):
                    run_foldseek_search(
                        input_path=None,
                        output_dir=output_dir,
                        demo=True,
                        db_cache=db_cache,
                    )

        assert (output_dir / "tables" / "hits.tsv").exists()
        assert (output_dir / "tables" / "hits_filtered.tsv").exists()

    def test_reproducibility_files_written(self, tmp_path):
        output_dir = tmp_path / "out"
        db_cache = tmp_path / "dbs"
        self._make_fake_db(db_cache, "pdb", "pdb")

        with patch("struct_predictor_foldseek.shutil.which", return_value="/usr/bin/foldseek"):
            with patch("struct_predictor_foldseek._foldseek_version", return_value="9.0"):
                with patch(
                    "struct_predictor_foldseek.subprocess.run",
                    side_effect=self._fake_run([_make_hit()]),
                ):
                    run_foldseek_search(
                        input_path=None,
                        output_dir=output_dir,
                        demo=True,
                        db_cache=db_cache,
                    )

        repro = output_dir / "reproducibility"
        assert (repro / "commands.sh").exists()
        assert (repro / "environment.txt").exists()
        env_text = (repro / "environment.txt").read_text()
        assert "9.0" in env_text

    def test_foldseek_nonzero_exit_raises(self, tmp_path):
        output_dir = tmp_path / "out"
        db_cache = tmp_path / "dbs"
        self._make_fake_db(db_cache, "pdb", "pdb")

        def _fail_run(cmd, **kwargs):
            mock_proc = MagicMock()
            mock_proc.returncode = 1
            mock_proc.stderr = "Segmentation fault"
            mock_proc.stdout = ""
            return mock_proc

        with patch("struct_predictor_foldseek.shutil.which", return_value="/usr/bin/foldseek"):
            with patch("struct_predictor_foldseek._foldseek_version", return_value="9.0"):
                with patch("struct_predictor_foldseek.subprocess.run", side_effect=_fail_run):
                    with pytest.raises(RuntimeError, match="foldseek exited with code 1"):
                        run_foldseek_search(
                            input_path=None,
                            output_dir=output_dir,
                            demo=True,
                            db_cache=db_cache,
                        )

    def test_no_hits_after_filtering(self, tmp_path):
        """When all hits are below threshold, report still generates cleanly."""
        output_dir = tmp_path / "out"
        db_cache = tmp_path / "dbs"
        self._make_fake_db(db_cache, "pdb", "pdb")

        low_hits = [_make_hit(f"H{i}", tmscore=0.1) for i in range(5)]

        with patch("struct_predictor_foldseek.shutil.which", return_value="/usr/bin/foldseek"):
            with patch("struct_predictor_foldseek._foldseek_version", return_value="9.0"):
                with patch(
                    "struct_predictor_foldseek.subprocess.run",
                    side_effect=self._fake_run(low_hits),
                ):
                    result = run_foldseek_search(
                        input_path=None,
                        output_dir=output_dir,
                        demo=True,
                        db_cache=db_cache,
                        min_tmscore=0.3,
                    )

        assert result["n_hits_raw"] == 5
        assert result["n_hits_filtered"] == 0
        assert result["top_hit"] is None
        assert (output_dir / "report.md").exists()
        report = (output_dir / "report.md").read_text()
        assert "No hits passed" in report

    def test_max_hits_limits_result_json(self, tmp_path):
        output_dir = tmp_path / "out"
        db_cache = tmp_path / "dbs"
        self._make_fake_db(db_cache, "pdb", "pdb")

        many_hits = [_make_hit(f"H{i:04d}", tmscore=0.9) for i in range(100)]

        with patch("struct_predictor_foldseek.shutil.which", return_value="/usr/bin/foldseek"):
            with patch("struct_predictor_foldseek._foldseek_version", return_value="9.0"):
                with patch(
                    "struct_predictor_foldseek.subprocess.run",
                    side_effect=self._fake_run(many_hits),
                ):
                    result = run_foldseek_search(
                        input_path=None,
                        output_dir=output_dir,
                        demo=True,
                        db_cache=db_cache,
                        max_hits=10,
                    )

        assert len(result["hits"]) == 10

    def test_input_path_passed_through(self, tmp_path):
        """When using --input, the query file is passed to foldseek."""
        output_dir = tmp_path / "out"
        db_cache = tmp_path / "dbs"
        self._make_fake_db(db_cache, "pdb", "pdb")

        query_cif = tmp_path / "my_protein.cif"
        query_cif.write_text(DEMO_STRUCTURE.read_text())

        captured_cmds: list[list] = []

        def _capture_run(cmd, **kwargs):
            captured_cmds.append(cmd)
            results_tsv = Path(cmd[4])
            results_tsv.parent.mkdir(parents=True, exist_ok=True)
            results_tsv.write_text("")  # empty = no hits
            mock_proc = MagicMock()
            mock_proc.returncode = 0
            mock_proc.stdout = mock_proc.stderr = ""
            return mock_proc

        with patch("struct_predictor_foldseek.shutil.which", return_value="/usr/bin/foldseek"):
            with patch("struct_predictor_foldseek._foldseek_version", return_value="9.0"):
                with patch("struct_predictor_foldseek.subprocess.run", side_effect=_capture_run):
                    run_foldseek_search(
                        input_path=query_cif,
                        output_dir=output_dir,
                        demo=False,
                        db_cache=db_cache,
                    )

        assert len(captured_cmds) == 1
        cmd = captured_cmds[0]
        assert str(query_cif) in cmd

    def test_multi_database_merges_hits(self, tmp_path):
        """Results from multiple databases should be merged and ranked together."""
        output_dir = tmp_path / "out"
        db_cache = tmp_path / "dbs"
        self._make_fake_db(db_cache, "pdb", "pdb")
        self._make_fake_db(db_cache, "afdb", "afdb50")

        call_count = [0]
        hits_per_db = [
            [_make_hit("PDB_A", tmscore=0.80)],
            [_make_hit("AFDB_B", tmscore=0.90)],
        ]

        def _multi_run(cmd, **kwargs):
            results_tsv = Path(cmd[4])
            results_tsv.parent.mkdir(parents=True, exist_ok=True)
            hits = hits_per_db[call_count[0] % 2]
            call_count[0] += 1
            with open(results_tsv, "w", encoding="utf-8") as fh:
                for h in hits:
                    fh.write("\t".join(str(h.get(c, "")) for c in _FOLDSEEK_COLS) + "\n")
            mock_proc = MagicMock()
            mock_proc.returncode = 0
            mock_proc.stdout = mock_proc.stderr = ""
            return mock_proc

        with patch("struct_predictor_foldseek.shutil.which", return_value="/usr/bin/foldseek"):
            with patch("struct_predictor_foldseek._foldseek_version", return_value="9.0"):
                with patch("struct_predictor_foldseek.subprocess.run", side_effect=_multi_run):
                    result = run_foldseek_search(
                        input_path=None,
                        output_dir=output_dir,
                        databases="pdb,afdb",
                        demo=True,
                        db_cache=db_cache,
                    )

        assert result["n_hits_raw"] == 2
        # Top hit should be AFDB_B (higher TM-score)
        assert result["top_hit"]["target"] == "AFDB_B"


# ---------------------------------------------------------------------------
# TestFoldseekVersion
# ---------------------------------------------------------------------------


class TestFoldseekVersion:
    def test_returns_string_on_success(self):
        mock_proc = MagicMock()
        mock_proc.stdout = "foldseek 9.427df8a"
        mock_proc.stderr = ""
        with patch("struct_predictor_foldseek.subprocess.run", return_value=mock_proc):
            v = _foldseek_version("/usr/bin/foldseek")
        assert "9" in v

    def test_returns_unknown_on_exception(self):
        with patch("struct_predictor_foldseek.subprocess.run", side_effect=Exception("boom")):
            v = _foldseek_version("/usr/bin/foldseek")
        assert v == "unknown"
