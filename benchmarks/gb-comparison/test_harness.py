"""Unit tests for the two harness decisions that read the record.

Both read `results.json` and say what ran, so both can claim more than
the record holds. Run them with:

    python3 -m pytest benchmarks/gb-comparison/test_harness.py

They read no result file and start no benchmark.
"""

import driver
import report

INSTANCES = ["cyclic-4-q", "katsura-4-q"]
TOOLS = [("singular", "-"), ("msolve", "-"), ("sylvester", "rational")]


def cell(status, basis=None):
    return {"status": status, "basis": basis}


def test_a_tool_on_path_with_no_cell_is_not_listed():
    results = {"cyclic-4-q|singular|-": cell("OK")}
    assert driver.tools_with_a_rational_cell(results, INSTANCES, TOOLS) == ["singular"]


def test_a_tool_whose_every_cell_failed_is_not_listed():
    results = {
        "cyclic-4-q|singular|-": cell("OK"),
        "cyclic-4-q|msolve|-": cell("DNF"),
        "katsura-4-q|msolve|-": cell("ERROR"),
    }
    assert driver.tools_with_a_rational_cell(results, INSTANCES, TOOLS) == ["singular"]


def test_every_tool_with_one_finished_cell_is_listed():
    results = {
        "cyclic-4-q|singular|-": cell("OK"),
        "katsura-4-q|msolve|-": cell("OK"),
        "cyclic-4-q|sylvester|rational": cell("OK"),
    }
    assert driver.tools_with_a_rational_cell(results, INSTANCES, TOOLS) == [
        "singular",
        "msolve",
        "sylvester",
    ]


def test_one_full_output_is_no_comparison():
    cells = [("singular", ["x"]), ("msolve", None)]
    assert report.rational_agreement(cells) == ("NO COMPARISON", [])


def test_no_full_output_is_no_comparison():
    cells = [("singular", None), ("msolve", None)]
    assert report.rational_agreement(cells) == ("NO COMPARISON", [])


def test_two_equal_full_outputs_agree():
    cells = [("singular", ["x"]), ("msolve", ["x"]), ("groebner.jl", None)]
    assert report.rational_agreement(cells) == ("AGREE", [])


def test_two_different_full_outputs_mismatch():
    cells = [("singular", ["x"]), ("msolve", ["y"])]
    assert report.rational_agreement(cells) == ("MISMATCH", [("singular", "msolve")])


def test_a_custom_json_record_keeps_its_csv_next_to_it():
    assert report.results_csv_path("/tmp/smoke.json") == "/tmp/smoke.csv"


def test_an_explicit_csv_path_wins():
    assert (
        report.results_csv_path("/tmp/smoke.json", "/tmp/table.csv") == "/tmp/table.csv"
    )
