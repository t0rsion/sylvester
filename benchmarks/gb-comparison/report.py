"""Build the CSV and Markdown fragments from one benchmark record."""

import csv
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_PATH = os.environ.get("GBBENCH_RESULTS", os.path.join(HERE, "results.json"))
RESULTS_CSV = os.environ.get("GBBENCH_RESULTS_CSV", os.path.join(HERE, "results.csv"))

FAMILIES = {
    "cyclic": [f"cyclic-{n}" for n in range(4, 9)],
    "katsura": [f"katsura-{n}" for n in range(4, 11)],
    "eco": [f"eco-{n}" for n in range(8, 12)],
    "noon": [f"noon-{n}" for n in range(3, 7)],
}
INSTANCES = [instance for family in FAMILIES.values() for instance in family]

TOOLS = [
    ("singular", "-"),
    ("msolve", "-"),
    ("m2", "-"),
    ("groebner.jl", "-"),
    ("sylvester-f4", "default"),
    ("sylvester-f4", "threads8"),
    ("sylvester-classic", "default"),
    ("sylvester-classic", "certified"),
    ("sylvester-f4", "certified"),
]
SYLVESTER = [tool for tool in TOOLS if tool[0].startswith("sylvester")]
SYLVESTER_RAW = [tool for tool in SYLVESTER if tool[1] == "default"]
SYLVESTER_RAW_ALL = [tool for tool in SYLVESTER if tool[1] != "certified"]
EXTERNAL = [tool for tool in TOOLS if not tool[0].startswith("sylvester")]

GATE_CELLS = ["cyclic-7", "katsura-9", "eco-9", "noon-6"]
GATE_GEOMEAN_MAX = 5.0
GATE_WORST_MAX = 10.0


def result(results, name, tool, config):
    """Return one cell, or `None` when the record has no cell."""
    return results.get(f"{name}|{tool}|{config}")


def fmt_time(cell):
    """Format one timing cell for a Markdown table."""
    if cell is None:
        return "?"
    if cell["status"] != "OK":
        return cell["status"]
    seconds = cell["seconds"]
    if seconds < 0.001:
        return "<0.001"
    return f"{seconds:.3f}" if seconds < 10 else f"{seconds:.1f}"


def fmt_optional(value, pattern="{:.6f}"):
    """Format an optional number as an empty CSV cell when absent."""
    return pattern.format(value) if value is not None else ""


def markdown_table(header, rows):
    """Print a Markdown table."""
    print("| " + " | ".join(header) + " |")
    print("|" + "---|" * len(header))
    for row in rows:
        print("| " + " | ".join(map(str, row)) + " |")


def csv_row(name, tool, config, cell):
    """Flatten one result cell."""
    spread = cell.get("spread_s") or [None, None]
    seconds = fmt_optional(cell.get("seconds")) if cell["status"] == "OK" else ""
    return [
        name,
        tool,
        config,
        seconds,
        cell.get("runs", ""),
        fmt_optional(spread[0]),
        fmt_optional(spread[1]),
        cell.get("size", ""),
        cell.get("cert_bytes", ""),
        fmt_optional(cell.get("verify_seconds")),
        cell.get("peak_rss_kb", ""),
        cell["status"],
    ]


def write_csv(results):
    """Write the scalar form of the record."""
    header = [
        "instance",
        "tool",
        "config",
        "median_seconds",
        "runs",
        "spread_min_s",
        "spread_max_s",
        "basis_size",
        "cert_bytes",
        "verify_seconds",
        "peak_rss_kb",
        "status",
    ]
    with open(RESULTS_CSV, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(header)
        for name in INSTANCES:
            for tool, config in TOOLS:
                cell = result(results, name, tool, config)
                if cell is not None:
                    writer.writerow(csv_row(name, tool, config, cell))


def finishers(results, name):
    """Return the successful cells of one instance."""
    rows = []
    for tool, config in TOOLS:
        cell = result(results, name, tool, config)
        if cell and cell["status"] == "OK":
            label = f"{tool}[{config}]" if config != "-" else tool
            rows.append((label, cell["size"], cell.get("lms"), cell.get("basis")))
    return rows


def weak_agreement(name, rows, mismatches):
    """Check basis size and leading monomials."""
    reference, ref_size, ref_lms, _ = rows[0]
    agreed = True
    for label, size, leading, _ in rows[1:]:
        same_leads = sorted(map(tuple, leading or [])) == sorted(
            map(tuple, ref_lms or [])
        )
        if size != ref_size or not same_leads:
            mismatches.append((name, reference, label, ref_size, size))
            agreed = False
    return agreed


def full_agreement(name, rows, mismatches, missing):
    """Check coefficients and every monomial of the basis."""
    comparable = [(label, basis) for label, _, _, basis in rows if basis is not None]
    absent = [label for label, _, _, basis in rows if basis is None]
    if absent:
        missing.append((name, absent))
    if len(comparable) < 2:
        return "NO FULL OUTPUT"
    reference, ref_basis = comparable[0]
    agreed = True
    for label, basis in comparable[1:]:
        if basis != ref_basis:
            mismatches.append((name, reference, label))
            agreed = False
    return "AGREE" if agreed else "MISMATCH"


def print_mismatches(title, rows):
    """Print one mismatch group when nonempty."""
    if rows:
        print(f"\n{title}:")
        for row in rows:
            print(row)


def report_correctness(results):
    """Check and print complete-basis agreement."""
    print("== CORRECTNESS CROSS-CHECK ==")
    print("The complete canonical basis is the primary check.")
    weak_mismatches = []
    full_mismatches = []
    missing = []
    table = []
    for name in INSTANCES:
        rows = finishers(results, name)
        if len(rows) < 2:
            print(f"{name}: fewer than 2 finishers")
            table.append((name, len(rows), "-", "-"))
            continue
        weak = weak_agreement(name, rows, weak_mismatches)
        full = full_agreement(name, rows, full_mismatches, missing)
        weak_text = "AGREE" if weak else "MISMATCH"
        print(f"{name}: {len(rows)} finishers, size/LM {weak_text}, full basis {full}")
        table.append((name, len(rows), weak_text, full))
    print_mismatches("SIZE/LM MISMATCH DETAILS", weak_mismatches)
    print_mismatches("FULL BASIS MISMATCH DETAILS", full_mismatches)
    print_mismatches("FINISHERS WITH NO FULL BASIS", missing)
    if not weak_mismatches and not full_mismatches:
        print("\nAll compared complete bases agree.")
    print("\n== CORRECTNESS TABLE (markdown) ==")
    markdown_table(
        ["instance", "finishers", "size/LM agreement", "full basis agreement"],
        table,
    )


def best_raw(results, name):
    """Return the fastest one-thread raw sylvester cell."""
    best = None
    label = None
    for tool, config in SYLVESTER_RAW:
        cell = result(results, name, tool, config)
        if cell and cell["status"] == "OK" and (best is None or cell["seconds"] < best):
            best = cell["seconds"]
            label = tool.split("-", 1)[1] + "/" + config
    return label, best


def report_main_table(results):
    """Print external and fastest raw timings."""
    rows = []
    for name in INSTANCES:
        row = [name]
        row.extend(fmt_time(result(results, name, *tool)) for tool in EXTERNAL)
        label, seconds = best_raw(results, name)
        cell = (
            "DNF"
            if seconds is None
            else f"{fmt_time({'status': 'OK', 'seconds': seconds})} ({label})"
        )
        rows.append(row + [cell])
    print("\n== MAIN TABLE ==")
    markdown_table(
        ["instance", "Singular", "msolve", "M2", "Groebner.jl", "sylvester"],
        rows,
    )


def report_raw_table(results):
    """Print all raw sylvester configurations."""
    rows = []
    for name in INSTANCES:
        times = [fmt_time(result(results, name, *tool)) for tool in SYLVESTER_RAW_ALL]
        rows.append([name, *times])
    print("\n== SYLVESTER APPENDIX: RAW CONFIGS ==")
    markdown_table(["instance", "f4/default", "f4/threads8", "classic/default"], rows)


def certified_row(results, name, backend, tool):
    """Build one certified-path table row."""
    raw = result(results, name, tool, "default")
    certified = result(results, name, tool, "certified")
    complete = (
        raw and certified and raw["status"] == "OK" and certified["status"] == "OK"
    )
    split = f"{certified['seconds'] - raw['seconds']:.6f}" if complete else "?"
    return [
        name,
        backend,
        fmt_time(raw),
        fmt_time(certified),
        split,
        certified.get("cert_bytes", "") if certified else "",
        fmt_optional(certified.get("verify_seconds")) if certified else "",
        raw.get("peak_rss_kb", "") if raw else "",
        certified.get("peak_rss_kb", "") if certified else "",
    ]


def report_certified_table(results):
    """Print certificate sizes, times, and peak memory."""
    rows = []
    for name in INSTANCES:
        rows.append(certified_row(results, name, "classic", "sylvester-classic"))
        rows.append(certified_row(results, name, "f4", "sylvester-f4"))
    print("\n== SYLVESTER APPENDIX: CERTIFIED ==")
    print("split = certified total - raw. verify is a later standalone check.\n")
    markdown_table(
        [
            "instance",
            "backend",
            "raw (s)",
            "certified total (s)",
            "split (s)",
            "cert bytes",
            "verify (s)",
            "peak RSS raw (KB)",
            "peak RSS certified (KB)",
        ],
        rows,
    )


def report_basis_sizes(results):
    """Print basis sizes grouped by finishing tool."""
    print("\n== BASIS SIZES ==")
    for name in INSTANCES:
        sizes = {}
        for tool, config in TOOLS:
            cell = result(results, name, tool, config)
            if cell and cell["status"] == "OK":
                sizes.setdefault(cell["size"], []).append(tool)
        print(name, sizes)


def report_versions(meta):
    """Print the recorded tool and machine metadata."""
    print("\n== ENGINE VERSIONS ==")
    versions = meta.get("versions", {})
    engines = ("singular", "msolve", "m2", "julia", "groebner.jl", "sylvester")
    for engine in engines:
        print(f"{engine}: {versions.get(engine) or '?'}")
    print(f"cores: {meta.get('cores', '?')}")
    print(f"memory cap: {meta.get('memmax', '?')}")
    print(f"field: F_{meta.get('p', '?')}")


def gate_row(results, name):
    """Build one F4-to-msolve gate row."""
    msolve = result(results, name, "msolve", "-")
    f4 = result(results, name, "sylvester-f4", "default")
    seconds = f4.get("seconds") if f4 and f4["status"] == "OK" else None
    label = "f4/default" if seconds is not None else None
    if not msolve or msolve["status"] != "OK" or seconds is None:
        msolve_seconds = msolve.get("seconds") if msolve else None
        return [
            name,
            label or "-",
            fmt_optional(seconds),
            fmt_optional(msolve_seconds),
            "not finished",
        ], None
    ratio = seconds / msolve["seconds"]
    return [
        name,
        label,
        fmt_optional(seconds),
        fmt_optional(msolve["seconds"]),
        f"{ratio:.2f}x",
    ], ratio


def report_gate(results):
    """Compute and print the release performance gate."""
    rows = []
    ratios = []
    for name in GATE_CELLS:
        row, ratio = gate_row(results, name)
        rows.append(row)
        if ratio is not None:
            ratios.append(ratio)
    print("\n== GATE (docs/v0.2-design.md section 14): sylvester F4 vs msolve ==")
    markdown_table(
        ["cell", "sylvester raw", "sylvester (s)", "msolve (s)", "ratio"], rows
    )
    if len(ratios) != len(GATE_CELLS):
        print(f"\nnot finished: {len(GATE_CELLS) - len(ratios)} gate cells")
        return
    geomean = math.exp(sum(math.log(ratio) for ratio in ratios) / len(ratios))
    worst = max(ratios)
    passed = geomean <= GATE_GEOMEAN_MAX and worst <= GATE_WORST_MAX
    print(f"\ngeometric mean: {geomean:.2f}x (gate: at most {GATE_GEOMEAN_MAX}x)")
    print(f"worst cell: {worst:.2f}x (gate: at most {GATE_WORST_MAX}x)")
    print(f"gate: {'PASS' if passed else 'FAIL'}")


def main():
    """Generate every report view."""
    with open(RESULTS_PATH) as file:
        results = json.load(file)
    write_csv(results)
    report_correctness(results)
    report_main_table(results)
    report_raw_table(results)
    report_certified_table(results)
    report_basis_sizes(results)
    report_versions(results.get("_meta", {}))
    report_gate(results)


if __name__ == "__main__":
    main()
