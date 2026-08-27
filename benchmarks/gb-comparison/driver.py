"""Benchmark driver. Runs every (instance, tool, config) cell sequentially.

Timing:
- Singular: internal rtimer (ms resolution) around std(I) only.
- M2: internal elapsedTiming around groebnerBasis only.
- msolve: median of repeated core_msolve calls in
  runner/msolve/msolve_inproc.c, not process wall clock. A separate
  untimed `msolve -g 2` call supplies the basis.
- Julia: @elapsed around groebner after a cyclic-3 warm-up in the same
  process.
- sylvester: Instant around the compute call only (parse excluded),
  measured inside the runner.

One thread except `sylvester-f4 threads8`. Singular and M2 use their
serial `std`/`groebnerBasis` paths.

core_msolve's `print_gb` path returns without cleanup. Measured directly
(5, 100, and 500 repetitions on cyclic-4 in one process): that leaks
close to a megabyte of heap per call. The runner forks one child per
repetition. See runner/msolve/msolve_inproc.c.

On success the driver runs twice more and takes the median of 3. msolve
takes the median of its in-process repetitions instead. A timeout (>120 s)
marks DNF and skips larger members of that family for the same tool and
config.
"""

import json
import os
import re
import statistics
import subprocess
import time

import canon

HERE = os.path.dirname(os.path.abspath(__file__))
INP = os.path.join(HERE, "inputs")
OUT = os.path.join(HERE, "outputs")
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
P = 1073741827
TIMEOUT = 120.0  # per computation
PROC_TIMEOUT = 150  # slack for startup/parse for external single-run procs
JULIA_PROC_TIMEOUT = 3 * 120 + 200

# Affinity string stored in the record. The caller pins the process with
# taskset; children inherit it.
CORES = os.environ.get("GBBENCH_CORES", "0-3,12-15")

FAMILIES = {
    "cyclic": [f"cyclic-{n}" for n in range(4, 9)],
    "katsura": [f"katsura-{n}" for n in range(4, 11)],
    "eco": [f"eco-{n}" for n in range(8, 12)],
    "noon": [f"noon-{n}" for n in range(3, 7)],
}

# Comma-separated instance allowlist, e.g. GBBENCH_INSTANCES=cyclic-4,katsura-5.
# Unset runs every instance.
_only = os.environ.get("GBBENCH_INSTANCES")
if _only:
    _allow = set(_only.split(","))
    FAMILIES = {
        fam: [n for n in members if n in _allow] for fam, members in FAMILIES.items()
    }
    FAMILIES = {fam: members for fam, members in FAMILIES.items() if members}

INSTANCES = [i for fam in FAMILIES.values() for i in fam]
NVARS = {}
for name in INSTANCES:
    with open(os.path.join(INP, name + ".syl")) as f:
        NVARS[name] = int(f.readline())

# (results label, config label, runner mode, thread count).
# A certified config is its own results row, not a variant of default.
SYLV_CONFIGS = [
    ("sylvester-f4", "default", "f4", 1),
    ("sylvester-f4", "threads8", "f4", 8),
    ("sylvester-classic", "default", "classic", 1),
    ("sylvester-classic", "certified", "certified", 1),
    ("sylvester-f4", "certified", "f4-certified", 1),
]


def parse_mono(s, nvars):
    return canon.parse_mono(s, nvars)


MEMMAX = os.environ.get("GBBENCH_MEMMAX", "16G")
RESULTS_PATH = os.environ.get("GBBENCH_RESULTS", os.path.join(HERE, "results.json"))

# Set when the last child was killed by the kernel rather than exiting.
LAST_OOM = [False]


def scoped(cmd):
    """Run a command in a memory-capped transient cgroup."""
    return [
        "systemd-run",
        "--user",
        "--scope",
        "-q",
        "-p",
        f"MemoryMax={MEMMAX}",
        "-p",
        "MemorySwapMax=0",
        "--",
    ] + cmd


def run_proc(cmd, timeout, extra_env=None):
    t0 = time.monotonic()
    LAST_OOM[0] = False
    env = None
    if extra_env:
        env = dict(os.environ)
        env.update(extra_env)
    try:
        cp = subprocess.run(
            scoped(cmd),
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
            check=False,
        )
        # -9/137 SIGKILL (cgroup OOM), -6/134 SIGABRT (Rust alloc failure).
        if cp.returncode in (-9, 137, -6, 134):
            LAST_OOM[0] = True
            return cp.returncode, cp.stdout, cp.stderr, time.monotonic() - t0, True
        return cp.returncode, cp.stdout, cp.stderr, time.monotonic() - t0, False
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        err = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or "")
        return -1, out, err, time.monotonic() - t0, True


def parse_marked_line(line, time_key, parsed):
    """Apply one recognized runner output line."""
    if line.startswith(time_key + " "):
        parsed["t"] = float(line.split()[1])
    elif line.startswith("SIZE "):
        parsed["size"] = int(line.split()[1])
    elif line.startswith("LM "):
        parsed["lms"].append(line[3:].strip())
    elif line.startswith("POLY "):
        body = line.removeprefix("POLY ").strip()
        if not body.isdigit():
            parsed["polys"].append(body)
    elif line.startswith("CERT_BYTES "):
        parsed["cert_bytes"] = int(line.split()[1])
    elif line.startswith("VERIFY_S "):
        parsed["verify_s"] = float(line.split()[1])
    elif line.startswith("PEAK_RSS_KB "):
        value = line.split()[1]
        parsed["peak_rss_kb"] = None if value == "null" else int(value)


def parse_marked_output(stdout, time_key):
    """Parse the common runner output fields."""
    parsed = {
        "t": None,
        "size": None,
        "lms": [],
        "polys": [],
        "cert_bytes": None,
        "verify_s": None,
        "peak_rss_kb": None,
    }
    for line in stdout.splitlines():
        parse_marked_line(line.strip(), time_key, parsed)
    return parsed


def parse_sylv_polys(stdout):
    """sylvester's "POLY <n>" blocks: n lines of "coeff e1 ... en" each.

    Returns a list of poly dicts (exponent tuple -> coefficient), already
    reduced mod P via `canon.poly_from_terms`.
    """
    lines = stdout.splitlines()
    polys = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("POLY "):
            n = int(line.split()[1])
            terms = []
            for j in range(1, n + 1):
                parts = lines[i + j].split()
                coeff = int(parts[0])
                exps = tuple(int(e) for e in parts[1:])
                terms.append((coeff, exps))
            polys.append(canon.poly_from_terms(terms, P))
            i += n + 1
        else:
            i += 1
    return polys


def text_polys_to_dicts(poly_exprs, nvars):
    return [canon.parse_full_poly(expr, nvars, P) for expr in poly_exprs]


def missing_basis(cmd, out, err):
    """An ERROR cell for a finisher whose full basis did not parse.

    report.py compares canonical bases. A size and leading-monomial list
    would drop out of that comparison without saying so.
    """
    detail = "no full basis in the output; " + (out + err)[-500:]
    return {"status": "ERROR", "detail": detail, "cmd": cmd}


def basis_json(poly_dicts):
    """canon.canon_basis(...) as plain nested lists, for JSON storage."""
    canonical = canon.canon_basis(poly_dicts, P)
    return [[[list(e), c] for e, c in poly] for poly in canonical]


def lms_from_dicts(poly_dicts):
    return [list(e) for e in canon.lms_of(poly_dicts)]


def canon_lms(lms):
    return sorted(lms)


def run_singular(name):
    cmd = ["Singular", "-q", os.path.join(INP, name + ".sing")]
    _rc, out, err, _wall, killed = run_proc(cmd, PROC_TIMEOUT)
    if killed:
        return {"status": "DNF", "cmd": cmd}
    parsed = parse_marked_output(out, "TIME_MS")
    if parsed["t"] is None or parsed["size"] is None:
        return {"status": "ERROR", "detail": (out + err)[-500:], "cmd": cmd}
    t = parsed["t"] / 1000.0
    if t > TIMEOUT:
        return {"status": "DNF", "cmd": cmd}
    nvars = NVARS[name]
    poly_dicts = text_polys_to_dicts(parsed["polys"], nvars)
    if not poly_dicts:
        return missing_basis(cmd, out, err)
    return {
        "status": "OK",
        "seconds": t,
        "size": parsed["size"],
        "lms": [parse_mono(s, nvars) for s in parsed["lms"]],
        "basis": basis_json(poly_dicts),
        "cmd": cmd,
    }


def run_m2(name):
    cmd = ["M2", "--script", os.path.join(INP, name + ".m2")]
    _rc, out, err, _wall, killed = run_proc(cmd, PROC_TIMEOUT)
    if killed:
        return {"status": "DNF", "cmd": cmd}
    parsed = parse_marked_output(out, "TIME_S")
    if parsed["t"] is None or parsed["size"] is None:
        return {"status": "ERROR", "detail": (out + err)[-500:], "cmd": cmd}
    if parsed["t"] > TIMEOUT:
        return {"status": "DNF", "cmd": cmd}
    nvars = NVARS[name]
    poly_dicts = text_polys_to_dicts(parsed["polys"], nvars)
    if not poly_dicts:
        return missing_basis(cmd, out, err)
    return {
        "status": "OK",
        "seconds": parsed["t"],
        "size": parsed["size"],
        "lms": [parse_mono(s, nvars) for s in parsed["lms"]],
        "basis": basis_json(poly_dicts),
        "cmd": cmd,
    }


def run_msolve_correctness(name):
    """One untimed `msolve -g 2` call. Supplies the basis. Not a timing cell."""
    outfile = os.path.join(OUT, name + ".msout")
    cmd = [
        "msolve",
        "-t",
        "1",
        "-g",
        "2",
        "-f",
        os.path.join(INP, name + ".ms"),
        "-o",
        outfile,
    ]
    rc, out, err, _wall, killed = run_proc(cmd, TIMEOUT + 10)
    if killed:
        return {"status": "DNF", "cmd": cmd}
    if rc != 0:
        return {"status": "ERROR", "detail": (out + err)[-500:], "cmd": cmd}
    with open(outfile) as f:
        text = f.read()
    body = "".join(line for line in text.splitlines() if not line.startswith("#"))
    body = body.strip().rstrip(":").strip()
    body = body.removeprefix("[")
    body = body.removesuffix("]")
    poly_exprs = [p for p in body.split(",") if p.strip()]
    nvars = NVARS[name]
    poly_dicts = text_polys_to_dicts(poly_exprs, nvars)
    if not poly_dicts:
        return {"status": "ERROR", "detail": "empty basis parse", "cmd": cmd}
    return {
        "status": "OK",
        "size": len(poly_dicts),
        "lms": lms_from_dicts(poly_dicts),
        "basis": basis_json(poly_dicts),
        "cmd": cmd,
    }


def run_msolve_timing(name):
    """Run the repeated in-process msolve timer."""
    binpath = os.path.join(HERE, "runner", "bin", "msolve-inproc")
    cmd = [binpath, os.path.join(INP, name + ".syl"), str(P), "1"]
    rc, out, err, _wall, killed = run_proc(cmd, TIMEOUT + 20)
    if killed:
        return {"status": "DNF", "cmd": cmd}
    if rc != 0:
        return {"status": "ERROR", "detail": (out + err)[-500:], "cmd": cmd}
    runs, count, peak_rss_kb = parse_msolve_runs(out)
    if not runs:
        return {"status": "ERROR", "detail": (out + err)[-500:], "cmd": cmd}
    if min(runs) > TIMEOUT:
        return {"status": "DNF", "cmd": cmd}
    return msolve_timing_result(cmd, runs, count, peak_rss_kb)


def parse_msolve_runs(stdout):
    """Parse the repeated msolve timer output."""
    runs = []
    count = None
    peak = None
    for line in stdout.splitlines():
        if line.startswith("RUN "):
            runs.append(float(line.split()[1]))
        elif line.startswith("COUNT "):
            count = int(line.split()[1])
        elif line.startswith("PEAK_RSS_KB "):
            value = line.split()[1]
            peak = None if value == "null" else int(value)
    return runs, count, peak


def msolve_timing_result(cmd, runs, count, peak_rss_kb):
    """Build one successful msolve timing cell."""
    return {
        "status": "OK",
        "seconds": statistics.median(runs),
        "runs": count if count is not None else len(runs),
        "spread_s": [min(runs), max(runs)],
        "peak_rss_kb": peak_rss_kb,
        "counters": counters(peak_rss_kb, None),
        "cmd": cmd,
    }


def run_julia_all(name):
    """One process: warmup plus 3 timed runs."""
    cmd = ["julia", "-t", "1", os.path.join(INP, name + ".jl")]
    _rc, out, err, _wall, killed = run_proc(cmd, JULIA_PROC_TIMEOUT)
    runs = [
        float(line.split()[1]) for line in out.splitlines() if line.startswith("RUN ")
    ]
    if julia_timed_out(killed, out, runs):
        return {"status": "DNF", "cmd": cmd}
    if not runs:
        return {"status": "ERROR", "detail": (out + err)[-500:], "cmd": cmd}
    parsed = parse_marked_output(out, "TIME_S")
    if parsed["size"] is None:
        return {"status": "ERROR", "detail": (out + err)[-500:], "cmd": cmd}
    nvars = NVARS[name]
    poly_dicts = text_polys_to_dicts(parsed["polys"], nvars)
    if not poly_dicts:
        return missing_basis(cmd, out, err)
    return {
        "status": "OK",
        "seconds": statistics.median(runs),
        "runs": len(runs),
        "size": parsed["size"],
        "lms": [parse_mono(s, nvars) for s in parsed["lms"]],
        "basis": basis_json(poly_dicts),
        "cmd": cmd,
    }


def julia_timed_out(killed, stdout, runs):
    """Return whether Julia crossed a process or computation limit."""
    return killed or "TIMEOUT" in stdout or (bool(runs) and max(runs) > TIMEOUT)


def run_sylvester(name, mode, threads):
    binpath = os.path.join(HERE, "runner", "bin", "sylv-runner")
    cmd = [binpath, os.path.join(INP, name + ".syl"), mode, str(threads), str(P)]
    _rc, out, err, _wall, killed = run_proc(
        cmd, PROC_TIMEOUT, extra_env={"RAYON_NUM_THREADS": str(threads)}
    )
    if killed or "STATUS TIMEOUT" in out:
        return {"status": "DNF", "cmd": cmd}
    parsed = parse_marked_output(out, "TIME_S")
    if parsed["t"] is None or parsed["size"] is None:
        return {"status": "ERROR", "detail": (out + err)[-500:], "cmd": cmd}
    if parsed["t"] > TIMEOUT:
        return {"status": "DNF", "cmd": cmd}
    poly_dicts = parse_sylv_polys(out)
    if not poly_dicts:
        return missing_basis(cmd, out, err)
    result = {
        "status": "OK",
        "seconds": parsed["t"],
        "size": parsed["size"],
        "lms": [[int(e) for e in s.split()] for s in parsed["lms"]],
        "basis": basis_json(poly_dicts),
        "peak_rss_kb": parsed["peak_rss_kb"],
        "counters": counters(
            parsed["peak_rss_kb"], parsed["size"], parse_counters(out)
        ),
        "cmd": cmd,
    }
    if mode.endswith("certified"):
        result["cert_bytes"] = parsed["cert_bytes"]
        result["verify_seconds"] = parsed["verify_s"]
    return result


# Field names of `sylvester::F4Counters`, plus `ComputeReport::threads_used`.
# A key stays null when the engine does not report it.
RESERVED_COUNTERS = [
    "pairs_generated",
    "pairs_discarded_product",
    "pairs_discarded_b",
    "pairs_discarded_m",
    "pairs_discarded_f",
    "batches",
    "matrix_rows",
    "matrix_columns",
    "matrix_nonzeros",
    "zero_rows",
    "new_pivots",
    "lane_restarts",
    "batch_retries",
    "basis_monomials",
    "threads_used",
]


def counters(peak_rss_kb, basis_size, engine=None):
    c = {"peak_rss_kb": peak_rss_kb, "basis_size": basis_size}
    for key in RESERVED_COUNTERS:
        c[key] = None
    for key, value in (engine or {}).items():
        c[key] = value
    return c


def parse_counters(stdout):
    """The runner's `COUNTERS <key>=<value> ...` line as a dict.

    Only numeric values go in, so `backend=` stays out. A cell with no
    such line yields an empty dict.
    """
    out = {}
    for line in stdout.splitlines():
        if not line.startswith("COUNTERS "):
            continue
        for field in line.split()[1:]:
            key, _, value = field.partition("=")
            if re.match(r"^\d+$", value):
                out[key] = int(value)
    return out


def tool_version(cmd, pattern, timeout=15):
    try:
        cp = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout, check=False
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    text = cp.stdout + cp.stderr
    m = re.search(pattern, text)
    return m.group(1) if m else None


def julia_pkg_version(pkg, timeout=60):
    script = (
        "import Pkg; d = Pkg.dependencies(); "
        f'for (_, p) in d; if p.name == "{pkg}"; println(p.version); end; end'
    )
    try:
        cp = subprocess.run(
            ["julia", "-e", script],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    out = cp.stdout.strip()
    return out or None


def sylvester_version():
    version = None
    cargo_toml = os.path.join(ROOT, "Cargo.toml")
    try:
        with open(cargo_toml) as f:
            for line in f:
                m = re.match(r'version\s*=\s*"([^"]+)"', line.strip())
                if m:
                    version = m.group(1)
                    break
    except OSError:
        pass
    commit = None
    try:
        cp = subprocess.run(
            ["git", "-C", ROOT, "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        if cp.returncode == 0:
            commit = cp.stdout.strip()
    except (OSError, subprocess.TimeoutExpired):
        pass
    if version and commit:
        return f"{version}+{commit}"
    return version or commit


def collect_versions():
    return {
        "singular": tool_version(["Singular", "--version"], r"version\s+([\d.]+)"),
        "msolve": tool_version(["msolve", "--help"], r"version\s+([\d.]+)"),
        "m2": tool_version(["M2", "--version"], r"([\d.]+)"),
        "julia": tool_version(["julia", "--version"], r"([\d.]+)"),
        "groebner.jl": julia_pkg_version("Groebner"),
        "sylvester": sylvester_version(),
    }


def median_cell(single_run, name):
    r1 = single_run(name)
    if r1["status"] != "OK":
        r1["runs"] = 1
        return r1
    times = [r1["seconds"]]
    for _ in range(2):
        r = single_run(name)
        if r["status"] != "OK":
            r["runs"] = len(times) + 1
            return r
        times.append(r["seconds"])
        if r["size"] != r1["size"] or canon_lms(r["lms"]) != canon_lms(r1["lms"]):
            return {
                "status": "ERROR",
                "detail": "nondeterministic output across runs",
                "cmd": r1.get("cmd"),
            }
        if r1.get("basis") is not None and r.get("basis") != r1.get("basis"):
            return {
                "status": "ERROR",
                "detail": "nondeterministic full basis across runs",
                "cmd": r1.get("cmd"),
            }
    r1["seconds"] = statistics.median(times)
    r1["runs"] = len(times)
    return r1


def msolve_cell(name):
    correctness = run_msolve_correctness(name)
    timing = run_msolve_timing(name)
    if timing["status"] != "OK":
        timing.setdefault("size", correctness.get("size"))
        timing.setdefault("lms", correctness.get("lms"))
        timing.setdefault("basis", correctness.get("basis"))
        return timing
    result = dict(timing)
    result["size"] = correctness.get("size")
    result["lms"] = correctness.get("lms")
    result["basis"] = correctness.get("basis")
    if "counters" in result:
        result["counters"]["basis_size"] = correctness.get("size")
    if correctness["status"] != "OK":
        result["status"] = correctness["status"]
        result["detail"] = correctness.get("detail")
    return result


def portable_command(command):
    """Remove machine-specific directory names from a stored command."""
    portable = []
    for argument in command:
        if not isinstance(argument, str) or not os.path.isabs(argument):
            portable.append(argument)
            continue
        try:
            inside_harness = os.path.commonpath((HERE, argument)) == HERE
        except ValueError:
            inside_harness = False
        if inside_harness:
            portable.append(os.path.relpath(argument, HERE))
        else:
            portable.append(os.path.basename(argument))
    return portable


def save_results(results):
    """Write the resumable record."""
    for cell in results.values():
        if isinstance(cell, dict) and isinstance(cell.get("cmd"), list):
            cell["cmd"] = portable_command(cell["cmd"])
    with open(RESULTS_PATH, "w") as file:
        json.dump(results, file, indent=1)


def load_results():
    """Load an existing record, or start an empty one."""
    os.makedirs(OUT, exist_ok=True)
    if not os.path.exists(RESULTS_PATH):
        return {}
    with open(RESULTS_PATH) as file:
        return json.load(file)


def update_metadata(results):
    """Record the tools and limits used by this run."""
    meta = results.setdefault("_meta", {})
    meta["versions"] = collect_versions()
    meta["cores"] = CORES
    meta["memmax"] = MEMMAX
    meta["p"] = P
    save_results(results)


def configured_tools():
    """Return every tool and configuration in record order."""
    tools = [
        ("singular", "-", lambda name: median_cell(run_singular, name)),
        ("msolve", "-", msolve_cell),
        ("m2", "-", lambda name: median_cell(run_m2, name)),
        ("groebner.jl", "-", run_julia_all),
    ]
    for tool, cfg, mode, threads in SYLV_CONFIGS:
        tools.append(
            (
                tool,
                cfg,
                lambda name, selected_mode=mode, selected_threads=threads: median_cell(
                    lambda instance: run_sylvester(
                        instance, selected_mode, selected_threads
                    ),
                    name,
                ),
            )
        )
    return tools


def run_cell(name, fn):
    """Run and normalize one new result cell."""
    print(f"RUNNING {name}", flush=True)
    started = time.monotonic()
    cell = fn(name.split("|", 1)[0])
    if cell.get("status") == "DNF" and LAST_OOM[0]:
        cell["status"] = "OOM"
        cell["detail"] = f"killed by the kernel at MemoryMax={MEMMAX}"
    cell["cell_wall_s"] = round(time.monotonic() - started, 3)
    cell["cores"] = CORES
    leading = cell.get("lms")
    if leading is not None:
        cell["lms"] = [list(term) for term in canon_lms(leading)]
    print(
        f"  -> {cell['status']} {cell.get('seconds', '')} size={cell.get('size', '')}",
        flush=True,
    )
    return cell


def run_family(results, members, tool, config, fn):
    """Run one tool through one monotone instance family."""
    stopped = False
    for name in members:
        key = f"{name}|{tool}|{config}"
        if key in results:
            stopped = stopped or results[key]["status"] in ("DNF", "OOM", "SKIP")
            continue
        if stopped:
            results[key] = {
                "status": "SKIP",
                "detail": "a smaller family member did not finish",
            }
        else:
            results[key] = run_cell(key, fn)
            stopped = results[key]["status"] in ("DNF", "OOM")
        save_results(results)


def main():
    """Run or resume the complete benchmark matrix."""
    results = load_results()
    update_metadata(results)
    for members in FAMILIES.values():
        for tool, config, fn in configured_tools():
            run_family(results, members, tool, config, fn)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
