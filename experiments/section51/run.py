# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Accuracy-matched full-plane benchmark using the unchanged four kernels."""
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import random
import statistics
import sys
import time

import mpmath as mp


OUT = Path(__file__).resolve().parent
from marcum._paths import ROOT
OUT = ROOT / "results/section51"
OUT.mkdir(parents=True, exist_ok=True)
sys.dont_write_bytecode = True
sys.path[:0] = [str(ROOT / "src/python"), str(ROOT / "experiments"), str(ROOT / "validation")]
from marcum.rebuilt import kernel
from fp_environment import load
from common import power_status

LANGUAGES = ("python", "cpp", "fortran")
METHODS = ("se", "af", "ni", "se_rr")
LIMITS = {"se": (0, 1000), "af": (1, 100), "ni": (0, 10), "se_rr": (1, 1000)}
GRID = [i / 2 for i in range(1, 41)]
WORKER_GUARD = None


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")
    temporary.replace(path)


def references():
    rows, provenance = [], {}
    for nu in (1, 10, 20):
        base = ROOT / 'data/section51/references' / ('nu%d' % nu)
        provenance[str(nu)] = {}
        for a in GRID:
            path = base / ('a_%03d.json' % (round(a*4)-1))
            column = json.loads(path.read_text(encoding='utf-8'))
            assert column['a'] == a and column.get('nu', 1) == nu
            assert float(column['max_relative_difference']) < 1e-25
            assert float(column['relative_omitted_tail_bound']) < 1e-90
            provenance[str(nu)][path.name] = sha(path)
            for b in GRID:
                rows.append((a, b, nu, column['q'][column['b'].index(b)]))
    return rows, provenance


def init_worker():
    global WORKER_GUARD
    WORKER_GUARD, _ = load()


def calibrate(job):
    language, method, a, b, nu, reference = job
    evaluate, _ = kernel(method, a, b, nu, language)
    first, last = LIMITS[method]
    selected, chosen_value, chosen_error, best, issue = None, None, None, mp.inf, ""
    with mp.workdps(100):
        target = mp.mpf(reference)
        for setting in range(first, last+1):
            try:
                WORKER_GUARD.million_set_binary64()
                value = evaluate(setting)
                if not math.isfinite(value):
                    continue
                error = abs(mp.mpf(value)-target)/target
                best = min(best, error)
                if error <= mp.mpf("1e-10"):
                    selected, chosen_value, chosen_error = setting, value, float(error)
                    break
            except (ArithmeticError, ValueError) as exc:
                issue = type(exc).__name__
    return dict(language=language, method=method, a=a, b=b, nu=nu, setting=selected,
                value=chosen_value, relative_error=chosen_error, reference=reference,
                best_error=float(best) if mp.isfinite(best) else None, issue=issue)


def measure(calibration, path):
    rows = []
    completed = set()
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            item = json.loads(line)
            rows.append(item)
            completed.add((item["language"], item["nu"], item["a"], item["b"], item["method"]))
    guard, _ = load()
    rng = random.Random(20261109)
    tasks = calibration[:]
    rng.shuffle(tasks)
    with path.open("a", encoding="utf-8", buffering=1) as handle:
        for index, row in enumerate(tasks):
            key = (row["language"], row["nu"], row["a"], row["b"], row["method"])
            if key in completed:
                continue
            result = dict(row)
            result.update(microseconds=None, samples=[], repeats=0)
            if row["setting"] is not None:
                evaluate, timer = kernel(row["method"], row["a"], row["b"], row["nu"], row["language"])
                guard.million_set_binary64()
                value = evaluate(row["setting"])
                with mp.workdps(100):
                    error = abs(mp.mpf(value)-mp.mpf(row["reference"]))/mp.mpf(row["reference"])
                    if not mp.isfinite(error) or error > mp.mpf("1e-10"):
                        raise ArithmeticError("Calibrated result failed the serial accuracy check")
                guard.million_set_binary64()
                pilot = timer(row["setting"], 10)
                repeats = max(10, min(30000, math.ceil(.0006/max(pilot, 1e-9))))
                for _ in range(5):
                    guard.million_set_binary64()
                    seconds = timer(row["setting"], repeats)
                    assert math.isfinite(seconds) and seconds > 0
                    result["samples"].append(seconds*1e6)
                result.update(microseconds=statistics.median(result["samples"]), repeats=repeats)
            handle.write(json.dumps(result, allow_nan=False)+"\n")
            rows.append(result)
            if (index+1) % 800 == 0:
                print("Timed", index+1, "/", len(tasks), flush=True)
    return rows


def main():
    before = power_status()
    assert before["ac_connected"], "AC power is required for publication timing."
    refs, provenance = references()
    jobs = [(language, method, a, b, nu, reference)
            for language in LANGUAGES for method in METHODS for a, b, nu, reference in refs]
    signature = dict(grid=GRID, nu=[1, 10, 20], relative_error="1e-10", limits=LIMITS,
                     raw_libraries={language: sha(ROOT / "build" / "rebuilt" / language / "rebuilt.dll")
                                    for language in ("cpp", "fortran")},
                     python_sources={name: sha(ROOT / "src/python" / "marcum" / name)
                                     for name in ("serr.py", "algorithm939.py")},
                     driver_sha256=sha(Path(__file__)), reference_provenance=provenance)
    write_json(OUT / "manifest.json", signature)
    calibration_path = OUT / "calibration.jsonl"
    done, calibration = set(), []
    if calibration_path.exists():
        for line in calibration_path.read_text(encoding="utf-8").splitlines():
            item = json.loads(line)
            calibration.append(item)
            done.add((item["language"], item["method"], item["a"], item["b"], item["nu"]))
    pending = [job for job in jobs if job[:5] not in done]
    print("Calibration remaining", len(pending), "/", len(jobs), flush=True)
    with calibration_path.open("a", encoding="utf-8", buffering=1) as handle:
        with ProcessPoolExecutor(max_workers=6, initializer=init_worker) as pool:
            for result in pool.map(calibrate, pending, chunksize=16):
                calibration.append(result)
                handle.write(json.dumps(result, allow_nan=False)+"\n")
                if len(calibration) % 400 == 0:
                    print("Calibrated", len(calibration), "/", len(jobs), flush=True)
    assert len(calibration) == len(jobs)
    print("Calibration complete. Starting serial kernel timing.", flush=True)
    rows = measure(calibration, OUT / "timings.jsonl")
    assert len(rows) == len(jobs)
    after = power_status()
    write_json(OUT / "status.json", dict(complete=True, evaluations=len(rows),
               qualified=sum(row["setting"] is not None for row in rows),
               power_before=before, power_after=after, publication_timing_allowed=bool(after["ac_connected"])))
    print("DONE", len(rows), "method/point/language combinations", flush=True)


if __name__ == "__main__":
    main()
