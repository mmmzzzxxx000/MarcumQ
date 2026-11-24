# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Build the unmodified official source plus C adapter; run the official driver.

Examples: python -B scripts/build_939.py
          python -B scripts/build_939.py --debug
No network access is performed. No outputs are written into upstream/.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
from source939 import ROOT, UPSTREAM, SHA256, verify_sources


def build(fc=None, debug=False):
    compiler = shutil.which(fc or os.environ.get("FC", "gfortran"))
    if not compiler:
        raise RuntimeError("gfortran not found; supply --fc with its executable path")
    hashes = verify_sources()
    directory = ROOT/"build"/"algorithm939"/("debug" if debug else "release")
    directory.mkdir(parents=True, exist_ok=True)
    src = UPSTREAM/"939"/"Fortran90"/"Src"
    driver = UPSTREAM/"939"/"Fortran90"/"Drivers"/"testmarq.f90"
    adapter = ROOT/"src"/"fortran"/"algorithm939_c.f90"
    flags = (["-O0", "-g", "-fcheck=all", "-fbacktrace"] if debug else ["-O3"])
    flags += ["-fPIC", "-frecursive", "-fno-fast-math", "-I", str(src), "-J", str(directory)]
    commands, log = [], []

    def run(arguments, timeout=120):
        command = [str(v) for v in arguments]
        commands.append(command)
        completed = subprocess.run(command, cwd=str(directory), stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, universal_newlines=True, timeout=timeout)
        log.append("COMMAND: "+repr(command)+"\n"+completed.stdout)
        (directory/"build.log").write_text("\n".join(log), encoding="utf-8")
        if completed.returncode:
            raise RuntimeError("Command failed; inspect "+str(directory/"build.log"))
        return completed.stdout

    version = run([compiler, "--version"]).splitlines()[0]
    objects = []
    for source in (src/"cons.f90", src/"incgam.f90", src/"marcumq.f90", adapter):
        obj = directory/(source.stem+".o")
        run([compiler]+flags+["-c", source, "-o", obj])
        objects.append(obj)
    suffix = ".dll" if os.name == "nt" else (".dylib" if sys.platform == "darwin" else ".so")
    library = directory/("algorithm939"+suffix)
    link = [compiler, "-shared"]+objects+["-o", library]
    if os.name == "nt":
        link += ["-static-libgfortran", "-static-libgcc"]
    run(link)
    executable = directory/("testmarq.exe" if os.name == "nt" else "testmarq")
    run([compiler]+flags+objects[:3]+[driver, "-o", executable])
    run([executable], timeout=120)
    report = directory/"testQP.dat"
    text = report.read_text(encoding="utf-8")
    numeric = re.findall(r"^\s*\d+\s+[-+.\d]+D[-+]\d+\s+[-+.\d]+D[-+]\d+\s+([-+.\d]+D[-+]\d+)\s+([-+.\d]+D[-+]\d+)", text, re.M|re.I)
    errors = [float(value.replace("D", "e").replace("d", "e")) for row in numeric for value in row]
    match = re.search(r"Maximum value of the recurrence check\s*=\s*([-+.\dEeDd]+)", text)
    recurrence = float(match.group(1).replace("D", "e")) if match else float("inf")
    if len(numeric) != 18 or not errors or max(errors) > 2e-11 or recurrence > 1e-10:
        raise RuntimeError("Official driver regression did not meet its published-precision checks: "+str(report))
    manifest = dict(schema_version=1, timestamp_utc=datetime.now(timezone.utc).isoformat(),
                    compiler=compiler, compiler_version=version, platform=platform.platform(),
                    python_pointer_bits=64 if sys.maxsize > 2**32 else 32,
                    profile="debug" if debug else "release", source_archive_sha256=SHA256,
                    upstream_files=hashes, upstream_unmodified=True,
                    adapter_sha256=hashlib.sha256(adapter.read_bytes()).hexdigest(),
                    library=library.name, library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
                    runtime_search_directory=str(Path(compiler).parent), commands=commands,
                    official_driver=dict(rows=18, max_reference_relative_error=max(errors),
                                         max_recurrence_residual=recurrence, passed=True))
    (directory/"build_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(dict(library=str(library), official_driver=manifest["official_driver"]), indent=2))
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fc")
    parser.add_argument("--debug", action="store_true")
    arguments = parser.parse_args()
    build(arguments.fc, arguments.debug)
