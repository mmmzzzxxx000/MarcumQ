# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Shared reproducibility and battery/development timing safeguards."""
import csv
from datetime import datetime,timezone,timedelta
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import ctypes


def power_status():
    """Read-only power query; unknown never qualifies as publication power."""
    if os.name!="nt":return dict(ac_connected=None,source="unavailable")
    class Power(ctypes.Structure):
        _fields_=[("ACLineStatus",ctypes.c_ubyte),("BatteryFlag",ctypes.c_ubyte),("BatteryLifePercent",ctypes.c_ubyte),("SystemStatusFlag",ctypes.c_ubyte),("BatteryLifeTime",ctypes.c_uint32),("BatteryFullLifeTime",ctypes.c_uint32)]
    value=Power()
    if not ctypes.windll.kernel32.GetSystemPowerStatus(ctypes.byref(value)):
        return dict(ac_connected=None,source="GetSystemPowerStatus failed")
    return dict(ac_connected=True if value.ACLineStatus==1 else (False if value.ACLineStatus==0 else None),
                battery_percent=int(value.BatteryLifePercent),source="GetSystemPowerStatus")


def timing_policy(requested=None):
    requested=requested or os.environ.get("MARCUM_RUN_PURPOSE","development")
    if requested not in ("development","publication"):raise ValueError("unknown timing purpose")
    power=power_status();allowed=requested=="publication" and power["ac_connected"] is True
    return dict(requested_purpose=requested,timing_purpose="publication" if allowed else "development",
                publication_timing_allowed=allowed,power=power,
                power_context="Initial unplugged-session timings are development only. Publication mode requires an explicit new run and verified AC power.")

from marcum._paths import ROOT
sys.path.insert(0,str(ROOT/"src"))
sys.path.insert(0,str(ROOT/"validation"))


def write_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def write_csv(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    keys=[]
    for row in rows:
        for key in row:
            if key not in keys:keys.append(key)
    with path.open("w",newline="",encoding="utf-8-sig") as handle:
        writer=csv.DictWriter(handle,fieldnames=keys);writer.writeheader()
        for row in rows:
            writer.writerow({key:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list,tuple)) else v for key,v in row.items()})


def environment():
    from marcum import __version__
    versions={}
    for name in ("numpy","scipy","mpmath","matplotlib"):
        try:versions[name]=__import__(name).__version__
        except ImportError:versions[name]="not installed"
    return dict(project_version=__version__,python=sys.version,platform=platform.platform(),versions=versions,
                timestamp=datetime.now(timezone(timedelta(hours=8))).isoformat(),
                **timing_policy(),
                thread_settings={n:os.environ.get(n,"not explicitly set") for n in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS")})


def load_config(path=None):
    path=Path(path or ROOT/"configs/section5_1.json")
    cfg=json.loads(path.read_text(encoding="utf-8"))
    if cfg.get("timing_purpose")!="development" or cfg.get("publication_timing_allowed") is not False:
        raise ValueError("This session permits development timings only; use a fresh approved setup for publication measurements")
    if set(cfg["languages"])!={"python","fortran","cpp"}:
        raise ValueError("The supplied complete suite must include all three language environments")
    if set(cfg["methods"])!={"se","af","ni","se_rr"}:
        raise ValueError("The suite must contain the four independent methods")
    if set(c["region"] for c in cfg["cases"])!={"I","II","III","IV"}:
        raise ValueError("All four regions must be represented")
    cfg["config_sha256"]=hashlib.sha256(path.read_bytes()).hexdigest()
    return cfg
