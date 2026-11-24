# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""ctypes binding to the unmodified official Algorithm 939 + thin C adapter.

No rtol is passed to Fortran. The upstream accuracy aim is used ONLY as a
conservative operational compatibility check, never as a certified error bound.
"""
import ctypes
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import threading
def validate_parameters(a, b, nu):
    values = []
    for name, value in (("a", a), ("b", b), ("nu", nu)):
        if isinstance(value, (bool, str, bytes)):
            raise ValueError(name + " must be a real scalar")
        try:
            value = float(value)
        except (TypeError, ValueError, OverflowError):
            raise ValueError(name + " must be a finite real scalar")
        if not math.isfinite(value):
            raise ValueError(name + " must be finite")
        if value < 0 or (name == "nu" and value == 0):
            raise ValueError("require a >= 0, b >= 0, nu > 0")
        values.append(value)
    return tuple(values)



UPSTREAM_RELATIVE_AIM = 1e-11  # upstream README; fixed, not a tunable tolerance
_CALL_LOCK = threading.RLock()


class Algorithm939Unavailable(RuntimeError):
    pass


def domain_issue(a, b, nu):
    x, y = (0.5*a)*a, (0.5*b)*b
    if not (math.isfinite(x) and math.isfinite(y)):
        return "a^2/2 or b^2/2 is outside the original Algorithm 939 domain"
    if not (1 <= nu <= 10000 and 0 <= x <= 10000 and 0 <= y <= 10000):
        return "Original Algorithm 939 requires 1<=nu<=10000 and 0<=a^2/2,b^2/2<=10000"
    if b > 0 and y == 0:
        return "Positive b^2/2 underflowed during the original-parameter conversion"
    return None


def default_library_path():
    override = os.environ.get("MARCUM939_LIBRARY")
    if override:
        return Path(override).expanduser().resolve()
    from .._paths import ROOT as root
    suffix = ".dll" if os.name == "nt" else (".dylib" if sys.platform == "darwin" else ".so")
    return root/"build"/"algorithm939"/"release"/("algorithm939"+suffix)


class _Library:
    def __init__(self, path):
        if not path.is_file():
            raise Algorithm939Unavailable("Native library is missing: {}. Run scripts/build_939.py explicitly.".format(path))
        self.manifest = {}
        self.directory_handles = []
        metadata = path.with_name("build_manifest.json")
        if metadata.exists():
            try:
                self.manifest = json.loads(metadata.read_text(encoding="utf-8"))
            except (OSError, ValueError) as exc:
                raise Algorithm939Unavailable("Invalid build manifest: "+str(exc))
            actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
            if actual_hash != self.manifest.get("library_sha256"):
                raise Algorithm939Unavailable("Native library does not match its build manifest")
        if os.name == "nt" and hasattr(os, "add_dll_directory"):
            for folder in (str(path.parent), self.manifest.get("runtime_search_directory")):
                if folder and Path(folder).is_dir():
                    self.directory_handles.append(os.add_dll_directory(folder))
        try:
            self.dll = ctypes.CDLL(str(path))
            version = self.dll.algorithm939_abi_version_c
            version.argtypes = []
            version.restype = ctypes.c_int
            if version() != 1:
                raise Algorithm939Unavailable("Unsupported Algorithm 939 C adapter ABI")
            self.function = self.dll.algorithm939_eval_c
            self.function.argtypes = [ctypes.c_double]*3 + [ctypes.POINTER(ctypes.c_double)]*2 + [ctypes.POINTER(ctypes.c_int)]
            self.function.restype = None
        except (OSError, AttributeError) as exc:
            raise Algorithm939Unavailable("Cannot load Algorithm 939: "+str(exc))
        self.path = str(path)
        self.build_id = self.manifest.get("library_sha256")


@lru_cache(maxsize=4)
def _load(path_string):
    return _Library(Path(path_string))


@dataclass(frozen=True)
class NativePair:
    p: float
    q: float
    ierr: int
    library: str
    build_id: str


def native_pair(a, b, nu):
    """Raw fixed-accuracy baseline call; caller must inspect ierr.

    All domain validation occurs before any native call. This entry is for
    comparisons against the released software, not a high-precision reference.
    """
    a, b, nu = validate_parameters(a, b, nu)
    issue = domain_issue(a, b, nu)
    if issue:
        raise ValueError(issue)
    library = _load(str(default_library_path()))
    p, q, ierr = ctypes.c_double(), ctypes.c_double(), ctypes.c_int()
    with _CALL_LOCK:
        library.function(nu, (0.5*a)*a, (0.5*b)*b,
                         ctypes.byref(p), ctypes.byref(q), ctypes.byref(ierr))
    return NativePair(p.value, q.value, ierr.value, library.path, library.build_id)

