# Project authors: Zhong-Xuan Mao and Jing-Feng Tian.
# Paper: Efficient numerical computation of the generalized Marcum
# Q-function with applications (2026).
# Preprint: https://arxiv.org/abs/2611.29209
# Citation: CITATION.cff / CITATION.bib at the repository root.
# Licensing and upstream attribution: LICENSE.md / third_party/NOTICE.md.

"""Pinned, official Algorithm 939 source acquisition and byte verification."""
import hashlib
import json
from pathlib import Path
import urllib.request
import zipfile

from marcum._paths import ROOT
SOURCE_DIR = ROOT/"third_party"/"algorithm939"
URL = "https://calgo.acm.org/939.zip"
SHA256 = "dbf77461a9c02c603f8eb2bfcfcef6e3d1cfbb7d628b54b18f79c4a7e62af9e4"
ARCHIVE = SOURCE_DIR/"939.zip"
UPSTREAM = SOURCE_DIR/"upstream"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def verify_sources():
    if not ARCHIVE.is_file():
        raise RuntimeError("Missing official archive; run scripts/fetch_939.py explicitly.")
    if digest(ARCHIVE.read_bytes()) != SHA256:
        raise RuntimeError("Algorithm 939 archive SHA256 differs from the pinned release.")
    hashes = {}
    with zipfile.ZipFile(str(ARCHIVE)) as archive:
        for info in archive.infolist():
            if info.is_dir():
                continue
            path = UPSTREAM/info.filename
            expected = archive.read(info.filename)
            if not path.is_file() or path.read_bytes() != expected:
                raise RuntimeError("Upstream file missing or modified: "+info.filename)
            hashes[info.filename] = digest(expected)
    return hashes


def fetch():
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    if not ARCHIVE.is_file():
        with urllib.request.urlopen(URL, timeout=45) as response:
            data = response.read()
        if digest(data) != SHA256:
            raise RuntimeError("Downloaded release does not match pinned SHA256; no files extracted.")
        ARCHIVE.write_bytes(data)
    if digest(ARCHIVE.read_bytes()) != SHA256:
        raise RuntimeError("Existing archive differs from pinned SHA256; refusing to overwrite it.")
    base = UPSTREAM.resolve()
    with zipfile.ZipFile(str(ARCHIVE)) as archive:
        for info in archive.infolist():
            if info.is_dir():
                continue
            target = (base/info.filename).resolve()
            if base not in target.parents:
                raise RuntimeError("Unsafe archive entry: "+info.filename)
            data = archive.read(info.filename)
            if target.exists() and target.read_bytes() != data:
                raise RuntimeError("Refusing to replace modified upstream file: "+info.filename)
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists():
                target.write_bytes(data)
    hashes = verify_sources()
    record = dict(algorithm=939, authors=["Amparo Gil", "Javier Segura", "Nico M. Temme"],
                  citation="ACM TOMS 40(3), Article 20 (2014), DOI 10.1145/2591004",
                  download_url=URL, index_url="https://calgo.acm.org/",
                  archive_sha256=SHA256, upstream_files=hashes,
                  license_url="https://www.acm.org/publications/policies/software-copyright-notice",
                  license_note="Upstream is subject to the ACM Software Copyright and License Agreement; no relicensing is applied.")
    (SOURCE_DIR/"source.json").write_text(json.dumps(record, indent=2)+"\n", encoding="utf-8")
    return record
