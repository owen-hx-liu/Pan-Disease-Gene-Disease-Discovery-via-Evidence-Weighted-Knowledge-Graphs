#!/usr/bin/env python3
"""Verify data artifacts against the frozen ``ARTIFACT_HASHES.txt`` manifest.

The manifest (repo root) freezes the sha256, byte count and line count of every
input/intermediate artifact so that every machine computes on identical bytes.
It is consumed by ``run_all.ps1`` / ``run_all.sh`` and can be run directly.

Manifest format, one artifact per non-comment line::

    <sha256>  <bytes>  <lines>  <relative/path>

Comment (``#``) and blank lines are ignored. Paths are repo-relative; this
script infers the repo root from its own location, so it runs from anywhere.

Modes
-----
``--require PATH [PATH ...]``
    Every listed path MUST appear in the manifest, exist on disk, and match.
    A missing file or a mismatch exits non-zero. Use this to gate a run on its
    *inputs* (edges_clean_integrated.csv + the Hetionet files) before compute.

``--check-existing`` (default when no mode is given)
    Verify every manifest entry that exists on disk. A present-but-mismatched
    file fails (exit 1); an absent file is reported ``SKIP`` and is not a
    failure. Use this *after* a run to confirm reproduced splits are
    byte-identical to the frozen ones.

Exit status: ``0`` all good, ``1`` any mismatch / required-missing / manifest
problem.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
MANIFEST = os.path.join(REPO_ROOT, "ARTIFACT_HASHES.txt")

_CHUNK = 1 << 20  # 1 MiB


def _load_manifest() -> "list[tuple[str, int, int, str]]":
    """Return [(sha256, bytes, lines, relpath), ...] from ARTIFACT_HASHES.txt."""
    if not os.path.exists(MANIFEST):
        sys.exit(f"[verify_hashes] manifest not found: {MANIFEST}")
    entries = []
    with open(MANIFEST, encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                sys.exit(f"[verify_hashes] malformed manifest line {lineno}: {raw!r}")
            sha, nbytes, nlines, path = parts[0], parts[1], parts[2], " ".join(parts[3:])
            try:
                entries.append((sha.lower(), int(nbytes), int(nlines), path))
            except ValueError:
                sys.exit(f"[verify_hashes] non-integer size/line count on line {lineno}: {raw!r}")
    return entries


def _sha256(path: str) -> "tuple[str, int]":
    """Streaming sha256 + byte count of a file."""
    h = hashlib.sha256()
    nbytes = 0
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(_CHUNK)
            if not chunk:
                break
            nbytes += len(chunk)
            h.update(chunk)
    return h.hexdigest(), nbytes


def _check_one(entry: "tuple[str, int, int, str]") -> "tuple[str, str]":
    """Return (status, detail) for one manifest entry. status in OK/FAIL/MISS."""
    sha_expect, bytes_expect, _lines_expect, relpath = entry
    abspath = os.path.join(REPO_ROOT, relpath)
    if not os.path.exists(abspath):
        return "MISS", f"{relpath}  (not on disk)"
    size = os.path.getsize(abspath)
    if size != bytes_expect:  # fast fail before hashing large files
        return "FAIL", f"{relpath}  size {size:,} != expected {bytes_expect:,}"
    sha_got, _ = _sha256(abspath)
    if sha_got.lower() != sha_expect:
        return "FAIL", f"{relpath}  sha256 {sha_got[:12]}... != expected {sha_expect[:12]}..."
    return "OK", f"{relpath}  ({size:,} bytes)"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--require", nargs="+", metavar="PATH",
                      help="these repo-relative paths must be listed, exist, and match")
    mode.add_argument("--check-existing", action="store_true",
                      help="verify every manifest entry present on disk (absent => SKIP)")
    args = ap.parse_args()

    entries = _load_manifest()
    by_path = {e[3]: e for e in entries}

    if args.require:
        targets = []
        for p in args.require:
            rel = os.path.relpath(os.path.abspath(p), REPO_ROOT).replace(os.sep, "/")
            if rel not in by_path:
                print(f"  FAIL  {rel}  (not listed in ARTIFACT_HASHES.txt)")
                return 1
            targets.append(by_path[rel])
        print(f"Verifying {len(targets)} required input artifact(s) against ARTIFACT_HASHES.txt:")
        failed = 0
        for entry in targets:
            status, detail = _check_one(entry)
            print(f"  {status:<4}  {detail}")
            if status != "OK":  # MISS is a failure in --require mode
                failed += 1
        if failed:
            print(f"\n[verify_hashes] {failed} required artifact(s) missing or mismatched -- ABORT.")
            return 1
        print("[verify_hashes] all required inputs verified.")
        return 0

    # default / --check-existing
    print("Verifying every ARTIFACT_HASHES.txt entry present on disk:")
    failed = skipped = ok = 0
    for entry in entries:
        status, detail = _check_one(entry)
        print(f"  {status:<4}  {detail}")
        if status == "OK":
            ok += 1
        elif status == "MISS":
            skipped += 1
        else:
            failed += 1
    print(f"\n[verify_hashes] {ok} OK, {skipped} not-yet-built (SKIP), {failed} MISMATCH.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
