#!/usr/bin/env python3
"""Resolve one read without prefix collisions or silent lane selection."""
import re
import sys
from pathlib import Path

def find_read(directory, sample, suffix):
    pattern = re.compile(re.escape(sample) + r"(?=[_.-])" + r".*" + suffix)
    matches = sorted(p for p in Path(directory).rglob("*")
                     if p.is_file() and pattern.fullmatch(p.name))
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one read for {sample!r} / {suffix!r}; "
                         f"found {len(matches)}: {matches}. Consolidate lanes explicitly.")
    return matches[0].resolve()

if __name__ == "__main__":
    try:
        print(find_read(*sys.argv[1:]))
    except (ValueError, re.error) as exc:
        sys.exit(f"ERROR: {exc}")
