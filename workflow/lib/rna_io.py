#!/usr/bin/env python3
"""Validate five-column RNA manifests and extract featureCounts matrices."""
import csv
import re
import sys
from pathlib import Path

def sample_from_r1(path):
    name = Path(path).name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    for suffix in ("_R1", "_1", "R1", "1"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name

def validate_manifest(path):
    references, samples, rows = {}, set(), 0
    with open(path, newline="") as handle:
        for line, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not row or not any(row) or row[0].lstrip().startswith("#"):
                continue
            if len(row) != 5 or not all(row):
                raise ValueError(f"Line {line}: expected five nonempty tab-separated columns")
            species, r1, r2, genome, gtf = row
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", species):
                raise ValueError(f"Line {line}: invalid species identifier {species!r}")
            for file in row[1:]:
                if not Path(file).is_file() or Path(file).stat().st_size == 0:
                    raise ValueError(f"Line {line}: missing/empty input {file}")
            pair = (str(Path(genome).resolve()), str(Path(gtf).resolve()))
            if references.setdefault(species, pair) != pair:
                raise ValueError(f"Line {line}: inconsistent genome/GTF for {species}")
            sample = sample_from_r1(r1)
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
                raise ValueError(f"Line {line}: invalid sample identifier {sample!r}")
            if (species, sample) in samples:
                raise ValueError(f"Line {line}: duplicate derived sample {species}/{sample}")
            samples.add((species, sample))
            if (r1.endswith(".gz") != r2.endswith(".gz")) or Path(r1).resolve() == Path(r2).resolve():
                raise ValueError(f"Line {line}: mates must be distinct with matching compression")
            rows += 1
    if not rows:
        raise ValueError("Manifest contains no data rows")

def extract_counts(source, destination):
    with open(source) as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t")
                if row and not row[0].startswith("#")]
    if not rows or rows[0][:6] != ["Geneid", "Chr", "Start", "End", "Strand", "Length"] or len(rows[0]) < 7:
        raise ValueError("Missing featureCounts header or sample columns")
    names = [Path(value).name.removesuffix(".sorted.bam") for value in rows[0][6:]]
    if len(set(names)) != len(names):
        raise ValueError("Duplicate sample names in featureCounts header")
    if any(len(row) != len(rows[0]) for row in rows[1:]):
        raise ValueError("Inconsistent featureCounts column count")
    with open(destination, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Geneid", *names])
        writer.writerows([row[0], *row[6:]] for row in rows[1:])

if __name__ == "__main__":
    try:
        if sys.argv[1] == "validate":
            validate_manifest(sys.argv[2])
        elif sys.argv[1] == "counts":
            extract_counts(*sys.argv[2:])
        else:
            raise ValueError("Expected validate or counts")
    except (ValueError, OSError) as exc:
        sys.exit(f"ERROR: {exc}")
