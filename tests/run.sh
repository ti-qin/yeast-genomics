#!/usr/bin/env bash
# Lightweight validation; does not launch real bioinformatics workloads.
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
for script in workflow/*.sh; do bash -n "$script"; done
python3 - <<'PY'
import ast
from pathlib import Path
for folder in ('workflow/lib', 'tests'):
    for file in Path(folder).glob('*.py'):
        ast.parse(file.read_text(), filename=str(file))
print('Bash and Python syntax checks passed')
PY
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
