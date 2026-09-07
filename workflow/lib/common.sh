#!/usr/bin/env bash
# Shared validation. Sourced by the numbered Bash entry points.
validate_value_options() {
  local value_options="$1"; shift
  while (( $# )); do
    if [[ " $value_options " == *" $1 "* ]]; then
      if (( $# < 2 )) || [[ -z "$2" || "$2" == --* || "$2" == -h ]]; then
        printf 'ERROR: %s requires a value\n' "$1" >&2
        exit 2
      fi
      shift 2
    else
      shift
    fi
  done
}
positive_int() {
  [[ "$2" =~ ^[1-9][0-9]*$ ]] || {
    printf 'ERROR: %s must be a positive integer (got %s)\n' "$1" "$2" >&2
    exit 2
  }
}
find_unique_read() {
  local helper_dir
  helper_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
  python3 "$helper_dir/find_read.py" "$@"
}
