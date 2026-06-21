#!/usr/bin/env bash
set -Eeuo pipefail

script_dir="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
exec python3 "${script_dir}/metamarker_profile.py" "$@"
