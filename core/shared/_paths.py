"""Shared paths for tracked code and the local, untracked Hong workspace."""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
LEGACY_HONG_ROOT = REPO_ROOT / "Hong"
CORE_RUNTIME = REPO_ROOT / "core" / "runtime"

if str(CORE_RUNTIME) not in sys.path:
    sys.path.insert(0, str(CORE_RUNTIME))
