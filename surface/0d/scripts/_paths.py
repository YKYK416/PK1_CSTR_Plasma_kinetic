"""Paths shared by the Hong surface closed-0D scripts."""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[3]
LEGACY_HONG_ROOT = REPO_ROOT / "Hong"
CORE_RUNTIME = REPO_ROOT / "core" / "runtime"

if str(CORE_RUNTIME) not in sys.path:
    sys.path.insert(0, str(CORE_RUNTIME))
