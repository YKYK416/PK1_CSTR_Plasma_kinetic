"""Paths shared by the Hong surface 0D-CSTR scripts."""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[3]
LEGACY_HONG_ROOT = REPO_ROOT / "Hong"

for path in (
    REPO_ROOT / "core" / "runtime",
    REPO_ROOT / "core" / "shared",
    REPO_ROOT / "surface" / "0d" / "scripts",
):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
