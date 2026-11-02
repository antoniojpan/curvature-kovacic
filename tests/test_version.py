"""The three version literals of the repository must agree.

`pyproject.toml`, `CITATION.cff` and `curvature_kovacic.__version__` are
independent literals and used to drift: v1.0.2 shipped the first two at 1.0.1
while the tag, `CITATION.cff` and the paper said 1.0.2.

Deriving `__version__` from the installed metadata is not an option here: with
an editable install, `importlib.metadata.version()` reports the `.dist-info`
generated at install time, not `pyproject.toml`, so the fix would replace a
transcription error with a stale-cache error. Hence the check below.

`tomllib` is not available on Python 3.10, so the values are read with regexes.
"""

import re
from pathlib import Path

import curvature_kovacic

ROOT = Path(__file__).resolve().parents[1]


def _first(text: str, pattern: str, origin: str) -> str:
    match = re.search(pattern, text, re.MULTILINE)
    assert match is not None, f"no version found in {origin}"
    return match.group(1)


def test_version_matches_pyproject():
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    declared = _first(text, r'^version\s*=\s*"([^"]+)"', "pyproject.toml")
    assert curvature_kovacic.__version__ == declared


def test_version_matches_citation():
    text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    cited = _first(text, r"^version:\s*(\S+)", "CITATION.cff")
    assert curvature_kovacic.__version__ == cited
