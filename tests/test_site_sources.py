from pathlib import Path
import re

import yaml


ROOT = Path(__file__).resolve().parents[1]
BOOK = ROOT / "book"


def test_toc_sources_exist():
    toc = yaml.safe_load((BOOK / "_toc.yml").read_text(encoding="utf-8"))
    files = [toc["root"]]
    for part in toc["parts"]:
        files.extend(chapter["file"] for chapter in part["chapters"])
    for stem in files:
        assert any((BOOK / f"{stem}{suffix}").exists() for suffix in (".md", ".ipynb")), stem


def test_no_ai_tool_residue_or_mainland_terms():
    forbidden = (
        "turn0search",
        "utm_source=chatgpt",
        "以下是修改後的版本",
        "可視化",
        "基於區域",
        "冷凍電顯",
        "模板比對",
    )
    sources = list(BOOK.glob("*.md")) + list(BOOK.glob("*.py"))
    for source in sources:
        text = source.read_text(encoding="utf-8")
        language_tokens = forbidden if source.name != "appendix_conventions.md" else forbidden[:3]
        for token in language_tokens:
            assert token not in text, f"{source.name}: forbidden token {token}"
        if source.name != "appendix_conventions.md":
            assert not re.search(r"(?<!演)算法", text), f"{source.name}: 中國用語「算法」"


def test_generated_outputs_are_not_in_book_tree():
    assert not (BOOK / "data" / "output").exists()
