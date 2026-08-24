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
        "單顆粒",
        "訊噪",
        "可視化",
        "基於區域",
        "冷凍電顯",
        "模板比對",
        "教材目前",
        "不會假裝",
        "老師維護",
        "provenance",
        "run_manifest",
        "sha256_file",
        "volume_sha256",
        "Jupyter Book",
        "_build/synthetic",
        "可辯護",
        "主張範圍",
        "下游",
        "短暂",
        "調制",
        "默認",
        "生成",
    )
    sources = list(BOOK.glob("*.md")) + list(BOOK.glob("*.py"))
    for source in sources:
        text = source.read_text(encoding="utf-8")
        for token in forbidden:
            assert token not in text, f"{source.name}: forbidden token {token}"
        assert not re.search(r"(?<!演)算法", text), f"{source.name}: 中國用語「算法」"
        assert not re.search(
            r"不是.{0,40}而是|不只是|不代表|不等於|不能只|不能單靠",
            text,
        ), f"{source.name}: 請改成直接說明成立條件與限制"


def test_student_facing_title_and_toc_terms():
    config = yaml.safe_load((BOOK / "_config.yml").read_text(encoding="utf-8"))
    toc = yaml.safe_load((BOOK / "_toc.yml").read_text(encoding="utf-8"))
    assert config["title"] == "從影像處理到 Cryo-EM 單粒子分析"
    assert any(part["caption"] == "Cryo-EM 單粒子分析" for part in toc["parts"])


def test_generated_outputs_are_not_in_book_tree():
    assert not (BOOK / "data" / "output").exists()
