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
        "模型與參數的參考來源",
        "這組合成資料的適用範圍",
        "原始 2SDR 合成實驗交代",
        "結論只適用",
        "這個範例涵蓋的條件",
        "這段示範包含的條件",
        "常見誤解是",
        "報告 FSC 時",
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
        assert not re.search(
            r"(?<![A-Za-z])(?:p|pp|Eq|Eqs|Fig)\.\s*\d|supplement\s+pp\.",
            text.split("## 延伸閱讀", maxsplit=1)[0],
        ), f"{source.name}: 精確頁碼與方程式定位應移到延伸閱讀"


def test_protected_chapters_keep_original_understanding_checks():
    expected = {
        "01_image_basics.py": 5,
        "02_filter_segment.py": 5,
        "03_fourier.py": 5,
        "04_wavelet.py": 4,
        "07_synthetic_data.py": 5,
    }
    for name, count in expected.items():
        text = (BOOK / name).read_text(encoding="utf-8")
        section = text.split("## 理解檢查", maxsplit=1)[1]
        questions = re.findall(r"^(?:# )?[1-9]\. ", section, flags=re.MULTILINE)
        assert len(questions) == count, name
        assert section.count("```{dropdown} 參考答案") == count, name


def test_student_facing_title_and_toc_terms():
    config = yaml.safe_load((BOOK / "_config.yml").read_text(encoding="utf-8"))
    toc = yaml.safe_load((BOOK / "_toc.yml").read_text(encoding="utf-8"))
    assert config["title"] == "從影像處理到 Cryo-EM 單粒子分析"
    assert any(part["caption"] == "Cryo-EM 單粒子分析" for part in toc["parts"])


def test_generated_outputs_are_not_in_book_tree():
    assert not (BOOK / "data" / "output").exists()


def test_revised_pages_integrate_checks_and_include_reading():
    toc = yaml.safe_load((BOOK / "_toc.yml").read_text(encoding="utf-8"))
    spa = next(part for part in toc["parts"] if part["caption"] == "Cryo-EM 單粒子分析")
    stems = [chapter["file"] for chapter in spa["chapters"]]
    for stem in stems + ["00_intro", "appendix_conventions"]:
        text = (BOOK / f"{stem}.md").read_text(encoding="utf-8")
        assert "理解檢查" not in text, stem
        assert "參考答案" not in text, stem
        assert "姿態" not in text, stem
        assert "## 延伸閱讀" in text, stem
        assert "{cite}" in text, stem


def test_protected_chapters_and_assets_are_byte_identical():
    import hashlib
    import json
    manifest = json.loads((ROOT / "notes/spa_revision_20260911/protected_files.json").read_text())
    for relative, expected in manifest.items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == expected, relative
