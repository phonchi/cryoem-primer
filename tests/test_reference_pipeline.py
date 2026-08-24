import argparse
import importlib.util
import json
import tempfile
import unittest
from unittest import mock
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "reference_pipeline.py"
SPEC = importlib.util.spec_from_file_location("reference_pipeline", MODULE_PATH)
rp = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(rp)


class ReferencePipelineTests(unittest.TestCase):
    def test_registry_is_public_and_marks_exact_duplicate(self):
        sha = "a" * 64
        rows = [
            {
                "title": "cache extract title",
                "source_path": "/private/machine/References/2016/file.pdf",
                "text_path": "/private/machine/cache/file.txt",
                "relative_path": "2016/Chapter-Five---Processing-of-Cryo-EM-Movie-Data_2016_Methods-in-Enzymology.pdf",
                "sha256": sha,
                "pages": 22,
                "size_bytes": 123,
                "status": "ok",
                "keyword_hits": ["ctf"],
                "summary_extract": "copyrighted extracted prose",
            },
            {
                "title": "cache extract title",
                "source_path": "/private/machine/References/file.pdf",
                "text_path": "/private/machine/cache/file2.txt",
                "relative_path": "Chapter-Five---Processing-of-Cryo-EM-Movie-Data_2016_Methods-in-Enzymology.pdf",
                "sha256": sha,
                "pages": 22,
                "size_bytes": 123,
                "status": "ok",
                "keyword_hits": ["ctf"],
                "summary_extract": "copyrighted extracted prose",
            },
        ]
        registry = rp.public_registry(rows)
        self.assertEqual(registry[0]["content_role"], "substantive")
        self.assertEqual(registry[0]["source_id"], "mie2016-ch05")
        self.assertEqual(registry[1]["content_role"], "duplicate_alias")
        serialized = json.dumps(registry)
        self.assertNotIn("/private/machine", serialized)
        self.assertNotIn("copyrighted extracted prose", serialized)

    def test_validate_marks_changed_claim_hash_needs_reverify(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            repo = root / "repo"
            sources = root / "References"
            cache = sources / "document_cache" / "references"
            catalog = repo / "notes" / "reference_catalog"
            digest = repo / "notes" / "ref_digest" / "digest.md"
            for directory in (cache, catalog, digest.parent, repo / "book"):
                directory.mkdir(parents=True, exist_ok=True)
            pdf = sources / "source.pdf"
            pdf.write_bytes(b"new PDF bytes")
            live_sha = rp.sha256_file(pdf)
            (cache / "document_index.jsonl").write_text("{}\n")
            rp.write_jsonl(catalog / "source_registry.jsonl", [{
                "source_id": "source", "relative_path": "source.pdf", "sha256": live_sha,
                "content_role": "substantive", "coverage_status": "full_digest",
                "digest_path": "notes/ref_digest/digest.md", "citation_key": "source",
            }])
            rp.write_jsonl(catalog / "claim_map.jsonl", [{
                "claim_id": "claim", "source_id": "source", "source_sha256": "0" * 64,
                "chapter": "01", "anchor": "a", "claim": "c", "locator": "p. 1", "scope": "s",
            }])
            digest.write_text("# digest\n")
            (catalog / "asset_manifest.yml").write_text("schema_version: 1\nassets:\n")
            (repo / "book" / "references.bib").write_text("@article{source, title={Source}}\n")
            args = argparse.Namespace(
                repo_root=repo, source_root=sources, cache_root=sources / "document_cache",
                catalog_dir=catalog, digest_dir=digest.parent, bib=repo / "book" / "references.bib",
            )
            report = rp.validate_paths(args, emit=False)
            self.assertEqual(report["errors"], [])
            self.assertEqual(report["needs_reverify"], ["claim"])

    def test_source_hash_change_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "paper.pdf"
            source.write_bytes(b"before")
            before = rp.live_pdf_state(Path(tmp))["paper.pdf"]["sha256"]
            source.write_bytes(b"after")
            after = rp.live_pdf_state(Path(tmp))["paper.pdf"]["sha256"]
            self.assertNotEqual(before, after)

    def test_claim_hash_is_not_silently_rebased(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "claims.jsonl"
            claims = [{"claim_id": "c", "source_id": "s", "chapter": "01", "anchor": "a", "claim": "claim", "locator": "p. 1", "scope": "scope"}]
            with mock.patch.object(rp, "CLAIMS", claims):
                rp.seed_claims([{"source_id": "s", "sha256": "a" * 64, "content_role": "substantive"}], path)
                rp.seed_claims([{"source_id": "s", "sha256": "b" * 64, "content_role": "substantive"}], path)
            row = rp.read_jsonl(path)[0]
            self.assertEqual(row["source_sha256"], "a" * 64)


if __name__ == "__main__":
    unittest.main()
