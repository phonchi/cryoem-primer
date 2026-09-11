#!/usr/bin/env python3
"""Validate local links, fragments, required pages, and image alt text."""

from __future__ import annotations

import argparse
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


REQUIRED_PAGES = {
    "index.html",
    "00_intro.html",
    "01_image_basics.html",
    "02_filter_segment.html",
    "03_fourier.html",
    "04_wavelet.html",
    "05_background.html",
    "05_image_formation.html",
    "05_statistical_inference.html",
    "06_alignment_classification.html",
    "06_dimension_reduction.html",
    "06_heterogeneity.html",
    "06_resolution_validation.html",
    "06_workflow.html",
    "06_reconstruction_validation.html",
    "07_synthetic_data.html",
    "08_resources.html",
    "appendix_conventions.html",
}


class PageParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []
        self.ids: set[str] = set()
        self.images_without_alt: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        if values.get("id"):
            self.ids.add(values["id"] or "")
        if tag == "a" and values.get("href"):
            self.links.append(values["href"] or "")
        if tag == "img" and not (values.get("alt") or "").strip():
            self.images_without_alt.append(values.get("src") or "<unknown>")


def parse_page(path: Path) -> PageParser:
    parser = PageParser()
    parser.feed(path.read_text(encoding="utf-8", errors="replace"))
    return parser


def local_target(root: Path, source: Path, href: str) -> tuple[Path, str] | None:
    parts = urlsplit(href)
    if parts.scheme or parts.netloc or href.startswith(("mailto:", "javascript:", "#")):
        return None
    raw_path = unquote(parts.path)
    if not raw_path:
        return source, parts.fragment
    target = (source.parent / raw_path).resolve()
    if raw_path.endswith("/"):
        target = target / "index.html"
    if target.suffix == "":
        html_target = target.with_suffix(".html")
        target = html_target if html_target.exists() else target
    try:
        target.relative_to(root.resolve())
    except ValueError:
        return target, parts.fragment
    return target, parts.fragment


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("html_dir", type=Path)
    args = parser.parse_args()
    root = args.html_dir.resolve()
    html_files = sorted(root.rglob("*.html"))
    parsed = {path: parse_page(path) for path in html_files}
    errors: list[str] = []

    missing_pages = sorted(REQUIRED_PAGES - {p.name for p in html_files if p.parent == root})
    errors.extend(f"missing required page: {name}" for name in missing_pages)

    for source, page in parsed.items():
        for image in page.images_without_alt:
            errors.append(f"{source.relative_to(root)}: image without alt text: {image}")
        for href in page.links:
            resolved = local_target(root, source, href)
            if resolved is None:
                continue
            target, fragment = resolved
            if not target.exists():
                errors.append(f"{source.relative_to(root)}: missing local target {href}")
                continue
            if fragment and target.suffix.lower() == ".html":
                target_page = parsed.get(target)
                if target_page is None:
                    target_page = parse_page(target)
                    parsed[target] = target_page
                if unquote(fragment) not in target_page.ids:
                    errors.append(
                        f"{source.relative_to(root)}: missing fragment #{fragment} in {target.relative_to(root)}"
                    )

    if errors:
        print("Site validation failed:")
        for error in errors:
            print(f"- {error}")
        return 1
    print(f"Site validation passed: {len(html_files)} HTML files, {len(REQUIRED_PAGES)} required pages")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
