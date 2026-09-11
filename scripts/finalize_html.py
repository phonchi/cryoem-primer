#!/usr/bin/env python3
"""Remove identical duplicate sphinx-thebe initializers from generated HTML.

sphinx-thebe 0.3.1 can register its inline const block more than once during a
Jupyter Book build. Re-declaring the same constants raises a browser SyntaxError.
Only byte-identical initializer bodies are deduplicated; distinct configurations
fail explicitly, and all chapter content and other scripts remain untouched.
"""
import argparse
from pathlib import Path
import re

THEBE_CONFIG = re.compile(r'<script>\s*(const THEBE_JS_URL\s*=[^<]*)</script>')


def deduplicate_thebe_config(html: str) -> tuple[str, int]:
    seen = None
    removed = 0

    def replace(match):
        nonlocal seen, removed
        body = match.group(1)
        if seen is None:
            seen = body
            return match.group(0)
        if body != seen:
            raise ValueError('Conflicting Thebe initializers in generated HTML')
        removed += 1
        return ''

    return THEBE_CONFIG.sub(replace, html), removed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('html_dir', type=Path)
    args = parser.parse_args()
    total = 0
    for path in sorted(args.html_dir.rglob('*.html')):
        original = path.read_text(encoding='utf-8')
        cleaned, count = deduplicate_thebe_config(original)
        if count:
            path.write_text(cleaned, encoding='utf-8')
            total += count
    print(f'Removed {total} identical duplicate Thebe initializer blocks')


if __name__ == '__main__':
    main()
