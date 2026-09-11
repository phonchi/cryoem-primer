"""Regression for duplicate inline constants emitted by sphinx-thebe."""
import importlib.util
from pathlib import Path

path = Path(__file__).resolve().parents[1] / 'scripts/finalize_html.py'
spec = importlib.util.spec_from_file_location('finalize_html', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
CONFIG = '<script>const THEBE_JS_URL = "example.js"; const thebe_selector = ".cell"</script>'


def test_identical_thebe_config_is_only_declared_once():
    source = '<head>' + CONFIG + '<script src="sphinx-thebe.js"></script>' + CONFIG + '</head><article>body</article>'
    cleaned, count = module.deduplicate_thebe_config(source)
    assert count == 1
    assert cleaned.count(CONFIG) == 1
    assert '<script src="sphinx-thebe.js"></script>' in cleaned
    assert cleaned.endswith('<article>body</article>')
    assert module.deduplicate_thebe_config(cleaned) == (cleaned, 0)


def test_other_scripts_and_conflicting_config_are_not_silently_removed():
    import pytest
    other = '<script>window.other = 1</script>'
    source = CONFIG + other + other
    assert module.deduplicate_thebe_config(source) == (source, 0)
    with pytest.raises(ValueError, match='Conflicting'):
        module.deduplicate_thebe_config(CONFIG + CONFIG.replace('example.js', 'different.js'))
