from pylint.reporters.text import TextReporter
import os.path
from pathlib import Path


class PyLintWriter(object):
    def __init__(self):
        self.content = []

    def write(self, st):
        self.content.append(st)

    def read(self):
        return self.content


def run_checker(fixture_name):
    writer = PyLintWriter()
    reporter = TextReporter(writer)
    from archimedes.checker import run

    success = run(
        [os.path.join(os.path.dirname(__file__), "fixtures", fixture_name)],
        reporter=reporter,
    )
    return writer.read(), success


def test_happy_path():
    output, success = run_checker("my_test_script.py.donotrun")
    assert output == []
    assert success


def test_sad_path():
    bad_output, success = run_checker("my_bad_script.py.donotrun")
    assert (
            "tests/fixtures/my_bad_script.py.donotrun:1:0: W9001: Use only safe imports. (non-safe-imports)"
            in bad_output
    )
    assert (
            "tests/fixtures/my_bad_script.py.donotrun:2:0: W9001: Use only safe imports. (non-safe-imports)"
            in bad_output
    )
    assert (
            "tests/fixtures/my_bad_script.py.donotrun:7:4: W9003: Use no unsafe globals or unsafe names (no-unsafe-globals)"
            in bad_output
    )
    assert (
            "tests/fixtures/my_bad_script.py.donotrun:8:13: W9003: Use no unsafe globals or unsafe names (no-unsafe-globals)"
            in bad_output
    )
    assert (
            "tests/fixtures/my_bad_script.py.donotrun:11:4: W9003: Use no unsafe globals or unsafe names (no-unsafe-globals)"
            in bad_output
    )
    assert (
            "tests/fixtures/my_bad_script.py.donotrun:13:4: W9003: Use no unsafe globals or unsafe names (no-unsafe-globals)"
            in bad_output
    )
    assert (
            "tests/fixtures/my_bad_script.py.donotrun:14:12: W9004: Use no private members (no-private-members)"
            in bad_output
    )
    assert not success
