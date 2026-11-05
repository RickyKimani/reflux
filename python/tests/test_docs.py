"""docs/python.md's examples run (M10a, exit criterion 11). Its python blocks run
in order, in one namespace, as a session would, and a block that prints is
followed by a text block of exactly what it prints. So the page can't drift
from the package: an example that stops working, or prints something else,
fails here."""

import contextlib
import io
import re

from support import ROOT

DOC = ROOT / "docs" / "python.md"

# A fenced block: its language, and its body. Fences aren't indented on the
# page, and none has another fence inside it.
BLOCK = re.compile(r"^```(\w*)\n(.*?)^```\n", re.M | re.S)


def blocks(text):
    """Each fenced block's language, body, and the line its body starts on."""
    return [(m[1], m[2], text.count("\n", 0, m.start()) + 2) for m in BLOCK.finditer(text)]


def test_every_example_runs_and_prints_what_the_page_shows():
    found = blocks(DOC.read_text(encoding="utf-8"))
    namespace = {"__name__": "python_md"}
    ran = 0
    for k, (lang, body, line) in enumerate(found):
        if lang != "python":
            continue
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            # Named by its place on the page, so a traceback points there.
            exec(compile(body, f"{DOC.relative_to(ROOT).as_posix()}:{line}", "exec"), namespace)
        ran += 1
        shown = found[k + 1] if k + 1 < len(found) and found[k + 1][0] == "text" else None
        if shown is None:
            assert out.getvalue() == "", f"the example at line {line} prints, but the page " \
                f"shows no output after it:\n{out.getvalue()}"
        else:
            assert out.getvalue() == shown[1], f"the example at line {line} printed\n" \
                f"{out.getvalue()}\nbut the page shows, at line {shown[2]}:\n{shown[1]}"
    assert ran >= 10, f"only {ran} examples: is the page still there?"


def test_the_runner_finds_blocks_and_their_lines():
    text = "intro\n```python\nx = 1\n```\nprose\n```text\n1\n```\n"
    assert blocks(text) == [("python", "x = 1\n", 3), ("text", "1\n", 7)]
