"""Threads (D5). Python threads call in at once, each registered with the
collector for the length of its call: without that, a collection from one
aborts with "Collecting from unknown thread", on Windows as on Linux. And the
library's first call, which starts V's runtime and mustn't race (reflux.h), is
the package's own, made by the import. Each test runs in a process of its
own, so that a crash is reported, not hidden."""

from support import isolated

CALLS = """
import threading
import reflux
from reflux import _native

wrong = []


def work():
    for _ in range(200):
        if reflux.collect() <= 0 or reflux.version() != reflux.__version__:
            wrong.append("collect or version")
        try:
            _native.free(10**15)
        except reflux.BadHandle as e:
            if str(e) != "1000000000000000 is not a handle":
                wrong.append(str(e))
        else:
            wrong.append("no BadHandle")


threads = [threading.Thread(target=work) for _ in range(8)]
for t in threads:
    t.start()
for t in threads:
    t.join()
print(len(wrong), wrong[:3])
"""


def test_threads_call_and_collect_at_once():
    assert isolated(CALLS).strip() == "0 []"


def test_an_import_from_a_thread_that_ends():
    # The import makes the first call, so V's runtime starts in that thread,
    # and here the thread ends before anything else calls in.
    code = ("import threading\n"
            "t = threading.Thread(target=lambda: __import__('reflux'))\n"
            "t.start()\n"
            "t.join()\n") + CALLS
    assert isolated(code).strip().endswith("0 []")
