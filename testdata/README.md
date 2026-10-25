# testdata

Fixtures that tests read, borrowed from Go's `testdata/` convention. V has
no convention of its own here.

Rules:

- Tests only **read** from here. Nothing in `v test` writes to it.
- Every file records where each value came from: the source (book and
  table, paper, database and version, or tool and version), and the date
  it was generated.
- Generated fixtures are produced by a script in [`../tools`](../tools) and
  committed. Tests never need Python, CoolProp or network access to run.
- JSON is the default format, read with V's `json2` module.
