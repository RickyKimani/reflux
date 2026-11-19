# Reflux CLI milestones

This is the implementation roadmap for Reflux's custom CLI module. The CLI
should provide a Cobra- and Fang-inspired experience without becoming a
general-purpose framework before Reflux needs one.

The migration is incremental: existing commands must remain usable after
every milestone.

Framework experiments may use the ignored root `cli_test/` harness. It is for
generic CLI behavior only; committed module tests belong in `cmdline/`, and
application integration tests remain in the repository root or their owning
module. No experiment in `cli_test/` may require changes to `main.v`.

## Scope and boundaries

The module lives at `cmdline/` and declares `module cmdline`, not `cli`: V
prefers a module in the project to vlib's, and `cli` would hide vlib's `cli`.
Code in the repository imports it as:

```v
import cmdline
```

The CLI module owns command metadata, parsing, dispatch, help rendering and
terminal presentation. It must not import `flowsheet`, `reactors`, `kinetics`
or other domain modules.

The root `main.v` remains the application adapter. It owns:

- process entry and `exit`;
- Reflux command handlers;
- mapping outcomes to exit codes;
- application-specific errors and reports;
- compatibility with `main_test.v`.

The existing exit-code contract is:

| Code | Meaning |
|---:|---|
| 0 | The command completed successfully |
| 1 | The command ran but a runtime, flowsheet or output operation failed |
| 2 | The invocation is invalid |

The parser and command framework must never call `exit()`. Only `main()` may
terminate the process.

## Design principles

1. Preserve existing command syntax and behavior before improving presentation.
2. Keep stdout, stderr, parsing and execution separately testable.
3. Keep JSON and CSV machine-readable and free of terminal styling.
4. Do not add a dependency until the standard library and existing code are
   insufficient.
5. Prefer explicit structs and interfaces over per-command closures. V does
   not free closures, so a command that runs repeatedly must not allocate a
   closure on every invocation.
6. Every completed milestone has focused tests and leaves the full CLI usable.

## C0 — Contract and baseline

Before changing command behavior:

- Treat `main_test.v` as the compatibility specification.
- Record stdout versus stderr behavior.
- Record the `0`, `1`, and `2` exit statuses.
- Preserve `--json FILE`, `--json=FILE`, `--csv FILE`, `--trace`, aliases and
  options after positional arguments.
- Identify which help text is temporary and which is user-facing.

**Exit criteria**

- The current CLI tests pass unchanged.
- The command list and syntax are written down in this file.
- No new parser behavior is introduced yet.

## C1 — Terminal and presentation foundation

Add the smallest presentation layer needed by help and diagnostics:

- terminal capability detection;
- `NO_COLOR` support;
- plain output for redirected streams;
- ANSI styling only when appropriate;
- optional Unicode symbols with ASCII fallbacks.

Keep the output abstraction injectable. The existing root `Console` remains
the compatibility layer while this work is developed.

A possible internal shape is:

```v
module cmdline

pub struct Terminal {
pub:
	interactive bool
	color       bool
	width       int
}

pub struct Style {
pub:
	terminal Terminal
}

pub fn (s Style) heading(text string) string {
	if !s.terminal.color {
		return text
	}
	return '\x1b[1;36m${text}\x1b[0m'
}
```

The exact palette is not part of the command contract. Keep it centralized
so it can change without touching command handlers.

**Exit criteria**

- Styled output is readable in an interactive terminal.
- `NO_COLOR=1` and redirected output contain no ANSI escapes.
- Windows output remains readable.
- Tests cover color enabled, color disabled and plain output.

**Status — done**

- `Theme` names an SGR code for each part of the output, and `reflux_theme`
  is reflux's, in xterm's 256 colours: copper for the heat (headings, drawn
  as labels, and the program's name), teal for the cooling water (commands),
  cornflower blue for flags, and greys for placeholders, comments and
  punctuation. Fang's layout, not Charm's colours.
- `terminal()` and `stderr_terminal()` detect each stream on its own, since
  stdout can be a file while stderr is the terminal. `NO_COLOR` turns colour
  off, and `CLICOLOR_FORCE` on into a pipe (`reflux help | less -R`). On
  Windows, V's runtime turns on escape processing only when stdout is a
  console, and `term` reads that mode, so stderr never gets escapes it can't
  show.
- A row's width is its plain label's, in visible characters, so escapes and
  non-ASCII labels don't misalign help.

## C2 — Command metadata and help

Define a command tree as data before implementing the parser. The first
version only needs to render help.

```v
module cmdline

pub struct Command {
pub:
	name        string
	aliases     []string
	summary     string
	description string
	usage       string
	examples    []string
	commands    []Command
	flags       []Flag
}

pub struct Flag {
pub:
	name        string
	alias       string
	description string
	value       bool
}
```

Use the tree for root help and command-specific help:

```text
reflux
├── solve
├── demo
├── kinetics
│   ├── check
│   └── explain
├── version
└── help
```

Help should have deterministic sections:

1. description;
2. usage;
3. commands;
4. arguments;
5. options;
6. examples;
7. the command's own sections, such as its exit status.

The description comes first, as in Fang, and has no heading; every other
section has one, in capitals (`USAGE`), drawn as a label in colour.

**Exit criteria**

- `reflux help` is rendered from command metadata.
- Command-specific help can be rendered without executing a command.
- Help works in colored and plain modes.
- Help tests do not depend on process exit.

**Status — done**

- `Command`, `Argument` and `Flag` metadata drive deterministic root and
  nested help rendering.
- All documented sections are covered: usage, description, commands,
  arguments, options and examples.
- Plain, styled and `--help`/`-h` command-specific rendering is covered by
  committed tests in `cmdline/cmdline_test.v`; visual output remains covered by the
  ignored `cli_test/` harness.
- No production files were changed; `main.v` and `main_test.v` remain the
  application migration boundary.
- Later: a `Section` is a titled block after the examples, an `Argument` can
  be `variadic` (the last only), and a `Flag` has a `value_name`. Usage
  lines and examples are highlighted: the program, then its commands, flags,
  and in a usage line its placeholders; an example's `#` line is a comment.
- Help follows the commands a request names and stops at the first word
  that isn't one, so `reflux solve plant.toml --help` is solve's help, not
  an unknown command `plant.toml`. `help` where a command is expected is a
  request too: `reflux help solve`, `reflux kinetics help check`.

## C3 — Parser core

Implement parsing independently from application execution.

Start with:

- root and nested command routing;
- positional arguments;
- long flags;
- boolean flags;
- flags with separate values;
- `--flag=value`;
- unknown command and flag errors;
- missing-value errors;
- duplicate flags.

The parser should return data or a typed error:

```v
module cmdline

pub struct Invocation {
pub:
	command    []string
	arguments  []string
	values     map[string]string
	switches   map[string]bool
}

pub struct ParseError {
pub:
	message string
	usage   string
}

pub fn parse(app Command, args []string) !Invocation {
	// Parsing only: no files are opened and no simulation is performed.
	return error('not implemented')
}
```

Use a result error rather than printing from the parser. The application
adapter decides whether an error is a usage error and how it is displayed.

**Exit criteria**

- Parser unit tests cover valid and invalid forms.
- Parsing never calls `exit()`, `println()` or `eprintln()`.
- Existing command syntax is accepted.
- Invalid invocations can be converted to exit code 2 by `main.v`.

**Status — done**

- Root and nested command routing, command aliases, positional arguments,
  long and short flags, separate and inline flag values, boolean switches,
  `--`, negative values, suggestions, missing values and duplicate detection
  are implemented.
- Empty inline values remain valid values and still participate in duplicate
  detection.
- Parser failures are typed `ParseError` values through the repository's
  `IError` convention; their messages remain stable for the application
  adapter and tests.
- Committed tests in `cmdline/cmdline_test.v` cover valid and invalid parser forms;
  the ignored `cli_test/` harness demonstrates the parsed invocation and
  errors without process termination.
- The parser returns result errors and performs no output or process
  termination. Production exit-code mapping remains intentionally deferred to
  the `main.v` migration boundary.

**Bugs fixed after it was marked done**, each with a test:

- An unknown option was an untyped error, not a `ParseError`; now every
  failure is one, and carries its command's usage line for a diagnostic.
- A negative number (`reflux demo 300 -1`) and a lone `-` were taken for
  unknown options, so "negative values" held only for a flag's value.
- A flag that takes a value swallowed the next flag (`--csv --trace` set
  csv to `--trace`); a flag there now means the value is missing, as
  `main.v` has it.
- `--trace=` was accepted: an empty inline value on a switch went unseen.
- A group with no command (`reflux kinetics`) parsed, and failed only at
  dispatch, as "no handler for kinetics": now it is a ParseError, and a
  handler can't be registered for a group.
- `-- solve` said "unknown command solve (did you mean solve?)".
- Aliases matched with their dashes stripped and names didn't, so `--run`
  ran the command aliased `run`, and `-version` ran `version`. Commands now
  match exactly; a long flag is `--name` and a short one `-alias` only.
- Suggestions took the one nearest name within three edits, so `zzz`
  suggested `fit` and any one-letter flag the first long one. Now up to
  three within two edits, as `main.v`'s did, against long flags only.
- Messages follow `main.v`'s: "no command", "no option".
- Validation allowed a command's dashed name (`--version`) only as its first
  alias, and named a nested command by its parent and index (`root
  reflux[0]`); now anywhere in the list, and by its path (`reflux kinetics
  check`).

## C4 — Dispatcher and application adapter

Add execution after parsing is stable. A handler receives a typed invocation
and an explicit context rather than raw process arguments.

```v
module cmdline

pub struct Context {
pub:
	invocation Invocation
}

pub type Handler = fn (Context) !int

pub struct RegisteredCommand {
pub:
	command Command
	handler Handler
}
```

If V's function-value requirements make this shape unsuitable, use a small
interface or explicit handler structs. Do not hide application state in
global variables.

The root application adapter should remain thin:

```v
fn main() {
	mut c := Console{live: true}
	exit(run(os.args[1..], mut c))
}

fn run(args []string, mut c Console) int {
	// Parse with cli, then map parse/runtime outcomes to 0, 1 or 2.
	// Reflux-specific handlers remain in this module.
	return 0
}
```

**Exit criteria**

- `version` and `help` use the new framework.
- `run()` remains callable from tests.
- Process termination remains confined to `main()`.
- The existing `Console` behavior is preserved.

**Status — done**

- `App` validates metadata, routes built-in help/version outcomes and
  dispatches registered commands without printing or terminating.
- Handlers receive a typed `Context` containing the parsed invocation rather
  than raw arguments; the interface form avoids per-invocation closures.
- Missing handlers and handler failures remain result errors for the
  application adapter to map to exit codes.
- Committed tests cover help/version outcomes, successful dispatch and
  missing-handler errors; the ignored harness demonstrates the same behavior
  with visible mock output.
- Production `main.v` and `main_test.v` remain unchanged, so the existing
  `Console` contract is preserved until the C5 migration.

## C5 — Incremental command migration

Migrate commands in increasing complexity:

1. `version`;
2. `help`;
3. `demo`;
4. `solve`;
5. `kinetics check`;
6. `kinetics explain`;
7. `study`, `optimize` and `fit`, which M8 added after this plan was written.

For each command:

- add command metadata;
- add parser tests;
- adapt the existing handler to typed arguments;
- retain application behavior;
- retain `main_test.v` coverage;
- remove only the old parsing code that is no longer reachable.

`solve` must retain:

- flowsheet loading and solving;
- trace output;
- JSON and CSV reports;
- file collision checks;
- failure reports;
- stale CSV removal after failure.

**Exit criteria**

- All existing command tests pass.
- No command silently changes its exit status.
- Machine-readable output contains no styling.
- Old parser code has been removed only after tests prove it is unused.

**Progress**

- Help is migrated. `commands.v` holds reflux's whole tree, every command
  with its summary, description, usage, arguments, flags and examples, and
  `run()` renders help from it for `reflux help solve`, `reflux solve
  --help`, `reflux kinetics help` and the rest: in colour on a terminal,
  plain in `main_test.v`. `reflux` alone prints it to stderr and exits 2, as
  before. The hand-written `usage` and `kinetics_usage` are gone: both are
  now the tree's help without colour, which errors still print after what
  was wrong, so `main_test.v`'s checks of them hold unchanged.
- `reflux solve --help` and the like were errors, exit 2; now they are help.
- `version` prints what it did; it doesn't need the framework yet.
- The commands still parse their own arguments, so every error message is
  `main.v`'s. `main_test.v` changed only where it read the old help's first
  line, `usage:`, and gained a test of every command's help.

## C6 — Diagnostics and polish

Once behavior is migrated, improve presentation:

- command and flag suggestions;
- consistent error prefixes;
- file and line locations;
- command-specific usage on errors;
- terminal-width-aware help;
- consistent headings, labels and spacing;
- optional progress indicators only for operations that benefit from them.

Diagnostics should remain useful when copied into an editor or redirected to
a file. Do not use color as the only way to distinguish severity.

Current framework progress:

- `Diagnostic` renders error, warning and info messages through the shared
  style layer without writing to a stream.
- Validation rejects malformed argument metadata and collisions between flag
  names and aliases before help, parsing or dispatch.
- Required positional arguments must precede optional positional arguments,
  matching the parser's arity model.
- Command groups cannot declare positional arguments or flags until inherited
  and global option parsing is implemented.
- Command aliases are checked symmetrically, while dashed aliases such as
  `--version` remain allowed for their canonical command.
- The ignored `cli_test/` harness covers visual behavior, and committed
  `cmdline/cmdline_test.v` covers deterministic help, parsing, validation and
  diagnostics without requiring the application entry point.
- The committed tests also cover the application boundary's help/version
  outcomes, typed handler dispatch and missing-handler errors.
- Diagnostics can include a command-specific usage line while remaining
  render-only and suitable for redirected output.

**Exit criteria**

- Invalid commands and flags provide useful suggestions.
- Diagnostics are readable with color disabled.
- Existing file and line locations are unchanged.
- Output is stable enough for deterministic tests.

## C7 — Executable integration and optional features

Add a small number of process-level tests for the actual executable. Keep most
tests in-process for speed and precise output capture.

Cover:

- exit code 0;
- runtime failure with exit code 1;
- invalid invocation with exit code 2;
- stdout and stderr separation;
- `--help` and `--version`;
- redirected output without ANSI escapes.

Optional later features:

- shell completion;
- generated Markdown or man-page help;
- command groups;
- environment-backed defaults;
- persistent global flags.

These are not prerequisites for the CLI migration.

**Exit criteria**

- The executable-level contract is tested on the supported platforms.
- The CLI framework remains smaller than the application logic it serves.
- Optional features have explicit tests before being considered complete.

## Initial implementation order

The first practical sequence is:

1. create `cmdline/cmdline_test.v`;
2. define minimal command metadata;
3. register the current command tree;
4. render root help;
5. add plain-output and color tests;
6. wire `help` through `run()`;
7. migrate `version`;
8. begin parser work.

Do not begin by rewriting `solve`. Proving the command tree, help renderer,
output policy and test seam with `help` and `version` keeps the first changes
small and reversible.

## Proposed files

Create files only when their responsibility has become real:

```text
cmdline/
├── milestones.md       # this roadmap
├── command.v           # command and metadata types
├── help.v              # help rendering
├── terminal.v          # TTY, color and width detection
├── style.v             # styling primitives
├── parser.v            # argument and flag parsing
├── error.v             # parser and framework errors
└── cmdline_test.v      # framework tests
```

`Console` should remain in `main.v` until there is a demonstrated need for a
shared terminal abstraction. If it later moves into `cmdline/`, make that a
separate migration with tests rather than combining it with parser work.
