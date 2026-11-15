module flowsheet

import flash

// Failure says why a stream or unit failed.
pub enum Failure {
	bad_specification // a value it can't use: a negative flow, a split of 1.2, …
	bad_flowsheet     // units and streams that don't join up: an open port, say
	flash_failed      // the flash it needed failed, for the reason in cause
	no_convergence    // a block's tears were still moving when the passes ran out
	unbalanced        // a balance didn't close: a bug, or a unit that's wrong
	bad_file          // a flowsheet file that can't be read: see line and key
	reaction_failed   // a reactor's equilibrium, integration or Newton failed (M9)
}

// FlowsheetError is every failure in this module, never a panic (D8). Test
// for it with `if err is FlowsheetError`.
pub struct FlowsheetError {
	Error
pub:
	kind Failure
	// unit is the unit at fault, by name, or empty for a feed.
	unit   string
	reason string
	// cause is the flash's own failure, such as .more_phases, for .flash_failed.
	cause ?flash.Failure
	// history is how the tears had converged, if it happened during a solve.
	history []TearHistory
	// line and key are where in a file it went wrong, if it was in one.
	line        int
	key         string
	suggestions []string // what a misspelt name may have meant
}

// msg is the reason, marked as the flowsheet's, with where it went wrong and
// what a misspelt name may have meant.
pub fn (e FlowsheetError) msg() string {
	mut place := ''
	if e.line > 0 {
		place += 'line ${e.line}: '
	}
	if e.key != '' {
		place += '${e.key}: '
	} else if e.unit != '' {
		place += '${e.unit}: '
	}
	hint := if e.suggestions.len > 0 { ' (did you mean ${e.suggestions.join(', ')}?)' } else { '' }
	return 'flowsheet: ${place}${e.reason}${hint}'
}

fn bad_spec(unit string, reason string) IError {
	return FlowsheetError{
		kind:   .bad_specification
		unit:   unit
		reason: reason
	}
}

fn bad_sheet(unit string, reason string) IError {
	return FlowsheetError{
		kind:   .bad_flowsheet
		unit:   unit
		reason: reason
	}
}

fn file_error(line int, key string, reason string, suggestions []string) IError {
	return FlowsheetError{
		kind:        .bad_file
		reason:      reason
		line:        line
		key:         key
		suggestions: suggestions
	}
}

// flash_failed passes on a flash's failure, naming the unit it happened in.
fn flash_failed(unit string, err IError) IError {
	return FlowsheetError{
		kind:   .flash_failed
		unit:   unit
		reason: err.msg()
		cause:  if err is flash.FlashError { err.kind } else { none }
	}
}
