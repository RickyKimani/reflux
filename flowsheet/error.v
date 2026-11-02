module flowsheet

import flash

// Failure says why a stream or unit failed.
pub enum Failure {
	bad_specification // a value it can't use: a negative flow, a split of 1.2, …
	flash_failed      // the flash it needed failed, for the reason in cause
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
}

// msg is the reason, marked as the flowsheet's and with the unit at fault.
pub fn (e FlowsheetError) msg() string {
	return if e.unit == '' { 'flowsheet: ${e.reason}' } else { 'flowsheet: ${e.unit}: ${e.reason}' }
}

fn bad_spec(unit string, reason string) IError {
	return FlowsheetError{
		kind:   .bad_specification
		unit:   unit
		reason: reason
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
