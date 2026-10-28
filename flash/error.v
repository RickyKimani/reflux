module flash

// Failure says why a flash or stability test failed.
pub enum Failure {
	bad_input      // a specification it can't use: a negative mole fraction, T ≤ 0, …
	no_convergence // the iteration ran out of steps without an answer it could prove
	no_volume_root // the model has no volume root at a state it needed
	// more_phases: a phase of the two-phase split is itself unstable. There
	// are three or more phases, or a better split, and M3 handles two (M6).
	more_phases
	// phases_swapped: a vapour-fraction flash converged, but the phase it
	// gave the fraction β is the denser one, so by D6 the answer's vapour
	// fraction is 1 − β: a dew point where a bubble point was asked for,
	// say. It's how the other kind of boundary shows near a critical point,
	// where bubble and dew curves meet.
	phases_swapped
	// outside_envelope: a vapour-fraction flash found no answer, and no PT
	// flash across a wide range of the unknown splits, so the point very
	// likely doesn't exist: a dew point above the cricondenbar, say. It's
	// evidence rather than proof, but it separates "no such point" from
	// .no_convergence.
	outside_envelope
}

// FlashError is every failure in this module: a kind and a reason, never a
// panic and never a result that wasn't shown to be right (D7). Test for it
// with `if err is FlashError`.
pub struct FlashError {
	Error
pub:
	kind   Failure
	reason string
}

pub fn (e FlashError) msg() string {
	return 'flash: ${e.reason}'
}

fn fail(kind Failure, reason string) IError {
	return FlashError{
		kind:   kind
		reason: reason
	}
}
