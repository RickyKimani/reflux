module flash

// Failure says why a flash or stability test failed.
pub enum Failure {
	bad_input      // a specification it can't use: a negative mole fraction, T ≤ 0, …
	no_convergence // the iteration ran out of steps without an answer it could prove
	no_volume_root // the model has no volume root at a state it needed
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
