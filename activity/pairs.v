module activity

import expr

// What the binary models (NRTL, UNIQUAC, Wilson) share: every pair of their
// components given exactly once.

fn expr_definitions() expr.Definitions {
	return expr.Definitions{}
}

fn pair_error(model string, reason string) IError {
	return ActivityError{
		model:  model
		reason: reason
	}
}

// check_pair refuses a pair that isn't two different components of names,
// or one given before, and marks it given.
fn check_pair(model string, names []string, i int, j int, mut given []bool) ! {
	nc := names.len
	if i < 0 || j < 0 || i >= nc || j >= nc {
		return pair_error(model, 'a pair (${i}, ${j}) of ${nc} components')
	}
	if i == j {
		return pair_error(model, '${names[i]} paired with itself')
	}
	if given[i * nc + j] {
		return pair_error(model, '${names[i]} with ${names[j]} given twice')
	}
	given[i * nc + j] = true
	given[j * nc + i] = true
}

// every_pair refuses a model missing a pair, naming them all. A caller
// that means a pair to have no energy of interaction says so, by giving it
// NRTL's τ = 0 or a τ or Λ of 1 (D7). That pair is ideal in NRTL and
// Wilson, and not in UNIQUAC, whose combinatorial term stays.
fn every_pair(model string, names []string, given []bool) ! {
	nc := names.len
	mut missing := []string{}
	for i in 0 .. nc {
		for j in i + 1 .. nc {
			if !given[i * nc + j] {
				missing << '${names[i]} with ${names[j]}'
			}
		}
	}
	if missing.len > 0 {
		return pair_error(model, 'no parameters for ${missing.join(', ')}')
	}
}
