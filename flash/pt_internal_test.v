module flash

import math
import thermo

// split's own guards, which a successful flash never reaches: a start that
// converges to a negative-flash solution, or onto the feed itself, is a
// failure, not a split.

fn methane_butane() !thermo.Cubic {
	methane := thermo.Component{
		name:  'methane'
		tc:    190.6
		pc:    45.99e5
		omega: 0.012
		mw:    0.016043
	}
	butane := thermo.Component{
		name:  'n-butane'
		tc:    425.1
		pc:    37.96e5
		omega: 0.2
		mw:    0.058123
	}
	return thermo.new_cubic(.pr, [methane, butane], kij: [0.0, 0.02, 0.02, 0.0])!
}

fn test_a_split_outside_zero_to_one_is_refused() {
	// 90% methane at 300 K and 20 bar is a single vapour phase, but Wilson's
	// K straddle 1 there, and Rachford–Rice puts the feed beyond the
	// tie-line: β > 1, a negative flash. Substitution from there converges
	// to a "split" that doesn't contain the feed, which split refuses.
	m := methane_butane()!
	z := [0.9, 0.1]
	assert stability(m, 300.0, 20e5, z)!.stable
	sp := Splitter{
		m:       m
		t:       300.0
		p:       20e5
		z:       z
		present: [0, 1]
	}
	if s := sp.split(wilson_k(m.components(), 300.0, 20e5)) {
		assert false, 'a split of a stable feed, β = ${s.beta}'
	} else {
		msg := err.msg()
		assert msg.contains('outside (0, 1)') || msg.contains('collapsed') || msg.contains('two-phase range'), msg
	}
}

fn test_a_split_that_collapses_onto_the_feed_is_refused() {
	// Started next to the trivial solution, K = e^±0.001, substitution on
	// these unstable feeds slides onto x = y = z, where the fugacities are
	// equal too. Without the collapse check that would pass for a split of
	// two identical phases.
	m := methane_butane()!
	for c in [[0.7, 330.0, 131.5e5], [0.6, 340.0, 122e5], [0.5, 300.0, 50e5]] {
		z := [c[0], 1.0 - c[0]]
		assert !stability(m, c[1], c[2], z)!.stable
		sp := Splitter{
			m:       m
			t:       c[1]
			p:       c[2]
			z:       z
			present: [0, 1]
		}
		if s := sp.split([math.exp(1e-3), math.exp(-1e-3)]) {
			assert false, '${c}: a split of β ${s.beta} with x ${s.x} and y ${s.y}'
		} else {
			assert err.msg().contains('collapsed'), err.msg()
		}
	}
}
