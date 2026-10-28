// The benchmarks' claims that can be checked (exit criterion 8): reflux,
// thermo and DWSIM ran the same flashes and got the same answers. The
// timings depend on the machine, so they aren't asserted;
// tools/benchmarks.vsh writes them to docs/benchmarks.md.
import benchmarks
import math

fn test_the_three_programs_ran_the_same_flashes() {
	fx := benchmarks.load()!
	assert fx.cases.len == benchmarks.kinds.len * fx.systems.len
	for c in fx.cases {
		s := fx.systems[c.system]
		fl := benchmarks.fluid(s)!
		calls := benchmarks.prepare(fl, s, c)!
		assert calls.len == fx.states
		r := calls[0]()!
		what := '${s.name}, ${c.kind}'
		others := {
			'thermo':        c.thermo
			'DWSIM':         c.dwsim_default
			'DWSIM at 1e-8': c.dwsim_tight
		}
		for name, t in others {
			assert t.error == '', '${what}: ${name} failed: ${t.error}'
			// thermo's FlashVL stops about 3e-8 from the answer in β, and
			// 3e-9 K in a bubble point. DWSIM is held to criterion 6's
			// targets, β to 1e-4 and T to 0.01 K: its rounded constants put
			// its bubble points 2–3 mK low, and its default tolerance its β
			// up to 8.5e-5 off.
			if c.kind in ['pt2', 'pt1'] {
				tol := if name == 'thermo' { 1e-7 } else { 1e-4 }
				assert math.abs(t.answer.beta - r.beta) <= tol, '${what}: ${name}'
			} else {
				tol := if name == 'thermo' { 1e-7 } else { 0.01 }
				assert math.abs(t.answer.t - r.t) <= tol, '${what}: ${name}'
			}
		}
		if c.kind == 'ph' {
			// reflux's H at the state, which its PH flash starts from, is
			// DWSIM's, from the same constants and Cp. DWSIM's rounded
			// constants move the departure by up to 1e-4 of itself (step 5
			// found); measured, DWSIM's H is 0.7–1.7 J/mol above reflux's,
			// 9e-5 of H. (thermo's H has its own Cp, so it isn't compared.)
			h := benchmarks.molar_h(fl, r)
			for name, t in others {
				if name != 'thermo' {
					assert math.abs(t.answer.h - h) <= 2e-4 * math.abs(h), '${what}: ${name} H ${t.answer.h} vs ${h}'
				}
			}
		}
	}
}

fn test_every_state_is_the_flash_it_says() {
	// Each state's call is the flash at that state: PT at (t[k], p[k]), PH
	// back at t[k] from its own H there, the bubble point at p[k].
	fx := benchmarks.load()!
	for c in fx.cases {
		s := fx.systems[c.system]
		fl := benchmarks.fluid(s)!
		calls := benchmarks.prepare(fl, s, c)!
		for k, f in calls {
			r := f()!
			what := '${s.name}, ${c.kind}, state ${k}'
			assert r.p == c.p[k], what
			match c.kind {
				'pt2', 'pt1' {
					assert r.t == c.t[k], what
					assert r.phases.len == if c.kind == 'pt2' { 2 } else { 1 }, what
				}
				'ph' {
					assert math.abs(r.t - c.t[k]) <= 1e-9 * c.t[k], what
				}
				else {
					assert r.beta == 0.0, what
				}
			}
		}
	}
}
