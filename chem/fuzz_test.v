module chem

// Rng is a small seeded generator for the fuzzing: splitmix64.
struct Rng {
mut:
	state u64
}

fn (mut r Rng) below(n int) int {
	r.state += 0x9e3779b97f4a7c15
	mut z := r.state
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return int((z ^ (z >> 31)) % u64(n))
}

fn text(mut rng Rng, pieces []string) string {
	mut t := ''
	for _ in 0 .. 1 + rng.below(12) {
		t += pieces[rng.below(pieces.len)]
	}
	return t
}

// 20,000 seeded strings of a formula's pieces and stray ones: each is
// refused with a FormulaError with a place and a reason, or reads to a
// formula that its own Hill form reads back to.
fn test_random_formulas_never_panic() {
	pieces := ['C', 'H', 'O', 'N', 'Cl', 'Co', 'Hg', 'Si', 'U', '(', ')', '[', ']', '2', '10',
		'0', '.5', '999999999', 'x', 'q', '-', ' ', '·', '₂', 'Uuo']
	mut rng := Rng{20261110}
	mut accepted := 0
	for _ in 0 .. 20000 {
		t := text(mut rng, pieces)
		f := parse(t) or {
			assert err is FormulaError, '"${t}": ${err.msg()}'
			e := err as FormulaError
			assert e.column >= 1 && e.column <= t.runes().len + 1 && e.reason != '', '"${t}": ${err.msg()}'
			continue
		}
		assert parse(f.str())! == f, '"${t}" wrote ${f}'
		assert f.counts.all(it.n > 0 && it.z >= 1 && it.z <= elements)
		accepted++
	}
	assert accepted > 1000
}

// 20,000 seeded strings of SMILES pieces: each is refused with a place and
// a reason, or gives a composition with no negative count.
fn test_random_smiles_never_panic() {
	pieces := ['C', 'c', 'N', 'n', 'O', 'o', 'S', 's', 'P', 'Cl', 'Br', 'I', 'F', 'B', '[', ']',
		'H', 'H2', '@', '@@', '@TH1', '+', '-', '++', '+2', '=', '#', '$', '/', '\\', ':', ':1',
		'(', ')', '1', '2', '%12', '%1', '.', '[nH]', '[13C]', '[O-]', '[Fe+2]', '[se]', '*', 'x',
		' ', '0']
	mut rng := Rng{20261111}
	mut accepted := 0
	for _ in 0 .. 20000 {
		t := text(mut rng, pieces)
		s := smiles(t) or {
			assert err is FormulaError, '"${t}": ${err.msg()}'
			e := err as FormulaError
			assert e.column >= 1 && e.column <= t.runes().len + 1 && e.reason != '', '"${t}": ${err.msg()}'
			continue
		}
		assert s.formula.counts.all(it.n > 0)
		accepted++
	}
	assert accepted > 1000
}
