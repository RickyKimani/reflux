module kinetics

import quantity

// Exact arithmetic on quantity's ratios, for stoichiometry: a balance is
// exact or it isn't, and a suggested equation's coefficients are whole.
// Every operation refuses to go out of range rather than round.

type Ratio = quantity.Ratio

const zero = Ratio{0, 1}
const one = Ratio{1, 1}

fn whole(n int) Ratio {
	return Ratio{i64(n), 1}
}

fn sub(a Ratio, b Ratio) !Ratio {
	return a.add(b.neg())
}

fn div(a Ratio, b Ratio) !Ratio {
	return a.mul(quantity.ratio(b.den, b.num)!)
}

fn is_zero(a Ratio) bool {
	return a.num == 0
}

fn gcd(a i64, b i64) i64 {
	mut x, mut y := if a < 0 { -a } else { a }, if b < 0 { -b } else { b }
	for y != 0 {
		x, y = y, x % y
	}
	return x
}

// decimal writes a ratio as a decimal where it is one, 2.5, and as p/q
// where it isn't, 1/3.
fn decimal(a Ratio) string {
	mut d := a.den
	mut twos, mut fives := 0, 0
	for d % 2 == 0 {
		d /= 2
		twos++
	}
	for d % 5 == 0 {
		d /= 5
		fives++
	}
	places := if twos > fives { twos } else { fives }
	if d != 1 || places > 6 {
		return a.str()
	}
	mut scale := i64(1)
	for _ in 0 .. places {
		scale *= 10
	}
	n := a.num * (scale / a.den)
	m := if n < 0 { -n } else { n }
	sign := if n < 0 { '-' } else { '' }
	if places == 0 {
		return '${sign}${m}'
	}
	frac := '${m % scale}'
	return '${sign}${m / scale}.${'0'.repeat(places - frac.len)}${frac}'
}

// reduce brings rows to reduced row echelon form in place and returns the
// pivot column of each nonzero row, in order.
fn reduce(mut rows [][]Ratio) ![]int {
	mut pivots := []int{}
	if rows.len == 0 {
		return pivots
	}
	n := rows[0].len
	mut r := 0
	for c in 0 .. n {
		if r == rows.len {
			break
		}
		mut p := -1
		for i in r .. rows.len {
			if !is_zero(rows[i][c]) {
				p = i
				break
			}
		}
		if p < 0 {
			continue
		}
		rows[r], rows[p] = rows[p], rows[r]
		lead := rows[r][c]
		for j in 0 .. n {
			rows[r][j] = div(rows[r][j], lead)!
		}
		for i in 0 .. rows.len {
			if i == r || is_zero(rows[i][c]) {
				continue
			}
			f := rows[i][c]
			for j in 0 .. n {
				rows[i][j] = sub(rows[i][j], f.mul(rows[r][j])!)!
			}
		}
		pivots << c
		r++
	}
	return pivots
}

// null_space is a basis of the vectors v with a·v = 0, one for each column
// without a pivot, that column's entry 1.
fn null_space(a [][]Ratio, n int) ![][]Ratio {
	mut rows := [][]Ratio{len: a.len, init: a[index].clone()}
	pivots := reduce(mut rows)!
	mut basis := [][]Ratio{}
	for f in 0 .. n {
		if f in pivots {
			continue
		}
		mut v := []Ratio{len: n, init: zero}
		v[f] = one
		for i, p in pivots {
			v[p] = rows[i][f].neg()
		}
		basis << v
	}
	return basis
}

// smallest_whole scales a null-space vector to the smallest whole numbers in
// its direction. Its free column's entry is 1, so multiplying by the least
// common multiple l of its denominators leaves the whole numbers with no
// common factor: a prime dividing l divides some denominator qⱼ to the
// power it divides l, so it doesn't divide l/qⱼ, nor pⱼ, which is prime to
// qⱼ, nor, therefore, pⱼ·l/qⱼ.
fn smallest_whole(v []Ratio) ![]i64 {
	mut l := i64(1)
	for x in v {
		g := gcd(l, x.den)
		l = l / g * x.den
		if l > 1_000_000 {
			return error('out of range')
		}
	}
	mut out := []i64{len: v.len}
	for i, x in v {
		out[i] = x.num * (l / x.den)
	}
	return out
}
