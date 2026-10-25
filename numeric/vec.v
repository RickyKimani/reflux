module numeric

// sum returns Σ xᵢ.
pub fn sum(x []f64) f64 {
	mut s := 0.0
	for v in x {
		s += v
	}
	return s
}

// dot returns Σ aᵢbᵢ. The arrays must be the same length.
pub fn dot(a []f64, b []f64) f64 {
	mut s := 0.0
	for i, v in a {
		s += v * b[i]
	}
	return s
}
