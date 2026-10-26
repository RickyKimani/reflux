module components

// check_cas checks a CAS registry number: 2–7 digits, 2 digits and a check
// digit, as in 106-97-8, with the check digit right. The CAS number becomes
// a file name, so this also keeps anything path-like out. ChemSep marks two
// isomers it lists apart, para- and ortho-hydrogen, with a letter suffix
// (1333-74-0p and 1333-74-0o), which is allowed.
pub fn check_cas(cas string) ! {
	not_cas := error('`${cas}` is not a CAS number')
	parts := cas.split('-')
	if parts.len != 3 || parts[0].len < 2 || parts[0].len > 7 || parts[1].len != 2 {
		return not_cas
	}
	last := parts[2]
	suffix_ok := last.len == 1 || (last.len == 2 && last[1] >= `a` && last[1] <= `z`)
	digits := parts[0] + parts[1] + last#[..1]
	if !suffix_ok || !digits.bytes().all(it.is_digit()) {
		return not_cas
	}
	body := digits#[..-1].reverse()
	mut sum := 0
	for i, d in body {
		sum += (i + 1) * int(d - `0`)
	}
	if sum % 10 != int(digits[digits.len - 1] - `0`) {
		return error('CAS number `${cas}` fails its check digit')
	}
}
