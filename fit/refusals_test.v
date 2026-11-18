// What a fit refuses (M8 steps 11 and 12, criterion 12): in the file,
// every key checked and every refusal at its line; in the data, each with
// its reason, the file and its line where it has one.
module fit

import flowsheet
import os

const scratch = os.join_path(os.temp_dir(), 'reflux_fit_test_${os.getpid()}')

fn testsuite_end() {
	os.rmdir_all(scratch) or {}
}

fn base_text() string {
	return os.read_file(os.join_path(@VMODROOT, 'testdata', 'm8_fit', 'ethylene.toml')) or {
		panic(err)
	}
}

fn vle_text() string {
	return os.read_file(os.join_path(@VMODROOT, 'testdata', 'm8_fit', 'ethane_propane.toml')) or {
		panic(err)
	}
}

// file_refused loads a fits file's text and expects its refusal at a line.
fn file_refused(text string, line int, says string) {
	flowsheet.load_fits_text(text, dir: os.join_path(@VMODROOT, 'testdata', 'm8_fit')) or {
		assert err is flowsheet.FlowsheetError, err.msg()
		e := err as flowsheet.FlowsheetError
		assert e.reason.contains(says), err.msg()
		assert e.line == line, 'line ${e.line}, not ${line}: ${err.msg()}'
		return
	}
	assert false, 'it loaded: ${says}'
}

// An empty list has no place in V's TOML tree, so its refusal is at its
// table's header.
fn test_the_file_refused() {
	text := base_text()
	at := 37 // parameters' line
	file_refused(text.replace('format = 3', 'format = 2'), 8, 'fits need format 3, not 2')
	file_refused(text.replace('parameters = ["k1.A", "k1.Ea", "K_E.A", "K_E.dH", "K_H", "k2.A", "k2.Ea"]',
		'parameters = []'), 35, 'a fit varies at least one parameter')
	file_refused(text.replace('"k2.Ea"]', '"k9.Ea"]'), at, 'no constant k9: the reactions have K_E, K_H, k1, k2')
	file_refused(text.replace('"k1.Ea"', '"k1.dH"'), at, "k1 is hydrogenation's rate constant, A·exp(−Ea/RT): fit k1.A or k1.Ea")
	file_refused(text.replace('"K_E.dH"', '"K_E.Ea"'), at, 'K_E is a term of site Ni, A·exp(−dH/RT): fit K_E.A or K_E.dH')
	file_refused(text.replace('"K_H", "k2.A"', '"K_H.A", "k2.A"'), at, 'a number that does not depend on T: fit K_H')
	file_refused(text.replace('"k2.Ea"]', '"k2.A"]'), at, 'k2.A is listed twice')
	file_refused(text.replace('"k2.Ea"]', '"kij[ethane, hydrogen]"]'), at, 'a fit is of rates or of kij, not both')
	file_refused(text.replace('{ name = "K_H", K = "0.3 bar^-1"', '{ name = "K_H", K = "0.3 bar^-1 * exp(500 K / T)"'),
		at, "written as an expression, which a fit can't take apart")
	file_refused(text + 'start = { "k3.A" = "1 mol/(kg*s*bar^2)" }\n', 39, 'k3.A is not one of the parameters')
	file_refused(text + 'start = { "k1.Ea" = "5 bar" }\n', 39, '"5 bar" is pressure, not')
	file_refused(text.replace('"r[hydrogenation]" = "5 %"', '"r[hydrogenation]" = "-5 %"'), 38,
		'must be above 0 %')
	file_refused(text.replace('"r[hydrogenation]" = "5 %"', '"r[hydrogenation]" = "five"'), 38,
		'neither a percentage nor a number with its unit')
	file_refused(text + 'method = "lm"\n', 39, 'no key method here')
	vle := vle_text()
	file_refused(vle.replace('"kij[ethane, propane]"', '"kij[ethane]"'), 18, 'a pair is written kij[a, b]')
	file_refused(vle.replace('"kij[ethane, propane]"', '"kij[ethane, ethane]"'), 18, 'a pair needs two components')
	file_refused(vle.replace('"kij[ethane, propane]"', '"kij[ethane, butane]"'), 18, 'butane is not one of the components')
	file_refused(vle.replace('[model]\nequation = "PR"\nkij = false\n', ''), 15, 'kij needs a [model]')
	file_refused(vle.replace('"kij[ethane, propane]"', '"k1.A"'), 18, 'the file has no reactions')
}

// data_refused writes data beside a fits file whose fit reads it, and
// expects the fit's refusal.
fn data_refused(text string, csv string, says string) {
	os.mkdir_all(scratch) or { panic(err) }
	os.write_file(os.join_path(scratch, 'data.csv'), csv) or { panic(err) }
	f := flowsheet.load_fits_text(text.replace('data = "ethylene.csv"', 'data = "data.csv"').replace('data = "ethane_propane.csv"',
		'data = "data.csv"'), dir: scratch) or { panic(err) }
	r := run(f, f.fits[0])
	assert !r.ok && r.error.contains(says), '${says}: ${r.error}'
}

const rate_head = 'T [K],P [bar],y[ethylene],y[hydrogen],y[ethane],y[methane],r[hydrogenation] [mol/(kg*s)],r[hydrogenolysis] [mol/(kg*s)]\n'

fn rate_rows(n int) string {
	mut out := ''
	for i in 0 .. n {
		t := 500 + 25 * (i % 5)
		out += '${t},${1 + i % 3},0.3,0.4,0.2,0.1,0.01,0.001\n'
	}
	return out
}

fn test_the_data_refused() {
	text := base_text()
	data_refused(text, 'P [bar],y[ethylene],r[hydrogenation] [mol/(kg*s)]\n2,0.5,0.01\n', 'rate data need columns T and P')
	data_refused(text, rate_head.replace('y[methane]', 'z[methane]') + rate_rows(10), 'no column z[methane] in rate data')
	data_refused(text, rate_head.replace('y[methane]', 'y[propane]') + rate_rows(10), 'column y[propane]')
	data_refused(text, rate_head + rate_rows(9) + '600,2,0.3,0.4,0.2,0.1,abc,0.001\n', 'data.csv:11, r[hydrogenation] [mol/(kg*s)]: "abc" is not a number')
	data_refused(text, rate_head + rate_rows(9) + '600,2,0.3\n', 'data.csv:11: it has 3 cells, not 8')
	data_refused(text, rate_head.replace('r[hydrogenation] [mol/(kg*s)]', 'r[hydrogenation] [mol/(m^3*s)]') +
		rate_rows(10), 'column r[hydrogenation] is ')
	data_refused(text, rate_head.replace('r[hydrogenolysis]', 'r[cracking]') + rate_rows(10), 'no reaction cracking with a rate')
	data_refused(text, rate_head + rate_rows(9) + '600,2,0.5,0.4,0.2,0.1,0.01,0.001\n', 'data.csv:11: the mole fractions add up to 1.2')
	data_refused(text, rate_head + rate_rows(9) + '600,2,0.3,0.4,0.2,0.1,0,0.001\n', 'data.csv:11: r[hydrogenation] is 0, which a relative weight makes an error of 0')
	data_refused(text.replace('"r[hydrogenation]" = "5 %"', '"r[hydrogenation]" = "0.1 bar"'), rate_head + rate_rows(10),
		'the weight of r[hydrogenation] is pressure')
	data_refused(text.replace('"r[hydrogenation]" = "5 %"', '"r[other]" = "5 %"'), rate_head + rate_rows(10),
		'a weight for r[other], which is not a measured rate')
	data_refused(text, rate_head + rate_rows(3), 'at least as many residuals (6) as parameters (7)')
	// At one T, an Arrhenius constant's A and Eₐ can't be told apart.
	mut one_t := ''
	for i in 0 .. 10 {
		one_t += '550,${1 + i % 4},${0.2 + 0.01 * i},0.4,0.2,0.1,${0.01 + 0.001 * i},${0.001 + 0.0001 * i}\n'
	}
	data_refused(text, rate_head + one_t, 'is not determined by the data')
	vle := vle_text()
	data_refused(vle.replace('weights = { P = "0.5 %", "y[ethane]" = "0.003" }', ''), 'T [K],P [Pa],x[ethane],y[ethane]\n250,1e6,0.5,0.7\n250,1.2e6,0.6,0.8\n',
		'a VLE fit needs the error of P')
	data_refused(vle.replace(', "y[ethane]" = "0.003"', ''), 'T [K],P [Pa],x[ethane],y[ethane]\n250,1e6,0.5,0.7\n250,1.2e6,0.6,0.8\n',
		'give y[ethane] its weight')
	data_refused(vle, 'T [K],P [Pa],x[methane],y[ethane]\n250,1e6,0.5,0.7\n', 'methane is not one of the components')
	data_refused(vle, 'T [K],P [Pa],x[ethane],y[ethane]\n250,1e6,1.5,0.7\n', 'x[ethane] must be a fraction')
}
