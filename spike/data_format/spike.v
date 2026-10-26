// M2 step 1 spike: can V's toml module carry reflux's component files?
//
//     v run spike/data_format
//
// Prints what works and what doesn't. Not part of the simulator.
module main

import os
import strconv
import time
import toml

const here = os.dir(@FILE)
const root = os.join_path(here, '..', '..')

fn main() {
	read_the_hand_written_file()!
	missing_keys_and_types()!
	bad_input_gives_a_useful_error()
	float_parsing_matches_strconv()!
	float_printing_round_trips()!
	parse_speed()!
}

// constant reads a constant that is either a plain number or an inline
// table { value, source }, and says where it came from.
fn constant(doc toml.Doc, key string, default_source string) !(f64, string) {
	v := doc.value_opt('constants.${key}')!
	if v is map[string]toml.Any {
		return v['value'] or { return error('${key}: inline table without `value`') }.f64(), v['source'] or {
			toml.Any(default_source)
		}.string()
	}
	return v.f64(), default_source
}

fn read_the_hand_written_file() ! {
	println('== reading 106-97-8.toml')
	doc := toml.parse_file(os.join_path(here, '106-97-8.toml'))!
	src := doc.value('source').string()
	println('   name ${doc.value('name').string()}, cas ${doc.value('ids.cas').string()}, pubchem ${doc.value('ids.pubchem').int()}')
	println('   aliases ${doc.value('ids.aliases').array().as_strings()}')
	for key in ['tc', 'pc', 'omega', 'mw', 'dhf', 'charge'] {
		v, from := constant(doc, key, src)!
		println('   ${key:-7} = ${v:-14} from ${from}')
	}
	mc := doc.value('parameters.mathias_copeman').array().map(it.f64())
	println('   mathias_copeman ${mc}')
	corr := doc.value('correlations').as_map()
	println('   correlations, in file order: ${corr.keys()}')
	cp := corr['ideal_gas_cp'] or { return error('no ideal_gas_cp') }
	println('   ideal_gas_cp: eq ${cp.value('eq').int()}, coeffs ${cp.value('coeffs').array().map(it.f64())}, units ${cp.value('units').string()}')
	groups := doc.value('groups.unifac_vle').as_map()
	println('   unifac_vle groups: ${groups.keys()} -> ${groups.values().map(it.int())}')
}

fn missing_keys_and_types() ! {
	println('== missing keys and types')
	doc := toml.parse_file(os.join_path(here, '106-97-8.toml'))!
	missing := doc.value('constants.no_such_key')
	println('   value() of a missing key: ${missing.type_name()}')
	if _ := doc.value_opt('constants.no_such_key') {
		println('   value_opt() of a missing key: no error (BAD)')
	} else {
		println('   value_opt() of a missing key: error, as it should be')
	}
	println('   int `charge = 0` read as f64: ${doc.value('constants.charge').f64()} (type ${doc.value('constants.charge').type_name()})')
	println('   string read as f64: ${doc.value('name').f64()} (silently 0 if so: must check types)')
}

fn bad_input_gives_a_useful_error() {
	println('== malformed input')
	bad := 'format = 1\n[constants]\ntc = 425.12\npc = "unterminated\n'
	toml.parse_text(bad) or {
		println('   error: ${err.msg()}')
		return
	}
	println('   no error (BAD)')
}

// Every distinct number string in the vendored ChemSep XML, parsed by the
// TOML parser and by strconv. They must give the same f64.
fn float_parsing_matches_strconv() ! {
	println('== TOML float parsing vs strconv, on every ChemSep number')
	mut seen := map[string]bool{}
	for f in ['chemsep1.xml', 'chemsep2.xml'] {
		text := os.read_file(os.join_path(root, 'data', 'vendor', 'chemsep', f))!
		for part in text.split('value="')[1..] {
			s := part.all_before('"')
			if s.len > 0 && (s[0].is_digit() || s[0] == `-`) && s.count('-') <= 1 + s.count('E-')
				&& !s.contains_any('/ ,') {
				strconv.atof64(s) or { continue }
				seen[s] = true
			}
		}
	}
	keys := seen.keys()
	mut doc_text := []string{}
	for i, s in keys {
		// TOML needs a decimal point or exponent to make a float.
		t := if s.contains_any('.eE') { s } else { s + '.0' }
		doc_text << 'x${i} = ${t}'
	}
	doc := toml.parse_text(doc_text.join('\n'))!
	mut mismatches := 0
	for i, s in keys {
		want := strconv.atof64(s)!
		got := doc.value('x${i}').f64()
		if got != want {
			if mismatches < 5 {
				println('   mismatch: ${s} -> toml ${got}, strconv ${want}')
			}
			mismatches++
		}
	}
	println('   ${keys.len} distinct numbers, ${mismatches} mismatches')
}

// The writer must print floats so that (a) they parse back to the same bits
// and (b) TOML sees a float, not an integer.
fn float_printing_round_trips() ! {
	println('== printing floats for the writer')
	samples := [425.12, 3.796e6, 0.0581222, -125790.0, 9.85743e-6, 2.0, 0.0, 1e21, 1e-300, -3.35913e-7,
		0.1 + 0.2]
	mut bad := 0
	for x in samples {
		s := x.str()
		looks_float := s.contains_any('.eE') || s.contains('inf') || s.contains('nan')
		back := toml.parse_text('x = ${s}')!.value('x').f64()
		ok := looks_float && back == x
		if !ok {
			bad++
		}
		println('   ${x:-22} -> "${s}" ${if ok { 'ok' } else { 'PROBLEM' }}')
	}
	println('   ${bad} problems')
}

fn parse_speed() ! {
	println('== parse speed')
	text := os.read_file(os.join_path(here, '106-97-8.toml'))!
	sw := time.new_stopwatch()
	for _ in 0 .. 642 {
		_ := toml.parse_text(text)!
	}
	println('   642 parses of one component file: ${sw.elapsed().milliseconds()} ms')
}
