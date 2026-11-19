// `reflux kinetics check` and `reflux kinetics explain` (M9, D12):
//
//     reflux kinetics check FILE [T P NAME=y ...]
//     reflux kinetics explain FILE T P NAME=y ...
module main

import equilibrium
import flowsheet
import kinetics
import math
import quantity

// kinetics_command runs a kinetics subcommand, and returns its exit status.
// Its help is run's, from the command tree (commands.v).
fn kinetics_command(args []string, mut c Console) int {
	sub := args[0] or { '' }
	if sub !in ['check', 'explain'] {
		what := if sub == '' {
			'it needs check or explain'
		} else {
			'no kinetics command ${sub}${did_you_mean(sub, ['check', 'explain'])}'
		}
		c.complain('reflux kinetics: ${what}\n\n${kinetics_usage}')
		return 2
	}
	rest := args[1..]
	if rest.len == 0 || rest[0].starts_with('-') {
		c.complain('reflux kinetics ${sub}: it needs a file\n\n${kinetics_usage}')
		return 2
	}
	file := rest[0]
	if sub == 'explain' && rest.len < 4 {
		c.complain('reflux kinetics explain: it needs a state: T, P and the mole fractions\n\n${kinetics_usage}')
		return 2
	}
	if sub == 'check' && rest.len in [2, 3] {
		c.complain('reflux kinetics check: a state is T, P and the mole fractions\n\n${kinetics_usage}')
		return 2
	}
	loaded := flowsheet.load_kinetics(file) or {
		c.complain(kinetics_failure(file, err))
		return 1
	}
	k := loaded.kinetics
	mut state := ?kinetics.State(none)
	if rest.len > 1 {
		state = read_state(k, rest[1..]) or {
			c.complain('reflux kinetics ${sub}: ${err.msg()}')
			return 1
		}
	}
	if sub == 'check' {
		mut findings := k.check(state)
		if s := state {
			findings << equilibrium.check(&k, none, s)
		}
		rated := '${k.rates.len} rate${if k.rates.len == 1 { '' } else { 's' }}'
		c.say('${file}: ${k.reactions.len} reaction${if k.reactions.len == 1 { '' } else { 's' }}, ${rated}, ${findings.len} finding${if findings.len == 1 {
			''
		} else {
			's'
		}}')
		for f in findings {
			c.say(f.str())
		}
		return 0
	}
	s := state or { return 2 }
	explained := k.explain(s) or {
		c.complain('reflux kinetics explain: ${err.msg()}')
		return 1
	}
	c.say(explained.map(it.str()).join('\n\n'))
	return 0
}

// kinetics_failure is a file's failure as a compiler puts one: the file and
// line first, then the key and why.
fn kinetics_failure(file string, err IError) string {
	if err is flowsheet.FlowsheetError {
		place := if err.line > 0 { '${file}:${err.line}' } else { file }
		key := if err.key != '' && err.key != file { '${err.key}: ' } else { '' }
		hint := if err.suggestions.len > 0 {
			' (did you mean ${err.suggestions.join(', ')}?)'
		} else {
			''
		}
		return '${place}: ${key}${err.reason}${hint}'
	}
	return '${file}: ${err.msg()}'
}

// read_state reads T and P with their units and each species' mole
// fraction, NAME=y, by any name the species answers to; a species left out
// is 0, and the fractions must sum to 1 within 1e-6.
fn read_state(k kinetics.Kinetics, args []string) !kinetics.State {
	t := quantity.parse(args[0], .temperature)!
	p := quantity.parse(args[1], .pressure)!
	mut y := []f64{len: k.species.len}
	mut named := []bool{len: k.species.len}
	for a in args[2..] {
		if !a.contains('=') {
			return error('${a}: a mole fraction is NAME=y, as CH4=0.2')
		}
		name := a.all_before('=').trim_space()
		i := kinetics.find(k.species, name)!
		if named[i] {
			return error('${name}: ${k.species[i].name} is given twice')
		}
		v := a.all_after('=').trim_space()
		x := v.f64()
		if v == '' || !(x >= 0 && x <= 1) || !v.contains_only('0123456789.eE+-') {
			return error('${a}: a mole fraction is a number from 0 to 1')
		}
		y[i] = x
		named[i] = true
	}
	mut total := 0.0
	for x in y {
		total += x
	}
	if math.abs(total - 1) > 1e-6 {
		return error('the mole fractions sum to ${kinetics.sig(total)}, not 1')
	}
	return kinetics.State{
		t: t
		p: p
		y: y
	}
}
