module reactors

import flowsheet
import kinetics
import quantity

// The reactors' tables in a flowsheet file (D11), each a unit type with a
// maker, which the command line gives flowsheet.load (D1). Reactors need
// format 2, and every reactor but the Gibbs reactor the file's reactions;
// a reactor's `reactions` names those it runs, and a conversion reactor's
// `conversions` its steps in order. docs/flowsheet-format.md describes the
// keys.

// makers is the reactors' unit types.
pub fn makers() map[string]flowsheet.Maker {
	return {
		'conversion_reactor':  flowsheet.Maker(make_conversion)
		'equilibrium_reactor': make_equilibrium
		'gibbs_reactor':       make_gibbs
		'cstr':                make_cstr
		'pfr':                 make_pfr
	}
}

const volume = quantity.dim(0, 3, 0, 0, 0, 0)
const length = quantity.dim(0, 1, 0, 0, 0, 0)
const density = quantity.dim(1, -3, 0, 0, 0, 0)

// format_2 refuses a reactor in a format-1 file.
fn format_2(s flowsheet.Section) ! {
	if s.context.format < 2 {
		return s.fail('type', 'a reactor needs format 2 (M9), and this file is format ${s.context.format}',
			[])
	}
}

// reactions_of is the file's reactions, which a reactor needs.
fn reactions_of(s flowsheet.Section) !kinetics.Kinetics {
	return s.context.kinetics or {
		return s.fail('type', "a reactor of this type needs the file's reactions, and it has none",
			[])
	}
}

// named_reactions is the reactions a reactor runs, by index: those its
// `reactions` names, or every one that keep allows.
fn named_reactions(mut s flowsheet.Section, k kinetics.Kinetics, keep fn (int) bool) ![]int {
	names := k.reactions.map(it.name)
	if !s.has('reactions') {
		out := []int{len: names.len, init: index}.filter(keep(it))
		if out.len == 0 {
			return s.fail('', "none of the file's reactions can run here: name them with reactions",
				[])
		}
		return out
	}
	mut out := []int{}
	for name in s.strings('reactions')! {
		i := names.index(name)
		if i < 0 {
			return s.fail('reactions', 'no reaction ${name}', quantity.closest(name, names))
		}
		if i in out {
			return s.fail('reactions', '${name} is listed twice', [])
		}
		out << i
	}
	if out.len == 0 {
		return s.fail('reactions', 'no reactions', [])
	}
	return out
}

// Heat is a reactor's thermal mode, as its table gives it.
struct Heat {
	thermal Thermal
	t       f64
	duty    f64
}

// heat reads T, duty or adiabatic = true, at most one: with none, the
// reactor is isothermal at its inlet's T.
fn heat(mut s flowsheet.Section) !Heat {
	given := ['T', 'duty', 'adiabatic'].filter(s.has(it))
	if given.len > 1 {
		return s.fail(given[1], 'give only one of T, duty and adiabatic, not ${given.join(' and ')}',
			[])
	}
	if s.has('T') {
		return Heat{
			thermal: .isothermal
			t:       s.quantity('T', .temperature)!
		}
	}
	if s.has('duty') {
		return Heat{
			thermal: .duty
			duty:    s.quantity('duty', .power)!
		}
	}
	if s.has('adiabatic') {
		if !s.boolean('adiabatic')! {
			return s.fail('adiabatic', "leave adiabatic out for a reactor at its inlet's T", [])
		}
		return Heat{
			thermal: .adiabatic
		}
	}
	return Heat{}
}

// drop is a reactor's pressure drop, 0 if not given.
fn drop(mut s flowsheet.Section) !f64 {
	return if s.has('pressure_drop') {
		s.quantity('pressure_drop', .pressure_difference)!
	} else {
		0.0
	}
}

// species_index finds a species by any name it answers to, failing at the
// section's key.
fn species_index(species []kinetics.Species, s flowsheet.Section, key string, name string) !int {
	return kinetics.find(species, name) or { return s.fail(key, err.msg(), []) }
}

fn make_gibbs(mut s flowsheet.Section) !flowsheet.Unit {
	format_2(s)!
	h := heat(mut s)!
	mut inert := []int{}
	if s.has('inert') {
		for name in s.strings('inert')! {
			inert << species_index(s.context.species, s, 'inert', name)!
		}
	}
	return new_gibbs(s.name,
		species: s.context.species
		inert:   inert
		thermal: h.thermal
		t:       h.t
		duty:    h.duty
		drop:    drop(mut s)!
	)!
}

fn make_equilibrium(mut s flowsheet.Section) !flowsheet.Unit {
	format_2(s)!
	k := reactions_of(s)!
	which := named_reactions(mut s, k, fn (_ int) bool {
		return true
	})!
	h := heat(mut s)!
	return new_equilibrium_reactor(s.name,
		species:   k.species
		reactions: which.map(k.reactions[it])
		thermal:   h.thermal
		t:         h.t
		duty:      h.duty
		drop:      drop(mut s)!
	)!
}

fn make_conversion(mut s flowsheet.Section) !flowsheet.Unit {
	format_2(s)!
	k := reactions_of(s)!
	names := k.reactions.map(it.name)
	mut steps := []Conversion{}
	for mut t in s.tables('conversions')! {
		name := t.text('reaction')!
		i := names.index(name)
		if i < 0 {
			return t.fail('reaction', 'no reaction ${name}', quantity.closest(name, names))
		}
		// A number, a fraction such as "27.4 %", or an expression in T such
		// as "0.85 - 0.0018*(T in degC)", with the file's definitions.
		e := t.expression('conversion', k.ctx, &k.defs, quantity.Dim{}, 'a conversion')!
		if what := beyond_t(e, &k.defs) {
			return t.fail('conversion', 'a conversion may depend on T alone, and this uses ${what}',
				[])
		}
		fixed := e.nodes.all(it.kind !in [.var, .member, .def])
		steps << Conversion{
			reaction:   k.reactions[i]
			key:        species_index(k.species, t, 'key', t.text('key')!)!
			conversion: if fixed { e.eval(&k.defs, []f64{len: k.ctx.slots()}) } else { 0 }
			expression: if fixed { none } else { e }
		}
		t.finish()!
	}
	h := heat(mut s)!
	return new_conversion_reactor(s.name,
		species: k.species
		steps:   steps
		thermal: h.thermal
		t:       h.t
		duty:    h.duty
		drop:    drop(mut s)!
		ctx:     k.ctx
		defs:    k.defs
	)!
}

fn make_cstr(mut s flowsheet.Section) !flowsheet.Unit {
	format_2(s)!
	k := reactions_of(s)!
	which := named_reactions(mut s, k, fn [k] (i int) bool {
		return i in k.rated
	})!
	h := heat(mut s)!
	return new_cstr(s.name,
		k:         k
		reactions: which
		volume:    s.measure('volume', volume, "a reactor's volume")!
		loading:   if s.has('catalyst_loading') {
			s.measure('catalyst_loading', density, "a catalyst's loading")!
		} else {
			0.0
		}
		thermal:   h.thermal
		t:         h.t
		duty:      h.duty
		drop:      drop(mut s)!
	)!
}

fn make_pfr(mut s flowsheet.Section) !flowsheet.Unit {
	format_2(s)!
	k := reactions_of(s)!
	which := named_reactions(mut s, k, fn [k] (i int) bool {
		return i in k.rated
	})!
	h := heat(mut s)!
	mut cat := Catalyst{}
	if s.has('catalyst') {
		mut c := s.child('catalyst')!
		cat = Catalyst{
			loading:  c.measure('loading', density, "a catalyst's loading")!
			void:     if c.has('void_fraction') { c.fraction('void_fraction')! } else { 0.0 }
			diameter: if c.has('particle_diameter') {
				c.measure('particle_diameter', length, "a particle's diameter")!
			} else {
				0.0
			}
		}
		c.finish()!
	}
	mut mode := PressureDrop.no_drop
	mut given := 0.0
	mut visc := ?Viscosity(none)
	if s.has('pressure_drop') {
		if s.text('pressure_drop')! == 'ergun' {
			if !s.has('length') {
				return s.fail('pressure_drop', "Ergun's pressure drop needs the reactor's length",
					[])
			}
			mode = .ergun
			visc = new_wilke(s.context.records) or { return s.fail('pressure_drop', err.msg(), []) }
		} else {
			mode = .given
			given = s.quantity('pressure_drop', .pressure_difference)!
		}
	}
	integrator := if s.has('integrator') {
		match s.text('integrator')! {
			'rodas4' { Integrator.rodas4 }
			'dopri5' { .dopri5 }
			else {
				return s.fail('integrator', 'no integrator ${s.text('integrator')!}: rodas4 or dopri5', quantity.closest(s.text('integrator')!, [
					'rodas4',
					'dopri5',
				]))
			}
		}
	} else {
		Integrator.rodas4
	}
	return new_pfr(s.name,
		k:             k
		reactions:     which
		volume:        s.measure('volume', volume, "a reactor's volume")!
		length:        if s.has('length') {
			s.measure('length', length, "a reactor's length")!
		} else {
			0.0
		}
		catalyst:      cat
		thermal:       h.thermal
		t:             h.t
		duty:          h.duty
		drop:          mode
		pressure_drop: given
		viscosity:     visc
		integrator:    integrator
		rtol:          if s.has('rtol') { s.number('rtol')! } else { 1e-8 }
	)!
}
