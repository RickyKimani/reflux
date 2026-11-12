module flowsheet

import encoding.csv
import flash
import json2
import numeric
import thermo
import v.vmod

// The report's own format (D11), which docs/report-format.md describes. What
// reads a report can rely on it: it changes only with a new number. Format 2
// (M9) adds a reactor's profile to its unit.
pub const report_format = 2

// Report is a solve's results, or its failure: what the command line writes as
// JSON and CSV, and prints as a stream table. Every number is SI, and every
// list by component is in the order of components. A failure has no streams or
// units, since nothing unproven is reported (D8), but it keeps the history.
pub struct Report {
pub:
	format     int
	reflux     string // the version that wrote it
	file       string // the flowsheet, as the command line named it
	status     string // 'solved' or 'failed'
	error      ?ErrorReport @[json_null]
	components []ComponentReport
	streams    []StreamReport // the feeds, then each unit's outlets, unit by unit
	units      []UnitReport
	solver     ?SolverReport @[json_null]
	blocks     []BlockReport // in the order they were solved
	notes      []string      // what the plan says about its choices
	balance    ?Balance @[json_null] // of the whole flowsheet
	worst_unit ?Balance @[json_null] // the worst unit's, for each balance
}

pub struct ComponentReport {
pub:
	name       string // as the file names it
	compound   string // as the database does
	molar_mass f64    // kg/mol
}

// StreamReport is a stream: its state and flows, on both bases (D9).
pub struct StreamReport {
pub:
	name            string
	from            ?string @[json_null] // the unit it leaves, or none for a feed
	to              ?string @[json_null] // the unit it enters, or none for a product
	t               f64     @[json: 'T'] // K
	p               f64     @[json: 'P'] // Pa
	vapour_fraction f64
	molar_flow      f64 // mol/s
	mass_flow       f64 // kg/s
	volume_flow     f64 // m³/s, at its T and P
	h               f64 @[json: 'H'] // J/mol
	s               f64 @[json: 'S'] // J/(mol·K)
	v               f64 @[json: 'V'] // m³/mol
	molar_flows     []f64 // mol/s
	mass_flows      []f64 // kg/s
	// An empty stream keeps the composition it was given, so it has fractions
	// though its flows are 0.
	mole_fractions []f64
	mass_fractions []f64
	phases         []PhaseReport // the less dense first
}

// PhaseReport is one of a stream's phases, per mole of the phase.
pub struct PhaseReport {
pub:
	phase          string // 'vapour' or 'liquid', as the flash labels it
	fraction       f64    // of the stream's moles
	mole_fractions []f64
	h              f64 @[json: 'H'] // J/mol
	s              f64 @[json: 'S'] // J/(mol·K)
	v              f64 @[json: 'V'] // m³/mol
}

pub struct UnitReport {
pub:
	name    string
	kind    string @[json: 'type'] // as the file gives it
	inlets  []string // stream names, in port order
	outlets []string
	energy  []EnergyStream // W, into the unit's material
	// balance is on the inlets it was calculated from, which for a tear's
	// consumer are a pass behind the stream reported. A reactor's is on its
	// elements.
	balance Balance
	// profile is a reactor's state along it, if it has one (format 2).
	profile ?[]ProfileReport @[json_null]
}

// ProfileReport is a point of a reactor's profile.
pub struct ProfileReport {
pub:
	volume      f64 // m³ from the inlet
	length      f64 // m from the inlet, for a reactor of known length; else 0
	t           f64 @[json: 'T'] // K
	p           f64 @[json: 'P'] // Pa
	molar_flows []f64 // mol/s
}

pub struct SolverReport {
pub:
	method         Method
	tol            f64
	max_iterations int
}

// BlockReport is a block of the plan, and how its tears converged.
pub struct BlockReport {
pub:
	units  []string // in calculation order
	tears  []string
	cycles int
	// residuals has a row for each pass, with each tear's largest change as a
	// fraction of its scale (D7). It's empty for a block with no tears, which is
	// calculated once, and for one in which a failure came before a pass finished.
	residuals [][]f64
}

pub struct ErrorReport {
pub:
	kind        string // a Failure, or 'other' for an error from outside this module
	message     string // all of it, as one line
	unit        string
	reason      string
	cause       ?string @[json_null] // the flash's own failure, such as more_phases
	line        int
	key         string
	suggestions []string
}

// version is reflux's own, from v.mod.
pub fn version() string {
	m := vmod.decode(@VMOD_FILE) or { return 'unknown' }
	return m.version
}

// report is a solved flowsheet's report.
pub fn (l Loaded) report(sol Solution) Report {
	from, to := l.sheet.ends()
	mut streams := []StreamReport{}
	for name in l.sheet.stream_order() {
		streams << stream_report(l.fluid, name, sol.streams[name], end(from, name), end(to,
			name))
	}
	mut units := []UnitReport{}
	for i, p in l.sheet.units {
		units << UnitReport{
			name:    p.unit.name
			kind:    l.types[i] or { '' }
			inlets:  p.inlets
			outlets: p.outlets
			energy:  sol.units[i].energy
			balance: sol.units[i].balance()
			profile: profile_report(sol.units[i].profile)
		}
	}
	return Report{
		format:     report_format
		reflux:     version()
		status:     'solved'
		components: l.component_reports()
		streams:    streams
		units:      units
		solver:     l.solver_report()
		blocks:     block_reports(l.sheet, sol.plan, sol.history)
		notes:      sol.plan.notes
		balance:    sol.balance
		worst_unit: sol.worst
	}
}

// profile_report is a unit's profile as the report gives it, none for a unit
// without one.
fn profile_report(p ?Profile) ?[]ProfileReport {
	pr := p or { return none }
	return pr.points.map(ProfileReport{
		volume:      it.volume
		length:      it.length
		t:           it.t
		p:           it.p
		molar_flows: it.flows
	})
}

// failure is the report of a flowsheet that loaded and failed in solving: the
// error, and the history as far as it got.
pub fn (l Loaded) failure(err IError) Report {
	history := if err is FlowsheetError { err.history } else { []TearHistory{} }
	return Report{
		format:     report_format
		reflux:     version()
		status:     'failed'
		error:      error_report(err)
		components: l.component_reports()
		solver:     l.solver_report()
		blocks:     block_reports(l.sheet, l.plan, history)
		notes:      l.plan.notes
	}
}

// load_failure is the report of a flowsheet that didn't load.
pub fn load_failure(err IError) Report {
	return Report{
		format: report_format
		reflux: version()
		status: 'failed'
		error:  error_report(err)
	}
}

fn (l Loaded) component_reports() []ComponentReport {
	comps := l.fluid.model.components()
	mut out := []ComponentReport{cap: comps.len}
	for i, c in comps {
		out << ComponentReport{
			name:       l.names[i] or { c.name }
			compound:   c.name
			molar_mass: c.mw
		}
	}
	return out
}

fn (l Loaded) solver_report() SolverReport {
	return SolverReport{
		method:         l.params.method
		tol:            l.params.tol
		max_iterations: l.params.max_iterations
	}
}

// stream_order is every stream once: the feeds, then each unit's outlets, unit
// by unit, which is the order a reader follows the process in.
fn (f Flowsheet) stream_order() []string {
	mut out := f.feeds.map(it.name)
	for p in f.units {
		out << p.outlets
	}
	return out
}

// ends is the unit each stream leaves, and the unit each enters.
fn (f Flowsheet) ends() (map[string]string, map[string]string) {
	mut from, mut to := map[string]string{}, map[string]string{}
	for p in f.units {
		for name in p.outlets {
			from[name] = p.unit.name
		}
		for name in p.inlets {
			to[name] = p.unit.name
		}
	}
	return from, to
}

// end is the unit at a stream's end, or none where it's a feed or a product.
fn end(units map[string]string, stream string) ?string {
	return units[stream] or { return none }
}

fn stream_report(fl thermo.Fluid, name string, st Stream, from ?string, to ?string) StreamReport {
	comps := fl.model.components()
	z := st.result.z
	mut mass := []f64{len: comps.len}
	mut per_mole := 0.0 // kg/mol of the mixture
	for i, c in comps {
		mass[i] = st.flows[i] * c.mw
		per_mole += z[i] * c.mw
	}
	mut w := []f64{len: comps.len}
	for i, c in comps {
		w[i] = z[i] * c.mw / per_mole
	}
	return StreamReport{
		name:            name
		from:            from
		to:              to
		t:               st.t
		p:               st.p
		vapour_fraction: st.beta
		molar_flow:      st.total()
		mass_flow:       numeric.sum(mass)
		volume_flow:     st.total() * st.v
		h:               st.h
		s:               st.s
		v:               st.v
		molar_flows:     st.flows.clone()
		mass_flows:      mass
		mole_fractions:  z.clone()
		mass_fractions:  w
		phases:          st.result.phases.map(phase_report(fl, it))
	}
}

// phase_report is a phase, whose state is one mole of it, as of_flash has it.
fn phase_report(fl thermo.Fluid, ph flash.Phase) PhaseReport {
	pr := fl.props(ph.state)
	return PhaseReport{
		phase:          if ph.label == .vapor { 'vapour' } else { 'liquid' }
		fraction:       ph.fraction
		mole_fractions: ph.state.n.clone()
		h:              pr.h
		s:              pr.s
		v:              ph.state.v
	}
}

// block_reports is the plan's blocks, each with the history of its tears. The
// history has an entry for each block with tears that solve reached, in order.
fn block_reports(f Flowsheet, plan Plan, history []TearHistory) []BlockReport {
	mut out := []BlockReport{cap: plan.blocks.len}
	mut next := 0
	for b in plan.blocks {
		mut rows := [][]f64{}
		if b.tears.len > 0 && next < history.len {
			rows = history[next].residuals.clone()
			next++
		}
		out << BlockReport{
			units:     b.units.map(f.units[it].unit.name)
			tears:     b.tears
			cycles:    b.cycles
			residuals: rows
		}
	}
	return out
}

fn error_report(err IError) ErrorReport {
	if err is FlowsheetError {
		return ErrorReport{
			kind:        err.kind.str()
			message:     err.msg()
			unit:        err.unit
			reason:      err.reason
			cause:       if c := err.cause { c.str() } else { none }
			line:        err.line
			key:         err.key
			suggestions: err.suggestions
		}
	}
	return ErrorReport{
		kind:    'other'
		message: err.msg()
		reason:  err.msg()
	}
}

// json is the report as JSON, a value to a line, so that two reports diff well.
// json2 writes each number as the shortest decimal that reads back as the same
// double.
pub fn (r Report) json() string {
	return json2.encode(r, prettify: true, indent_string: '  ') + '\n'
}

// csv is the streams as CSV: a row for each, and a column for each of a
// stream's values in the JSON but its phases, named by its JSON key with its
// unit. Each number reads back as the same double.
pub fn (r Report) csv() string {
	mut head := ['stream', 'from', 'to', 'T [K]', 'P [Pa]', 'vapour_fraction', 'molar_flow [mol/s]',
		'mass_flow [kg/s]', 'volume_flow [m3/s]', 'H [J/mol]', 'S [J/(mol K)]', 'V [m3/mol]']
	for key in ['molar_flows', 'mass_flows', 'mole_fractions', 'mass_fractions'] {
		unit := match key {
			'molar_flows' { ' [mol/s]' }
			'mass_flows' { ' [kg/s]' }
			else { '' }
		}
		for c in r.components {
			head << '${key}.${c.name}${unit}'
		}
	}
	mut w := csv.new_writer()
	w.write(head) or {}
	for s in r.streams {
		mut row := [s.name, s.from or { '' }, s.to or { '' }]
		for x in [s.t, s.p, s.vapour_fraction, s.molar_flow, s.mass_flow, s.volume_flow, s.h, s.s,
			s.v] {
			row << x.str()
		}
		for values in [s.molar_flows, s.mass_flows, s.mole_fractions, s.mass_fractions] {
			for x in values {
				row << x.str()
			}
		}
		w.write(row) or {}
	}
	return w.str()
}
