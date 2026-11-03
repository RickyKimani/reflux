// reflux against DWSIM 9 (exit criteria 2 and 5), on testdata/m4_dwsim.json,
// written by hand on Windows by tools/gen_m4_dwsim.py: each unit alone from the
// same inlet, and the reference flowsheet on DWSIM's constants and kᵢⱼ
// (testdata/flowsheets/reference_dwsim.toml).
//
// Criterion 2's targets are T to 0.01 K, P to 0.01%, β and compositions to
// 1e-4, and duty and power to 0.01%. A deviation may exceed them only when the
// fixture reproduces DWSIM's number from something in DWSIM itself, as M3's
// criterion 6 was amended: here DWSIM's rounded constants, its pump's liquid
// density and its ideal mixing entropy. reflux must then match the exact model
// with DWSIM's constants. Every deviation has a bound that can't grow.
import flash
import flowsheet { EnergyStream, Outcome, Placed, Stream, Unit }
import json2
import math
import os
import thermo

struct DwStream {
	t         f64
	p         f64
	f         f64
	flows     []f64
	z         []f64
	beta      f64
	h         f64
	s         f64
	v         f64
	y         []f64
	x         []f64
	h_exact   f64
	h_rounded f64
}

struct DwSpec {
	kind          string @[json: 'type']
	mode          string
	value         f64
	pressure_drop f64
	efficiency    f64
	fractions     []f64
}

struct DwExact {
	t     f64
	h     f64
	beta  f64
	v     f64
	t_s   f64
	h_in  f64
	h_out f64
}

struct DwMixing {
	inlet  f64
	vapour f64
	liquid f64
}

struct DwUnit {
	name         string
	dwsim        string
	spec         DwSpec
	inlets       []DwStream
	outlets      []DwStream
	energy       f64
	from_streams f64
	exact        DwExact
	mixing_error DwMixing
}

struct DwComponent {
	name  string
	tc    f64
	pc    f64
	omega f64
	mw    f64
}

struct DwBalance {
	component f64
	energy    f64
	total     f64
}

struct DwReference {
	iterations int
	streams    map[string]DwStream
	energy     map[string]f64
	balance    DwBalance
}

struct DwOpenLoop {
	handed_back DwStream
	flows       f64
	t           f64
	p           f64
}

struct DwFixture {
	components []DwComponent
	kij        [][]f64
	units      []DwUnit
	reference  DwReference
	open_loop  DwOpenLoop
}

fn testdata(name string) string {
	return os.join_path(@VMODROOT, 'testdata', name)
}

const fixture = json2.decode[DwFixture](os.read_file(testdata('m4_dwsim.json')) or { panic(err) }) or {
	panic(err)
}

// loaded is the reference flowsheet on DWSIM's constants, whose fluid every
// case uses.
const loaded = flowsheet.load(testdata('flowsheets/reference_dwsim.toml')) or { panic(err) }

// DWSIM's gas constant, in its ideal mixing entropy.
const r_dwsim = 8.314

fn test_the_flowsheet_has_dwsims_constants() {
	// The generated component files and kᵢⱼ, as the fluid has them.
	comps := loaded.fluid.model.components()
	for i, c in fixture.components {
		assert comps[i].tc == c.tc && comps[i].pc == c.pc && comps[i].omega == c.omega && comps[i].mw == c.mw
	}
	m := loaded.fluid.model
	if m is thermo.Cubic {
		for i, row in fixture.kij {
			for j, k in row {
				assert m.kij[i * row.len + j] == k, '${i} ${j}'
			}
		}
		return
	}
	assert false, 'not a cubic'
}

fn heater_mode(mode string) flowsheet.HeaterMode {
	return match mode {
		'outlet_t' { .outlet_t }
		'duty' { .duty }
		'outlet_vapour_fraction' { .outlet_vapour_fraction }
		else { .delta_t }
	}
}

// unit_of is reflux's unit for a case, by its spec.
fn unit_of(c DwUnit) !Unit {
	s := c.spec
	if s.kind == 'heater' {
		return flowsheet.new_heater('unit',
			mode:          heater_mode(s.mode)
			value:         s.value
			pressure_drop: s.pressure_drop
		)!
	}
	if s.kind == 'valve' {
		return flowsheet.new_valve('unit', mode: .outlet_p, value: s.value)!
	}
	if s.kind == 'pump' {
		return flowsheet.new_pump('unit', mode: .outlet_p, value: s.value, efficiency: s.efficiency)!
	}
	if s.kind == 'compressor' {
		return flowsheet.new_compressor('unit',
			mode:       .outlet_p
			value:      s.value
			efficiency: s.efficiency
		)!
	}
	if s.kind == 'expander' {
		return flowsheet.new_expander('unit',
			mode:       .outlet_p
			value:      s.value
			efficiency: s.efficiency
		)!
	}
	if s.kind == 'mixer' {
		return flowsheet.new_mixer('unit', c.inlets.len)!
	}
	if s.kind == 'splitter' {
		return flowsheet.new_splitter('unit', s.fractions)!
	}
	if s.kind == 'flash_drum' {
		mode := if s.mode == 'duty' { flowsheet.DrumMode.duty } else { flowsheet.DrumMode.outlet_t }
		return flowsheet.new_flash_drum('unit', mode: mode, value: s.value)!
	}
	// DWSIM's separator leaves both outlets at the inlet's T and P.
	return flowsheet.new_component_separator('unit',
		fractions: s.fractions
		t:         c.outlets[0].t
		p:         c.outlets[0].p
	)!
}

// Deviation is how far a reflux stream is from DWSIM's: T in K and P as a
// fraction; the flows as a fraction of DWSIM's total; and β and the mole
// fractions, overall and in each phase, as differences, the worst of them.
struct Deviation {
	t           f64
	p           f64
	composition f64
}

fn worst(a []f64, b []f64) f64 {
	mut w := 0.0
	for i, x in a {
		w = math.max(w, math.abs(x - b[i]))
	}
	return w
}

fn deviation(s Stream, d DwStream) Deviation {
	mut c := math.max(worst(s.flows, d.flows) / d.f, math.abs(s.beta - d.beta))
	c = math.max(c, worst(s.result.z, d.z))
	if s.result.phases.len == 2 && d.y.len > 0 && d.x.len > 0 {
		c = math.max(c, worst(s.result.phases[0].state.n, d.y))
		c = math.max(c, worst(s.result.phases[1].state.n, d.x))
	}
	return Deviation{
		t:           math.abs(s.t - d.t)
		p:           math.abs(s.p - d.p) / d.p
		composition: c
	}
}

fn inlet_of(fl thermo.Fluid, d DwStream) !Stream {
	return flowsheet.feed_tp(fl, d.flows, d.t, d.p)!
}

fn energy_of(e []EnergyStream) f64 {
	mut q := 0.0
	for x in e {
		q += x.rate
	}
	return q
}

// relative is a duty's deviation from DWSIM's, as a fraction of DWSIM's, or of
// F·RT for a unit that should take none.
fn relative(q f64, dwsim f64, c DwUnit) f64 {
	mut flow := 0.0
	for s in c.inlets {
		flow += s.f
	}
	scale := if dwsim == 0.0 { flow * thermo.r_gas * c.inlets[0].t } else { math.abs(dwsim) }
	return math.abs(q - dwsim) / scale
}

// Run is a case as reflux calculates it.
struct Run {
	inlets  []Stream
	outcome Outcome
}

fn run(c DwUnit) !Run {
	fl := loaded.fluid
	inlets := c.inlets.map(inlet_of(fl, it) or { panic(err) })
	return Run{inlets, unit_of(c)!.calculate(fl, inlets)!}
}

// The bounds that can't grow, as measured on 2026-11-03: T in K, the worst of
// the flows, β and compositions, and the duty or work. Each within criterion 2's
// targets, but for the three cases a reproduction below explains.
const bounds = {
	'cooler to a T':                          [1e-9, 1.5e-5, 3.4e-5]
	'heater by its duty':                     [1.8e-3, 2.5e-5, 1e-9]
	'heater to a vapour fraction':            [5.8e-3, 1e-5, 1.3e-4]
	'heater by a rise, with a pressure drop': [1e-9, 1.6e-5, 4.6e-5]
	'valve':                                  [8.6e-5, 2.1e-5, 1e-9]
	'pump':                                   [0.3, 1e-9, 0.42]
	'compressor':                             [3.5e-4, 1e-9, 2.1e-6]
	'expander':                               [0.43, 2.8e-3, 2.6e-2]
	'mixer':                                  [1.3e-4, 2e-5, 1e-9]
	'splitter':                               [1e-9, 1.5e-5, 1e-9]
	'adiabatic drum':                         [1e-9, 1.6e-5, 1e-9]
	'drum at a T':                            [1e-9, 2.5e-5, 5.6e-5]
	'component separator':                    [1e-9, 3.3e-5, 7.8e-6]
}

const explained = ['heater to a vapour fraction', 'pump', 'expander']

fn test_each_unit_is_as_close_to_dwsim_as_it_was() {
	assert bounds.len == fixture.units.len
	for c in fixture.units {
		r := run(c)!
		bound := bounds[c.name] or { panic('no bound for ${c.name}') }
		mut worst_t, mut worst_c := 0.0, 0.0
		for k, s in r.outcome.outlets {
			d := deviation(s, c.outlets[k])
			assert d.p <= 1e-12, '${c.name}: P'
			worst_t = math.max(worst_t, d.t)
			worst_c = math.max(worst_c, d.composition)
		}
		// The separator's own number for its duty isn't its streams': that's
		// what goes in less what comes out, which the fixture's sign undoes.
		e := relative(energy_of(r.outcome.energy), c.energy, c)
		assert worst_t <= bound[0] && worst_c <= bound[1] && e <= bound[2], '${c.name}: T ${worst_t}, composition ${worst_c}, energy ${e}'
		if c.name !in explained {
			assert worst_t <= 0.01 && worst_c <= 1e-4 && e <= 1e-4, c.name
		}
	}
}

fn case(name string) DwUnit {
	return fixture.units.filter(it.name == name)[0]
}

// duty is F·H out less F·H in, by one of the fixture's enthalpies.
fn duty(c DwUnit, h fn (DwStream) f64) f64 {
	mut q := 0.0
	for s in c.outlets {
		q += s.f * h(s)
	}
	for s in c.inlets {
		q -= s.f * h(s)
	}
	return q
}

fn test_dwsims_duties_are_its_rounded_constants() {
	// Each of DWSIM's duties is its streams' enthalpies by the PR with its own
	// rounded Ωa, Ωb and R, to its convergence; the same states by the exact PR
	// give reflux's. That's the vapour-fraction heater's 1.3e-4, over the
	// target, and the rest of the duties' smaller differences.
	for c in fixture.units.filter(it.energy != 0.0 && it.name != 'pump' && it.name != 'expander') {
		rounded := duty(c, fn (s DwStream) f64 {
			return s.h_rounded
		})
		assert math.abs(rounded - c.from_streams) <= 1e-6 * math.abs(c.from_streams), '${c.name}: ${rounded} against ${c.from_streams}'
	}
	// Where reflux's states are DWSIM's, its duty is the exact PR's at them: one
	// outlet at a given T, or a separator's, whose compositions are given. A
	// drum's two outlets carry each program's own split.
	for name in ['cooler to a T', 'component separator'] {
		c := case(name)
		r := run(c)!
		exact := duty(c, fn (s DwStream) f64 {
			return s.h_exact
		})
		q := energy_of(r.outcome.energy)
		assert math.abs(q - exact) <= 1e-9 * math.abs(exact), '${name}: ${q} against ${exact}'
	}
	// At a vapour fraction the outlet's T moves too: reflux's is the exact PR's.
	c := case('heater to a vapour fraction')
	r := run(c)!
	out := r.outcome.outlets[0]
	assert math.abs(out.t - c.exact.t) <= 1e-7, '${out.t} against ${c.exact.t}'
	exact := c.outlets[0].f * c.exact.h - c.inlets[0].f * c.inlets[0].h_exact
	q := energy_of(r.outcome.energy)
	assert math.abs(q - exact) <= 1e-9 * math.abs(exact), '${q} against ${exact}'
}

fn test_dwsims_pump_is_its_own_density() {
	// DWSIM's pump work is F·v·ΔP/η with the v of its own reported density, not
	// its PR's. reflux's is the exact PR's v. Given DWSIM's work, reflux's energy
	// balance puts the outlet where DWSIM's is.
	c := case('pump')
	i, o := c.inlets[0], c.outlets[0]
	dp := c.spec.value - i.p
	dwsim := i.f * i.v * dp / c.spec.efficiency
	assert math.abs(dwsim - c.energy) <= 2e-5 * c.energy, '${dwsim} against ${c.energy}'
	r := run(c)!
	exact := i.f * c.exact.v * dp / c.spec.efficiency
	q := energy_of(r.outcome.energy)
	assert math.abs(q - exact) <= 1e-9 * exact, '${q} against ${exact}'
	inlet := r.inlets[0]
	back := flash.ph(loaded.fluid, c.spec.value, inlet.h + c.energy / i.f, inlet.flows,
		estimate: inlet.t
	)!
	assert math.abs(back.t - o.t) <= 1e-3, '${back.t} against ${o.t}'
}

// dwsim_mixing is DWSIM's ideal mixing entropy less reflux's, J/(mol·K): DWSIM
// has −R Σ xᵢ ln xᵢ · M/Mᵢ with R 8.314, where −R Σ xᵢ ln xᵢ belongs, M being
// the phase's molar mass (M3, docs/spike-dwsim.md).
fn dwsim_mixing(x []f64, mw []f64, r f64) f64 {
	mut m := 0.0
	for i, xi in x {
		m += xi * mw[i]
	}
	mut s := 0.0
	for i, xi in x {
		if xi > 0.0 {
			s += -r_dwsim * xi * math.log(xi) * m / mw[i] + r * xi * math.log(xi)
		}
	}
	return s
}

// dwsim_s is a state's molar S with DWSIM's mixing term in each phase.
fn dwsim_s(r flash.FlashResult, mw []f64) f64 {
	mut s := 0.0
	for ph in r.phases {
		s += ph.fraction * (loaded.fluid.props(ph.state).s + dwsim_mixing(ph.state.n, mw, thermo.r_gas))
	}
	return s
}

fn test_dwsims_expander_is_its_mixing_entropy() {
	// The outlet is two phases of different compositions, so DWSIM's wrong
	// mixing entropy no longer cancels: its isentropic outlet comes 0.142
	// J/(mol·K) off. reflux's isentropic path, by DWSIM's S, ends where DWSIM's
	// expander does.
	c := case('expander')
	i, o := c.inlets[0], c.outlets[0]
	mw := fixture.components.map(it.mw)
	// The formula is DWSIM's own term, as M3 measured it for these phases.
	for x, got in {
		'inlet':  [dwsim_mixing(i.z, mw, r_dwsim), c.mixing_error.inlet]
		'vapour': [dwsim_mixing(o.y, mw, r_dwsim), c.mixing_error.vapour]
		'liquid': [dwsim_mixing(o.x, mw, r_dwsim), c.mixing_error.liquid]
	} {
		assert math.abs(got[0] - got[1]) <= 1e-9, '${x}: ${got[0]} against ${got[1]}'
	}
	fl := loaded.fluid
	inlet := inlet_of(fl, i)!
	s_in := dwsim_s(inlet.result, mw)
	p := c.spec.value
	mut lo, mut hi := 250.0, i.t
	for hi - lo > 1e-9 {
		mid := 0.5 * (lo + hi)
		if dwsim_s(flash.pt(fl.model, mid, p, inlet.flows)!, mw) > s_in {
			hi = mid
		} else {
			lo = mid
		}
	}
	isentropic := flowsheet.feed_tp(fl, inlet.flows, 0.5 * (lo + hi), p)!
	h_out := inlet.h + c.spec.efficiency * (isentropic.h - inlet.h)
	out := flowsheet.feed_ph(fl, inlet.flows, p, h_out)!
	w := i.f * (h_out - inlet.h)
	// From 0.42 K, 2.7e-3 and 2.6% to these, inside the targets: what's left is
	// DWSIM's rounded constants.
	assert math.abs(out.t - o.t) <= 2.5e-3 && math.abs(out.beta - o.beta) <= 1.4e-5
	assert math.abs(w - c.energy) <= 5.7e-5 * math.abs(c.energy), '${(w - c.energy) / c.energy}'
	// And without it, reflux's expander is the exact PR's.
	r := run(c)!
	got := r.outcome.outlets[0]
	assert math.abs(got.t - c.exact.t) <= 1e-6 && math.abs(got.beta - c.exact.beta) <= 1e-9
	q := energy_of(r.outcome.energy)
	exact := i.f * (c.exact.h_out - c.exact.h_in)
	assert math.abs(q - exact) <= 1e-9 * math.abs(exact), '${q} against ${exact}'
}

// DwsimPump is reflux's pump with DWSIM's liquid volume, which DWSIM takes from
// its own density: the one difference in the loop beyond its rounded constants.
struct DwsimPump {
	name       string
	p_out      f64
	efficiency f64
	v          f64 // m³/mol
}

fn (p DwsimPump) inlet_ports() []string {
	return ['in']
}

fn (p DwsimPump) outlet_ports() []string {
	return ['out']
}

fn (p DwsimPump) calculate(fl thermo.Fluid, inlets []Stream) !Outcome {
	s := inlets[0]
	w := p.v * (p.p_out - s.p) / p.efficiency
	return Outcome{
		outlets: [flowsheet.feed_ph(fl, s.flows, p.p_out, s.h + w)!]
		energy:  [EnergyStream{'work', s.total() * w}]
	}
}

// The reference flowsheet's bounds: T in K and the composition's worst, for
// each stream. The recycle leaves the pump, whose work is DWSIM's density's.
const reference_bounds = {
	'feed':      [1e-9, 2.2e-5]
	'mixed':     [5.5e-3, 1.1e-5]
	'cold':      [1e-9, 1.6e-5]
	'gas':       [1e-9, 1.1e-5]
	'liquid':    [1e-9, 1.9e-5]
	'back':      [7.7e-5, 1.9e-5]
	'on':        [7.7e-5, 1.9e-5]
	'recycle':   [2.92e-2, 1.9e-5]
	'let down':  [7.8e-4, 2.4e-5]
	'flash gas': [7.8e-4, 4.5e-5]
	'product':   [7.8e-4, 2.4e-5]
}

fn test_the_reference_flowsheet_is_as_close_to_dwsim_as_it_was() {
	sol := loaded.solve()!
	ref := fixture.reference
	assert reference_bounds.len == ref.streams.len
	for name, d in ref.streams {
		dev := deviation(sol.streams[name], d)
		bound := reference_bounds[name] or { panic(name) }
		assert dev.p <= 1e-12 && dev.t <= bound[0] && dev.composition <= bound[1], '${name}: ${dev}'
		if name != 'recycle' {
			assert dev.t <= 0.01 && dev.composition <= 1e-4, name
		}
	}
	energy := sol.units.map(energy_of(it.energy))
	names := sol.units.map(it.name)
	cooler := energy[names.index('cooler')]
	assert math.abs(cooler - ref.energy['cooler']) <= 1.71e-4 * math.abs(ref.energy['cooler'])
	pump := energy[names.index('pump')]
	assert math.abs(pump - ref.energy['pump']) <= 0.43 * ref.energy['pump']
	// DWSIM's own balances, closed to its floor (docs/spike-dwsim-flowsheet.md).
	assert ref.balance.component <= 1e-7 && ref.balance.energy <= 1e-6 && ref.balance.total <= 1e-10
}

fn test_with_dwsims_pump_the_loop_is_dwsims() {
	// With DWSIM's liquid volume in its pump, reflux's loop is within criterion
	// 2's targets of DWSIM's everywhere, the recycle included. The cooler's duty
	// is then DWSIM's by DWSIM's rounded constants, and the exact PR's on reflux's.
	ref := fixture.reference
	mut sheet := flowsheet.Flowsheet{
		feeds: loaded.sheet.feeds
	}
	for p in loaded.sheet.units {
		if p.unit.name == 'pump' {
			sheet.units << Placed{
				...p
				unit: DwsimPump{
					name:       'pump'
					p_out:      31e5
					efficiency: 0.75
					v:          ref.streams['back'].v
				}
			}
		} else {
			sheet.units << p
		}
	}
	sol := sheet.solve(loaded.fluid, loaded.params)!
	for name, d in ref.streams {
		dev := deviation(sol.streams[name], d)
		// The recycle from 0.029 K to 1.7e-6 K; the rest as before.
		bound := if name == 'recycle' { 2e-6 } else { reference_bounds[name][0] }
		assert dev.t <= bound && dev.composition <= 1e-4, '${name}: ${dev}'
	}
	names := sol.units.map(it.name)
	pump := energy_of(sol.units[names.index('pump')].energy)
	assert math.abs(pump - ref.energy['pump']) <= 1e-4 * ref.energy['pump'], '${pump}'
	// The cooler's 1.7e-4 was mostly the pump's: now 3.3e-5, inside the target.
	// DWSIM's own number is its rounded constants', and reflux's the exact PR's at
	// DWSIM's states, to what the states' small differences allow.
	mixed, cold := ref.streams['mixed'], ref.streams['cold']
	rounded := cold.f * cold.h_rounded - mixed.f * mixed.h_rounded
	assert math.abs(rounded - ref.energy['cooler']) <= 1e-6 * math.abs(rounded)
	exact := cold.f * cold.h_exact - mixed.f * mixed.h_exact
	cooler := energy_of(sol.units[names.index('cooler')].energy)
	assert math.abs(cooler - exact) <= 2.1e-6 * math.abs(exact), '${(cooler - exact) / exact}'
	assert math.abs(cooler - ref.energy['cooler']) <= 3.4e-5 * math.abs(exact)
}

fn test_dwsim_hands_back_refluxs_recycle() {
	// The open loop (D12): DWSIM's own loop, fed reflux's converged recycle,
	// hands it back to within the same differences as the closed loop has, the
	// pump's 0.029 K among them.
	ol := fixture.open_loop
	assert ol.flows <= 1.5e-5 && math.abs(ol.t) <= 2.92e-2 && ol.p == 0.0, '${ol}'
}
