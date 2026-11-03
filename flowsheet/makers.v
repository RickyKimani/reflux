module flowsheet

// Maker builds a unit from its table in a file. Unit families of later
// milestones add theirs through LoadParams.makers, so that the loader needn't
// import them (D1). It reads the specification only: the ports are the
// loader's.
pub type Maker = fn (mut s Section) !Unit

// built_in_makers is the unit types of M4.
pub fn built_in_makers() map[string]Maker {
	return {
		'mixer':               Maker(make_mixer)
		'splitter':            make_splitter
		'valve':               make_valve
		'heater':              make_heater
		'pump':                make_pump
		'compressor':          make_compressor
		'expander':            make_expander
		'flash_drum':          make_flash_drum
		'component_separator': make_component_separator
	}
}

// pressure_spec reads which of outlet_P and delta_P is given.
fn pressure_spec(mut s Section) !(PressureMode, f64) {
	key := s.one_of(['outlet_P', 'delta_P'])!
	if key == 'outlet_P' {
		return PressureMode.outlet_p, s.quantity(key, .pressure)!
	}
	return PressureMode.delta_p, s.quantity(key, .pressure_difference)!
}

// machine_spec reads a pump's, compressor's or expander's specification.
fn machine_spec(mut s Section) !MachineParams {
	mode, value := pressure_spec(mut s)!
	return MachineParams{
		mode:       mode
		value:      value
		efficiency: s.fraction('efficiency')!
	}
}

fn make_mixer(mut s Section) !Unit {
	inlets := s.strings('in')!
	p := if s.has('P') { s.quantity('P', .pressure)! } else { 0.0 }
	return new_mixer(s.name, inlets.len, outlet_p: p)!
}

fn make_splitter(mut s Section) !Unit {
	return new_splitter(s.name, s.fractions('fractions')!)!
}

fn make_valve(mut s Section) !Unit {
	mode, value := pressure_spec(mut s)!
	return new_valve(s.name, mode: mode, value: value)!
}

fn make_heater(mut s Section) !Unit {
	key := s.one_of(['outlet_T', 'delta_T', 'duty', 'outlet_vapour_fraction'])!
	mode, value := match key {
		'outlet_T' { HeaterMode.outlet_t, s.quantity(key, .temperature)! }
		'delta_T' { HeaterMode.delta_t, s.quantity(key, .temperature_difference)! }
		'duty' { HeaterMode.duty, s.quantity(key, .power)! }
		else { HeaterMode.outlet_vapour_fraction, s.fraction(key)! }
	}
	drop := if s.has('pressure_drop') {
		s.quantity('pressure_drop', .pressure_difference)!
	} else {
		0.0
	}
	return new_heater(s.name, mode: mode, value: value, pressure_drop: drop)!
}

fn make_pump(mut s Section) !Unit {
	return new_pump(s.name, machine_spec(mut s)!)!
}

fn make_compressor(mut s Section) !Unit {
	return new_compressor(s.name, machine_spec(mut s)!)!
}

fn make_expander(mut s Section) !Unit {
	return new_expander(s.name, machine_spec(mut s)!)!
}

fn make_flash_drum(mut s Section) !Unit {
	inlets := s.strings('in')!
	key := s.one_of(['T', 'duty'])!
	mode, value := if key == 'T' {
		DrumMode.outlet_t, s.quantity(key, .temperature)!
	} else {
		DrumMode.duty, s.quantity(key, .power)!
	}
	p := if s.has('P') { s.quantity('P', .pressure)! } else { 0.0 }
	return new_flash_drum(s.name, inlets: inlets.len, mode: mode, value: value, outlet_p: p)!
}

fn make_component_separator(mut s Section) !Unit {
	return new_component_separator(s.name,
		fractions: s.by_component('fractions', .fraction)!
		t:         s.quantity('T', .temperature)!
		p:         s.quantity('P', .pressure)!
	)!
}
