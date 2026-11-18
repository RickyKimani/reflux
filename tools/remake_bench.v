// remake_bench times what remaking a feed or a unit costs (M8 step 3, D4):
// for each unit type, Loaded.with on one of its numbers, against the unit's
// own calculation on its solved inlets, which is what a study's point then
// spends. Every M4 unit is in every_unit.toml or reference.toml; the
// reactors are a chain of each type.
//
//     v -prod -cc gcc -o remake_bench.exe tools/remake_bench.v
//     remake_bench.exe
module main

import flowsheet
import os
import reactors
import time

const chain = 'format = 2

[components]
names = ["methane", "water", "hydrogen", "carbon monoxide", "carbon dioxide"]

[model]
equation = "PR"

[feed.feed]
T = "900 K"
P = "5 bar"
flow = "10 mol/s"
mole_fractions = { methane = 0.25, water = 0.7, hydrogen = 0.05, "carbon monoxide" = 0, "carbon dioxide" = 0 }

[unit.convert]
type = "conversion_reactor"
in = "feed"
out = "s1"
conversions = [{ reaction = "reforming", key = "methane", conversion = 0.3 }]

[unit.tube]
type = "pfr"
in = "s1"
out = "s2"
volume = "0.5 m^3"
length = "2 m"
catalyst = { loading = "10 kg/m^3" }
pressure_drop = "0.2 bar"

[unit.tank]
type = "cstr"
in = "s2"
out = "s3"
volume = "0.2 m^3"
catalyst_loading = "10 kg/m^3"
T = "850 K"

[unit.shift]
type = "equilibrium_reactor"
in = "s3"
out = "s4"
reactions = ["shift"]
T = "600 K"

[unit.gibbs]
type = "gibbs_reactor"
in = "s4"
out = "product"
T = "1000 K"
'

fn us_since(t0 u64) f64 {
	return f64(time.sys_mono_now() - t0) / 1e3
}

fn main() {
	root := @VMODROOT
	kinetics := os.read_file(os.join_path(root, 'testdata', 'kinetics', 'xu_froment.toml'))!
	files := {
		'every_unit.toml': os.read_file(os.join_path(root, 'testdata', 'flowsheets', 'every_unit.toml'))!
		'reference.toml':  os.read_file(os.join_path(root, 'testdata', 'flowsheets', 'reference.toml'))!
		'chain':           chain + kinetics[kinetics.index('[site.Ni]') or { 0 }..]
	}
	println('${'unit':10} ${'type':20} ${'key':34} ${'remake, µs':>11} ${'calculate, µs':>14}')
	mut seen := []string{}
	for _, text in files {
		l := flowsheet.load_text(text, makers: reactors.makers())!
		sol := l.solve()!
		for i, p in l.sheet.units {
			kind := l.types[i]
			if kind in seen {
				continue
			}
			readings := l.readings_of('unit', p.unit.name)
			if readings.len == 0 {
				continue
			}
			seen << kind
			input := l.input(readings[0].key)!
			reps := 200
			mut t0 := time.sys_mono_now()
			for _ in 0 .. reps {
				l.with([input], [input.value])!
			}
			remake := us_since(t0) / f64(reps)
			t0 = time.sys_mono_now()
			for _ in 0 .. reps {
				p.unit.calculate(l.fluid, sol.units[i].inlets)!
			}
			calculate := us_since(t0) / f64(reps)
			println('${p.unit.name:10} ${kind:20} ${input.key:34} ${remake:11.1f} ${calculate:14.1f}')
		}
		if l.sheet.feeds.len > 0 {
			name := l.sheet.feeds[0].name
			readings := l.readings_of('feed', name)
			input := l.input(readings.filter(it.key.ends_with('.T') || it.key.ends_with('.P'))[0].key)!
			reps := 200
			t0 := time.sys_mono_now()
			for _ in 0 .. reps {
				l.with([input], [input.value])!
			}
			println('${name:10} ${'feed':20} ${input.key:34} ${us_since(t0) / f64(reps):11.1f} ${'(a flash)':>14}')
		}
	}
}
