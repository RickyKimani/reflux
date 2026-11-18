module main

import fit
import flowsheet
import time

// fit_command runs `reflux fit FILE [NAME] [--json OUT] [--csv OUT]` (M8's
// D14, D15, D17): a fit of the file's, named or its only one, printed as a
// table, with its JSON and CSV if asked for. It exits 0 with a fit, 1 when
// it couldn't run or found none, and 2 for a wrong command.
fn fit_command(args []string, mut c Console) int {
	if args.any(it == '--workers' || it.starts_with('--workers=')) {
		c.complain("reflux fit: --workers is study's and optimize's\n\n${usage}")
		return 2
	}
	a := parse_named('fit', args, mut c) or { return 2 }
	o := options(a.file_args) or {
		c.complain('reflux fit: ${err.msg()}\n\n${usage}')
		return 2
	}
	fits := flowsheet.load_fits(o.file) or {
		r := fit.load_failure(o.file, err)
		e := r.error or { flowsheet.ErrorReport{} }
		c.complain('reflux: ${located(o.file, e)}')
		if o.json != '' {
			write(o.json, r.json(), mut c)
		}
		return 1
	}
	k := pick('fit', o.file, fits.fits.map(it.name), a.name, mut c) or { return err.msg().int() }
	t0 := time.sys_mono_now()
	fitted := fit.run(fits, fits.fits[k])
	seconds := f64(time.sys_mono_now() - t0) / 1e9
	r := fitted.report(o.file, fits)
	taken := if seconds < 1 {
		'${flowsheet.figures(1e3 * seconds)} ms'
	} else {
		'${flowsheet.figures(seconds)} s'
	}
	c.say(r.table() + '\nIn ${taken}.')
	if o.json != '' && !write(o.json, r.json(), mut c) {
		return 1
	}
	if o.csv != '' && !write(o.csv, r.csv(), mut c) {
		return 1
	}
	return if fitted.ok { 0 } else { 1 }
}
