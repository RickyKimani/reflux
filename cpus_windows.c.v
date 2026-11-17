module main

import os

// cpus is the number of processors, a study's workers by default (M8's
// D8), as Windows gives it every process. vlib's runtime counts them too,
// but tcc can't build it on Windows.
fn cpus() int {
	n := os.getenv('NUMBER_OF_PROCESSORS').int()
	return if n > 0 { n } else { 1 }
}
