module main

fn C.sysconf(name i32) i64

// cpus is the number of processors online, a study's workers by default
// (M8's D8). vlib's runtime counts them too, but tcc can't build it on
// Windows.
fn cpus() int {
	n := int(C.sysconf(C._SC_NPROCESSORS_ONLN))
	return if n > 0 { n } else { 1 }
}
