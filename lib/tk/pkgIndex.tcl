if {[package vcompare [package provide Tcl] 8.4] != 0} { return }
package ifneeded Tk 8.4 "if {{$tcl_platform(os)} == {Darwin}} {[list load [file join $dir .. $env(MACHINE)-binaries libtk8.4.dylib] Tk]} else {[list load [file join $dir .. $env(MACHINE)-binaries libtk8.4.so] Tk]}"
