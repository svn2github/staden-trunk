#
# Assert function. If "expr" is false we generate a Tcl error to be dealth
# with by the normal error handler.
#
proc assert {expr} {
    if {![uplevel 1 [list expr $expr]]} {
	return -code error "Assertion failed: $expr"
    }
} 