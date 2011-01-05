# A wrapper around spinbox to add a -label parameter.

package provide Xspinbox 1.0

# Create this to make sure there are registered in auto_mkindex
# these must come before the [widget create ...]
proc Xspinbox args {}
proc xspinbox args {}

# Make combobox a top-level namespace item. The reason is that widget.tcl
# doesn't cope well when using -base on a widget within it's own namespace.

widget create Xspinbox -type frame -base spinbox -components {
    xlabel
} -options {
    {-label             label           Label           {}}
    {-state		state		State		normal}
    {-from              from            From            1}
    {-to                to              To              1000}
}

namespace eval ::Widget::Xspinbox {

;proc construct {w} {
    variable $w
    upvar 0 $w data

    pack $data(xlabel)  -side left -fill both
    pack $data(spinbox) -side right
}

;proc configure {w args} {
    variable $w
    upvar 0 $w data

    set alist {}

    foreach {key val} $args {
	switch -- $key {
	    -label {
                $data(xlabel) configure -text $val
	    }
	    -state {
		$data(xlabel) configure -state $val
		$data(spinbox) configure -state $val
	    }
	    -from {
		if {[$data(spinbox) cget -to] < $val} {
		    $data(spinbox) configure -to $val
		}
		$data(spinbox) configure -from $val
	    }
	    -to {
		if {[$data(spinbox) cget -from] > $val} {
		    $data(spinbox) configure -from $val
		}
		$data(spinbox) configure -to $val
	    }
	    default {
		lappend alist $key $val
	    }
	}
    }

    if {$alist != {}} {
	eval $data(spinbox) configure $alist
    }
}

}; # matches namespace eval

