package provide Xscrolledtext 1.0

# Create this to make sure there are registered in auto_mkindex
# these must come before the [widget create ...]
proc Xscrolledtext args {}
proc xscrolledtext args {}

# Make combobox a top-level namespace item. The reason is that widget.tcl
# doesn't cope well when using -base on a widget within it's own namespace.

widget create Xscrolledtext -type frame -base text -components {
    {scrollbar yscrollbar y {-takefocus 0 -bd 1 -orient v \
				  -command [list $data(text) yview]}}
    {scrollbar xscrollbar x {-takefocus 0 -bd 1 -orient h \
				  -command [list $data(text) xview]}}
    label
} -options {
    {-label             label           Label           {}}
}

namespace eval ::Widget::Xscrolledtext {

;proc construct {w} {
    variable $w
    upvar 0 $w data
    variable x_disp 0
    variable y_disp 0

    $data(text) configure -xscrollcommand [list [namespace current]::xs $w]
    $data(text) configure -yscrollcommand [list [namespace current]::ys $w]

    grid columnconfig $w 0 -weight 1
    grid rowconfigure $w 1 -weight 1
    grid $data(text) -row 1 -sticky nsew
}

;proc configure {w args} {
    variable $w
    upvar 0 $w data

    foreach {key val} $args {
	switch -- $key {
	    -label {
                $data(label) configure -text $val
		if {$val != ""} {
		    grid $data(label) -row 0 -sticky nsew
		} else {
		    grid forget $data(label)
		}
	    }
	    * {
                eval $data(text) configure $args
	    }
	}
    }
}

# Callbacks when scrolling in X/Y. These are triggered when the data changes
# or also when the user scrolls. We check here to see if the scrollbar now
# needs displaying or hiding.
;proc xs {w start end} {
    variable $w
    variable x_disp
    upvar 0 $w data

    $data(xscrollbar) set $start $end

    if {$x_disp && $start == 0 && $end == 1} {
	grid forget $data(xscrollbar)
	set x_disp 0
    } elseif {!$x_disp && !($start == 0 && $end == 1)} {
	grid $data(xscrollbar) -row 2 -column 0 -sticky nsew
	set x_disp 1
    }
}

;proc ys {w start end} {
    variable $w
    variable y_disp
    upvar 0 $w data

    $data(yscrollbar) set $start $end

    if {$y_disp && $start == 0 && $end == 1} {
	grid forget $data(yscrollbar)
	set y_disp 0
    } elseif {!$y_disp && !($start == 0 && $end == 1)} {
	grid $data(yscrollbar) -row 1 -column 1 -sticky nsew
	set y_disp 1
    }
}

}; # matches namespace eval