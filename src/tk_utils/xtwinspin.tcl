# A pair of xspinbox widgets

package provide Xtwinspin 1.0

# Create this to make sure there are registered in auto_mkindex
# these must come before the [widget create ...]
proc Xtwinspin args {}
proc xtwinspin args {}

# Make combobox a top-level namespace item. The reason is that widget.tcl
# doesn't cope well when using -base on a widget within it's own namespace.


widget create Xtwinspin -type frame -base frame -components {
    {xspinbox xspinbox1 xspinbox1 {-label {Start position}}}
    {xspinbox xspinbox2 xspinbox2 {-label {End position}}}
} -options {
    {-label1             label           Label           {}}
    {-label2             label           Label           {}}
    {-from               from            From            {}}
    {-to                 to              To              {}}
    {-state		 state		 State           normal}
    {-single             single          Single          0}
    {-bg                -background}
    {-background         background      Background      {}}
}

namespace eval ::Widget::Xtwinspin {

;proc construct {w} {
    variable $w
    upvar 0 $w data

    pack $data(xspinbox1) -fill both -expand 1
    pack $data(xspinbox2) -fill both -expand 1
}

;proc configure {w args} {
    variable $w
    upvar 0 $w data

    set alist {}
    foreach {key val} $args {
	switch -- $key {
	    -label1 {
                $data(xspinbox1) configure -label $val
	    }
	    -label2 {
                $data(xspinbox2) configure -label $val
	    }

	    -from {
                $data(xspinbox1) configure -from $val
                $data(xspinbox2) configure -from $val
		if {[$data(xspinbox1) get] > $val} {
		     $data(xspinbox1) set $val
		}
	    }

	    -to {
                $data(xspinbox1) configure -to $val
                $data(xspinbox2) configure -to $val
		if {[$data(xspinbox2) get] < $val} {
		     $data(xspinbox2) set $val
		}
	    }

	    -single {
		puts val=$val
		if {$val == 1 || $val == "yes"} {
		    puts "pack forget $data(xspinbox2)"
		    pack forget $data(xspinbox2)
		} else {
		    puts "pack $data(xspinbox2) -fill both -expand 1"
		    pack $data(xspinbox2) -fill both -expand 1
		}
	    }

	    default {
		lappend alist $key $val
	    }
	}
    }

    if {$alist != {}} {
	eval $data(xspinbox1) configure $alist
	eval $data(xspinbox2) configure $alist
    }
}

;proc _get {w} {
    variable $w
    upvar 0 $w data

    set start [$data(xspinbox1) get]
    set end   [$data(xspinbox2) get]

    if {$start > $end} {
	return {}
    } else {
	return [list $start $end]
    }
}

;proc _get_start {w} {
    variable $w
    upvar 0 $w data

    return [$data(xspinbox1) get]
}

;proc _get_end {w} {
    variable $w
    upvar 0 $w data

    return [$data(xspinbox2) get]
}

;proc _set {w start end} {
    variable $w
    upvar 0 $w data

    $data(xspinbox1) set $start
    $data(xspinbox2) set $end

    return ""
}

;proc _set_start {w val} {
    variable $w
    upvar 0 $w data

    $data(xspinbox1) set $val
}

;proc _set_end {w val} {
    variable $w
    upvar 0 $w data

    $data(xspinbox2) set $val
}

}; # matches namespace eval

