package provide Xradiobox 1.0

# Create this to make sure there are registered in auto_mkindex
# these must come before the [widget create ...]
proc Xradiobox args {}
proc xradiobox args {}

widget create Xradiobox -type frame -base labelframe -components {
} -options {
    {-orient		orient		Orient		{horizontal}}
    {-labeltext		ALIAS labelframe -text}
}

namespace eval ::Widget::Xradiobox {;

;proc construct {w} {
    global tcl_platform
    variable $w
    upvar 0 $w data

    set data(tag_num) 0
    pack $data(labelframe) -side top -fill both -padx 5 -pady 5
    $data(labelframe) configure -text "test"
}

;proc configure {w args} {
    variable $w
    upvar 0 $w data

    foreach {key val} $args {
	set data($key) $val
	switch -- $key {
	    -labeltext {
		$data(labelframe) configure -text $val
	    }
	}
    }    
}

;proc _add {w tag args} {
    variable $w
    upvar 0 $w data

    set _arg(-text) "Tag $data(tag_num)"
    foreach {k v} $args {
	set _arg($k) $v
    }

    set r [radiobutton $data(labelframe).b_$data(tag_num) \
	       -text $_arg(-text) \
	       -variable [namespace current]::${w}(tag) \
	       -value $tag]

    if {[string match h* $data(-orient)]} {
	pack $r -side left -expand 1 -anchor w
    } else {
	pack $r -side top -expand 1 -anchor w
    }

    incr data(tag_num)
}

;proc _get {w} {
    variable $w
    upvar 0 $w data

    return $data(tag)
}

;proc _select {w tag} {
    variable $w
    upvar 0 $w data

    set data(tag) $tag
}

}; # end namespace eval ::Widget::Xradiobox

