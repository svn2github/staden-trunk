package provide Xlistnotebook 1.0

# Create this to make sure there are registered in auto_mkindex
# these must come before the [widget create ...]
proc Xlistnotebook args {}
proc xlistnotebook args {}

widget create Xlistnotebook -type frame -base xmclistbox -components {
    {frame	frame		frame		{}}
} -options {
    {-raisecommand	raiseCommand	RaiseCommand	{}}
    {-bg		-background}
    {-background	background	Background	{#d9d9d9}}
    {-fg		-foreground}
    {-foreground	foreground	Foreground	{black}}
    {-borderwidth	borderwidth	Borderwidth	{}}
    {-bd		-borderwidth}
    {-relief		relief		Relief		flat}
}

namespace eval ::Widget::Xlistnotebook {;

;proc construct {w} {
    variable $w
    upvar 0 $w data

    set data(pane_num) 0
    set data(raised) ""

    pack $data(xmclistbox) -side left -fill both
    pack $data(frame) -side right -fill both -expand 1

    # We want this triggered after the standard mclistbox event
    bind $data(xmclistbox) <ButtonRelease-1> \
	"after idle {[namespace current]::select_pane $w}"
}

;proc configure {w args} {
    variable $w
    upvar 0 $w data

    set truth {^(1|yes|true|on)$}

    foreach {key val} $args {
	switch -- $key {
	    -background {
		$data(frame)     configure -bg $val
		$data(base)      configure -bg $val
		$data(container) configure -bg $val
	    }
	    -borderwidth {
		$data(frame)	configure -bd $val
	    }
	    -relief {
		$data(frame)	configure -relief $val
	    }
	}
	set data($key) $val
    }
}

# Add a new pane item.
#
# This puts the text into the listbox and creates a frame dedicated to
# this pane. The frame pathname is returned.
;proc _add {w name item args} {
    variable $w
    upvar 0 $w data
    
    set pnum [incr data(pane_num)]

    foreach {opt arg} $args {
	switch -- $opt {
	    -raisecommand {
		set data(raisecommand_$pnum) $arg
	    }
	}
    }

    $data(xmclistbox) insert end $item
    set pane $data(frame).pane_$pnum
    set data(translate_$name) $pnum
    set data(translate_$pnum) $pnum
    set data(back_translate_$pnum) $name

    frame $pane
    return $pane
}

;proc _raise {w {name {}}} {
    variable $w
    upvar 0 $w data

    if {$name == {}} {
	return $data(raised)
    } else {
	raise_pane $w $data(translate_$name)
    }
}

;proc _update {w row_name column text} {
    variable $w
    upvar 0 $w data

    $data(xmclistbox) update $data(translate_$row_name) $column $text
}

;proc _row_state {w name state} {
    variable $w
    upvar 0 $w data

    $data(xmclistbox) row_state $data(translate_$name) $state
}

;proc _subpanel {w name} {
    variable $w
    upvar 0 $w data
    
    set pnum $data(translate_$name)
    return $data(frame).pane_$pnum    
}

;proc _number_to_name {w num} {
    variable $w
    upvar 0 $w data

    return $data(back_translate_$num)
}

;proc raise_pane {w pane_num} {
    variable $w
    upvar 0 $w data

    # Remove old pane
    set slaves [pack slaves $data(frame)]
    pack forget $slaves

    # Pack new pane
    pack $data(frame).pane_$pane_num -fill both -expand 1
    set data(raised) $data(back_translate_$pane_num)

    $data(xmclistbox) select $pane_num

    # Callback -raisecommand procedure
    if {$data(-raisecommand) != ""} {
	uplevel #0 $data(-raisecommand)
    }
    if {[info exists data(raisecommand_$pane_num)]} {
	uplevel #0 $data(raisecommand_$pane_num)
    }
}

;proc select_pane {w} {
    variable $w
    upvar 0 $w data

    set sel [$data(xmclistbox) curselection]
    if {$sel == ""} {
	return
    }

    raise_pane $w $sel
}

}; # end namespace eval ::Widget::Xlistnotebook


