# A modified label widget which supports the -state option. Useful for when
# building megawidgets so that when disabling the widget we can change the
# colour of the associated label.

proc xlabel {w args} {
    global tcl_version

    if {$tcl_version == "8.3"} {
	eval button $w -takefocus 0 -bd 0 -padx 1 -pady 1 -relief flat $args
	bindtags $w "$w Label . all"
	return $w
    } else {
	return [eval label $w $args]
    }
}
