package provide Xcheck 1.0

# Create this to make sure there are registered in auto_mkindex
# these must come before the [widget create ...]
proc Xcheck args {}
proc xcheck args {}


widget create Xcheck -type frame -base checkbutton -components {
    {checkbutton check check {-image ::Images::cross \
			      -indicatoron 0 \
			      -selectimage ::Images::tick \
			      -variable [namespace current]::${w}(Config)}}
} -options {
    {-bg		-background}
    {-background	ALIAS checkbutton -background}
    {-highlightbackground ALIAS checkbutton -highlightbackground}
    {-config		config		Config		0}
}


namespace eval ::Widget::Xcheck {;

;proc construct {w} {
    variable $w
    upvar 0 $w data

    pack $data(checkbutton) -side left -fill both
}

;proc configure {w args} {
    variable $w
    upvar 0 $w data

    set truth {^(1|yes|true|on)$}

    foreach {key val} $args {
	switch -- $key {
	    -config {
		if {[set val [regexp $truth $val]]} {
		    catch {pack $data(check) -side left -before $data(checkbutton)}
		} else {
		    catch {pack forget $data(check)}
		}
	    }
	    -background {
		$data(check)       configure -bg $val
		$data(checkbutton) configure -bg $val
		$data(container)   configure -bg $val
	    }
	}
	set data($key) $val
    }
}

;proc _check_config {w args} {
    variable $w
    upvar 0 $w data

    return $data(Config)
}

}; # end namespace eval ::Widget::Xcheck
