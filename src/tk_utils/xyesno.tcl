package provide Xyn 1.0

# Create this to make sure there are registered in auto_mkindex
# these must come before the [widget create ...]
proc Xyn args {}
proc xyn args {}


widget create Xyn -type frame -base frame -components {
    {checkbutton check check {-image ::Images::cross \
			      -indicatoron 0 \
			      -selectimage ::Images::tick \
			      -variable [namespace current]::${w}(Config)}}
    xlabel
    {radiobutton yes yes {-text Yes \
			  -variable [namespace current]::${w}(Var) \
			  -value 1}}
    {radiobutton no  no  {-text No \
			  -variable [namespace current]::${w}(Var) \
			  -value 0}}
} -options {
    {-bg		-background}
    {-background	ALIAS xlabel -background}
    {-fg		-foreground}
    {-foreground	ALIAS xlabel -foreground}
    {-config		config		Config		0}
    {-state		state		State		normal}
    {-label		label		Text		{}}
    {-orientation	orientation	Orientation	vertical}
    {-ycommand		yCommand	YCommand	{}}
    {-ncommand		nCommand	NCommand	{}}
    {-default		default		Default		{}}
    {-variable		variabel	Variable	{}}
}


namespace eval ::Widget::Xyn {;

;proc construct {w} {
    variable $w
    upvar 0 $w data

    pack $data(xlabel) -side left  -fill both
    pack $data(no)     -side right -fill both
    pack $data(yes)    -side right -fill both
    set data(Var) 1
    set data(VarName) ${w}(Var)
}

;proc configure {w args} {
    variable $w
    upvar 0 $w data

    set truth {^(1|yes|true|on)$}

    foreach {key val} $args {
	switch -- $key {
	    -background {
		$data(xlabel)    configure -bg $val
		$data(yes)       configure -bg $val
		$data(no)        configure -bg $val
		$data(base)      configure -bg $val
		$data(container) configure -bg $val
	    }
	    -foreground {
		$data(xlabel)    configure -fg $val
		$data(yes)       configure -fg $val
		$data(no)        configure -fg $val
	    }
	    -label {
		$data(xlabel)    configure -text $val
	    }
	    -orientation {
		catch {pack forget $data(xlabel) $data(yes) $data(no)}
		if {[string match h* $val]} {
		    pack $data(xlabel) -side left
	       	    pack $data(no) $data(yes) -side right
		} else {
		    pack $data(xlabel) -side top
	       	    pack $data(no) $data(yes) -side top -anchor w
		}
	    }
	    -state {
		$data(xlabel) configure -state $val
		$data(yes)    configure -state $val
		$data(no)     configure -state $val
	    }
	    -config {
		if {[set val [regexp $truth $val]]} {
		    catch {pack $data(check) -side left -before $data(xlabel)}
		} else {
		    catch {pack forget $data(check)}
		}
	    }
	    -ycommand {
		$data(yes) configure -command $val
	    }
	    -ncommand {
		$data(no) configure -command $val
	    }
	    -default {
		set $data(VarName) $val
	    }
	    -variable {
		$data(yes) configure -variable $val
		$data(no) configure -variable $val
		set data(VarName) $val
	    }
	}
	set data($key) $val
    }
}

;proc _get {w args} {
    variable $w
    upvar 0 $w data

    return [set $data(VarName)]
}

;proc _set {w val} {
    variable $w
    upvar 0 $w data

    set v [set $data(VarName) $val]
    if {$v == 0 && $data(-ncommand) != ""} {
	uplevel #0 $data(-ncommand)
    } elseif {$data(-ycommand) != ""} {
	uplevel #0 $data(-ycommand)
    }
}

;proc _check_config {w args} {
    variable $w
    upvar 0 $w data

    return $data(Config)
}

}; # end namespace eval ::Widget::Xyn
