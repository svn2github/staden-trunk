package provide Xget_fname 1.0

# Create this to make sure there are registered in auto_mkindex
# these must come before the [widget create ...]
proc Xget_fname args {}
proc xget_fname args {}

widget create Xget_fname -type frame -base xentry -components {
    {button button button {-text Browse}}
} -options {
    {-bg		-background}
    {-background	ALIAS xentry -background}
    {-state		ALIAS xentry -state}
    {-type		type		Type		{load}}
    {-initialdir	directory	Directory	{}}
    {-filetypes		fileTypes	FileTypes	{}}
    {-text		text		Text		{}}
}

namespace eval ::Widget::Xget_fname {;
;proc construct {w} {
    variable $w
    upvar 0 $w data

    pack $data(xentry) -side left -fill both -expand 1
    pack $data(button) -side right -fill x

    set data(Cmd1) ""
    set data(Cmd2) ""
    set data(Cmd3) ""

    configure $w -type load
}

;proc configure {w args} {
    variable $w
    upvar 0 $w data

    set truth {^(1|yes|true|on)$}

    foreach {key val} $args {
	switch -- $key {
	    -background {
		$data(xentry)    configure -bg $val
		$data(button)    configure -bg $val
		$data(container) configure -bg $val
	    }
	    -type {
		switch -glob -- $val {
		    load {
			set data(Cmd1) \
			    "XInvokeFileBrowser $data(xentry) open"
			$data(xentry) configure -type fileinput
		    }
		    load_multiple {
			set data(Cmd1) \
			    "XInvokeFileBrowser $data(xentry) open_multiple"
			$data(xentry) configure -type "fileinput 0 1"
		    }
		    load_optional {
			set data(Cmd1) \
			    "XInvokeFileBrowser $data(xentry) open"
			$data(xentry) configure -type "fileinput 1"
		    }
		    save {
			set data(Cmd1) \
			    "XInvokeFileBrowser $data(xentry) save"
			$data(xentry) configure -type fileoutput
		    }
		    save_optional {
			set data(Cmd1) \
			    "XInvokeFileBrowser $data(xentry) save"
			$data(xentry) configure -type "fileoutput 1"
		    }
		}
		$data(button) configure \
		    -command "$data(Cmd1) $data(Cmd2) $data(Cmd3)"
	    }
	    -initialdir {
		set data(Cmd2) [list -initialdir $val]
		if {$data(Cmd1) != ""} {
		    $data(button) configure \
			-command "$data(Cmd1) $data(Cmd2) $data(Cmd3)"
		}
	    }
	    -filetypes {
		set data(Cmd3) [list -filetypes $val]
		if {$data(Cmd1) != ""} {
		    $data(button) configure \
			-command "$data(Cmd1) $data(Cmd2) $data(Cmd3)"
		}
	    }
	    -text {
		$data(xentry) configure -label $val
	    }
	    -state {
		$data(xentry) configure -state $val
		$data(button) configure -state $val
	    }
	}
	set data($key) $val
    }
}

}; # end namespace eval ::Widget::Xget_fname


proc get_fname {args} {
    return [eval xget_fname $args]
}