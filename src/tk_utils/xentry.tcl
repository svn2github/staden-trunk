package provide Xentry 1.0

# Create this to make sure there are registered in auto_mkindex
# these must come before the [widget create ...]
proc Xentry args {}
proc xentry args {}

widget create Xentry -type frame -base entry -components {
    {checkbutton check check {-image ::Images::cross \
			      -indicatoron 0 \
			      -selectimage ::Images::tick \
			      -variable [namespace current]::${w}(Config)}}
    xlabel
} -options {
    {-label		label		Label		{}}
    {-config		config		Config		0}
    {-default		defaultText	DefaultText	0}
    {-state		state		State		normal}
    {-checkcommand	checkCommand	CheckCommand	{}}
    {-type		type		Type		{}}
    {-textvariable	textVariable	TextVariable	{}}
    {-bg		-background}
    {-background	background	Background	{#d9d9d9}}
    {-entrybg		ALIAS entry -background}
    {-fg		-foreground}
    {-foreground	ALIAS entry -foreground}
}

namespace eval ::Widget::Xentry {;

;proc construct {w} {
    global tcl_platform
    variable $w
    upvar 0 $w data

    pack $data(xlabel) -side left -fill both
    pack $data(entry) -side right

    bind $data(entry) <Return> "$w get"
    catch {bind $data(entry) <KP_Enter> "$w get"}

    # Hack for windows
    if {$tcl_platform(platform) == "windows"} {
	catch {[winfo parent $data(entry)] configure -bg SystemButtonFace}
    }
}

;proc configure {w args} {
    variable $w
    upvar 0 $w data

    set truth {^(1|yes|true|on)$}

    foreach {key val} $args {
	switch -- $key {
	    -label {
		$data(xlabel) configure -text $val
	    }
	    -default {
		$data(entry) delete 0 end
		$data(entry) insert 0 $val
	    }
	    -state {
		$data(xlabel) configure -state $val
		$data(entry) configure -state $val
	    }
	    -config {
		if {[set val [regexp $truth $val]]} {
		    catch {pack $data(check) -side left -before $data(xlabel)}
		} else {
		    catch {pack forget $data(check)}
		}
	    }
	    -textvariable {
		$data(entry) configure -textvariable $val
	    }
	    -background {
		$data(xlabel)     configure -bg $val
		#$data(base)      configure -bg $val
		$data(container) configure -bg $val
	    }
	    -entrybg {
		$data(entry)     configure -bg $val
	    }
	    -foreground {
		$data(xlabel) configure -fg $val
		$data(entry) configure -fg $val
	    }
	    -type {
		if {$val != ""} {
		    set data(-checkcommand) [namespace current]::check_$val
		} else {
		    set data(-checkcommand) ""
		}
	    }
	}
	set data($key) $val
    }
}

;proc _get {w args} {
    variable $w
    upvar 0 $w data

    if {$data(-checkcommand) != ""} {
	set comm [lindex $data(-checkcommand) 0]
	set cargs [lrange $data(-checkcommand) 1 end]
	if {[eval [list $comm] $w [list [$data(entry) get]] $cargs] == 0} {
	    raise [winfo toplevel $w]
	    focus $data(entry)
	    $data(entry) icursor end
	    return ""
	}
    }

    return [$data(entry) get]
}

;proc _get2 {w args} {
    variable $w
    upvar 0 $w data

    return [$data(entry) get]
}

;proc _check_config {w args} {
    variable $w
    upvar 0 $w data

    return $data(Config)
}

# -----------------------------------------------------------------------------
# Predefined type checking commands.

# -type int ?from ?to??
;proc check_int {w value {from -} {to -}} {
    if {[regexp {^[+-]?[0-9]+$} $value] == 0} {
	return 0
    }
    if {[string compare $from "-"]} {
	if {$value < $from} {
	    return 0
	}
    }
    if {[string compare $to "-"]} {
	if {$value > $to} {
	    return 0
	}
    }
    return 1
}

# -type float ?from ?to??
;proc check_float {w value {from -} {to -}} {
    set r {[+-]?([0-9]+|[0-9]+\.[0-9]*|[0-9]*\.[0-9]+)([Ee][+-]?[0-9]+)?$}
    if {[regexp $r $value] == 0} {
	return 0
    }
    if {[string compare $from "-"]} {
	if {$value < $from} {
	    return 0
	}
    }
    if {[string compare $to "-"]} {
	if {$value > $to} {
	    return 0
	}
    }
    return 1
}

# -type fileinput ?optional?
;proc check_fileinput {w path {optional 0} {multiple 0}} {
    if {$optional && $path == ""} {
	return 1
    } else {
        if {$multiple} {
	    set response 1
	    foreach f $path {
		if {[XCheckOpenFile [expandpath $f] $w] == 0} {
		    set response 0
		    break
		}
	    }
	    return $response
	} else {   
	    set response [XCheckOpenFile [expandpath $path] $w]
	    return $response
	}
    }
}

# -type fileoutput ?optional?
;proc check_fileoutput {w path {optional 0}} {
    set filename [expandpath $path]
    if {$optional && $filename == ""} {
	return 1;
    } else {
	set response [XCheckSaveFile "$filename" $w]
	if {$response && [file exists "$filename"]} {
	    DeleteFile $filename
	}
	return $response
    }
}

# -type directoryoutput ?optional?
;proc check_directoryoutput {w path {optional 0}} {
    set filename [expandpath $path]
    if {$optional && $filename == ""} {
	return 1;
    }
    if {[file exists $filename] && [file isdirectory $filename] == 0} {
	tk_messageBox -icon error -type ok -title "Not a directory" \
	    -message "'$filename' exists, but is not a directory."
	return 0;
    }
    if {[file exists $filename]} {
	set ret [tk_messageBox -icon warning -type yesnocancel \
		     -default no -title "Directory Exists" \
		     -message "Do you wish to overwrite $filename?"]
	if {$ret == yes} {
	    return 1
	} else {
	    return 0
	}
    }

    return 1
}

}; # end namespace eval ::Widget::Xentry

