#
# Seqentry
# -----------------------------------------------------------------------------
# Implements a "query for sequence file" widget. It provides both
# local access and remote access (via EMBOSS).
#
# Public methods: 
#	get		Gets the name of a sequence selected.
#
#	sequence	Returns a chosen sequence, or "" if none found.
#


# -----------------------------------------------------------------------------
# seqentry class
# -----------------------------------------------------------------------------
class Seqentry {
    inherit itk::Widget

    constructor {args} {}

    public method get {}
    public method sequence {{format plain}}

    protected method _id_databases {}
    protected method _browse {}
    protected method _dbname_changed {args}

    variable _personal "Personal file"

    common _dbname
    common _ename
}

#
# Provide a lowercase access method for the class.
# 
proc ::seqentry {pathName args} {
    uplevel ::Seqentry $pathName $args
}

# -----------------------------------------------------------------------------
# Constructor
# -----------------------------------------------------------------------------
body Seqentry::constructor {args} {

    # database name selector
    itk_component add dbname {
	iwidgets::combobox $itk_component(hull).dbname \
	    -labeltext "Database" \
	    -labelpos n \
	    -textvariable [scope _dbname($this)] \
	    -editable false
    } {
	usual
    }
    pack $itk_component(dbname) -fill both -side left
    trace variable [scope _dbname] w [code $this _dbname_changed]

    # Filename / entryname box
    itk_component add entry {
	iwidgets::entryfield $itk_component(hull).entry \
	    -labelpos n \
	    -labeltext "Entry/Filename"
    } {
	usual
    }
    pack $itk_component(entry) -fill both -expand 1 -side left -padx 5

    # Browse button
    itk_component add browse {
	button $itk_component(hull).browse \
	    -text Browse \
	    -command [code $this _browse]
    }
    pack $itk_component(browse) -fill x -side right

    # Initialize the widget based on the command line options.
    eval itk_initialize $args

    # Read the databases
    foreach db [_id_databases] {
	$itk_component(dbname) insert list end $db
	set _ename($this,$db) ""
    }
    $itk_component(dbname) insert list 0 $_personal
    set _dbname($this) $_personal
    $itk_component(entry) delete 0 end

    #set ffgg [$itk_component(entry) cget -foreground]
    #[$itk_component(dbname) component entry] configure -disabledforeground $ffgg
    #set bbgg [$itk_component(entry) cget -textbackground]
    #[$itk_component(dbname) component entry] configure -disabledbackground $bbgg
}

# -----------------------------------------------------------------------------
# PUBLIC METHOD: get
#
# Returns a sequence USA (emboss format) or just plain filename
# -----------------------------------------------------------------------------
body Seqentry::get {} {
    if {$_dbname($this) == $_personal} {
	set prefix ""
    } else {
	set prefix $_dbname($this):
    }
    set e [$itk_component(entry) get]

    if {$e != ""} {
	return $prefix$e
    } else {
	return ""
    }
}


# -----------------------------------------------------------------------------
# PUBLIC METHOD: sequence
#
# Returns the actual sequence the user has selected.
# Format may be "entry" for the original database entry, or one of the
# supported EMBOSS formats. Defaults to "plain"
# -----------------------------------------------------------------------------
body Seqentry::sequence {{format plain}} {
    if {[set sname [get]] == ""} {
	return ""
    }
    set seq ""
    if {$format == "entry"} {
	if {$_dbname($this) == $_personal} {
	    # Read directly file disk, as entret seems to fail on some 
	    # disk files. This also means that we can read local files 
	    # without needing to have EMBOSS installed.
	    set sname_list $sname
	    set fd ""
	    foreach sname $sname_list {
		if {[catch "lappend fd [open $sname r]" err]} {
		    tk_messageBox \
			    -icon error \
			    -message "Failed to open file '$sname'.\
			    Error message is:\n\n$err" \
			    -parent $itk_component(hull)
		    return ""
		}
	    }
	} else {
	    if {[llength $sname] > 1} {
		tk_messageBox \
			-icon error \
			-message "Unable to open multiple entries."\
			-parent $itk_component(hull)
		return ""
	    }
	    set fd [open "|entret -outfile stdout $sname"]
	}
    } else {
	set fd [open "|seqret -outseq $format:stdout $sname"]
    }

    foreach i $fd {
	lappend seq [read $i]
    }
    
    #set seq [read $fd]
    global errorCode errorInfo
    set errorCode NONE
    catch {close $fd}
    if {$errorCode != "NONE"} {
	tk_messageBox \
	    -icon error \
	    -message "Failed to read file '$sname'\nError message is:\n\n\n \
		      $errorInfo" \
	    -parent $itk_component(hull)
	return ""
    }

    return $seq
}

# -----------------------------------------------------------------------------
# PRIVATE METHOD: _id_databases
#
# Uses EMBOSS showdb to determine which databases we can search by IDs
# Returns a list of database names
# -----------------------------------------------------------------------------
body Seqentry::_id_databases {} {
    if {[catch {set fd [open "|showdb -only -id" r]}]} {
	# Showdb command failed
	return
    }
    set dbs [read $fd]
    # Ignore stderr stuff.
    catch {close $fd}

    set r ""
    foreach {db status} $dbs {
	if {$status == "OK"} {
	    lappend r $db
	}
    }

    return $r
}

# -----------------------------------------------------------------------------
# PRIVATE METHOD: _browse
#
# Invokes a filebrowser for picking a sequence file. Only used for input
# type 'filename'.
# -----------------------------------------------------------------------------
body Seqentry::_browse {} {
    set f [tk_getOpenFile -parent $itk_component(hull) -multiple 65000]
    if {$f != ""} {
	set _ename($this,$_dbname($this)) $f
	$itk_component(entry) xview [string last / $f]
    }
}

# -----------------------------------------------------------------------------
# PRIVATE METHOD: _dbname_changed
#
# Trace callback. Invoked when a new dbname is picked
# -----------------------------------------------------------------------------
body Seqentry::_dbname_changed {args} {
    $itk_component(entry) configure \
	-textvariable [scope _ename($this,$_dbname($this))]
    if {$_dbname($this) == $_personal} {
	$itk_component(browse) configure -state normal
    } else {
	$itk_component(browse) configure -state disabled
    }
}
