#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#-----------------------------------------------------------------------------
# Input/Output File only dialogue
#-----------------------------------------------------------------------------

proc getFname {t title type {command {}} args} {
    frame $t

    if {"$type" == "load"} {
	set ty "CheckInput"
        button $t.browse -text "Browse" \
	    -command "InvokeFileBrowser $t.entry open;
                      focus \[entrybox_path $t.entry\]"
    } elseif {"$type" == "load_optional"} {
	set ty "CheckInputOptional"
        button $t.browse -text "Browse" \
	    -command "InvokeFileBrowser $t.entry open"
    } elseif {"$type" == "load_optional_default"} {
	set ty "CheckInputOptionalDefault $args"
        button $t.browse -text "Browse" \
	    -command "InvokeFileBrowser $t.entry open"
    } elseif {"$type" == "save_optional"} {
	set ty "CheckOutputOptional"
        button $t.browse -text "Browse" \
	    -command "InvokeFileBrowser $t.entry save"
    } else {
	set ty "CheckOutput"
        button $t.browse -text "Browse" \
	    -command "InvokeFileBrowser $t.entry save"
    }

    entrybox $t.entry \
	-title $title \
	-width 15 \
	-type $ty \
	-default $args

    raise $t.browse

    pack $t.entry -side left -expand 1 -fill x
    pack $t.browse -side right

    return $t.entry
}

#
# Returns the entrybox list component
#
proc getFname_in_name {f} {
    return [expandpath [entrybox_get $f.entry]]
}

#
# Returns the temp. list name for a given file name
#
proc FileListName {file} {
    return "#$file"
}

#
# Configures a getFname box
#
proc getFname_configure {path args} {
    eval entrybox_configure $path.entry $args
    eval $path.browse configure $args
}

#-----------------------------------------------------------------------------
# Input/Output List only dialogue
#-----------------------------------------------------------------------------

proc getLname {t title type} {
    frame $t

    entrybox $t.entry \
	-title $title \
	-width 15
    
    button $t.browse -text "Browse" -command "ListBrowse $t $t.entry $type"
    
    pack $t.entry -side left -expand 1 -fill x
    pack $t.browse -side right

    return $t.entry
}

#
# Returns the entrybox list component
#
proc getLname_in_name {f} {
    return [entrybox_get $f.entry]
}

#
# Configures a getLname box
#
proc getLname_configure {path args} {
    getFname_configure "$args"
}


#-----------------------------------------------------------------------------
# Input only List/File dialogue
#
# Apologies for the mess!
# In the past, the position on the list of radiobuttons was the key to knowing
# which button meant what. the rightmost (1) is list, 2nd from right is file,
# 3rd from right is 'single file', and anything else disables the browse
# button.
#
# This is inflexible, but the old syntax still holds. Now if we specify a
# command for the radiolist buttons, this command may end with "#type" where
# "type" is one of "list", "fofn", "selection" or "biolims". (Selection allows
# for multiple files to be selected, but works for just a single file too and
# so it replaces single mode.)
#-----------------------------------------------------------------------------

#
# The "Input list or file?" dialogues
#
proc lorf_in {f l cmds args} {
    eval frame $f $args
    global $f.Format
    
    set val 0
    foreach b [keylget l WHICH.BUTTONS] {
	set cmd [lindex $cmds $val]
	if {![regexp {.*#([^#]*)} $cmd dummy type]} {
	    set type [lindex {list fofn} $val]
	}
	lappend buts "$b -command {lorf_in_set $f [list $type] [list $cmd];
    		      global $f.Format; set $f.Format [list $type]}"
	incr val
    }

    # The name entrybox
    frame $f.name
    entrybox $f.name.entry \
	-title [keylget l NAME.NAME] \
	-default [keylget l NAME.VALUE] \
	-width 15
    button $f.name.browse \
	-text [keylget l NAME.BROWSE]

    # The format selector
    radiolist $f.format \
	-title [keylget l WHICH.NAME] \
	-buttons [set buts] \
	-default [keylget l WHICH.VALUE] \
	-orient horizontal

    pack $f.name.browse -side right -fill x
    pack $f.name.entry -side left -fill x -expand 1
    pack $f.format $f.name -side top -fill both -expand 1
    raise $f.name

    return $f
}

#
# configures an list of file box
#
proc lorf_in_configure {path args} {
    set in_arg 0
    set arglist ""
    set state ""

    # Process command line args
    foreach i $args {
	if {$in_arg} {
	     if {$option == "-state"} {
		set state "-state $i"
	     } else {
		lappend arglist $option $i
	     }
	     set in_arg 0
	} else {
	     set option $i
	     set in_arg 1
	}
    }

    eval radiolist_configure $path.format -state $i
    
    eval entrybox_configure $path.name.entry -state $i
    eval $path.name.browse configure -state $i
}


#
# Called from (or called to) set the radiolist component of the file/list
# selection dialogues. If command is set, this will be called.
#
# 'type' controls the operation of the Browse command. It may be one of:
#   selection		Direct multiple selection of files
#   fofn		A filebrowser to choose a file of filenames
#   list		A list browser to choose a list of filenames
#   other		Disables the entry box.
proc lorf_in_set {w type command} {
    set entry_w  $w.name.entry
    set browse_w $w.name.browse
    switch $type {
	"selection" {
	    $browse_w configure \
		    -state normal \
		    -command "InvokeFileBrowser $entry_w openmulti"
	    entrybox_configure $entry_w \
		    -state normal -type ""
	}
        "biolims" {
	    $browse_w configure \
		    -state normal \
		    -command "InvokeBiolimsBrowser $entry_w openmulti"
	    entrybox_configure $entry_w \
		    -state normal -type ""
	}
	"fofn" {
	    $browse_w configure \
		    -state normal \
		    -command "InvokeFileBrowser $entry_w open"
	    entrybox_configure $entry_w -type CheckInput \
		    -state normal
	}
	"list" {
	    $browse_w configure \
		    -state normal \
		    -command "ListBrowse $w $entry_w read"
	    entrybox_configure $entry_w -type CheckList \
		    -state normal
	}
	default {
	    $browse_w configure \
		    -state disabled
	    entrybox_configure $entry_w \
		    -state disabled
	}
    }

    if {$command != ""} {
	eval $command
    }
}


#
# Returns the radio list component
#
# Too many functions are now using this to make it easy to change. They expect
# a numerical number, rather than a proper string type. So we have to convert
# to the appropriate value.
#
proc lorf_in_get {w} {
    global $w.Format

    if {[set $w.Format] == "selection" || [set $w.Format] == "biolims"} {
	if {[string match "#*" [entrybox_get $w.name.entry]]} {
	    return 1; # List type
	} else {
	    if {[CheckInput $w.name.entry] == 0} {
		tkwait variable forever
	    }
	    return 3; # Single file
	}
    } else {
        return [radiolist_get $w.format]
    }
}


#
# Returns the entrybox list component
#
proc lorf_in_name {f} {
    return [expandpath [entrybox_get $f.name.entry]]
}

#
# Returns a list for a specific 'lorf' box. Loads the file to a list if
# appropriate.
# Uses temporary lists named "#name" for files.
#
proc lorf_get_list {f} {
    set name [lorf_in_name $f]

    global $f.Format
    if {[set form [set $f.Format]] == ""} {
	set form [lorf_in_get $f]
    }

    if {"$name" == ""} { return "" }

    # Find or generate a list name
    switch -regexp $form {
	list|2 {
	    set lname $name
	}
	fofn|1 {
	    set lname "#$name"
	    if {![ListWritable "$lname"]} {return ""}
	    if {[ListLoad $name "$lname"] == ""} {
		ListDelete $lname
		tk_messageBox \
			-icon error \
			-title "Empty list" \
			-message "'$name' contains zero names. Aborting." \
			-type ok \
			-parent $f
		#wait forever...
		set re_enter 0
		tkwait variable re_enter
		return ""
	    }
	}
	selection {
	    set lname $name
	}
	biolims {
	    set lname $name
	}
	default {
	    return ""
	}
    }

    # Load the list
    set items ""
    foreach i [ListToItems "$lname"] {
	if {"[file dirname $i]" == "[pwd]"} {
	    lappend items [file tail $i]
	} else {
	    lappend items $i
	}
    }
    if {$items == ""} {
	tk_messageBox \
		-icon error \
		-title "Empty list" \
		-message "'$name' contains zero names. Aborting." \
		-type ok \
		-parent $f
	#wait forever...
	set re_enter 0
	tkwait variable re_enter
	return ""
    }

    # For temporary lists - delete them
    if {$form != "list"} {
	ListDelete $lname
    }

    return $items
}

#-----------------------------------------------------------------------------
# Output only List/File dialogue
#-----------------------------------------------------------------------------

#
# The "Output list or file?" dialogues
#
proc lorf_out {f l cmds args} {
    eval frame $f $args
    
    set val 1
    foreach b [keylget l WHICH.BUTTONS] {
	lappend buts \
	    "$b -command {lorf_out_set $f $val \
		{[lindex $cmds [expr $val-1]]}}"
	incr val
    }

    # The name entrybox
    frame $f.name
    entrybox $f.name.entry \
	-title [keylget l NAME.NAME] \
	-default [keylget l NAME.VALUE] \
	-width 15
    button $f.name.browse \
	-text [keylget l NAME.BROWSE]

    # The format selector
    radiolist $f.format \
	-title [keylget l WHICH.NAME] \
	-buttons [set buts] \
	-default [keylget l WHICH.VALUE] \
	-orient horizontal

    pack $f.name.browse -side right -fill x
    pack $f.name.entry -side left -fill x -expand 1
    pack $f.format $f.name -side top -fill both -expand 1
    raise $f.name

    return $f
}


#
# Sets the radiolist component (invoking commands as appropriate)
#
# 1 for list
# 2 for file
# 3 for something else (disables the entrybox)
#
proc lorf_out_set {f value ocmd} {
    if {"$value" == 2} {
        $f.name.browse configure -state normal \
	    -command "InvokeFileBrowser $f.name.entry save"
	entrybox_configure $f.name.entry -state normal -type CheckOutput
    } elseif {"$value" == 1} {
	$f.name.browse configure -state normal \
	    -command "ListBrowse $f $f.name.entry create"
	entrybox_configure $f.name.entry -state normal \
	    -type CheckOutList
    } else {
	entrybox_configure $f.name.entry -state disabled
	$f.name.browse configure -state disabled
    }

    eval $ocmd
}


#
# Returns the radio list component
#
proc lorf_out_get {f} {
    return [radiolist_get $f.format]
}


#
# Returns the entrybox list component
# Uses temporary lists named "#name" for files.
# Also ensures that the list is clear before writing to it.
#
proc lorf_out_name {f} {

    set form [radiolist_get $f.format]
    set name [entrybox_get $f.name.entry]

    if {$name == ""} {
	return ""
    }

    if { "$form" == 2} { #file
	set lname "#$name"
    } else {
	set lname "$name"
    }

    #if [ListWritable $lname] {return $lname}
    ListCreate2 $lname ""
    return $lname

    #return ""
}

#
# Saves a list to a file. Also removes the temporary list name.
#
proc lorf_out_save {name} {
    scan $name #%s f

    ListSave $f $name
    ListDelete $name
}


#-----------------------------------------------------------------------------
# Misc utils
#-----------------------------------------------------------------------------
proc lorf_focus {f} {
    entrybox_focus $f.name.entry
}


