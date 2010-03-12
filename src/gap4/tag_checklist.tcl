#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995-1997. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

# Catch as we need to use this in non-X mode.
catch {package require tablelist}
catch {namespace import tablelist::tablelist}

##############################################################################
#set the active tag defaults
proc TagDialog {name path command {version {}}} {
    global NGTag_$name$version

    if {[xtoplevel $path -resizable 1] == ""} {return}
    wm title $path "Active tags"

    set arr [array get NGTag_$name$version]
    tag_checklist $path.check $arr -selectmode extended -width 40 -height 20

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $path.but \
	-bd 2 \
	-relief groove \
	-ok_command "TagDialog_ok $path $name [list $command] \
		      [list $version]" \
	-asdefault_command "TagDialog_default $path $name" \
	-cancel_command "destroy $path $name" \
	-help_command "show_help gap4 {Conf-Tag}"
    
    pack $path.check -fill both -expand 1
    pack $path.but -fill both
}

proc TagDialog_ok {path name command version} {
    global NGTag_$name$version

    array set NGTag_$name$version [$path.check list2array]
    destroy $path

    if {$command != {}} {
	eval $command
    }
}    

proc TagDialog_default {path name} {
    global NGTag gap_defs env

    array set tmp [$path.check list2array]
    set active_tags ""
    for {set i 0} {$i < $NGTag(num_tags)} {incr i} {
        if {$tmp($i)} {
            lappend active_tags $NGTag($i,tagid)
        }
    }

    keylset gap_defs $name $active_tags
    update_defs gap_defs $env(HOME)/.gaprc $name
}

##############################################################################
#sets up a local version of the active tag list
proc SetDefaultTags {name {version {}}} {
    global gap_defs NGTag_$name$version NGTag_num_tags${name}${version} NGTag
    
    # Initialise all tags to off by default
    set list ""
    for {set i 0} {$i < $NGTag(num_tags)} {incr i} {
	set NGTag_${name}${version}($i) 0
	lappend list $NGTag($i,tagid)
    }

    set NGTag_num_tags${name}${version} $NGTag(num_tags)

    foreach i [keylget gap_defs $name] {
	if {$i == {*}} {
	    # Add all tags
	    for {set j 0} {$j < $NGTag(num_tags)} {incr j} {
		set NGTag_${name}${version}($j) 1
	    }
	} else {
	    if {[string match -* $i]} {
		set c [string trimleft $i -]
		if {[set ind [lsearch -exact $list $c]] != -1} {
		    set NGTag_${name}${version}($ind) 0
		}
	    } else {
		set c [string trimleft $i +]
		if {[set ind [lsearch -exact $list $c]] != -1} {
		    set NGTag_${name}${version}($ind) 1
		}
	    }
	}
    }
}


##############################################################################
#Gets the default tag list
proc GetDefaultTags {name {version {}}} {
    global NGTag NGTag_$name$version NGTag_num_tags${name}${version} gap_defs

    set active_tags ""

    for {set i 0} {$i < [set NGTag_num_tags${name}${version}]} {incr i} {
        if {[set NGTag_${name}${version}($i)]} {
            lappend active_tags $NGTag($i,tagid)
        }
    }

    return $active_tags
}


##############################################################################
#
# tag_checklist:
#
# Creates a tag checklist frame containing a listbox of all tags types.
# Extra command line arguments are as per the frame widget.
#
proc tag_checklist {path array args} {
    global NGTag
    global gap_defs

    # Create the frame
    eval frame $path -class TagCheckList

    rename $path $path.widget
    proc $path {command args} {
	set path [lindex [info level [info level]] 0]
	if {[info procs tag_checklist_$command] != ""} {
	    return [eval tag_checklist_$command $path $args]
	} else {
	    return [eval $path.list [list $command] $args]
	}
    }

    bind $path <Destroy> {
	catch {rename $path {}}
    }

    # Create and populate our tablelist
    set l [eval [list tablelist $path.list \
			       -columns {0 "Index"
				   0 "Code"
				   0 "Tag name"} \
			       -labelcommand tablelist::sortByColumn \
			       -exportselection 0 \
			       -stretch 1 \
			       -yscrollcommand [list $path.yscroll set]
			  ] \
	       $args]
    $l columnconfigure 0 -sortmode integer
    for {set i 0} {$i < $NGTag(num_tags)} {incr i} {
	$l insert end [list $i $NGTag($i,tagid) $NGTag($i,tagtype)]
    }

    tag_checklist_array2list $l $array

    # Add a scrollbar    
    scrollbar $path.yscroll -command "$path.list yview"

    grid columnconfigure $path 0 -weight 1
    grid rowconfigure $path 0 -weight 1

    grid $path.list $path.yscroll -sticky nsew

    set mode [$l cget -selectmode]
    if {$mode == "extended" || $mode == "extended2" || $mode == "multiple"} {
	# Add select all and deselect all messages.
	frame $path.buttons
	button $path.buttons.select_all \
	    -text "Select all" \
	    -command "$path.list selection clear 0 end;
		     $path.list selection set 0 end"
	button $path.buttons.clear_all \
	    -text "Clear all" \
	    -command "$path.list selection clear 0 end"
	pack $path.buttons.select_all $path.buttons.clear_all \
	    -side left -expand 1
	grid $path.buttons -columnspan 2 -sticky nsew
    }

    return $path
}

proc tag_checklist_component {w comp} {
    switch $comp {
	frame {
	    return $w.widget
	}
	list {
	    return $w.list
	}
	default {
	    return -code error "No component '$comp'"
	}
    }
}

proc tag_checklist_selection {w} {
    set result ""
    foreach item [$w.list curselection] {
	lappend result [$w.list get $item]
    }

    return $result
}

proc tag_checklist_list2array {w} {
    set w $w.list

    set nitems [$w index end]
    set arr ""

    for {set i 0} {$i < $nitems} {incr i} {
	set v [$w get $i]
	set used([lindex $v 0]) [$w selection includes $i]
    }

    for {set i 0} {$i < $nitems} {incr i} {
	lappend arr $i $used($i)
    }

    return $arr
}

proc tag_checklist_array2list {w arr} {
    $w selection clear 0 end
    foreach {ind val} $arr {
	if {$val} {
	    $w selection set $ind
	}
    }
}


##############################################################################
#convert C tag ids and types into tcl ids and types
proc InitTagArray { } {
    global NGTag

    set taglist [get_tag_array]
    set NGTag(num_tags) [llength $taglist]
    for {set i 0}  {$i < $NGTag(num_tags)} {incr i} {
        set t [lindex $taglist $i]
	set NGTag($i,tagtype) [lindex $t 0]
	set NGTag($i,tagid)   [lindex $t 1]
	set NGTag($i,comment) [lindex $t 2]
    }
}

##############################################################################
# Convert a tag type into an ID
proc tag_type2id {type {default 0}} {
    global NGTag

    for {set i 0}  {$i < $NGTag(num_tags)} {incr i} {
	if {$NGTag($i,tagid) == $type} {
	    return $i
	}
    }
    return $default;
}