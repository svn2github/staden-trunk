#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

# jkb: 9/2/00
# Why does this widget pack items with the last button left-most for
# horizontal packing and first button top-most for vertical packing!?
# Too late to change...

#
# Configures a radiolist dialogue
#
proc radiolist_configure {path args} {
    global $path.Radio
    global $path.num_buts

    set in_arg 0
    set buttons ""
    set arglist ""
    set title ""
    set orient ""
    set font ""
    set state ""
    set foreground ""

    # Process command line args
    foreach i $args {
	if {$in_arg} {
	    if {$option == "-buttons"} {
		set buttons $i
	    } elseif {$option == "-default"} {
		set $path.Radio $i
	    } elseif {$option == "-title"} {
		set title "-text {$i}"
	    } elseif {$option == "-orient"} {
		set orient "$i"
	     } elseif {$option == "-font"} {
		set font "-font {$i}"
	     } elseif {$option == "-state"} {
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

    eval $path configure $arglist

    if {"$title" != "" || "$foreground" != "" || "$font" != "" || $state != ""} {
        eval $path.label configure $title $foreground $font $state
    }

    if {"$state" != ""} {
	for {set i 1} {$i < [set $path.num_buts]} {incr i} {
	    eval $path.r$i configure $state
	}
    }

    # Add the items
    set count 1
    set comm ""
    foreach i $buttons {
	set args [lrange $i 1 end]
	eval radiobutton $path.r$count -text {[lindex $i 0]} -value $count \
	    -variable $path.Radio $args 

	if {"$orient" == "horizontal"} {
	    pack $path.r$count -side right
	    set comm "raise [list $path.r$count]; $comm"
	} else {
	    pack $path.r$count -side top -anchor w
	}
	incr count
	set $path.num_buts $count
    }
    if {$comm != ""} {eval $comm}

    return $orient
}

#
# radiolist:
#
# Creates a radiolist frame containing a title and a collection of buttons.
# Command line arguments are as per the frame widget, with the addition of:
#	-buttons {{button ?args?} {button ?args?} ...}
#	-default default_value
#	-title   title
#
proc radiolist {path args} {
    global $path.Radio

    # Create the frame and label
    frame $path -class RadioList
    xlabel $path.label

    # Configure, which also adds the buttons
    set orient [eval radiolist_configure $path $args]

    if {"$orient" == "horizontal"} {
	pack $path.label -side left
    } else {
	pack $path.label -side top -before $path.r1
    }

    # Invoke the default operations
    eval set r \$\{$path.Radio\}
    eval $path.r$r invoke

    return $path
}


#
# radiolist_get:
#
# Given a radiolist path we return the currently selected item
#
proc radiolist_get {path} {
    global $path.Radio
    eval return \$\{$path.Radio\}
}


#
# radiolist_destroy:
#
# Destroys a radiolist path
#
proc radiolist_destroy {path} {
    global $path.Radio
    unset $path.Radio

    destroy $path
}

