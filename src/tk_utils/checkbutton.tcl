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
# Configures a checklist dialogue
#
proc checklist_configure {path args} {
    set in_arg 0
    set buttons ""
    set arglist ""
    set title ""
    set orient ""

    # Process command line args
    foreach i $args {
	if {$in_arg} {
	    if {$option == "-buttons"} {
		set buttons $i
	    } elseif {$option == "-default"} {
		set count 1
		foreach c $i {
		    global $path.Check_$count
		    set $path.Check_$count $c
		    incr count
		}
	    } elseif {$option == "-title"} {
		set title "-text {$i}"
	    } elseif {$option == "-orient"} {
		set orient "$i"
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

    if {"$title" != ""} {
        eval $path.label configure $title
    }

    # Add the items
    set count 1
    set comm ""
    foreach i $buttons {
	set args [lrange $i 1 end]
	eval checkbutton $path.r$count -text {[lindex $i 0]} \
	    -variable $path.Check_$count $args

	if {"$orient" == "horizontal"} {
	    pack $path.r$count -side right
	    set comm "raise [list $path.r$count]; $comm"
	} else {
	    pack $path.r$count -side top -anchor w
	}
        #pack $path.r$count -side top -anchor w
	incr count
    }
    if {$comm != ""} {eval $comm}

    return $orient
}

#
# checklist:
#
# Creates a checklist frame containing a title and a collection of buttons.
# Command line arguments are as per the frame widget, with the addition of:
#	-buttons {{button ?args?} {button ?args?} ...}
#	-default default_value
#	-title   title
#
proc checklist {path args} {
    # Create the frame and label
    frame $path -class CheckList
    label $path.label
    pack $path.label

    # Configure, which also adds the buttons
    set orient [eval checklist_configure $path $args]

    if {"$orient" == "horizontal"} {
	pack $path.label -side left
    } else {
	pack $path.label -side top -before $path.r1
    }
    return $path
}


#
# checklist_get:
#
# Given a checklist path we return the currently selected item
#
proc checklist_get {path inum} {
    global $path.Check_$inum
    eval return \$\{$path.Check_$inum\}
}


#
# checklist_destroy:
#
# Destroys a checklist path
#
proc checklist_destroy {path} {
    destroy $path
}
