#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#
# Configures a yesno dialog
#
proc yes_no_configure {path args } {
    global $path.YesNo

    set in_arg 0
    set arglist ""
    set title ""
    set orient ""
    set ycommand ""
    set ncommand ""
    set state ""

    # Process command line args
    foreach i $args {
	if {$in_arg} {
	    if {$option == "-title"} {
		set title "-text {$i}"
	    } elseif {$option == "-state"} {
		set state $i
	    } elseif {$option == "-default"} {
		set $path.YesNo $i
	    } elseif {$option == "-orient"} {
		set orient "$i"
	    } elseif {$option == "-ycommand"} {
		set ycommand "-command {$i}"
	    } elseif {$option == "-ncommand"} {
		set ncommand "-command {$i}"
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

    if {$state != ""} {
	eval $path.yes configure -state $state
	eval $path.no configure -state $state
	eval $path.label configure -state $state
    }
    if {"$ycommand" != ""} {
	eval $path.yes configure $ycommand
    }

    if {"$ncommand" != ""} {
	eval $path.no configure $ncommand
    }

    if {"$title" != ""} {
	eval $path.label configure $title
    }
    return $orient
}

#
# Creates a yes no dialogue.
# Command line arguments are as per the frame widget, with the addition of:
#       -title     title
#	-ycommand  script
#	-ncommand  script
#	-default   0/1
#       -state     normal/disabled
#
proc yes_no {path args} {
    global $path.YesNo

    # Create the frame, label and radiobuttons
    frame $path -class YesNo
    xlabel $path.label
    radiobutton $path.no -text No -variable $path.YesNo -value 0
    radiobutton $path.yes -text Yes -variable $path.YesNo -value 1

    # Configure
    eval yes_no_configure $path $args

    # Pack and return
    # Configure, which also adds the buttons
    set orient [eval yes_no_configure $path $args]

    if {"$orient" == "horizontal"} {
	pack $path.label -side left
	pack $path.no $path.yes -side right
    } else {
	pack $path.label -side top
	pack $path.yes $path.no -side top -anchor w
    }
    raise $path.no

    #pack $path.label -side top
    #pack $path.no $path.yes -side left

    # Invoke the default operations
    eval set y \$\{$path.YesNo\}
    if {"$y" == 1} {
        $path.yes invoke
    } elseif {"$y" == 0} {
        $path.no invoke
    }

    return $path
}

#
# yes_no_get:
#
# Given a yes no path we return the currently selected item
#
proc yes_no_get { path } {
    global $path.YesNo

    eval return \$\{$path.YesNo\}
}


#
# yes_no_destroy:
#
# Destroys a yes_no path
#
proc yes_no_destroy {path} {
    global $path.YesNo
    unset $path.YesNo

    destroy $path
}
