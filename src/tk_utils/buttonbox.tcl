#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#
# Configures a okcancelhelp dialog
#
proc okcancelhelp_configure {path args } {
    set in_arg 0
    set arglist ""
    set ok_command ""
    set apply_command ""
    set cancel_command ""
    set help_command ""
    set perm_command ""
    set asdefault_command ""

    # Process command line args
    foreach i $args {
	if {$in_arg} {
	    if {$option == "-ok_command"} {
		set ok_command "-command {$i}"
	    } elseif {$option == "-apply_command"} {
		set apply_command "-command {$i}"
	    } elseif {$option == "-cancel_command"} {
		set cancel_command "-command {$i}"
	    } elseif {$option == "-help_command"} {
		set help_command "-command {$i}"
	    } elseif {$option == "-perm_command"} {
		set perm_command "-command {$i}"
	    } elseif {$option == "-asdefault_command"} {
		set asdefault_command "-command {$i}"
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
    if {"$ok_command" != ""} {
	eval $path.ok configure $ok_command
    }
    if {"apply_command" != ""} {
	eval $path.apply configure $apply_command
    }
    if {"$cancel_command" != ""} {
	eval $path.cancel configure $cancel_command
    }
    if {"$help_command" != ""} {
	eval $path.help configure $help_command
    }
    if {"perm_command" != ""} {
	eval $path.perm configure $perm_command
    }
    if {"asdefault_command" != ""} {
	eval $path.asdefault configure $asdefault_command
    }
}

#
# Creates an OK Cancel Help buttons
# Command line arguments are as per the frame widget, with the addition of:
#	-ok_command  script
#	-asdefault_command script
#	-apply_command script
#	-perm_command script
#	-cancel_command  script
#	-help_command script
#
proc okcancelhelp {path args } {
    set in_arg 0
    set arglist ""

    # Create the frame with our arguments
    frame $path -class OKCancelHelp 

    button $path.ok -text OK
    button $path.apply -text Apply
    button $path.asdefault -text "As Default"
    button $path.cancel -text Cancel
    button $path.help -text Help
    button $path.perm -text "OK Permanent"
    # Configure
    eval okcancelhelp_configure $path $args

    pack $path.ok -side left -expand 1 -anchor c
    if {[$path.perm cget -command] != ""} {
	pack $path.perm -side left -expand 1 -anchor c
    }

    if {[$path.apply cget -command] != ""} {
	pack $path.apply -side left -expand 1 -anchor c
    }

    if {[$path.asdefault cget -command] != ""} {
	pack $path.asdefault -side left -expand 1 -anchor c
    }

    pack $path.cancel $path.help -side left -expand 1 -anchor c

    # Add return binding to parent window.
    bind [winfo toplevel $path] <Return> "catch {$path.ok invoke}"
}

proc okapplycancelhelp {args} {
    return [eval okcancelhelp $args]
}

