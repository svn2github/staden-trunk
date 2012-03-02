#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
##############################################################################
#convert entry values to scale values
proc scalevalue { path } {
    global entry_

    set entryvalue [$path.entry get]
    if {[info exists entry_($path,type)] } {

	#check entry is of correct type
	if { ![eval $entry_($path,type) $path] } {
	    bell
	} else {
	    set v [$path.scale get]
	    #HACK - to deal simply with the problem that once the scale is
	    #at it's limits, it doesn't update and therefore the 
	    #scalebox_update_entry function is never called and the entrybox
	    #isn't set to the scale limit
	    $path.scale set [expr $v-1]
	    $path.scale set [expr $v+1]
	    $path.scale set $entryvalue
	}
	#$path.entry delete 0 end
    }
}

#
# Configures a scalebox
#
proc scalebox_configure {path args} {
    global entry_

    set in_arg 0
    set arglist ""
    set title ""
    set command ""
    set type CheckInt
    set orient ""
    set from ""
    set to ""
    set state ""
    set sarg 0
    set earg 0
    set default ""
    set resolution ""
    set digits ""
    set varg ""

    #set default type to be CheckInt
    set entry_($path,type) $type

    # Process command line args
    foreach i $args {
	if {$in_arg} {
	    if {$option == "-title"} {
		set title "-text {$i}"
	     } elseif {$option == "-type"} {
		 set entry_($path,type) "$i"
	     } elseif {$option == "-orient" || $option == "orientation"} {
		set orient "-orient $i"
		set sarg 1
	     } elseif {$option == "-resolution"} {
		set from "-resolution $i"
		set sarg 1
	     } elseif {$option == "-digits"} {
		set from "-digits $i"
		set sarg 1
	     } elseif {$option == "-from"} {
		set from "-from $i"
		set sarg 1
	     } elseif {$option == "-to"} {
		set to "-to $i"
		set sarg 1
	     } elseif {$option == "-default"} {
		 set default $i
	     } elseif {$option == "-state"} {
		set state "$i"
	     } elseif {$option == "-width"} {
		set width "-width $i"
		set earg 1
	     } elseif {$option == "-variable"} {
		 set varg $i
	     } elseif {$option == "-command"} {
		 #add new command before original command
		 set comm [lindex [$path.scale configure -command] 4]
		 set command [format \
			 {-command {%s [%s get]; %s} } \
			 $i [list $path.scale] $comm]
		 #set command "-command {$i}"
		 set sarg 1
	     } else {
		lappend arglist $option $i
	     }

	     set in_arg 0
	} else {
	     set option $i
	     set in_arg 1
	}
    }

    eval $path.scale configure $arglist

    if {"$title" != ""} {
	eval $path.label configure $title
    }
    if {$sarg} {
	eval $path.scale configure $from $to $orient $command $digits $resolution
    }
    if {$earg} {
	eval $path.entry configure $width
    }
    if {$varg != ""} {
	eval $path.entry configure -textvariable $varg
    }
    if {$default != ""} {
	eval $path.scale set $default
	scalebox_update_entry $path.scale $path.entry $default
    }

    if {"$state" == "normal"} {
	$path.label configure -state normal
	$path.scale configure -state normal
	$path.entry configure -state normal
    } elseif {"$state" == "disabled"} {
	$path.label configure -state disabled
	$path.scale configure -state disabled
	$path.entry configure -state disabled
    }
}

proc scalebox_update_entry {path entry value } {
    
    $entry delete 0 end
    #must use [$path get] here and NOT value
    $entry insert 0 [$path get]
}

#
# Creates a scalebox
# Command line arguments are as per the frame widget, with the addition of:
#       -title     title
#	-command   script
#	-default   value
#       -resolution value
#       -digits    value
#	-from	   value
#	-to	   value
#	-orient	   vert/horiz
#	-state	   normal/disabled
#	-type	   CheckInt/CheckFloat
#
proc scalebox { path args } {
    global entry_

    # Create the frame, label, entry and scale
    frame $path -class ScaleBox

    #create label and entry frame
    frame $path.le

    xlabel $path.label -anchor w
    entry $path.entry
    scale $path.scale -showvalue 0 -command "scalebox_update_entry $path.scale $path.entry"

    # Configure
    eval scalebox_configure $path $args

    pack $path.label -side left -in $path.le
    pack $path.entry -side right -in $path.le 
    pack $path.scale -side right

    pack $path.le -side bottom -pady 2 -fill x

    #update slider value from entrybox value
    bind $path.entry <Return> "scalevalue $path"
    bind $path.entry <Any-Leave> "CheckLeave $path"
    bind $path.entry <Any-FocusOut> "CheckLeave $path"

    return $path
}

proc CheckLeave { path } {

    if {[$path.entry get] != ""} {
	scalevalue $path
    }

}

#
# scalebox_get:
#
# Given a scalebox path we return the current value
#
proc scalebox_get {path} {
    CheckLeave $path
    return [$path.scale get]
}

proc scalebox_set {path val} {
    $path.scale set $val
    scalebox_update_entry $path.scale $path.entry $val
}

#
# scalebox_destroy:
#
# Destroys a scalebox path
#
proc scalebox_destroy {path} {
    global $path.Scale
    unset $path.Scale

    destroy $path
}

