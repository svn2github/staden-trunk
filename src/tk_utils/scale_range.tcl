#
# Copyrightquer (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

##############################################################################
#ensure start scale bar does not become larger than the end scale bar
proc CheckStartLimits { start end min_value value} {
    if {$value > [expr {[$end.scale get] - $min_value}]} {
	$end.scale set [expr {[$start.scale get] + $min_value}]
	if {[$end.scale get] < [expr {[$start.scale get] + $min_value}]} {
	    $start.scale set [expr {[$end.scale get] - $min_value}]
	}
    }
}

##############################################################################
#ensure end scale bar does not become smaller than the start scale bar
proc CheckEndLimits { start end min_value value} {
    if {$value < [expr {[$start.scale get] + $min_value}]} {
	$start.scale set [expr {[$end.scale get] - $min_value}]
	if {[$start.scale get] == 1} {
	    $end.scale set [expr {$min_value+1}]
	}
    }
}

#
#Configure a scale_range widget
#
proc scale_range_configureOLD {path args} {

    if {[winfo exists $path.from] } {
	eval scalebox_configure $path.from $args
	eval scalebox_configure $path.to $args
    }
}

proc scale_range_configure {path args} {

    set title ""
    set in_arg 0
    set arglist ""
    set command ""
    set checked_ok 0
    set end_value ""
    set start_value ""
    set start_name ""
    set end_name   ""
    set min_value ""

    foreach i $args {
	if {$in_arg} {
	    if {$option == "-title"} {
		set title $i
	    } elseif {$option == "-end_value"} {
		set end_value $i
	    } elseif {$option == "-start_value"} {
		set start_value $i
	    } elseif {$option == "-end_name"} {
		set end_name $i
	    } elseif {$option == "-start_name"} {
		set start_name $i
	    } elseif {$option == "-min_value"} {
		set min_value $i
	    } else {
		lappend arglist $option $i
		set in_arg 0
	    }
	     set in_arg 0
	} else {
	    set option $i
	    set in_arg 1
	}
    }

    if {[winfo exists $path.from] } {
	eval scalebox_configure $path.from $arglist
	eval scalebox_configure $path.to $arglist
    }

    if {$title != ""} {
	$path.l.label configure -text $title
    }
    if {$start_value != "" && [winfo exists $path.from]} {
	scalebox_configure $path.from -default $start_value
    }
    if {$end_value != "" && [winfo exists $path.to]} {
	scalebox_configure $path.to -default $end_value
    }
    if {$start_name != "" && [winfo exists $path.from]} {
	scalebox_configure $path.from -title $start_name
    }
    if {$end_name != "" && [winfo exists $path.to]} {
	scalebox_configure $path.to -title $end_name
    }

    if {$min_value != "" && [winfo exists $path.from]} {
	scalebox_configure $path.from -command "CheckStartLimits $path.from $path.to $min_value"
	scalebox_configure $path.from -command "CheckEndLimits $path.from $path.to $min_value"
    }
}

#
#Creates a scale_range widget
#from is the minimum value
#to is the maximum value
#start_value is the starting default for the from scale
#end_value is the starting default for the to scale
#min_value is the minimum difference allowed between the from & to scales
#
proc scale_range {path args} {

    set title ""
    set in_arg 0
    set arglist ""
    set command ""
    set checked_ok 0
    set from 1
    set to 1
    set end_value 0 
    set start_value 0
    set start_name "Start position"
    set end_name   "End position"
    set min_value 0

    foreach i $args {
	if {$in_arg} {
	    if {$option == "-title"} {
		set title $i
	    } elseif {$option == "-from"} {
		set from $i
	    } elseif {$option == "-to"} {
		set to $i
	    } elseif {$option == "-end_value"} {
		set end_value $i
	    } elseif {$option == "-start_value"} {
		set start_value $i
	    } elseif {$option == "-end_name"} {
		set end_name $i
	    } elseif {$option == "-start_name"} {
		set start_name $i
	    } elseif {$option == "-min_value"} {
		set min_value $i
	    } elseif {$option == "-command"} {
		set command "$i"
	    } else {
		lappend arglist $option $i
		set in_arg 0
	    }
	     set in_arg 0
	} else {
	    set option $i
	    set in_arg 1
	}
    }
    
    eval frame $path $arglist
    
    if {$title != ""} {
	frame $path.l
	label $path.l.label -text $title
	pack $path.l.label -side left
	pack $path.l -fill x
    }

    scalebox $path.from \
	-title $start_name \
	-orient horizontal \
	-from $from \
	-to $to \
	-default $start_value\
	-width 7 \
	-type CheckInt \
	-command "CheckStartLimits $path.from $path.to $min_value"
	
    scalebox $path.to \
	-title $end_name \
	-orient horizontal \
	-from $from \
	-to $to \
	-default $end_value \
	-width 7 \
	-type CheckInt \
	-command "CheckEndLimits $path.from $path.to $min_value"

    eval scale_range_configure $path $arglist

    pack $path.from -side top -fill x
    pack $path.to -side top -fill x

}

#
#Return the from value
#
proc scale_range_from {path} {

    return [scalebox_get $path.from]

}

#
#Return the to value
#
proc scale_range_to {path} {

    return [scalebox_get $path.to]

}

