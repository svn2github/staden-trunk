#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1998. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

proc get_def {var} {
    global $var
    if {[info commands ${var}_com] != ""}  {
	return ""
    } elseif {[info exists $var]} {
	return [set $var]
    } else {
	return ""
    }
}

proc is_command {var} {
    if {[info commands ${var}_com] != ""} {
	return 1
    } else {
	return 0
    }
}

#proc error {args} {
#    puts "ERROR: $args"
#}

# Checks that all paramaters listed in 'args' are global variables that both
# exists and have non blank contents. Alternatively if a ${param}_com funtion
# exists this will also be allowed.
proc global_param_exists {args} {
    foreach param $args {
	global $param
	if {!(([info exists $param] && [set $param] != "") ||
	      [info commands ${param}_com] != "")} { return $param }
    }

    return ""
}

# Checks that all the parameters listed in 'args' are local variables of the
# namespace 'ns' and that they have non blank contents.
proc local_param_exists {ns args} {
    foreach param $args {
	if {!([info exists ${ns}::$param] && [set ${ns}::$param] != "")} {
	    return $param
	}
    }

    return ""
}

#-----------------------------------------------------------------------------
# A simple drag and drop interface for moving objects between a pair of
# connected listboxes. The initial object is in the listbox with path $from,
# and the 'other' box is in $to.
# This also allows movement of objects within a single listbox
#
# Handle the initial grabbing of an object
proc dnd_grab {to from x y {delete_to 1} {delete_from 1} {drag_to 1} {drag_from 1}} {
    global grabbed_data
    global grabbed_from
    global grabbed_ind
    global grabbed_to
    global grabbed_delete
    global grabbed_drag
    
    set grabbed_delete($to) $delete_to
    set grabbed_delete($from) $delete_from
    set grabbed_drag($to) $drag_to
    set grabbed_drag($from) $drag_from

    set grabbed_ind [$from nearest $y]

    if {[$from get $grabbed_ind] != ""} {
	set grabbed_data [$from get $grabbed_ind]
	set grabbed_from $from
	set grabbed_to $to
    } else {
	set grabbed_data ""
    }
}

# Handle movement of grabbed_data around the windows
proc dnd_motion {x y} {
    global grabbed_data
    global grabbed_from
    global grabbed_to
    global grabbed_ind
    global grabbed_delete
    global grabbed_drag

    if {$grabbed_data == ""} {
	return
    }

    if {[set new_win [winfo containing $x $y]] == "" ||
	($new_win != $grabbed_to && $new_win != $grabbed_from)} {
	return
    }

    set new_x [expr $x-[winfo rootx $new_win]]
    set new_y [expr $y-[winfo rooty $new_win]]
    set nearest [$new_win nearest $new_y]
    if {[set bbox [$new_win bbox $nearest]] != ""} {
        set top_y [lindex [$new_win bbox $nearest] 1]
        set bot_y [expr $top_y+[lindex [$new_win bbox $nearest] 3]]
        if {[expr ($new_y-$top_y) > ($bot_y-$new_y)]} {
	    incr nearest
        }
    }

    if {$grabbed_from == $new_win} {
	if {$grabbed_drag($grabbed_from)} {
	    # Moving within a window
	    if {$grabbed_ind >= $nearest} {
	        $grabbed_from delete $grabbed_ind
	        $grabbed_from insert $nearest $grabbed_data
	        set grabbed_ind $nearest
	    } else {
	        $grabbed_from insert $nearest $grabbed_data
	        $grabbed_from delete $grabbed_ind
	        set grabbed_ind [expr $nearest-1]
	    }
	    $grabbed_from select clear 0 end
	    $grabbed_from select set $grabbed_ind
	}
    } else {
	# Swapping windows - so swap grabbed_from and grabbed_to
	if {$grabbed_delete($grabbed_from)} {
	    $grabbed_from delete $grabbed_ind
	}
	if {$grabbed_delete($grabbed_to)} {
	    $grabbed_to insert $nearest $grabbed_data
	} else {
	    set ind [lsearch [$grabbed_to get 0 end] $grabbed_data]
	    if {$ind != -1} {
	        set nearest $ind
	    }
	}
	$grabbed_from select clear 0 end
	$grabbed_to select clear 0 end
	set f $grabbed_to
	set grabbed_to $grabbed_from
	set grabbed_from $f
	set grabbed_ind $nearest
	$grabbed_from select clear 0 end
	$grabbed_from select set $grabbed_ind
	$grabbed_from see $grabbed_ind
    }
}

# tidy up
proc dnd_release {} {
    global grabbed_data
    set grabbed_data ""
}

# Replace newlines with '\n'
proc strip_nl {str} {
   regsub -all "\n" $str {\n} str
   return $str
}

# Force a filename to be a full pathname
proc fullpath {file {dir {}} {only_if_not_present 0}} {
    if {$dir == ""} {
	set dir [pwd]
    }
    if {[file pathtype $file] != "absolute"} {
	if {$only_if_not_present} {
	    if {[file exists $file]} {
		return $file
	    }
	}
	return [file join $dir $file]
    } else {
	return $file
    }
}
