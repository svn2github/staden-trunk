#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#
#Configure a contig_id widget
#
proc contig_id_configure {path args} {

    eval entrybox_configure $path.ent $args

    if {[winfo exists $path.lreg] } {
	eval scalebox_configure $path.lreg $args
	eval scalebox_configure $path.rreg $args
    }
}

#
#Creates a contig_id widget
#
#-io io handle - must be supplied
#-range boolean - whether to draw the start and end scalebars (default 1)
#-trace boolean - whether to update with the contig selector (defalult 1)
#-default value
proc contig_id {path args} {
    global LREG
    global RREG
    global CurContig
    global $path.Io
    global db_namelen

    set in_arg 0
    set arglist ""
    set io ""
    set command ""
    set checked_ok 0
    set range 1
    set end_value $RREG
    set start_value $LREG
    set trace 1
    set default "$CurContig"

    frame $path -relief groove -bd 2

    #extract io info from argument list
    foreach i $args {
	if {$in_arg} {
	    if {$option == "-io"} {
		set $path.Io $i
		set io $i
		set in_arg 0
	    } elseif {$option == "-range"} {
		set range $i
	    } elseif {$option == "-end_value"} {
		set end_value $i
	    } elseif {$option == "-start_value"} {
		set start_value $i
	    } elseif {$option == "-command"} {
		set command "$i"
      	    } elseif {$option == "-trace"} {
		set trace "$i"
      	    } elseif {$option == "-default"} {
		set default "$i"
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

    if {$range} {
        set c_num [db_info get_contig_num $io $CurContig]
       	if {$c_num < 0} {
	    InitContigGlobals $io
            set c_num [db_info get_contig_num $io $CurContig]
        }
        set length [c_length $io $c_num]
	entrybox $path.ent \
		-title "Contig identifier" \
		-type "CheckContig $io"\
		-default $default \
	        -command "update_limits $io $path" \
		-width $db_namelen

	scalebox $path.lreg \
		-title "Start position" \
		-orient horizontal \
		-from 1 \
		-to $length \
		-default $start_value\
		-width 7 \
		-type CheckInt \
		-command "CheckStartLimits $path.lreg $path.rreg 0"
	
	scalebox $path.rreg \
		-title "End position" \
		-orient horizontal \
		-from 1 \
		-to $length \
		-default $end_value \
		-width 7 \
		-type CheckInt \
		-command "CheckEndLimits $path.lreg $path.rreg 0"
    } else {
	entrybox $path.ent \
		-title "Contig identifier" \
		-type "CheckContigName $io"\
		-default $default \
		-width $db_namelen \
		-command "contig_id_callback2 [list $command]"
    }
    eval contig_id_configure $path $arglist
    entrybox_configure $path.ent -exportselection 0
    
    [entrybox_path $path.ent] selection range 0 end
    focus [entrybox_path $path.ent]

    if {$range} {
	#bindings
	bind $path.ent.entry <Return> "contig_id_callback $io %W $path; $command"
	bind $path.ent.entry <Any-FocusOut> "contig_id_callback $io %W $path"
	bind $path.ent.entry <Any-Leave> "contig_id_callback $io %W $path"
    }

    pack $path.ent -side top -fill x 
    if {$range} {
	pack $path.lreg -side top -fill x
	pack $path.rreg -side top -fill x
    }

    if {"$trace" == 1} {
        # Monitor the contig_selector for changes in the selected contig.
        global c_id_contig
        trace variable c_id_contig w "contig_id_trace $path $io"
	
        # Make sure we remove the trace when this widget is destroyed.
        bind $path <Any-Destroy> "contig_id_destroy $path $io"
    }
    if {"$trace" == 2} {
        # Monitor the contig_selector for changes in the selected contig.
        global c_id_contig

        trace variable c_id_contig w "contig_id_trace \[GetCurFrame \[winfo parent $path\]\] $io"
	
        # Make sure we remove the trace when this widget is destroyed.
        bind $path <Any-Destroy> "contig_id_destroy $path $io"
    }
}

proc update_limits {io path c_name } {

if {[winfo exists $path.lreg]} {
        if {[set c_num [db_info get_contig_num $io $c_name]] != -1} {
            set len [c_length $io $c_num]
            scalebox_configure $path.lreg -to $len -default 1
            scalebox_configure $path.rreg -to $len -default $len
        }
}
}

#
# The trace callback
#
proc contig_id_trace {path io name element op} {
    global $name

    entrybox_delete $path.ent 0 end
    entrybox_insert $path.ent 0 [set $name]
    if {[winfo exists $path.lreg]} {
	if {[string compare normal \
		[lindex [$path.lreg.scale configure -state] 4]] == 0} {

	    if {[set c_num [db_info get_contig_num $io [set $name]]] != -1} {
		set len [c_length $io $c_num]
		scalebox_configure $path.lreg -to $len -default 1
		scalebox_configure $path.rreg -to $len -default $len
	    }
	}
    }
}


#
# Shutdown the contig_id to remove the trace
#
proc contig_id_destroy {path io} {
    global c_id_contig
    global $path.Io

    trace vdelete c_id_contig w "contig_id_trace $path $io"
    trace vdelete c_id_contig w "contig_id_trace \[GetCurFrame \[winfo parent $path\]\] $io"
    unset $path.Io
}


#
#Return the gel name in the entry box
#Updates the global CurContig, LREG and RREG too.
#
proc contig_id_gel {path} {
    global $path.Io

    set r [entrybox_get $path.ent]
    if {"$r" != ""} {
	if [winfo exists $path.lreg] {
	    SetContigGlobals [set $path.Io] $r \
		[contig_id_lreg $path] \
		[contig_id_rreg $path]
	} else {
	    SetContigGlobals [set $path.Io] $r
	}
    }

    return $r
}

#
#Return the lreg value
#
proc contig_id_lreg {path} {

    return [scalebox_get $path.lreg]

}
#
#Return the rreg value
#
proc contig_id_rreg {path} {

    return [scalebox_get $path.rreg]

}

#
# Focus on the name region
#
proc contig_id_focus {path} {
    entrybox_focus $path.ent
}

#
# Only called when "-range 1" is used.
#
# This keeps check of the last name entered and only resets the values if the
# reading name has changed since then. Note that without this it is possible
# to get event loops due to this callback being called from the
# Any-Leave binding. If the "last" bit is removed, make sure that the
# bindings are also adjusted.
#
proc contig_id_callback {io entry path} {
    global $entry.Last

    if {![info exists $entry.Last]} {
	set $entry.Last ""
    }
    if {[string compare [set $entry.Last] [$entry get]] != 0} {
	set $entry.Last [$entry get]
	UpdateContigLimits $io $path.lreg.scale $path.rreg.scale $entry
    }
}

# Only called when "-range 0" is used.
proc contig_id_callback2 {command args} {
    uplevel #0 $command
}
