#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#
#Configure a seq_id widget
#
proc seq_id_configure {path args} {

    set in_arg 0
    set arglist ""
    set command ""
    set checked_ok 0
    set range 1
    set direction 0
    set browse 1
    set default ""
    set end ""
    set start ""
    set from ""
    set to ""
    set title ""
    set min_value 1
    set update_cmd {""}
    set browse_cmd {""}
    set feature ""
    set single ""

    #extract io info from argument list
    foreach i $args {
	if {$in_arg} {
	    if {$option == "-range"} {
		set range $i
	    } elseif {$option == "-title"} {
		set title $i
	    } elseif {$option == "-from"} {
		set from $i
	    } elseif {$option == "-to"} {
		set to $i
	    } elseif {$option == "-end_value"} {
		set end $i
	    } elseif {$option == "-start_value"} {
		set start $i
	    } elseif {$option == "-min_value"} {
		set min_value $i
	    } elseif {$option == "-command"} {
		set command "$i"
	    } elseif {$option == "-update_cmd"} {
		set update_cmd "$i"
      	    } elseif {$option == "-browse"} {
		set browse "$i"
      	    } elseif {$option == "-browse_cmd"} {
		set browse_cmd "$i"
	    } elseif {$option == "-default"} {
		set default "$i"	
            } elseif {$option == "-feature"} {
		set feature "$i"
	    } elseif {$option == "-single"} {
		set single "$i"
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

    eval $path.ent configure $arglist
}

proc CheckSequenceName {name} {

    set seq_id [name_to_seq_id $name]
    if {$seq_id == -1} {
 	tk_messageBox -icon error -type ok -title "Sequence identifer" \
		-message "Sequence name does not exist"
	return 1
    }
    return 0
}

proc CheckSequence {name end} {

    if {[CheckSequenceName $name] == 0} {
	set seq_id [name_to_seq_id $name]
	#now check if ranges are correct
	set length [seq_info $seq_id length]
	if {$end > $length} {
	    tk_messageBox -icon error -type ok -title "Sequence identifer"\
		-message "You have chosen an end \
		    point greater than the length of the sequence. Try \
		    pressing return after entering the seq identifier \
		    to update the sequence limits"
	    return 1
	} 
	return 0
    }
    return 1
}

proc CheckSeq {path } {

#    set name [$path.ent.entry get]
#    set end [scale_range_to $path.range]

    set name [$path.ent get_seqname]
    set end [$path.ent get_e]
    eval set response \{[CheckSequence $name $end]\}
    if { $response == 1} {
	return 0
    }
    return 1
}

proc CheckSeqName {path } {

    set name [$path.entry get]

    [entrybox_path $path] xview [string last / $name]
    eval set response \{[CheckSequenceName $name]\}
    if { $response == 1} {
	return 0
    }
    return 1
}

#
#Creates a seq_id widget
#
#-range boolean - whether to draw the start and end scalebars (default 1)
#-direction boolean - whether to allow selection of seq for horizontal or 
#                     vertical
#-browse boolean - whether to create a sequence browse button
#-default value


proc seq_id {path args} {
    global HORIZONTAL $path.update_cmd SeqName CurFrame spin_defs

    set in_arg 0
    set arglist ""
    set command ""
    set checked_ok 0
    set range 1
    set direction 0
    set browse 1
    set default ""
    set end 100
    set start 1
    set from 1
    set to 100
    set title ""
    set min_value 1
    set update_cmd {""}
    set browse_cmd {""}
    set feature ""
    set single ""

    frame $path
    bind $path <Destroy> "seq_id_destroy $path"

    #extract io info from argument list
    foreach i $args {
	if {$in_arg} {
	    if {$option == "-range"} {
		set range $i
	    } elseif {$option == "-title"} {
		set title $i
	    } elseif {$option == "-default"} {
		set default $i
	    } elseif {$option == "-from"} {
		set from $i
	    } elseif {$option == "-to"} {
		set to $i
	    } elseif {$option == "-end_value"} {
		set end $i
	    } elseif {$option == "-start_value"} {
		set start $i
	    } elseif {$option == "-min_value"} {
		set min_value $i
	    } elseif {$option == "-command"} {
		set command "$i"
	    } elseif {$option == "-update_cmd"} {
		set update_cmd "$i"
      	    } elseif {$option == "-browse"} {
		set browse "$i"
      	    } elseif {$option == "-browse_cmd"} {
		set browse_cmd "$i"
            } elseif {$option == "-feature"} {
		set feature "$i"
	    } elseif {$option == "-single"} {
		set single "$i"
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
    set $path.update_cmd $update_cmd

    frame $path.id
#    if {$title != ""} {
#	label $path.title -text $title
#	pack $path.title -anchor w -side top
#    }
###################################################################
    if {$range} {
	if {$feature == "yes" && $single == "yes"} {	     
	    xcomborange $path.ent -labeltext "$title" \
		-feature "$feature" \
		-single "$single" \
		-default "$default" \
		-start_value $start \
		-end_value $end
	} elseif {$feature == "yes"} {
	    xcomborange $path.ent \
		-labeltext "$title" \
		-feature "$feature" \
		-default "$default" \
		-start_value $start \
		-end_value $end
	} elseif {$single == "yes"} {
	    xcomborange $path.ent \
		-labeltext "$title" \
		-single "$single" \
		-default "$default" \
		-start_value $start \
		-end_value $end 
	} else {
	    xcomborange $path.ent \
		-labeltext "$title" \
		-default "$default" \
		-start_value $start \
		-end_value $end
	}
    } else {
	xcomborange $path.ent \
	    -trange no \
	    -default "$default"\
	    -start_value $start \
	    -end_value $end
    }

    eval seq_id_configure $path $arglist
  
    global seq_id_loop
    set seq_id_loop 0

    if {$range} {
	#bindings
	bind $path.ent <Return> "set seq_id_loop 0; seq_id_callback $path $update_cmd; $command"
	bind $path.ent <Any-FocusOut> "seq_id_callback $path $update_cmd"
	bind $path.ent <Any-Leave> "seq_id_callback $path $update_cmd"
    } else {
	bind $path.ent <Return> "$command"
	bind $path.ent <Return> "set seq_id_loop 0;seq_id_callback $path $update_cmd; $command"
	bind $path.ent <Any-FocusOut> "seq_id_callback $path $update_cmd"
	bind $path.ent <Any-Leave> "seq_id_callback $path $update_cmd"    }
    
	pack $path.id -side top -fill x
	pack $path.ent -in $path.id -expand 1 -fill x 

	set CurFrame $path
	trace variable SeqName w "seq_id_trace $path $update_cmd"
}

proc seq_id_destroy {path} {
    global $path.update_cmd SeqName CurFrame

    if {[info exists CurFrame]} {
	unset CurFrame
    }
    trace vdelete SeqName w "seq_id_trace $path [set $path.update_cmd]"
}

#
#Return the sequence name in the entry box
#
proc seq_id_name {path {avoid {}}} {

    set re_enter 0
    set name [$path.ent get_seqname]
    
    if {$name == ""} {
	bell

	#awful hack to avoid tk_messageBox invoking the leave binding 
	#specifically for sip_similar_spans.tcl win len updates
	if {$avoid != ""} {
	    set b [bind $avoid <Any-Leave>]
	    bind $avoid <Any-Leave> {}
	}
	tk_messageBox -icon error -type ok -title "Sequence identifer" \
		-message "Sequence name does not exist" -parent $path
	if {$avoid != ""} {
	    bind $avoid <Any-Leave> $b
	}
	raise [winfo toplevel $path]
	    
	tkwait variable re_enter
	return $name
    }
    return $name
}

#
#Return the lreg value
#
proc seq_id_from {path} {
    set re_enter 0
    set start [$path.ent get_s]
    set end [$path.ent get_e]

    if { $start > $end  } {
	bell
	catch ClearBusy
	tk_messageBox -icon error -type ok -title "Sequence identifer"\
		-message "You have chosen a start position \
		greater than the end position "

	raise [winfo toplevel $path]

	tkwait variable re_enter
    }

    return [$path.ent get_s]
}

#
#Return the rreg value
#
proc seq_id_to {path} {
    
    set re_enter 0
    set start [$path.ent get_s]
    set end [$path.ent get_e]
    if { $start > $end } {
	bell
	catch ClearBusy
	tk_messageBox -icon error -type ok -title "Sequence identifer"\
		-message "You have chosen a start position\
	        greater than the end position "

	raise [winfo toplevel $path]

	tkwait variable re_enter 
    }
    return [$path.ent get_e]  
}

#
#Return the single entry value
#

proc seq_id_single {path} {

    return [$path.ent get_single]  
}

#
#Return the tag of the currently selected radiobutton. 
#
proc seq_id_feat {path} {

    return [$path.ent get_selected]
}

#
#Returns the contents of the listbox element indicated by the current 
#selection indexes.
#
proc seq_id_cds {path} {
    return [$path.ent get_curselection]
}

#
# Focus on the name region
#
proc seq_id_focus {path} {
    entrybox_focus $path.ent
}

proc seq_id_trace {path update_cmd name element op} { 
    global CurFrame

    if {[info exists CurFrame]} {
	set path $CurFrame
    }

    upvar $name x
    seq_id_configure $path -default $x
    seq_id_callback $path $update_cmd
}

# Callback from focus-out/leave binding on the sequence identifier.
# This function has been crafted more out of trial and error than
# understanding! Basically, due to interactions between tk_dialog and the
# contig_id leave/focus-out bindings, it's possible to get continuous error
# loops. I suspect that this is because tk_dialog does some very odd things
# (ie create, withdraw, update, redisplay) that forces multiple binding
# callbacks for a single error.
# For this reason we need to check whether we've already entered this routine
# once to block these requests. It's ugly and is likely to break at the
# slightest change. You have been warned...
proc seq_id_callback {path update_cmd} {
    global seq_id_loop

    if {$seq_id_loop == 0} {
	set seq_id_loop 1
#    	focus $path
    	update
    	UpdateSeqLimits $path $update_cmd
       	set seq_id_loop 0
    }
}

########### FIXME #################

proc UpdateSeqLimits {path update_cmd} {
    global $path.old_id

    #for case of sequence name only
    if {![winfo exists $path.range.from]} {
	if {$update_cmd != ""} {
	    eval $update_cmd
	}
	return
    }
    
    set name [$path.ent get_seqname]
#    set name [entrybox_get $path.ent]
    set seq_id [name_to_seq_id $name]
    set start [seq_info $seq_id start]
    set end [seq_info $seq_id end]
    set length [seq_info $seq_id length]

    if {![info exists $path.old_id]} {
	set $path.old_id $seq_id
    }

    if {$update_cmd != ""} {
	eval $update_cmd
    }

    set $path.old_id $seq_id
}

