#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
global InCommunication
set InCommunication 0

proc SendTo {} {
    global HORIZONTAL sip_defs ssend

    set t .sendto
    if {[xtoplevel $t -resizable 0] == ""} return
    
    seq_id $t.seq_id -range 0 -default [get_active_seq_name $HORIZONTAL]\
	-browse_cmd sip_seq_browser

    frame $t.interps -bd 2 -relief raised
    label $t.interps.title -text "Program name"
    frame $t.interps.lb
    listbox $t.interps.lb.l -height 5 -yscrollcommand "$t.interps.lb.s set"
    scrollbar $t.interps.lb.s -orient vert -command "$t.interps.lb.l yview"

    set interps ""
    set interps_done ""
    set thisapp [tk appname]
    catch {unset ssend};set ssend() ""; unset ssend(); # create a blank array
    foreach i [keylget sip_defs SIP_SEND_TO] {
	set ssend([lindex $i 0]) [lindex $i 1]
    }
    
    foreach i [winfo interps] {
	if {$i == $thisapp} {
	    continue
	}
	if {[info exists gsend($i)]} {
	    lappend interps_done $i
	}

	#puts "checking $i"
	if {[catch "send $i {info commands EventHandler}" ret] == 1} {
	    #puts "failing $i"
	    continue
	}
	
	#puts "ret=$ret"
	if {$ret != ""} {
	    $t.interps.lb.l insert end $i
	    lappend interps $i
	}
    }
    foreach i [array names ssend] {
	if {[lsearch -exact $interps_done $i] == -1} {
	    $t.interps.lb.l insert end $i
	    lappend interps $i
        }
    }

    pack $t.seq_id
    pack $t.interps.lb.l -side left -fill both -expand 1
    pack $t.interps.lb.s -side left -fill both
    pack $t.interps.title $t.interps.lb -side top -fill both

    okcancelhelp $t.but \
	-ok_command "SendTo2 $t $t.seq_id $t.interps.lb.l ssend" \
	-cancel_command "destroy $t" \
	-help_command "show_help interface {UI-SendTo}"

    pack $t.interps -side top -fill both
    pack $t.but -side bottom -fill both
}

proc SendTo2 {t s_id lb ssend_p} {
    upvar #0 $ssend_p ssend
    global seq sip_defs cur_dir
    global HORIZONTAL VERTICAL

    if {[set name [$lb curselection]] == ""} {
	focus $lb
	bell
	return
    }
    set name [$lb get $name]

    # If necessary, start up child process and wait for it to be ready.
    if {[lsearch -exact [winfo interps] $name] == -1} {
	eval exec $ssend($name) &
    }
    set i 0
    while {[lsearch -exact [winfo interps] $name] == -1} {
	if {[incr i] == 10} {
	    verror ERR_WARN sendto "couldn't start $name"
	    bell
	    return
        }
	exec sleep 1
    }

    #get seq name from seq_id entrybox
    set seq_name [seq_id_name $s_id]
    set seq_id [name_to_seq_id $seq_name]
    set seq_num [get_seq_num $seq_id]

    destroy $t
    set l ""
    keylset l sequence [seq_info $seq_id sequence]

    keylset l command [list send [tk appname]]
    keylset l from [tk appname]
    keylset l seq_id $seq_name

    set rid [send $name EventHandler [list [tk appname]=$seq_name] \
	INIT_SEQUENCE $l]
    #puts ">>>>>>>> SENT [tk appname]=$seq_name    RID = $rid <<<<<<<<<<<<"

    global commn_${rid}
    upvar #0 commn_${rid} commn
    set commn(command) "send {$name}"
    set commn(seq_id) $seq_id

    #if loaded sequence via sip file menu
    if {![info exists cur_dir] || $cur_dir == $VERTICAL} {
	set cur_dir $HORIZONTAL
    } else {
	set cur_dir $VERTICAL
    }

    #create cursors only when necessary, not on sending
    #set commn(cursor_id) [create_cursor -seq_num $seq_num -line_width [keylget sip_defs SIP.CURSOR.LINE_WIDTH] -direction $cur_dir]

    seq_sender -rid $rid -seq_id $seq_id
}

proc LoadSequence {seq name} {
    global HORIZONTAL VERTICAL cur_dir
    set w .tmp
    global $w.success_h $w.success_v

    if {[info exists cur_dir]} {
	#alternate between horizontal and vertical sequence
	if {$cur_dir == $HORIZONTAL} {
	    set cur_dir $VERTICAL
	} else {
	    set cur_dir $HORIZONTAL
	}
    } else {
	set cur_dir $HORIZONTAL
    }

    #loading seq from sequence
    set seq_id [read_sequence -sequence $seq -entry $name \
		    -direction $cur_dir]

    set $w.success_h $seq_id

    if {[set $w.success_h] != -1} {
	if {[llength [sequence_names]] < 1} {
	    SeqActivateMenu_Initial
	} elseif {[llength [sequence_names]] == 1} {
	    #SeqActivateMenu_Open1
	    SeqActivateMenu_Open
	} else {
	    SeqActivateMenu_Open
	}
    }

    return $seq_id
}

proc EventHandler {id type args} {
    #puts SIP_EH
    #puts [info level [info level]]
    upvar #0 commn_$id commn
    upvar #0 cursor_h_$id cursor_lookup_h
    upvar #0 cursor_v_$id cursor_lookup_v

    global HORIZONTAL VERTICAL

    global communicating
    if {[info exists communicating($id)]} {

	if {$communicating($id) == 1} { return }
    }

    #puts "id=$id type=$type args=$args"

    global sip_defs
    if {$type == "INIT_SEQUENCE"} {
	global num_event_ids

	#puts INIT_SEQUENCE
        if {![info exists num_event_ids]} {
           set num_event_ids 0
        }

        set idname ${id}_${num_event_ids}
	set seq_id [LoadSequence [keylget args sequence] [keylget args seq_id]]
	set seq_num [get_seq_num $seq_id]

	#puts "seq_num $seq_num seq-id $seq_id"

        global commn_$idname cur_dir
       	upvar #0 commn_$idname commn
#	set commn(cursor_id) [create_cursor -seq_num $seq_num -line_width [keylget sip_defs SIP.CURSOR.LINE_WIDTH] -direction $cur_dir]
	set commn(seq_id) $seq_id
	set commn(command) [keylget args command]

	set commn(reg_id) [seq_sender -rid $idname -seq_id $seq_id]

        incr num_event_ids
        return $idname
    }

    if {$type == "STOP_SEQUENCE"} {
	#puts "Cancelling communication with $id"

	#need to save seq_id before do QUIT because this unsets commn
	set seq_id $commn(seq_id)
	seq_result_update -index $commn(reg_id) -job QUIT

	foreach c [array names cursor_lookup_h] {
	    delete_cursor -seq_num [get_seq_num $seq_id] \
		-id $cursor_lookup_h($c) -private 1
	}

	foreach c [array names cursor_lookup_v] {
	    delete_cursor -seq_num [get_seq_num $seq_id] \
		-id $cursor_lookup_v($c) -private 1
	}

	#check if no more sequences - shutdown sip
	if {[num_sequences] == 0} {
	    ExitSip
	}
	unset commn
        return
    }
    
    if {![info exists commn]} {
	#puts "Attempt to communicate with cancelled id $id"
	return
    }

    if {$type == "CURSOR_NOTIFY"} {

	global dir
	
	set created ""
	set apos [keylget args abspos]
	set job  [keylget args job]
	set cid  [keylget args id]
	set seq_id $commn(seq_id)
	set seq_num [get_seq_num $seq_id]

	#puts "cid=$cid [info exists cursor_lookup_h($cid)] [info exists cursor_lookup_v($cid)]"

	#alternate which cursor is created
	if {![info exists cursor_lookup_h($cid)] && \
		![info exists cursor_lookup_v($cid)]} {

	    if {![info exists dir] || $dir == $VERTICAL} {
		set dir $HORIZONTAL
		set cursor_lookup_h($cid) \
		    [create_cursor \
			 -seq_num $seq_num \
			 -line_width [keylget sip_defs SIP.CURSOR.LINE_WIDTH] \
			 -direction $HORIZONTAL \
			 -private 1]
		#puts "sip cursor H $cid $cursor_lookup_h($cid)"

		lappend created $cursor_lookup_h($cid)
	    } else {
		set dir $VERTICAL
		set cursor_lookup_v($cid) \
		    [create_cursor \
			 -seq_num $seq_num \
			 -line_width [keylget sip_defs SIP.CURSOR.LINE_WIDTH] \
			 -direction $VERTICAL \
			 -private 1]

		#puts "sip cursor V $cid $cursor_lookup_v($cid)"
		lappend created $cursor_lookup_v($cid)
	    }
	} else {
	    #don't create new cursor AND increment
            if {[lsearch $job INCREMENT] != -1} {
		if {[info exists cursor_lookup_h($cid)]} {
		    #puts INCREMENT_H
		    cursor_ref -seq_num $seq_num -id $cursor_lookup_h($cid)\
			-direction $HORIZONTAL -ref 1
		} 
		if {[info exists cursor_lookup_v($cid)]} {
		    #puts INCREMENT_V
		    cursor_ref -seq_num $seq_num -id $cursor_lookup_v($cid)\
			-direction $VERTICAL -ref 1
		}
            }
	}	    

# 	if {![info exists cursor_lookup_h($cid)]} {
# 	    set cursor_lookup_h($cid) \
# 		[create_cursor \
# 		    -seq_num $seq_num \
# 		    -line_width [keylget sip_defs SIP.CURSOR.LINE_WIDTH] \
# 		    -direction $HORIZONTAL \
# 		    -private 1]
# 	    lappend created $cursor_lookup_h($cid)
# 	}
# 	if {![info exists cursor_lookup_v($cid)]} {
# 	    set cursor_lookup_v($cid) \
# 		[create_cursor \
# 		    -seq_num $seq_num \
# 		    -line_width [keylget sip_defs SIP.CURSOR.LINE_WIDTH] \
# 		    -direction $VERTICAL \
# 		    -private 1]
# 	    lappend created $cursor_lookup_v($cid)
# 	}

	if {[lsearch $job DELETE] != -1} {
	    if {[info exists cursor_lookup_h($cid)]} {
		delete_cursor -seq_num $seq_num -id $cursor_lookup_h($cid) \
		    -private 1
		
		#puts "DELETE H -id $cid $cursor_lookup_h($cid)"

	    }
	    if {[info exists cursor_lookup_v($cid)]} { 
		delete_cursor -seq_num $seq_num -id $cursor_lookup_v($cid) \
		    -private 1

		#puts "DELETE V -id $cid $cursor_lookup_v($cid)"

	    }
	} else {
	    #don't delete a cursor AND decrement
	    if {[lsearch $job DECREMENT] != -1} {
		if {[info exists cursor_lookup_h($cid)]} {
		    #puts DECREMENT_H
		    cursor_ref -seq_num $seq_num -id $cursor_lookup_h($cid)\
			-direction $HORIZONTAL -ref -1
		}
		if {[info exists cursor_lookup_v($cid)]} {
		    #puts DECREMENT_V
		    cursor_ref -seq_num $seq_num -id $cursor_lookup_v($cid)\
			-direction $VERTICAL -ref -1
		}
	    }		
	}

	if {[lsearch $job MOVE] != -1} {
	    if {[info exists cursor_lookup_h($cid)]} {
		#puts "cursor_notify -seq_num $seq_num -id $cursor_lookup_h($cid) \
							  -pos $apos -direction $HORIZONTAL"
		cursor_notify -seq_num $seq_num -id $cursor_lookup_h($cid) \
		    -pos $apos -direction $HORIZONTAL
	    }
	    if {[info exists cursor_lookup_v($cid)]} {  
		#puts "cursor_notify -seq_num $seq_num -id $cursor_lookup_v($cid) \
							  -pos $apos -direction $VERTICAL"
		cursor_notify -seq_num $seq_num -id $cursor_lookup_v($cid) \
		    -pos $apos -direction $VERTICAL
}
	}
	if {$created != ""} {
	    set clist [query_cursor -cursorid $created -seq_num $seq_num]
	    
	    #puts "SIP CURSOR $clist [keylget clist refs]"
	    return "$created [keylget clist refs]"
	} else {
	    return $created
	}
    }
}
