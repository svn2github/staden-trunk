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
    global nip_defs nsend

    set t .sendto
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Send to"

    seq_id $t.seq_id -range 0 -default [get_active_seq_name] -browse_cmd nip_seq_browser

    frame $t.interps -bd 2 -relief raised
    label $t.interps.title -text "Program name"
    frame $t.interps.lb
    listbox $t.interps.lb.l -height 5 -yscrollcommand "$t.interps.lb.s set"
    scrollbar $t.interps.lb.s -orient vert -command "$t.interps.lb.l yview"

    set interps ""
    set interps_done ""
    set thisapp [tk appname]
    catch {unset nsend};set nsend() ""; unset nsend(); # create a blank array
    foreach i [keylget nip_defs NIP_SEND_TO] {
	set nsend([lindex $i 0]) [lindex $i 1]
    }
    
    foreach i [winfo interps] {
	if {$i == $thisapp} {
	    continue
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
    foreach i [array names nsend] {
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
	-ok_command "SendTo2 $t $t.seq_id $t.interps.lb.l nsend" \
	-cancel_command "destroy $t" \
	-help_command "show_help interface {UI-SendTo}"

    pack $t.interps -side top -fill both
    pack $t.but -side bottom -fill both
}

proc SendTo2 {t s_id lb nsend_p} {
    upvar #0 $nsend_p nsend
    global contig nip_defs HORIZONTAL

    if {[set name [$lb curselection]] == ""} {
	focus $lb
	bell
	return
    }
    set name [$lb get $name]

    # If necessary, start up child process and wait for it to be ready.
    if {[lsearch -exact [winfo interps] $name] == -1} {
	eval exec $nsend($name) &
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

    #create cursors only when necessary, not on sending
    #set commn(cursor_id) [create_cursor -seq_num $seq_num -line_width \
			      [keylget nip_defs NIP.CURSOR.LINE_WIDTH] \
			     -direction $HORIZONTAL]
    seq_sender -rid $rid -seq_id $seq_id

}


proc NipLoadSequence {seq name} {

    #loading seq from sequence
    set seq_id [nip_set_sequence -sequence $seq -entry $name]
    menu_state_set nip_menu 4 .mainwin.menus

    return $seq_id
}

proc EventHandler {id type args} {
    #puts "NIP_EH"

    #puts [info level [info level]]
    upvar #0 commn_$id commn
    upvar #0 cursor_$id cursor_lookup
    
    global communicating HORIZONTAL
    if {[info exists communicating($id)]} {
	if {$communicating($id) == 1} { return }
    }
    
    #puts "NIP id $id type $type args $args"

    global nip_defs
    if {$type == "INIT_SEQUENCE"} {
	global num_event_ids 
	
        if {![info exists num_event_ids]} {
	    set num_event_ids 0
        }

        set idname ${id}_${num_event_ids}
	set seq_id [NipLoadSequence [keylget args sequence] [keylget args seq_id]]
	
        global commn_$idname
	
       	upvar #0 commn_$idname commn
	#set commn(cursor_id) [create_cursor \
	#    -seq_num [get_seq_num $seq_id] \
	#   -line_width \
	#  [keylget nip_defs NIP.CURSOR.LINE_WIDTH]\
	# -direction $HORIZONTAL]

        set commn(seq_id) $seq_id
	set commn(command) [keylget args command]
    
        set commn(reg_id) [seq_sender -rid $idname -seq_id $seq_id]
    
        incr num_event_ids
        return $idname
    }
    
    if {$type == "STOP_SEQUENCE"} {
	#puts "Cancelling communication with $id"
	#delete_cursor -seq_id $commn(seq_id) -id $commn(cursor_id)
	#nip_delete_sequence -seq_id $commn(seq_id)

	#need to save seq_id before do QUIT because this unsets commn
	set seq_id $commn(seq_id)
	set seq_num [get_seq_num $seq_id]
	seq_result_update -index $commn(reg_id) -job QUIT

	foreach c [array names cursor_lookup] {
	    #puts "delete_cursor -seq_num $seq_num \
	    #-id $cursor_lookup($c) -private 1"
	    delete_cursor -seq_num [get_seq_num $seq_id] \
		-id $cursor_lookup($c) -private 1
	}
	#check if no more sequences - shutdown nip
	if {[num_sequences] == 0} {
	    ExitNip
	}
	unset commn
        return
    }

    if {![info exists commn]} {
	#puts "Attempt to communicate with cancelled id $id"
	return
    }

    if {$type == "CURSOR_NOTIFY"} {

	set apos [keylget args abspos]
	set job  [keylget args job]
	set cid  [keylget args id]

	set seq_id $commn(seq_id)
	set seq_num [get_seq_num $seq_id]
	set created ""

        #puts "cid $cid"
	if {![info exists cursor_lookup($cid)]} {
	    set cursor_lookup($cid) \
		[create_cursor \
		     -seq_num $seq_num \
		     -line_width [keylget nip_defs NIP.CURSOR.LINE_WIDTH] \
		     -direction $HORIZONTAL \
		     -private 0]
	    set created $cursor_lookup($cid)
            #puts "CREATED NEW CURSOR $created $cid"
	} else {
	    #don't create new cursor and increment
            if {[lsearch $job INCREMENT] != -1} {
		cursor_ref -seq_num $seq_num -id $cursor_lookup($cid)\
		    -direction $HORIZONTAL -ref 1
            }
	}
	
	if {[lsearch $job DELETE] != -1} {
	    if {[info exists cursor_lookup($cid)]} {
		delete_cursor -seq_num $seq_num -id $cursor_lookup($cid) \
		    -private 1
		#puts "DELETE -id $cid $cursor_lookup($cid)"
	    }
	} elseif {[lsearch $job DECREMENT] != -1} {
	    #don't delete a cursor AND decrement
	    cursor_ref -seq_num $seq_num -id $cursor_lookup($cid)\
		-direction $HORIZONTAL -ref -1
	}
	
	if {[lsearch $job MOVE] != -1} {
	    if {[info exists cursor_lookup($cid)]} {
		#puts "cursor_notify -seq_num $seq_num -id $cursor_lookup($cid) \
							-pos $apos"

		cursor_notify -seq_num $seq_num -id $cursor_lookup($cid)\
		    -pos $apos -direction $HORIZONTAL
	    }
	}
	if {$created != ""} {
	    set clist [query_cursor -cursorid $created -seq_num $seq_num]
	    
	    #puts "NIP CURSOR $clist [keylget clist refs]"
	    return "$created [keylget clist refs]"
	} else {
	    return $created
	}
	#return $created
    }
}
