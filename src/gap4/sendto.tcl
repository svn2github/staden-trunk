global InCommunication
set InCommunication 0

proc SendTo {io} {
    global gap_defs gsend

    set t .sendto
    if {[xtoplevel $t -resizable 0] == ""} return

    wm title $t "Send to"

    contig_id $t.id -io $io -range 0 -command "SendTo2 $io $t $t.id"

    frame $t.interps -bd 2 -relief raised
    label $t.interps.title -text "Program name"
    frame $t.interps.lb
    listbox $t.interps.lb.l -height 5 -yscrollcommand "$t.interps.lb.s set"
    scrollbar $t.interps.lb.s -orient vert -command "$t.interps.lb.l yview"

    set interps ""
    set interps_done ""
    set thisapp [tk appname]
    catch {unset gsend};set gsend() ""; unset gsend(); # create a blank array
    foreach i [keylget gap_defs GAP_SEND_TO] {
	set gsend([lindex $i 0]) [lindex $i 1]
    }
    foreach i [winfo interps] {
	if {$i == $thisapp} {
	    continue
	}
	if {[info exists gsend($i)]} {
	    lappend interps_done $i
	}
	set cval [catch {send $i {info commands EventHandler}} send_ret]
	if {$cval == 0 && $send_ret == "EventHandler"} {
	    $t.interps.lb.l insert end $i
	    lappend interps $i
	}
    }
    foreach i [array names gsend] {
	if {[lsearch -exact $interps_done $i] == -1} {
	    $t.interps.lb.l insert end $i
	    lappend interps $i
        }
    }

    pack $t.interps.lb.l -side left -fill both -expand 1
    pack $t.interps.lb.s -side left -fill both
    pack $t.interps.title $t.interps.lb -side top -fill both

    okcancelhelp $t.but \
	-ok_command "SendTo2 $io $t $t.id $t.interps.lb.l gsend" \
	-cancel_command "destroy $t" \
	-help_command "show_help interface {UI-SendTo}"

    pack $t.id $t.interps -side top -fill both
    pack $t.but -side bottom -fill both
}

proc SendTo2 {io t id lb gsend_p} {
    upvar #0 $gsend_p gsend

    if {[set reading [contig_id_gel $id]] == ""} {
	focus $id
	bell
	return
    }

    if {[set name [$lb curselection]] == ""} {
	focus $lb
	bell
	return
    }
    set name [$lb get $name]

    # If necessary, start up child process and wait for it to be ready.
    if {[lsearch -exact [winfo interps] $name] == -1} {
	eval exec $gsend($name) &
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

    destroy $t
    SetContigGlobals $io $reading

    set contig [db_info get_contig_num $io $reading]
    # puts "Reading = $reading, contig = $contig"

    set seq [calc_consensus -io $io -contigs "{=$contig}"]

    set l ""
    keylset l sequence $seq
    keylset l command [list send [tk appname]]
    keylset l from [tk appname]
    #keylset l seq_id $contig
    keylset l seq_id [left_gel $io $contig]

    set rid [send $name EventHandler [list [tk appname]=$contig] \
	INIT_SEQUENCE $l]
    # puts ">>>>>>>> SENT [tk appname]=$contig    RID = $rid <<<<<<<<<<<<"

    #set cid [create_cursor -io $io -cnum $contig]

    global commn_${rid}
    upvar #0 commn_${rid} commn
    set commn(remote_id) $rid

    #set commn(cursor_id) $cid

    set commn(contig) $contig
    set commn(command) "send {$name}"
    set commn(remote_name) $name
    set commn(from) [tk appname]

    set reg [contig_register \
	-io $io \
	-contig $contig \
	-command [list SendToCallback $io $rid] \
	-flags {REQUIRED OPS CURSOR_NOTIFY}]
    set commn(reg_id) $reg
}

proc SendToCallback {io rid type id args} {
    global InCommunication
    upvar #0 commn_${rid} commn
    upvar #0 cursor_$rid cursor

    #puts "GAP callback id=$id type=$type args=$args"

    #puts "[clock format [clock clicks]]:START [info level [info level]]"

    if {$type == "REGISTER"} {
	set c [keylget args contig]
    }

    if {$type == "CURSOR_NOTIFY"} {
	set InCommunication 1

	if {[catch {set cursors [eval $commn(command) \
		[list EventHandler [list $rid] CURSOR_NOTIFY $args]]}] == 0} {

	    set gap_curs [keylget args id]
	    set cursor([lindex $cursors 0]) $gap_curs
	    
	    #puts "gap cursor [lindex $cursors 0] $cursor([lindex $cursors 0])"
	}
        set InCommunication 0
    }

    if {$type == "QUERY_NAME"} {
	return "Send to $commn(remote_name), [left_gel $io $commn(contig)]"
    }

    if {$type == "QUERY_PARAMS"} {
	return "unknown"
    }

    if {$type == "GET_OPS"} {
	return "Information Remove"
    }

    if {$type == "INVOKE_OP"} {
	SendToInvokeOp $io $rid $type $id $args
    }

    if {$type == "DELETE" || $type == "QUIT"} {
	set InCommunication 1
	catch {eval $commn(command) [list EventHandler [list $rid] \
		STOP_SEQUENCE]}
	contig_deregister -io $io -id $commn(reg_id)
	#delete_cursor -io $io -cnum $commn(contig) -id $commn(cursor_id)
	unset commn
	set InCommunication 0
    }
    # puts "[clock format [clock clicks]]:END [info level [info level]]"
}

proc SendToInvokeOp {io rid type id ops} {
    # puts "[clock format [clock clicks]]:START [info level [info level]]"

    global InCommunication
    upvar #0 commn_${rid} commn

    set op [keylget ops op]

    if {$op == 0} {
	# Information
	vmessage "SendTo information"
	vmessage "    Contig [left_gel $io $commn(contig)]"
	vmessage "    Command \"$commn(command)\""

    } elseif {$op == 1} {
	# Remove
	set InCommunication 1
	catch [list eval $commn(command) [list EventHandler [list $rid] \
		STOP_SEQUENCE]]
	contig_deregister -io $io -id $commn(reg_id)
	#delete_cursor -io $io -cnum $commn(contig) -id $commn(cursor_id)
	unset commn
	set InCommunication 0
   }

   # puts "[clock format [clock clicks]]:END [info level [info level]]"
}

proc EventHandler {id type args} {
    #puts "gap EH"
    global InCommunication
    global io num_event_ids
    upvar #0 commn_$id commn
    upvar #0 cursor_$id cursor

    #puts "[clock format [clock clicks]]:START [info level [info level]]"

    # Loop avoidance
    if {$InCommunication} {
	return
    }

    #puts "id=$id type=$type args=$args"

    if {$type == "INIT_SEQUENCE"} {
	if {![info exists num_event_ids]} {
	   set num_event_ids 0
	}

	set idname ${id}/[tk appname]=${num_event_ids}
        global commn_$idname
       	upvar #0 commn_$idname commn

	set commn(contig) [db_info get_contig_num $io [keylget args seq_id]]
	
	if {$commn(contig) == -1} {
	    verror ERR_WARN sendto "invalid contig name"
	    bell 
	    return
	}

	#set commn(contig) [keylget args seq_id]; #FIXME
#	set cid [create_cursor -io $io -cnum $commn(contig)]
	set commn(remote_id) $id
#	set commn(cursor_id) $cid
	set commn(command) [keylget args command]
	set commn(from) [tk appname]
	set commn(remote_name) [keylget args from]
	set reg [contig_register \
	    -io $io \
	    -contig $commn(contig) \
	    -command [list SendToCallback $io $idname] \
	    -flags {REQUIRED OPS CURSOR_NOTIFY}]
        set commn(reg_id) $reg

	incr num_event_ids
	#puts "returning $idname"
	return $idname
    }

    if {$type == "STOP_SEQUENCE"} {
#	delete_cursor \
	    -io $io \
	    -cnum [set commn(contig)] \
	    -id [set commn(cursor_id)]
	contig_deregister -io $io -id $commn(reg_id)
	unset commn
	return
    }

    if {![info exists commn]} {
	#puts "Attempt to communicate with cancelled id $id"
	return
    }

    if {$type == "CURSOR_NOTIFY"} {
	set contig $commn(contig)
	set pos  [keylget args pos]
	set seq  [keylget args seq]
	set apos [keylget args abspos]
	set job  [keylget args job]
	set cid  [keylget args id]
	set created ""

	#puts "cid $cid"
	if {![info exists cursor($cid)]} {
	    set cursor($cid) [create_cursor -io $io -cnum $contig -private 1]
	    #puts "GAP create cursor $cid $cursor($cid)"
	    set created $cursor($cid)
	} else {
	    if {[lsearch $job INCREMENT] != -1} {
		cursor_ref -io $io -cnum $contig -id $cursor($cid) -ref 1
	    }
	    
	    if {[lsearch $job DECREMENT] != -1} {
		cursor_ref -io $io -cnum $contig -id $cursor($cid) -ref -1
	    }
	}
	keylset args id $cursor($cid)

	#puts "$job: contig_notify -io $io -type $type -cnum $contig -args $args"
	if {[lsearch $job MOVE] != -1} {
	    #puts "contig_notify -io $io -type $type -cnum $contig -args $args"
	    contig_notify -io $io -type $type -cnum $contig -args $args
	}

	if {$created != ""} {
	    set clist [query_cursor -io $io -cursorid $cursor($cid)]
	    
	    #puts "GAP CURSOR $clist [keylget clist refs]"
	    return "$created [keylget clist refs]"
	} else {
	    return $created
	}
    }
    # puts "[clock format [clock clicks]]:END [info level [info level]]"
}




