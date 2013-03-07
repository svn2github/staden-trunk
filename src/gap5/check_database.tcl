# Check Database
proc CheckDatabase {io} {
    global gap5_defs

    set t .checkdb
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Check Database"

    contig_id $t.id -io $io -range 0
    lorf_in $t.infile [keylget gap5_defs CHECK_DATABASE.INFILE] \
	"{contig_id_configure $t.id -state disabled} \
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state normal}" -bd 2 -relief groove
    
#    xradiobox $t.speed \
#	-labeltext "Attention to detail" \
#	-orient horiztonal
#    $t.speed add 1 -text "Minimal (fast)"
#    $t.speed add 2 -text "Thorough (slow)"
#    $t.speed select [keylget gap5_defs CHECK_DATABASE.DETAIL]
    radiolist $t.speed \
	-title "Attention to detail" \
	-bd 0 \
	-relief groove \
	-orient horizontal \
	-default [keylget gap5_defs CHECK_DATABASE.DETAIL] \
	-buttons {{{Minimal (fast)}} {{Thorough (slow)}}}

    if {![$io read_only]} {
	xyn $t.fix \
	    -label "Fix problems?" \
	    -default [keylget gap5_defs CHECK_DATABASE.FIX]
    } else {
	xyn $t.fix \
	    -label "Fix problems?" \
	    -default 0 \
	    -state disabled
    }

    okcancelhelp $t.ok \
	-ok_command "CheckDatabase2 $io $t" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 CheckDatabase" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.speed $t.fix $t.ok -side top -fill x
}

proc CheckDatabase2 {io w} {
    set fix [$w.fix get]
    set level [radiolist_get $w.speed]

    if {[lorf_in_get $w.infile] == 4} {
	set list [list [contig_id_gel $w.id]]
    } elseif {[lorf_in_get $w.infile] == 3} {
	set list "*all*"
    } else {
	set list [lorf_get_list $w.infile]
    }

    destroy $w

    $io flush

    SetBusy

    # Store messages in the log file if we're writing one.
    set prev_vmlog [log_vmessage 1]

    if {$list == "*all*"} {
	$io check $fix $level
    } else {
	vfuncheader "Check Database"
	
	set nc [llength $list]
	set idx 1
	set err 0
	set fixed 0
	foreach c $list {
	    set crec [db_info get_contig_num $io $c]
	    if {$crec == 0} {
		verror "Unknown contig name '$c'"
	    } else {
		vmessage "--Checking contig $c/#$crec ($idx of $nc)"
		set c [$io get_contig $crec]
		foreach {err_inc fixed_inc} [$c check $fix] break
		incr err $err_inc
		incr fixed $fixed_inc
		$c delete
		update
		update idletasks
	    }
	    incr idx
	}

	vmessage "\n*** Total number of errors: $err ***"
	if {$fix} {
	    vmessage "*** Attempted to fix      : $fixed ***"
	}
    }

    log_vmessage $prev_vmlog

    $io flush

    ClearBusy
}

