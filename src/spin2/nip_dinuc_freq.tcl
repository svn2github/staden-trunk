proc CountDinucFreq { } {
    global nip_defs

    set t .dinuc_freq
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Count dinucleotide frequencies"
    
    set seq_id [get_active_seq_id] 
    global $seq_id.start $seq_id.end

    set seq_length [seq_info $seq_id length] 
    set seq_start [seq_info $seq_id start] 
    set seq_end [seq_info $seq_id end] 

    if {[info exists $seq_id.start]} {
	set seq_start [set $seq_id.start]
    }    
    if {[info exists $seq_id.end]} {
	set seq_end [set $seq_id.end]
    }
    keylset us RANGE [keylget nip_defs NIP.DINUC_FREQ.RANGE]
    seq_id $t.range -range 1 -browse 1 -from 1 -to $seq_length \
	-start_value $seq_start -end_value $seq_end -min_value 1 \
	-default [seq_info $seq_id name] \
	-update_cmd [list [list seq_range_updates $t.range]]\
	-browse_cmd nip_seq_browser

    #########################################################################
    #OK and Cancel buttons
    okcancelhelp $t.ok_cancel \
	    -ok_command "CountDinucFreq2 $t.range; destroy $t" \
	    -cancel_command "seq_id_destroy $t.range; destroy $t" \
	    -help_command "show_help spin {SPIN-Dinucleotide-Freq}" \
	    -bd 2 \
	    -relief groove

    #final packing
    pack $t.range -fill both
    pack $t.ok_cancel -side left -fill both -expand 1
}

proc CountDinucFreq2 {range } {
    global PROTEIN

    set seq_id [name_to_seq_id [seq_id_name $range]]
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "dinucleotide frequencies" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    SetBusy

    count_dinuc_freq \
	-seq_id $seq_id \
	-start [seq_id_from $range] \
	-end [seq_id_to $range]

    ClearBusy

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
    seq_id_destroy $range
}
