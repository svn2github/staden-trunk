proc TranslateFT { } {
    global nip_defs 

    set t .translateft
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Find ORF: write as feature table"
    
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
    keylset us RANGE [keylget nip_defs NIP.TRANSLATE_FT.RANGE]
    seq_id $t.range -range 1 -browse 1 -from 1 -to $seq_length \
	-start_value $seq_start -end_value $seq_end -min_value 1 \
	-default [seq_info $seq_id name] \
	-update_cmd [list [list seq_range_updates $t.range]]\
	-browse_cmd nip_seq_browser

    #########################################################################
    #strand selection
    strand_both $t.strand

    #########################################################################
    #minimum orf
    keylset orf MIN_ORF [keylget nip_defs NIP.TRANSLATE_FT.MIN_ORF]
    scalebox $t.min_orf -title [keylget orf MIN_ORF.NAME] \
	-orient horizontal \
	-from [keylget orf MIN_ORF.MIN] \
	-to [keylget orf MIN_ORF.MAX] \
	-default [keylget orf MIN_ORF.VALUE] \
	-width 5 \
	-type CheckInt

    #########################################################################
    #OK and Cancel buttons
    okcancelhelp $t.ok_cancel \
	    -ok_command "TranslateFT2 $t.range $t.strand $t.min_orf; destroy $t" \
	    -cancel_command "seq_id_destroy $t.range; destroy $t" \
	    -help_command "show_help spin {SPIN-Open-Reading-Frames}" \
	    -bd 2 \
	    -relief groove

    #final packing
    pack $t.range -fill both
    pack $t.strand -fill both
    pack $t.min_orf -fill both
    pack $t.ok_cancel -fill x
}

proc TranslateFT2 {range strand min_orf} {
    global PROTEIN

    set seq_id [name_to_seq_id [seq_id_name $range]]
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "open reading frames" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    SetBusy

    translate_orf_to_feature_table \
	-seq_id $seq_id \
	-start [seq_id_from $range] \
	-end [seq_id_to $range] \
	-strand [strand_get $strand]\
	-min_orf [scalebox_get $min_orf]

    ClearBusy
    seq_id_destroy $range
    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
}


proc TranslateFasta { } {
    global nip_defs 

    set t .plot_base_comp
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Find ORF: write protein as fasta file"
    
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
    keylset us RANGE [keylget nip_defs NIP.TRANSLATE_FASTA.RANGE]
    seq_id $t.range -range 1 -browse 1 -from 1 -to $seq_length \
	-start_value $seq_start -end_value $seq_end -min_value 1 \
	-default [seq_info $seq_id name] \
	-update_cmd [list [list seq_range_updates $t.range]]\
	-browse_cmd nip_seq_browser

    #########################################################################
    #strand selection
    strand_both $t.strand

    #########################################################################
    #minimum orf
    keylset orf MIN_ORF [keylget nip_defs NIP.TRANSLATE_FASTA.MIN_ORF]
    scalebox $t.min_orf -title [keylget orf MIN_ORF.NAME] \
	-orient horizontal \
	-from [keylget orf MIN_ORF.MIN] \
	-to [keylget orf MIN_ORF.MAX] \
	-default [keylget orf MIN_ORF.VALUE] \
	-width 5 \
	-type CheckInt

    #########################################################################
    #output filename
    
    keylset fn FILENAME [keylget nip_defs NIP.TRANSLATE_FASTA.FILENAME]
    entrybox $t.filename \
	-title   [keylget fn FILENAME.NAME] \
	-default [keylget fn FILENAME.VALUE] \
	-type CheckOutput

    #########################################################################
    #OK and Cancel buttons
    okcancelhelp $t.ok_cancel \
	    -ok_command "TranslateFasta2 $t $t.range $t.strand $t.min_orf $t.filename" \
	    -cancel_command "seq_id_destroy $t.range; destroy $t" \
	    -help_command "show_help spin {SPIN-Open-Reading-Frames}" \
	    -bd 2 \
	    -relief groove

    #final packing
    pack $t.range -fill both
    pack $t.strand -fill both
    pack $t.min_orf -fill both
    pack $t.filename -fill both
    pack $t.ok_cancel -fill x

}

proc TranslateFasta2 {t range strand min_orf filename} {
    global PROTEIN

    set seq_id [name_to_seq_id [seq_id_name $range]]
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "open reading frames" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }


    set fn [entrybox_get $filename]
    if {$fn == ""} {
	raise $t
	return
    }

    SetBusy

    translate_orf_to_fasta \
	-seq_id $seq_id \
	-start [seq_id_from $range] \
	-end [seq_id_to $range] \
	-strand [strand_get $strand]\
	-min_orf [scalebox_get $min_orf]\
	-filename $fn

    ClearBusy
    seq_id_destroy $range
    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
    destroy $t
}
