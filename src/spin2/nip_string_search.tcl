#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc NipStringSearch { } {
    global nip_defs

    set w .nip_search
    if {[SeqedSearchDialog $w "NipStringSearch2 $w $w.range; seq_id_destroy $w.range; destroy $w" "destroy $w" "show_help spin {SPIN-String-Search}"] == -1} {
	return
    }
    
    keylset us RANGE [keylget nip_defs NIP.STRING_SEARCH.RANGE]

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

    seq_id $w.range -range 1 -browse 1 -from 1 -to $seq_length \
	-start_value $seq_start -end_value $seq_end -min_value 1 \
	-default [seq_info $seq_id name]\
	-update_cmd [list [list seq_range_updates $w.range]]\
	-browse_cmd nip_seq_browser

    pack $w.range -before $w.strand -fill both
    frame $w.separator -bd 2 -relief raised -height 2
    pack $w.separator -before $w.strand -side top -fill x -padx 10 -pady 5
}

proc plot_string {seq_id strand result_id results} {
    global nip_defs tk_utils_defs spin_defs HORIZONTAL

    if {$result_id == -1} {
	return
    }

    set type [keylget nip_defs STRINGSEARCH]
    set frame 0

    set element_info [create_seq_element $seq_id -1 $type $strand $frame \
	    $HORIZONTAL X "CANVAS" [keylget spin_defs CONTAINER.TITLE] \
	    [keylget tk_utils_defs ELEMENT.PLOT_WIDTH] \
	    [keylget tk_utils_defs ELEMENT.SINGLE.PLOT_HEIGHT]]
    
    set c_win [keylget element_info container_win]
    set c_id [keylget element_info container_id]
    set e_win [keylget element_info element_win]
    set e_id [keylget element_info element_id]
    set orientation [keylget element_info orientation]
    
    update idletasks
    nip_string_search plot -element $c_win$e_win\
	    -container $c_win\
	    -seq_id $seq_id \
	    -result_id $result_id\
	    -results $results\
	    -container_id [keylget element_info container_id]\
	    -element_id [keylget element_info element_id]\
	    -element_type "CANVAS"\
	    -fill [keylget nip_defs NIP.STRING_SEARCH.COLOUR] \
	    -width 2 \
	    -tick_ht [keylget nip_defs NIP.STRING_SEARCH.TICK_HT]\
	    -orientation $orientation

    
    seqed_element_bindings $c_id $c_win$e_win $e_id

    #update result list
    result_list_update $c_win 
}

proc NipStringSearch2 {t range} {
    global nip_defs tk_utils_defs
    global HORIZONTAL SCALE_BAR PROTEIN TOP_S BOTTOM_S

    set strand [SeqedSearch_strand $t]
    set match [SeqedSearch_match $t]
    set string [SeqedSearch_string $t]
    set use_iub [SeqedSearch_use_iub $t]

    set seq_id [name_to_seq_id [seq_id_name $range]]
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "string search" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    SetBusy
    if {[expr $strand & $TOP_S]} {
	set res [nip_string_search create -seq_id $seq_id \
		-start [seq_id_from $range] \
		-end [seq_id_to $range]\
		-strand $TOP_S -min_pmatch $match -string $string\
		-use_iub $use_iub]

	#remove id from start of list
	set result_id [lindex $res 0]
	set results [lindex $res 1]
	
	# stop windows from hiding the plot
	wm withdraw $t
	
	plot_string $seq_id $TOP_S $result_id $results
	ClearBusy
	
    }
    if {[expr $strand & $BOTTOM_S]} {
	set res [nip_string_search create -seq_id $seq_id \
		-start [seq_id_from $range] \
		-end [seq_id_to $range]\
		-strand $BOTTOM_S -min_pmatch $match -string $string\
		-use_iub $use_iub]
	#remove id from start of list
	set result_id [lindex $res 0]
	set results [lindex $res 1]
	
	# stop windows from hiding the plot
	wm withdraw $t
	
	plot_string $seq_id $BOTTOM_S $result_id $results
	ClearBusy
    }

    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
}

