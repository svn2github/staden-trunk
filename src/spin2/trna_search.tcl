#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc trnaSearch {} {
    global nip_defs 

    set w .trna_search
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "tRNA gene search"
 
    keylset us RANGE [keylget nip_defs NIP.TRNA.RANGE]

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
	-title [keylget us RANGE.NAME] \
	-update_cmd [list [list seq_range_updates $w.range]]\
	-browse_cmd nip_seq_browser

    strand_both $w.strand
    
    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $w.button -bd 2 -relief groove \
	-ok_command "trnaSearch2 $w $w.range \[strand_get $w.strand\]; destroy $w"\
	-cancel_command "seq_id_destroy $w.range; destroy $w" \
	-help_command "show_help spin {SPIN-TRNA-Search}"

    pack $w.range $w.strand $w.button -fill x
}

proc plot_trna_search {seq_id strand result_id result} {
    global nip_defs tk_utils_defs spin_defs
    global HORIZONTAL SCALE_BAR

    set type [keylget nip_defs TRNA]

    set element_info [create_seq_element $seq_id -1 $type $strand 0 \
	    $HORIZONTAL XY "CANVAS" [keylget spin_defs CONTAINER.TITLE] \
	    [keylget tk_utils_defs ELEMENT.PLOT_WIDTH] \
	    [keylget tk_utils_defs ELEMENT.SINGLE.PLOT_HEIGHT]]
    
    set c_win [keylget element_info container_win]
    set c_id [keylget element_info container_id]
    set e_win [keylget element_info element_win]
    set e_id [keylget element_info element_id]
    set orientation [keylget element_info orientation]
    
    nip_trna_search plot -element $c_win$e_win\
	    -container $c_win\
	    -seq_id $seq_id \
	    -result_id $result_id\
	    -results $result\
	    -container_id [keylget element_info container_id]\
	    -element_id [keylget element_info element_id]\
	    -element_type "CANVAS"\
	    -fill [keylget nip_defs NIP.TRNA.COLOUR]\
	    -tick_ht [keylget nip_defs NIP.TRNA.TICK_HT]\
	    -orientation $orientation

    seqed_element_bindings $c_id $c_win$e_win $e_id
    
    #update result list
    result_list_update $c_win
}

proc trnaSearch2 {t range strand} {
    global nip_defs tk_utils_defs
    global HORIZONTAL SCALE_BAR PROTEIN TOP_S BOTTOM_S

    set seq_id [name_to_seq_id [seq_id_name $range]]
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "trna search" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    SetBusy
    
    if {[expr $strand & $TOP_S]} {
	
	set res [nip_trna_search create -start [seq_id_from $range] \
		       -end [seq_id_to $range] -seq_id $seq_id -strand $TOP_S]
	
	set result_id [lindex $res 0]
	set result [lindex $res 1]
	if {$result_id != -1} {
	    # stop windows from hiding the plot
	    wm withdraw $t
	
	    plot_trna_search $seq_id $TOP_S $result_id $result
	}
    }
    if {[expr $strand & $BOTTOM_S]} {
	
	set res [nip_trna_search create -start [seq_id_from $range] \
		       -end [seq_id_to $range] -seq_id $seq_id -strand $BOTTOM_S]
	
	set result_id [lindex $res 0]
	set result [lindex $res 1]
	if {$result_id != -1} {
	    # stop windows from hiding the plot
	    wm withdraw $t
	    
	    plot_trna_search $seq_id $BOTTOM_S $result_id $result
	}
    }

    ClearBusy
    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
}
