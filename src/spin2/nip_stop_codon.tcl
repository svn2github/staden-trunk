#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc NipStopCodons { } {
    global nip_defs 

    set t .stop_codons

    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Plot stop codons"
  
    #select plot range
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

    keylset us RANGE [keylget nip_defs NIP.STOP_CODON.RANGE]
    seq_id $t.range -range 1 -browse 1 -from 1 -to $seq_length \
	-start_value $seq_start -end_value $seq_end -min_value 1 \
	-default [seq_info $seq_id name]\
 	-update_cmd [list [list seq_range_updates $t.range]]\
	-browse_cmd nip_seq_browser

    #########################################################################
    #strand selection
    strand_both $t.strand

    #########################################################################
    #OK and Cancel buttons
    okcancelhelp $t.ok_cancel \
	    -ok_command "NipStopCodon2 $t $t.range \[strand_get $t.strand\]; destroy $t" \
	    -cancel_command "seq_id_destroy $t.range; destroy $t" \
	    -help_command "show_help spin {SPIN-Stop-Codon-Search}" \
	    -bd 2 \
	    -relief groove

    #final packing
    pack $t.range -fill both
    pack $t.strand -fill both
    pack $t.ok_cancel -fill x
}

proc plot_codons {seq_id strand r_id res col_list} {
    global nip_defs spin_defs tk_utils_defs HORIZONTAL TOP_S BOTTOM_S
    global HORIZONTAL SCALE_BAR

    upvar $r_id result_id $res results

    set type [keylget nip_defs STOPCODON]

    set cnt [array size result_id]

    for {set i 0} {$i < $cnt} {incr i} {

	if {$result_id($i) != -1} {
	    set frame [expr $i + 1]
	    set element_info [create_seq_element $seq_id -1 $type $strand $frame \
		    $HORIZONTAL X "CANVAS" [keylget spin_defs CONTAINER.TITLE] \
		    [keylget tk_utils_defs ELEMENT.PLOT_WIDTH] \
		    [keylget tk_utils_defs ELEMENT.SINGLE.PLOT_HEIGHT]]

	    set c_win [keylget element_info container_win]
	    set c_id [keylget element_info container_id]
	    set e_win [keylget element_info element_win]
	    set e_id [keylget element_info element_id]
	    set row [keylget element_info row]
	    set orientation [keylget element_info orientation]

	    update idletasks

	    nip_stop_codons plot -element $c_win$e_win\
		    -container $c_win\
		    -seq_id $seq_id \
		    -result_id $result_id($i)\
		    -results $results($i)\
		    -container_id [keylget element_info container_id]\
		    -element_id [keylget element_info element_id]\
		    -element_type "CANVAS"\
		    -fill [lindex $col_list $i]\
		    -tick_ht [keylget nip_defs NIP.STOP_CODON.TICK_HT]\
		    -orientation $orientation

	    seqed_element_bindings $c_id $c_win$e_win $e_id
	}
    }
    #update result list
    result_list_update $c_win 

}

proc NipStopCodon2 {t range strand } {
    global PROTEIN nip_defs TOP_S BOTTOM_S

    set seq_id [name_to_seq_id [seq_id_name $range]]
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "stop codons" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    SetBusy

    if {$t != ""} {
	wm withdraw $t
    }

    set col_list "[keylget nip_defs NIP.STOP_CODON.COLOUR.F1] [keylget nip_defs NIP.STOP_CODON.COLOUR.F2] [keylget nip_defs NIP.STOP_CODON.COLOUR.F3]"

    if {[expr $strand & $TOP_S]} {
	set res [nip_stop_codons create\
		-start [seq_id_from $range] \
		-end [seq_id_to $range] \
		-strand $TOP_S -seq_id $seq_id]

	if {[llength $res] == 0} {
	    ClearBusy
	    seq_id_destroy $range
	    return
	}
	set cnt [expr [llength $res] / 2]
	set n 0 
	for {set i 0} {$i < $cnt} {incr i} {
	    set result_id($i) [lindex $res $n]
	    incr n
	    set results($i) [lindex $res $n]
	    incr n
	}

	if {-1 == [plot_codons $seq_id $TOP_S result_id results $col_list]} {
	    ClearBusy
	    seq_id_destroy $range
	    return
	}
    }

    if {[expr $strand & $BOTTOM_S]} {
	set res [nip_stop_codons create\
		-start [seq_id_from $range] \
		-end [seq_id_to $range] \
		-strand $BOTTOM_S -seq_id $seq_id]

	if {[llength $res] == 0} {
	    ClearBusy
	    seq_id_destroy $range
	    return
	}
	set cnt [expr [llength $res] / 2]
	set n 0 
	for {set i 0} {$i < $cnt} {incr i} {
	    set result_id($i) [lindex $res $n]
	    incr n
	    set results($i) [lindex $res $n]
	    incr n
	}
	
	if {-1 == [plot_codons $seq_id $BOTTOM_S result_id results $col_list]} {
	    ClearBusy
	    seq_id_destroy $range
	    return
	}
 
    }
    ClearBusy

    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
}

proc NipStartCodons { } {
    global nip_defs 

    set t .start_codons

    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Plot start codons"
  
    #select plot range
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

    keylset us RANGE [keylget nip_defs NIP.START_CODON.RANGE]
    seq_id $t.range -range 1 -browse 1 -from 1 -to $seq_length \
	-start_value $seq_start -end_value $seq_end -min_value 1 \
	-default [seq_info $seq_id name]\
 	-update_cmd [list [list seq_range_updates $t.range]]\
	-browse_cmd nip_seq_browser

    #########################################################################
    #strand selection
    strand_both $t.strand 

    #########################################################################
    #OK and Cancel buttons
    okcancelhelp $t.ok_cancel \
	    -ok_command "NipStartCodon2 $t $t.range \[strand_get $t.strand\]; destroy $t" \
	    -cancel_command "seq_id_destroy $t.range; destroy $t" \
	    -help_command "show_help spin {SPIN-Start-Codon-Search}" \
	    -bd 2 \
	    -relief groove

    #final packing
    pack $t.range -fill both
    pack $t.strand -fill both
    pack $t.ok_cancel -fill x
}

proc NipStartCodon2 {t range strand} {
    global PROTEIN nip_defs TOP_S BOTTOM_S

    set seq_id [name_to_seq_id [seq_id_name $range]]
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "start codons" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    set col_list "[keylget nip_defs NIP.START_CODON.COLOUR.F1] [keylget nip_defs NIP.START_CODON.COLOUR.F2] [keylget nip_defs NIP.START_CODON.COLOUR.F3]"


    SetBusy
    if {$t != ""} {
	wm withdraw $t
    }

    if {[expr $strand & $TOP_S]} {
	set res [nip_start_codons create\
	    -start [seq_id_from $range] \
	    -end [seq_id_to $range] \
	    -strand $TOP_S -seq_id $seq_id]
    
	if {[llength $res] == 0} {
	    ClearBusy
	    seq_id_destroy $range
	    return
	}
	set cnt [expr [llength $res] / 2]
	set n 0 
	for {set i 0} {$i < $cnt} {incr i} {
	    set result_id($i) [lindex $res $n]
	    incr n
	    set results($i) [lindex $res $n]
	    incr n
	}
	if {-1 == [plot_codons $seq_id $TOP_S result_id results $col_list]} {
	    ClearBusy
	    seq_id_destroy $range
	    return
	}
    }
    if {[expr $strand & $BOTTOM_S]} {
	set res [nip_start_codons create\
	    -start [seq_id_from $range] \
	    -end [seq_id_to $range] \
	    -strand $BOTTOM_S -seq_id $seq_id]
    
	if {[llength $res] == 0} {
	    ClearBusy
	    seq_id_destroy $range
	    return
	}
	set cnt [expr [llength $res] / 2]
	set n 0 
	for {set i 0} {$i < $cnt} {incr i} {
	    set result_id($i) [lindex $res $n]
	    incr n
	    set results($i) [lindex $res $n]
	    incr n
	}
	if {-1 == [plot_codons $seq_id $BOTTOM_S result_id results $col_list]} {
	    ClearBusy
	    seq_id_destroy $range
	    return
	}
    }
    ClearBusy
    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
}


