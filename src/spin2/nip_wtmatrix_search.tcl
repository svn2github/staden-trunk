#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc NipWtMatrixSearch { } {
    global nip_defs

    set w .nip_matrix
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Search using weight matrix"

    keylset us RANGE [keylget nip_defs NIP.WTMATRIX_SEARCH.RANGE]

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

    #strand
    strand_both $w.strand

    #get weight matrix file
    keylset d MATRIX [keylget nip_defs NIP.WTMATRIX_SEARCH.MATRIX]
    eval getFname $w.matrix [list [keylget d MATRIX.NAME]] load {}

    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $w.button -bd 2 -relief groove \
	-ok_command "NipWtMatrixSearch2 $w $w.range $w.matrix \[strand_get $w.strand\]; destroy $w"\
	-cancel_command "seq_id_destroy $w.range; destroy $w" \
	-help_command "show_help spin {SPIN-Weight-Matrix-Search}"

    pack $w.range $w.matrix $w.strand $w.button -fill both
}

proc plot_wtmatrix {seq_id strand result_id results} {
    global nip_defs tk_utils_defs spin_defs
    global HORIZONTAL SCALE_BAR

    if {$result_id == "-1"} {
	return
    }
    set type [keylget nip_defs WTMATRIXSEARCH]
    set element_info [create_seq_element $seq_id -1 $type $strand 0 \
	    $HORIZONTAL XY "CANVAS" [keylget spin_defs CONTAINER.TITLE] \
	    [keylget tk_utils_defs ELEMENT.PLOT_WIDTH] \
	    [keylget tk_utils_defs ELEMENT.SINGLE.PLOT_HEIGHT]]
    
    set c_win [keylget element_info container_win]
    set c_id [keylget element_info container_id]
    set e_win [keylget element_info element_win]
    set e_id [keylget element_info element_id]
    set orientation [keylget element_info orientation]

    nip_wtmatrix_search plot -element $c_win$e_win\
		    -container $c_win\
		    -seq_id $seq_id \
		    -result_id $result_id\
		    -results $results\
		    -container_id [keylget element_info container_id]\
		    -element_id [keylget element_info element_id]\
		    -element_type "CANVAS"\
		    -fill [keylget nip_defs NIP.WTMATRIX_SEARCH.COLOUR]\
		    -tick_ht [keylget nip_defs NIP.WTMATRIX_SEARCH.TICK_HT]\
		    -orientation $orientation
    
    seqed_element_bindings $c_id $c_win$e_win $e_id

    set brief $c_win[keylget tk_utils_defs CONTAINER.BRIEF.WIN]
    #$c_win$e_win bind S <Any-Motion> "highlight_line %W $brief"

    #update result list
    result_list_update $c_win
}

proc NipWtMatrixSearch2 {t range matrix strand} {
    global nip_defs HORIZONTAL PROTEIN TOP_S BOTTOM_S

    set matrix [getFname_in_name $matrix]
    if {$matrix == ""} {
	bell
	raise $t
	return
    }

    set seq_id [name_to_seq_id [seq_id_name $range]]
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "weight matrix search" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    set type [keylget nip_defs WTMATRIXSEARCH]

    SetBusy
    if {[expr $strand & $TOP_S]} {
	set res [nip_wtmatrix_search create -seq_id $seq_id \
		-start [seq_id_from $range] \
		-end [seq_id_to $range]\
		-wt_matrix $matrix\
		-strand $TOP_S]

	set result_id [lindex $res 0]
	set result [lindex $res 1]
	# stop windows from hiding the plot
	wm withdraw $t
	
	plot_wtmatrix $seq_id $TOP_S $result_id $result
    }
    
    if {[expr $strand & $BOTTOM_S]} {
	set res [nip_wtmatrix_search create -seq_id $seq_id \
		-start [seq_id_from $range] \
		-end [seq_id_to $range]\
		-wt_matrix $matrix\
		-strand $BOTTOM_S]

	set result_id [lindex $res 0]
	set result [lindex $res 1]
	# stop windows from hiding the plot
	wm withdraw $t
	
	plot_wtmatrix $seq_id $BOTTOM_S $result_id $result
    }

    ClearBusy

    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
}

