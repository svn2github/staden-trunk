#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc SpliceSearch {} {
    global nip_defs

    set w .splice_search
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Search for splice junctions"

    keylset us RANGE [keylget nip_defs NIP.SPLICE.RANGE]
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
	-title [keylget us RANGE.NAME]\
	-update_cmd [list [list seq_range_updates $w.range]]\
	-browse_cmd nip_seq_browser

    strand_both $w.strand

    #get donor weight matrix file
    keylset d DONOR [keylget nip_defs NIP.SPLICE.DONOR]
    eval getFname $w.donor [list [keylget d DONOR.NAME]] load {} \
	    "-default [keylget d DONOR.VALUE]"
    [entrybox_path $w.donor.entry] xview [string last / [keylget d DONOR.VALUE]]

    #get acceptor weight matrix file
    keylset d ACCEPTOR [keylget nip_defs NIP.SPLICE.ACCEPTOR]
    eval getFname $w.acceptor [list [keylget d ACCEPTOR.NAME]] load {} \
	    "-default [keylget d ACCEPTOR.VALUE]"
    [entrybox_path $w.acceptor.entry] xview [string last / [keylget d ACCEPTOR.VALUE]]

    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $w.button -bd 2 -relief groove \
	-ok_command "SpliceSearch2 $w $w.range \[strand_get $w.strand\] $w.donor $w.acceptor; destroy $w"\
	-cancel_command "seq_id_destroy $w.range; destroy $w" \
	-help_command "show_help spin {SPIN-Splice-Site-Search}"

    pack $w.range $w.strand $w.donor $w.acceptor $w.button -fill x
}

proc plot_splice {seq_id strand r_id res} {
    global nip_defs tk_utils_defs spin_defs
    global HORIZONTAL SCALE_BAR

    upvar $r_id result_id $res results
    set type [keylget nip_defs SPLICE]

    set col_list "[keylget nip_defs NIP.SPLICE.COLOUR.F1] [keylget nip_defs NIP.SPLICE.COLOUR.F2] [keylget nip_defs NIP.SPLICE.COLOUR.F3] [keylget nip_defs NIP.SPLICE.COLOUR.F1] [keylget nip_defs NIP.SPLICE.COLOUR.F2] [keylget nip_defs NIP.SPLICE.COLOUR.F3]"

    set element_info [create_seq_element $seq_id -1 $type $strand 0 \
	    $HORIZONTAL XY "CANVAS" [keylget spin_defs CONTAINER.TITLE] \
	    [keylget tk_utils_defs ELEMENT.PLOT_WIDTH] \
	    [keylget tk_utils_defs ELEMENT.SINGLE.PLOT_HEIGHT]]
    
    set c_win [keylget element_info container_win]
    set c_id [keylget element_info container_id]
    set e_win [keylget element_info element_win]
    set e_id [keylget element_info element_id]
    set orientation [keylget element_info orientation]

    set cnt [array size result_id]
    for {set i 0} {$i < $cnt} {incr i} {
	if {$result_id($i) != -1} {
	    set frame [expr $i + 1]

	    update idletasks

	    splice_search plot -element $c_win$e_win\
		    -container $c_win\
		    -seq_id $seq_id \
		    -result_id $result_id($i)\
		    -results $results($i)\
		    -container_id [keylget element_info container_id]\
		    -element_id [keylget element_info element_id]\
		    -element_type "CANVAS"\
		    -fill [lindex $col_list $i]\
		    -tick_ht [keylget nip_defs NIP.SPLICE.TICK_HT]\
		    -orientation $orientation

	    seqed_element_bindings $c_id $c_win$e_win $e_id
	}
    }
    set brief $c_win[keylget tk_utils_defs CONTAINER.BRIEF.WIN]
    #$c_win$e_win bind S <Any-Motion> "highlight_line %W $brief"

    #update result list
    result_list_update $c_win
}

proc SpliceSearch2 {t range strand donor acceptor} {
    global PROTEIN TOP_S BOTTOM_S

    set seq_id [name_to_seq_id [seq_id_name $range]]
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "splice search" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    set donor [getFname_in_name $donor]
    if {$donor == ""} {
	bell
	raise $w
	return
    }
    set acceptor [getFname_in_name $acceptor]
    if {$acceptor == ""} {
	bell
	raise $w
	return
    }

    SetBusy

    if {[expr $strand & $TOP_S]} {
	set res [splice_search create \
		-seq_id $seq_id -start [seq_id_from $range] \
		-end [seq_id_to $range]\
		-donor $donor\
		-acceptor $acceptor\
		-strand $TOP_S]

	if {[llength $res] == 0} {
	    ClearBusy
	    seq_id_destroy $range
	    return
	} else {
	    set cnt [expr [llength $res] / 2]
	}
	set n 0 
	for {set i 0} {$i < $cnt} {incr i} {
	    set result_id($i) [lindex $res $n]
	    incr n
	    set results($i) [lindex $res $n]
	    incr n
	}
	# stop windows from hiding the plot
	wm withdraw $t
	
	plot_splice $seq_id $TOP_S result_id results
    }

    if {[expr $strand & $BOTTOM_S]} {
	set res [splice_search create \
		-seq_id $seq_id -start [seq_id_from $range] \
		-end [seq_id_to $range]\
		-donor $donor\
		-acceptor $acceptor\
		-strand $BOTTOM_S]
	
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

	# stop windows from hiding the plot
	wm withdraw $t
	
	plot_splice $seq_id $BOTTOM_S result_id results
    }

    ClearBusy
    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
}
