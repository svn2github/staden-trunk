#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc plot_global_align {seq_id_h seq_id_v res_id s_array num_seq_array res strand} {
    global sip_defs tk_utils_defs spin_defs
    global DNA PROTEIN HORIZONTAL VERTICAL ZOOM_SCALE TOP_S BOTH

    upvar $res_id result_id $s_array seq_array $res results
    
    set type [keylget sip_defs SIMILARSPANS]
    set frame 0

    set element_info [create_seq_element $seq_id_h $seq_id_v $type $strand $frame \
	    $BOTH XY "CANVAS" [keylget spin_defs CONTAINER.TITLE] \
	    [keylget tk_utils_defs DOT.PLOT_WIDTH] \
 	    [keylget tk_utils_defs DOT.PLOT_HEIGHT]]

    set c_win [keylget element_info container_win]
    set e_win [keylget element_info element_win]
    set e_id [keylget element_info element_id]
    set row [keylget element_info row]

    update idletasks

    for {set i 0} {$i < $num_seq_array} {incr i} { 

	#don't want to plot results of no matches
	if {$result_id($i) == -1} {
	    continue
	}

	sip_global_align plot -element $c_win$e_win\
		-container $c_win\
		-seq_id_h $seq_array($i,h) \
		-seq_id_v $seq_array($i,v) \
		-result_id $result_id($i)\
		-results $results($i)\
		-container_id [keylget element_info container_id]\
		-element_id [keylget element_info element_id]\
		-element_type "CANVAS"
    }

    seqpair_element_bindings $c_win$e_win $e_id
    
    #update result list
    result_list_update $c_win 
}

proc AlignSeqs2 {t range_h range_v strand match mismatch start_gap cont_gap} {
    global HORIZONTAL VERTICAL DNA PROTEIN TOP_S BOTTOM_S
    global sip_defs tk_utils_defs

    #must check here if sequence name and range are OK
    if {![CheckSeq $range_h]} {
	raise $t
	tkwait variable wait_forever
    }
    if {![CheckSeq $range_v]} {
	raise $t
	tkwait variable wait_forever
    }
    
    set line_width [keylget sip_defs SIP.CURSOR.LINE_WIDTH]
 
    set seq_id_h [name_to_seq_id [seq_id_name $range_h]]
    set seq_id_v [name_to_seq_id [seq_id_name $range_v]]

    set start_h [seq_id_from $range_h]
    set end_h [seq_id_to $range_h]
    set start_v [seq_id_from $range_v]
    set end_v [seq_id_to $range_v]

    CheckSequenceTypes seq_id_h seq_id_v type_h type_v $start_h $end_h $start_v $end_v

    CreateSeqArray $seq_id_h $seq_id_v $start_h $end_h $start_v $end_v seq_array num_seq_array

    if {$type_h == $PROTEIN || $type_v == $PROTEIN} {
	set m_score 0
	set mm_score 0
    } else {
	set m_score [entrybox_get $match]
	set mm_score [entrybox_get $mismatch]
    }
    set st_gap [entrybox_get $start_gap]
    set ct_gap [entrybox_get $cont_gap]

    SetBusy
    set final_id -1

    if {[expr $strand & $TOP_S]} {
	for {set i 0} {$i < $num_seq_array} {incr i} { 
	    set res [sip_global_align create \
		    -seq_id_h $seq_array($i,h) \
		    -seq_id_v $seq_array($i,v) \
		    -start_h $seq_array($i,start_h) \
		    -end_h $seq_array($i,end_h) \
		    -start_v $seq_array($i,start_v) \
		    -end_v $seq_array($i,end_v) \
		    -strand $TOP_S\
		    -match $m_score \
		    -mismatch $mm_score \
		    -start_gap $st_gap \
		    -cont_gap $ct_gap]   
	
	    set result_id($i) [lindex $res 0]
	    set results($i) [lindex $res 1]
	    
	    if {$result_id($i) != -1} {
		set final_id 0
	    }
	}
	#if no results have been found, return
	if {$final_id == -1} {
	    return
	}
    
	
	# stop windows from hiding the plot
	wm withdraw $t
	
	plot_global_align $seq_id_h $seq_id_v result_id seq_array $num_seq_array results $TOP_S
    }
     if {[expr $strand & $BOTTOM_S]} {
	for {set i 0} {$i < $num_seq_array} {incr i} { 
	    set res [sip_global_align create \
		    -seq_id_h $seq_array($i,h) \
		    -seq_id_v $seq_array($i,v) \
		    -start_h $seq_array($i,start_h) \
		    -end_h $seq_array($i,end_h) \
		    -start_v $seq_array($i,start_v) \
		    -end_v $seq_array($i,end_v) \
		    -strand $BOTTOM_S\
		    -match $m_score \
		    -mismatch $mm_score \
		    -start_gap $st_gap \
		    -cont_gap $ct_gap]   
	
	    set result_id($i) [lindex $res 0]
	    set results($i) [lindex $res 1]
	    
	    if {$result_id($i) != -1} {
		set final_id 0
	    }
	}
	#if no results have been found, return
	if {$final_id == -1} {
	    return
	}
    
	
	# stop windows from hiding the plot
	wm withdraw $t
	
	plot_global_align $seq_id_h $seq_id_v result_id seq_array $num_seq_array results $BOTTOM_S
    }
    

    ClearBusy
    sequence_list_update
    seq_id_destroy $range_h
    seq_id_destroy $range_v
    
     for {set i 0} {$i < $num_seq_array} {incr i} { 
	global $seq_array($i,h).start $seq_array($i,h).end
	global $seq_array($i,v).start $seq_array($i,v).end
	set $seq_array($i,h).start $seq_array($i,start_h)
	set $seq_array($i,h).end $seq_array($i,end_h)
	set $seq_array($i,v).start $seq_array($i,start_v) 
	set $seq_array($i,v).end $seq_array($i,end_v)
    }
}

##############################################################################
proc AlignUpdates {range_h range_v match mismatch start_gap cont_gap job} {
    global sip_defs DNA PROTEIN old_id_h old_id_v

    #puts AlignUpdates
    set seq_id_h [name_to_seq_id [seq_id_name $range_h]]
    set seq_id_v [name_to_seq_id [seq_id_name $range_v]]

    if {![info exists old_id_h]} {
	set old_id_h $seq_id_h
    }
    if {![info exists old_id_v]} {
	set old_id_v $seq_id_v
    }

    #update range of new sequence
    if {$old_id_h != $seq_id_h} {
	seq_range_updates $range_h
    }
    if {$old_id_v != $seq_id_v} {
	seq_range_updates $range_v
    }

    set tmp 0
    GetSeqTypes seq_id_h seq_id_v type_h type_v tmp tmp tmp tmp use_av_comp

    keylset sg START_GAP [keylget sip_defs SIP.ALIGN.START_GAP]
    keylset cg CONT_GAP [keylget sip_defs SIP.ALIGN.CONT_GAP]

    if {$type_h == $PROTEIN || $type_v == $PROTEIN} {
	entrybox_configure $match -state disabled
	entrybox_configure $mismatch -state disabled
	entrybox_configure $start_gap -default [keylget sg START_GAP.PROTEIN.VALUE]
	entrybox_configure $cont_gap -default [keylget cg CONT_GAP.PROTEIN.VALUE]

    } else {
	entrybox_configure $match -state normal
	entrybox_configure $mismatch -state normal
	entrybox_configure $start_gap -default [keylget sg START_GAP.DNA.VALUE]
	entrybox_configure $cont_gap -default [keylget cg CONT_GAP.DNA.VALUE]
    }
    set old_id_h $seq_id_h
    set old_id_v $seq_id_v
}

proc AlignSeqs { } {
    global sip_defs DNA PROTEIN HORIZONTAL VERTICAL

    set seq_id_h [get_active_seq_id $HORIZONTAL] 
    set seq_id_v [get_active_seq_id $VERTICAL] 

    if {$seq_id_h == -1 && $seq_id_v == -1} {
	verror ERR_WARN "Align sequences" "Horizontal and vertical sequence has not been set in the sequence manager"
	return
    }

    if {$seq_id_h == -1} {
	set seq_id_h $seq_id_v
    }

    if {$seq_id_v == -1} {
	set seq_id_v $seq_id_h
    }

    set s .align_sequences
    if {[xtoplevel $s -resizable 0] == ""} return
    wm title $s "align sequences"
    
    global $seq_id_h.start $seq_id_h.end
    global $seq_id_v.start $seq_id_v.end

    #########################################################################
    #plot ranges
    set seq_length_h [seq_info $seq_id_h length] 
    set seq_length_v [seq_info $seq_id_v length] 
    set seq_start_h [seq_info $seq_id_h start] 
    set seq_start_v [seq_info $seq_id_v start] 
    set seq_end_h [seq_info $seq_id_h end] 
    set seq_end_v [seq_info $seq_id_v end] 

    if {[info exists $seq_id_h.start]} {
	set seq_start_h [set $seq_id_h.start]
    }    
    if {[info exists $seq_id_h.end]} {
	set seq_end_h [set $seq_id_h.end]
    }
    if {[info exists $seq_id_v.start]} {
	set seq_start_v [set $seq_id_v.start]
    }    
    if {[info exists $seq_id_v.end]} {
	set seq_end_v [set $seq_id_v.end]
    }

    keylset us RANGE [keylget sip_defs SIP.CS.RANGE_H]
    seq_id $s.range_h -range 1 -browse 1 -from 1 -to $seq_length_h \
	-start_value $seq_start_h -end_value $seq_end_h -min_value 1 \
	-default [seq_info $seq_id_h name]\
	-title [keylget us RANGE.NAME] \
	-update_cmd [list [list AlignUpdates $s.range_h $s.range_v $s.match $s.mismatch $s.start_gap $s.cont_gap 0]]\
	-browse_cmd seq_browser

    keylset us RANGE [keylget sip_defs SIP.CS.RANGE_V]
    seq_id $s.range_v -range 1 -browse 1 -from 1 -to $seq_length_v \
	-start_value $seq_start_v -end_value $seq_end_v -min_value 1 \
	-default [seq_info $seq_id_v name]\
	-title [keylget us RANGE.NAME] \
	-update_cmd [list [list AlignUpdates $s.range_h $s.range_v $s.match $s.mismatch $s.start_gap $s.cont_gap 1]]\
	-browse_cmd seq_browser

    pack $s.range_h -anchor w -fill x
    pack $s.range_v -anchor w -fill x

    #########################################################################
    #strand selection
    strand_both $s.strand    

    pack $s.strand -fill x
    #read in score matrix match and mismatch values for DNA only
    keylset ms MATCH [keylget sip_defs SIP.ALIGN.MATCH]
    entrybox $s.match \
	    -title [keylget ms MATCH.NAME] \
	    -default [keylget ms MATCH.VALUE] \
	    -width 5 \
	    -type "CheckIntRange [keylget ms MATCH.MIN] \
	                         [keylget ms MATCH.MAX]"

    keylset mm MISMATCH [keylget sip_defs SIP.ALIGN.MISMATCH]
    entrybox $s.mismatch \
	    -title [keylget mm MISMATCH.NAME] \
	    -default [keylget mm MISMATCH.VALUE] \
	    -width 5 \
	    -type "CheckIntRange [keylget mm MISMATCH.MIN] \
	                         [keylget mm MISMATCH.MAX]"
    
    set tmp 0
    GetSeqTypes seq_id_h seq_id_v type_h type_v tmp tmp tmp tmp use_av_comp
    if {$type_h == $PROTEIN || $type_v == $PROTEIN} {
	entrybox_configure $s.match -state disabled
	entrybox_configure $s.mismatch -state disabled
    }
    pack $s.match $s.mismatch -fill x

    keylset sg START_GAP [keylget sip_defs SIP.ALIGN.START_GAP]
    if {$type_h == $DNA} {
	set value [keylget sg START_GAP.DNA.VALUE]
    } else {
	set value [keylget sg START_GAP.PROTEIN.VALUE]
    }

    entrybox $s.start_gap \
	    -title [keylget sg START_GAP.NAME] \
	    -default $value \
	    -width 5 \
	    -type "CheckIntRange [keylget sg START_GAP.MIN] \
	                         [keylget sg START_GAP.MAX]"
	    
    keylset cg CONT_GAP [keylget sip_defs SIP.ALIGN.CONT_GAP]
    if {$type_h == $DNA} {
	set value [keylget cg CONT_GAP.DNA.VALUE]
    } else {
	set value [keylget cg CONT_GAP.PROTEIN.VALUE]
    }
    entrybox $s.cont_gap \
	    -title [keylget cg CONT_GAP.NAME] \
	    -default $value \
	    -width 5 \
	    -type "CheckIntRange [keylget cg CONT_GAP.MIN] \
	                         [keylget cg CONT_GAP.MAX]"
	    
    pack $s.start_gap $s.cont_gap -fill x 

    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $s.button -bd 2 -relief groove \
	    -ok_command "AlignSeqs2 $s $s.range_h $s.range_v \
	    \[strand_get $s.strand\] $s.match $s.mismatch $s.start_gap $s.cont_gap; destroy $s" \
	    -cancel_command "seq_id_destroy $s.range_h; seq_id_destroy $s.range_v; destroy $s" \
	    -help_command "show_help spin {SPIN-Align Sequences}"
    pack $s.button -fill x
}
