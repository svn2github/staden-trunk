#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
##############################################################################
proc plot_best_diagonals {seq_id_h seq_id_v res_id s_array num_seq_array} {
    global sip_defs tk_utils_defs
    global DNA PROTEIN HORIZONTAL VERTICAL ZOOM_SCALE

    upvar $res_id result_id $s_array seq_array 

    set r_id [CreateRasterDot raster [lindex $seq_id_h 0] \
	    [lindex $seq_id_v 0] \
	    [keylget sip_defs RASTER.TITLE]\
	    [keylget sip_defs RASTER.PLOT_HEIGHT] \
	    [keylget sip_defs RASTER.PLOT_WIDTH] \
	    [keylget sip_defs RULER.PLOT_HEIGHT] \
	    [keylget sip_defs RULER.PLOT_WIDTH]]

    set line_width [keylget sip_defs SIP.QS.LINE_WIDTH]

    #must register each of translated sequences with raster.
    #Note the first sequence was used in CreateRasterDot
    if {[llength $seq_id_h] > 1} {
	add_seq_to_raster -raster_id $r_id \
	    -seq_id [lindex $seq_id_h 1] \
	    -direction $HORIZONTAL -line_width $line_width	    
	add_seq_to_raster -raster_id $r_id \
	    -seq_id [lindex $seq_id_h 2] \
	    -direction $HORIZONTAL -line_width $line_width
    }
    if {[llength $seq_id_v] > 1} {
	add_seq_to_raster -raster_id $r_id \
	    -seq_id [lindex $seq_id_v 1] \
	    -direction $VERTICAL -line_width $line_width	    
	add_seq_to_raster -raster_id $r_id \
	    -seq_id [lindex $seq_id_v 2] \
	    -direction $VERTICAL -line_width $line_width
    }

    set r_win [winfo parent $raster]

    for {set i 0} {$i < $num_seq_array} {incr i} { 

	#don't want to plot results of no matches
	if {$result_id($i) == -1} {
	    continue
	}
	sip_best_diagonals plot -result_id $result_id($i) -window $raster \
		-raster_id $r_id \
		-seq_id_h $seq_array($i,h) \
		-seq_id_v $seq_array($i,v) \
		-fill [keylget sip_defs SIP.QS.COLOUR] \
		-width $line_width
	
	keybox_add $r_win.key$r_id \
		-text "[seq_result_key_name -index $result_id($i)]" \
		-background [keylget sip_defs SIP.QS.COLOUR]\
		-enter "EnterKey $raster $result_id($i)" -motion MotionRaster \
		-leave "LeaveKey $raster" \
		-drop "DropResult $result_id($i) $ZOOM_SCALE" \
		-menu "seq_result_keybox_update $r_win $result_id($i) \[seq_result_names -result_id $result_id($i)\]"
    }
    seq_result_list_update [keylget tk_utils_defs RASTER.RESULTS.WIN]
}

proc QuickScan2 {t range_h range_v win_len word_len min_match sd} {
    global sip_defs

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

     #do this before SetBusy in case of invalid entries
    set winlen [entrybox_get $win_len]
    set minmatch [entrybox_get $min_match]
    set wordlen [entrybox_get $word_len]
    set SD [entrybox_get $sd]

    CheckSequenceTypes seq_id_h seq_id_v type_h type_v $start_h $end_h $start_v $end_v

    CreateSeqArray $seq_id_h $seq_id_v $start_h $end_h $start_v $end_v seq_array num_seq_array

    SetBusy
    set final_id -1
    for {set i 0} {$i < $num_seq_array} {incr i} { 
	set result_id($i) [sip_best_diagonals create \
		-seq_id_h $seq_array($i,h) \
		-seq_id_v $seq_array($i,v) \
		-start_h $seq_array($i,start_h) \
		-end_h $seq_array($i,end_h) \
		-start_v $seq_array($i,start_v) \
		-end_v $seq_array($i,end_v)\
		-win_len $winlen \
		-min_match $minmatch \
		-word_len $wordlen\
		-sd $SD]
	if {$result_id($i) != -1} {
	    set final_id 0
	}
    }

    #check if any results have been found
    if {$final_id != -1} {

	# stop windows from hiding the plot
	wm withdraw $t

	plot_best_diagonals $seq_id_h $seq_id_v result_id seq_array $num_seq_array
	
	sequence_list_update

	for {set i 0} {$i < $num_seq_array} {incr i} { 
	    global $seq_array($i,h).start $seq_array($i,h).end
	    global $seq_array($i,v).start $seq_array($i,v).end
	    set $seq_array($i,h).start $seq_array($i,start_h)
	    set $seq_array($i,h).end $seq_array($i,end_h)
	    set $seq_array($i,v).start $seq_array($i,start_v) 
	    set $seq_array($i,v).end $seq_array($i,end_v)
	}
    }

    ClearBusy
    seq_id_destroy $range_h
    seq_id_destroy $range_v
    destroy $t
	
}

##############################################################################
proc QuickScanUpdates {range_h range_v win_len min_score word_len} {
    global sip_defs DNA PROTEIN

    SeqNameUpdates $range_h $range_v $win_len $min_score
    
    set seq_id_h [name_to_seq_id [seq_id_name $range_h]]
    set seq_id_v [name_to_seq_id [seq_id_name $range_v]]

    set tmp 0
    GetSeqTypes seq_id_h seq_id_v type_h type_v tmp tmp tmp tmp use_av_comp

    keylset wl WORD_LEN [keylget sip_defs SIP.QS.WORD_LEN]
    if {$type_h == $PROTEIN || $type_v == $PROTEIN} {
	entrybox_configure $word_len \
	    -title "[keylget wl WORD_LEN.NAME] ([keylget wl WORD_LEN.PROTEIN.MIN] \
	    to [keylget wl WORD_LEN.PROTEIN.MAX])"\
	    -default [keylget wl WORD_LEN.PROTEIN.VALUE] \
	    -type "CheckIntRange [keylget wl WORD_LEN.PROTEIN.MIN]\
	                         [keylget wl WORD_LEN.PROTEIN.MAX]"
    } else {
	entrybox_configure $word_len \
		-title "[keylget wl WORD_LEN.NAME] ([keylget wl WORD_LEN.DNA.MIN] \
		to [keylget wl WORD_LEN.DNA.MAX])"\
		-default [keylget wl WORD_LEN.DNA.VALUE] \
		-type "CheckIntRange [keylget wl WORD_LEN.DNA.MIN]\
		[keylget wl WORD_LEN.DNA.MAX]"
    }

}

##############################################################################
proc QuickScan { } {
    global sip_defs PROTEIN DNA HORIZONTAL VERTICAL

    set seq_id_h [get_active_seq_id $HORIZONTAL] 
    set seq_id_v [get_active_seq_id $VERTICAL] 

    #check to see if horizontal and vertical sequences have been set
    #if {$seq_id_h == -1 || $seq_id_v == -1} {
	#verror ERR_WARN "Find best diagonals" "Horizontal or vertical sequence has not been set in the sequence manager"
	#return
    #}

    if {$seq_id_h == -1 && $seq_id_v == -1} {
	verror ERR_WARN "Find best diagonals" "Horizontal and vertical sequence has not been set in the sequence manager"
	return
    }

    if {$seq_id_h == -1} {
	set seq_id_h $seq_id_v
    }

    if {$seq_id_v == -1} {
	set seq_id_v $seq_id_h
    }
    

    set s .find_best_diagonals
    if {[xtoplevel $s -resizable 0] == ""} return
    wm title $s "find best diagonals"

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
	-update_cmd [list [list QuickScanUpdates $s.range_h $s.range_v $s.win_len $s.min_score $s.word_len]]\
	-browse_cmd seq_browser

    keylset us RANGE [keylget sip_defs SIP.CS.RANGE_V]
    seq_id $s.range_v -range 1 -browse 1 -from 1 -to $seq_length_v \
	-start_value $seq_start_v -end_value $seq_end_v -min_value 1 \
	-default [seq_info $seq_id_v name]\
	-title [keylget us RANGE.NAME] \
	-update_cmd [list [list QuickScanUpdates $s.range_h $s.range_v $s.win_len $s.min_score $s.word_len]]\
	-browse_cmd seq_browser

    #########################################################################
    #window length and min score widgets
    if {[winlen_minscore $s.range_h $s.range_v $s.win_len $s.min_score $seq_id_h $seq_id_v] == -1} {
	destroy $s
	return 
    }

    pack $s.range_h -anchor w -fill x
    pack $s.range_v -anchor w -fill x
    pack $s.win_len -anchor w -fill x
    pack $s.min_score -anchor w -fill x

    #########################################################################
    #word length
    keylset wl WORD_LEN [keylget sip_defs SIP.QS.WORD_LEN]

    set tmp 0
    GetSeqTypes seq_id_h seq_id_v type_h type_v tmp tmp tmp tmp use_av_comp
    
    if {$type_h == $PROTEIN || $type_v == $PROTEIN} {
	entrybox $s.word_len \
	    -title "[keylget wl WORD_LEN.NAME] ([keylget wl WORD_LEN.PROTEIN.MIN] \
	    to [keylget wl WORD_LEN.PROTEIN.MAX])"\
	    -default [keylget wl WORD_LEN.PROTEIN.VALUE] \
	    -width 5 \
	    -type "CheckIntRange [keylget wl WORD_LEN.PROTEIN.MIN]\
	                         [keylget wl WORD_LEN.PROTEIN.MAX]"
    } else {
	entrybox $s.word_len \
		-title "[keylget wl WORD_LEN.NAME] ([keylget wl WORD_LEN.DNA.MIN] \
		to [keylget wl WORD_LEN.DNA.MAX])"\
		-default [keylget wl WORD_LEN.DNA.VALUE] \
		-width 5 \
		-type "CheckIntRange [keylget wl WORD_LEN.DNA.MIN]\
		[keylget wl WORD_LEN.DNA.MAX]"
    }
    pack $s.word_len -anchor w -fill x
    
    #########################################################################
    #number of standard deviations
    keylset sd SD [keylget sip_defs SIP.QS.SD]
    entrybox $s.sd \
	    -title "[keylget sd SD.NAME] ([keylget sd SD.MIN] to \
	                                  [keylget sd SD.MAX])"\
	    -default [keylget sd SD.VALUE] \
	    -width 5 \
	    -type "CheckFloatRange [keylget sd SD.MIN] [keylget sd SD.MAX]"
    pack $s.sd -anchor w -fill x
    

    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $s.button -bd 2 -relief groove \
	-ok_command "QuickScan2 $s $s.range_h $s.range_v $s.win_len $s.word_len \
	    $s.min_score $s.sd" \
	-cancel_command "seq_id_destroy $s.range_h; seq_id_destroy $s.range_v; unset $s.min_score.score; destroy $s" \
	-help_command "show_help spin {SPIN-Find Best Diagonals}"
    
    pack $s.button -fill x
}
