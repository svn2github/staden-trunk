#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#############################################################################
#do all necessary translations to make sequence types consistent
proc CheckSequenceTypes {id_h id_v t_h t_v start_h end_h start_v end_v} {
    upvar $id_h seq_id_h $id_v seq_id_v $t_h type_h $t_v type_v
    global DNA PROTEIN

    #check here the types of sequence (dna or protein). If have dna vs protein
    #then automatically translate the dna into three reading frames
    set type_h [seq_info $seq_id_h type]
    set type_v [seq_info $seq_id_v type]

    #check for _rf123 ie dna seq to be treated as protein
    set name_h [string first _rf123 [seq_info $seq_id_h name]]
    set name_v [string first _rf123 [seq_info $seq_id_v name]]

    if {$name_h > -1} {
 	set seq_id_h [seq_translate_seq -seq_id $seq_id_h -f1 1 -f2 1 -f3 1 -start $start_h -end $end_h]
	set type_h $PROTEIN
    }
    if {$name_v > -1} {
	set seq_id_v [seq_translate_seq -seq_id $seq_id_v -f1 1 -f2 1 -f3 1 -start $start_v -end $end_v]
	set type_v $PROTEIN
    }

    if {$type_h == $type_v} {
	set type $type_h
    } else {
	set type $PROTEIN
	if {$type_h == $DNA} {
	    set seq_id_h [seq_translate_seq -seq_id $seq_id_h -f1 1 -f2 1 -f3 1 -start $start_h -end $end_h]
	} else {
	    set seq_id_v [seq_translate_seq -seq_id $seq_id_v -f1 1 -f2 1 -f3 1 -start $start_v -end $end_v]
	}
    }
}

proc CreateSeqArray {seq_id_h seq_id_v start_h end_h start_v end_v s_array num_ele } {
    upvar $s_array seq_array $num_ele cnt

    set cnt 0

    if {[llength $seq_id_h] > 1} {
	for {set i 0} {$i < 3} {incr i} {
	    if {[llength $seq_id_v] > 1} {
		for {set j 0} {$j < 3} {incr j} {
		    set seq_array($cnt,h) [lindex $seq_id_h $i]
		    #set seq_array($cnt,start_h) [expr ($start_h+2)/3]
		    #set seq_array($cnt,end_h) [expr ($end_h+2)/3]
		    set seq_array($cnt,start_h) 1
		    set seq_array($cnt,end_h) [expr ($end_h - $start_h + 1 - $i)/3]
		    set seq_array($cnt,v) [lindex $seq_id_v $j]
		    set seq_array($cnt,start_v) 1
		    set seq_array($cnt,end_v) [expr ($end_v - $start_v + 1 - $j)/3]
		    #set seq_array($cnt,start_v) [expr ($start_v+2)/3]
		    #set seq_array($cnt,end_v) [expr ($end_v+2)/3]
		    incr cnt
		}
	    } else {
		set seq_array($cnt,h) [lindex $seq_id_h $i]
		#set seq_array($cnt,start_h) [expr ($start_h+2)/3]
		set seq_array($cnt,start_h) 1
		set seq_array($cnt,end_h) [expr ($end_h - $start_h + 1 - $i)/3]
		set seq_array($cnt,v) $seq_id_v
		set seq_array($cnt,start_v) $start_v
		set seq_array($cnt,end_v) $end_v
		incr cnt
	    }
	}
    } elseif {[llength $seq_id_v] > 1} {
	for {set i 0} {$i < 3} {incr i} {
	    set seq_array($cnt,h) $seq_id_h
	    #set seq_array($cnt,start_h) $start_h
	    set seq_array($cnt,start_h) 1
	    set seq_array($cnt,end_h) $end_h
	    set seq_array($cnt,v)  [lindex $seq_id_v $i]
	    set seq_array($cnt,start_v) [expr ($start_v+2)/3]
	    set seq_array($cnt,end_v) [expr ($end_v - $start_v + 1 - $i)/3]
	    incr cnt
	}
    } else {
	set seq_array(0,h) $seq_id_h
	set seq_array(0,v) $seq_id_v
	set seq_array(0,start_h) $start_h
	set seq_array(0,end_h) $end_h
	set seq_array(0,start_v) $start_v
	set seq_array(0,end_v) $end_v
	set cnt 1
    }
}

proc plot_similar_spans {seq_id_h seq_id_v res_id s_array num_seq_array} {
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

    set line_width [keylget sip_defs SIP.CS.LINE_WIDTH]

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

	#puts "PLOT $seq_array($i,h) $seq_array($i,v)"
	#don't want to plot results of no matches
	if {$result_id($i) == -1} {
	    continue
	}

	sip_similar_spans plot -result_id $result_id($i) -window $raster \
		-raster_id $r_id \
		-seq_id_h $seq_array($i,h) \
		-seq_id_v $seq_array($i,v) \
		-fill [keylget sip_defs SIP.CS.COLOUR] \
		-width $line_width

	keybox_add $r_win.key$r_id \
		-text "[seq_result_key_name -index $result_id($i)]" \
		-background [keylget sip_defs SIP.CS.COLOUR] \
		-enter "EnterKey $raster $result_id($i)" -motion MotionRaster \
		-leave "LeaveKey $raster" \
		-drop "DropResult $result_id($i) $ZOOM_SCALE" \
		-menu "seq_result_keybox_update $r_win $result_id($i) \[seq_result_names -result_id $result_id($i)\]"
    }
    seq_result_list_update [keylget tk_utils_defs RASTER.RESULTS.WIN]

}

#
#remove direction option, only compare in forward direction only.
proc CompSpans2 {t range_h range_v win_len_win min_score_win} {
    global sip_defs tk_utils_defs
    global DNA PROTEIN HORIZONTAL VERTICAL
    global $min_score_win.score

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

    # win_len needs to be odd. If it is even the algorithm will not work,
    # and infact will corrupt memory!
    set win_len [expr {[entrybox_get $win_len_win]|1}]
    set min_score [entrybox_get $min_score_win]

    set min_score_val [SipFindScore $range_h $range_v $win_len_win]

    #move this check from CheckMinScore since this caused awful nested loops
    if {$min_score < $min_score_val} {
	set res [tk_messageBox -icon warning -type yesnocancel \
		-title "too many matches" \
		-message "This score could produce more than MAXMATCHES of matches and therefore could take a considerable amount of time. Do you really want to do this?"]
	switch $res {
	    yes {}
	    no return
	    cancel return
	}
    }

   
    CheckSequenceTypes seq_id_h seq_id_v type_h type_v $start_h $end_h $start_v $end_v

    CreateSeqArray $seq_id_h $seq_id_v $start_h $end_h $start_v $end_v seq_array num_seq_array

    SetBusy

    set final_id -1
    for {set i 0} {$i < $num_seq_array} {incr i} { 
	set result_id($i) [sip_similar_spans create \
		-seq_id_h $seq_array($i,h) \
		-seq_id_v $seq_array($i,v) \
		-win_len $win_len \
		-min_match $min_score \
		-start_h $seq_array($i,start_h) \
		-end_h $seq_array($i,end_h) \
		-start_v $seq_array($i,start_v) \
		-end_v $seq_array($i,end_v)]
	if {$result_id($i) != -1} {
	    set final_id 0
	}
    }

    #only delete the raster if no other results are plotted in it
    #if {[raster_results -id $r_id -option number] == 0} {
	#seq_result_update -index $r_id -job QUIT
	#ClearBusy
	#return
    #}

    #check if any results have been found
    if {$final_id != -1} {

	# stop windows from hiding the plot
	wm withdraw $t

	plot_similar_spans $seq_id_h $seq_id_v result_id seq_array $num_seq_array

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
    unset $min_score_win.score
    destroy $t
}

proc GetType {seq_id_h seq_id_v} {

    set type_h [seq_info $seq_id_h type]
    set type_v [seq_info $seq_id_v type]

    if {$type_h == $type_v} {
	return $type_h
    } else {
	tk_messageBox -icon error -type ok \
		-title "sequences of different types" \
		-message "both sequences must either be DNA or protein"
	return ""
    }

}

#takes into account dna vs protein, _rf123 etc
proc GetSeqTypes {id_h id_v t_h t_v s_h e_h s_v e_v comp} {
    upvar $id_h seq_id_h $id_v seq_id_v $t_h type_h $t_v type_v $s_h start_h $e_h end_h $s_v start_v $e_v end_v $comp use_av_comp
    global PROTEIN DNA

    set type_h [seq_info $seq_id_h type]
    set type_v [seq_info $seq_id_v type]

    set use_av_comp 0

    #check for _rf123 ie dna seq to be treated as protein
    set name_h [string first _rf123 [seq_info $seq_id_h name]]
    set name_v [string first _rf123 [seq_info $seq_id_v name]]

    set old_type_h $type_h
    set old_type_v $type_v

    if {$type_h != $type_v} {
	if {$type_h == $DNA} {
	    #use only protein seq to estimate score
	    set seq_id_h $seq_id_v
	    set type_h $PROTEIN
	    set start_h $start_v
	    set end_h $end_v
	} else {
	    set seq_id_v $seq_id_h
	    set type_v $PROTEIN
	    set start_v $start_h
	    set end_v $end_h
	}
    }

    if {$name_h > -1} {
	set type_h $PROTEIN
	#only need to use use_av_comp if neither type is protein
	if {$old_type_v == $DNA} {
	    set use_av_comp 1
	}
    }
    if {$name_v > -1} {
	set type_v $PROTEIN
	if {$old_type_h == $DNA} {
	    set use_av_comp 1
	}
    }
}

#updates to be done when the sequence name changes
proc SeqNameUpdates {range_h range_v win_len min_score} {
    global old_id_h old_id_v

    #puts "SeqNameUpdates"
    set seq_id_h [name_to_seq_id [seq_id_name $range_h]]
    set seq_id_v [name_to_seq_id [seq_id_name $range_v]]

    if {![info exists old_id_h]} {
	set old_id_h $seq_id_h
    }
    if {![info exists old_id_v]} {
	set old_id_v $seq_id_v
    }

    #puts "old_id_h $old_id_h $seq_id_h"

    #update range of new sequence
    if {$old_id_h != $seq_id_h} {
	seq_range_updates $range_h
    }
    if {$old_id_v != $seq_id_v} {
	seq_range_updates $range_v
    }
    set start_h [seq_id_from $range_h]
    set end_h [seq_id_to $range_h]
    set start_v [seq_id_from $range_v]
    set end_v [seq_id_to $range_v]

    UpdateMinScore $seq_id_h $seq_id_v $start_h $end_h $start_v $end_v $win_len $min_score
    set old_id_h $seq_id_h
    set old_id_v $seq_id_v
}

#updates to be done when the win length changes
proc WinLenUpdates {range_h range_v win_len min_score args} {

    #puts WinLenUpdates

    set seq_id_h [name_to_seq_id [seq_id_name $range_h $win_len.entry]]
    set seq_id_v [name_to_seq_id [seq_id_name $range_v $win_len.entry]]
    set start_h [seq_id_from $range_h]
    set end_h [seq_id_to $range_h]
    set start_v [seq_id_from $range_v]
    set end_v [seq_id_to $range_v]

    UpdateMinScore $seq_id_h $seq_id_v $start_h $end_h $start_v $end_v $win_len $min_score

}

##############################################################################
proc CompareSpans { } {
    global sip_defs PROTEIN DNA HORIZONTAL VERTICAL

    set seq_id_h [get_active_seq_id $HORIZONTAL] 
    set seq_id_v [get_active_seq_id $VERTICAL] 

    #check to see if horizontal and vertical sequences have been set
    #if {$seq_id_h == -1 || $seq_id_v == -1} {
	#verror ERR_WARN "Find similar spans" "Horizontal or vertical sequence has not been set in the sequence manager"
	#return
    #}

    if {$seq_id_h == -1 && $seq_id_v == -1} {
	verror ERR_WARN "Find similar spans" "Horizontal and vertical sequence has not been set in the sequence manager"
	return
    }

    if {$seq_id_h == -1} {
	set seq_id_h $seq_id_v
    }

    if {$seq_id_v == -1} {
	set seq_id_v $seq_id_h
    }

    set s .find_similar_spans
    if {[xtoplevel $s -resizable 0] == ""} return
    wm title $s "find similar spans"

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
	-update_cmd [list [list SeqNameUpdates $s.range_h $s.range_v $s.win_len $s.min_score]] \
	-browse_cmd seq_browser

    keylset us RANGE [keylget sip_defs SIP.CS.RANGE_V]
    seq_id $s.range_v -range 1 -browse 1 -from 1 -to $seq_length_v \
	-start_value $seq_start_v -end_value $seq_end_v -min_value 1 \
	-default [seq_info $seq_id_v name]\
	-title [keylget us RANGE.NAME] \
	-update_cmd [list [list SeqNameUpdates $s.range_h $s.range_v $s.win_len $s.min_score]]\
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
    #ok cancel help buttons 
    okcancelhelp $s.button -bd 2 -relief groove \
	-ok_command "CompSpans2 $s $s.range_h $s.range_v $s.win_len \
	$s.min_score"\
	-cancel_command "seq_id_destroy $s.range_h; seq_id_destroy $s.range_v; unset $s.min_score.score; destroy $s" \
	-help_command "show_help spin {SPIN-Find similar spans}"

    pack $s.button -fill x
}

proc update_score {min_score} {
    global $min_score.score

    set $min_score.score [entrybox_get $min_score]
}

##############################################################################
proc winlen_minscore {range_h range_v win_len min_score seq_id_h seq_id_v} {
    global sip_defs 
    global $min_score.score
    global $win_len.win_len

    #########################################################################
    #window length
    keylset wl WIN_LEN [keylget sip_defs SIP.CS.WIN_LEN]
    set $win_len.win_len ""

    entrybox $win_len \
	    -title "[keylget wl WIN_LEN.NAME] ([keylget wl WIN_LEN.MIN]\
	             to [keylget wl WIN_LEN.MAX])" \
	    -default [keylget wl WIN_LEN.VALUE]\
	    -width 5 \
	    -type "CheckIntRange [keylget wl WIN_LEN.MIN] \
	                         [keylget wl WIN_LEN.MAX] " \
	    -command "WinLenUpdates $range_h $range_v $win_len $min_score"

    #update min score title when leave win len entry box
    bind $win_len.entry <Any-Leave> "WinLenUpdates $range_h $range_v $win_len $min_score"

    #########################################################################
    #minimum score
    keylset ms MIN_SCORE [keylget sip_defs SIP.CS.MIN_SCORE]

    #entrybox $min_score \
#	-width 5 \
#	-type "CheckIntMin [SipFindScore $range_h $range_v $win_len]"

    entrybox $min_score \
	    -width 5 \
	    -type "CheckMinScore [SipFindScore $range_h $range_v $win_len]"

    global $min_score.score
    #bind $min_score.entry <Any-Leave> "puts here"
    bind $min_score.entry <Any-Leave> "update_score $min_score"

    WinLenUpdates $range_h $range_v $win_len $min_score

    return 0
}

proc SipFindScore {range_h range_v win_len} {
    global DNA

    set max_matches [get_max_matches]

    update idletasks

    set seq_id_h [name_to_seq_id [seq_id_name $range_h]]
    set seq_id_v [name_to_seq_id [seq_id_name $range_v]]
    set start_h [seq_id_from $range_h]
    set end_h [seq_id_to $range_h]
    set start_v [seq_id_from $range_v]
    set end_v [seq_id_to $range_v]


    GetSeqTypes seq_id_h seq_id_v type_h type_v start_h end_h start_v end_v use_av_comp

    sip_find_probs -win_len [entrybox_get $win_len] -seq_id_h $seq_id_h \
	-seq_id_v $seq_id_v -start_h $start_h -end_h $end_h \
	-start_v $start_v -end_v $end_v -use_av_comp $use_av_comp
    set score [sip_find_score -win_len [entrybox_get $win_len] \
		   -num_matches $max_matches\
		   -seq_id_h $seq_id_h -seq_id_v $seq_id_v \
		   -start_h $start_h -end_h $end_h \
		   -start_v $start_v -end_v $end_v -use_av_comp $use_av_comp]
    return $score
}

##############################################################################
proc UpdateMinScore {seq_id_h seq_id_v start_h end_h start_v end_v win_len min_score} {
    global sip_defs DNA PROTEIN $min_score.score
    
    GetSeqTypes seq_id_h seq_id_v type_h type_v start_h end_h start_v end_v use_av_comp

    set max_matches [get_max_matches]

    set type_h [seq_info $seq_id_h type]
    set type_v [seq_info $seq_id_v type]

    sip_find_probs -win_len [entrybox_get $win_len] -seq_id_h $seq_id_h \
	-seq_id_v $seq_id_v -start_h $start_h -end_h $end_h \
	-start_v $start_v -end_v $end_v -use_av_comp $use_av_comp \
	-type_h $type_h -type_v $type_v

    #update title on minscore entrybox
    keylset ms MIN_SCORE [keylget sip_defs SIP.CS.MIN_SCORE]
    set min_score_val [sip_find_score -win_len [entrybox_get $win_len] \
	    -num_matches $max_matches -seq_id_h $seq_id_h -seq_id_v $seq_id_v \
            -start_h $start_h -end_h $end_h \
	    -start_v $start_v -end_v $end_v -use_av_comp $use_av_comp]
    set max_score_val [sip_find_score -win_len [entrybox_get $win_len] \
	         -num_matches 1 -seq_id_h $seq_id_h \
                 -seq_id_v $seq_id_v -start_h $start_h -end_h $end_h \
	         -start_v $start_v -end_v $end_v -use_av_comp $use_av_comp]
    entrybox_configure $min_score \
	    -title "[keylget ms MIN_SCORE.NAME]\
	    ($min_score_val to $max_score_val)" \
	    -type "CheckMinScore $min_score_val"

    #update default entry in minscore entrybox
    set score [sip_find_score -win_len [entrybox_get $win_len] \
		   -num_matches [get_def_matches]\
		   -seq_id_h $seq_id_h -seq_id_v $seq_id_v \
		   -start_h $start_h -end_h $end_h \
		   -start_v $start_v -end_v $end_v -use_av_comp $use_av_comp]
    if {$score < $min_score_val} {
	set score $min_score_val
    }
    if {$score > $max_score_val} {
	set score $max_score_val
    }

    if {![info exists $min_score.score] } {
	set $min_score.score $score
	entrybox_configure $min_score -default $score
    }


    set $min_score.score $score
    entrybox_configure $min_score -default [set $min_score.score]
}

##############################################################################
proc CheckMinScore {min path} {

    set current [$path.entry get]
    if {![isinteger $current]} {
	return 0
    }

    #moved this check into the final press OK button check.
    if {$current < $min} {
	return 1
    }
	
    return 1
}
