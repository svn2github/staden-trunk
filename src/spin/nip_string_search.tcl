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

proc plot_string {seq_id result_id} {
    global nip_defs tk_utils_defs
    global HORIZONTAL SCALE_BAR

    if {$result_id == -1} {
	return
    }
    set type [keylget nip_defs STRINGSEARCH]
    set r_id [CreateRasterGraph raster [list [list $seq_id $HORIZONTAL]] $type 0\
		  [keylget nip_defs RASTER.TITLE]\
		  [keylget nip_defs NIP.STRING_SEARCH.PLOT_HEIGHT] \
		  [keylget nip_defs RASTER.PLOT_WIDTH] \
		  [keylget nip_defs RULER.PLOT_HEIGHT] \
		  [keylget nip_defs RULER.PLOT_WIDTH]]

    nip_string_search plot -window $raster -window_id $r_id\
	     -result_id $result_id -seq_id $seq_id \
	     -fill [keylget nip_defs NIP.STRING_SEARCH.COLOUR] \
	     -width 2 -tick_ht [keylget nip_defs NIP.STRING_SEARCH.TICK_HT]

    set r_win [winfo parent $raster]
    keybox_add $r_win.key$r_id \
	-text "[seq_result_key_name -index [lindex $result_id 0]]" \
	-background  [keylget nip_defs NIP.STRING_SEARCH.COLOUR]\
	-enter "EnterKey $raster [lindex $result_id 0]" -motion MotionRaster \
	-leave "LeaveKey $raster" \
	-drop "DropResult [lindex $result_id 0] $SCALE_BAR"\
	-menu "seq_result_keybox_update $r_win $result_id \[seq_result_names -result_id $result_id\]"

    #update result list
    seq_result_list_update [keylget tk_utils_defs RASTER.RESULTS.WIN]
}

proc NipStringSearch2 {t range} {
    global nip_defs tk_utils_defs
    global HORIZONTAL SCALE_BAR PROTEIN

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

    set result_id [nip_string_search create -seq_id $seq_id \
	    -start [seq_id_from $range] \
	    -end [seq_id_to $range]\
	    -strand $strand -min_pmatch $match -string $string\
	    -use_iub $use_iub]

    # stop windows from hiding the plot
    wm withdraw $t
    
    plot_string $seq_id $result_id
    ClearBusy

    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
}

