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

    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $w.button -bd 2 -relief groove \
	-ok_command "trnaSearch2 $w $w.range; destroy $w"\
	-cancel_command "seq_id_destroy $w.range; destroy $w" \
	-help_command "show_help spin {SPIN-TRNA-Search}"

    pack $w.range $w.button -fill x
}

proc plot_trna_search {seq_id result_id} {
    global nip_defs tk_utils_defs
    global HORIZONTAL SCALE_BAR

    set type [keylget nip_defs TRNA]
    set r_id [CreateRasterGraph raster [list [list $seq_id $HORIZONTAL]] $type 0\
		  [keylget nip_defs RASTER.TITLE]\
		  [keylget nip_defs NIP.TRNA.PLOT_HEIGHT] \
		  [keylget nip_defs RASTER.PLOT_WIDTH] \
		  [keylget nip_defs RULER.PLOT_HEIGHT] \
		  [keylget nip_defs RULER.PLOT_WIDTH]]

    nip_trna_search plot -window $raster \
	    -fill [keylget nip_defs NIP.TRNA.COLOUR] -width 2\
	    -tick_ht [keylget nip_defs NIP.TRNA.TICK_HT] \
	    -window_id $r_id -seq_id $seq_id -result_id $result_id]

    
    set r_win [winfo parent $raster]
    keybox_add $r_win.key$r_id \
	-text "[seq_result_key_name -index [lindex $result_id 0]]" \
	-background  [keylget nip_defs NIP.TRNA.COLOUR]\
	-enter "EnterKey $raster [lindex $result_id 0]" -motion MotionRaster \
	-leave "LeaveKey $raster" \
	-drop "DropResult [lindex $result_id 0] $SCALE_BAR" \
	-menu "seq_result_keybox_update $r_win $result_id \[seq_result_names -result_id $result_id\]"
 
    #update result list
    seq_result_list_update [keylget tk_utils_defs RASTER.RESULTS.WIN]
}

proc trnaSearch2 {t range } {
    global nip_defs tk_utils_defs
    global HORIZONTAL SCALE_BAR PROTEIN

    set seq_id [name_to_seq_id [seq_id_name $range]]
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "trna search" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    SetBusy

    set result_id [nip_trna_search create -start [seq_id_from $range] \
		       -end [seq_id_to $range] -seq_id $seq_id]

    if {$result_id == -1} {
	ClearBusy
	seq_id_destroy $range
	return
    }

    # stop windows from hiding the plot
    wm withdraw $t
    
    plot_trna_search $seq_id $result_id
    ClearBusy
    
    #update result list
    seq_result_list_update [keylget tk_utils_defs RASTER.RESULTS.WIN]
    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
}
