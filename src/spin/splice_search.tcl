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
	-ok_command "SpliceSearch2 $w $w.range $w.donor $w.acceptor; destroy $w"\
	-cancel_command "seq_id_destroy $w.range; destroy $w" \
	-help_command "show_help spin {SPIN-Splice-Site-Search}"

    pack $w.range $w.donor $w.acceptor $w.button -fill x
}

proc plot_splice {seq_id result_id} {
    global nip_defs tk_utils_defs
    global HORIZONTAL SCALE_BAR

    #if no results were found
    if {$result_id == "-1 -1 -1"} {
	return
    }
    set type [keylget nip_defs SPLICE]

    set r_id [CreateRasterGraph raster [list [list $seq_id $HORIZONTAL]] $type 1\
		  [keylget nip_defs RASTER.TITLE] \
		  [keylget nip_defs NIP.SPLICE.PLOT_HEIGHT] \
		  [keylget nip_defs RASTER.PLOT_WIDTH] \
		  [keylget nip_defs RULER.PLOT_HEIGHT] \
		  [keylget nip_defs RULER.PLOT_WIDTH]]

    set col_list "[keylget nip_defs NIP.SPLICE.COLOUR.F1] [keylget nip_defs NIP.SPLICE.COLOUR.F2] [keylget nip_defs NIP.SPLICE.COLOUR.F3]"

    splice_search plot -window $raster \
		       -tick_ht [keylget nip_defs NIP.SPLICE.TICK_HT] \
		       -fill $col_list -width 2\
		       -window_id $r_id -seq_id $seq_id -result_id $result_id

    set r_win [winfo parent $raster]
    keybox_add $r_win.key$r_id \
	-text "[seq_result_key_name -index [lindex $result_id 0]]" \
	-background [lindex $col_list 0] \
	-enter "EnterKey $raster [lindex $result_id 0]" -motion MotionRaster \
	-leave "LeaveKey $raster" \
	-drop "DropResult [lindex $result_id 0] $SCALE_BAR"\
	-menu "seq_result_keybox_update $r_win [lindex $result_id 0] \[seq_result_names -result_id [lindex $result_id 0]\]"

    keybox_add $r_win.key$r_id \
	-text "[seq_result_key_name -index [lindex $result_id 1]]" \
	-background  [lindex $col_list 1] \
	-enter "EnterKey $raster [lindex $result_id 1]" -motion MotionRaster \
	-leave "LeaveKey $raster" \
	-drop "DropResult [lindex $result_id 1] $SCALE_BAR"\
	-menu "seq_result_keybox_update $r_win [lindex $result_id 1] \[seq_result_names -result_id [lindex $result_id 1]\]"

    keybox_add $r_win.key$r_id \
	-text "[seq_result_key_name -index [lindex $result_id 2]]" \
	-background [lindex $col_list 2] \
	-enter "EnterKey $raster [lindex $result_id 2]" -motion MotionRaster\
	-leave "LeaveKey $raster" \
	-drop "DropResult [lindex $result_id 2] $SCALE_BAR"\
	-menu "seq_result_keybox_update $r_win [lindex $result_id 2] \[seq_result_names -result_id [lindex $result_id 2]\]"

    #update result list
    seq_result_list_update [keylget tk_utils_defs RASTER.RESULTS.WIN]
}

proc SpliceSearch2 {t range donor acceptor} {
    global PROTEIN

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
    set result_id [splice_search create \
	    -seq_id $seq_id -start [seq_id_from $range] \
	    -end [seq_id_to $range]\
	    -donor $donor\
	    -acceptor $acceptor]

    # stop windows from hiding the plot
    wm withdraw $t
    
    plot_splice $seq_id $result_id
    ClearBusy

    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
}
