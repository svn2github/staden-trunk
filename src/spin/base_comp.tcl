#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc PlotBaseComp { } {
    global nip_defs 

    set t .plot_base_comp
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Plot base composition"
    
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
    keylset us RANGE [keylget nip_defs NIP.PBC.RANGE]
    seq_id $t.range -range 1 -browse 1 -from 1 -to $seq_length \
	-start_value $seq_start -end_value $seq_end -min_value 1 \
	-default [seq_info $seq_id name] \
	-update_cmd [list [list seq_range_updates $t.range]]\
	-browse_cmd seq_browser

    global $t.t $t.c $t.a $t.g
    set $t.t [keylget nip_defs NIP.PBC.T]
    set $t.c [keylget nip_defs NIP.PBC.C]
    set $t.a [keylget nip_defs NIP.PBC.A]
    set $t.g [keylget nip_defs NIP.PBC.G]
    frame $t.nt
    checkbutton $t.nt.cb_t -text T -variable $t.t
    checkbutton $t.nt.cb_c -text C -variable $t.c
    checkbutton $t.nt.cb_a -text A -variable $t.a
    checkbutton $t.nt.cb_g -text G -variable $t.g

    keylset wl WIN_LEN [keylget nip_defs NIP.PBC.WIN_LEN]
    set win_length [keylget wl WIN_LEN.VALUE]

    set seq_length [expr [seq_id_to $t.range] - [seq_id_from $t.range] + 1]

    if {$win_length > $seq_length} {
	set win_length [expr $seq_length / 2]
    }
    set max_length [keylget wl WIN_LEN.MAX]

    entrybox $t.win_len \
	-title "[keylget wl WIN_LEN.NAME] ([keylget wl WIN_LEN.MIN]\
	             to $max_length)" \
	-default $win_length \
	-width 5 \
	-type "CheckIntRange [keylget wl WIN_LEN.MIN] $max_length "
	    		 
    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "PlotBaseComp2 $t $t.range $t.win_len"\
	-cancel_command "seq_id_destroy $t.range; destroy $t" \
	-help_command "show_help spin {SPIN-Plot-Base-Composition}"

    pack $t.range
    pack $t.nt.cb_t $t.nt.cb_c $t.nt.cb_a $t.nt.cb_g -side left -fill x -expand 1
    pack $t.nt -fill x 
    pack $t.win_len -fill x
    pack $t.button -side bottom -fill x
}

proc plot_base_comp {seq_id result_id} {
    global nip_defs tk_utils_defs
    global HORIZONTAL SCALE_BAR

    set type [keylget nip_defs BASECOMP]
    set r_id [CreateRasterGraph raster [list [list $seq_id $HORIZONTAL]] $type 0\
		  [keylget nip_defs RASTER.TITLE]\
		  [keylget nip_defs RASTER.PLOT_HEIGHT] \
		  [keylget nip_defs RASTER.PLOT_WIDTH] \
		  [keylget nip_defs RULER.PLOT_HEIGHT] \
		  [keylget nip_defs RULER.PLOT_WIDTH]]

    nip_base_comp plot -window $raster \
	    -window_id $r_id \
	    -seq_id $seq_id \
	    -result_id $result_id \
	    -fill [keylget nip_defs NIP.PBC.COLOUR] \
	    -width [keylget nip_defs NIP.PBC.L_WIDTH]
    
    set r_win [winfo parent $raster]
    #add key for all results 
    keybox_add $r_win.key$r_id \
	-text "[seq_result_key_name -index $result_id]" \
	-background  [keylget nip_defs NIP.PBC.COLOUR]\
	-enter "EnterKey $raster $result_id" -motion MotionRaster \
	-leave "LeaveKey $raster" \
	-drop "DropResult $result_id $SCALE_BAR" \
	-menu "seq_result_keybox_update $r_win $result_id \[seq_result_names -result_id $result_id\]"
    fit_on_screen $r_win

    #update result list
    seq_result_list_update [keylget tk_utils_defs RASTER.RESULTS.WIN]
}

proc PlotBaseComp2 { t range win_len} {
    global nip_defs tk_utils_defs
    global $t.t $t.c $t.a $t.g
    global HORIZONTAL SCALE_BAR PROTEIN

    if {![set $t.a] && ![set $t.c] && ![set $t.g] && ![set $t.t]} {
	bell
	verror ERR_WARN "plot base composition" "no base types selected"
	return
    } 
    set window_length [entrybox_get $win_len]
    if {[expr $window_length % 2] == 0} {
	bell
	verror ERR_WARN "plot base composition" "window length must be odd"
	return
    }	

    if {$window_length > [expr [seq_id_to $range] - [seq_id_from $range] + 1]} {
	bell 
	verror ERR_WARN "plot base composition" "window length must be less than sequence length"
	return
    }

    set seq_id [name_to_seq_id [seq_id_name $range]]
    
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "base composition" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    SetBusy

    set result_id [nip_base_comp create -seq_id $seq_id \
		       -start [seq_id_from $range] \
		       -end [seq_id_to $range] \
		       -win_len $window_length \
		       -a [set $t.a] -c [set $t.c] \
		       -g [set $t.g] -t [set $t.t]]

    # stop windows from hiding the plot
    wm withdraw $t

    plot_base_comp $seq_id $result_id
    ClearBusy

    #failed to do base composition
    if {$result_id == "-1"} {
	seq_result_update -index $result_id -job QUIT
	return
    }

    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]

    destroy $t
}

proc CountBaseComp { } {
    global nip_defs

    set t .count_base_comp
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Count sequence composition"
    
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

    keylset us RANGE [keylget nip_defs NIP.PBC.RANGE]
    seq_id $t.range -range 1 -browse 1 -from 1 -to $seq_length \
	-start_value $seq_start -end_value $seq_end -min_value 1 \
	-default [seq_info $seq_id name] \
	-update_cmd [list [list seq_range_updates $t.range]]\
	-browse_cmd nip_seq_browser
    
    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "CountBaseComp2 $t.range; destroy $t"\
	-cancel_command "seq_id_destroy $t.range; destroy $t" \
	-help_command "show_help spin {SPIN-Base-Composition}"

    pack $t.range
    pack $t.button -side bottom -fill x
}

proc CountBaseComp2 {range} {

    set seq_id [name_to_seq_id [seq_id_name $range]] 
    count_base_comp -seq_id $seq_id \
	-start [seq_id_from $range] \
	-end [seq_id_to $range]
    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
}
