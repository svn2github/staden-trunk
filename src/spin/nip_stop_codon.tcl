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

proc plot_codons {seq_id result_id col_list} {
    global nip_defs tk_utils_defs
    global HORIZONTAL SCALE_BAR

    #if no results were found
    if {$result_id == "-1 -1 -1"} {
	return
    }
    set type [keylget nip_defs STOPCODON]

    set cnt 0
    set r_id_list ""
    set raster_list ""
    set result_id_list ""
    for {set i 0} {$i < 3} {incr i} {
	if {[lindex $result_id $i] != -1} {
	    set frame [expr $i + 1]
	    set r_id [CreateRasterGraph raster [list [list $seq_id $HORIZONTAL]] $type $frame\
		    [keylget nip_defs RASTER.TITLE]\
		    [keylget nip_defs RASTER.SINGLE.PLOT_HEIGHT] \
		    [keylget nip_defs RASTER.PLOT_WIDTH] \
		    [keylget nip_defs RULER.PLOT_HEIGHT] \
		    [keylget nip_defs RULER.PLOT_WIDTH]]
	    lappend r_id_list $r_id
	    lappend raster_list $raster
	    lappend result_id_list [lindex $result_id $i]

	}
    }

    nip_stop_codons plot -window $raster_list\
		       -tick_ht [keylget nip_defs NIP.STOP_CODON.TICK_HT] \
		       -fill $col_list\
		       -width 0\
		       -window_id $r_id_list\
		       -seq_id $seq_id -result_id $result_id_list

    set cnt [llength $result_id_list]
    for {set i 0} {$i < $cnt} {incr i} {
	set r_win [winfo parent [lindex $raster_list $i]]
	keybox_add $r_win.key[lindex $r_id_list $i] \
		-text "[seq_result_key_name -index [lindex $result_id_list $i]]" \
		-background [lindex $col_list $i] \
		-enter "EnterKey [lindex $raster_list $i] [lindex $result_id_list $i]" \
		-motion MotionRaster \
		-leave "LeaveKey [lindex $raster_list $i]" \
		-drop "DropResult [lindex $result_id_list $i] $SCALE_BAR"\
		-menu "seq_result_keybox_update $r_win [lindex $result_id_list $i] \[seq_result_names -result_id [lindex $result_id_list $i]\]"
    }
}

proc NipStopCodon2 {t range strand } {
    global PROTEIN nip_defs

    set seq_id [name_to_seq_id [seq_id_name $range]]
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "stop codons" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    SetBusy

    set result_id [nip_stop_codons create\
	    -start [seq_id_from $range] \
	    -end [seq_id_to $range] \
	    -strand $strand -seq_id $seq_id]
    
    if {$t != ""} {
	wm withdraw $t
    }


    set col_list "[keylget nip_defs NIP.STOP_CODON.COLOUR.F1] [keylget nip_defs NIP.STOP_CODON.COLOUR.F2] [keylget nip_defs NIP.STOP_CODON.COLOUR.F3]"

    plot_codons $seq_id $result_id $col_list
 
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
    #strand_both $t.strand -bd 2 -relief groove

    #########################################################################
    #OK and Cancel buttons
    okcancelhelp $t.ok_cancel \
	    -ok_command "NipStartCodon2 $t $t.range; destroy $t" \
	    -cancel_command "seq_id_destroy $t.range; destroy $t" \
	    -help_command "show_help spin {SPIN-Start-Codon-Search}" \
	    -bd 2 \
	    -relief groove

    #final packing
    pack $t.range -fill both
    #pack $t.strand -fill both
    pack $t.ok_cancel -fill x
}

proc NipStartCodon2 {t range} {
    global PROTEIN nip_defs

    set seq_id [name_to_seq_id [seq_id_name $range]]
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "start codons" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    SetBusy

    set result_id [nip_start_codons create\
	    -start [seq_id_from $range] \
	    -end [seq_id_to $range] \
	    -seq_id $seq_id]
    
    wm withdraw $t

    set col_list "[keylget nip_defs NIP.START_CODON.COLOUR.F1] [keylget nip_defs NIP.START_CODON.COLOUR.F2] [keylget nip_defs NIP.START_CODON.COLOUR.F3]"

    plot_codons $seq_id $result_id $col_list
 
    ClearBusy

    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
}


