#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc FindWinLen { default } {

    set seq_length [expr [s_length] /3]

    if {$seq_length < $default} {
	
	set middle [expr $seq_length / 2]
	while {[WinLenCheck $seq_length $middle 0] != 1} {
	    incr middle -1
	}
	return $middle
    }
    return $default
}

proc WinLenCheck {seq_length win_length show_error } {
    global nip_defs

    keylset wl WIN_LEN [keylget nip_defs NIP.PGS.WIN_LEN]

    #must be odd
    if {$win_length % 2 == 0} {
	if {$show_error} {
	    verror ERR_WARN "gene search" "window length must be odd"
	}
	return 0
    } 

    #must be less than the sequence length
    if {[expr $win_length * 3] > $seq_length} {
	if {$show_error} {
	    verror ERR_WARN "gene search" "3 * window length must be less than the sequence length"
	}
	return 0
    }

    #must be greater than min win length
    if {$win_length < [keylget wl WIN_LEN.MIN]} {
	return 0
    }

    #must be less than max win length
    if {$win_length > [keylget wl WIN_LEN.MAX]} {
	return 0
    }
    return 1
}

proc CheckWinLen {seq_length path} {

    set win_length [[entrybox_path $path] get]
    if {[WinLenCheck $seq_length $win_length 1] == 0} {
	return 0
    }
    return 1
}

#check if valid values have been used
proc GeneSearchInputCheck {range win_len} {
    
    set window_length [expr [entrybox_get $win_len] * 3]
    if {[expr $window_length % 2] == 0} {
	bell
	verror ERR_WARN "plot gene search" "window length must be odd"
	return -1
    }	
    if {[expr $window_length % 3] != 0} {
	bell
	verror ERR_WARN "plot gene search" "window length must be a muliple of 3"
	return -1
    }	
    if {$window_length > [expr [seq_id_to $range] - [seq_id_from $range] - 1]} {
	bell 
	verror ERR_WARN "plot gene search" "window length must be less than sequence length"
	return -1
    }
    return 0
}

proc BaseBias { } {
    global nip_defs

    set w .base_bias
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Base bias"

    keylset wl WIN_LEN [keylget nip_defs NIP.PGS.WIN_LEN]

    set win_length [FindWinLen [keylget wl WIN_LEN.VALUE]]
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

    entrybox $w.win_len \
	-title "[keylget wl WIN_LEN.NAME] ([keylget wl WIN_LEN.MIN]\
	             to [keylget wl WIN_LEN.MAX])" \
	-default $win_length \
	-width 5 \
	-type "CheckWinLen $seq_length"
	    		 
    keylset us RANGE [keylget nip_defs NIP.PGS.RANGE]
    seq_id $w.range -range 1 -browse 1 -from 1 -to $seq_length \
	-start_value $seq_start -end_value $seq_end -min_value 1 \
	-default [seq_info $seq_id name]\
	-update_cmd [list [list seq_range_updates $w.range]]\
	-browse_cmd nip_seq_browser

    #########################################################################
    #strand selection
    strand_both $w.strand    

    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $w.button -bd 2 -relief groove \
	-ok_command "BaseBias2 $w $w.win_len $w.range \[strand_get $w.strand\]"\
	-cancel_command "seq_id_destroy $w.range; destroy $w" \
	-help_command "show_help spin {SPIN-Uneven-Positional-Base-Freqs}"

    pack $w.range $w.win_len -anchor w -fill x
    pack $w.strand -fill x 
    pack $w.button -side bottom -fill x
}


proc BaseBias2 { w win_len range strand} {
    global nip_defs spin_defs tk_utils_defs PROTEIN TOP_S BOTTOM_S
    global HORIZONTAL SCALE_BAR

    if {-1 == [GeneSearchInputCheck $range $win_len]} {
	return 
    }
    set win_length [entrybox_get $win_len]
    set seq_id [name_to_seq_id [seq_id_name $range]]

    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "base bias" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    SetBusy
    if {[expr $strand & $TOP_S]} {
	set res [nip_base_bias create -win_len $win_length \
		-start [seq_id_from $range] \
		-end [seq_id_to $range]\
		-seq_id $seq_id \
		-strand $TOP_S]

	set result_id [lindex $res 0]
	set results [lindex $res 1]

	# stop windows from hiding the plot
	wm withdraw $w
	if {-1 == [plot_base_bias $seq_id $result_id $results $TOP_S]} {
	    seq_id_destroy $range
	    ClearBusy
	    return
	}
    }

    if {[expr $strand & $BOTTOM_S]} {
	set res [nip_base_bias create -win_len $win_length \
		-start [seq_id_from $range] \
		-end [seq_id_to $range]\
		-seq_id $seq_id \
		-strand $BOTTOM_S]

	set result_id [lindex $res 0]
	set results [lindex $res 1]

	# stop windows from hiding the plot
	wm withdraw $w
	if {-1 == [plot_base_bias $seq_id $result_id $results $BOTTOM_S]} {
	    #seq_result_update -index $r_id1 -job QUIT
	    seq_id_destroy $range
	    ClearBusy
	    return
	}
    }
    ClearBusy

    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
    destroy $w
}

proc plot_base_bias {seq_id result_id results strand} {
    global nip_defs spin_defs tk_utils_defs
    global HORIZONTAL SCALE_BAR
    
    if {$result_id == -1} {
	return -1
    }
    set type [keylget nip_defs BASEBIAS]

    set frame 0

    set element_info [create_seq_element $seq_id -1 $type $strand $frame \
	    $HORIZONTAL XY "CANVAS" [keylget spin_defs CONTAINER.TITLE] \
	    [keylget tk_utils_defs ELEMENT.PLOT_WIDTH] \
	    [keylget tk_utils_defs ELEMENT.PLOT_HEIGHT]]
    set c_win [keylget element_info container_win]
    set c_id [keylget element_info container_id]
    set e_win [keylget element_info element_win]
    set e_id [keylget element_info element_id]
    set orientation [keylget element_info orientation]

    update idletasks
    
    nip_base_bias plot -element $c_win$e_win\
	    -container $c_win\
	    -seq_id $seq_id \
	    -result_id $result_id\
	    -results $results\
	    -container_id [keylget element_info container_id]\
	    -element_id [keylget element_info element_id]\
	    -element_type "CANVAS"\
	    -fill [keylget nip_defs NIP.BASEBIAS.COLOUR] \
	    -orientation $orientation

    seqed_element_bindings $c_id $c_win$e_win $e_id
    fit_on_screen $c_win

    #update result list
    result_list_update $c_win 
    return 0
}

#use by codon pref and author test algorithms
proc GeneSearch_CTable { w type } {
    global nip_defs

    keylset wl WIN_LEN [keylget nip_defs NIP.PGS.WIN_LEN]

    set win_length [FindWinLen [keylget wl WIN_LEN.VALUE]]
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

    entrybox $w.win_len \
	-title "[keylget wl WIN_LEN.NAME] ([keylget wl WIN_LEN.MIN]\
	             to [keylget wl WIN_LEN.MAX])" \
	-default $win_length \
	-width 5 \
	-type "CheckWinLen $seq_length"
	    		 
    keylset us RANGE [keylget nip_defs NIP.PGS.RANGE]
    seq_id $w.range -range 1 -browse 1 -from 1 -to $seq_length \
	-start_value $seq_start -end_value $seq_end -min_value 1 \
	-default [seq_info $seq_id name]\
	-update_cmd [list [list seq_range_updates $w.range]]\
	-browse_cmd nip_seq_browser

    #get codon table file name

    if {[string compare $type "load"] == 0} {
	keylset ct C_TABLE [keylget nip_defs NIP.PGS.C_TABLE]
    } else {
	keylset ct C_TABLE [keylget nip_defs NIP.PGS.C_TABLE_O]
    }
    eval getFname $w.fname [list [keylget ct C_TABLE.NAME]] $type {} \
	    [keylget ct C_TABLE.VALUE]

    pack $w.range $w.win_len $w.fname -anchor w -fill x
}

proc plot_author_test {seq_id r_id res strand} {
    global nip_defs spin_defs tk_utils_defs HORIZONTAL
    global SCALE_BAR

    upvar $r_id result_id $res results
    
    set type [keylget nip_defs CODONPREF]

    set cnt 0

    set col_list "[keylget nip_defs NIP.AUTHOR.COLOUR.F1] [keylget nip_defs NIP.AUTHOR.COLOUR.F2] [keylget nip_defs NIP.AUTHOR.COLOUR.F3]"

    for {set i 0} {$i < 3} {incr i} {
	if {$result_id($i) != -1} {
	    set frame [expr $i + 1]
	    set element_info [create_seq_element $seq_id -1 $type $strand $frame \
		    $HORIZONTAL XY "CANVAS" [keylget spin_defs CONTAINER.TITLE] \
		    [keylget tk_utils_defs ELEMENT.PLOT_WIDTH] \
		    [keylget tk_utils_defs ELEMENT.PLOT_HEIGHT]]
	    set c_win [keylget element_info container_win]
	    set c_id [keylget element_info container_id]
	    set e_win [keylget element_info element_win]
	    set e_id [keylget element_info element_id]
	    set orientation [keylget element_info orientation]

	    update idletasks

	    nip_author_test plot -element $c_win$e_win\
		    -container $c_win\
		    -seq_id $seq_id \
		    -result_id $result_id($i)\
		    -results $results($i)\
		    -container_id [keylget element_info container_id]\
		    -element_id [keylget element_info element_id]\
		    -element_type "CANVAS"\
		    -fill [lindex $col_list $i]\
		    -orientation $orientation

	    seqed_element_bindings $c_id $c_win$e_win $e_id
	}
    }
    fit_on_screen $c_win

    #update result list
    result_list_update $c_win 
    return 0
}

proc plot_codon_pref {seq_id r_id res strand} {
    global nip_defs spin_defs tk_utils_defs HORIZONTAL TOP_S

    upvar $r_id result_id $res results

    set type [keylget nip_defs CODONPREF]

    set col_list "[keylget nip_defs NIP.CODONPREF.COLOUR.F1] [keylget nip_defs NIP.CODONPREF.COLOUR.F2] [keylget nip_defs NIP.CODONPREF.COLOUR.F3]"
    
    for {set i 0} {$i < 3} {incr i} {
	if {$result_id($i) != -1} {
	    set frame [expr $i + 1]
	    set element_info [create_seq_element $seq_id -1 $type $strand $frame \
		    $HORIZONTAL XY "CANVAS" [keylget spin_defs CONTAINER.TITLE] \
		    [keylget tk_utils_defs ELEMENT.PLOT_WIDTH] \
		    [keylget tk_utils_defs ELEMENT.PLOT_HEIGHT]]
	    set c_win [keylget element_info container_win]
	    set c_id [keylget element_info container_id]
	    set e_win [keylget element_info element_win]
	    set e_id [keylget element_info element_id]
	    set orientation [keylget element_info orientation]

	    update idletasks

	    nip_codon_pref plot -element $c_win$e_win\
		    -container $c_win\
		    -seq_id $seq_id \
		    -result_id $result_id($i)\
		    -results $results($i)\
		    -container_id [keylget element_info container_id]\
		    -element_id [keylget element_info element_id]\
		    -element_type "CANVAS"\
		    -fill [lindex $col_list $i]\
		    -orientation $orientation

	    seqed_element_bindings $c_id $c_win$e_win $e_id
	}
    }
    fit_on_screen $c_win
    #update result list
    result_list_update $c_win 
}


proc CodonPref { } {
    global nip_defs

    set w .codon_pref
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Codon preferences"
    
    GeneSearch_CTable $w load_optional

    keylset m MODE [keylget nip_defs NIP.PGS.MODE]
    set b1 [keylget m MODE.BUTTON.1]
    set b2 [keylget m MODE.BUTTON.2]

    frame $w.m

    global $w.mode1 $w.mode2
    set $w.mode1 [keylget m MODE.VALUE1]
    set $w.mode2 [keylget m MODE.VALUE2]
    checkbutton $w.m.type1 \
	    -text [keylget m MODE.BUTTON.3] \
	    -variable $w.mode1
    checkbutton $w.m.type2 \
	    -text [keylget m MODE.BUTTON.4] \
	    -variable $w.mode2

    pack $w.m.type1 $w.m.type2 -anchor w
    pack $w.m -fill x

    #########################################################################
    #strand selection
    strand_both $w.strand    
    pack $w.strand -fill x 

    frame $w.separator -bd 2 -relief raised -height 2

    global $w.sc
    keylset s STOP [keylget nip_defs NIP.PGS.STOP]
    set $w.sc [keylget s STOP.VALUE]
    frame $w.s
    checkbutton $w.s.stop \
	    -text [keylget s STOP.NAME] \
	    -variable $w.sc

    pack $w.s.stop -side left

    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $w.button -bd 2 -relief groove \
	-ok_command "CodonPref2 $w $w.win_len $w.range $w.fname \[strand_get $w.strand\]"\
	-cancel_command "seq_id_destroy $w.range; destroy $w" \
	-help_command "show_help spin {SPIN-Codon-Usage-Method}"
    
    pack $w.separator -side top -fill x -padx 10 -pady 5
    pack $w.s $w.button -fill x 
}

proc CodonPref2 { w win_len range fname strand} {
    global nip_defs tk_utils_defs PROTEIN TOP_S BOTTOM_S
    global $w.sc $w.mode1 $w.mode2

    if {-1 == [GeneSearchInputCheck $range $win_len]} {
	return 
    }
    set win_length [entrybox_get $win_len]
    set codon_table [getFname_in_name $fname]

    #check if valid codon table
    if {$codon_table != ""} {
	if {![valid_codon_table -codon_table $codon_table]} {
	    verror ERR_WARN "codon preference" "Invalid codon table"
	    raise $w
	    return
	}
    }

    set mode1 0
    set mode2 0
    if {[set $w.mode1]} {
	set mode1 2
    }
    if {[set $w.mode2]} {
	set mode2 4
    }
    set option [expr $mode1 + $mode2]

    set seq_id [name_to_seq_id [seq_id_name $range]]

    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "codon preference" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    SetBusy

    if {[expr $strand & $TOP_S]} {
	
	set res [nip_codon_pref create -win_len $win_length \
		-start [seq_id_from $range] \
		-end [seq_id_to $range] \
		-codon_table $codon_table \
		-option $option \
		-seq_id $seq_id\
		-strand $TOP_S]
	set n 0 
	for {set i 0} {$i < 3} {incr i} {
	    set result_id($i) [lindex $res $n]
	    incr n
	    set results($i) [lindex $res $n]
	    incr n
	}

	if {-1 == [plot_codon_pref $seq_id result_id results $TOP_S]} {
	    ClearBusy
	    seq_id_destroy $range
	    return
	}

	if {[set $w.sc]} {
	    NipStopCodon2 "" $range $TOP_S
	}
    }
    if {[expr $strand & $BOTTOM_S]} {
	
	set res [nip_codon_pref create -win_len $win_length \
		-start [seq_id_from $range] \
		-end [seq_id_to $range] \
		-codon_table $codon_table \
		-option $option \
		-seq_id $seq_id\
		-strand $BOTTOM_S]
	set n 0 
	for {set i 0} {$i < 3} {incr i} {
	    set result_id($i) [lindex $res $n]
	    incr n
	    set results($i) [lindex $res $n]
	    incr n
	}

	if {-1 == [plot_codon_pref $seq_id result_id results $BOTTOM_S]} {
	    ClearBusy
	    seq_id_destroy $range
	    return
	}

	if {[set $w.sc]} {
	    NipStopCodon2 "" $range $BOTTOM_S
	}
    }

    # stop windows from hiding the plot
    wm withdraw $w
	
    ClearBusy

    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
    destroy $w
}

proc AuthorTest { } {
    global nip_defs

    set w .author_test
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Author test"
    
    #GeneSearch_CTable $w load

    keylset wl ERROR [keylget nip_defs NIP.PGS.ERROR]

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

    entrybox $w.error \
	-title "[keylget wl ERROR.NAME] ([keylget wl ERROR.MIN] to [keylget wl ERROR.MAX])" \
	-default [keylget wl ERROR.VALUE] \
	-width 5 \
	-type "CheckFloatRange [keylget wl ERROR.MIN] [keylget wl ERROR.MAX]"
	    		 
    keylset us RANGE [keylget nip_defs NIP.PGS.RANGE]
    seq_id $w.range -range 1 -browse 1 -from 1 -to $seq_length \
	-start_value $seq_start -end_value $seq_end -min_value 1 \
	-default [seq_info $seq_id name]\
	-update_cmd [list [list seq_range_updates $w.range]]\
	-browse_cmd nip_seq_browser

    #get codon table file name
    keylset ct C_TABLE [keylget nip_defs NIP.PGS.C_TABLE]
    eval getFname $w.fname [list [keylget ct C_TABLE.NAME]] load {} \
	    [keylget ct C_TABLE.VALUE]

    #########################################################################
    #strand selection
    strand_both $w.strand    

    pack $w.range $w.error $w.fname $w.strand -anchor w -fill x

    frame $w.separator -bd 2 -relief raised -height 2

    global $w.sc
    keylset s STOP [keylget nip_defs NIP.PGS.STOP]
    set $w.sc [keylget s STOP.VALUE]
    frame $w.s
    checkbutton $w.s.stop \
	    -text [keylget s STOP.NAME] \
	    -variable $w.sc

    pack $w.s.stop -side left

    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $w.button -bd 2 -relief groove \
	-ok_command "AuthorTest2 $w $w.error $w.range $w.fname \[strand_get $w.strand\]" \
	-cancel_command "seq_id_destroy $w.range; destroy $w" \
	-help_command "show_help spin {SPIN-Author-Test}"

    pack $w.separator -side top -fill x -padx 10 -pady 5
    pack $w.s $w.button -fill x
}

proc AuthorTest2 { w error range fname strand} {
    global tk_utils_defs nip_defs $w.sc PROTEIN TOP_S BOTTOM_S

    set error [entrybox_get $error]
    set codon_table [getFname_in_name $fname]
    #check if valid codon table
    if {![valid_codon_table -codon_table $codon_table]} {
	verror ERR_WARN "author test" "Invalid codon table"
	return
    }

    set seq_id [name_to_seq_id [seq_id_name $range]]

    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "author test" "unable to process protein sequences"
	seq_id_destroy $range
	return
    }

    SetBusy

    if {[expr $strand & $TOP_S]} {
	set res [nip_author_test create -error $error \
		-start [seq_id_from $range] \
		-end [seq_id_to $range] \
		-codon_table $codon_table \
		-seq_id $seq_id\
		-strand $TOP_S]

	set n 0 
	for {set i 0} {$i < 3} {incr i} {
	    set result_id($i) [lindex $res $n]
	    incr n
	    set results($i) [lindex $res $n]
	    incr n
	}

	# stop windows from hiding the plot
	wm withdraw $w
	
	if {-1 == [plot_author_test $seq_id result_id results $TOP_S]} {
	    ClearBusy
	    seq_id_destroy $range
	    return
	}
    
	if {[set $w.sc]} {
	    NipStopCodon2 "" $range $TOP_S
	}
    }
    if {[expr $strand & $BOTTOM_S]} {
	set res [nip_author_test create -error $error \
		-start [seq_id_from $range] \
		-end [seq_id_to $range] \
		-codon_table $codon_table \
		-seq_id $seq_id\
		-strand $BOTTOM_S]

	set n 0 
	for {set i 0} {$i < 3} {incr i} {
	    set result_id($i) [lindex $res $n]
	    incr n
	    set results($i) [lindex $res $n]
	    incr n
	}

	# stop windows from hiding the plot
	wm withdraw $w
	
	if {-1 == [plot_author_test $seq_id result_id results $BOTTOM_S]} {
	    ClearBusy
	    seq_id_destroy $range
	    return
	}
    
	if {[set $w.sc]} {
	    NipStopCodon2 "" $range $BOTTOM_S
	}
    }

    ClearBusy

    seq_id_destroy $range

    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range]
    set $seq_id.end [seq_id_to $range]
    
    destroy $w
}

