#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc NipRestrictionEnzymeDialogue { } {
    global tk_utils_defs nip_defs

    set w [keylget tk_utils_defs R_ENZ.WIN]
    if {[xtoplevel $w] == ""} return
    wm title $w "Select restriction enzymes"

    global $w.list

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

    if {![seq_renzbox $w $seq_id 1 $seq_length $seq_start $seq_end]} {
	#pressed cancel
	return
    }

    set list [set $w.list]
    set filename [renzbox_filename $w]
    set num_items [llength $list]

    set results [renzbox_get_seqid $w]

    NipCreateREnzDisplay $list $filename \
	[name_to_seq_id [lindex $results 0]] [lindex $results 1] \
	[lindex $results 2]
}

proc nip_scroll_x {id command x args} {
    nip_scroll_canvas -id $id -xscrollcommand "$command $x $args"
}

proc nip_scroll_y {id command y args} {
    nip_scroll_canvas -id $id -yscrollcommand "$command $y $args"
}

proc nip_zoom {seq_id id args } {
    eval nip_zoom_canvas -seq_id $seq_id -id $id $args
}

proc nip_zoomback {seq_id id } {
    nip_zoom_canvas -seq_id $seq_id -id $id
}

proc nip_resize {id} {
    if {$id != -1} {
	nip_resize_canvas -id $id
    }
}

proc NipCreateREnzDisplay {list filename seq_id from to} {
    global tk_utils_defs PROTEIN

    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "restriction enzymes" "unable to process protein sequences"
	return
    }

    set w .nip_r_enzyme_map
    #generate a new display
    set num_display [next_renz_display]
    set w $w$num_display
   
    global $w.renz_id
    set $w.renz_id -1

    set x_scroll "nip_scroll_x \[set $w.renz_id\]"
    set y_scroll "nip_scroll_y \[set $w.renz_id\]"
    set zoomback_command "nip_zoomback $seq_id \[set $w.renz_id\]"
    set zoom_command "nip_zoom $seq_id \[set $w.renz_id\]"
    set cursor_command "nip_canvas_cursor_x -id \[set $w.renz_id\]"

    set invoke_command "nip_canvas_cursor_editor $seq_id \[set $w.renz_id\]"
    #set invoke_command "NipCreateCanvasEditorCursor $w \[set $w.renz_id\]"

    set renz_name_command "nip_get_renz_name -id \[set $w.renz_id\]"
    set renz_info_command "nip_get_renz_info -id \[set $w.renz_id\]"

    set config_command "nip_resize \[set $w.renz_id\]"

    set min_win_ht [keylget tk_utils_defs R_ENZ.PLOT_HEIGHT]
    set num_enz [llength $list]
    set tick_ht [keylget tk_utils_defs R_ENZ.TICK_HEIGHT]

    #HACK - want a better way of deciding on the text offset
    set text_offset [expr $tick_ht * 1.5]

    #set optimal window height for names and plot
    set height [expr ($num_enz + 2) * $tick_ht]

    #HACK to ensure window does not get too big!
    if { $height > $min_win_ht} {
	set height $min_win_ht
    }

    # stop windows from hiding the plot
    # wm withdraw $w FIXME
    renz_map $w -zoom_command $zoom_command \
	-zoomback_command $zoomback_command\
	-scrollbar_x_cmd $x_scroll -scrollbar_y_cmd $y_scroll \
	-cursor_cmd $cursor_command \
	-renz_name_cmd $renz_name_command \
	-renz_info_cmd $renz_info_command \
	-width [keylget tk_utils_defs R_ENZ.PLOT_WIDTH]\
	-height $height \
	-ruler_height [keylget tk_utils_defs R_ENZ.RULER.PLOT_HEIGHT]\
	-names_width [keylget tk_utils_defs R_ENZ.NAME_WIDTH] \
	-selectbackground [keylget tk_utils_defs R_ENZ.SELECT_COLOUR] \
	-config_cmd $config_command \
	-invoke_cmd $invoke_command \
	-tick_ht [keylget tk_utils_defs R_ENZ.TICK_HEIGHT] \
	-text_offset $text_offset \
	-text_fill [keylget tk_utils_defs R_ENZ.TEXT_COLOUR]
 
    # Main Menu Bar
    NipCreateREnzMenu $w

    wm protocol $w WM_DELETE_WINDOW "NipREnzStartShutdown $w"

#    update idletasks
#    19/1/99 johnt - need just update for WINNT
    update
    SetBusy

    set $w.renz_id [renz_map_plot $w "nip_plot_renz -enzymes {$list} -num_enzymes [llength $list] -file {$filename} -seq_id $seq_id -start $from -end $to]"]

    ClearBusy
    if {[set $w.renz_id] == -1} {
	return
    }

    wm title $w "SPIN Restriction Enzyme Map (\#[set $w.renz_id])"
    #update result list
    seq_result_list_update [keylget tk_utils_defs RASTER.RESULTS.WIN]

    global $seq_id.start $seq_id.end
    set $seq_id.start $from
    set $seq_id.end $to

    wm minsize $w [winfo width $w] [winfo height $w]
    wm maxsize $w [winfo screenwidth $w] [winfo height $w]
   
}

##############################################################################
#stand alone restriction enzyme display
#menu display
proc NipCreateREnzMenu {w } {
    global $w.renz_id

    #set up a File menu
    menubutton $w.menubar.file -text "File" -menu $w.menubar.file.opts
    menu $w.menubar.file.opts
    $w.menubar.file.opts add command -label "Exit" \
	    -command "NipREnzStartShutdown $w"

    #set up a View menu
    menubutton $w.menubar.view -text "View" -menu $w.menubar.view.opts
    menu $w.menubar.view.opts
    $w.menubar.view.opts add command -label "Results manager" \
	    -command "rasterResultsManager {show_help spin {SPIN-Result-Manager}}"

    #set up a Results menu
    menubutton $w.menubar.results -text "Results" -menu $w.menubar.results.opts
    menu $w.menubar.results.opts -postcommand "NipGetREnzResults $w"

    #$w.menubar.results.opts add command -label "Output enzyme by enzyme" \
	-command "nip_renz_info -result_id \[set $w.renz_id\] -option 0"
    #$w.menubar.results.opts add command -label "Output ordered on position" \
	-command "nip_renz_info -result_id \[set $w.renz_id\] -option 1"

    menubutton $w.menubar.help -text "Help" -menu [set m $w.menubar.help.opts]
    menu $m
    $m add command -label "Introduction" \
	-command "show_help spin {SPIN-Restrict-Introduction}"
    $m add command -label "Selecting Enzymes" \
	-command "show_help spin {SPIN-Restrict-Selecting}"
    $m add command -label "Examining the Plot" \
	-command "show_help spin {SPIN-Restrict-Examining}"
    $m add command -label "Reconfiguring the Plot" \
	-command "show_help spin {SPIN-Restrict-Reconfig}"

    #do the packing
    pack $w.menubar.file $w.menubar.view $w.menubar.results -side left
    pack $w.menubar.help -side right

}

##############################################################################
proc NipGetREnzResults {w } {
    global $w.renz_id
    
    $w.menubar.results.opts delete 1 end
    
    result_list_popup_single [set $w.renz_id]\
	[seq_get_ops -index [set $w.renz_id]] \
	$w.menubar.results.opts

}

##############################################################################
#executed when the exit command is chosen from the File menu of stand alone
#restriction enzyme plot
proc NipREnzStartShutdown {w } {
    global $w.renz_id

    if {[info exists $w.renz_id]} {
	#nip_result_delete -result_id [set $w.renz_id]
	seq_result_update -index [set $w.renz_id] -job QUIT
    }
}
