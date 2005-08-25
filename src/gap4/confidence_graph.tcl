#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

##############################################################################
#called from C
#delete confidence graph
proc DeleteConfidenceGraph {io c_win conf_win cons_id} {
    global $conf_win.conf_id $c_win.row

    consistency_result_list_update $io $c_win $cons_id

    set new_height [expr [winfo height $c_win] - [winfo height $conf_win]]

    unset $conf_win.conf_id

    set id [get_consistency_window_id $conf_win]

    array set info [grid info $c_win.qg$id]
    grid_delete $c_win row $info(-row) 1

    #need to set the weight of the deleted row to 0
    grid rowconfigure $c_win [set $c_win.row] -weight 0 -minsize 0
    incr $c_win.row -1

    $conf_win delete all
    destroy $conf_win
    destroy $c_win.qg$id
    destroy $c_win.vscroll$id
    destroy $c_win.vruler$id

    #shrink window
    wm geometry $c_win [winfo width $c_win]x$new_height

    #removed last row
    if {[set $c_win.row] == 1} {
	wm resizable $c_win 1 0
    }
    update idletasks
}

##############################################################################
#user interface dialog box for confidence values graph
proc ConfidenceGraph { io mode } {
    global gap_defs

    set f [keylget gap_defs CONFIDENCE_GRAPH.WIN]$mode
    if {[xtoplevel $f -resizable 0] == ""} return
    if {$mode == "conf"} {
	wm title $f "confidence graph"
    } elseif {$mode == "second"} {
	wm title $f "2nd-highest confidence graph"
    } else {
	wm title $f "diploid graph"
    }

    contig_id $f.id \
	    -io $io \
	    -range 1
    
    lorf_in $f.infile [keylget gap_defs CONFIDENCE_GRAPH.INFILE] \
	    "{contig_id_configure $f.id -state disabled} \
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state normal}" -bd 2 -relief groove

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "ConfGraph_OK_Pressed $io $mode $f $f.id $f.infile" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Consistency-Confidence}" \
	    -bd 2 \
	    -relief groove

    pack $f.infile $f.id $f.ok_cancel -side top -fill both
}

##############################################################################
#stand alone quality display
proc ConfGraph_OK_Pressed { io mode f id infile} {
    
    set contig_list ""
    if {[lorf_in_get $infile] == 4} {
	#single contig
        if {[set contign [contig_id_gel $id]] == ""} {bell; return}
       	if {[set lreg  [contig_id_lreg $id]] == ""} {bell; return}
	if {[set rreg  [contig_id_rreg $id]] == ""} {bell; return}
	set contig_list "{$contign $lreg $rreg}"
	SetContigGlobals $io $contign
    } elseif {[lorf_in_get $infile] == 3 } {
	#all contigs
	set contig_list [CreateAllContigList $io]
    } else {
	#list or file
	set contig_list [lorf_get_list $infile]
	set contig_list [remove_contig_duplicates -io $io -contigs $contig_list]
    }

    if {$contig_list == ""} {
	raise $f
	return
    }

    #set contig [contig_id_gel $id]
    #set lreg [contig_id_lreg $id]
    #set rreg [contig_id_rreg $id]
    #SetContigGlobals $io $contig $lreg $rreg

    # stop windows from hiding the plot
    destroy $f

    CreateNewConfidenceGraph $io $mode $contig_list

}

##############################################################################
proc AddConfidenceCrossHair {io conf_id c_win conf_win y} {
    global $c_win.cursor

    if {[set $c_win.cursor]} {
	draw_canvas_cursor_y -io $io -id $conf_id -y [$conf_win canvasy $y]
    } else {
	#delete_canvas_cursor -io $io -id $conf_id
	$conf_win delete cursor_y
    }
}

##############################################################################
proc SetConfidenceBindings { io c_win conf_win} {
    global $c_win.cons_id $conf_win.conf_id

    bind $conf_win <Any-Leave> "+delete_canvas_cursor -io $io -id [set $c_win.cons_id]; $c_win.brief configure -text \"\""

    bind $conf_win <Any-Leave> "+$conf_win delete cursor_y"

    bind $conf_win <Any-Motion> "AddConsistencyCrossHair $io [set $c_win.cons_id] $c_win $conf_win %x"

    bind $conf_win <Any-Motion> "+AddConfidenceCrossHair $io [set $conf_win.conf_id] $c_win $conf_win %y"

    bind $conf_win <Any-Motion> "+ $c_win.brief configure -text \"Confidence graph (#[set $conf_win.conf_id]) \""

    # Double button 1 and 2 to move or create an editor
    bind $conf_win <<move-create>> "
	consistency_cursor_editor $io $conf_win [set $c_win.cons_id] [set $conf_win.conf_id] %x
    "
    bind $conf_win <<use>> "
	consistency_cursor_editor $io $conf_win [set $c_win.cons_id] [set $conf_win.conf_id] %x
    "
}

##############################################################################
#called from main gap4 menu. Create new display
proc CreateNewConfidenceGraph {io mode contig_list} {
    global gap_defs 

    set result [CreateConsistencyDisplay $io $contig_list]
    set c_win [keylget result cons_win]
    set cons_id [keylget result cons_id]

    CreateConfidenceGraph $io $mode $cons_id $c_win
}

##############################################################################
proc CreateConfidenceGraph {io mode cons_id c_win} {
    global gap_defs $c_win.row

    set orig_height [winfo height $c_win]

    set id [next_consistency_window]

    set conf_win $c_win[keylget gap_defs CONFIDENCE_GRAPH.WIN]$id
    global $conf_win.conf_id 

    set scroll x
    set height [keylget gap_defs CONFIDENCE_GRAPH.PLOT_HEIGHT]
    set borderwidth [keylget gap_defs CONFIDENCE_GRAPH.BORDERWIDTH]
    set width [keylget gap_defs CONSISTENCY_DISPLAY.PLOT_WIDTH]

    set zoom_cmd [list "consistency_zoom $io $cons_id \[set $conf_win.conf_id\] $scroll"]

    set row_num [c_get_next_row $c_win]

    ##########################################################################
    #create confidence graph
    frame $c_win.qg$id -bd $borderwidth -relief groove
    canvasbox $conf_win -width $width -height $height \
	-bd 0 -highlightthickness 0 \
	-xscrollcommand "$c_win.hscroll set" \
	-yscrollcommand "$c_win.vscroll$id set" \
	-zoom_command $zoom_cmd 
    
    set width [keylget gap_defs CONSISTENCY_DISPLAY.RULER.PLOT_HEIGHT]
    canvasbox $c_win.vruler$id -width $width -height $height

    scrollbar $c_win.vscroll$id -relief sunken \
	-command "gc_scroll_y $io \[set $conf_win.conf_id\]"
    
    grid_insert $c_win row $row_num 1
    grid columnconfig $c_win 1 -weight 1
    grid rowconfig $c_win $row_num -weight 100

    #need to pack conf_win in a frame to allow a border around it
    grid $c_win.qg$id -row $row_num -column 1 -sticky nesw

    pack $conf_win -in $c_win.qg$id -padx [keylget gap_defs CONSISTENCY_DISPLAY.PADX] -pady [keylget gap_defs CONFIDENCE_GRAPH.PADY] -fill both -expand yes
    #grid $c_win.vscroll$id -row $row_num -column 2 -sticky ns

    grid $c_win.vruler$id -row $row_num -column 0 -sticky ns -pady [keylget gap_defs CONFIDENCE_GRAPH.PADY]

    #update geometry so the toplevel window can't grow larger than the screen
    update_geom $c_win $conf_win

    #register quality graph and do first display
    switch $mode {
	"conf" {
	    set $conf_win.conf_id [confidence_graph \
				       -io $io \
				       -frame $c_win \
				       -window $conf_win \
				       -id $cons_id \
				       -win_ruler $c_win.vruler$id]
	}
	"second" {
	    set $conf_win.conf_id [second_confidence_graph \
				       -io $io \
				       -frame $c_win \
				       -window $conf_win \
				       -id $cons_id \
				       -win_ruler $c_win.vruler$id]
	}
	"discrep" {
	    set $conf_win.conf_id [discrepancy_graph \
				       -io $io \
				       -frame $c_win \
				       -window $conf_win \
				       -id $cons_id \
				       -win_ruler $c_win.vruler$id \
				       -two_alleles 0]
	}
	"discrep2" {
	    set $conf_win.conf_id [discrepancy_graph \
				       -io $io \
				       -frame $c_win \
				       -window $conf_win \
				       -id $cons_id \
				       -win_ruler $c_win.vruler$id \
				       -two_alleles 1]
	}
    }

    
    if {[set $conf_win.conf_id] == -1} {
	verror ERR_WARN "Gap4" "Too many windows"

	grid_delete $c_win row $row_num 1

	#need to set the weight of the deleted row to 0
	grid rowconfigure $c_win $row_num -weight 0 -minsize 0
	incr $c_win.row -1

	$conf_win delete all
	destroy $conf_win
	destroy $c_win.qg$id
	destroy $c_win.vscroll$id
	destroy $c_win.vruler$id

	#shrink window
	wm geometry $c_win [winfo width $c_win]x$orig_height
	update idletasks
	return
    }

    consistency_result_list_update $io $c_win $cons_id

    #bind the configure actions to the toplevel
    bind $c_win <Any-Configure> "+
    if {\[winfo toplevel %W\] == \"%W\"} {
	update idletasks
	resize_canvas -io $io -id [set $conf_win.conf_id]
    }
    "

    SetCanvasBindings $conf_win $zoom_cmd    
    SetConfidenceBindings $io $c_win $conf_win
}
