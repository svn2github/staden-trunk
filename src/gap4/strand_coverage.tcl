##############################################################################
#called from C
#delete strand coverage histogram
proc DeleteStrandCoverage {io c_win s_win cons_id} {
    global $s_win.strand_id $c_win.row

    consistency_result_list_update $io $c_win $cons_id

    set new_height [expr [winfo height $c_win] - [winfo height $s_win]]

    unset $s_win.strand_id
    set id [get_consistency_window_id $s_win]

    array set info [grid info $c_win.s$id]
    grid_delete $c_win row $info(-row) 1

    #need to set the weight of the deleted row to 0
    grid rowconfigure $c_win [set $c_win.row] -weight 0 -minsize 0
    incr $c_win.row -1

    $s_win delete all
    destroy $s_win
    destroy $c_win.s$id

    #shrink window
    wm geometry $c_win [winfo width $c_win]x$new_height

    #removed last row
    if {[set $c_win.row] == 1} {
	wm resizable $c_win 1 0
    }
    update idletasks
}

##############################################################################
#user interface dialog box for strand coverage histogram
proc StrandCoverage { io } {
    global gap_defs

    set f [keylget gap_defs STRAND_COVERAGE.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "strand coverage"

    contig_id $f.id \
	    -io $io \
	    -range 1
    
    lorf_in $f.infile [keylget gap_defs STRAND_COVERAGE.INFILE] \
	    "{contig_id_configure $f.id -state disabled} \
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state normal}" -bd 2 -relief groove

    ###########################################################################
    #strand selection
    keylset st STRAND [keylget gap_defs STRAND_COVERAGE.STRAND]
    set b1 [keylget st STRAND.BUTTON.1]
    set b2 [keylget st STRAND.BUTTON.2]
    set b3 [keylget st STRAND.BUTTON.3]

    radiolist $f.strand \
	    -title [keylget st STRAND.NAME]\
	    -bd 2 \
	    -relief groove \
	    -default [keylget st STRAND.VALUE] \
	    -buttons [format { \
	    { %s } { %s } { %s } } \
	    [list $b1] [list $b2] [list $b3]]

    ###########################################################################
    #problem selection
    keylset st PROBLEM [keylget gap_defs STRAND_COVERAGE.PROBLEM]
    set b1 [keylget st PROBLEM.BUTTON.1]
    set b2 [keylget st PROBLEM.BUTTON.2]

    radiolist $f.problem \
	    -title [keylget st PROBLEM.NAME]\
	    -bd 2 \
	    -relief groove \
	    -default [keylget st PROBLEM.VALUE] \
	    -buttons [format { \
	    { %s } { %s } } \
	    [list $b1] [list $b2]]
    

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "StrandCoverage2 $io $f $f.id $f.infile $f.strand $f.problem" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Consistency-StrandCov}" \
	    -bd 2 \
	    -relief groove

    pack $f.infile $f.id $f.strand $f.problem $f.ok_cancel -side top -fill both
}

##############################################################################
#stand alone quality display
proc StrandCoverage2 { io f id infile strand problem} {
    
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

    set str [radiolist_get $strand]
    set prob [radiolist_get $problem]

    # stop windows from hiding the plot
    destroy $f

    CreateNewStrandCoverage $io $contig_list $str $prob
}

##############################################################################
proc AddStrandCoverageCrossHair {io strand_id c_win s_win y} {
    global $c_win.cursor

    if {[set $c_win.cursor]} {
	draw_canvas_cursor_y -io $io -id $strand_id -y [$s_win canvasy $y]
    } else {
	$s_win delete cursor_y
    }
}

##############################################################################
proc SetStrandCoverageBindings { io c_win s_win} {
    global $c_win.cons_id $s_win.strand_id

    bind $s_win <Any-Leave> "+delete_canvas_cursor -io $io -id [set $c_win.cons_id]; $c_win.brief configure -text \"\""

    bind $s_win <Any-Motion> "AddConsistencyCrossHair $io [set $c_win.cons_id] $c_win $s_win %x"

    bind $s_win <Any-Motion> "+AddStrandCoverageCrossHair $io [set $s_win.strand_id] $c_win $s_win %y"

    bind $s_win <Any-Motion> "+ $c_win.brief configure -text \"Strand coverage  (#[set $s_win.strand_id])\""

    # Double button 1 or 2 to move or create an editor
    bind $s_win <<move-create>> "
	consistency_cursor_editor $io $s_win [set $c_win.cons_id] [set $s_win.strand_id] %x
    "
    bind $s_win <<use>> "
	consistency_cursor_editor $io $s_win [set $c_win.cons_id] [set $s_win.strand_id] %x
    "
}

##############################################################################
#called from main gap4 menu. Create new display
proc CreateNewStrandCoverage {io contig_list strand problem} {
    global gap_defs 

    set result [CreateConsistencyDisplay $io $contig_list]
    set c_win [keylget result cons_win]
    set cons_id [keylget result cons_id]

    CreateStrandCoverage $io $cons_id $c_win $strand $problem
}

##############################################################################
proc ConsistencyStrandCoverage {io cons_id c_win } {
    global gap_defs

    set f [keylget gap_defs STRAND_COVERAGE.WIN]2
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "strand coverage"

    ###########################################################################
    #strand selection
    keylset st STRAND [keylget gap_defs STRAND_COVERAGE.STRAND]
    set b1 [keylget st STRAND.BUTTON.1]
    set b2 [keylget st STRAND.BUTTON.2]
    set b3 [keylget st STRAND.BUTTON.3]

    radiolist $f.strand \
	    -title [keylget st STRAND.NAME]\
	    -bd 2 \
	    -relief groove \
	    -default [keylget st STRAND.VALUE] \
	    -buttons [format { \
	    { %s } { %s } { %s } } \
	    [list $b1] [list $b2] [list $b3]]

    ###########################################################################
    #problem selection
    keylset st PROBLEM [keylget gap_defs STRAND_COVERAGE.PROBLEM]
    set b1 [keylget st PROBLEM.BUTTON.1]
    set b2 [keylget st PROBLEM.BUTTON.2]

    radiolist $f.problem \
	    -title [keylget st PROBLEM.NAME]\
	    -bd 2 \
	    -relief groove \
	    -default [keylget st PROBLEM.VALUE] \
	    -buttons [format { \
	    { %s } { %s } } \
	    [list $b1] [list $b2]]
    

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "ConsistencyStrandCoverage2 $f $io $cons_id $c_win $f.strand $f.problem" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Strand Coverage:Introduction}" \
	    -bd 2 \
	    -relief groove

    pack $f.strand $f.problem $f.ok_cancel -side top -fill both
}

proc ConsistencyStrandCoverage2 {f io cons_id c_win strand problem} {

    set str [radiolist_get $strand]
    set prob [radiolist_get $problem]

    # stop windows from hiding the plot
    destroy $f    

    CreateStrandCoverage $io $cons_id $c_win $str $prob
}

##############################################################################
proc CreateStrandCoverage {io cons_id c_win strand problem} {
    global gap_defs $c_win.row

    set orig_height [winfo height $c_win]

    set id [next_consistency_window]

    set s_win $c_win[keylget gap_defs STRAND_COVERAGE.WIN]$id

    global $s_win.strand_id

    set scroll x
    set height [keylget gap_defs STRAND_COVERAGE.PLOT_HEIGHT]
    set borderwidth [keylget gap_defs STRAND_COVERAGE.BORDERWIDTH]
    set width [keylget gap_defs CONSISTENCY_DISPLAY.PLOT_WIDTH]

    set zoom_cmd [list "consistency_zoom $io $cons_id \[set $s_win.strand_id\] $scroll"]

    set row_num [c_get_next_row $c_win]

    ##########################################################################
    #create strand coverage histogram
    frame $c_win.s$id -bd $borderwidth -relief groove
    canvasbox $s_win -width $width -height $height \
	-bd 0 -highlightthickness 0 \
	-xscrollcommand "$c_win.hscroll set" \
	-zoom_command $zoom_cmd 
    
    grid_insert $c_win row $row_num 1
    grid columnconfig $c_win 1 -weight 1
    grid rowconfig $c_win $row_num -weight 1

    #need to pack s_win in a frame to allow a border around it
    grid $c_win.s$id -row $row_num -column 1 -sticky nesw

    pack $s_win -in $c_win.s$id -padx [keylget gap_defs CONSISTENCY_DISPLAY.PADX] -pady [keylget gap_defs STRAND_COVERAGE.PADY] -fill both -expand yes

    #update geometry so the toplevel window can't grow larger than the screen
    update_geom $c_win $s_win

    #register strand coverage histogram and do first display
    set $s_win.strand_id \
	    [strand_coverage -io $io -frame $c_win -window $s_win -id $cons_id -strand $strand -problems $problem]

    if {[set $s_win.strand_id] == -1} {
	verror ERR_WARN "Gap4" "Too many windows"

	grid_delete $c_win row $row_num 1

	#need to set the weight of the deleted row to 0
	grid rowconfigure $c_win $row_num -weight 0 -minsize 0
	incr $c_win.row -1

	$s_win delete all
	destroy $s_win
	destroy $c_win.s$id

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
	resize_canvas -io $io -id [set $s_win.strand_id]
    }
    "

    SetCanvasBindings $s_win $zoom_cmd    
    SetStrandCoverageBindings $io $c_win $s_win
}
