##############################################################################
#called from C
#delete reading coverage histogram
proc DeleteReadingCoverage {io c_win rc_win cons_id} {
    global $rc_win.rcov_id $c_win.row

    consistency_result_list_update $io $c_win $cons_id

    set new_height [expr [winfo height $c_win] - [winfo height $rc_win]]

    unset $rc_win.rcov_id
    set id [get_consistency_window_id $rc_win]

    array set info [grid info $c_win.rc$id]
    grid_delete $c_win row $info(-row) 1

    #need to set the weight of the deleted row to 0
    grid rowconfigure $c_win [set $c_win.row] -weight 0 -minsize 0
    incr $c_win.row -1

    $rc_win delete all
    destroy $rc_win
    destroy $c_win.rc$id
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
#user interface dialog box for reading coverage histogram
proc ReadingCoverage { io } {
    global gap_defs

    set f [keylget gap_defs READING_COVERAGE.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "reading coverage"

    contig_id $f.id \
	    -io $io \
	    -range 1
    
    lorf_in $f.infile [keylget gap_defs READING_COVERAGE.INFILE] \
	    "{contig_id_configure $f.id -state disabled} \
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state normal}" -bd 2 -relief groove

    ###########################################################################
    #strand selection
    keylset st STRAND [keylget gap_defs READING_COVERAGE.STRAND]
    set b1 [keylget st STRAND.BUTTON.1]
    set b2 [keylget st STRAND.BUTTON.2]
    set b3 [keylget st STRAND.BUTTON.3]
    set b4 [keylget st STRAND.BUTTON.4]

    radiolist $f.strand \
	    -title [keylget st STRAND.NAME]\
	    -bd 2 \
	    -relief groove \
	    -default [keylget st STRAND.VALUE] \
	    -buttons [format { \
	    { %s } { %s } { %s } { %s } } \
	    [list $b1] [list $b2] [list $b3] [list $b4]]

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "ReadingCoverage2 $io $f $f.id $f.infile $f.strand" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Consistency-ReadingCov}" \
	    -bd 2 \
	    -relief groove

    pack $f.infile $f.id $f.strand $f.ok_cancel -side top -fill both
}

##############################################################################
#stand alone quality display
proc ReadingCoverage2 { io f id infile strand} {
    
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

    # stop windows from hiding the plot
    destroy $f

    CreateNewReadingCoverage $io $contig_list $str
}

##############################################################################
proc AddReadingCoverageCrossHair {io rcov_id c_win rc_win y} {
    global $c_win.cursor

    if {[set $c_win.cursor]} {
	draw_canvas_cursor_y -io $io -id $rcov_id -y [$rc_win canvasy $y]
    } else {
	$rc_win delete cursor_y
    }
}

##############################################################################
proc SetReadingCoverageBindings { io c_win rc_win} {
    global $c_win.cons_id $rc_win.rcov_id

    bind $rc_win <Any-Leave> "+delete_canvas_cursor -io $io -id [set $c_win.cons_id]; $c_win.brief configure -text \"\""

    bind $rc_win <Any-Motion> "AddConsistencyCrossHair $io [set $c_win.cons_id] $c_win $rc_win %x"

    bind $rc_win <Any-Motion> "+AddReadingCoverageCrossHair $io [set $rc_win.rcov_id] $c_win $rc_win %y"

    bind $rc_win <Any-Motion> "+ $c_win.brief configure -text \"Reading coverage (#[set $rc_win.rcov_id])\""

    # Double button 1 or 2 to move or create an editor
    bind $rc_win <<move-create>> "
	consistency_cursor_editor $io $rc_win [set $c_win.cons_id] [set $rc_win.rcov_id] %x
    "
    bind $rc_win <<use>> "
	consistency_cursor_editor $io $rc_win [set $c_win.cons_id] [set $rc_win.rcov_id] %x
    "
}

##############################################################################
#called from main gap4 menu. Create new display
proc CreateNewReadingCoverage {io contig_list strand} {
    global gap_defs 

    set result [CreateConsistencyDisplay $io $contig_list]
    set c_win [keylget result cons_win]
    set cons_id [keylget result cons_id]

    CreateReadingCoverage $io $cons_id $c_win $strand
}

##############################################################################
proc ConsistencyReadingCoverage {io cons_id c_win } {
    global gap_defs

    set f [keylget gap_defs READING_COVERAGE.WIN]2
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "reading coverage"

    ###########################################################################
    #strand selection
    keylset st STRAND [keylget gap_defs READING_COVERAGE.STRAND]
    set b1 [keylget st STRAND.BUTTON.1]
    set b2 [keylget st STRAND.BUTTON.2]
    set b3 [keylget st STRAND.BUTTON.3]
    set b4 [keylget st STRAND.BUTTON.4]

    radiolist $f.strand \
	    -title [keylget st STRAND.NAME]\
	    -bd 2 \
	    -relief groove \
	    -default [keylget st STRAND.VALUE] \
	    -buttons [format { \
	    { %s } { %s } { %s } { %s } } \
	    [list $b1] [list $b2] [list $b3] [list $b4]]

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "ConsistencyReadingCoverage2 $f $io $cons_id $c_win $f.strand" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Reading Coverage:Introduction}" \
	    -bd 2 \
	    -relief groove

    pack $f.strand $f.ok_cancel -side top -fill both

}

proc ConsistencyReadingCoverage2 {f io cons_id c_win strand } {

    
    set str [radiolist_get $strand]

    # stop windows from hiding the plot
    destroy $f

    CreateReadingCoverage $io $cons_id $c_win $str
}

##############################################################################
proc CreateReadingCoverage {io cons_id c_win strand} {
    global gap_defs $c_win.row

    set orig_height [winfo height $c_win]

    set id [next_consistency_window]

    set rc_win $c_win[keylget gap_defs READING_COVERAGE.WIN]$id
    global $rc_win.rcov_id

    set scroll x
    set height [keylget gap_defs READING_COVERAGE.PLOT_HEIGHT]
    set borderwidth [keylget gap_defs READING_COVERAGE.BORDERWIDTH]
    set width [keylget gap_defs CONSISTENCY_DISPLAY.PLOT_WIDTH]

    set zoom_cmd [list "consistency_zoom $io $cons_id \[set $rc_win.rcov_id\] $scroll"]

    #allow -height and -width to have affect
    wm geometry $c_win {}

    set row_num [c_get_next_row $c_win]

    ##########################################################################
    #create reading coverage histogram
    frame $c_win.rc$id -bd $borderwidth -relief groove
    canvasbox $rc_win -width $width -height $height \
	-bd 0 -highlightthickness 0 \
	-xscrollcommand "$c_win.hscroll set" \
	-yscrollcommand "$c_win.vscroll$id set" \
	-zoom_command $zoom_cmd 
    
    set width [keylget gap_defs CONSISTENCY_DISPLAY.RULER.PLOT_HEIGHT]
    canvasbox $c_win.vruler$id -width $width -height $height

    scrollbar $c_win.vscroll$id -relief sunken -command "gc_scroll_y $io \[set $rc_win.rcov_id\]"

    #set toplevel geometry when resized window
    set new_height [expr [winfo height $c_win] + [winfo reqheight $rc_win] + \
	    2 * $borderwidth]

    grid_insert $c_win row $row_num 1
    grid columnconfig $c_win 1 -weight 1
    grid rowconfig $c_win $row_num -weight 100

    #need to pack rc_win in a frame to allow a border around it
    grid $c_win.rc$id -row $row_num -column 1 -sticky nesw

    pack $rc_win -in $c_win.rc$id -padx [keylget gap_defs CONSISTENCY_DISPLAY.PADX] -pady [keylget gap_defs READING_COVERAGE.PADY] -fill both -expand yes

    #grid $c_win.vscroll$id -row $row_num -column 2 -sticky ns
    grid $c_win.vruler$id -row $row_num -column 0 -sticky ns -pady [keylget gap_defs READING_COVERAGE.PADY]

    #update geometry so the toplevel window can't grow larger than the screen
    update_geom $c_win $rc_win

    #register reading coverage histogram and do first display
    set $rc_win.rcov_id \
	    [reading_coverage -io $io -frame $c_win -window $rc_win -id $cons_id -win_ruler $c_win.vruler$id -strand $strand]

    if {[set $rc_win.rcov_id] == -1} {
	verror ERR_WARN "Gap4" "Too many windows"

	grid_delete $c_win row $row_num 1

	#need to set the weight of the deleted row to 0
	grid rowconfigure $c_win $row_num -weight 0 -minsize 0
	incr $c_win.row -1

	$rc_win delete all
	destroy $rc_win
	destroy $c_win.rc$id
	destroy $c_win.vscroll$id
	destroy $c_win.vruler$id

	#shrink window
	wm geometry $c_win [winfo width $c_win]x$orig_height
	#update_geom $c_win $c_win
	update idletasks
	return
    }


    consistency_result_list_update $io $c_win $cons_id

    #bind the configure actions to the toplevel
    bind $c_win <Any-Configure> "+
    if {\[winfo toplevel %W\] == \"%W\"} {
	update idletasks
	resize_canvas -io $io -id [set $rc_win.rcov_id]
    }
    "

    SetCanvasBindings $rc_win $zoom_cmd    
    SetReadingCoverageBindings $io $c_win $rc_win
}
