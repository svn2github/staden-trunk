##############################################################################
#called from C
#delete readpair coverage histogram
proc DeleteReadPairCoverage {io c_win rp_win cons_id} {
    global $rp_win.rpair_id $c_win.row

    consistency_result_list_update $io $c_win $cons_id

    set new_height [expr [winfo height $c_win] - [winfo height $rp_win]]

    if {[info exists $rp_win.rpair_id]} {
	unset $rp_win.rpair_id
    }

    set id [get_consistency_window_id $rp_win]

    array set info [grid info $c_win.rp$id]
    grid_delete $c_win row $info(-row) 1

    #need to set the weight of the deleted row to 0
    grid rowconfigure $c_win [set $c_win.row] -weight 0 -minsize 0
    incr $c_win.row -1

    $rp_win delete all
    destroy $rp_win
    destroy $c_win.rp$id
    destroy $c_win.vscroll$id
    destroy $c_win.vruler$id

    #shrink window
    wm geometry $c_win [winfo width $c_win]x$new_height
    update idletasks
}

##############################################################################
#user interface dialog box for readpair coverage histogram
proc ReadPairCoverage { io } {
    global gap_defs

    set f [keylget gap_defs READPAIR_COVERAGE.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "readpair coverage"

    contig_id $f.id \
	    -io $io \
	    -range 1
    
    lorf_in $f.infile [keylget gap_defs READPAIR_COVERAGE.INFILE] \
	    "{contig_id_configure $f.id -state disabled} \
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state normal}" -bd 2 -relief groove

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "ReadPairCoverage2 $io $f $f.id $f.infile" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Consistency-ReadPairCov}" \
	    -bd 2 \
	    -relief groove

    pack $f.infile $f.id $f.ok_cancel -side top -fill both
}

##############################################################################
#readpair coverage
proc ReadPairCoverage2 { io f id infile} {
    
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

    # stop windows from hiding the plot
    destroy $f

    CreateNewReadPairCoverage $io $contig_list
}

##############################################################################
proc AddReadPairCoverageCrossHair {io rpair_id c_win rp_win y} {
    global $c_win.cursor

    if {[set $c_win.cursor]} {
	draw_canvas_cursor_y -io $io -id $rpair_id -y [$rp_win canvasy $y]
    } else {
	$rp_win delete cursor_y
    }
}

##############################################################################
proc SetReadPairCoverageBindings { io c_win rp_win} {
    global $c_win.cons_id $rp_win.rpair_id

    bind $rp_win <Any-Leave> "+delete_canvas_cursor -io $io -id [set $c_win.cons_id]; $c_win.brief configure -text \"\""

    bind $rp_win <Any-Motion> "AddConsistencyCrossHair $io [set $c_win.cons_id] $c_win $rp_win %x"

    bind $rp_win <Any-Motion> "+AddReadPairCoverageCrossHair $io [set $rp_win.rpair_id] $c_win $rp_win %y"

    bind $rp_win <Any-Motion> "+ $c_win.brief configure -text \"Readpair coverage (#[set $rp_win.rpair_id])\""

    # Double button 1 or 2 to move or create an editor
    bind $rp_win <<move-create>> "
	consistency_cursor_editor $io $rp_win [set $c_win.cons_id] [set $rp_win.rpair_id] %x
    "
    bind $rp_win <<use>> "
	consistency_cursor_editor $io $rp_win [set $c_win.cons_id] [set $rp_win.rpair_id] %x
    "
}

##############################################################################
#called from main gap4 menu. Create new display
proc CreateNewReadPairCoverage {io contig_list} {
    global gap_defs 

    set result [CreateConsistencyDisplay $io $contig_list]
    set c_win [keylget result cons_win]
    set cons_id [keylget result cons_id]

    CreateReadPairCoverage $io $cons_id $c_win
}

##############################################################################
proc CreateReadPairCoverage {io cons_id c_win} {
    global gap_defs $c_win.row

    #set toplevel geometry when resized window
    set orig_height [winfo height $c_win]

    set id [next_consistency_window]

    set rp_win $c_win[keylget gap_defs READPAIR_COVERAGE.WIN]$id
    global $rp_win.rpair_id 

    set scroll b
    set height [keylget gap_defs READPAIR_COVERAGE.PLOT_HEIGHT]
    set borderwidth [keylget gap_defs READPAIR_COVERAGE.BORDERWIDTH]
    set width [keylget gap_defs CONSISTENCY_DISPLAY.PLOT_WIDTH]

    set zoom_cmd [list "consistency_zoom $io $cons_id \[set $rp_win.rpair_id\] $scroll"]

    set row_num [c_get_next_row $c_win]

    ##########################################################################
    #create readpair coverage histogram
    frame $c_win.rp$id -bd $borderwidth -relief groove
    canvasbox $rp_win -width $width -height $height \
	-bd 0 -highlightthickness 0 \
	-xscrollcommand "$c_win.hscroll set" \
	-yscrollcommand "$c_win.vscroll$id set" \
	-zoom_command $zoom_cmd 
    
    set width [keylget gap_defs CONSISTENCY_DISPLAY.RULER.PLOT_HEIGHT]
    canvasbox $c_win.vruler$id -width $width -height $height

    scrollbar $c_win.vscroll$id -relief sunken -command "gc_scroll_y $io \[set $rp_win.rpair_id\]"

    grid_insert $c_win row $row_num 1

    grid columnconfig $c_win 1 -weight 1
    grid rowconfig $c_win $row_num -weight 100

    #need to pack rp_win in a frame to allow a border around it
    grid $c_win.rp$id -row $row_num -column 1 -sticky nesw

    pack $rp_win -in $c_win.rp$id -padx [keylget gap_defs CONSISTENCY_DISPLAY.PADX] -pady [keylget gap_defs READPAIR_COVERAGE.PADY] -fill both -expand yes


    grid $c_win.vscroll$id -row $row_num -column 2 -sticky ns
    grid $c_win.vruler$id -row $row_num -column 0 -sticky ns -pady [keylget gap_defs READPAIR_COVERAGE.PADY]

    #update geometry so the toplevel window can't grow larger than the screen
    update_geom $c_win $rp_win

    #register readpair coverage histogram and do first display
    set $rp_win.rpair_id \
	    [readpair_coverage -io $io -frame $c_win -window $rp_win -id $cons_id -win_ruler $c_win.vruler$id]


    #happens when no readpairs are present or there are too many windows up
    if {[set $rp_win.rpair_id] == -2} {
	grid rowconfig $c_win $row_num -weight 1
	return
    } 
    if {[set $rp_win.rpair_id] == -1} {
	verror ERR_WARN "Gap4" "Too many windows"

	grid_delete $c_win row $row_num 1

	#need to set the weight of the deleted row to 0
	grid rowconfigure $c_win $row_num -weight 0 -minsize 0
	incr $c_win.row -1

	$rp_win delete all
	destroy $rp_win
	destroy $c_win.rp$id
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
	resize_canvas -io $io -id [set $rp_win.rpair_id]
    }
    "
    SetCanvasBindings $rp_win $zoom_cmd    
    SetReadPairCoverageBindings $io $c_win $rp_win
}
