#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc consistency_zoom {io cons_id id scroll args } {

    eval zoom_canvas -io $io -id $cons_id -r_id $id -direction $scroll $args
}

proc CreateConsistencyRuler {io t r_win hscroll} {
    global gap_defs config$t.ruler $t.cons_id config$t.ticks

    
    set config$t.ruler 1
    set config$t.ticks 1

    set NGContig(min_x) 100
    set scroll x
    set borderwidth [keylget gap_defs CONSISTENCY_DISPLAY.BORDERWIDTH]
    set height [keylget gap_defs CONSISTENCY_DISPLAY.RULER.PLOT_HEIGHT]
    set width [keylget gap_defs CONSISTENCY_DISPLAY.PLOT_WIDTH]

    #allow -height and -width to have affect
    wm geometry $t {}
    
    set zoom_cmd [list "gap_zoom $io \[set $t.cons_id\] x"]

    frame $t.r -bd 2 -relief groove
    canvasbox $r_win -width $width \
	    -height $height \
	    -xscrollcommand "$hscroll set" \
	    -highlightthickness 0 \
	    -bd 0 -zoom_command $zoom_cmd
    
    #set toplevel geometry when resized window
    set new_height [expr [winfo height $t] + [winfo reqheight $r_win] + \
	    2 * $borderwidth]
    
    set row_num [expr [GetConsistencyRow $t.hscroll] - 1]

    grid_insert $t row $row_num 1
    grid rowconfig $t $row_num
    grid $t.r -row $row_num -column 1 -sticky nsew
    pack $r_win -in $t.r -padx [keylget gap_defs CONSISTENCY_DISPLAY.PADX] -fill both -expand yes
    update idletasks

    #set toplevel geometry when resized window
    #check that things have been packed!
    if {[winfo height $t] != 1} {
	wm geometry $t [winfo width $t]x$new_height
    }
    SetCanvasBindings $r_win $zoom_cmd
}

proc DisplayConsistencyRuler {io t cons_id r_win hscroll} {
    global config$t.ruler
    global gap_defs

     if { [set config$t.ruler] } {
	#plot ruler

	#HACK to work round the horrible problem in zooming, turning ruler
	#off, scrolling and turning ruler on again.
	set scroll [lindex [$hscroll get] 0]
	if { ![winfo exists $r_win] } {
	    CreateConsistencyRuler $io $t $r_win $hscroll
	    SetConsistencyRulerBindings $io $t $r_win
	}
	#display_ruler -io $io -id $cons_id -win_ruler $r_win
	create_consistency_ruler -io $io -id $cons_id
	uplevel #0 eval [$hscroll cget -command] moveto $scroll
    } else {
	#delete ruler
	set new_height [expr [winfo height $t] - [winfo height $t.r]]
	
	#delete_window -io $io -id $cons_id -window $r_win
	delete_consistency_ruler -io $io -id $cons_id -window $r_win

	#check to see if ruler was last window in consistency display and
	#therefore the entire consistency display may have gone by this time
	if ![winfo exists $t] {
	    return
	}
	set row_num [GetConsistencyRow $t.hscroll]
	grid_delete $t row [expr $row_num - 1] 1

	#need to set the weight of the deleted row to 0
	grid rowconfigure $t [expr $row_num - 1] -weight 0 -minsize 0
	$r_win delete all
	destroy $r_win 
	destroy $t.r
	wm geometry $t [winfo width $t]x$new_height
    }
}

proc DisplayConsistencyRulerTicks {io t id hscroll } {
    global config$t.ticks
    global gap_defs

    set r_win $t[keylget gap_defs CONSISTENCY_DISPLAY.RULER.WIN]

    if {[set config$t.ticks]} {
	display_ruler_ticks -io $io -id $id -ticks 1
    } else {
	if {[winfo exists $r_win]} {
	    $r_win delete tick
	    ConsistencyRulerWindowSize 0 $t $r_win
	}
    }
}

proc ConsistencyRulerWindowSize {disp_ticks t r_win } {
    global gap_defs

    if {$disp_ticks} {
	set bbox [$r_win bbox tick]
    } else {
	set bbox [$r_win bbox contig]
    }

    if {[lindex $bbox 3] < [winfo height $r_win]} {
	
	if {[lindex $bbox 3] > [keylget gap_defs CONSISTENCY_DISPLAY.RULER.PLOT_HEIGHT]} {
	    set height [lindex $bbox 3]
	} else {
	    set height [keylget gap_defs CONSISTENCY_DISPLAY.RULER.PLOT_HEIGHT]
	}
	set new_height [expr abs([winfo height $r_win] - $height)]

	$r_win config -height $height
	set new_height [expr [winfo height $t] - $new_height]
	wm geometry $t [winfo width $t]x$new_height
	#update idletasks
    }

    if {[lindex $bbox 3] > [winfo height $r_win]} {

	set height [lindex $bbox 3]
	set new_height [expr abs([winfo height $r_win] - $height)]
	$r_win config -height $height
	set new_height [expr [winfo height $t] + $new_height]
	wm geometry $t [winfo width $t]x$new_height
	#update idletasks
    }
}

#zoom in buttons
proc ZoomInConsistency {amount zoom_cmd} {

    set l $amount
    set r [expr 1.0-$amount]

    #set w [winfo width $canvas]
    #set h [winfo height $canvas]
    #set $canvas.areaX1 [$canvas canvasx [expr $w*$l]]
    #set $canvas.areaX2 [$canvas canvasx [expr $w*$r]]
    #set $canvas.areaY1 [$canvas canvasx [expr $h*$l]]
    #set $canvas.areaY2 [$canvas canvasx [expr $h*$r]]
    #set $canvas.zoom 1

     #need to run the command at a global level
    uplevel #0 eval $zoom_cmd -amount $amount
   
}

##############################################################################
#consistency display
proc CreateConsistencyDisplay {io contig_list} {
    global gap_defs

    set num_display [next_consistency_display]

    set t [keylget gap_defs CONSISTENCY_DISPLAY.WIN]$num_display
    if {[xtoplevel $t] == ""} return
    fix_maxsize $t
    wm protocol $t WM_DELETE_WINDOW "ConsistencyStartShutdown $io $t"
    wm title $t "Consistency display: [lindex $contig_list 0] #[db_info get_read_num $io [lindex [lindex $contig_list 0] 0]]"

    global $t.cons_id

    set r_win $t[keylget gap_defs CONSISTENCY_DISPLAY.RULER.WIN]
    set borderwidth [keylget gap_defs CONSISTENCY_DISPLAY.BORDERWIDTH]
    set width [keylget gap_defs CONSISTENCY_DISPLAY.PLOT_WIDTH]

    #set scroll to be "both" because of zoom buttons
    set zoom_cmd [list "gap_zoom $io \[set $t.cons_id\] x"]
    
    scrollbar $t.hscroll -orient horizontal -relief sunken -command \
	    "gc_scroll_x $io \[set $t.cons_id\]"

    ##########################################################################
    # Main Menu Bar
    #frame $t.menubar -relief raised -borderwidth 2
    CreateConsistencyMenu $io $t $r_win

    ##########################################################################
    #button bar
    frame $t.buttons

    #zoom back button
    button $t.buttons.back -text "zoom out" -command "ZoomBackCanvas $io \[set $t.cons_id\]"
    button $t.buttons.zoomin10 -text "+10%" \
	-command "eval $zoom_cmd -amount 0.05"
    button $t.buttons.zoomin50 -text "+50%" \
	-command "eval $zoom_cmd -amount 0.1666" 

    #cursor checkbutton
    global $t.cursor
    checkbutton $t.buttons.cursor -text crosshairs -variable $t.cursor

    #cursor local, total and y position labels
    set cursor_t [keylget gap_defs CONSISTENCY_DISPLAY.CURSOR1]
    set cursor_l [keylget gap_defs CONSISTENCY_DISPLAY.CURSOR2]
    set cursor_y [keylget gap_defs CONSISTENCY_DISPLAY.CURSORY]
    label $t$cursor_t -bd 2 -relief sunken -width 6
    label $t$cursor_l -bd 2 -relief sunken -width 6
    label $t$cursor_y -bd 2 -relief sunken -width 6

    #information line
    label $t.brief

    pack $t.buttons.zoomin10 $t.buttons.zoomin50 $t.buttons.back \
	 $t.buttons.cursor -expand yes -side left

    pack $t$cursor_l $t$cursor_t $t$cursor_y -in $t.buttons -side left -expand yes 
   
    #grid columnconfig $t 0 -weight 1
    #grid rowconfig $t 2 -weight 2
   
    #grid $t.menubar -row 0 -column 0 -sticky ew -columnspan 3
    grid $t.buttons -row 1 -column 0 -sticky ew -columnspan 3
    grid $t.hscroll -row 4 -column 1 -sticky ew
    grid $t.brief   -row 5 -column 1 -sticky ew

    global $t.first_row
    set $t.first_row 2

    #must ensure that the packing is complete before calling consistency_reg
    #which interrogates the canvas width and height
    tkwait visibility $t.buttons

    CreateConsistencyRuler $io $t $r_win $t.hscroll

    #if user tries to destroy window 
    wm protocol $t WM_DELETE_WINDOW "ConsistencyStartShutdown $io $t"
    
    #must ensure that the packing is complete before calling 
    #consistency_display which interrogates the canvas width and height
    tkwait visibility $r_win

    set $t.cons_id [consistency_display -io $io -contigs $contig_list -frame $t -win_ruler $r_win -cursor_width [keylget gap_defs CONSISTENCY_DISPLAY.CURSOR_WIDTH] -cursor_fill [keylget gap_defs CONSISTENCY_DISPLAY.CURSOR_COLOUR]]

    #bind the configure actions to the toplevel
    bind $t <Any-Configure> "
    if {\[winfo toplevel %W\] == \"%W\"} {
	update idletasks
	resize_canvas -io $io -id [set $t.cons_id]
    }
    "

    keylset result cons_win $t
    keylset result cons_id [set $t.cons_id]

    #SetCanvasBindings $r_win $zoom_cmd
    SetConsistencyRulerBindings $io $t $r_win
    return $result
}

##############################################################################
#unique identifier to allow several consistency displays
proc next_consistency_display { } {
    global next_consist_display

    if {![info exists next_consist_display]} {
	set next_consist_display 0
	return $next_consist_display
    }
    incr next_consist_display
    return $next_consist_display
}

##############################################################################
#unique identifier to allow several windows of the same type to be displayed
#within the same consistency display
proc next_consistency_window { } {
    global next_consist_window

    if {![info exists next_consist_window]} {
	set next_consist_window 0
	return $next_consist_window
    }
    incr next_consist_window
    return $next_consist_window
}

##############################################################################
#find the id which must be of the format $win$id
proc get_consistency_window_id {win} {

    regexp {[0-9]+$} $win id
    return $id
}

##############################################################################
proc ConsistencyStartShutdown {io t } {
    global $t.cons_id

    if {[info exists $t.cons_id]} {
	result_quit -io $io -id [set $t.cons_id]
    }
}

##############################################################################
proc CreateConsistencyMenu {io t r_win} {
    global consistency_menu

    $t configure -menu $t.menubar
    menu $t.menubar
    create_menus $consistency_menu $t.menubar

    menu_state_set consistency_menu 12 $t.menubar
 
}

proc c_get_next_row {t} {
    global $t.row

    if {![info exists $t.row]} {
	set $t.row 2
	return [set $t.row]
    } 

    incr $t.row
    return [set $t.row]
}

#find the row of a consistency window
proc GetConsistencyRow { win } {
    array set info [grid info $win]
    return $info(-row)
}

#delete entire consistency display and unset global variables
proc DeleteConsistencyDisplay {t} {
    global config$t.ruler config$t.ticks
    global $t.contig_list

    #unset all global variables
    if [info exists config$t.ruler] {
	unset config$t.ruler     
    }
    if [info exists config$t.ticks] {
	unset config$t.ticks
    }
    #HACK - do I need this?
    #unset $t.contig_list

    destroy $t
}

##############################################################################
#display the consistency crosshair
proc AddConsistencyCrossHair {io id t canvas x} {
    global $t.cursor

    if {[set $t.cursor]} {
	draw_canvas_cursor_x -io $io -id $id -x [$canvas canvasx $x]
    } else {
	delete_canvas_cursor -io $io -id $id
    }
}

##############################################################################
#Creates and/or moves an editor cursor from within the template display
proc consistency_cursor_editor {io win c_id id x} {

    set cnum [consistency_contig -io $io -id $c_id -x [$win canvasx $x]]
    canvas_cursor_editor $io $cnum $id $x $win
}

##############################################################################
proc SetConsistencyRulerBindings { io t r_win } {
    global $t.cons_id

    bind $r_win <Any-Leave> "+delete_canvas_cursor -io $io -id [set $t.cons_id]"
    bind $r_win <Any-Motion> "AddConsistencyCrossHair $io [set $t.cons_id] $t $r_win %x"

    # Double button 1 and 2 to move or create an editor
    bind $r_win <<move-create>> "
	consistency_cursor_editor $io $r_win [set $t.cons_id] [set $t.cons_id] %x
    "
    bind $r_win <<use>> "
	consistency_cursor_editor $io $r_win [set $t.cons_id] [set $t.cons_id] %x
    "
    bind $r_win <<menu>> "set cnum \[consistency_contig \
		  -io $io \
		  -id [set $t.cons_id] \
		  -x \[$r_win canvasx %x\]\];
                  PopUpSingleContigMenu $io %W =\$cnum %X %Y"
}

##############################################################################
proc PopUpSingleContigMenu {io r_win c_name X Y} {
    global read_only

    if [winfo exists $r_win.m] {destroy $r_win.m}

    #otherwise, popup the contig menu
    create_popup $r_win.m "Contig Commands"

    set c_num [db_info get_contig_num $io $c_name]

    $r_win.m add command -label Information \
	    -command "PrintContigDetails $r_win $io $c_num"
    $r_win.m add command -label "Edit contig" \
	-command "edit_contig -io $io -contig [left_gel $io $c_num]"
    if {!$read_only} {
	$r_win.m add command -label "Complement contig" \
	    -command "complement_contig -io $io -contigs \"=$c_num\";
                      SetContigGlobals $io [left_gel $io $c_num]"

    }
    $r_win.m add command -label "List notes" \
	    -command "NoteSelector $io contig $c_name"
    tk_popup $r_win.m [expr $X-20] [expr $Y-10]
}

##############################################################################
#delete individual window from consistency display
proc DeleteConsistencyWindow {t win row_num} {
    global $t.row $t.first_row

    set new_height [expr [winfo height $t] - [winfo height $win]]

    grid_delete $t row $row_num 1

    #need to set the weight of the deleted row to 0
    grid rowconfigure $win $row_num -weight 0 -minsize 0

    incr $t.row -1

    #remove toplevel if no more slaves left
    if {[set $t.row] == [set $t.first_row]} {
	destroy $t
    }
    wm geometry $t [winfo width $t]x$new_height
}

#update geometry so the toplevel window can't grow larger than the screen
proc update_geom { t w } {

    wm geometry $t {}
    tkwait visibility $w

    set height [winfo height $t]
    if {$height > [lindex [wm maxsize $t] 1]} {
	set height [lindex [wm maxsize $t] 1]
    }
    wm geometry $t [winfo width $t]x$height
}
