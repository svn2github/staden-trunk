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
#delete stand alone restriction enzyme display
proc DeleteQualDisplay { f q_win} {
    global $f.qual_id

    unset $f.qual_id

    $q_win delete all
    destroy $q_win
    destroy $f
}


##############################################################################
#stand alone quality display
#user interface dialog box for quality plot
proc QualityPlotDisplay { io } {
    global gap_defs

    set f [keylget gap_defs QUALITY.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Quality plot"

    contig_id $f.id \
	    -state normal \
	    -io $io \
	    -range 1
    
    pack $f.id -side top
    
    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "Qual_OK_Pressed $io $f $f.id" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Quality}" \
	    -bd 2 \
	    -relief groove

    pack $f.ok_cancel -side top -fill both
   
}

##############################################################################
#stand alone quality display
proc Qual_OK_Pressed { io f id} {
    
    set contig [contig_id_gel $id]
    set lreg [contig_id_lreg $id]
    set rreg [contig_id_rreg $id]
    SetContigGlobals $io $contig $lreg $rreg

     # stop windows from hiding the plot
    destroy $f

    CreateQualityDisplay $io $contig $lreg $rreg

}

##############################################################################
proc AddQualityCursor {io id f q_win r_win x} {
    global $f.cursor

    if {[set $f.cursor]} {
	draw_canvas_cursor_x -io $io -id $id -x [$q_win canvasx $x]
    } else {
	DeleteQualityCursor $q_win $r_win
    }
}

##############################################################################
proc DeleteQualityCursor { q_win r_win} {

    $q_win delete cursor_x
    $r_win delete cursor_x

}

##############################################################################
proc SetQualityBindings { io f q_win r_win contig} {
    global $f.qual_id

    bind $q_win <Any-Leave> "DeleteQualityCursor $q_win $r_win"

    bind $q_win <Any-Motion> "AddQualityCursor $io [set $f.qual_id] $f $q_win $r_win %x"

    bind $q_win <<menu>> "PopUpQualityMenu $io $f $q_win [set $f.qual_id] %X %Y"

    # Double button 1 and 2 to move or create an editor
    bind $q_win <<move-create>> "
	canvas_cursor_editor $io [db_info get_contig_num $io $contig] \
		[set $f.qual_id] %x $q_win
    "
    bind $q_win <<use>> "
	canvas_cursor_editor $io [db_info get_contig_num $io $contig] \
		[set $f.qual_id] %x $q_win
    "
}

##############################################################################
proc SetQualityRulerBindings { io f q_win r_win contig} {
    global $f.qual_id

    bind $r_win <Any-Leave> "DeleteQualityCursor $q_win $r_win"

    bind $r_win <Any-Motion> "AddQualityCursor $io [set $f.qual_id] $f $q_win $r_win %x"

    bind $r_win <<menu>> "PopUpSingleContigMenu $io %W $contig %X %Y"

    
    # Double button 1 and 2 to move or create an editor
    bind $r_win <<move-create>> "
	canvas_cursor_editor $io [db_info get_contig_num $io $contig] \
		[set $f.qual_id] %x $q_win
    "
    bind $r_win <<use>> "
	canvas_cursor_editor $io [db_info get_contig_num $io $contig] \
		[set $f.qual_id] %x $q_win
    "
}

##############################################################################
proc PopUpQualityMenu {io f q_win qual_id X Y} {

    if [winfo exists $q_win.m] {destroy $q_win.m}

    #otherwise, popup the quality menu
    create_popup $q_win.m "Quality Commands"

    result_list_popup_single $io $qual_id [reg_get_ops -io $io -id $qual_id] $q_win.m

    tk_popup $q_win.m [expr $X-20] [expr $Y-10]
}

##############################################################################
#stand alone quality display
proc next_quality_display { } {
    global next_q_display

    if {![info exists next_q_display]} {
	set next_q_display 0
	return $next_q_display
    }
    incr next_q_display
    return $next_q_display

}

##############################################################################
#stand alone quality display
proc CreateQualityMenu {io f} {

    #set up a File menu
    menubutton $f.menubar.file -text "File" -menu $f.menubar.file.opts
    menu $f.menubar.file.opts
    $f.menubar.file.opts add command -label "Exit" \
	    -command "QualityStartShutdown $io $f"

    #HACK - to do
    # Help menu
    menubutton $f.menubar.help -text "Help" -menu [set m $f.menubar.help.opts]
    menu $m
    $m add command -label "Introduction" \
	-command "show_help gap4 {Quality}"
    $m add command -label "Examining the Plot" \
	-command "show_help gap4 {Quality-Examining}"

    #do the packing
    pack $f.menubar.file -side left
    pack $f.menubar.help -side right

}


##############################################################################
#stand alone quality display
proc CreateQualityDisplay {io contig lreg rreg } {
    global gap_defs

    set num_display [next_quality_display]

    set f [keylget gap_defs QUALITY.WIN]$num_display
    if {[xtoplevel $f] == ""} return
    wm protocol $f WM_DELETE_WINDOW "QualityStartShutdown $io $f"
    wm title $f "Quality plot: $contig #[db_info get_read_num $io $contig]"

    global $f.qual_id

    set q_win $f[keylget gap_defs QUALITY.WIN]
    set r_win $f[keylget gap_defs QUALITY.RULER.WIN]
    set scroll x
    set borderwidth [keylget gap_defs QUALITY.BORDERWIDTH]
    set width [keylget gap_defs QUALITY.PLOT_WIDTH]
    set height [keylget gap_defs QUALITY.PLOT_HEIGHT]

    set zoom_cmd [list "gap_zoom $io \[set $f.qual_id\] $scroll"]

    ##########################################################################
    #create quality display
    frame $f.q -bd $borderwidth -relief groove
    canvasbox $q_win -width $width -height $height \
	-bd 0 -highlightthickness 0 \
	-zoom_command $zoom_cmd 
    
    scrollbar $f.hscroll -orient horizontal -relief sunken -command \
	    "gc_scroll_x $io \[set $f.qual_id\]"

    ##########################################################################
    # Main Menu Bar
    frame $f.menubar -relief raised -borderwidth 2
    CreateQualityMenu $io $f

    ##########################################################################
    #button bar
    frame $f.buttons

    #zoom back button
    button $f.buttons.back -text "zoom out" -command "ZoomBackCanvas $io \[set $f.qual_id\]"
    #button $f.buttons.back -text "zoom out" -command 
    button $f.buttons.zoomin10 -text "+10%" \
	-command "ZoomInCanvas $q_win 0.05" 
    button $f.buttons.zoomin50 -text "+50%" \
	-command "ZoomInCanvas $q_win 0.1666" 
    
    #cursor checkbutton
    global $f.cursor
    checkbutton $f.buttons.cursor -text crosshairs -variable $f.cursor

    set cursor_t [keylget gap_defs QUALITY.CURSOR]
    label $f$cursor_t -bd 2 -relief sunken -width 6

    pack $f.buttons.zoomin10 $f.buttons.zoomin50 $f.buttons.back \
	 $f.buttons.cursor -expand yes -side left
    pack $f$cursor_t -in $f.buttons -side left -expand yes 

    ##########################################################################
    #ruler
    frame $f.r -bd $borderwidth -relief groove
    canvasbox $r_win -width $width \
	    -height [keylget gap_defs QUALITY.RULER.PLOT_HEIGHT] \
	    -xscrollcommand "$f.hscroll set" \
	    -bd 0 -highlightthickness 0 -zoom_command $zoom_cmd 

    pack $f.menubar -side top -fill x
    pack $f.buttons -side top -fill x
    pack $f.q -side top -fill both -expand yes
    pack $q_win -in $f.q -padx 5 -fill both -expand yes
    pack $f.hscroll -side bottom -fill x
    pack $f.r -side bottom -fill x
    pack $r_win -in $f.r -padx 5 -fill x -expand yes

    #must ensure that the packing is complete before calling quality_reg
    #which interrogates the canvas width and height
    tkwait visibility $q_win

    set c_name $contig
    set contig "{$contig $lreg $rreg}"
    #register quality plot and do first display
    set $f.qual_id \
	    [display_quality -io $io -contigs $contig -frame $f -window $q_win]

    #if user tries to destroy window 
    wm protocol $f WM_DELETE_WINDOW "QualityStartShutdown $io $f"

    #bind the configure actions to the toplevel
    bind $f <Any-Configure> "
	    if {\[winfo toplevel %W\] == \"%W\"} {
		update idletasks
		resize_canvas -io $io -id [set $f.qual_id]
	    }
	"

    SetCanvasBindings $q_win $zoom_cmd
    SetCanvasBindings $r_win $zoom_cmd

    SetQualityBindings $io $f $q_win $r_win $c_name
    SetQualityRulerBindings $io $f $q_win $r_win $c_name

    wm maxsize $f [winfo screenwidth $f] [winfo height $f]
    #set contig_num [db_info get_contig_num $io $contig]
    #set contig_length [c_length $io $contig_num]
}

##############################################################################
#start stand alone quality shutdown
proc QualityStartShutdown {io f } {
    global $f.qual_id

    if {[info exists $f.qual_id]} {
	result_delete -io $io -id [set $f.qual_id]
    }
}

##############################################################################
#template display quality plot functions
##############################################################################

##############################################################################
#bindings specific to the quality canvas
proc SetTemplateQualityBindings { io f canvas contig num} {
    global $f.template_id
    global $f.qual_t_id$num

    bind $canvas <Any-Leave> "delete_canvas_cursor -io $io -id [set $f.template_id]"
    bind $canvas <Any-Motion> "AddTemplateCrossHair $io [set $f.template_id] $f $canvas %x"

    bind $canvas <<menu>> "PopUpQualityMenu $io $f $canvas \[set $f.qual_t_id$num\] %X %Y"

    # Double button 1 and 2 to move or create an editor
    bind $canvas <<move-create>> "
	    canvas_cursor_editor $io [db_info get_contig_num $io $contig] \
	        [set $f.template_id] %x $canvas
    "
    bind $canvas <<use>> "
	    canvas_cursor_editor $io [db_info get_contig_num $io $contig] \
	        [set $f.template_id] %x $canvas
    "
}

##############################################################################
#Delete a quality plot canvas from template display
proc DeleteTemplateQualityPlot {f q_win} {
    global gap_defs
    global $f.num_q_win

    incr $f.num_q_win -1
    set q $f[keylget gap_defs TEMPLATE.QUALITY.WIN]
    set num [string range $q_win [string length $q] end]

    global config$f.quality$num $f.qual_t_id$num
    unset $f.qual_t_id$num
    unset config$f.quality$num

    #if deleting exitting gap4, $f may not exist
    if {![winfo exists $f]} {
	return
    }

    #set new_height [expr [winfo height $f] - [winfo height $q_win]]
    set new_height [expr [winfo height $f] - [winfo height $f.q$num]]

    $q_win delete all
    destroy $q_win
    destroy $f.q$num
    if {[set $f.num_q_win] == 0} {
	destroy $f.qual
	unset $f.num_q_win
    }

    wm geometry $f [winfo width $f]x$new_height
}

##############################################################################
proc CreateTemplateQualityPlot {io f q_win hscroll scale_label num contig} {
    global gap_defs
    global $f.template_id
    global $f.num_q_win

    if {![info exists $f.num_q_win]} {
	set $f.num_q_win 0
    }
    incr $f.num_q_win

    set borderwidth [keylget gap_defs TEMPLATE.BORDERWIDTH]

    #only scroll in x direction
    set scroll x
    set height [keylget gap_defs TEMPLATE.QUALITY.PLOT_HEIGHT]
    set width [winfo width $f[keylget gap_defs TEMPLATE.TEMPLATE.WIN]]

    #allow -height and -width to have affect
    wm geometry $f {}
   
    if {![winfo exists $f.qual]} {
	frame $f.qual
    }
    frame $f.q$num -bd $borderwidth -relief groove
    set zoom_cmd [list "gap_zoom $io \[set $f.template_id\] $scroll"]
    canvasbox $q_win -width $width -height $height \
	-borderwidth 0 \
	-highlightthickness 0 \
	-xscrollcommand "$hscroll set" \
	-zoom_command $zoom_cmd
    
    #set toplevel geometry when resized window
    set new_height [expr [winfo height $f] + [winfo reqheight $q_win] + \
	    2 * $borderwidth]

    
    grid $f.qual -row 4 -column 0 -sticky ew
    pack $f.q$num -in $f.qual -fill x -expand yes
    pack $q_win -in $f.q$num -padx 5 -fill x -expand yes
    update idletasks

    if {[winfo height $f] != 1} {
	wm geometry $f [winfo width $f]x$new_height
	#wm geometry $f [winfo width $f]x[winfo height $f]
    }
    #SetCanvasBindings $io [set $f.template_id] $q_win $scroll 
    SetCanvasBindings $q_win $zoom_cmd
    SetTemplateQualityBindings $io $f $q_win $contig $num

    $q_win delete all
}

##############################################################################
#template quality display
proc next_template_quality_display { } {
    global next_t_q_display

    if {![info exists next_t_q_display]} {
	set next_t_q_display 0
	return $next_t_q_display
    }
    incr next_t_q_display
    return $next_t_q_display
}

proc unset_next_template_quality_display { } {
    global next_t_q_display

    unset next_t_q_display
}

##############################################################################
#template display menu option to display quality plot
proc TemplateQualityPlot {io f template_id hscroll scale_label contig num} {
    global gap_defs
    global $f.contig_list
    global $f.qual_t_id$num
    global config$f.quality$num
    
    set q_win $f[keylget gap_defs TEMPLATE.QUALITY.WIN]$num
    if {[set config$f.quality$num]} {
	if {![template_win_free -io $io -id $template_id]} {
	    verror ERR_WARN template_display "no more windows available"
	    bell
	    set config$f.quality$num 0
	    return
	}

	#HACK to work round the horrible problem in zooming, turning ruler
	#off, scrolling and turning ruler on again.
	set scroll [lindex [$hscroll get] 0]
	
	if {![winfo exists $q_win]} {
	    CreateTemplateQualityPlot $io $f $q_win $hscroll $scale_label $num $contig
	}
	set $f.qual_t_id$num \
		[display_template_quality -io $io \
		-contigs $contig \
		-frame $f -win_quality $q_win -id $template_id]

	uplevel #0 eval [$hscroll cget -command] moveto $scroll

	if {[set $f.qual_t_id$num] == -1} {
	    DeleteTemplateQualityPlot $f $q_win
	}

    } else {
	delete_window -io $io -id $template_id -window $q_win
	result_delete -io $io -id [set $f.qual_t_id$num]
    }
}

