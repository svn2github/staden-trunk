#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

##############################################################################
#return the name of a stop codon
proc GetCodonName { f plot item io} {
    global $f.codon_id

    set index [GetREnzIndex $plot $item]
    set name [get_codon_name -io $io -id [set $f.codon_id] -codon $index\
	    -cnum 0]
    return $name
}

##############################################################################
proc HighlightCodon {f sc_win label io} {
    global restoreCmd
    global $sc_win.prev_item $sc_win.item_num

    set nearest [$sc_win find withtag current]
    #only do this code if the nearest item is different from the previous one
    if {$nearest != [set $sc_win.prev_item]} {
 
	#unhighlight object
    	if {[set $sc_win.item_num] != 0} {
	    eval $restoreCmd($sc_win,r_enz)
	} 

	if {$nearest != 0} {
	    InitialSettings $sc_win $sc_win $nearest r_enz
	    $sc_win itemconfig $nearest -fill white
	    $sc_win raise $nearest
	    set $sc_win.item_num $nearest
	    
	    $label configure -text "[GetCodonName $f $sc_win $nearest $io] \
		    position [GetREnzPos $sc_win $nearest]"
	}
    }
    #set previous item
    set $sc_win.prev_item $nearest

}
##############################################################################
proc AddCodonCursor {io id f sc_win r_win x} {
    global $f.cursor

    if {[set $f.cursor]} {
	draw_canvas_cursor_x -io $io -id $id -x [$sc_win canvasx $x]
    } else {
	DeleteCodonCursor $sc_win $r_win
    }
}

##############################################################################
proc DeleteCodonCursor { sc_win r_win} {

    $sc_win delete cursor_x
    $r_win delete cursor_x

}

##############################################################################
#bindings specific to the codon map
proc SetCodonBindings { f sc_win names r_win label dist io contig} {
    global $f.codon_id

    bind $sc_win <Any-Leave> "DeleteCodonCursor $sc_win $r_win"

    bind $sc_win <Any-Motion> "AddCodonCursor $io [set $f.codon_id] $f $sc_win $r_win %x"

    #any mouse motion - highlight nearest cut line
    $sc_win bind S <Any-Motion> [format {
	HighlightCodon %s %%W %s %s
    } [list $f] [list $label] [list $io]]
    $sc_win bind S <Shift-Motion> {;}

    #button-1 in plot canvas find the distance between 2 cut lines
    $sc_win bind S <<select>> "FindDistance $f %W $label $dist codon"

     # Double button 2 to move or create an editor
    bind $sc_win <<move-create>> "
	canvas_cursor_editor $io [db_info get_contig_num $io $contig] \
		[set $f.codon_id] %x $sc_win 
    "
    bind $sc_win <<use>> "
	canvas_cursor_editor $io [db_info get_contig_num $io $contig] \
		[set $f.codon_id] %x $sc_win 
    "
}

##############################################################################
proc SetCodonRulerBindings { io f sc_win r_win contig} {
    global $f.codon_id

    bind $r_win <Any-Leave> "DeleteCodonCursor $sc_win $r_win"
    bind $r_win <Any-Motion> "AddCodonCursor $io [set $f.codon_id] $f $sc_win $r_win %x"

    bind $r_win <<menu>> "PopUpSingleContigMenu $io %W $contig %X %Y"

    # Double button 2 to move or create an editor
    #bind $r_win <<move-create>> "CreateCanvasEditorCursor $io $f $sc_win [set $f.codon_id] %x $contig"
    
    bind $r_win <<move-create>> "
	canvas_cursor_editor $io [db_info get_contig_num $io $contig] \
		[set $f.codon_id] %x $r_win
    "
    bind $r_win <<use>> "
	canvas_cursor_editor $io [db_info get_contig_num $io $contig] \
		[set $f.codon_id] %x $r_win
    "
}

##############################################################################
#called from gaprc menu
proc StopCodonDisplay { io } {
    global gap_defs

    set f [keylget gap_defs STOP_CODON_DIAL.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Stop codon map"

    #########################################################################
    #contig id 
    contig_id $f.sc \
	    -io $io 

    #########################################################################
    #strand selection
    keylset st STRAND [keylget gap_defs STOP_CODON.STRAND]
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

    #########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "StopCodonOKPressed $io $f" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Stops}" \
	    -bd 2 \
	    -relief groove

    #final packing
    pack $f.sc -fill both
    pack $f.strand -fill both
    pack $f.ok_cancel -fill x
}

##############################################################################
proc StopCodonOKPressed { io d } {

    set gel_name [contig_id_gel $d.sc]
    set lreg [contig_id_lreg $d.sc]
    set rreg [contig_id_rreg $d.sc]
    SetContigGlobals $io $gel_name $lreg $rreg
    set contig_num [db_info get_contig_num $io [contig_id_gel $d.sc]]

    #strand = 1 for plus strand
    #strand = 2 for complementary strand
    #strand = 3 for both strands
    set strand [radiolist_get $d.strand]

    destroy $d
   
    CreateStopCodonDisplay $io $contig_num $lreg $rreg $strand

}
##############################################################################
#stand alone stop codon display
proc next_stop_codon_display { } {
    global next_c_display

    if {![info exists next_c_display]} {
	set next_c_display 0
	return $next_c_display
    }
    incr next_c_display
    return $next_c_display
}

##############################################################################
proc CreateStopCodonDisplay {io contig_num lreg rreg strand } {
    global CurContig
    global gap_defs

    #generate a new display
    set num_display [next_stop_codon_display]
    set f [keylget gap_defs STOP_CODON.WIN]$num_display
    xtoplevel $f

    wm protocol $f WM_DELETE_WINDOW "CodonStartShutdown $io $f"
    wm title $f "Stop codons: $CurContig #[db_info get_read_num $io $CurContig]"

    global $f.codon_id
    global $f.prev_strand
    
    #initialise the strand 
    set $f.prev_strand $strand
    set sc_win $f[keylget gap_defs CODON.WIN]
    set r_win $f[keylget gap_defs CODON.RULER.WIN]
    set scroll x
    set borderwidth [keylget gap_defs CODON.BORDERWIDTH]
    set width [keylget gap_defs CODON.PLOT_WIDTH]
    set height [keylget gap_defs CODON.PLOT_HEIGHT]

    global $sc_win.prev_item $sc_win.item_num
    set $sc_win.prev_item 0
    set $sc_win.item_num 0

    set tick_ht [keylget gap_defs CODON.TICK_HEIGHT]
    set offset $tick_ht
    
    #create stop codon canvas
    frame $f.m
    frame $f.m.top 
    frame $f.m.top.left
    frame $f.m.top.right
    frame $f.m.bot
    frame $f.m.bot.right

    #if both strands are to be displayed, double canvas height
    if {$strand == 3} {
	set height [expr $height * 9 / 5]
    }

    set zoom_cmd [list "gap_zoom $io \[set $f.codon_id\] $scroll"]
    canvasbox $sc_win -width $width -height $height \
	    -bd 0 -highlightthickness 0\
	    -xscrollcommand "$f.hscroll set" \
	    -closeenough 4 -zoom_command $zoom_cmd
    
    scrollbar $f.hscroll -orient horizontal -relief sunken -command \
	    "gc_scroll_x $io \[set $f.codon_id\]"

    ##########################################################################
    # Main Menu Bar
    frame $f.menubar -relief raised -borderwidth 2
    CreateCodonMenu $io $f $sc_win $contig_num $sc_win $strand

    ##########################################################################
    #create canvas of frame names
    canvasbox $f.names -width [keylget gap_defs CODON.NAMES.PLOT_WIDTH]\
	-height $height -bd 0 -highlightthickness 0

    ##########################################################################
    #button bar
    frame $f.buttons
    #zoom back button
    button $f.buttons.back -text "zoom out" -command "ZoomBackCanvas $io \[set $f.codon_id\]"
    button $f.buttons.zoomin10 -text "+10%" \
	-command "ZoomInCanvas $sc_win 0.05" 
    button $f.buttons.zoomin50 -text "+50%" \
	-command "ZoomInCanvas $sc_win 0.1666" 
    
    #refresh button
    button $f.buttons.refresh -text Refresh -command "RefreshCodon $io $f \
	    $sc_win $contig_num \[set $f.Strand\] 1"

    #cursor checkbutton
    global $f.cursor
    checkbutton $f.buttons.cursor -text crosshairs -variable $f.cursor

    #cursor position label
    set cursor_t [keylget gap_defs CODON.CURSOR]
    label $f$cursor_t -bd 2 -relief sunken -width 6

    #stop codon distance label
    label $f.buttons.dist -relief sunken -bd 2 -width 6

    pack $f.buttons.zoomin10 $f.buttons.zoomin50 $f.buttons.back \
	 $f.buttons.refresh $f.buttons.cursor -expand yes -side left
    pack $f$cursor_t -in $f.buttons -side left -expand yes 
    pack $f.buttons.dist -side left -expand yes 

    ##########################################################################
    #create ruler canvas
    set zoom_cmd [list "gap_zoom $io \[set $f.codon_id\] $scroll"]
    canvasbox $r_win -width $width \
	-height [keylget gap_defs CODON.RULER.PLOT_HEIGHT] \
	-xscrollcommand "$f.hscroll set" \
	-bd 0 -highlightthickness 0 \
	-zoom_command $zoom_cmd

    #filler frame for ruler
    frame $f.m.bot.left -width [winfo reqwidth $f.names] \
	    -height [winfo reqheight $r_win]

    ##########################################################################
    label $f.brief
    pack $f.brief -side bottom -fill both

    pack $f.menubar -side top -fill x
    pack $f.buttons -side top -fill x
    pack $f.names -in $f.m.top.left -side left -fill both -expand yes
    pack $sc_win -in $f.m.top.right -side left -fill both -expand yes
    pack $f.m -side left -fill both -expand yes
    pack $f.m.top -side top -fill both -expand yes
    pack $f.m.top.left -side left -fill both
    pack $f.m.top.right -side right -fill both -expand yes
    pack $f.m.bot -side bottom -fill x -expand yes
    pack $f.m.bot.left -side left -fill x
    pack $f.m.bot.right -side right -fill x -expand yes
    pack $f.hscroll -in $f.m.bot.right -side bottom -fill x
    pack $r_win -in $f.m.bot.right -side bottom -fill x -expand yes

    #must ensure that the packing is complete before calling stopcodon_reg
    #which interrogates the canvas width and height
    tkwait visibility $sc_win
    
    #do plotting
    set $f.codon_id [plot_stop_codons -strand $strand\
	    -tick_height $tick_ht -win_names $f.names -window $sc_win \
	    -yoffset $offset -io $io -contigs "{=$contig_num $lreg $rreg}" \
	    -frame $f]

    #bind the configure actions to the toplevel
    bind $f <Any-Configure> "
	    if {\[winfo toplevel %W\] == \"%W\"} {
		update idletasks
		resize_canvas -io $io -id [set $f.codon_id]
	    }
	"
    #set bindings
    set contig_name [left_gel $io $contig_num]
    SetCanvasBindings $sc_win $zoom_cmd
    SetCanvasBindings $r_win $zoom_cmd
    SetCodonBindings $f $sc_win $f.names $r_win $f.brief $f.buttons.dist $io $contig_name
    SetCodonRulerBindings $io $f $sc_win $r_win $contig_name

    update idletasks
    wm minsize $f [winfo reqwidth $f.buttons] [winfo height $f]
    wm maxsize $f [winfo screenwidth $f] [winfo height $f]

}

##############################################################################
#menu display
proc CreateCodonMenu { io f sc_win contig_num plot strand} {
    global $f.Strand
    global gap_defs

    set $f.Strand $strand
    #set up a File menu
    menubutton $f.menubar.file -text "File" -menu $f.menubar.file.opts
    menu $f.menubar.file.opts
    $f.menubar.file.opts add command -label "Exit" \
	    -command "CodonStartShutdown $io $f"

    #set up View menu
    menubutton $f.menubar.view -text "View" -menu $f.menubar.view.opts
    menu $f.menubar.view.opts

    keylset st STRAND [keylget gap_defs STOP_CODON.STRAND]
    set b1 [keylget st STRAND.BUTTON.1]
    set b2 [keylget st STRAND.BUTTON.2]
    set b3 [keylget st STRAND.BUTTON.3]

    $f.menubar.view.opts add radiobutton -label $b1 \
	    -command "RefreshCodon $io $f $sc_win $contig_num 1 0" -value 1 \
	    -variable $f.Strand
    $f.menubar.view.opts add radiobutton -label $b2 \
	    -command "RefreshCodon $io $f $sc_win $contig_num 2 0" -value 2 \
	    -variable $f.Strand
    $f.menubar.view.opts add radiobutton -label $b3 \
	    -command "RefreshCodon $io $f $sc_win $contig_num 3 0" -value 3 \
	    -variable $f.Strand

    # Help menu
    menubutton $f.menubar.help -text "Help" -menu [set m $f.menubar.help.opts]
    menu $m
    $m add command -label "Introduction" \
	-command "show_help gap4 {Stops}"
    $m add command -label "Examining the Plot" \
	-command "show_help gap4 {Stops-Examining}"
    $m add command -label "Updating the Plot" \
	-command "show_help gap4 {Stops-Updating}"

    #do the packing
    pack $f.menubar.file $f.menubar.view -side left
    pack $f.menubar.help -side right

}

##############################################################################
proc RefreshCodon {io f sc_win contig_num strand update} {
    global NGRec
    global gap_defs
    global $f.codon_id
    global $f.prev_strand

    #save current window height minus the stop codon canvas height
    set old_height [expr [winfo height $f] - [winfo height $sc_win]]
    set sc_height [winfo height $sc_win]

    #change plot height if necessary
    if {$strand == 3 && [set $f.prev_strand] != 3} {
	#need to double display
	set height [expr $sc_height * 9 / 5]
    } elseif {$strand != 3 && [set $f.prev_strand] == 3} {
	#need to halve display
	set height [expr $sc_height * 5 / 9]
    } else {
	#stay the same
	set height $sc_height
    }
    set new_height [expr $old_height + $height]

    # Shrink windows to min height to let the "wm geometry" resize them.
    $sc_win configure -height 1
    [winfo parent $sc_win].names configure -height 1

    #only way I could "disable" minsize and maxsize - very HACKy
    wm minsize $f 1 1
    wm maxsize $f [winfo screenwidth $f] [winfo screenheight $f]

    #resize toplevel and canvas
    wm geometry $f [winfo width $f]x$new_height
    $sc_win configure -height $height

    #set min and maxsize back again
    wm minsize $f [winfo reqwidth $f.buttons] $new_height
    wm maxsize $f [winfo screenwidth $f] $new_height

    refresh_codon_map -id [set $f.codon_id] -io $io -cnum $contig_num \
	    -strand $strand -update $update
    set $f.prev_strand $strand
}

##############################################################################
#executed when the exit command is chosen from the File menu of stop codon
proc CodonStartShutdown { io f } {
    global $f.codon_id

    if {[info exists $f.codon_id]} {
	result_delete -io $io -id [set $f.codon_id]
    }
}

proc DeleteCodonPlot {f sc_win} {
    global $f.codon_id
    global $f.cursor
    
    unset $f.codon_id
    unset $f.cursor

    $sc_win delete all
    destroy $sc_win
    destroy $f

}
