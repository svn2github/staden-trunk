#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

proc renz_map_configure {w args} {
    global $w.zoom_cmd $w.zoomback_cmd $w.scrollbar_x_cmd $w.scrollbar_y_cmd
    global $w.cursor_cmd $w.renz_name_cmd $w.renz_info_cmd $w.invoke_cmd
    global $w.tick_ht $w.text_offset $w.text_fill $w.config_cmd

    set in_arg 0
    set arglist ""
    set height ""
    set width ""
    set ruler_height ""
    set names_width ""
    set selectbackground ""
    set tick_ht ""

    set $w.zoom_cmd ""
    set $w.zoomback_cmd ""
    set $w.scrollbar_x_cmd ""
    set $w.scrollbar_y_cmd ""
    set $w.cursor_cmd ""
    set $w.renz_name_cmd ""
    set $w.renz_info_cmd ""
    set $w.invoke_cmd ""
    set $w.tick_ht ""
    set $w.text_offset ""
    set $w.text_fill ""
    set $w.config_cmd ""

    # Process command line args
    foreach i $args {
	if {$in_arg} {
	     if {$option == "-zoom_command"} {
		set $w.zoom_cmd "{$i}"
	     } elseif {$option == "-zoomback_command"} {
		 set $w.zoomback_cmd "{$i}"
	     } elseif {$option == "-scrollbar_x_cmd"} {
		 set $w.scrollbar_x_cmd "{$i}"
	     } elseif {$option == "-scrollbar_y_cmd"} {
		 set $w.scrollbar_y_cmd "{$i}"
	     } elseif {$option == "-cursor_cmd"} {
		 set $w.cursor_cmd "{$i}"
	     } elseif {$option == "-renz_name_cmd"} {
		 set $w.renz_name_cmd "{$i}"
	     } elseif {$option == "-renz_info_cmd"} {
		 set $w.renz_info_cmd "{$i}"
	     } elseif {$option == "-width"} {
		set width "-width $i"
	     } elseif {$option == "-height"} {
		set height "-height $i"
	     } elseif {$option == "-ruler_height"} {
		set ruler_height "-height $i"
	     } elseif {$option == "-names_width"} {
		set names_width "-width $i"
	     } elseif {$option == "-tick_ht"} {
		 set $w.tick_ht "{$i}"
	     } elseif {$option == "-selectbackground"} {
		 set selectbackground "-selectbackground $i"
	     } elseif {$option == "-text_offset"} {
		 set $w.text_offset "{$i}"
	     } elseif {$option == "-text_fill"} {
		 set $w.text_fill "{$i}"
	     } elseif {$option == "-invoke_cmd"} {
		 set $w.invoke_cmd "{$i}"
	     } elseif {$option == "-config_cmd"} {
		 set $w.config_cmd "{$i}"
	     } else {
		lappend arglist $option $i
	     }
	     set in_arg 0
	} else {
	     set option $i
	     set in_arg 1
	}
    }
    eval $w.renz configure $arglist

    if {$height != ""} {
	eval $w.renz config $height
	eval $w.names config $height
    } 

    if {$ruler_height != ""} {
	eval $w.ruler config $ruler_height
    } 
    if {$width != ""} {
	eval $w.renz config $width
	eval $w.ruler config $width
    }

    if {$names_width != ""} {
	eval $w.names config $names_width
    }

    if {$selectbackground != ""} {
	eval $w.names config $selectbackground
    }

    if {[set $w.scrollbar_x_cmd] != ""} {
	eval $w.hscroll config -command [set $w.scrollbar_x_cmd]
    }
    if {[set $w.scrollbar_y_cmd] != ""} {
	eval $w.vscroll config -command [set $w.scrollbar_y_cmd]
    }
    if {[set $w.zoom_cmd] != ""} {
	set zoom_cmd "[set $w.zoom_cmd] -direction b"
	eval canvasbox_configure $w.renz -zoom_command [list $zoom_cmd]
	set zoom_cmd "[set $w.zoom_cmd] -direction x"
	eval canvasbox_configure $w.ruler -zoom_command [list $zoom_cmd]
    }
    if {[set $w.zoomback_cmd] != ""} {
	eval $w.buttons.back config -command [set $w.zoomback_cmd]
    }
}

#create restriction enzyme canvas
proc renz_map {w args} {
    global tk_utils_defs
    global $w.cursor_cmd $w.renz_name_cmd $w.renz_info_cmd $w.invoke_cmd
    global $w.config_cmd
    xtoplevel $w

    canvasbox $w.renz -yscrollcommand "$w.vscroll set" \
	    -bd 0 -highlightthickness 0

    scrollbar $w.hscroll -orient horizontal -relief sunken
    scrollbar $w.vscroll -orient vertical -relief sunken 

    frame $w.menubar -relief raised -borderwidth 2
 
    frame $w.buttons
    #zoom back button
    button $w.buttons.back -text "zoom out" 
    button $w.buttons.zoomin10 -text "+10%" \
	-command "ZoomInCanvas $w.renz 0.05" 
    button $w.buttons.zoomin50 -text "+50%" \
	-command "ZoomInCanvas $w.renz 0.1666" 
  
    #cursor checkbutton
    global $w.cursor
    checkbutton $w.buttons.cursor -text crosshairs -variable $w.cursor

    #cursor position label
    set cursor_t [keylget tk_utils_defs R_ENZ.CURSOR]
    label $w$cursor_t -bd 2 -relief sunken -width 6

    #restriction enzyme distance label
    label $w.buttons.dist -relief sunken -bd 2 -width 6

    pack $w.buttons.zoomin10 $w.buttons.zoomin50 $w.buttons.back \
	 $w.buttons.cursor -expand yes -side left
    pack $w$cursor_t -in $w.buttons -side left -expand yes 
    pack $w.buttons.dist -side left -expand yes 

    ##########################################################################
    #create canvas of names of restriction enzymes
    canvasbox $w.names \
	    -selectbackground [keylget tk_utils_defs R_ENZ.SELECT_COLOUR] \
	    -bd 0 -highlightthickness 0

    ##########################################################################
    #create ruler canvas
    canvasbox $w.ruler -xscrollcommand "$w.hscroll set" \
	-bd 0 -highlightthickness 0

    ##########################################################################
    label $w.brief 
  
    grid columnconfig $w 1 -weight 1
    grid rowconfig $w 2 -weight 1

    grid $w.menubar -row 0 -column 0 -sticky ew -columnspan 3
    grid $w.buttons -row 1 -column 0 -sticky ew -columnspan 3
    grid $w.names -row 2 -column 0 -sticky ns 
    grid $w.renz -row 2 -column 1 -sticky nsew
    grid $w.ruler -row 3 -column 1 -sticky ew
    grid $w.hscroll -row 4 -column 1 -sticky ew
    grid $w.brief -row 5 -column 1 -sticky ew
    grid $w.vscroll -row 2 -column 2 -sticky ns
    
    eval renz_map_configure $w $args

    global $w.renz.prev_item $w.renz.item_num
    set $w.renz.prev_item 0
    set $w.renz.item_num 0
    set $w.sel_enz ""

    SetREnzBindings $w $w.renz $w.ruler $w.names $w.brief $w.buttons.dist \
	[set $w.cursor_cmd] [set $w.renz_name_cmd] [set $w.renz_info_cmd] \
       [set $w.invoke_cmd] [set $w.config_cmd]

    SetREnzRulerBindings $w $w.renz $w.ruler [set $w.cursor_cmd] [set $w.invoke_cmd] 
}


##############################################################################
proc next_renz_display { } {
    global next_r_display

    if {![info exists next_r_display]} {
	set next_r_display 0
	return $next_r_display
    }
    incr next_r_display
    return $next_r_display
}

##############################################################################
proc renz_map_plot {w cmd} {
    global $w.scrollbar_x_cmd $w.scrollbar_y_cmd
    global $w.tick_ht $w.text_offset $w.text_fill

    return [eval $cmd -frame $w -win_names $w.names -window $w.renz \
		-win_ruler $w.ruler \
		-text_offset [set $w.text_offset] \
		-text_fill [set $w.text_fill]  \
		-tick_height [set $w.tick_ht] -yoffset [set $w.tick_ht]]
}

##############################################################################
#called from C
#delete stand alone restriction enzyme display
proc DeleteREnzPlot { w re_win} {
    global REnzyme
    global $w.renz_id
    global $w.cursor
    global $w.pos
    global $w.plot

    if {[info exists REnzyme]} {
	unset REnzyme
    }
    unset $w.renz_id
    unset $w.cursor

    if {[info exists $w.pos]} {
	unset $w.pos
    } 
    if {[info exists $w.plot]} {
	unset $w.plot
    }

    $re_win delete all
    destroy $re_win
    destroy $w
}

#create selectable rectangle
proc CreateSelRect {names item_num } {
    global tk_utils_defs

    set colour [keylget tk_utils_defs R_ENZ.SELECT_COLOUR]
    set box [$names bbox re_$item_num]
    set x1 [lindex $box 0]
    set y1 [lindex $box 1]
    set x2 [winfo width $names]
    set y2 [lindex $box 3]
    $names lower [$names create rectangle $x1 $y1 $x2 $y2 -fill $colour \
	    -outline "" -tag r_$item_num]
}

#after a redraw, add selections back again
#called from C
proc ReSelectRect {f names} {
    global $f.sel_enz

    if {[info exists $f.sel_enz]} {
	foreach tag [set $f.sel_enz] {
	    CreateSelRect $names [string trim $tag re_]
	}
    }
}


##############################################################################
#
proc AddREnzCursor {f re_win r_win x cursor_cmd} {
    global $f.cursor

    if {[set $f.cursor]} {
	uplevel #0 eval $cursor_cmd -x [$re_win canvasx $x]
    } else {
	DeleteREnzCursor $re_win $r_win
    }
}
##############################################################################
#
proc DeleteREnzCursor {re_win r_win} {

    $re_win delete cursor_x
    $r_win delete cursor_x
}

##############################################################################
#invoked by moving the mouse over a restriction enzyme cut line
proc HighlightREnz {f plot label name_cmd} {
    global restoreCmd
    global $plot.prev_item $plot.item_num

    set nearest [$plot find withtag current]

    #only do this code if the nearest item is different from the previous one

    if {$nearest != [set $plot.prev_item]} {
 
	#unhighlight object
    	if {[set $plot.item_num] != 0} {
	    eval $restoreCmd($plot,r_enz)
	} 

	if {$nearest != 0} {
	    InitialSettings $plot $plot $nearest r_enz
	    $plot itemconfig $nearest -fill white
	    $plot raise $nearest
	    set $plot.item_num $nearest
	    $label configure -text "[GetREnzName $f $plot $nearest $name_cmd] cut position [GetREnzPos $plot $nearest]"
	}
    }
    #set previous item
    set $plot.prev_item $nearest

}

proc RenzInvokeCmd {re_win x invoke_cmd} {
    uplevel #0 eval $invoke_cmd $x $re_win
}

##############################################################################
#bindings specific to the stand alone restriction enzyme 
proc SetREnzBindings {f re_win r_win names_win label dist cursor_cmd renz_name_cmd renz_info_cmd invoke_cmd config_cmd} {
    global REnzY1
    global $f.renz_id $names_win.Move

    bind $re_win <Any-Leave> "DeleteREnzCursor $re_win $r_win"

    bind $re_win <Any-Motion> "AddREnzCursor $f $re_win $r_win %x $cursor_cmd"

    #any mouse motion - highlight nearest cut line
    $re_win bind S <Any-Motion> \
	"HighlightREnz $f $re_win $label $renz_name_cmd"
    $re_win bind S <Shift-Motion> {;}

    #button-1 in plot canvas find the distance between 2 cut lines
    $re_win bind S <<select>> "FindDistance $f %W $label $dist cut"

    #button-1 in names canvas - select a name
    bind $names_win <<select>> "PickREnz $names_win $f"

    #button-2 in names canvas - move a name and assoc cut lines
    set $names_win.Move 0
    $names_win bind S <<move>> \
	    "if {\[SelectREnzPlot $f $names_win $re_win %y\] == -1} {set $names_win.Move 0} else {set $names_win.Move 1} "
    $names_win bind S <<move-drag>> \
	    "if {\[set $names_win.Move\] == 1} {
		MoveREnzPlot $names_win $re_win %y
	    }"
    $names_win bind S <<move-release>> \
	    "if {\[set $names_win.Move\] == 1} {
		DropREnzPlot $f $names_win $re_win %y
		set $names_win.Move 0
	    }"

    #button-2 in plot canvas - disable binding which scrolls the entire canvas
    #bind $re_win <2> {}
    #bind $re_win <B2-Motion> {}

    #auto scrolling
    bind $names_win <<move-autoscroll>> "set $f.auto_scroll 1; REnzAutoScroll $names_win $re_win %y $f y"

    bind $names_win <<stop-autoscroll>> "+REnzStopScroll $f"

    #button-3 in plot canvas - invoke a popup menu
    $re_win bind S <<menu>> "PopUpREnzMenu $f $re_win %X %Y $renz_info_cmd"

    #Double button 1 or 2 to move or create an editor
    bind $re_win <<move-create>> "RenzInvokeCmd $re_win %x $invoke_cmd"
    bind $re_win <<use>> "RenzInvokeCmd $re_win %x $invoke_cmd"

    #must bind names_win here as it is the last window to be updated by the
    #Configure callback. If use re_win, the scrollregion is calculated
    #incorrectly for the names win as it hasn't been updated
    
    #07.03.00
    #changed this back again because it introduces another bug because
    #increasing the re_win in x only doesn't cause the binding to be activated
    #solution for now: change back again but alter nip so that the canvas
    #canvas can't be resized in y. This should be fixed at a later date to
    #show more enzymes rather than just increase their size

    bind $re_win <Any-Configure> "uplevel #0 eval $config_cmd"
    #bind $names_win <Any-Configure> "uplevel #0 eval $config_cmd"
}

##############################################################################
proc SetREnzRulerBindings {f re_win r_win cursor_cmd invoke_cmd} {
    global $f.renz_id

    bind $r_win <Any-Leave> "DeleteREnzCursor $re_win $r_win"
    bind $r_win <Any-Motion> "AddREnzCursor $f $re_win $r_win %x $cursor_cmd"

    # Double button 2 to move or create an editor
    bind $r_win <<move-create>> "RenzInvokeCmd $re_win %x $invoke_cmd"
}

##############################################################################
#return the index of a restriction enzyme
proc GetREnzIndex {plot item } {

    foreach tag [$plot gettags $item] {	
	if {[string compare [string range $tag 0 2] re_] == 0} {
	    return [string trim $tag re_]
	} 
    }
    #no tags
    return -1
}

##############################################################################
#return the name of a restriction enzyme
proc GetREnzName { f plot item name_cmd} {

    set index [GetREnzIndex $plot $item]

    set name [uplevel #0 eval $name_cmd -enzyme $index]
    return $name
}

##############################################################################
#return the restriction enzyme cut position
proc GetREnzPos { plot id } {

    foreach tag [$plot gettags $id] {
	if {[string compare [string range $tag 0 3] pos_] == 0} {
	    return [string trim $tag pos_]
	}
    }
}

##############################################################################
#return the entire restriction enzyme tag
proc GetREnzTag { c i} {

    foreach tag [$c gettags $i] {
	if {[string compare [string range $tag 0 2] re_] ==0} {
	    return $tag
	}
    }
}

##############################################################################
#find the distance between two restriction enzyme cuts
proc FindDistance { f re_win label dist text} {
    global $f.pos
    global $f.plot

    if {[info exists $f.plot]} {
	if {[set $f.plot] != $re_win} {
	    return
	}
    }
    set x [GetREnzPos $re_win current]

    if { [info exists $f.pos]} {
	$dist configure -text [expr abs($x - [set $f.pos])]
	$label configure -text "Distance: [expr abs($x - [set $f.pos])]"
	unset $f.pos
	unset $f.plot
	bell
    } else {
	set $f.pos $x
	set $f.plot $re_win
	$label configure -text "Select another $text"
	bell
    }
}

##############################################################################
#highlight selected enzyme name with a selectable rectangle
proc PickREnz { names f} {
    global $f.sel_enz

    #check the current selection is valid
    if {[$names find withtag current] == ""} {
	return
    }
    #check tag is selectable (ie exclude the rectangles)
    if {[lsearch [$names gettag [$names find withtag current]] S] == -1} {
	return
    }
    set tag [GetREnzTag $names current]
    set item_num [GetREnzIndex $names current]

    #remove rectangle if it already exists - ie toggle off
    if {[info exists $f.sel_enz]} {
	set pos [lsearch [set $f.sel_enz] $tag]
    } else {
	set pos -1
    }
    if { $pos > -1} {
	$names delete r_$item_num
	set $f.sel_enz [lreplace [set $f.sel_enz] $pos $pos]
	return
    }
    
    CreateSelRect $names $item_num
    lappend $f.sel_enz $tag

}

##############################################################################
#stand alone restriction enzyme display
#selected enzyme name and plot to move - set up extents of movement and first
#y position
proc SelectREnzPlot {f names r_enz y} {
    global REnzyme
    global REnzY1
    global NGConst

    #puts "start SELECTRENZPLOT"
    #if no cuts, don't allow user to move
    if {[$r_enz coords S] == ""} {
	return -1
    }

    set REnzyme($f,smallest) [winfo height $names]
    set REnzyme($f,largest) 0

    set REnzY1 [$names canvasy $y]
    set REnzyme($f,orig_pos) [lindex [$names coord current] 1]

    foreach i [$names find withtag S] {
	set item_pos [lindex [$names coords $i] 1]
	if {$item_pos < $REnzyme($f,smallest)} {
	    set REnzyme($f,smallest) $item_pos
	}
	#find largest coord
	if {$item_pos > $REnzyme($f,largest)} {
	    set REnzyme($f,largest) $item_pos
	}
    }
    return 0
}

##############################################################################
#move a restriction enzyme plot in the stand alone display
#must move the enzyme name, its associated 'selecting' rectangle and its plot
proc MoveREnzPlot { names r_enz y} {
    global REnzY1

    #puts "start MOVERENZPLOT"

    #set min_x $NGConst(ymargin)
    #set max_x [expr $NGConst($names,win_len) - $NGConst(ymargin)]
    #set min_y $NGConst(xmargin)
    #set max_y [expr $NGConst($names,win_ht) - $NGConst(xmargin)]
    
    set min_x 1
    set max_x [winfo width $names]
    set min_y 1
    set max_y [winfo height $names]

    if {($y > $min_y) && ($y < $max_y)} {
	set y [$names canvasy $y]
	set dx 0
	set dy [expr $y - $REnzY1]

	set tag [GetREnzTag $names current]
	set item_num [GetREnzIndex $names current]
	#move both name and cuts
	$names move $tag $dx $dy
	$names move r_$item_num $dx $dy
	$r_enz move $tag $dx $dy
	set REnzY1 $y
    }
}

##############################################################################
#drop a restriction enzyme plot in the stand alone display
#must move the enzyme name, its associated 'selecting' rectangle and its plot
proc DropREnzPlot {f names r_enz y} {
    global REnzyme
    global gap_defs

    #puts "start DROPRENZPLOT [$r_enz coords S]"
    REnzStopScroll $f

    set current [$names find withtag current]; #moved item
    set cur_pos [lindex [$names coord $current] 1]; #y coord of moved item

    #set offset to be tick height
    set offset [expr [lindex [$r_enz coords S] 3] - \
	    [lindex [$r_enz coords S] 1]]

    #set text_offset to be the distance to the first text item
    set text_offset $REnzyme($f,smallest)

    #fudge factor necessary due to precision errors
    set multiplier [expr int((($cur_pos - $text_offset)/$offset) + 0.0001)]
    set final_pos [expr $text_offset + ($multiplier * $offset)]

    if {$final_pos == $REnzyme($f,orig_pos)} {
	set tag [GetREnzTag $names $current]; #corresponding cut plot tag
	set item_num [GetREnzIndex $names $current]
	$names move $current 0 [expr $final_pos - $cur_pos]
	$names move r_$item_num 0 [expr $final_pos - $cur_pos]
	$r_enz move $tag 0 [expr $final_pos - $cur_pos]
	return
    }

    if {$cur_pos >= $REnzyme($f,orig_pos)} {
	foreach i [$names find withtag S] {
	    set item_pos [lindex [$names coords $i] 1]
	    if {($item_pos > $REnzyme($f,orig_pos)) && \
		    ($item_pos <= $cur_pos) && ($i != $current)} {
		set tag [GetREnzTag $names $i]
		set item_num [GetREnzIndex $names $i]
		$names move $i 0 [expr $offset * -1]
		$names move r_$item_num 0 [expr $offset * -1]
		$r_enz move $tag 0 [expr $offset * -1]
	    }
	}
	if {$final_pos > $REnzyme($f,largest)} {
	    set final_pos $REnzyme($f,largest)
	}
	set tag [GetREnzTag $names $current]; #corresponding cut plot tag
	set item_num [GetREnzIndex $names $current]
	$names move $current 0 [expr $final_pos - $cur_pos]
	$names move r_$item_num 0 [expr $final_pos - $cur_pos]
	$r_enz move $tag 0 [expr $final_pos - $cur_pos]

    } elseif {$cur_pos < $REnzyme($f,orig_pos)} {
	#must round up here
	set multiplier [expr ceil(($cur_pos - $text_offset)/$offset)]
	set final_pos [expr $text_offset + ($multiplier * $offset)]
	foreach i [$names find withtag S] {
	    set item_pos [lindex [$names coords $i] 1]
	    if {($item_pos < $REnzyme($f,orig_pos)) && \
		    ($item_pos >= $cur_pos) && ($i != $current)} {
		set tag [GetREnzTag $names $i]
		set item_num [GetREnzIndex $names $i]
		$names move $i 0 $offset
		$names move r_$item_num 0 $offset
		$r_enz move $tag 0 $offset
	    }
	}
	if {$final_pos < $REnzyme($f,smallest)} {
	    set final_pos $REnzyme($f,smallest)
	}
	set tag [GetREnzTag $names $current]; #corresponding cut plot tag
	set item_num [GetREnzIndex $names $current]; #corresponding rectangle
	$names move $current 0 [expr $final_pos - $cur_pos]
	$names move r_$item_num 0 [expr $final_pos - $cur_pos]
	$r_enz move $tag 0 [expr $final_pos - $cur_pos]
    }
}

##############################################################################
#stand alone restriction enzyme display
#auto scrolls whilst the mouse pointer is outside the canvas
proc REnzDoScroll { f names_win re_win direction scroll } {
    global CanvasConst
    global $f.auto_scroll

    if {[set $f.auto_scroll] && [winfo exists $names_win]} {
	set bbox [$names_win bbox S]
	set min_y [lindex $bbox 1]
	set max_y [lindex $bbox 3]
	set first_unit [lindex [$f.vscroll get] 0]
	set unit [expr double($CanvasConst(auto_incr))/($max_y - $min_y)]
	$re_win yview moveto [expr $first_unit + ($direction * $unit)]
	$names_win yview moveto [expr $first_unit + ($direction * $unit)]

	after 100 "REnzDoScroll $f $names_win $re_win $direction $scroll" 
    }
}

##############################################################################
#stand alone restriction enzyme display
proc REnzAutoScroll {names_win re_win y f scroll} {
    global CanvasConst

    #puts "start AUTOSCROLL y $y"

    set min_y 1
    set max_y [winfo height $names_win]

    if {$y < $min_y } {
	after $CanvasConst(auto_time) \
	    "REnzDoScroll $f $names_win $re_win -1 $scroll"
    } elseif {$y > $max_y} {
	after $CanvasConst(auto_time) \
	    "REnzDoScroll $f $names_win $re_win +1 $scroll"
	    
    }
} 

##############################################################################
#
proc REnzStopScroll {f} {
    global $f.auto_scroll

    set $f.auto_scroll 0
}

##############################################################################
#popup menu invoked when click on cut line with button-3
proc PopUpREnzMenu { f re_win X Y renz_info_cmd} {
    global initialCol
    global $re_win.prev_item

    if [winfo exists $re_win.m] {destroy $re_win.m}

    set item [GetREnzTag $re_win current]
    set index [GetREnzIndex $re_win $item]

    #starting colour
    set colour $initialCol($re_win,r_enz)

    #cmd to execute when ok button on colourbox pressed
    set ok_cmd "ConfigureREnzPlot $re_win $item"
    #cmd to execute when changing colours on colourbox
    set update_cmd "UpdateREnzPlot $re_win $item"
    #cmd to execute when cancel button on colourbox pressed
    set cancel_cmd "$re_win itemconfigure $item -fill $colour; InitialSettings $re_win $re_win $item r_enz"

    create_popup $re_win.m Commands
    #$re_win.m add command -label information \
	    -command "destroy $re_win.m; \
                      get_r_enz_info -id $id -enzyme $index -io $io -cnum $contig"
    $re_win.m add command -label information \
	    -command "destroy $re_win.m; \
                      uplevel #0 eval $renz_info_cmd -enzyme $index"
    $re_win.m add command -label configure \
	    -command "destroy $re_win.m; if {[winfo exists $re_win.cbox]} {return}; \
	    xtoplevel $re_win.cbox -resizable 0;\
	    wm title $re_win.cbox {Configure colours}; \
	    ColourBox $re_win.cbox $colour {$ok_cmd} {$update_cmd} {$cancel_cmd};\
	    set $re_win.prev_item 0"
     
    tk_popup $re_win.m [expr $X-20] [expr $Y-10]
}

##############################################################################
proc UpdateREnzPlot {canvas item colour } {

    #puts "update item $item colour $colour"
    $canvas itemconfigure $item -fill $colour
    InitialSettings $canvas $canvas $item r_enz
}

#############################################################################
proc ConfigureREnzPlot {canvas item colour } {
    global restoreCmd
    
    #puts "configure item $item colour $colour"
    eval $restoreCmd($canvas,r_enz)
    UpdateREnzPlot $canvas $item $colour
    InitialSettings $canvas $canvas $item r_enz
    eval $restoreCmd($canvas,r_enz)
}
