#job can be NEW_ROW,
proc create_canvas_element {c_win c_id e_win e_id width height row column orientation scale job} {
    global tk_utils_defs HORIZONTAL VERTICAL BOTH $e_win.column $e_win.row TOP_S BOTTOM_S canvas_menu env

    #puts "**************** create_canvas_element $c_win $e_win o=$orientation r=$row c=$column job=$job scale=$scale width=$width height=$height"
    set RULER_ROW [keylget tk_utils_defs CONTAINER.RULER_ROW]
    set RULER_COL [keylget tk_utils_defs CONTAINER.RULER_COL]
    set MAX_ROW [keylget tk_utils_defs CONTAINER.MAX_ROW]

    container_menubar $c_win $c_id $canvas_menu
    set_element_specs $e_win "menu [list $canvas_menu]"

    #puts "RULER_ROW $RULER_ROW $RULER_COL"

    #FIXME
    set e_specs [get_element_specs $e_win]

    #wm geometry $c_win {}

    set strand [keylget e_specs strand]

    #location of length ruler
    if {$orientation == $HORIZONTAL || $orientation == $BOTH} {
	set sb_h [find_widget_path $c_win[keylget tk_utils_defs CONTAINER.HSCROLL.WIN] column $column]
	set hscroll_cmd "$sb_h set "
    } 

    if {$orientation == $VERTICAL || $orientation == $BOTH} {
	set sb_v [find_widget_path $c_win[keylget tk_utils_defs CONTAINER.VSCROLL.WIN] row $row]

	set vscroll_cmd "$sb_v set "
    }
    if {$orientation == $HORIZONTAL || $orientation == $BOTH} {
	#grid_insert $c_win row $row 1

	if {[expr $strand & $TOP_S]} {
	    #always ensure that RULER_ROW is always RULER_ROW!

	    #puts "RULER_ROW $RULER_ROW"
	    #grid_delete $c_win row [expr $RULER_ROW - 2] 1 

	} else {
	    #always ensure that MAX_ROW is always MAX_ROW!
	    #grid_delete $c_win row [expr $MAX_ROW - 2] 1 
	}
    } else {
	#grid_insert $c_win column $column 1
    }

    if {$scale == "XY"} {
	#vertical length ruler scrollbar
	if {[info exists vscroll_cmd]} {
	    if {![winfo exists $sb_v]} {
		scrollbar $sb_v -orient vertical \
			-command "scroll_y $c_id \[get_element_row $sb_v\]"
		#grid $sb_v -row $row -column [expr $column+1] -sticky ns
		grid $sb_v -row $row -column [keylget tk_utils_defs CONTAINER.MAX_COL] -sticky ns
	    }
	}
	#horizontal length ruler scrollbar
	if {[info exists hscroll_cmd]} {
	    if {![winfo exists $sb_h]} {
		scrollbar $sb_h -orient horizontal \
			-command "scroll_x $c_id \[get_element_column $sb_h\]"
		grid $sb_h -row [expr $MAX_ROW - 1] -column $column -sticky ew
	    }
	}

	if {![info exists hscroll_cmd]} {

	    #puts "HORIZONTAL AMPLITUDE "
	    set sb_h [find_widget_path $c_win[keylget tk_utils_defs CONTAINER.HSCROLL.WIN] column $column]
	    set hscroll_cmd "$sb_h set "
	    scrollbar $sb_h -orient horizontal \
		    -command "scroll_x $c_id \[get_element_column $sb_h\]"
	    grid $sb_h -row [expr $MAX_ROW - 1] -column $column -sticky ew

	}

	#horizontal length, therefore vertical amplitude
	if {![info exists vscroll_cmd]} {
	    set sb_v [find_widget_path $c_win[keylget tk_utils_defs CONTAINER.VSCROLL.WIN] row $row]
	    
	    #puts "VERTICAL $sb_v $row [expr $column+1]"

	    set vscroll_cmd "$sb_v set "
	    if {![winfo exists $sb_v]} {
		scrollbar $sb_v -orient vertical \
			-command "scroll_y $c_id \[get_element_row $sb_v\]"
		#grid $sb_v -row $row -column [expr $column+1] -sticky ns
		grid $sb_v -row $row -column [keylget tk_utils_defs CONTAINER.MAX_COL] -sticky ns
	    }
	}
    }

    if {$scale == "X"} {
	#horizontal length ruler scrollbar
	if {[info exists hscroll_cmd]} {
	    if {![winfo exists $sb_h]} {
		scrollbar $sb_h -orient horizontal \
			-command "scroll_x $c_id \[get_element_column $sb_h\]"
		grid $sb_h -row [expr $MAX_ROW - 1] -column $column -sticky ew
	    }
	}
    }

    if {$scale == "Y"} {

    }
    #puts "canvasbox $e_win -width $width -height $height -bg green \
	    -bd 0 -highlightthickness 0 -xscrollcommand $hscroll_cmd\
	    -zoom_command $zoom_cmd"

    set zoom_cmd "container_zoom -container_id $c_id"

    if {[info exists hscroll_cmd] && [info exists vscroll_cmd]} {
	frame $c_win.f$e_id -bd 2 -relief groove
	canvasbox $e_win -width $width -height $height\
		-bd 0 -highlightthickness 0 -xscrollcommand $hscroll_cmd\
		-yscrollcommand $vscroll_cmd -zoom_command $zoom_cmd

    } elseif {[info exists hscroll_cmd]} {
	frame $c_win.f$e_id -bd 2 -relief groove
	canvasbox $e_win -width $width -height $height\
		-bd 0 -highlightthickness 0 -xscrollcommand $hscroll_cmd\
		-zoom_command $zoom_cmd
    } elseif {[info exists vscroll_cmd]} {
	frame $c_win.f$e_id -bd 2 -relief groove
	canvasbox $e_win -width $width -height $height \
		-bd 0 -highlightthickness 0 -yscrollcommand $vscroll_cmd\
		-zoom_command $zoom_cmd
    }

    canvas $c_win.f$e_id.l -width 8 -height 16
    $c_win.f$e_id.l create bitmap 5 10 -bitmap @[file join $env(STADTABL) 2bars_v.bmp]

    grid rowconfig $c_win $row -weight 1
    grid columnconfig $c_win $column -weight 1

    grid $c_win.f$e_id -row $row -column $column -sticky ewns
    pack $c_win.f$e_id.l -side left -anchor n
    pack $e_win -in $c_win.f$e_id -fill both -expand yes

    set $e_win.row $row
    set $e_win.column $column

    #set keypath [find_widget_path $c_win[keylget tk_utils_defs CONTAINER.KEYBOX.WIN] row $row]

    #if {![winfo exists $keypath]} {
	#keybox $keypath
	#add key for all results 
	
	#keybox_add $keypath $row -text all \
		-background white \
		-motion motion_canvas \
		-enter "enter_key $e_win -1" \
		-leave "leave_key $e_win" \
		-drop "drop_element $e_win" \
		-menu "seq_result_keybox_update -1 \[seq_result_names -element_id $e_id\]"

	#grid $keypath -row $row -column [expr $column+2] -sticky nsew
    #}

    canvas_element_bindings $c_win $e_id $e_win $c_win.f$e_id.l
    update idletasks
 
    if {[winfo height $c_win] != 1} {
       #wm geometry $c_win [winfo width $c_win]x$new_height
    }

    #puts "GRID [grid slaves $c_win -row 251]"

    #foreach win [grid slaves $c_win ] {
	#puts "GRID $win [grid info $win]"
    #}

}

#certain widgets (scrollbars, keyboxes) are defined per row/column (orient)
#but are given a unique number (w_cnt) so need to search if widget already
#exists for a row/column.
proc find_widget_path {widget orient num {job ""}} {
    global w_cnt 

    #puts "find_widget $widget $orient $num"

    if {![info exists w_cnt]} {
	set w_cnt 0
	#puts "here1 $widget$w_cnt"
	return $widget$w_cnt
    } elseif {[string compare $job "NEW"] == 0} {
	incr w_cnt
	#puts "here2  $widget$w_cnt"
	return $widget$w_cnt	
    } else {
	regexp {^(\.[^.]+)(\.[^.]+)$} $widget dummy parent win
	set children [winfo children $parent]
	foreach child $children {
	    if {[regexp "($win)(\[0-9\]*)$" $child dummy dummy val]} {
		if {[string compare $orient "row"] == 0} {
		    set widget_orient [get_element_row $widget$val]
		} else {
		    set widget_orient [get_element_column $widget$val]
		}
		if {$num == $widget_orient} {
		    #puts "here3 $widget$val"
		    return $widget$val
		}
	    }
	}
	incr w_cnt
	#puts "here4 $widget$w_cnt"
	return $widget$w_cnt
    }
}

proc get_next_element_frame { } {
    global id

    if {![info exists id]} {
	set id 0
    } else {
	incr id
    }
    return $id

}

proc create_canvas_ruler {c_win c_id orientation row column ruler_width ruler_height type} {
    global HORIZONTAL VERTICAL tk_utils_defs env

    #puts "***************create_canvas_ruler $c_win $c_id $orientation r=$row c=$column"

    set zoom_cmd "container_zoom -container_id $c_id"

    set element_info [get_new_element]

    set e_id [keylget element_info element_id]
    set e_win [keylget element_info element_win]

    set args ""

    if {$orientation == $HORIZONTAL} {

	set sb_h [find_widget_path $c_win[keylget tk_utils_defs CONTAINER.HSCROLL.WIN] column $column]
	set hscroll_cmd "$sb_h set "

	#puts "HORIZONTAL RULER $c_win$e_win $ruler_height $ruler_width $row $column"
	frame $c_win.f$e_id -bd 2 -relief groove
	canvasbox $c_win$e_win -width $ruler_width -height $ruler_height\
		-bd 0 -highlightthickness 0\
		-xscrollcommand $hscroll_cmd -zoom_command $zoom_cmd

	if {[string compare $type LENGTH] == 0} {
	    canvas $c_win.f$e_id.l -width 8 -height 16
	    $c_win.f$e_id.l create bitmap 5 10 -bitmap @[file join $env(STADTABL) 2bars_v.bmp]
	    pack $c_win.f$e_id.l -side left -anchor n
	    set args $c_win.f$e_id.l
	}

	#need to ensure the row weight is 0
	grid rowconfig $c_win $row -weight 0

	grid $c_win.f$e_id -row $row -column $column -sticky ew
	pack $c_win$e_win -in $c_win.f$e_id -fill both -expand yes
    }

    if {$orientation == $VERTICAL} {
	set sb_v [find_widget_path $c_win[keylget tk_utils_defs CONTAINER.VSCROLL.WIN] row $row]
	set vscroll_cmd "$sb_v set "

	#puts "VERTICAL RULER $e_win $ruler_height $ruler_width $row $column"
	frame $c_win.f$e_id -bd 2 -relief groove
	canvasbox $c_win$e_win -width $ruler_width -height $ruler_height\
		-bd 0 -highlightthickness 0 \
		-yscrollcommand $vscroll_cmd -zoom_command $zoom_cmd

	if {[string compare $type LENGTH] == 0} {
	    canvas $c_win.f$e_id.l -width 8 -height 16
	    $c_win.f$e_id.l create bitmap 5 10 -bitmap @[file join $env(STADTABL) 2bars_v.bmp]
	    pack $c_win.f$e_id.l -side left -anchor n
	    set args $c_win.f$e_id.l
	}

	grid $c_win.f$e_id -row $row -column $column -sticky ns
	pack $c_win$e_win -in $c_win.f$e_id -fill both -expand yes
    }

    canvas_element_bindings $c_win $e_id $c_win$e_win $args

    global $c_win$e_win.row $c_win$e_win.column
    set $c_win$e_win.row $row
    set $c_win$e_win.column $column

    update idletasks

    #foreach win [grid slaves $c_win ] {
	#puts "GRID $win [grid info $win]"
    #}

    return [list $e_id $c_win$e_win];
}

proc handle_leave {handle e_win} {
    global e_from e_to cursor_win

    set e_from $e_win
    set e_to ""
    set cursor_win $handle
    #puts "LEAVE $e_from"
}

proc handle_popup_menu {h_win X Y e_id} {

    if [winfo exists $h_win.m] {destroy $h_win.m}

    set m [create_popup $h_win.m Commands]
    seq_result_keybox_update -1 [seq_result_names -element_id $e_id] $m
    tk_popup $m $X $Y
}

proc canvas_element_bindings {c_win e_id e_win handle} {

    if {$handle != ""} {
	bind $handle <<select>> "collapse_e_win $c_win $e_win $e_id $handle"
	bind $handle <B2-Leave> "handle_leave $handle $e_win"
	bind $handle <Alt-Button-1> "handle_leave $handle $e_win"
	bind $handle <<move-release>> "drop_element $e_win %x %y"
	bind $handle <<move-drag>> "element_move_result $e_win %x %y %X %Y"

	bind $handle <<menu>> "handle_popup_menu $handle %X %Y $e_id"
    }

    bind $e_win <Any-Leave> "+delete_container_crosshair -element_id $e_id"
    bind $e_win <Any-Motion> "add_container_crosshair $c_win $e_id %x %y"
    bind $e_win <Any-Configure> "element_resize -element_id $e_id"
    bind $e_win <Any-Leave> "+unhighlight_graph %W"

    bind $e_win <Any-Motion> "+highlight_tag $c_win $e_win %x %y"

    bind $e_win <<menu>> "element_popup_menu %W %x %y %X %Y"
    bind $e_win <<select>> "+element_select_result %W %x %y"
    bind $e_win <<move>> "+element_select_result %W %x %y"
    bind $e_win <<move-drag>> "+element_move_result %W %x %y %X %Y "
    bind $e_win <<move-release>> "+element_drop_result %X %Y"
}

proc collapse_e_win {c_win e_win e_id handle} {
    global $e_win.pack_info tk_utils_defs HORIZONTAL env

    foreach i [$handle find all] {
        set bmap [lindex [$handle itemconfig $i -bitmap] 4]
    }
    $handle delete all

    if {[regexp {2bars_h.bmp} $bmap]} {
	#expanding element

	$handle configure -width 8 -height 16
	$handle create bitmap 5 10 -bitmap @[file join $env(STADTABL) 2bars_v.bmp]

	set e_struct [element_struct $e_id]
	set ruler_id [keylget e_struct ruler_id]
	if {[keylget e_struct ruler_id] != -1} {
	    set r_struct [element_struct [keylget e_struct ruler_id]]
	    eval pack [keylget r_struct win] [keylget $e_win.pack_info ruler]
	    grid $c_win.f$ruler_id
	    if {[keylget e_struct orientation] == $HORIZONTAL} {
		
		grid [keylget $e_win.pack_info sb_v]
	    }
	}
	eval pack $e_win [keylget $e_win.pack_info element]
    } else {
	#collapsing element

	$handle configure -width 16 -height 8
	$handle create bitmap 10 5 -bitmap @[file join $env(STADTABL) 2bars_h.bmp]

	keylset $e_win.pack_info row [get_element_row $e_win]

	set e_struct [element_struct $e_id]

	set ruler_id [keylget e_struct ruler_id]
	if {[keylget e_struct ruler_id] != -1} {
	    set r_struct [element_struct [keylget e_struct ruler_id]]
	    set ruler_win [keylget r_struct win]

	    keylset $e_win.pack_info ruler [pack info $ruler_win]
	    pack forget $ruler_win
	    grid remove $c_win.f$ruler_id
 
	    if {[keylget e_struct orientation] == $HORIZONTAL} {
		set sb_v [find_widget_path $c_win[keylget tk_utils_defs CONTAINER.VSCROLL.WIN] row [keylget $e_win.pack_info row]]
		
		keylset $e_win.pack_info sb_v $sb_v
		grid remove $sb_v
	    }
	}
	keylset $e_win.pack_info element [pack info $e_win]
	pack forget $e_win
    }
}

#called from C
proc draw_canvas_crosshairX {c_win e_win cx wx} {
    global tk_utils_defs 

    set colour [keylget tk_utils_defs CROSSHAIR.COLOUR]
    set line_width [keylget tk_utils_defs CROSSHAIR.LINE_WIDTH]
    
    set min_y [$e_win canvasy 0]
    set max_y [$e_win canvasy [winfo height $e_win]]

    if {[$e_win find withtag crosshair_x] == ""} {
	$e_win create line $cx $min_y $cx $max_y -tag "crosshair_x NS"\
		-fill $colour -width $line_width
	$e_win lower crosshair_x
    } else {
	$e_win coords crosshair_x $cx $min_y $cx $max_y
    }

    set buttons [keylget tk_utils_defs CONTAINER.BUTTONS.WIN]
    set pos1 [keylget tk_utils_defs CONTAINER.POS1.NAME]
    set pos2 [keylget tk_utils_defs CONTAINER.POS2.NAME]

    set x_format %d
    container_position_label $c_win$buttons$pos1 $wx $x_format

}

##############################################################################
#called from C
proc draw_canvas_crosshairY {c_win e_win cy wy} {
    global tk_utils_defs 

    set colour [keylget tk_utils_defs CROSSHAIR.COLOUR]
    set line_width [keylget tk_utils_defs CROSSHAIR.LINE_WIDTH]
    
    set min_x [$e_win canvasx 0]
    set max_x [$e_win canvasx [winfo width $e_win]]

    if {[$e_win find withtag crosshair_y] == ""} {
	$e_win create line $min_x $cy $max_x $cy -tag crosshair_y\
		-fill $colour -width $line_width
	$e_win lower crosshair_y
    } else {
	$e_win coords crosshair_y $min_x $cy $max_x $cy
    }
    set buttons [keylget tk_utils_defs CONTAINER.BUTTONS.WIN]
    set pos1 [keylget tk_utils_defs CONTAINER.POS1.NAME]
    set pos2 [keylget tk_utils_defs CONTAINER.POS2.NAME]

    set y_format %f
    container_position_label $c_win$buttons$pos2 $wy $y_format
}

proc motion_element {x y} {
    global e_to spin_defs e_from

    set win [winfo containing $x $y]
    
    if {$win == ""} {
	return
    }

    if {[string compare [winfo class $win] Canvas] != 0 && [string compare [winfo class $win] Seqed] != 0 } {
	return
    }
    #HACK to eliminate handle canvas
    set type [$win type 1]
    if {[string compare $type bitmap] == 0} {
	return
    }
    set e_to $win

    if {[string compare [winfo class $win] Canvas] == 0} {
	highlight_position_canvas $e_to $x $y
    } elseif {[string compare [winfo class $win] Seqed] == 0} {
	highlight_position_seqed $e_to $x $y
    }
}

proc enter_key {e_win result_id } {
    global tk_utils_defs

    if {$result_id == -1} {
	#if all - no status line
	return 
    }

    set c_win [winfo parent $e_win]

    set brief [keylget tk_utils_defs CONTAINER.BRIEF.WIN]
    $c_win$brief config -text [seq_get_brief -index $result_id]

}

proc leave_key {e_win} {
    global e_from e_to 

    set e_from $e_win
    set e_to ""
    #puts "LEAVE $e_from"
}

##############################################################################
proc drop_result {id x y} {
    after 50 [list where_result $id $x $y]
}

##############################################################################
proc drop_element {e_win x y} {

    after 50 [list where_element $e_win $x $y]
}

proc where_result {result_id x y} {
    global TOP MIDDLE BOTTOM e_to e_from prev_element HORIZONTAL

    if {![info exists e_from] || $e_from == ""} {
	return
    }
    global $e_from.location

    if {$e_to != ""} {
	set gr_type_from [get_element_gr_type [get_element_id $e_from]]
	set gr_type_to [get_element_gr_type [get_element_id $e_to]]

	set pos [find_result_position $e_from $e_to [set $e_from.location]]

	set row [lindex $pos 0]
	set column [lindex $pos 1]

	#puts "ROW $row COLUMN $column"

	if {$row == -2 } {
	    #puts "DO NOTHING"
	} elseif {$row == -1} {
	    #superimpose two elements
	    #puts superimpose
	    if {$gr_type_from == $gr_type_to} {
		move_result_superimpose $e_from $e_to $result_id
	    }
	} else {
	    #create new element for result
	    #puts "new element $e_to"
	    set c_win [winfo parent $e_to]
	    #move_result_new $e_from $e_to $result_id $row [get_current_column $c_win] 
	    set seq_id_h [find_seq_id -result_id $result_id -orient $HORIZONTAL]
	    #set column [seq_find_column -seq_id $seq_id_h]
	    move_result_new $e_from $e_to $result_id $row $column 
	}
    } else {
	#create new container window for result
	#puts "new container"
	init_result_container $e_from $result_id
    }
    
    #foreach win [grid slaves .container_win0] {
	#puts "$win [get_element_row $win]"
	#puts "$win [grid info $win]"
    #}

    remove_highlight_position_canvas
    catch {unset prev_element}
    set e_from ""
    set e_to ""

}

proc where_element {e_win x y} {
    global TOP MIDDLE BOTTOM e_to e_from prev_element

    if {![info exists e_from] || $e_from == ""} {
	return
    }
    global $e_from.location
    #puts "where_element e_from $e_from"

    if {$e_to != ""} {
	set gr_type_from [get_element_gr_type [get_element_id $e_from]]
	set gr_type_to [get_element_gr_type [get_element_id $e_to]]

	set pos [find_result_position $e_from $e_to [set $e_from.location]]

	set row [lindex $pos 0]
	set column [lindex $pos 1]

	#puts "ROW $row COLUMN $column"
	if {$row == -1} {
	    #superimpose two elements
	    #puts superimpose
	    if {$gr_type_from == $gr_type_to} {
		move_element_superimpose $e_from $e_to
	    }
	} elseif {[winfo parent $e_from] == [winfo parent $e_to]} {
	    #move element in same container
	   # move_element_same $e_from $e_to $row $column

	    #seems easier to just recreate element rather than figure out what
	    #to ungrid and grid again
	    move_result_new $e_from $e_to -1 $row $column
	} else {
	    #move element to a different container
	    move_element_different $e_from $e_to $row
	}

    } else {
	#create new container window for result
	#puts "new container"
	init_element_container $e_from
    }
    
    remove_highlight_position_canvas
    catch {unset prev_element}
    set e_from ""
    set e_to ""
}

proc remove_highlight_position_canvas { } {
    global prev_element e_from cursor_win

    catch {$cursor_win configure -cursor left_ptr}

    #if {![info exists prev_element]} {
	#return
    #}
    #catch {$prev_element delete highlight_rect}
}

proc highlight_position_canvas {canvas x y} {
    global TOP MIDDLE BOTTOM LEFT RIGHT prev_element e_from cursor_win
    global $e_from.location moving_cursor

    #puts "highlight_position_canvas $canvas $e_from moving_cursor [info exists moving_cursor]"

    if {[info exists moving_cursor]} {
	return
    }

    if {$canvas == ""} {
	catch {unset prev_element}
 	return
    }

    set posy [winfo rooty $canvas]
    set posx [winfo rootx $canvas]


    set win_ht [winfo height $canvas]
    set win_wd [winfo width $canvas]
    set ht3 [expr $win_ht / 3]
    set wd3 [expr $win_wd / 3]

    set top [expr $posy + $ht3]
    set bottom [expr $top + $ht3]

    set left [expr $posx + $wd3]
    set right [expr $left + $wd3]
	
    #puts "y $y posy $posy top $top bottom $bottom"
	
    set y1 [expr $posy - ($posy + $win_ht)]
    set x1 [expr $posx - ($posx + $win_wd)]
    set d1_m [expr double ($y1) / $x1]
    set d1_c [expr $posy - ($posx * $d1_m)]
    
    set y1 [expr ($posy + $win_ht) - $posy]
    set x1 [expr $posx - ($posx + $win_wd)]
    set d2_m [expr double ($y1) / $x1]
    set d2_c [expr ($posy + $win_ht) - ($posx * $d2_m)]

    set d1_x [expr double($y - $d1_c) / $d1_m]
    set d2_x [expr double($y - $d2_c) / $d2_m]

    set d1_y [expr $d1_m * $x + $d1_c]
    set d2_y [expr $d2_m * $x + $d2_c]
	
    #above d1 && above d2
    if {[expr $y > $d1_y] && [expr $y > $d2_y]} {
	#puts BOTTOM
	set location $BOTTOM
	$cursor_win configure -cursor bottom_side
    }

    #above d1 && below d2
    if {[expr $y > $d1_y] && [expr $y <= $d2_y]} {
	#puts LEFT
	set location $LEFT
	$cursor_win configure -cursor left_side
    }

    #below d1 && above d2
    if {[expr $y <= $d1_y] && [expr $y > $d2_y]} {
	#puts RIGHT
	set location $RIGHT
	$cursor_win configure -cursor right_side
    }
    #below d1 && below d2
    if {[expr $y <= $d1_y] && [expr $y <= $d2_y]} {
	#puts TOP
	set location $TOP
	$cursor_win configure -cursor top_side
    }

    if {$y >= $top && $y < $bottom && $x >= $left && $x < $right} {
	#puts MIDDLE
	set location $MIDDLE
	$cursor_win configure -cursor dotbox
    }

    set $e_from.location $location
}

#NB assume ruler to left and scrollbar to right
proc get_element_num_columns {e_win job} {
    global RULER_AMP LEFT RIGHT

    #puts get_element_num_columns

    set num_column 0

    set row [get_element_row $e_win]
    set column [get_element_column $e_win]

    set c_win [winfo parent $e_win]
    set e_struct [element_struct [get_element_id $e_win]]
    set children_pos [keylget e_struct children_position]
    
    if {$job == $LEFT} {
	#set win [grid slaves $c_win -row $row -column [expr $column - 1]]
	#check if next column is a ruler
	#set e_struct [element_struct [get_element_id $win]]
	
	if {$children_pos < 0} {
	    set num_column [keylget e_struct num_children_left]
	}
	#puts "NUM_COL $num_column"
    }    

    if {$job == $RIGHT} {
	#set win [grid slaves $c_win -row $row -column [expr $column + 1]]
	#if {[string compare [winfo class $win] scrollbar] == 0} {
	 #   incr num_column
	#}
	if {$children_pos > 0} {
	    set num_column [keylget e_struct num_children_right]
	}
    }
    return $num_column
}

#automatic element positioning when create element
proc get_element_position {c_win e_win_to strand orientation} {
    global tk_utils_defs HORIZONTAL VERTICAL TOP MIDDLE BOTTOM LEFT RIGHT

    if {$orientation & $HORIZONTAL} {
	if {[string compare $e_win_to ""] != 0 } {
	    set row [expr [get_next_row $c_win $strand] - 1]
	    set column [get_element_column $e_win_to]
	    for {set i $row} {$i < [keylget tk_utils_defs CONTAINER.MAX_ROW]} {incr i} {
		if {[llength [grid slaves $c_win -row $i -column $column]] == 0} {
		    set next_row $i
		    break
		}
	    }
	    if {[llength [grid slaves $c_win -row [expr $row + 1] -column $column]] > 0} {
		grid_insert $c_win row $next_row
	    }
	    set row $next_row
	} else {
	    set row [get_next_row $c_win $strand]
	    set column [keylget tk_utils_defs CONTAINER.RULER_COL]
	}
	
    } else {
	if {[string compare $e_win_to ""] != 0 } {
	    set list [find_result_position "" $e_win_to $LEFT]
	    set row [lindex $list 0]
	    set column [lindex $list 1]
	} else {
	    set row [get_next_row $c_win $strand]
	    set column [keylget tk_utils_defs CONTAINER.RULER_COL]
	}
    }
 
    return "$row $column"
}

##############################################################################
#find where to position the result when let go of the mouse
proc find_result_position {e_from e_to location} {
    global tk_utils_defs TOP MIDDLE BOTTOM LEFT RIGHT $e_from.location RULER_AMP

    set c_win [winfo parent $e_to]
    global $c_win.info

    set row [get_element_row $e_to]
    set row_from [get_element_row $e_from]
    set first [keylget $c_win.info row_first]

    set column [get_element_column $e_to]
    set column_from [get_element_column $e_from]

    #puts "find_result_position $canvas $row_from $row $first [element_info [get_element_id $e_from] num_results]"

    if {$row == $row_from && [element_info [get_element_id $e_from] num_results] == 1 && $e_from == $e_to} {
	    return -2
    }

    if {$column == $column_from && [element_info [get_element_id $e_from] num_results] == 1 && $e_from == $e_to} {
	return -2
    }

    if {$location == $MIDDLE} {
	return -1
    }

    if {$location == $TOP} {

	if {[llength [grid slaves $c_win -row [expr $row - 1] -column $column]] > 0} {
	    grid_insert $c_win row $row
	} else {
	    set row [expr $row - 1]
	}
	#grid $canvas -row $row -column $column -sticky nesw
	#puts "TOP $row $column"
	return "$row $column"
    }

    if {$location == $BOTTOM} {

	if {[string compare $e_from ""] == 0} {
	    for {set i $row} {$i < [keylget tk_utils_defs CONTAINER.MAX_ROW]} {incr i} {
		if {[llength [grid slaves $c_win -row $i -column $column]] == 0} {
		    set next_row $i
		    break
		}
	    }
	} else {
	    set next_row [expr $row + 1]
	}

	if {[llength [grid slaves $c_win -row [expr $row + 1] -column $column]] > 0} {
	    grid_insert $c_win row $next_row
	}
	set row $next_row
	#puts "BOTTOM $row $column"
	return "$row $column"
   }

    if {$location == $RIGHT} {

	set next_column [find_next_column $e_to $RIGHT]

	if {[string compare $e_from ""] != 0} {
	    set num_columns [get_element_num_columns $e_from $LEFT]
	} else {
	    set num_columns 0
	}
	if {[llength [grid slaves $c_win -row $row -column $next_column]] > 0} {
	    grid_insert $c_win column $next_column [expr $num_columns + 1]
	}
	set column [expr $next_column + $num_columns]
	#puts "RIGHT $row $column $num_columns"
	return "$row $column"
    }

    if {$location == $LEFT} {

	set next_column [find_next_column $e_to $LEFT]

	#are there any windows to the left?
	if {[llength [grid slaves $c_win -row $row -column $next_column]] > 0} {
	    #how many columns does e_from have?
	    if {[string compare $e_from ""] != 0} {
		set num_columns [get_element_num_columns $e_from $LEFT]
	    } else {
		set num_columns 0
	    }
	    grid_insert $c_win column $next_column [expr $num_columns + 1]
	}
	#puts "LEFT $row $column"
	return "$row $next_column"
    }
}

proc find_result_positionOLD {e_from canvas} {
    global TOP MIDDLE BOTTOM LEFT RIGHT $e_from.location

    set c_win [winfo parent $canvas]
    global $c_win.info

    set row [get_element_row $canvas]
    set row_from [get_element_row $e_from]
    set first [keylget $c_win.info row_first]

    set column [get_element_column $canvas]
    set column_from [get_element_column $e_from]

    #puts "find_result_position $canvas $row_from $row $first [element_info [get_element_id $e_from] num_results]"

    if {$row == $row_from && [element_info [get_element_id $e_from] num_results] == 1 && $e_from == $canvas} {
	    return -2
    }

    if {$column == $column_from && [element_info [get_element_id $e_from] num_results] == 1 && $e_from == $canvas} {
	return -2
    }

    #foreach tag [$canvas gettags highlight_rect] {
	#if {[string compare [string range $tag 0 1] p_] == 0} {
	  #  set pos [string trim $tag p_]
	 #   if {$pos == $TOP} {
	 #}
     #}
 #}

    set r $row
    set c $column

    if {[set $e_from.location] == $TOP} {
	#set row [expr $row - 1]
	if {$row < $first} {
	    set row $first
	}
	#puts "ROW1 $row"
	set r [expr $row - 1]
    } elseif {[set $e_from.location] == $MIDDLE} {
	return -1;
    } elseif {[set $e_from.location] == $BOTTOM} {
	#puts "ROW2 [expr $row + 1]"
	set r [expr $row + 1]
    } elseif {[set $e_from.location] == $LEFT} {
	set c [expr $column - 1]
    } else {
	set c [expr $column + 1]
    }

    return "$r $c"
}

#make new container by moving single result
proc init_result_container {e_from result_id} {
    global TOP_S tk_utils_defs HORIZONTAL VERTICAL

    #puts init_result_container

    set c_info [get_element_info $e_from]

    set c_id_old [keylget c_info container_id]
    set c_win_old [keylget c_info container_win]$c_id_old
    set e_id_old [keylget c_info element_id]
    set e_win_old [keylget c_info element_win]$e_id_old

    #check if e_from container contains only 1 result
    set num [num_container_results -container_id $c_id_old]

    if {$num == 1} {
	return
    }

    set element_info [get_new_container]
    set c_id_new [keylget element_info container_id]
    set c_win_new [keylget element_info container_win]

    set element_info [get_new_element]
    set e_id_new [keylget element_info element_id]
    set e_win_new [keylget element_info element_win]

    set e_specs [get_element_specs $e_from]
    set orientation [keylget e_specs orientation]

    #if orientation is currently vertical, swap to horizontal for new container
    if {$orientation == $VERTICAL} {
	set orientation $HORIZONTAL
    }
    set_element_specs $e_from "orientation $orientation"

    create_container $c_win_new $c_id_new [keylget e_specs title]

    if {$result_id != -1} {
	set win_size [seq_result_info -index $result_id -option win_size]
	#puts "win_size $win_size $result_id"
	set width [lindex $win_size 0]
	set height [lindex $win_size 1]
    } else {
	set height [winfo height $e_from]
	set width [winfo width $e_from]
    }

    set_element_specs $c_win_new$e_win_new [join $e_specs]

    set seq_id_h [find_seq_id -result_id $result_id -orient $HORIZONTAL]
    set column [seq_find_column -seq_id $seq_id_h]

    create_canvas_element $c_win_new $c_id_new $c_win_new$e_win_new $e_id_new $width $height [get_next_row $c_win_new $TOP_S] $column $orientation [keylget e_specs scale] NEW_ROW

    #set row_old [element_info $e_id_old row]
    set row_old [get_element_row $c_win_old$e_win_old]

    update_container -new_container_id $c_id_new -new_container_win $c_win_new -old_container_id $c_id_old -new_element_id $e_id_new -new_element_win $c_win_new$e_win_new -old_element_id $e_id_old -job NEW -result_id $result_id -new_orientation $orientation

    set row_new [element_info $e_id_new row]

    result_list_update $c_win_old
    result_list_update $c_win_new
} 

proc move_result_superimpose {e_win_old e_win_new result_id} {

    #puts move_result_superimpose

    #Don't bother to superimpose if not moving the window
    if {$e_win_old == $e_win_new} {
	return
    }

    set c_info [get_element_info $e_win_old]
    set c_id_old [keylget c_info container_id]
    set c_win_old [keylget c_info container_win]$c_id_old
    set e_id_old [keylget c_info element_id]

    set c_info [get_element_info $e_win_new]
    set c_win_new [keylget c_info container_win]
    set e_id_new [keylget c_info element_id]

    set row_old [get_element_row $e_win_old]
    set row_new [get_element_row $e_win_new]

    set e_specs [get_element_specs $e_win_old]
    set orientation [keylget e_specs orientation]

    update_container -old_container_id $c_id_old -old_element_id $e_id_old -new_element_id $e_id_new -job SUPERIMPOSE -result_id $result_id -new_orientation $orientation
}

#move result to new element
proc move_result_new {e_win_old e_to result_id row column } {
    global TOP_S tk_utils_defs HORIZONTAL VERTICAL $e_win_old.location LEFT RIGHT TOP BOTTOM

    #puts "move_result_new $e_to $e_win_old $row $column $result_id"

    set c_info [get_element_info $e_to]
    set c_id [keylget c_info container_id]
    set c_win [keylget c_info container_win]$c_id

    set c_info [get_element_info $e_win_old]
    set e_id_old [keylget c_info element_id]

    set element_info [get_new_element]
    set e_id_new [keylget element_info element_id]
    set e_win_new [keylget element_info element_win]

    set e_specs [get_element_specs $e_win_old]
    #set from_orientation [keylget e_specs orientation]

    set e_struct [element_struct [get_element_id $e_win_old]]
    set from_orientation [keylget e_struct orientation]
    set from_row [get_element_row $e_win_old]
    set from_column [get_element_column $e_win_old]

    set e_struct [element_struct [get_element_id $e_to]]
    set to_orientation [keylget e_struct orientation]

    #puts "from_row $from_row $row from_column $from_column $column location [set $e_win_old.location] $from_orientation"

    if {($to_orientation & $VERTICAL && $from_orientation & $HORIZONTAL) ||
    ($to_orientation & $HORIZONTAL && $from_orientation & $VERTICAL)} {
	#check if need to change orientation
	if {($from_orientation == $HORIZONTAL && ([set $e_win_old.location] == $LEFT || [set $e_win_old.location] == $RIGHT)) ||  
	($from_orientation == $VERTICAL && ([set $e_win_old.location] == $TOP || [set $e_win_old.location] == $BOTTOM))} {  
	    if {$from_orientation & $HORIZONTAL} {
		set from_orientation $VERTICAL
	    } else {
		set from_orientation $HORIZONTAL
	    }
	}
    }

    if {$result_id == -1} {
	if {$from_orientation & $HORIZONTAL} {
	    set width [winfo width $e_win_old]
	    set height [winfo height $e_win_old]
	} 
	if {$from_orientation & $VERTICAL} {
	    set height [winfo width $e_win_old]
	    set width [winfo height $e_win_old]
	} 
    } else {

	if {$from_orientation & $HORIZONTAL} {
	    set win_size [seq_result_info -index $result_id -option win_size]
	    set width [lindex $win_size 0]
	    set height [lindex $win_size 1]
	}
	if {$from_orientation == $VERTICAL} {
	    set win_size [seq_result_info -index $result_id -option win_size]
	    set height [lindex $win_size 0]
	    set width [lindex $win_size 1]
	} 
    }

    #puts "width $width height $height result_id $result_id [keylget e_specs scale]"

    #puts "e_specs $e_specs"
    set_element_specs $c_win$e_win_new [join $e_specs]
    set_element_specs $c_win$e_win_new "orientation $from_orientation"

    create_canvas_element $c_win $c_id $c_win$e_win_new $e_id_new $width $height $row $column $from_orientation [keylget e_specs scale] NEW_ROW

    set row_old [get_element_row $e_win_old]

    update_container -old_container_id $c_id -new_element_id $e_id_new -new_element_win $c_win$e_win_new -old_element_id $e_id_old -job NEW -result_id $result_id -new_container_id $c_id -new_orientation $from_orientation

    set row_new [element_info $e_id_new row]

    result_list_update $c_win
    #puts "END move_result_new"
}

#move element onto a existing element
proc move_element_superimpose {e_from e_to} {

    #puts move_element_superimpose

    #Don't bother to superimpose if not moving the window
    if {$e_from == $e_to} {
	return
    }

    set einfo [get_element_info $e_from]
    set c_id_old [keylget einfo container_id]
    set c_win_old [keylget einfo container_win]$c_id_old
    set e_id_old [keylget einfo element_id]
    set e_win_old [keylget einfo element_win]$e_id_old

    set einfo [get_element_info $e_to]
    set c_id_new [keylget einfo container_id]
    set c_win_new [keylget einfo container_win]$c_id_new
    set e_id_new [keylget einfo element_id]

    #must move_key before update_container otherwise I delete the key before
    #had chance to move it!
    #move_key $c_win_old $c_win_new $e_id_old $e_id_new $c_win_old$e_win_old $e_to

    #puts "NEW $e_id_new OLD $e_id_old"
    set e_specs [get_element_specs $c_win_old$e_win_old]
    set orientation [keylget e_specs orientation]

   update_container -old_container_id $c_id_old -old_element_id $e_id_old -new_element_id $e_id_new -new_element_win $e_to -new_container_id $c_id_new -new_container_win $c_win_new -job SUPERIMPOSE -new_orientation $orientation

}

#move element to new element in same container
proc move_element_same {e_from e_to row column} {
    global tk_utils_defs

    #puts move_element_same

    #don't bother to move if not moving element
    if {$e_from == $e_to} {
	return
    }
    set einfo [get_element_info $e_from]
    set e_id_old [keylget einfo element_id]
    set c_win [winfo parent $e_from]

    set frame_from [keylget einfo container_win][keylget einfo container_id].f[keylget einfo element_id]

    set einfo [get_element_info $e_to]
    set frame_to [keylget einfo container_win][keylget einfo container_id].f[keylget einfo element_id]

    array set info_from [grid info $frame_from]
    array set info_to [grid info $frame_to]

    set row_from [get_element_row $e_from]
    set column_from [get_element_column $e_from]

    #puts "row=$row column=$column from=$row_from $column_from"

    set win_list [grid_find $c_win row [get_element_row $e_from]]

    #set win_list [grid slaves $c_win]

    foreach win $win_list {
	array set info_$win [grid info $win]
	grid forget $win
    }

    #remove element row from grid
    grid_delete $c_win row $row_from $info_from(-rowspan)

    #insert new element row into grid
    grid_insert $c_win row $row $info_to(-rowspan)

    #repack the frames
    grid rowconfig $c_win $row -weight 1 -minsize 40
    grid columnconfig $c_win $column -weight 1 -minsize 40

    foreach win $win_list {
	set info_${win}(-row) $row
	#set info_${win}(-column) $column
	#puts "grid $win [array get info_$win]"
	eval grid $win [array get info_$win]
    }
}

#move element to new element in different container
proc move_element_different {e_from e_to row} {
    global tk_utils_defs

    #puts move_element_different

    set einfo [get_element_info $e_from]
    set c_id_old [keylget einfo container_id]
    set c_win_old [keylget einfo container_win]$c_id_old
    set e_id_old [keylget einfo element_id]
    set e_win_old [keylget einfo element_win]$e_id_old

    set einfo [get_element_info $e_to]
    set c_id_new [keylget einfo container_id]
    set c_win_new [keylget einfo container_win]$c_id_new

    set element_info [get_new_element]
    set e_id_new [keylget element_info element_id]
    set e_win_new [keylget element_info element_win]

    set height [winfo height $e_from]
    set width [winfo width $e_from]

    set e_specs [get_element_specs $e_from]
    set orientation [keylget e_specs orientation]

    set_element_specs $c_win_new$e_win_new [join $e_specs]

    #puts "CURRENT COLUMN [get_current_column $c_win_new]"
    

    create_canvas_element $c_win_new $c_id_new $c_win_new$e_win_new $e_id_new $width $height $row [get_current_column $c_win_new] [keylget e_specs orientation] [keylget e_specs scale] NEW_ROW

    #must move_key before update_container otherwise I delete the key before
    #had chance to move it!
    #move_key $c_win_old $c_win_new $e_id_old $e_id_new $c_win_old$e_win_old $c_win_new$e_win_new


    #puts "move_element_different $c_win_new$e_win_new"

    update_container -new_container_id $c_id_new -new_container_win $c_win_new -old_container_id $c_id_old -new_element_id $e_id_new -new_element_win $c_win_new$e_win_new -old_element_id $e_id_old -job NEW -new_orientation $orientation

}

#create a new container for element
proc init_element_container {e_from } {
    global TOP_S tk_utils_defs

    set seqed 0
    set c_info [get_element_info $e_from]

    set c_id_old [keylget c_info container_id]
    set c_win_old [keylget c_info container_win]$c_id_old
    set e_id_old [keylget c_info element_id]
    set e_win_old [keylget c_info element_win]$e_id_old

    #check if e_from container contains only 1 result
    set num [num_container_elements -container_id $c_id_old]

    if {$num == 1} {
	return
    }

    set element_info [get_new_container]
    set c_id_new [keylget element_info container_id]
    set c_win_new [keylget element_info container_win]

    set e_specs [get_element_specs $e_from]
    set orientation [keylget e_specs orientation]

    create_container $c_win_new $c_id_new [keylget e_specs title]

    set_element_specs $c_win_new$e_win_old [join $e_specs]
    
    if {[string compare [winfo class $e_from] Seqed] == 0} {

	
    } elseif {[string compare [winfo class $e_from] Canvas] == 0} {
	set height [winfo height $e_from]
	set width [winfo width $e_from]
    
	create_canvas_element $c_win_new $c_id_new $c_win_new$e_win_old $e_id_old $width $height [get_next_row $c_win_new $TOP_S] [get_next_column $c_win_new] [keylget e_specs orientation] [keylget e_specs scale] NEW_ROW
    }

    #must move_key before update_container otherwise I delete the key before
    #had chance to move it!
    #move_key $c_win_old $c_win_new $e_id_old $e_id_old $c_win_old$e_win_old $c_win_new$e_win_old

    update_container -new_container_id $c_id_new -new_container_win $c_win_new -old_container_id $c_id_old -new_element_id $e_id_old -new_element_win $c_win_new$e_win_old -old_element_id $e_id_old -job NEW -new_orientation $orientation

    result_list_update $c_win_old
    result_list_update $c_win_new

}

proc configure_canvas_result {e_id e_win result_id line_width colour} {

    config_result -element_id $e_id -result_id $result_id -width $line_width -fill $colour

    #puts "configure_canvas_result id$result_id"
    $e_win itemconfigure id$result_id -width $line_width -fill $colour
}

proc update_canvas_result {e_id e_win result_id line_width colour} {

    $e_win itemconfigure id$result_id -width $line_width -fill $colour
}

proc get_result_id_tag {win item} {

    foreach tag [$win gettags $item] {
	if {[regexp {id([0-9]+)$} $tag dummy id]} {
	    return $id
	}
    }
    return -1
}

proc brief_tag {brief win item result_id} {

    foreach tag [$win gettags $item] {
	if {[regexp {d([0-9]+).([0-9]+)$} $tag dummy darrays darray]} {
	    set line [seq_get_brief_tag -array_type "d" -arrays $darrays -array $darray -result_id $result_id]
	    $brief configure -text $line
	}
    }
}

proc unhighlight_graph {e_win} {
    global $e_win.prev

    #puts "unhighlight_graph $e_win"
    if {[info exists $e_win.prev]} {
	if {[keylget $e_win.prev id] != -1} {
	    set item [keylget $e_win.prev id]
	    set colour [keylget $e_win.prev colour]
	    set line_width [keylget $e_win.prev line_width]

	   # puts " UNHIGHLIGHT [set $e_win.prev]"
	   # puts [$e_win find all]

	    #puts "UNHIGHLIGHT $item"
	    $e_win itemconfig $item -fill $colour -width $line_width
	    unset $e_win.prev
	}
    }
}


proc select_result {e_win x y cmd} {
    global $e_win.prev tk_utils_defs element_item current_item 

    if {![info exists $e_win.prev]} {
	keylset $e_win.prev id -1
    }

    set tags [$e_win find withtag NS]
    foreach tag $tags {
	$e_win itemconfigure $tag -state hidden
    }

    set item [$e_win find closest [$e_win canvasx $x] [$e_win canvasy $y]]
    set current_item $item

    #puts "item $current_item"
    foreach tag $tags {
	$e_win itemconfigure $tag -state normal
    }

    set tags [$e_win gettags $item]
    foreach tag $tags {
	if {[regexp {id([0-9]+)$} $tag dummy result_id]} {

	    if {$item != [keylget $e_win.prev id]} {
		unhighlight_graph $e_win

		keylset $e_win.prev id $item
		keylset $e_win.prev colour [lindex [$e_win itemconfig $item -fill] 4]
		keylset $e_win.prev line_width [lindex [$e_win itemconfig $item -width] 4]
		#puts "PREV $e_win [set $e_win.prev]"

	    }
	    #puts "HIGHLIGHT $item $cmd"
	    $e_win itemconfig $item -fill [keylget tk_utils_defs ELEMENT.SELECT_COLOUR] -width 0.0
	    eval $cmd
	}
    }
}

proc highlight_tag {c_win e_win x y} {
    global tk_utils_defs

    set brief $c_win[keylget tk_utils_defs CONTAINER.BRIEF.WIN]

    select_result $e_win $x $y "$brief configure -text \[seq_get_brief -index \$result_id\]"
}

proc element_popup_menu_cmd {e_win X Y result_id} {

    if [winfo exists $e_win.m] {destroy $e_win.m}
    
    set name [seq_result_key_name -index $result_id]
    set m [create_popup $e_win.m $name]
    tk_popup $m $X $Y
    seq_result_keybox_popup $m $result_id
}

proc element_popup_menu {e_win x y X Y} {

    select_result $e_win $x $y "element_popup_menu_cmd $e_win $X $Y \$result_id"
} 

proc element_select_result {e_win x y} {
    global element_item e_from e_to cursor_win 

    select_result $e_win $x $y "set element_item \$result_id"

    set e_from $e_win
    set e_to ""
    set cursor_win $e_win
}

proc element_set_result {e_win result_id} {
    global element_item e_from e_to

    set element_item $result_id

    set e_from $e_win
    set e_to ""
    set cursor_win $e_win
}

proc element_move_result {e_win x y X Y} {
    global e_to moving_cursor current_item moving_result

    set tags [$e_win gettags $current_item]
    foreach tag $tags {
	if {$tag == "cursor"} {
	    return
	}
    }

    if {[info exists moving_cursor]} {
	return
    }
    set moving_result 1

    set e_to ""
    motion_element $X $Y

}

proc element_drop_result { X Y } {
    global element_item moving_cursor moving_result
    
    if {[info exists moving_cursor]} {
	return
    }
    if {![info exists moving_result] && ![info exists moving_cursor]} {
	return
    }

    if {[info exists element_item]} {
	drop_result $element_item $X $Y
    } else {
	remove_highlight_position_canvas 
    }
    if {[info exists moving_result]} {
	unset moving_result
    }
}



