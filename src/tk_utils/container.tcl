#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

proc container_start_shutdown {c_win c_id} {
    global $c_win.crosshair

    set $c_win.crosshair 0
    container_shutdown -container_id $c_id

    destroy $c_win
}

proc container_results_manager {help_cmd args} {
     global tk_utils_defs

    #puts "***************** container_results_manager ******************"

    set t [keylget tk_utils_defs CONTAINER.RESULTS.WIN]

    if {[xtoplevel $t] == ""} return
    wm title $t "Results manager"

    result_list_create $t $help_cmd
    eval result_list_update $args
  
}

proc container_menubar {c_win c_id menu} {
    global container_menu tk_utils_defs

    set mb [keylget tk_utils_defs CONTAINER.MENUBAR]

    if {![winfo exists $c_win$mb]} {
	$c_win configure -menu $c_win$mb
	menu $c_win$mb
    }
    merge_menus $menu $c_win$mb
}

proc container_buttonbar {t buttons zoom_cmd} {
    global tk_utils_defs $t.crosshair 

    #crosshair checkbutton
    checkbutton $buttons.cross -text crosshairs -variable $t.crosshair

    #position labels
    set pos1 [keylget tk_utils_defs CONTAINER.POS1.NAME]
    set pos2 [keylget tk_utils_defs CONTAINER.POS2.NAME]
    label $buttons$pos1 -bd 2 -relief sunken -width [keylget tk_utils_defs CONTAINER.POS1.WIDTH]
    label $buttons$pos2 -bd 2 -relief sunken -width [keylget tk_utils_defs CONTAINER.POS1.WIDTH]

    #zoom out button
    button $buttons.zoom_out -text "zoom out" -command "eval $zoom_cmd"

    #zoom in buttons
    button $buttons.zoomin10 -text "+10%" \
	-command "eval $zoom_cmd -amount 0.05"
    button $buttons.zoomin50 -text "+50%" \
	-command "eval $zoom_cmd -amount 0.1666" 
     
    pack $buttons.zoomin10 $buttons.zoomin50 $buttons.zoom_out $buttons.cross \
	    $buttons$pos1 $buttons$pos2 -side left -fill x -expand yes
}

proc add_container_column {c_win container_id column} {
    global tk_utils_defs

    set MAX_ROW [keylget tk_utils_defs CONTAINER.MAX_ROW]

    #Horizontal scrollbar
    set sb_h [keylget tk_utils_defs CONTAINER.HSCROLL.WIN]$column
    scrollbar $c_win$sb_h -orient horizontal \
	    -command "scroll_x $container_id $column"

    grid $c_win$sb_h -row [expr $MAX_ROW - 1] -column $column -sticky ew

}

proc create_container {t container_id title} {
    global tk_utils_defs
    global $t.info HORIZONTAL

    set MAX_ROW [keylget tk_utils_defs CONTAINER.MAX_ROW]
    set MAX_COL [keylget tk_utils_defs CONTAINER.MAX_COL]

    set RULER_ROW [keylget tk_utils_defs CONTAINER.RULER_ROW]
    set RULER_COL [keylget tk_utils_defs CONTAINER.RULER_COL]

    if {[xtoplevel $t] == ""} return
    #fix_maxsize $t

    wm protocol $t WM_DELETE_WINDOW "container_start_shutdown $t $container_id"
    wm title $t $title

    #initialise info structure
    keylset $t.info row_top [expr $RULER_ROW/2]
    keylset $t.info row_bottom $RULER_ROW
    keylset $t.info row_first 2
    keylset $t.info column $RULER_COL
    keylset $t.info column_left [expr $RULER_COL + ($MAX_COL/2)]

    ##########################################################################
    # Main Menu Bar
    #container_menubar $t $container_id

    ##########################################################################
    #Button bar
    set buttons [keylget tk_utils_defs CONTAINER.BUTTONS.WIN]
    frame $t$buttons

    set zoom_cmd "container_zoom -container_id $container_id"
    container_buttonbar $t $t$buttons $zoom_cmd

    #add_container_column $t $container_id [expr $RULER_COL + 1]

    ##########################################################################
    #Brief line
    set brief [keylget tk_utils_defs CONTAINER.BRIEF.WIN]
    label $t$brief

    #puts "BRIEF [expr $RULER_COL + 0]"
    grid $t$buttons -row 2 -column 0 -sticky ew -columnspan $MAX_COL
    grid $t$brief -row $MAX_ROW -column [expr $RULER_COL + 0] -sticky ew -columnspan $MAX_COL

}

proc get_prev_row {c_win strand} {
    global $c_win.info TOP_S BOTTOM_S

    if {[expr $strand & $TOP_S]} {
	set top [keylget $c_win.info row_top]
	incr top -1
	keylset $c_win.info row_top $top
	return $top
    } elseif {[expr $strand & $BOTTOM_S]} {
	set bottom [keylget $c_win.info row_bottom]
	incr bottom -1
	keylset $c_win.info row_bottom $bottom
	return $bottom
    }
}

proc get_next_row {c_win strand} {
    global $c_win.info TOP_S BOTTOM_S

    #puts "get_next_row $strand"

    if {[expr $strand & $TOP_S]} {
	set top [keylget $c_win.info row_top]
	incr top
	keylset $c_win.info row_top $top
	return $top
    } elseif [expr $strand & $BOTTOM_S] {
	set bottom [keylget $c_win.info row_bottom]
	incr bottom
	keylset $c_win.info row_bottom $bottom
	return $bottom
    }
}

proc get_current_column {c_win} {
    global $c_win.info 

    set column [keylget $c_win.info column]
    return $column
}

proc get_next_column {c_win} {
    global $c_win.info 

    set column [keylget $c_win.info column]

    #puts "get_next_column $c_win $column"

    incr column
    keylset $c_win.info column $column
    return $column

}

#the element id is always the last number of the element_win
proc get_element_id {element_win} {
    regexp {[0-9]+$} $element_win id
    return $id

}

proc get_element_row { e_win} {

    if {[winfo exists $e_win]} {
	set einfo [get_element_info $e_win]
	if {$einfo == ""} {
	    set e_frame $e_win
	} else {
	    set e_frame [keylget einfo container_win][keylget einfo container_id].f[keylget einfo element_id]
	}
	array set info [grid info $e_frame]
	#parray info
	#puts "array size [array size info]"

	if {[array size info] == 0} {
	    return -1
	}
	#puts "get_element_row $e_win $info(-row)"
	return $info(-row)
    } else {
	return -1
    }
}

proc get_element_column { e_win} {

    #puts "get_element_column $e_win"
    if {[winfo exists $e_win]} {
	set einfo [get_element_info $e_win]
	if {$einfo == ""} {
	    set e_frame $e_win
	} else {
	    set e_frame [keylget einfo container_win][keylget einfo container_id].f[keylget einfo element_id]
	}
	array set info [grid info $e_frame]
	#puts "get_element_column $e_win $info(-row)"
	return $info(-column)
    } else {
	return -1
    }
}

proc extract_element_name {e_frame} {
    global tk_utils_defs
    
    if {[winfo exists $e_frame]} {

	set c [keylget tk_utils_defs CONTAINER.WIN]

	regexp "($c)(\[0-9\]*).f(\[0-9\]*)" $e_frame win c_win c_id e_id
	return $c_win$c_id[keylget tk_utils_defs ELEMENT.WIN]$e_id

    } else {
	return $e_frame
    }


}

proc scroll_x {container_id column args } {

    #puts "container_scroll_x $container_id $column $args"
    container_scroll_x -container_id $container_id -column $column -command $args
    

}

proc scroll_y {container_id row args } {

    #puts "container_scroll_y $container_id $column $args"
    container_scroll_y -container_id $container_id -row $row -command $args
}

proc container_position_label {label position format } {

    eval scan $position $format pos
    $label configure -text $pos
}

#draw a crosshair
proc add_container_crosshair {c_win e_id x y} {
    global $c_win.crosshair

    if {[set $c_win.crosshair]} {
	draw_container_crosshair -element_id $e_id -x $x -y $y
    } else {
	delete_container_crosshair -element_id $e_id 
    }
}

proc get_container_info {element} {
    global tk_utils_defs

    set c [keylget tk_utils_defs CONTAINER.WIN]

    regexp "($c)(\[0-9\]*)" $element win c_win c_id

    #puts "get_container_id $element $win $c_win $c_id"

    keylset info container_id $c_id
    keylset info container_win $c_win

    return $info
}

proc get_element_info {element} {
    global tk_utils_defs

    set c [keylget tk_utils_defs CONTAINER.WIN]
    set e [keylget tk_utils_defs ELEMENT.WIN]

    if {[regexp "($c)(\[0-9\]*)($e)(\[0-9\]*)" $element win c_win c_id e_win e_id] == 0} {
	return ""
    }

    #puts "get_container_id $element $win $c_win $c_id $e_win $e_id"

    keylset info container_id $c_id
    keylset info container_win $c_win
    keylset info element_id $e_id
    keylset info element_win $e_win

    return $info
}

proc set_element_specs {e_win list} {
    global $e_win.specs

    #puts "set_element_specs $list"
    foreach {label value} $list {
	keylset $e_win.specs $label $value
    }
}

proc get_element_specs {e_win} {
    global $e_win.specs

    if {[info exists $e_win.specs]} {
	return [set $e_win.specs]
    } else {
	return ""
    }
}

proc result_config {result_id line_width colour e_id e_win} {

    if {[xtoplevel [set f .cbox] -resizable 0] == ""} return

    frame $f.lw
    # Line width
    scale $f.lw.scale \
	    -label "Line width" \
	    -from 0 \
	    -to 10 \
	    -orient horiz
    $f.lw.scale set $line_width

    pack $f.lw -side top -fill both
    pack $f.lw.scale -side top -fill x

    #cmd to execute when ok button on colourbox pressed
    set e_specs [get_element_specs $e_win]
    set e_type [keylget e_specs type]

    if {$e_type == "CANVAS"} {
	set ok_cmd "configure_canvas_result $e_id $e_win $result_id \[$f.lw.scale get\]"

	#cmd to execute when changing colours on colourbox
	set update_cmd "update_canvas_result $e_id $e_win $result_id \[$f.lw.scale get\]"
    } elseif {$e_type == "RASTER"} {
	set ok_cmd "ConfigureRasterPlot $result_id \[$f.lw.scale get\]"

	#cmd to execute when changing colours on colourbox
	set update_cmd "UpdateRasterPlot $result_id"
    }


    #cmd to execute when cancel button on colourbox pressed
    set cancel_cmd "set dummy 0"

    ColourBox $f $colour $ok_cmd $update_cmd $cancel_cmd
}

proc delete_sb_h {c_win column} {
    global tk_utils_defs

    #puts "&&&&&&&&& delete_sb_h $c_win $column"

    set sb_h [keylget tk_utils_defs CONTAINER.HSCROLL.WIN]$column
    destroy $c_win$sb_h
}

proc delete_sb_v {c_win row} {
    global tk_utils_defs

    #puts "&&&&&&&&& delete_sb_v $c_win $row"
    set sb_v [keylget tk_utils_defs CONTAINER.VSCROLL.WIN]$row
    destroy $c_win$sb_v
}

#if only have rulers 
proc shuffle_rulers {ruler_win } {
    global tk_utils_defs

    set row [get_element_row $ruler_win]
    set len_row [keylget tk_utils_defs CONTAINER.RULER_ROW]
    puts "shuffle_rulers $row $len_row $ruler_win [winfo parent $ruler_win]"

    foreach win [grid slaves [winfo parent $ruler_win]] {
	#puts "win $win [grid info $win]"
	array set info [grid info $win]
	#puts "$win [winfo class $win] $info(-row)"
	#puts "[pack slaves $win]"
	
	set e_win [pack slaves $win]
	if {[string compare $e_win ""] != 0} {
	    puts "e_win $e_win"
	    set e_specs [get_element_specs $e_win]
	    if {[llength $e_specs] > 0} {
		set e_type [keylget e_specs type]
		if {[string compare $e_type CANVAS] != 0} {
		    puts "ruler $e_win"
		}
	    }
	}
    }
}

proc rotate_element {e_win ruler_win seq_id result_id orientation row} {
    global tk_utils_defs 

    set_element_specs $e_win "orientation $orientation"

    #puts ROTATE
    #set column [expr [get_current_column [winfo parent $e_win]] - 1]
    set column [seq_find_column -seq_id $seq_id]

    #find all ruler elements and move 1 to left
    #set wins [keylget tk_utils_defs ELEMENT.WIN]
    #foreach win [grid slaves [winfo parent $e_win] -column $column] {
	#puts "win $win"

	#if {[string first $wins $win] != -1} {
	    #array set info_$win [grid info $win]
	    #grid forget $win
	    #set info_${win}(-column) [expr $column - 1]
	    #eval grid $win [array get info_$win]
	#}
    #}

    #puts "rotate_element row = $row"
    move_result_new $e_win $e_win $result_id $row $column
}

proc tcl_delete_row {c_win row} {
    
    set list [grid_find $c_win row $row]

    #puts "&&&&&&&&& tcl_delete_row $row $list"
    foreach win $list {
	#need to destroy element windows packed inside packing frame
	foreach i [pack slaves $win] {
	    destroy $i
	}
	grid rowconfig $c_win $row -weight 0
	destroy $win
    }
}

proc tcl_delete_column {c_win column} {
    
    set list [grid_find $c_win column $column]

    #puts "&&&&&&&&& tcl_delete_column $column $list"
    foreach win $list {
	#need to destroy element windows packed inside packing frame
	foreach i [pack slaves $win] {
	    destroy $i
	}
	grid columnconfig $c_win $column -weight 0
	destroy $win
    }
}

proc update_container_menu {c_win c_id e_win} {
    global container_menu tk_utils_defs

    if {![catch [set e_specs [get_element_specs $e_win]]]} {
	return
    }

    if {![catch [set menu [keylget e_specs menu]]]} {
	return
    }

    global $menu
    container_menubar $c_win $c_id $menu

    result_list_update $c_win
}

proc delete_menubar {c_win} {
    global tk_utils_defs menulist

    set mb [keylget tk_utils_defs CONTAINER.MENUBAR]
    set mn $c_win$mb
    global $mn.MS

    destroy $mn

    foreach menu $menulist {
	global $menu
	if {[info exists $menu]} {
	    unset $menu
	}
    }
    set menulist {}
}

#check if column is empty except for win
proc is_empty_column {win} {

    set parent [winfo parent $win]
    set column [get_element_column $win]

    set list [grid_find $parent column $column]

    if {[llength $list] > 1} {
	return 0
    } 
    return 1
}

proc is_empty_row {win} {

    set parent [winfo parent $win]
    set row [get_element_row $win]

    set list [grid_find $parent row $row]

    if {[llength $list] > 1} {
	return 0
    } 
    return 1
}

proc tcl_delete_element {e_win} {
    
    if {[winfo exists $e_win]} {
	set einfo [get_element_info $e_win]
	if {$einfo == ""} {
	    set e_frame $e_win
	} else {
	    set e_frame [keylget einfo container_win][keylget einfo container_id].f[keylget einfo element_id]
	}
	destroy $e_win
	destroy $e_frame
    }
}   

#find the next available column when moving e_to to LEFT or RIGHT of e_from
proc find_next_column {e_to pos} {
    global LEFT RIGHT

    set found 0
    set column [get_element_column $e_to]
    set c_win [winfo parent $e_to]

    if {$pos == $LEFT} {
	incr column -1
	while {!$found} {
	    
	    set list [grid_find $c_win column $column]
	    if {[llength $list] == 0} {
		set found 1
	    } else {
		set win [lindex $list 0]
		
		set name [extract_element_name $win]
		set e_struct [element_struct [get_element_id $name]]
		if {[keylget e_struct parent] == -1} {
		    set found 1
		} else {
		    incr column -1
		}
	    }
	}
    }
    if {$pos == $RIGHT} {
	incr column 1
	while {!$found} {
	    
	    set list [grid_find $c_win column $column]
	    if {[llength $list] == 0} {
		set found 1
	    } else {
		set win [lindex $list 0]
		
		set name [extract_element_name $win]
		set e_struct [element_struct [get_element_id $name]]
		if {[keylget e_struct parent] == -1} {
		    set found 1
		} else {
		    incr column 1
		}
	    }
	}
    }
    return $column
}

