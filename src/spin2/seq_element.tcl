proc create_seq_element {seq_id_h seq_id_v plot_type strand frame orientation scale element_type title width height {container_id -1}} {
    global HORIZONTAL tk_utils_defs TOP_S LEFT BOTTOM

    set element_info [get_element_name -seq_id_h $seq_id_h -seq_id_v $seq_id_v -plot_type $plot_type -frame $frame -container_id $container_id]

    set c_win [keylget element_info container_win]
    set c_id [keylget element_info container_id]

    #puts "*************create_seq_element c_win $c_win $c_id container_id $container_id"

    if {![winfo exists $c_win]} {
	create_container $c_win [keylget element_info container_id] $title
    } else {
	raise $c_win
	focus $c_win
    }

    set e_win [keylget element_info element_win]
    set e_id [keylget element_info element_id]
    set e_win_to [keylget element_info e_win_to]

    set new_orientation [keylget element_info orientation]

    #if orientation is changed, swap seq_id_h & seq_id_v, width & height
    if {$new_orientation != $orientation} {
	set tmp $seq_id_v 
	set seq_id_v $seq_id_h
	set seq_id_h $tmp
	set tmp $width
	set width $height
	set height $tmp
    }

    #puts "*************create_seq_element e_win  $e_win $e_id"

    if {![winfo exists $c_win$e_win]} {
	set_element_specs $c_win$e_win "strand $strand orientation $orientation scale $scale type $element_type title $title"

	#puts "set_element_specs strand $strand orientation $new_orientation scale $scale type $element_type title $title"

	set list [get_element_position $c_win $e_win_to $strand $new_orientation]

	set row [lindex $list 0]
	set column [lindex $list 1]

	#puts "ROW $row COLUMN $column"
	if {[string compare $element_type "CANVAS"] == 0} {
	    create_canvas_element $c_win $c_id $c_win$e_win $e_id $width $height $row $column $new_orientation $scale NEW_ROW
	}
	keylset element_info row $row
	keylset element_info column $column 
    } else {
	keylset element_info row [get_element_row $c_win$e_win]
	keylset element_info column [get_element_column $c_win$e_win]
    }

    return $element_info
}

proc element_create_seqed {c_id e_win e_id px py} {
    #puts element_create_seqed 

    #only want to create a new seqed if one doesn't already exist
    set result [seq_find_cursor -element_id $e_id]

    set cursor_id [lindex $result 0]
    set seq_id [lindex $result 1]
    
    set result [pixel_to_world -element_id $e_id -x [$e_win canvasx $px] -y [$e_win canvasy $py]]
    set wx [lindex $result 0]
    set wy [lindex $result 1]

    if {$cursor_id < 0} {
	#no cursors 
	set row 1
	#puts "POS $wx"

	#SeqedDisplayPacked [expr int($wx)] $seq_id [winfo parent $e_win] $row [get_element_column $e_win]
	SeqedDisplay $c_id [expr int($wx)] $seq_id

    } else {
	canvas_move_cursor -element_id $e_id -pos [$e_win canvasx $px] -cursor_id $cursor_id

	#puts "element_create_seqed [$e_win canvasx $px]"
    }
}

proc seqed_element_bindings {c_id e_win e_id } {

    # Double button 1 and 2 to move or create an editor
    bind $e_win <<move-create>> "element_create_seqed $c_id $e_win $e_id %x %y"
    bind $e_win <<use>> "+element_create_seqed $c_id $e_win $e_id %x %y"
}


proc element_create_seqpair {e_win e_id px py} {
    global HORIZONTAL VERTICAL

    #only want to create a new seqpair if one doesn't already exist
    set result [seq_find_cursor -element_id $e_id -direction $HORIZONTAL]
    set cursor_id_h [lindex $result 0]
    set seq_id_h [lindex $result 1]
    
    #only want to create a new seqpair if one doesn't already exist
    set result [seq_find_cursor -element_id $e_id -direction $VERTICAL]
    set cursor_id_v [lindex $result 0]
    set seq_id_v [lindex $result 1]

    set result [pixel_to_world -element_id $e_id -x [$e_win canvasx $px] -y [$e_win canvasy $py]]
    set wx [lindex $result 0]
    set wy [lindex $result 1]

    set wy [invert_wy -element_id $e_id -y $wy]

    if {$cursor_id_h < 0 && $cursor_id_v < 0} {
	#no cursors 
	#puts "canvas_move_cursor_h pos [$e_win canvasx $px] [expr int($wx)]"
	
	set result_id [seq_find_result_id -element_id $e_id -seq_id_h $seq_id_h -seq_id_v $seq_id_v]

	SequencePairDisplay $wx $wy $seq_id_h $seq_id_v $cursor_id_h $cursor_id_v $result_id
	
    } else {
	#puts "element_create_seqpair $wx $px [$e_win canvasx $px]"

	canvas_move_cursor -element_id $e_id -pos [$e_win canvasx $px] -cursor_id $cursor_id_h -direction $HORIZONTAL
	canvas_move_cursor -element_id $e_id -pos [$e_win canvasy $py] -cursor_id $cursor_id_v -direction $VERTICAL
    }
}

proc seqpair_element_bindings {e_win e_id } {

    # Double button 1 and 2 to move or create an editor
    bind $e_win <<move-create>> "element_create_seqpair $e_win $e_id %x %y"
    bind $e_win <<use>> "+element_create_seqpair $e_win $e_id %x %y"
}

proc seq_canvas_move_cursor {e_win e_id cursor_id direction x y} {
    global HORIZONTAL moving_cursor moving_result

    #puts seq_canvas_move_cursor

    if {[info exists moving_result]} {
	return
    }

    set moving_cursor 1

    if {$direction == $HORIZONTAL} {

	#puts "seq_canvas_move_cursor $e_win x $x [$e_win canvasx $x]"
	canvas_move_cursor -element_id $e_id -pos [$e_win canvasx $x] -cursor_id $cursor_id -direction $direction
    } else {
	canvas_move_cursor -element_id $e_id -pos [$e_win canvasy $y] -cursor_id $cursor_id -direction $direction

    }
}

proc unset_moving_cursor { } {
    global moving_cursor moving_result

    if {[info exists moving_result]} {
	return
    }

    if {[info exists moving_cursor]} {
	unset moving_cursor
    }
}

# Create a new visible cursor
proc seq_canvas_cursor_create {seq_id e_win cursor_id e_id direction colour} {
    global HORIZONTAL
    
    #puts seq_canvas_cursor_create

    if {$direction == $HORIZONTAL} {
	set height [$e_win canvasy [winfo height $e_win]]
	
	#HACK - why 10?
	set line   [$e_win create line 10 [$e_win canvasy 0] 10 $height \
		-fill $colour \
		-tag "cursor_$cursor_id cursor"]
    } else {
	set width [$e_win canvasx [winfo width $e_win]]
	
	#HACK - why 10?
	set line   [$e_win create line [$e_win canvasx 0] 10 $width 10\
		-fill $colour \
		-tag "cursor_$cursor_id cursor"]

    }
    #bind $e_win <<move-drag>> "canvas_move_cursor -element_id $e_id -pos \[$e_win canvasx %x\] -cursor_id $cursor_id -direction $direction"

    #FIXME - not sure this is the effect we want - should we click on the
    #cursor to select it?

    bind $e_win <<move-drag>> "+seq_canvas_move_cursor $e_win $e_id $cursor_id $direction %x %y"

    bind $e_win <<move-release>> "+unset_moving_cursor"
}

# Move an existing cursor to 'x'; called from C
proc seq_canvas_cursor_move_x {seq_id e_win cursor_id e_id colour x wx} {
    global tk_utils_defs HORIZONTAL

   #puts "seq_canvas_cursor_move_x $seq_id $e_win $cursor_id $e_id $colour $x $wx"
    if {![winfo exists $e_win]} {
	return
    }

    # Create the cursor if it doesn't currently exist
    if {$seq_id != -1} {
	if {[$e_win find withtag cursor_$cursor_id] == ""} {
	    seq_canvas_cursor_create $seq_id $e_win $cursor_id $e_id $HORIZONTAL $colour
	}
    }
    # And now move it
    set height [$e_win canvasy [winfo height $e_win]]

    set y1 [lindex [$e_win cget -scrollregion] 1]
    set y2 [lindex [$e_win cget -scrollregion] 3]

    #no scrolling allowed eg ruler
    if {$y2 == 0} {
	set y2 [winfo height $e_win]
    }

    $e_win coords cursor_$cursor_id $x $y1 $x $y2
    #puts "seq_canvas_cursor_move_x x $x wx $wx"

    set buttons [keylget tk_utils_defs CONTAINER.BUTTONS.WIN]
    set pos1 [keylget tk_utils_defs CONTAINER.POS1.NAME]
    set x_format %d
    set c_win [winfo parent $e_win]
    container_position_label $c_win$buttons$pos1 $wx $x_format
}

# Move an existing cursor to 'y'; called from C
proc seq_canvas_cursor_move_y {seq_id e_win cursor_id e_id colour y wy} {
    global tk_utils_defs VERTICAL

    #puts "seq_canvas_cursor_move_y $seq_id $e_win $cursor_id $e_id $colour $y $wy"

    # Create the cursor if it doesn't currently exist
    if {$seq_id != -1} {
	if {[$e_win find withtag cursor_$cursor_id] == ""} {
	    seq_canvas_cursor_create $seq_id $e_win $cursor_id $e_id $VERTICAL $colour
	}
    }
    # And now move it
    set width [$e_win canvasx [winfo width $e_win]]

    set x1 [lindex [$e_win cget -scrollregion] 0]
    set x2 [lindex [$e_win cget -scrollregion] 2]

    #no scrolling allowed eg ruler
    if {$x2 == 0} {
	set x2 [winfo width $e_win]
    }

    $e_win coords cursor_$cursor_id $x1 $y $x2 $y

    set buttons [keylget tk_utils_defs CONTAINER.BUTTONS.WIN]
    set pos2 [keylget tk_utils_defs CONTAINER.POS2.NAME]
    set y_format %d
    set c_win [winfo parent $e_win]
    container_position_label $c_win$buttons$pos2 $wy $y_format

}



