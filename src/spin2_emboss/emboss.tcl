proc plot_emboss {seq_id result_id results name } {
    global spin_defs nip_defs tk_utils_defs
    global colour col_index NUM_COLOURS HORIZONTAL SCALE_BAR TOP_S

    set type [keylget spin_defs EMBOSS.TYPE]
    set frame 0
    set strand $TOP_S
    
    set element_info [create_seq_element $seq_id -1 $type $strand $frame \
	    $HORIZONTAL XY "CANVAS" [keylget spin_defs CONTAINER.TITLE] \
	    [keylget tk_utils_defs ELEMENT.PLOT_WIDTH] \
 	    [keylget tk_utils_defs ELEMENT.PLOT_HEIGHT]]

    set c_win [keylget element_info container_win]
    set c_id [keylget element_info container_id]
    set e_win [keylget element_info element_win]
    set e_id [keylget element_info element_id]
    set row [keylget element_info row]

    update idletasks

    emboss plot -element $c_win$e_win\
	    -container $c_win\
	    -seq_id_h $seq_id \
	    -result_id $result_id\
	    -results $results\
	    -container_id [keylget element_info container_id]\
	    -element_id [keylget element_info element_id]\
	    -element_type "CANVAS"\
	    -graph 0\
	    -name $name\
	    -fill red \
	    -width [keylget spin_defs EMBOSS.L_WIDTH]
    
    seqed_element_bindings $c_id $c_win$e_win $e_id

    #update result list
    result_list_update $c_win 
}

proc plot_emboss_dot {seq_id_h seq_id_v result_id results name} {
    global spin_defs nip_defs tk_utils_defs
    global colour col_index NUM_COLOURS HORIZONTAL SCALE_BAR TOP_S BOTH

    set type [keylget spin_defs EMBOSS.TYPE]
    set frame 0
    set strand $TOP_S

    set element_info [create_seq_element $seq_id_h $seq_id_v $type $strand $frame \
	    $BOTH XY "CANVAS" [keylget spin_defs CONTAINER.TITLE] \
	    [keylget tk_utils_defs DOT.PLOT_WIDTH] \
 	    [keylget tk_utils_defs DOT.PLOT_HEIGHT]]

    set c_win [keylget element_info container_win]
    set c_id [keylget element_info container_id]
    set e_win [keylget element_info element_win]
    set e_id [keylget element_info element_id]
    set row [keylget element_info row]

    update idletasks

    emboss plot -element $c_win$e_win\
	    -container $c_win\
	    -seq_id_h $seq_id_h \
	    -seq_id_v $seq_id_v \
	    -result_id $result_id\
	    -results $results\
	    -container_id [keylget element_info container_id]\
	    -element_id [keylget element_info element_id]\
	    -element_type "CANVAS"\
	    -graph 1\
	    -name $name\
	    -fill red \
	    -width [keylget spin_defs EMBOSS.L_WIDTH]
    
    seqed_element_bindings $c_id $c_win$e_win $e_id

    #update result list
    result_list_update $c_win 
}

proc calc_zoom_origin {min1 min2 max1 max2 } {

    if {[expr ($max2 + $min1 - $max1 - $min2)] == 0.0} {
	return [expr ($min1 * $max2) - ($max1 * $min2) / 0.1]
    } else {
	return [expr ((($min1 * $max2) - ($max1 * $min2)) / ($max2 + $min1 - $max1 - $min2))]
    }

}

proc calc_zoom_sf {min1 min2 max1 max2} {
    if {[expr $max2 - $min2] == 0.0} {
	return 1.0
    } else {
	return [expr (($max1 - $min1) / ($max2 - $min2))]
    }
}


proc plot_emboss_graphic {program fs} {
    global e_incr

    if {![info exist e_incr]} {
	set e_incr 0
    }

    array set emboss_colour {0 "BLACK" 1 "RED" 2 "YELLOW" 3 "GREEN" \
	    4 "AQUAMARINE" 5 "PINK" 6 "WHEAT" 7 "GREY" 8 "BROWN" 9 "BLUE" \
	    10 "BLUEVIOLET" 11 "CYAN" 12 "TURQUOISE" 13 "MAGENTA" 14 "SALMON" \
	    15 "WHITE"}

    set t .{$program}_g$e_incr
    xtoplevel $t
    wm title $t "graphic $e_incr"

    set width 700
    set height 700
    canvasbox $t.c -bd 0 -highlightthickness 0 -width $width -height $height
    pack $t.c

    set scale [gets $fs]
    #puts "scale $scale"
    scan $scale "##Screen x1 %f y1 %f x2 %f y2%f" s_x1 s_y1 s_x2 s_y2

    while {![eof $fs]} {
	set line [gets $fs]
	if {[regexp {^Line} $line]} {
	    scan $line "Line x1 %f y1 %f x2 %f y2 %f colour %d" x1 y1 x2 y2 colour
	    $t.c create line $x1 [expr $s_y2 - $y1 + $s_y1] $x2 [expr $s_y2 - $y2 + $s_y1] -fill $emboss_colour($colour) -width 0 -capstyle round
	    
	}
	if {[regexp {^Text1} $line]} {
	    scan $line "Text1 x1 %f y1 %f colour %d size %f %s" x1 y1 colour size text
	    #text might contain white spaces
	    regexp {size [0-9.]* (.*)} $line dummy text
	    #puts "text1=$text"
	    $t.c create text $x1 [expr $s_y2 - $y1 + $s_y1] -text $text -fill $emboss_colour($colour) -anchor w
	}
	if {[regexp {^Text2} $line]} {
	    scan $line "Text2 x1 %f y1 %f colour %d size %f %s" x1 y1 colour size text
	    #text might contain white spaces
	    regexp {size [0-9.]* (.*)} $line dummy text
	    #puts "text2=$text"
	    $t.c create text $x1 [expr $s_y2 - $y1 + $s_y1] -text $text -fill $emboss_colour($colour)
	}
	if {[regexp {^Text3} $line]} {
	    scan $line "Text3 x1 %f y1 %f colour %d size %f %s" x1 y1 colour size text
	    #text might contain white spaces
	    regexp {size [0-9.]* (.*)} $line dummy text
	    #puts "text3=$text"
	    $t.c create text $x1 [expr $s_y2 - $y1 + $s_y1] -text $text -fill $emboss_colour($colour) -anchor e
	}
	if {[regexp {^Textline} $line]} {
	    scan $line "Textline x1 %f y1 %f x2 %f y2 %f colour %d size %f %s" x1 y1 x2 y2 colour size text
	    #text might contain white spaces
	    regexp {size [0-9.]* (.*)} $line dummy text
	    $t.c create text $x1 [expr $s_y2 - $y1 + $s_y1] -text $text -fill $emboss_colour($colour)
	}
	if {[regexp {^Rectangle} $line]} {
	    scan $line "Rectangle x1 %f y1 %f x2 %f y2 %f colour %d" x1 y1 x2 y2 colour
	    $t.c create rectangle $x1 [expr $s_y2 - $y1 + $s_y1] $x2 [expr $s_y2 - $y2 + $s_y1] -outline $emboss_colour($colour)
	}
	if {[regexp {^Shaded Rectangle} $line]} {
	    scan $line "Shaded Rectangle x1 %f y1 %f x2 %f y2 %f colour %d" x1 y1 x2 y2 colour
	    $t.c create rectangle $x1 [expr $s_y2 - $y1 + $s_y1] $x2 [expr $s_y2 - $y2 + $s_y1] -fill $emboss_colour($colour)
	}
    }

    set x_origin [calc_zoom_origin $s_x1 0 $s_x2 $width] 
    set y_origin [calc_zoom_origin $s_y1 0 $s_y2 $height]
    set sf_x [calc_zoom_sf 0 $s_x1 $width $s_x2]
    set sf_y [calc_zoom_sf 0 $s_y1 $height $s_y2]

    #puts "SCALE $x_origin $y_origin $sf_x $sf_y"

    $t.c scale all $x_origin $y_origin $sf_x $sf_y
    incr e_incr
}
