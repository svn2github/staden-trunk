#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

##############################################################################
proc SequencePairDisplay {wx wy seq_id_h seq_id_v cursor_id_h cursor_id_v result_id} {
    global spin_defs 

    #puts "SequencePairDisplay ID $result_id $seq_id_h $seq_id_v"
    #FIXME - more than 1 seq disp per result?
    set seq_win [keylget spin_defs SEQ_DISP.WIN]$result_id
    if {[winfo exists $seq_win]} {
	wm deiconify $seq_win
	raise $seq_win
	return
    }

    global $seq_win.seqdisp_id

    CreateSequenceDisplay $seq_win $result_id $wx $wy
    set $seq_win.seqdisp_id [seq_pair_display -window $seq_win \
				 -result_id $result_id \
				 -seq_id_h $seq_id_h -seq_id_v $seq_id_v \
				 -cursor_id_h $cursor_id_h \
				 -cursor_id_v $cursor_id_v \
				 -x $wx -y $wy]
}

##############################################################################
proc seq_replot { f diff seq1 seq2 left1 left2 result_index} {
    global $f.width $f.p_tag1 $f.p_tag2

    set middle [expr [set $f.width] / 2]

    $diff config -state normal
    $seq1 config -state normal
    $seq2 config -state normal
    
    $diff delete 0.0 end

    $seq1 delete 1.0 end
    $seq2 delete 1.0 end
    $seq1 delete 2.0 end
    $seq2 delete 2.0 end
    $seq1 insert 1.end \n
    $seq1 insert 2.end \n
    $seq2 insert 1.end \n
    $seq2 insert 2.end \n
    
    #puts "seq_replot $left1 $left2"
    update_seq_pair -win_diff $diff -win_1 $seq1 -win_2 $seq2 -left1 $left1 -left2 $left2 -win_len [set $f.width] -result_index $result_index

    $seq1 tag remove cursor1 2.[set $f.p_tag1]
    $seq2 tag remove cursor2 1.[set $f.p_tag2]

    $seq1 tag add cursor1 2.$middle
    $seq2 tag add cursor2 1.$middle
    $seq1 tag config cursor1 -background black -foreground white
    $seq2 tag config cursor2 -background black -foreground white

    $seq1 config -state disabled
    $seq2 config -state disabled
    $diff config -state disabled

    set $f.p_tag1 $middle
    set $f.p_tag2 $middle

}

##############################################################################
proc seq_scroll {f seq1 seq2 number diff result_index seq_dir} {
    global $f.lock
    
    if {[set $f.lock]} {
	seq_scroll12 $f $seq1 $seq2 $number $diff $result_index
    } elseif {$seq_dir == 1} {
	seq_scroll1 $f $seq1 $seq2 $number $diff $result_index
    } elseif {$seq_dir == 2} {
	seq_scroll2 $f $seq1 $seq2 $number $diff $result_index
    }

}

##############################################################################
#start1 is the window central BASE in seq1 ie tagged
#left1 is the left most BASE in seq1
#p1 is the previous start1
proc seq_scroll1 {f seq1 seq2 number diff result_index} {
    global $f.p1 $f.p2 $f.max1 $f.max2 $f.width $f.lock $f.seqdisp_id
    global HORIZONTAL VERTICAL   
    
    if {[expr $number + [set $f.p1]] > [set $f.max1]} {
	set number [expr [set $f.max1] - [set $f.p1]]
    }

    set middle [expr [set $f.width] / 2]
    set start1 [expr [set $f.p1] + $number]
    if {$start1 < 1} {set start1 1}
    #if {$start1 > [set $f.max1]} {set start1 [set $f.max1]}
    set left1 [expr $start1 - $middle]
    set left2 [expr [set $f.p2] - $middle]
 
    #puts "seq_scroll1 $left1 [set $f.p1]"
    seq_replot $f $diff $seq1 $seq2 $left1 $left2 $result_index
    set $f.p1 $start1

    #also need to send callback to raster cursors since scrolling effectively
    #moves the cursor
    seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $HORIZONTAL -pos [set $f.p1]
    seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $VERTICAL -pos [set $f.p2]
}

##############################################################################
proc seq_scroll2 {f seq1 seq2 number diff result_index} {
    global $f.p1 $f.p2 $f.max1 $f.max2 $f.width $f.lock $f.seqdisp_id
    global HORIZONTAL VERTICAL
    
    #puts seq_scroll2

    if {[expr $number + [set $f.p2]] > [set $f.max2]} {
	set number [expr [set $f.max2] - [set $f.p2]]
    }

    set middle [expr [set $f.width] / 2]
    set start2 [expr [set $f.p2] + $number]
    if {$start2 < 1} {set start2 1}

    set left2 [expr $start2 - $middle]
    set left1 [expr [set $f.p1] - $middle]

    seq_replot $f $diff $seq1 $seq2 $left1 $left2 $result_index
    set $f.p2 $start2

    #also need to send callback to raster cursors since scrolling effectively
    #moves the cursor
    seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $HORIZONTAL -pos [set $f.p1]
    seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $VERTICAL -pos [set $f.p2]
}

##############################################################################
proc seq_scroll12 {f seq1 seq2 number diff result_index} {
    global $f.p1 $f.p2 $f.max1 $f.max2 $f.width $f.lock $f.seqdisp_id
    global HORIZONTAL VERTICAL 
    
    #puts seq_scroll12

    set middle [expr [set $f.width] / 2]

    set start1 [expr [set $f.p1] + $number]
    set start2 [expr [set $f.p2] + $number]
    #if {$start1 < 1} {set start1 1}
    #if {$start2 < 1} {set start2 1}

    #puts "start1 $start1 start2 $start2 max1 [set $f.max1] max2 [set $f.max2] p1 [set $f.p1] p2 [set $f.p2]"

   if { ($start1 < 1) && ($start2 < 1)} {
	return
    }

    if { ($start2 > [set $f.max2]) && ($start1 > [set $f.max1])} {
	return
    }

    set left1 [expr $start1 - $middle]
    set left2 [expr $start2 - $middle]

    #puts "seq_scroll12 $left1 [set $f.p1] $left2 [set $f.p2]"
    seq_replot $f $diff $seq1 $seq2 $left1 $left2 $result_index
    set $f.p1 $start1
    set $f.p2 $start2

    #also need to send callback to raster cursors since scrolling effectively
    #moves the cursor
    seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $HORIZONTAL -pos [set $f.p1]
    seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $VERTICAL -pos [set $f.p2]
}

##############################################################################
proc SeqPairRescaleSeqs {f seq1 seq2 diff result_index char_width} {
    global $f.p1 $f.p2 $f.max1 $f.width

    #puts "cursor pos [set $f.p1] [set $f.p2]"
    #puts "[winfo width $seq1] $char_width"
    set $f.width [expr round(ceil([winfo width $seq1] / $char_width))]
    set middle [expr [set $f.width] / 2]
    set left1 [expr [set $f.p1] - $middle]
    set left2 [expr [set $f.p2] - $middle]
 
    #puts "rescale width [set $f.width]"
    seq_replot $f $diff $seq1 $seq2 $left1 $left2 $result_index
}

##############################################################################
proc FindNearestMatch {f result_index diff seq1 seq2 match} {
    global $f.p1 $f.p2 $f.max1 $f.max2 $f.width $f.seqdisp_id
    global HORIZONTAL VERTICAL

    set nearest [nearest_match -x [set $f.p1] -y [set $f.p2] -result_index $result_index -match $match]
    #puts "NEAREST $nearest"
    
    set start1 [lindex $nearest 0]
    set start2 [lindex $nearest 1]
    set middle [expr [set $f.width] / 2]
    if {$start1 < 1} {set start1 1}
    if {$start1 > [set $f.max1]} {set start1 [set $f.max1]}
    if {$start2 < 1} {set start2 1}
    if {$start2 > [set $f.max2]} {set start2 [set $f.max2]}
    set left1 [expr $start1 - $middle]
    set left2 [expr $start2 - $middle]

    seq_replot $f $diff $seq1 $seq2 $left1 $left2 $result_index

    set $f.p1 $start1
    set $f.p2 $start2

    seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $HORIZONTAL -pos [set $f.p1]
    seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $VERTICAL -pos [set $f.p2]
    
}

##############################################################################
proc scrollbar_view1 {f sb result_index args} {
    global $f.width $f.p1 $f.prev_pos1 spin_defs HORIZONTAL VERTICAL

    
    #puts "scrollbar_view1 $args"
    set length [seq_result_info -index $result_index -direction $HORIZONTAL -option length]
    if {![info exists $f.width]} {
	set width [keylget spin_defs SEQ_DISP.PLOT_WIDTH]
	set end [expr (double($width)/$length)]
	$sb set 0 $end
	return
    }
    set length [expr $length + [set $f.width]]

    set middle [expr double([set $f.width]) / 2]
    set left1 [expr double([set $f.p1])/$length]
    set right [expr $left1 + (double([set $f.width])/$length)]

    #puts "p1 [set $f.p1]  left [expr [set $f.p1] - $middle]"
    #puts "left1 $left1 right $right width [set $f.width]"

    $sb set $left1 $right
    set $f.prev_pos1 [set $f.p1]
}

##############################################################################
proc scrollbar_view2 {f sb result_index args} {
    global $f.width $f.p2 $f.prev_pos2 spin_defs VERTICAL

    
    #puts "scrollbar_view2 $args"
    set length [seq_result_info -index $result_index -direction $VERTICAL -option length]

    if {![info exists $f.width]} {
	set width [keylget spin_defs SEQ_DISP.PLOT_WIDTH]
	set end [expr (double($width)/$length)]
	$sb set 0 $end
	return
    }
    set length [expr $length + [set $f.width]]

    set middle [expr double([set $f.width]) / 2]
    set left1 [expr double([set $f.p2])/$length]
    set right [expr $left1 + (double([set $f.width])/$length)]

    #puts "p2 [set $f.p2]  left [expr [set $f.p2] - $middle]"
    #puts "left1 $left1 right $right width [set $f.width]"

    $sb set $left1 $right
    set $f.prev_pos2 [set $f.p2]
}

##############################################################################
proc scrollbar_move1 {f seq1 seq2 diff result_index sb args } {
    global $f.lock $f.prev_pos1 $f.p1 $f.p2 $f.width $f.seqdisp_id HORIZONTAL VERTICAL

    set len [seq_result_info -index $result_index -direction $HORIZONTAL -option length] 

    if {[llength $args] == 2} {
	#moveto fraction
	set num [expr ($len + [set $f.width]) * [lindex $args 1]]
	if {![info exists $f.prev_pos1]} {
	    set $f.prev_pos1 0
	}
	set number [expr $num - [set $f.prev_pos1]]
    } elseif {[llength $args] == 3} {
	#scroll number unit/page
	set number [lindex $args 1]
	if {[string compare [lindex $args 2] pages] == 0} {
	    set number [expr 0.9 * $number * [set $f.width]]
	}
    }

    set number [expr round($number)]
    #puts "number $number"

    seq_scroll $f $seq1 $seq2 $number $diff $result_index 1

    #also need to send callback to raster cursors since scrolling effectively
    #moves the cursor
    seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $HORIZONTAL -pos [set $f.p1]
    if {[set $f.lock]} {
	seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $VERTICAL -pos [set $f.p2]
    }

}

##############################################################################
proc scrollbar_move2 {f seq1 seq2 diff result_index sb args } {
    global $f.lock $f.prev_pos2 $f.p1 $f.p2 $f.width $f.seqdisp_id HORIZONTAL VERTICAL
    
    set len [seq_result_info -index $result_index -direction $VERTICAL -option length] 

    if {[llength $args] == 2} {
	#moveto fraction
	set num [expr ($len + [set $f.width]) * [lindex $args 1]]
	if {![info exists $f.prev_pos2]} {
	    set $f.prev_pos2 0
	}
	set number [expr $num - [set $f.prev_pos2]]
    } elseif {[llength $args] == 3} {
	#scroll number unit
	set number [lindex $args 1]
	if {[string compare [lindex $args 2] pages] == 0} {
	    set number [expr 0.9 * $number * [set $f.width]]
	}
    }
    set number [expr round($number)]
    #puts "number $number"

    seq_scroll $f $seq1 $seq2 $number $diff $result_index 2
 
    #also need to send callback to raster cursors since scrolling effectively
    #moves the cursor
    if {[set $f.lock]} {
	seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $HORIZONTAL -pos [set $f.p1]
    }
    seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $VERTICAL -pos [set $f.p2]
}

##############################################################################
#called from C (result_index IS result_index - no conversion necessary)
proc CreateSequenceDisplay { f result_index pos1 pos2} {
    global $f.p_tag1 $f.p_tag2 $f.p1 $f.p2 spin_defs HORIZONTAL
    global VERTICAL
    global $f.max1 $f.max2
    global $f.seq_disp $f.width
    global $f.lock

    xtoplevel $f
    wm title $f "Sequence Comparison Display"
    fix_maxsize $f

    #puts "start CreateSequenceDisplay $f"
    set width [keylget spin_defs SEQ_DISP.PLOT_WIDTH]
    set middle [expr $width / 2]
    set font sheet_font

    text [set diff $f.diff] -height 1 -width $width -wrap none -font $font

    text [set seq1 $f.seq1] -height 2 -width $width -wrap none \
	    -xscrollcommand "scrollbar_view1 $f $f.sb1 $result_index" \
	    -font $font
    text [set seq2 $f.seq2] -height 2 -width $width -wrap none \
	    -xscrollcommand "scrollbar_view2 $f $f.sb2 $result_index" \
	    -font $font
    scrollbar $f.sb1 -orient horizontal \
	    -command "scrollbar_move1 $f $seq1 $seq2 $diff $result_index $f.sb1"
    scrollbar $f.sb2 -orient horizontal \
	    -command "scrollbar_move2 $f $seq1 $seq2 $diff $result_index $f.sb2"
    label [set label1 $f.label1] -width 16 -font title_font
    label [set label2 $f.label2] -width 16 -font title_font

    frame $f.menubar -bd 2 -relief raised
    #set up File menu
    menubutton $f.menubar.file -text "File" -menu $f.menubar.file.opts
    set m [menu $f.menubar.file.opts]
    $m add command -label "Exit" \
	    -command "SeqDispStartShutdown $f"

    # Trap window quit so that clean-up is done properly
    wm protocol $f WM_DELETE_WINDOW "SeqDispStartShutdown $f"

    #nearest match button - not applicable for temporary results
    set type [seq_result_info -index $result_index -option type]
    if {$type == 0} {
	#permenant
	button $f.menubar.nearest -text "Nearest match" \
		-command "FindNearestMatch $f $result_index $diff $seq1 $seq2 1"\
		-state normal

	button $f.menubar.nearest_d -text "Nearest dot" \
		-command "FindNearestMatch $f $result_index $diff $seq1 $seq2 0"\
		-state normal
    } else {
	#temporary
	button $f.menubar.nearest -text "Nearest match" \
		-command "FindNearestMatch $f $result_index $diff $seq1 $seq2 1"\
		-state disabled
	button $f.menubar.nearest_d -text "Nearest dot" \
		-command "FindNearestMatch $f $result_index $diff $seq1 $seq2 0"\
		-state disabled
    }

    #lock button
    checkbutton $f.menubar.lock -text "Lock"\
	    -command "" -variable $f.lock -bd 2 -relief raised 

    #help button
    button $f.menubar.help -text "Help" \
	    -command "show_help spin {SPIN-Sequence-Comparison Display}"

    pack $f.menubar.file $f.menubar.nearest $f.menubar.nearest_d $f.menubar.lock -side left -fill both
    pack $f.menubar.help -side right
    #pack $f.menubar -side top -fill x

    $seq1 insert 1.end \n
    $seq1 insert 2.end \n
    $seq2 insert 1.end \n
    $seq2 insert 2.end \n

    $label1 configure -anchor w -text [seq_result_info -index $result_index -direction $HORIZONTAL -option basename]
    $label2 configure -anchor w -text [seq_result_info -index $result_index -direction $VERTICAL -option basename]

    #starting sequence positions
    #set left [expr -1 * $middle + 1]
    set left1 [expr $pos1 - $middle]
    set left2 [expr $pos2 - $middle]

    update_seq_pair -win_diff $diff -win_1 $seq1 -win_2 $seq2 -left1 $left1 -left2 $left2 -win_len $width -result_index $result_index

    set $f.max1 [seq_result_info -index $result_index -direction $HORIZONTAL -option length] 
    set $f.max2 [seq_result_info -index $result_index -direction $VERTICAL -option length] 

    $diff config -state disabled
    $seq1 config -state normal
    $seq2 config -state disabled


    repeater $f.ll1 "seq_scroll $f $seq1 $seq2 [expr -1 * $width] $diff $result_index 1" -text << -pady 0
    repeater $f.l1  "seq_scroll $f $seq1 $seq2 -1 $diff $result_index 1" -text < -pady 0
    repeater $f.r1  "seq_scroll $f $seq1 $seq2 1 $diff $result_index 1" -text > -pady 0
    repeater $f.rr1 "seq_scroll $f $seq1 $seq2 $width $diff $result_index 1" -text >> -pady 0 

    repeater $f.ll2 "seq_scroll $f $seq1 $seq2 [expr -1 * $width] $diff $result_index 2" -text << -pady 0
    repeater $f.l2  "seq_scroll $f $seq1 $seq2 -1 $diff $result_index 2" -text < -pady 0
    repeater $f.r2  "seq_scroll $f $seq1 $seq2 1 $diff $result_index 2" -text > -pady 0
    repeater $f.rr2 "seq_scroll $f $seq1 $seq2 $width $diff $result_index 2" -text >> -pady 0

    set $f.p_tag1 $middle
    set $f.p_tag2 $middle
    set $f.p1 $pos1
    set $f.p2 $pos2
    $seq1 tag add cursor1 2.$middle
    $seq2 tag add cursor2 1.$middle
    $seq1 tag config cursor1 -background black -foreground white
    $seq2 tag config cursor2 -background black -foreground white
    
    grid columnconfig $f 4 -weight 1

    grid $f.menubar -row 0 -column 0 -sticky ew -columnspan 5
    grid $f.ll1 -row 1 -column 0 -sticky ew
    grid $f.l1  -row 1 -column 1 -sticky ew
    grid $f.r1  -row 1 -column 2 -sticky ew
    grid $f.rr1 -row 1 -column 3 -sticky ew

    grid $f.sb1 -row 1 -column 4 -sticky ew
    grid $f.label1 -row 2 -column 0 -sticky w -columnspan 4

    grid $f.seq1 -row 2 -column 4 -sticky ew
    grid $f.diff -row 3 -column 4 -sticky ew
    grid $f.seq2 -row 4 -column 4 -sticky ew
    
    grid $f.label2 -row 4 -column 0 -sticky w -columnspan 4
    grid $f.sb2 -row 5 -column 4 -sticky ew

    grid $f.ll2 -row 5 -column 0 -sticky ew
    grid $f.l2  -row 5 -column 1 -sticky ew
    grid $f.r2  -row 5 -column 2 -sticky ew
    grid $f.rr2 -row 5 -column 3 -sticky ew

    update

    set char_width [font measure [$seq1 cget -font] 0]
    set $f.width [expr round(ceil([winfo width $seq1] / $char_width))]

    bind $f <Configure> "SeqPairRescaleSeqs $f $seq1 $seq2 $diff $result_index $char_width"

    bind $seq1 <Any-Enter> {focus %W}
    bind $seq1 <<select>> "SeqPairCursorSet1 $f $seq1 $seq2 $diff $result_index %x %y"
    bind $seq1 <Key-Left> "SeqPairCursor $f $seq1 $seq2 -1 $diff $result_index 1"
    bind $seq1 <Key-Right> "SeqPairCursor $f $seq1 $seq2 1 $diff $result_index 1"

    bind $seq2 <Any-Enter> {focus %W}
    bind $seq2 <<select>> "SeqPairCursorSet2 $f $seq1 $seq2 $diff $result_index %x %y"
    bind $seq2 <Key-Left> "SeqPairCursor $f $seq1 $seq2 -1 $diff $result_index 2"
    bind $seq2 <Key-Right> "SeqPairCursor $f $seq1 $seq2 1 $diff $result_index 2"

    #set only my bindings ie not text widget
    bindtags $seq1 "$seq1 all"
    bindtags $seq2 "$seq2 all"
}

##############################################################################
#position cursor with mouse on horizontal seq
proc SeqPairCursorSet1 {f seq1 seq2 diff result_index x y} {
    global $f.p1 HORIZONTAL $f.seqdisp_id $f.width
    
    #FIXME - better way of getting char from line.char format 
    set n [$seq1 index @$x,$y]
    regexp {\.[0-9]+$} $n char
    set char [string range $char 1 [string length $char]]
    #puts "char $char"
    set middle [expr [set $f.width] / 2]

    #puts "SeqPairCursorSet1 [set $f.p1]"
    set $f.p1 [expr [set $f.p1] + $char - $middle]
    seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $HORIZONTAL -pos [set $f.p1]
    
}

##############################################################################
#position cursor with mouse on vertical seq
proc SeqPairCursorSet2 {f seq1 seq2 diff result_index x y} {
    global $f.p2 VERTICAL $f.seqdisp_id $f.width
   
    #FIXME - better way of getting char from line.char format 
    set n [$seq2 index @$x,$y]
    regexp {\.[0-9]+$} $n char
    set char [string range $char 1 [string length $char]]
    #puts "char $char"
    set middle [expr [set $f.width] / 2]

    set $f.p2 [expr [set $f.p2] + $char - $middle]
    seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $VERTICAL -pos [set $f.p2]
    
}

##############################################################################
#move cursor
proc SeqPairCursor {f seq1 seq2 number diff result_index seq_dir} {
    global $f.p1 $f.p2 HORIZONTAL VERTICAL $f.seqdisp_id $f.lock

    #puts "POSITION [$seq1 get 1.1] [$seq1 get 2.1]"

    if {[set $f.lock]} {
	seq_scroll12 $f $seq1 $seq2 $number $diff $result_index
	seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $HORIZONTAL -pos [set $f.p1]
	seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $VERTICAL -pos [set $f.p2]
 	
    } elseif {$seq_dir == 1} {
	seq_scroll1 $f $seq1 $seq2 $number $diff $result_index
	seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $HORIZONTAL -pos [set $f.p1]
    } elseif {$seq_dir == 2} {
	seq_scroll2 $f $seq1 $seq2 $number $diff $result_index
	seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $VERTICAL -pos [set $f.p2]

    }
}

##############################################################################
proc SeqDispStartShutdown {f } {
    global $f.p_tag1 $f.p_tag2 $f.p1 $f.p2 $f.seqdisp_id

    if {![info exists $f.seqdisp_id]} {
	return
    }
    #puts "SeqDispStartShutdown $f [set $f.seqdisp_id]"

    seq_result_update -index [set $f.seqdisp_id] -job QUIT
    
    if {[info exists $f.p1]} {
	unset $f.p1 
    }
    if {[info exists $f.p2]} {
	unset $f.p2 
    }
    if {[info exists $f.p_tag1]} {
	unset $f.p_tag1 
    }
    if {[info exists $f.p_tag2]} {
	unset $f.p_tag2 
    }
    destroy $f
}

##############################################################################
#called from C - callback to move seq disp cursor 
proc seq_disp_show_cursor {f result_index pos direction} {
    global $f.p1 $f.p2 $f.max1 $f.max2 $f.width HORIZONTAL VERTICAL $f.lock
    global $f.seqdisp_id $f.been_here
    
    #puts "seq_disp_show_cursor $f $result_index $pos $direction"

    set middle [expr [set $f.width] / 2]

    if {$direction == $HORIZONTAL} {
	set start1 $pos
	set start2 [set $f.p2]
    } else {
	set start1 [set $f.p1]
	set start2 $pos
    }
    if { ($start1 < 1) && ($start2 < 1)} {
	return
    }
    if { ($start2 > [set $f.max2]) && ($start1 > [set $f.max1])} {
	return
    }

    if {[set $f.lock]} {
	#need to stop callback ie if move horizontal cursor, then I need a
	#callback to similarly move the vertical cursor but unless I stop the
	#callback loop here, the vertical cursor movement will cause movement
	#of the horizontal cursor etc etc
	if {![info exists $f.been_here] || [set $f.been_here] == 0} {
	    set $f.been_here 1
	}
	if {$direction == $HORIZONTAL} {
	    if {[set $f.been_here] == 1} {
		set left1 [expr $pos - $middle]
		set w_diff [expr $start1 - [set $f.p1]]
		set left2 [expr [set $f.p2] + $w_diff - $middle]
		set start2 [expr [set $f.p2] + $w_diff]

		set $f.been_here 2
		seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $VERTICAL -pos [expr [set $f.p2] + $w_diff]
		set $f.been_here 0
	    } else {
		return
	    }

	} else {
	    if {[set $f.been_here] == 1} {
		set left2 [expr $pos - $middle]
		set w_diff [expr $start2 - [set $f.p2]]
		set left1 [expr [set $f.p1] + $w_diff - $middle]
		set start1 [expr [set $f.p1] + $w_diff]

		set $f.been_here 2
		seq_pair_move_cursor -seqdisp_id [set $f.seqdisp_id] -direction $HORIZONTAL -pos [expr [set $f.p1] + $w_diff]
		set $f.been_here 0
	    } else {
		return
	    }
	}
    } else {
	if {$direction == $HORIZONTAL} {
	    set left1 [expr $pos - $middle]
	    set left2 [expr [set $f.p2] - $middle]
	} else {
	    set left1 [expr [set $f.p1] - $middle]
	    set left2 [expr $pos - $middle]
	}
    }
    #puts "left1 $left1 $left2"
    seq_replot $f $f.diff $f.seq1 $f.seq2 $left1 $left2 $result_index

    #puts "start $start1 $start2"
    set $f.p1 $start1
    set $f.p2 $start2
}

