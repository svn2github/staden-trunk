#-----------------------------------------------------------------------------
# Cursor commands
#-----------------------------------------------------------------------------

# Create a new visible cursor
proc ed_nip_canvas_cursor_create {seq_id canvwin cursorid regid colour} {
    global nip_defs
    set height [$canvwin canvasy [winfo height $canvwin]]

    #HACK - why 10?
    set line   [$canvwin create line 10 [$canvwin canvasy 0] 10 $height \
	-fill $colour \
	-tag "cursor_$cursorid cursor"]
    $canvwin bind $line <<move-drag>> \
	"ed_nip_canvas_cursor_drag $seq_id $canvwin $cursorid $regid $colour %x"
}

# Delete an existing cursor
proc nip_canvas_cursor_delete {canvwin cursorid} {
    $canvwin delete cursor_$cursorid
}

# Move an existing cursor to 'x'
proc ed_nip_canvas_cursor_move {seq_id canvwin cursorid regid colour x} {
    # Create the cursor if it doesn't currently exist

    if {$seq_id != -1} {
	if {[$canvwin find withtag cursor_$cursorid] == ""} {
	    ed_nip_canvas_cursor_create $seq_id $canvwin $cursorid $regid $colour
	}
    }

    # And now move it
    set height [$canvwin canvasy [winfo height $canvwin]]

    set y1 [lindex [$canvwin cget -scrollregion] 1]
    set y2 [lindex [$canvwin cget -scrollregion] 3]
   
    # no scrolling allowed eg ruler
    if {$y2 == 0} {
	set y2 [winfo height $canvwin]
    }
    $canvwin coords cursor_$cursorid $x $y1 $x $y2

    # delete selection when cursor moved
    $canvwin delete dummy

}

# Drag an existing cursor to pixel coord 'px'
proc ed_nip_canvas_cursor_drag {seq_id canvwin cursorid regid colour px} {
    # Convert the 'px' to a base x.
    set x [ed_nip_canvas_to_world -id $regid -x [$canvwin canvasx $px]]

    # Move the cursor
    ed_nip_canvas_cursor_move -1 $canvwin $cursorid 0 $colour $px

    # Generate an event
    if {$x < 1} {
	set x 1
    } elseif {$x > [expr [ed_s_length -seq_id $seq_id]+1]} {
	set x [expr [ed_s_length -seq_id $seq_id] + 1]
    }

    #keylset l id $cursorid abspos $x sent_by $regid
    ed_cursor_notify -seq_num [get_editor_reg_num $seq_id] -id $cursorid -pos $x
}

# Use for creating a new editor at pixel position 'px', or for moving
# an existing editor (and cursor) to 'px'.
proc ed_nip_canvas_cursor_editor {seq_id regid px canvwin} {

    #puts "ed_nip_canvas_cursor_editor:seq_id=$seq_id"
    # Get the new cursor position in base coordinates
    set x [ed_nip_canvas_to_world -id $regid -x [$canvwin canvasx $px]]
    if {$x < 1} {
	set x 1
    } elseif {$x > [expr [ed_s_length -seq_id $seq_id]+1]} {
	set x [expr [ed_s_length -seq_id $seq_id]+1]
    }

    # Find a cursor id
    set found 0
    foreach item [$canvwin find withtag cursor] {
	#puts "item=$item"
	regsub {.*cursor_([^ ]*|$f).*} [$canvwin gettags $item] {\1} c
	
	# We've now got a cursor id, check if the editor is using it
	# This is done somewhat crudely by checking if the cursor
	# is 'private'.
	#set seq_num [get_seq_num $seq_id]
	#puts "ed_nip_canvas_cursor_editor: c=$c"
	set seq_num [get_editor_reg_num $seq_id]
	set clist [ed_query_cursor -cursorid $c -seq_num $seq_num]
	if {[keylget clist private] != 0 &&
	    [keylget clist id] == $seq_id} {
	    set found 1
	    break
	}
    }
    #if {$found} {
	 #If we've found an editor, then simply move it
	#ed_cursor_notify -seq_num $seq_num -id $c -pos $x
    #} else {
	# Otherwise create a new editor
	set seq_name [get_name_from_id $seq_id]
	CreateSequenceEditor $seq_name $x
    #}
}

#return the restriction enzyme cut position
proc EdGetREnzPosY { plot id } {

    foreach tag [$plot gettags $id] {
	if {[string compare [string range $tag 0 2] re_] == 0} {
	    return [string trim $tag re_]
	}
    }
}

proc GetRenzSeqId {plot id} {

    foreach tag [$plot gettags $id] {
	if {[string compare [string range $tag 0 5] seqid_] == 0} {
	    return [string trim $tag seqid_]
	}
    }
}

proc HighLightFragment {f re_win seq_id y x1 x2 xw1 xw2 name1 name2} {
    
    global tk_utils_defs
    global $f.fragment
    
    set tick_ht [keylget tk_utils_defs R_ENZ.TICK_HEIGHT]
    set yoffset [expr ($y + 2)*$tick_ht]
    
    if {[info exists $f.fragment]} {
	 $re_win delete fragment
    } 
    set $f.fragment [$re_win create rectangle \
	    [expr $x1+1] [expr $yoffset-$tick_ht+2] [expr $x2-1] [expr $yoffset-2] \
	    -fill green \
	    -outline beige \
	    -width 1 \
	    -tag "fragment seqid_$seq_id xL_$xw1 xR_$xw2 name1_$name1 name2_$name2"]

    #Added this binding in case the inseration will be made 
    #in the end of the fragment 
    #$re_win bind fragment <1> "$re_win delete fragment"
}

#find the distance between two restriction enzyme cuts
proc EdFindDistance { f re_win label dist text x name_cmd} {

    global $f.pos 
    global $f.plot
    global $f.posY $f.posw $f.name

    if {[info exists $f.plot]} {
	if {[set $f.plot] != $re_win} {
	    return
	}
    }
    set xw [GetREnzPos $re_win current]

    if { [info exists $f.pos]} {
	$dist configure -text [expr abs($xw - [set $f.pos])] 
	$label configure -text "Distance: [expr abs($xw - [set $f.posw])]"
	set nearest [$re_win find withtag current]
	set name [GetREnzName $f $re_win $nearest $name_cmd]
	set seq_id [GetRenzSeqId $re_win current]

	#HighLightFragment $f $re_win $seq_id [set $f.posY] [set $f.pos] $x [set $f.posw] $xw [set $f.name] $name
	
	sequences_redisplay_graphic select $seq_id [set $f.posw] $xw [set $f.name] $name 
	unset $f.pos
	unset $f.posw
	unset $f.posY
	unset $f.plot
	unset $f.name
	bell
    } else {	
	set $f.pos $x
	set $f.posw $xw
	set $f.posY [EdGetREnzPosY $re_win current]
	set $f.plot $re_win
	set nearest [$re_win find withtag current]
	set $f.name [GetREnzName $f $re_win $nearest $name_cmd]
	$label configure -text "Select another $text"
	bell
    }
}

#create restriction enzyme canvas
proc ed_renz_map {w args} {
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

    EdSetREnzBindings $w $w.renz $w.ruler $w.names $w.brief $w.buttons.dist \
	[set $w.cursor_cmd] [set $w.renz_name_cmd] [set $w.renz_info_cmd] \
       [set $w.invoke_cmd] [set $w.config_cmd]

    SetREnzRulerBindings $w $w.renz $w.ruler [set $w.cursor_cmd] [set $w.invoke_cmd] 
}

proc PopUpEditMenu {f re_win X Y name_cmd} {

    set tags [$re_win gettags current]

    set seq_id 0
    set x1 0 
    set x2 0
    set pos 0
    set name1 ""
    set name2 ""
    set name ""

    foreach tag $tags {
	if {[string compare [string range $tag 0 5] seqid_] == 0} {
	    set seq_id  [string trim $tag seqid_]
	}
	if {[string compare [string range $tag 0 2] xL_] == 0} {
	    set x1  [string trim $tag xL_]
	}
	if {[string compare [string range $tag 0 2] xR_] == 0} {
	    set x2  [string trim $tag xR_]
	}
	if {[string compare [string range $tag 0 5] name1_] == 0} {
	    set name1  [string trim $tag name1_]
	}
	if {[string compare [string range $tag 0 5] name2_] == 0} {
	    set name2  [string trim $tag name2_]
	}
	if {[string compare [string range $tag 0 3] pos_] == 0} {
	    set pos  [string trim $tag pos_]
	}
    }


    if {[winfo exists $re_win.m]} {destroy $re_win.m}
    menu $re_win.m -tearoff no
    foreach tag $tags {
	if {[string compare $tag fragment] == 0} {
	    $re_win.m add command -label cut \
		    -command "destroy $re_win.m; \
		    sequences_redisplay_graphic cut $seq_id $x1 $x2 $name1 $name2"
		    
	    $re_win.m add command -label copy \
		    -command "sequences_redisplay_graphic copy $seq_id $x1 $x2 $name1 $name2"
		    
	    $re_win.m add command -label paste \
		    -state disabled 
	    $re_win.m add command -label replace \
		    -command "sequences_redisplay_graphic replace $seq_id $x1 $x2 $name1 $name2"
		    
	}
	if {[string compare $tag S] == 0} {
	    set nearest [$re_win find withtag current]
	    set name [GetREnzName $f $re_win $nearest $name_cmd]
	    $re_win.m add command -label cut \
		    -state disabled
	    $re_win.m add command -label copy \
		    -state disabled
	    $re_win.m add command -label paste \
		    -command "sequences_redisplay_graphic paste $seq_id $pos $name"
		    
	    $re_win.m add command -label replace \
		    -state disabled  
	}
	 tk_popup $re_win.m [expr $X-20] [expr $Y+5] 
    }
}

#bindings specific to the stand alone restriction enzyme 
proc EdSetREnzBindings {f re_win r_win names_win label dist cursor_cmd renz_name_cmd renz_info_cmd invoke_cmd config_cmd} {
    global REnzY1
    global $f.renz_id $names_win.Move

    bind $re_win <Any-Leave> "DeleteREnzCursor $re_win $r_win"

    bind $re_win <Any-Motion> "AddREnzCursor $f $re_win $r_win %x $cursor_cmd"

    #any mouse motion - highlight nearest cut line
    $re_win bind S <Any-Motion> \
	"HighlightREnz $f $re_win $label $renz_name_cmd"
    $re_win bind S <Shift-Motion> {;}

    #button-1 in plot canvas find the distance between 2 cut lines
    $re_win bind S <<select>> "EdFindDistance $f %W $label $dist cut %x $renz_name_cmd"

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

    # added for sequence editing 
    #button-2 in plot canvas - invoke a popup menu to do edit
    bind $re_win <2> "PopUpEditMenu $f $re_win %X %Y $renz_name_cmd"
    # added for sequence editing

    #auto scrolling
    bind $names_win <<move-autoscroll>> "set $f.auto_scroll 1; REnzAutoScroll $names_win $re_win %y $f y"

    bind $names_win <<stop-autoscroll>> "+REnzStopScroll $f"

    #button-3 in plot canvas - invoke a popup menu
    #$re_win bind S <<menu>> "PopUpREnzMenu $f $re_win %X %Y $renz_info_cmd"

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


