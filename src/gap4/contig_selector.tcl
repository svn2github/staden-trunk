#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

##############################################################################
#hack to get auto scrolling to stop if the mouse re-enters the canvas 
proc StopScroll {f} {
    global NGRec

    set NGRec($f,auto_scroll) 0
}

##############################################################################
#auto scrolls whilst the mouse pointer is outside the canvas
proc DoScroll { f csh_win csp_win direction scroll } {
    global NGRec
    global CanvasConst

    if {$NGRec($f,auto_scroll)} {
	set bbox [$csh_win bbox S]
	set min_x [lindex $bbox 0]
	set max_x [lindex $bbox 2]
	set min_y [lindex $bbox 1]
	set max_y [lindex $bbox 3]
	
	set first_unit [lindex [$f.hscroll get] 0]
	set unit [expr double($CanvasConst(auto_incr))/($max_x - $min_x)]
	$csh_win xview moveto [expr $first_unit + ($direction * $unit)]

	if {[winfo exists $csp_win] } {
	    $csp_win xview moveto [expr $first_unit + ($direction * $unit)]
	}
	after 100 "DoScroll $f $csh_win $csp_win $direction $scroll" 
    }
}

##############################################################################
proc AutoScroll {canvas csp_win x y f scroll} {
    global CanvasConst

    #puts "start AUTOSCROLL y $y"

    set min_x 1
    set max_x [winfo width $canvas]
    if {$x < $min_x } {
	after $CanvasConst(auto_time) "DoScroll $f $canvas $csp_win -1 $scroll"
    } elseif {$x > $max_x} {
	after $CanvasConst(auto_time) "DoScroll $f $canvas $csp_win +1 $scroll"
    }
} 

##############################################################################
#move a selected item in a canvas
proc itemMove {io csh_win csp_win x y} {
    global $csh_win.areaX1 $csh_win.areaY1

    set min_x 1
    set max_x [winfo width $csh_win]
    set min_y 1
    set max_y [winfo height $csh_win]

    set list [GetContigList]
    set list [remove_contig_duplicates -io $io -contigs $list]

    if {($x > $min_x) && ($x < $max_x) && ($y > $min_y) && ($y < $max_y)} {
	set x [$csh_win canvasx $x]
	set y [$csh_win canvasy $y]
	set dx [expr $x - [set $csh_win.areaX1]]
	#set dy [expr $y - [set $csh_win.areaY1]]

	#disallow movement in the y direction
	set dy 0
	foreach i $list {
	    set c_num [db_info get_contig_num $io $i]
	    $csh_win move num_$c_num $dx $dy
	    if {[winfo exists $csp_win]} {
		$csp_win move num_$c_num $dx $dy
	    }
	}
	set $csh_win.areaX1 $x
	set $csh_win.areaY1 $y
    }
}

##############################################################################
#drop selected contig list and reorder the contiguous line of contigs
proc itemDrop {io f csh_win x} {
    global $f.cs_id

    #unmark the canvas to avoid accidental rubber-banding by <<zoom-drag>>
    itemUnMark $csh_win

    set x [$csh_win canvasx $x]
    set order_list [GetContigList]

    update_contig_order -io $io -id [set $f.cs_id] -contigs "$order_list" -x $x

    #for the moment, clear the "contigs" list upon completion of the drag
    ClearContigSelection $io $csh_win
}

##############################################################################
#local
#double the height of id
#height of the left hand separator
proc HighlightSeparator { csh_win pos} {
    global gap_defs

    set id [$csh_win find withtag sep_$pos]

    set tick_ht [keylget gap_defs CONTIG_SEL.TICK_HEIGHT]
    set x [lindex [$csh_win coords $id] 0]
    set y1 [lindex [$csh_win coords $id] 1]
    set y2 [lindex [$csh_win coords $id] 3]
    set y1 [expr $y1 - $tick_ht]
    set y2 [expr $y2 + $tick_ht]
    $csh_win coords $id $x $y1 $x $y2
}

##############################################################################
proc ActivateCSMenu {m state } {

    $m entryconfigure "Clear all" -state $state
    $m entryconfigure "Display diagonal" -state $state
}

##############################################################################
#display a vertical view of the contig selector but without buttons and
#selection facilities and the dot plot
proc ContigComparator { io } {
    global NGRec
    global gap_defs

    set f [keylget gap_defs CONTIG_SEL.WIN]

    global $f.cs_id
    global $f.diagonal

    #puts "%%%%%%%%%%%%%%%%%%start CSPLOT%%%%%%%%%%%%%%%%%%"

    #if contig selector does not exist, create it
    if {![winfo exists $f]} {
	ContigSelector $io
    }
    raise $f
    wm deiconify $f
    wm resizable $f 1 1
    wm title $f "Contig Comparator"

    #if contig comparator does exist, return, else create it
    if {[winfo exists $f[keylget gap_defs CONTIG_CSP.WIN]]} {
	return
    }

    set NGRec($f,curr_item) 0; #used in highlightnearest
    set NGRec($f,last_result) 0; #used in InvokeNext

    if {![info exists $f.diagonal]} {
	set $f.diagonal 0; #diagonal not to be displayed
    }
    ActivateCSMenu $f.menubar.view normal

    set numcontigs [db_info num_contigs $io]

    #allow -height and -width to have affect
    wm geometry $f {}; #needed to grow window if resized cs_h

#    update idletasks; #necessary for sizing of canvases correctly!
    
    set csh_win $f[keylget gap_defs CONTIG_SEL.WIN]
    set csv_win $f[keylget gap_defs CONTIG_CSV.WIN]
    set csp_win $f[keylget gap_defs CONTIG_CSP.WIN]
    set borderwidth [keylget gap_defs CONTIG_SEL.BORDERWIDTH]

    set width [winfo width $csh_win]
    set height [winfo height $csh_win]

    $csh_win configure -height $height -width $width

    #create cs_v and cs_plot windows
    frame $f.csv -bd $borderwidth -relief groove
    set zoom_cmdv [list "gap_zoom $io \[set $f.cs_id\] y"]
    canvasbox $csv_win -height $width \
	-width $height\
	-yscrollcommand "$f.vscroll set" \
	-bd 0 -highlightthickness 0\
	-zoom_command $zoom_cmdv 
    
    frame $f.csp -bd $borderwidth -relief groove
    set zoom_cmdp [list "gap_zoom $io \[set $f.cs_id\] b"]
    canvasbox $csp_win -height $width -width $width \
	-bd 0 -highlightthickness 0\
	-cursor crosshair -closeenough 4 \
	-zoom_command $zoom_cmdp 
    
    scrollbar $f.vscroll -relief sunken -command "gc_scroll_y $io \[set $f.cs_id\]"

    grid rowconfig $f 3 -weight 1

    grid $f.csv -row 3 -column 0 -sticky ns
    pack $csv_win -in $f.csv -pady 5 -fill y -expand yes
    grid $f.vscroll -row 3 -column 2 -sticky ns
    grid $f.csp -row 3 -column 1 -sticky nsew
    pack $csp_win -in $f.csp -padx 5 -pady 5 -fill both -expand yes

    #HACK - cp_id is never used
    update idletasks

    catch {tkwait visibility $csv_win}

    set cp_id [display_contig_comparator -io $io \
	    -window $csp_win -win_vertical $csv_win -id [set $f.cs_id]]

    #need this to ensure that a resize event happens after call to
    #display_contig_comparator which sets up the new canvas structure values
    resize_canvas -io $io -id [set $f.cs_id]

    #bind the configure actions to the toplevel
    bind $f <Any-Configure> "
	    if {\[winfo toplevel %W\] == \"%W\"} {
		update idletasks
		resize_canvas -io $io -id [set $f.cs_id]
	    }
	"
    SetCanvasBindings $csp_win $zoom_cmdp
    SetCanvasBindings $csv_win $zoom_cmdv

    #SetCanvasBindings $io [set $f.cs_id]  $csp_win "b"
    #SetCanvasBindings $io [set $f.cs_id]  $csv_win "y"
    SetCSPBindings $io $f $csp_win $csh_win $csv_win $f.brief    
    SetCSVBindings $io $f $csv_win
}

##############################################################################
#plot tags on contig selector
proc DisplayCSTags { io f csh_win} {
    global $f.cs_id

    #clear all tags
    $csh_win delete tag

    display_cs_tags -io $io -id [set $f.cs_id]
}

##############################################################################
proc itemBind {io canvas x y} {

    set c_list [GetContigNum $canvas current]
    #do not toggle contig selection where selecting contigs with the "move"
    #mouse button
    UpdateContigListNoToggle $io $c_list
    UpdateContigId $io $canvas current

    #SelectSingleContig $io $canvas
    itemMark $canvas $x $y
}

##############################################################################
#highlight contigs as move mouse over
proc HighlightItems {io f plot brief} {
    global restoreCmd
    global gap_defs
    global $f.prev_item

    if {![info exists $f.prev_item]} {
	set $f.prev_item ""
    }
    set nearest [$plot find withtag current]

    #puts "NEAREST $nearest prev [set $f.prev_item]"
    #only do this code if the nearest item is different from the previous one
    if {$nearest != [set $f.prev_item]} {
 
	#unhighlight object
    	if {$restoreCmd($f,item_num) != 0} {
	    eval $restoreCmd($f,reading)
	    #keep readings always raised
	    set item [lindex $restoreCmd($f,reading) 2]
	    if {[string compare [lindex [GetItemInfo $io $plot $item] 0] \
		    Reading:] == 0} {
		$plot raise $item
	    }
	} 

	if {$nearest != 0} {
	    InitialSettings $f $plot $nearest reading
	    $plot itemconfig $nearest \
		    -fill [keylget gap_defs CONTIG_SEL.HIGHLIGHT_COLOUR]
	    $plot raise $nearest
	    set restoreCmd($f,item_num) $nearest
	    $brief configure -text [GetItemInfo $io $plot $nearest]
	}
    }
    #set previous item
    set $f.prev_item $nearest
}

##############################################################################
proc UnHighlight {f plot } {
    global restoreCmd

    #unhighlight object
    if {$restoreCmd($f,match_num) != 0} {
	foreach tag [$plot gettags $restoreCmd($f,match_num)] {
	    if {[string compare [string range $tag 0 3] num_] == 0} {
		eval $restoreCmd($f,c_h)
		eval $restoreCmd($f,c_v)
	    }
	}
	eval $restoreCmd($f,item)
    } 

}

##############################################################################
#Highlights an match
proc HighlightMatch {f csh_win csv_win plot brief item} {
    global restoreCmd
    global NGRec
    global gap_defs

    if {$restoreCmd($f,match_num) != 0} {
        foreach tag [$plot gettags $restoreCmd($f,match_num)] {
    	    if {[string compare [string range $tag 0 3] num_] == 0} {
    	        eval $restoreCmd($f,c_h)
    	        eval $restoreCmd($f,c_v)
    	    }
        }
        eval $restoreCmd($f,item)
    } 

    if {"[$plot type $item]" == ""} {
	set item 0
    }

    if {$item != 0} {
        #highlight object
        InitialSettings $f $plot $item item
        #find which contigs are involved
        set i 0
        foreach tag [$plot gettags $item] {
    	    if {[string compare [string range $tag 0 3] num_] == 0} {
    	        set contig($i) $tag
    	        incr i
    	    }
        }
        
        #highlight the item
        $plot itemconfig $item \
		-fill [keylget gap_defs CONTIG_SEL.HIGHLIGHT_COLOUR]
        $plot raise $item
        
        #find left hand coords of contigs
        set pos1 [lindex [$csh_win coords $contig(0)] 0]
        set pos2 [lindex [$csh_win coords $contig(1)] 0]
        #highlight the contigs aswell, horizontal contigs > vertical 
        if {$pos1 > $pos2} {
    	    InitialSettings $f $csh_win $contig(0) c_h
    	    InitialSettings $f $csv_win $contig(1) c_v
    	    $csh_win itemconfig $contig(0) -fill [keylget gap_defs CONTIG_SEL.HIGHLIGHT_COLOUR]
    	    $csv_win itemconfig $contig(1) -fill [keylget gap_defs CONTIG_SEL.HIGHLIGHT_COLOUR]
        } else {
    	    InitialSettings $f $csh_win $contig(1) c_h
    	    InitialSettings $f $csv_win $contig(0) c_v
    	    $csh_win itemconfig $contig(1) -fill [keylget gap_defs CONTIG_SEL.HIGHLIGHT_COLOUR]
    	    $csv_win itemconfig $contig(0) -fill [keylget gap_defs CONTIG_SEL.HIGHLIGHT_COLOUR]
        }
        set restoreCmd($f,match_num) $item
        $brief configure -text [obj_get_brief $item]
    }

    #set previous item
    set NGRec($f,curr_item) $item
}

##############################################################################
#highlight the nearest line and the associated contigs
proc HighlightNearest {f hor ver plot brief} {
    global restoreCmd
    global NGRec

    set nearest [$plot find withtag current]

    #puts "nearest $nearest"
    #only do this code if the nearest item is different from the previous one
    if {$nearest != $NGRec($f,curr_item)} {
	HighlightMatch $f $hor $ver $plot $brief $nearest
    }
}

##############################################################################
#bindings specific to vertical contig selector display
proc SetCSVBindings { io f csv_win} {
    global $f.cs_id
    global restoreCmd

    bind $csv_win <Any-Leave> "delete_canvas_cursor -io $io -id [set $f.cs_id]"
    bind $csv_win <Any-Motion> "AddCSVCrossHair $io [set $f.cs_id] $f $csv_win %y"
}

##############################################################################
proc SetCSPBindings {io f csp_win csh_win csv_win brief} {
    global $f.cs_id
    global restoreCmd
   
    bind $csp_win <Any-Leave> "delete_canvas_cursor -io $io -id [set $f.cs_id]"
    bind $csp_win <Any-Motion> "AddCSHCrossHair $io [set $f.cs_id] $f $csp_win %x"
    bind $csp_win <Any-Motion> "+AddCSVCrossHair $io [set $f.cs_id] $f $csp_win %y"
    #HACK - redo?
    set restoreCmd($f,match_num) 0
    #highlighting the nearest object to the cursor
    #if leave the cs_plot canvas, unhighlight current selection
    bind $csp_win <Any-Leave> "+UnHighlight $f $csp_win"
    bind $csp_win <Any-Enter> "+\
	    InitialSettings $f $csp_win $restoreCmd($f,match_num) item; \
	    HighlightMatch $f $csh_win $csv_win %W $brief \$NGRec($f,curr_item)"

    $csp_win bind S <Any-Motion> [format {
	HighlightNearest %s %s %s %%W %s
    } [list $f] [list $csh_win] [list $csv_win] [list $brief] ]
    $csp_win bind S <Shift-Motion> {;}

    bind $csp_win <<menu>> [format {+ 
        PopUpMenu %s %%W %%x %%y %%X %%Y %s
    } [list $f] [list $csp_win] ]

    bind $csp_win <<use>> "DefaultOp $f $csp_win %x %y"
}

##############################################################################
# Call the default operation
proc DefaultOp {f canvas x y} {
    global NGRec

    if {$NGRec($f,curr_item) == 0} {
	return
    }

    # Find object number
    set nearest $NGRec($f,curr_item)
    set result 0
    foreach tag [$canvas gettags $nearest] {
	if {[string compare $tag ignore] == 0} {
	    return
	}
	if {[string match p_* $tag]} {
	    set result $tag
	}
    }

    set count 0
    obj_get_ops list $nearest
    obj_invoke_op $nearest -2; # -2 is an alias for the default operation
}

##############################################################################
# Sets the last used result to 'result'
proc CSLastUsed {result} {
    global NGRec gap_defs

    set f [keylget gap_defs CONTIG_SEL.WIN]
    set NGRec($f,last_result) $result
    $f.buttons.next configure -state normal
}

proc CSLastUsedFree {result} {
    global NGRec gap_defs

    set f [keylget gap_defs CONTIG_SEL.WIN]
    if {"$result" == "$NGRec($f,last_result)"} {
	set NGRec($f,last_result) 0
	$f.buttons.next configure -state disabled
	UnHighlight $f $f.csp
    }

}

##############################################################################
# Call the default operation on the next match
proc InvokeNext {f} {
    global NGRec

    if {$NGRec($f,last_result) == 0} {
	return
    }

    obj_invoke_next $NGRec($f,last_result)
}

proc ScanNext {f brief horiz ver canvas} {
    global NGRec

    if {$NGRec($f,last_result) == 0} {
	return
    }

    set next [obj_get_next $NGRec($f,last_result)]
    if {$next != -1} {
       	HighlightMatch $f $horiz $ver $canvas $brief $next
    }
}

##############################################################################
#toggle diagonal line on and off
proc DisplayDiagonal {f plot io} {
    global $f.diagonal
    global $f.cs_id

    if {[set $f.diagonal]} {
	$plot delete diagonal
	display_cs_diagonal -io $io -id [set $f.cs_id]
    } else {
	$plot delete diagonal
    }
}

##############################################################################
#display the contig selector crosshair
proc AddCSHCrossHair {io id f cs_win x} {
    global $f.cursor

    if {[set $f.cursor]} {
	draw_canvas_cursor_x -io $io -id $id -x [$cs_win canvasx $x]
    } else {
	delete_canvas_cursor -io $io -id $id
    }
}

##############################################################################
#display the contig selector crosshair
proc AddCSVCrossHair {io id f cs_win y} {
    global $f.cursor

    if {[set $f.cursor]} {
	draw_canvas_cursor_y -io $io -id $id -y [$cs_win canvasy $y]
    } else {
	delete_canvas_cursor -io $io -id $id
    }
}

##############################################################################
#bindings specific to contig selector display
proc SetCSHBindings { io f cs_h csp_win select_item  scroll brief label_x1 label_x2} {
    global NGRec
    global NGWinType
    global restoreCmd
    global $f.cs_id
    global $cs_h.Move $cs_h.Select

    set $cs_h.Move 0
    set $cs_h.Select 0
    
    set restoreCmd($f,item_num) 0

    #puts "start SETCSBINDINGS"
    set csh_win $cs_h
    bind $csh_win <Any-Leave> "delete_canvas_cursor -io $io -id [set $f.cs_id]"
    bind $csh_win <Any-Motion> "AddCSHCrossHair $io [set $f.cs_id] $f $cs_h %x"

    #button 2 for moving - select and move
    $cs_h bind contig <<move>> \
	"if {\[set $cs_h.Select\] == 0} {
            set $cs_h.Move 1
    	    itemBind $io $csh_win %x %y
	}"
	
    #moving a contig
    $cs_h bind contig <<move-drag>> \
	"if {\[set $cs_h.Move\] == 1} {
    	    itemMove $io %W $csp_win %x %y
	}"

    $cs_h bind contig <<move-release>> \
	"if {\[set $cs_h.Move\] == 1} {
    	   itemDrop $io $f %W %x
	   set $cs_h.Move 0
	}"
	   
    #auto scrolling
    bind $cs_h <<move-autoscroll>> [format {+ 
	set NGRec(%s,auto_scroll) 1; 
	AutoScroll %%W %s %%x %%y %s %s
    } [list $f] [list $csp_win] [list $f] [list $scroll] ]

    bind $cs_h <<stop-autoscroll>> "+StopScroll $f"


    bind $cs_h <Key-c> "ClearContigSelection $io $cs_h"
    #button 1 for selection
    $cs_h bind contig <<select>> "SelectSingleContig $io $cs_h"
    bind $cs_h <<select>> \
	"if {\[set $cs_h.Move\] == 0} {
	    set $cs_h.Select 1
    	    itemMark $cs_h %x %y
	}"
    bind $cs_h <<select-drag>> \
	"if {\[set $cs_h.Select\] == 1} {
    	    itemStroke $cs_h %x %y
	}"
    bind $cs_h <<select-release>> \
	"if {\[set $cs_h.Select\] == 1} {
	    SelectContigs $f $cs_h $io
	    set $cs_h.Select 0
	}"

    #bind button 3 to pop-up menu
    $cs_h bind tag <<menu>> \
	"PopUpCSTagMenu $io %W \[%W find withtag current\] %X %Y"
    $cs_h bind contig <<menu>> \
	"PopUpCSContigMenu $io %W \[%W find withtag current\] %X %Y %x %y"

    $cs_h bind contig <Any-Motion> "HighlightItems $io $f %W $brief"
    $cs_h bind contig <Shift-Motion> {;}

    bind $cs_h <Any-Leave> "+UnHighlightItems $f"
}

##############################################################################
#called from C: rehighlight contigs after redrawing of contig selector
proc ReHighlightContigSelection { io csh_win} {
    global gap_defs

    set list [GetContigList]
    foreach c $list {
	set c_num [db_info get_contig_num $io [lindex $c 0]]
	$csh_win itemconfig hl_$c_num \
		    -width [keylget gap_defs "CONTIG_SEL.LINE_BOLD"]
    }

}

##############################################################################
proc ClearContigSelection { io csh_win} {
    global gap_defs

    set list [GetContigList]
    foreach c $list {
	set c_num [db_info get_contig_num $io [lindex $c 0]]
	$csh_win itemconfig hl_$c_num \
		    -width [keylget gap_defs "CONTIG_SEL.LINE_WIDTH"]
    }
    ListClear contigs
}

##############################################################################
proc DrawSelection { io csh_win h_list } {
    global gap_defs ${csh_win}.Selected

    if {![info exists ${csh_win}.Selected]} {
	array set ${csh_win}.Selected {}
    }

    # Turn the list into a hash for fast lookup
    array set items ""
    foreach i $h_list {
	set items($i) 1
    }

    # Unselect any items in Selected but not in h_list
    set fine_width [keylget gap_defs "CONTIG_SEL.LINE_WIDTH"]
    foreach c [array names ${csh_win}.Selected] {
	if {![info exists items($c)]} {
	    set c_num [db_info get_contig_num $io $c]
	    $csh_win itemconfig hl_$c_num -width $fine_width
	    unset ${csh_win}.Selected($c)
	}
    }

    # Select any items in h_list but not in Selected
    set bold_width [keylget gap_defs "CONTIG_SEL.LINE_BOLD"]
    foreach c $h_list {
	if {![info exists ${csh_win}.Selected($c)]} {
	    set c_num [db_info get_contig_num $io $c]
	    $csh_win itemconfig hl_$c_num -width $bold_width
	    set ${csh_win}.Selected($c) 1
	}
    }
}

##############################################################################
#updates the list of contigs if select from contig selector
#HACK not implemented yet - should be called from DrawSelection
proc HighlightListContig {csh_win item highlight} {

    set f [winfo parent $csh_win]
    set t $f.l

    set index ""
    foreach tag [$csh_win gettags $item] {
	if {[string compare [string range $tag 0 1] c_] == 0} {
	    set index [string trim $tag c_]
	    $t.list selection set $index
	}
    }

}

##############################################################################
#click on single contig
proc SelectSingleContig {io csh_win } {

    lappend c_list [GetContigNum $csh_win current]
    UpdateContigList $io $c_list
    UpdateContigId $io $csh_win current
}

##############################################################################
#when drag out region in contig selector, want to toggle between highlighted
#and non-highlighted contigs
proc SelectContigs {f csh_win io } {

    set c_list ""

    set list [itemsUnderArea $f $csh_win contig]

    #unmark the canvas to avoid accidental rubber-banding by <<zoom-drag>>
    itemUnMark $csh_win

    if {[llength $list] == 0} {
	return
    }

    #must ensure that contigs are in left to right order. This is not
    #guarenteed with "$canvas find enclosed"
    set order_list ""
    foreach i $list {
	foreach tag [$csh_win gettags $i] {
	    if {[string compare [string range $tag 0 1] c_] == 0} {
		set order [string trim $tag c_]
		lappend order_list $order
	    }
	}
    }

    set order_list [lsort -integer -increasing $order_list]

    global $csh_win.Cnum
    foreach i $order_list {
	lappend c_list [GetContigNum $csh_win [set ${csh_win}.Cnum($i)]]
    }    

    UpdateContigList $io $c_list
}


##############################################################################
#print contig selector tag info
proc PrintCSTagDetails { io canvas current } {

    foreach tag [$canvas gettags $current] {
	if {[string compare [string range $tag 0 1] t_] == 0} {
	    set t_num [string trim $tag t_]
	}
    }
    set a [io_read_annotation $io $t_num]

    set str ""
    append str  "position [keylget a position] \n"
    append str "length [keylget a length] \n"
    append str "type [keylget a type] \n"
    if {[keylget a annotation] != 0} {
	append str "comment '[io_read_text $io [keylget a annotation]]' \n"
    }

    
    vfuncgroup 3 "Contig selector"
    #vfuncheader "Contig selector"
    vmessage $str
    messagebox $canvas $str
}

##############################################################################
#Edit the contig at the location of a tag
proc EditCSTagDetails { io canvas current } {
    set c_num [GetContigNum $canvas $current]
    if {$c_num == 0} {
        bell;
        return
    }

    set t_num 0
    set r_num 0

    foreach tag [$canvas gettags $current] {
	if {[string match t_* $tag]} {
	    set t_num [string trim $tag t_]
	}
	if {[string match rnum_* $tag]} {
	    set r_num [string trim $tag rnum_]
	}
    }
    set a [io_read_annotation $io $t_num]

    set t_pos [keylget a position]
    if {$r_num != 0} {
	set r [io_read_reading $io $r_num]
	if {[keylget r sense] == 0} {
	    set t_pos [expr {[keylget r position]+$t_pos-[keylget r start]-1}]
	} else {
	    set t_pos [expr {[keylget r position]+[keylget r length]-
			     ($t_pos + [keylget a length]-1) -
			     [keylget r start]}]
	}
    }
    
    edit_contig -io $io -contig [left_gel $io $c_num] -pos $t_pos
}

##############################################################################
#update the contig identifier box from the current contig in the contig
#selector
proc UpdateContigId { io cs_h current } {
    global c_id_contig
    
    set c_id_contig [left_gel $io [GetContigNum $cs_h $current]]
}

##############################################################################
proc PopUpCSTagMenu {io canvas current X Y} {

    if {[winfo exists $canvas.m]} {destroy $canvas.m}

    create_popup $canvas.m "Tag Commands"
    $canvas.m add command -label information \
	    -command "destroy $canvas.m; \
	    PrintCSTagDetails $io $canvas $current"
    $canvas.m add command -label "Edit contig at tag" \
	    -command "destroy $canvas.m; \
	    EditCSTagDetails $io $canvas $current"
    tk_popup $canvas.m [expr $X-20] [expr $Y-10]

}


##############################################################################
proc popup_cs_contig_1 {io canvas obj cpos} {

    set c_num [GetContigNum $canvas $obj]
    if {$c_num == 0} {
        bell;
        return
    }

    destroy $canvas.m
    edit_contig -io $io -contig [left_gel $io $c_num] -pos $cpos
}

proc popup_cs_contig_2 {io canvas obj} {
    global gap_defs

    set c_num [GetContigNum $canvas $obj]
    if {$c_num == 0} {
        bell;
        return
    }

    destroy $canvas.m
    set c_name [left_gel $io $c_num]
    CreateTemplateDisplay $io $c_name
}

proc popup_cs_contig_3 {io canvas obj} {
    global gap_defs

    set c_num [GetContigNum $canvas $obj]
    if {$c_num == 0} {
        bell;
        return
    }

    destroy $canvas.m

    complement_contig -io $io -contigs "=$c_num"
    SetContigGlobals $io [left_gel $io $c_num]
}

proc popup_cs_contig_cnotes {io canvas obj} {
    global gap_defs
    set c_num 0
    
    set c_num [GetContigNum $canvas $obj]
    if {$c_num == 0} {
        bell;
        return
    }

    destroy $canvas.m

    NoteSelector $io contig "=$c_num"
    SetContigGlobals $io [left_gel $io $c_num]
}

proc PopUpCSContigMenu {io canvas current X Y x y} {
    global read_only gap_defs
    if {[winfo exists $canvas.m]} {destroy $canvas.m}

    set cnum [GetContigNum $canvas $current]
    # Work out position clicked within the contig
    if {($x == 0 && $y == 0) || [keylget gap_defs CONTIG_SEL.EDITOR_POS_1]} {
	set cpos 1
    } else {
	set c [io_read_contig $io $cnum]
	set clen [keylget c length]
	foreach {x1 y1 x2 y2} [$canvas coords current] {}
	set x [$canvas canvasx $x]
	set cpos [expr {int(((double($x-$x1))/($x2-$x1)) * $clen)}]
	if {$cpos < 1} {set cpos 1}
	if {$cpos > $clen} {set cpos $clen}
    }

    set name [left_gel $io $cnum]
    create_popup $canvas.m "Contig Commands ($name)"
    $canvas.m add command -label "Edit contig" \
	    -command "popup_cs_contig_1 $io $canvas $current $cpos"
    $canvas.m add command -label "Template display" \
	    -command "popup_cs_contig_2 $io $canvas $current"
    if {!$read_only} {
	$canvas.m add command -label "Complement contig" \
	    -command "popup_cs_contig_3 $io $canvas $current"
    }
    $canvas.m add command -label "List notes" \
	-command "popup_cs_contig_cnotes $io $canvas $current"
    tk_popup $canvas.m [expr $X-20] [expr $Y-10]
}


##############################################################################
#create a pop-up menu on an item "near" to the cursor
proc PopUpMenu { f canvas x y X Y parent} {
    global NGRec

    set nearest $NGRec($f,curr_item)

    set ignore 0; #false
    foreach tag [$canvas gettags $nearest] {
	if {[string compare $tag ignore] == 0} {
	    set ignore 1; #true
	}
    }

    if {($nearest != 0) && !$ignore} {
	# Create menu from operations list
	if {[winfo exists $parent.m]} {destroy $parent.m}
	obj_get_ops list $nearest
	if {$list == ""} { return }

	set count 0
        create_popup $parent.m Commands
	foreach i $list {
	    if {$i == "SEPARATOR"} {
		$parent.m add separator
	    } elseif {$i == "IGNORE"} {
		incr count
      	    } else {
       	        $parent.m add command -label $i \
		    -command "destroy $parent.m; obj_invoke_op $nearest $count"
	        incr count
       	    }
	}

	tk_popup $parent.m [expr $X-20] [expr $Y-10]
    }
}

##############################################################################
#get the contig number from the contig selector
proc GetContigNum {cs id } {
    set contig_num ""
    foreach tag [$cs gettags $id] {
	if {[string compare [string range $tag 0 3] num_] == 0} {
	    set contig_num [string trim $tag num_]
	}
    }
    return $contig_num
}


##############################################################################
proc cs_config_colour {cs result colour} {
    global restoreCmd

    set f [winfo parent $cs]
    $cs itemconfigure $result -fill $colour    
    InitialSettings $f $cs $result item
    eval $restoreCmd($f,item)
    matchresult_configure -result $result -colour $colour -csplot $cs
}

##############################################################################
proc cs_config_width {cs result width} {

    $cs itemconfigure $result -width $width
    matchresult_configure -result $result -width $width -csplot $cs
}

##############################################################################
proc cs_config_quit {cs result} {
    set t .csconfig_$result

    if {[winfo exists $t]} {destroy $t}
}

##############################################################################
proc UpdateCSPlotWidth {canvas item width } {

    $canvas itemconfigure $item -width $width
}

##############################################################################
proc UpdateCSPlot {f canvas item colour } {

    $canvas itemconfigure $item -fill $colour
    InitialSettings $f $canvas $item item
}

##############################################################################
proc get_width {lw } {

    set width [$lw.scale get]
    return $width
}

##############################################################################
proc cs_config {cs result} {
    set lw [lindex [$cs itemconfigure $result -width] 4]
    if {"$lw" == ""} {
	bell
       	return
    }

    set t .csconfig_$result
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Configure result"

    set f [winfo parent $cs]
    # Line width
    frame $t.lw -bd 3 -relief raised
    scale $t.lw.scale \
	    -label "Line width" \
	    -from 0 \
	    -to 10 \
	    -orient horiz \
	    -command "UpdateCSPlotWidth $cs $result"
    $t.lw.scale set $lw

    pack $t.lw -side top -fill both
    pack $t.lw.scale -side top -fill x

    # Colour
    set colour [lindex [$cs itemconfigure $result -fill] 4]
    set width [lindex [$cs itemconfigure $result -width] 4]

    #cmd to execute when ok button on colourbox pressed
    set ok_cmd "cs_config_width $cs $result \[get_width $t.lw]; cs_config_colour $cs $result"
    #cmd to execute when changing colours on colourbox
    set update_cmd "UpdateCSPlot $f $cs $result"
    #set update_cmd "$cs itemconfigure $result -fill $colour"
    #cmd to execute when cancel button on colourbox pressed
    set cancel_cmd "$cs itemconfigure $result -fill $colour -width $width; \
	        InitialSettings $f $cs $result item"

    ColourBox $t $colour $ok_cmd $update_cmd $cancel_cmd
    
    pack $t.col -side top -fill both
}

##############################################################################
proc CSTagCheckList {io parent path canvas} {
    TagDialog CONTIG_SEL.TAGS $path "DisplayCSTags $io $parent $canvas" {}
}

##############################################################################
#
proc DeleteCSPlotCanvas {f } {
    global gap_defs

    set new_height [expr [winfo height $f] - [winfo height $f[keylget gap_defs CONTIG_CSP.WIN]]]
    destroy $f[keylget gap_defs CONTIG_CSP.WIN]
    destroy $f.vscroll
    destroy $f[keylget gap_defs CONTIG_CSV.WIN]
    destroy $f.csp
    destroy $f.csv
    wm geometry $f [winfo width $f]x$new_height

}

##############################################################################
proc CSClearAll {io f csp_win csv_win } {
    global $f.cs_id
    global gap_defs

    wm title $f "Contig Selector"

    set cursor_ty [keylget gap_defs CONTIG_SEL.CURSOR1_Y]
    set cursor_ly [keylget gap_defs CONTIG_SEL.CURSOR2_Y]
    $f$cursor_ty configure -text ""
    $f$cursor_ly configure -text ""

    ActivateCSMenu $f.menubar.view disabled
    set new_height [expr [winfo height $f] - [winfo height $csp_win]]

    destroy $f.vscroll
    delete_window -io $io -id [set $f.cs_id] -window $csp_win
    delete_window -io $io -id [set $f.cs_id] -window $csv_win
    destroy $f[keylget gap_defs CONTIG_CSP.WIN]
    destroy $f[keylget gap_defs CONTIG_CSV.WIN]
    destroy $f.csp
    destroy $f.csv

    #MUST be done after delete_window
    clear_cp -io $io -id [set $f.cs_id]
    wm resizable $f 1 0
}

##############################################################################
#HACK to complete
proc DeleteContigSelector {f} {

    if {[winfo exists $f]} {
	destroy $f
    }
}

##############################################################################
proc CSStartShutdown {io f} {
    global $f.cs_id NGList gap_defs

    set csh_win $f[keylget gap_defs CONTIG_SEL.WIN]
    set trace_cmd "ContigSelector_ContigsList $io $csh_win"
    catch {trace vdelete NGList(contigs) w [list $trace_cmd]}
    global $csh_win.Selected
    catch {unset ${csh_win}.Selected}

    if {[info exists $f.cs_id]} {
	result_quit -io $io -id [set $f.cs_id]
    }
}

##############################################################################
proc sort_c_num {x y } {

    regexp {.*#([0-9]+)\)$} $x dummy a
    regexp {.*#([0-9]+)\)$} $y dummy b

    if {$a < $b} {return -1}
    if {$a == $b} {return 0}
    if {$a > $b} {return 1}
}

##############################################################################
proc sort_c_len {x y } {

    set a [lindex $x 2]
    set b [lindex $y 2]

    if {$a < $b} {return 1}
    if {$a == $b} {return 0}
    if {$a > $b} {return -1}
}

##############################################################################
proc CreateCSMenu {io f csh_win csp_win csv_win} {
    global gap_defs selector_menu
    global $f.diagonal

    SetDefaultTags CONTIG_SEL.TAGS
    $f configure -menu $f.menubar
    menu $f.menubar
    create_menus $selector_menu $f.menubar

    menu_state_set selector_menu 12 $f.menubar
    ActivateCSMenu $f.menubar.view disabled
}

##############################################################################
proc SetContigRid { cs rid } {
    global $cs.rid

    set $cs.rid $rid
}

##############################################################################
proc GetContigRid { cs } {
    global $cs.rid

    return [set $cs.rid]
}

###############################################################################
# Callback when the global "contigs" list is modified
proc ContigSelector_ContigsList {io csh_win name1 name2 op} {
    global NGList

    if {![winfo exists $csh_win]} return

    DrawSelection $io $csh_win $NGList(contigs)
}

##############################################################################
#called from gaprc menu
proc ContigSelector { io } {
    global NGRec
    global gap_defs

    set f [keylget gap_defs CONTIG_SEL.WIN] 
    if {[xtoplevel $f] == ""} return
    fix_maxsize $f
    global $f.cs_id
    wm resizable $f 1 0

    #wm protocol $f WM_DELETE_WINDOW {puts "Please use the quit button!"}
    wm protocol $f WM_DELETE_WINDOW "CSStartShutdown $io $f"
    wm title $f "Contig Selector"

    set csh_win $f[keylget gap_defs CONTIG_SEL.WIN]
    set csp_win $f[keylget gap_defs CONTIG_CSP.WIN]
    set csv_win $f[keylget gap_defs CONTIG_CSV.WIN]
    set scroll x
    set borderwidth [keylget gap_defs CONTIG_SEL.BORDERWIDTH]
    set width [keylget gap_defs CONTIG_SEL.PLOT_WIDTH]
    set height [keylget gap_defs CONTIG_SEL.PLOT_HEIGHT]

    wm minsize $f 0 0

    set numcontigs [db_info num_contigs $io]
    set min_x 0
    set max_x [expr [db_info t_contig_length $io] + $numcontigs]

    ##########################################################################
    #create contig selector
    frame $f.csh -bd $borderwidth -relief groove
    #canvas $csh_win -width $width -height $height \
	    -bd 0 -highlightthickness 0\
	    -xscrollcommand "$f.hscroll set" \
	    -closeenough 4
    set zoom_cmd [list "gap_zoom $io \[set $f.cs_id\] $scroll"]
    canvasbox $csh_win -width $width -height $height \
	    -bd 0 -highlightthickness 0\
	    -xscrollcommand "$f.hscroll set" \
	    -closeenough 4 -zoom_command $zoom_cmd
	    
    scrollbar $f.hscroll -orient horizontal -relief sunken -command \
	    "gc_scroll_x $io \[set $f.cs_id\]"
    
    ##########################################################################
    # Main Menu Bar
    #frame $f.menubar -relief raised -borderwidth 2
    CreateCSMenu $io $f $csh_win $csp_win $csv_win

    ##########################################################################
    #button bar
    frame $f.buttons

    # Next button
    button $f.buttons.next -text "Next" -command "InvokeNext $f" -state disabled
    bind $f.buttons.next <Any-Enter> "+
	if {\[%W cget -state\] != \"disabled\"} {
	    ScanNext $f $f.brief $csh_win $csv_win $csp_win
	}"
    bind $f.buttons.next <Any-ButtonRelease> "+
	if {\[%W cget -state\] != \"disabled\"} {
	    ScanNext $f $f.brief $csh_win $csv_win $csp_win
	}"
    bind $f.buttons.next <Any-Leave> "+
	if {\[%W cget -state\] != \"disabled\"} {
	    UnHighlight $f $csp_win
	}"
     #zoom back button
    button $f.buttons.back -text "zoom out" -command "ZoomBackCanvas $io \[set $f.cs_id\]"
    button $f.buttons.zoomin10 -text "+10%" \
	-command "if {\[winfo exists $csp_win\]} {
			ZoomInCanvas $csp_win 0.05
		  } else {
			ZoomInCanvas $csh_win 0.05
		  }"
    button $f.buttons.zoomin50 -text "+50%" \
	-command "if {\[winfo exists $csp_win\]} {
			ZoomInCanvas $csp_win 0.1666
		  } else {
			ZoomInCanvas $csh_win 0.1666
		  }"
    
    #cursor checkbutton
    global $f.cursor
    checkbutton $f.buttons.cursor -text crosshairs -variable $f.cursor

    #cursor local and total position labels
    set cursor_tx [keylget gap_defs CONTIG_SEL.CURSOR1_X]
    set cursor_lx [keylget gap_defs CONTIG_SEL.CURSOR2_X]
    set cursor_ty [keylget gap_defs CONTIG_SEL.CURSOR1_Y]
    set cursor_ly [keylget gap_defs CONTIG_SEL.CURSOR2_Y]
    label $f$cursor_tx -bd 2 -relief sunken -width 6
    label $f$cursor_lx -bd 2 -relief sunken -width 6
    label $f$cursor_ty -bd 2 -relief sunken -width 6
    label $f$cursor_ly -bd 2 -relief sunken -width 6

    ##########################################################################
    label $f.brief_dummy
    label $f.brief

    pack $f.buttons.next $f.buttons.zoomin10 $f.buttons.zoomin50 \
	 $f.buttons.back $f.buttons.cursor -expand yes -side left
    pack $f$cursor_lx $f$cursor_tx $f$cursor_ly $f$cursor_ty -in $f.buttons -side left -expand yes 

    grid columnconfig $f 1 -weight 1

    #grid $f.menubar -row 0 -column 0 -sticky ew -columnspan 3
    grid $f.buttons -row 1 -column 0 -sticky ew -columnspan 3
    grid $f.csh     -row 2 -column 1 -sticky ew
    pack $csh_win -in $f.csh -padx 5 -fill x -expand yes
    grid $f.hscroll -row 4 -column 1 -sticky ew
    grid $f.brief_dummy -row 5 -column 0 -sticky ew -columnspan 3
    place $f.brief -in $f.brief_dummy -relx 0

    #need to ensure the windows are packed before doing the plot which uses
    #the canvas width and height to do the scaling
    tkwait visibility $csh_win

    set NGRec($f,prev_length) [db_info t_contig_length $io]
    set NGRec($f,last_result) 0

    set $f.cs_id [display_contig_selector -io $io -frame $f -window $csh_win]

    #bind the configure actions to the toplevel
    bind $f <Any-Configure> "
	    if {\[winfo toplevel %W\] == \"%W\"} {
		update idletasks
		resize_canvas -io $io -id [set $f.cs_id]
	    }
	"

    #SetCanvasBindings $io [set $f.cs_id]  $csh_win $scroll
    SetCanvasBindings $csh_win $zoom_cmd
    SetCSHBindings $io $f $csh_win $csp_win contig $scroll $f.brief \
	    $f$cursor_lx $f$cursor_tx

    # Trace changes in the global "contigs" list
    global NGList
    set trace_cmd "ContigSelector_ContigsList $io $csh_win"
    trace variable NGList(contigs) w $trace_cmd
}

#HACK - to check
proc ContigInitReg { io } {
    global do_csel

    if {[db_info num_contigs $io] > 0 && $do_csel} {
	ContigSelector $io
    }
}

#set contig globals after call to CS_reg may have changed contig lengths or
#deleted CurContig
proc ContigParams { io } {
    global CurContig
    global LREG
    global RREG

    if {[db_info get_contig_num $io $CurContig] != -1} {
	set LREG 1
	set RREG [c_length $io [db_info get_contig_num $io $CurContig]]
    } else {
	set longest [db_info longest_contig $io]
	set CurContig [left_gel $io $longest]
	set LREG
	set RREG [c_length $io $longest]
    }
}
