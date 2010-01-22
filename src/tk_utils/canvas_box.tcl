#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

proc itemNearest {plot x y} {
    set halo 200

    #convert into canvas coords
    set x [$plot canvasx $x]
    set y [$plot canvasy $y]

    set nearest [$plot find closest $x $y]

    if { $nearest != ""} {
	return $nearest
    }
    return 0
}

##############################################################################
proc FindZoomArea {canvas left right top bottom } {
    global $canvas.areaX1 $canvas.areaY1 $canvas.areaX2 $canvas.areaY2

    upvar $left z_left $right z_right $top z_top $bottom z_bottom

    #find left, top, right, bottom pixel coords of zoom rectangle
    if {[set $canvas.areaX1] < [set $canvas.areaX2]} {
	set z_left [set $canvas.areaX1]
	set z_right [set $canvas.areaX2]
    } else {
	set z_left [set $canvas.areaX2]
    	set z_right [set $canvas.areaX1]
    }
    if {[set $canvas.areaY1] < [set $canvas.areaY2]} {
	set z_top [set $canvas.areaY1]
	set z_bottom [set $canvas.areaY2]
    } else {
	set z_top [set $canvas.areaY2]
	set z_bottom [set $canvas.areaY1]
    }

    if {$z_left < [$canvas canvasx 0]} {
	set z_left [$canvas canvasx 0]
    }
    if {$z_right > [$canvas canvasx [winfo width $canvas]]} {
	set z_right [$canvas canvasx [winfo width $canvas]]
    }
}

##############################################################################
#tk_utils
proc ZoomCanvas {canvas zoom_cmd} {
    global $canvas.areaX1 $canvas.areaY1 $canvas.areaX2 $canvas.areaY2
    global $canvas.button $canvas.zoom

    #unset canvas.button upon button release
    if {![info exists $canvas.zoom]} {
	#delete any rectangle occurs with:
	#<B1-Press> <Motion> <B3-Press> <B3-Release> <B1-Release> => box left on screen.
	$canvas delete area
	return
    }

    if {[info exists $canvas.button]} {
	unset $canvas.button
    }

    if {![info exists $canvas.areaX1]} {
	return
    }

    if {![info exists $canvas.zoom]} {
	return
    }
    unset $canvas.zoom

   #minimum size of zooming rectangle
    set delta 10

    #if an area has not been marked out, ie areaX2 & areaY2 are 0
    if { (([set $canvas.areaX2] == 0) && ([set $canvas.areaY2] == 0)) || \
	  ([expr abs([set $canvas.areaX1]-[set $canvas.areaX2])] < $delta) ||\
	    ([expr abs([set $canvas.areaY1]-[set $canvas.areaY2])] < $delta)} {
	$canvas delete area
	ReSetZoomRect $canvas
	bell
	return
    }

    FindZoomArea $canvas z_left z_right z_top z_bottom
    
    #need to run the command at a global level
    uplevel #0 eval $zoom_cmd -x1 $z_left -y1 $z_top -x2 $z_right \
	-y2 $z_bottom

    ReSetZoomRect $canvas
}

proc ZoomInCanvas {canvas amount} {
    global $canvas.areaX1 $canvas.areaX2 $canvas.areaY1 $canvas.areaY2
    global $canvas.zoom

    set l $amount
    set r [expr 1.0-$amount]

    set w [winfo width $canvas]
    set h [winfo height $canvas]
    set $canvas.areaX1 [$canvas canvasx [expr ($w*$l) - 0]]
    set $canvas.areaX2 [$canvas canvasx [expr ($w*$r) - 0]]
    set $canvas.areaY1 [$canvas canvasy [expr $h*$l]]
    set $canvas.areaY2 [$canvas canvasy [expr $h*$r]]
    set $canvas.zoom 1
    eval ZoomCanvas $canvas [canvasbox_zoom_command $canvas]
}

##############################################################################
#set NGConstants
proc InitCanvasConstants { } {
    global CanvasConst

    set CanvasConst(xmargin) 10
    set CanvasConst(ymargin) 10

    #tick height
    set CanvasConst(tick_ht) 10

    #increment for autoscrolling. Larger num increase speed
    set CanvasConst(auto_incr) 10
    #time delay in ms for autoscrolling. Larger num decrease speed
    set CanvasConst(auto_time) 100
}

##############################################################################
proc RestoreCmd { c obj fill} {

    #puts "RestoreCmd $c $obj $fill"
    if {![winfo exists $c]} {
	return
    }

    $c itemconfig $obj -fill $fill
    $c lower $obj

    #raise item above the cursor 
    if {[$c find withtag cursor_x] != "" } {
	$c lower cursor_x
    }  
    if {[$c find withtag cursor_x1] != "" } {
	$c lower cursor_x1
    }  
    if {[$c find withtag cursor_y] != "" } {
	$c lower cursor_y
    }  
}

##############################################################################
#save the initial settings of an object
proc InitialSettings {f c obj type } {
    global restoreCmd initialCol

    set fill [$c itemcget $obj -fill] 

    #check that obj is valid
    if {$fill == ""} {
	set fill black
    }
    #set restoreCmd($f,$type) "$c itemconfig $obj -fill $fill; $c lower $obj"  
    set restoreCmd($f,$type) "RestoreCmd $c $obj $fill"
    set initialCol($f,$type) $fill
}

##############################################################################
proc UnHighlightItems { f } {
    global $f.prev_item
    global restoreCmd
    
    if {[info exists $f.prev_item]} {
	unset $f.prev_item
    }

    if {[info exists restoreCmd($f,reading)]} {
	eval $restoreCmd($f,reading)
    }
    if {[info exists restoreCmd($f,contig)]} {
	eval $restoreCmd($f,contig)
    }

}

##############################################################################
proc ReSetZoomRect {canvas } {
    global $canvas.areaX1 $canvas.areaY1 $canvas.areaX2 $canvas.areaY2

    $canvas delete area

    if [info exists $canvas.areaX1] {
        unset $canvas.areaX1 $canvas.areaX2 $canvas.areaY1 $canvas.areaY2
    }
}
#end ReSetZoomRect

##############################################################################
#tk_utils
proc SetCanvasBindings {canvas zoom_cmd} {
    global $canvas.zoom

#Disabled as it breaks things.
#    #scrolling using button 2
#    if {([string compare $scroll y] == 0) || \
#	    ([string compare $scroll b] == 0)} {
#	bind $canvas <<move>> "$canvas scan mark %x %y"
#	bind $canvas <<move-drag>> "$canvas scan dragto %x %y"
#    } elseif {([string compare $scroll x] == 0)} {
#	bind $canvas <<move>> "$canvas scan mark %x %y"
#	bind $canvas <<move-drag>> "$canvas scan dragto %x 0"
#    }

    #zoom area attached to mouse control button 3
    bind $canvas <<zoom>> "set $canvas.zoom 1; itemMark $canvas %x %y"
bind $canvas <<zoom-drag>> "if {\[info exists $canvas.zoom\]} {itemStroke $canvas %x %y}"
bind $canvas <<zoom-release>> "if {\[info exists $canvas.zoom\]} {ZoomCanvas $canvas $zoom_cmd}"

    bind $canvas <Escape> "%W delete area; ReSetZoomRect $canvas; bell"
}

proc itemsUnderArea {f canvas tag } {
    global $canvas.areaX1 $canvas.areaY1 $canvas.areaX2 $canvas.areaY2

    set list ""
    if {![info exists $canvas.areaX1]} {
	return
    }

    set area [$canvas find withtag area]

    #remove selecting rectangle
    $canvas delete area

    #foreach item in the enclosing box
    foreach i [$canvas find enclosed [set $canvas.areaX1] [set $canvas.areaY1]\
	    [set $canvas.areaX2] [set $canvas.areaY2]] {
	if {[lsearch [$canvas gettags $i] $tag] != -1} {
	    lappend list $i
	}
    }
    return $list
}

##############################################################################
proc itemMark {canvas x y} {
    global $canvas.areaX1 $canvas.areaY1 $canvas.areaX2 $canvas.areaY2
    global $canvas.button

    set $canvas.areaX1 [$canvas canvasx $x]
    set $canvas.areaY1 [$canvas canvasy $y]
    set $canvas.areaX2 [set $canvas.areaX1]
    set $canvas.areaY2 [set $canvas.areaY1]

    set $canvas.button 1
    $canvas delete area
}

##############################################################################
proc itemUnMark {canvas} {
    global $canvas.areaX1 $canvas.areaY1 $canvas.areaX2 $canvas.areaY2
    global $canvas.button

    if {[info exists $canvas.areaX1]} {
	unset $canvas.areaX1
    }
    if {[info exists $canvas.areaY1]} {
	unset $canvas.areaY1
    }
    if {[info exists $canvas.areaX2]} {
	unset $canvas.areaX2
    }
    if {[info exists $canvas.areaY2]} {
	unset $canvas.areaY2
    }
    if {[info exists $canvas.button]} {
	unset $canvas.button
    }
}

##############################################################################
proc itemStroke {canvas x y} {
    global $canvas.areaX1 $canvas.areaY1 $canvas.areaX2 $canvas.areaY2
    global $canvas.button
    global NGRec
    global NGConst

    if {![info exists $canvas.areaX1]} {
	return
    }
    #only stroke if previously set canvas.button. This is unset upon release
    if {![info exists $canvas.button]} {
	unset $canvas.areaX1 $canvas.areaY1 $canvas.areaX2 $canvas.areaY2
	return
    }

    set x [$canvas canvasx $x]
    set y [$canvas canvasy $y]

    if {([set $canvas.areaX1] != $x) && ([set $canvas.areaY1] != $y)} {
	$canvas delete area
	$canvas addtag area withtag [$canvas create rect [set $canvas.areaX1] \
		[set $canvas.areaY1] $x $y\
		-outline black]

	set $canvas.areaX2 $x
	set $canvas.areaY2 $y
    }
}

##############################################################################
#remove borders and highlight thickness from winfo 
proc GetCanvasWidth { canvas } {

    set width [winfo width $canvas]
    set highlight [$canvas cget -highlightthickness]
    set bd [$canvas cget -bd]

    #puts "GetCanvasWidth $canvas $width $highlight $bd"
    return [expr $width - (2 * $highlight) - (2 * $bd)]
}
##############################################################################
#remove borders and highlight thickness from winfo 
proc GetCanvasHeight { canvas } {

    set height [winfo height $canvas]
    set highlight [$canvas cget -highlightthickness]
    set bd [$canvas cget -bd]
    #puts "GetCanvasHeight $canvas $height $highlight $bd"
    return [expr $height - (2 * $highlight) - (2 * $bd)]
}


##############################################################################
proc ZoomBackCanvas {io id} {

    zoom_canvas -io $io -id $id

}

##############################################################################
#called from C
proc DrawCanvasCursorXOLD {f canvas cx colour line_width} {
    global gap_defs $f.cursor_obj

    set height [$canvas canvasy [winfo height $canvas]]
    if {![info exists $f.cursor_obj]} {
	
	set $f.cursor_obj [$canvas create line $cx [$canvas canvasy 0] $cx $height -tag cursor_x\
		-fill $colour -width $line_width]
	$canvas lower [set $f.cursor_obj]
    } else {
	$canvas coords [set $f.cursor_obj] $cx [$canvas canvasy 0] $cx $height
    }
}

##############################################################################
#called from C
proc DrawCanvasCursorX {f canvas cx colour line_width} {
    global gap_defs 

    set height [$canvas canvasy [winfo height $canvas]]
    if {[$canvas find withtag cursor_x] == ""} {
	$canvas create line $cx [$canvas canvasy 0] $cx $height -tag cursor_x\
		-fill $colour -width $line_width
	$canvas lower cursor_x
    } else {
	$canvas coords cursor_x $cx [$canvas canvasy 0] $cx $height
    }
}

##############################################################################
#called from C
proc DrawCanvasCursorY {f canvas cy colour line_width} {
    global gap_defs

    set width [$canvas canvasx [winfo width $canvas]]
    if {[$canvas find withtag cursor_y] == ""} {
	$canvas create line [$canvas canvasx 0] $cy $width $cy -tag cursor_y\
		-fill $colour -width $line_width
	$canvas lower cursor_y
    } else {
	$canvas coords cursor_y [$canvas canvasx 0] $cy $width $cy
    }
}

##############################################################################
#called from C - contig selector - display y cursor as x cursor
proc DrawCanvasCursorX1 {f canvas cx colour line_width} {
    global gap_defs 

    set height [$canvas canvasy [winfo height $canvas]]
    if {[$canvas find withtag cursor_x1] == ""} {
	$canvas create line $cx [$canvas canvasy 0] $cx $height -tag cursor_x1\
		-fill $colour -width $line_width
	$canvas lower cursor_x1
    } else {
	$canvas coords cursor_x1 $cx [$canvas canvasy 0] $cx $height
    }
}


proc canvasbox_configure {path args} {
    global $path.zoom_cmd $path.zoomback_cmd
    set in_arg 0
    set arglist ""
    set $path.zoom_cmd ""
    set $path.zoomback_cmd ""

    # Process command line args
    foreach i $args {
	if {$in_arg} {

	     if {$option == "-zoom_command"} {
		set $path.zoom_cmd "{$i}"
	     } elseif {$option == "-zoomback_command"} {
		 set $path.zoomback_cmd "{$i}"
	     } else {
		lappend arglist $option $i
	     }
	     set in_arg 0
	} else {
	     set option $i
	     set in_arg 1
	}
    }
    eval $path configure $arglist

    if {[set $path.zoom_cmd] != ""} {
	#set bindings
	SetCanvasBindings $path [set $path.zoom_cmd]
    }
}

proc canvasbox {path args } {
    global $path.zoom_cmd
    # Create the canvas
    canvas $path
    
    # Configure	
    eval canvasbox_configure $path $args

    InitCanvasConstants
    return $path
}

proc canvasbox_zoom_command {path} {
    global $path.zoom_cmd

    return [set $path.zoom_cmd]
}

proc canvasbox_zoomback_command {path} {
    global $path.zoomback_cmd

    return [set $path.zoomback_cmd]
}
