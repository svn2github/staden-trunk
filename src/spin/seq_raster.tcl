#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
##############################################################################
##############################################################################
#essential global variables for use with the raster
proc RasterGlobals { } {
    global ZOOM_BOX SCALE_BAR ZOOM_SCALE

    #raster types
    #zoom box only
    set ZOOM_BOX  1   

    #scale bars only
    set SCALE_BAR 2

    #both zoom and scale
    set ZOOM_SCALE 3
}

proc RasterEnter {raster x y x_format y_format args} {
    #puts ENTER 

    set r_win [winfo parent $raster]
    global $r_win.cross $r_win.onscreen

    set $r_win.onscreen 1
    if {[set $r_win.cross]} {
	draw_raster_xh $raster $x $y $x_format $y_format
    }
    
}

proc RasterLeave {raster} {

    #puts LEAVE 
    set r_win [winfo parent $raster]
    global $r_win.cross $r_win.onscreen

    if {[set $r_win.cross]} {
	remove_raster_xh $raster
    }
    set $r_win.onscreen 0
}

proc RasterMove {raster x y x_format y_format} {

    #puts MOVE
    set r_win [winfo parent $raster]
    global $r_win.cross $r_win.onscreen
    if {[set $r_win.cross] && [set $r_win.onscreen]} {
	remove_raster_xh $raster
	draw_raster_xh $raster $x $y $x_format $y_format
    }
}

#Raster class bindings
proc RasterBindings {x_format y_format} {

    bind Raster <Any-Motion> "RasterMove %W %x %y \$x_format \$y_format"
    bind Raster <Any-Enter>  "RasterEnter %W %x %y \$x_format \$y_format"
    bind Raster <Any-Leave>  "RasterLeave %W"

    #ability to move crosshairs off plot but maintain their position
    bind Raster <Control-Any-Motion> {break;}
}

#draw cross-hairs on raster display
proc raster_display_xh { r_win stem raster} {
    global $r_win.cross

    if {[set $r_win.cross]} {
	raster_xhinit $r_win $raster $stem
    } else {
	raster_xhdel $r_win $raster $stem
    }
}

##############################################################################
#set up bindings for cross-hairs on raster display
proc raster_xhinit {r_win raster stem} {
    global $raster.Lastx $r_win.Lasty
    global $r_win.win_list
    global $r_win.cross $raster.visible

    set $raster.visible 0
    set cursor_env [$raster envcreate -function xor -fg #d9d9d9]
    set cnt 0

    foreach i [set $r_win.win_list] {
	if {[string compare $stem$i $raster] == 0} {
	    lappend i $cursor_env
	    set $r_win.win_list [lreplace [set $r_win.win_list] $cnt $cnt $i]
	}
	incr cnt
    }
    set $raster.Lastx ""
    #set $r_win.Lasty ""
}

##############################################################################
#delete cross-hairs and bindings
proc raster_xhdel {r_win raster stem} {
    global $raster.Lastx $r_win.Lasty $r_win.win_list $raster.visible

    if {[winfo exists $r_win.buttons.pos1]&&
	[winfo exists $r_win.buttons.pos2]} {
	$r_win.buttons.pos1 configure -text ""
	$r_win.buttons.pos2 configure -text ""
    }
    if {[info exists $raster.Lastx]} {unset $raster.Lastx}
    if {[info exists $r_win.Lasty]} {unset $r_win.Lasty}
    if {[info exists $raster.visible]} {unset $raster.visible}
}


##############################################################################
#draw cross-hairs on raster display
proc draw_raster_xh {raster x y x_format y_format} {

    set r_win [winfo parent $raster]
    set stem [GetRasterStem $r_win]

    global $raster.visible $r_win.win_list $r_win.cross
    global $raster.Lastx $r_win.Lasty

    if {[set $r_win.cross] == 0} {
	return
    }

    #if cursor is still visible, don't draw
    if {[set $raster.visible] == 1} {
	return
    }
    #puts "raster_xh $r_win $stem x $x lastx [set $raster.Lastx] winlist [set $r_win.win_list]"
   
    #puts "config [$raster envconfigure]"
    set b [$raster cget -bd]
    set wid [expr [winfo width $raster]-2*$b]
    set x [expr $x-$b]
    set y [expr $y-$b]

    set r_id [GetRasterId $raster]

    if {[winfo exists $r_win.buttons.pos1]&&[winfo exists $r_win.buttons.pos2]} {
	rasterPositionLabel $r_win.buttons.pos1 $r_win.buttons.pos2 $raster $x $y $x_format $y_format
    }
    foreach win [set $r_win.win_list] {
	    
	set i [lindex $win 0]
	set s [lindex $win 1]
	
	if {$r_id == $i} {
	    set style $s
	}
	set hei [expr [winfo height $stem$i]-2*$b]
	
	eval $stem$i draw_line [$stem$i toworld $x 0] [$stem$i toworld $x $hei] -style $s
    }
    eval $raster draw_line [$raster toworld 0 $y] [$raster toworld $wid $y] -style $style
	
    update idletasks
    set $raster.Lastx $x
    set $r_win.Lasty $y
    set $raster.visible 1
}

##############################################################################
#remove cross-hairs on raster display
proc remove_raster_xh {raster} {

    set r_win [winfo parent $raster]
    set stem [GetRasterStem $r_win]

    global $raster.visible $raster.Lastx $r_win.Lasty $r_win.win_list $r_win.cross

    if {[set $r_win.cross] == 0} {
	return
    }

    #if raster is not visible then don't remove it
    if {[set $raster.visible] == 0} {
	return
    }
    #puts "raster_xh $r_win $stem lastx [set $raster.Lastx] winlist [set $r_win.win_list]"
   
    #puts "config [$raster envconfigure]"
    set b [$raster cget -bd]
    set wid [expr [winfo width $raster]-2*$b]
    #set x [expr $x-$b]
    #set y [expr $y-$b]

    set r_id [GetRasterId $raster]

    foreach win [set $r_win.win_list] {

	set i [lindex $win 0]
	set s [lindex $win 1]
	
	if {$r_id == $i} {
	    set style $s
	    #puts "r_id $r_id i $i style $s"
	}
	set hei [expr [winfo height $stem$i]-2*$b]
	
	eval $stem$i draw_line [$stem$i toworld [set $raster.Lastx] 0]\
	    [$stem$i toworld [set $raster.Lastx] $hei] -style $s
    }
    eval $raster draw_line [$raster toworld 0 [set $r_win.Lasty]] \
	[$raster toworld $wid [set $r_win.Lasty]] -style $style
    set $raster.visible 0
}

##############################################################################
proc rasterPositionLabel {labelx labely w x y x_format y_format} {

    set position [rasterGetPos $w $x $y]
    #$labelx configure -text [expr int([lindex $position 0])]
    eval scan [expr [lindex $position 0] + 0.5] $x_format pos
    $labelx configure -text $pos
    #eval scan [expr [lindex $position 1] + 0.5] $y_format pos
    eval scan [lindex $position 1] $y_format pos
    $labely configure -text $pos
}

##############################################################################
#convert from raster top = low bottom = high to top = high bottom = low
proc rasterY { w y } {

    set min_y [lindex [$w world_size] 1]
    set max_y [lindex [$w world_size] 3]

    return [expr $max_y - $y + $min_y]

}

##############################################################################
#convert from raster to world for entering into the crosshair label
proc rasterGetPos {w rx ry} {

    set world [$w toworld $rx $ry]
    set wx [lindex $world 0]
    set wy [rasterY $w [lindex $world 1]]

    set position "$wx $wy"
    #puts "************rasterGetPos $position"

    return $position
}

##############################################################################
proc RasterStartShutdown {r_win} {
    global $r_win.win_list $r_win.cross 

    set $r_win.cross 0
    foreach win [set $r_win.win_list] {
	#puts "RasterStartShutdown win $win"
	
	set id [lindex $win 0]
	seq_result_update -index $id -job QUIT
    }
    destroy $r_win
}

##############################################################################
#callback from DELETE or QUIT from c
proc removeRaster { raster manager} {

    set r_win [winfo parent $raster]
    global $r_win.win_list $r_win.row $r_win.first_row

    #puts "REMOVE_RASTER $raster $manager"

    set id [GetRasterId $raster]
    seq_result_list_update $manager

    set cnt 0
    #remove raster from win_list
    foreach win [set $r_win.win_list] {
	set w_id [lindex $win 0]
	if {$id == $w_id} {
	    set $r_win.win_list [lreplace [set $r_win.win_list] $cnt $cnt]
	}
	incr cnt
    }

    #must get info on anything except $raster which is packed -in
    array set info [grid info $r_win.b$id]
    grid_delete [winfo parent $raster] row $info(-row) 1
    incr $r_win.row -1

    #need to set the weight of the deleted row to 0
    grid rowconfigure $r_win [set $r_win.row] -weight 0 -minsize 0

    destroy $r_win.label$id
    destroy $r_win.ruler_v$id
    destroy $r_win.sb_v$id
    destroy $r_win.key$id
    destroy $raster
    destroy $r_win.b$id

    #remove toplevel if no more slaves left
    if {[set $r_win.row] == [set $r_win.first_row]} {
	destroy $r_win
    }
}

##############################################################################
proc rasterResultsManager {help_cmd } {
    global tk_utils_defs

    set t [keylget tk_utils_defs RASTER.RESULTS.WIN]

    if {[xtoplevel $t] == ""} return
    wm title $t "Results manager"

    seq_result_list_create $t $help_cmd
    seq_result_list_update $t

}

##############################################################################
#hide all results in dot plot
proc rasterHide {raster } {

    $raster clear
    set id [GetRasterId $raster]
    seq_result_update -index $id -job HIDE
}

##############################################################################
#reveal all results in a raster
proc rasterReveal {raster } {

   # $raster clear
    set id [GetRasterId $raster]

    #puts "ID $id"
    seq_result_update -index $id -job REVEAL

    raster_replot_all -raster $raster
}

proc RemoveWin {r_win raster_id} {
    global $r_win.win_list

    #puts "REMOVEWIN"

    set cnt 0
    #remove raster from win_list
    foreach win [set $r_win.win_list] {
	set w_id [lindex $win 0]
	if {$raster_id == $w_id} {
	    set $r_win.win_list [lreplace [set $r_win.win_list] $cnt $cnt]
	}
	incr cnt
    }
}

##############################################################################
proc init_crosshair {r_win stem} {
    global $r_win.win_list 
    
    foreach win [set $r_win.win_list] {
	set i [lindex $win 0]
	raster_display_xh $r_win $stem $stem$i
    }
}

proc GetRasterWinList {r_win} {
    #set r_win [winfo parent $raster]
    global $r_win.win_list

    set list ""
    set raster_stem [GetRasterStem $r_win]

    foreach win [set $r_win.win_list] {
	set id [lindex $win 0]
	lappend list $raster_stem$id
    }
    return $list
}

proc GetRasterIdList {r_win} {

    #set r_win [winfo parent $raster]
    global $r_win.win_list

    set list ""
    foreach win [set $r_win.win_list] {
	set id [lindex $win 0]
	lappend list $id
    }
    return $list
}

proc GetRasterWindowSize { raster} {

    set r_win [winfo parent $raster]
    set raster_stem [GetRasterStem $r_win]

    global $r_win.win_list

    set x0 100000000
    set y0 100000000
    set x1 0
    set y1 0

    foreach win [set $r_win.win_list] {
	set id [lindex $win 0]
	
	set world_size [$raster_stem$id world_size]
	set xx0 [lindex $world_size 0]
	set yy0 [lindex $world_size 1]
	set xx1 [lindex $world_size 2]
	set yy1 [lindex $world_size 3]
	if {$xx0 < $x0} {
	    set x0 $xx0
	}
	if {$yy0 < $y0} {
	    set y0 $yy0
	}
	if {$xx1 > $x1} {
	    set x1 $xx1
	}
	if {$yy1 > $y1} {
	    set y1 $yy1
	}
    }
    return "$x0 $y0 $x1 $y1"
}


##############################################################################
#plot horizontal ruler and ticks
proc rasterHRuler {raster min_x max_x} {

    #puts "rasterHRuler $raster $min_x $max_x"

    #HACK do I need this?
    #update idletasks

    set r_win [winfo parent $raster]
    set r_width [winfo width $raster]
    set p0 0
    set p1 [expr $p0 + $r_width]
    set ruler_offset 10

    $r_win.ruler_h delete all
    $r_win.ruler_h create line $p0 $ruler_offset $p1 $ruler_offset \
	-tags ruler_h

    set num_ticks 4

    set ticks [ruler_ticks -min $min_x -max $max_x -num_ticks $num_ticks]
    eval PlotTicks_h $raster $r_win.ruler_h $ticks
}

##############################################################################
#plot vertical ruler and ticks
#called from C and tcl
proc rasterVRuler {raster min_y max_y} {

    #puts "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^rasterVRuler $raster $min_y $max_y"

    set zoom [raster_results -id [GetRasterId $raster] -option zoom]

    if {$zoom <= 0} {
	return
    }

    set r_win [winfo parent $raster]
    set id [GetRasterId $raster]

    set r_height [winfo height $raster]
    set p0 0
    set p1 [expr $p0 + $r_height]

    set ruler_offset [expr [winfo width $r_win.ruler_v$id] - 10]

    $r_win.ruler_v$id delete all
    $r_win.ruler_v$id create line $ruler_offset $p0 $ruler_offset $p1

    set num_ticks 4

    #puts "-min $min_y -max $max_y"
    set ticks [ruler_ticks -min $min_y -max $max_y -num_ticks $num_ticks]
    #puts "tick $ticks"

    eval PlotTicks_v $raster $r_win.ruler_v$id $ticks
}

##############################################################################
proc raster_setvars {r_win stem raster} {
    global $r_win.cross

    if {![info exists $r_win.cross]} {
	set $r_win.cross 0
    }
    raster_display_xh $r_win $stem $raster

}

##############################################################################
#called by interactive change in window size
proc RescaleRaster {raster} {

    #puts "start RescaleRaster $raster"
    rasterHRuler $raster [lindex [$raster world_size] 0] [lindex [$raster world_size] 2]
    rasterVRuler $raster [lindex [$raster world_size] 1] [lindex [$raster world_size] 3]

    #regrid the keys
    set r_win [winfo parent $raster]
    set r_id [GetRasterId $raster]

    keybox_grid $r_win.key$r_id
}

##############################################################################
proc scrollXCmd { r_win raster ruler_h args} {
    
    #puts "start scrollXCmd $raster $args"

    global $r_win.win_list
    
    #check if win_list contains any items
    if {[llength [set $r_win.win_list]] == 0} {
	return
    }
    foreach win [set $r_win.win_list] {
	set i [lindex $win 0]
	eval $raster$i xview $args
    }
    #HACK - the above command doesn't update "world" unless I do this
    update idletasks
    rulerTick_h $raster[lindex [lindex [set $r_win.win_list] 0] 0] $ruler_h
}

##############################################################################
proc scrollYCmd {raster ruler_v args} {

    #puts "start scrollYCmd"
    eval $raster yview $args
    
    #HACK - the above command doesn't update "world" unless I do this
    update idletasks
    rulerTick_v $raster $ruler_v
}

##############################################################################
#plot horizontal ruler ticks on th fly ie when scrolling or zooming
proc rulerTick_h {raster ruler_h} {

    #puts "start rulerTick_h $raster"

    set world [$raster world]
    set wx0 [lindex $world 0]
    set wx1 [lindex $world 2]
    #puts "-min $wx0 -max $wx1"

    set num_ticks 4
    set ticks [ruler_ticks -min $wx0 -max $wx1 -num_ticks $num_ticks]    
    #puts "ticks $ticks"

    eval PlotTicks_h $raster $ruler_h $ticks
}

##############################################################################
#plot vertical ruler ticks on the fly ie when scrolling or zooming
proc rulerTick_v {raster ruler_v} {

    set zoom [raster_results -id [GetRasterId $raster] -option zoom]
    #if no zooming or unknown
    if {$zoom <= 0} {
	return
    }

    set world [$raster world]
    set wy0 [rasterY $raster [lindex $world 3]]
    set wy1 [rasterY $raster [lindex $world 1]]
    set num_ticks 4
    #puts "^^^^^^^^^^^^^^^^rulerTick_v $raster -min $wy0 -max $wy1"

    set ticks [ruler_ticks -min $wy0 -max $wy1 -num_ticks $num_ticks]
    #puts "***********ticks $ticks"

    eval PlotTicks_v $raster $ruler_v $ticks
}

##############################################################################
#draw horizontal ruler ticks
proc PlotTicks_h {raster ruler firstTick step numTicks } {
    global $ruler.step

    set $ruler.step $step
    set ruler_offset 10
    set text_offset 25
    set tick_ht 10
    
    $ruler delete tick
    
    set pos [lindex [$raster topixmap $firstTick 0] 0]
    set tick $firstTick
    
    $ruler create line $pos $ruler_offset $pos \
	    [expr $ruler_offset + $tick_ht] -tag tick
    $ruler create text $pos $text_offset\
	    -text [format %g $tick] -tag tick
    
    set cstep [expr [lindex [$raster topixmap $step 0] 0] - \
	    [lindex [$raster topixmap 0 0] 0]]
    
    for {set i 1} {$i <= $numTicks} {incr i} {
	set pos [expr $pos + $cstep]
	set tick [expr $tick + $step]
	$ruler create line $pos $ruler_offset $pos \
		[expr $ruler_offset + $tick_ht] -tag tick
	$ruler create text $pos $text_offset\
		-text [format %g $tick] -tag tick
    }
}

##############################################################################
#draw vertical ruler ticks
proc PlotTicks_v {raster ruler firstTick step numTicks } {
    global $ruler.step 

    #puts PlotTicks_v 

    set $ruler.step $step
    set ruler_offset [expr [winfo width $ruler] - 10]
    set text_offset [expr [winfo width $ruler] - 40]

    set tick_ht 10

    $ruler delete tick

    set min_y [lindex [$raster world_size] 1]
    set max_y [lindex [$raster world_size] 3]

    #puts "raster $raster $min_y $max_y"
    #convert into raster world ie top = low, bottom = high
    set pos [lindex [$raster topixmap 0 [expr $max_y - $firstTick + $min_y]] 1]

    #puts "PlotTicks_v firsttick $firstTick pos $pos max $max_y"

    set tick $firstTick
    $ruler create line $ruler_offset $pos \
	    [expr $ruler_offset - $tick_ht] $pos -tag tick
    $ruler create text $text_offset $pos -text $tick -tag tick

    set cstep [expr [lindex [$raster topixmap 0 $step] 1] - \
	    [lindex [$raster topixmap 0 0] 1]]

    for {set i 1} {$i <= $numTicks} {incr i} {
	set pos [expr $pos - $cstep]
	set tick [format "%g" [expr $tick + $step]]
	$ruler create line $ruler_offset $pos \
		[expr $ruler_offset - $tick_ht] $pos -tag tick
	$ruler create text $text_offset $pos -text $tick -tag tick
    }
}

##############################################################################
proc RasterScaleX {r_win stem ruler_h value} {
    global $r_win.win_list

    #puts "RasterScaleX $r_win $value"

    foreach win [set $r_win.win_list] {
	set i [lindex $win 0]
	$stem$i xmag $value
    }

    rulerTick_h $stem[lindex [lindex [set $r_win.win_list] 0] 0] $ruler_h 
}

##############################################################################
proc RasterScaleY {r_win stem ruler_v value} {
    global $r_win.win_list

    #puts "RasterScaleY $value"

    foreach win [set $r_win.win_list] {
	set i [lindex $win 0]
	$stem$i ymag $value
	rulerTick_v $stem$i $ruler_v$i
    }
}

proc InitRasterColours {raster } {
    global spin_defs

    $raster envcreate -fg [keylget spin_defs RASTER.COLOUR.0]
    $raster envcreate -fg [keylget spin_defs RASTER.COLOUR.1]
    $raster envcreate -fg [keylget spin_defs RASTER.COLOUR.2]
    $raster envcreate -fg [keylget spin_defs RASTER.COLOUR.3]
    $raster envcreate -fg [keylget spin_defs RASTER.COLOUR.4]
    $raster envcreate -fg [keylget spin_defs RASTER.COLOUR.5]
    $raster envcreate -fg [keylget spin_defs RASTER.COLOUR.6]
    $raster envcreate -fg [keylget spin_defs RASTER.COLOUR.7]
    $raster envcreate -fg [keylget spin_defs RASTER.COLOUR.8]
    $raster envcreate -fg [keylget spin_defs RASTER.COLOUR.9]

}

##############################################################################
#need to regrid the x and y scalebars if inserted a new row at first position
proc MoveScaleBars {r_win row} {
    global $r_win.first_row

puts MoveScaleBars
return

    if {[set $r_win.first_row] == $row} {
	array set info [grid info $r_win.scale_x]
	set info(-row) [set $r_win.first_row]
	eval grid $r_win.scale_x [array get info]

	array set info [grid info $r_win.scale_y]
	set info(-row) [set $r_win.first_row]
	eval grid $r_win.scale_y [array get info]
    }
}

##############################################################################
#create raster with associated scrollbar and ruler
proc RasterFrame {r_win raster id raster_height raster_width ruler_width row r_type} {
    global $r_win.row $r_win.win_list $r_win.first_row 
    global $r_win.borderwidth
    global SCALE_BAR ZOOM_SCALE

    #puts "***************RasterFrame $r_win $raster row $row width $ruler_width"
    
    wm geometry $r_win {}

    set new_height [expr [winfo height $r_win] + $raster_height + \
			2 * [set $r_win.borderwidth]]

    #puts "height $raster_height new $new_height"
    #puts "width [winfo width $r_win]"

    frame $r_win.b$id -bd [set $r_win.borderwidth] -relief groove
    raster $raster -height $raster_height -width $raster_width \
	-bd 0 -dbl 0\
	-xscrollcommand "$r_win.sb_h set " \
	-yscrollcommand "$r_win.sb_v$id set " -scrollincrement 10

    label $r_win.label$id -text "p \#$id"

    InitRasterColours $raster

    canvas $r_win.ruler_v$id -width $ruler_width -height $raster_height \
	-highlightthickness 0 

    scrollbar $r_win.sb_v$id -orient vertical \
	-command "scrollYCmd $raster $r_win.ruler_v$id"

    keybox $r_win.key$id

    #add key for all results 
    keybox_add $r_win.key$id -text all \
	-background white \
	-motion MotionRaster \
	-enter "EnterKey $raster -1" \
	-leave "LeaveKey $raster" \
	-drop "DropRaster $raster $r_type" \
	-menu "seq_result_keybox_update $r_win -1 \[seq_result_names -raster_id $id\]"
    fit_on_screen $r_win

    bind $raster <Any-Configure> "RescaleRaster $raster"

    SetRasterBinding $raster 

    set raster_stem [GetRasterStem $r_win]

    grid_insert $r_win row $row 1
    grid rowconfig $r_win $row -weight 1 -minsize 40

    #grid rowconfig $r_win.b$id $row -weight 1 -minsize 40
    #grid columnconfig $r_win.b$id 3 -weight 1 -minsize 40

    grid $r_win.b$id -row $row -column 4 -sticky nsew

    grid $r_win.sb_v$id -row $row -column 5 -sticky ns 
    #grid $raster -in $r_win.b$id -row $row -column 3 -sticky nsew
    pack $raster -in $r_win.b$id -fill both -expand 1

    grid $r_win.label$id -row $row -column 2 -sticky n
    grid $r_win.ruler_v$id -row $row -column 3 -sticky ns
    grid $r_win.key$id -row $row -column 6 -sticky nsew

    #update idletasks
    #wm geometry $r_win [winfo width $r_win]x${new_height}

    #need to regrid the x and y scalebars if inserted a new row at first position
    #if {$r_type == $SCALE_BAR || $r_type == $ZOOM_SCALE} {
	#MoveScaleBars $r_win $row
    #}
    incr $r_win.row
    lappend $r_win.win_list $id
    raster_setvars $r_win $raster_stem $raster
}

##############################################################################
proc AddRaster {r_win raster raster_height raster_width ruler_width row seq_id r_type} {
    global raster_cmds raster_reg_cmd tk_utils_defs
    global $r_win.row $r_win.win_list $r_win.first_row
    global ZOOM_BOX ZOOM_SCALE

    #register raster with each sequence, irrespective of whether a new
    #raster widget is created
    #puts [info level [info level]]
    #puts "ADDRASTER $seq_id"
    
    if {[info exists raster_reg_cmd]} {
	set id [eval $raster_reg_cmd]
    } else {
	verror ERR_WARN "AddRaster" "A registration command must be set"
    }

    set raster_stem $raster
    set raster $raster$id
    if [winfo exists $raster] { raise $raster; return $id}

    frame $r_win.right.$id
    RasterFrame $r_win $raster $id $raster_height $raster_width $ruler_width $row $r_type

    #module specific commands eg nip4 bindings
    if {[info exists raster_cmds]} {
	eval $raster_cmds
    }
    seq_result_list_update [keylget tk_utils_defs RASTER.RESULTS.WIN]

    #only do this if zooming is applicable using a zoom box eg sip4 
    if {$r_type == $ZOOM_BOX || $r_type == $ZOOM_SCALE} {
	raster_zoominit $r_win $raster_stem $raster
	$r_win.buttons.back configure -command "rasterZoomUndo $r_win $raster_stem $raster"
    }
    fit_on_screen $r_win
    return $id
}

proc UpdateRasterPlot {index colour} {

} 

proc ConfigureRasterPlot {result_id line_width colour} {
    raster_config -index $result_id -width $line_width -fill $colour
    ConfigKey $result_id $colour
    
} 

proc ConfigureRasterWidth {index line_width} {
    raster_config -index $index -width $line_width
} 

##############################################################################
#called from C
proc RasterConfig {result_id} {

    if {[xtoplevel [set f .cbox] -resizable 0] == ""} return

    set line_width [lindex [lindex [raster_getconfig -index $result_id] 1] 1]

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

    #starting colour
    set colour [lindex [lindex [raster_getconfig -index $result_id] 0] 1]

    #cmd to execute when ok button on colourbox pressed
    set ok_cmd "ConfigureRasterPlot $result_id \[$f.lw.scale get\]"

    #cmd to execute when changing colours on colourbox
    set update_cmd "UpdateRasterPlot $result_id"

    #cmd to execute when cancel button on colourbox pressed
    set cancel_cmd "set dummy 0"

    ColourBox $f $colour $ok_cmd $update_cmd $cancel_cmd
}

##############################################################################
proc DropResult {id r_type x y} {

    after 50 [list WhereResult $id $r_type $x $y]
}

##############################################################################
proc DropRaster {id r_type x y} {

    after 50 [list WhereRaster $id $r_type $x $y]
}

##############################################################################
#find position after dropped a result
proc WhereResult {id r_type x y} {
    global r_from r_to w_to prev_raster

    if {![info exists r_from] || $r_from == ""} {
	return
    }

    # only allow placement onto other Rasters or new window
    if {$r_to != "" && [string compare [winfo class $r_to] "Raster"] != 0} {
	return
    }

    #puts "WhereResult from $r_from to $r_to win $w_to"
    RemoveHighlightPosition
    catch {unset prev_raster}

    if {$r_to != ""} {
	set row [FindResultPosition [winfo parent $r_to] $x $y 0]
	if {$row == -999} {
	    return
	} elseif {$row == -1} {
	    #put on top
	    RasterLeave $r_to
	    MoveRasterResult_Superimpose $id $r_type $r_from $r_to $x $y SUPER
	    RasterEnter $r_to $x $y %d %6f
	} else {
	    #create new raster for result
	    RasterLeave $r_to
	    set r_to [MoveRasterResult_New $id $r_type $r_from $w_to $x $y $row]
	    RasterEnter $r_to $x $y %d %6f
	}

    } else {
	#create a new window for result
 	InitResultWindow $r_from $id $r_type
     }
    set r_from ""
    set r_to ""
    set w_to ""
}

##############################################################################
#find position after dropped a raster
proc WhereRaster {id r_type x y} {
    global r_from r_to w_to prev_raster

    if {![info exists r_from] || $r_from == ""} {
	return
    }

    # only allow placement onto other Rasters or new window
    if {$r_to != "" && [string compare [winfo class $r_to] "Raster"] != 0} {
	return
    }

    RemoveHighlightPosition
    catch {unset prev_raster}

    if {$r_to != ""} {
	set row [FindResultPosition [winfo parent $r_to] $x $y 1]
	if {$row == -999} {
	    return
	} elseif {$row == -1} {
	    #superimpose two rasters
	    RasterLeave $r_to
	    MoveRaster_Superimpose $r_from $r_to $x $y
	    RasterEnter $r_to $x $y %d %6f
	} elseif {[winfo parent $r_from] == [winfo parent $r_to]} {
	    #move raster within the same window
	    RasterLeave $r_to
	    MoveRaster_Same $r_from $r_to $x $y $row $r_type
	    RasterEnter $r_to $x $y %d %6f
	} else {
	    #move raster to a different window
	    RasterLeave $r_to
	    MoveRaster_Different $r_type $r_from $r_to $x $y $row
	    RasterEnter $r_to $x $y %d %6f
	}
    } else {
	#create new raster window
	InitRasterWindow $r_from $r_type
    }
    set r_from ""
    set r_to ""
}

##############################################################################
proc EnterKey {raster result_id} {

    if {$result_id == -1} {
	#if all - no status line
	return 
    }

    set r_win [winfo parent $raster]
    $r_win.status_line config -text [seq_get_brief -index $result_id]

}


##############################################################################
proc LeaveKey {raster} {
    global r_from r_to 

    set r_from $raster
    set r_to ""
    #puts "LEAVE $r_from"

}

##############################################################################
proc EnterRaster {raster} {
    global r_to
    set r_to $raster

    #puts "ENTER $r_to"
}

##############################################################################
proc MotionRaster {X Y} {
    global r_to
    set r_to [winfo containing $X $Y]
    HighlightPosition [winfo containing $X $Y] $Y
}

##############################################################################
#if leave window frame, set w_to to ""
proc LeaveRWin {r_win} {
    global w_to
    
    set w_to ""
    #puts "LEAVE RWIN $w_to"
}

##############################################################################
proc EnterRWin {r_win} {
    global w_to

    set w_to $r_win
    #puts "ENTER RWIN $w_to"
}

##############################################################################
proc RemoveHighlightPosition {} {
    global prev_style prev_pos prev_raster

    if {![info exists prev_raster]} {
	return
    }

    #puts "RemoveHighlightPosition $prev_raster draw_rectangle $prev_pos -style $prev_style"

    eval $prev_raster draw_rectangle $prev_pos -style $prev_style
}

##############################################################################
proc HighlightPosition {raster Y} {
    global prev_style prev_pos prev_raster $raster.style

    #puts "HighlightPosition $raster"
    RemoveHighlightPosition

    if {$raster == ""} {
	#RemoveHighlightPosition $raster
	catch {unset prev_raster}
 	return
    }
    if {[string compare [winfo class $raster] Raster] != 0} {
	#RemoveHighlightPosition $raster
	catch {unset prev_raster}
	return
    }

    set r_win [winfo parent $raster]
    if {![info exists prev_pos]} {
	set prev_pos "0 0 0 0"
    }

    if {![info exists $raster.style]} {
	set $raster.style [$raster envcreate -function xor -fg #d9d9d9 -linewidth 3]
    }

    set pos [winfo rooty $raster]
    set win_ht [winfo height $raster]
    set ht3 [expr $win_ht / 3]

    set top [expr $pos + $ht3]
    set bottom [expr $top + $ht3]
	
    #puts "y $Y pos $pos top $top bottom $bottom win_ht $win_ht ht3 $ht3"
	
    set x0 0
    set x1 [winfo width $raster]
    if {$Y >= $pos && $Y < $top} {
	set y0 0
	set y1 $ht3
    } elseif {$Y >= $top && $Y < $bottom} {
	set y0 $ht3
	set y1 [expr $ht3 + $ht3]
    } else {
	set y0 [expr $ht3 + $ht3]
	set y1 $win_ht
    }

    set wx0 [lindex [$raster toworld $x0 $y0] 0]
    set wy0 [lindex [$raster toworld $x0 $y0] 1]
    set wx1 [lindex [$raster toworld $x1 $y1] 0]
    set wy1 [lindex [$raster toworld $x1 $y1] 1]

    set position "$wx0 $wy0 $wx1 $wy1"

    #if {$position != $prev_pos} {
	#eval $raster draw_rectangle [set $raster.prev_pos] -style [set $raster.style]
    #}

    eval $raster draw_rectangle $position -style [set $raster.style]

    set prev_style [set $raster.style]
    set prev_pos $position
    set prev_raster $raster
}

##############################################################################
#find where to position the result when let go of the mouse
proc FindResultPosition {r_win x y job} {
    global $r_win.win_list $r_win.first_row
    global r_from r_to

    set row 0
    set ht [winfo rooty $r_win]
    set min $ht
    set ht_list ""

    #puts "FindResultPosition"

    set raster_stem [GetRasterStem $r_win]
    set max_y 0
    
    if {![info exists $r_win.win_list]} {
	return -999
    }

    foreach win [set $r_win.win_list] {
	set id [lindex $win 0]
	set raster $raster_stem$id
	set win_y [winfo rooty $raster]
	if {$win_y > $max_y} {
	    set max_r $raster
	}
	set max_y $win_y
	lappend ht_list "[format %10d $win_y] [winfo height $raster]"
    }	

    #awful cheat to sort ht_list based on ascii - but need to format with
    #leading spaces
    set ht_list [lsort $ht_list]

    #add bottom window height
    lappend ht_list "[expr [format %10d [lindex [lindex $ht_list end] 0] + [winfo height $max_r]]] [winfo height $raster]"
    
    set cnt 1
    foreach l $ht_list {
        foreach {ht win_ht} $l {break}

	set ht3 [expr $win_ht / 3]
	set top [expr $ht + $ht3]
	set bottom [expr $ht + $ht3 + $ht3]
	
	#puts "y $y ht $ht top $top bottom $bottom end [expr $ht + $win_ht] win_ht $win_ht"
	if {$y > $ht && $y < $top} {
	    #puts "PLACE AT ABOVE"
	    break
	} elseif {$y >= $top && $y < $bottom} {
	    #puts "PLACE ON TOP"
	    return -1
	} elseif {$y >= $bottom && $y < [expr $ht + $win_ht]} {
	    #puts "PLACE AT BOTTOM $job"
	    #if moving a result, incr row cos later I delete a raster
	    #if moving a raster, don't incr row
	    if {$job == 0} {
		incr row
	    }

	    #8/1/98 seem to need this now, eg if do compare spans twice, move
	    #one plot on its own and then move to bottom using move raster

	    #9/3/00 don't seem to need this - can't get previous comment to
	    #go wrong. However, the bug this causes is if you move a plot
	    #using the all tag to the bottom of the bottom plot, an extra
	    #frame is inserted

	    #17/10/02 all depends on whether you are dropping a raster from
	    #within a window or from a different window as to whether you
	    #delete and also create a raster
	    if {[winfo parent $r_from] != [winfo parent $r_to]} {
		incr row
	    }
	    if {$job == 1} {
		#incr row
	    }
	    break
	} else {
	    #puts "next row"
	    incr cnt
	    incr row
	}
    }
    #puts "ROW [expr $row + [set $r_win.first_row]]"

    return [expr $row + [set $r_win.first_row]]
}


##############################################################################
#HACK to generalise
proc UpdateScaleBars {r_win old_size new_size} {

    puts UpdateScaleBars

    set from [$r_win.scale_x cget -from]
    set to [$r_win.scale_x cget -to]
    set from [$r_win.scale_y cget -from]
    set to [$r_win.scale_y cget -to]

    set o_x0 [lindex $old_size 0]
    set o_x1 [lindex $old_size 2]
    set o_range [expr $o_x1 - $o_x0]

    set n_x0 [lindex $new_size 0]
    set n_x1 [lindex $new_size 2]
    set n_range [expr $n_x1 - $n_x0]

    #set value [expr double($n_range) / $o_range]

    #set x [expr double($n_range) * 100 / $o_range]
    #set value [expr 100 - $x + 1]

    set value [raster_get_xmag -new_range $n_range -old_range $o_range]

    $r_win.scale_x set $value

    set o_y0 [lindex $old_size 1]
    set o_y1 [lindex $old_size 3]
    set o_range [expr $o_y1 - $o_y0]

    set n_y0 [lindex $new_size 1]
    set n_y1 [lindex $new_size 3]
    set n_range [expr $n_y1 - $n_y0]

    set value [expr double ($n_range / $o_range)]

    $r_win.scale_y set $value
}

##############################################################################
#superimpose a result on another raster
proc MoveRasterResult_Superimpose {seq_id r_type from to x y job} {
    
    #puts "MoveRasterResult_Superimpose $from $to"
    #Don't bother to superimpose if not moving the window
    if {$from == $to} {
	return
    }

    set from_p [winfo parent $from]

    #update new raster name in C structure and replot
    set old_scroll [$to world_size]
    
    #puts "!!!!!!!!!!!!!!!MoveRasterResult_Superimpose $from $to"
    set r_id [GetRasterId $from]

    global $from_p.win_list
    #puts "************************win list [set $from_p.win_list]"

    #want to remove old raster from list if there was only a single result in
    #the raster and the old raster is not the new raster before doing updates
    #which looks through all the windows in the win_list
    if {([GetRasterId $to] != $r_id) && 
	[raster_results -id $r_id -option number] == 1} {
	RemoveWin $from_p $r_id
    }

    #puts "************************win list [set $from_p.win_list]"

    update_raster_window -old $from -new $to -result_id $seq_id -new_id [GetRasterId $to] -old_id $r_id -job $job

    set new_scroll [$to world_size]

    #puts "OLD_SCROLL $old_scroll"
    #puts "NEW_SCROLL $new_scroll"
    UpdateKey $seq_id $r_id $from_p $to $r_type

    #find number of results are in raster
    #puts "num results [raster_results -id $r_id -option number]"
    #if {[raster_results -id $r_id -option number] == 0} {
	#seq_result_update -index [GetRasterId $from] -job QUIT
    #}
    #puts "END MoveRasterResult_Superimpose"
}

##############################################################################
#move result to a new raster
proc MoveRasterResult_New {result_id r_type from w_to x y row} {
    global tk_utils_defs

    #puts "********MoveRasterResult_New to $w_to FROM $from $r_type row $row"
    set id [GetRasterId $from]

    set from_p [winfo parent $from]

    set size [GetRasterResultSize $result_id]

    set raster_width [lindex $size 0]
    set raster_height [lindex $size 1]
    set ruler_width [winfo width $from_p.ruler_v$id]
    set end [string range $from [string length $from_p] end]

    #create new raster
    set raster $w_to[keylget tk_utils_defs RASTER.R.WIN]

    set seq_id_array [get_result_seq_id -result_id $result_id]

    #puts "result_id $result_id seq_id_array $seq_id_array"
    set id [AddRaster $w_to $raster $raster_height $raster_width \
		$ruler_width $row $seq_id_array $r_type]

    set raster $raster$id
    MoveRasterResult_Superimpose $result_id $r_type $from $raster $x $y ADD
    return $raster
}

##############################################################################
proc MoveRaster_Superimpose {from to x y} {

    #puts MoveRaster_Superimpose
    #don't bother to superimpose if not moving window
    if {$from == $to} {
	return
    }

    #update new raster name in C structure
    set old_scroll [$to world_size]

    set from_id [GetRasterId $from]
    set to_id [GetRasterId $to]
    set from_p [winfo parent $from]
    set to_p [winfo parent $to]

    RemoveWin $from_p $from_id
    #puts "MoveRaster_Superimpose $from $to"

    MoveKey $from_p $to_p $from_id $to_id $to

    update_raster_window -old $from -new $to -new_id [GetRasterId $to] -old_id [GetRasterId $from] -job SUPER

    set new_scroll [$to world_size]


    #destroy from raster
    #seq_result_update -index [GetRasterId $from] -job QUIT

    #need to update plot
    #raster_replot_all -raster $to

    #remove toplevel if no more slaves left
    global $from_p.row $from_p.first_row
    if {[set $from_p.row] == [set $from_p.first_row]} {
	destroy $from_p
    }
}

##############################################################################
#move raster within the same window
proc MoveRaster_Same {from to x y row r_type} {
    global SCALE_BAR ZOOM_SCALE
    #puts "^^^^^^^MoveRaster_Same $from $to $x $y row $row"
    

    #don't bother to superimpose if not moving window
    if {$from == $to} {
	return
    }

    set r_win [winfo parent $to]
    global $r_win.win_list

    set id [GetRasterId $from]
    set label $r_win.label$id
    set sb_v $r_win.sb_v$id
    set ruler_v $r_win.ruler_v$id
    set key $r_win.key$id
    set frame_from $r_win.b$id

    set id_to [GetRasterId $to]
    set frame_to $r_win.b$id_to

    #remember all grid info for all parts of the raster
    array set info_from [grid info $frame_from]
    array set info_to [grid info $frame_to]
    array set info_label [grid info $label]
    array set info_sb_v [grid info $sb_v]
    array set info_ruler_v [grid info $ruler_v]
    array set info_key [grid info $key]
    
    #remove all parts of raster
    grid forget $frame_from
    grid forget $label
    grid forget $sb_v
    grid forget $ruler_v
    grid forget $key

    #remove raster row from grid
    grid_delete $r_win row $info_from(-row) $info_from(-rowspan)
    #insert new raster row into grid
    grid_insert $r_win row $row $info_to(-rowspan)

    #set the new row value for each frame
    set info_from(-row) $row
    set info_label(-row) $row
    set info_sb_v(-row) $row
    set info_ruler_v(-row) $row
    set info_key(-row) $row

    #repack the frames
    grid rowconfig $r_win $row -weight 1 -minsize 40

    eval grid $frame_from [array get info_from]
    eval grid $label [array get info_label]
    eval grid $sb_v [array get info_sb_v]
    eval grid $ruler_v [array get info_ruler_v]
    eval grid $key [array get info_key]

     #foreach slave [grid slaves $r_win] {
 	#array set info [grid info $slave]
 	#puts "$slave column $info(-column) row $info(-row) "
     #}

    #need to regrid the x and y scalebars if inserted a new row at first position
    #if {$r_type == $SCALE_BAR || $r_type == $ZOOM_SCALE} {
        #MoveScaleBars $r_win $row
    #}
}

##############################################################################
#move raster to another window
proc MoveRaster_Different {r_type from to x y row} {
    global tk_utils_defs

    #puts "MoveRaster_Different row $row"

    #store info about from raster
    set from_p [winfo parent $from]
    set end [string range $from [string length $from_p] end]

    set to_p [winfo parent $to]
    set raster $to_p$end

    set from_id [GetRasterId $from]
    
    set raster_height [winfo height $from]
    set raster_width [winfo width $from]
    set ruler_width [winfo width $from_p.ruler_v$from_id]

    #create new raster
    set raster $to_p[keylget tk_utils_defs RASTER.R.WIN]

    set seq_id_array [get_result_seq_id -result_id [GetRasterId $from]]

    set id [AddRaster $to_p $raster $raster_height $raster_width \
		$ruler_width $row $seq_id_array $r_type]

    set raster $raster$id

    RemoveWin $from_p $from_id

    MoveKey $from_p $to_p [GetRasterId $from] $id $raster

    #update new raster name in C structure
    update_raster_window -old $from -new $raster -new_id [GetRasterId $raster] -old_id [GetRasterId $from] -job ADD

    #UpdateKey $seq_id $from $to

    #destroy from raster
    #seq_result_update -index [GetRasterId $from] -job QUIT

    #remove toplevel if no more slaves left
    #global $from_p.row $from_p.first_row
    #if {[set $from_p.row] == [set $from_p.first_row]} {
	#destroy $from_p
    #}
}

##############################################################################
#create a new window for a single result
proc InitResultWindow {from result_id r_type} {

    #puts "InitResultWindow $result_id"
    #check if the from raster contains 1 raster
    set from_p [winfo parent $from]

    global $from_p.row $from_p.first_row $from_p.title $from_p.borderwidth

    set from_id [GetRasterId $from]

    #if raster only contains one plot, no point moving result 
    if {([expr [set $from_p.row] - 1] == [set $from_p.first_row]) && 
	[raster_results -id $from_id -option number] == 1} {
	#puts "NO POINT"
	return
    }
    #store info about from raster
    set from_p [winfo parent $from]

    set size [GetRasterResultSize $result_id]
    set raster_width [lindex $size 0]
    set raster_height [lindex $size 1]

    set id [GetRasterId $from]
    #set title [GetRasterTitle]
    set title [set $from_p.title]

    #create new toplevel window
    set r_win [GetRasterWindow]

    set ruler_height [winfo height $from_p.ruler_h]
    set ruler_width [winfo width $from_p.ruler_v[GetRasterId $from]]

    CreateRasterWindow $r_type $r_win $title $ruler_height \
	$raster_width [set $from_p.borderwidth]

    global $r_win.row
    set row [set $r_win.row]

    #create new raster
    set raster_stem [GetRasterStem $r_win]

    set seq_id [get_result_seq_id -result_id $result_id]

    set r_id [AddRaster $r_win $raster_stem $raster_height $raster_width \
		$ruler_width $row $seq_id $r_type]

    set raster $raster_stem$r_id

    if {[raster_results -id [GetRasterId $from] -option number] == 0} {
	RemoveWin $from_p $from_id
    }
    #update new raster name in C structure and replot
    update_raster_window -old $from -new $raster -result_id $result_id \
	-new_id [GetRasterId $raster] -old_id [GetRasterId $from] -job NEW

    #puts "result $result_id from $from_id"
    UpdateKey $result_id $from_id $from_p $raster $r_type

    #need to destroy this last because I use the old raster in 
    #update_raster_window
    #if {[raster_results -id [GetRasterId $from] -option number] == 0} {
	#seq_result_update -index [GetRasterId $from] -job QUIT
    #}

}

##############################################################################
#create a new window for a whole raster
proc InitRasterWindow {from r_type} {

    #puts InitRasterWindow
    #check if the from raster contains 1 raster
    set from_p [winfo parent $from]

    global $from_p.row $from_p.first_row $from_p.title $from_p.borderwidth

    #if raster only contains one plot
    if {[expr [set $from_p.row] - 1] == [set $from_p.first_row]} {
	#puts "NO POINT"
	return
    }

    #store info about from raster
    set from_p [winfo parent $from]
    set raster_height [winfo height $from]
    set raster_width [winfo width $from]
    set ruler_height [winfo height $from_p.ruler_h]
    set ruler_width [winfo width $from_p.ruler_v[GetRasterId $from]]

    set from_id [GetRasterId $from]
    #set title [GetRasterTitle]
    set title [set $from_p.title]

    #create new toplevel window
    set r_win [GetRasterWindow]

    CreateRasterWindow $r_type $r_win $title $ruler_height $raster_width\
	[set $from_p.borderwidth]

    #create new raster
    global $r_win.row
    set row [set $r_win.row]

    #create new raster
    set raster_stem [GetRasterStem $r_win]

    set seq_id_array [get_result_seq_id -result_id [GetRasterId $from]]
    set r_id [AddRaster $r_win $raster_stem $raster_height $raster_width \
		$ruler_width $row $seq_id_array $r_type]

    set raster $raster_stem$r_id

    RemoveWin $from_p $from_id

    MoveKey $from_p $r_win $from_id [GetRasterId $raster] $raster

    #update new raster name in C structure and replot
    update_raster_window -old $from -new $raster -new_id [GetRasterId $raster] -old_id $from_id -job NEW

    #UpdateKey $seq_id $from $raster

    #need to destroy this last because I use the old raster in 
    #update_raster_window
    #seq_result_update -index [GetRasterId $from] -job QUIT
}

proc GetRasterWindow { } {
    global num_raster_win tk_utils_defs

    if {![info exists num_raster_win]} {
	set num_raster_win 0
    }

    #puts "*****************GetRasterWindow $num_raster_win"

    set r_win [keylget tk_utils_defs RASTER.WIN]$num_raster_win
    incr num_raster_win
    return $r_win
}

##############################################################################
proc SetRasterBinding {raster} {

    #bind $raster <B2-ButtonRelease> "+DropRaster %W %X %Y"
    #bind $raster <B2-Leave> "+LeaveRaster %W"
    bind $raster <Enter> "+EnterRaster %W"
}

proc ReSetRasterSize {raster width height} {

    $raster configure -width $width -height $height
    set raster_id [GetRasterId $raster]
    set r_win [winfo parent $raster]

    #delete contents of ruler if not needed
    set zoom [raster_results -id $raster_id -option zoom]
    if {$zoom <= 0} {
        $r_win.ruler_v$raster_id delete all
    }

    #need to configure ruler canvas aswell
    $r_win.ruler_v$raster_id configure -height $height
}

proc GetRasterId {raster} {
    regexp {[0-9]+$} $raster id
    return $id
}

proc GetRasterStem {r_win} {
    global tk_utils_defs

    #regsub {[0-9]+$} $raster "" raster_stem

    return $r_win[keylget tk_utils_defs RASTER.R.WIN]
}

#method of getting the parent of raster without the raster having to exist
proc GetRasterParent {raster} {

    regexp {^\.[a-z]+[0-9]+} $raster r_win
    return $r_win
}

proc RemoveRasterResultKey {raster text } {

    #puts "RemoveRasterResultKey $raster $text"

    set raster_id [GetRasterId $raster]
    set r_win [winfo parent $raster]

    keybox_delete $r_win.key$raster_id $text
}

proc rasterInitZoom {raster} {
    set r_win [winfo parent $raster]
    global $r_win.zoom_list 

    set world [$raster world]
    if {![info exists $r_win.zoom_list]} {
	set $r_win.zoom_list "{$world}"
    }
}

##############################################################################
#called from C when replot all
proc rasterRescaleZoom {w} {
    
    rasterInitZoom $w

}

#HACK - can't remember why I needed this and can't get things to go wrong
#without it
proc rasterRescaleZoomOLD { w } {

    set r_win [winfo parent $w]
    global $r_win.zoom_list $w.max_y

    puts rasterRescaleZoom
    set world [$w world]
    puts "world $world"
    set $w.max_y [lindex $world 3]
    puts "max_y [set $w.max_y]"
    if {![info exists $r_win.zoom_list]} {
	set $r_win.zoom_list "{$world}"
    }
    puts "zoomlist [set $r_win.zoom_list]"
    set length [expr [llength [set $r_win.zoom_list]] - 1]

    set orig_world [lindex [set $r_win.zoom_list] $length]
    puts "orig $orig_world"
    set dx [expr [lindex $world 2] - [lindex $orig_world 2]]
    set dy [expr [lindex $world 3] - [lindex $orig_world 3]]
    puts "dx $dx dy $dy"
 
    #necessary cos plotting upside down therefore must subtract new height
    #from zoom list!!!
    set l ""
    foreach i [set $r_win.zoom_list] {
	lappend l [list [lindex $i 0] \
		[expr [lindex $i 1]+$dy] \
		[expr [lindex $i 2]+$dx] \
		[expr [lindex $i 3]+$dy]]
    }
    #the last item in the zoom list ie the original zoom, is a special case
    #because dy only wants to be deleted from y2 not y1 (= 1)
    set i [lindex [set $r_win.zoom_list] end]
    set item [list [lindex $i 0] \
		[lindex $i 1] \
		[expr [lindex $i 2]+$dx] \
		[expr [lindex $i 3]+$dy]]

    puts "old zoomlist [set $r_win.zoom_list]"
    puts "item $item"

    set $r_win.zoom_list $l
    set $r_win.zoom_list [lreplace [set $r_win.zoom_list] end end $item]

    #puts "new zoomlist [set $r_win.zoom_list]"

    #HACK to dismiss first case when do want to rescale!
    if {[llength [set $r_win.zoom_list]] > 1} {
	eval $w world [lindex [set $r_win.zoom_list] 0]
    }
}


##############################################################################
#zooming commands
##############################################################################

##############################################################################
proc rasterZoom {r_win stem w s} {
    global $w.X1Y1 $w.X1Y1X2Y2 $w.Lastx $w.Lasty 
    global $r_win.zoom_list $r_win.win_list 

    #puts "RASTERZOOM $r_win $stem $w $s"
    if {![info exists $w.X1Y1]} {
	return
    }
    remove_raster_xh $w

    set newworld [set $w.X1Y1X2Y2]

    #avoid "coordinates must define a rectangle" error
    if {([lindex $newworld 0] == [lindex $newworld 2]) || \
	    ([lindex $newworld 1] == [lindex $newworld 3])} {
	raster_replot_zoom -raster $w
	return
    }

    set cur_world "[expr int([lindex $newworld 0])] [expr int([lindex $newworld 1])] [expr int([lindex $newworld 2])] [expr int([lindex $newworld 3])]"

    #set cur_world $newworld
    if {[lindex $newworld 0] > [lindex $newworld 2]} {
	set cur_world "[lindex $newworld 2] [lindex $newworld 1] [lindex $newworld 0] [lindex $newworld 3]"
    }
    if {[lindex $newworld 1] > [lindex $newworld 3]} {
	set cur_world "[lindex $cur_world 0] [lindex $newworld 3] [lindex $cur_world 2] [lindex $newworld 1]"
    }

    #puts "rasterZoom world $cur_world"

    set $r_win.zoom_list [linsert [set $r_win.zoom_list] 0 $cur_world]
    #puts "zoom [set $r_win.zoom_list]"
    set $w.X1Y1X2Y2 ""
    unset $w.X1Y1

    #zoom each raster in window
    foreach win [set $r_win.win_list] {
	set i [lindex $win 0]
	set raster $stem$i
	eval $stem$i world $cur_world
	raster_replot_zoom -raster $stem$i
    }
}

##############################################################################
proc rasterZoomUndo {r_win stem w} {
     global $r_win.zoom_list $w.Lastx $w.Lasty $r_win.win_list


    set $w.Lastx ""
    set $w.Lasty ""

    #puts "zoom_list [set $r_win.zoom_list]"
    if { [llength [set $r_win.zoom_list]] == 1} {
	return
    }
    #remove first element in zoom_list (ie last zoom)
    set $r_win.zoom_list [lrange [set $r_win.zoom_list] 1 [llength [set $r_win.zoom_list]]]
    #set world to the new first element in zoom_list (ie previous zoom)
    set newworld [lindex [set $r_win.zoom_list] 0]

    foreach win [set $r_win.win_list] {
	set i [lindex $win 0]
	#puts "undo $newworld"
	eval $stem$i world $newworld
	#puts "zoom [set $r_win.zoom_list]"
	raster_replot_zoom -raster $stem$i
    }
}

##############################################################################
proc rasterMark {w x y} {
    global $w.X1Y1 

    #before do zoom, ensure that don't have a "plot only" result 
    #ie num_matches is -1
    #set results [sip_num_result_type -type 2d]
    #if {$results < 1} {
	#bell
	#return
    #}
    set $w.X1Y1 [$w toworld $x $y]
}

##############################################################################
proc rasterStroke {w old new x y s} {
    global $w.X1Y1 $w.X1Y1X2Y2

    set cursor_env [$w envcreate -function xor -fg #d9d9d9]

    #$w.X1Y1 contains the top left corner of new box
    #$w.X1Y1X2Y2 contains x1y1 and x2y2 of previous box
    if {![info exists $w.X1Y1]} {
	return
    }

    if {$old  && [set $w.X1Y1X2Y2] != ""} {
	eval $w draw_rectangle [set $w.X1Y1X2Y2] -style $s
    }
    if {$new } {
	set X2Y2 [$w toworld $x $y]
	eval $w draw_rectangle [set $w.X1Y1] $X2Y2 -style $s
	set $w.X1Y1X2Y2 "[set $w.X1Y1] $X2Y2"
    }
}

##############################################################################
proc raster_zoominit {r_win stem w} {
    global $w.X1Y1X2Y2

    set zoom_env [$w envcreate -function xor -fg #d9d9d9]

    #bind $w <Button-1> "rasterGetPos %W %x %y"

    bind $w <Control-Button-3> "rasterMark %W %x %y"
    bind $w <Any-B3-Motion>  "rasterStroke %W 1 1 %x %y $zoom_env"
    bind $w <B3-ButtonRelease>  "rasterZoom $r_win $stem $w $zoom_env"
    bind $w <Control-B3-ButtonRelease> "rasterZoom $r_win $stem $w $zoom_env"

    #bind $w <B3-ButtonRelease>  "rasterZoom %W $zoom_env"
    #bind $w <Control-B3-ButtonRelease> "rasterZoom %W $zoom_env"
		  
    set $w.X1Y1X2Y2 ""
}

#called from C
#want to change the original zoom extents when add or remove sequences of
#different lengths
proc update_zoom_list {r_win job zoom} {
    global $r_win.zoom_list
    #puts "update_zoom_list1 [set $r_win.zoom_list]"

    if {$job == 0} {
	#add new extents to zoom list
	set $r_win.zoom_list [linsert [set $r_win.zoom_list] end $zoom]
    } else {
	#remove previous extents
	set $r_win.zoom_list [lreplace [set $r_win.zoom_list] end end]

	set cur_extents [lindex [set $r_win.zoom_list] end]
	#puts "cur_extents $cur_extents"

	#if current last element in zoom list is NOT equal to the current
	#extents (held in zoom) then must add current extents to end of 
	#zoom list
	if {([llength $cur_extents] == 0) || 
	    ([lindex $cur_extents 0] != [lindex $zoom 0] || 
	    [lindex $cur_extents 1] != [lindex $zoom 1] ||
	    [lindex $cur_extents 2] != [lindex $zoom 2] || 
	    [lindex $cur_extents 3] != [lindex $zoom 3])} {
	    set $r_win.zoom_list [linsert [set $r_win.zoom_list] end $zoom]
	}

    }

    #puts "update_zoom_list2 [set $r_win.zoom_list]"
}

#called from C: ReSetRasterWindowWorld
#if need to change world because scroll region has shrunk, also need to 
#change the associated zoom element in the zoom_list
proc change_zoom_list {r_win old_ele new_ele} {
    global $r_win.zoom_list

    #puts "old_ele $old_ele"
    #puts "new_ele $new_ele"

    if {$old_ele == $new_ele} {
	#puts SAME
	#return
    }
    set cnt 0
    foreach z [set $r_win.zoom_list] {
	if {[lindex $z 0] == [lindex $old_ele 0] &&
	    [lindex $z 1] == [lindex $old_ele 1] &&
	    [lindex $z 2] == [lindex $old_ele 2] &&
	    [lindex $z 3] == [lindex $old_ele 3]} {
	} 
	if {$z == $old_ele} {
	    break
	}
	incr cnt
    }
    #puts "change_zoom_list1 [set $r_win.zoom_list]"
    #set $r_win.zoom_list [lreplace [set $r_win.zoom_list] $cnt $cnt $new_ele]
    set l ""
    lappend l $new_ele
    lappend l [lindex [set $r_win.zoom_list] end]
    set $r_win.zoom_list $l
    #puts "change_zoom_list2 [set $r_win.zoom_list]"
}

##############################################################################
proc UpdateKey {result_id from_id from_win r_to r_type} {
    
    set r_win [winfo parent $r_to]
    #set from_id [GetRasterId $r_from]
    set to_id [GetRasterId $r_to]

    set text [seq_result_key_name -index $result_id]
    set colour [seq_result_info -index $result_id -option colour]

    #delete key in previous raster
    if {[winfo exists $from_win.key$from_id]} {
	keybox_delete $from_win.key$from_id $text
    }
    set r_win [winfo parent $r_to]
    #add key in new raster
    keybox_add $r_win.key$to_id -text $text \
	-background $colour \
	-enter "EnterKey $r_to $result_id" \
	-leave "LeaveKey $r_to" \
	-drop "DropResult $result_id $r_type" \
	-motion MotionRaster \
	-menu "seq_result_keybox_update $r_win $result_id \[seq_result_names -raster_id $to_id\]"
    fit_on_screen $r_win
}


##############################################################################
#move an entire key from 1 raster to another
proc MoveKey {from_p to_p f_id t_id raster} {

    set oldkey $from_p.key$f_id
    set newkey $to_p.key$t_id

    for {set i 1} {$i < [keybox_size $from_p.key$f_id]} {incr i} {

	set enter [keybox_entrycget $oldkey $i -enter]
	#puts "^^^^^^^^^^^^^^^^^^^^^ENTER $enter"
	scan $enter "%s %s %s" cmd r result_id
	#puts "cmd $cmd raster $r result_id $result_id"

	set r_win [winfo parent $raster]
	keybox_add $newkey -text [keybox_entrycget $oldkey $i -text] \
	    -background [keybox_entrycget $oldkey $i -background]\
	    -enter "$cmd $raster $result_id "\
	    -leave "LeaveKey $raster" \
	    -drop [keybox_entrycget $oldkey $i -drop] \
	    -motion MotionRaster \
	    -menu "seq_result_keybox_update $r_win $result_id \[seq_result_names -raster_id $t_id\]"
	fit_on_screen $r_win
    }

    for {set i [expr [keybox_size $from_p.key$f_id] - 1]} {$i >= 0} {incr i -1} {
	keybox_delete $oldkey $i
    }

}

proc ConfigKey {result_id colour} {

    set raster [seq_result_info -index $result_id -option raster]
    set r_win [winfo parent $raster]


    set text [seq_result_key_name -index $result_id]
    keybox_entryconfigure $r_win.key[GetRasterId $raster] $text -background $colour


}

#return the default raster window size for a particular result
proc GetRasterResultSize {result_id} {

    return [seq_result_info -index $result_id -option win_size]
}
    
##############################################################################
#create a raster window
proc CreateRasterWindow {r_type r_win title ruler_height raster_width borderwidth} {
    global colour col_index $r_win.row $r_win.first_row $r_win.borderwidth
    global $r_win.title tk_utils_defs NUM_COLOURS spin_defs
    global ZOOM_BOX ZOOM_SCALE

    #puts "*****************CreateRasterWindow $r_win"
    set MAX_ROW 999

    #create raster window
    xtoplevel $r_win
    fix_maxsize $r_win

    frame $r_win.left
    frame $r_win.right
    wm protocol $r_win WM_DELETE_WINDOW "RasterStartShutdown $r_win"
    wm title $r_win $title
    
    frame $r_win.r_main
    grid $r_win.r_main -row 2 -column 0 
    set $r_win.row 2
    set $r_win.first_row 2

    #set some global info about all rasters
    set $r_win.borderwidth $borderwidth
    set $r_win.title $title

    set raster_stem [GetRasterStem $r_win]
    #puts "$r_win $raster_stem"

    frame $r_win.menubar -relief raised -bd 2
    seq_rasterMenu $r_win $title $r_win.menubar
    frame $r_win.buttons
    canvas $r_win.ruler_h -height $ruler_height -highlightthickness 0\
	    -width $raster_width
    
    if {$r_type == $ZOOM_BOX || $r_type == $ZOOM_SCALE} {
	rasterButtonsDot $r_win $raster_stem $r_win.buttons 
    } else { 
	rasterButtonsGraph $r_win $raster_stem $r_win.buttons 
    }

    scrollbar $r_win.sb_h -orient horizontal \
	-command "scrollXCmd $r_win $raster_stem $r_win.ruler_h" 

    frame $r_win.scale

    scale $r_win[keylget tk_utils_defs RASTER.SCALEX.WIN]\
	-from [keylget tk_utils_defs RASTER.SCALEX.MIN] \
	-to [keylget tk_utils_defs RASTER.SCALEX.MAX] \
	-showvalue 0 -orient horizontal -length $ruler_height\
	-command "RasterScaleX $r_win $raster_stem $r_win.ruler_h"

    label $r_win.label_x -text X

    scale $r_win[keylget tk_utils_defs RASTER.SCALEY.WIN]\
	-from [keylget tk_utils_defs RASTER.SCALEY.MIN] \
	-to [keylget tk_utils_defs RASTER.SCALEY.MAX] \
	-showvalue 0 -orient horizontal -length $ruler_height\
	-command "RasterScaleY $r_win $raster_stem $r_win.ruler_v"

    label $r_win.label_y -text Y

    label $r_win.status_line

    grid columnconfig $r_win 4 -weight 1
    grid rowconfig $r_win 2 -weight 1

    grid $r_win.menubar -row 0 -column 0 -sticky ew -columnspan 7
    grid $r_win.buttons -row 1 -column 0 -sticky ew -columnspan 7

    #grid $r_win.scale_x -row 2 -column 0 -sticky ns
    #grid $r_win.scale_y -row 2 -column 1 -sticky ns

    grid $r_win.ruler_h -row [expr $MAX_ROW - 3] -column 4 -sticky ew
    grid $r_win.sb_h -row [expr $MAX_ROW - 2] -column 4 -sticky ew

    grid $r_win.scale -row [expr $MAX_ROW - 1] -column 4 -sticky ew

    grid columnconfig $r_win.scale 1 -weight 1
    grid columnconfig $r_win.scale 3 -weight 1

    grid $r_win.label_x -in $r_win.scale -row 0 -column 0
    grid $r_win.scale_x -in $r_win.scale -row 0 -column 1 -sticky ew
    grid $r_win.label_y -in $r_win.scale -row 0 -column 2
    grid $r_win.scale_y -in $r_win.scale -row 0 -column 3 -sticky ew

    grid $r_win.status_line -row $MAX_ROW -column 4 -sticky ew

    #necessary to ensure the window is packed before pack any rasters
    update idletasks

    bind $r_win <Enter> "EnterRWin $r_win"
    bind $r_win <B2-Leave> "LeaveRWin $r_win"
    bind $r_win <Alt-Button-1> "LeaveRWin $r_win"

    set col_index 0
    set colour(0) [keylget spin_defs RASTER.COLOUR.0]
    set colour(1) [keylget spin_defs RASTER.COLOUR.1]
    set colour(2) [keylget spin_defs RASTER.COLOUR.2]
    set colour(3) [keylget spin_defs RASTER.COLOUR.3]
    set colour(4) [keylget spin_defs RASTER.COLOUR.4]
    set colour(5) [keylget spin_defs RASTER.COLOUR.5]
    set colour(6) [keylget spin_defs RASTER.COLOUR.6]
    set colour(7) [keylget spin_defs RASTER.COLOUR.7]
    set colour(8) [keylget spin_defs RASTER.COLOUR.8]
    set colour(9) [keylget spin_defs RASTER.COLOUR.9]
    set NUM_COLOURS 10
}

proc seq_rasterMenu {r_win title menubar} {
    global $r_win.replot

    #set up File menu
    menubutton $menubar.file -text "File" -menu $menubar.file.opts
    set m [menu $menubar.file.opts]
    $m add command -label "Exit" \
	    -command "RasterStartShutdown $r_win"
    
    #set up a View menu
    menubutton $menubar.view -text "View" -menu $menubar.view.opts
    set m [menu $menubar.view.opts]
    $m add command -label "Results manager" \
	    -command "rasterResultsManager {show_help spin {SPIN-Result-Manager}}"

    #set up a Results menu
    menubutton $menubar.results -text "Results" -menu $menubar.results.opts 
    set m [menu $menubar.results.opts -tearoff 0]
    bind $menubar.results <<select-menu>> "+ seq_result_menu_update $r_win \[seq_result_names\] "

    #set up Help menu
    menubutton $menubar.help -text "Help" -menu $menubar.help.opts
    set m [menu $menubar.help.opts]
    $m add command -label "Spin Plot" \
	-command "show_help spin {SPIN-Spin-Plot}" -state normal

    #do the packing
    pack $r_win.menubar.file $r_win.menubar.view $r_win.menubar.results \
	-side left -fill x
    pack $r_win.menubar.help -side right
}
