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

proc RasterCreateSeqPairDisp {raster_id raster rx ry} {
    global HORIZONTAL VERTICAL

    #puts "RasterCreateSeqPairDisp $raster_id rx $rx ry $ry"

    #only want to create a new seqed if one doesn't already exist
    set result_h [raster_find_edcursor -raster $raster -id $raster_id -pos $rx -direction $HORIZONTAL]
    set result_v [raster_find_edcursor -raster $raster -id $raster_id -pos $rx -direction $VERTICAL]

    set cursor_id_h [lindex $result_h 0]
    set seq_id_h [lindex $result_h 1]
    set cursor_id_v [lindex $result_v 0]
    set seq_id_v [lindex $result_v 1]
    set b [$raster cget -bd]
    
    #set rx [expr $rx -$b]
    #set ry [expr $ry -$b]

    if {$cursor_id_h < 0 || $cursor_id_v < 0} {
	#no cursors 

	set result_id [seq_find_result -raster $raster -seq_id_h $seq_id_h -seq_id_v $seq_id_v]

	set wx [expr int([lindex [$raster toworld $rx 0] 0])]
	set wy [expr int([lindex [$raster toworld 0 $ry] 1])]

	#turn y upside down
	set wy [rasterY $raster $wy]

	SequencePairDisplay $wx $wy $seq_id_h $seq_id_v $cursor_id_h $cursor_id_v $result_id

    } else {

	raster_move_cursor -id $raster_id -raster $raster -pos $rx -cursor_id $cursor_id_h -direction $HORIZONTAL

	raster_move_cursor -id $raster_id -raster $raster -pos $ry -cursor_id $cursor_id_v -direction $VERTICAL
    }
}

##############################################################################
proc RasterMoveHCursor {raster_id raster x cid } {
    global HORIZONTAL

    #puts "RasterMoveHCursor x $x"
    set b [$raster cget -bd]
    set x [expr $x -$b]

    raster_move_cursor -id $raster_id -raster $raster -pos $x -cursor_id $cid -direction $HORIZONTAL
}

##############################################################################
proc RasterMoveVCursor {raster_id raster y cid } {
    global VERTICAL 

    #puts "RasterMoveVCursor y $y"
    set b [$raster cget -bd]
    set y [expr $y -$b]

    raster_move_cursor -id $raster_id -raster $raster -pos $y -cursor_id $cid -direction $VERTICAL
}

##############################################################################
proc RasterMoveCursorDot {raster_id raster rx ry} {
    global HORIZONTAL

    #puts RasterMoveCursorDot

    set result [raster_select_cursor_dot -id $raster_id -raster $raster -x $rx -y $ry]

    set cid_h [lindex $result 0]
    set cid_v [lindex $result 1]

    #puts "******************** $cid_h $cid_v"

    if {$cid_h == -1 && $cid_v == -1} {
	#move neither
	bind $raster <<move-drag>> {}
    } elseif {$cid_h != -1 && $cid_v != -1} {
	#move both together
	bind $raster <<move-drag>> "RasterMoveHCursor $raster_id $raster %x $cid_h; RasterMoveVCursor $raster_id $raster %y $cid_v"
    } elseif {$cid_v == -1} {
	#move horizontal only
	bind $raster <<move-drag>> "RasterMoveHCursor $raster_id $raster %x $cid_h"
    } elseif {$cid_h == -1} {
	#move vertical only
	bind $raster <<move-drag>> "RasterMoveVCursor $raster_id $raster %y $cid_v"
    }
}

##############################################################################
proc DotBindingCommands {id raster} {

    bind $raster <<move>> "+RasterMoveCursorDot $id %W %x %y"
    bind $raster <<move-create>> "RasterCreateSeqPairDisp $id %W %x %y"

    bind $raster <<select>> "+RasterMoveCursorDot $id %W %x %y"
    bind $raster <<use>> "RasterCreateSeqPairDisp $id %W %x %y"

}

##############################################################################
proc CreateRasterDot {r seq_id_h seq_id_v title raster_height raster_width \
			  ruler_height ruler_width} {
    global tk_utils_defs spin_defs raster_cmds raster_reg_cmd
    global SCALE_BAR ZOOM_SCALE
    upvar $r raster

    #puts CreateRasterDot 

    set raster [get_raster_frame_dot -seq_id_h $seq_id_h -seq_id_v $seq_id_v]
    set id [GetRasterId $raster]

    #set r_win [winfo parent $raster]
    #cannot used winfo parent here because r_win may not exist yet
    regexp {^\.[a-z]+[0-9]+} $raster r_win

    global $r_win.row
    #set raster_stem $r_win[keylget tk_utils_defs RASTER.R.WIN]

    #create new raster window if necessary
    if {![winfo exists $r_win]} {
	CreateRasterWindow $ZOOM_SCALE $r_win $title $ruler_height $raster_width [keylget spin_defs DOT.RASTER.BORDERWIDTH]
    }

    set row [set $r_win.row]

    #HACK - but I can't see an easier way of ensuring that the bindings are
    #set when I create rasters within raster.tcl
    set raster_reg_cmd "seq_raster_reg -window \$raster -seq_id \$seq_id"
    set raster_cmds "DotBindingCommands \$id \$raster"    

    if {![winfo exists $raster]} {
	frame $r_win.right.$id
	RasterFrame $r_win $raster $id $raster_height $raster_width $ruler_width $row $ZOOM_SCALE

	DotBindingCommands $id $raster
	seq_result_list_update [keylget tk_utils_defs RASTER.RESULTS.WIN]
	
	set raster_stem [GetRasterStem $r_win]
	raster_zoominit $r_win $raster_stem $raster
	$r_win.buttons.back configure -command "rasterZoomUndo $r_win $raster_stem $raster"
    }

    #raster_zoominit $r_win $raster_stem $raster
    #HACK - terrible - but I only know raster at this point
    #$r_win.buttons.back configure -command "raster_zoom_undo_win $r_win $raster_stem $raster"
    #$r_win.buttons.back configure -command "rasterZoomUndo $raster"

    #raise window above main window - necessary on windows
    raise $r_win
    return $id
}

##############################################################################
proc rasterButtonsDot {r_win stem buttons} {
    global $r_win.cross

    button $buttons.back -text "zoom out"

    #crosshair checkbutton
    checkbutton $buttons.cross -text crosshairs -variable $r_win.cross\
	    -command "init_crosshair $r_win $stem"
    
    #position label
    label $buttons.pos1 -bd 2 -relief sunken -width 6
    label $buttons.pos2 -bd 2 -relief sunken -width 6

    pack $buttons.back $buttons.cross $buttons.pos1 $buttons.pos2 \
	    -expand yes -side left
}
