#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
##############################################################################
proc RasterCreateSeqed { raster_id raster rx } {

    #only want to create a new seqed if one doesn't already exist
    set result [raster_find_edcursor -raster $raster -id $raster_id -pos $rx]
    set cursor_id [lindex $result 0]
    set seq_id [lindex $result 1]

    #puts "RasterCreateSeqed $cursor_id $seq_id"

    if {$cursor_id < 0} {
	#no cursors 
	set world [$raster toworld $rx 0]
	SeqedDisplay [expr int([lindex $world 0])] $seq_id
    } else {
	raster_move_cursor -id $raster_id -raster $raster -pos $rx -cursor_id $cursor_id
    }
}

##############################################################################
proc BindingCommandsGraph {id raster} {

    bind $raster <<move>> "+RasterMoveCursorGraph $id %W %x"
    bind $raster <<move-create>> "RasterCreateSeqed $id %W %x"

    bind $raster <<select>> "+RasterMoveCursorGraph $id %W %x"
    bind $raster <<use>> "RasterCreateSeqed $id %W %x"

}

##############################################################################
proc RasterMoveCursorGraph {raster_id raster rx} {

    set cid [raster_select_cursor_graph -id $raster_id -raster $raster -x $rx]
    if {$cid == -1} {
	bind $raster <<move-drag>> {}
    } else {
	bind $raster <<move-drag>> "raster_move_cursor -id $raster_id -raster $raster -pos %x -cursor_id $cid"
    }
}

##############################################################################
proc CreateRasterGraph {r seq_id type frame title raster_height raster_width \
	ruler_height ruler_width} {
    global tk_utils_defs spin_defs raster_cmds raster_reg_cmd
    global SCALE_BAR
    upvar $r raster
    
    #find the appropriate raster window for the current active sequence
    #puts CreateRasterGraph
    set raster [get_raster_frame_graph -seq_id [lindex $seq_id 0] -type $type -frame $frame]
   # puts "raster=$raster"

    set id [GetRasterId $raster]

    #cannot used winfo parent here because r_win may not exist yet
    regexp {^\.[a-z]+[0-9]+} $raster r_win

    global $r_win.row

    #create new raster window if necessary
    if {![winfo exists $r_win]} {
	CreateRasterWindow $SCALE_BAR $r_win $title $ruler_height $raster_width [keylget spin_defs GRAPH.RASTER.BORDERWIDTH]
    } else {
	raise $r_win
	focus $r_win	
    }

    set row [set $r_win.row]

    #HACK - but I can't see an easier way of ensuring that the bindings are
    #set when I create rasters within raster.tcl
    set raster_reg_cmd "seq_raster_reg -window \$raster -seq_id \$seq_id"
    set raster_cmds "BindingCommandsGraph \$id \$raster"    

    if {![winfo exists $raster]} {
	frame $r_win.right.$id
	RasterFrame $r_win $raster $id $raster_height $raster_width $ruler_width $row $SCALE_BAR

	BindingCommandsGraph $id $raster
	seq_result_list_update [keylget tk_utils_defs RASTER.RESULTS.WIN]
    }
    return $id
}

##############################################################################
proc rasterButtonsGraph {r_win stem buttons} {
    global $r_win.cross 

    #crosshair checkbutton
    checkbutton $buttons.cross -text crosshairs -variable $r_win.cross\
	-command "init_crosshair $r_win $stem"
    
    #position label
    label $buttons.pos1 -bd 2 -relief sunken -width 6
    label $buttons.pos2 -bd 2 -relief sunken -width 6

    pack $buttons.cross $buttons.pos1 $buttons.pos2 \
	    -expand yes -side left

}
