#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

proc gap_zoom {io id scroll args } {

    eval zoom_canvas -io $io -id $id -direction $scroll $args
}

#-----------------------------------------------------------------------------
# Cursor commands
#-----------------------------------------------------------------------------
# Create a new visible cursor
proc canvas_cursor_create {io cnum canvwin cursorid regid} {
    global gap_defs

    set height [$canvwin canvasy [winfo height $canvwin]]
    set colour [keylget gap_defs \
	TEMPLATE.EDITOR.CROSSHAIR_COLOUR[expr $cursorid%5]]
    set line   [$canvwin create line 10 [$canvwin canvasy 0] 10 $height \
	-fill $colour \
	-tag "cursor_$cursorid cursor"]
    $canvwin bind $line <<move-drag>> \
	"canvas_cursor_drag $io $cnum $canvwin $cursorid $regid %x"
}

# Delete an existing cursor
proc canvas_cursor_delete {io canvwin cursorid} {
    $canvwin delete cursor_$cursorid
}

# Move an existing cursor to 'x'
proc canvas_cursor_move {io cnum canvwin cursorid regid x} {
    # Create the cursor if it doesn't currently exist
    if {$cnum != 0} {
	if {[$canvwin find withtag cursor_$cursorid] == ""} {
	    canvas_cursor_create $io $cnum $canvwin $cursorid $regid
	}
    }

    # And now move it
    set height [$canvwin canvasy [winfo height $canvwin]]

    set y1 [lindex [$canvwin cget -scrollregion] 1]
    set y2 [lindex [$canvwin cget -scrollregion] 3]
   
    #no scrolling allowed eg ruler
    if {$y2 == 0} {
	set y2 [winfo height $canvwin]
    }

    #puts "canvas_cursor_move $canvwin $x"
    $canvwin coords cursor_$cursorid $x $y1 $x $y2
}

# Drag an existing cursor to pixel coord 'px'
proc canvas_cursor_drag {io cnum canvwin cursorid regid px} {

    #puts "canvas_cursor_drag $canvwin $px [$canvwin canvasx $px] [$canvwin canvasx 0]"

    #need to avoid confusions with selecting readings
    itemUnMark $canvwin
    
    # Convert the 'px' to a base x.
    set x [canvas_to_world -io $io -id $regid -x [$canvwin canvasx $px] \
	-cnum $cnum]

    # Move the cursor
    #canvas_cursor_move $io 0 $canvwin $cursorid 0 $px
    canvas_cursor_move $io 0 $canvwin $cursorid 0 [$canvwin canvasx $px]

    # Generate an event
    if {$x < 1} {
	set x 1
    } elseif {$x > [expr [c_length $io $cnum]+1]} {
	set x [expr [c_length $io $cnum]+1]
    }
    keylset l id $cursorid seq -1 pos -1 abspos $x sent_by $regid job MOVE
    contig_notify -io $io -cnum $cnum -type CURSOR_NOTIFY -args "$l"
}

# Use for creating a new editor at pixel position 'px', or for moving
# an existing editor (and cursor) to 'px'.
proc canvas_cursor_editor {io cnum regid px canvwin} {

    # Get the new cursor position in base coordinates
    set x [canvas_to_world -io $io -id $regid -x [$canvwin canvasx $px] \
	-cnum $cnum]
    if {$x < 1} {
	set x 1
    } elseif {$x > [expr [c_length $io $cnum]+1]} {
	set x [expr [c_length $io $cnum]+1]
    }

    # Find a cursor id
    set found 0
    foreach item [$canvwin find withtag cursor] {
	regsub {.*cursor_([^ ]*|$f).*} [$canvwin gettags $item] {\1} c

	# We've now got a cursor id, check if the editor is using it
	# This is done somewhat crudely by checking if the cursor
	# is 'private'.
	set clist [query_cursor -io $io -cursorid $c]
	if {[keylget clist private] != 0 &&
	    [keylget clist contig] == $cnum} {
	    set found 1
	    break
	}
    }

    if {$found} {
	# If we've found an editor, then simply move it
	keylset l id $c seq -1 pos -1 abspos $x sent_by $regid job MOVE

        contig_notify -io $io -cnum $cnum -type CURSOR_NOTIFY -args "$l"
    } else {
	# Otherwise create a new editor
	edit_contig -io $io -contig =$cnum -pos $x
    }
}


##############################################################################
#scroll all active windows in the x direction
proc gc_scroll_x { io id command x args} {
    scroll_canvas -io $io -id $id -xscrollcommand "$command $x $args"
} 
#end gc_scroll_x

##############################################################################
#scroll all active windows in the y direction
proc gc_scroll_y { io id command y args} {
    scroll_canvas -io $io -id $id -yscrollcommand "$command $y $args"
} 
#end gc_scroll_y
