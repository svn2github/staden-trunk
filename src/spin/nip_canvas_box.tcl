#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#-----------------------------------------------------------------------------
# Cursor commands
#-----------------------------------------------------------------------------

# Create a new visible cursor
proc nip_canvas_cursor_create {seq_id canvwin cursorid regid colour} {
    global nip_defs

    set height [$canvwin canvasy [winfo height $canvwin]]

    #HACK - why 10?
    set line   [$canvwin create line 10 [$canvwin canvasy 0] 10 $height \
	-fill $colour \
	-tag "cursor_$cursorid cursor"]
    $canvwin bind $line <<move-drag>> \
	"nip_canvas_cursor_drag $seq_id $canvwin $cursorid $regid $colour %x"
}

# Delete an existing cursor
proc nip_canvas_cursor_delete {canvwin cursorid} {
    $canvwin delete cursor_$cursorid
}

# Move an existing cursor to 'x'
proc nip_canvas_cursor_move {seq_id canvwin cursorid regid colour x} {
    # Create the cursor if it doesn't currently exist
    if {$seq_id != -1} {
	if {[$canvwin find withtag cursor_$cursorid] == ""} {
	    nip_canvas_cursor_create $seq_id $canvwin $cursorid $regid $colour
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

    $canvwin coords cursor_$cursorid $x $y1 $x $y2
}

# Drag an existing cursor to pixel coord 'px'
proc nip_canvas_cursor_drag {seq_id canvwin cursorid regid colour px} {
    # Convert the 'px' to a base x.
    set x [nip_canvas_to_world -id $regid -x [$canvwin canvasx $px]]

    # Move the cursor
    nip_canvas_cursor_move -1 $canvwin $cursorid 0 $colour $px

    # Generate an event
    if {$x < 1} {
	set x 1
    } elseif {$x > [expr [s_length -seq_id $seq_id]+1]} {
	set x [expr [s_length -seq_id $seq_id]+1]
    }
    #keylset l id $cursorid abspos $x sent_by $regid
    cursor_notify -seq_num [get_seq_num $seq_id] -id $cursorid -pos $x
}

# Use for creating a new editor at pixel position 'px', or for moving
# an existing editor (and cursor) to 'px'.
proc nip_canvas_cursor_editor {seq_id regid px canvwin} {
    # Get the new cursor position in base coordinates

    set x [nip_canvas_to_world -id $regid -x [$canvwin canvasx $px]]
    if {$x < 1} {
	set x 1
    } elseif {$x > [expr [s_length -seq_id $seq_id]+1]} {
	set x [expr [s_length -seq_id $seq_id]+1]
    }

    # Find a cursor id
    set found 0
    foreach item [$canvwin find withtag cursor] {
	regsub {.*cursor_([^ ]*|$f).*} [$canvwin gettags $item] {\1} c

	# We've now got a cursor id, check if the editor is using it
	# This is done somewhat crudely by checking if the cursor
	# is 'private'.
	set seq_num [get_seq_num $seq_id]

	set clist [query_cursor -cursorid $c -seq_num $seq_num]
	if {[keylget clist private] != 0 &&
	    [keylget clist id] == $seq_id} {
	    set found 1
	    break
	}
    }

    if {$found} {
	# If we've found an editor, then simply move it
	cursor_notify -seq_num $seq_num -id $c -pos $x
	
    } else {
	# Otherwise create a new editor
	SeqedDisplay $x $seq_id
    }
}
