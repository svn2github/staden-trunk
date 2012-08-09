#
# The contig editor consists of multiple components, following the standard
# model, view, controller pattern.
#
# 1) Controller: A top-level window, w, containing one or two editors plus
#    all the usual GUI buttons and menus.
#    Display settings are local to this and are held within the global tcl
#    array named after $w. In the code I typically upvar this to a
#    local array named opt().
#    Button/checkbutton GUI textvariables are stored within opt
#    starting with capital letters.
#        opt(PackSequences)
#        opt(HideAnno)
#
# 2) Model: The contig itself (with an associated io). The gapio ($io) is
#    typically a child io so we can edit (copy on write).
#    The contig has a global array named contigIO_$crec, usually
#    upvared to cio. This is the level at which we handle undo and
#    contig registration events.
#        cio(Undo) = <commands>
#        cio(Redo) = <commands>
#        cio(io)   = io=0x69ec60
#        cio(base) = io=0x68afa0
#        cio(crec) = 298459
#        cio(ref)  = 2
#
# 3) View: A names and editor tk widget (ednames and editor), attached to the
#    io/contig. This will be something like $w.seqs.sheet.
#    It doesn't have much Tcl data locally as the settings are mainly
#    held in the C widgets.
#        .e1.ed1.pane.seq.sheet(displayPos) = 658122
#        .e1.ed1.pane.seq.sheet(parent)     = .e1.ed1.pane
#        .e1.ed1.pane.seq.sheet(reg)        = 2
#        .e1.ed1.pane.seq.sheet(top)        = .e1

# The join editor is a stack of two contig editors. As such we have
# toplevel window with settings that govern both visible editors.
# $opt(all_editors) is a list of editors visible for window.

catch {package require tablelist}
catch {namespace import tablelist::tablelist}

#-----------------------------------------------------------------------------
# IO/contig specific components. (child IOs, undo history)
#-----------------------------------------------------------------------------

# Conceptually each edit had an equal and opposite edit, although it may
# require some compound edits to achieve these.
#
# When we make a change we store on the undo list the command(s) required
# to reverse that change. (We actually just store an opcode & operands as
# this is more space efficient - see io_undo_exec.) 
#
# Note that in order for undo to restore the editing cursor position we
# typically store a cursor move event plus the edit event. These are grouped
# as a tcl list so they can be undone with as a single event. Cursor
# positioning is not a data-change, rather a view change. Hence the editor
# that has the cursor moved is the editor the user clicked Undo in (while
# the edit itself may have been made in another editor window on the same
# contig).
#
# Redo, for now, has been disabled until we get a realiable undo
# implementation fully tested.

proc io_child {io crec} {
    upvar \#0 contigIO_$crec cio

    if {![info exists cio(io)]} {
	set cio(base) $io
	set child     [$io child]
	set cio(io)   $child
	set cio(crec) $crec
	set cio(ref)  0
    }

    incr cio(ref)

    return $cio(io)
}

proc io_detach {crec} {
    upvar \#0 contigIO_$crec cio

    if {![info exists cio]} return

    incr cio(ref) -1

    if {[set cio(ref)] == 0} {
	$cio(io) close
	unset cio
    }
}

proc io_undo_exec {w crec cmdu} {
    upvar \#0 contigIO_$crec cio
    global $w

    set io [$w io]

    foreach cmd $cmdu {
	foreach {code op1 op2 op3 op4 op5} $cmd break
	switch -- $code {
	    C_SET {
		$w set_cursor $op1 $op2 $op3 1
		# This may change the Cutoffs status, so update button
		set top [set ${w}(top)]
		global $top
		set ${top}(Cutoffs) [lindex [$w configure -display_cutoffs] 4]
	    }

	    B_REP {
		set seq [$io get_sequence $op1]
		$seq replace_base $op2 $op3 $op4
		$seq delete
	    }

	    B_INS {
		set seq [$io get_sequence $op1]
		$seq insert_base $op2 $op3 $op4
		$seq delete
	    }

	    B_DEL {
		set seq [$io get_sequence $op1]
		$seq delete_base $op2
		$seq delete
	    }

	    B_MOVE {
		set c [$io get_contig $cio(crec)]
		foreach {p f} [$c remove_sequence $op1] break;
		$c add_sequence $op1 $op2 $p $f
		$c delete

		eval $w set_cursor [$w get_cursor relative]
	    }

	    B_CUT {
		set seq [$io get_sequence $op1]
		$seq set_clips $op2 $op3
		$seq delete
	    }

	    C_CUT {
		set contig [$io get_contig $op1]
		foreach {rec l r} $op4 {
		    set s [$io get_seq $rec]
		    $s set_clips_no_invalidate $l $r
		    $s delete
		}
		$contig invalidate_consensus $op2 $op3
		$contig delete
	    }

	    C_INS {
		set contig [$io get_contig $op1]
		$contig insert_base $op2
		foreach seq $op3 {
		    foreach {rec pos base val cut start} $seq break;
		    set seq [$io get_sequence $rec]

		    if {$pos >= abs([$seq get_length])} {
			# Add to end as "$contig insert_base" hasn't
			$seq insert_base $pos $base $val

			# Whether or not the insert also moved depends on
			# orientation, so we may also need to fix this
			set s_pos [$seq get_position]
			if {$s_pos != $start} {
			    $contig move_seq $rec [expr {$start-$s_pos}]
			}
		    } else {
			set s_pos [$seq get_position]
			if {$start != $s_pos} {
			    # $contig insert_base moved the sequence instead
			    # of inserting to it, probably because we're
			    # adding to the first base. Manually insert instead
			    $seq insert_base $pos $base $val

			    # Whether or not the insert also moved depends on
			    # orientation, so we may also need to fix this
			    set s_pos [$seq get_position]
			    if {$s_pos != $start} {
				$contig move_seq $rec [expr {$start-$s_pos}]
			    }
			} else {
			    # Otherwise it did insert a *, but it may not have
			    # been the original base we removed and the quality
			    # is almost certainly wrong. Replace it with the
			    # correct value.
			    $seq replace_base $pos $base $val
			}
		    }
			
		    foreach {b q c} [$seq get_base $pos] break
		    if {$c != $cut} {
			# Base was deleted from cutoff, but now in
			# used portion. Adjust clips to compensate
			set orient [$seq get_orient]
			set len [$seq get_length]
			foreach {l r} [$seq get_clips] break;
			if {$orient != 0} {
			    set pos [expr {abs($len)-$pos-1}]
			}

			if {$pos == [expr {$l-1}]} {
			    incr l
			} elseif {$pos == [expr {$r-1}]} {
			    incr r -1
			}

			$seq set_clips $l $r
		    }
		    $seq delete
		}
		$contig delete
	    }
	    
	    C_DEL {
		set contig [$io get_contig $op1]
		$contig delete_base $op2
		$contig delete
	    }

	    C_MOVE {
		set contig [$io get_contig $op1]
		$contig set_visible_start $op2
		$contig delete

		# Adjust scrollbar & cursor
		eval $w set_cursor [$w get_cursor relative]
	    }

	    C_RENAME {
		set iob [$io base]
		set contig [$iob get_contig $op1]
		$contig set_name $op2
		$contig delete

		contig_notify -io $iob -type LENGTH \
		    -cnum [$w contig_rec] -args {}

		contig_notify \
		    -io $iob \
		    -type RENAME \
		    -cnum [$w contig_rec] \
		    -args [list name $op2]

		$iob flush
	    }

	    C_SHIFT {
		set contig [$io get_contig $op1]
		$contig shift_base $op2 $op3
		set seqs [$contig seqs_in_range [expr {$op2-1}] [expr {$op2+2}]]
		set anno [$contig anno_in_range [expr {$op2-1}] [expr {$op2+2}]]
		$contig delete

		if {$seqs != $op4} {
		    array set old {}
		    array set new {}
		    foreach _ $op4 {
			set old([lindex $_ 2]) [lindex $_ 0]
		    }
		    foreach _ $seqs {
			set new([lindex $_ 2]) [lindex $_ 0]
		    }

		    foreach seq [array names new] {
			if {![info exists old($seq)]} {
			    continue
			}
			if {$new($seq) == $old($seq)} {
			    continue
			}

			set contig [$io get_contig $op1]
			foreach {p f} [$contig remove_sequence $seq] break;
			$contig add_sequence $seq $old($seq) $p $f
			$contig delete
			set s [$io get_sequence $seq]
			$s move_annos [expr {$old($seq)-$new($seq)}]
			$s delete
		    }
		}

		# Shift back consensus tags which shouldn't have moved
		set contig [$io get_contig $op1]
		foreach anno $op5 {
		    $contig move_anno $anno -1
		}
		$contig delete
	    }

	    T_DEL {
		set tag [$io get_anno_ele $op1]
		$tag remove
	    }

	    T_NEW {
		array set d $op1

		# If we have a record already, we never actually deallocated it
		# so reuse that record and add it back. This fixes bugs where
		# higher up in the undo list we may still have that tag record
		# number referred to.
		if {[info exists d(rec)] && $d(rec) > 0} {
		    set c [$io get_contig $cio(crec)]
		    # The anno isn't actually in $c any more as we're adding
		    # it back, but "$c move_anno" doesn't care: if the remove
		    # fails it will assume it wasn't there and so defaults
		    # to being an add_anno method instead (this doesn't
		    # exist, hence the hacky approach for quickness).
		    $c move_anno $d(rec) $d(start) $d(end) $d(otype) $d(orec)
		    $c delete
		    set rec $d(rec)
		} else {
		    set rec [$io new_anno_ele $d(otype) $d(orec) $d(start) $d(end)]
		}
		set t [$io get_anno_ele $rec]
		$t set_comment $d(anno)
		$t set_type $d(type)
		$t delete
	    }

	    T_MOVE {
		set s [$io get_sequence $op1]
		$s move_annos $op2
		$s delete
	    }

	    T_MOD {
		array set d $op2
		set tag [$io get_anno_ele $op1]
		if {[$tag get_comment] != $d(anno)} {
		    $tag set_comment $d(anno)
		}
		if {[$tag get_type] != $d(type)} {
		    $tag set_type $d(type)
		}
		if {[lsearch -exact {+ - . ?} [$tag get_direction]] != $d(strand)} {
		    $tag set_direction [string index "+-.?" $d(strand)]
		}


		if {[$tag get_obj_type] != $d(otype) || 
		    [$tag get_obj_rec]  != $d(orec) || 
		    [$tag get_position] != "$d(start) $d(end)"} {
		    $tag move $d(otype) $d(orec) $d(start) $d(end)
		}
		$tag delete
	    }
	    
	    default {
		puts stderr "Unknown undo command: $cmd"
	    }
	}
    }
}

proc io_store_undo {crec cmdu cmdr} {
    upvar \#0 contigIO_$crec cio

    lappend cio(Undo) [list $cmdu $cmdr]
    set cio(Redo) ""

    contig_notify -io $cio(base) -cnum $crec -type GENERIC \
	-args [list TASK_GENERIC "" data {undo normal redo disable}]
    contig_notify -io $cio(base) -cnum $crec -type CHILD_EDIT -args ""
}

proc io_undo {ed crec} {
    upvar \#0 contigIO_$crec cio
    
    foreach {cmdu cmdr} [lindex [set cio(Undo)] end] break
    io_undo_exec $ed $crec $cmdu
    #eval $cmdu
    #lappend cio(Redo) [list $cmdu $cmdr]
    set cio(Undo) [lrange $cio(Undo) 0 end-1]

    if {[llength $cio(Undo)] == 0} {
	contig_notify -io $cio(base) -cnum $crec -type GENERIC \
	    -args [list TASK_GENERIC "" data {undo disable}]
    }
    contig_notify -io $cio(base) -cnum $crec -type GENERIC \
	-args [list TASK_GENERIC "" data {redo normal}]
    contig_notify -io $cio(base) -cnum $crec -type CHILD_EDIT -args ""
}

proc io_undo_state {crec} {
    upvar \#0 contigIO_$crec cio
    if {[info exists cio(Undo)] && [llength $cio(Undo)] != 0} {
	return normal
    }
    return disabled
}

proc io_redo {crec} {
    upvar \#0 contigIO_$crec cio
    
    foreach {cmdu cmdr} [lindex [set cio(Redo)] end] break
    eval $cmdr
    lappend cio(Undo) [list $cmdu $cmdr]
    set cio(Redo) [lrange $cio(Redo) 0 end-1]

    if {[llength $cio(Redo)] == 0} {
	contig_notify -io $cio(base) -cnum $crec -type GENERIC \
	    -args [list TASK_GENERIC "" data {redo disable}]
    }
    contig_notify -io $cio(base) -cnum $crec -type GENERIC \
	-args [list TASK_GENERIC "" data {undo normal}]
    contig_notify -io $cio(base) -cnum $crec -type CHILD_EDIT -args ""
}

proc io_redo_state {crec} {
    upvar \#0 contigIO_$crec cio
    if {[info exists cio(Redo)] && [llength $cio(Redo)] != 0} {
	return normal
    }
    return disabled
}

#-----------------------------------------------------------------------------
# Contig registration hookups. This data is obviously held per view rather
# than per contig or per io.
#-----------------------------------------------------------------------------
proc contig_register_callback {ed type id args} {
    global $ed
    set w [set ${ed}(top)]
    global $w

    #puts [info level [info level]]

    switch $type {
	QUERY_NAME {
	    return "Contig Editor"
	}

	CHILD_EDIT -
	LENGTH {
	    $ed clear_visibility_cache
	    $ed redraw
	    # A bit obscure, but it ensures edSetApos() is called in C,
	    # keeping cached absolute and relative positions in sync
	    # incase the edit was moving a sequence.
	    eval $ed set_cursor [$ed get_cursor relative] 0
	}

	RENAME {
	    set io [set ${w}(io)]

	    if {[llength [set ${w}(all_editors)]] == 1} {
		set c1 [$io get_contig [set ${w}(contig)]]
		wm title $w "Edit: [$c1 get_name]"
		$c1 delete
	    } else {
		set c1 [$io get_contig [set ${w}(contig)]]
		set c2 [$io get_contig [set ${w}(contig2)]]
		wm title $w "Join: [$c1 get_name] / [$c2 get_name]"
		$c1 delete
		$c2 delete
	    }
	}

	GENERIC {
	    if {$ed != [set ${w}(curr_editor)]} return
	    foreach {component arg1 arg2 arg3} [lindex $args 2] {
		switch $component {
		    "undo" {
			$w.toolbar.undo configure -state $arg1
		    }
		    "redo" {
			#$w.toolbar.redo configure -state $arg1
		    }
		    "set_cursor" {
			if {[$ed contig_rec] == $arg1} {
			    $ed set_cursor 17 $arg1 $arg2
			} else {
			    $ed set_cursor 18 $arg1 $arg2
			}

			raise [set ${w}(win)]
		    }
		}
	    }
	}
	
	CURSOR_NOTIFY {
	    foreach a $args {
		foreach {k v} $a break;
		set arg($k) $v
	    }

	    # Only move the cursor if it's not sent by oursleves and
	    # it's our primary cursor.
	    # EDIT: Removed " || $arg(id) == 0"
	    if {[set ${ed}(reg)] != $arg(sent_by) && \
		    ($arg(id) == [$ed cursor_id])} {
		if {[$ed contig_rec] == $arg(seq) || $arg(seq) == -1} {
		    $ed set_cursor 17 [$ed contig_rec] $arg(abspos)
		} else {
		    $ed set_cursor 18 $arg(seq) $arg(pos)
		}

		raise [set ${w}(win)]
	    }
	}

	RAISE {
	    
	}

	JOIN_TO {
	    foreach a $args {
		foreach {k v} $a break;
		set arg($k) $v
	    }

	    # NB: What happens when we have contig 1 and contig 2 open, plus
	    # a join editor on 1+2 which we then join. We end up with two
	    # $io objects both for the same contig.
	    #
	    # We can ignore this whole problem as we can only join
	    # after saving.
	    #io_detach $arg(contig_num)
	    set ${w}(io) [io_child [set ${w}(io_base)] $arg(contig)]
	    $ed io [set ${w}(io)]

	    # Update the editor cached contig record details
	    set ${w}(-contig) $arg(contig)
	    $ed change_contig $arg(contig)
	    upvar \#0 contigIO_$arg(contig) cio
	    set cio(Undo) ""
	    set cio(Redo) ""
	    set cio(io) [$ed io]
	    set cio(crec) $arg(contig)
	    incr cio(ref)
	    
	    # If cursor was on the old contig record, move it to the
	    # new contig instead.
	    foreach {type rec pos} [$ed get_cursor relative] break
	    if {$rec == $arg(contig_num)} {
		set rec $arg(contig)
		incr pos $arg(offset)
		$ed set_cursor $type $rec $pos 0
	    }

	    # Correct the window title
	    upvar \#0 $w opt

	    set c [$opt(io) get_contig $opt(-contig)]
	    if {[info exists opt(-contig2)]} {
		set c2 [$opt(io2) get_contig $opt(-contig2)]
		wm title $w "Join: [$c get_name] / [$c2 get_name]"
		$c2 delete
	    } else {
		wm title $w "Edit: [$c get_name]"
	    }
	    $c delete

	    # Finally adjust the display.
	    set pos [$ed xview]
	    incr pos $arg(offset)
	    $ed xview $pos
	}

	GET_LOCK {
	    foreach a $args {
		foreach {k v} $a break;
		set arg($k) $v
	    }

	    # Check for lock 2 => WRITE
	    if {[expr {$arg(lock)&2}] == 2} {
		# Attempt to exit => save dialogue box
		if {[info exists ${w}(-contig2)] || ![editor_exit $w 1]} {
		    return [expr {$arg(lock) & ~2}]  ;# Disallow lock
		}
		update idletask                      ;# Ensure it's flushed
	    }

	    return $arg(lock)                        ;# Permit lock
	}

	REGISTER -
	DEREGISTER {
	    # nothing to do
	}

	HIGHLIGHT_READ {
	    set ${w}(ListSize) [ListSize [set ${w}(OutputList)]]
	    $ed redraw
	}

	QUIT {
	    return [editor_exit $w]
	}
	
	default {
	    puts "Event '$type $id $args' not handled"
	}
    }

    return
}

#-----------------------------------------------------------------------------
# User dialogue for starting up the editors
#-----------------------------------------------------------------------------
proc EditContig2 {io t id} {
    if {[set reading [contig_id_gel $id]] == ""} return
    if {[set crec [contig_id_rec $id]] == ""} return

    destroy $t
    catch {edit_contig -io $io -contig $crec -reading $reading}
    SetContigGlobals $io $reading
}

proc EditContig {io} {
    set t .cedialog
    if {[xtoplevel $t -resizable 0] == ""} { return }
    wm title $t "Edit contig"

    contig_id $t.id -io $io -range 0 -command "EditContig2 $io $t $t.id"
    okcancelhelp $t.but \
	-ok_command "EditContig2 $io $t $t.id" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 {Editor}"

    pack $t.id -side top -fill both
    pack $t.but -side bottom -fill both
}

proc JoinContig2 {io t id1 id2} {
    if {[set read1 [contig_id_gel $id1]] == ""} return
    if {[set read2 [contig_id_gel $id2]] == ""} return
    if {[set crec1 [contig_id_rec $id1]] == ""} return
    if {[set crec2 [contig_id_rec $id2]] == ""} return

    destroy $t
    join_contig -io $io \
	-contig  $crec1 -reading  $read1 -pos  1 \
	-contig2 $crec2 -reading2 $read2 -pos2 1
    SetContigGlobals $io $read1
}

proc JoinContig {io} {
    set t .jedialog
    if {[xtoplevel $t -resizable 0] == ""} { return }
    wm title $t "Join contigs"

    contig_id $t.id1 -io $io -range 0 -default "" -trace 2
    contig_id $t.id2 -io $io -range 0 -default "" -trace 0
    okcancelhelp $t.but \
	-ok_command "JoinContig2 $io $t $t.id1 $t.id2" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 {Editor-Joining}"

    #initialise current frame
    SetCurFrame $t [list $t.id1 $t.id2]
    bind [entrybox_path $t.id1.ent] <<select>> "SetCurFrame $t {$t.id1 $t.id2}"
    bind [entrybox_path $t.id2.ent] <<select>> "SetCurFrame $t {$t.id2 $t.id1}"
    pack $t.id1 $t.id2 -side top -fill both
    pack $t.but -side bottom -fill both
}

#-----------------------------------------------------------------------------
# Externally usable functions.
#-----------------------------------------------------------------------------

# The top-level interface called from TCL
proc edit_contig {args} {
    eval contig_editor [next_editor] $args
}

# The top-level interface called from TCL
proc join_contig {args} {
    eval contig_editor [next_editor] $args
}

# Move or create an editor at contig rec, cursor rec and position.
proc create_or_move_editor {io contig cursor_rec cursor_pos} {
    # Look for an existing editor
    set found 0
    foreach result [result_names -io $io] {
	foreach {crec id name} $result break;
	if {$name == "Contig Editor" && $crec == $contig} {
	    set found $id
	    break
	}
    }

    # If found, send a cursor movement event
    if {$found} {
	#puts "1 contig_notify -io $io -cnum $contig -type GENERIC -args [list TASK_GENERIC "" data "set_cursor $cursor_rec $cursor_pos"]"
	contig_notify -io $io -cnum $contig -type GENERIC \
	    -args [list TASK_GENERIC "" data "set_cursor $cursor_rec $cursor_pos"]
    } else {
	# Not found, so launch a new editor
	edit_contig \
	    -io      $io \
	    -contig  $contig \
	    -reading "#$cursor_rec" \
	    -pos     $cursor_pos
    }
}

#-----------------------------------------------------------------------------
# Internally usable functions.
#-----------------------------------------------------------------------------

# Allocates a window pathname for the contig editor
set editor_num 0
proc next_editor {} {
    global editor_num
    incr editor_num
    return ".e$editor_num"
}

#
# Creates a contig editor mega-widget. The layout top to bottom is:
# Menu bar
# Tool/button bar
# Paned window (left/right split)
# Status line
#
# The paned window contains the actual editor itself with the left side being
# Name + scrollbar underneath and the right side being editor + 2 scrollbars
# below and right.
#
# Returns the pathname of the newly created object: '$path'
#
proc contig_editor {w args} {
    global gap5_defs
    upvar \#0 $w opt

    # Initialise the $path global array - the instance data for this widget
    foreach {arg val} $args {
	set opt($arg) $val
    }

    set opt(win) $w
    set opt(Disagreements)    [keylget gap5_defs CONTIG_EDITOR.DISAGREEMENTS]
    set opt(DisagreeMode)     [keylget gap5_defs CONTIG_EDITOR.DISAGREE_MODE]
    set opt(DisagreeCase)     [keylget gap5_defs CONTIG_EDITOR.DISAGREE_CASE]
    set opt(DisagreeQuality)  [keylget gap5_defs CONTIG_EDITOR.DISAGREE_QUAL]
    set opt(PackSequences)    [keylget gap5_defs CONTIG_EDITOR.PACK_SEQUENCES]
    set opt(StripeMode)       [keylget gap5_defs CONTIG_EDITOR.STRIPE_SIZE]
    set opt(Quality)          [keylget gap5_defs CONTIG_EDITOR.SHOW_QUALITY]
    set opt(MappingQuality)   [keylget gap5_defs CONTIG_EDITOR.SHOW_MAPPING_QUALITY]
    set opt(Cutoffs)          [keylget gap5_defs CONTIG_EDITOR.SHOW_CUTOFFS]
    set opt(PosType)	      [keylget gap5_defs CONTIG_EDITOR.POS_TYPE]
    set opt(HideAnno)         0
    set opt(Status)           "--- Status info here ---"
    set opt(GroupByPrimary)   [keylget gap5_defs CONTIG_EDITOR.GROUP_BY_PRIMARY]
    set opt(GroupBySecondary) [keylget gap5_defs CONTIG_EDITOR.GROUP_BY_SECONDARY]
    set opt(GroupPrimary)     [keylget gap5_defs CONTIG_EDITOR.GROUP_PRIMARY]
    set opt(GroupSecondary)   [keylget gap5_defs CONTIG_EDITOR.GROUP_SECONDARY]

    set opt(io_base) $opt(-io)
    set opt(io) [io_child $opt(-io) $opt(-contig)]

    set join [info exists opt(-contig2)]

    #set opt(contig) [contig_order_to_number -io $opt(-io) -order 0]
    set opt(contig) $opt(-contig)
    if {[info exists opt(-reading)]} {
	if {[regexp {^\#([0-9]+)$} $opt(-reading) _dummy rec]} {
	    set opt(-reading) $rec
	} else {
	    set opt(-reading) [$opt(io) seq_name2rec $opt(-reading)]
	    if {$opt(-reading) == -1} { set opt(-reading) 0 }
	}
    } else {
	set opt(-reading) 0
    }
    if {![info exists opt(-pos)]} {
	if {$opt(-reading) != 0 && $opt(-reading) != $opt(-contig)} {
	    set s [$opt(-io) get_sequence $opt(-reading)]
	    foreach {cl cr} [$s get_clips] break;
	    if {[$s get_orient]} {
		set cl [expr {abs([$s get_length])-$cr}]
	    }
	    $s delete
	    set opt(-pos) $cl
	} else {
	    set opt(-pos) 1
	}
    }
    if {$join} {
	set opt(contig2) $opt(-contig2)
	set opt(io2) [io_child $opt(-io) $opt(-contig2)]
	if {[info exists opt(-reading2)]} {
	    if {[regexp {^\#([0-9]+)$} $opt(-reading2) _dummy rec]} {
		set opt(-reading2) $rec
	    } else {
		set opt(-reading2) [$opt(io) seq_name2rec $opt(-reading2)]
		if {$opt(-reading2) == -1} { set opt(-reading2) 0 }
	    }
	} else {
	    set opt(-reading2) 0
	}
	if {![info exists opt(-pos2)]}     {
	    if {$opt(-reading2) != 0 && $opt(-reading2) != $opt(-contig2)} {
		set s [$opt(-io) get_sequence $opt(-reading2)]
		foreach {cl cr} [$s get_clips] break;
		if {[$s get_orient]} {
		    set cl [expr {abs([$s get_length])-$cr}]
		}
		$s delete
		set opt(-pos2) [expr {$cl-1}]
	    } else {
		set opt(-pos2) 1
	    }
	}
    }

    # Create the window layout
    if {![winfo exists $w]} {
	toplevel $w
	set c [$opt(io) get_contig $opt(contig)]
	#$c dump_ps /tmp/tree.ps
	if {[$opt(io) read_only]} {
	    set extra "   *** READ-ONLY ***"
	} else {
	    set extra ""
	}
	if {$join} {
	    set c2 [$opt(io2) get_contig $opt(contig2)]
	    wm title $w "Join: [$c get_name] / [$c2 get_name]$extra"
	} else {
	    wm title $w "Edit: [$c get_name]$extra"
	}
    }
    wm resizable $w 1 1

    # The toolbar 
    set tool [frame $w.toolbar -bd 0]
    checkbutton $tool.cutoffs \
	-variable ${w}(Cutoffs) \
	-text Cutoffs \
	-command "editor_cutoffs $w"
    checkbutton $tool.quality \
	-variable ${w}(Quality) \
	-text Quality \
	-command "editor_quality $w"
    
    xcombobox $tool.list_xc \
	-label "          List:" \
	-textvariable ${w}(OutputList) \
	-command editor_olist_switch \
 	-postcommand "editor_olist_populate $tool.list_xc" \
	-default readings \
	-width 15
    bind $tool.list_xc <Any-Leave> {
	editor_olist_switch %W [%W get]
	focus [winfo toplevel %W]
    }

    label $tool.list_cl -text "  Size:"
    menubutton $tool.list_c \
	-textvariable ${w}(ListSize) \
	-menu $tool.list_c.opts \
	-indicatoron 1
    menu $tool.list_c.opts -tearoff 0
    $tool.list_c.opts add command \
	-label Save \
	-command "SaveListDialog $tool.list_save \[set ${w}(OutputList)\]"
    $tool.list_c.opts add command \
	-label Load \
	-command "LoadListDialog $tool.list_load \[set ${w}(OutputList)\]"
    $tool.list_c.opts add command \
	-label Edit \
	-command "ListEdit \[set ${w}(OutputList)\]"
    $tool.list_c.opts add separator
    $tool.list_c.opts add command \
	-label Clear \
	-command "ListClear \[set ${w}(OutputList)\]"
    $tool.list_c.opts add command \
	-label Delete \
	-command "ListDelete \[set ${w}(OutputList)\]; \
                  editor_olist_switch $w readings"

    global $w
    set ${w}(OutputList) readings
    set ${w}(ListSize) [ListSize [set ${w}(OutputList)]]

    button $tool.undo    -text Undo -command "editor_undo $w" \
	-state [io_undo_state $opt(contig)]
#    button $tool.redo    -text Redo -command "editor_redo $w" \
	-state [io_redo_state $opt(contig)]
    button $tool.search  -text Search \
	-command "create_search_win $w.search \"editor_search $w\" 0"
    button $tool.save -text Save -command "editor_save $w"
    wm protocol $w WM_DELETE_WINDOW "editor_exit $w"
    pack $tool.undo $tool.search $tool.cutoffs $tool.quality \
	-side left
    pack $tool.save -side right

    # Highlights of the current editor so we know what window the button
    # applies to.
    bind $tool.undo <Any-Enter> "editor_hl \[set ${w}(curr_editor)\] red;
                                 editor_undo_info $w"
    bind $tool.undo <Any-Leave> "editor_hl \[set ${w}(curr_editor)\] \#d9d9d9;
                                 editor_undo_info $w 1"
#    bind $tool.redo <Any-Enter> "editor_hl \[set ${w}(curr_editor)\] red"
#    bind $tool.redo <Any-Leave> "editor_hl \[set ${w}(curr_editor)\] \#d9d9d9"

    if {$join} {
	set opt(Lock) 1; # See default in tkEditor.c link_to command
	checkbutton $tool.lock \
	    -variable ${w}(Lock) \
	    -text Lock \
	    -command "editor_lock $w"
	pack $tool.lock -side left

	button $tool.align \
	    -text Align -padx 2 \
	    -command "\[set ${w}(curr_editor)\] join_align"
	button $tool.alignL \
	    -text "<" -padx 2 \
	    -command "\[set ${w}(curr_editor)\] join_align 0 1"
	button $tool.alignR \
	    -text ">" -padx 2 \
	    -command "\[set ${w}(curr_editor)\] join_align 1 0"
	pack $tool.alignL $tool.align $tool.alignR \
	    -side left -fill both -padx 0

	button $tool.join -text Join -command "editor_join $w"
	pack $tool.join -side right

	if {[$c get_rec] == [$c2 get_rec]} {
	    $tool.join configure -state disabled
	}
    }

    pack $tool.list_xc $tool.list_cl $tool.list_c -side left

    if {[$opt(io) read_only]} {
	$tool.save configure -state disabled
	catch {$tool.join configure -state disabled}
    }

    # The editor(s) itself
    if {$join} {
	set pane0 $w.ed0
	set e [editor_pane $w $pane0 1 "" opt]
	set opt(editor2) $e
	lappend opt(all_editors) $e

	set pane1 $w.ed1
	set e [editor_pane $w $pane1 0 2 opt]
	set opt(editor1) $e
	lappend opt(all_editors) $e
    } else {
	set pane1 $w.ed1
	set e [editor_pane $w $pane1 0 "" opt]
	set opt(editor1) $e
	lappend opt(all_editors) $e
    }

    # Difference bar for the join editor
    if {$join} {
	set diffs [diff_pane $w]
	$e link_to $opt(editor2) $diffs.pane.seq.sheet
    }

    # The bottom status line
    set status $w.status
    frame $status -bd 2 -relief groove
    label $status.dummy
    label $status.l -textvariable ${w}(Status)
    pack  $status.dummy -fill both
    place $status.l -relx 0

    # Menu bar
    global contig_editor_main_menu
    $w configure -menu $w.menubar
    menu $w.menubar 
    create_menus $contig_editor_main_menu $w.menubar
    bind  $w.menubar <ButtonPress> "tag_repopulate_menu \[set ${w}(curr_editor)\]"

    bind $w.menubar <1> {+
	global editor_right_click
	set editor_right_click 0
    }

    # Packing
    grid rowconfigure $w 3 -weight 1
    grid columnconfigure $w 0 -weight 1
    grid $tool   -sticky nsew -row 0
    if {$join} {
	grid rowconfigure $w 1 -weight 1
	grid $pane0  -sticky nsew -row 1
	grid $diffs  -sticky nsew -row 2
    }
    grid $pane1  -sticky nsew -row 3
    grid $status -sticky nsew -row 4

    # Synchronised pane movement
    set opt(panes) $pane1.pane
    if {$join} {
	lappend opt(panes) $diffs.pane $pane0.pane
    }

    foreach p $opt(panes) {
	bind $p <ButtonRelease-1> "+sync_panes %W $w 1 1"
	bind $p <ButtonRelease-2> "+sync_panes %W $w 0 1"

	bind $p <Any-B2-Motion> "+sync_panes %W $w 0 0"
    }
    sync_panes $pane1.pane opt 0 1

    # Grid control
    # In theory this should work, but in practice getting the panedwindow
    # widget to correctly honour the side of the internal components is
    # a massive exercise in frustration, let alone making the pane "sash"
    # only move in increments of one font element.

    update idletasks
#    array set font [font metrics sheet_font]
#    wm grid $w [winfo width .e1] [winfo height .e1] \
#	[font measure sheet_font A] $font(-linespace)

    $e redraw
    $e show_cursor

    $c delete
    if {$join} {
	$c2 delete
    }

    after idle "focus \[set ${w}(curr_editor)\]"
}

proc editor_search {w args} {
    global $w
    set ed [set ${w}(curr_editor)]
    eval $ed search $args
}

proc editor_lock {w} {
    global $w
    set ed [set ${w}(curr_editor)]
    $ed lock [set ${w}(Lock)]
}

proc editor_save {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	if {[$ed save] != 0} {
	    bell
	}
    }
}

proc editor_join {w} {
    upvar \#0 $w opt

    # Gather location information about the join
    set ed  [lindex $opt(all_editors) 0]
    set ed2 [lindex $opt(all_editors) 1]

    if {[$ed contig_rec] == [$ed2 contig_rec]} {
	tk_messageBox \
	    -icon error \
	    -message "Cannot join a contig to itself." \
	    -title Error \
	    -type ok \
	    -parent $w
	return
    }

    if {[catch {foreach {len mis} [$ed join_mismatch] break}]} {
	set ret [tk_messageBox \
		     -icon error \
		     -message "Contigs do not overlap." \
		     -title "Error" \
		     -type ok \
		     -parent $w]
	return
    }

    set ret [tk_messageBox \
		 -icon question \
		 -message [format "Overlap length:\t\t%d\nPercentage mismatch:\t%5.2f%%\n\nMake join?" $len [expr {(100.0*$mis)/$len}]] \
		 -title "Join and quit Editor?" \
		 -type yesnocancel \
		 -parent $w]

    if {$ret == "cancel"} return

    if {$ret == "yes"} {
	foreach ed $opt(all_editors) {
	    if {[$ed save] != 0} {
		bell
		return
	    }
	}
	set ed [lindex $opt(all_editors) 0]

	# Force trace shutdown first as these get updated during join,
	# accessing now freed memory.
	global gap5_defs
	catch {destroy $opt(editor1)[keylget gap5_defs TRACE_DISPLAY.WIN]}
	catch {destroy $opt(editor2)[keylget gap5_defs TRACE_DISPLAY.WIN]}
	update idletasks

	$ed join
    }

    editor_exit $w
}

# Returns true if we really want to exit
#         false if we cannot exit (user hit cancel, or failed to save).
proc editor_exit {w {get_lock 0}} {
    global $w

    if {[winfo exists $w.save_dialog]} return

    set ret ""
    foreach ed [set ${w}(all_editors)] {
	if {![[set ${w}(io)] read_only] && [$ed edits_made]} {
	    if {$ret == ""} {
		# Ask once only and apply to both in join editor
		set ret [tk_messageBox \
			     -icon question \
			     -title "Save changes" \
			     -message "Edits have been made. Save changes?" \
			     -default yes \
			     -type yesnocancel \
			     -parent $w]
	    }
	    
	    if {$ret == "cancel"} {
		return 0
	    } elseif {$ret == "yes"} {
		if {[$ed save] != 0} {
		    bell
		    return 0
		}
	    }
	}
    }

    set detach ""
    foreach ed [set ${w}(all_editors)] {
	global $ed
	set id [set ${ed}(reg)]
	contig_deregister -io [set ${w}(io_base)] -id $id
	set id [set ${ed}(reg_all)]
	contig_deregister -io [set ${w}(io_base)] -id $id
	lappend detach [$ed contig_rec]
    }

    destroy $w

    foreach crec $detach {
	io_detach $crec
    }

    return 1
}

proc display_pos_set {ed pos} {
    global $ed
    upvar \#0 [set ${ed}(top)] opt

    if {$opt(PosType) == "P"} {
	set ${ed}(displayPos) $pos
    } else {
	set ${ed}(displayPos) [lindex [$ed reference_pos $pos] 0]
    }
}

proc jog_editor {cmd dist} {
    $cmd xview scroll $dist units
}

# Lays out the editor pane, with options of the scrollbars being above or
# below the editor (used in the join editor).
# Above indicates whether the scrollbars are above the text panels or below.
proc editor_pane {top w above ind arg_array} {
    upvar $arg_array opt
    global gap5_defs

    set hei [keylget gap5_defs CONTIG_EDITOR.MAX_HEIGHT]
    set fn_height [font metrics sheet_font -linespace]
    set lines [expr {int([winfo screenheight $top]/$fn_height)}]

    if {[info exists opt(Lock)]} {
	incr lines -9; # guess
	if {$hei > $lines/2} {
	    set hei [expr {int($lines/2)}]
	}
    } else {
	incr lines -6; # guess
	if {$hei > $lines} {
	    set hei $lines
	}
    }

    if {$above != 0} {
	set above 1
	set jogrow 0
	set scrollrow 1
	set textrow 2
	set cattop 0
    } else {
	set textrow   0
	set scrollrow 1
	set jogrow 2
	set cattop 1
    }

    set f $w
    frame $f -bd 3

    set w $f.pane
    panedwindow $w -orient horiz -bd 1 -relief sunken \
	-showhandle 1 -sashrelief raised

    frame $w.name -bd 0 -highlightthickness 0
    frame $w.seq -bd 0 -highlightthickness 0

    $w add $w.name $w.seq

    # Seqs panel
    set ed $w.seq.sheet
    if {$opt(Disagreements)} {
	set dis_mode $opt(DisagreeMode)
    } else {
	set dis_mode 0
    }
    set ed [editor $ed \
		-width  [keylget gap5_defs CONTIG_EDITOR.SEQ_WIDTH] \
		-height $hei \
		-xscrollcommand "display_pos_set $ed \[$ed xview\];
                                 $w.seq.x set" \
		-yscrollcommand "$w.seq.y set" \
		-qual_fg     [keylget gap5_defs CONTIG_EDITOR.QUAL_IGNORE] \
		-diff1_bg    [keylget gap5_defs CONTIG_EDITOR.DIFF1_BG] \
		-diff2_bg    [keylget gap5_defs CONTIG_EDITOR.DIFF2_BG] \
		-diff1_fg    [keylget gap5_defs CONTIG_EDITOR.DIFF1_FG] \
		-diff2_fg    [keylget gap5_defs CONTIG_EDITOR.DIFF2_FG] \
		-stripe_bg   [keylget gap5_defs CONTIG_EDITOR.STRIPE_BG] \
		-stripe_mode [keylget gap5_defs CONTIG_EDITOR.STRIPE_SIZE] \
		-qualcolour0 [keylget gap5_defs CONTIG_EDITOR.QUAL0_COLOUR] \
		-qualcolour1 [keylget gap5_defs CONTIG_EDITOR.QUAL1_COLOUR] \
		-qualcolour2 [keylget gap5_defs CONTIG_EDITOR.QUAL2_COLOUR] \
		-qualcolour3 [keylget gap5_defs CONTIG_EDITOR.QUAL3_COLOUR] \
		-qualcolour4 [keylget gap5_defs CONTIG_EDITOR.QUAL4_COLOUR] \
		-qualcolour5 [keylget gap5_defs CONTIG_EDITOR.QUAL5_COLOUR] \
		-qualcolour6 [keylget gap5_defs CONTIG_EDITOR.QUAL6_COLOUR] \
		-qualcolour7 [keylget gap5_defs CONTIG_EDITOR.QUAL7_COLOUR] \
		-qualcolour8 [keylget gap5_defs CONTIG_EDITOR.QUAL8_COLOUR] \
		-qualcolour9 [keylget gap5_defs CONTIG_EDITOR.QUAL9_COLOUR] \
		-bd 0 \
	        -consensus_at_top            $cattop \
	        -stack_mode                  $opt(PackSequences) \
		-display_differences         $dis_mode \
		-differences_case_sensitive  $opt(DisagreeCase) \
		-display_differences_quality $opt(DisagreeQuality) \
		-display_quality             $opt(Quality) \
		-display_mapping_quality     $opt(MappingQuality) \
		-display_cutoffs             $opt(Cutoffs) \
	        -hide_anno                   $opt(HideAnno) \
		-group_by_primary   	     $opt(GroupPrimary) \
		-group_by_secondary 	     $opt(GroupSecondary) \
		-pos_type		     [scan $opt(PosType) %c] \
		-fg black \
		-indelcolour [keylget gap5_defs CONTIG_EDITOR.INDEL_COLOUR] \
	        -bg [tk::Darken [. cget -bg] 115]]
    if {[info exists opt(curr_editor)]} {
	$opt(curr_editor) configure -hollow_cursor 1
    }
    set opt(curr_editor) $ed
    $ed configure -hollow_cursor 0

    # X and y scrollbars
    scrollbar $w.seq.x -orient horiz -repeatinterval 10
    scrollbar $w.seq.y -orient vert  -repeatinterval 10

    # Names panel
    set edname $w.name.sheet
    ednames $edname \
	-width 15 \
	-height $hei \
	-xscrollcommand "$w.name.x set" \
	-bd 0 \
	-fg black \
	-bg [tk::Darken [. cget -bg] 115]
    scrollbar $w.name.x -orient horiz

    entry $w.name.pos \
	-textvariable ${ed}(displayPos) \
	-font {Helvetica -12} \
	-width 15
    bind $w.name.pos <Return> "editor_goto $ed $w.seq.sheet"
    catch {bind $w.name.pos <KP_Enter> "editor_goto $ed $w.seq.sheet"}
    button $w.name.pos_type \
	-textvariable ${top}(PosType) \
	-padx 2 \
	-pady 0 \
	-highlightthickness 0 \
	-command "toggle_editor_pos_type $top"
#    xmenubutton $w.name.pos_type \
#	-textvariable $w.name.PosType \
#	-menu $w.name.pos_type.m
#    set m [menu $w.name.pos_type.m]
#    $m add radiobutton \
#	-label Padded \
#	-value P \
#	-variable $w.name.PosType \
#	-command "puts \[set $w.name.PosType\]"
#    $m add radiobutton \
#	-label Reference \
#	-value R \
#	-variable $w.name.PosType \
#	-command "puts \[set $w.name.PosType\]"
#    global $w.name.PosType
#    set $w.name.PosType P

    # The jog control for scrolling the editor
    set posh1 [font metrics [$w.name.pos cget -font] -linespace]
    incr posh1 [expr {2*([$w.name.pos cget -borderwidth]+
			[$w.name.pos cget -highlightthickness]+
			1)}]
    set posh2 [font metrics [$w.name.pos_type cget -font] -linespace]
    incr posh2 [expr {2*([$w.name.pos_type cget -borderwidth]+
			[$w.name.pos_type cget -highlightthickness]+
			[$w.name.pos_type cget -pady]+
			1)}]
    set posh [expr {$posh1>$posh2?$posh1:$posh2}]

    incr posh -8;# 4 lots of jog borderwidth
    jog $w.seq.jog \
	-orient horiz \
	-command "jog_editor $w.seq.sheet" \
	-repeatinterval 50 -width $posh \
	-borderwidth 2 \
	-highlightthickness 0

    # Pack it all in the paned window
    grid rowconfigure $w.name $textrow -weight 1
    grid columnconfigure $w.name 0 -weight 1
    grid $w.name.sheet - -row $textrow   -sticky nsew
    grid $w.name.x     - -row $scrollrow -sticky nsew
    grid $w.name.pos $w.name.pos_type -row $jogrow    -sticky nsew

    $w.name configure -width [expr {[keylget gap5_defs CONTIG_EDITOR.NAMES_WIDTH]*[font measure sheet_font A]}]

    grid rowconfigure $w.seq $textrow -weight 1
    grid columnconfigure $w.seq 0 -weight 1
    grid $w.seq.sheet $w.seq.y -row $textrow   -sticky nsew
    grid $w.seq.x              -row $scrollrow -sticky nsew
    grid $w.seq.jog            -row $jogrow -sticky nsew

    focus $w.seq.sheet

    # Initialise with an IO and link name/seq panel together
    global $ed $edname
    set ${ed}(parent) $w
    set ${ed}(top) $top 
    set ${ed}(reg) [contig_register \
			-io $opt(io_base) \
			-contig $opt(contig$ind) \
			-command "contig_register_callback $ed" \
			-flags [list ALL GENERIC CHILD_EDIT]]
    set ${ed}(reg_all) [contig_register \
			-io $opt(io_base) \
			-contig 0 \
			-command "contig_register_callback $ed" \
			-flags [list HIGHLIGHT_READ]]
    $ed init $opt(io$ind) $opt(contig$ind) $opt(-reading$ind) $opt(-pos$ind) $w.name.sheet

    if {$ind == 2} {
	set ${ed}(side) top
    } else {
	set ${ed}(side) bottom
    }
    set ${edname}(ed) $ed

    # Force new style mode
    $w.name.x set 0.0 0.1
    $w.seq.x set 0.0 0.1
    $w.seq.y set 0.0 0.1

    $w.seq.x configure -command "$w.seq.sheet xview"
    $w.seq.y configure -command "$w.seq.sheet yview"
    $w.name.x configure -command "$w.name.sheet xview"

    grid rowconfigure $f 0 -weight 1
    grid columnconfigure $f 1 -weight 1
    grid $f.pane           -row 0 -columnspan 2 -sticky nsew

    # Force scrollbar to be set to the correct size.
    $w.name.sheet xview moveto 0.0

    # Force the editing cursor to be visible
    eval $ed set_cursor [$ed get_cursor relative] 1

    return $ed
}

# Returns the current active editor for a join editor
proc curr_ed {e} {
    upvar \#0 $e ea
    upvar \#0 $ea(top) opt
    return $opt(curr_editor)
}

# The "differences" bar that separates a pair of join editors
proc diff_pane {w} {
    upvar \#0 $w opt
    set ed $opt(curr_editor)

    set w $w.diffs
    frame $w -bd 0 -padx 3
    set p [panedwindow $w.pane -orient horiz -bd 2 -relief sunken \
	      -showhandle 1 -sashrelief raised]
    pack $p -fill both -expand 1

    frame $p.name -bd 0 -highlightthickness 0
    label $p.name.diff -text "Differences" -anchor w
    button $p.name.dleft  -text "<" -pady 0 -padx 2 -highlightthickness 0 \
	-command "$ed prev_difference"
    button $p.name.dright -text ">" -pady 0 -padx 2 -highlightthickness 0 \
	-command "$ed next_difference"
    pack $p.name.dleft -side left
    pack $p.name.dright -side right
    pack $p.name.diff -expand 1 -side left

    frame $p.seq -bd 0 -highlightthickness 0
    sheet $p.seq.sheet
    pack $p.seq.sheet -fill both -expand 1

    $p add $p.name $p.seq
    
    return $w
}

# Synchronises the pane position between a list of panes
proc sync_panes {w arg_array proxy round} {
    upvar $arg_array opt

    if {[winfo class $w] != "Panedwindow"} {
	set w $w.pane
    }

    if {$proxy} {
	foreach {px y} [$w proxy coord] {}
    } else {
	foreach {px y} [$w sash coord 0] {}
    }

    if {$round} {
	set fs [font measure sheet_font A]
	set x [expr {int($px / $fs)*$fs+[$w cget -bd]+[$w cget -sashpad]+3}]
    } else {
	set x $px
    }

    foreach pane $opt(panes) {
	if {[winfo class $pane] != "Panedwindow"} {
	    set pane $pane.pane
	}

	if {!$round && $pane == $w} continue

	after idle "$pane sash place 0 $x $y"
    }
}

# Highlights the curr_editor window
proc editor_hl {ed col} {
    set w [winfo parent [winfo parent [winfo parent $ed]]]
    $w configure -bg $col
}

# A coordinate jump
proc editor_goto {ed w} {
    upvar \#0 $ed eopt
    upvar \#0 $eopt(top) opt

    set unpadded 0
    set pos $eopt(displayPos)

    if {[regexp {^([-+]?[0-9]+)[rR]$} $pos _ _p]} {
	set opt(PosType) R
	set pos $_p
    } elseif {[regexp {^([-+]?[0-9]+)[pP]$} $pos _ _p]} {
	set opt(PosType) P
	set pos $_p
    } elseif {[regexp {^([-+]?[0-9]+)[uU]$} $pos _ _p]} {
	set unpadded 1
	set pos $_p
    } elseif {![regexp {^[-+]?[0-9]+$} $pos]} {
	set opt(Status) "**ERROR**: $pos is not a number"
	bell
	return
    }

    set opt(Status) ""

    if {$unpadded} {
	set ppos [consensus_padded_pos \
		      -io     [$ed io]  \
		      -contig [$ed contig_rec] \
		      -pos    $pos]
	$w xview $ppos
    } elseif {$opt(PosType) == "P"} {
	$w xview $pos
    } else {
	set crec [$ed contig_rec]
	set c [[$ed io] get_contig $crec]
	catch {set ppos [$c ref_to_padded $pos [$ed xview]]}
	$c delete
	
	if {[info exists ppos] && $ppos != ""} {
	    $w xview $ppos
	} else {
	    set opt(Status) "**WARNING**: Unable to find reference position $pos"
	    bell
	}
    }
}

# Callback for cutoffs button
proc editor_cutoffs {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	$ed configure -display_cutoffs $opt(Cutoffs)
	# Move cursor to the consensus if disabling the cutoffs has now
	# hidden it.
	if {$opt(Cutoffs) == 0} {
	    set curr [$ed get_cursor relative]
	    eval $ed set_cursor [$ed get_cursor absolute] 0
	    eval $ed set_cursor $curr 0
	}
	$ed redraw
    }
}

# Callback for quality button
proc editor_quality {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	$ed configure -display_quality $opt(Quality)
	$ed redraw
    }
}

proc editor_disagreements {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	if {$opt(Disagreements)} {
	    $ed configure -display_differences $opt(DisagreeMode)
	} else {
	    $ed configure -display_differences 0
	}
	$ed configure -differences_case_sensitive $opt(DisagreeCase)
	$ed redraw
    }
}

proc editor_toggle_annos {w} {
    global $w

    if {[info exists ${w}(all_editors)]} {
	# Called either via the menu, after setting ther new value
	upvar \#0 $w opt
    } else {
	# Or via control-Q on an editor, to toggle the value
	upvar \#0 [set ${w}(top)] opt
	set opt(HideAnno) [expr {1-$opt(HideAnno)}]
    }

    foreach ed $opt(all_editors) {
	$ed configure -hide_anno $opt(HideAnno)
	$ed redraw
    }
}

proc ed2name {w} {
    return [regsub {\.seq\.sheet} $w .name.sheet]
}

proc name2ed {w} {
    return [regsub {\.name\.sheet} $w .seq.sheet]
}

proc set_editor_pack_sequences {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	$ed configure -stack_mode $opt(PackSequences)
	$ed xview scroll 1 units
	$ed xview scroll -1 units

	# Force redraw of name scrollbar
	[ed2name $ed] xview scroll 0 units
    }
}

proc set_editor_stripe_mode {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	$ed configure -stripe_mode [expr {2*$opt(StripeMode)}]
	$ed redraw
    }
}

proc set_editor_mapping_quality {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	$ed configure -display_mapping_quality $opt(MappingQuality)
	$ed redraw
    }
}

proc set_differences_quality_callback {w val} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	$ed configure -display_differences_quality $val
	$ed redraw
    }
}

proc set_differences_quality {w} {
    upvar \#0 $w opt
    global gap5_defs

    set ed $opt(curr_editor)

    set t $ed.qual_diff
    if {[xtoplevel $t -resizable 0] == ""} {return}
    wm title $t "Set differences quality"

    set start $opt(DisagreeQuality)

    scalebox $t.qual \
	-title "Quality" \
	-orient horizontal \
	-from 0 \
	-to 100 \
	-width 5 \
	-variable ${w}(DisagreeQuality) \
	-type CheckInt \
	-command "set_differences_quality_callback $w"
    $t.qual.scale configure -length 150

    okcancelhelp $t.ok \
	-ok_command "destroy $t" \
	-cancel_command "set_differences_quality_callback $w $start;
                         set ${w}(DisagreeQuality) $start;
                         destroy $t" \
        -help_command "show_help gap5 {Editor-Disagree}"

    pack $t.qual $t.ok -side top -fill both
}

proc toggle_editor_pos_type {w} {
    upvar \#0 $w opt

    switch -- $opt(PosType) {
	P {
	    set opt(PosType) R
	}
	R {
	    set opt(PosType) P
	}
    }

    set_editor_pos_type $w
}

proc set_editor_pos_type {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	$ed configure -pos_type [scan $opt(PosType) %c]
	$ed redraw
    }
}


proc editor_group_by {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	if {$opt(GroupByPrimary)} {
	    $ed configure -group_by_primary $opt(GroupPrimary)
	} else {
	    $ed configure -group_by_primary 0
	}
	
	if {$opt(GroupBySecondary)} {
	    $ed configure -group_by_secondary $opt(GroupSecondary)
	} else {
	    $ed configure -group_by_secondary 0
	}
	
	if {$opt(GroupByPrimary) && $opt(GroupPrimary) == 7} {
	    editor_sort_by_sequence $ed
	} else {
	    $ed set_sort_order 
	    $ed redraw
	}
    }
}



# sort_by_base depends on a base in the consensus
# (type 17) being set.  Code in set_sort_order forces
# a redraw.
proc order_update_on_consensus {w} {
    global $w
    
    foreach {type rec pos} [$w get_number] break
    
    set tl [winfo toplevel $w]
    
    if {$type == 17} {
    	upvar \#0 $tl opt
	
	foreach ed $opt(all_editors) {
	    $ed set_base_sort_point
	    $ed set_sort_order 
    	    $ed redraw
	}
    }
}


# if there is a new selection and group by sequence
# is chosen then sort.
proc editor_new_select_sort {w} {
    global $w
    
    set tl [winfo toplevel $w]
    upvar \#0 $tl opt

    if {$opt(GroupByPrimary) && $opt(GroupPrimary) == 7} {
    	editor_sort_by_sequence $w
    }
}
    


proc editor_sort_by_sequence {w} {
    set w2 [selection own]
    
    if {$w2 == ""} {
    	set w2 $w
    } else {
    	if {[winfo class $w2] != "Editor"} {
	    bell
	    return
	}
    }
    
    foreach {otype orec start end} [$w2 select get] break
    
    if {$start == $end} {
    	return
    }
    
    set tl [winfo toplevel $w2]
    
    upvar \#0 $tl opt
    
    foreach ed $opt(all_editors) {
    	set opt(GroupByPrimary) 1
	set opt(GroupPrimary) 7
	$ed configure -group_by_primary $opt(GroupPrimary)
	
    	$ed set_sequence_sort
	$ed set_sort_order
	
	# Force the editing cursor to be visible
    	eval $w set_cursor [$w get_cursor relative] 1
	$ed redraw
    }
}

#-----------------------------------------------------------------------------
# Output list handling
proc editor_olist_populate {w} {
    global NGLists NGSpecial

    # Filter some.
    set valid_lists ""
    foreach l $NGLists {
	if {[lsearch -exact $NGSpecial $l] == -1 || $l == "readings"} {
	    lappend valid_lists $l
	}
    }
    $w configure -values $valid_lists
}

proc editor_olist_switch {w l} {
    set w [winfo toplevel $w]
    upvar #0 $w opt

    global NGList
    if {![ListExists2 $l]} {
	ListCreate2 $l {} SEQID
    }
    InitListTrace $l
    set opt(OutputList) $l
    set opt(ListSize) [ListSize $opt(OutputList)]

    foreach ed $opt(all_editors) {
	$ed configure -output_list $l 
	$ed redraw
   }
}

#-----------------------------------------------------------------------------
# Undo support
proc store_undo {w cmdu cmdr} {
    io_store_undo [$w contig_rec] $cmdu $cmdr
}

proc editor_undo {top} {
    upvar \#0 $top opt
    
    set w $opt(curr_editor)
    upvar \#0 contigIO_[$w contig_rec] cio
    if {[llength $cio(Undo)] == 0} {
	bell
	return
    }
    io_undo $w [$w contig_rec]
    editor_undo_info $top
}

# Updates the information line with the top-most item on the undo stack
proc editor_undo_info {top {clear 0}} {
    upvar \#0 $top opt
    
    set w $opt(curr_editor)
    set crec [$w contig_rec]

    upvar \#0 contigIO_$crec cio

    if {$clear || ![info exists cio(Undo)] || $cio(Undo) == ""} {
	set opt(Status) ""
	return
    }
    foreach {cmdu cmdr} [lindex [set cio(Undo)] end] break

    set io [$w io]

    set c [$io get_contig $crec]

    set msg ""
    foreach cmd $cmdu {
	foreach {code op1 op2 op3 op4} $cmd break
	switch -- $code {
	    C_SET { }

	    B_REP {
		set s [$io get_sequence $op1]
		lappend msg "Change base in seq [$s get_name] at $op2 to base $op3, qual $op4"
		$s delete
	    }

	    B_INS {
		set s [$io get_sequence $op1]
		lappend msg "Insert base in seq [$s get_name] at $op2 with base $op3, qual $op4"
		$s delete
	    }

	    B_DEL {
		set s [$io get_sequence $op1]
		lappend msg "Delete base in seq [$s get_name] at $op2"
		$s delete
	    }

	    B_MOVE {
		set s [$io get_sequence $op1]
		lappend msg "Move seq [$s get_name] to position $op2"
		$s delete
	    }

	    B_CUT {
		set s [$io get_sequence $op1]
		lappend msg "Set clip points for seq [$s get_name] to left $op2, right $op3"
		$s delete
	    }

	    C_INS {
		foreach seq $op3 {
		    foreach {rec pos base val cut} $seq break;
		    append b $base
		}
		lappend msg "Insert column into contig at position $op2, bases $b"
	    }
	    
	    C_DEL {
		lappend msg "Delete column from contig at position $op2"
	    }

	    C_MOVE {
		lappend msg "Set contig start to $op2"
	    }

	    C_RENAME {
		lappend msg "Rename contig to $op2"
	    }

	    C_SHIFT {
		lappend msg "Shift column from contig at position $op2, dir $op3"
	    }
	    
	    T_DEL {
		set tag [$io get_anno_ele $op1]
		set otype [$tag get_obj_type]
		if {$otype == 18} {
		    set s [$io get_sequence [$tag get_obj_rec]]
		    set obj [$s get_name]
		    $s delete
		} else {
		    set obj "<consensus>"
		}
		foreach {start end contig} [$tag get_position] break;
		lappend msg "Remove annotation \#$op1: type=[$tag get_type], text=\"[$tag get_comment]\", object=$obj, position=$start..$end"
		$tag delete
	    }

	    T_NEW {
		array set d $op1
		if {$d(otype) == 18} {
		    set s [$io get_sequence $d(orec)]
		    set obj [$s get_name]
		    $s delete
		} else {
		    set obj "<consensus>"
		}
		lappend msg "Create annotation: type=$d(type), text=\"$d(anno)\", object=$obj, position=$d(start)..$d(end)"
	    }

	    T_MOVE {
		lappend msg "Move annotations on seq \#$op1 to $op2"
	    }

	    T_MOD {
		array set d $op2
		if {$d(otype) == 18} {
		    set s [$io get_sequence $d(orec)]
		    set obj [$s get_name]
		    $s delete
		} else {
		    set obj "<consensus>"
		}
		lappend msg "Modify annotation \#$op1: type=$d(type), text=\"$d(anno)\", object=$obj, position=$d(start)..$d(end)"
	    }
	    
	    default {
		lappend msg "Unknown undo command: $cmd"
	    }
	}
    }

    $c delete

    regsub -all "\n" $msg "\\n" msg
    set opt(Status) [string range [join $msg " / "] 0 1000]
}

# proc editor_redo {top} {
#     upvar \#0 $top opt
#     set w $opt(curr_editor)
#     io_redo [$w contig_rec]
# } 

#-----------------------------------------------------------------------------
# Basic sequence editing function

# Returns original read pos as it is needed for undo
proc editor_shift_seq {io rec dir {move_anno 1}} {
    # Dir 1 => right, -1 => left
    set seq [$io get_sequence $rec]
    set rpos [$seq get_position]
    set cnum [$seq get_contig]
    $seq delete
    set c [$io get_contig $cnum]
    foreach {pair_rec flags} [$c remove_sequence $rec] break;
    $c add_sequence $rec [expr {$rpos+$dir}] $pair_rec $flags
    $c delete
    if {$move_anno} {
	set seq [$io get_sequence $rec]
	$seq move_annos -1
	$seq delete
    }

    # FIXME: try $c move_seq instead

    return $rpos
}

proc editor_edit_base {w call where} {
    upvar $w opt

    set io [$w io]
    upvar $opt(top) top
    
    if {$where == "" || [$io read_only]} {
	bell
	return
    }

    foreach {type rec pos} $where break;
    if {$type == 18} {
	set seq [$io get_sequence $rec]
	foreach {old_call old_conf} [$seq get_base $pos] break
	$seq replace_base $pos $call 100
	$seq delete

	foreach {type rec _pos} [$w get_cursor relative] break
	$w cursor_right
	store_undo $w \
	    [list \
		 [list C_SET $type $rec $pos] \
		 [list B_REP $rec $pos $old_call $old_conf] ] {}
    }

    $w redraw
}

# End 0 means insert to left side, keeping right alignment static.
# End 1 is an insertion to right, keeping left side static.
proc editor_insert_gap {w where {end 1}} {
    upvar $w opt

    set io [$w io]

    if {$where == "" || [$io read_only]} {
	bell
	return
    }

    foreach {type rec pos} $where break;
    if {$type == 18} {
	set seq [$io get_sequence $rec]
	$seq insert_base $pos * -1
	set cmp [expr {([$seq get_length] < 0) ^ [$seq get_orient]}] 
	$seq delete
	set undo ""

	if {$cmp} {
	    set rpos [editor_shift_seq $io $rec 1 0]
	    lappend undo [list B_MOVE $rec $rpos]
	}

	lappend undo \
	    [list C_SET $type $rec $pos] \
	    [list B_DEL $rec $pos]
	eval lprepend undo [editor_handle_tag_insertion $io $rec $pos]

	if {$end == 0} {
	    # Also shift sequence left one base.
	    set rpos [editor_shift_seq $io $rec -1]
	    lprepend undo \
		[list T_MOVE $rec 1] \
		[list B_MOVE $rec $rpos]
	}
	
	store_undo $w $undo {}
    } else {
	set contig [$io get_contig $rec]
	$contig insert_base $pos
	$contig delete

	store_undo $w \
	    [list \
		 [list C_SET $type $rec $pos] \
		 [list C_DEL $rec $pos] ] {}
    }
    #$w cursor_right

    # Force the editing cursor to be visible
    eval $w set_cursor [$w get_cursor relative] 1
    
    $w redraw
}

# Handles the shifting of tags when we remove from a sequence srec at
# position pos.
#
# Returns a portion of Undo list to reverse the process
proc editor_handle_tag_deletion {io srec pos} {
    set s [$io get_sequence $srec]
    set crec [$s get_contig]
    set lpos [$s get_position]
    set rpos [expr {$lpos + abs([$s get_length])}]
    $s delete

    incr pos $lpos
    #puts "DEL $lpos..$rpos $pos in $crec"

    set c [$io get_contig $crec]
    set annos [$c anno_in_range $pos $rpos]
    $c delete

    set undo ""

    foreach anno $annos {
	if {[lindex $anno 8] != $srec} continue
	foreach {astart aend arec} $anno break
	#puts -nonewline "#$arec\t$astart..$aend"

	set a [$io get_anno_ele $arec]
	set d(type)   [$a get_type]
	foreach {d(start) d(end)} [$a get_position] break;
	set d(otype)  [$a get_obj_type]
	set d(orec)   [$a get_obj_rec]
	set d(anno)   [$a get_comment]
	set d(strand) [lsearch -exact {+ - . ?} [$a get_direction]]
	set d(rec)    $arec

	if {$astart > $pos} {
	    #puts "Move left"
	    incr astart -1
	    incr aend   -1
	    incr astart [expr {-1*$lpos}]
	    incr aend   [expr {-1*$lpos}]
	    $a move $astart $aend
	    $a delete

	    lappend undo [list T_MOD $arec [array get d]]
	} elseif {$aend >= $pos} {
	    incr aend   -1
	    incr astart [expr {-1*$lpos}]
	    incr aend   [expr {-1*$lpos}]
	    if {$aend - $astart < 0} {
		#puts "Delete"
		$a remove
		lappend undo [list T_NEW [array get d]]
	    } else {
		#puts "Shrink"
		$a move $astart $aend
		$a delete
		lappend undo [list T_MOD $arec [array get d]]
	    }
	} else {
	    #puts "To left => do nothing"
	    $a delete
	}
    }

    return $undo
}

# Handles the insertion to and movement of tags when we insert to a sequence
# at position 'pos'.
#
# Returns a portion of Undo list to reverse the process
proc editor_handle_tag_insertion {io srec pos} {
    set s [$io get_sequence $srec]
    set crec [$s get_contig]
    set lpos [$s get_position]
    set rpos [expr {$lpos + abs([$s get_length])}]
    $s delete

    incr pos $lpos
    #puts "INS $lpos..$rpos $pos in $crec"

    set c [$io get_contig $crec]
    set annos [$c anno_in_range $pos $rpos]
    $c delete

    set undo ""

    foreach anno $annos {
	if {[lindex $anno 8] != $srec} continue
	foreach {astart aend arec} $anno break
	#puts -nonewline "#$arec\t$astart..$aend"

	set a [$io get_anno_ele $arec]
	set d(type)   [$a get_type]
	foreach {d(start) d(end)} [$a get_position] break;
	set d(otype)  [$a get_obj_type]
	set d(orec)   [$a get_obj_rec]
	set d(anno)   [$a get_comment]
	set d(strand) [lsearch -exact {+ - . ?} [$a get_direction]]
	set d(rec)    $arec

	if {$astart >= $pos} {
	    #puts "Move right"
	    incr astart 1
	    incr aend   1
	    incr astart [expr {-1*$lpos}]
	    incr aend   [expr {-1*$lpos}]
	    $a move $astart $aend
	    lappend undo [list T_MOD $arec [array get d]]
	} elseif {$aend >= $pos} {
	    incr aend   1
	    incr astart [expr {-1*$lpos}]
	    incr aend   [expr {-1*$lpos}]
	    #puts "Grow"
	    $a move $astart $aend
	    lappend undo [list T_MOD $arec [array get d]]
	} else {
	    #puts "To left => do nothing"
	}

	$a delete
    }

    return $undo
}

# Dir is 0 for delete base left of cursor (standard Backspace functionality)
# Dir is 1 for delete base under cursor (standard Delete key)
#
# Do not confuse with end 0/1, which controls whether the left or right
# portion of the alignment is fixed. Also see editor_insert_gap.
proc editor_delete_base {w where {end 1} {dir 0} {powerup 0}} {
    upvar $w opt

    set io [$w io]

    if {$where == "" || [$io read_only]} {
	bell
	return
    }

    foreach {type rec pos} $where break;

    # Delete base to left and move cursor, vs delete under cursor and move seq
    if {$dir == 0} {
	$w cursor_left
	incr pos -1
	#editor_move_seq $w $where 1
    }

    if {$type == 18} {
	set seq [$io get_sequence $rec]
	foreach {old_base old_conf} [$seq get_base $pos] break
	foreach {l0 r0} [$seq get_clips] break;
	set cmp [expr {([$seq get_length] < 0) ^ [$seq get_orient]}] 

	if {$old_base != "*" && !$powerup} {
	    bell
	    if {$dir == 0} {
		$w cursor_right
	    }
	    return
	}

	set undo ""
	if {$cmp} {
	    $seq move_annos 1
	    lappend undo [list T_MOVE $rec -1]
	}

	$seq delete_base $pos


	# Identify if we're in a situation where we need to undo the clip
	# points too. This occurs when, for example, we delete the last base
	# of the left-cutoff data. Undoing this would be an insertion to the
	# start of the used portion, making that base now visible instead.
	# We use a belt and braces method by temporarily adding a base back
	# to check our clip points are consistent
	$seq insert_base $pos A 0
	foreach {l1 r1} [$seq get_clips] break;
	$seq delete_base $pos

	$seq delete

	lprepend undo [list C_SET $type $rec [expr {$pos+1}]] \
	              [list B_INS $rec $pos $old_base $old_conf]

	eval lprepend undo [editor_handle_tag_deletion $io $rec $pos]

	if {$end == 0} {
	    # Move annos
	    set seq [$io get_sequence $rec]
	    $seq move_annos 2
	    $seq delete
	    lprepend undo [list T_MOVE $rec -2]

	    # Also shift sequence right one base.
	    set rpos [editor_shift_seq $io $rec 1]
	    
	    if {$l0 != $l1 || $r0 != $r1} {
		lprepend undo \
		    [list B_CUT $rec $l0 $r0] \
		    [list T_MOVE $rec 1] \
		    [list B_MOVE $rec $rpos]
	    } else {
		lprepend undo \
		    [list T_MOVE $rec 1] \
		    [list B_MOVE $rec $rpos]
	    }
	} else {
	    if {$l0 != $l1 || $r0 != $r1} {
		lprepend undo [list B_CUT $rec $l0 $r0]
	    }
	}

	# FIXME: to do...
	#
	# Also check if the deletion causes the contig extents to change.
	# If so, we need to check consensus annotation locations and amend.

	store_undo $w $undo {}
    } else {
	set contig [$io get_contig $rec]

	set cons [calc_consensus -io $io -contigs "{=$rec $pos $pos}"]
	if {$cons != "*" && !$powerup} {
	    bell
	    $w cursor_right
	    $contig delete
	    return
	}

	# get pileup at consensus position so we can restore it
	set pileup [$contig get_pileup $pos]

	$contig delete_base $pos
	$contig delete

	store_undo $w \
	    [list \
		 [list C_SET $type $rec [expr {$pos+1}]] \
		 [list C_INS $rec $pos $pileup]] {}
    }
    
    $w redraw
}

proc editor_shift {w where dir} {
    upvar $w opt

    set io [$w io]

    if {$where == "" || [$io read_only]} {
	bell
	return
    }

    foreach {type rec pos} $where break;
    if {$type != 17} return

    # Undo of shift-left is *mostly* shift-right, and vice versa.
    # Where it fails is if a read starts at that point. Eg:
    #
    #     cursor pos
    #     v
    #  --------        A
    #     --------     B
    #        --------  C
    #
    # control-left here moves seq C but not seq B. control-right
    # (assumed undo) moves both B+C.
    # Solution to shift-left at pos is to shift-right at pos+1 and also
    # to manually shift-right any sequences which weren't at pos before
    # the shift-left but are now. Similarly for opposite shifts.

    set contig [$io get_contig $rec]
    set seqs [$contig seqs_in_range [expr {$pos-1}] [expr {$pos+1}]]

    if {$dir == 1} {
	set anno {}
    } else {
	set anno [$contig anno_in_range [expr {$pos-1}] [expr {$pos+1}]]

	set cons_anno {}
	foreach a $anno {
	    # Consensus anno only
	    if {[lindex $a 12] != 0} continue

	    # And only those starting at this base, which won't be shifted.
	    if {[lindex $a 0] != $pos} continue

	    lappend cons_anno [lindex $a 2]
	}
	set anno $cons_anno
    }

    set shifted [$contig shift_base $pos $dir]
    $contig delete

    if {$shifted > 0} {
	store_undo $w \
	    [list \
		 [list C_SET $type $rec $pos] \
		 [list C_SHIFT $rec $pos [expr {-$dir}] $seqs $anno] ] {}

	#$w cursor_right
	$w redraw
    } else {
	bell
    }
}

proc editor_set_confidence {w where qual} {
    upvar $w opt

    set io [$w io]
    upvar $opt(top) top
    
    if {$where == "" || [$io read_only]} {
	bell
	return
    }

    foreach {type rec pos} $where break;
    if {$type == 18} {
	set seq [$io get_sequence $rec]
	foreach {old_call old_conf} [$seq get_base $pos] break
	$seq replace_base $pos $old_call $qual
	$seq delete

	foreach {type rec _pos} [$w get_cursor relative] break
	$w cursor_right
	store_undo $w \
	    [list \
		 [list C_SET $type $rec $pos] \
		 [list B_REP $rec $pos $old_call $old_conf] ] {}
    }

    $w redraw
}

proc editor_increment_confidence {w where amount} {
    upvar $w opt

    set io [$w io]
    upvar $opt(top) top
    
    if {$where == "" || [$io read_only]} {
	bell
	return
    }

    foreach {type rec pos} $where break;
    if {$type == 18} {
	set seq [$io get_sequence $rec]
	foreach {old_call old_conf} [$seq get_base $pos] break
	set new_conf [expr {$old_conf+$amount}]
	if {$new_conf > 100} {set new_conf 100}
	if {$new_conf <   0} {set new_conf   0}
	$seq replace_base $pos $old_call $new_conf
	$seq delete

	foreach {type rec _pos} [$w get_cursor relative] break
	store_undo $w \
	    [list \
		 [list C_SET $type $rec $pos] \
		 [list B_REP $rec $pos $old_call $old_conf] ] {}
    }

    $w redraw
    update_brief $w
}

proc editor_move_seq {w where direction} {
    upvar $w opt

    set io [$w io]

    if {$where == "" || [$io read_only]} {
	bell
	return
    }

    foreach {type rec _pos} $where break;
    if {$type == 17} {
	# Consensus
	return [editor_shift $w $where $direction]
    }

    if {$type != 18} {
	# sequences only
	bell
	return
    }

    set seq  [$io get_sequence $rec]
    set cnum [$seq get_contig]
    set pos  [$seq get_position]
    $seq delete

    set upos $pos; # copy for undo

    set c [$io get_contig $cnum]
    foreach {pair_rec flags} [$c remove_sequence $rec] break;
    incr pos $direction
    $c add_sequence $rec $pos $pair_rec $flags
    $c delete

    set seq [$io get_sequence $rec]
    $seq move_annos $direction
    $seq delete

    store_undo $w \
	[list \
	     [list C_SET $type $rec $_pos] \
	     [list T_MOVE $rec [expr {-$direction}]] \
	     [list B_MOVE $rec $upos] ] {}

    $w redraw
}

proc editor_clip_seq {w where end} {
    upvar $w opt

    set io [$w io]

    if {$where == "" || [$io read_only]} {
	bell
	return
    }

    foreach {type rec pos} $where break;
    if {$type == 17} {
	return [editor_clip_contig $w $where $end]
    }

    if {$type != 18} {
	# sequences only
	bell
	return
    }

    set seq  [$io get_sequence $rec]
    set orient [$seq get_orient]

    foreach {l r} [$seq get_clips] break;
    set ul $l
    set ur $r

    switch $end {
	l {
	    if {$orient} {
		set r [expr {abs([$seq get_length])-$pos}]
	    } else {
		set l [expr {$pos+1}]
	    }
	}
	r {
	    if {$orient} {
		set l [expr {abs([$seq get_length])-$pos+1}]
	    } else {
		set r $pos
	    }
	}
    }

    if {$l <= $r} {
	$seq set_clips $l $r
    } else {
	bell
    }
    $seq delete

    store_undo $w \
	[list \
	     [list B_CUT $rec $ul $ur] \
	     [list C_SET $type $rec $pos]] {}
    $w redraw
}

proc editor_clip_contig {w where end} {
    upvar $w opt

    set io [$w io]

    foreach {type crec pos} $where break;
    set c [$io get_contig $crec]
    set seqs [$c seqs_in_range $pos $pos]
    set undo {}

    set min $pos
    set max $pos
    foreach x $seqs {
	set dir [lindex $x 14]
	set st  [lindex $x 0]
	set en  [lindex $x 1]
	set rec [lindex $x 2]

	set min [expr {$st<$min?$st:$min}]
	set max [expr {$en>$max?$en:$max}]

	set s [$io get_seq $rec]
	foreach {l r} [$s get_clips] break;
	lappend undo $rec $l $r

	if {([$s get_length] < 0) ^ $dir} {
	    set l2 [expr {abs([$s get_length])-$r+1}]
	    set r  [expr {abs([$s get_length])-$l+1}]
	    set l  $l2
	}

	set l [expr {$l+$st-1}]
	set r [expr {$r+$st-1}]

	if {$l < $pos && $r >= $pos} {
	    if {$end == "l"} {
		set l $pos
	    } else {
		set r [expr {$pos-1}]
	    }

	    set l [expr {$l + 1 - ($st)}]
	    set r [expr {$r + 1 - ($st)}]
	    if {([$s get_length] < 0) ^ $dir} {
		set l2 [expr {abs([$s get_length])-$r+1}]
		set r  [expr {abs([$s get_length])-$l+1}]
		set l  $l2
	    }

	    #puts "$s set_clips_no_invalidate $l $r: extents $st $en"
	    $s set_clips_no_invalidate $l $r
	}
	$s delete
    }

    #puts "invalidate from $min to $max"
    $c invalidate_consensus $min $max
    $c delete

    store_undo $w \
	[list \
	     [list C_CUT $crec $min $max $undo] \
	     [list C_SET $type $crec $pos]] {}

    $w redraw
}

# Updates the editor status line for editor $w.
# x and y are optional, but if set specify the location to display (eg
# for mouse motion events).
proc update_brief {w {name 0} {x {}} {y {}}} {
    global gap5_defs
    global $w

    if {$x != "" && $y != ""} {
	foreach {type rec pos} [$w get_number $x $y] break
    } else {
	foreach {type rec pos} [$w get_number] break
    }

    if {$name} {
	set w [name2ed $w]
	global $w
    }

    if {![info exists type]} {
	set w [set ${w}(top)]
	global $w
	set ${w}(Status) ""
	return
    }

    if {$name} {
	switch $type {
	    18 {
		set msg [$w get_seq_status $type $rec $pos \
			     [keylget gap5_defs READ_BRIEF_FORMAT]]
	    }
	    17 {
		set msg [$w get_seq_status $type $rec $pos \
			     [keylget gap5_defs CONTIG_BRIEF_FORMAT]]
	    }
	    21 {
		set msg "tag in name?"
	    }
	    default {
		set msg "Unknown data type $type"
	    }
	}
    } else {
	switch $type {
	    18 {
		set msg [$w get_seq_status $type $rec $pos \
			     [keylget gap5_defs BASE_BRIEF_FORMAT1]]
	    }
	    17 {
		set msg [$w get_seq_status $type $rec $pos \
			     [keylget gap5_defs BASE_BRIEF_FORMAT2]]
	    }
	    21 {
		set msg [$w get_seq_status $type $rec $pos \
			     [keylget gap5_defs TAG_BRIEF_FORMAT]]
		regsub -all "\n" $msg { // } msg
	    }
	    default {
		set msg "Unknown data type $type"
	    }
	}
    }

    # Squash white-spacing for ease of display
    regsub -all {[\n\t ]+} $msg " " msg

    set w [set ${w}(top)]
    global $w
    set ${w}(Status) $msg
}

proc editor_name_select {w where {num 0}} {
    global $w

    if {$where == ""} return

    # Add name to XA_PRIMARY selection
    foreach {type rec pos} $where break;
    if {$num} {
	set name "#$rec"
    } else {
	switch $type {
	    18 {
		set seq [[$w io] get_sequence $rec]
		set name [$seq get_name]
		$seq delete
	    }
	    17 {
		set c [[$w io] get_contig $rec]
		set name [$c get_name]
		$c delete
	    }
	}
    }

    # FIXME: underline name too, via $ed call?
#    [set ${w}(ed)] ...

    # Ignore parameters for now. Assume reading name length <= maxbytes.
    catch {rename editor_selection_handler {}}
    proc editor_selection_handler {offset maxbytes} [list return $name]

    selection own $w
    selection handle $w editor_selection_handler

    # For Windows...
    clipboard clear
    clipboard append $name
}

proc editor_set_start {ed} {
    set io [$ed io]
    set c [$io get_contig [$ed contig_rec]]

    set w $ed.move
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Set Contig Start"

    xentry $w.pos \
	-label "Visible start position" \
	-default [$c get_visible_start] \
	-type int
    okcancelhelp $w.ok \
	-ok_command "editor_set_start2 $ed $w" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap5 {Editor-Move}"

    pack $w.pos -side top -fill both
    pack $w.ok -side bottom -fill both
}

proc editor_set_start2 {ed w} {
    set pos [$w.pos get]
    if {[regexp {^[-+]?[0-9]+$} $pos] == 0} { 
	bell
	return
    }

    destroy $w

    set io [$ed io]
    set c [$io get_contig [$ed contig_rec]]
    set old_pos [$c get_visible_start]
    $c set_visible_start $pos

    # Adjust scrollbar and editor cursor
    $ed xview [expr {[$ed xview]-($old_pos - $pos)}]

    foreach {type rec apos} [$ed get_cursor relative] break
    if {$rec == [$ed contig_rec]} {
	$ed set_cursor $type $rec [expr {$apos - ($old_pos - $pos)}]
    }

    # Force relative to absolute cursor map to update
    eval $ed set_cursor [$ed get_cursor relative] 0

    store_undo $ed \
	[list \
	     [list C_MOVE [$ed contig_rec] $old_pos] \
	     [list C_SET $type $rec $apos]] {}
}

proc editor_set_name {ed} {
    set io [[$ed io] base]
    set c [$io get_contig [$ed contig_rec]]

    set w $ed.name
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Set Contig Name"

    xentry $w.name \
	-label "Contig name" \
	-default [$c get_name]
    okcancelhelp $w.ok \
	-ok_command "editor_set_name2 $ed $w" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap5 {Editor-Rename}"

    pack $w.name -side top -fill both
    pack $w.ok -side bottom -fill both
}

proc editor_set_name2 {ed w} {
    set io [[$ed io] base]

    set nm [$w.name get]
    if {$nm == ""} {
	bell
	return
    }

    set c [$io get_contig [$ed contig_rec]]
    set old_name [$c get_name]
    $c delete

    if {![contig_rename $io [$ed contig_rec] $nm $w]} {
	return
    }

    store_undo $ed \
	[list \
	     [list C_RENAME [$ed contig_rec] $old_name $nm]] {}
}

proc editor_template_display {ed} {
    set io [[$ed io] base]
    set xpos [expr {[$ed xview]-13000}]; # a hack as T.disp pos is left edge
    CreateTemplateDisplay $io [$ed contig_rec] $xpos

    # Force cursor to be visible
    after idle "after 100 {$ed set_cursor [$ed get_cursor relative] 0}"
}

#-----------------------------------------------------------------------------
# Align cutoff data

# Exists as "string reverse" in tcl 8.5
proc string_reverse {str} {
    set rstr ""
    for {set i [expr {[string length $str]-1}]} {$i >=0} {incr i -1} {
	append rstr [string index $str $i]
    }
    return $rstr
}

proc editor_align_cutoff {ed} {
    global gap5_defs

    set io [$ed io]

    if {[$io read_only]} {
	bell
	return
    }

    # Fetch start/end relative to sequence start
    foreach {type rec start end} [$ed select get] break
    if {$type != 18} {
	# Only work on sequences
	bell
	return
    }

    if {$end < $start} {
	set tmp $end
	set end $start
	set start $tmp
    }

    # Get absolute position in consensus
    set s [$io get_sequence $rec]
    set pos [$s get_position]
    $s delete
    set cstart [expr {$start+$pos}]
    set cend   [expr {$end  +$pos}]
    set cons [calc_consensus \
		  -io $io \
		  -contigs [list [list =[$ed contig_rec] $cstart $cend]]]
    
    #puts "Sequence selected = '$seq', pos=$start..$end"
    #puts "Cons =              '$cons', pos=$cstart..$cend"

    # Get sequence, in contig orientation, with clip points
    set s [$io get_sequence $rec]
    set seq [$s get_seq]
    foreach {sleft sright} [$s get_clips] break
    if {[$s get_orient] != 0} {
	set orient 1
	set seq [string map -nocase {A T C G G C T A} [string_reverse $seq]]
	set oleft  [expr {abs([$s get_length])-$sright+1}]
	set oright [expr {abs([$s get_length])-$sleft+1}]
    } else {
	# Contig orient left/right
	set orient 0
	set oleft  $sleft
	set oright $sright
    }

    incr oleft -1
    incr oright -1
    #puts "$start..$end vs $oleft..$oright (DIR $sleft,$sright)"

    # Undo stack
    foreach {_type _rec _pos} [$ed get_cursor relative] break
    set undo {}
    lappend undo [list C_SET $_type $_rec $_pos]
    lappend undo [list B_CUT $rec $sleft $sright]

    # Extend read to include the new underlined region.
    # Somewhat fiddly due to the way complemented reads work.
    if {$start < $oleft && $end > $oright} {
	# extended both
	set lmode 1
	set rmode 1
	if {$orient} {
	    set sleft [expr {abs([$s get_length])-$end}]
	    set sright [expr {abs([$s get_length])-$start}]
	} else {
	    set sleft [expr {$start+1}]
	    set sright [expr {$end+1}]
	}
	#puts "new clips $start..$end, $sleft..$sright"
	$s set_clips $sleft $sright
    } elseif {$start < $oleft} {
	# extended left
	set lmode 1
	set rmode 0
	if {$orient} {
	    set sright [expr {abs([$s get_length])-$start}]
	} else {
	    set sleft [expr {$start+1}]
	}
	#puts "new clips $start..$oright, $sleft..$sright"
	$s set_clips $sleft $sright
    } elseif {$end > $oright} {
	# extended right
	set lmode 0
	set rmode 1
	if {$orient} {
	    set sleft [expr {abs([$s get_length])-$end}]
	} else {
	    set sright [expr {$end+1}]
	}
	#puts "new clips $oleft..$end, $sleft..$sright"
	$s set_clips $sleft $sright
    } else {
	# Internal alignment
	set lmode 0
	set rmode 0
	#puts "no new clips"
    }
    $s delete

    #puts "Clipped region=[expr $oleft+[$s get_position]]..[expr $oright+[$s get_position]]"
    set seq [string range $seq $start $end]

    # Perform alignment
    foreach {cons_a seq_a} [align_seqs -seq1 $cons -seq2 $seq] break
    #puts "Orig $cons"
    #puts "Cons $cons_a"
    #puts "Seq  $seq_a"
    #puts "Orig $seq"

    set orig_len  [string length $cons]
    set align_len [string length $cons_a]

    # Edit alignment.
    #
    # Loop: i  = pos in aligned strings.
    #       rc = relative pos in contig (relative to start of $cons)
    #       rs = relative pos in sequence
    #       np = number of additional pads in consensus
    set c [$io get_contig [$ed contig_rec]]
    for {set np 0;set i 0; set rc 0; set rs 0} {$i < $align_len} {} {
	set c1 [string index $cons $rc]
	set c2 [string index $cons_a $i]
	set s1 [string index $seq $rs]
	set s2 [string index $seq_a $i]

	#puts "$i $rc $rs $c1/$c2 $s2/$s1"

	# We can get ins/del to read, and ins (no del) to cons.
	if {$s1 == $s2} {
	    incr rs
	    incr i
	} elseif {$s2 == "*"} {
	    # Ins to read
	    #puts "INS read $rs/[expr {$rs+$start}]"
	    set s [$io get_sequence $rec]
	    $s insert_base [expr {$rs+$start}] * -1
	    $s delete

	    lappend undo [list B_DEL $rec [expr {$rs+$start}]]

	    incr start
	    incr i
	} elseif {$s1 == "*"} {
	    # Del from read
	    #puts "DEL read $rs/[expr {$rs+$start}]"
	    set s [$io get_sequence $rec]
	    $s delete_base [expr {$rs+$start}]
	    $s delete

	    lappend undo [list B_INS $rec [expr {$rs+$start}] * -1]

	    incr start -1
	    incr rs
	    continue
	} else {
	    puts "ERROR $s1/$s2"
	    break
	}

	if {$c1 == $c2} {
	    incr rc
	} elseif {$c2 == "*"} {
	    # Ins to cons
	    #puts "INS cons $rc/[expr {$rc+$cstart+$np}]"
	    $c insert_base [expr {$rc+$cstart+$np}] * -1
	    set s [$io get_sequence $rec]
	    $s delete_base [expr {$rs+$start-1}]
	    $s delete

	    lappend undo [list C_DEL [$ed contig_rec] [expr {$rc+$cstart+$np}] * -1]
	    lappend undo [list B_INS $rec [expr {$rs+$start-1}] * -1]

	    incr np
	} else {
	    puts "ERROR $c1/$c2"
	    break
	}
    }

    $c delete

    $ed clear_visibility_cache
    $ed redraw

    set i [llength $undo]
    set rev {}
    while {[incr i -1] >= 0} {
	lappend rev [lindex $undo $i]
    }
    #puts $undo
    #puts $rev

    store_undo $ed $rev {}
}

#-----------------------------------------------------------------------------
# Break contig

proc editor_break_contig {ed} {
    global gap5_defs

    set io [$ed io]

    if {[$io read_only]} {
	bell
	return
    }

    # show cursor
    eval $ed set_cursor [$ed get_cursor relative] 1

    # Get reading position
    foreach {type rec pos} [$ed get_cursor relative] break
    if {$type == 18} {
	# Coord is clipped start point of this read
	set arec [$ed contig_rec]
	set s [$io get_sequence $rec]
	#set apos [$s get_position]
	set apos [$s get_clipped_position]
	$s delete
    } else {
	# On consensus, coord is this pos itself
	foreach {atype arec apos} [$ed get_cursor absolute] break
    }

    upvar contigIO_$arec cio
    
    set ret [tk_messageBox \
		 -icon  question \
		 -title "Break contig" \
		 -message "This will save changes and break the contig in two, with all reads starting on or to the right of consensus position $apos being in the second contig. Continue?" \
		 -type yesno \
		 -parent $ed]

    if {$ret == "no"} {
	return
    }

    $ed save

    set io $cio(base)

    if {![quit_displays -io $io -msg "break contig"]} {
	# Someone's too busy to shutdown?
	ClearBusy
	return
    }

    set cr [break_contig \
		-io $io \
		-contig $arec \
		-pos $apos \
		-break_holes [keylget gap5_defs CONTIG_EDITOR.REMOVE_HOLES]]
    ContigSelector $io
    ContigInitReg $io

    edit_contig -io $io -contig $arec -pos $apos
    edit_contig -io $io -contig $cr   -pos 1
}

#-----------------------------------------------------------------------------
# Tag editor windows

# Functions to make tag edits or to be called by the undo/redo stack.
proc U_tag_change {w rec new_a} {
    set io [$w io]

    #-- Get existing tag
    set old_a ""
    if {$rec != -1} {
	set tag [$io get_anno_ele $rec]
	set d(type)    [$tag get_type]
	foreach {d(start) d(end)} [$tag get_position] break;
	set d(otype)   [$tag get_obj_type]
	set d(orec)    [$tag get_obj_rec]
	set d(anno)    [$tag get_comment]
	set d(default) "?"
	set d(strand)  [lsearch -exact {+ - . ?} [$tag get_direction]]
	set d(rec)     $rec

	set old_a [array get d]
	unset d
    }

    #-- Create, Modify or Delete existing tag
    array set d $new_a

    if {$new_a == ""} {
	# Delete
	$tag remove

	store_undo $w \
	    [list \
		 [list T_NEW $old_a]] {}

#	    [list U_tag_change $w -1 $old_a] \
	    [list U_tag_change $w $rec ""]

    } elseif {$rec == -1} {
	if {$d(start) > $d(end)} {
	    set t $d(start)
	    set d(start) $d(end)
	    set d(end) $t

	}
	# Clip coords to consensus
	if {$d(otype) == 17} {
	    set c [$io get_contig [$w contig_rec]]
	    set cstart [$c get_visible_start]
	    set cend   [$c get_visible_end]
	    $c delete

	    if {$d(start)  < $cstart} {
		set d(start) $cstart
	    }
	    if {$d(end)  > $cend} {
		set d(end) $cend
	    }
	}

	# Create
	if {$d(start) <= $d(end)} {
	    set rec [$io new_anno_ele $d(otype) $d(orec) $d(start) $d(end)]
	    set t [$io get_anno_ele $rec]
	    $t set_comment $d(anno)
	    $t set_type $d(type)
	    $t set_direction [string index "+-.?" $d(strand)]
	    $t delete
	
	    store_undo $w \
		[list \
		     [list T_DEL $rec]] {}
	} else {
	    # Clipped and still negative size => entirely off contig
	    bell;
	}

    } else {
	array set ot $old_a

	# Modify
	if {[$tag get_comment] != $d(anno)} {
	    $tag set_comment $d(anno)
	}
	if {[$tag get_type] != $d(type)} {
	    $tag set_type $d(type)
	}
	if {[lsearch -exact {+ - . ?} [$tag get_direction]] != $d(strand)} {
	    $tag set_direction [string index "+-.?" $d(strand)]
	}

	if {[info exists d(start)] && [info exists d(end)]} {
	    if {$ot(start) != $d(start) ||
		$ot(end)   != $d(end)   ||
		$ot(otype) != $d(otype) ||
		$ot(orec)  != $d(orec)} {
		$tag set_obj_rec  $d(orec)
		$tag set_obj_type $d(otype)
		puts ""
		puts "Undo move tag"
		puts "$tag move $d(otype) $d(orec) $d(start) $d(end)"
		$tag move $d(otype) $d(orec) $d(start) $d(end)
	    }
	}

	$tag delete

	store_undo $w \
	    [list \
		 [list T_MOD $rec $old_a]] {}
    }

    unset d
}

proc tag_repopulate_menu {w} {
    foreach {rtype rrec rpos} [$w get_cursor relative] break
    foreach {atype arec apos} [$w get_cursor absolute] break

    upvar contigIO_$arec cio
    upvar $w wa

    set me $wa(top).menubar.commands.edit_tag
    $me delete 0 end
    $me configure -tearoff 0

    set md $wa(top).menubar.commands.delete_tag
    $md delete 0 end
    $md configure -tearoff 0

    # Perhaps not the most efficient approach, but it works
    set c [$cio(io) get_contig $cio(crec)]
    foreach anno [$c anno_in_range $apos $apos] {
	if {$rrec == [lindex $anno 8] || \
		($cio(crec) == $rrec && [lindex $anno 12] == 0)} {
	    foreach {start end rec itype} $anno break
	    set type ""
	    while {$itype > 0} {
		set type "[format "%c" [expr {$itype & 0xff}]]$type"
		set itype [expr {$itype >> 8}]
	    }
	    $me add command -label "$type $start..$end \#$rec" \
		-command "tag_editor_launch $w $rec"
	    $md add command -label "$type $start..$end \#$rec" \
		-command "tag_editor_delete $w $rec"
	}
    }
    $c delete
}

proc tag_editor_launch {w where} {
    if {$where == ""} {
	bell
	return
    }

    if {[llength $where] != 1} {
	foreach {type rec pos} $where break;
	if {$type != 21} return
    } else {
	set rec $where
    }

    set tag [[$w io] get_anno_ele $rec]
    global .Tag.$rec
    upvar \#0 .Tag.$rec d

    set d(strand)  [lsearch -exact {+ - . ?} [$tag get_direction]]
    set d(type)    [$tag get_type]
    set d(otype)   [$tag get_obj_type]
    set d(orec)    [$tag get_obj_rec]
    set d(anno)    [$tag get_comment]
    set d(default) "?" 
    set d(rec)     $rec

    $tag delete

    create_tag_editor $w.tag_$rec "tag_editor_callback $w $rec" .Tag.$rec
}

proc tag_editor_delete {w where} {
    if {$where == ""} {
	bell
	return
    }

    if {[llength $where] != 1} {
	foreach {type rec pos} $where break;
	if {$type != 21} return
    } else {
	set rec $where
    }

    U_tag_change $w $rec ""
#    set tag [[$w io] get_anno_ele $rec]
#    $tag remove

    $w redraw
}

proc tag_editor_create {w} {
    global $w
    set rec -1
    
    foreach {otype orec start end} [$w select get] break;

    global .Tag.$rec
    upvar \#0 .Tag.$rec d

    set d(strand) 0
    set d(type)    "COMM"
    set d(otype)   $otype
    set d(orec)    $orec
    set d(start)   $start
    set d(end)     $end
    set d(anno)    ""
    set d(default) "?"
    set d(rec)     $rec


    create_tag_editor $w.tag_$rec "tag_editor_callback $w $rec" .Tag.$rec
}

proc tag_editor_callback {w rec cmd args} {
    upvar \#0 .Tag.$rec d
    set f $w.tag_$rec
    set io [$w io]

    switch $cmd {
	"save" {
#	    if {$rec == -1} {
#		# Allocate a new item
#		set rec [$io new_anno_ele $d(otype) $d(orec) $d(start) $d(end)]
#		puts "New tag with rec $rec"
#		set t [$io get_anno_ele $rec]
#		$t set_comment $d(anno)
#		$t set_type $d(type)
#		$t delete
#	    } else {
#		set t [$io get_anno_ele $rec]
#
#		if {[$t get_comment] != $d(anno)} {
#		    $t set_comment $d(anno)
#		}
#		if {[$t get_type] != $d(type)} {
#		    $t set_type $d(type)
#		}
#		$t delete
#	    }
	    U_tag_change $w $rec [array get d]

	    $w redraw
	    destroy $f
	    unset d
	}

	"quit" {
	    destroy $f
	    unset d
	}

	"move" {
	    set w2 [selection own]
	    if {$w2 == ""} {
		set w2 $w
	    } else {
		if {[winfo class $w2] != "Editor"} {
		    bell
		    return
		}
	    }
	    if {$w2 != $w} {
		# Move to a different contig, so delete from here
		# and add to there (add => rec -1)
		U_tag_change $w $rec ""
		set rec -1
	    }
	    foreach {otype orec start end} [$w2 select get] break;
	    set d(otype) $otype
	    set d(orec)  $orec
	    set d(start) $start
	    set d(end)   $end
	    U_tag_change $w2 $rec [array get d]
	    $w2 redraw
	    if {$w != $w2} {
		$w redraw
	    }
	    destroy $f
	}

	"copy" {
	    set w2 [selection own]
	    if {$w2 == ""} {
		set w2 $w
	    } else {
		if {[winfo class $w2] != "Editor"} {
		    bell
		    return
		}
	    }
	    foreach {otype orec start end} [$w2 select get] break;
	    set d(otype) $otype
	    set d(orec)  $orec
	    set d(start) $start
	    set d(end)   $end
	    set d(rec)   -1
	    U_tag_change $w2 -1 [array get d]
	    $w2 redraw
	    destroy $f
	}

	default {
	    puts "ERROR: Unknown tag callback command: $cmd"
	}
    }
}

proc editor_menu {w x y} {
    upvar $w opt

    tk_popup $opt(top).menubar.commands \
	[expr $x+[winfo rootx $w]] [expr $y+[winfo rooty $w]]
}

#-----------------------------------------------------------------------------
# Trace display
proc show_trace {w loc} {
    set io [$w io]
    
    if {$loc == ""} return
    foreach {type rec pos} $loc break;

    if {$type == 18} {
	set s [$io get_seq $rec]
	set name [$s get_name]
	set t [trace_add $w.traces $name $w $name]
	after 100 "$t xview 0"

	global $w
	lappend ${w}(Traces) $t
    }
}


#-----------------------------------------------------------------------------
# Auto-scrolling when dragging selections to outside the editor window.
proc editor_select_scroll {e x} {
    global $e.AutoScroll
    set wid [winfo width $e]
    set dist 0
    if {$x < 0} {
	set dist [expr {$x / 10}]
    } elseif {$x > $wid} {
	set dist [expr {($x-$wid)/10}]
    }

    if {$dist} {
	if {![info exists $e.AutoScroll]} {
	    set $e.AutoScroll $dist
	    editor_autoscroll $e
	} else {
	    set $e.AutoScroll $dist
	}
    } else {
	catch {unset $e.AutoScroll}
    }
}

proc editor_autoscroll {e} {
    global $e.AutoScroll $e.AutoScrollEvent
    if {[catch {jog_editor $e [set $e.AutoScroll]}] == 0} {
	if {[set $e.AutoScroll] > 0} {
	    $e select to [lindex [$e configure -width] end]
	} else {
	    $e select to 0
	}
	set $e.AutoScrollEvent [after 50 "editor_autoscroll $e"]
    }
}

proc editor_autoscroll_cancel {e} {
    global $e.AutoScrollEvent
    catch {after cancel [set $e.AutoScrollEvent]; unset $e.AutoScrollEvent}
}


#-----------------------------------------------------------------------------
# Selecting and deselecting reads

proc editor_select_dialog {ed sel} {
    global gap5_defs
    global $ed
    global editor_right_click

    set pos [lindex [$ed get_cursor absolute] 2]

    if {$editor_right_click} {
	return [editor_select_reads $ed $sel 0 $pos $pos]
    }

    set t $ed.select_dialog
    if {[xtoplevel $t] == ""} return
    wm title $t "Select Reads"
    
    # Positional parameters
    entrybox $t.start \
	-title "From consensus base" \
	-default $pos \
	-width 10 \
	-type CheckInt

    entrybox $t.end \
	-title "To consensus base" \
	-default $pos \
	-width 10 \
	-type CheckInt

    yes_no $t.contained \
	-title "Entirely contained in range" \
	-orient horizontal \
	-bd 0 \
	-default 0

    okcancelhelp $t.ok \
	-bd 2 -relief groove \
	-ok_command "editor_select_reads $ed $sel    \
                         \[yes_no_get $t.contained\] \
                         \[entrybox_get $t.start\]   \
                         \[entrybox_get $t.end\];    \
                         destroy $t" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 {Editor-Select Seq}"

    pack $t.start $t.end $t.contained $t.ok -side top -fill both
}

proc editor_select_reads {ed sel contained start end} {
    SetBusy

    global $ed
    upvar #0 [set ${ed}(top)] opt

    set io [$ed io]
    set c [$io get_contig [$ed contig_rec]]
    set rlist {}
    foreach x [$c seqs_in_range $start $end] {
	foreach {x_st x_en x_rec} $x break

	# Adjust for clipped region
	set r [$io get_seq $x_rec]
	foreach {cl cr} [$r get_clips] break

	if {[$r get_orient] == 0} {
	    set c_st [expr {$x_st+$cl-1}]
	    set c_en [expr {$x_st+$cr-1}]
	} else {
	    set c_st [expr {$x_en-$cr+1}]
	    set c_en [expr {$x_en-$cl+1}]
	}
	$r delete

	if {$contained} {
	    if {$c_st < $start || $c_en > $end} continue
	} else {
	    if {$c_en < $start || $c_st > $end} continue
	}
	lappend rlist "#$x_rec"
    }
    $c delete

    if {$rlist != {}} {
	if {![ListExists2 $opt(OutputList)]} {
	    ListCreate2 $opt(OutputList) {} SEQID
	    InitListTrace $opt(OutputList)
	}

	UpdateReadingListItem $opt(OutputList) $rlist $sel
    }

    ClearBusy
}

#-----------------------------------------------------------------------------
# Oligo selection

# Creates the Find Primer dialogue window.
proc editor_oligo_dialog {ed} {
    global gap5_defs

    set direction    [keylget gap5_defs SELECT_OLIGOS.DIRECTION]
    set search_ahead [keylget gap5_defs SELECT_OLIGOS.SEARCH_AHEAD]
    set search_back  [keylget gap5_defs SELECT_OLIGOS.SEARCH_BACK]
    set read_length  [keylget gap5_defs SELECT_OLIGOS.READ_LENGTH]

    set t $ed.oligo_dialog
    if {[xtoplevel $t] == ""} return
    wm title $t "Find Primer-Walk"

    # Positional parameters
    frame $t.pos -bd 1 -relief groove
    radiolist $t.pos.dir \
	-title Direction \
	-bd 0 \
	-orient horizontal \
	-default [expr {3-$direction}] \
	-buttons {Backwards Forwards}

    entrybox $t.pos.e1 \
	-title "Search window bases ahead" \
	-default $search_ahead \
	-width 5 \
	-type {CheckIntRange 0 1000}

    entrybox $t.pos.e2 \
	-title "Search window bases back" \
	-default $search_back \
	-width 5 \
	-type {CheckIntRange 0 1000}

    entrybox $t.pos.e4 \
	-title "Average read length" \
	-default $read_length \
	-width 5 \
	-type {CheckIntRange 1 5000}
    pack $t.pos.e4 $t.pos.e2 $t.pos.e1 $t.pos.dir -side bottom -fill x

    # Primer3 parameters
    foreach i [keylget gap5_defs PRIMER] {
	set pdefs([lindex $i 0]) [lindex $i 1]
    }

    frame $t.p3 -bd 1 -relief groove
    frame $t.p3.tm -bd 0 -relief groove
    label $t.p3.tm.label -text "Melting temperature"
    entrybox $t.p3.tm.min \
	-default $pdefs(min_tm) \
	-title "Min" \
	-width 5 \
	-type {CheckFloatRange 0 100}

    entrybox $t.p3.tm.opt \
	-default $pdefs(opt_tm) \
	-title "Opt" \
	-width 5 \
	-type {CheckFloatRange 0 100}

    entrybox $t.p3.tm.max \
	-default $pdefs(max_tm) \
	-title "Max" \
	-width 5 \
	-type {CheckFloatRange 0 100}

    pack $t.p3.tm.label -side left 
    pack $t.p3.tm.max $t.p3.tm.opt $t.p3.tm.min -side right 

    frame $t.p3.length -bd 0 -relief groove
    label $t.p3.length.label -text "Primer length"
    entrybox $t.p3.length.min \
	-default $pdefs(min_len) \
	-title "Min"\
	-width 5 \
	-type {CheckIntRange 1 100}

    entrybox $t.p3.length.opt \
	-default $pdefs(opt_len) \
	-title "Opt"\
	-width 5 \
	-type {CheckIntRange 1 100}

    entrybox $t.p3.length.max \
	-default $pdefs(max_len) \
	-title "Max" \
	-width 5 \
	-type {CheckIntRange 1 100}

    pack $t.p3.length.label -side left
    pack $t.p3.length.max $t.p3.length.opt $t.p3.length.min -side right 

    frame $t.p3.gc -bd 0 -relief groove
    label $t.p3.gc.label -text "GC content (%)"
    entrybox $t.p3.gc.min \
	-default $pdefs(min_gc) \
	-title "Min"\
	-width 5 \
	-type {CheckIntRange 1 100}

    entrybox $t.p3.gc.opt \
	-default $pdefs(opt_gc) \
	-title "Opt"\
	-width 5 \
	-type {CheckIntRange 1 100}

    entrybox $t.p3.gc.max \
	-default $pdefs(max_gc) \
	-title "Max" \
	-width 5 \
	-type {CheckIntRange 1 100}

    yes_no $t.p3.gc_clamp \
    	    -title "GC Clamp" \
	    -orient horizontal \
	    -bd 0 \
	    -default $pdefs(gc_clamp)

    frame $t.p3.conc -bd 0
    entrybox $t.p3.conc.dna \
	-title "DNA concentration (nM)" \
	-width 5 \
	-default $pdefs(dna_conc) \
	-type {CheckFloatRange 0.01 1000000}

    entrybox $t.p3.conc.salt \
	-title "Salt concentration (mM)" \
	-width 5 \
	-default $pdefs(salt_conc) \
	-type {CheckFloatRange 10 1000}

    entrybox $t.p3.conc.mg \
	-title "Magnesium concentration (mM)" \
	-width 5 \
	-default $pdefs(mg_conc) \
	-type {CheckFloatRange 0 1000000}

    entrybox $t.p3.conc.dntp \
	-title "Total dNTP concentration (mM)" \
	-width 5 \
	-default $pdefs(dntp_conc) \
	-type {CheckFloatRange 0 1000000}

    pack $t.p3.conc.dna $t.p3.conc.salt $t.p3.conc.mg $t.p3.conc.dntp \
	-fill x -side bottom

    pack $t.p3.gc.label -side left 
    pack $t.p3.gc.max $t.p3.gc.opt $t.p3.gc.min -side right 
    pack $t.p3.conc $t.p3.tm $t.p3.length $t.p3.gc $t.p3.gc_clamp \
	-side bottom -fill x

    okcancelhelp $t.but \
	-bd 2 -relief groove \
	-ok_command "if {\[editor_oligo_report $ed $t\] == 0} {destroy $t}" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 {Editor-Primer Selection}"

    pack $t.but $t.p3 $t.pos -side bottom -fill x
}

# Called from the Find Primer dialogue.
# Returns 0 on success
#        -1 on failure (eg invalid parameters)
proc editor_oligo_report {ed t} {
    global gap5_defs

    # First convert the dialogue params into a P3 format string.
    foreach i [keylget gap5_defs PRIMER] {
	set pdefs([lindex $i 0]) [lindex $i 1]
    }

    set pdefs(min_tm) [entrybox_get $t.p3.tm.min] 
    set pdefs(opt_tm) [entrybox_get $t.p3.tm.opt]
    set pdefs(max_tm) [entrybox_get $t.p3.tm.max]
    
    set pdefs(min_len) [entrybox_get $t.p3.length.min]
    set pdefs(opt_len) [entrybox_get $t.p3.length.opt]
    set pdefs(max_len) [entrybox_get $t.p3.length.max]
    
    set pdefs(min_gc) [entrybox_get $t.p3.gc.min]
    set pdefs(opt_gc) [entrybox_get $t.p3.gc.opt]
    set pdefs(max_gc) [entrybox_get $t.p3.gc.max]

    set pdefs(gc_clamp) [yes_no_get $t.p3.gc_clamp]

    set pdefs(dna_conc)  [entrybox_get $t.p3.conc.dna]
    set pdefs(salt_conc) [entrybox_get $t.p3.conc.salt]
    set pdefs(mg_conc)   [entrybox_get $t.p3.conc.mg]
    set pdefs(dntp_conc) [entrybox_get $t.p3.conc.dntp]

    set p3_params [array get pdefs]
    eval keylset primer_defs $p3_params
    keylset gap5_defs PRIMER $primer_defs


    # Other params too
    if {[set search_ahead_val [entrybox_get $t.pos.e1]] == ""} {
	entrybox_focus $t.pos.e1; return -1
    }

    if {[set search_back_val [entrybox_get $t.pos.e2]] == ""} {
	entrybox_focus $t.pos.e2; return -1
    }

    if {[set read_length_val [entrybox_get $t.pos.e4]] == ""} {
	entrybox_focus $t.pos.e4; return -1
    }
    set direction [expr {[radiolist_get $t.pos.dir]-1}]

    keylset gap5_defs SELECT_OLIGOS.SEARCH_AHEAD $search_ahead_val
    keylset gap5_defs SELECT_OLIGOS.SEARCH_BACK  $search_back_val
    keylset gap5_defs SELECT_OLIGOS.READ_LENGTH  $read_length_val


    # Now actually do-it
    set w $ed.oligos
    xtoplevel $w
    global $w
    catch {unset $w}

    if {![winfo exists $w.list]} {
	wm geometry $w 640x400
	wm title $w "Oligos"

	# The main list
	tablelist $w.list \
	    -columns { 8 "Score" \
		      10 "Start" \
		      10 "End" \
		       8 "GC %" \
		       8 "Temperature" \
		      50 "Sequence"} \
	    -labelcommand tablelist::sortByColumn \
	    -exportselection 0 \
	    -stretch 0 \
	    -yscrollcommand [list $w.yscroll set]
	scrollbar $w.yscroll -command "$w.list yview"

	grid rowconfigure $w 0 -weight 1
	grid columnconfigure $w 0 -weight 1
	grid $w.list $w.yscroll -sticky nsew

	# Further details on the active primer
	set d [frame $w.details]
	label $d.seq_l           -text "Sequence:"                 -anchor w
	label $d.seq_r           -textvariable ${w}(sequence)      -anchor w
	label $d.temp_l          -text "Melting temp.:"            -anchor w
	label $d.temp_r          -textvariable ${w}(temperature)   -anchor e
	label $d.self_any_l      -text "Self-any:"                 -anchor w
	label $d.self_any_r      -textvariable ${w}(self_any)      -anchor e
	label $d.self_end_l      -text "Self-end:"                 -anchor w
	label $d.self_end_r      -textvariable ${w}(self_end)      -anchor e
	label $d.end_stability_l -text "End stability:"            -anchor w
	label $d.end_stability_r -textvariable ${w}(end_stability) -anchor e

	grid columnconfigure $d 0 -weight 0
	grid columnconfigure $d 1 -weight 1
	grid columnconfigure $d 2 -minsize 20
	grid columnconfigure $d 3 -weight 0
	grid columnconfigure $d 4 -weight 1

	grid $d.seq_l           -row 0 -column 0 -sticky nsew
	grid $d.seq_r           -row 0 -column 3 -columnspan 2 -sticky nsew
	grid $d.self_any_l      $d.self_any_r      x \
	     $d.temp_l          $d.temp_r          -sticky nsew
	grid $d.self_end_l      $d.self_end_r      x \
	     $d.end_stability_l $d.end_stability_r -sticky nsew

	xcombobox $d.name \
	    -text "Seq.name to tag" \
	    -textvariable ${w}(read) \
	    -valuesvariable ${w}(read_values)

	xcombobox $d.template \
	    -text "Template name" \
	    -textvariable ${w}(template) \
	    -valuesvariable ${w}(template_values)

	grid $d.name - x $d.template - -sticky nsew
	grid $d - -sticky nsew


	# Add and close buttons
	set d [frame $w.buttons -bd 2 -relief groove]
	button $d.add \
	    -text "Add annotation" \
	    -command "editor_oligo_add $ed $w $direction"
	button $d.close \
	    -text "Close" \
	    -command "unset $w; destroy $w"

	pack $d.add $d.close -side left -expand 1
	grid $d - -sticky nsew
    } else {
	$w.list delete 0 end
    }

    set oligos [$ed select_oligo $direction \
		    $search_ahead_val $search_back_val \
		    $read_length_val $p3_params]

    foreach oligo $oligos {
	array set o $oligo
	$w.list insert end [list [format %6.2f $o(quality)] \
				$o(start) \
				$o(end) \
				[format %5.1f $o(GC)] \
				[format %5.1f $o(temperature)] \
				$o(sequence)]
	set ${w}(line_$o(sequence)) $oligo
    }

    bind $w.list <<TablelistSelect>> "+editor_oligo_select $ed %W $w"

    return 0
}

proc editor_oligo_select {ed tl w} {
    global $w

    set line [$tl get [$tl curselection]]
    foreach {score start end gc temp seq} $line {}
    array set $w [set ${w}(line_$seq)]

    $ed select set $start $end

    # Obtain a list of overlapping reads. Ideally we need to find
    # overlapping templates instead - more work to do here.
    set io [$ed io]
    set c [$io get_contig [$ed contig_rec]]
    set name_list "(consensus)"
    set temp_list ""

    foreach x [$c seqs_in_range $start $end] {
	foreach {x_st x_en x_rec} $x break
	if {$x_st <= $start && $x_en >= $end} {

	    # Adjust for clipped region
	    set r [$io get_seq $x_rec]
	    foreach {cl cr} [$r get_clips] break

	    if {[$r get_orient] == 0} {
		set c_st [expr {$x_st+$cl-1}]
		set c_en [expr {$x_st+$cr-1}]
	    } else {
		set c_st [expr {$x_en-$cr+1}]
		set c_en [expr {$x_en-$cl+1}]
	    }

	    lappend ${w}(rname_[$r get_name]) $x_rec

	    if {$c_st <= $start && $c_en >= $end} {
		lappend name_list [$r get_name]
	    }
	    lappend temp_list [regsub {\.[^.]*$} [$r get_name] {}]
	    
	    $r delete
	}
    }
    $c delete

    # Update GUI
    set ${w}(read) [lindex $name_list 0]
    set ${w}(read_values)  $name_list

    set ${w}(template) [lindex $temp_list 0]
    set ${w}(template_values)  $temp_list


    global .Selection
    set .Selection $seq
    selection own .
    selection handle . editor_oligo_handle_selection

    # For Windows...
    clipboard clear
    clipboard append $seq
}

proc editor_oligo_handle_selection {offset maxbytes} {
    global .Selection
    return [string range ${.Selection} $offset [expr {$offset+$maxbytes-1}]]
}

proc editor_oligo_add {ed w direction} {
    global $w

    # Get record number(s) for name (maybe multiple seqs with this name)
    set read [set ${w}(read)]

    if {[info exists ${w}(rname_$read)]} {
	set recs [set ${w}(rname_$read)]
    } else {
	# User typed in another name not listed. We'll need to verify it.
	# FIXME: shouldn't be necessary?

	set recs {}
    }

    set use_rec ""
    set offset 0
    foreach rec $recs {
	set r [[$ed io] get_seq $rec]

	set pos [$r get_position]
	set len [$r get_length]
	foreach {cl cr} [$r get_clips] break
	$r delete

	if {$len >= 0} {
	    set r_st [expr {$pos+$cl-1}]
	    set r_en [expr {$pos+$cl-1}]
	} else {
	    set r_st [expr {$pos-$len-$cr}]
	    set r_en [expr {$pos-$len-$cl}]
	}

	# Ideal case
	if {[set ${w}(start)] >= $r_st && [set ${w}(end)] <= $r_en} {
	    set use_rec $rec
	    set offset [expr {-$pos}]
	    break;
	}

	# Next-best in cutoff data, but keep searching
	set r_st $pos
	set r_en [expr {$pos+abs($len)-1}]
	if {[set ${w}(start)] >= $r_st && [set ${w}(end)] <= $r_en} {
	    set use_rec $rec
	    set offset [expr {-$pos}]
	}
    }

    if {$use_rec == ""} {
	# Consensus
	set d(orec) [$ed contig_rec]
	set d(otype) 17
	set d(start)   [set ${w}(start)]
	set d(end)     [set ${w}(end)]
    } else {
	# Sequences
	set d(orec) $use_rec
	set d(otype) 18
	set d(start)   [expr {[set ${w}(start)]+$offset}]
	set d(end)     [expr {[set ${w}(end)]+$offset}]
    }

    set d(type)    "OLIG"
    set d(strand)  [expr {1-$direction}]

    set d(anno)    "Template	[set ${w}(template)]
Oligoname	??
GC		[set ${w}(GC)]
Temperature	[set ${w}(temperature)]
Score		[set ${w}(quality)]
Date_picked	[clock format [clock seconds]]
Sequence	[set ${w}(sequence)]"

    U_tag_change $ed -1 [array get d]
}

proc editor_move_pad {w dir {powerup 0}} {
    upvar $w opt

    set io [$w io]
    upvar $opt(top) top

    set where [$w get_number]
    if {$where == "" || [$io read_only]} {
	bell
	return
    }

    foreach {type rec pos} $where break;

    if {$type != 18} { 
	bell
	return
    }

    if {$dir == -1 && $pos == 0} {
	bell
	return
    }
    set s [$io get_sequence $rec]
    if {$dir == 1 && $pos >= [expr {abs([$s get_length])-1}]} {
	$s delete
	bell
	return
    }
    
    if {$dir == 1} {
	set p0 $pos
	set p1 [expr {$pos+1}]
    } else {
	set p0 $pos
	set p1 [expr {$pos-1}]
    }
    foreach {call0 conf0} [$s get_base $p0] break
    foreach {call1 conf1} [$s get_base $p1] break

    if {$call0 != "*" && $powerup == 0} {
	$s delete
	bell
	return
    }

    $s replace_base $p0 $call1 $conf1
    $s replace_base $p1 $call0 $conf0
    $s delete

    if {$dir == 1} {
	$w cursor_right
    } else {
	$w cursor_left
    }

    store_undo $w \
	[list \
	     [list B_REP $rec $p0 $call0 $conf0] \
	     [list B_REP $rec $p1 $call1 $conf1] \
	     [list C_SET $type $rec $pos]] {}
	    
    $w redraw
}

proc save_editor_settings {w} {
    upvar \#0 $w opt
    global gap5_defs env

    set name $w.ed1.pane.name.sheet
    set name_wid [lindex [$name configure -width] 4]

    set ed $opt(curr_editor)
    set ed_wid [lindex [$ed configure -width] 4]
    set ed_hei [lindex [$ed configure -height] 4]

    set C CONTIG_EDITOR
    keylset gap5_defs $C.DISAGREEMENTS      $opt(Disagreements)
    keylset gap5_defs $C.DISAGREE_MODE      $opt(DisagreeMode)
    keylset gap5_defs $C.DISAGREE_CASE      $opt(DisagreeCase)
    keylset gap5_defs $C.DISAGREE_QUAL      $opt(DisagreeQuality)
    keylset gap5_defs $C.PACK_SEQUENCES     $opt(PackSequences)
    keylset gap5_defs $C.SHOW_QUALITY       $opt(Quality)
    keylset gap5_defs $C.SHOW_CUTOFFS       $opt(Cutoffs)
    keylset gap5_defs $C.NAMES_WIDTH        $name_wid
    keylset gap5_defs $C.MAX_HEIGHT	    $ed_hei
    keylset gap5_defs $C.SEQ_WIDTH	    $ed_wid
    keylset gap5_defs $C.GROUP_BY_PRIMARY   $opt(GroupByPrimary)
    keylset gap5_defs $C.GROUP_BY_SECONDARY $opt(GroupBySecondary)
    keylset gap5_defs $C.GROUP_PRIMARY      $opt(GroupPrimary)
    keylset gap5_defs $C.GROUP_SECONDARY    $opt(GroupSecondary)

    # Write to the .gaprc file
    update_defs gap5_defs $env(HOME)/.gap5rc \
	$C.DISAGREEMENTS    \
	$C.DISAGREE_MODE    \
	$C.DISAGREE_CASE    \
	$C.DISAGREE_QUAL    \
	$C.PACK_SEQUENCES   \
	$C.SHOW_QUALITY     \
	$C.SHOW_CUTOFFS     \
	$C.NAMES_WIDTH	    \
	$C.MAX_HEIGHT       \
	$C.SEQ_WIDTH	    \
    	$C.GROUP_BY_PRIMARY  	\
    	$C.GROUP_BY_SECONDARY 	\
    	$C.GROUP_PRIMARY    	\
    	$C.GROUP_SECONDARY
}

#-----------------------------------------------------------------------------
# Generic bindings
bind EdNames <Any-Motion> {update_brief %W 1 @%x @%y}

bind Editor <Any-Motion> {update_brief %W 0 @%x @%y}

# Jump to read-pair
bind EdNames <<menu>> {
    global %W

    set ed [set %W(ed)]

    foreach {type rec pos} [%W get_number @%x @%y] break
    if {![info exists type]} {
	return
    }
    
    switch $type {
	17 {
	    create_popup %W.m "Commands for contig \#$rec"
	}
	18 {
	    create_popup %W.m "Commands for sequence \#$rec"
	}
	default {
	    return
	}
    }

    %W.m add command \
	-label "Copy name to clipboard" \
	-command "editor_name_select $ed {$type $rec $pos}"
    %W.m add command \
	-label "Copy #number to clipboard" \
	-command "editor_name_select $ed {$type $rec $pos} 1"
    
    if {$type != 18} {
	%W.m add separator
	%W.m add command \
	    -label "Change contig start" \
	    -command "editor_set_start $ed"
	%W.m add command \
	    -label "Change contig name" \
	    -command "editor_set_name $ed"
	tk_popup %W.m [expr %X-20] [expr %Y-10]
	return
    }

    # Create the goto... menu
    #%W.m add cascade -label "Goto..." -menu %W.m.goto
    #menu %W.m.goto

    set other_end [$ed get_template_seqs $rec]
    if {[llength $other_end] >= 1} {
	# Compute distance from this sequence to the end of this contig.
	set s [[$ed io] get_seq $rec]
	set sc [$s get_contig]
	set c [[$ed io] get_contig $sc]
	set pos [$s get_position]
	set lsize [expr {$pos-[$c get_start]}]
	set rsize [expr {[$c get_end] - ($pos+abs([$s get_length])-1)}]
	set dist [expr {$lsize < $rsize ? $lsize : $rsize}]
	$c delete
	$s delete

	%W.m add separator
	set to_join {}
	foreach orec $other_end {
	    set s [[$ed io] get_seq $orec]
	    set pos [$s get_position]

	    # Get distance of other sequence from the end of its contig
	    set sc [$s get_contig]
	    set c [[$ed io] get_contig $sc]
	    set lsize [expr {$pos-[$c get_start]}]
	    set rsize [expr {[$c get_end] - ($pos+abs([$s get_length])-1)}]
	    incr dist [expr {$lsize < $rsize ? $lsize : $rsize}]

	    # Get the base contig IO
	    set crec [$ed contig_rec]
	    upvar \#0 contigIO_$crec cio
	    set base_io $cio(base)

	    if {$sc == [$ed contig_rec]} {
		%W.m add command \
		    -label "Goto [$s get_name] (@ $pos)" \
		    -command "$ed set_cursor 18 $orec 0"
	    } else {
		set cname [$c get_name]

		%W.m add command \
		    -label "Goto [$s get_name] (Contig '$cname' @ $pos, size ~ $dist)" \
		    -command "create_or_move_editor $base_io $sc $orec 0"

		set ts [$io get_seq $rec]

		if {[$ts get_template_orient] != [$s get_template_orient]} {
		    lappend to_join \
			[list %W.m add command \
			     -label "Join to [$s get_name] (Contig '$cname' @ $pos, complemented)" \
			     -command "
                                 complement_contig \
                                  -io $base_io \
                                  -contigs =[$s get_contig];
                                 join_contig \
                                  -io $base_io \
                                  -contig   [$ed contig_rec] \
                                  -reading  \#$rec \
                                  -pos      0 \
                                  -contig2  $sc \
                                  -reading2 \#$orec \
                                  -pos2     0"]
		} else {
		    lappend to_join \
			[list %W.m add command \
			     -label "Join to [$s get_name] (Contig '$cname' @ $pos)" \
			     -command "join_contig \
                                  -io $base_io \
                                  -contig   [$ed contig_rec] \
                                  -reading  \#$rec \
                                  -pos      0 \
                                  -contig2  $sc \
                                  -reading2 \#$orec \
                                  -pos2     0"]
		}
		$ts delete
	    }

	    $s delete
	    $c delete
	}

	if {[llength $to_join] != 0} {
	    %W.m add separator
	    foreach j $to_join {
		eval $j
	    }
	}
    }

    tk_popup %W.m [expr %X-20] [expr %Y-10]

    # $ed set_cursor 18 $other_end 1
}

bind Editor <<select>> {
    focus %W

    set w [winfo toplevel %W]
    if {![string match [set ${w}(curr_editor)] %W]} {
	[set ${w}(curr_editor)] configure -hollow_cursor 1
	set ${w}(curr_editor) %W
	[set ${w}(curr_editor)] configure -hollow_cursor 0
	$w.toolbar.undo configure -state [io_undo_state [%W contig_rec]]
#	$w.toolbar.redo configure -state [io_redo_state [%W contig_rec]]
    }
    set _sel [%W get_number @%x @%y 0 1]
    if {$_sel == ""} {
	unset _sel
	return
    } else {
	eval %W set_cursor $_sel 1
	update_brief %W
    }
    
    if {[string match "17*" $_sel]} {
    	upvar \#0 $w opt
	
	foreach ed $opt(all_editors) {
	    $ed set_base_sort_point
	    $ed set_sort_order 
    	    $ed redraw
	}
    }
    
    unset _sel
    %W select clear
}

proc editor_return {w} {
    global $w
    foreach {type rec pos} [$w get_number] break;

    if {$type != 17} return

    set upos [consensus_unpadded_pos \
		  -io     [$w io] \
		  -contig $rec \
		  -pos    $pos]
    
    set w [set ${w}(top)]
    global $w
    set ${w}(Status) "Padded position $pos, unpadded $upos"
}

bind Editor <Key-Return> {editor_return %W}
catch {bind Editor <Key-KP_Enter> {editor_return %W}}

bind Editor <Key-Left>		{%W cursor_left; update_brief %W; order_update_on_consensus %W}
bind Editor <Control-Key-b>	{%W cursor_left; update_brief %W; order_update_on_consensus %W}

bind Editor <Key-Right>		{%W cursor_right; update_brief %W; order_update_on_consensus %W}
bind Editor <Control-Key-f>	{%W cursor_right; update_brief %W; order_update_on_consensus %W}

bind Editor <Key-Up>		{%W cursor_up;    update_brief %W; order_update_on_consensus %W}
bind Editor <Control-Key-p>	{%W cursor_up;    update_brief %W; order_update_on_consensus %W}

bind Editor <Key-Down>		{%W cursor_down;  update_brief %W; order_update_on_consensus %W}
bind Editor <Control-Key-n>	{%W cursor_down;  update_brief %W; order_update_on_consensus %W}

# Not all X11 servers have these keysyms
catch {
    bind Editor <Key-KP_Left>	{%W cursor_left;  update_brief %W}
    bind Editor <Key-KP_Right>	{%W cursor_right; update_brief %W}
    bind Editor <Key-KP_Up>	{%W cursor_up;    update_brief %W}
    bind Editor <Key-KP_Down>	{%W cursor_down;  update_brief %W}
}

bind Editor <Control-Key-a>	{%W read_start;   update_brief %W}
bind Editor <Alt-Key-a>		{%W read_start2;  update_brief %W}
bind Editor <Meta-Key-a>	{%W read_start2;  update_brief %W}
bind Editor <Escape><Key-a>	{%W read_start2;  update_brief %W}

bind Editor <Control-Key-e>	{%W read_end;     update_brief %W}
bind Editor <Alt-Key-e>		{%W read_end2;    update_brief %W}
bind Editor <Meta-Key-e>	{%W read_end2;    update_brief %W}
bind Editor <Escape><Key-e>	{%W read_end2;    update_brief %W}

bind Editor <Alt-Key-comma>	{%W contig_start; update_brief %W}
bind Editor <Meta-Key-comma>	{%W contig_start; update_brief %W}
bind Editor <Escape><Key-comma>	{%W contig_start; update_brief %W}

bind Editor <Alt-Key-period>	{%W contig_end;   update_brief %W}
bind Editor <Meta-Key-period>	{%W contig_end;   update_brief %W}
bind Editor <Escape><Key-period> {%W contig_end;  update_brief %W}

bind Editor <Double-1> {%W display_trace}
bind Editor <Control-Key-t> {%W display_trace}

bind Editor <Control-Key-l> {editor_move_pad %W -1}
bind Editor <Control-Key-r> {editor_move_pad %W +1}
bind Editor <Alt-Key-l> {editor_move_pad %W -1 1}
bind Editor <Alt-Key-r> {editor_move_pad %W +1 1}

# Editing commands
bind Editor <Key-a> {editor_edit_base %W a [%W get_number]}
bind Editor <Key-c> {editor_edit_base %W c [%W get_number]}
bind Editor <Key-g> {editor_edit_base %W g [%W get_number]}
bind Editor <Key-t> {editor_edit_base %W t [%W get_number]}
bind Editor <Key-u> {editor_edit_base %W t [%W get_number]}
bind Editor <Key-n> {editor_edit_base %W n [%W get_number]}
bind Editor <Key-minus> {editor_edit_base %W - [%W get_number]}
catch {bind Editor <Key-KP_Subtract> {editor_edit_base %W - [%W get_number]}}
bind Editor <Key-asterisk> {editor_edit_base %W * [%W get_number]}
catch {bind Editor <Key-KP_Multiply> {editor_edit_base %W * [%W get_number]}}
bind Editor <Key-i> {editor_insert_gap %W [%W get_number] 1}
bind Editor <Key-I> {editor_insert_gap %W [%W get_number] 0}
bind Editor <Key-Insert> {editor_insert_gap %W [%W get_number] 1}
bind Editor <Shift-Key-Insert> {editor_insert_gap %W [%W get_number] 0}
bind Editor <Key-Delete> {editor_delete_base %W [%W get_number] 1 1}
bind Editor <Key-BackSpace> {editor_delete_base %W [%W get_number] 1 0}
bind Editor <Shift-Key-Delete> {editor_delete_base %W [%W get_number] 0 1}
bind Editor <Shift-Key-BackSpace> {editor_delete_base %W [%W get_number] 0 0}
bind Editor <Control-Key-Delete> {editor_delete_base %W [%W get_number] 1 1 1}
bind Editor <Control-Key-BackSpace> {editor_delete_base %W [%W get_number] 1 0 1}
bind Editor <Shift-Control-Key-Delete> {editor_delete_base %W [%W get_number] 0 1 1}
bind Editor <Shift-Control-Key-BackSpace> {editor_delete_base %W [%W get_number] 0 0 1}

bind Editor <Key-bracketleft>  {editor_set_confidence %W [%W get_number] 0}
bind Editor <Key-bracketright> {editor_set_confidence %W [%W get_number] 100}
bind Editor <Shift-Key-Up>     {editor_increment_confidence %W [%W get_number] 1}
bind Editor <Control-Key-Up>   {editor_increment_confidence %W [%W get_number] 10}
bind Editor <Shift-Key-Down>   {editor_increment_confidence %W [%W get_number] -1}
bind Editor <Control-Key-Down> {editor_increment_confidence %W [%W get_number] -10}

bind Editor <Control-Key-Left>  {editor_move_seq %W [%W get_number] -1}
bind Editor <Control-Key-Right> {editor_move_seq %W [%W get_number]  1}

bind Editor <Control-Key-q>     {editor_toggle_annos %W}

bind Editor <Key-less>          {editor_clip_seq %W [%W get_number] l}
bind Editor <Key-greater>       {editor_clip_seq %W [%W get_number] r}

bind Editor <Control-Key-z>     {editor_undo [winfo toplevel %W]}

# MouseWheel scrolling
bind Editor  <MouseWheel> {%W yview scroll [expr {-(%D)}] units}
bind EdNames <MouseWheel> {%W yview scroll [expr {-(%D)}] units}
if {[tk windowingsystem] eq "x11"} {
    bind Editor <4>               {%W yview scroll  -1 units}
    bind Editor <5>               {%W yview scroll  +1 units}
    bind Editor <Control-4>       {%W yview scroll -10 units}
    bind Editor <Control-5>       {%W yview scroll +10 units}

    bind Editor <Shift-4>         {%W xview scroll  -1 units}
    bind Editor <Shift-5>         {%W xview scroll  +1 units}
    bind Editor <Shift-Control-4> {%W xview scroll -10 units}
    bind Editor <Shift-Control-5> {%W xview scroll +10 units}

    bind EdNames <4>               {%W yview scroll  -1 units}
    bind EdNames <5>               {%W yview scroll  +1 units}
    bind EdNames <Control-4>       {%W yview scroll -10 units}
    bind EdNames <Control-5>       {%W yview scroll +10 units}

    bind EdNames <Shift-4>         {%W xview scroll  -1 units}
    bind EdNames <Shift-5>         {%W xview scroll  +1 units}
}

# These keysyms have different names, so try both and ignore errors
catch {
bind Editor <Key-Page_Down> {%W xview scroll  +1000 units}
bind Editor <Key-Page_Up>   {%W xview scroll  -1000 units}
bind Editor <Shift-Key-Page_Down> {%W xview scroll  +10000 units}
bind Editor <Shift-Key-Page_Up>   {%W xview scroll  -10000 units}
bind Editor <Control-Key-Page_Down> {%W xview scroll  +100000 units}
bind Editor <Control-Key-Page_Up>   {%W xview scroll  -100000 units}
bind Editor <Shift-Control-Key-Page_Down> {%W xview scroll  +1000000 units}
bind Editor <Shift-Control-Key-Page_Up>   {%W xview scroll  -1000000 units}
}

catch {
bind Editor <Key-Home>			{%W read_start; update_brief %W}
bind Editor <Key-End>			{%W read_end;   update_brief %W}
bind Editor <Key-Next>                  {%W xview scroll +1 pages}
bind Editor <Key-Prior>                 {%W xview scroll -1 pages}
bind Editor <Shift-Key-Next>            {%W xview scroll +1000 units}
bind Editor <Shift-Key-Prior>           {%W xview scroll -1000 units}
bind Editor <Control-Key-Next>          {%W xview scroll +10000 units}
bind Editor <Control-Key-Prior>         {%W xview scroll -10000 units}
bind Editor <Shift-Control-Key-Next>    {%W xview scroll +100000 units}
bind Editor <Shift-Control-Key-Prior>   {%W xview scroll -100000 units}
}

# Selection control for adding tags
bind Editor <<select-drag>> {%W select to @%x; editor_select_scroll %W %x}
bind Editor <<select-to>>   {%W select to @%x}
bind Editor <<select-release>>	{editor_autoscroll_cancel %W; editor_new_select_sort %W}
bind EdNames <2> {editor_name_select %W [%W get_number @%x @%y] 1}

bind EdNames <<select>> {
    upvar #0 [winfo toplevel %W] opt

    global EdNames_select EdNames_select_last
    global EdNames_select_x EdNames_select_ctg
    set EdNames_select 1

    set where [%W get_number @%x @%y]
    if {$where == ""} return

    if {![ListExists2 $opt(OutputList)]} {
	ListCreate2 $opt(OutputList) {} SEQID
	InitListTrace $opt(OutputList)
    }

    foreach {type rec pos} $where break
#    if {$type != 18} return

    set EdNames_select [UpdateReadingListItem $opt(OutputList) [list "\#$rec"] -1]
    editor_name_select %W [%W get_number @%x @%y]

    if {$type == 18} {
	set EdNames_select_last $rec
    } else {
	set EdNames_select_last -1
    }

    set ed [name2ed %W]
    set EdNames_select_ctg [$ed contig_rec]
    set EdNames_select_x   [$ed xview]
}

bind EdNames <<select-to>> {
    upvar #0 [winfo toplevel %W] opt

    global EdNames_select EdNames_select_last
    global EdNames_select_x EdNames_select_ctg

    if {![info exists EdNames_select_last]} return
    if {$EdNames_select_last == -1} return

    # Disallow scrolling or attempts to select between two editors
    set ed [name2ed %W]
    if {$EdNames_select_ctg != [$ed contig_rec]} {bell; return}
    if {$EdNames_select_x   != [$ed xview]}      {bell; return}

    set where [%W get_number @%x @%y]
    if {$where == ""} return

    foreach {type rec pos} $where break
    if {$type != 18} return

    set recs [%W recs_between $EdNames_select_last $rec]
    lappend recs $rec; #include "to" reading

    if {![ListExists2 $opt(OutputList)]} {
	ListCreate2 $opt(OutputList) {} SEQID
	InitListTrace $opt(OutputList)
    }

    set reads ""
    foreach r $recs { lappend reads "#$r" }
    UpdateReadingListItem $opt(OutputList) $reads $EdNames_select
}

bind EdNames <<select-drag>> {
    upvar #0 [winfo toplevel %W] opt

    set where [%W get_number @%x @%y]
    if {$where == ""} return

    foreach {type rec pos} $where break
    if {$type != 18} return

    if {![ListExists2 $opt(OutputList)]} {
	ListCreate2 $opt(OutputList) {} SEQID
	InitListTrace $opt(OutputList)
    }

    global EdNames_select
    UpdateReadingListItem $opt(OutputList) [list "\#$rec"] $EdNames_select
    editor_name_select %W [%W get_number @%x @%y]
}

# Searching
bind Editor <<search>>		{
    create_search_win [set %W(top)].search "%W search" 1
}
bind Editor <<rsearch>>		{
    create_search_win [set %W(top)].search "%W search" -1
}

# Tag editing
bind Editor <Key-F11> {
    set w [%W get_xy]
    if {$w == ""} {
	%W show_cursor

	set w [%W get_xy]
	if {$w == ""} {
	    bell
	    return
	}
    }
    foreach {x y} $w break
    tag_editor_launch %W [%W get_number $x $y]
}

bind Editor <Key-F12> {
    set w [%W get_xy]
    if {$w == ""} {
	%W show_cursor
	bell
	return
    }
    foreach {x y} $w break
    tag_editor_delete %W [%W get_number $x $y]
}

set editor_right_click 0
bind Editor <<menu>> {
    set _sel [%W get_number @%x @%y 0 1]
    if {$_sel == ""} {
	unset _sel
	return
    } else {
	eval %W set_cursor $_sel
    }
    unset _sel

    global editor_right_click
    set editor_right_click 1
    tag_repopulate_menu %W
    editor_menu %W %x %y
}

bind Editor <<menu-release>> {
    puts menu_release
    global editor_right_click
    set editor_right_click 0
}

bind Editor <Control-Key-c> {
    upvar \#0 [winfo toplevel %W] opt
 
    set w $opt(curr_editor)
    set crec [$w contig_rec]
    set c [[$w io] get_contig $crec]
    vmessage "Check contig: [$c check] errs"
    $c delete
}

bind Editor <Key-F9> {editor_oligo_dialog %W}

# Tag macros
set ed_macro_keys ""
for {set i 1} {$i <= 10} {incr i} {
    bind Editor <Control-Key-F$i> "tag_macro_copy %W F$i; break"
    bind Editor <Shift-Key-F$i> "tag_macro_create %W F$i;break"
    bind Editor <F$i> "tag_macro_invoke %W F$i;%W select clear;break"
    bind EdNames <Shift-Key-F$i> \
	"tag_macro_create \[edname_to_editor %W\] F$i;break"
    bind EdNames <F$i> "tag_macro_invoke \[edname_to_editor %W\] F$i; \
                       \[edname_to_editor %W\] select clear; \
                       break"
    # If XF86_Switch_VT_1 (etc) keysyms exist then they are likely to be
    # replacing Shift-F1. This seems to be the default in some newer
    # XFree86 installations
    catch {
	bind Editor <XF86_Switch_VT_$i> "tag_macro_create %W F$i;break"
	bind EdNames <XF86_Switch_VT_$i> \
	    "tag_macro_create \[edname_to_editor %W\] F$i;break"
    }
    lappend ed_macro_keys F$i
}
tag_macro_load
