#-----------------------------------------------------------------------------
# User dialogue for starting up the editors
#-----------------------------------------------------------------------------
proc EditContig2 {io t id} {
    if {[set reading [contig_id_gel $id]] == ""} return
    if {[set crec [contig_id_rec $id]] == ""} return

    destroy $t
    edit_contig -io $io -contig $crec -reading $reading
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
	-help_command "show_help gap4 {Editor}"

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
	-contig $crec1 -reading $read1 -pos 1 -contig2 $crec2 -pos2 1
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
	-help_command "show_help gap4 {Editor-Joining}"

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
    upvar \#0 $w opt

    # Initialise the $path global array - the instance data for this widget
    foreach {arg val} $args {
	set opt($arg) $val
    }
    set opt(win) $w
    set opt(Disagreements) 0
    set opt(DisagreeMode)  1
    set opt(DisagreeCase)  1
    set opt(Status)        "--- Status info here ---"

    set opt(-io) [$opt(-io) child]

    set join [info exists opt(-contig2)]

    #set opt(contig) [contig_order_to_number -io $opt(-io) -order 0]
    set opt(contig) $opt(-contig)
    if {![info exists opt(-reading)]} { set opt(-reading) 0 }
    if {![info exists opt(-pos)]}     { set opt(-pos) 1 }
    if {$join} {
	set opt(contig2) $opt(-contig2)
	if {![info exists opt(-reading2)]} { set opt(-reading2) 0 }
	if {![info exists opt(-pos2)]}     { set opt(-pos2) 1 }
    }

    # Create the window layout
    if {![winfo exists $w]} {
	toplevel $w
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
    button $tool.undo    -text Undo -command "editor_undo $w" -state disabled
    button $tool.redo    -text Redo -command "editor_redo $w" -state disabled
    button $tool.search  -text Search
    button $tool.save -text Save -command "editor_save $w"
    wm protocol $w WM_DELETE_WINDOW "editor_exit $w"
    pack $tool.undo $tool.redo $tool.search $tool.cutoffs $tool.quality \
	-side left
    pack $tool.save -side right

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
    }

    # The editor(s) itself
    if {$join} {
	set pane0 $w.ed0
	set e [editor_pane $w $pane0 1 2 opt]
	set opt(editor2) $e
    }
    set pane1 $w.ed1
    set e [editor_pane $w $pane1 0 "" opt]
    set opt(editor1) $e

    # Difference bar for the join editor
    if {$join} {
	set diffs [diff_pane $w.diffs]
	$e link_to $opt(editor2) $diffs.seq.sheet
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
	lappend opt(panes) $diffs $pane0.pane
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
    array set font [font metrics sheet_font]
#    wm grid $w [winfo width .e1] [winfo height .e1] \
#	[font measure sheet_font A] $font(-linespace)


    $e redraw
}

proc editor_lock {w} {
    global $w
    set ed [set ${w}(curr_editor)]
    $ed lock [set ${w}(Lock)]
}

proc editor_save {w} {
    global $w
    set ed [set ${w}(curr_editor)]
    if {[$ed save] != 0} {
	bell
    }
}

proc editor_exit {w} {
    global $w
    set ed [set ${w}(curr_editor)]
    if {[$ed edits_made]} {
	set ret [tk_messageBox \
		     -icon question \
		     -title "Save changes" \
		     -message "Edits have been made. Save changes?" \
		     -type yesnocancel \
		     -parent $w]
	if {$ret == "cancel"} {
	    return
	} elseif {$ret == "yes"} {
	    if {[$ed save] != 0} {
		bell
		return
	    }
	}
    }

    destroy $w
}

proc display_pos_set {ed pos} {
    global $ed
    set ${ed}(displayPos) $pos
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
    frame $f -bd 0

    set w $f.pane
    panedwindow $w -orient horiz -bd 2 -relief sunken

    frame $w.name -bd 0 -highlightthickness 0
    frame $w.seq -bd 0 -highlightthickness 0

    $w add $w.name $w.seq

    # Seqs panel
    set ed $w.seq.sheet
    set ed [editor $ed \
		-width 80 \
		-height 16 \
		-xscrollcommand "display_pos_set $ed \[$ed xview\];
                                 $w.seq.x set" \
		-yscrollcommand "$w.seq.y set" \
		-qual_fg     [keylget gap5_defs CONTIG_EDITOR.QUAL_IGNORE] \
		-diff1_bg    [keylget gap5_defs CONTIG_EDITOR.DIFF1_BG] \
		-diff2_bg    [keylget gap5_defs CONTIG_EDITOR.DIFF2_BG] \
		-diff1_fg    [keylget gap5_defs CONTIG_EDITOR.DIFF1_FG] \
		-diff2_fg    [keylget gap5_defs CONTIG_EDITOR.DIFF2_FG] \
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
	        -consensus_at_top $cattop]
    set opt(curr_editor) $ed

    # X and y scrollbars
    scrollbar $w.seq.x -orient horiz -repeatinterval 30
    scrollbar $w.seq.y -orient vert 

    # Names panel
    set edname $w.name.sheet
    ednames $edname \
	-width 15 \
	-height 16 \
	-xscrollcommand "$w.name.x set" \
	-bd 0

    scrollbar $w.name.x -orient horiz

    entry $w.name.pos \
	-textvariable ${ed}(displayPos) \
	-font {Helvetica -12} \
	-width 15
    bind $w.name.pos <Return> "editor_goto $ed $w.seq.sheet"

    # The jog control for scrolling the editor
    set posh [font metrics [$w.name.pos cget -font] -linespace]
    incr posh -2
    jog $w.seq.jog \
	-orient horiz \
	-command "jog_editor $w.seq.sheet" \
	-repeatinterval 50 -width $posh

    # Pack it all in the paned window
    grid rowconfigure $w.name $textrow -weight 1
    grid columnconfigure $w.name 0 -weight 1
    grid $w.name.sheet -row $textrow   -sticky nsew
    grid $w.name.x     -row $scrollrow -sticky nsew
    grid $w.name.pos   -row $jogrow    -sticky nsew

    $w.name configure -width [expr {16*[font measure sheet_font A]}]

    grid rowconfigure $w.seq $textrow -weight 1
    grid columnconfigure $w.seq 0 -weight 1
    grid $w.seq.sheet $w.seq.y -row $textrow   -sticky nsew
    grid $w.seq.x              -row $scrollrow -sticky nsew
    grid $w.seq.jog            -row $jogrow -sticky nsew

    focus $w.seq.sheet

    # Initialise with an IO and link name/seq panel together
    puts "$ed init $opt(-io) $opt(contig$ind) $opt(-reading$ind) $opt(-pos$ind) $w.name.sheet"
#    $ed init $opt(-io) $opt(contig$ind) $opt(-reading$ind) $opt(-pos$ind) $w.name.sheet
    $ed init $opt(-io) $opt(contig$ind) 0 $opt(-pos$ind) $w.name.sheet
    global $ed $edname
    set ${ed}(parent) $w
    set ${ed}(top) $top
    set ${ed}(Undo) ""
    set ${ed}(Redo) ""
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

    return $ed
}

# The "differences" bar that separates a pair of join editors
proc diff_pane {w} {
    panedwindow $w -orient horiz -bd 2 -relief sunken

    frame $w.name -bd 0 -highlightthickness 0
    label $w.name.diff -text "Differences"
    pack $w.name.diff

    frame $w.seq -bd 0 -highlightthickness 0
    sheet $w.seq.sheet
    pack $w.seq.sheet -fill both -expand 1

    $w add $w.name $w.seq
    
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

# A coordinate jump
proc editor_goto {ed w} {
    upvar \#0 $ed eopt

    set pos $eopt(displayPos)

    if {[regexp {^[-+]?[0-9]+$} $pos] == 1} {
	$w xview $pos
    }
}

# Callback for cutoffs button
proc editor_cutoffs {w} {
    upvar \#0 $w opt

    set ed $opt(curr_editor)
    $ed configure -display_cutoffs $opt(Cutoffs)
    $ed redraw
}

# Callback for quality button
proc editor_quality {w} {
    upvar \#0 $w opt

    set ed $opt(curr_editor)
    $ed configure -display_quality $opt(Quality)
    $ed configure -display_mapping_quality $opt(Quality)
    $ed redraw
}

proc editor_disagreements {w} {
    upvar \#0 $w opt

    set ed $opt(curr_editor)
    if {$opt(Disagreements)} {
	$ed configure -display_differences $opt(DisagreeMode)
    } else {
	$ed configure -display_differences 0
    }
    $ed configure -differences_case_sensitive $opt(DisagreeCase)
    $ed redraw
}

proc set_differences_quality_callback {ed val} {
    $ed configure -display_differences_quality $val
    $ed redraw
}

proc set_differences_quality {w} {
    upvar \#0 $w opt
    global gap5_defs

    set ed $opt(curr_editor)

    set t $ed.qual_diff
    if {[xtoplevel $t -resizable 0] == ""} {return}
    wm title $t "Set differences quality"

    set start [lindex [$ed configure -display_differences_quality] 3]

    scalebox $t.qual \
	-title "Quality" \
	-orient horizontal \
	-from 0 \
	-to 100 \
	-width 5 \
	-default $start \
	-type CheckInt \
	-command "set_differences_quality_callback $ed"
    $t.qual.scale configure -length 150

    okcancelhelp $t.ok \
	-ok_command "keylset gap5_defs CONTIG_EDITOR.DISAGREE_QUAL \[scalebox_get $t.qual\]; destroy $t" \
	-cancel_command "set_differences_quality_callback $ed $start; destroy $t" \
        -help_command "show_help gap4 {Editor-Differences Quality}"

    pack $t.qual $t.ok -side top -fill both
}

proc editor_disagreements_quality {w} {
    $ed configure  $opt(DisagreeQuality)
    $ed redraw
}


#-----------------------------------------------------------------------------
# Undo support
proc store_undo {w cmdu cmdr} {
    global $w

    lappend ${w}(Undo) [list $cmdu $cmdr]
    set ${w}(Redo) ""
    set top [set ${w}(top)]
    $top.toolbar.undo configure -state normal
    $top.toolbar.redo configure -state disabled
}

proc editor_undo {top} {
    upvar \#0 $top opt
    set w $opt(curr_editor)
    global $w

    foreach {cmdu cmdr} [lindex [set ${w}(Undo)] end] break
    eval $cmdu
    lappend ${w}(Redo) [list $cmdu $cmdr]
    set ${w}(Undo) [lrange [set ${w}(Undo)] 0 end-1]

    if {[llength [set ${w}(Undo)]] == 0} {
	$top.toolbar.undo configure -state disabled
    }
    $top.toolbar.redo configure -state normal

    $w redraw
}

proc editor_redo {top} {
    upvar \#0 $top opt
    set w $opt(curr_editor)
    global $w

    foreach {cmdu cmdr} [lindex [set ${w}(Redo)] end] break
    eval $cmdr
    lappend ${w}(Undo) [list $cmdu $cmdr]
    set ${w}(Redo) [lrange [set ${w}(Redo)] 0 end-1]

    if {[llength [set ${w}(Redo)]] == 0} {
	$top.toolbar.redo configure -state disabled
    }
    $top.toolbar.undo configure -state normal

    $w redraw
} 


#-----------------------------------------------------------------------------
# Basic sequence editing function
proc editor_edit_base {w call where} {
    upvar $w opt

    #set io $opt(-io)
    set io [$w io]

    if {$where == ""} {
	bell
	return
    }

    foreach {type rec pos} $where break;
    if {$type == 18} {
	set seq [$io get_sequence $rec]
	foreach {old_call old_conf} [$seq get_base $pos] break
	$seq replace_base $pos $call 30
	$w cursor_right
	store_undo $w \
	    "$seq replace_base $pos $old_call $old_conf; $w cursor_left" \
	    "$seq replace_base $pos $call 30; $w cursor_right"
    }

    $w redraw
}

proc editor_insert_gap {w where} {
    upvar $w opt

    set io [$w io]

    if {$where == ""} {
	bell
	return
    }

    foreach {type rec pos} $where break;
    if {$type == 18} {
	set seq [$io get_sequence $rec]
	$seq insert_base $pos * 20

	store_undo $w \
	    "$seq delete_base $pos; $w cursor_left" \
	    "$seq insert_base $pos * 20; $w cursor_right"
    } else {
	set contig [$io get_contig $rec]

	$contig insert_base $pos

	store_undo $w \
	    "$contig delete_base $pos; $w cursor_left" \
	    "$contig insert_base $pos; $w cursor_right"
    }
    $w cursor_right
    
    $w redraw
}

proc editor_move_seq {w where direction} {
    upvar $w opt

    set io [$w io]

    if {$where == ""} {
	bell
	return
    }

    foreach {type rec pos} $where break;
    if {$type != 18} {
	# sequences only
	bell
	return
    }

    set seq  [$io get_sequence $rec]
    set cnum [$seq get_contig]
    set pos  [$seq get_position]

    puts "Moving $rec in $pos@$cnum by $direction"
    set c [$io get_contig $cnum]
    $c remove_sequence $rec
    incr pos $direction
    $c add_sequence $rec $pos

    $w redraw
}

# Updates the editor status line for editor $w.
# x and y are optional, but if set specify the location to display (eg
# for mouse motion events).
proc update_brief {w {x {}} {y {}}} {
    global gap5_defs
    global $w

    if {$x != "" && $y != ""} {
	foreach {type rec pos} [$w get_number $x $y] break
    } else {
	foreach {type rec pos} [$w get_number] break
    }

    if {![info exists type]} {
	return
    }

    if {$type == 18} {
	set msg [$w get_seq_status $type $rec $pos \
		     [keylget gap5_defs BASE_BRIEF_FORMAT1]]
    } else {
	set msg [$w get_seq_status $type $rec $pos \
		     [keylget gap5_defs BASE_BRIEF_FORMAT2]]
    }

    set w [set ${w}(top)]
    global $w
    set ${w}(Status) $msg
}

#-----------------------------------------------------------------------------
# Trace display
proc show_trace {w loc} {
    set io [$w io]
    
    if {$loc == ""} return
    foreach {type rec pos} $loc break;

    puts $rec/$pos
    if {$type == 18} {
	set s [$io get_seq $rec]
	set name [$s get_name]
	set t [trace_add $w.traces $name $w $name]
	after 100 "$t xview 0"

	global $w
	lappend ${w}(Traces) $t
    }

    parray $w
}

#-----------------------------------------------------------------------------
# Generic bindings
bind Editor <Any-Enter> {
    focus %W
}

bind EdNames <Any-Motion> {update_brief [set %W(ed)] @%x @%y}

bind Editor <Any-Motion> {
    global gap5_defs
    global %W

    foreach {type rec pos} [%W get_number @%x @%y] break
    if {![info exists type]} {
	return
    }

    if {$type == 18} {
	set msg [%W get_seq_status $type $rec $pos \
		     [keylget gap5_defs BASE_BRIEF_FORMAT1]]
    } else {
	set msg [%W get_seq_status $type $rec $pos \
		     [keylget gap5_defs BASE_BRIEF_FORMAT2]]
    }

    set w [set %W(top)]
    global $w
    set ${w}(Status) $msg
}

# Jump to read-pair
bind EdNames <3> {
    global %W

    set ed [set %W(ed)]
    foreach {type rec pos} [$ed get_number @%x @%y] break
    if {![info exists type]} {
	return
    }

    if {$type == 18} {
	set other_end [$ed get_template_seqs $rec]
	puts rec=$rec,other_end=$other_end
	if {[llength $other_end] != 1} return
	set s [[$ed io] get_seq $other_end]
	$ed set_cursor 18 $other_end 1
    }
}

bind Editor <1> {
    focus %W
    set d [%W get_number @%x @%y]
    if {$d == ""} {
	return
    } else {
	eval %W set_cursor $d
    }
}

bind Editor <Key-Left>		{%W cursor_left; update_brief %W}
bind Editor <Control-Key-b>	{%W cursor_left; update_brief %W}

bind Editor <Key-Right>		{%W cursor_right; update_brief %W}
bind Editor <Control-Key-f>	{%W cursor_right; update_brief %W}

bind Editor <Key-Up>		{%W cursor_up;    update_brief %W}
bind Editor <Control-Key-p>	{%W cursor_up;    update_brief %W}

bind Editor <Key-Down>		{%W cursor_down;  update_brief %W}
bind Editor <Control-Key-n>	{%W cursor_down;  update_brief %W}

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

# Editing commanda
bind Editor <Key-a> {editor_edit_base %W a [%W get_number]}
bind Editor <Key-c> {editor_edit_base %W c [%W get_number]}
bind Editor <Key-g> {editor_edit_base %W g [%W get_number]}
bind Editor <Key-t> {editor_edit_base %W t [%W get_number]}
bind Editor <Key-i> {editor_insert_gap %W [%W get_number]}

bind Editor <Control-Key-Left>  {editor_move_seq %W [%W get_number] -1}
bind Editor <Control-Key-Right> {editor_move_seq %W [%W get_number]  1}

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