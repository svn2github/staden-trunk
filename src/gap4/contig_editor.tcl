#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#
# This file creates the contig editor window and adds on the appropriate
# generic bindings (which are written in C).
# 

global gap_defs
set f [keylget gap_defs CONTIG_EDITOR.FONT]
option add *Sheet.font $f
option add *Editor.font $f
option add *EdNames.font $f
option add *diff_label.l.font $f
unset f

if {![info exists .cedit.SE_fast_delete_anno]} {
    set .cedit.SE_fast_delete_anno \
	[keylget gap_defs CONTIG_EDITOR.SE_FAST_DELETE_ANNOTATION]
}

# Remember last location the user scrolled the names display to
set ednames_scrollpos [keylget gap_defs CONTIG_EDITOR.NAMES_XPOS]

#-----------------------------------------------------------------------------
# User dialogue for starting up the editors
#-----------------------------------------------------------------------------
proc EditContig2 {io t id} {
    if {[set reading [contig_id_gel $id]] == ""} {
	return
    }

    destroy $t
    edit_contig -io $io -contig $reading -reading $reading
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

    destroy $t
    join_contig -io $io \
	-contig1 $read1 -reading1 $read1 -pos1 1 \
	-contig2 $read2 -reading2 $read2 -pos2 1
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

#
# Creates a contig editor. If 'join' is true then this editor is assumed to be
# part of a join editor (hence subtle changes such as the addition of a lock
# button occur).
#
proc create_editor {w edname join reveal ccut qcut dbptr {sets {}}} {
    set ed     $w.$edname
    global $ed.Repeater $ed.Reveal $ed.Tlist gap_defs db_namelen
    global licence
    global tk_utils_defs
    global $w.LREG
    global $w.RREG
    global tcl_version

    if {![winfo exists $w]} {
	xtoplevel $w
	fix_maxsize $w
    }
    wm resizable $w 1 0

    catch {set $w.Grab [grab current]}
    frame $w.grab
    while {[catch {grab $w.grab} var]} { exec sleep 1 }

    wm protocol $w WM_DELETE_WINDOW {;}

    wm geometry $w {}

    # Create frames - a complete editor complex
    # An editor is split into a top (buttons), middle (buttons + scroll) and
    # bottom (names + sequences).
    # See editor_quit_internal if you change this layout.
    set top	$ed.top
    set middle	$ed.middle
    set bottom	$ed.bottom
    set status  $w.status
    set tabs    $ed.tabs

    set editor	$ed.seqs
    set names	$ed.names
    set fine 	$ed.fine
    set scrollx	$ed.scrollx
    set scrolly	$ed.scrolly
    set namesx	$ed.namesx
    set cutoffs	$ed.cutoffs
    set buttons	$ed.buttons

    frame $ed     -bd 0
    frame $top    -bd 0
    frame $middle -bd 0
    frame $bottom -bd 0

    # The very top of the editor is now a tabset to give the feel of a tabbed
    # notebook. We only have one editor, but the tabs switch between
    # sequence sets
    global $ed.Sets
    set $ed.Sets $sets
    if {[llength $sets]} {
	set $ed.Tlist [editor_tabs_create $tabs $editor $sets]
    } else {
	set $ed.Tlist {}
    }

    # TOP: On the left we have the confidence and quality cutoff adjusters.
    # On the right are a series of buttons
    frame $cutoffs -bd 2 -relief raised
    frame $buttons -bd 2 -relief raised
    set font [keylget gap_defs CONTIG_EDITOR.FONT_MB]

    option add *buttons*font $font
    option add *cutoffs*font $font

    frame $cutoffs.c -bd 2 -relief raised
    frame $cutoffs.q -bd 2 -relief raised
    if {$tcl_version == "8.3"} {
	repeater $cutoffs.c.down "c_down $editor $cutoffs.c.val" -text "<" -padx 0
	repeater $cutoffs.c.up   "c_up   $editor $cutoffs.c.val" -text ">" -padx 0
	repeater $cutoffs.q.down "q_down $editor $cutoffs.q.val" -text "<" -padx 0
	repeater $cutoffs.q.up   "q_up   $editor $cutoffs.q.val" -text ">" -padx 0
	label $cutoffs.c.label -text C:
	label $cutoffs.q.label -text Q:
	entry $cutoffs.c.val -width 3
	entry $cutoffs.q.val -width 3
	eval $cutoffs.c.val insert 0 $ccut
	eval $cutoffs.q.val insert 0 $qcut

	pack $cutoffs.c.down -fill x -side left -expand 1
	pack $cutoffs.c.label $cutoffs.c.val -side left
	pack $cutoffs.c.up -fill x -side left -expand 1
	pack $cutoffs.q.down -fill x -side left -expand 1
	pack $cutoffs.q.label $cutoffs.q.val -side left
	pack $cutoffs.q.up -fill x -side left -expand 1
    } else {
	# Tk 8.4 added the Spinbox
	label $cutoffs.c.label -text Cons
	label $cutoffs.q.label -text Qual
	spinbox $cutoffs.c.val \
		-command "c_set $editor $cutoffs.c.val" \
		-from -1 \
		-to 100 \
		-width 3
	spinbox $cutoffs.q.val \
		-command "q_set $editor $cutoffs.q.val" \
		-from -1 \
		-to 100 \
		-width 3
	$cutoffs.c.val set $ccut
	$cutoffs.q.val set $qcut

	pack $cutoffs.c.label $cutoffs.c.val -side left
	pack $cutoffs.q.label $cutoffs.q.val -side left
    }
    bind $cutoffs.c.val <Key-Return> "c_set $editor $cutoffs.c.val"
    bind $cutoffs.q.val <Key-Return> "q_set $editor $cutoffs.q.val"
    bind $cutoffs.c.val <Any-Leave> "c_set $editor $cutoffs.c.val"
    bind $cutoffs.q.val <Any-Leave> "q_set $editor $cutoffs.q.val"
    bind $cutoffs.c.val <Any-FocusOut> "c_set $editor $cutoffs.c.val"
    bind $cutoffs.q.val <Any-FocusOut> "q_set $editor $cutoffs.q.val"

    if {$join} {
	button $buttons.quit \
	    -command "editor_quit $w $ed $editor $buttons.quit" \
	    -text Join/Quit -padx 2
    } else {
	button $buttons.quit \
	    -command "editor_quit $w $ed $editor $buttons.quit" \
	    -text Quit -padx 2
    }
    if {$licence(type) != "v"} {
        checkbutton $buttons.insert -text Insert -variable $ed.Insert \
	    -command "$editor set_insert" -bd 2 -relief raised -padx 2
        xmenubutton $buttons.superedit -text "Edit Modes >>" \
	    -menu $buttons.superedit.menu.edit_modes -padx 2
    }
    checkbutton $buttons.reveal -text "Cutoffs"  -bd 2 -relief raised\
	-variable $ed.Reveal -command "$editor set_reveal" -padx 2
    set $ed.Reveal $reveal
    if {$licence(type) != "v"} {
        button $buttons.undo -text "Undo" -command "$editor undo" -padx 2
    }
    button $buttons.next -text "Next Search" -padx 2\
	-command "create_search_win $editor.search \"$editor search\" 0"
    if {$join} {
        checkbutton $buttons.lock -text "Lock" -command "$editor join_lock" \
		-variable $w.Lock  -bd 2 -relief raised -padx 2
	$buttons.lock select
        button $buttons.align \
	    -text "Align" \
	    -command "SetBusy; catch {$editor join_align}; ClearBusy" \
	    -padx 2
    }
    xmenubutton $buttons.commands -text "Commands >>" \
	-menu $buttons.commands.menu.commands -padx 2
    xmenubutton $buttons.settings -text "Settings >>" \
	-menu $buttons.settings.menu.settings -padx 2
    xmenubutton $buttons.help -text "Help >>" \
	-menu $buttons.help.menu.help -padx 2

    create_editor_menus $dbptr $join $ed $editor $names \
	$buttons.commands.menu $buttons.settings.menu $buttons.help.menu \
	$buttons.superedit.menu 

    # MIDDLE: Movement control. We have fine movement by the buttons on the
    # left and larger movement by the scrollbar on the right.
    frame $fine
    option add *fine*font $font
    frame $namesx

    repeater $fine.ll "scroll_ll $editor $scrollx" -text <<
    repeater $fine.l  "scroll_l  $editor $scrollx" -text < 
    repeater $fine.r  "scroll_r  $editor $scrollx" -text > 
    repeater $fine.rr "scroll_rr $editor $scrollx" -text >>
    scrollbar $namesx.s -orient horiz -command "scroll_names  $names"
    scrollbar $scrollx -orient horiz -command "scroll_editor $editor"
    global tcl_platform
    if {$tcl_platform(os) == "Darwin"} {
	# FIXME - This works around a bug in Aqua-tk
	scrollbar $scrolly -orient vert -command "$editor yview" -bd 1
    } else {
	scrollbar $scrolly -orient vert -command "$editor yview"
    }
 
    # Compute maximum width and height for editor, based on screen size
    set max [expr {([winfo screenwidth $w] - 42 - \
	            [keylget tk_utils_defs X_BORDER_SIZE]) / \
	    [font measure sheet_font "A"] - \
	    [keylget gap_defs CONTIG_EDITOR.NAMES_WIDTH]}]
    if {[keylget gap_defs CONTIG_EDITOR.SEQ_WIDTH] > $max} {
	keylset gap_defs CONTIG_EDITOR.SEQ_WIDTH $max
    }

    set fh [font metrics sheet_font -linespace]
    set max [expr {([winfo screenheight $w] - 115 - \
	            [keylget tk_utils_defs Y_BORDER_SIZE]) / $fh - 1}]

    if {$max > [keylget gap_defs CONTIG_EDITOR.MAX_HEIGHT]} {
	set max [keylget gap_defs CONTIG_EDITOR.MAX_HEIGHT]
    }

    if {$join} {
	set max_height [expr int(($max-4)/2)]
    } else {
	set max_height $max
    }

    # BOTTOM: Left is the names display, right are the sequences.
    editor $editor \
	-width [keylget gap_defs CONTIG_EDITOR.SEQ_WIDTH] \
        -height 1 -bd 2 -relief raised \
	-xscrollcommand "$scrollx set" \
	-yscrollcommand "$scrolly set" \
	-highlightcommand "$status.l configure -text" \
	-lightcolour [keylget gap_defs CONTIG_EDITOR.CUTOFF_COLOUR] \
	-fg [keylget gap_defs CONTIG_EDITOR.BASE_COLOUR] \
        -qualcolour0 [keylget gap_defs CONTIG_EDITOR.QUAL0_COLOUR] \
        -qualcolour1 [keylget gap_defs CONTIG_EDITOR.QUAL1_COLOUR] \
        -qualcolour2 [keylget gap_defs CONTIG_EDITOR.QUAL2_COLOUR] \
        -qualcolour3 [keylget gap_defs CONTIG_EDITOR.QUAL3_COLOUR] \
        -qualcolour4 [keylget gap_defs CONTIG_EDITOR.QUAL4_COLOUR] \
        -qualcolour5 [keylget gap_defs CONTIG_EDITOR.QUAL5_COLOUR] \
        -qualcolour6 [keylget gap_defs CONTIG_EDITOR.QUAL6_COLOUR] \
        -qualcolour7 [keylget gap_defs CONTIG_EDITOR.QUAL7_COLOUR] \
        -qualcolour8 [keylget gap_defs CONTIG_EDITOR.QUAL8_COLOUR] \
        -qualcolour9 [keylget gap_defs CONTIG_EDITOR.QUAL9_COLOUR] \
	-qual_fg     [keylget gap_defs CONTIG_EDITOR.QUAL_IGNORE] \
	-diff_bg     [keylget gap_defs CONTIG_EDITOR.DIFF_COLOUR] \
	-editcolour0 [keylget gap_defs CONTIG_EDITOR.EDIT_DEL_COLOUR] \
	-editcolour1 [keylget gap_defs CONTIG_EDITOR.EDIT_BASE_COLOUR] \
	-editcolour2 [keylget gap_defs CONTIG_EDITOR.EDIT_PAD_COLOUR] \
	-editcolour3 [keylget gap_defs CONTIG_EDITOR.EDIT_CONF_COLOUR] \
	-tmplcolour0 [keylget gap_defs CONTIG_EDITOR.TEMP_DIST_COLOUR] \
	-tmplcolour1 [keylget gap_defs CONTIG_EDITOR.TEMP_STRAND_COLOUR] \
	-tmplcolour2 [keylget gap_defs CONTIG_EDITOR.TEMP_PRIM_COLOUR] \
	-tmplcolour3 [keylget gap_defs CONTIG_EDITOR.TEMP_OTHER_COLOUR] \
	-tmplcolour4 [keylget gap_defs CONTIG_EDITOR.TEMP_END_COLOUR] \
	-tmplcolour5 [keylget gap_defs CONTIG_EDITOR.TEMP_ENDSPAN_COLOUR] \
	-max_height  $max_height \
	-bg [tk::Darken [. cget -bg] 115]

    ednames $names \
        -width [keylget gap_defs CONTIG_EDITOR.NAMES_WIDTH] \
	-height 1 -bd 2 -relief raised \
	-lightcolour [keylget gap_defs CONTIG_EDITOR.CUTOFF_COLOUR] \
	-fg [keylget gap_defs CONTIG_EDITOR.BASE_COLOUR] \
	-xscrollcommand "$namesx.s set" \
	-bg [tk::Darken [. cget -bg] 115]
    
    # STATUS: a brief status line at the bottom to report tag and reading info
    if {![winfo exists $status]} {
        frame $status -bd 2 -relief groove
        label $status.dummy
        label $status.l
    }

    # Packing and placing
    pack $ed -side top -fill both -expand 1
    pack $cutoffs.c $cutoffs.q -side left -expand 1

    pack $buttons.help -side right -fill both
    pack $buttons.quit -side right -fill both
    if {$licence(type) != "v"} {
        pack $buttons.insert $buttons.superedit $buttons.reveal \
	    $buttons.undo $buttons.next -side left -fill both;
    } else {
        pack $buttons.reveal $buttons.next -side left -fill both;
    }
    if "$join" {
        pack $buttons.lock $buttons.align -side left -fill both -padx 1
    }
    pack $buttons.commands $buttons.settings -side left -fill both

    pack $fine.ll $fine.l $fine.r $fine.rr -side left -fill both -expand 1

    pack $namesx.s -side left -fill both -expand 1    

    grid $cutoffs -row 0 -column 0 -sticky ew
    grid $buttons -row 0 -column 1 -sticky ew
    grid $fine    -row 1 -column 0 -sticky ew
    grid $scrollx -row 1 -column 1 -sticky nsew
    grid $namesx  -row 2 -column 0 -sticky ew
    grid $names   -row 3 -column 0 -sticky ew
    grid $editor  -row 2 -rowspan 2 -column 1 -sticky ew
    grid $scrolly -row 2 -rowspan 2 -column 2 -sticky ns
    grid columnconfigure $ed 1 -weight 1
    tkwait visibility $names

    # If the width of the names display is less than the width of the fine
    # control or the cutoffs, then we need to round up names to the next width
    # before proceeding.
    set wid [winfo width $names]
    set fw [font measure [lindex [$names configure -font] 4] 0]
    set fh [font metrics [lindex [$names configure -font] 4] -linespace]
    set maxwid [expr [winfo width $fine] > [winfo width $cutoffs] \
		? [winfo width $cutoffs] : [winfo width $fine]]
    set bdw [lindex [$names configure -bd] 4]
    set nwid [expr ($maxwid - 2 * $bdw + $fw - 1) / $fw]
    $names configure -width $nwid
#    $fine configure -height $fh
    $namesx configure -height $fh
    set wid [winfo width $names]
    grid forget $cutoffs $buttons $namesx $fine $scrollx $names $editor $scrolly

    pack $top $middle -side top -fill x
    pack $bottom -side top -fill both -expand 1

    pack $cutoffs -in $top -side left -fill both
    pack $buttons -in $top -side left -fill both -expand 1
    pack $fine -in $middle -side left -fill both
    frame $middle.pad -width [winfo width $scrolly]
    pack $middle.pad -side right -fill both
    pack $scrollx -in $middle -side left -fill both -expand 1
    frame $bottom.l -bd 0
    pack $bottom.l -in $bottom -side left -fill both
    pack $namesx $names -in $bottom.l -side top -fill x
    pack $editor  -in $bottom -side left -anchor n -fill x -expand 1
    pack $scrolly -in $bottom -side right -fill both

    pack $status -side bottom -fill both -expand 1
    pack $status.dummy -fill both
    place $status.l -relx 0

    # Shuffle things to align correctly - yuk.
    pack propagate $cutoffs 0
    pack propagate $fine 0
    pack propagate $namesx 0
    $cutoffs configure -width $wid
    $fine    configure -width $wid
#    $namesx  configure -width $wid

    # Eeek. When the user resizes the window, we need to reset the geometry to
    # {} to allow the editor to specify its own dimensions (it uses the user X
    # size, and its own Y size), otherwise user requested ones are honoured
    # forever more. If we set wm geometry now then loops ensue.
    # Doing it when idle still causes loops; why?
    # So we hack it by waiting a second. Hopefully this doesn't
    # cause any problems, but it still means it can loop on slow
    # systems.
    #
    # 11th March 2003:
    # Since upgrading to Tk8.4.0 the editor resizing has broken (it always
    # jumps back to its original size). However commenting out this "fix"
    # now fixes the resize problem and it seems that we no longer need this
    # change anyway.
    #
    # 4th May 2004:
    # It's about time this is totally rewritten! However, I put the hack back
    # again. Since the change to disable 'SetGrid' on the sheet widget to fix
    # a bug in the MetaCity wm, this broke twm.
    #
    bind $w <Any-Configure> {
   	if {[winfo toplevel %W] == "%W"} {
   	    after 1000 {if {[winfo exists %W]} {wm geometry %W {}}}
   	}
    }

    SetDefaultTags CONTIG_EDITOR.TAGS $editor
    wm protocol $w WM_DELETE_WINDOW \
	"editor_quit $w $ed $editor $buttons.quit"

    if {[set $w.Grab] != ""} {
	catch {grab [set $w.Grab]}
    } else {
	catch {grab release $w.grab}
    }
    unset $w.Grab
    destroy $w.grab

    focus $editor
    bindtags $editor "$editor allfocus Editor all"
    bindtags $names  "$names allfocus EdNames all"

    bind $editor <Alt-Key-Up> "$tabs prev"
    bind $editor <Alt-Key-Down> "$tabs next"

    # Put a trace on variable NAMED ${editor}(tags) and use this to update
    # tag menus.
    global $editor.Tags
    trace variable $editor.Tags w update_tag_menus

    global ednames_scrollpos
    after idle "catch {$names xview $ednames_scrollpos}"

    return "$editor $names"
}

proc create_editor_tabs_old {tabs sets editor} {
    if {[xtoplevel $tabs -resizable 1] == ""} {return}
    set tabs $tabs.t
    iwidgets::tabset $tabs \
	-tabpos n \
	-command "$editor view_set"

    set gnum 1
    $tabs add -label All
    foreach g $sets {
	$tabs add -label $gnum
	incr gnum
    }

    pack $tabs -side bottom -fill both -expand 1
}

proc editor_tabs_create {t editor {sets NONE}} {
    if {[xtoplevel $t -resizable 1] == ""} {return}
    set l $t.list
    tablelist $t.list \
	-columns {0 "Set" 0 "No. seqs"} \
	-labelcommand tablelist::sortByColumn \
	-exportselection 0 \
	-stretch 0 \
	-selectmode extended \
	-yscrollcommand [list $t.yscroll set]
    $t.list columnconfigure 1 -sortmode integer

    frame $t.buttons -bd 0
    button $t.buttons.cancel \
	-text "Cancel" \
	-command "destroy $t"

    button $t.buttons.help \
	-text "Help" \
	-command "show_help gap4 FIXME"

    pack $t.buttons.cancel $t.buttons.help -side left -expand 1

    # Add a scrollbar    
    scrollbar $t.yscroll -command "$t.list yview"

    grid columnconfigure $t 0 -weight 1
    grid rowconfigure $t 0 -weight 1

    grid $t.list $t.yscroll -sticky nsew
    grid $t.buttons -sticky nsew

    # Bindings
    bind [$t.list bodypath] <<menu>> \
	"puts foo"
    bind [$t.list bodypath] <<select>> \
	"editor_tabs_select $editor %W %x %y"

    # Populate
    if {$sets != "NONE"} {
	editor_tabs_repopulate $t.list $sets
    }

    return $t.list
}

proc editor_tabs_repopulate {w sets} {
    $w selection clear 0 end
    $w delete 0 end

    set gnum 1
    $w insert end [list All 0]
    foreach g $sets {
	$w insert end [list $gnum [llength $g]]
	incr gnum
    }
}

proc editor_tabs_select {editor w x y} {
    foreach {line char} [split [$w index @$x,$y] .] {}
    set grp [lindex [$w get $line.0 $line.end] 0]
    if {$grp == "All"} {set grp 0}
    $editor view_set $grp
}

# Updates the commands menu Edit Tags and Delete Tags cascades so that the
# tags under the current cursor are listed.
proc update_tag_menus {name1 name2 op} {
    global $name1 gap_defs
    regexp {(.*)\.[^.]*\.[^.]*$} $name1 _ e

    set menu $e.buttons.commands.menu.commands.edit_tag
    $menu configure -tearoff 0
    $menu delete 0 end
    set count 0
    foreach tag [set $name1] {
	foreach {ptr type st len} $tag {}
	set end [expr {$st+$len-1}]
	if {$count >= [keylget gap_defs MAX_MENU_ITEMS]} {
	    set count 0
	    $menu add cascade \
		-label "More..." \
		-menu $menu.more
	    set menu [menu $menu.more -tearoff 0]
	}
	$menu add command \
	    -label "$type $st-$end" \
	    -command "$e.seqs edit_anno $ptr"
	incr count
    }

    set menu $e.buttons.commands.menu.commands.delete_tag
    $menu configure -tearoff 0
    $menu delete 0 end
    set count 0
    foreach tag [set $name1] {
	foreach {ptr type st len} $tag {}
	set end [expr {$st+$len-1}]
	if {$count >= [keylget gap_defs MAX_MENU_ITEMS]} {
	    set count 0
	    $menu add cascade \
		-label "More..." \
		-menu $menu.more
	    set menu [menu $menu.more -tearoff 0]
	}
	$menu add command \
	    -label "$type $st-$end" \
	    -command "$e.seqs delete_anno $ptr"
	incr count
    }

    if {[keylget gap_defs CONTIG_EDITOR.TAG_POPUPS]} {
	if {[catch {handle_tag_popups $name1 $e} err]} {
	    global errorInfo; puts $err,$errorInfo
	}
    }
}

catch {font create tooltip_font    -family Helvetica -size -12}
proc tag_popup {w msg} {
    # Work out coordinates
    set x [expr {[winfo pointerx $w]+5}]
    set y [expr {[winfo pointery $w]+15}]
    set padx 5
    set bd 1
    set xdim 0
    set lines 0
    foreach line [split $msg "\n"] {
	set dim [font measure tooltip_font $line]
	if {$dim > $xdim} {
	    set xdim $dim
	}
	incr lines
    }
    incr xdim [expr {2*($padx+$bd)}]
    array set metrics [font metrics tooltip_font]
    set ydim [expr {$lines*$metrics(-linespace) + 2*$bd+2}]

    # Make sure it'llbe on the screen.
    if {[expr {$x+$xdim}] > [winfo screenwidth $w]} {
	set x [expr {[winfo screenwidth $w]-$xdim}]
    }
    if {[expr {$y+$ydim}] > [winfo screenheight $w]} {
	set y [expr {[winfo pointery $w]-$ydim-10}]
    }

    # Create the window itself
    set p $w.tag_tooltip
    if {[winfo exists $p]} {destroy $p}
    toplevel $p
    wm overrideredirect $p 1
    wm geometry $p +$x+$y

    label $p.l \
	-bg lightyellow \
	-justify left \
	-padx $padx \
	-bd $bd -relief solid \
	-font tooltip_font \
	-text $msg
    pack $p.l -fill both -expand 1

    # Bindings to remove it
    bind $p <1> {destroy [winfo toplevel %W]}
}

proc handle_tag_popups {t e} {
    global $t gap_defs

    set msg ""
    set msg2 ""
    foreach tag [set $t] {
	foreach {ptr type st len} $tag {}
	if {$msg != ""} {append msg "\n"}
	set txt [$e.seqs get_anno $ptr]
	append msg "$type:$txt"
	append msg2 "$type\("
	set first 1
	foreach line [split $txt "\n"] {
	    foreach {var val} [split $line =] {break}
	    if {$var == "TABNAME"} {
		if {!$first} {append msg2 " "}
		set first 0
		append msg2 "$val:"
	    }
	    if {$var == "type"} {
		append msg2 "$val"
	    }
	}
	append msg2 "\) "
    }

    if {[keylget gap_defs CONTIG_EDITOR.TAG_POPUPS] == 2} {
	if {$msg == {}} {
	    catch {destroy $e.tag_tooltip}
	} else {
	    tag_popup $e $msg
	}
    }
    if {$msg2 != ""} {
	after idle "[lindex [$e.seqs configure -highlightcommand] 4] [list $msg2]"
    }
}
 
# Adds a 'diff' bar between a couple of editors. Assumes at least one editor
# is already packed in 'w'.
proc create_editor_diff {w dname edname} {
    set e $w.$edname.seqs
    frame $w.$dname -bd 0 -relief raised
    frame $w.$dname.diff_label -bd 2 -relief raised
    button $w.$dname.diff_label.prev -text "<" -command "$e prev_difference"
    label $w.$dname.diff_label.l -text "Differences"
    button $w.$dname.diff_label.next -text ">" -command "$e next_difference"
    sheet $w.$dname.diffs -width 60 -height 1 -bd 2 -relief raised

    pack $w.$dname -side top -fill both
    pack $w.$dname.diff_label -side left -fill both
    pack $w.$dname.diff_label.prev -side left
    pack $w.$dname.diff_label.l -side left -fill both -expand 1
    pack $w.$dname.diff_label.next -side left
    pack $w.$dname.diffs -fill both -expand 1

    pack propagate $w.$dname.diff_label 0
    $w.$dname.diff_label configure -width [winfo width $w.$edname.names]

    return $w.$dname.diffs
}

proc init_editor_states {w e dbptr} {
    global gap_defs consensus_mode licence read_only
    global $w.ShowQuality $w.AminoMode $w.DisplayTraces $w.AutoSave \
	   $w.GroupBy $w.ShowEdits $w.Disagreements $w.CompareStrands \
	   $w.ShowCQuality $w.DisagreeMode $w.ShowUnpadded \
	   $w.ConsensusAlgorithm _$dbptr.StoreUndo $w.DisplayTraces \
	   $w.DiffTraces $w.ReadPairTraces $w.DisplayMiniTraces
    global $w.Status0 $w.Status1 $w.Status2 $w.Status3 $w.Status4 \
	   $w.Status5 $w.Status6 $w.Status7
    global $w.LREG $w.RREG $w.TemplateNames

    # Initialise contig editor state
    foreach [list $w.LREG $w.RREG] [$e get_extents] {}
    if {[set $w.DisplayTraces]} {$e autodisplay_traces}
    $e show_mini_traces [set $w.DisplayMiniTraces]
    if {[set $w.AutoSave] && !$read_only}      {$e auto_save}
    if {[set $w.DiffTraces]}    {$e autodiff_traces}
    if {[set $w.ReadPairTraces]}    {$e read_pair_traces}
    set _$dbptr.StoreUndo [$e store_undo]
    if {[$e number_of_views] == 1} {
	if {[$e join_mode]} {
	    set su [keylget gap_defs JOIN_EDITOR.STORE_UNDO]
	} else {
	    set su [keylget gap_defs CONTIG_EDITOR.STORE_UNDO]
	}
	if {$su != [$e store_undo]} {
	    set _$dbptr.StoreUndo $su
	    $e store_undo $su
	}
    }
    if {[set $w.Disagreements]} {
      $e show_differences [set $w.DisagreeMode]
    }
    set $w.ConsensusAlgorithm $consensus_mode
    $e compare_strands [set $w.CompareStrands]
    $e translation_mode [lindex {1 3} [set $w.AminoMode]]
    $e show_quality [set $w.ShowQuality]
    $e show_consensus_quality [set $w.ShowCQuality]
    $e show_edits [set $w.ShowEdits]
    $e show_template_names [set $w.TemplateNames]
    editor_set_superedit $e $w
    $e set_grouping [set $w.GroupBy]
    $e insert_confidence [keylget gap_defs CONTIG_EDITOR.INSERTION_CONFIDENCE]
    $e replace_confidence [keylget gap_defs CONTIG_EDITOR.REPLACE_CONFIDENCE]
    $e set_trace_lock [keylget gap_defs CONTIG_EDITOR.TRACE_LOCK]
    $e set_unpadded_ruler [set $w.ShowUnpadded]
    $e set_consensus_mode [set $w.ConsensusAlgorithm]
    $e set_disagreement_cutoff [keylget gap_defs CONTIG_EDITOR.DISAGREE_QUAL]
    for {set i 0} {$i <= 6} {incr i} {
        if {[set $w.Status$i]} {$e status add $i}
    }
    set_editor_write_mode $e [$e write_mode]
    if {[$e join_mode]} {
	menu_state_set contig_editor_settings_menu  -4 $w.buttons.settings.menu
	menu_state_set contig_editor_commands_menu  -4 $w.buttons.commands.menu
	menu_state_set contig_editor_help_menu      -4 $w.buttons.help.menu
	if {$licence(type) != "v"} {
	    menu_state_set contig_editor_editmodes_menu \
		-4 $w.buttons.editmodes.menu
	}
    }

    editor_trace_config $e $w
}

#-----------------------------------------------------------------------------
# Internally usable functions.
#-----------------------------------------------------------------------------

# ====superedit============= #
# o Allow insert in read     #
# o Allow del in read        #
# o Allow insert any in cons #
# o Allow del dash in cons   #
# o Allow del any in cons    #
# o Allow reading shift      #
# o Allow transpose any      #
# o Allow uppercase          #
# -------------------------- #
# x Edit by base type        #
# x Edit by confidence       #
# ========================== #
#
# ====settings============== #
# > Status Line		     #
# > Trace Display            #
# -------------------------- #
# o Highlight Disagreements  #
# x   By dots                #
# x   By foreground colour   #
# x   By background colour   #
# o   Case sensitive         #
# -------------------------- #
# o Compare Strands	     #
# o Toggle auto-save	     #
# o 3 Character Amino Acids  #
# o Group readings by temp.  #
#   Set Active Tags	     #
#   Set Output List	     #
#   Set Default Confidence   #
# ========================== #
#
# ====Status Line=========== #
# o Show Strands             #
# o translate frame 1+       #
# o translate frame 2+       #
# o translate frame 3+       #
# o translate frame 1-       #
# o translate frame 2-       #
# o translate frame 3-       #
# -------------------------- #
#   Translate + frames       #
#   Translate - frames       #
#   Translate all frames     #
# -------------------------- #
#   Remove all               #
# ========================== #
#
# ====Trace Display========= #
# o Auto-display Traces      #
# o Auto-diff Traces         #
# -------------------------- #
# x No differencing          #
# x Diff against consensus   #
# x Diff against specific    #
# -------------------------- #
# o Only matching reads      #
# o Ignore selected read     #
# -------------------------- #
# o Show positive differences#
# o Y scale differences      #
# ========================== #

# ====commands============== #
#   Search		     #
# -------------------------- #
#   Create Tag		     #
#   Edit Tag		     #
#   Delete Tag		     #
# -------------------------- #
#   Save Contig		     #
#   Dump Contig to File	     #
#   Save Consensus Trace     #
# -------------------------- #
#   Select Primer	     #
#   Align		     #
#   Shuffle Pads	     #
#   Remove reading	     #
#   Break contig	     #
# ========================== #
#
# ====help================== #
#   Help                     #
#   ....                     #
proc create_editor_menus {dbptr join w e n m1 m2 m3 m4} {
    global gap_defs $w.ShowQuality $w.AminoMode $w.DisplayTraces $w.AutoSave
    global $w.SE_ins_read $w.SE_del_read $w.SE_ins_cons $w.SE_del_dash_cons
    global $w.SE_del_any_cons $w.SE_read_shift $w.SE_trans_any
    global $w.SE_uppercase $w.SE_edit_mode $w.SE_replace_cons
    global $w.TraceDiff $w.TraceConsMatch $w.TraceConsSelect $w.DisagreeMode
    global $w.TraceDiffAlgorithm $w.TraceDiffScale $w.GroupBy
    global $w.ShowEdits $w.Disagreements $w.CompareStrands
    global $w.ShowCQuality $w.DisagreeCase $w.TemplateNames
    global $w.Status0 $w.Status1 $w.Status2 $w.Status3
    global $w.Status4 $w.Status5 $w.Status6 $w.Status7
    global $w.ShowUnpadded $w.DiffTraces $w.ReadPairTraces $w.DisplayMiniTraces
    global licence

    set $w.Disagreements  [keylget gap_defs CONTIG_EDITOR.DISAGREEMENTS]
    set $w.DisagreeMode	  [keylget gap_defs CONTIG_EDITOR.DISAGREE_MODE]
    set $w.DisagreeCase	  [keylget gap_defs CONTIG_EDITOR.DISAGREE_CASE]
    set $w.CompareStrands [keylget gap_defs CONTIG_EDITOR.COMPARE_STRANDS]
    set $w.DisplayTraces  [keylget gap_defs CONTIG_EDITOR.AUTO_DISPLAY_TRACES]
    set $w.DisplayMiniTraces \
	                  [keylget gap_defs CONTIG_EDITOR.DISPLAY_MINI_TRACES]
    set $w.DiffTraces     [keylget gap_defs CONTIG_EDITOR.AUTO_DIFF_TRACES]
    set $w.ReadPairTraces [keylget gap_defs CONTIG_EDITOR.READ_PAIR_TRACES]
    set $w.AutoSave       [keylget gap_defs CONTIG_EDITOR.AUTO_SAVE]
    set $w.Status0	  [keylget gap_defs CONTIG_EDITOR.STATUS_STRAND]
    set $w.Status1	  [keylget gap_defs CONTIG_EDITOR.STATUS_FRAME1P]
    set $w.Status2	  [keylget gap_defs CONTIG_EDITOR.STATUS_FRAME2P]
    set $w.Status3	  [keylget gap_defs CONTIG_EDITOR.STATUS_FRAME3P]
    set $w.Status4	  [keylget gap_defs CONTIG_EDITOR.STATUS_FRAME1M]
    set $w.Status5	  [keylget gap_defs CONTIG_EDITOR.STATUS_FRAME2M]
    set $w.Status6	  [keylget gap_defs CONTIG_EDITOR.STATUS_FRAME3M]
    set $w.Status7	  [keylget gap_defs CONTIG_EDITOR.STATUS_AUTO_TRANSLATE]
    set $w.AminoMode	  [keylget gap_defs CONTIG_EDITOR.AMINO_ACID_MODE]
    set $w.GroupBy  [keylget gap_defs CONTIG_EDITOR.GROUP_BY]
    set $w.ShowQuality	  [keylget gap_defs CONTIG_EDITOR.SHOW_QUALITY]
    set $w.ShowCQuality	  [keylget gap_defs CONTIG_EDITOR.SHOW_CONSENSUS_QUALITY]

    set $w.ShowUnpadded   [keylget gap_defs CONTIG_EDITOR.SHOW_UNPADDED]
    set $w.ShowEdits	  [keylget gap_defs CONTIG_EDITOR.SHOW_EDITS]
    set $w.TemplateNames  [keylget gap_defs CONTIG_EDITOR.TEMPLATE_NAMES]
    set $w.TraceDiff	  [keylget gap_defs CONTIG_EDITOR.TRACE_DIFF]
    set $w.TraceConsMatch [keylget gap_defs CONTIG_EDITOR.TRACE_CONS_MATCH]
    set $w.TraceConsSelect [keylget gap_defs CONTIG_EDITOR.TRACE_CONS_SELECT]
    set $w.TraceDiffAlgorithm \
			  [keylget gap_defs CONTIG_EDITOR.TRACE_DIFF_ALGORITHM]
    set $w.TraceDiffScale [keylget gap_defs CONTIG_EDITOR.TRACE_DIFF_SCALE]

    set $w.SE_ins_read	     [keylget gap_defs CONTIG_EDITOR.SE_INS_ANY]
    set $w.SE_del_read	     [keylget gap_defs CONTIG_EDITOR.SE_DEL_READ]
    set $w.SE_ins_cons	     [keylget gap_defs CONTIG_EDITOR.SE_INS_CONS]
    set $w.SE_del_dash_cons  [keylget gap_defs CONTIG_EDITOR.SE_DEL_DASH_CONS]
    set $w.SE_del_any_cons   [keylget gap_defs CONTIG_EDITOR.SE_DEL_ANY_CONS]
    set $w.SE_replace_cons   [keylget gap_defs CONTIG_EDITOR.SE_REPLACE_CONS]
    set $w.SE_read_shift     [keylget gap_defs CONTIG_EDITOR.SE_READ_SHIFT]
    set $w.SE_trans_any	     [keylget gap_defs CONTIG_EDITOR.SE_TRANS_ANY]
    set $w.SE_uppercase	     [keylget gap_defs CONTIG_EDITOR.SE_UPPERCASE]
    set $w.SE_edit_mode	     [keylget gap_defs CONTIG_EDITOR.SE_EDIT_MODE]

    # Commands menu
    global contig_editor_commands_menu
    menu $m1
    create_menus $contig_editor_commands_menu $m1 [keylget gap_defs MENU_LEVEL]

    # Settings menu
    global contig_editor_settings_menu
    menu $m2
    create_menus $contig_editor_settings_menu $m2 [keylget gap_defs MENU_LEVEL]
    catch {$m2.settings.trace_display configure -disabledforeground blue}

    # Help menu
    global contig_editor_help_menu
    menu $m3
    create_menus $contig_editor_help_menu $m3 [keylget gap_defs MENU_LEVEL]

    if {$licence(type) != "v"} {
        # Edit Modes menu
        global contig_editor_editmodes_menu
        menu $m4
        create_menus $contig_editor_editmodes_menu $m4 [keylget gap_defs MENU_LEVEL]
    }

    # Load the tag macros
    tag_macro_load
}

proc editor_to_ed {editor} {
    regsub "\.seqs" $editor "" ed
    return $ed
}

proc edname_to_editor {edname} {
    regsub "names$" $edname "seqs" ed
    return $ed
}

proc editor_to_edname {ed} {
    regsub "seqs$" $ed "names" edname
    return $edname
}

proc popup_editor_menu {w x y} {
    regsub "\.seqs" $w ".buttons.commands.menu.commands" m

    tk_popup $m [expr $x+[winfo rootx $w]] [expr $y+[winfo rooty $w]]
}

proc destroy_editor_menu {m} {
    if {[winfo exists $m]} {destroy $m}
}

proc editor_quit {top ed e object} {
    global $ed.Repeater gap_defs read_only

    if {$read_only} {
	$e quit
	return
    }

    # Object details where the dialog will appear.
    set x [winfo rootx $object]
    set y [winfo rooty $object]

    set list [$e get_hidden_reads]
    set list [eval get_read_names -io [$e io] $list]
    ListCreate2 disassemble $list SEQID

    if {[$e join_mode]} {
	foreach {perc tgood tbad} [$e join_percentage] break
	set nt [expr {$tgood+$tbad}]
	set message "$nt spanning template[expr {$nt==1?{}:{s}}]\n"
	if {$nt} {
	    append message "of which $tgood [expr {$tgood==1?{is}:{are}}] good\n"
	    append message "and $tbad [expr {$tbad==1?{is}:{are}}] bad.\n\n"
	} else {
	    append message "\n"
	}
	    
	if {$perc == -1} {
	    set ret [tk_messageBox \
		    -icon warning \
		    -message "${message}Contigs do not overlap." \
		    -title "Quit Editor" \
		    -type okcancel \
		    -parent $top]
        } else {
	    set ret [tk_messageBox \
		    -icon question \
		    -message "[format "${message}Percentage Mismatch: %5.2f%%\nMake join?" $perc]" \
		    -title "Quit Editor" \
		    -type yesnocancel \
		    -parent $top]
       	}
	if {$ret == "cancel"} {
	    # Cancel
	    return
	} elseif {$ret == "yes"} {
	    set mperc [keylget gap_defs DIFF_WARNING]
	    if {$perc >= $mperc} {
		bell
		set ret [tk_messageBox \
			-icon warning \
			-message "WARNING! Percentage mismatch is above $mperc%. Cancel join?"\
			-title "Quit Editor" \
			-type yesno \
			-parent $top]

		if {$ret == "yes"} {
		    # Cancel
		    return
		}
	    }
	    # join
	    toplevel .grab_me
	    wm withdraw .grab_me
	    while {[catch {grab .grab_me} var]} { exec sleep 1 }
	    update
	    $e join
	}
    } elseif {[$e edits_made]} {
	set ret [tk_messageBox \
		-icon question \
		-title "Quit Editor" \
		-message "Save changes" \
		-type yesnocancel \
		-parent $top]
	update
        if {$ret == "cancel"} {
	    # Cancel
	    return
        } elseif {$ret == "yes"} {
	    # Save data
	    $e save
       	}
    }

    # Disassemble if needed.
    if {$list != ""} {
	set ret [tk_messageBox \
		-icon question \
		-title "Disassemble readings" \
		-message "Reading(s) have been selected for disassembly. Do this now?" \
		-type yesnocancel \
		-parent $top]
	if {$ret == "cancel"} {
	    return
	} elseif {$ret == "yes"} {
	    update
	    editor_disassembly $list [$e io] $top
	    return
	}
    }

#    destroy $top
	$e quit

	catch {grab release .grab_me}
	catch {destroy .grab_me}
}

proc set_differences_mode {e w type} {
    global $w.DisagreeMode $w.DisagreeCase $w.Disagreements

    set mode [expr {[set $w.DisagreeMode]+4*[set $w.DisagreeCase]}]

    if {$type == "toggle"} {
	if {[$e show_differences] != 0} {
	    $e show_differences 0
	} else {
	    $e show_differences $mode
	}
    } else {
	$e show_differences $mode
	set $w.Disagreements 1
    }
}

proc set_differences_quality {e} {
    set t $e.qual_diff
    if {[xtoplevel $t -resizable 0] == ""} {return}
    wm title $t "Set differences quality"

    set start [$e set_disagreement_cutoff];

    scalebox $t.qual \
	-title "Quality" \
	-orient horizontal \
	-from 0 \
	-to 100 \
	-width 5 \
	-default [$e set_disagreement_cutoff] \
	-type CheckInt \
	-command "$e set_disagreement_cutoff"
    $t.qual.scale configure -length 150

    okcancelhelp $t.ok \
	-ok_command "keylset gap_defs CONTIG_EDITOR.DISAGREE_QUAL \[$e set_disagreement_cutoff\]; destroy $t" \
	-cancel_command "$e set_disagreement_cutoff $start; destroy $t" \
        -help_command "show_help gap4 {Editor-Differences QUality}"

    pack $t.qual $t.ok -side top -fill both
}

proc editor_disassembly {list io top} {
    DisEditorReadingsDialog $io $list $top.dis
}

proc editor_break_contig {w e} {
    set read [$e get_read_number]
    set io [$e io]

    global read_only
    if {$read_only} {
	bell
	return
    }

    if {$read < 1} {
        verror ERR_WARN break_contig \
		"Editing cursor must be placed on a reading"
	bell
	return
    }

    set r [io_read_reading $io $read]
    if {[keylget r left] == 0} {
        verror ERR_WARN break_contig \
		"No 'left part' of contig to break off. Please move the editing cursor"
	bell
	return
    }

    set ret [tk_messageBox \
	    -icon question \
	    -title "Break contig" \
	    -message "This will save changes and break the contig into two, with reading $read starting the second contig. Continue?" \
	    -type yesno \
	    -parent $w]
    if {$ret == "no"} {
	# No
	return
    }
    $e save

    if {![quit_displays $io "break contig"]} {
	# Someone's too busy to shutdown?
	ClearBusy
	return
    }

    break_contig -io $io -readings [r_name $io $read]
    ContigInitReg $io

    set left [keylget r left]
    edit_contig -io $io \
	-contig "#$left" \
	-reading "#$left"
    edit_contig -io $io \
	-contig "#$read" \
	-reading "#$read"
}

#
# We've been forced to shutdown for some reasons by the internals of the
# editor C code. We know the 'editor' pathname, but little else. We rely on
# knowing how to construct the top level again for this... yuk
#
# Also note that we shut down the trace display first, as doing so updates
# the names display. Due to organisation of the code this generates errors
# if we let the destroy command destroy things in it's own order. (It's easier
# to fix this than add all the checks.)
#
proc editor_quit_internal {w} {
    if {[winfo exists $w.traces]} {destroy $w.traces}
    regsub {\.[^.]*.\.seqs} $w "" top
    update
    destroy $top
}

# Scrollbar callbacks
proc scroll_editor {seqs args} {
    eval $seqs xview $args
}

proc scroll_names {names args} {
    global ednames_scrollpos
    eval $names xview $args
    set ednames_scrollpos [$names xview]
}

proc scroll_ll {seqs sbar} {
    scroll_editor $seqs scroll -40 units
}

proc scroll_l {seqs sbar} {
    scroll_editor $seqs scroll -1 unit
}

proc scroll_r {seqs sbar} {
    scroll_editor $seqs scroll 1 unit
}

proc scroll_rr {seqs sbar} {
    scroll_editor $seqs scroll 40 units
}


# Consenus and quality cutoff percentages.
proc c_down {ed entry} {
    eval set c \"[$entry get]\"

    if ![regexp {^[+-]?[0-9]+$} "$c"] {
	set c 0
    }

    if {$c > 100} {
        $entry delete 0 end
        $entry insert 0 100
        eval $ed set_ccutoff 100
    } else {
	incr c -1
	if {$c < 1} {set c 1}
        eval $ed set_ccutoff $c
        $entry delete 0 end
        $entry insert 0 $c
    }
}

proc c_up {ed entry} {
    eval set c \"[$entry get]\"

    if ![regexp {^[+-]?[0-9]+$} "$c"] {
	set c 0
    }

    if {$c < 1} {
        $entry delete 0 end
        $entry insert 0 1
        eval $ed set_ccutoff 1
    } else {
	incr c
	if {$c > 100} {set c 100}
        eval $ed set_ccutoff $c
        $entry delete 0 end
        $entry insert 0 $c
    }
}

proc c_set {ed entry} {
    eval set c \"[$entry get]\"

    if ![regexp {^[+-]?[0-9]+$} "$c"] {
	set c 0
    }

    if {$c < 1} {set c 1}
    if {$c > 100} {set c 100}

    eval $ed set_ccutoff $c
    $entry delete 0 end
    $entry insert 0 $c
}

proc q_down {ed entry} {
    eval set c \"[$entry get]\"

    if ![regexp {^[+-]?[0-9]+$} "$c"] {
	set c 0
    }

    if {$c > 100} {
        $entry delete 0 end
        $entry insert 0 100
        eval $ed set_qcutoff 100
    } else {
	incr c -1
	if {$c < -1} {set c -1}
        eval $ed set_qcutoff $c
        $entry delete 0 end
        $entry insert 0 $c
    }
}

proc q_up {ed entry} {
    eval set c \"[$entry get]\"

    if ![regexp {^[+-]?[0-9]+$} "$c"] {
	set c 0
    }

    if {$c < -1} {
        $entry delete 0 end
        $entry insert 0 -1
        eval $ed set_qcutoff -1
    } else {
	incr c
	if {$c > 100} {set c 100}
        eval $ed set_qcutoff $c
        $entry delete 0 end
        $entry insert 0 $c
    }
}

proc q_set {ed entry} {
    eval set c \"[$entry get]\"

    if ![regexp {^[+-]?[0-9]+$} "$c"] {
	set c 0
    }

    if {$c < -1} {set c -1}
    if {$c > 100} {set c 100}

    eval $ed set_qcutoff $c
    $entry delete 0 end
    $entry insert 0 $c
}

#
# Sets the default list for the middle button to output to
#
proc editor_setlist {t n} {
    global gap_defs $n.List

    if {[xtoplevel $t -resizable 0] == ""} {return}
    wm title $t "Set output list"

    getLname $t.list "List name" create

    okcancelhelp $t.ok \
       -ok_command "editor_setlist2 $t $n.List \[entrybox_get $t.list.entry\]"\
        -cancel_command "destroy $t" \
	-help_command "show_help gap4 {Editor-Output List}" \
        -bd 2 -relief groove

    pack $t.list $t.ok -side top -fill x
}

proc editor_setlist2 {t l_var l_name} {
    global $l_var

    set $l_var $l_name
    if {![ListExists2 $l_name]} {
        ListCreate2 $l_name {} SEQID
    }
    ListEdit $l_name
    destroy $t
}

# Adds to the active list the reading under a specific x,y coordinate on
# the ednames (n) window
proc editor_addlist {n x y} {
    global $n.List

    set highlighted [$n highlight -1 $y]

    # Get reading name
    if {[set num [$n get_number $x $y]] != "" && "$num" != 0} {
        set name [$n get_name $num]
    } else {
	bell
	return
    }


    # Selection to PRIMARY STRING
    copy_name $n $x $y

    # Add name to list
    if {[info exist $n.List]} {
	if {$highlighted} {
	    editor_addlist_read $n [lindex $name 1]
	} else {
	    editor_dellist_read $n [lindex $name 1]
	}
    }
}

# Adds to the active list one or more sequence names
proc editor_addlist_read {n names} {
    global $n.List
    UpdateReadingDisplays_disable
    if {[info exist $n.List]} {
	ListAppend [set $n.List] $names
	foreach name $names {
	    $n highlight 1 =$name
	}
    } else {
	foreach name $names {
	    $n highlight 1 =$name
	}
    }
    UpdateReadingDisplays_enable
    UpdateReadingDisplays
}

proc editor_dellist_read {n names} {
    global $n.List
    UpdateReadingDisplays_disable
    if {[info exist $n.List]} {
	ListSubtract [set $n.List] $names
	foreach name $names {
	    $n highlight 0 =$name
	}
    } else {
	foreach name $names {
	    $n highlight 0 =$name
	}
    }
    UpdateReadingDisplays_enable
    UpdateReadingDisplays
}

# As editor_addlist_read, but it takes a set of {name contig} pairs as the
# items to add to the list. We just produce a new list consisting of the names
# and use editor_addlist_read
proc editor_addlist_template {n tseqs} {
    set names ""
    foreach {name cnum pos} $tseqs {
	lappend names $name
    }
    editor_addlist_read $n $names
}

proc editor_dellist_template {n tseqs} {
    set names ""
    foreach {name cnum pos} $tseqs {
	lappend names $name
    }
    editor_dellist_read $n $names
}

proc editor_clearlist {n} {
    global $n.List
    if {[info exist $n.List]} {
	catch {
	    foreach name [ListGet [set $n.List]] {
		$n highlight 0 =$name
	    }
	}
	ListClear [set $n.List]
    } else {
	ListClear readings
    }
}

proc copy_name {n x y} {
    # Get reading name
    if {[set num [$n get_number $x $y]] != "" && "$num" != 0} {
        set name [$n get_name $num]
    } else {
	return
    }

    # Selection to PRIMARY STRING
    catch {rename editor_selection_handler {}}
    # Ignore parameters for now. Assume reading name length <= maxbytes.
    proc editor_selection_handler {offset maxbytes} \
	[list return [lindex $name 1]]
    selection own $n
    selection handle $n editor_selection_handler

    # For Windows...
    clipboard clear
    clipboard append $name
}

#
# Sets the default confidence for new and replaced bases
#
proc set_default_confidence {e t} {
    if {[xtoplevel $t -resizable 0] == ""} {return}
    wm title $t "Set default confidences"

    scalebox $t.insert \
	-title "Confidence for inserted bases"\
	-orient horizontal \
	-from -1 \
	-to 100 \
	-width 5 \
	-default [$e insert_confidence]\
	-type CheckInt

    scalebox $t.replace \
	-title "Confidence for replaced bases"\
	-orient horizontal \
	-from -1 \
	-to 100 \
	-width 5 \
	-default [$e replace_confidence]\
	-type CheckInt

    okcancelhelp $t.ok \
        -ok_command "$e replace_confidence \[scalebox_get $t.replace\];
		     $e insert_confidence \[scalebox_get $t.insert\];
		     destroy $t" \
	-help_command "show_help gap4 {Editor-Default Confidence}" \
        -cancel_command "destroy $t" \
        -bd 2 -relief groove

    pack $t.replace $t.insert $t.ok -side top -fill both
}

proc ed_list_confidence {e t} {
    set w [editor_to_ed $e]

    global gap_defs $w.ConsensusAlgorithm
    global $w.LREG
    global $w.RREG

    if {[set $w.ConsensusAlgorithm] != 2} {
	tk_messageBox \
		-icon info \
		-title "Error" \
		-message "Error rates can only be computed when the\
		      `confidence' consensus algorithm is in use." \
		-type ok \
		-parent $w
	return
    }

    if {[xtoplevel $t -resizable 0] == ""} {return}
    wm title $t "List confidence"

    set end_value [set $w.RREG]
    set start_value [set $w.LREG]

    set extents [$e get_extents]
    set left [lindex $extents 0]
    set right [lindex $extents 1]

    scalebox $t.lreg \
	    -title "Start position" \
	    -orient horizontal \
	    -from $left \
	    -to $right \
	    -default $start_value\
	    -width 7 \
	    -type CheckInt \
	    -command "CheckStartLimits $t.lreg $t.rreg 0"
    
    scalebox $t.rreg \
	    -title "End position" \
	    -orient horizontal \
	    -from $left \
	    -to $right\
	    -default $end_value \
	    -width 7 \
	    -type CheckInt \
	    -command "CheckEndLimits $t.lreg $t.rreg 0"

    yes_no $t.summary \
	-title "Only update information line" \
	-orient horizontal \
	-default [keylget gap_defs CONTIG_EDITOR.LIST_CONFIDENCE.ONLY_NUM_ERRS]

    okapplycancelhelp $t.ok \
        -ok_command "ed_list_confidence2 $e $t $t.lreg $t.rreg $t.summary; \
		     destroy $t" \
	-apply_command "ed_list_confidence2 $e $t $t.lreg $t.rreg $t.summary" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap4 {Editor-Comm-List Confidence}" \
        -bd 2 -relief groove

    pack $t.lreg $t.rreg $t.summary $t.ok -side top -fill both
}

proc ed_list_confidence2 {ed t lreg_w rreg_w summary_w} {
    set w [editor_to_ed $ed]
    global $w.LREG
    global $w.RREG
    global gap_defs

    set lreg [scalebox_get $lreg_w]
    set rreg [scalebox_get $rreg_w]
    set summary [yes_no_get $summary_w]
    $ed list_confidence $lreg $rreg $summary

    keylset gap_defs CONTIG_EDITOR.LIST_CONFIDENCE.ONLY_NUM_ERRS $summary

    set $w.LREG $lreg
    set $w.RREG $rreg
}

#
# Sets the tag types that are to be displayed
#
proc editor_setannos {e t} {
    TagDialog CONTIG_EDITOR.TAGS $t \
	"eval $e set_displayed_annos \[GetDefaultTags CONTIG_EDITOR.TAGS $e\];
	 global $e.VisAnnos;
	 if {\[info exists $e.VisAnnos\]} {unset $e.VisAnnos}"\
	$e
}

#
# Toggles visibility of annotations
#
proc toggle_annos {e} {
    global $e.VisAnnos
    if {[info exists $e.VisAnnos]} {
	eval $e set_displayed_annos [set $e.VisAnnos]
	unset $e.VisAnnos
    } else {
	set $e.VisAnnos [$e get_displayed_annos]
	$e set_displayed_annos
    }
}

#
# Turns on and off the status lines
#
proc editor_set_status {e w value args} {
    global $w.Status0 $w.Status7
    global $w.Status1 $w.Status2 $w.Status3 $w.Status4 $w.Status5 $w.Status6

    if {$value == 0} {
	foreach i $args {
    	    $e status delete $i
	    set $w.Status$i 0
	}
    } else {
	foreach i $args {
    	    $e status add $i
	    set $w.Status$i 1
	}
    }
}

#
# Sets the superedit mode
#
proc editor_set_superedit {e w} {
    global $w
    global $w.SE_ins_read $w.SE_del_read $w.SE_ins_cons $w.SE_del_dash_cons
    global $w.SE_del_any_cons $w.SE_read_shift $w.SE_trans_any
    global $w.SE_uppercase $w.SE_edit_mode $w.SE_replace_cons

    lappend mode \
	[set $w.SE_ins_read] \
	[set $w.SE_del_read] \
	[set $w.SE_ins_cons] \
	[set $w.SE_del_dash_cons] \
	[set $w.SE_del_any_cons] \
	[set $w.SE_replace_cons] \
	[set $w.SE_read_shift] \
	[set $w.SE_trans_any] \
	[set $w.SE_uppercase] \
	[set $w.SE_edit_mode]

    $e superedit $mode
}

proc editor_set_superedit_set {e w n} {
    global $w gap_defs
    global $w.SE_ins_read $w.SE_del_read $w.SE_ins_cons $w.SE_del_dash_cons
    global $w.SE_del_any_cons $w.SE_read_shift $w.SE_trans_any
    global $w.SE_uppercase $w.SE_edit_mode $w.SE_replace_cons

    set sets [keylget gap_defs CONTIG_EDITOR.SE_SET]
    set x [keylget sets $n]

    set $w.SE_ins_read		[lindex $x 0]
    set $w.SE_del_read		[lindex $x 1]
    set $w.SE_ins_cons		[lindex $x 2]
    set $w.SE_del_dash_cons	[lindex $x 3]
    set $w.SE_del_any_cons	[lindex $x 4]
    set $w.SE_replace_cons	[lindex $x 5]
    set $w.SE_read_shift	[lindex $x 6]
    set $w.SE_trans_any		[lindex $x 7]
    set $w.SE_uppercase		[lindex $x 8]
    set $w.SE_edit_mode		[lindex $x 9]

    $e superedit $x
}

#
# Sets the trace differencing mode
#
# Mode 0 == no diffs
# Mode 1 == diff against computed consensus
# Mode 2 == diff against specified trace
#
proc editor_trace_diff {e w} {
    global $w.TraceDiff

    if {[set $w.TraceDiff] == 2} {
	trace_diff_specific $w $e
    } elseif {[set $w.TraceDiff] == 1} {
	$e trace_comparator 0
    } else {
        $e trace_comparator
    }
}

#
# Configures the consensus trace modes
#
proc editor_trace_config {e w} {
    global $w.TraceConsMatch $w.TraceConsSelect
    global $w.TraceDiffAlgorithm $w.TraceDiffScale

    $e trace_config [set $w.TraceConsMatch] [set $w.TraceConsSelect] \
	[set $w.TraceDiffAlgorithm] [set $w.TraceDiffScale]
}

#
# Given a sequence number, this returns a list of sequence names starting
# at that sequence number and chaining right.
#
proc ednames_to_right {w seq_num} {
    set names [$w get_names_to_right $seq_num]
    set rnames ""
    foreach n $names {
	lappend rnames [lindex $n 1]
    }

    return $rnames
}


#
# Pops up a list of commands to run on a status line
#
proc ednames_menu {w x y X Y} {
    global licence
    set e [edname_to_editor $w]
    set ed [winfo parent $e]
    global $ed.Tlist
    set tlist [set $ed.Tlist]

    # Sequence/consensus names
    set seq_num [$w get_number @$x @$y]
    if {$seq_num == ""} {
	return
    }
    set name [$w get_name $seq_num]

    if {[winfo exists $w.m]} {destroy $w.m}
    if {[lindex $name 0] == "CONSENSUS"} {
	create_popup $w.m "Commands for consensus"
	set cnum [$w get_contig_number]
	$w.m add command -label "List notes" \
	    -command "NoteSelector [$e io] contig =$cnum"
    } else {
	set flags [$e get_flags $seq_num]
	set refseq [lsearch -exact $flags REFSEQ]
	set refseq [expr {$refseq == -1 ? 0 : 1}]
	set refneg [lsearch -exact $flags REFTRACE_NEG]
	set refpos [lsearch -exact $flags REFTRACE_POS]
	set reftrace [expr {($refneg == -1 && $refpos == -1) ? 0 : 1}]
	set tseqs [$e get_template_seqs $seq_num]

	create_popup $w.m "Commands for [lindex $name 1]"
	set rnum [$w get_read_number $seq_num]

	if {[llength $tseqs] != 0} {
	    # Goto menu
	    $w.m add cascade -label "Goto..." -menu $w.m.goto
	    menu $w.m.goto
	    set this_contig [$e get_contig_num]
	    set diff_contig_count 0
	    foreach {seq dummy_cnum pos} $tseqs {
		# Recalculate cnum. When joining contigs a contig may get
		# renumbered. It may or may not be this contig. If it is not
		# this contig and this contig is not getting joined, then this
		# contig will not receive any events and hence cannot update
		# its internal data structures to recalculate the contig
		# numbers.
		set cnum [db_info get_contig_num [$e io] $seq]
		if {$cnum == -1} {
		    verror ERR_WARN get_contig_num \
			"Failed to identify contig for sequence $seq"
		    set cnum $dummy_cnum
		}
		if {$cnum != $this_contig} {
		    set cname " (Contig [left_gel [$e io] $cnum] @ $pos)"
		    incr diff_contig_count
		} else {
		    set cname " @ $pos"
		}
		$w.m.goto add command -label "Goto $seq$cname" \
		    -command "ed_goto_seq $e $seq $cnum"
	    }

	    # Join to menu
	    if {$diff_contig_count > 0} {
		$w.m add cascade -label "Join to..." -menu $w.m.jointo
		menu $w.m.jointo
		foreach {seq dummy_cnum pos} $tseqs {
		    set cnum [db_info get_contig_num [$e io] $seq]
		    if {$cnum == -1 || $cnum == $this_contig} {
			continue
		    }
		    $w.m.jointo add command \
			-command "ed_jointo_seq $e $cnum $seq \
                                      [$w get_contig_number] \
                                      [lindex $name 1]" \
			-label "Join to $seq ([left_gel [$e io] $cnum] @ $pos)"
		}
	    }
	}

	$w.m add command -label "Select this reading" \
	    -command "editor_addlist_read [editor_to_edname $e] \
                      [lindex $name 1]"
	$w.m add command -label "Select this reading and all to right" \
	    -command "editor_addlist_read [editor_to_edname $e] \
                      \[ednames_to_right $w $seq_num\]"
	$w.m add command -label "Deselect this reading" \
	    -command "editor_dellist_read [editor_to_edname $e] \
                      [lindex $name 1]"
	$w.m add command -label "Deselect this reading and all to right" \
	    -command "editor_dellist_read [editor_to_edname $e] \
                      \[ednames_to_right $w $seq_num\]"
	$w.m add command -label "Select readings on this template" \
	    -command "editor_addlist_template [editor_to_edname $e] \
                      [list $tseqs]"
	$w.m add command -label "Deselect readings on this template" \
	    -command "editor_dellist_template [editor_to_edname $e] \
                      [list $tseqs]"
	$w.m add command -label "List notes" \
	    -command "NoteSelector [$e io] reading #$rnum"

	$w.m add separator
	if {!$refseq} {
	    $w.m add command -label "Set as reference sequence" \
		-command "ednames_menu_refseq $e $seq_num"
	} else {
	    $w.m add command -label "Clear as reference sequence" \
		-command "$e set_reference_seq $seq_num"
	}
	if {!$reftrace} {
	    $w.m add command -label "Set as reference trace" \
		-command "ednames_menu_reftrace $e $seq_num"
	} else {
	    $w.m add command -label "Clear as reference trace" \
		-command "$e set_reference_trace $seq_num 0"
	}

	if {$licence(type) == "f"} {
	    $w.m add separator
	    $w.m add command -label "Remove reading (this only)" \
		-command "$e hide_read $seq_num"
	    $w.m add command -label "Remove reading and all to right" \
		-command "$e hide_read -$seq_num"
	}
	$w.m add separator
	$w.m add command -label "Clear selection" \
	    -command "editor_clearlist [editor_to_edname $e]"
	$w.m add cascade -label "Move selection to..." \
	    -menu $w.m.moveto
	menu $w.m.moveto
	$w.m.moveto add command -label All \
	    -command "ed_move_set $e {$tlist} 0"
	$w.m.moveto add command -label New \
	    -command "ed_move_set $e {$tlist} +"
	set gnum 1
	foreach grp [set ::[editor_to_ed $e].Sets] {
	    $w.m.moveto add command -label "Set $gnum" \
		-command "ed_move_set $e {$tlist} $gnum"
	    incr gnum
	}
    }

    tk_popup $w.m [expr $X-20] [expr $Y-10]
}

# Function to jump to a specific sequence in the editor.
proc ed_goto_seq {e name cnum} {
    edit_contig -io [$e io] -contig =$cnum -reading $name -reuse 1
}

# Function to invoke join editor between a specific pair of sequences.
proc ed_jointo_seq {e cnum1 name1 cnum2 name2} {
    # check if contigs need complementing.
    set io [$e io]
    set rn1 [db_info get_read_num $io $name1]
    set rn2 [db_info get_read_num $io $name2]
    if {$rn1 < 1 || $rn2 < 1} {
	bell
	return
    }
    set r1 [io_read_reading $io $rn1]
    set d1 [string map {3 1 4 2} [keylget r1 primer]]
    set r2 [io_read_reading $io $rn2]
    set d2 [string map {3 1 4 2} [keylget r2 primer]]
    if {$d1 != 0 && $d2 != 0} {
	if {[keylget r1 sense] == [keylget r2 sense]} {
	    if {$d1 == $d2} {
		set comp 0
	    } else {
		set comp 1
	    }
	} else {
	    if {$d1 == $d2} {
		set comp 1
	    } else {
		set comp 0
	    }
	}
    } else {
	set comp 0
    }

    # Complement if needed
    global read_only
    if {!$read_only && $comp} {
	complement_contig -io $io -contigs $name1
    }

    # Bring up the join editor
    join_contig \
	-io $io \
	-contig1 =$cnum1 \
	-reading1 $name1 \
	-contig2 =$cnum2 \
	-reading2 $name2 
}

# Function to move sequences between sets
proc ed_move_set {e list grp} {
    set ed [editor_to_ed $e]
    set n [editor_to_edname $e]
    global $ed.Sets $ed.Tlist
    global $n.List

    if {$list == ""} {
	set $ed.Tlist [editor_tabs_create $ed.tabs $e ""]
	set list [set $ed.Tlist]
    } elseif {![winfo exists [set $ed.Tlist]]} {
	set $ed.Tlist [editor_tabs_create $ed.tabs $e]
    }

    if {[info exists $n.List]} {
	set sel [ListGet [set $n.List]]
    } else {
	set sel [ListGet readings]
    }

    set sets [set $ed.Sets]

    # Remove from current sets
    set newsets ""
    foreach g $sets {
	foreach r $sel {
	    set ind [lsearch $g $r]
	    if {$ind >= 0} {
		set g [lreplace $g $ind $ind]
	    }
	}
	lappend newsets $g
    }

    # Move to new sets
    if {$grp == "+"} {
	lappend newsets $sel
	set grp [llength $newsets]
    } elseif {$grp != 0} {
	set g [lindex $newsets [expr {$grp-1}]]
	foreach r $sel {
	    lappend g $r
	}
	set newsets [lreplace $newsets [expr {$grp-1}] [expr {$grp-1}] $g]
    }

    set $ed.Sets $newsets

    # Update displays; both editor structure and tablelist
    $e move_to_set $grp $sel
    editor_tabs_repopulate $list $newsets
}

# Dialogue for setting a reference sequence
proc ednames_menu_refseq {e seq_num} {
    set w $e.win
    if {[xtoplevel $w -resizable 0] == ""} {return}
    wm title $w "[string trim [$e get_name $seq_num]]"

    xentry $w.number \
	-label "First base number" \
	-default 1 \
	-width 10 \
	-type "int 0" \

    xentry $w.length \
	-label "Sequence length" \
	-width 10 \
	-type "int 1" \

    xyn $w.xyn \
	-label "Circular sequence?" \
	-orient horiz \
	-ycommand "$w.length configure -state normal" \
	-ncommand "$w.length configure -state disabled"

    $w.xyn set 0

    okcancelhelp $w.ok \
	    -ok_command "ednames_menu_refseq_ok $e $seq_num $w" \
	    -cancel_command "destroy $w" \
	    -help_command "show_help gap4 {Editor-Reference-Sequences}" \
	    -bd 2 \
	    -relief groove

    pack $w.number $w.xyn $w.length $w.ok -fill both
}

proc ednames_menu_refseq_ok {e seq_num w} {
    set offset [$w.number get]
    if {[$w.xyn get]} {
	set length [$w.length get]
    } else {
	set length 0
    }
    
    if {$offset == "" || $length == ""} {
	bell
	return
    }
    
    destroy $w
    $e set_reference_seq $seq_num $length $offset
}

# Dialogue for setting a reference trace
proc ednames_menu_reftrace {e seq_num} {
    set w $e.win
    if {[xtoplevel $w -resizable 0] == ""} {return}
    wm title $w "[string trim [$e get_name $seq_num]]"

    global $w.control
    if {![info exists $w.control]} {
	set $w.control -
    }

    set f [frame $w.controlf]
    radiobutton $f.p \
	-text "Positive control"\
	-variable $w.control \
	-value +
    radiobutton $f.m \
	-text "Negative control"\
	-variable $w.control \
	-value -
    pack $f.m $f.p -side top -anchor w

    okcancelhelp $w.ok \
	    -ok_command "ednames_menu_reftrace_ok $e $seq_num $w" \
	    -cancel_command "destroy $w" \
	    -help_command "show_help gap4 {Editor-Reference-Traces}" \
	    -bd 2 \
	    -relief groove

    pack $w.controlf $w.ok -fill both
}

proc ednames_menu_reftrace_ok {e seq_num w} {
    global $w.control
    if {[set $w.control] == "+"} {
	set control 1
    } elseif {[set $w.control] == "-"} {
	set control -1
    } else {
	set control 0
    }
    $e set_reference_trace $seq_num $control
    destroy $w
}

#
# Creates a dialogue for the "Report Mutations" command
#
proc report_mutations_dialog {e} {
    global gap_defs

    set w [editor_to_ed $e].mutations
    if {[xtoplevel $w -resizable 0] == ""} {return}
    wm title $w "Report Mutations"

    radiolist $w.tags \
	-title "Find mutations by" \
	-bd 2 -relief groove -orient horizontal \
	-default [keylget gap_defs CONTIG_EDITOR.MUTATIONS_TAGGED] \
	-buttons {{Tags} {Differences}}

    radiolist $w.position \
	-title "Sort by" \
	-bd 2 -relief groove -orient horizontal \
	-default [keylget gap_defs CONTIG_EDITOR.MUTATIONS_SORT] \
	-buttons {{Position} {Sequence}}

    xentry $w.directory \
        -label "HTML report directory" \
	-default "mutation_report" \
	-type directoryoutput

    checkbutton $w.detailed \
        -text "Detailed HTML report" \
	-variable $w.Detailed
    global $w.Detailed
    set $w.Detailed 1

    okcancelhelp $w.ok \
	-ok_command "
	    report_mutations_dialog2 $w $e \
		\[expr {2-\[radiolist_get $w.tags\]}\] \
		\[expr {2-\[radiolist_get $w.position\]}\]" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap4 {Editor-Comm-Report-Mutations}" \
	-bd 2 \
	-relief groove

    pack $w.tags $w.position $w.directory -fill both -expand 1
    pack $w.detailed -anchor w
    pack $w.ok -fill both -expand 1
}

proc report_mutations_dialog2 {w e tags pos} {
    global $w.Detailed
    set detailed [set $w.Detailed]
    set dir [$w.directory get]
    if {$dir == "" && [$w.directory get2] != ""} {
	bell
	return
    }

    if {$dir == ""} {
	set detailed 0
    } else {
	incr detailed
    }

    if {$dir != ""} {
	catch {file mkdir $dir}
    }

    destroy $w
    if {[catch {set html [$e report_mutations $tags $pos $dir $detailed]} err]} {
	tk_messageBox \
	    -icon error \
	    -message $err \
	    -title "Report Mutations" \
	    -type ok \
	    -parent $e
	return
    }

    if {$dir != ""} {
	set fd [open $dir/index.html w]
	puts $fd $html
	close $fd
    }
}

#
# Dump contig user interface
#
proc dump_contig_dialog { e } {
    set w [editor_to_ed $e]
    global gap_defs db_namelen
    global $w.LREG
    global $w.RREG

    set f [keylget gap_defs DUMP_CONTIG.WIN]

    if {[xtoplevel $f -resizable 0] == ""} {return}
    wm title $f "Dump contig to file"

    set end_value [set $w.RREG]
    set start_value [set $w.LREG]

    ###########################################################################    #output file
    getFname $f.output [keylget gap_defs DUMP_CONTIG.NAME] save {} \
	[keylget gap_defs DUMP_CONTIG.VALUE]
    
    set extents [$e get_extents]
    set left [lindex $extents 0]
    set right [lindex $extents 1]
    scalebox $f.lreg \
	    -title "Start position" \
	    -orient horizontal \
	    -from $left \
	    -to $right \
	    -default $start_value\
	    -width 7 \
	    -type CheckInt \
	    -command "CheckStartLimits $f.lreg $f.rreg 0"
    
    scalebox $f.rreg \
	    -title "End position" \
	    -orient horizontal \
	    -from $left \
	    -to $right\
	    -default $end_value \
	    -width 7 \
	    -type CheckInt \
	    -command "CheckEndLimits $f.lreg $f.rreg 0"

    scalebox $f.llength \
	    -title "Line length" \
	    -orient horizontal \
	    -from 20 \
	    -to 1000 \
	    -default [keylget gap_defs DUMP_CONTIG.LINE_LENGTH] \
    	    -width 7 \
	    -type CheckInt

    scalebox $f.nwidth \
	    -title "Sequence name width" \
	    -orient horizontal \
	    -from 0 \
	    -to $db_namelen \
	    -default [keylget gap_defs DUMP_CONTIG.NAME_WIDTH] \
	    -width 7 \
	    -type CheckInt

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "Dump_OK_Pressed $e $f $f.output $f.lreg $f.rreg \
			 $f.llength $f.nwidth" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Editor-Comm-Dump}" \
	    -bd 2 \
	    -relief groove

    ###########################################################################
    #packing
    pack $f.output $f.lreg $f.rreg $f.llength $f.nwidth $f.ok_cancel -fill x
}

proc Dump_OK_Pressed {e f output lreg rreg llength nwidth} {
    set w [editor_to_ed $e]
    global $w.LREG
    global $w.RREG
    global gap_defs

    set out_file [getFname_in_name $output]
    if {$out_file == ""} return

    set $w.LREG [scalebox_get $lreg] 
    set $w.RREG [scalebox_get $rreg]
    set llen [scalebox_get $llength]
    set nwid [scalebox_get $nwidth]
    keylset gap_defs DUMP_CONTIG.LINE_LENGTH $llen
    keylset gap_defs DUMP_CONTIG.NAME_WIDTH $nwid
    destroy $f

    $e dump_contig $out_file [set $w.LREG] [set $w.RREG] $llen $nwid
}

#
# Save Consensus Trace dialog
#
proc consensus_trace_dialog { e } {
    set w [editor_to_ed $e]
    global gap_defs
    global $w.LREG
    global $w.RREG

    set f [keylget gap_defs SAVE_CON_TRACE.WIN]

    if {[xtoplevel $f -resizable 0] == ""} {return}
    wm title $f "Save consensus trace"

    set end_value [set $w.RREG]
    set start_value [set $w.LREG]

    ###########################################################################    #output file
    getFname $f.output [keylget gap_defs SAVE_CON_TRACE.NAME] save {} \
	[keylget gap_defs SAVE_CON_TRACE.VALUE]
    
    set extents [$e get_extents]
    set left [lindex $extents 0]
    set right [lindex $extents 1]
    scalebox $f.lreg \
	    -title "Start position" \
	    -orient horizontal \
	    -from $left \
	    -to $right \
	    -default $start_value\
	    -width 7 \
	    -type CheckInt \
	    -command "CheckStartLimits $f.lreg $f.rreg 0"
    
    scalebox $f.rreg \
	    -title "End position" \
	    -orient horizontal \
	    -from $left \
	    -to $right\
	    -default $end_value \
	    -width 7 \
	    -type CheckInt \
	    -command "CheckEndLimits $f.lreg $f.rreg 0"

    radiolist $f.strand \
	-title "strand" \
	-bd 2 -relief groove -orient horizontal \
	-default [keylget gap_defs SAVE_CON_TRACE.STRAND] \
	-buttons {{Forward} {Reverse}}

    yes_no $f.matching \
	-title "Use only matching reads" \
	-bd 2 -relief groove \
	-orient horizontal \
	-default [keylget gap_defs SAVE_CON_TRACE.MATCHING] \
	

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "consensus_trace_dialog2 $e $f $f.output $f.lreg \
		$f.rreg $f.strand $f.matching" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Editor-Comm-Consensus Trace}" \
	    -bd 2 \
	    -relief groove

    ###########################################################################
    #packing
    pack $f.output $f.lreg $f.rreg $f.strand $f.matching $f.ok_cancel -fill x
}

proc consensus_trace_dialog2 {e f output lreg rreg strand matching} {
    set w [editor_to_ed $e]
    global $w.LREG
    global $w.RREG

    set out_file [getFname_in_name $output]
    if {$out_file == ""} return

    set $w.LREG [scalebox_get $lreg] 
    set $w.RREG [scalebox_get $rreg]
    set matching [yes_no_get $matching]
    set strand [radiolist_get $strand]
    if {$strand == 2} {set strand -1}
    destroy $f

    $e consensus_trace $out_file [set $w.LREG] [set $w.RREG] $strand $matching
}

proc set_editor_write_mode {e mode} {
    global licence

    set w [editor_to_ed $e]
    if {$mode == 0} {
	menu_state_set contig_editor_settings_menu  -2 $w.buttons.settings.menu
	menu_state_set contig_editor_commands_menu  -2 $w.buttons.commands.menu
	menu_state_set contig_editor_help_menu      -2 $w.buttons.help.menu
	if {$licence(type) != "v"} {
	    menu_state_set contig_editor_editmodes_menu \
		-2 $w.buttons.editmodes.menu
	}
    } else {
        menu_state_set contig_editor_settings_menu  2 $w.buttons.settings.menu
        menu_state_set contig_editor_commands_menu  2 $w.buttons.commands.menu
        menu_state_set contig_editor_help_menu      2 $w.buttons.help.menu
	if {$licence(type) != "v"} {
	    menu_state_set contig_editor_editmodes_menu \
		2 $w.buttons.editmodes.menu
	}
    }
}

proc editor_shuffle_pads {e} {
    global gap_defs

    $e shuffle_pads \
	    [keylget gap_defs CONTIG_EDITOR.SHUFFLE_PADS.CONSENSUS_MODE] \
	    [keylget gap_defs CONTIG_EDITOR.SHUFFLE_PADS.CONSENSUS_CUTOFF]
}

proc editor_set_status_line {e text {delay 0}} {
    [winfo parent [winfo parent $e]].status.l configure -text $text
    if {$delay != 0} {
	after $delay "catch {editor_set_status_line $e {}}"
    }
}

proc editor_change_consensus_algorithm {e w} {
    global $w.ConsensusAlgorithm

    $e set_consensus_mode [set $w.ConsensusAlgorithm]
}

proc save_editor_settings {e w} {
    puts w=$w
    global gap_defs env
    global $w.ShowQuality $w.AminoMode $w.DisplayTraces $w.AutoSave
    global $w.SE_ins_read $w.SE_del_read $w.SE_ins_cons $w.SE_del_dash_cons
    global $w.SE_del_any_cons $w.SE_read_shift $w.SE_trans_any
    global $w.SE_uppercase $w.SE_edit_mode $w.SE_replace_cons
    global $w.TraceDiff $w.TraceConsMatch $w.TraceConsSelect $w.DisagreeMode
    global $w.TraceDiffAlgorithm $w.TraceDiffScale $w.GroupBy
    global $w.ShowEdits $w.Disagreements $w.CompareStrands
    global $w.ShowCQuality $w.DisagreeCase $w.TemplateNames
    global $w.Status0 $w.Status1 $w.Status2 $w.Status3
    global $w.Status4 $w.Status5 $w.Status6 $w.Status7
    global $w.ShowUnpadded $w.DiffTraces $w.ReadPairTraces $w.DisplayMiniTraces

    # Update the in-memory definitions so that subsequent editors keep these
    # settings.
    set C CONTIG_EDITOR
    keylset gap_defs $C.DISAGREEMENTS          [set $w.Disagreements]
    keylset gap_defs $C.DISAGREE_MODE          [set $w.DisagreeMode]
    keylset gap_defs $C.DISAGREE_CASE          [set $w.DisagreeCase]
    keylset gap_defs $C.COMPARE_STRANDS        [set $w.CompareStrands]
    keylset gap_defs $C.AUTO_DISPLAY_TRACES    [set $w.DisplayTraces]
    keylset gap_defs $C.DISPLAY_MINI_TRACES    [set $w.DisplayMiniTraces]
    keylset gap_defs $C.AUTO_DIFF_TRACES       [set $w.DiffTraces]
    keylset gap_defs $C.READ_PAIR_TRACES       [set $w.ReadPairTraces]
    keylset gap_defs $C.AUTO_SAVE              [set $w.AutoSave]
    keylset gap_defs $C.STATUS_STRAND          [set $w.Status0]
    keylset gap_defs $C.STATUS_FRAME1P         [set $w.Status1]
    keylset gap_defs $C.STATUS_FRAME2P         [set $w.Status2]
    keylset gap_defs $C.STATUS_FRAME3P         [set $w.Status3]
    keylset gap_defs $C.STATUS_FRAME1M         [set $w.Status4]
    keylset gap_defs $C.STATUS_FRAME2M         [set $w.Status5]
    keylset gap_defs $C.STATUS_FRAME3M         [set $w.Status6]
    keylset gap_defs $C.STATUS_AUTO_TRANSLATE  [set $w.Status7]
    keylset gap_defs $C.AMINO_ACID_MODE        [set $w.AminoMode]
    keylset gap_defs $C.GROUP_BY	       [set $w.GroupBy]
    keylset gap_defs $C.SHOW_QUALITY           [set $w.ShowQuality]
    keylset gap_defs $C.SHOW_CONSENSUS_QUALITY [set $w.ShowCQuality]
    keylset gap_defs $C.SHOW_UNPADDED          [set $w.ShowUnpadded]
    keylset gap_defs $C.SHOW_EDITS             [set $w.ShowEdits]
    keylset gap_defs $C.TRACE_DIFF             [set $w.TraceDiff]
    keylset gap_defs $C.TRACE_CONS_MATCH       [set $w.TraceConsMatch]
    keylset gap_defs $C.TRACE_CONS_SELECT      [set $w.TraceConsSelect]
    keylset gap_defs $C.TRACE_DIFF_ALGORITHM   [set $w.TraceDiffAlgorithm]
    keylset gap_defs $C.TRACE_DIFF_SCALE       [set $w.TraceDiffScale]
    keylset gap_defs $C.SE_INS_ANY             [set $w.SE_ins_read]
    keylset gap_defs $C.SE_DEL_READ            [set $w.SE_del_read]
    keylset gap_defs $C.SE_INS_CONS            [set $w.SE_ins_cons]
    keylset gap_defs $C.SE_DEL_DASH_CONS       [set $w.SE_del_dash_cons]
    keylset gap_defs $C.SE_DEL_ANY_CONS        [set $w.SE_del_any_cons]
    keylset gap_defs $C.SE_REPLACE_CONS        [set $w.SE_replace_cons]
    keylset gap_defs $C.SE_READ_SHIFT          [set $w.SE_read_shift]
    keylset gap_defs $C.SE_TRANS_ANY           [set $w.SE_trans_any]
    keylset gap_defs $C.SE_UPPERCASE           [set $w.SE_uppercase]
    keylset gap_defs $C.SE_EDIT_MODE           [set $w.SE_edit_mode]
    keylset gap_defs $C.TEMPLATE_NAMES         [set $w.TemplateNames]

    # Write to the .gaprc file
    update_defs gap_defs $env(HOME)/.gaprc \
	$C.DISAGREEMENTS \
	$C.DISAGREE_MODE \
	$C.DISAGREE_CASE \
	$C.DISAGREE_QUAL \
	$C.COMPARE_STRANDS \
	$C.AUTO_DISPLAY_TRACES \
	$C.AUTO_DIFF_TRACES \
	$C.READ_PAIR_TRACES \
	$C.AUTO_SAVE \
	$C.STATUS_STRAND \
	$C.STATUS_FRAME1P \
	$C.STATUS_FRAME2P \
	$C.STATUS_FRAME3P \
	$C.STATUS_FRAME1M \
	$C.STATUS_FRAME2M \
	$C.STATUS_FRAME3M \
	$C.AMINO_ACID_MODE \
	$C.GROUP_BY \
	$C.SHOW_QUALITY \
	$C.SHOW_CONSENSUS_QUALITY \
	$C.SHOW_UNPADDED \
	$C.SHOW_EDITS \
	$C.TRACE_DIFF \
	$C.TRACE_CONS_MATCH \
	$C.TRACE_CONS_SELECT \
	$C.TRACE_DIFF_ALGORITHM \
	$C.TRACE_DIFF_SCALE \
	$C.DISPLAY_MINI_TRACES \
	$C.SE_INS_ANY \
	$C.SE_DEL_READ \
	$C.SE_INS_CONS \
	$C.SE_DEL_DASH_CONS \
	$C.SE_DEL_ANY_CONS \
	$C.SE_REPLACE_CONS \
	$C.SE_READ_SHIFT \
	$C.SE_TRANS_ANY \
	$C.SE_UPPERCASE \
	$C.SE_EDIT_MODE \
	$C.TEMPLATE_NAMES

    # Also save the tag macros (tag_editor.tcl)
    tag_macro_save
}

proc set_mini_traces {ed val} {
     set w [winfo parent [edname_to_editor $ed]]
     global $w.DisplayMiniTraces
     set $w.DisplayMiniTraces $val
     $ed show_mini_traces $val
}

# Mousewheel callback
proc editor_mouse_wheel {e dist speed} {
    global gap_defs
    switch $speed {
	"slow" {
	    set step 1
	}
	"fast" {
	    set step 20
	}
	default {
	    set step [keylget gap_defs CONTIG_EDITOR.MOUSEWHEEL_SPEED]
	}
    }

    if {$dist > 0} {
	set dist $step
    } elseif {$dist < 0} {
	set dist -$step
    }

    $e yview scroll $dist units
}

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
    if {[catch {scroll_editor $e scroll [set $e.AutoScroll] units}] == 0} {
	if {[set $e.AutoScroll] > 0} {
	    $e select adjust [lindex [$e configure -width] end]
	} else {
	    $e select adjust 0
	}
	set $e.AutoScrollEvent [after 50 "editor_autoscroll $e"]
    }
}

proc editor_autoscroll_cancel {e} {
    global $e.AutoScrollEvent
    catch {after cancel [set $e.AutoScrollEvent]; unset $e.AutoScrollEvent}
}

#
# Add bindings
#
bind Editor <<select>>	{
    focus %W
    if {[%W cursor_set @%x @%y] == 0} {
	%W select from @%x
    } else {
	%W select clear
    }
    %W update_brief_base @%x @%y
}

if {[keylget gap_defs CONTIG_EDITOR.AUTO_FOCUS] == 1} {
    bind Editor <Any-Enter> {focus %W}
}

bind Editor <Key-Return>	{
    %W cursor_set @%x @%y
    %W update_brief_base -mode2 @%x @%y
}
bind Editor <Key-KP_Enter>	{%W update_brief_base -mode2}
bind Editor <Control-1>		{popup_editor_menu %W %x %y}
#bind Editor <<menu>>		"[bind Editor <<select>>]
#				 popup_editor_menu %W %x %y"
bind Editor <<menu>>		{
    %W cursor_set @%x @%y
    popup_editor_menu %W %x %y
}
bind Editor <<select-drag>>	{%W select adjust @%x;
				 editor_select_scroll %W %x}
bind Editor <<select-release>>	{editor_autoscroll_cancel %W}
bind Editor <<move>>		{
    %W cursor_set @%x @%y
    %W update_brief_base
}
bind Editor <<trace>> {
    if {[%W cursor_set @%x @%y] == 0} {
	%W invoke_trace
    }
}
bind Editor <Any-Motion>	{%W update_brief @%x @%y}
bind Editor <Shift-Motion>	{;}
bind Editor <Control-Key-t>	{%W invoke_trace}
bind Editor <Control-Key-q>	{toggle_annos %W}

bind Editor <<select-to>>	{%W select to @%x}

bind Editor <Key-Left>		{%W cursor_left;  %W update_brief_base}
bind Editor <Control-Key-b>	{%W cursor_left;  %W update_brief_base}

bind Editor <Key-Right>		{%W cursor_right; %W update_brief_base}
bind Editor <Control-Key-f>	{%W cursor_right; %W update_brief_base}

bind Editor <Key-Up>		{%W cursor_up;    %W update_brief_base}
bind Editor <Control-Key-p>	{%W cursor_up;    %W update_brief_base}

bind Editor <Key-Down>		{%W cursor_down;  %W update_brief_base}
bind Editor <Control-Key-n>	{%W cursor_down;  %W update_brief_base}

# 22/1/99 johnt - KP_Down/UP/Right/Left keysyms are not defined on some systems
catch {
  bind Editor <Key-KP_Down>	{%W cursor_down;  %W update_brief_base}
  bind Editor <Key-KP_Up>	{%W cursor_up;    %W update_brief_base}
  bind Editor <Key-KP_Right>	{%W cursor_right; %W update_brief_base}
  bind Editor <Key-KP_Left>	{%W cursor_left;  %W update_brief_base}
}

bind Editor <Control-Key-a>	{%W read_start}
bind Editor <Control-Key-e>	{%W read_end}
bind Editor <Meta-Key-a>	{%W read_start2}
bind Editor <Alt-Key-a>		{%W read_start2}
bind Editor <Escape><Key-a>	{%W read_start2}
bind Editor <Meta-Key-e>	{%W read_end2}
bind Editor <Alt-Key-e>		{%W read_end2}
bind Editor <Escape><Key-e>	{%W read_end2}
bind Editor <Meta-Key-comma>	{%W contig_start}
bind Editor <Alt-Key-comma>	{%W contig_start}
bind Editor <Escape><Key-comma>	{%W contig_start}
bind Editor <Meta-Key-period>	{%W contig_end}
bind Editor <Alt-Key-period>	{%W contig_end}
bind Editor <Escape><Key-period> {%W contig_end}

if {$licence(type) != "v"} {
  bind Editor <Control-Key-l>	{%W transpose_left}
  bind Editor <Control-Key-r>	{%W transpose_right}

  bind Editor <Key-bracketleft>	{%W set_confidence 0;
				 %W update_brief_base}
  bind Editor <Key-bracketright> {%W set_confidence 100;
				 %W update_brief_base}
  bind Editor <Shift-Key-Up>	{%W increment_confidence 1;
				 %W update_brief_base}
  bind Editor <Control-Key-Up>	{%W increment_confidence 10;
				 %W update_brief_base}
  bind Editor <Shift-Key-Down>	{%W increment_confidence -1;
				 %W update_brief_base}
  bind Editor <Control-Key-Down> {%W increment_confidence -10;
				 %W update_brief_base}

  # 22/1/99 johnt - KP_Down/UP/Right/Left keysyms are not defined on some
  # systems
  catch {
    bind Editor <Shift-Key-KP_Up>	{%W increment_confidence 1;
				 	%W update_brief_base}
    bind Editor <Control-Key-KP_Up>	{%W increment_confidence 10;
				 	%W update_brief_base}
    bind Editor <Shift-Key-KP_Down>	{%W increment_confidence -1;
				 	%W update_brief_base}
    bind Editor <Control-Key-KP_Down>	{%W increment_confidence -10;
				 	%W update_brief_base}
  }

  bind Editor <BackSpace>		{%W delete_key}
  bind Editor <Control-Key-Delete>	{%W delete_left_key}
  catch {bind Editor <Control-Key-KP_Delete>	{%W delete_left_key}}
  bind Editor <Control-Key-d>		{%W cursor_right;
  					 %W delete_key;
					 %W update_brief_base}
  if {[keylget gap_defs CONTIG_EDITOR.DELETE_RIGHT] == 1} {
    bind Editor <Delete>		{%W cursor_right;
  					 %W delete_key;
					 %W update_brief_base}
    catch {bind Editor <KP_Delete>	{%W cursor_right;
  					 %W delete_key;
					 %W update_brief_base}}
  } else {
    bind Editor <Delete>		{%W delete_key}
    catch {bind Editor <KP_Delete>	{%W delete_key}}
  }

  bind Editor <Meta-Key-Left>		{%W extend_left}
  bind Editor <Control-Key-Left>	{%W extend_left}
  bind Editor <Alt-Key-Left>		{%W extend_left}
  bind Editor <Meta-Key-Right>		{%W extend_right}
  bind Editor <Control-Key-Right>	{%W extend_right}
  bind Editor <Alt-Key-Right>		{%W extend_right}
  bind Editor <Key-less>		{%W zap_left}
  bind Editor <Key-greater>		{%W zap_right}

  bind Editor <<save>>		{%W save}

  # 'Named key' bindings. Use catch incase they do not exist on all systems.
  if {[keylget gap_defs CONTIG_EDITOR.MOVE_ON_EDITS]} {
      # The [string length] avoids having Control, Shift, etc sent over as
      # edit commands. We only send an edit request if it generates an ascii
      # symbol.
      bind Editor <Key>		{if {[string length %A]} {%W edit_key %A}}
      catch {bind Editor <Key-Insert>	{%W edit_key *}}
  } else {
      bind Editor <Key>		{if {[string length %A]} {%W edit_key %A}}
      catch {bind Editor <Key-Insert>	{%W edit_key -nomove *}}
  }
  catch {bind Editor <Key-Undo>		{%W undo}}
  catch {bind Editor <Key-DRemove>	{%W cursor_right;
					 %W delete_key;
					 %W update_brief_base}}

  bind Editor <<toggle-read-only>> \
	{set_editor_write_mode %W [%W write_mode -1]}
}

proc editor_collapse {en x y} {
    set ed [edname_to_editor $en]
    incr y
    set set [$ed get_set [$ed get_number $x $y]]
    if {$set < 0} {
	set set [expr {-$set}]
	$ed collapse_set $set -1
    }
}

bind Editor <Meta-Key-v>	{scroll_ll %W [editor_to_ed %W].scrollx;
				 scroll_ll %W [editor_to_ed %W].scrollx}
bind Editor <Alt-Key-v>		{scroll_ll %W [editor_to_ed %W].scrollx;
				 scroll_ll %W [editor_to_ed %W].scrollx}
bind Editor <Escape><Key-v>	{scroll_ll %W [editor_to_ed %W].scrollx;
				 scroll_ll %W [editor_to_ed %W].scrollx}
bind Editor <Control-Key-v>	{scroll_rr %W [editor_to_ed %W].scrollx;
				 scroll_rr %W [editor_to_ed %W].scrollx}

bind Editor <Control-Key-i>	{
	%W set_insert
	set [editor_to_ed %W].Insert [lindex {1 0} \
	    [set [editor_to_ed %W].Insert]]
}

bind Editor <Control-Key-underscore> {%W undo}
bind Editor <Control-Key-h>	{%W hide_read}

bind Editor <Control-Key-x>	{#noop;}
bind Editor <Control-Key-x><Control-Key-g> {#noop;}
bind Editor <<search>>		{create_search_win %W.search "%W search" 1}
bind Editor <<rsearch>>		{create_search_win %W.search "%W search" -1}
bind Editor <<copy>>		{
    if {[selection own] == "%W"} {
	clipboard clear
	clipboard append [selection get]
    }
}
bind Editor <<paste>>		{break;}

#bind EdNames <<select>>		{%W highlight -1 @%x @%y}
bind EdNames <<select>>		{
    foreach {x y} [%W coords @%x @%y] break;
    if {$x == 0} {
	editor_collapse %W $x $y
    } else {
	editor_addlist %W @%x @%y
    }
}
bind EdNames <<move>>		{editor_addlist %W @%x @%y}
bind EdNames <<copy>>		{copy_name %W @%x @%y}
bind EdNames <<menu>>		{ednames_menu %W %x %y %X %Y}
bind EdNames <Any-Motion>	{%W update_brief @%x @%y}
bind EdNames <Shift-Motion>	{;}
bind EdNames <Control-Key-h>	{[edname_to_editor %W] hide_read}

# 'Named key' bindings. Use catch incase they do not exist on all systems.
catch {bind Editor <Key-Next>	{create_search_win %W.search "%W search"  1 0}}
catch {bind Editor <Key-Prior>	{create_search_win %W.search "%W search" -1 0}}
catch {bind Editor <Key-Find>	{create_search_win %W.search "%W search"  1 0}}
catch {bind Editor <Key-Begin>	{%W read_start}}
catch {bind Editor <Key-Home>	{%W read_start}}
catch {bind Editor <Key-End>	{%W read_end}}
catch {bind Editor <Control-Key-Begin>	{%W read_start2}}
catch {bind Editor <Control-Key-End>	{%W read_end2}}

# Mouse wheel support
bind Editor <MouseWheel> {editor_mouse_wheel %W %D normal}
if {[tk windowingsystem] eq "x11"} {
    bind Editor <Shift-4>   {editor_mouse_wheel %W -200 slow}
    bind Editor <4>         {editor_mouse_wheel %W -200 normal}
    bind Editor <Control-4> {editor_mouse_wheel %W -200 fast}
    bind Editor <Shift-5>   {editor_mouse_wheel %W +200 slow}
    bind Editor <5>         {editor_mouse_wheel %W +200 normal}
    bind Editor <Control-5> {editor_mouse_wheel %W +200 fast}
}

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
bind Editor <F11> "%W select clear; %W edit_anno"
bind Editor <Shift-Key-F11> {
    if {[%W cursor_set @%x @%y] == 0} {
	%W select clear
	%W edit_anno
    }
}
bind Editor <F12> {
    if {[set .cedit.SE_fast_delete_anno]} {
	%W select clear
	%W delete_anno
    }
}
bind Editor <Shift-Key-F12> {
    if {[set .cedit.SE_fast_delete_anno]} {
	if {[%W cursor_set @%x @%y] == 0} {
	    %W select clear
	    %W delete_anno
	}
    }
}

bind Editor <Control-Key-0> {set_mini_traces %W 0}
bind Editor <Control-Key-1> {set_mini_traces %W 1}
bind Editor <Control-Key-2> {set_mini_traces %W 2}
bind Editor <Control-Key-3> {set_mini_traces %W 3}
bind Editor <Control-Key-4> {set_mini_traces %W 4}
bind Editor <Control-Key-5> {set_mini_traces %W 5}

