#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#==============================================================================
# Externally usable functions
#==============================================================================

# FIXME: direction of top fails at present due to bugs in the tk grid command.
# We override this and force users to add traces to the bottom of the display.
# It's not really a limitation now as we have the Y scrollbar.
keylset gap5_defs TRACE_DISPLAY.DIRECTION bottom

# Height of the trace widget including scrollbar
set trace_height [keylget gap5_defs TRACE_DISPLAY.TRACE_HEIGHT]

# The constant width of the trace display canvas (totalling all columns)
set trace_display_width [keylget gap5_defs TRACE_DISPLAY.WINDOW_WIDTH]
if {$trace_display_width > [expr [winfo screenwidth .]-30]} {
    set trace_display_width [expr [winfo screenwidth .]-30]
}

# Compact mode, where the command buttons are replaced by popup menus
set trace_compact_mode [keylget gap5_defs TRACE_DISPLAY.COMPACT_MODE]

# for multi-column trace displays
if {![info exists trace_columns]} {
    set trace_columns [keylget gap5_defs TRACE_DISPLAY.COLUMNS]
}
set trace_rows    [keylget gap5_defs TRACE_DISPLAY.ROWS]

# Whether to show the confidence values
set trace_confidence [keylget gap5_defs TRACE_DISPLAY.SHOW_CONFIDENCE]

#
# Adds a new trace 'file', to the toplevel 'w' for editor 'e'
#
# This returns the pathname of the trace widget created.
#
proc trace_add {w file e title} {
    global trace_columns

    if {[catch {set t [new_trace_create2 $w $e $title 1]} err]} {
	puts $err
    }

    if {[catch {$t.trace load $file}] == 1} {
	# Trace not found, or corrupted
	destroy $t
	error ""
    }

    set format [$t.trace format]
    if {$format == "PLN" || $format == "EXP"} {
	destroy $t
	error ""
    }

    if {[wm state [winfo toplevel $t.trace]] == "iconic"} {
        wm deiconify [winfo toplevel $t.trace]
    }

    return $t.trace
}

#
# Creates a new trace display, but with no trace loaded.
#
proc trace_create {w e title} {
    set t [new_trace_create2 $w $e $title 0]

    return $t.trace
}

#==============================================================================
# Internally usable functions
#==============================================================================

proc recurse_trace_bind {w} {
    foreach i [winfo children $w] {
	bindtags $i "tdiff [bindtags $i]"
	if {"[winfo children $i]" != ""} {
	    recurse_trace_bind $i
	}
    }
}

#
# Destroy the trace with dnatrace pathname of $t
# Called from C when (for example) double clicking on consensus to get all
# traces; we first remove the present ones
#
proc dnatrace_remove {t} {
    destroy [winfo parent $t]
}

#
# Remove the trace with frame $f from the top level display $top
# 
proc trace_remove {f w e} {
    global $w.after_id

    # This removes from editor and this display
    $e delete_trace $f.trace

    if {![winfo exists $w]} {return}

    global $w.NTraces trace_height $w.Order

    if {[llength [winfo children $w.traces.c.f]] == 1} {
	destroy $w
	unset $w.NTraces
	# Ugh - get rid of any pending after events. If we do not do this, the
	# regrid function may be called and this caused the tk grid command
	# to crash!
	foreach i [set $w.after_id] {
	    catch {after cancel $i} err
	}
	set $w.after_id ""
	update idletasks
	return
    }

    incr $w.NTraces -1
    set l {}
    foreach i [set $w.Order] {
	if {$i != $f} {
	    lappend l $i
	}
    }
    set $w.Order $l
    regsub -all {\.[^#]*} [set $w.Order] {} x

    append $w.after_id [after idle "trace_regrid $w $w.traces.c $w.traces.s"] " "
}


#
# Quits the trace display and unsets any global variables used
#
proc trace_display_quit {w e} {
    global $w.after_id

    bind $w <Destroy> {}
    if {[winfo exists $w]} {destroy $w}
    foreach i [set $w.after_id] {
	catch {after cancel $i} err
    }
    set $w.after_id ""

    global trace_columns
    catch {trace vdelete trace_columns w "trace_columns_trace $w"} err
    catch {trace vdelete trace_rows w "trace_rows_trace $w"} err
}


#
# Scale the y magnification of trace 't' to a value 'val'/10. Called from the
# ymag scale.
#
proc trace_yscale {t val} {
    $t ymag [expr double($val) / 10]
}

#
# Scale the x magnification of trace 't' to a value 'val'. Called from the
# xmag scale. This sets the global definition too so that the value is
# permanent (for the lifetime of the application).
#
proc trace_xscale {t val} {
    global gap5_defs
    $t xmag $val
    keylset gap5_defs TRACE_DISPLAY.XMAG $val
}

#
# Displays the information about a trace 't' for tracedisplay 'w', editor 'e'
# 
proc trace_display_info {t w e title} {
    # If it already exists then reset the text
    if {[xtoplevel $w.t_info -resizable 0] == ""} {
	$w.t_info.text configure -text [$t info]
        wm title $w.t_info "Info: $title"
	return
    }

    wm title $w.t_info "Info: $title"

    message $w.t_info.text -text "[$t info]" -aspect 200 -bd 1 -relief groove
    button $w.t_info.quit -text "Close" -command "destroy $w.t_info"

    pack $w.t_info.text -side top -fill both -expand 1
    pack $w.t_info.quit -side bottom
}

#
# Sets difference mode where two traces can be picked to display the
# trace differences.
# WARNING: This is all a bit hideous; especially the cursor bindings. If only
# we could use bindtags for cursors too.
#
proc trace_display_diff {but t w e} {
    $but configure -cursor iron_cross
    [winfo parent $but] configure -cursor iron_cross
    bind $but <Destroy> "trace_display_diff_cancel $but %W"

    bind tdiff <<select>> "trace_display_diff_callback $but %W $t $e; break"
    bind tdiff <<not-select>> "trace_display_diff_cancel $but %W; break"
    bind tdiff <Any-Enter> {
	catch {%W configure -cursor iron_cross}
	for {set w %W} {$w != "." && [winfo class $w] != "TraceRow"} \
	    {set w [winfo parent $w]} {}
	if {$w != "."} {
	    global $w.Diff
	    set $w.Diff blue
	    trace_highlight $w
	    unset $w.Diff
	}
	break}
    bind tdiff <Any-Leave> {catch {%W configure -cursor top_left_arrow}; break}
}

proc trace_display_diff_callback {but w t e} {
    if {[regsub {(.*\.traces.c.f.[^.]*)\..*} $w {\1.trace} t2] == 1} {
        $e diff_trace $t $t2
    } else {
	bell
    }

    trace_display_diff_cancel $but $w
}

proc trace_display_diff_cancel {but w} {
    bind tdiff <<select>> {}
    bind tdiff <<not-select>> {}
    bind tdiff <Any-Enter> {}
    bind tdiff <Any-Leave> {}
    catch {$w configure -cursor top_left_arrow} z
    catch {[winfo parent $w] configure -cursor top_left_arrow} z
}

#
# GUI for the 'trace differences against specific sequence' mode.
#
proc trace_diff_specific {w e} {
    set t $w.diff_spec
    if {[xtoplevel $t -resizable 0] == ""} {return}
    wm title $t "Diff against specific trace"
    wm resizable $t 0 0

    entrybox $t.name \
	-title "Compare which reading?" \
	-type "CheckEditorGel $e" \
	-command "trace_diff_specific2 $e $t"

    okcancelhelp $t.ok \
        -ok_command "trace_diff_specific2 $e $t \[$t.name.entry get\]"\
        -cancel_command "destroy $t" \
	-help_command "show_help gap5 {Editor-Trace Display}" \
        -bd 2 -relief groove

    pack $t.name $t.ok -side top -fill x
}

proc trace_diff_specific2 {e w name} {
    set num [$e find_read $name]
    if {$num == -1} {bell;return;}

    destroy $w
    $e trace_comparator $name
}

proc CheckEditorGel {e path} {
    set n 
    if {"[$e find_read [$path.entry get]]" == -1} {
	return 1
    } else {
	return 0
    }
}


#
# Does the actual work of creating a trace display; the common component
# of trace_add and trace_create.
#
proc new_trace_create2 {w e title allow_diff} {
    global $w.Incr
    global gap5_defs tk_utils_defs db_namelen
    global $w.WinPos
    global $w.NTraces
    global $w.Order
    global $w.Trace_width
    global $w.Highlight $w.HighlightOld $w.after_id
    global tcl_platform
    global trace_height trace_columns trace_rows
    global trace_confidence
    global trace_display_width
    global $w.HasConf
    global trace_compact_mode $w.Trace_Compact_Mode
    global $w.Trace_Label_Width


    set dt [keylget gap5_defs TRACE_DISPLAY]

    # Create master trace window
    if {![winfo exists $w]} {
	xtoplevel $w -resizable 0 -focus 0
	wm title $w "Trace display"
	if {[info exists $w.WinPos]} {
	    wm geometry $w [set $w.WinPos]
	}

        frame $w.bar -bd 2 -relief raised
	frame $w.traces -bd 2 -relief raised
	frame $w.info
	set $w.NTraces 0
	set $w.Order {}
	set $w.Highlight {}
	set $w.HighlightOld [$w cget -bg]
	set $w.Incr 0
	set $w.Trace_Compact_Mode $trace_compact_mode
	trace_set_compact $w

	# Canvas with scrollbar

	canvas $w.traces.c -bd 0 -yscrollcommand "$w.traces.s set" \
	    -yscrollincrement $trace_height -takefocus 0 \
	    -width $trace_display_width \
	     -highlightthickness 0
	frame $w.traces.c.f -bd 1 -takefocus 0 -highlightthickness 0
	$w.traces.c create window 0 0 -window $w.traces.c.f -anchor nw
	scrollbar $w.traces.s -command "$w.traces.c yview" -orient vert
	bind $w <Any-Configure> {trace_window_moved %W %x %y}

	# Create menu bar
	checkbutton $w.bar.lock \
	     -text Lock -variable $w.Lock \
	     -command "$e set_trace_lock"
	global $w.Lock
       	set $w.Lock [keylget gap5_defs CONTIG_EDITOR.TRACE_LOCK]
	pack $w.bar.lock -side left -padx 2m

	label $w.bar.cols -text Columns:
	label $w.bar.c1 -text 1 -relief raised
	label $w.bar.c2 -text 2 -relief raised
	label $w.bar.c3 -text 3 -relief raised
	label $w.bar.c4 -text 4 -relief raised
	bind $w.bar.c1 <1> {set trace_columns 1}
	bind $w.bar.c2 <1> {set trace_columns 2}
	bind $w.bar.c3 <1> {set trace_columns 3}
	bind $w.bar.c4 <1> {set trace_columns 4}
	pack $w.bar.cols $w.bar.c1 $w.bar.c2 $w.bar.c3 $w.bar.c4 -side left

	label $w.bar.rows -text "   Rows:"
	label $w.bar.r1 -text 1 -relief raised
	label $w.bar.r2 -text 2 -relief raised
	label $w.bar.r3 -text 3 -relief raised
	label $w.bar.r4 -text 4 -relief raised
	label $w.bar.r5 -text 5 -relief raised
	label $w.bar.r6 -text 6 -relief raised
	bind $w.bar.r1 <1> {set trace_rows 1}
	bind $w.bar.r2 <1> {set trace_rows 2}
	bind $w.bar.r3 <1> {set trace_rows 3}
	bind $w.bar.r4 <1> {set trace_rows 4}
	bind $w.bar.r5 <1> {set trace_rows 5}
	bind $w.bar.r6 <1> {set trace_rows 6}
	pack $w.bar.rows $w.bar.r1 $w.bar.r2 $w.bar.r3 $w.bar.r4 $w.bar.r5 \
	    $w.bar.r6 -side left

	# Show confidence
	checkbutton $w.bar.conf -text "Show confidence" \
	    -variable trace_confidence \
	    -command "trace_show_confidence $w"
	pack $w.bar.conf -side left -padx 2m

       	# Save configuration
	button $w.bar.save -text "Save settings" \
	    -command "trace_save_layout"
	pack $w.bar.save -side left -padx 2m

	# Compact mode - disable for now
#	checkbutton $w.bar.compact -text "Compact" \
#	    -variable $w.Trace_Compact_Mode \
#	    -command "trace_set_compact $w"
#	pack $w.bar.compact -side left -padx 2m

	# Help menu
	button $w.bar.quit \
	     -text Close \
	     -command "trace_display_quit $w $e"
	xmenubutton $w.bar.help \
	    -text "Help >>" \
	    -menu $w.bar.help.m -padx 2
	menu $w.bar.help.m
	$w.bar.help.m add command \
	    -label "Trace Display" \
	    -command "show_help gap5 {Editor-Traces}"
	$w.bar.help.m add command \
	    -label "Settings" \
	    -command "show_help gap5 {Editor-Trace Display}"
	$w.bar.help.m add command \
	    -label "Consensus Trace" \
	    -command "show_help gap5 {Editor-Comm-Consensus Trace}"
	pack $w.bar.help $w.bar.quit -side right -fill both
	
	# Information line
	label $w.info.dummy
	label $w.info.l
	pack $w.info -side bottom -fill both
	pack $w.info.dummy -fill both
	place $w.info.l -relx 0

	pack $w.bar -side top -fill x
	pack $w.traces   -side bottom -fill both -expand 1
	grid $w.traces.c $w.traces.s -sticky nsew
	grid columnconfigure $w.traces 0 -weight 1
	grid rowconfigure $w.traces 0 -weight 1

	trace variable trace_columns w "trace_columns_trace $w"
	trace variable trace_rows w "trace_rows_trace $w"

	bind $w <Destroy> "if {\"%W\" == \"$w\"} {trace_display_quit $w $e}; break"
	reshape_trace_window $w $w.traces.c $w.traces.s $trace_columns
	set trace_rows $trace_rows
    }

    # Add an individual trace to the window
    regsub { *$} $title {} title
    set ymag 10
    if {[string match "diff*" $title]} {
        set ymag 22
	regexp {diffs: (.)([^ ]*) \#([^ ]*)} $title dummy type t1 t2
	if {$type == "="} {
	    set type ""
	    set t1 C
	    set title "Difference between Consenus and trace [r_name [$e io] $t2]"
	} else {
	    set title "Difference between trace [r_name [$e io] $t1] and trace [r_name [$e io] $t2]"
	}
	set tit_num "$type$t1/#$t2"
	set tname trace${t1}_${t2}
	set tmp $tname
	set copy 0
	while {[winfo exists $w.traces.c.f.${tname}_r]} {
	    incr copy
	    set tname ${tmp}#$copy
	}
	if {$t1 == "C"} {
	    set tit_name "consensus/#$t2"
	} else {
	    set tit_name ""
	}
	set info_state disabled
        set tf_r $w.traces.c.f.${tname}_r[set $w.Incr]
	set $w.HasConf($tf_r) 0
    } elseif {[string match "consensus*" $title]} {
	set tit_num $title
	set tname tracecons
	set tit_name "consensus"
	set info_state disabled
        set tf_r $w.traces.c.f.${tname}_r[set $w.Incr]
	set $w.HasConf($tf_r) 0
    } else {
        set tname [string tolower $title]
	regsub -all {\.} $tname _ tname
	set tit_name ""
	set tit_num $title
	set info_state normal
	set tf_r $w.traces.c.f.${tname}_r[set $w.Incr]
	set $w.HasConf($tf_r) 1
    }
    incr $w.Incr

    set tnum [set $w.NTraces]
    incr $w.NTraces
    if {[keylget gap5_defs TRACE_DISPLAY.DIRECTION] == "top"} {
	set tnum [expr 9983-$tnum]
    }

    grid rowconfigure $w.traces.c.f [expr $tnum/$trace_columns] -weight 1
    if {[keylget gap5_defs TRACE_DISPLAY.DIRECTION] == "top"} {
        set $w.Order "$tf_r [set $w.Order]"
    } else {
        lappend $w.Order $tf_r
    }
    regsub -all {\.[^\#]*} [set $w.Order] {} x

    # Main frame, with <Destroy> binding.
    frame $tf_r -bd 2 -relief sunken -highlightthickness 2 -class TraceRow
    set tf_l [frame $tf_r.l -width [set $w.Trace_Label_Width]]
    grid propagate $tf_l 0
    bind $tf_r <Destroy> "trace_remove $tf_r $w $e"

    # Labels and buttons / menu
    label $tf_l.spacer -text ""
    menu $tf_r.m
    frame $tf_l.labels
    frame $tf_l.labels.b
    button $tf_l.labels.b.info -text "Info" -pady 2 \
	-command "trace_display_info $tf_r.trace $w $e {$title}" \
        -state $info_state
    $tf_r.m add command \
	-label "Information" \
	-command "trace_display_info $tf_r.trace $w $e {$title}" \
	-state $info_state
    # Disable difference trace for now as not implemented in gap5 yet
    if {0 && $allow_diff} {
	button $tf_l.labels.b.diff -text "Diff" -pady 2 \
	    -command "trace_display_diff $tf_l.labels.b.diff $tf_r.trace $w $e"
	$tf_r.m add command \
	    -label "Trace difference" \
	    -command "trace_display_diff $tf_l.labels.b.diff $tf_r.trace $w $e"
    } else {
	button $tf_l.labels.b.diff -text "Save" -pady 2 \
	    -command "trace_save $tf_r.trace"
	$tf_r.m add command \
	    -label "Save" \
	    -command "trace_save $tf_r.trace"
    }
    button $tf_l.labels.b.complement -text "Comp." -pady 2 \
	-command "trace_complement $e $tf_r.trace"
    $tf_r.m add command \
	-label "Complement" \
	-command "trace_complement $e $tf_r.trace"
    button $tf_l.labels.b.quit -text "Cancel" -pady 2 \
	-command "destroy $tf_r"
    $tf_r.m add command \
	-label "Quit" \
	-command "destroy $tf_r"

    pack $tf_l.labels.b -side left -fill both
    pack $tf_l.labels.b.info $tf_l.labels.b.diff \
	$tf_l.labels.b.complement $tf_l.labels.b.quit -side top -fill x
    bind $tf_r <Any-Enter> "$w.info.l configure -text [list $title]"

    # Trace + scrollbar
    if {[set $w.HasConf($tf_r)]} {
        set tc $trace_confidence
    } else {
        set tc 0
    }

    # Height minus borderwidth and highlightthickess
    set remaining_height [expr {$trace_height-2*(2+2+2)}]
	
    dnatrace $tf_r.trace \
	-bd 2 \
	-relief ridge \
	-xscrollcommand "$tf_r.scroll set" \
	-showedits 0 \
	-showconf 	$tc \
	-bg 	        [keylget dt BACKGROUND] \
	-xmag	        [expr double([keylget dt XMAG])/100] \
	-width		[set $w.Trace_width] \
	-height 	$remaining_height \
	-line_width     [keylget gap5_defs TRACE_DISPLAY.LINE_WIDTH]
    if {$allow_diff} {
	recurse_trace_bind $tf_r
    }
    bind $tf_r.trace <<menu>> "tk_popup $tf_r.m %X %Y"

    if {$tcl_platform(platform) == "unix"} {
	scrollbar $tf_r.scroll \
		-command "trace_scroll $e $tf_r.trace" \
		-orient horizontal \
		-width 10 \
		-bd 2 \
		-highlightthickness 0
    } else {
	scrollbar $tf_r.scroll \
		-command "trace_scroll $e $tf_r.trace" \
		-orient horizontal \
		-bd 2 \
		-highlightthickness 0
    }

    # X and Y scales
    frame $tf_l.xmag
    label $tf_l.xmag.l -text "X"
    scale $tf_l.xmag.s -from 10 -to 1000 -showvalue 0 \
	-command "trace_xscale $tf_r.trace"
    $tf_l.xmag.s set [keylget dt XMAG]
    frame $tf_l.ymag
    label $tf_l.ymag.l -text "Y"
    scale $tf_l.ymag.s -from 10 -to 100 -showvalue 0 \
	-command "trace_yscale $tf_r.trace"
    $tf_l.ymag.s set $ymag
    if {[set $w.Trace_Compact_Mode]} {
	$tf_l.xmag.s configure -width 12
	$tf_l.ymag.s configure -width 12
    } else {
	$tf_l.xmag.s configure -width 15
	$tf_l.ymag.s configure -width 15
    }
    pack $tf_l.xmag.l $tf_l.ymag.l -side top -fill x
    pack $tf_l.xmag.s $tf_l.ymag.s -side bottom -fill both -expand 1

    grid columnconfigure $tf_r 0 -weight 1
    grid rowconfigure $tf_r 0 -weight 1

    grid $tf_l.spacer -column 0 -row 0 -sticky nsew
    if {![set $w.Trace_Compact_Mode]} {
	grid $tf_l.labels -column 0 -row 1 -sticky nsew
	grid $tf_l.xmag   -column 1 -row 1 -sticky nsew
	grid $tf_l.ymag   -column 2 -row 1 -sticky nsew
    } else {
	grid $tf_l.xmag   -column 0 -row 1 -sticky nsew
	grid $tf_l.ymag   -column 1 -row 1 -sticky nsew
    }

    grid $tf_r.l      -column 0 -row 0 -sticky nsew -rowspan 2
    grid $tf_r.trace  -column 1 -row 0 -sticky nsew
    grid $tf_r.scroll -column 1 -row 1 -sticky nsew

    grid $tf_r -row [expr {$tnum/$trace_columns}] \
	-column [expr {$tnum%$trace_columns}] -sticky nsew

    trace_highlight $tf_r

    global $w.PendingReconfigure
    if {![info exists $w.PendingReconfigure]} {
	append $w.after_id [after idle "trace_reconfigure $w $w.traces.c $w.traces.s \[set $w.PendingReconfigure\]; unset $w.PendingReconfigure"] " "
    }
    set $w.PendingReconfigure $tf_r


    if {[keylget gap5_defs TRACE_DISPLAY.FULL_NAME]} {
	if {$tit_name != ""} {
	    set tit_name " / $tit_name"
	}
	label $tf_r.name -text "$tit_num$tit_name"
        place $tf_r.name -relx 0
        raise $tf_r.name
	frame $tf_r.trace.flash \
	    -height [expr {$trace_height-52}] \
	    -width 0 \
            -bg [keylget tk_utils_defs TRACE.COLOUR_CURSOR]

        raise $tf_r.trace.flash
    }

    bind $tf_r.trace <Any-Motion> \
	"$w.info.l configure -text \[%W base_info @%x\]"

    # Mouse-wheel support
    focus $tf_r.trace
    bind $tf_r.trace <MouseWheel> \
	"$w.traces.c yview scroll \[expr {-(%D)/120}\] units"
    if {[tk windowingsystem] eq "x11"} {
	foreach c {trace name} {
	    bind $tf_r.$c <4> "$w.traces.c yview scroll -1 units"
	    bind $tf_r.$c <5> "$w.traces.c yview scroll +1 units"
	    bind $tf_r.$c <Control-4> \
		"$w.traces.c yview scroll -\$trace_rows units"
	    bind $tf_r.$c <Control-5> \
		"$w.traces.c yview scroll +\$trace_rows units"
	}
    }


    # Without this use "auto-diff traces" gives blank windows. Everything
    # has been gridded and sized correctly, but tk seems to have calculated
    # the geometry information wrongly. At this stage (when idle) changing
    # the frame border width forces a recomputation of the geometry and the
    # traces magically reappear. Hence our "after idle" request.
    #
    # It should be noted that this bug does not occur for other ways of
    # bringing up traces. I believe it is due to the way this is done in the
    # C code for auto-diff, which invokes all of the traces and then sets the
    # number of trace columns, without doing an intervening update. Adding an
    # update idletasks (in the C trace_columns function) does correct the
    # missing traces bug, but bizarrely causes another bug where the number of
    # columns does not change.
    after idle "catch {$w.traces.c.f configure -bd 0}"

    reshape_trace_window $w $w.traces.c $w.traces.s $trace_columns

    return $tf_r
}

proc trace_reconfigure {w c ys {show_trace {}}} {
    global gap5_defs
    global trace_height $w.NTraces
    global trace_columns trace_rows trace_display_width

    if {$trace_rows < 1} {
	set trace_rows 1
    }
    set trace_display_height [expr $trace_rows * $trace_height]

    set h [expr {$trace_height * [expr {([set $w.NTraces]+($trace_columns-1))/$trace_columns}]}]
    $c itemconfigure 1 -height $h
    $c.f configure -height $h

    $c configure -scrollregion "0 0 [winfo width $c.f] $h"
    if {$h >= $trace_display_height} {set h $trace_display_height}
    $c configure -height $h -width $trace_display_width

    if {$show_trace != {}} {
	trace_visible $w $ys $show_trace
    }
}

proc trace_regrid {w c ys} {
    global $w.Order trace_columns gap5_defs $w.Trace_Label_Width

    foreach r [set $w.Order] {
	grid forget $r
    }

    if {[keylget gap5_defs TRACE_DISPLAY.DIRECTION] == "top"} {
	set n [expr 9983+[llength [set $w.Order]]]
	set inc 1
    } else {
        set n 0
        set inc 1
    }

    foreach r [set $w.Order] {
        grid $r -row [expr {$n/$trace_columns}] \
	    -column [expr {$n%$trace_columns}] -sticky nsew
	incr n $inc
    }

    trace_reconfigure $w $c $ys
}

proc trace_scroll {e tr args} {
    $e trace_scroll $tr $args
}


proc reshape_trace_window {w c ys cols} {
    global $w.Trace_width trace_display_width trace_columns
    global $w.Order $w.Trace_Label_Width

    $w.bar.c1 configure -relief raised
    $w.bar.c2 configure -relief raised
    $w.bar.c3 configure -relief raised
    $w.bar.c4 configure -relief raised
    $w.bar.c$cols configure -relief sunken

    set trace_columns $cols
    # 8 = 2*(tf_r borderwidth + highlightthickness)
    # trace_display_width is canvas width, which has size 0 borders.
    set $w.Trace_width [expr ($trace_display_width)/$cols - \
			    [set $w.Trace_Label_Width] - 8]

    if {[set $w.Order] == {}} {
	return
    }

    foreach r [set $w.Order] {
	$r.trace configure -width [set $w.Trace_width]
    }
    trace_regrid $w $c $ys
}

proc trace_changed_rows {w c ys args} {
    global trace_rows
    for {set i 1} {$i <= 6} {incr i} {
	$w.bar.r$i configure -relief raised
    }
    $w.bar.r$trace_rows configure -relief sunken
    trace_reconfigure $w $c $ys
}

proc trace_window_moved {w x y} {
    if {$w != [winfo toplevel $w]} { return }

    global $w.WinPos
    # set $w.WinPos +$x+$y
    regsub {^[^+]*\+} [wm geometry $w] {+} z
    set $w.WinPos $z
}

proc trace_visible {w ys trace} {
    global gap5_defs $w.Order $w.NTraces trace_columns trace_rows

    regsub -all {\.[^#]*} [set $w.Order] {} x
    update idletasks

    # See if the trace is visible
    set st [expr round([lindex [$ys get] 0]*[set $w.NTraces])]
    set en [expr round([lindex [$ys get] 1]*[set $w.NTraces])]
    set vis 0
    set order [set $w.Order]
    for {set i $st} {$i < $en} {incr i} {
	if {"[lindex $order $i]" == "$trace"} {
	     set vis 1
	     break
	}
    }
    if {$vis} return

    # Find trace index
    set en [set $w.NTraces]
    for {set i 0} {$i < $en} {incr i} {
	if {"[lindex $order $i]" == "$trace"} {
	     break
	}
    }
    if {$i == $en} {
	puts "Trace not in display!?"
	return
    }

    # Make this index visible
    if {$i > [expr [set $w.NTraces]-$trace_rows]} {
	set i [expr [set $w.NTraces]-$trace_rows]
    }

    set i [expr $i-($i%$trace_columns)]
    set cmd [$ys cget -command]
    uplevel #0 $cmd moveto [expr $i.0/[set $w.NTraces]]

    return
}

# Highlights a specific trace. Handy for keep track of which trace was
# last displayed
proc trace_highlight {trace {make_visible 0}} {
    set w [winfo toplevel $trace]
    global $w.Highlight $w.HighlightOld

    if {[set $w.Highlight] != "" && [winfo exists [set $w.Highlight]]} {
	[set $w.Highlight] configure -highlightbackground [set $w.HighlightOld]
    }

    if {$trace != ""} {
	global $trace.Diff
	if {[info exist $trace.Diff]} {
	    set colour [set $trace.Diff]
	} else {
	    set colour red
	}
        $trace configure -highlightbackground $colour
	set $w.Highlight $trace
	if {$make_visible} {
	    trace_visible $w $w.traces.s $trace
	    trace_flash $trace.trace.flash
	}
    }
}

proc trace_complement {e t} {
     $t complement
}

proc trace_save {t} {
    set fname [tk_getSaveFile -parent $t]
    if {$fname == ""} {
	return
    }

    catch {$t save $fname ZTR} err
    if {$err != ""} {
	tk_messageBox -type ok -icon error -parent $t \
	    -message "Failed to save trace file"
    }
}

proc trace_save_layout {} {
    global trace_rows trace_columns trace_display_width gap5_defs env
    global trace_confidence trace_compact_mode

    keylset gap5_defs TRACE_DISPLAY.ROWS [list $trace_rows]
    keylset gap5_defs TRACE_DISPLAY.COLUMNS $trace_columns
    keylset gap5_defs TRACE_DISPLAY.WINDOW_WIDTH $trace_display_width
    keylset gap5_defs TRACE_DISPLAY.DIRECTION \
	[keylget gap5_defs TRACE_DISPLAY.DIRECTION]
    keylset gap5_defs TRACE_DISPLAY.SHOW_CONFIDENCE $trace_confidence
    keylset gap5_defs TRACE_DISPLAY.COMPACT_MODE $trace_compact_mode

    update_defs gap5_defs $env(HOME)/.gap5rc \
	TRACE_DISPLAY.ROWS \
	TRACE_DISPLAY.COLUMNS \
	TRACE_DISPLAY.WINDOW_WIDTH \
	TRACE_DISPLAY.DIRECTION \
	TRACE_DISPLAY.SHOW_CONFIDENCE \
	TRACE_DISPLAY.COMPACT_MODE
}

proc trace_flash {w {ind 5} {wid 20}} {
    place $w -relx 0.505 -y 0 -anchor n
    for {set t 0} {$t < $ind} {incr t} {
	after 30
	if {[catch {$w configure -width $wid} err]} {
	    return
        }
        set wid [expr $wid-5]
	update idletasks
    }
    place forget $w
}   


proc trace_show_confidence {w} {
    global trace_confidence $w.Order $w.HasConf

    foreach r [set $w.Order] {
	if {[set $w.HasConf($r)]} {
	    $r.trace configure -showconf $trace_confidence
	} else {
	    $r.trace configure -showconf 0
	}
    }
}

proc trace_set_compact {w} {
    global $w.Trace_Compact_Mode $w.Trace_Label_Width trace_compact_mode
    global $w.Order

    set trace_compact_mode [set $w.Trace_Compact_Mode]
    if {[set $w.Trace_Compact_Mode]} {
	set $w.Trace_Label_Width 44
    } else {
	set $w.Trace_Label_Width 123
    }

    foreach t [set $w.Order] {
	$t.l configure -width [set $w.Trace_Label_Width]
	if {[set $w.Trace_Compact_Mode]} {
	    $t.l.xmag.s configure -width 12
	    $t.l.ymag.s configure -width 12
	    grid forget $t.l.labels
	    grid $t.l.xmag    -column 0 -row 1 -sticky nsew
	    grid $t.l.ymag    -column 1 -row 1 -sticky nsew
	} else {
	    $t.l.xmag.s configure -width 15
	    $t.l.ymag.s configure -width 15
	    grid $t.l.labels  -column 0 -row 1 -sticky nsew
	    grid $t.l.xmag    -column 1 -row 1 -sticky nsew
	    grid $t.l.ymag    -column 2 -row 1 -sticky nsew
	}
    }

    if {[set $w.Order] != ""} {
	global trace_columns
	reshape_trace_window $w $w.traces.c $w.traces.s $trace_columns
    }
}

proc trace_columns_trace {w args} {
    global trace_columns
    if {[winfo exists $w]} {
	reshape_trace_window $w $w.traces.c $w.traces.s $trace_columns
    } else {
	trace vdelete trace_columns w "trace_columns_trace $w"
    }
}

proc trace_rows_trace {w args} {
    if {[winfo exists $w]} {
	trace_changed_rows $w $w.traces.c $w.traces.s
    } else {
	trace vdelete trace_rows w "trace_rows_trace $w"
    }
}


# -----------------------------------------------------------------------------
# Height is in lines of text
proc trace_small_add {w file e seqnum height} {
    array set f [font metrics sheet_font]
    # trace_small_add .cedit0.0.seqs.traces xb54f3.s1ta.ztr .cedit0.0.seqs 12
    set t $e.trace_$seqnum
    dnatrace $t \
	-height [expr {$height * $f(-linespace)}] \
	-showsequence 0 \
	-shownumbers 0 \
	-showedits 0 \
	-showconf 0 \
	-xmag 1 \
	-width [winfo width $e] \
	-background white \
	-bd 0 \
	-showends 1

    $t load $file
    $t resample [font measure sheet_font A]

    bind $t <<trace>> " 
        set relx \[expr {%x+\[winfo rootx $t\]-\[winfo rootx $e\]}\] 
        set rely \[expr {%y+\[winfo rooty $t\]-\[winfo rooty $e\]}\] 
        if {\[$e cursor_set @\$relx @\$rely\] == 0} { 
            $e invoke_trace 
        } 
    " 

    return $t
}

