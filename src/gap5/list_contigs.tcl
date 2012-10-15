#
# Implements the "List Contigs" panel - a textual alternative to the contig
# selector window.
#

##############################################################################
# Main entry point - produces a tablelist containing a list of contig details.
# This may be used for selecting and displaying the "contigs" list.
#
# io		a GapIO handle
# parent	from which window should this be a child
# csh_win	unused (and optional) - will be removed later

proc InitListContigs {io parent {csh_win {}}} {
    global NGList

    set t $parent.list_contigs
    if {[xtoplevel $t] == ""} return
    wm title $t "Contig List"

    # Create our tablelist
    if {[$io db_version] >= 6} {
	tablelist $t.list \
	    -columns {0 "Name" 10 "Length" 10 "Clipped len" 14 "# sequences" 14 "# annotations" 20 "Scaffold"} \
	    -labelcommand tablelist::sortByColumn \
	    -exportselection 0 \
	    -stretch 0 \
	    -selectmode extended \
	    -yscrollcommand [list $t.yscroll set]
	$t.list columnconfigure 0 -sortmode dictionary
	$t.list columnconfigure 1 -sortmode integer
	$t.list columnconfigure 2 -sortmode integer
	$t.list columnconfigure 3 -sortmode integer
	$t.list columnconfigure 4 -sortmode integer
	$t.list columnconfigure 5 -sortmode command \
	    -formatcommand ListContigsScaffoldFormat \
	    -sortcommand [list ListContigsScaffoldSort $t.list]
	#$t.list columnconfigure 6 -sortmode integer
    } elseif {[$io db_version] >= 5} {
	tablelist $t.list \
	    -columns {30 "Name" 10 "Length" 14 "# sequences" 14 "# annotations" 20 "Scaffold"} \
	    -labelcommand tablelist::sortByColumn \
	    -exportselection 0 \
	    -stretch 0 \
	    -selectmode extended \
	    -yscrollcommand [list $t.yscroll set]
	$t.list columnconfigure 0 -sortmode dictionary
	$t.list columnconfigure 1 -sortmode integer
	$t.list columnconfigure 2 -sortmode integer
	$t.list columnconfigure 3 -sortmode integer
	$t.list columnconfigure 4 -sortmode command \
	    -formatcommand ListContigsScaffoldFormat \
	    -sortcommand [list ListContigsScaffoldSort $t.list]
    } else {
	tablelist $t.list \
	    -columns {30 "Name" 10 "Length" 14 "# sequences" 14 "# annotations"} \
	    -labelcommand tablelist::sortByColumn \
	    -exportselection 0 \
	    -stretch 0 \
	    -selectmode extended \
	    -yscrollcommand [list $t.yscroll set]
	$t.list columnconfigure 0 -sortmode dictionary
	$t.list columnconfigure 1 -sortmode integer
	$t.list columnconfigure 2 -sortmode integer
	$t.list columnconfigure 3 -sortmode integer
    }
    
    frame $t.buttons -bd 0
    button $t.buttons.cancel \
	-text "Cancel" \
	-command "destroy $t"

    button $t.buttons.save \
    	-text "Save order" \
	-command "ListContigsSave $io $t.list"
    if {[$io read_only]} {
	$t.buttons.save configure -state disabled
    }

    button $t.buttons.copy \
    	-text Copy \
	-command "ListContigsCopy $t.list"

    button $t.buttons.help \
	-text "Help" \
	-command "show_help gap5 Contig-Selector-Contigs"

    repeater $t.buttons.up "ListContigsUp $io $t.list" \
	-text "\u2191" \

    repeater $t.buttons.down "ListContigsDown $io $t.list" \
	-text "\u2193" 

    pack $t.buttons.cancel $t.buttons.save $t.buttons.copy $t.buttons.help \
	$t.buttons.up $t.buttons.down -side left -expand 1

    # Add a scrollbar    
    scrollbar $t.yscroll -command "$t.list yview"

    grid columnconfigure $t 0 -weight 1
    grid rowconfigure $t 0 -weight 1

    grid $t.list $t.yscroll -sticky nsew
    grid $t.buttons -sticky nsew

    bind [$t.list bodypath] <<menu>> \
	"ListContigsMenuBinding $io $t.list %x %y %X %Y"
    bind [$t.list bodypath] <<select>> \
	"+ListContigsSelectPressBinding $io %W %x %y"
    bind [$t.list bodypath] <<select-release>> \
	"+ListContigsSelectReleaseBinding $io $t.list"
    bind [$t.list bodypath] <Control-Key-a> \
	"ListContigsSelectAll $io $t.list"

    focus [$t.list bodypath]

    wm geometry $t 800x300

    set trace_cmd "ListContigsUpdate $io $t.list"
    trace variable NGList(contigs) w $trace_cmd
    bind $t <Destroy> \
	"ListContigsExit $io $t [list $trace_cmd]"
    wm protocol $t WM_DELETE_WINDOW \
	"ListContigsExit $io $t [list $trace_cmd]"

    # Register with all contigs so we can update when contigs are joined
    # or broken apart. Contig 0 implies all contigs.
    global $t.Reg
    set $t.Reg [contig_register \
		    -io $io \
		    -contig 0 \
		    -command "ListContigsCallback $io $t.list" \
		    -flags {REQUIRED LENGTH JOIN_TO DELETE COMPLEMENT RENAME ORDER}]
    
    # Populate the list
    ListContigsRepopulate $io $t.list
}

##############################################################################
# Formatting of scaffold name & index into just scaffold name
; proc ListContigsScaffoldFormat {name} {
    return $name
    #return [lindex [regexp -inline {(.*)/(-?[0-9]+)} $name] 1]
}

##############################################################################
# Sort function for scaffolds. If we reverse scaffolds we want to keep the
# members numerically sorted even though the scaffolds themselves are
# reverse sorted.
; proc ListContigsScaffoldSort {w n1 n2} {
    regexp {(.*)/(-?[0-9]+)} $n1 _ s1 c1
    regexp {(.*)/(-?[0-9]+)} $n2 _ s2 c2
    if {$s1 != $s2} {
	if {[lsort -dictionary [list $s1 $s2]] == "$s1 $s2"} {
	    return -1
	} else {
	    return 1
	}
    }

    # Identical scaffolds, so sort by index into scaffold. We have to
    # reverse this if we're doing a reverse search though. To get the
    # sort order we sneakily peek into the internal structures for the
    # tablelist widget via  {::tablelist::ns${w}::data}(4-sortOrder)
    set order [set ::tablelist::ns${w}::data(4-sortOrder)]
    if {$order == "decreasing"} {
	return [expr {$c2 < $c1 ? -1 : ($c2 > $c1 ? 1 : 0)}]
    } else {
	return [expr {$c1 < $c2 ? -1 : ($c1 > $c2 ? 1 : 0)}]
    }
}

##############################################################################
# Creates a menu for the list contigs window. This is a duplicate of the
# contig selector menu, but the contig selector menu code requires knowledge
# of canvas ids and canvas tag values right down to the final menu callbacks,
# which makes it nigh on impossible to use from something other than the
# csh_win canvas.
#
# FIXME: rewrite the csh_win canvas stuff to call this menu code instead.
; proc ListContigsPopupMenu {io crec w x y X Y} {
    global read_only gap5_defs

    create_popup $w.m "Contig Commands (\#$crec)"
    $w.m add command -label "Edit contig" \
	-command "edit_contig -io $io -contig $crec"
    $w.m add command -label "Template Display" \
	-command "CreateTemplateDisplay $io $crec"
    if {!$read_only} {
	$w.m add command -label "Complement contig" \
	    -command "complement_contig -io $io -contigs =$crec"
	set scaf [lindex [$w get @$x,$y] 4]
	set scaf [lindex [regexp -inline {(.*)/.*} $scaf] 1]
	if {$scaf != "(none)"} {
	    $w.m add command -label "Complement scaffold" \
		-command "complement_scaffold -io $io -scaffolds [list $scaf]"
	}
	$w.m add separator
	$w.m add command -label "Rename contig" \
	    -command "ListContigsRename $w $io $crec"
	if {[$io db_version] >= 5} {
	    $w.m add command -label "Change scaffold" \
		-command "ListContigsChangeScaffold $w $io $crec"
	}
    }
#    $w.m add command -label "List notes" \
#	-command "NoteSelector $io contig $crec"
    tk_popup $w.m [expr $X-20] [expr $Y-10]
}

##############################################################################
# Rename Contig callback
; proc ListContigsRename {t io crec} {
    set w [toplevel $t.rename]
    set c [$io get_contig $crec]

    wm title $w "Rename Contig: [$c get_name]"
    
    xentry $w.name \
	-label "New contig name" \
	-default [$c get_name]
    $c delete
    
    okcancelhelp $w.ok \
	-ok_command "contig_rename $io $crec \[$w.name get\] $w; ListContigsRepopulate $io $t" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap5 Rename-Contig"

    pack $w.name $w.ok -side top -fill both
}

##############################################################################
# Change Scaffold callback
; proc ListContigsChangeScaffold {t io crec} {
    set w [xtoplevel $t.scaffold]

    set c [$io get_contig $crec]
    wm title $w "Change Scaffold: [$c get_name]"
    
    if {[$c get_scaffold] > 0} {
	set f [$io get_scaffold [$c get_scaffold]]
	set fname [$f get_name]
	$f delete
    } else {
	set fname ""
    }
    $c delete

    xentry $w.name \
	-label "New scaffold name (blank for none)" \
	-default $fname
    
    okcancelhelp $w.ok \
	-ok_command "contig_add_to_scaffold $io $crec \[$w.name get\] $w; ListContigsRepopulate $io $t" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap5 Scaffolds"

    pack $w.name $w.ok -side top -fill both
}

##############################################################################
# Button-1 callback from the listcontig tablelist widget
; proc ListContigsMenuBinding {io w x y X Y} {
    # Compute row and column clicked on
    set l [$w bodypath] 
    set x [expr {$x + [winfo x $l]}]
    set y [expr {$y + [winfo y $l]}]
    foreach {row col} [split [$w nearestcell $x $y] ,] {break}

    # Find the contig identifier
    set ident [lindex [lindex [$w get $row] 0] 1]
    set ident [string range $ident 2 end-1]

    ListContigsPopupMenu $io $ident $w $x $y $X $Y
}		

; proc ListContigsSelectPressBinding {io w x y} {
    # Set the global contig id
    global c_id_contig
    foreach {line char} [split [$w index @$x,$y] .] {}
    set c_id_contig [lindex [$w get $line.0 $line.end] 0]
}

; proc ListContigsSelectReleaseBinding {io w} {
    global NGListTag

    set contigs ""
    foreach item [$w curselection] {
	set ident [lindex [lindex [$w get $item] 0] 0]
	#set cnum [db_info get_contig_num $io $ident]
	lappend contigs $ident
    }
    ListCreate2 contigs $contigs $NGListTag(contigs)
}

; proc ListContigsSelectAll {io w} {
    global NGListTag
    set contigs {}
    set end [$w index end]
    for {set item 0} {$item <= $end} {incr item} {
	set ident [lindex [lindex [$w get $item] 0] 0]
	lappend contigs $ident
    }
    ListCreate2 contigs $contigs NGListTag(contigs)
}

##############################################################################
# Called when the global "contigs" list gets updated
; proc ListContigsUpdate {io w {name1 {}} {name2 {}} {op {}}} {
    global NGList

    # Incase we get an update after the window was destroyed.
    if {![winfo exists $w]} {
	return
    }

    # Map names to list indices
    set index 0
    foreach item [$w get 0 end] {
	set name [lindex [lindex $item 0] 0]
	set ind($name) $index
	incr index
    }

    # Loop through names in NGList adding to the tablelist selection
    $w selection clear 0 end
    foreach contig $NGList(contigs) {
	if {[info exists ind($contig)]} {
	    $w selection set $ind($contig)
	}
    }
}

##############################################################################
# Called by the registration scheme when contigs are joined or broken apart.
# We need to repopulate the list.
; proc ListContigsCallback {io w type id args} {
    switch $type {
	"QUERY_NAME" {
	    return "Contig list"
	}

	"QUERY_PARAMS" {
	    return ""
	}


	"DELETE" -
	"QUIT" {
	    destroy [winfo toplevel $w]
	    return 1
	}

	"LENGTH" -
	"COMPLEMENT" -
	"JOIN_TO" -
	"ORDER" -
	"RENAME" {
	    global $w.redraw_pending
	    if {![info exist $w.redraw_pending]} {
		set $w.redraw_pending 1
		after idle "
		    ListContigsRepopulate $io $w
		    unset $w.redraw_pending
		"
	    }
	}
    }

    return ""
}

; proc ListContigsRepopulate {io w} {
    set vers [$io db_version]

    if {![winfo exists $w]} return
    set top [expr {int([lindex [$w yview] 0]*[$w index end])}]

    # Clear current setup
    $w selection clear 0 end
    $w delete 0 end

    set name_list [CreateAllContigList $io]
    set num_list [CreateAllContigListNumbers $io]
    foreach name $name_list num $num_list {
	set cstruct [$io get_contig $num]
	set clen [$cstruct get_length]
	set nreads [$cstruct nseqs]
	set nanno [$cstruct nanno]
	#set time [$cstruct get_timestamp]
	set scaffold [$cstruct get_scaffold]
	if {$vers >= 6} {
	    set clipped_len [$cstruct get_clipped_length]
	}
	if {$scaffold != 0} {
	    set fstruct [$io get_scaffold $scaffold]
	    set ctgs [$fstruct get_contigs]
	    set fname "[$fstruct get_name]/[lsearch $ctgs $num]"
	    $fstruct delete
	} else {
	    set fname "(none)/0"
	}
	$cstruct delete
	if {$vers >= 6} {
	    $w insert end [list "$name (=$num)" $clen $clipped_len $nreads $nanno $fname]
	} else {
	    $w insert end [list "$name (=$num)" $clen $nreads $nanno $fname]
	}
    }
    if {[$w sortcolumn] != -1} {
	$w sortbycolumn [$w sortcolumn] -[$w sortorder]
    }

    # Redraw the selection
    ListContigsUpdate $io $w

    $w yview $top
}

##############################################################################
# Exit handler - removes variable trace and deregisters the window
; proc ListContigsExit {io w trace_cmd} {
    global NGList
    global $w.Reg

    bind $w <Destroy> {}

    catch {trace vdelete NGList(contigs) w [list $trace_cmd]}
    catch {contig_deregister -io $io -id [set $w.Reg]}
    destroy $w
}

##############################################################################
# The "Save" button callback. This saves the order of the contig list back
# to the database.
proc ListContigsSave {io t} {
    # Select contig order
    set order ""
    foreach line [$t get 0 end] {
	lappend order [lindex $line 0]
    }

    save_contig_order -io $io -contigs $order
    update_scaffold_order -io $io
}

##############################################################################
# The "Copy" button callback. This copies the current list selection to the
# cut-n-paste selection, in a tab-formatted text form.
proc ListContigsCopy {t} {
    # Generate data string
    set data ""
    foreach line [$t curselection] {
	set text [$t get $line]
	append data \
	    [format "%s\t%s\t%s\t%s\n" \
		 [lindex $text 0] \
		 [lindex $text 1] \
		 [lindex $text 2] \
		 [lindex $text 3]]
    }
    global .Selection
    set .Selection $data

    # Selection to PRIMARY STRING
    selection own .
    selection handle . "ListContigsSelection"

    # For Windows...
    clipboard clear
    clipboard append $data
}

proc ListContigsSelection {offset maxbytes} {
    global .Selection
    return [string range ${.Selection} $offset [expr {$offset+$maxbytes-1}]]
}

##############################################################################
# Moving contig Up or Down
proc ListContigsUp {io w} {
    # If we have 8 entries, then valid indices are 0 1 2 3 4 5 6 7
    # We could have selected 2, 3 and 5.
    # "Up" then should move 2, 3 and 5 one higher, effectively reordering:
    #
    # 0 1 2 3 4 5 6 7 start
    # 0 2 1 3 4 5 6 7 <- 2
    # 0 2 3 1 4 5 6 7 <- 3
    # 0 2 3 1 5 4 6 7 <- 5

    set end [$w index end]
    set y [$w index @0,0]
    set l [$w get 0 end]

    set n 0
    foreach i [$w curselection] {
	if {$i == $n} {incr n; continue}
	set l [lreplace $l [expr {$i-1}] $i [lindex $l $i] [lindex $l [expr {$i-1}]]]
    }

    $w selection clear 0 end
    $w delete 0 end
    eval [list $w] insert end $l

#    if {[$w sortcolumn] != -1} {
#	$w sortbycolumn [$w sortcolumn] -[$w sortorder]
#    }

    # Redraw the selection
    ListContigsUpdate $io $w

    if {$y > 0} {incr y -1}
    $w yview $y
}

proc ListContigsDown {io w} {
    # Move 2, 3, 5 down =>
    #
    # 0 1 2 3 4 5 6 7 start
    # 0 1 2 3 4 6 5 7 -> 5
    # 0 1 2 4 3 6 5 7 -> 3
    # 0 1 4 2 3 6 5 7 -> 2
    
    set end [$w index end]
    set y [$w index @0,0]
    set l [$w get 0 end]

    set n [expr {[llength $l]-1}]
    foreach i [lreverse [$w curselection]] {
	if {$i == $n} {incr n -1; continue}
	set l [lreplace $l $i [expr {$i+1}] [lindex $l [expr {$i+1}]] [lindex $l $i]]
    }

    $w selection clear 0 end
    $w delete 0 end
    eval [list $w] insert end $l

#    if {[$w sortcolumn] != -1} {
#	$w sortbycolumn [$w sortcolumn] -[$w sortorder]
#    }

    # Redraw the selection
    ListContigsUpdate $io $w

    # Force the current selection to be visible
    incr y
    $w yview $y
}
