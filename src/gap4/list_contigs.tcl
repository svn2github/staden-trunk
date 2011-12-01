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
    global NGList read_only

    set t $parent.list_contigs
    if {[xtoplevel $t] == ""} return
    wm title $t "Contig List"

    # Create our tablelist
    tablelist $t.list \
	-columns {0 "Name" 10 "Length" 14 "# sequences"} \
	-labelcommand tablelist::sortByColumn \
	-exportselection 0 \
	-stretch 0 \
	-selectmode extended \
	-yscrollcommand [list $t.yscroll set]
    $t.list columnconfigure 1 -sortmode integer
    $t.list columnconfigure 2 -sortmode integer
    
    frame $t.buttons -bd 0
    button $t.buttons.cancel \
	-text "Cancel" \
	-command "destroy $t"

    button $t.buttons.save \
    	-text "Save order" \
	-command "ListContigsSave $io $t.list"
    if {$read_only} {
	$t.buttons.save configure -state disabled
    }

    button $t.buttons.copy \
    	-text Copy \
	-command "ListContigsCopy $t.list"

    button $t.buttons.help \
	-text "Help" \
	-command "show_help gap4 Contig-Selector-Contigs"

    pack $t.buttons.cancel $t.buttons.save $t.buttons.copy $t.buttons.help \
	-side left -expand 1

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

    wm geometry $t 400x200

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
		    -flags {REQUIRED LENGTH JOIN_TO DELETE COMPLEMENT}]
    
    # Populate the list
    ListContigsRepopulate $io $t.list
}

##############################################################################
# Creates a menu for the list contigs window. This is a duplicate of the
# contig selector menu, but the contig selector menu code requires knowledge
# of canvas ids and canvas tag values right down to the final menu callbacks,
# which makes it nigh on impossible to use from something other than the
# csh_win canvas.
#
# FIXME: rewrite the csh_win canvas stuff to call this menu code instead.
; proc ListContigsPopupMenu {io rnum w X Y} {
    global read_only gap_defs

    set name [io_read_reading_name $io $rnum]
    create_popup $w.m "Contig Commands ($name)"
    $w.m add command -label "Edit contig" \
	-command "edit_contig -io $io -contig {\#$rnum}"
    $w.m add command -label "Template Display" \
	-command "CreateTemplateDisplay $io {\#$rnum}"
    if {!$read_only} {
	$w.m add command -label "Complement contig" \
	    -command "complement_contig -io $io -contigs {\#$rnum}"
    }
    $w.m add command -label "List notes" \
	-command "NoteSelector $io contig {\#$rnum}"
    tk_popup $w.m [expr $X-20] [expr $Y-10]
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
    set rnum [lindex [regexp -inline  {\(\#([0-9]+)\)} [lindex [$w get $row] 0]] 1]

    ListContigsPopupMenu $io $rnum $w $X $Y
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
; proc ListContigsCallback {io w id args} {
    switch $id {
	"QUERY_NAME" {
	    return "Contig list"
	}

	"QUERY_PARAMS" {
	    return ""
	}


	"DELETE" -
	"QUIT" {
	    destroy [winfo toplevel $w]
	}

	"LENGTH" -
	"COMPLEMENT" -
	"JOIN_TO" {
	    ListContigsRepopulate $io $w
	}
    }

    return ""
}

; proc ListContigsRepopulate {io w} {
    set top [expr {int([lindex [$w yview] 0]*[$w index end])}]

    # Clear current setup
    $w selection clear 0 end
    $w delete 0 end

    set db [io_read_database $io]
    set ncontigs [keylget db num_contigs]
    for {set cnum 1} {$cnum <= $ncontigs} {incr cnum} {
        set order [contig_order_to_number -io $io -order $cnum]
        set cstruct [io_read_contig $io $order]
	set clen [keylget cstruct length]
	set num [keylget cstruct left]
	set name [io_read_reading_name $io $num]

	# Slow - should be cached!
	set nreads 0
	for {set rnum $num} {$rnum != 0} {set rnum [keylget rstruct right]} {
	    set rstruct [io_read_reading $io $rnum]
	    incr nreads
	}

	$w insert end [list "$name (#$num)" $clen $nreads]
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
	    [format "%s\t%s\t%s\n" \
		 [lindex $text 0] \
		 [lindex $text 1] \
		 [lindex $text 2]]
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
    return [string range ${.Selection} $offset $maxbytes]
}
