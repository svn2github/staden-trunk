#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#
# TODO:
# *	Add note type to experiment file and enter them at assembly
#
# *	Add actions to note - eg bring up editor at pos X

proc get_temp_file {} {
    global temp_file_index
       

    if {![info exists temp_file_index]} {set temp_file_index 0}
    set prefix "/tmp/tcl_[pid]_"
    while {[file exists $prefix$temp_file_index]} {
	incr temp_file_index	
    }

    return $prefix$temp_file_index
}

proc NoteSelectorCreate {io} {
    global gap_defs read_only
    global db_namelen
    set w [keylget gap_defs NOTE_SELECTOR.WIN]
    global $w.Type
    global $w.Ident
    global $w.Reg
    global $w.Buffer
    global $w.Redisp

    xtoplevel $w
    wm title $w "Edit notebooks"
    set $w.Type ""
    set $w.Ident ""
    set $w.Buffer 0
    set $w.Redisp 0

    # Register with the 1st contig. REG_NOTE requests are sent to all contigs
    # anyway.
    set $w.Reg [contig_register \
	-io $io \
	-contig 0 \
	-command [list NoteSelectorCallback $io $w] \
	-flags {REQUIRED NOTE BUFFER_START BUFFER_END JOIN_TO}]

    wm protocol $w WM_DELETE_WINDOW "NoteSelectorExit $io $w"

    # The note menus
    global select_notes_menu
    $w configure -menu $w.menubar
    menu $w.menubar
    create_menus $select_notes_menu $w.menubar
    if {$read_only} {
        menu_state_set select_notes_menu -16 $w.menubar
    }

    # The 'select read' window. Initially ungridded
    entrybox $w.which_read \
	-title "Please enter reading identifier" \
	-type "CheckContigName $io" \
	-width $db_namelen \
	-command "SelectReadNote2 $io $w $w.which_read"

    # The 'select contig' window. Initially ungridded
    contig_id $w.which_contig \
	-io $io \
	-range 0 \
	-command "SelectContigNote2 $io $w $w.which_contig"
    entrybox_configure $w.which_contig.ent \
	    -title "Please enter contig identifier"

    # The note selector window
    text $w.selection \
	-width 100 \
	-height 10 \
	-yscrollcommand "$w.sel_ys set" \
	-xscrollcommand "$w.sel_xs set" \
	-wrap none
    scrollbar $w.sel_ys \
	-orient vert \
	-command "$w.selection yview"
    scrollbar $w.sel_xs \
	-orient horiz \
	-command "$w.selection xview"

    # Text tags for the note selector
    $w.selection tag configure selected -background grey70 -relief raised -borderwidth 1
    $w.selection tag configure selectable
    $w.selection tag configure highlight -relief raised -borderwidth 2
    catch {eval font create bold_text_font [font configure text_font] -weight bold}
    $w.selection tag configure title -font bold_text_font

    # Bindings in the note selector
    bind $w.selection <1> {
        tk::TextButton1 %W %x %y
        set ::tk::Priv(selectMode) line
        %W tag remove sel 0.0 end
        tk::TextSelectTo %W %x %y
        catch {%W mark set insert sel.first}
	%W tag remove sel 1.0 3.0
    }
    bind $w.selection <B1-Motion> {
        set ::tk::Priv(x) %x
        set ::tk::Priv(y) %y
        tk::TextSelectTo %W %x %y
	%W tag remove sel 1.0 3.0
    }
    $w.selection tag bind selectable <<use>> "\
	set lastLine \[%W index {@%x,%y linestart}\]
	%W tag remove selected \$lastLine \"\$lastLine lineend\"
	EditNote $io \[NoteSelectorNoteNum %W %x %y\]"
    $w.selection tag bind selectable <Any-Enter> {
	set lastLine [%W index {@%x,%y linestart}]
	%W tag add highlight $lastLine "$lastLine lineend"
    }
    $w.selection tag bind selectable <Any-Motion> {
	set lastLine [%W index {@%x,%y linestart}]
	%W tag remove highlight 1.0 end
	%W tag add highlight $lastLine "$lastLine lineend"
    }
    $w.selection tag bind selectable <Any-Leave> {
	%W tag remove highlight 1.0 end
    }

    # Remove most default Text bindings
    bindtags $w.selection "$w.selection . all"
    bind $w.selection <2> [bind Text <2>]
    bind $w.selection <B2-Motion> [bind Text <B2-Motion>]
    bind $w.selection <ButtonRelease-2> [bind Text <ButtonRelease-2>]

    # Grid them
    NoteSelectorGridNotes $w notes
}

# Grid and ungrid the note display.
proc NoteSelectorGridNotes {w mode} {
    if {$mode == "notes"} {
	grid columnconfigure $w 0 -weight 1
	grid rowconfigure    $w 1 -weight 1
	grid $w.selection -row 1 -column 0 -sticky nsew
	grid $w.sel_ys    -row 1 -column 1 -sticky ns
	grid $w.sel_xs    -row 2 -column 0 -sticky ew
	catch {grid forget $w.which_contig}
	catch {grid forget $w.which_read}
    } elseif {$mode == "contigs"} {
	grid forget $w.selection $w.sel_ys $w.sel_xs
	catch {grid forget $w.which_read}
    } else {
	grid forget $w.selection $w.sel_ys $w.sel_xs
	catch {grid forget $w.which_contig}
    }
}


#
# The 'note' selector.
#
proc NoteSelector {io {type {}} {ident {}}} {
    global gap_defs
    set w [keylget gap_defs NOTE_SELECTOR.WIN]

    if {![winfo exists $w]} {
	NoteSelectorCreate $io
    } else {
	raise $w
	wm deiconify $w
    }

    if {$type == "database" } {
	SelectDBNote $io $w
    } elseif {$type == "reading"} {
	SelectReadingNote3 $io $w $ident
    } elseif {$type == "contig"} {
	SelectContigNote3 $io $w $ident
    }
}

proc NoteSelectorExit {io w} {
    global $w.Reg

    contig_deregister -io $io -id [set $w.Reg]
    destroy $w
}

proc NoteSelectorNoteNum {w x y} {
    # Extract from text window and parse out the note number
    set lastLine [$w index "@$x,$y linestart"]
    set n [$w get $lastLine "$lastLine lineend"]
    set tags [$w tag names $lastLine]
    regexp {#([0-9]+)} $tags dummy nn
    return $nn
}

proc NoteSelectorCallback {io w type id args} {
    global $w.Buffer $w.Redisp

    if {$type == "QUERY_NAME"} {
	return "Note selector"
    }
    if {$type == "QUERY_PARAMS"} {
	return "none"
    }
    if {$type == "DELETE" || $type == "QUIT"} {
	NoteSelectorExit $io $w
    }

    if {$type == "BUFFER_START"} {
	incr $w.Buffer
	return
    }

    if {$type == "BUFFER_END"} {
	incr $w.Buffer -1
	if {[set $w.Buffer] == 0 && [set $w.Redisp]} {
	    set $w.Redisp 0
	    SelectNoteRedisplay $io $w
	}
    }
	
    if {$type == "NOTE" || $type == "JOIN_TO"} {
	if {[set $w.Buffer] > 0} {
	    set $w.Redisp 1
	} else {
	    SelectNoteRedisplay $io $w
	}
    }
}

proc SelectNoteRedisplay {io w} {
    global $w.Type $w.Ident

    if {[set $w.Type] == "database"} {
	SelectDBNote $io $w
    } elseif {[set $w.Type] == "reading"} {
	SelectReadingNote3 $io $w [set $w.Ident]
    } else {
	SelectContigNote3 $io $w [set $w.Ident]
    }
}

# Tcl interface to C new_note function
proc NewNote {io w} {
    global $w.Type $w.Ident

    if {[set $w.Type] == ""} {
	bell
	return
    }

    if {[set $w.Type] == "database"} {
	set num 0
    } elseif {[set $w.Type] == "reading"} {
	set num [db_info get_read_num $io [set $w.Ident]]
    } else {
	set num [db_info get_contig_num $io [set $w.Ident]]
    }

    set nn [new_note \
		-io $io \
		-type COMM \
		-to [set $w.Type] \
		-number $num]
    if {$nn != -1} {
	EditNote $io $nn
    }
}

# Tcl interface to C delete_note function
proc DeleteNote {io w} {
    set range [$w.selection tag ranges sel]

    if {$range == ""} {
	return
    }

    set start [expr int([lindex $range 0])]
    set end   [expr int([lindex $range 1])]
    set del_line ""
    set note_list ""
    for {set i $start} {$i < $end} {incr i} {
	set tags [$w.selection tag names $i.0]
	if {0 == [regexp {#([0-9]+)} $tags dummy nn]} {
	    continue
	}
	lappend note_list $nn
	set del_line "$i $del_line"
    }

    contig_notify -io $io -type BUFFER_START -cnum 0 -args {}
    foreach nn $note_list {
	delete_note -io $io -note $nn
    }
    contig_notify -io $io -type BUFFER_END -cnum 0 -args {}

#    # Remove from display
#    foreach l $del_line {
#	$w.selection delete $l.0 [expr {$l+1}].0
#    }
}

proc ListSelectedNotes {io w first clear {from {}}} {
    if {$clear == 1} {
	$w.selection delete 1.0 end
	# Header
	$w.selection insert end "Notes from: "
	$w.selection insert end "$from\n" title
	$w.selection insert end \
	  "Location of note              Number Type  Creation date        Modifcation date     1st line" title
    }

    if {$first == 0} {
	$w.selection insert end "\n -- No notes found --"
	return
    }

    for {set i $first} {$i != 0} {set i [keylget n next]} {
	set n [io_read_note $io $i]
	if {[keylget n annotation]} {
	    set title [io_read_text $io [keylget n annotation]]
	    set title [lindex [split $title \n] 0]
	} else {
	    set title ""
	}
	$w.selection insert end \
	    [format "\n%-30s  %-4d %.4s  %s  %s  %s" \
		$from $i [keylget n type] \
		[clock format [keylget n ctime] -format {%Y/%m/%d %T}] \
		[clock format [keylget n mtime] -format {%Y/%m/%d %T}] \
		$title] \
	    [list #$i selectable]
    }
}

proc SelectDBNote {io w} {
    global $w.Type
    set $w.Type database
    NoteSelectorGridNotes $w notes
    set db [io_read_database $io]
    ListSelectedNotes $io $w [keylget db notes] 1 Database
}

proc SelectReadingNote {io w} {
    NoteSelectorGridNotes $w readings
    grid $w.which_read -row 0 -column 0 -sticky nsw
}

proc SelectReadNote2 {io w which_read item} {
    set rname [entrybox_get $which_read]
    NoteSelectorGridNotes $w notes
    SelectReadingNote3 $io $w $rname
}

proc SelectReadingNote3 {io w ident} {
    set rnum [db_info get_read_num $io $ident]
    set r [io_read_reading $io $rnum]
    set rname [io_read_text $io [keylget r name]]
    ListSelectedNotes $io $w [keylget r notes] 1 "Read $rname"

    global $w.Type $w.Ident
    set $w.Type reading
    set $w.Ident $ident
}

proc SelectContigNote {io w} {
    NoteSelectorGridNotes $w contigs
    grid $w.which_contig -row 0 -column 0 -sticky nsw
}

proc SelectContigNote2 {io w which_contig} {
    set cname [contig_id_gel $which_contig]
    NoteSelectorGridNotes $w notes
    SelectContigNote3 $io $w $cname
}

proc SelectContigNote3 {io w ident} {
    set cnum [db_info get_contig_num $io $ident]
    set c [io_read_contig $io $cnum]
    set cname [left_gel $io $cnum]
    ListSelectedNotes $io $w [keylget c notes] 1 "Contig $cname"

    global $w.Type $w.Ident
    set $w.Type contig
    set $w.Ident $ident
}

proc EditNoteCheckFd {fd n} {
    global done_view_note_$n
    read $fd
    if {[eof $fd]} {close $fd; set done_view_note_$n 1}
}

# Returns the index into NoteDB for note with type 'id'
proc FindNoteInd {id} {
    global NoteDB

    for {set i 0} {$i < $NoteDB(num_notes)} {incr i} {
	if {[string compare $NoteDB($i,id) $id] == 0} {
	    return $i
	}
    }

    return -1
}

proc EditNoteType {w menubar item type} {
    global $w.Type
    global NoteDB

    if {[winfo exists $w.note.t]} {
        set anno [$w.note.t get 1.0 end]
        regsub "\n$" $anno {} anno
       	if {$anno == " -- No text attached to this note --"} {
	    set anno ""
       	}

        set ind [FindNoteInd [set $w.Type]]
        if {[string compare $anno $NoteDB($ind,dt)] == 0} {
	    $w.note.t delete 1.0 end
            set ind [FindNoteInd $type]
	    if {$NoteDB($ind,dt) != ""} {
	        $w.note.t insert 1.0 $NoteDB($ind,dt)
	    } else {
		$w.note.t insert 1.0 " -- No text attached to this note --"
	    }
        }
    }

    set $w.Type $type
    $menubar entryconfigure $item -label Type:$type
}

proc EditNoteSave {io w nn} {
    global $w.Type $w.SavedType $w.SavedText
    set anno [$w.note.t get 1.0 end]
    regsub "\n$" $anno {} anno
    if {$anno == " -- No text attached to this note --"} {
        set anno ""
    }
    set type [set $w.Type]

    edit_note -io $io -note $nn -type $type -comment $anno
    set $w.SavedType $type
    set $w.SavedText $anno
}

proc EditNoteDelete {io w nn} {
    if {[tk_messageBox \
	    -icon question \
	    -title "Delete Note" \
	    -message "Do you wish to delete this note?" \
	    -type yesno \
	    -parent $w] == "yes"} {
	delete_note -io $io -note $nn
    }
}

proc EditNoteExit {io w nn} {
    global $w.Reg $w.Type $w.SavedType $w.SavedText read_only
    set type [set $w.Type]
    regsub "\n*\$" [$w.note.t get 1.0 end] "" text
    regsub "\n*\$" [set $w.SavedText] "" savedtext

    if {$read_only == 0 &&
	([string compare [set $w.SavedType] $type] != 0 || \
	 [string compare $savedtext $text] != 0)} {
        set ans [tk_messageBox \
		-icon question \
		-title "Save Note" \
		-message "There are unsaved changes to this note. Do you wish to save them?" \
		-type yesnocancel \
		-parent $w]
	if {$ans == "yes"} {
	    # Yes
	    EditNoteSave $io $w $nn
        } elseif {$ans == "cancel"} {
	    # Cancel
	    return
	}
    }

    if {[info exists $w.Reg]} {
        contig_deregister -io $io -id [set $w.Reg]
    }
    destroy $w
}

proc EditNote {io nn} {
    global gap_defs read_only
    global edit_note_menu

    set w [keylget gap_defs NOTE_EDITOR.WIN]$nn

    # Create a note view/editor window
    if {[xtoplevel $w] == ""} return
    wm title $w "Edit note: #$nn"
    global $w.Reg $w.SavedText $w.SavedType

    wm protocol $w WM_DELETE_WINDOW "EditNoteExit $io $w $nn"
	
    # Any contig will do as the callbacks are the standard ones plus
    # REG_NOTE ones (sent to all contigs).
    set $w.Reg [contig_register \
	-io $io \
	-contig 0 \
	-command [list EditNoteCallback $io $w $nn] \
	-flags {REQUIRED NOTE}]

    set n [io_read_note $io $nn]
    if {[keylget n annotation]} {
        set note_text [io_read_text $io [keylget n annotation]]
    } else {
        set note_text " -- No text attached to this note --"
    }

    # The menus
    $w configure -menu $w.menubar
    menu $w.menubar
    create_menus $edit_note_menu $w.menubar
    if {$read_only} {
        menu_state_set edit_note_menu -16 $w.menubar
    }

    # Type menu is dynamically created
    set t_ind [$w.menubar index Type]
    add_note_menus $w.menubar.type "EditNoteType $w $w.menubar $t_ind"
    EditNoteType $w $w.menubar $t_ind [keylget n type]
    set $w.SavedType [keylget n type]
    set $w.SavedText $note_text

    # The creation/modification times
    frame $w.times
    frame $w.times.c
    label $w.times.c.l -text "Creation date \t:"
    label $w.times.c.r \
	-text "[clock format [keylget n ctime] -format {%Y/%m/%d %T}]" \
	-font text_font
    pack $w.times.c.l $w.times.c.r -side left

    frame $w.times.m
    label $w.times.m.l -text "Modification date \t:"
    label $w.times.m.r \
	-text "[clock format [keylget n mtime] -format {%Y/%m/%d %T}]" \
	-font text_font
    pack $w.times.m.l $w.times.m.r -side left

    pack $w.times.c $w.times.m -side left -expand 1
    pack $w.times -side top -fill both


    # The editor panel itself
    if {[set com [keylget gap_defs NOTE_EDITOR.PROGRAM]] != ""} {
	# Start up NOTE_EDITOR to perform the editing
	set tmp_file [get_temp_file]
	set fd [open $tmp_file w]
	puts $fd $note_text
	close $fd

	global done_view_note_$nn
	set fd [open "|[keylget gap_defs NOTE_EDITOR.PROGRAM] $tmp_file	" r]
	fconfigure $fd -blocking 0
	fileevent $fd readable "EditNoteCheckFd $fd $nn"
	vwait done_view_note_$nn

	set fd [open $tmp_file r]
	set note_text [read $fd]
	close $fd
    } else {
	# Create our own simple editor window
	# The scrolled text widget
	frame $w.note
	pack $w.note -side bottom -fill both -expand 1
        text $w.note.t \
	    -width 80 \
            -yscrollcommand "$w.note.ys set" \
            -xscrollcommand "$w.note.xs set"
        scrollbar $w.note.ys \
            -orient vert \
            -command "$w.note.t yview"
        scrollbar $w.note.xs \
            -orient horiz \
            -command "$w.note.t xview"
    
        grid columnconfigure $w 0 -weight 1    
        grid rowconfigure    $w 0 -weight 1 
        grid columnconfigure $w.note 0 -weight 1    
        grid rowconfigure    $w.note 0 -weight 1 
        grid $w.note.t  -row 0 -column 0 -sticky nsew
        grid $w.note.ys -row 0 -column 1 -sticky ns
        grid $w.note.xs -row 1 -column 0 -sticky ew
    
	$w.note.t insert 1.0 $note_text
    }
}

proc EditNoteCallback {io w nn type id args} {
    if {$type == "QUERY_NAME"} {
	return "Note editor"
    }

    if {$type == "QUERY_PARAMS"} {
	return "none"
    }

    if {$type == "DELETE" || $type == "QUIT"} {
	global $w.Reg
	contig_deregister -io $io -id [set $w.Reg]
	destroy $w
	return
    }

    if {$type == "NOTE"} {
	if {[keylget args note] != $nn} {
	    return
	}
	if {[keylget args task] == "DELETE"} {
	    global $w.Reg
	    contig_deregister -io $io -id [set $w.Reg]
	    destroy $w
	    return
	}
    }
}

proc add_note_menus {mpath command} {
    global NoteDB gap_defs

    for {set i 0; set j 1} {$i < $NoteDB(num_notes)} {incr i; incr j} {
        $mpath add command \
	    -label "$NoteDB($i,id): $NoteDB($i,type)" \
 	    -command "$command $NoteDB($i,id)"
	if {$j >= [keylget gap_defs MAX_MENU_ITEMS]} {
	    $mpath add cascade -label "More..." -menu $mpath.more
	    set mpath [menu $mpath.more -tearoff 0]
	    set j 0
	}
    }
}
