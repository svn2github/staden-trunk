#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc CreateModeButton { mode } {
    global $mode.button 

    set m_button(0) "EntryName"
    set m_button(1) "AccessionNumber"

    set $mode.button 0
    menubutton $mode -text $m_button([set $mode.button]) -indicatoron 1 \
	    -menu $mode.m -bd 2 -relief raised -width 15
    menu $mode.m -tearoff 0

    for {set i 0} {$i < 2} {incr i} {
	$mode.m add command -label $m_button($i) \
		-command "$mode configure -text $m_button($i); \
		set $mode.button $i"
    }
}

proc GetMode { mode } {
    global $mode.button
    return [set $mode.button]
}

#find entrynames within a personal library file
proc InvokeArchiveBrowser { l fn entry} {

    set list [get_archive_list -file [entrybox_get $fn]]

    if {[llength $list] == 0} {
	return
    }
    set lb [SipListBox $l $list]
    bind $lb <Double-1> "destroy $l; break"
    bind $lb <1> "
    	entrybox_delete $entry 0 end;
	entrybox_insert $entry 0 \[%W get \[%W index @%x,%y\]\]
    "
}

#library browser
proc InvokeLibBrowser { s f num_lib mode lib_list} {
    global $s.options
    global sip_defs

    set library [GetCheckLib $f.lib]
    set mode [GetMode $mode]

    if {[xtoplevel $s -resizable 0] == ""} return
    SetCurFrame $s $f
    set $s.options "$library $mode "
    trace variable $s.options w "UpdateEntry $s {$lib_list}"
    SeqLibraries $s [set $s.options]
    tkwait variable $s.options
}

proc SipListBox {l list } {

    if {[xtoplevel $l -resizable 0] == ""} return
    listbox $l.lists -yscrollcommand "$l.scrolly set"
    scrollbar $l.scrolly -command "$l.lists yview" -orient vertical
    button $l.cancel -text Cancel -command "destroy $l"

    pack $l.cancel -side bottom -fill x
    pack $l.scrolly -side right -fill y
    pack $l.lists -fill both
    foreach i $list {
	$l.lists insert end $i
    }
    return $l.lists

}

proc UpdateEntry {s lib_list name element op} {

    set f [GetCurFrame $s]

    #check that the sequence browser has not been destroyed!
    if {![winfo exists $f]} {
	bell
	tk_messageBox -icon error -type ok -title "Error" \
		-message "No sequence browser" \
		-parent $f
	return
    }
    #update text entry box
    upvar $name x
    entrybox_delete $f.word 0 end
    entrybox_insert $f.word 0 [lindex $x 2]
    #update library
    SetCheckLib $f.lib [lindex $x 0] $lib_list
}

#unused?
proc seq_entry { f } {
    global seqlib_defs

    set s [keylget seqlib_defs SEQ.WIN]
    frame $f

    SetCurFrame $s $f

    #create list of available libraries
    set lib_list [list_seq_libs]
    lappend lib_list {personal file}
    set num_lib [llength $lib_list]

    #search word
    entrybox $f.word -width 15 -type "CheckStringExists"

    #filename
    entrybox $f.fn -title Filename -width 15 \
	    -state disabled

    button $f.browse1 -text Browse -command "InvokeFileBrowser $f.fn open" \
	    -state disabled
    button $f.browse2 -text Browse -command "InvokeLibBrowser $s $f $num_lib \
	    $f.mode {$lib_list}"

    CreateLibListPersonal $f.lib $lib_list -1 $f.fn $f.browse1 $f.browse2 $f.word $s $f $f.mode

    CreateModeButton $f.mode
    pack $f.lib $f.fn $f.browse1 $f.mode $f.word $f.browse2 -side left
    bindtags $f.word "Text . all $f.word"
#    bind $f.word <Any-FocusIn> "puts here1; SetCurFrame $s $f"
    bind [entrybox_path $f.word] <1> "SetCurFrame $s $f"

#HACK!!!!
    #entrybox_insert $f.word 0 hsproperd
}

proc get_seq_entry { f } {

    set lib [GetCheckLib $f.lib]
    set mode [GetMode $f.mode]
    set file [entrybox_get $f.fn]
    set word [entrybox_get $f.word]

    return "{$lib} {$mode} {$word} {$file}"
}

