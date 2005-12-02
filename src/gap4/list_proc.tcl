#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
set NGList() ""
set NGLists ""
set _NEXT_LIST_EDITOR 0

proc null_func {} {
}

proc ListGet {name} {
    global NGList io

    if {![info exists io]} {
	return $NGList($name)
    }

    if {$io >= 0 && [regexp {^\[(.*)\]$} $name dummy id] == 1} {
	set l {}
	for {set rn [db_info chain_left $io $id]} {$rn > 0} {} {
	    set r [io_read_reading $io $rn]
	    lappend l [io_read_text $io [keylget r name]]
	    set rn [keylget r right]
	}
	return $l
    } elseif {[regexp {^{(.*)}$} $name dummy id] == 1} {
	set l {}
	foreach i [ListGet $id] {
	    set l "$l [ListGet \[$i\]]"
	}
	return $l
    } elseif {$io >= 0 && $name == "allcontigs"} {
        return [CreateAllContigList $io]
    } elseif {$io >= 0 && $name == "allreadings"} {
	return [ListGet "{allcontigs}"]
    } else {
        return $NGList($name)
    }
}

#############################################################################
# List file IO
proc ListLoad {file name {tag {}}} {
    global NGList NGLists

    if {[catch {set f [open "$file"]}]} {
	bell
	return
    }

    set list ""
    while {[gets $f line] >= 0} {
	lappend list $line
    }
    close $f

    # Special case for readings list
    if {$name == "readings"} {
	SetReadingList $list
    } else {
	ListCreate $name $list
	if {$tag != ""} {
	    global NGListTag
	    set NGListTag($name) $tag
	}
    }

    return $list
}

proc ListSave {file name} {
    global NGList

    if {![ListExists $name]} {
	bell
	return 0
    }

    if {[catch {set f [open "$file" w]}]} {
	bell
	return 1
    }

    foreach i [ListGet $name] {
	puts $f $i
    }

    close $f
    return 1
}

#############################################################################
# Creating, appending, deleting, copying and checking
proc ListCreate2 {name list {tag {}}} {
    global NGList NGLists NGListTag

    if {![info exists NGList($name)]} {
        lappend NGLists $name
    } 
    set NGList($name) $list
    set NGListTag($name) $tag
}

proc ListCreate {name list} {
    global NGList NGLists

    if {![ListWritable $name]} {return 0}
    ListCreate2 $name $list
    return 1
}

proc ListAppend {name list} {
    global NGList NGListTag

    ListEditSave $name
    foreach item $list {
	lappend NGList($name) $item
    }
    if {[info exists NGListTag($name)] && $NGListTag($name) == "SEQID"} {
	ListRemoveDuplicates $name
    }
}

proc ListRemoveDuplicates {name} {
    global NGList

    # Convert list into array
    set l ""
    set ind 0
    foreach item $NGList($name) {
	lappend l $item $ind
	incr ind
    }
    array set dummy $l

    # Reverse the array to order via index
    set l ""
    foreach {item val} [array get dummy] {
	lappend l $val $item
    }
    unset dummy
    array set dummy $l

    # Convert the array back into a list
    set l ""
    foreach ind [lsort -integer [array names dummy]] {
	lappend l $dummy($ind)
    }
    set NGList($name) $l
}

proc ListSubtract {name list} {
    global NGList NGListTag

    ListEditSave $name
    set l $NGList($name)
    foreach item $list {
	while {[set ind [lsearch -exact $l $item]] != -1} {
	    set l [lreplace $l $ind $ind]
	}
    }
    set NGList($name) $l
}

proc ListClear {name} {
    global NGList NGLists

    set NGList($name) ""
}

proc ListPrint {name} {
    global NGList tk_utils_defs

    if {![ListExists $name]} {
	bell
	return 0
    }

    ListEditSave $name

    if {[catch {set f [open "|[keylget tk_utils_defs PRINT_COMMAND]" w]}]} {
	bell
	return 0
    }

    # Title
    puts $f "List:$name       Date:[clock format [clock seconds]]\n"

    # List
    foreach i [ListGet $name] {
	puts $f $i
    }

    close $f
    return 1
}

proc ListDelete {name} {
    global NGList NGLists NGSpecial NGListTag

    if {[ListExists2 $name]} {
	if {[info exists NGSpecial] && [lsearch -exact $NGSpecial $name] != -1 } {
	    tk_messageBox -icon error -type ok -title "Special list" \
		    -message "This list cannot be deleted"
	    return 0
	}

	unset NGList($name)
	unset NGListTag($name)
	set NGLists [lreplace $NGLists [lsearch -exact $NGLists $name] [lsearch -exact $NGLists $name]]
	return 1
    } else {
	tk_messageBox -icon error -type ok -title "No such list" \
		-message "List does not exist"
	return 0
    }
}

proc ListCopy {from to} {
    global NGList NGLists

    if {$from == $to} {
	tk_messageBox -icon error -type ok -title "Copy list" \
	    -message "Source and destination list names need to be different"
	return 0
    }

    if {![ListExists $from]} {return 0}
    return [ListCreate $to [ListGet $from]]
}

proc ListExists2 {name} {
    global NGList

    return [info exists NGList($name)]
}

proc ListExists {name {gap 1}} {
    global NGList io

    if {![ListNameValid $name]} {
	return 0
    }

    if {[info exists io] && $gap && $io >= 0 && \
	[regexp {^\[(.*)\]$} $name dummy id] == 1 && \
	[db_info get_read_num $io $id] != -1} {
	return 1
    }

    if {$gap && [regexp {^{(.*)}$} $name dummy id] == 1 && \
	[info exists NGList($id)]} {
	return 1
    }
    
    if {[info exists NGList($name)]} {
    	return 1
    } else {
	tk_messageBox -icon error -type ok -title "No such list" \
		-message "List does not exist"
	return 0
    }
}

proc ListWritable {name} {
    global NGSpecial

    if {![ListNameValid $name]} {
	return 0
    }

    if {[ListExists2 $name]} {

	if {[lsearch -exact $NGSpecial $name] != -1 } {
	    tk_messageBox -icon error -type ok -title "Special list" \
		    -message "This list cannot be replaced"
	    return 0
	}

	set answer [tk_messageBox -icon question -type yesno \
		-title "List exists" -message "List already exists. Replace?"]
	case $answer {
	    no {return 0} 
	    yes {ListDelete $name}
	}
    }

    return 1
}

proc ListNameValid {name} {
    if {"$name" == ""} {
	tk_messageBox -icon error -type ok -title "No list name" \
		-message "You have not entered a list name"
	return 0
    }
    return 1
}

proc ListBrowse_ok {entry listbox} {
    entrybox_delete $entry 0 end
    if {[$listbox curselection] != ""} {
	set item [$listbox get [$listbox curselection]]
	entrybox_insert $entry 0 $item
	return $item
   } else {
	return ""
    }
}

proc ListBrowse {path entry {type {read}}} {
    set l [ListListbox $path.browser "ListBrowse_ok $entry" $type] 

    bind $l <<use>> "
    	entrybox_delete $entry 0 end;
	entrybox_insert $entry 0 \[%W get \[%W index @%x,%y\]\]
	destroy $path.browser; break
    "
}

proc CheckOutList {path} {
    return [ListWritable [$path.entry get]]
}

#############################################################################
# Conversion of of string-list-item types
proc ListToItems {name} {
    global NGList

    if {![ListExists $name]} {return ""}

    ListEditSave $name

    set items ""
    foreach i [ListGet $name] {
	set l [lindex $i 0]
	if {"$l" != ""} {
	    lappend items $l
       	}
    }

    return $items
}

proc StringToList {string} {
    return [split $string \n]
}


#############################################################################
# The list editing display
proc ListNextEditor {} {
    global _NEXT_LIST_EDITOR

    return [incr _NEXT_LIST_EDITOR]
}

# Does an editor for this list exist already?
proc ListEditExists {name} {
    global gap_defs

    if {![winfo exists [keylget gap_defs EDIT_LIST.WIN]]} {
	return ""
    }
	
    foreach i [winfo children [keylget gap_defs EDIT_LIST.WIN]] {
	global $i.CurList
	if {[info exists $i.CurList] && [set $i.CurList] == "$name"} {
	    return $i	    
	}
    }

    return ""
}

# Saves a list editor - if needed
proc ListEditSave {name} {
    global NGList

    if {[set e [ListEditExists $name]] != ""} {
	#YUK!
	regsub \n*$ [$e.list get 1.0 end] "" t
	set NGList($name) [StringToList $t]
    }
}

# Updates a list editor
proc ListEditUpdate {t name args} {
    global NGList NGListTag

    if {![winfo exists $t]} {
        global $t.CurList
        trace vdelete NGLists w "UpdateListListbox $t"
	if {[info exists $t.CurList]} {
	    unset $t.CurList
	}
	return
    }

    $t delete 1.0 end
    if {$NGListTag($name) != ""} {
	if {$NGListTag($name) == "SEQID"} {
	    foreach i $NGList($name) {
		regexp {([[:space:]]*)([^[:space:]]*)(.*)} $i _ l m r
		$t insert end $l "" $m $NGListTag($name) "$r\n" ""
	    }
	} else {
	    foreach i $NGList($name) {
		$t insert end "$i" $NGListTag($name) "\n"
	    }
	}
    } else {
	foreach i $NGList($name) {
	    $t insert end "$i\n"
	}
    }
}

#invoke a filebrowser from a browse button attached to EditList
#displays ALL files, no filetype has been implemented 
proc ListEditBrowse { t name } {
    global NGList

    set filelist [tk_getOpenFile -multiple 65000 -parent $t]
    set fcount [llength $filelist]

    #Strip off cwd
    for {set i 0 } { $i < $fcount } { incr i } {
	if {[string match [pwd]* [lindex $filelist $i]]} {
          set filelist \
	  [lreplace $filelist $i $i [string range [lindex $filelist $i] [expr [string length [pwd]]+1] end]]
      }
    }

    #Add to the list
    ListEditSave $name
    set NGList($name) [concat $NGList($name) $filelist]
    ListEditUpdate $t.list $name
}

proc ListEdit {name {read_only 0}} {
    global gap_defs NGList NGSpecial tcl_platform
    
    if {![ListExists $name 0]} {bell; return 0}

    if {![winfo exists [set w [keylget gap_defs EDIT_LIST.WIN]]]} {
	frame $w
    }

    # We cannot edit the special lists
    if {[lsearch -exact $NGSpecial $name] != -1} {
	set read_only 1
    }

    # Create the new editor if it doesn't exist
    if {[set t [ListEditExists $name]] == ""} {
        set t "$w.[ListNextEditor]"
        xtoplevel $t
	if {$read_only} {
	    wm title $t "Viewing list: \"$name\""
	} else {
	    wm title $t "Editing list: \"$name\""
	}
    
        global $t.CurList
        set $t.CurList $name
        
        # The menu bar
        frame $t.mbar -relief raised -bd 2
	button $t.mbar.quit \
	    -text Ok \
	    -command "ListEditSave [list $name];destroy $t"
	button $t.mbar.clear \
	    -text Clear \
	    -command "ListClear [list $name]"
	button $t.mbar.search \
	    -text Search \
	    -command "tout_search $t.list"
	if { $tcl_platform(platform) != "windows" } {
	    button $t.mbar.print \
		-text Print \
		-command "ListPrint [list $name]"
	}
	button $t.mbar.save \
	    -text Save \
	    -command "SaveThisListDialog [list $name] $t.save"
	button $t.mbar.browse \
	    -text Browse \
	    -command "ListEditBrowse $t [list $name]"
	if {$read_only} {
	    if {$name != "readings" && $name != "contigs"} {
		$t.mbar.clear configure -state disabled
	    }
	    $t.mbar.browse configure -state disabled
	}

       	wm protocol $t WM_DELETE_WINDOW "ListEditSave [list $name]; destroy $t"
    
        # The text display and it's scrollbars.
        text $t.list -width 30 -wrap none\
    	    -xscrollcommand "$t.f.scrollx set" \
    	    -yscrollcommand "$t.scrolly set"
	if {$read_only} {
	    bindtags $t.list "$t.list TextRO . all"
	}
	gap4_text_init $t.list
        frame $t.f
        frame $t.f.padding
        scrollbar $t.f.scrollx -orient horiz -command "$t.list xview"
        scrollbar $t.scrolly   -orient vert  -command "$t.list yview"

	if { $tcl_platform(platform) != "windows" } {
	    pack $t.mbar.quit $t.mbar.clear $t.mbar.search $t.mbar.print \
		$t.mbar.save $t.mbar.browse -side left
	} else {
	    pack $t.mbar.quit $t.mbar.clear $t.mbar.search $t.mbar.browse \
		$t.mbar.save -side left
	}
        pack $t.mbar -side top -fill x
    
        pack $t.f -side bottom -fill x
        pack $t.scrolly -side right -fill y
        pack $t.list -expand 1 -fill both
        pack $t.f.scrollx -side left -fill x -expand 1
        pack $t.f.padding -side right
    
        # What a lovely hack
        pack propagate $t.f.padding 0
        $t.f.padding configure -width 18 -height 18
        update idletasks
        $t.f.padding configure \
            -width  [winfo width $t.scrolly] \
            -height [winfo width $t.scrolly]

	# Monitor changes to this list
	trace variable NGList($name) w "ListEditUpdate $t.list [list $name]"
	trace variable NGList($name) u "destroy $t"
    } else {
        raise $t
    }

    # But we always update the list
    ListEditUpdate $t.list $name

    return 1
}

#############################################################################
# List dialogues operations

#
# A scrolled listbox containing a list of available "List"s.
#
proc ListListboxQuit {t} {
    global NGLists
    global $t.CurList

    trace vdelete NGLists w "UpdateListListbox $t"
    if {[info exists $t.CurList]} {
        unset $t.CurList
    }
    if {[winfo exists $t]} {destroy $t}
}

proc ListListbox_ok {lists} {
    return [$lists selection get]
}

proc ListListbox {t {ok_command ListListbox_ok} {type read}} {
    global NGList NGLists NGSpecial

    if {[xtoplevel $t] == ""} return
    wm title $t "List browser"
    global $t.Type
    set $t.Type $type

    listbox $t.lists -yscrollcommand "$t.scrolly set"
    scrollbar $t.scrolly -command "$t.lists yview" -orient vertical

    foreach i $NGLists {
	if {$type == "read" || [lsearch -exact $NGSpecial $i] == -1} {
	    $t.lists insert 0 $i
	}
	if {$type == "load" && $i == "readings"} {
	    $t.lists insert 0 $i
	}
	if {$type == "edit" && ($i == "readings" || $i == "contigs")} {
	    $t.lists insert 0 $i
	}
    }

    trace variable NGLists w "UpdateListListbox $t"
    wm protocol $t WM_DELETE_WINDOW "null_func"

    frame $t.f
    okcancelhelp $t.f.ok \
	-ok_command "if {\[eval $ok_command $t.lists\] != \"\"} {
    		ListListboxQuit $t
        }" \
	-cancel_command "ListListboxQuit $t" \
	-help_command "show_help gap4 {List-Commands}"
 
    pack $t.f -side bottom -fill x
    pack $t.f.ok -side left -fill x
    pack $t.scrolly -side right -fill y
    pack $t.lists -fill both -expand 1

    wm protocol $t WM_DELETE_WINDOW "ListListboxQuit $t"

    return $t.lists
}

proc UpdateListListbox {t name element op} {
    global NGList NGLists NGSpecial
    global $t.Type

    if {![winfo exists $t]} {
	ListListboxQuit $t
	return
    }

    $t.lists delete 0 end
    foreach i $NGLists {
	if {[set $t.Type] == "read" || [lsearch -exact $NGSpecial $i] == -1} {
	    $t.lists insert end $i
	}
    }
}


#############################################################################
# Main List menu dialogues

#
# Edit a list
#
proc EditList2 {t name} {
    if {![ListExists $name 0]} {bell; return 0}
    destroy $t
    ListEdit $name
}

proc EditListDialog {} {
    global gap_defs

    if {![winfo exists [set w [keylget gap_defs EDIT_LIST.WIN]]]} {
	frame $w
    }

    set t [keylget gap_defs EDIT_LIST.WIN].dialog

    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Edit list"

    entrybox_configure [getLname $t.list \
		 	    [keylget gap_defs EDIT_LIST.LISTNAME] edit] \
	-command "EditList2 $t"

    okcancelhelp $t.ok \
        -ok_command "entrybox_command $t.list.entry" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap4 {List-Commands}" \
        -bd 2 -relief groove

    pack $t.list $t.ok -side top -fill x
} 

#
# Print a list
# 
proc PrintList2 {t list} {
    ListPrint $list
    destroy $t
}

proc PrintListDialog {} {
    global gap_defs

    set t [keylget gap_defs PRINT_LIST.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Print list"

    getLname $t.list [keylget gap_defs PRINT_LIST.LISTNAME] read

    okcancelhelp $t.ok \
        -ok_command "PrintList2 $t \[entrybox_get $t.list.entry\]" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap4 {List-Commands}" \
        -bd 2 -relief groove

    pack $t.list $t.ok -side top -fill both
} 

#
# Load a list from a file
# 
proc LoadList2 {t file name {tag {}}} {
    if {"$file" == ""} {
	bell
	return
    }

    if {"$name" == ""} {set name $file}
    if {$name != "readings"} {
	if {![ListWritable $name]} return
    }
    ListLoad $file $name $tag
    destroy $t
    ListEdit $name
}

proc LoadListDialog {} {
    global gap_defs

    set t [keylget gap_defs LOAD_LIST.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Load list"

    getFname $t.file [keylget gap_defs LOAD_LIST.FILENAME] load
    getLname $t.list [keylget gap_defs LOAD_LIST.LISTNAME] load
    xyn $t.seqid -label "Reading list?"

    okcancelhelp $t.ok \
        -ok_command "LoadList2 $t \[entrybox_get $t.file.entry\] \
		               \[entrybox_get $t.list.entry\] \
			       \[lindex {{} SEQID} \[$t.seqid get\]\]" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap4 {List-Commands}" \
        -bd 2 -relief groove

    pack $t.file $t.list $t.seqid $t.ok -side top -fill both
} 

#
# Save a list to a file
# 
proc SaveList2 {t file name} {
    if {"$file" == ""} {
	bell
	return
    }

    if {"$name" == ""} {set name $file}
    if {[ListSave $file $name]} {destroy $t}
}

proc SaveListDialog {{t {}} {listname {}}} {
    global gap_defs

    if {$t == ""} {
	set t [keylget gap_defs SAVE_LIST.WIN]
    }
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Save list"

    getLname $t.list [keylget gap_defs SAVE_LIST.LISTNAME] read
    getFname $t.file [keylget gap_defs SAVE_LIST.FILENAME] save

    if {$listname != ""} {
	$t.list.entry.entry insert end $listname
    }

    okcancelhelp $t.ok \
        -ok_command "SaveList2 $t \[entrybox_get $t.file.entry\] \
		               \[entrybox_get $t.list.entry\]" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap4 {List-Commands}" \
        -bd 2 -relief groove

    pack $t.list $t.file $t.ok -side top -fill both
} 

proc SaveThisListDialog {listname {t {}}} {
    global gap_defs

    if {$t == ""} {
	set t [keylget gap_defs SAVE_LIST.WIN]
    }
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Save list"

    getFname $t.file [keylget gap_defs SAVE_LIST.FILENAME] save

    okcancelhelp $t.ok \
        -ok_command "SaveList2 $t \[entrybox_get $t.file.entry\] \
		              [list $listname]" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap4 {List-Commands}" \
        -bd 2 -relief groove

    pack $t.file $t.ok -side top -fill both
} 

#
# Delete a list
#
proc DeleteList2 {t list} {
    if {[ListDelete $list]} {
	destroy $t
    }
}

proc DeleteListDialog {} {
    global gap_defs

    set t [keylget gap_defs DELETE_LIST.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Delete list"

    entrybox_configure [getLname $t.list \
		 	    [keylget gap_defs DELETE_LIST.LISTNAME] delete] \
	-command "DeleteList2 $t"

    okcancelhelp $t.ok \
        -ok_command "entrybox_command $t.list.entry" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap4 {List-Commands}" \
        -bd 2 -relief groove

    pack $t.list $t.ok -side top -fill x
} 

#
# Create a list
#
proc CreateList2 {t list} {
    destroy $t
    if {[ListCreate $list ""]} {ListEdit $list}
}

proc CreateListDialog {} {
    global gap_defs

    set t [keylget gap_defs CREATE_LIST.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Create list"

    entrybox_configure [getLname $t.list \
			    [keylget gap_defs CREATE_LIST.LISTNAME] create] \
	-command "CreateList2 $t"


    okcancelhelp $t.ok \
        -ok_command "entrybox_command $t.list.entry" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap4 {List-Commands}" \
        -bd 2 -relief groove

    pack $t.list $t.ok -side top -fill x
}


#
# Copy a list
#
proc CopyList2 {t from to} {
    if {$from == "" || $to == ""} {
	bell
	return
    }

    if {[ListCopy $from $to]} {destroy $t}
}

proc CopyListDialog {} {
    global gap_defs

    set t [keylget gap_defs COPY_LIST.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Copy list"

    getLname $t.from [keylget gap_defs COPY_LIST.FROMLIST] read
    getLname $t.to   [keylget gap_defs COPY_LIST.TOLIST]   create

    okcancelhelp $t.ok \
        -ok_command "CopyList2 $t \[entrybox_get $t.from.entry\] \
		       		  \[entrybox_get $t.to.entry\]" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap4 {List-Commands}" \
        -bd 2 -relief groove

    pack $t.from $t.to $t.ok -side top -fill x
}

#############################################################################
# Some specific list creation functions

#list of left reading names of all contigs
proc CreateAllContigListNumbers { io } {
    global gap_defs

    set clist ""
    set num_contigs [db_info num_contigs $io]
    for {set i 1} {$i <= $num_contigs} {incr i} {
	set order [contig_order_to_number -io $io -order $i]
       	set c [io_read_contig $io $order]
	#set r [io_read_reading $io [keylget c left]]
	set r [keylget c left]
	lappend clist $r
    }
    #puts "LIST $clist"
    return $clist
}

#list of left reading names of all contigs
proc CreateAllContigList { io } {
    global gap_defs

    set clist ""

    set num_contigs [db_info num_contigs $io]
    for {set i 1} {$i <= $num_contigs} {incr i} {
	set order [contig_order_to_number -io $io -order $i]
       	set c [io_read_contig $io $order]
	set r [io_read_reading $io [keylget c left]]
	set contig_name [io_read_text $io [keylget r name]]
	lappend clist $contig_name
    }
    #puts "LIST $clist"
    return $clist
}

#############################################################################
proc InitLists {} {
    global NGSpecial

    set NGSpecial {contigs readings allcontigs allreadings}

    InitReadingList

    ListCreate2 contigs {} SEQID
    ListCreate2 readings {} SEQID
    ListCreate2 allcontigs {} SEQID
    ListCreate2 allreadings {} SEQID

}

proc InitListTrace { } {
    global NGList

    #update displays when list is changed eg editted or deleted
    trace variable NGList(readings) w "UpdateReadingDisplays"
}

#takes a list of contigs numbers created from the contig selector and toggles 
#the selection in NGList(contigs). Returns a list of contig names and toggle
#values (0 to unhighlight, 1 to highlight)
proc UpdateContigList { io list } {
    global c_list NGList NGListTag

    set c_list $NGList(contigs)

    #convert "contigs" list into an array indexed by contig number and
    #containing the position of that contig in the list
    array set cnums {}
    set pos 0
    foreach item $c_list {
	set cnums([db_info get_contig_num $io $item]) $pos
	incr pos
    }

    #search and toggle contig numbers in our array
    foreach item $list {
        if {[info exists cnums($item)]} {
	    unset cnums($item)
        } else {
	    incr pos
	    set cnums($item) -$pos
        }
    }

    # Reverse X(a)=b => X(b)=a
    array set rev {}
    foreach i [array names cnums] {
	set rev($cnums($i)) $i
    }

    # Reconstruct c_list
    set n_list {}
    foreach i [array names rev] {
	set n $rev($i)
	if {$i < 0} {
	    lappend n_list [left_gel $io $rev($i)]
	} else {
	    lappend n_list [lindex $c_list $i]
	}
    }

    ListCreate2 contigs $n_list $NGListTag(contigs)
}

proc UpdateContigListOld { io list } {
    global c_list NGList NGListTag

    set h_list ""
    set n_list ""
    set c_list $NGList(contigs)

    #convert contig list to numbers
    foreach item $c_list {
        lappend n_list [db_info get_contig_num $io $item]
    }
    #search and toggle
    foreach item $list {
        set pos [lsearch $n_list $item]
        set name [left_gel $io $item]
        if {$pos == -1} {
            #not found item in list therefore add it
            lappend n_list $item
            lappend c_list $name
            lappend h_list "$name 1"
        } else {
            #if find item, remove it from the list
            set n_list [lreplace $n_list $pos $pos]
            set c_list [lreplace $c_list $pos $pos]
            lappend h_list "$name 0"
        }
    }
    ListCreate2 contigs $c_list $NGListTag(contigs)
    return $h_list
}


#takes a list of contigs numbers created from the contig selector. Returns a 
#list of contig names and toggle values (0 to unhighlight, 1 to highlight)
proc UpdateContigListNoToggle { io list } {
    global c_list NGList NGListTag

    set h_list ""    
    set n_list ""
    set c_list $NGList(contigs)

    #convert contig list to numbers
    foreach item $c_list {
	lappend n_list [db_info get_contig_num $io $item]
    }
    #search and toggle 
    foreach item $list {
	set pos [lsearch $n_list $item]
	set name [left_gel $io $item]
	if {$pos == -1} {
	    #not found item in list therefore add it
	    lappend n_list $item
	    lappend c_list $name
	    lappend h_list "$name 1"
	} else {
	    #if find item, remove it from the list
	    #set n_list [lreplace $n_list $pos $pos]
	    #set c_list [lreplace $c_list $pos $pos]
	    #lappend h_list "$name 0"
	}
    }
    ListCreate2 contigs $c_list $NGListTag(contigs)
    return $h_list
}

proc GetContigList { } {
    global NGList

    return $NGList(contigs)

}

proc InitReadingList { } {
    global r_list

    set r_list {}

}

#result of trace on NGList(readings)
#HACK with global io but I can't see an easier way of solving this yet
proc UpdateReadingDisplays_disable {} {
    global NGList_urd
    set NGList_urd 0
}

proc UpdateReadingDisplays_enable {} {
    global NGList_urd
    set NGList_urd 1
}

set NGList_urd 1

proc UpdateReadingDisplays {args} {
    global NGList r_list io NGList_urd

    # Optimisation for when this gets triggered LOTS of times, eg
    # when adding many many items to the readings list.
    if {!$NGList_urd} {
	return
    }

    PruneReadingList

    if {![info exists r_list]} {
	set r_list $NGList(readings)
    }

    set h_list ""
    #look for new readings
    foreach r_name $NGList(readings) {
	if {[lsearch -exact $r_list $r_name] == -1} {
	    if {[db_info get_read_num $io $r_name] != -1} {
		lappend h_list "$r_name 1"
	    }
	}
    }
    #look for deleted readings
    foreach r_name $r_list {
	if {[lsearch -exact $NGList(readings) $r_name ] == -1} {
	    if {[db_info get_read_num $io $r_name] != -1} {
	        lappend h_list "$r_name 0"
	    }
	}
    }

    foreach r $h_list {
	reg_notify_highlight -io $io -reading [lindex $r 0] \
		-highlight [lindex $r 1]
	
    }
    set r_list $NGList(readings)

}

# Requires a global io variable
proc PruneReadingList {} {
    global NGList io
    set l $NGList(readings)
    set l2 ""
    foreach i $l {
	if {[db_info get_read_num $io $i] != -1} {
	    lappend l2 $i
	}
    }
    set NGList(readings) $l2
}

proc GetReadingList { } {
    global NGList
    PruneReadingList
    return $NGList(readings)
}

#given a reading name and a highlight status, add or delete from list
#called from contig editor
proc UpdateReadingListItem { r_name highlight} {
    global r_list NGList NGListTag

    if {$highlight} {
	lappend NGList(readings) $r_name
    } else {
	set r_list $NGList(readings)
	set pos [lsearch -exact $r_list $r_name]
	set r_list [lreplace $r_list $pos $pos]
	ListCreate2 readings $r_list $NGListTag(readings)
    }
}

proc SetReadingList { list } {
    global r_list NGList NGListTag

    set r_list $NGList(readings)
    foreach item $list {
	set pos [lsearch -exact $r_list $item]
	if {$pos == -1} {
	    #not found item in list therefore add it
	    lappend r_list $item
	}
    }
    ListCreate2 readings $r_list $NGListTag(readings)
}

#takes a list of readings created from the template display and toggles the
#selection in NGList(readings). Returns a list of reading names and toggle
#values (0 to unhighlight, 1 to highlight)
proc ToggleReadingList { list } {
    global r_list NGList NGListTag

    set h_list ""    
    set r_list $NGList(readings)
    foreach item $list {
	set pos [lsearch -exact $r_list $item]
	if {$pos == -1} {
	    #not found item in list therefore add it
	    lappend r_list $item
	    lappend h_list "$item 1"
	} else {
	    #if find item, remove it from the list
	    set r_list [lreplace $r_list $pos $pos]
	    lappend h_list "$item 0"
	}
    }
    ListCreate2 readings $r_list $NGListTag(readings)
    return $h_list
}

#------------------------------------------------------------------------------

proc ContigsToReadings {io} {
    global gap_defs
    set l [keylget gap_defs CONTIGS_TO_LIST]

    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Contigs to readings"

    contig_id $t.id \
	-state normal\
	-io $io \
	-range 0

    lorf_in $t.infile [keylget l INFILE] \
	"{contig_id_configure $t.id -state disabled} \
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state normal}" -bd 2 -relief groove

    lorf_out $t.outfile [keylget l OUTFILE] \
	{} -bd 2 -relief groove

    okcancelhelp $t.ok \
	-ok_command "ContigsToReadings2 $io $t $t.id $t.infile $t.outfile" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {List-ContigToRead}" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.outfile $t.ok -side top -fill both
}

proc ContigsToReadings2 {io t id infile outfile} {
    if {[lorf_in_get $infile] == 4} {
	set ltmp [contig_id_gel $id]
	SetContigGlobals $io $ltmp
	set list [ListGet "\[$ltmp\]"]
    } elseif {[lorf_in_get $infile] == 3} {
	set list [ListGet "allreadings"]
    } else {
	set ltmp [lorf_get_list $infile]
	set list {}
	foreach i $ltmp {
	    set list "$list [ListGet \[$i\]]"
	}
    }
    if {[set out    [lorf_out_name $outfile]] == ""} {bell; return}
    if {[set format [lorf_out_get  $outfile]] == ""} {bell; return}
    destroy $t
    update idletasks

    ListCreate2 $out $list
    if {"$format" == 2} {
	lorf_out_save $out
    } else {
        ListEdit $out
    }
}

proc MinimalCoverage {io} {
    global gap_defs
    set l [keylget gap_defs MINIMAL_COVERAGE]

    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Minimal coverage"

    contig_id $t.id \
	-state normal\
	-io $io \
	-range 0

    lorf_in $t.infile [keylget gap_defs MINIMAL_COVERAGE.INFILE] \
	"{contig_id_configure $t.id -state disabled} \
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state normal}" -bd 2 -relief groove

    lorf_out $t.outfile [keylget gap_defs MINIMAL_COVERAGE.OUTFILE] \
	{} -bd 2 -relief groove

    okcancelhelp $t.ok \
	-ok_command "MinimalCoverage2 $io $t $t.id $t.infile $t.outfile" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {List-MinCoverage}" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.outfile $t.ok -side top -fill both
}

proc MinimalCoverage2 {io t id infile outfile} {
    if {[lorf_in_get $infile] == 4} {
	set list [contig_id_gel $id]
	SetContigGlobals $io $list
    } elseif {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $infile]
    }
    if {[set out    [lorf_out_name $outfile]] == ""} {bell; return}
    if {[set format [lorf_out_get  $outfile]] == ""} {bell; return}
    destroy $t
    update idletasks

    set r [minimal_coverage -io $io -contigs $list]

    ListCreate2 $out $r
    if {"$format" == 2} {
	lorf_out_save $out
    } else {
        ListEdit $out
    }
}

proc UnattachedReadings {io} {
    global gap_defs
    set l [keylget gap_defs UNATTACHED_READINGS]

    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Unattached readings"

    lorf_out $t.outfile [keylget gap_defs UNATTACHED_READINGS.OUTFILE] \
	{} -bd 2 -relief groove

    okcancelhelp $t.ok \
	-ok_command "UnattachedReadings2 $io $t $t.outfile" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {List-Unattached}" \
	-bd 2 -relief groove

    pack $t.outfile $t.ok -side top -fill both
}

proc UnattachedReadings2 {io t outfile} {
    if {[set out    [lorf_out_name $outfile]]   == ""} {bell; return}
    if {[set format [lorf_out_get  $outfile]] == ""} {bell; return}

    destroy $t
    update idletasks

    set r [unattached_readings -io $io]

    ListCreate2 $out $r
    if {"$format" == 2} {
	lorf_out_save $out
    } else {
        ListEdit $out
    }
}
