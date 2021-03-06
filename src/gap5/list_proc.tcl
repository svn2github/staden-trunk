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
	set crec [db_info get_contig_num $io $id]
	if {$crec == "" || $crec <= 0} {
	    return ""
	}
	set c [$io get_contig $crec]
	foreach x [$c seqs_in_range [$c get_start] [$c get_end]] {
	    lappend l "#[lindex $x 2]"
	}
	$c delete
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

proc ListDelete {name {parent .}} {
    global NGList NGLists NGSpecial NGListTag NGList_traced

    if {[ListExists2 $name]} {
	if {[info exists NGSpecial] && [lsearch -exact $NGSpecial $name] != -1 } {
	    tk_messageBox -icon error -type ok -title "Special list" \
		-message "This list cannot be deleted" \
		-parent $parent
	    return 0
	}

	unset NGList($name)
	unset NGListTag($name)
	global NGList_read_hash_$name
	catch {unset NGList_read_hash_$name}
	set NGLists [lreplace $NGLists [lsearch -exact $NGLists $name] [lsearch -exact $NGLists $name]]

	if {[info exists NGList_traced($name)]} {
	    unset NGList_traced($name)
	    trace remove variable NGList_traced($name) write [list UpdateReadingDisplays $name]
	}
	return 1
    } else {
	tk_messageBox -icon error -type ok -title "No such list" \
	    -message "List does not exist" \
	    -parent $parent
	return 0
    }
}

proc ListCopy {from to {parent .}} {
    global NGList NGLists

    if {$from == $to} {
	tk_messageBox -icon error -type ok -title "Copy list" \
	    -message "Source and destination list names need to be different" \
	    -parent $parent
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

proc ListWritable {name {parent .}} {
    global NGSpecial

    if {![ListNameValid $name $parent]} {
	return 0
    }

    if {[ListExists2 $name]} {

	if {[lsearch -exact $NGSpecial $name] != -1 } {
	    tk_messageBox -icon error -type ok -title "Special list" \
		    -message "This list cannot be replaced" \
		    -parent $parent
	    return 0
	}

	set answer [tk_messageBox -icon question -type yesno \
			-title "List exists" \
			-message "List already exists. Replace?" \
			-parent $parent]
	case $answer {
	    no {return 0} 
	    yes {ListDelete $name $parent}
	}
    }

    return 1
}

proc ListNameValid {name {parent .}} {
    if {"$name" == ""} {
	tk_messageBox -icon error -type ok -title "No list name" \
	    -message "You have not entered a list name"\
	    -parent $parent
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
    return [ListWritable [$path.entry get] $path]
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
    global gap5_defs
    set ed_list [keylget gap5_defs EDIT_LIST.WIN]
    global $ed_list
    upvar \#0 $ed_list el_array
    
    if {![array exists el_array] || ![info exists el_array($name)]} {
	return ""
    }
    
    if {![winfo exists $el_array($name)]} {
	# Someone closed the window behind our backs
	unset el_array($name)
	return ""
    }

    return $el_array($name)
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
	#foreach i $NGList($name) {
	#    $t insert end "$i\n"
	#}
	$t insert end [join $NGList($name) "\n"] "" "\n"
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

proc ListEditClose { t name { save 1 } } {
    global gap5_defs NGList
    set ed_list [keylget gap5_defs EDIT_LIST.WIN]
    global $ed_list
    upvar \#0 $ed_list el_array

    trace remove variable NGList($name) write "ListEditUpdate $t.list [list $name]"

    if { $save } { ListEditSave [list $name] }
    unset el_array($name)

    destroy $t
}

proc ListEdit {name {read_only 0}} {
    global gap5_defs NGList NGSpecial tcl_platform
    
    if {![ListExists $name 0]} {bell; return 0}

    set ed_list [keylget gap5_defs EDIT_LIST.WIN]
    global $ed_list
    upvar \#0 $ed_list el_array

    # We cannot edit the special lists
    if {[lsearch -exact $NGSpecial $name] != -1} {
	set read_only 1
    }

    # Create the new editor if it doesn't exist
    if {[set t [ListEditExists $name]] == ""} {
        set t "${ed_list}_[ListNextEditor]"
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
	    -command "ListEditClose $t [list $name]"
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

	wm protocol $t WM_DELETE_WINDOW "ListEditClose $t [list $name]"
    
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
	trace variable NGList($name) u "ListEditClose $t [list $name] 0"

	# Put the window in our array of list windows
	set el_array($name) $t
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
	-help_command "show_help gap5 {List-Commands}"
 
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
    global NGListTag

    if {![ListExists $name 0]} {bell; return 0}
    destroy $t

    if {[yes_no_get $t.extended] \
	    && $name ne "contigs" \
	    && $name ne "allcontigs" \
	    && $NGListTag($name) != ""} {
	ListEditMulti $name
    } else {
	ListEdit $name
    }
}

proc EditListDialog {} {
    global gap5_defs

    set t [keylget gap5_defs EDIT_LIST.WIN]_dialog

    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Edit list"

    entrybox_configure [getLname $t.list \
		 	    [keylget gap5_defs EDIT_LIST.LISTNAME] edit] \
	-command "EditList2 $t"

    yes_no $t.extended \
	-title "Extended reading list view" \
	-orient horizontal \
	-bd 0 \
	-default 0

    okcancelhelp $t.ok \
        -ok_command "entrybox_command $t.list.entry" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap5 {List-Commands}" \
        -bd 2 -relief groove

    pack $t.list $t.extended $t.ok -side top -fill x
} 

#
# Print a list
# 
proc PrintList2 {t list} {
    ListPrint $list
    destroy $t
}

proc PrintListDialog {} {
    global gap5_defs

    set t [keylget gap5_defs PRINT_LIST.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Print list"

    getLname $t.list [keylget gap5_defs PRINT_LIST.LISTNAME] read

    okcancelhelp $t.ok \
        -ok_command "PrintList2 $t \[entrybox_get $t.list.entry\]" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap5 {List-Commands}" \
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
	if {![ListWritable $name $t]} return
    }
    ListLoad $file $name $tag
    destroy $t
    ListEdit $name
}

proc LoadListDialog {{t {}} {listname {}}} {
    global gap5_defs

    if {$t == ""} {
	set t [keylget gap5_defs LOAD_LIST.WIN]
    }
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Load list"

    getFname $t.file [keylget gap5_defs LOAD_LIST.FILENAME] load
    getLname $t.list [keylget gap5_defs LOAD_LIST.LISTNAME] load
    xyn $t.seqid -label "Reading list?"

    if {$listname != ""} {
	$t.list.entry.entry insert end $listname
    }

    okcancelhelp $t.ok \
        -ok_command "LoadList2 $t \[entrybox_get $t.file.entry\] \
		               \[entrybox_get $t.list.entry\] \
			       \[lindex {{} SEQID} \[$t.seqid get\]\]" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap5 {List-Commands}" \
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
    global gap5_defs

    if {$t == ""} {
	set t [keylget gap5_defs SAVE_LIST.WIN]
    }
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Save list"

    getLname $t.list [keylget gap5_defs SAVE_LIST.LISTNAME] read
    getFname $t.file [keylget gap5_defs SAVE_LIST.FILENAME] save

    if {$listname != ""} {
	$t.list.entry.entry insert end $listname
    }

    okcancelhelp $t.ok \
        -ok_command "SaveList2 $t \[entrybox_get $t.file.entry\] \
		               \[entrybox_get $t.list.entry\]" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap5 {List-Commands}" \
        -bd 2 -relief groove

    pack $t.list $t.file $t.ok -side top -fill both
} 

proc SaveThisListDialog {listname {t {}}} {
    global gap5_defs

    if {$t == ""} {
	set t [keylget gap5_defs SAVE_LIST.WIN]
    }
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Save list"

    getFname $t.file [keylget gap5_defs SAVE_LIST.FILENAME] save

    okcancelhelp $t.ok \
        -ok_command "SaveList2 $t \[entrybox_get $t.file.entry\] \
		              [list $listname]" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap5 {List-Commands}" \
        -bd 2 -relief groove

    pack $t.file $t.ok -side top -fill both
} 

#
# Delete a list
#
proc DeleteList2 {t list} {
    if {[ListDelete $list $t]} {
	destroy $t
    }
}

proc DeleteListDialog {} {
    global gap5_defs

    set t [keylget gap5_defs DELETE_LIST.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Delete list"

    entrybox_configure [getLname $t.list \
		 	    [keylget gap5_defs DELETE_LIST.LISTNAME] delete] \
	-command "DeleteList2 $t"

    okcancelhelp $t.ok \
        -ok_command "entrybox_command $t.list.entry" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap5 {List-Commands}" \
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

proc CreateListDialog {{wait 0}} {
    global gap5_defs

    set t [keylget gap5_defs CREATE_LIST.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Create list"

    entrybox_configure [getLname $t.list \
			    [keylget gap5_defs CREATE_LIST.LISTNAME] create] \
	-command "CreateListDialog2 $t"


    okcancelhelp $t.ok \
        -ok_command "entrybox_command $t.list.entry" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap5 {List-Commands}" \
        -bd 2 -relief groove

    pack $t.list $t.ok -side top -fill x
    if {$wait} {
	global $t._list
	vwait $t._list
	return [set $t._list]
    }
}

proc CreateListDialog2 {t list} {
    global $t._list
    set $t._list $list
    CreateList2 $t $list
}


#
# Copy a list
#
proc CopyList2 {t from to} {
    if {$from == "" || $to == ""} {
	bell
	return
    }

    if {[ListCopy $from $to $t]} {destroy $t}
}

proc CopyListDialog {} {
    global gap5_defs

    set t [keylget gap5_defs COPY_LIST.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Copy list"

    getLname $t.from [keylget gap5_defs COPY_LIST.FROMLIST] read
    getLname $t.to   [keylget gap5_defs COPY_LIST.TOLIST]   create

    okcancelhelp $t.ok \
        -ok_command "CopyList2 $t \[entrybox_get $t.from.entry\] \
		       		  \[entrybox_get $t.to.entry\]" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap5 {List-Commands}" \
        -bd 2 -relief groove

    pack $t.from $t.to $t.ok -side top -fill x
}

#############################################################################
# Some specific list creation functions

#list of record numbers of all contigs
proc CreateAllContigListNumbers { io } {
    global gap5_defs

    set clist ""
    set num_contigs [$io num_contigs]
    for {set i 0} {$i < $num_contigs} {incr i} {
	set order [contig_order_to_number -io $io -order $i]
	lappend clist $order
    }
    return $clist
}

#list of record numbers in "=num" syntax of all contigs
proc CreateAllContigList=Numbers { io } {
    global gap5_defs

    set clist ""
    set num_contigs [$io num_contigs]
    for {set i 0} {$i < $num_contigs} {incr i} {
	set order [contig_order_to_number -io $io -order $i]
	lappend clist =$order
    }
    return $clist
}

#list of names of all contigs
proc CreateAllContigList { io } {
    global gap5_defs

    set clist ""

    set num_contigs [$io num_contigs]
    for {set i 0} {$i < $num_contigs} {incr i} {
	set order [contig_order_to_number -io $io -order $i]
	set c [$io get_contig $order]
	lappend clist [$c get_name]
	$c delete
    }
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

proc InitListTrace { {list readings} } {
    global NGList NGList_traced

    if {![info exists NGList_traced]} {array set NGList_traced {}}

    if {![info exists NGList_traced($list)]} {
	#update displays when list is changed eg editted or deleted
	trace variable NGList($list) w [list UpdateReadingDisplays $list]
	set NGList_traced($list) 1
    }

    # While fun, tracing each independent set and unset within the array
    # has the potential to spam events if we load an entire new list from 
    # a file or hit the Clear button.
    # So for now we just disable this feature until we can find an
    # efficient route, if any exists.

    #global NGList_read_hash
    #trace variable NGList_read_hash wu "UpdateReadingDisplaysHash"
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
	    set c [$io get_contig $rev($i)]
	    lappend n_list [$c get_name]
	    $c delete
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
	set c [$io get_contig $item]
	set name [$c get_name]
	$c delete
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
    global NGList_read_hash_readings
    array set NGList_read_hash_readings ""

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

# Converts from NGList(readings) text contents to NGList_read_hash.
# This automatically generates reg_notify_highlight events for changes
# to the hash.
proc UpdateReadingDisplays {list args} {
    global NGList io NGList_urd NGList_read_hash_$list  NGList_size
    upvar #0 NGList_read_hash_$list hash

    # Optimisation for when this gets triggered LOTS of times, eg
    # when adding many many items to the readings list.
    if {!$NGList_urd} {
	return
    }

    PruneReadingList $list

    array set tmp ""
    foreach r_name $NGList($list) {
	set tmp($r_name) 1
    }

    contig_notify -io $io -type BUFFER_START -cnum 0 -args {}

    #look for new readings
    foreach r_name [array names tmp] {
	if {![info exists hash($r_name)]} {
	    if {[regexp {^\#(\d+)} $r_name _ rec]} {
		if {[$io rec_exists 18 $rec] == 0} {
		    verror ERR_WARN UpdateReadingDisplays \
			"Sequence $r_name is invalid"
		    continue
		}
	    }
	    set hash($r_name) ""
	}
    }

    #look for deleted readings
    foreach r_name [array names hash] {
	if {![info exists tmp($r_name)]} {
	    unset hash($r_name)
	}
    }

    set NGList($list) [array names hash]
    set NGList_size($list) [llength $NGList($list)]

    contig_notify -io $io -type BUFFER_END -cnum 0 -args {}

    # Send a catch-all notification to indicate that the reading list
    # has changed, but not precisely what changed. This avoids an endless
    # stream of events for windows that only just care about whether they
    # should redraw.
    contig_notify -io $io -type HIGHLIGHT_READ -cnum 0 -args {}
}

# Returns the number of elements in a list
proc ListSize {list} {
    global NGList
    return [llength $NGList($list)]
}

proc UpdateReadingDisplaysHash {var subvar op} {
    global io

    if {$op == "w"} {
	reg_notify_highlight -io $io -reading $subvar -highlight 1
    } elseif {$op == "u"} {
	reg_notify_highlight -io $io -reading $subvar -highlight 0
    }
}

# Requires a global io variable
proc PruneReadingList { {list readings} } {
    global NGList io
    set l $NGList($list)
    set l2 ""
    foreach i $l {
	if {[db_info get_read_num $io $i] != -1} {
	    lappend l2 $i
	}
    }
    set NGList($list) $l2
}

proc GetReadingList { } {
    global NGList
    PruneReadingList
    return $NGList(readings)
}

# Given a reading name and a highlight status, add or delete from list
# called from contig editor.
# If $pair == 1 then we also add the reading mate-pair to the list too
#
#Returns 1 if last item is selected, 0 if unselected
proc UpdateReadingListItem { list r_names highlight {pair 0}} {
    global r_list NGList NGListTag NGList_read_hash_$list
    upvar #0 NGList_read_hash_$list hash
    global io

    set plist {}
    if {$pair} {
	foreach r_name $r_names {
	    set rec [db_info get_read_num $io $r_name]
	    if {$rec > 0} {
		set s [$io get_seq $rec]
		set rec [$s get_pair]
		if {$rec > 0} {
		    lappend plist "#$rec"
		}
		$s delete
	    }
	}
    }

    foreach r_name "$plist $r_names" {
	if {$highlight == -1} {
	    # Toggle
	    if {[info exists hash($r_name)]} {
		unset hash($r_name)
		set r 0
	    } else {
		set hash($r_name) ""
		set r 1
	    }
	} else {
	    if {$highlight == 1} {
		if {![info exists hash($r_name)]} {
		    set hash($r_name) ""
		}
		set r 1
	    } else {
		if {[info exists hash($r_name)]} {
		    unset hash($r_name)
		}
		set r 0
	    }
	}
    }

    ListCreate2 $list [array names hash] $NGListTag($list)

#    if {$highlight} {
#	lappend NGList(readings) $r_name
#    } else {
#	set r_list $NGList(readings)
#	set pos [lsearch -exact $r_list $r_name]
#	set r_list [lreplace $r_list $pos $pos]
#	ListCreate2 readings $r_list $NGListTag(readings)
#    }

    return $r
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
    global gap5_defs
    set l [keylget gap5_defs CONTIGS_TO_LIST]

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
	-help_command "show_help gap5 {List-ContigToRead}" \
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

    ListCreate2 $out $list SEQID
    if {"$format" == 2} {
	lorf_out_save $out
    } else {
        ListEdit $out
    }
}

#------------------------------------------------------------------------------

proc MinimalCoverage {io} {
    global gap5_defs
    set l [keylget gap5_defs MINIMAL_COVERAGE]

    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Minimal coverage"

    contig_id $t.id \
	-state normal\
	-io $io \
	-range 0

    lorf_in $t.infile [keylget gap5_defs MINIMAL_COVERAGE.INFILE] \
	"{contig_id_configure $t.id -state disabled} \
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state normal}" -bd 2 -relief groove

    lorf_out $t.outfile [keylget gap5_defs MINIMAL_COVERAGE.OUTFILE] \
	{} -bd 2 -relief groove

    okcancelhelp $t.ok \
	-ok_command "MinimalCoverage2 $io $t $t.id $t.infile $t.outfile" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 {List-MinCoverage}" \
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

#------------------------------------------------------------------------------

proc UnattachedReadings {io} {
    global gap5_defs
    set l [keylget gap5_defs UNATTACHED_READINGS]

    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Unattached readings"

    lorf_out $t.outfile [keylget gap5_defs UNATTACHED_READINGS.OUTFILE] \
	{} -bd 2 -relief groove

    okcancelhelp $t.ok \
	-ok_command "UnattachedReadings2 $io $t $t.outfile" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 {List-Unattached}" \
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

#-----------------------------------------------------------------------------

proc PairReadings {io} {
    global gap5_defs
    set l [keylget gap5_defs PAIR_READINGS]

    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Pair readings"

    lorf_in $t.infile [keylget l INFILE] {} -bd 2 -relief groove

    lorf_out $t.outfile [keylget l OUTFILE] \
	{} -bd 2 -relief groove

    okcancelhelp $t.ok \
	-ok_command "PairReadings2 $io $t $t.infile $t.outfile" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 {List-Pair}" \
	-bd 2 -relief groove

    pack $t.infile $t.outfile $t.ok -side top -fill both
}

proc PairReadings2 {io t infile outfile} {
    set list [lorf_get_list $infile]
    if {$list == ""} {bell; return}

    if {[set out    [lorf_out_name $outfile]]   == ""} {bell; return}
    if {[set format [lorf_out_get  $outfile]] == ""} {bell; return}

    destroy $t
    update idletasks

    set rlist {}
    foreach rec [pair_readings -io $io -readings $list] {
	lappend rlist "#$rec"
    }

    ListCreate2 $out $rlist SEQID
    if {"$format" == 2} {
	lorf_out_save $out
    } else {
        ListEdit $out
    }
}

#-----------------------------------------------------------------------------
# Convert a list of #num to name
proc NumToName {io} {
    global gap5_defs
    set l [keylget gap5_defs NUMLIST_TO_NAME]

    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Reading numbers to names"

    getLname $t.from [keylget gap5_defs NUMLIST_TO_NAME.FROMLIST] read
    getLname $t.to   [keylget gap5_defs NUMLIST_TO_NAME.TOLIST]   create

    okcancelhelp $t.ok \
	-ok_command "NumToName2 $io $t $t.from $t.to" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 {List-NumToName}" \
	-bd 2 -relief groove

    pack $t.from $t.to $t.ok -side top -fill both
}


proc NumToName2 {io t from to} {
    set from [entrybox_get $from.entry]
    set to   [entrybox_get $to.entry]

    if {$from == "" || $to == ""} {
	bell
	return
    }
    set parent [winfo parent $t]
    destroy $t

    set l {}
    set nerrs 0
    foreach n [ListGet $from] {
	set num [lindex [regexp -inline {^\#?([0-9]+)$} $n] 1]
	if {$num} {
	    if {[catch {set s [$io get_seq $num]}]} {
		set s ""
	    }

	    if {$s != ""} {
		set name [$s get_name]
		$s delete
	    } else {
		set num ""
	    }
	}

	if {$num == ""} {
	    tk_messageBox -icon error -type ok -title "List conversion" \
		-message "Failure to convert seq number '$n' to a name" \
		-parent $parent
	    return
	} 

	lappend l $name
    }

    ListCreate2 $to $l SEQID
    ListEdit $to
}

#-----------------------------------------------------------------------------
# Upgraded multi-column reading list
# Does an editor for this list exist already?
proc ListEditMultiExists {name} {
    global gap5_defs
    set ed_list_multi [keylget gap5_defs EDIT_LIST.WIN]_multi
    global $ed_list_multi
    upvar \#0 $ed_list_multi elm_array

    if {![array exists elm_array] || ![info exists elm_array($name)]} {
	return ""
    }

    if {![winfo exists $elm_array($name)]} {
	# Someone closed the window behind our backs
	unset elm_array($name)
	return ""
    }

    return $elm_array($name)
}

proc ListEditMultiClose { t name } {
    global gap5_defs
    set ed_list_multi [keylget gap5_defs EDIT_LIST.WIN]_multi
    global $ed_list_multi
    upvar \#0 $ed_list_multi elm_array

    unset elm_array($name)

    destroy $t
}

proc ListEditMulti {name {read_only 0}} {
    global gap5_defs NGList NGSpecial tcl_platform io
    
    if {![ListExists $name 0]} {bell; return 0}

    set ed_list_multi [keylget gap5_defs EDIT_LIST.WIN]_multi
    global $ed_list_multi
    upvar \#0 $ed_list_multi elm_array

#    # We cannot edit the special lists
#    if {[lsearch -exact $NGSpecial $name] != -1} {
#	set read_only 1
#    }
    set read_only 1

    # Create the new editor if it doesn't exist
    if {[set t [ListEditMultiExists $name]] == ""} {
        set t "${ed_list_multi}_[ListNextEditor]"
        xtoplevel $t
	wm geometry $t 990x400
	if {$read_only} {
	    wm title $t "Viewing list: \"$name\""
	} else {
	    wm title $t "Editing list: \"$name\""
	}
    
        global $t.CurList
        set $t.CurList $name
        
       	#wm protocol $t WM_DELETE_WINDOW "ListEditSave [list $name]; destroy $t"
	wm protocol $t WM_DELETE_WINDOW "ListEditMultiClose $t [list $name]"

        # The text display and it's scrollbars.
	tablelist $t.list \
	    -columns {20 Name 10 Number 15 Contig 10 Position 15 Scaffold 10 {Pair Num} 15 {Pair Contig} 10 {Pair Pos} 15 {Pair Scaf}} \
	    -labelcommand tablelist::sortByColumn \
	    -exportselection 0 \
	    -selectmode extended \
	    -yscrollcommand [list $t.scrolly set]
	
	$t.list columnconfigure 0 -bg #e0e0e0
	$t.list columnconfigure 1 -sortmode integer
	$t.list columnconfigure 3 -sortmode integer
	$t.list columnconfigure 4 -sortmode command \
	    -formatcommand ListContigsScaffoldFormat \
	    -sortcommand [list ListContigsScaffoldSort $t.list]
	$t.list columnconfigure 5 -bg #e0e0e0 -sortmode command \
	    -sortcommand SortOptionalInteger
	$t.list columnconfigure 6 -bg #e0e0e0
	$t.list columnconfigure 7 -bg #e0e0e0 -sortmode command \
	    -sortcommand SortOptionalInteger
	$t.list columnconfigure 8 -bg #e0e0e0 -sortmode command \
	    -formatcommand ListContigsScaffoldFormat \
	    -sortcommand [list ListContigsScaffoldSort $t.list]

        scrollbar $t.scrolly   -orient vert  -command "$t.list yview"

        pack $t.scrolly -side right -fill y
        pack $t.list -expand 1 -fill both
    
	# Monitor changes to this list
	trace variable NGList($name) w "ListEditMultiUpdate $t.list [list $name]"
	trace variable NGList($name) u "ListEditMultiClose $t [list $name]"

	# The command bar.
	set bar [frame $t.bar]
	pack $bar -side bottom -fill both

	# Output list name + size/menu
	set lf [frame $bar.name -bd 0]
	xcombobox $lf.entry \
	    -textvariable ${t}(OutputList) \
	    -command ListEditMultiList \
	    -postcommand "editor_olist_populate $lf.entry" \
	    -label "Output list:" \
	    -default $name \
	    -width 15
	bind $lf.entry <Any-Leave> {
	    ListEditMultiList %W [%W get]
	    #focus [winfo toplevel %W]
	}

	label $lf.size -text "  Size:"
	menubutton $lf.menu \
	    -textvariable ${t}(ListSize) \
	    -menu $lf.menu.opts \
	    -indicatoron 1
	menu $lf.menu.opts -tearoff 0
	$lf.menu.opts add command \
	    -label Save \
	    -command "SaveListDialog $bar.list_sv \[set ${t}(OutputList)\]"
	$lf.menu.opts add command \
	    -label Load \
	    -command "LoadListDialog $bar.list_ld \[set ${t}(OutputList)\]"
	$lf.menu.opts add command \
	    -label Clear \
	    -command "ListClear \[set ${t}(OutputList)\]"
	$lf.menu.opts add command \
	    -label Delete \
	    -command "ListDelete \[set ${t}(OutputList)\]"
	
	pack $lf.entry $lf.size $lf.menu -side left -fill both

	# Which end to output
	set ef [frame $bar.end -bd 0]
	radiobutton $ef.read \
	    -text "Read" \
	    -variable ${t}(End) \
	    -value 1
	radiobutton $ef.pair \
	    -text "Pair" \
	    -variable ${t}(End) \
	    -value 2
	radiobutton $ef.both \
	    -text "Both" \
	    -variable ${t}(End) \
	    -value 3
	global $t
	set ${t}(End) 3
	pack $ef.read $ef.pair $ef.both -side left -fill both

	# Add, remove, set commands
	button $bar.add -text "Add selection" \
	    -command "ListEditMultiSave $t add"
	button $bar.del -text "Remove selection" \
	    -command "ListEditMultiSave $t del"
	button $bar.set -text "Set to selection" \
	    -command "ListEditMultiSave $t set"
	
	pack $lf $ef $bar.add $bar.del $bar.set -side left -expand 1
	
	# Put the window in our array of multi-list windows
	set elm_array($name) $t
    } else {
        raise $t
    }

    # But we always update the list
    ListEditMultiUpdate $t.list $name

    bind [$t.list bodypath] <<menu>> \
	"ListEditMultiMenu $io $t.list %x %y %X %Y"

    ListEditMultiList $t.list $name

    return 1
}

proc SortOptionalInteger {n1 n2} {
    if {$n1 == "" && $n2 == ""} {
	return 0
    } elseif {$n1 == ""} {
	return -1
    } elseif {$n2 == ""} {
	return 1
    } else {
	return [expr {$n1 - $n2}]
    }
}

# Tablelist update
proc ListEditMultiUpdate {t name args} {
    global NGList NGListTag io

    if {![winfo exists $t]} {
        global $t.CurList
        trace vdelete NGLists w "UpdateListListbox $t"
	if {[info exists $t.CurList]} {
	    unset $t.CurList
	}
	return
    }

    $t delete 0 end
    if {$NGListTag($name) != ""} {
	if {$NGListTag($name) == "SEQID"} {
	    foreach i $NGList($name) {
		regexp {([[:space:]]*)([^[:space:]]*)(.*)} $i _ l m r
		set rec [db_info get_read_num $io $m]
		if {[$io rec_exists 18 $rec] == 0} {
		    $t insert end "$i"
		    continue
		}
		set s [$io get_sequence $rec]
		set sname [$s get_name]
		set pos [$s get_position]
		set crec [$s get_contig]
		set c [$io get_contig $crec]
		set cname [$c get_name]
		set frec [$c get_scaffold]
		if {$frec != 0} {
		    set f [$io get_scaffold $frec]
		    set ctgs [$f get_contigs]
		    set fname "[$f get_name]/[lsearch $ctgs $crec]"
		    $f delete
		} else {
		    set fname "(none)/0"
		}
		$c delete

		foreach {pair_rec pair_contig pair_start pair_end} \
		    [$s get_pair_pos] break
		if {$pair_rec} {
		    set c [$io get_contig $pair_contig]
		    set pcname [$c get_name]
		    set frec [$c get_scaffold]
		    if {$frec != 0} {
			set f [$io get_scaffold $frec]
			set ctgs [$f get_contigs]
			set pfname "[$f get_name]/[lsearch $ctgs $crec]"
			$f delete
		    } else {
			set pfname "(none)/0"
		    }
		    $c delete
		    
		    $t insert end [list $sname $rec $cname $pos $fname \
				      $pair_rec $pcname $pair_start $pfname]
		} else {
		    $t insert end [list $sname $rec $cname $pos $fname]
		}
		$s delete
	    }
	} else {
	    foreach i $NGList($name) {
		$t insert end "$i"
	    }
	}
    } else {
	foreach i $NGList($name) {
	    $t insert end "$i\n"
	}
    }

    # Incase same list is both input and output
    #    ListEditMultiList $t $name
}

proc ListEditMultiMenu {io w x y X Y} {
    # Compute row and column clicked on
    set l [$w bodypath] 
    set x [expr {$x + [winfo x $l]}]
    set y [expr {$y + [winfo y $l]}]
    foreach {row col} [split [$w nearestcell $x $y] ,] {break}

    # puts row=$row,col=$col

    # Find the contig identifier
    if {$col >= 5} {
	set rec [lindex [$w get $row] 5]
	set pos [lindex [$w get $row] 7]
    } else {
	set rec [lindex [$w get $row] 1]
	set pos [lindex [$w get $row] 3]
    }
    set contig [db_info get_contig_num $io #$rec]
    # puts $rec/[$w get $row]

    if { $contig <= 0 } {
	# Not a valid read
	return
    }

    create_popup $w.m "Commands (\#$rec)"
    $w.m add command -label "Edit Contig" \
	-command "edit_contig -io $io -contig $contig -reading #$rec"
    $w.m add command -label "Template Display" \
	-command "CreateTemplateDisplay $io $contig $pos"

    tk_popup $w.m [expr $X-20] [expr $Y-10]
}

proc ListEditMultiSave {w op} {
    global NGList
    set t [winfo toplevel $w]
    upvar #0 $t opt

    set name [$w.bar.name.entry get]
    set reads ""
    if {[llength [$w.list curselection]] <= 1} {
	set l [list [$w.list get [$w.list curselection]]]
    } else {
	set l [$w.list get [$w.list curselection]]
    }
    foreach r $l {
	if {$opt(End) & 1} {
	    set val [lindex $r 1]
	    if { $val != "" } {
		lappend reads "#$val"
	    }
	}
	if {$opt(End) & 2} {
	    set val [lindex $r 5]
	    if { $val != ""} {
		lappend reads "#$val"
	    }
	}
    }

    switch $op {
	"add" {
	    append NGList($name) " $reads"
	    ListRemoveDuplicates $name
	}

	"del" {
	    set n $NGList($name)
	    foreach r $reads {
		if {[set ind [lsearch -exact $n $r]] != -1} {
		    set n [lreplace $n $ind $ind]
		}
	    }
	    set NGList($name) $n
	}

	"set" {
	    set NGList($name) $reads
	}
    }
}

proc ListEditMultiList {w l} {
    set w [winfo toplevel $w]
    upvar #0 $w opt

    global NGList
    if {![ListExists2 $l]} {
	ListCreate2 $l {} SEQID
    }
    InitListTrace $l

    global NGList_size
    set opt(OutputList) $l
    set opt(ListSize) [ListSize $opt(OutputList)]
    catch {set NGList_size($opt(OutputList)) $opt(ListSize)}

    $w.bar.name.menu configure -textvariable NGList_size($opt(OutputList))
}
