#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995-1997. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

# Catch as we need to use this in non-X mode.
catch {package require tablelist}
catch {namespace import tablelist::tablelist}

##############################################################################
#set the active tag defaults
proc TagDialog {name path command {version {}}} {
    global NGTag_$name$version

    if {[xtoplevel $path -resizable 1] == ""} {return}
    wm title $path "Active tags"

    set arr [array get NGTag_$name$version]
    tag_checklist $path.check $arr -selectmode extended -width 40 -height 20

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $path.but \
	-bd 2 \
	-relief groove \
	-ok_command "TagDialog_ok $path $name [list $command] \
		      [list $version]" \
	-asdefault_command "TagDialog_default $path $name" \
	-cancel_command "destroy $path $name" \
	-help_command "show_help gap5 {Conf-Tag}"
    
    pack $path.check -fill both -expand 1
    pack $path.but -fill both
}

proc TagDialog_ok {path name command version} {
    global NGTag_$name$version

    array set NGTag_$name$version [$path.check list2array]
    destroy $path

    if {$command != {}} {
	eval $command
    }
}    

proc TagDialog_default {path name} {
    global NGTag gap5_defs env

    array set tmp [$path.check list2array]
    set active_tags ""
    for {set i 0} {$i < $NGTag(num_tags)} {incr i} {
        if {$tmp($i)} {
            lappend active_tags $NGTag($i,tagid)
        }
    }

    keylset gap5_defs $name $active_tags
    update_defs gap5_defs $env(HOME)/.gap5rc $name
}

##############################################################################
#sets up a local version of the active tag list
proc SetDefaultTags {name {version {}}} {
    global gap5_defs NGTag_$name$version NGTag_num_tags${name}${version} NGTag
    
    # Initialise all tags to off by default
    set list ""
    for {set i 0} {$i < $NGTag(num_tags)} {incr i} {
	set NGTag_${name}${version}($i) 0
	lappend list $NGTag($i,tagid)
    }

    set NGTag_num_tags${name}${version} $NGTag(num_tags)

    foreach i [keylget gap5_defs $name] {
	if {$i == {*}} {
	    # Add all tags
	    for {set j 0} {$j < $NGTag(num_tags)} {incr j} {
		set NGTag_${name}${version}($j) 1
	    }
	} else {
	    if {[string match -* $i]} {
		set c [string trimleft $i -]
		if {[set ind [lsearch -exact $list $c]] != -1} {
		    set NGTag_${name}${version}($ind) 0
		}
	    } else {
		set c [string trimleft $i +]
		if {[set ind [lsearch -exact $list $c]] != -1} {
		    set NGTag_${name}${version}($ind) 1
		}
	    }
	}
    }
}


##############################################################################
#Gets the default tag list
proc GetDefaultTags {name {version {}}} {
    global NGTag NGTag_$name$version NGTag_num_tags${name}${version} gap5_defs

    set active_tags ""

    for {set i 0} {$i < [set NGTag_num_tags${name}${version}]} {incr i} {
        if {[set NGTag_${name}${version}($i)]} {
            lappend active_tags $NGTag($i,tagid)
        }
    }

    return $active_tags
}


##############################################################################
#
# tag_checklist:
#
# Creates a tag checklist frame containing a listbox of all tags types.
# Extra command line arguments are as per the frame widget.
#
proc tag_checklist {path array args} {
    global NGTag
    global gap5_defs

    # Create the frame
    eval frame $path -class TagCheckList

    rename $path $path.widget
    proc $path {command args} {
	set path [lindex [info level [info level]] 0]
	if {[info procs tag_checklist_$command] != ""} {
	    return [eval tag_checklist_$command $path $args]
	} else {
	    return [eval $path.list [list $command] $args]
	}
    }

    bind $path <Destroy> {
	catch {rename $path {}}
    }

    # Create and populate our tablelist
    set l [eval [list tablelist $path.list \
			       -columns {0 "Index"
				   0 "Code"
				   0 "Tag name"} \
			       -labelcommand tablelist::sortByColumn \
			       -exportselection 0 \
			       -stretch 1 \
			       -yscrollcommand [list $path.yscroll set]
			  ] \
	       $args]
    $l columnconfigure 0 -sortmode integer
    for {set i 0} {$i < $NGTag(num_tags)} {incr i} {
	$l insert end [list $i $NGTag($i,tagid) $NGTag($i,tagtype)]
    }

    tag_checklist_array2list $l $array

    # Add a scrollbar    
    scrollbar $path.yscroll -command "$path.list yview"

    grid columnconfigure $path 0 -weight 1
    grid rowconfigure $path 0 -weight 1

    grid $path.list $path.yscroll -sticky nsew

    set mode [$l cget -selectmode]
    if {$mode == "extended" || $mode == "extended2" || $mode == "multiple"} {
	# Add select all and deselect all messages.
	frame $path.buttons
	button $path.buttons.select_all \
	    -text "Select all" \
	    -command "$path.list selection clear 0 end;
		     $path.list selection set 0 end"
	button $path.buttons.clear_all \
	    -text "Clear all" \
	    -command "$path.list selection clear 0 end"
	pack $path.buttons.select_all $path.buttons.clear_all \
	    -side left -expand 1
	grid $path.buttons -columnspan 2 -sticky nsew
    }

    return $path
}

proc tag_checklist_component {w comp} {
    switch $comp {
	frame {
	    return $w.widget
	}
	list {
	    return $w.list
	}
	default {
	    return -code error "No component '$comp'"
	}
    }
}

proc tag_checklist_selection {w} {
    set result ""
    foreach item [$w.list curselection] {
	lappend result [$w.list get $item]
    }

    return $result
}

proc tag_checklist_list2array {w} {
    set w $w.list

    set nitems [$w index end]
    set arr ""

    for {set i 0} {$i < $nitems} {incr i} {
	set v [$w get $i]
	set used([lindex $v 0]) [$w selection includes $i]
    }

    for {set i 0} {$i < $nitems} {incr i} {
	lappend arr $i $used($i)
    }

    return $arr
}

proc tag_checklist_array2list {w arr} {
    $w selection clear 0 end
    foreach {ind val} $arr {
	if {$val} {
	    $w selection set $ind
	}
    }
}


##############################################################################
#convert C tag ids and types into tcl ids and types
proc InitTagArray { } {
    global NGTag

    set taglist [get_tag_array]
    set NGTag(num_tags) [llength $taglist]
    for {set i 0}  {$i < $NGTag(num_tags)} {incr i} {
        set t [lindex $taglist $i]
	set NGTag($i,tagtype) [lindex $t 0]
	set NGTag($i,tagid)   [lindex $t 1]
	set NGTag($i,comment) [lindex $t 2]
    }
}

##############################################################################
# Convert a tag type into an ID
proc tag_type2id {type {default 0}} {
    global NGTag

    for {set i 0}  {$i < $NGTag(num_tags)} {incr i} {
	if {$NGTag($i,tagid) == $type} {
	    return $i
	}
    }
    return $default;
}

##############################################################################
# Removes a set of tags
proc DeleteTags {io} {
    global gap5_defs
    set t .delete_tags

    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Delete Tags"

    contig_id $t.id -io $io -range 0

    lorf_in $t.infile [keylget gap5_defs REMOVE_PAD_COLUMNS.INFILE] \
	"{contig_id_configure $t.id -state disabled} \
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state normal}" -bd 2 -relief groove

    frame $t.tags
    xentry $t.tags.l \
	-label "Tags to remove" \
	-default "(pick; * for all)" \
	-type str \
	-textvariable $t.Tags
    button $t.tags.b \
	-text "Select types" \
	-command "DeleteTags_pick $t $t.tags.b"
    pack $t.tags.l -side left -expand 1 -fill both
    pack $t.tags.b -side right

    #--- OK/cancel/help
    okcancelhelp $t.ok \
	-ok_command "DeleteTags2 $io $t" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 DeleteTags" \
	-bd 2 \
	-relief groove

    pack $t.infile $t.id $t.tags $t.ok -side top -fill both -expand 1
}

;proc DeleteTags2 {io t} {
    if {[lorf_in_get $t.infile] == 4} {
	set list [list [contig_id_gel $t.id]]
    } elseif {[lorf_in_get $t.infile] == 3} {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $t.infile]
    }

    global $t.Tags
    set tag_types [set $t.Tags]
    if {$tag_types == "*"} {
	set tag_types ""; # All
    }

    log_call delete_tags -io $io -contigs $list -tag_types $tag_types -verbose 1
    destroy $t
}

# Called when Select Types button is pressed
;proc DeleteTags_pick {top t} {
    # Create a modal window
    set w $t.select
    if {[winfo exists $w]} {
	return
    }
    modal $w -resizable 1 -width 200 -height 200
    wm transient $w [winfo toplevel $t]
    #wm overrideredirect $w 1

    # Move it relative to the button
    set rootx [winfo rootx $t]
    set rooty [winfo rooty $t]
    incr rootx 20 
    incr rooty 20 

    on_screen $w $rootx $rooty 300 300

    # Add the tag selector
    tag_checklist $w.type {} -selectmode extended
    pack $w.type -side top -fill both -expand 1
    #$w.type see 0

    okcancelhelp $w.ok \
	-ok_command "DeleteTags_picked $top $w.type" \
	-cancel_command "destroy $w"

    pack $w.ok -fill both
}

# Called when a tag selection has been made and we've hit OK.
;proc DeleteTags_picked {t w} {
    set list ""
    foreach x [$w selection] {
	lappend list [lindex $x 1]
    }

    global $t.Tags
    if {$list != ""} {
	set $t.Tags $list
    } else {
	set $t.Tags "*"
    }

    destroy [winfo toplevel $w]
}
