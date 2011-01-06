#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#
# Displays a tag editor.
#
# We fill out a tcl array for the tag which is the converted to a tag structure
# in C.
#
# Our data ($d) array contains elements "type", "strand", and "anno".
# Additionally we have the name of a command ($c) passed to us to call for
# certain actions (save, quit, etc).
#
proc create_tag_editor {w c data} {
    upvar #0 $data d
    global NGTag gap_defs default_tag_type licence read_only

    if {[xtoplevel $w] == ""} { return }
    # wm minsize $w 280 100
    wm protocol $w WM_DELETE_WINDOW "$c quit"
    bind $w <Destroy> "if {\"%W\" == \"$w\"} {tag_editor_destroy %W}"

    if {![info exists default_tag_type]} {
	set default_tag_type COMM
    }

    # Buttons at the top
    frame $w.bar

    # A GUI frame for gui based tags
    frame $w.gui

    # Scrollable text panel
    set v $w.disp
    frame $v
    frame $v.b
    frame $v.b.fudge -bd 2 -relief sunken
 
    button $w.bar.cancel -text " Cancel " -command "$c quit"

    if {[info exists d(macro)]} {
	button $w.bar.save   -text " Save Macro" \
	 -command "tag_editor_save [list $c] $data $w"
	button $w.bar.delete   -text " Delete Macro" \
	 -command "tag_macro_delete $w $data"
    } else {
	button $w.bar.save   -text " Save " \
	 -command "tag_editor_save [list $c] $data $w"
	button $w.bar.move -text " Move " \
	    -command "tag_editor_moveorcopy move \
                      [list $c] $data \[$v.t get 1.0 end-1chars\]"
	button $w.bar.copy -text " Copy " \
	    -command "tag_editor_moveorcopy copy \
                      [list $c] $data \[$v.t get 1.0 end-1chars\]"
    }

    if {$d(type) == "NONE"} {
	set type $default_tag_type
    } else {
	set type $d(type)
    }

    set tid [tag_type2id $type -1]
    if {$tid != -1} {
	set long_type $NGTag($tid,tagtype)
    } else {
	# Unknown tag type - Make one up
	set tid $NGTag(num_tags)
	set long_type $d(type)
	set NGTag($tid,tagtype) $type
	set NGTag($tid,tagid) $type
	set NGTag($tid,comment) ""
	incr NGTag(num_tags)
    }

    button $w.bar.type -text "Type:$long_type" \
	-command "tag_editor_select_type $w.bar.type $type \
		  {show_help gap4 Editor-Annotations} \
		  {tag_editor_set_type $w.bar.type $w $v.t $data}"

    button $w.bar.strand -text [tag_editor_strand $d(strand)] \
	-command "tag_editor_set_strand $w.bar.strand $data"
    button $w.bar.help -text "Help" \
	-command "show_help gap4 {Editor-Annotations}"
 
    scrollbar $v.sy -orient vert -command "$v.t yview" -bd 2 -relief sunken
    scrollbar $v.b.sx -orient horiz -command "$v.t xview" -bd 2 -relief sunken
    text $v.t -height 8 -width 20 -yscrollcommand "$v.sy set" \
	-xscrollcommand "$v.b.sx set" -bd 2 -relief sunken

    if {$d(type) == "NONE"} {
	tag_editor_set_type $w.bar.type $w $v.t $data \
	    [tag_type2id $default_tag_type]
    } else {
        $v.t insert 1.0 $d(anno)
    }
    bind $v.t <Any-KeyPress> "upvar #0 $data d; set d(default) 0"
 
    if {$read_only && ![info exists d(macro)]} {
	$w.bar.save configure -state disabled
	$w.bar.move configure -state disabled
	$w.bar.copy configure -state disabled
	$w.bar.type configure -state disabled
	$w.bar.strand configure -state disabled
	$v.t configure -state disabled
    }

    # Button bar
    pack $w.bar -side top -fill both
    pack $w.bar.cancel -side left
    if {$licence(type) != "v" || [info exists d(macro)]} {
	pack $w.bar.save -side left
    }
    if {$licence(type) != "v" && ![info exists d(macro)]} {
	pack $w.bar.move -side left
	pack $w.bar.copy -side left
    }
    if {[info exists d(macro)]} {
	pack $w.bar.delete -side left
    }
    pack $w.bar.type -side left -fill both
    pack $w.bar.strand -side left
    pack $w.bar.help -side right
 
    # A scrolled view of data
    pack $v.b -side bottom -fill both
    pack $v.t -fill both -side left -expand 1
    pack $v.sy -side right -fill y 
    pack $v.b.sx -side left -fill both -expand 1

    # Hideous trickery to get the scrollbars to position correctly.
    pack $v.b.fudge -side right
    pack propagate $v.b.fudge 0
    $v.b.fudge configure -width 18 -height 18
    update idletasks
    pack propagate $v.b.fudge 0
    $v.b.fudge configure -width [winfo width $v.sy] -height [winfo width $v.sy]

    # Pack the appropriate dialogue component; text frame or gui
    tag_editor_set_type $w.bar.type $w $v.t $data $tid
}

proc on_screen {win x y w h} {
    set win [winfo toplevel $win]
    set sw [winfo screenwidth $win]
    set sh [winfo screenheight $win]

    # Compensate for window borders
    set wp [winfo toplevel [winfo parent $win]]
    set y_geom [lindex [split [wm geometry $wp] x+] 3]
    set y_real [winfo y $wp]
    set y_diff [expr {$y_real-$y_geom}]
    incr sh -$y_diff

    if {[expr {($x+$w) > $sw}]} {
	set x [expr {$sw-$w}]
    }

    if {[expr {($y+$h) > $sh}]} {
	set y [expr {$sh-$h}]
    }

    if {$x < 0} {set x 0}
    if {$y < 0} {set y 0}

    wm geometry $win ${w}x$h+$x+$y
}

proc tag_editor_destroy {w} {
    # Destroy the GUI details array so that variable traces are removed and
    # so that we have no old values being used the next time we open this
    # dialogue.
    global ::$w.GUI
    catch {unset ::$w.GUI}
}

proc tag_editor_select_type {button type help command} {
    global default_tag_type NGTag

    # Create a modal window
    set w $button.select
    if {[winfo exists $w]} {
	return
    }
    modal $w -resizable 1
    wm transient $w [winfo toplevel $button]
    # wm overrideredirect $w 1

    # Move it relative to the button
    set rootx [winfo rootx $button]
    set rooty [winfo rooty $button]
    incr rootx 20 
    incr rooty 20 

    on_screen $w $rootx $rooty 300 300

    # Initialise all tags to off by default
    set list ""
    for {set i 0} {$i < $NGTag(num_tags)} {incr i} {
	lappend list $NGTag($i,tagid)
    }
    if {[set ind [lsearch -exact $list $type]] != -1} {
	set arr "$ind 1"
    } else {
	set ind 0
	set arr ""
    }

    # Add the tag selector
    tag_checklist $w.type $arr -selectmode browse -width 35 -height 10
    pack $w.type -side top -fill both -expand 1
    $w.type see $ind

    okcancelhelp $w.ok \
	-ok_command "set type \[lindex \[lindex \[$w.type selection\] 0\] 1\];
		     $command \[tag_type2id \$type\];
                     after idle {destroy $w}" \
	-cancel_command "destroy $w" \
	-help_command $help

    set body [[$w.type component list] bodypath]
    bind $body <Double-1> "$w.ok.ok invoke"

    pack $w.ok -fill both
}

#
# Updates the tag type in 'data' and changes the menu name too
#
proc tag_editor_set_type {button w textwin data typeid} {
    upvar #0 $data d
    global default_tag_type NGTag

    set type $NGTag($typeid,tagtype)
    $button configure -text "Type:$type"

    if {[string match {#!acdtag:*} $NGTag($typeid,comment)]} {
	pack forget $w.disp
	catch {destroy $w.gui.f}
	pack [frame $w.gui.f] -fill both -expand 1
	set tid $NGTag($typeid,tagid)
	if {[info commands ::tag_gui::${tid}::create_dialogue] == {}} {
	    set code [acd2tag::parse $NGTag($typeid,comment) ::tag_gui::$tid]
	    #catch {set fd [open /tmp/jkb.out w]; puts $fd $code; close $fd}
	    eval $code
	}
	set d(namespace) ::$w.GUI
	catch {array set $d(namespace) [::acd2tag::str2array $d(anno)]} err
	#set d(namespace) ::tag_gui::${type}::$w
	set $w.GUI $d(namespace)
	::tag_gui::${tid}::create_dialogue $w.gui.f $d(namespace)
	pack $w.gui -side bottom -fill both -expand 1
    } else {
	catch {destroy $w.gui.f}
	pack forget $w.gui
	pack $w.disp -side bottom -fill both -expand 1
    }

    if {$d(default) == 1} {
	if {[string match {#!acdtag:*} $NGTag($typeid,comment)]} {
	    set d(anno) ""
	} else {
	    set d(anno) $NGTag($typeid,comment)
	    $textwin delete 1.0 end
	    $textwin insert 1.0 $d(anno)
	}
    }
    set default_tag_type $NGTag($typeid,tagid)
    set d(type) $default_tag_type
}

#
# Saves a tag (which exits the editor).
#
proc tag_editor_save {com data w} {
    global NGTag
    upvar #0 $data d

    set typeid [tag_type2id $d(type)]
    if {[string match {#!acdtag:*} $NGTag($typeid,comment)]} {
	set anno [::acd2tag::get_values ::tag_gui::$d(type) $d(namespace)]
    } else {
	set anno [$w.disp.t get 1.0 end-1chars]
    }

    set d(anno) $anno
    eval $com save
}

proc tag_editor_moveorcopy {method com data anno} {
    upvar #0 $data d

    set owner [selection own]
    if {$owner == ""} {
	bell
	return
    }
    if {[winfo class $owner] != "Editor"} {
	bell
	return
    }

    set d(anno) $anno
    eval $com $method $owner
}

proc tag_editor_strand {strand} {
    return [lindex {"----->" "<-----" "<---->"} $strand]
}

#
# Set the strand for a tag and update the strand button text.
#
proc tag_editor_set_strand {w data} {
    upvar #0 $data d

    if {$d(strand) == 0} {
	set d(strand) 1
    } elseif {$d(strand) == 1} {
	set d(strand) 2
    } elseif {$d(strand) == 2} {
	set d(strand) 0
    }

    $w configure -text [tag_editor_strand $d(strand)]
}

#
# Brings up the tag creation dialogue, which is used to define a tag macro
#
proc tag_macro_create {ed key} {
    set d ed_macro_$key
    global $d
    set w $ed.macro_$key
    if {![info exists $d]} {
	set ${d}(type) COMM
	set ${d}(anno) ""
	set ${d}(strand) 2
	set ${d}(default) 0
	set ${d}(set) 0
    }
    set ${d}(macro) 1
    set ${d}(backup) [array get $d]

    if {![winfo exists $w]} {
	create_tag_editor $w [list tag_macro_callback $w $d] $d
    } else {
	wm deiconify $w
	raise $w
    }
}

# Handles the "quit" and "save" callbacks from the macro dialogue
proc tag_macro_callback {win data command} {
    upvar #0 $data d
    
    if {$command == "quit"} {
	array set d $d(backup)
    } else {
	set d(set) 1
    }

    destroy $win
}

# Invokes a macro - ie adds a tag to the current cursor position
proc tag_macro_invoke {ed key} {
    set d ed_macro_$key
    global $d gap_defs

    if {![info exists $d] || [set ${d}(set)] == 0} {
	bell
	return
    }

    $ed create_anno [set ${d}(type)] [set ${d}(anno)] [set ${d}(strand)]
    if {[keylget gap_defs CONTIG_EDITOR.MACRO_AUTOEDIT]} {
	$ed edit_anno
    }
}

# Called by control-F* keys. Copies the tag under the cursor to the associated
# macro.
proc tag_macro_copy {ed key} {
    set d ed_macro_$key
    global $d
    set w $ed.macro_$key

    # only top-most tag.
    set tag [lindex [$ed list_anno] 0]
    foreach {ptr type st len} $tag {}
    set ${d}(type)    $type
    set ${d}(anno)    [$ed get_anno $ptr]
    set ${d}(strand)  0
    set ${d}(default) 0
    set ${d}(set)     1
    set ${d}(macro)   1
}

# Loads macro defintions from the gap_defs (and hence from .gaprc)
proc tag_macro_load {} {
    global ed_macro_keys gap_defs

    set c CONTIG_EDITOR
    foreach key $ed_macro_keys {
	if {[catch {set macro [keylget gap_defs $c.TAG_MACRO_$key]}]} {
	    continue
	}
	if {$macro == ""} continue
	global ed_macro_$key
	# Expand \n and \\ back - study this carefully if it needs fixing.
	set macro [string map "\\\\\\\\ \\\\ \\\\n \\n" $macro]
	array set ed_macro_$key $macro
	array set ed_macro_$key {set 1}
    }
}

# Saves macro definitions to the .gaprc file
proc tag_macro_save {} {
    global ed_macro_keys gap_defs env

    set c CONTIG_EDITOR
    foreach key $ed_macro_keys {
	global ed_macro_$key
	if {[info exists ed_macro_$key]} {
	    array set copy [array get ed_macro_$key]
	    catch {unset copy(backup)}
	    catch {unset copy(macro)}
	    catch {unset copy(set)}
	    # Escape newlines and backslashes
	    set str [array get copy]
	    set str [string map "\\n \\\\n \\\\ \\\\\\\\" $str]
	    keylset gap_defs $c.TAG_MACRO_$key $str
	} else {
	    keylset gap_defs $c.TAG_MACRO_$key {}
	}
	update_defs gap_defs $env(HOME)/.gaprc $c.TAG_MACRO_$key
    }
}

proc tag_macro_delete {win data} {
    global $data
    unset $data
    destroy $win
}
