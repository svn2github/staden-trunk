#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#
# Brings up a search window
#
# search_setup array contains use_text use_strand use_tag use_where value1_name value2_name
set search_setup(position)		"0 0 0 0 {padded position}"
set search_setup(uposition)		"0 0 0 0 {unpadded position}"
set search_setup(problem)		"0 0 0 0"
set search_setup(anno)			"0 0 0 0 string"
set search_setup(sequence)		"0 1 0 1 sequence {num. mismatches}"
set search_setup(quality)		"0 0 0 0"
set search_setup(consquality)		"0 0 0 0 value"
set search_setup(file)			"1 0 0 0 filename"
set search_setup(name)			"0 0 0 0 name"
set search_setup(edit)			"0 0 0 0"
set search_setup(verifyand)		"0 0 0 0"
set search_setup(verifyor)		"0 0 0 0"
set search_setup(discrepancy)		"0 0 0 0 value"
set search_setup(consdiscrepancy)	"0 0 0 0 value"
set search_setup(tag)			"0 0 1 0"
set search_setup(difference)		"0 0 0 0"

proc create_search_win {w com {dir 0} {raise 1}} {
     global NGTag gap_defs search_setup

    # Remembered from previous uses of this search window.
    # FIXME: should we allow this?
    global $w.Type
    global $w.Direction
    global $w.Tag
    global $w.Where
    global $w.File $w.FileName
    global $w.Strand

    set df $w.direction_frame

    # FIXME: not an ideal assumption, but currently valid.
    set e [lindex $com 0]
	
    if {[winfo exists $w]} {
	if {[wm state $w] == "iconic"} {
	    wm deiconify $w
	    return
	} else {
	    if {$dir == 0} {
		set dir [set $w.Direction]
	    }
	    if {$dir == 1 || $dir == "forward"} {
		$df.forward invoke
	    } else {
		$df.backward invoke
	    }
	}
	if {$raise} {
	    raise $w
	}
	do_search $w $com
	return
    }

    if {[winfo exists $w]} {
	if {$raise} {
	    raise $w
	}
	return
    }
    catch {unset $w.File}; # Just incase
    set $w.FileName ""

    if {![info exists $w.Tag]} {
	set $w.Tag [keylget gap_defs CONTIG_EDITOR.SEARCH.TAG_DEF]
    }

    if {![info exists $w.Direction]} {
	set init_d 1
    } else {
	set init_d 0
    }

    if {![info exists $w.Strand]} {
	set init_s 1
    } else {
	set init_s 0
    }

    if {![info exists $w.Type]} {
	set init_t 1
    } else {
	set init_t 0
    }

    if {![info exists $w.Where]} {
	set init_w 1
    } else {
	set init_w 0
    }

    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Search"

    # Also see search_setup and do_file_search - keep in sync
    set tf $w.type_frame
    set df $w.direction_frame
    set wf $w.where_frame
    set of $w.option_frame
    set bf $w.button_frame

    # Direction frame
    frame $df -bd 2 -relief sunken
    label $df.label -text "Direction"
    radiobutton $df.backward -text "backward" -variable $w.Direction \
	-relief flat -value backward
    radiobutton $df.forward -text "forward" -variable $w.Direction \
	-relief flat -value forward
    pack $df.label -side top
    pack $df.backward -side left
    pack $df.forward -side right

    # Type frame
    frame $tf -bd 2 -relief sunken
    label $tf.label -text "Search by"
    frame $tf.l
    frame $tf.r
    radiobutton $tf.l.position  -text "padded position" -variable $w.Type \
	-relief flat -command "search_setup $w position" \
	-value position
    radiobutton $tf.l.problem   -text "problem"      -variable $w.Type \
	-relief flat -command "search_setup $w problem" \
	-value problem
    radiobutton $tf.l.anno      -text "annotation"   -variable $w.Type \
	-relief flat -command "search_setup $w anno" \
	-value anno
    radiobutton $tf.l.sequence  -text "sequence"     -variable $w.Type \
	-relief flat -command "search_setup $w sequence" \
	-value sequence
    radiobutton $tf.l.quality   -text "quality"      -variable $w.Type \
	-relief flat -command "search_setup $w quality" \
	-value quality
    radiobutton $tf.l.conqual  -text "consensus quality" -variable $w.Type \
	-relief flat -command "search_setup $w consquality" \
	-value consquality
    radiobutton $tf.l.file     -text "from file" -variable $w.Type \
	-relief flat -command "search_setup $w file" \
	-value file
    radiobutton $tf.r.uposition  -text "unpadded position" -variable $w.Type \
	-relief flat -command "search_setup $w uposition" \
	-value uposition
    radiobutton $tf.r.name      -text "reading name" -variable $w.Type \
	-relief flat -command "search_setup $w name" \
	-value name
    radiobutton $tf.r.edit      -text "edit"         -variable $w.Type \
	-relief flat -command "search_setup $w edit" \
	-value edit
    radiobutton $tf.r.verifyand -text "evidence for edits (1)" \
	-variable $w.Type \
	-relief flat -command "search_setup $w verifyand" \
	-value verifyand
    radiobutton $tf.r.verifyor  -text "evidence for edits (2)   " \
	-variable $w.Type \
	-relief flat -command "search_setup $w verifyor" \
	-value verifyor
    radiobutton $tf.r.tag   -text "tag type" -variable $w.Type \
	-relief flat -command "search_setup $w tag" \
	-value tag
    radiobutton $tf.r.discrepancies -text "discrepancies" -variable $w.Type \
	-relief flat -command "search_setup $w discrepancy" \
	-value discrepancy
    radiobutton $tf.r.consdiscrep -text "Consensus discrepancies" -variable $w.Type \
	-relief flat -command "search_setup $w consdiscrepancy" \
	-value consdiscrepancy

    pack $tf.label -side top
    pack $tf.l -side left
    pack $tf.r -side right
    pack $tf.l.position $tf.l.problem $tf.l.anno $tf.l.sequence $tf.l.quality \
	$tf.l.conqual $tf.l.file -anchor w
    pack $tf.r.uposition $tf.r.name $tf.r.edit $tf.r.verifyand $tf.r.verifyor \
	$tf.r.tag $tf.r.discrepancies $tf.r.consdiscrep -anchor w

    # For join editor only
    if {[$e join_mode]} {
	radiobutton $tf.l.difference -text "consensus diff." \
	    -variable $w.Type \
	    -relief flat -command "search_setup $w difference" \
	    -value difference
	pack $tf.l.difference -anchor w
    }


    # Options: values
    frame $of
    entrybox $of.val1 \
	 -title "Value:" -relief sunken -command "do_search $w \"$com\"" \
	 -exportselection 0
    entrybox $of.val2 \
	 -title "Value:" -relief sunken -command "do_search $w \"$com\"" \
	 -exportselection 0
    

    # Options: strand selector - only needed for sequence search
    frame $of.strand
    label $of.strand.label -text "strand"
    radiobutton $of.strand.top    -text top    -variable $w.Strand -value +
    radiobutton $of.strand.bottom -text bottom -variable $w.Strand -value -
    radiobutton $of.strand.both   -text both   -variable $w.Strand -value =
    pack $of.strand.label -side left
    pack $of.strand.both $of.strand.bottom $of.strand.top -side right


    # Options: text box - for 'from file' search
    frame $of.text -bd 2 -relief groove
    label $of.text.label -text "File comments"
    text $of.text.text -width 30 -height 5 -takefocus 0 -bg [. cget -bg]
    bindtags $of.text.text "$of.text.text . all"
    pack $of.text.label -side top -fill x
    pack $of.text.text -side top -fill both
    bind $of.text.text <<Clear>> {break;}
    bind $of.text.text <<Paste>> {break;}
    bind $of.text.text <2> {break;}
    bind $of.text.text <KeyPress> {break;}

    # Options: tag selector
    frame $of.tag -bd 2 -relief groove
    label $of.tag.label -text "Tag type"
    button $of.tag.button -text "" \
	 -command "tag_editor_select_type $tf.r.tag \[set $w.Tag\] \
		{show_help gap4 Editor-Searching} \
		{search_set_tag_type $of.tag.button $w.Tag}"
    search_set_tag_type $of.tag.button $w.Tag [tag_type2id [set $w.Tag]]
	
    pack $of.tag.label -side left
    pack $of.tag.button -side right -fill both -expand 1

    # Options: where
    set wf $of.where
    frame $wf -bd 0
    label $wf.label -text "Where"
    radiobutton $wf.seq  -text "Sequence"  -variable $w.Where \
	-relief flat -value 1
    radiobutton $wf.cons -text "Consensus" -variable $w.Where \
	-relief flat -value 2
    radiobutton $wf.both -text "Both"      -variable $w.Where \
	-relief flat -value 3
    pack $wf.label -side left
    pack $wf.both $wf.cons $wf.seq -side right

    # Buttons
    frame $bf -bd 2 -relief sunken
    button $bf.search -text "Search" -command "do_search $w \"$com\""
    button $bf.quit -text "Cancel" -command "destroy $w"
    button $bf.help -text "Help" \
	-command "show_help gap4 {Editor-Searching}"
    pack $bf.search $bf.quit -side left -expand 1
    pack $bf.help -side right -expand 1

    # key bindings
    bind $w <<search>> "do_search $w \"$com\""

    # Packing/Gridding
    grid columnconfigure $w 0 -weight 1
    grid $tf -sticky nsew
    grid $df -sticky nsew
    grid $of -sticky nsew
    grid $bf -sticky nsew
    grid columnconfigure $of 0 -weight 1

    # Init
    if {$init_d} {
	set $w.Direction [keylget gap_defs CONTIG_EDITOR.SEARCH.DEFAULT_DIRECTION]
    }
    if {$init_s} {
	set $w.Strand [keylget gap_defs CONTIG_EDITOR.SEARCH.DEFAULT_STRAND]
    }
    if {$init_t} {
	set $w.Type [keylget gap_defs CONTIG_EDITOR.SEARCH.DEFAULT_TYPE]
    }
    if {$init_w} {
	set $w.Where [keylget gap_defs CONTIG_EDITOR.SEARCH.DEFAULT_WHERE]
    }
    eval search_setup $w [set $w.Type]
}

proc search_set_tag_type {button var id} {
    global $var NGTag
    set type $NGTag($id,tagid)
    set $var $type
    $button configure -text "$type ($NGTag($id,tagtype))"
}

proc do_search {w com args} {
    global $w.Direction
    global $w.Type
    global $w.Tag
    global $w.Where
    global $w.Strand
    global $w.LastTagWindow
    global .cedit.Defaults

    set dir	[set $w.Direction]
    set type	[set $w.Type]
    set tag	[set $w.Tag]
    set where	[set $w.Where]
    set strand	[set $w.Strand]

    # Remember last search as the defaults for new editors.
    global gap_defs 
    keylset gap_defs CONTIG_EDITOR.SEARCH.DEFAULT_DIRECTION 	$dir
    keylset gap_defs CONTIG_EDITOR.SEARCH.DEFAULT_STRAND 	$strand
    keylset gap_defs CONTIG_EDITOR.SEARCH.DEFAULT_TYPE 		$type
    keylset gap_defs CONTIG_EDITOR.SEARCH.TAG_DEF 		$tag

    set value1 [entrybox_get $w.option_frame.val1]
    set value2 [entrybox_get $w.option_frame.val2]

    set .cedit.Defaults(${type}1) $value1
    set .cedit.Defaults(${type}2) $value2

    if {$type == "position" \
	    || $type == "anno" \
	    || $type == "name" \
	    || $type == "consquality" \
	    || $type == "discrepancy"
	    || $type == "consdiscrepancy"} {
        set r [eval $com $dir $strand $type [list $value1]]
    } elseif {$type == "uposition"} {
	set r [eval $com $dir $strand position [list u$value1]]
    } elseif {$type == "sequence"} {
	set r [eval $com $dir $strand $type [list $value1#$value2#$where]]
    } elseif {$type == "tag"} {
        set r [eval $com $dir $strand $type $tag]
    } elseif {$type == "file"} {
	do_file_search $w $com $dir $type $value1
	set r 1
    } elseif {$type == "difference"} {
	set e [lindex $com 0]; #editor name
	if {$dir == "forward"} {
	    $e next_difference
	} else {
	    $e prev_difference
	}
	set r 1
    } else {
	set r [eval $com $dir $strand $type]
    }

    if {[keylget gap_defs CONTIG_EDITOR.SEARCH.AUTO_EDIT_TAG] &&
	$type == "anno" && $r != 0} {
	if {[info exists $w.LastTagWindow] &&
	    [winfo exist [set $w.LastTagWindow]]} {
	    catch {destroy [set $w.LastTagWindow]}
	}
	set $w.LastTagWindow [[lindex $com 0] edit_anno]
	after idle "focus [lindex $com 0]"
    }

    if {$r == 0} {bell}
}

#
# Handles the "from file" search type
#
proc do_file_search {w com dir type fname} {
    # Keep these in sync with create_search_win
    set text_win $w.option_frame.text

    if {$fname == ""} {
	bell
	return
    }

    global $w.File $w.FileName $w.FileIndex

    # Handle file opening and closing, read all into $w.File array
    if {[set $w.FileName] != $fname} {
	catch {unset $w.File}
	set $w.File() ""; unset $w.File(); # blank array
	if {[catch {set fd [open $fname r]}]} {
	    bell
	    tk_messageBox -message "Couldn't open $fname" \
		-type ok \
		-parent $w
	    return
	}
	set $w.FileName $fname
	for {set i 0} {[gets $fd l] != -1} {} {
	    if {$l != ""} {
	        set $w.File($i) $l
		incr i
	    }
	}
	if {$dir == "forward"} {
	    set $w.FileIndex -1
	} else {
	    set $w.FileIndex $i
	}
        close $fd
    }

    # Determine direction
    if {$dir == "forward"} {
	set dir 1
    } else {
	set dir -1
    }

    # Loop through items until we find a valid one
    while {1} {
	set ind [expr [set $w.FileIndex]+$dir]

 	if {![info exists $w.File($ind)]} {
	    bell
	    return
	}

	set $w.FileIndex $ind
	set l [set $w.File($ind)]

	# Move editor
	if {[eval $com $dir = position "@[lindex $l 1]/[lindex $l 0]"] == 1} {
	    break;
	}
    }

    # We've found one, so update the text box
    $text_win.text delete 1.0 end
    $text_win.text insert end "Sequence: [lindex $l 0]\n"
    $text_win.text insert end "Position: [lindex $l 1]\n"
    $text_win.text insert end [subst [lrange $l 2 end]]
}


proc search_setup {w name} {
    global search_setup gap_defs

    foreach {use_text use_strand use_tag use_where val1_name val2_name} \
	$search_setup($name) {}

    # Keep these in sync with create_search_win
    set val1_win $w.option_frame.val1
    set val2_win $w.option_frame.val2
    set label1_win $w.option_frame.val1.label
    set label2_win $w.option_frame.val2.label
    set text_win $w.option_frame.text
    set strand_win $w.option_frame.strand
    set where_win $w.option_frame.where
    set tag_win $w.option_frame.tag
 
    global $w.Defaults .cedit.Defaults
    if {![info exists $w.Defaults(${name}1)]} {
	if {[info exists .cedit.Defaults(${name}1)]} {
	    set $w.Defaults(${name}1) [set .cedit.Defaults(${name}1)]
	    set $w.Defaults(${name}2) [set .cedit.Defaults(${name}2)]
	} elseif {$val1_name != ""} {
	    set u [keylget gap_defs \
		    CONTIG_EDITOR.SEARCH.[string toupper $name]_DEF]
	    set $w.Defaults(${name}1) [lindex $u 0]
	    set $w.Defaults(${name}2) [lindex $u 1]
	} else {
	    set $w.Defaults(${name}1) ""
	    set $w.Defaults(${name}2) ""
	}
    }

    # Special case for sequence search 
    if {$val1_name != ""} {
	entrybox_configure $val1_win -textvariable $w.Defaults(${name}1)
	grid $val1_win -row 0 -sticky nsew
    } else {
	grid forget $val1_win
    }

    if {$val2_name != ""} {
	entrybox_configure $val2_win -textvariable $w.Defaults(${name}2)
	grid $val2_win -row 1 -sticky nsew
    } else {
	grid forget $val2_win
    }
 
    if {$use_text} {
	grid $text_win -row 3 -sticky nsew
    } else {
	grid forget $text_win
	$text_win.text delete 1.0 end
	global $w.FileName $w.File
	set $w.FileName ""
	catch {unset $w.File}
    }

    if {$use_strand} {
	grid $strand_win -row 2 -sticky nsew
    } else {
	grid forget $strand_win
    }

    if {$use_tag} {
	grid $tag_win -row 4 -sticky nsew
    } else {
	grid forget $tag_win
    }

    if {$use_where} {
	grid $where_win -row 5 -sticky nsew
    } else {
	grid forget $where_win
    }

    if {$val1_name != ""} {
        $label1_win configure -text $val1_name
    }

    if {$val2_name != ""} {
        $label2_win configure -text $val2_name
    }

    [entrybox_path $val1_win] selection range 0 end
    [entrybox_path $val2_win] selection range 0 end
    focus [entrybox_path $val1_win]

    $w.option_frame configure -height 1
}
