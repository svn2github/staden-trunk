#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 2001. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

# ----------------------------------------------------------------------------
# Text display bindings
proc gap4_text_init {t} {
    set bg [tk::Darken [$t cget -background] 95]
    $t tag configure HLIGHT -relief raised -borderwidth 2 -foreground blue

    # Sequence SEQIDs
    $t tag configure SEQID -underline 1; #-background $bg
    $t tag bind SEQID <<use>> "gap4_text_SEQID_selected \$io %W SEQID"
    $t tag bind SEQID <<menu>> "gap4_text_SEQID_popup \$io %W SEQID %X %Y; break"
    $t tag bind SEQID <Any-Enter> "gap4_text_enter \$io %W SEQID"
    $t tag bind SEQID <Any-Motion> "gap4_text_motion \$io %W SEQID"
    $t tag bind SEQID <Any-Leave> "gap4_text_leave \$io %W SEQID"

    # TEMPLATE names
    # What can we do here?

    # TAGIDs
    $t tag configure TAGID -underline 1; # -background $bg
    $t tag bind TAGID <<use>> "gap4_text_TAGID_selected \$io %W TAGID"
    $t tag bind TAGID <<menu>> "gap4_text_TAGID_popup \$io %W TAGID %X %Y; break"
    $t tag bind TAGID <Any-Enter> "gap4_text_enter \$io %W TAGID"
    $t tag bind TAGID <Any-Motion> "gap4_text_motion \$io %W TAGID"
    $t tag bind TAGID <Any-Leave> "gap4_text_leave \$io %W TAGID"
}

;proc gap4_text_enter {io w tag} {
    set range [$w tag prevrange $tag current]
    if {[llength $range] != 2} return 
    foreach {s e} [$w tag ranges HLIGHT] {
	$w tag remove HLIGHT $s $e
    }
    $w tag add HLIGHT [lindex $range 0] [lindex $range 1]
}

;proc gap4_text_motion {io w tag} {
    set range [$w tag prevrange $tag current]
    if {[llength $range] != 2} return 
    foreach {s e} [$w tag ranges HLIGHT] {
	$w tag remove HLIGHT $s $e
    }
    $w tag add HLIGHT [lindex $range 0] [lindex $range 1]
}

;proc gap4_text_leave {io w tag} {
    foreach {s e} [$w tag ranges HLIGHT] {
	$w tag remove HLIGHT $s $e
    }
}

;proc gap4_text_SEQID_selected {io w tag} {
    set range [$w tag prevrange $tag current]
    set text [$w get [lindex $range 0] [lindex $range 1]]

    if {[db_info get_read_num $io $text] <= 0} {
	verror ERR_WARN "Reading '$text' not found in database"
	bell
	return
    }

    edit_contig -io $io -contig $text -reading $text -reuse 1
}

;proc gap4_text_SEQID_popup {io w tag X Y} {
    set range [$w tag prevrange $tag current]
    set text [$w get [lindex $range 0] [lindex $range 1]]

    if {[db_info get_read_num $io $text] <= 0} {
	verror ERR_WARN "Reading '$text' not found in database"
	bell
	return
    }

    if {[winfo exists $w.m]} {destroy $w.m}
    create_popup $w.m "Reading/Contig Commands ($text)"

    $w.m add command -label "Edit contig" \
	-command "edit_contig -io $io -contig $text -reading $text -reuse 1"

    $w.m add command -label "Template display" \
	-command "CreateTemplateDisplay $io $text"

    $w.m add command -label "List reading notes" \
	-command "NoteSelector $io reading $text"

    $w.m add command -label "List contig notes" \
	-command "NoteSelector $io contig $text"

    tk_popup $w.m [expr $X-20] [expr $Y-10]
}

;proc gap4_text_TAGID_parse {io w tag} {
    set range [$w tag prevrange $tag current]
    set text [$w get [lindex $range 0] [lindex $range 1]]

    # Parse the tagged text, eg "tagnum (reading readid, pos position (cpos))"
    set l [regexp -inline \
	       {.*([0-9])+\s+\(reading\s+(.*), pos ([0-9]+) \(.*\)\).*} \
	       $text]
    if {$l != ""} {
	return [lreplace $l 0 0 1]
    } else {
	set l [regexp -inline \
		   {.*([0-9])+\s+\(contig\s+(.*), pos ([0-9]+)\).*} $text]
	return [lreplace $l 0 0 0]
    }
}

;proc gap4_text_TAGID_selected {io w tag} {
    foreach {is_read tnum id pos} [gap4_text_TAGID_parse $io $w $tag] {}

    if {[db_info get_read_num $io $id] <= 0} {
	verror ERR_WARN "Reading '$id' not found in database"
	bell
	return
    }

    # Invoke the editor
    if {$is_read} {
	edit_contig -io $io -contig $id -reading $id -pos $pos -reuse 1
    } else {
	edit_contig -io $io -contig $id -pos $pos -reuse 1
    }
}


;proc gap4_text_TAGID_popup {io w tag X Y} {
    foreach {is_read tnum id pos} [gap4_text_TAGID_parse $io $w $tag] {}

    if {[db_info get_read_num $io $id] <= 0} {
	verror ERR_WARN "Reading '$id' not found in database"
	bell
	return
    }

    if {[winfo exists $w.m]} {destroy $w.m}
    create_popup $w.m "Reading/Contig Commands (tag #$tnum)"

    if {$is_read} {
	$w.m add command -label "Edit contig" \
	    -command "edit_contig -io $io -contig $id -pos $pos -reuse 1"
    } else {
	$w.m add command -label "Edit contig" \
	    -command "edit_contig -io $io -contig $id -reading $id -pos $pos -reuse 1"
    }

    $w.m add command -label "Template display" \
	-command "CreateTemplateDisplay $io $id"

    $w.m add command -label "List reading notes" \
	-command "NoteSelector $io reading $id"

    $w.m add command -label "List contig notes" \
	-command "NoteSelector $io contig $id"

    tk_popup $w.m [expr $X-20] [expr $Y-10]
}

# ----------------------------------------------------------------------------
# Common search dialog

;proc SearchDialog {io type tag title help_node} {
    global gap_defs
    
    set w .search_seq_$type
    if {[xtoplevel $w -resizable 0] == ""} return
    if {$title == ""} {
	set title "Search $type names"
    }
    wm title $w $title

    entrybox $w.pattern \
	-title "Sequence $type pattern"

    yes_no $w.case \
	-title "Case insensitive" \
	-bd 0 \
	-orient horizontal \
	-default [keylget gap_defs SEARCH_CASE]

    radiolist $w.ptype \
	-title "Pattern type" \
	-bd 0 \
	-relief groove \
	-orient horizontal \
	-default [keylget gap_defs SEARCH_MODE] \
	-buttons {
	    {{regular expression}}
	    {{wild-cards}}
	    {{sub-string}}
	}

    entrybox $w.list \
	-title "Save to list named" \
	-default [keylget gap_defs SEARCH_LIST_[string toupper $type]]

    okcancelhelp $w.ok \
	-ok_command "SearchDialog2 $io $type $w 1 $tag" \
	-apply_command "SearchDialog2 $io $type $w 0 $tag" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap4 [list $help_node]" \
	-bd 2 -relief groove

    pack $w.pattern $w.case $w.ptype $w.list $w.ok -side top -fill both
}

;proc SearchDialog2 {io type w destroy {tag {}}} {
    if {[set pattern [entrybox_get $w.pattern]] == ""} {
	return
    }
    set lpattern [string tolower $pattern]
    set ptype [radiolist_get $w.ptype]
    set case [yes_no_get $w.case]
    set lname [entrybox_get $w.list]

    if {$destroy} {
	destroy $w
    }

    set list [search_${type}_list $io]

    set db [io_read_database $io]
    set nr [keylget db num_readings]

    set l ""
    if {$ptype == 1 && $case} {
	foreach {id name} $list {
	    if {[regexp -nocase ".*$pattern.*" $name]} {lappend l $id $name}
	}
    } elseif {$ptype == 1} {
	foreach {id name} $list {
	    if {[regexp ".*$pattern.*" $name]} {lappend l $id $name}
	}
    } elseif {$ptype == 2 && $case} {
	foreach {id name} $list {
	    if {[string match -nocase "$pattern" $name]} {
		lappend l $id $name
	    }
	}
    } elseif {$ptype == 2} {
	foreach {id name} $list {
	    if {[string match "$pattern" $name]} {lappend l $id $name}
	}
    } elseif {$ptype == 3 && $case} {
	foreach {id name} $list {
	    if {[string first $lpattern [string tolower $name]] >= 0} {
		lappend l $id $name
	    }
	}
    } elseif {$ptype == 3} {
	foreach {id name} $list {
	    if {[string first $lpattern $name] >= 0} {
		lappend l $id $name
	    }
	}
    } else {
	verror ERR_WARN "Unknown search pattern type: $ptype"
	return
    }

    if {$lname != ""} {
	set l2 ""
	if {[info commands search_${type}_format] != ""} {
	    foreach {id name} $l {
		lappend l2 [search_${type}_format $io $id $name]
	    }
	} else {
	    foreach {id name} $l {
		lappend l2 $name
	    }
	}
	ListCreate2 $lname $l2 $tag
	ListEdit $lname
    }

    foreach {id name} $l {
	search_${type}_disp $io $id $name
    }
}

# ----------------------------------------------------------------------------
# Search by sequence name

proc SearchSeqDialog {io} {
    vfuncheader "Search sequence names"
    SearchDialog $io sequence SEQID "Search Sequence Names" List-SearchSequenceNames
}

;proc search_sequence_list {io} {
    set db [io_read_database $io]
    set nr [keylget db num_readings]

    set l ""

    for {set i 1} {$i <= $nr} {incr i} {
	lappend l $i [r_name $io $i]
    }

    return $l
}

;proc search_sequence_disp {io id name} {
    vmessage_tagged $name SEQID
}


# ----------------------------------------------------------------------------
# Search by template name

proc SearchTemplateDialog {io} {
    global _tmp_templates

    vfuncheader "Search template names"
    SearchDialog $io template {}  "Search Template Names" List-SearchTemplateNames

    catch {unset _tmp_templates}
}

;proc search_template_list {io} {
    global _tmp_templates

    set db [io_read_database $io]
    set nt [keylget db Ntemplates]
    set nr [keylget db num_readings]

    set l ""

    # Search for template names
    for {set i 1} {$i <= $nt} {incr i} {
	set t [io_read_template $io $i]
	set tname($i) [io_read_text $io [keylget t name]]
	lappend l $i $tname($i)
	set _tmp_templates($tname($i)) ""
    }

    # Precompute lookup table of template to sequence names
    for {set i 1} {$i <= $nr} {incr i} {
	set r [io_read_reading $io $i]
	set name [io_read_text $io [keylget r name]]
	set temp [keylget r template]
	lappend _tmp_templates($tname($temp)) $name
    }

    return $l
}

;proc search_template_disp {io id name} {
    global _tmp_templates
    vmessage -nonewline "$name "
    set seqs $_tmp_templates($name)
    if {$seqs != ""} {
	vmessage -nonewline "(Sequence names: "
	foreach seq [lrange $seqs 0 end-1] {
	    vmessage_tagged -nonewline "$seq" SEQID
	    vmessage -nonewline " "
	}
	vmessage_tagged  -nonewline [lindex $seqs end] SEQID
	vmessage ")"
    } else {
	vmessage "(no sequences using this template)"
    }
}

# ----------------------------------------------------------------------------
# Search by annotation contents

proc SearchAnnoDialog {io} {
    vfuncheader "Search annotation contents"
    SearchDialog $io annotation TAGID "Search Annotation Contents" List-SearchAnnotations
}

;proc search_annotation_list {io} {
    set db [io_read_database $io]
    set nr [keylget db num_readings]
    set nc [keylget db num_contigs]

    set l ""

    # Search readings
    for {set i 1} {$i <= $nr} {incr i} {
	set r [io_read_reading $io $i]
	for {set x [keylget r annotations]} {$x!=0} {set x [keylget t next]} {
	    set t [io_read_annotation $io $x]
	    set a [keylget t annotation]
	    if {$a != 0} {
		lappend l [list $x $i] [io_read_text $io $a]
	    }
	}
    }

    # Search contigs
    for {set i 1} {$i <= $nc} {incr i} {
	set c [io_read_contig $io $i]
	for {set x [keylget c annotations]} {$x!=0} {set x [keylget t next]} {
	    set t [io_read_annotation $io $x]
	    set a [keylget t annotation]
	    if {$a != 0} {
		lappend l [list $x -$i] [io_read_text $io $a]
	    }
	}
    }

    return $l
}

;proc search_annotation_format {io id name} {
    foreach {id rcnum} $id {}
    set t [io_read_annotation $io $id]
    set type [keylget t type]

    if {[keylget t annotation] != 0} {
	set anno [io_read_text $io [keylget t annotation]]
	regsub -all "\n" $anno " " anno
    } else {
	set anno ""
    }
    regsub "\n$" $name {} name

    if {$rcnum > 0} {
	set r [io_read_reading $io $rcnum]
	set rpos [keylget t position]
	set rstart [keylget r start]
	if {[keylget r sense] == 0} {
	    set rpos [expr {$rpos-$rstart}]
	} else {
	    set rlen [keylget r length]
	    set len  [keylget t length]
	    set rpos [expr {$rlen-($rpos+$len-1+$rstart)+1}]
	}
	set cpos [expr {$rpos + [keylget r position]-1}]
	return "Tag $id (reading [r_name $io $rcnum], pos $rpos ($cpos)) $anno"
    } else {
	set rcnum [db_info chain_left $io =[expr {-$rcnum}]]
	set pos [keylget t position]
	return "Tag $id (contig [r_name $io $rcnum], pos $pos) $anno"
    }
}

;proc search_annotation_disp {io id name} {
    foreach {id rcnum} $id {}
    set t [io_read_annotation $io $id]
    set type [keylget t type]

    regsub "\n$" $name {} name

    if {$rcnum > 0} {
	set r [io_read_reading $io $rcnum]
	set rpos [keylget t position]
	set rstart [keylget r start]
	if {[keylget r sense] == 0} {
	    set rpos [expr {$rpos-$rstart}]
	} else {
	    set rlen [keylget r length]
	    set len  [keylget t length]
	    set rpos [expr {$rlen-($rpos+$len-1+$rstart)+1}]
	}
	set cpos [expr {$rpos + [keylget r position]-1}]
	vmessage_tagged \
	    "Tag id " "" \
	    "$id (reading [r_name $io $rcnum], pos $rpos ($cpos))" TAGID \
	    ", type $type, contents=\"$name\"\n" ""
    } else {
	set rcnum [db_info chain_left $io =[expr {-$rcnum}]]
	set pos [keylget t position]
	vmessage_tagged \
	    "Tag id " "" \
	    "$id (contig [r_name $io $rcnum], pos $pos)" TAGID \
	    ", type $type, contents=\"$name\"\n" ""
    }
}
