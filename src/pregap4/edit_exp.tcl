proc edit_exp_types {} {
    global line_types line_comval line_com line_val

    set w .edit_exp_types

    if {[modal $w -resizable 1] == ""} {
	return
    }
    wm title $w "Experiment file line types"
    wm geometry $w 500x500

    # Create vertically scrolled frame
    frame $w.sf -bd 2 -relief groove
    pack $w.sf -side top -fill both -expand 1

    set c $w.sf.canvas
    scrollbar $w.sf.yscroll -command "$c yview"
    pack $w.sf.yscroll -side right -fill y

    canvas $c -yscrollcommand "$w.sf.yscroll set"
    pack $c -side right -fill both -expand 1

    set f [frame $c.sframe]
    $c create window 0 0 -window $c.sframe -anchor nw
    bind $c <Any-Configure> "
	%W itemconfigure 1 \
	    -width  \[winfo width %W\]"

    # Create line type windows
    grid columnconfigure $f 0 -weight 0
    grid columnconfigure $f 1 -weight 0
    grid columnconfigure $f 2 -weight 1
    foreach lt [lsort $line_types] {
	if {[catch {set com [info body ${lt}_com]}] == 0} {
	    set line_com($lt) $com
	} else {
	    set line_com($lt) ""
	}
	global $lt
	if {[info exists $lt]} {
	    set line_val($lt) [set $lt]
	} else {
	    set line_val($lt) ""
	}

	label $f._${lt}_label \
	    -text "Type: $lt   "
	global $f._${lt}_Com_Val
	tk_optionMenu $f._${lt}_com_val line_comval($lt) Value Command
	trace variable line_comval($lt) w "edit_exp_type_trace $f"
	entry $f._${lt}_entry -width 20
	set line_comval($lt) [lindex {Value Command} [is_command $lt]]

	grid $f._${lt}_label $f._${lt}_com_val $f._${lt}_entry \
	    -sticky nsew -padx 5

	bind $f._${lt}_label   <Any-Enter> "+edit_exp_type_description $w $lt"
	bind $f._${lt}_com_val <Any-Enter> "+edit_exp_type_description $w $lt"
	bind $f._${lt}_entry   <Any-Enter> "+edit_exp_type_description $w $lt"
    }

    # Add other windows
    label $w.description
    pack $w.description -fill both -side top -pady 5 -padx 5

    frame $w.buttons
    pack $w.buttons -fill both -side top

    button $w.buttons.okm \
	-text "OK (in memory)" \
	-command "edit_exp_type_ok $w 0"
    button $w.buttons.oks \
	-text "OK (and save)" \
	-command "edit_exp_type_ok $w 1"
    button $w.buttons.cancel \
	-text Cancel \
	-command "edit_exp_type_cancel $w"
    button $w.buttons.help \
	-text Help \
	-command "show_help pregap4 {Pregap4-Database-LineTypes}"
    pack $w.buttons.okm $w.buttons.oks $w.buttons.cancel -side left
    pack $w.buttons.help -side right

    update idletasks
    $c configure -scrollregion "0 0 [winfo width $f] [winfo height $f]"
}

proc edit_exp_type_trace {f aname lt op} {
    global line_comval line_com line_val

    $f._${lt}_entry configure -textvariable {}
    $f._${lt}_entry delete 0 end
    if {$line_comval($lt) == "Value"} {
	$f._${lt}_entry insert end $line_val($lt)
	$f._${lt}_entry configure -textvariable line_val($lt)
    } else {
	$f._${lt}_entry insert end $line_com($lt)
	$f._${lt}_entry configure -textvariable line_com($lt)
    }
}

proc edit_exp_type_ok {w save} {
    global line_comval line_com line_val line_types

    set global_conf ""

    foreach lt [lsort $line_types] {
	if {$line_comval($lt) == "Value"} {
	    catch {rename ${lt}_com {}}
	    global $lt
	    if {$line_val($lt) != ""} {
		append global_conf "set $lt [list $line_val($lt)]\n"
	    } else {
		catch {unset $lt}
	    }
	} else {
	    if {$line_com($lt) != ""} {
		append global_conf "proc ${lt}_com {} [list $line_com($lt)]\n"
	    } else {
		catch {rename ${lt}_com {}}
	    }
	}
    }

    if {$global_conf != ""} {
        uplevel #0 eval [list $global_conf]
    }

    if {$save} {
        set_conf_data_section global_variables $global_conf
        write_conf_file
    }

    edit_exp_type_cancel $w
}

proc edit_exp_type_cancel {w} {
    global line_comval line_com line_val

    unset line_comval; # should remove traces
    unset line_com
    unset line_val
    destroy $w    
}

proc edit_exp_type_description {w lt} {
    switch $lt {
      CF	{set desc "cloning vector sequence file"; set eg "lorist2.vector"}
      CN	{set desc "clone name"; set eg "B0334"}
      CS	{set desc "cloning vector sequence present in sequence"; set eg ""}
      CV	{set desc "cloning vector type"; set eg "Lorist2"}
      DR	{set desc "direction of read"; set eg ""}
      DT	{set desc "date of experiment"; set eg ""}
      EN	{set desc "experiment name"; set eg ""}
      EX	{set desc "experimental notes"; set eg ""}
      FM	{set desc "sequencing vector fragmentation method"; set eg ""}
      LN	{set desc "local format trace file name"; set eg ""}
      LT	{set desc "local format trace file type"; set eg ""}
      MC	{set desc "machine on which experiment ran"; set eg ""}
      MN	{set desc "machine generated trace file name"; set eg ""}
      MT	{set desc "machine generated trace file type"; set eg ""}
      OP	{set desc "operator"; set eg "fbloggs@xyzzy"}
      PN	{set desc "primer name"; set eg ""}
      QR	{set desc "poor quality sequence present at right (3') end"; set eg ""}
      SC	{set desc "sequencing vector cloning site"; set eg "6249"}
      SF	{set desc "sequencing vector sequence file"; set eg "m13mp18.vector"}
      SI	{set desc "sequencing vector insertion length"; set eg "1400..2000"}
      SL	{set desc "sequencing vector present at left (5') end"; set eg ""}
      SP	{set desc "sequencing vector primer site (relative to cloning site)"; set eg ""}
      SQ	{set desc "sequence"; set eg ""}
      SR	{set desc "sequencing vector present at right (3') end"; set eg ""}
      ST	{set desc "strands"; set eg "2"}
      SV	{set desc "sequencing vector type"; set eg "M13mp18"}
      TN	{set desc "template name"; set eg ""}
      QL	{set desc "poor quality sequence present at left (5') end"; set eg ""}
      PS	{set desc "processing status"; set eg ""}
      CC	{set desc "comments"; set eg ""}
      SS	{set desc "sequence to screen against"; set eg ""}
      TG	{set desc "gel tag line"; set eg ""}
      ID	{set desc "identifier"; set eg ""}
      AQ	{set desc "average quality measure"; set eg ""}
      PR	{set desc "primer type"; set eg "1"}
      LI	{set desc "subclone library (mtd)"; set eg ""}
      LE	{set desc "subclone library entry (well)"; set eg ""}
      TC	{set desc "contig tag line"; set eg ""}
      AC	{set desc "accession number"; set eg ""}
      BC	{set desc "base calling software"; set eg ""}
      ON	{set desc "original base numbers (positions)"; set eg ""}
      AV	{set desc "accuracy (quality) values"; set eg ""}
      PC	{set desc "position in contig"; set eg ""}
      SE	{set desc "sense, whether it is complemented"; set eg ""}
      CL	{set desc "cloning vector left en"; set eg ""}
      CR	{set desc "cloning vector right en"; set eg ""}
      AP	{set desc "assembly position"; set eg ""}
      CH	{set desc "special chemistry used (eg taq)"; set eg "1"}
      PD	{set desc "primer data - the sequence of a primer"; set eg ""}
      WT	{set desc "wild type trace"; set eg "wild_type.scf"}
      default	{set desc ""; set eg ""}
    }
    if {$eg != ""} {
    $w.description configure \
	-text "Line type $lt : $desc. Example value : \"$eg\""
    } else {
	$w.description configure -text "Line type $lt : $desc"
    }
}
