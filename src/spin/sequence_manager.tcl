#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.


proc ReplaceSeq { listbox direction } {

    set index [$listbox curselection]
    seq_set_active_seq -index $index -direction $direction
}

#
# Creates a sequence list window
# 
proc sequence_list_create { } {

    set t .seq_sequence
    if {[xtoplevel $t -resizable 1] == ""} return
    wm title $t "Sequence manager"

    set font text_font
    listbox $t.l -yscrollcommand "$t.s set" -width 80 -selectmode single \
	    -font $font
    scrollbar $t.s -orient vertical -command "$t.l yview"

    frame $t.button
    button $t.button.h -text Horizontal -command "ReplaceSeq $t.l 0"
    button $t.button.v -text Vertical -command "ReplaceSeq $t.l 1"
    frame $t.bb
    button $t.bb.q -text "OK" -command "destroy $t"
    button $t.bb.h -text "Help" -command "show_help spin {SPIN-Sequence Manager}"

    bind $t.l <Motion> {
	%W selection clear 0 end
	%W selection set [%W index @%x,%y]
	%W selection anchor [%W index @%x,%y]
    }

    pack $t.bb -side bottom -fill x
    pack $t.bb.q -side left -expand 1
    pack $t.bb.h -side right -expand 1
    pack $t.s -side right -fill y
    pack $t.l -side left -expand 1 -fill both
}

#destroy a results list window
proc sequence_list_destroy { } {
    set t .seq_sequence
    if {[winfo exists $t]} {
	destroy $t
    }
}

#
# Displays a menu from a given sequence list item
#
proc sequence_list_popup {lbox seq_num x y} {
    if [winfo exists $lbox.m] {destroy $lbox.m}

    #header for popup menu is 1st or 2nd item in list box ie entryname
    #depending on whether it is active or not

    set sequences [sequence_names]
    set m [create_popup $lbox.m [lindex [lindex $sequences $seq_num] 1]] 

    set ops [seq_get_seq_ops -seq_num $seq_num]
    #puts "OPS $ops"
    set count 0
    foreach i $ops {
	if {$i == "SEPARATOR"} {
	    $m add separator
	} else {
	    if {$i != "PLACEHOLDER"} {
		if {$i == "Type"} {
		    SequenceTypeCascadeMenu $m $seq_num $count
		} else {
		    $m add command \
			-label $i \
			-command "destroy $m; \
			    seq_invoke_seq_op -seq_num $seq_num -job $count"
		}
	    }
	    incr count
	}
    }
    $lbox selection clear 0 end
    $lbox selection set $seq_num
    tk_popup $m $x $y
}

proc SequenceTypeCascadeMenu { m seq_num count } {
    #circle if 1; linear if 0 
    set seq_id [get_seq_id $seq_num]
    global $m.type$seq_id

    set $m.type$seq_id [seq_info $seq_id structure]

    $m add cascade -label "Sequence type" -menu $m.seq_type
    menu $m.seq_type

    $m.seq_type add radiobutton -label "linear" \
	-command "destroy $m; seq_invoke_seq_op -seq_num $seq_num -job $count -data 0" -variable $m.type$seq_id -value 0

    $m.seq_type add radiobutton -label "circular" \
	-command "destroy $m; seq_invoke_seq_op -seq_num $seq_num -job $count -data 1" -variable $m.type$seq_id -value 1
}


#
# Updates the results list window (always .seq_sequence)
# 
proc sequence_list_update { } {
    set t .seq_sequence; # toplevel window name

    if {![winfo exists $t]} {
	# Grab the registration list. Of format:
	set sequences [sequence_names]
	# Even if sequence manager window does not exist,
	# still need to check if menu state needs resetting.
	if {[llength $sequences] < 1} {
	    SeqActivateMenu_Initial
	}
	return
    }
    wm deiconify $t
    raise $t

    # Clear and remove old binding
    $t.l delete 0 end
    bind $t.l <<menu>> {}

    # Grab the registration list. Of format:
    set sequences [sequence_names]
    if {[llength $sequences] < 1} {
	SeqActivateMenu_Initial
    } elseif {[llength $sequences] == 1} {
	#ActivateMenu_Open1
    }

    foreach i $sequences {
	$t.l insert end [format "%s %-55s %s %6s %s %s" \
                            [lindex $i 0] \
                            [lindex $i 1] \
                            [lindex $i 2] \
                            [lindex $i 3] \
                            [lindex $i 4] \
			    [lindex $i 5]]
    }
    # Reenable our binding with our new list if > 0 sequences
    if {[llength $sequences] > 0} {
	bind $t.l <<menu>> {sequence_list_popup %W [%W index @%x,%y] %X %Y}
	#bind $t.l <<select>> {sequence_browse_seq %W 0}
	#bind $t.l <<use>> {sequence_browse_seq %W 1}
    }
}

#
# Selection binding - set global variable SeqName which can be used by other
# dialogues
#
proc sequence_browse_seq {path delete} {
    global SeqName

    set line [lindex [$path curselection] 0]
    set sequences [sequence_names]
    set name [lindex [lindex $sequences $line] 1]

    set SeqName $name
    #if called from <<use>> binding and has associated seq_id box then
    #destroy list box.

    if {$delete && [trace vinfo SeqName] != ""} {
	destroy [winfo parent $path]
    }
}

proc set_range_d {args} {
    global HORIZONTAL VERTICAL

    set t .set_range

    #invoked from Sequences menu
    if {$args != ""} {
	if {[modal $t -resizable 0] == ""} return
	#invoked from sequence manager
	set active_id $args
    } else {
	if {[xtoplevel $t -resizable 0] == ""} return
	set active_id [get_active_seq_id $HORIZONTAL]
	if {$active_id == -1} {
	    set active_id [get_active_seq_id $VERTICAL]
	}
	
	#if no active sequence is set, take the first sequence in list
	if {$active_id == -1} {
	    set active_id [get_seq_id -seq_num 0]
	} 
    }
    wm title $t "Set range"
    
    set seq_length [seq_info $active_id length]
    set seq_start [seq_info $active_id start] 
    set seq_end [seq_info $active_id end] 

    if {$args == ""} {
	seq_id $t.range -range 1 \
	    -default [seq_info $active_id name] -from 1 -to $seq_length \
	    -start_value $seq_start -end_value $seq_end\
	    -browse_cmd seq_browser
    } else {
	set cs [labelframe $t.range -text "Base range"]
	xtwinspin $cs.t -from $seq_start -to $seq_end
	pack $cs.t -fill both
	
#	scale_range $t.range -title "Base range" \
#	    -from 1 -to $seq_length -start_value $seq_start \
#	    -end_value $seq_end -min_value 1
    }
    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "set_range_d2 $t.range [list $args]; if {[list $args]==\"\"} {seq_id_destroy $t.range}; destroy $t" \
	-cancel_command "if {[list $args]==\"\"} {seq_id_destroy $t.range}; destroy $t" \
	-help_command "show_help spin {SPIN-Set Range}"
    
    pack $t.range -fill x
    pack $t.button -side bottom -fill x
}

proc set_range_d2 {range id} {

    if {$id == ""} {
	set id [name_to_seq_id [seq_id_name $range]]

	seq_set_range -seq_id $id -start [seq_id_from $range] \
	    -end [seq_id_to $range]
    } else {
	set re_enter 0
	set st_t [$range.t get_s]
	set en_t [$range.t get_e]
	if { $st_t > $en_t } {
	    tk_messageBox -icon error -type ok -title "Sequence identifer"\
		-message "You have chosen a start position \
		greater than the end position  "

	raise [winfo toplevel $range]

	tkwait variable re_enter   
	}
	seq_set_range -seq_id $id -start $st_t \
	    -end $en_t
    }
    sequence_list_update
}

##############################################################################
#sequence manager functions accessed from the main menu

proc file_save_d {args} {
    global HORIZONTAL VERTICAL

    set t .file_save

    if {$args == ""} {
	#invoked from Sequence menu
	if {[xtoplevel $t -resizable 0] == ""} return
	set active_id [get_active_seq_id $HORIZONTAL]
	
	if {$active_id == -1} {
	    set active_id [get_active_seq_id $VERTICAL]
	}
	
	#if no active sequence is set, take the first sequence in list
	if {$active_id == -1} {
	    set active_id [get_seq_id -seq_num 0]
	} 

	set seq_length [seq_info $active_id length] 
	set seq_start [seq_info $active_id start] 
	set seq_end [seq_info $active_id end] 

	seq_id $t.seq_id -range 1  \
	    -default [seq_info $active_id name] -from 1 -to $seq_length \
	    -start_value $seq_start -end_value $seq_end\
	    -browse_cmd sip_seq_browser
##############################################################################

    } else {
	#invoked from sequence manager
	if {[modal $t -resizable 0] == ""} return
	set active_id $args

	set seq_length [seq_info $active_id length] 
	set seq_start [seq_info $active_id start] 
	set seq_end [seq_info $active_id end] 
#	scale_range $t.seq_id -title "range" \
#	    -from 1 -to $seq_length -start_value $seq_start \
#	    -end_value $seq_end -min_value 1
#############################################################################

	set cs [labelframe $t.seq_id -text "Range"]
	xtwinspin $cs.t -from $seq_start -to $seq_end
	pack $cs.t -fill both
    }
    wm title $t Save

    set entry [getFname $t.file "Output filename" save]
    entrybox_configure $entry -width 20

    #####################################
    set b1 Fasta
    set b2 EMBL
    radiolist $t.format \
	-title "Save as" \
	-default 2 \
	-orient horizontal \
	-buttons [format {{%s} {%s}} \
		      [list $b1] [list $b2]]
    ######################################

    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "file_save_d2 $t $t.seq_id $t.file $t.format [list $args];\
                     if {[list $args]==\"\"} {seq_id_destroy $t.seq_id}"\
	-cancel_command "if {[list $args]==\"\"} {seq_id_destroy $t.seq_id}; destroy $t" \
	-help_command "show_help spin {SPIN-Save Sequence}"

    pack $t.seq_id -fill x
    pack $t.file -fill x
    pack $t.format -fill x
    pack $t.button -fill x
}

proc file_save_d2 {t seq_id file format id} {

    if {$id == ""} {
	set id [name_to_seq_id [seq_id_name $seq_id]]
	set from [seq_id_from $seq_id]
	set to [seq_id_to $seq_id]
    } else {
	set from [$seq_id.t get_s]
	set to [$seq_id.t get_e]
    }
    set file [getFname_in_name $file]
    set format [radiolist_get $format]

    seq_file_save -seq_id $id -start $from -end $to -format $format -file $file 
    destroy $t
}

proc file_delete_d { } {
    global HORIZONTAL VERTICAL

    set t .file_delete
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t Delete

    set active_id [get_active_seq_id $HORIZONTAL]
    if {$active_id == -1} {
	set active_id [get_active_seq_id $VERTICAL]
    }

    #if no active sequence is set, take the first sequence in list
    if {$active_id == -1} {
	set active_id [get_seq_id -seq_num 0]
    } 

    seq_id $t.seq_id -range 0 -default [seq_info $active_id name]\
	    -browse_cmd seq_browser

    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "file_delete_d2 $t $t.seq_id"\
	-cancel_command "seq_id_destroy $t.seq_id; destroy $t" \
	-help_command "show_help spin {SPIN-Delete Sequence}"

    pack $t.seq_id -fill x
    pack $t.button -fill x
}

proc file_delete_d2 {t seq_id} {

    set id [name_to_seq_id [seq_id_name $seq_id]]
#    puts $id
    seq_file_delete -seq_id $id
    seq_id_destroy $seq_id
    sequence_list_update
    destroy $t
}

proc complement_d { } {
    global HORIZONTAL VERTICAL

    set t .complement_seq
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t Complement

    set active_id [get_active_seq_id $HORIZONTAL]
    if {$active_id == -1} {
	set active_id [get_active_seq_id $VERTICAL]
    }

    #if no active sequence is set, take the first sequence in list
    if {$active_id == -1} {
	set active_id [get_seq_id -seq_num 0]
    } 

    seq_id $t.seq_id -range 0 -default [seq_info $active_id name]\
	    -browse_cmd seq_browser

    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "complement_d2 $t $t.seq_id"\
	-cancel_command "seq_id_destroy $t.seq_id; destroy $t" \
	-help_command "show_help spin {SPIN-Complement Sequence}"

    pack $t.seq_id -fill x
    pack $t.button -fill x
}

proc complement_d2 {t seq_id} {
    global DNA PROTEIN

    set id [name_to_seq_id [seq_id_name $seq_id]]
    if {[seq_info $id type] == $PROTEIN} {
	verror ERR_WARN "Complement sequence" "Unable to complement a protein sequence"
	return
    }

    seq_complement -seq_id $id
    sequence_list_update
    seq_id_destroy $seq_id
    destroy $t
}

proc interconvert_d { } {
    global HORIZONTAL VERTICAL

    set t .interconvert_seq
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Interconvert t and u"

    set active_id [get_active_seq_id $HORIZONTAL]
    if {$active_id == -1} {
	set active_id [get_active_seq_id $VERTICAL]
    }

    #if no active sequence is set, take the first sequence in list
    if {$active_id == -1} {
	set active_id [get_seq_id -seq_num 0]
    } 

    seq_id $t.seq_id -range 0 -default [seq_info $active_id name]\
	    -browse_cmd seq_browser

    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "interconvert_d2 $t $t.seq_id"\
	-cancel_command "seq_id_destroy $t.seq_id; destroy $t" \
	-help_command "show_help spin {SPIN-Interconvert t and u}"

    pack $t.seq_id -fill x
    pack $t.button -fill x
}

proc interconvert_d2 {t seq_id} {
    global PROTEIN
   
    set id [name_to_seq_id [seq_id_name $seq_id]]

    if {[seq_info $id type] == $PROTEIN} {
	verror ERR_WARN "Interconvert sequence" "Unable to interconvert t and u for a protein sequence"
	return
    }

    seq_interconvert -seq_id $id
    sequence_list_update
    seq_id_destroy $seq_id
    destroy $t
}

proc scramble_d { } {
    global HORIZONTAL VERTICAL

    set t .scramble_seq
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Scramble"

    set active_id [get_active_seq_id $HORIZONTAL]
    if {$active_id == -1} {
	set active_id [get_active_seq_id $VERTICAL]
    }

    #if no active sequence is set, take the first sequence in list
    if {$active_id == -1} {
	set active_id [get_seq_id -seq_num 0]
    } 

    seq_id $t.seq_id -range 0 -default [seq_info $active_id name]\
	    -browse_cmd seq_browser

    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "scramble_d2 $t $t.seq_id"\
	-cancel_command "seq_id_destroy $t.seq_id; destroy $t" \
	-help_command "show_help spin {SPIN-Scramble Sequence}"

    pack $t.seq_id -fill x
    pack $t.button -fill x
}

proc scramble_d2 {t seq_id} {
    
    set id [name_to_seq_id [seq_id_name $seq_id]]

    seq_scramble -seq_id $id
    sequence_list_update
    seq_id_destroy $seq_id
    destroy $t
}

proc set_horizontal_d { } {

    set t .set_horizontal
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Set horizontal sequence"

    #set default as first sequence in list
    set active_id [get_seq_id -seq_num 0]

    seq_id $t.seq_id -range 0 -default [seq_info $active_id name]\
	    -browse_cmd seq_browser

    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "set_horizontal_d2 $t $t.seq_id"\
	-cancel_command "seq_id_destroy $t.seq_id; destroy $t" \
	-help_command "show_help spin {SPIN-Change Active Sequence}"

    pack $t.seq_id -fill x
    pack $t.button -fill x
}

proc set_horizontal_d2 {t seq_id} {
    global HORIZONTAL

    set id [name_to_seq_id [seq_id_name $seq_id]]

    seq_set_active_seq -seq_id $id -direction $HORIZONTAL
    sequence_list_update
    seq_id_destroy $seq_id
    destroy $t
}

proc set_vertical_d { } {

    set t .set_vertical
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Set vertical sequence"

    #set default as first sequence in list
    set active_id [get_seq_id -seq_num 0]
    #set default as second sequence in list
    #set active_id [get_seq_id -seq_num 1]

    seq_id $t.seq_id -range 0 -default [seq_info $active_id name] \
	    -browse_cmd seq_browser

    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "set_vertical_d2 $t $t.seq_id"\
	-cancel_command "seq_id_destroy $t.seq_id; destroy $t" \
	-help_command "show_help spin {SPIN-Change Active Sequence}"

    pack $t.seq_id -fill x
    pack $t.button -fill x
}

proc set_vertical_d2 {t seq_id} {
    global VERTICAL

    set id [name_to_seq_id [seq_id_name $seq_id]]
    #puts $id

    #set direction [seq_info $seq_id direction]

    seq_set_active_seq -seq_id $id -direction $VERTICAL
    sequence_list_update
    seq_id_destroy $seq_id
    destroy $t
}

proc set_structure_d {structure } {

    set t .set_strucutre
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Set sequence structure"

    #set default as first sequence in list
    set active_id [get_seq_id -seq_num 0]

    seq_id $t.seq_id -range 0 -default [seq_info $active_id name]\
	    -browse_cmd seq_browser

    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "set_structure_d2 $t $t.seq_id $structure"\
	-cancel_command "seq_id_destroy $t.seq_id; destroy $t" \
	-help_command "show_help spin {SPIN-Change Active Sequence}"

    pack $t.seq_id -fill x
    pack $t.button -fill x
}

proc set_structure_d2 {t seq_id structure} {
    global HORIZONTAL

    set id [name_to_seq_id [seq_id_name $seq_id]]

    set_seq_structure -seq_id $id -structure $structure
    sequence_list_update
    seq_id_destroy $seq_id
    destroy $t
}

proc translate_d {args} {
    global HORIZONTAL VERTICAL

    set t .translate
#    if {[xtoplevel $t -resizable 0] == ""} return
#    wm title $t Translate

    if {$args != ""} {
	if {[modal $t -resizable 0] == ""} return
	set active_id $args
	frame $t.seq_id
    } else {
	if {[xtoplevel $t -resizable 0] == ""} return
	set active_id [get_active_seq_id $HORIZONTAL]

	if {$active_id == -1} {
	    set active_id [get_active_seq_id $VERTICAL]
	}

	#if no active sequence is set, take the first sequence in list
	if {$active_id == -1} {
	    set active_id [get_seq_id -seq_num 0]
	} 
	
	seq_id $t.seq_id -range 0 \
		-default [seq_info $active_id name] \
		-browse_cmd seq_browser
    }
    wm title $t Translate
    
    global $t.f1 $t.f2 $t.f3 $t.all
    set $t.f1 0
    set $t.f2 0
    set $t.f3 0
    set $t.all 0

    frame $t.frame
    checkbutton $t.frame.1 -text "frame 1" -variable $t.f1
    checkbutton $t.frame.2 -text "frame 2" -variable $t.f2
    checkbutton $t.frame.3 -text "frame 3" -variable $t.f3
    checkbutton $t.frame.all -text "all together" -variable $t.all
    
    pack $t.frame.1 $t.frame.2 $t.frame.3 $t.frame.all -side left -fill x

    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "translate_d2 $t $t.seq_id [list $args]; if {[list $args]==\"\"} {seq_id_destroy $t.seq_id}"\
	-cancel_command "if {[list $args]==\"\"} {seq_id_destroy $t.seq_id}; destroy $t" \
	-help_command "show_help spin {SPIN-Translate Sequence}"

    pack $t.seq_id $t.frame -fill x
    pack $t.button -fill x
}

proc translate_d2 {t seq_id id} {
    global $t.f1 $t.f2 $t.f3 $t.all PROTEIN

    if {$id == ""} {
	set id [name_to_seq_id [seq_id_name $seq_id]]
	seq_id_destroy $seq_id
    }

    if {[seq_info $id type] == $PROTEIN} {
	verror ERR_WARN "Translate sequence" "Unable to translate a protein sequence"
	return
    }

    seq_translate_seq -seq_id $id -f1 [set $t.f1] -f2 [set $t.f2] -f3 [set $t.f3] -all [set $t.all]
    sequence_list_update
    destroy $t
}

proc rotate_d {args} {
    global HORIZONTAL VERTICAL sip_defs

    set t .rotate

    if {$args != ""} {
	if {[modal $t -resizable 0] == ""} return
	set active_id $args
	#frame $t.seq_id
	set start [seq_info $active_id start]
	set end [seq_info $active_id end]

	#puts new-spinbox	
	xspinbox $t.seq_id \
	    -label "Origin" \
	    -width 7 \
	    -from $start \
	    -to $end

	$t.seq_id delete 0 end
	$t.seq_id insert end "$start"
    } else {
	if {[xtoplevel $t -resizable 0] == ""} return
	set active_id [get_active_seq_id $HORIZONTAL]
	if {$active_id == -1} {
	    set active_id [get_active_seq_id $VERTICAL]
	}

	#if no active sequence is set, take the first sequence in list
	if {$active_id == -1} {
	    set active_id [get_seq_id -seq_num 0]
	}
	set start [seq_info $active_id start]
	set end [seq_info $active_id end]
	seq_id $t.seq_id -range 1 -default [seq_info $active_id name]\
	    -update_cmd [list [list RotateUpdate $t.seq_id $t.origin]] \
	    -browse_cmd seq_browser \
	    -start_value $start \
	    -end_value $end \
	    -single yes \
	    -label_single "Origin"
    }
#	seq_id $t.seq_id -range 0 -default [seq_info $active_id name]\
#	    -update_cmd [list [list RotateUpdate $t.seq_id $t.origin]] \
#	    -browse_cmd seq_browser \
#	    -start_value $start \
#	    -end_value $end

    wm title $t "Rotate"

  
#    scalebox $t.origin -title [keylget sip_defs SIP.ROTATE.ORIGIN.NAME] \
#	-default [keylget sip_defs SIP.ROTATE.ORIGIN.VALUE] \
#	-width 6 -from $start -to $end -type CheckInt -orient horizontal\
    
    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "rotate_d2 $t $t.seq_id [list $args]; if {[list $args]==\"\"} {seq_id_destroy $t.seq_id}"\
	-cancel_command "if {[list $args]==\"\"} {seq_id_destroy $t.seq_id}; destroy $t" \
	-help_command "show_help spin {SPIN-Rotate Sequence}"

    pack $t.seq_id -fill x
    pack $t.button -fill x
}

proc rotate_d2 {t seq_id id} {
    
    if {$id == ""} {
	set id [name_to_seq_id [seq_id_name $seq_id]]
	set ori [seq_id_single $seq_id]
	seq_id_destroy $seq_id
    } else {
	set ori [$seq_id get]
    }

    if {[seq_info $id is_sub_seq] == 1} {
	verror ERR_WARN Rotate "Unable to rotate a sub-sequence"
	return
    }
      
    seq_rotate -seq_id $id -origin $ori

    sequence_list_update
    destroy $t
}

proc RotateUpdate {seq_id origin} {

    set id [name_to_seq_id [seq_id_name $seq_id]]

#    set start [seq_info $id start]
#    set end [seq_info $id end]
#    scalebox_configure $origin -from [seq_info $id start]\
#	-to [seq_info $id end]
}
proc copy_range_d {args} {
    global HORIZONTAL VERTICAL

    set t .copy_range

    #invoked from Sequences menu
    if {$args != ""} {
	if {[modal $t -resizable 0] == ""} return
	#invoked from sequence manager
	set active_id $args
    } else {
	if {[xtoplevel $t -resizable 0] == ""} return
	set active_id [get_active_seq_id $HORIZONTAL]
	if {$active_id == -1} {
	    set active_id [get_active_seq_id $VERTICAL]
	}
	
	#if no active sequence is set, take the first sequence in list
	if {$active_id == -1} {
	    set active_id [get_seq_id -seq_num 0]
	} 
    }
    wm title $t "Copy range"
    
    set seq_length [seq_info $active_id length]
    set seq_start [seq_info $active_id start] 
    set seq_end [seq_info $active_id end] 

    if {$args == ""} {
	seq_id $t.range -range 1 \
	    -default [seq_info $active_id name] -from 1 -to $seq_length \
	    -start_value $seq_start -end_value $seq_end\
	    -browse_cmd seq_browser
    } else {
	set cs [labelframe $t.range -text "Base range"]
	xtwinspin $cs.t -from $seq_start -to $seq_end
	pack $cs.t -fill both
	
#	scale_range $t.range -title "Base range" \
#	    -from 1 -to $seq_length -start_value $seq_start \
#	    -end_value $seq_end -min_value 1
    }
    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "copy_range_d2 $t.range [list $args]; if {[list $args]==\"\"} {seq_id_destroy $t.range}; destroy $t" \
	-cancel_command "if {[list $args]==\"\"} {seq_id_destroy $t.range}; destroy $t" \
	-help_command "show_help spin {SPIN-Set Range}"
    
    pack $t.range -fill x
    pack $t.button -side bottom -fill x
}

proc copy_range_d2 {range id} {

    if {$id == ""} {
	set id [name_to_seq_id [seq_id_name $range]]

	seq_copy_range -seq_id $id -start [seq_id_from $range] \
	    -end [seq_id_to $range]
    } else {
	set re_enter 0
	set st_t [$range.t get_start]
	set en_t [$range.t get_end]
	if { $st_t > $en_t } {
	    tk_messageBox -icon error -type ok -title "Sequence identifer"\
		-message "You have chosen a start position \
		greater than the end position  "

	raise [winfo toplevel $range]

	tkwait variable re_enter   
	}
	seq_copy_range -seq_id $id -start $st_t -end $en_t
    }
    sequence_list_update
}

proc seq_range_updates {range} {
    global old_id

    set seq_id [name_to_seq_id [seq_id_name $range]]
    global $seq_id.start $seq_id.end

    if {![info exists old_id]} {
	set old_id $seq_id
    }

    if {$old_id != $seq_id} {
	if {[info exists $seq_id.start]} {
	    seq_id_configure $range -start_value [set $seq_id.start]
	}    
	if {[info exists $seq_id.end]} {
	    seq_id_configure $range -end_value [set $seq_id.end]
	}
    }
    set old_id $seq_id
}

#browse command for seq_id to invoke sequence manager
proc seq_browser {path} {
    global CurFrame

    set CurFrame $path
    sequence_list_create 
    sequence_list_update 
}

#browse command for seq_id to invoke sip sequence manager
proc sip_seq_browser {path} {
    global CurFrame

    set CurFrame $path
  #  sip_sequence_list_create 
  #  sip_sequence_list_update
  #  yy
     sequence_list_create 
     sequence_list_update
      }

#browse command for seq_id to invoke nip sequence manager
proc nip_seq_browser {path} {

 #   nip_sequence_list_create 
 #   nip_sequence_list_update
 # yy 
     sequence_list_create 
     sequence_list_update
}


