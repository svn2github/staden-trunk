#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc Cons_SelMask { parent f } {
    global gap_defs

    keylset sm SELMASK [keylget gap_defs CONSENSUS.SELMASK]
    set b1 [keylget sm SELMASK.BUTTON.1]
    set b2 [keylget sm SELMASK.BUTTON.2]
    set b3 [keylget sm SELMASK.BUTTON.3]

    frame $f -bd 2 -relief groove
    button $f.but \
	    -text "Select tags" \
	    -command "TagDialog CONSENSUS.TAGS \
			$parent[keylget gap_defs SELECT_TAGS.WIN] {}"

    radiolist $f.rl \
	    -title [keylget sm SELMASK.NAME] \
	    -default [keylget sm SELMASK.VALUE]\
	    -buttons [format { \
	    { %s -command { %s configure -state normal; \
	    SetDefaultTags %s }} \
	    { %s -command { %s configure -state normal; \
	    SetDefaultTags %s }} \
	    { %s -command { %s configure -state disabled}} } \
	    [list $b1] [list $f.but] CONSENSUS.TAGS \
	    [list $b2] [list $f.but] CONSENSUS.TAGS \
	    [list $b3] [list $f.but] ]
    pack $f.rl -side left
    pack $f.but -side right
}

proc NormalDialog { io } {
    global gap_defs

    set f [keylget gap_defs CONSENSUS.1.WIN]
    global $f.format.expt.Radio
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Save consensus: normal"

    ###########################################################################
    #contig identifier widget
    contig_id $f.id -io $io 

    lorf_in $f.infile [keylget gap_defs CONSENSUS.INFILE] \
	"{contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state normal}
	" -bd 2 -relief groove

 
    ###########################################################################
    #select masking
    Cons_SelMask $f $f.sel_mask

    yes_no $f.pads \
	    -title "Strip pads" \
	    -relief groove -bd 2 -orient horizontal\
	    -default [keylget gap_defs CONSENSUS.NORMAL.STRIP_PADS]

    ###########################################################################
    #Reading annotations
    radiolist $f.annos \
	-title "Output reading annotations  " \
	-default [keylget gap_defs CONSENSUS.READ_ANNOTATIONS] \
	-orient horizontal \
	-buttons {{all} {{non-cutoff}} {none}}

    ###########################################################################
    #Reading notes
    yes_no $f.notes \
	    -title "Output reading notes" \
	    -relief groove -bd 2 -orient horizontal\
	    -default [keylget gap_defs CONSENSUS.READ_NOTES]

    radiolist $f.template \
	-title "Name consensus by" \
	-default [keylget gap_defs CONSENSUS.NAME_BY] \
	-orient horizontal \
	-buttons {{{left-most reading}} {{left-most template}}}

    ###########################################################################
    #format

    frame $f.format

    keylset ex EXPT [keylget gap_defs CONSENSUS.EXPT]
    set $f.format.expt.Radio [keylget ex EXPT.VALUE]
    keylset fo FORMAT [keylget gap_defs CONSENSUS.NORMAL.FORMAT]
    set b1 [keylget fo FORMAT.BUTTON.1]
    set b2 [keylget fo FORMAT.BUTTON.2]
    set b3 [keylget fo FORMAT.BUTTON.3]
    set b4 [keylget fo FORMAT.BUTTON.4]
 
    frame $f.format.dummy
    radiolist $f.format.main \
	    -title [keylget fo FORMAT.NAME] \
	    -default [keylget fo FORMAT.VALUE] \
	    -orient horizontal \
	    -buttons [format { \
	    { %s -command {radiolist_configure %s -state disabled;\
			   yes_no_configure %s -state disabled} } \
	    { %s -command {radiolist_configure %s -state disabled;\
			   yes_no_configure %s -state disabled} } \
	    { %s -command {radiolist_configure %s -state disabled;\
			   yes_no_configure %s -state disabled} } \
	    { %s -command {radiolist_configure %s -state normal;\
			   yes_no_configure %s -state normal} } }\
	    [list $b1] [list $f.annos] [list $f.notes] \
	    [list $b2] [list $f.annos] [list $f.notes] \
	    [list $b3] [list $f.annos] [list $f.notes] \
	    [list $b4] [list $f.annos] [list $f.notes] ] 

    pack $f.format.main -fill x
    pack $f.format.dummy -fill x

    ###########################################################################    #output file
    keylset op OUTPUT [keylget gap_defs CONSENSUS.OUTPUT]
    getFname $f.output [keylget op OUTPUT.NAME] save {} \
	[keylget op OUTPUT.VALUE]

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "Normal_OK_Pressed $f $io $f.infile $f.id \
	    $f.sel_mask.rl $f.pads $f.notes $f.template $f.annos \
	    $f.format.main $f.output" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Con-Normal}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################
    pack $f.infile -fill x
    pack $f.id -fill x
    pack $f.sel_mask -fill x
    pack $f.pads -fill x
    pack $f.template -fill x
    pack $f.format -fill x
    pack $f.annos -fill x
    pack $f.notes -fill x
    pack $f.output -fill x
    pack $f.ok_cancel -fill x

}

proc Normal_OK_Pressed {f io infile id sel_mask strippads notes template annos format output} {
    global gap_defs

    set gel_anno 0; #no gel annotations with expt file
    set truncate 1; #no gel annotations in hidden data with expt file
    set note_val 0
    set active_tags {}

    #single or list or file
    if {[lorf_in_get $infile] == 4} {
	set gel_name [contig_id_gel $id]
	set lreg [contig_id_lreg $id]
	set rreg [contig_id_rreg $id]
	SetContigGlobals $io $gel_name $lreg $rreg
	set list "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $infile]
    }

    if {$list  == ""} {
	raise $f
	return
    }

    set masking [radiolist_get $sel_mask]

    #if masking mode is 1 or 2 (mark or mask active tags)
    if {($masking == 1) || ($masking == 2)} {
	set active_tags [GetDefaultTags CONSENSUS.TAGS]
    }
    set out_format [radiolist_get $format]

    set strip [yes_no_get $strippads]

    #expt format chosen
    if { $out_format == 4 } {
	set expt [radiolist_get $annos]
	if {$expt == 1} {
	    set gel_anno 1
	    set truncate 0
	} elseif {$expt == 2 } {
	    set gel_anno 1
	    set truncate 1
	}
	set note_val [yes_no_get $notes]
    }
    #set out_file [entrybox_get $output.entry]
    set out_file [getFname_in_name $output]

    if {$out_file == ""} return

    SetBusy
    get_consensus -io $io \
	    -contigs $list \
	    -type normal\
	    -mask [lindex {mask mark none} [expr $masking-1]] \
	    -format $out_format\
	    -annotations $gel_anno \
            -truncate $truncate \
	    -notes $note_val \
	    -outfile $out_file \
	    -tag_types $active_tags \
	    -strip_pads $strip \
	    -name_format [radiolist_get $template]
    ClearBusy
    destroy $f
}

proc ExtendedDialog { io } {
    global gap_defs

    set f [keylget gap_defs CONSENSUS.2.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Save consensus: extended"

    ###########################################################################
    #contig identifier widget
    contig_id $f.id -io $io 

    lorf_in $f.infile [keylget gap_defs CONSENSUS.INFILE] \
	"{contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state normal}
	" -bd 2 -relief groove

    
    ###########################################################################
    #hidden data
    frame $f.hidden -bd 2 -relief groove

    HiddenParametersDialogInit $f.ops
    button $f.hidden.options \
	-text "Edit hidden data paramaters" \
	-command "HiddenParametersDialog $f $f.ops"

    pack $f.hidden.options -side top -anchor w


    ###########################################################################
    #select masking
    Cons_SelMask $f $f.sel_mask

    yes_no $f.pads \
	    -title "Strip pads" \
	    -relief groove -bd 2 -orient horizontal\
	    -default [keylget gap_defs CONSENSUS.EXTENDED.STRIP_PADS]

    ###########################################################################
    #format
    keylset fo FORMAT [keylget gap_defs CONSENSUS.EXTENDED.FORMAT]
    set b1 [keylget fo FORMAT.BUTTON.1]
    set b2 [keylget fo FORMAT.BUTTON.2]
 
    radiolist $f.format \
	    -title [keylget fo FORMAT.NAME] \
	    -default [keylget fo FORMAT.VALUE] \
	    -orient horizontal \
	    -buttons [format { \
	    { %s } {%s } }\
	    [list $b1] [list $b2] ] 

    ###########################################################################    #output file
    keylset op OUTPUT [keylget gap_defs CONSENSUS.OUTPUT]
    getFname $f.output [keylget op OUTPUT.NAME] save {} \
	[keylget op OUTPUT.VALUE]

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "Extended_OK_Pressed $f $io $f.infile $f.id \
	    $f.sel_mask.rl $f.pads $f.format $f.output $f.ops" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Con-Extended}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################
    pack $f.infile -fill x
    pack $f.id -fill x
    pack $f.hidden -fill x
    pack $f.sel_mask -fill x
    pack $f.pads -fill x
    pack $f.format -fill x
    pack $f.output -fill x
    pack $f.ok_cancel -fill x
}

proc Extended_OK_Pressed {f io infile id sel_mask strippads format output \
			  hidden_ops} {
    global gap_defs

    set gel_anno 0; #no gel annotations with expt file
    set truncate 0; #no gel annotations in hidden data with expt file
    set active_tags {}

    #list or file or single
    if {[lorf_in_get $infile] == 1} {
	set list [lorf_get_list $infile]
    } elseif {[lorf_in_get $infile] == 2} {
	set list [lorf_get_list $infile]
    } elseif {[lorf_in_get $infile] == 4} {
	set gel_name [contig_id_gel $id]
	set lreg [contig_id_lreg $id]
	set rreg [contig_id_rreg $id]
	SetContigGlobals $io $gel_name $lreg $rreg
	set list "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    }

    if {$list  == ""} {
	raise $f
	return
    }

    #Hidden data parameters
    upvar #0 $hidden_ops data
    set win_size $data($data(which_mode)_win_size)
    set max_dash $data(base_max_dash)
    set min_conf $data(conf_min_conf)
    set use_conf [lsearch {base conf} $data(which_mode)]

    set masking [radiolist_get $sel_mask]
    #if masking mode is 1 or 2 (mark or mask active tags)
    if {($masking == 1) || ($masking == 2)} {
	set active_tags [GetDefaultTags CONSENSUS.TAGS]
    }
    set out_format [radiolist_get $format]

    set strip [yes_no_get $strippads]

    #set out_file [entrybox_get $output]
    set out_file [getFname_in_name $output]

    if {$out_file == ""} return

    SetBusy
    get_consensus -io $io \
	    -contigs $list \
	    -type extended\
	    -mask [lindex {mask mark none} [expr $masking-1]] \
	    -format $out_format\
	    -win_size $win_size \
	    -max_dashes $max_dash \
	    -min_conf $min_conf \
	    -use_conf $use_conf \
	    -outfile $out_file \
	    -tag_types $active_tags \
	    -strip_pads $strip
    ClearBusy
    destroy $f
}

proc UnfinishedDialog { io } {
    global gap_defs

    set f [keylget gap_defs CONSENSUS.3.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Save consensus: unfinished"

    ###########################################################################
    #contig identifier widget
    contig_id $f.id -io $io 

    lorf_in $f.infile [keylget gap_defs CONSENSUS.INFILE] \
	"{contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state normal}
	" -bd 2 -relief groove

  
    yes_no $f.pads \
	    -title "Strip pads" \
	    -relief groove -bd 2 -orient horizontal\
	    -default [keylget gap_defs CONSENSUS.UNFINISHED.STRIP_PADS]

    ###########################################################################
    #format
    keylset fo FORMAT [keylget gap_defs CONSENSUS.UNFINISHED.FORMAT]
    set b1 [keylget fo FORMAT.BUTTON.1]
    set b2 [keylget fo FORMAT.BUTTON.2]
 
    radiolist $f.format \
	    -title [keylget fo FORMAT.NAME] \
	    -default [keylget fo FORMAT.VALUE] \
	    -orient horizontal \
	    -buttons [format { \
	    { %s } {%s } }\
	    [list $b1] [list $b2] ] 

    ###########################################################################    #output file
    keylset op OUTPUT [keylget gap_defs CONSENSUS.OUTPUT]
    getFname $f.output [keylget op OUTPUT.NAME] save {} \
	[keylget op OUTPUT.VALUE]

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "Unfinished_OK_Pressed $f $io $f.infile $f.id \
	    $f.pads $f.format $f.output" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Con-Unfinished}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################

    pack $f.infile -fill x
    pack $f.id -fill x
    pack $f.pads -fill x
    pack $f.format -fill x
    pack $f.output -fill x
    pack $f.ok_cancel -fill x

}
proc Unfinished_OK_Pressed {f io infile id strippads format output} {
    #list or file or single
    if {[lorf_in_get $infile] == 1} {
	set list [lorf_get_list $infile]
    } elseif {[lorf_in_get $infile] == 2} {
	set list [lorf_get_list $infile]
    } elseif {[lorf_in_get $infile] == 4} {
	set gel_name [contig_id_gel $id]
	set lreg [contig_id_lreg $id]
	set rreg [contig_id_rreg $id]
	SetContigGlobals $io $gel_name $lreg $rreg
	set list "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    }

    if {$list  == ""} {
	raise $f
	return
    }

    set out_format [radiolist_get $format]

    set strip [yes_no_get $strippads]

    #set out_file [entrybox_get $output]
    set out_file [getFname_in_name $output]

    if {$out_file == ""} return

    SetBusy
    get_consensus -io $io \
	    -contigs $list \
	    -type unfinished\
	    -format $out_format\
	    -outfile $out_file \
	    -strip_pads $strip
    ClearBusy
    destroy $f
}

proc QualityDialog { io } {
    global gap_defs

    set f [keylget gap_defs CONSENSUS.4.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Save consensus: quality"

    ###########################################################################
    #contig identifier widget
    contig_id $f.id -io $io 

    lorf_in $f.infile [keylget gap_defs CONSENSUS.INFILE] \
	"{contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state disabled}
	 {contig_id_configure $f.id -state normal}
	" -bd 2 -relief groove

    ###########################################################################
    #format
    keylset fo FORMAT [keylget gap_defs CONSENSUS.QUALITY.FORMAT]
    set b1 [keylget fo FORMAT.BUTTON.1]
 
    radiolist $f.format \
	    -title [keylget fo FORMAT.NAME] \
	    -default [keylget fo FORMAT.VALUE] \
	    -orient horizontal \
	    -buttons [format { { %s } } [list $b1] ] 

    ###########################################################################    #output file
    keylset op OUTPUT [keylget gap_defs CONSENSUS.OUTPUT]
    getFname $f.output [keylget op OUTPUT.NAME] save {} \
	[keylget op OUTPUT.VALUE]

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "Quality_OK_Pressed $f $io $f.id $f.infile $f.format\
	    $f.output" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Con-Quality}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################
    pack $f.infile -fill x
    pack $f.id -fill x
    pack $f.format -fill x
    pack $f.output -fill x
    pack $f.ok_cancel -fill x


}

proc Quality_OK_Pressed {f io id infile format output} {
    #list or file or single
    if {[lorf_in_get $infile] == 1} {
	set list [lorf_get_list $infile]
    } elseif {[lorf_in_get $infile] == 2} {
	set list [lorf_get_list $infile]
    } elseif {[lorf_in_get $infile] == 4} {
	set gel_name [contig_id_gel $id]
	set lreg [contig_id_lreg $id]
	set rreg [contig_id_rreg $id]
	SetContigGlobals $io $gel_name $lreg $rreg
	set list "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    }

    if {$list  == ""} {
	raise $f
	return
    }

    set out_format [radiolist_get $format]
    set out_file [getFname_in_name $output]

    if {$out_file == ""} return

    SetBusy
    get_consensus -io $io \
	    -contigs $list \
	    -type quality\
	    -format $out_format\
	    -outfile $out_file
    ClearBusy
    destroy $f
}
