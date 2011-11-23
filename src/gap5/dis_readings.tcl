#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#-----------------------------------------------------------------------------
# Disassemble Readings

proc DisReadings { io } {
    global gap5_defs

    set f [keylget gap5_defs DIS_READINGS.WIN]
    set cs [keylget gap5_defs CONTIG_SEL.WIN]

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Disassemble readings"

    contig_id $f.id -title "Reading identifier" -io $io -range 0

    #create file of filenames entry windows
    lorf_in $f.infile [keylget gap5_defs DIS_READINGS.INFILE] \
	"{contig_id_configure $f.id -state disabled} \
	 {contig_id_configure $f.id -state disabled}\
	 {contig_id_configure $f.id -state normal}" -bd 2 -relief groove


    ###########################################################################
    frame $f.break -bd 2 -relief groove
    checkbutton $f.break.b \
	-text "Break contig at holes" \
	-variable $f.Break
    global $f.Break
    set $f.Break [keylget gap5_defs DIS_READINGS.BREAK_CONTIG]
    pack $f.break.b -anchor w

    ###########################################################################
    keylset st SELTASK [keylget gap5_defs DIS_READINGS.SELTASK]
    set b1 [keylget st SELTASK.BUTTON.1]
    set b2 [keylget st SELTASK.BUTTON.2]
    set b3 [keylget st SELTASK.BUTTON.3]

    radiolist $f.sel_task \
	    -title [keylget st SELTASK.NAME] \
	    -bd 2 \
	    -relief groove \
	    -default [keylget st SELTASK.VALUE] \
	-buttons [format { {%s} {%s} {%s}} \
		      [list $b1] [list $b2] [list $b3]]

    ###########################################################################
    # Whether to duplicate tags
    xyn $f.constags \
	-label "Duplicate consensus tags?" \
	-orient horiz \
	-default [keylget gap5_defs DIS_READINGS.DUPCONSTAGS.VALUE]


    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
        -ok_command "OK_Pressed_DisReading $io $f $cs $f.infile $f.id $f.sel_task $f.constags"\
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 {Disassemble}" \
	-bd 2 \
	-relief groove
    ###########################################################################

    pack $f.infile -fill x
    pack $f.id -fill x
    pack $f.break -fill x
    pack $f.sel_task -fill x
    # still sometimes duplicates tags when requested not to, so option
    # is disabled for now.
    # pack $f.constags -fill x
    pack $f.ok_cancel -fill x

}

proc OK_Pressed_DisReading { io f cs infile id sel_task constags } {
    global gap5_defs $f.Break

    set iopt [expr {[radiolist_get $sel_task]-1}]
    set dup_tags [$constags get]
    
    #special case for a single reading
    if {[lorf_in_get $infile] == 3} {
	set list [contig_id_gel $id]
    } else {
	if {[set list [lorf_get_list $infile]] == ""} {bell; return}
    }

    #convert reading numbers into reading names
#    set list [eval get_read_names -io $io $list]
    destroy $f
    update idletasks

    if {![quit_displays -io $io -msg "disassemble_readings"]} {
	# Someone's too busy to shutdown?
	return
    }

    set result [disassemble_readings \
		    -io $io \
		    -readings $list \
		    -move $iopt \
		    -duplicate_tags $dup_tags \
		    -remove_holes [set $f.Break]]
    #if database is empty, destroy contig selector and set menus back to
    #as if opened new database
    if {[db_info num_contigs $io] == 0} {
	set cs_win [keylget gap5_defs CONTIG_SEL.WIN]
	destroy $cs_win
	DisableMenu_Open
	ActivateMenu_New
	return
    }

    ContigInitReg $io
    InitContigGlobals $io

    # Called to increase maxseq
    PrintDatabaseInfo $io
}

#
# For disassembly when called from the editor - we already have a list in
# this case so there's fewer dialogue bits to set
#
proc DisEditorReadingsDialog { io list f } {
    global gap5_defs

    if {[modal $f] == ""} return
    wm title $f "Disassemble readings"
    

    ###########################################################################
    keylset st SELTASK [keylget gap5_defs DIS_READINGS.SELTASK]
    set b1 [keylget st SELTASK.BUTTON.1]
    set b2 [keylget st SELTASK.BUTTON.2]
    set b3 [keylget st SELTASK.BUTTON.3]

    radiolist $f.sel_task \
	    -title [keylget st SELTASK.NAME] \
	    -bd 2 \
	    -relief groove \
	    -default [keylget st SELTASK.VALUE] \
	    -buttons [format { {%s} {%s} {%s} } \
			  [list $b1] [list $b2] [list $b3]]

    ###########################################################################
    # Whether to duplicate tags
    xyn $f.constags \
	-label "Duplicate consensus tags?" \
	-orient horiz \
	-default [keylget gap5_defs DIS_READINGS.DUPCONSTAGS.VALUE]

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
        -ok_command "OK_Pressed_EdDisReading $io {$list} $f $f.sel_task $f.constags"\
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 {Disassemble}" \
	-bd 2 \
	-relief groove
    ###########################################################################

    pack $f.sel_task -fill x
    # still sometimes duplicates tags when requested not to, so option
    # is disabled for now.
    # pack $f.constags -fill x; 
    pack $f.ok_cancel -fill x

}

proc OK_Pressed_EdDisReading { io list f sel_task constags } {
    global gap5_defs

    set iopt [expr {[radiolist_get $sel_task]-1}]
    set dup_tags [$constags get]
 
    #convert reading numbers into reading names
    set list [eval get_read_names -io $io $list]
    destroy $f
    update idletasks

    if {![quit_displays -io $io -msg "disassemble_readings"]} {
	# Someone's too busy to shutdown?
	return
    }
    set result [disassemble_readings \
		    -io $io \
		    -readings $list \
		    -move $iopt \
		    -duplicate_tags $dup_tags]
    if { [string compare $result "OK"] == 0} {

	#if database is empty, destroy contig selector and set menus back to
	#as if opened new database
	if {[db_info num_contigs $io] == 0} {
	    set cs_win [keylget gap5_defs CONTIG_SEL.WIN]
	    destroy $cs_win
	    DisableMenu_Open
	    ActivateMenu_New
	    return
	}
    }

    ContigInitReg $io
}

#-----------------------------------------------------------------------------
# Disassemble Contigs

proc DisContigs { io } {
    global gap5_defs

    set f [keylget gap5_defs DIS_CONTIGS.WIN]
    set cs [keylget gap5_defs CONTIG_SEL.WIN]

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Delete contigs"

    contig_id $f.id -title "Contig identifier" -io $io -range 0

    #create file of filenames entry windows
    lorf_in $f.infile [keylget gap5_defs DIS_READINGS.INFILE] \
	"{contig_id_configure $f.id -state disabled} \
	 {contig_id_configure $f.id -state disabled}\
	 {contig_id_configure $f.id -state normal}" -bd 2 -relief groove

    label $f.warn -text "WARNING: this removes entire contigs. To remove individual sequences please use Disassemble Readings" -fg red -wrap 400


    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
        -ok_command "OK_Pressed_DisContig $io $f $cs $f.infile $f.id"\
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 {Delete Contigs}" \
	-bd 2 \
	-relief groove
    ###########################################################################

    pack $f.infile -fill x
    pack $f.id -fill x
    pack $f.warn
    pack $f.ok_cancel -fill x

}

proc OK_Pressed_DisContig { io f cs infile id } {
    global gap5_defs

    #special case for a single reading
    if {[lorf_in_get $infile] == 3} {
	set list [contig_id_gel $id]
    } else {
	if {[set list [lorf_get_list $infile]] == ""} {bell; return}
    }

    destroy $f
    update idletasks

    if {![quit_displays -io $io -msg "disassemble_contigs"]} {
	# Someone's too busy to shutdown?
	return
    }

    set result [disassemble_contigs \
		    -io $io \
		    -contigs $list]
    #if database is empty, destroy contig selector and set menus back to
    #as if opened new database
    if {[db_info num_contigs $io] == 0} {
	set cs_win [keylget gap5_defs CONTIG_SEL.WIN]
	destroy $cs_win
	DisableMenu_Open
	ActivateMenu_New
	return
    }

    ContigInitReg $io
    InitContigGlobals $io

    # Called to increase maxseq
    PrintDatabaseInfo $io
}

#-----------------------------------------------------------------------------
# Remove Contig Holes
proc RemoveContigHoles {io} {
    global gap5_defs

    set t .remove_contig_holes
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Remove Contig Holes"

    contig_id $t.id -io $io
    lorf_in $t.infile [keylget gap5_defs CONSENSUS.INFILE] \
        "{contig_id_configure $t.id -state disabled}
         {contig_id_configure $t.id -state disabled}
         {contig_id_configure $t.id -state disabled}
         {contig_id_configure $t.id -state normal}" \
	-bd 2 -relief groove

    okcancelhelp $t.ok_cancel \
	    -ok_command "RemoveContigHoles2 $io $t"\
	    -cancel_command "destroy $t" \
	    -help_command "show_help gap5 {Remove Contig Holes}" \
	    -bd 2 \
	    -relief groove

    pack $t.infile $t.id $t.ok_cancel -side top -fill x
}

proc RemoveContigHoles2 {io t} {
    if {[lorf_in_get $t.infile] == 4} {
	set gel_name [contig_id_gel $t.id]
	set lreg [contig_id_lreg $t.id]
	set rreg [contig_id_rreg $t.id]
	set list [list [list $gel_name $lreg $rreg]]
    } elseif {[lorf_in_get $t.infile] == 3} {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $t.infile]
    }

    if {![quit_displays -io $io -msg "break_contig_holes"]} {
	# Someone's too busy to shutdown?
	return
    }
    destroy $t

    break_contig_holes -io $io -contigs $list

    if {[db_info num_contigs $io] == 0} {
	set cs_win [keylget gap5_defs CONTIG_SEL.WIN]
	destroy $cs_win
	DisableMenu_Open
	ActivateMenu_New
	return
    }

    ContigInitReg $io
}