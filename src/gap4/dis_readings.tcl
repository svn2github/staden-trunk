#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc DisReadingsDialog { io f cs} {
    global gap_defs

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Disassemble readings"

    contig_id $f.id -title "Reading identifier" -io $io -range 0

    #create file of filenames entry windows
    lorf_in $f.infile [keylget gap_defs DIS_READINGS.INFILE] \
	"{contig_id_configure $f.id -state disabled} \
	 {contig_id_configure $f.id -state disabled}\
	 {contig_id_configure $f.id -state normal}" -bd 2 -relief groove


    ###########################################################################
    keylset st SELTASK [keylget gap_defs DIS_READINGS.SELTASK]
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
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
        -ok_command "OK_Pressed_DisReading $io $f $cs $f.infile $f.id $f.sel_task"\
	-cancel_command "destroy $f" \
	-help_command "show_help gap4 {Disassemble}" \
	-bd 2 \
	-relief groove
    ###########################################################################

    pack $f.infile -fill x
    pack $f.id -fill x
    pack $f.sel_task -fill x
    pack $f.ok_cancel -fill x

}

proc OK_Pressed_DisReading { io f cs infile id sel_task } {
    global gap_defs

    set iopt [expr {[radiolist_get $sel_task]-1}]

    #special case for a single reading
    if {[lorf_in_get $infile] == 3} {
	set list [contig_id_gel $id]
    } else {
	if {[set list [lorf_get_list $infile]] == ""} {bell; return}
    }

    #convert reading numbers into reading names
    set list [eval get_read_names -io $io $list]
    destroy $f
    update idletasks

    if {![quit_displays $io "disassemble_readings"]} {
	# Someone's too busy to shutdown?
	return
    }
    set result [disassemble_readings -io $io -readings $list -move $iopt]
    #if database is empty, destroy contig selector and set menus back to
    #as if opened new database
    if {[db_info num_contigs $io] == 0} {
	set cs_win [keylget gap_defs CONTIG_SEL.WIN]
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
    global gap_defs

    if {[modal $f] == ""} return
    wm title $f "Disassemble readings"
    

    ###########################################################################
    keylset st SELTASK [keylget gap_defs DIS_READINGS.SELTASK]
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
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
        -ok_command "OK_Pressed_EdDisReading $io {$list} $f $f.sel_task"\
	-cancel_command "destroy $f" \
	-help_command "show_help gap4 {Disassemble}" \
	-bd 2 \
	-relief groove
    ###########################################################################

    pack $f.sel_task -fill x
    pack $f.ok_cancel -fill x

}

proc OK_Pressed_EdDisReading { io list f sel_task } {
    global gap_defs

    set iopt [expr {[radiolist_get $sel_task]-1}]
 
    #convert reading numbers into reading names
    set list [eval get_read_names -io $io $list]
    destroy $f
    update idletasks

    if {![quit_displays $io "disassemble_readings"]} {
	# Someone's too busy to shutdown?
	return
    }
    set result [disassemble_readings -io $io -readings $list -move $iopt]
    if { [string compare $result "OK"] == 0} {

	#if database is empty, destroy contig selector and set menus back to
	#as if opened new database
	if {[db_info num_contigs $io] == 0} {
	    set cs_win [keylget gap_defs CONTIG_SEL.WIN]
	    destroy $cs_win
	    DisableMenu_Open
	    ActivateMenu_New
	    return
	}
    }
    ContigInitReg $io
    InitContigGlobals $io

    # Called to increase maxseq
    PrintDatabaseInfo $io
}
