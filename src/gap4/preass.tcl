#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc Preassemble {io} {
    global gap_defs

    set f [keylget gap_defs PREASS.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Enter pre-assembled data"

    ###########################################################################
    #file of reading names
    getFname $f.reads [keylget gap_defs PREASS.INPUT.NAME] load

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "PreAssemble_OkPressed $io $f \
	    \[getFname_in_name $f.reads \]" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Assembly-Pre}" \
	    -bd 2 \
	    -relief groove

    pack $f.reads -fill x
    pack $f.ok_cancel -fill x
}

proc PreAssemble_OkPressed { io f reads} {
    global FBRec
    global gap_defs

    if {![quit_displays $io "auto_assemble"]} {
        # Someone's too busy to shutdown?
        return
    }

    set lname [FileListName $reads]
    ListLoad $reads $lname
    SetBusy
    pre_assemble -io $io -files [ListToItems $lname]
    ClearBusy
    ListDelete $lname

    destroy $f
    #check that the database has contigs in it!
    if {[db_info num_contigs $io] > 0} {
	ActivateMenu_Open
	InitContigGlobals $io
	
	#draw the contig selector if it does not already exist
	set cs_win [keylget gap_defs CONTIG_SEL.WIN]
	if {![winfo exists $cs_win]} {
	    ContigSelector $io
	} else {
	    ContigInitReg $io
	    raise $cs_win
	}
    }
}
