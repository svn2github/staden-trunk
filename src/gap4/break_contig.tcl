#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc BreakContig {io} {
    global gap_defs

    set t [keylget gap_defs BREAK_CONTIG.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Break contig"

    label $t.l -text "Reading to be at the left end of the contig:"
    contig_id $t.id \
	    -command "BC_OK_Pressed $io $t $t.id"\
	    -io $io \
	    -range 0

    okcancelhelp $t.ok_cancel \
	    -ok_command "BC_OK_Pressed $io $t $t.id"\
	    -cancel_command "destroy $t" \
	    -help_command "show_help gap4 {Break Contig}" \
	    -bd 2 \
	    -relief groove
    
    ###########################################################################
    pack $t.l $t.id $t.ok_cancel -side top -fill x
}

proc BC_OK_Pressed {io t id} {
    if {[set g [contig_id_gel $id]] == ""} {bell; return}

    destroy $t
    SetBusy

    if {![quit_displays $io "break contig"]} {
	# Someone's too busy to shutdown?
	ClearBusy
	return
    }
    break_contig -io $io -readings $g
    ContigInitReg $io

    ClearBusy
}
