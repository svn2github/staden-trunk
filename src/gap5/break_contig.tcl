#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc BreakContig {io} {
    global gap5_defs

    set t [keylget gap5_defs BREAK_CONTIG.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Break contig"

    contig_id $t.id -io $io -range 0 -command "BC_OK_Pressed $io $t $t.id $t.pos"
    xentry $t.pos -label "Base position"

    frame $t.holes -bd 2 -relief flat
    checkbutton $t.holes.b \
	-text "Remove contig holes" \
	-variable $t.RemoveHoles
    pack $t.holes.b -anchor w
    global $t.RemoveHoles
    set $t.RemoveHoles [keylget gap5_defs BREAK_CONTIG.REMOVE_HOLES]

    okcancelhelp $t.ok_cancel \
	    -ok_command "BC_OK_Pressed $io $t $t.id $t.holes $t.pos"\
	    -cancel_command "destroy $t" \
	    -help_command "show_help gap5 {Break Contig}" \
	    -bd 2 \
	    -relief groove
    
    ###########################################################################
    pack $t.id $t.holes $t.pos $t.ok_cancel -side top -fill x
}

proc BC_OK_Pressed {io t id holes pos} {
    global $t.RemoveHoles

    if {[set crec [contig_id_rec $id]] == ""} return
    if {[set pos [$pos get]] == ""} return

    destroy $t
    SetBusy

    if {![quit_displays -io $io -msg "break contig"]} {
	# Someone's too busy to shutdown?
	ClearBusy
	return
    }

    if {[catch {log_call break_contig \
		    -io $io \
		    -contig $crec \
		    -pos $pos \
		    -break_holes [set $t.RemoveHoles]} err]} {
	verror ERR_WARN break_contig $err
    }
    ContigSelector $io
    ContigInitReg $io

    ClearBusy
}
