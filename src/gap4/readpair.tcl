#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc ReadPairDialog { io f} {
    global gap_defs

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Find read pairs"

    ###########################################################################
    #input 
    lorf_in $f.infile [keylget gap_defs READPAIR.INFILE] \
	"" -bd 2 -relief groove

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "ReadPairs_OK_Pressed $io $f $f.infile"\
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Read Pairs}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################
    #final packing

    pack $f.infile -fill x
    pack $f.ok_cancel -fill x

}

proc ReadPairs_OK_Pressed {io f infile} {
    global gap_defs

    if {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $infile]
    }

    destroy $f

    ContigComparator $io

    SetBusy
    find_read_pairs -io $io -contigs $list
    ClearBusy
}
