#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc RescanOKPressed { score } {

    return [entrybox_get $score]
}


proc sip_rescan_matches { } {
    set f .rescan_matches
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f Rescan 
    global $f.ret

    entrybox $f.score \
	    -title Score \
	    -default 1\
	    -width 5\
	    -type "CheckIntMin 0"
    pack $f.score -anchor w

    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $f.button -bd 2 -relief groove \
	    -ok_command "global $f.ret; \
	                 set $f.ret \[entrybox_get $f.score\]; \
			 destroy $f"\
	    -cancel_command "destroy $f" \
	    -help_command "show_help spin {SPIN-Find similar spans}"
    pack $f.button -fill x

    tkwait variable $f.ret
    return [set $f.ret]
}
