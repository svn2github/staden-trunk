#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc ComplementContig {io} {
    global gap_defs

    set l [keylget gap_defs COMPLEMENT_CONTIG]
    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Complement contig"

    contig_id $t.id \
	    -command "ComplementContig2 $io $t $t.id" \
	    -io $io \
	    -range 0

    okcancelhelp $t.ok \
	-ok_command "ComplementContig2 $io $t $t.id" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Complement}" \
	-bd 2 -relief groove

    pack $t.id $t.ok -side top -fill x
}

proc ComplementContig2 {io t id} {
    if {[set c [contig_id_gel $id]] == ""} {bell; return}

    destroy $t
    update idletasks

    set cnum [db_info get_contig_num $io $c]
    complement_contig -io $io -contigs "$c"
    SetContigGlobals $io [left_gel $io $cnum]
}
