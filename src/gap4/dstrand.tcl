#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc DoubleStrand {io} {
    global gap_defs

    set l [keylget gap_defs DOUBLE_STRAND]
    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Double strand"

    contig_id $t.id -io $io -range 0

    lorf_in $t.infile [keylget gap_defs DOUBLE_STRAND.INFILE] \
	"{contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state normal}
	" -bd 2 -relief groove

    entrybox $t.maxmis \
	-title   "[keylget l MAXMIS.NAME]" \
	-default "[keylget l MAXMIS.VALUE]" \
	-width 5 \
	-type "CheckIntMin 0"

    scalebox $t.maxperc \
	-orient horizontal -width 5 \
	-title   "[keylget l MAXPERC.NAME]" \
	-from    "[keylget l MAXPERC.MIN]" \
	-to      "[keylget l MAXPERC.MAX]" \
	-default "[keylget l MAXPERC.VALUE]"

    okcancelhelp $t.ok \
	-ok_command "DoubleStrand2 $io $t $t.maxmis $t.maxperc \
			$t.infile $t.id" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Double Strand}" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.maxmis $t.maxperc $t.ok \
	-side top -fill x
}

proc DoubleStrand2 {io t maxmis maxperc lorf id} {
    if {[set maxm [entrybox_get $maxmis]] == ""} {bell; return}
    if {[set maxp [scalebox_get $maxperc]]  == ""} {bell; return}

    if {[lorf_in_get $lorf] == 4} {
        if {[set contign [contig_id_gel $id]] == ""} {bell; return}
#      	if {[set lreg  [contig_id_lreg $id]] == ""} {bell; return}
#	if {[set rreg  [contig_id_rreg $id]] == ""} {bell; return}
#	set list "{$contign $lreg $rreg}"
	set list "{$contign}"
	SetContigGlobals $io $contign
    } elseif {[lorf_in_get $lorf] == 3} {
	set list [CreateAllContigList $io]
    } else {
	if {[set list [lorf_get_list $lorf]] == ""} {bell; return}
    }

    destroy $t
#    update idletasks

    SetBusy
    double_strand \
	-io $io \
	-max_nmismatch $maxm \
	-max_pmismatch $maxp \
	-contigs $list
    ClearBusy
}
