#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc CheckAssembly {io} {
    global gap5_defs
 
    set l [keylget gap5_defs CHECK_ASSEMBLY]
    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Check assembly"

    contig_id $t.id -io $io -range 0

    lorf_in $t.infile [keylget gap5_defs CHECK_ASSEMBLY.INFILE] \
	"{contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state normal}
	" -bd 2 -relief groove

    frame $t.mp -bd 2 -relief groove
    scalebox $t.mp.max_perc \
	-bd 2 -relief groove \
	-orient horizontal -width 5 \
	-title   "[keylget l MAXPERC.NAME]" \
	-from    "[keylget l MAXPERC.MIN]" \
	-to      "[keylget l MAXPERC.MAX]" \
	-default "[keylget l MAXPERC.VALUE]"

    keylset ws WINSIZE [keylget gap5_defs CHECK_ASSEMBLY.USED.WINSIZE]
    scalebox $t.mp.win_size \
            -title [keylget ws WINSIZE.NAME] \
            -orient horizontal \
            -to [keylget ws WINSIZE.MAX] \
            -from [keylget ws WINSIZE.MIN] \
            -default [keylget ws WINSIZE.VALUE]\
	    -bd 2 -relief groove \
            -width 5 \
            -type CheckInt

    xyn $t.mp.ignore_N \
	-label "Ignore N bases" \
	-default 1

    pack $t.mp.max_perc $t.mp.win_size $t.mp.ignore_N -side top -fill x

    okcancelhelp $t.ok \
	-ok_command "CheckAssembly2 $io $t $t.infile $t.id \
			$t.mp.max_perc $t.mp.win_size $t.mp.ignore_N" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Check Assembly}" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.mp $t.ok -side top -fill x
}

proc CheckAssembly2 {io t lorf id maxperc u_winsize ignore_N} {
    global gap5_defs
    global NGRec

    set args ""

    if {[set mperc   [scalebox_get $maxperc]] == ""} {bell; return}
    if {[set wsize   [scalebox_get $u_winsize]] == ""} {bell; return}
    set ignore_N [$ignore_N get]

    if {[lorf_in_get $lorf] == 4} {
        if {[set contign [contig_id_gel $id]]  == ""} {bell; return}
	set list $contign
	SetContigGlobals $io $contign
    } elseif {[lorf_in_get $lorf] == 3} {
	set list [CreateAllContigList $io]
    } else {
	if {[set list [lorf_get_list $lorf]] == ""} {bell; return}
    }

    destroy $t
    SetBusy
    update idletasks
    ContigComparator $io

    check_assembly \
	-io $io \
	-max_pmismatch $mperc \
	-contigs $list \
	-win_size $wsize \
	-ignore_N $ignore_N

    ClearBusy
}
