#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc CheckAssembly {io} {
    global gap_defs
 
    set l [keylget gap_defs CHECK_ASSEMBLY]
    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Check assembly"

    contig_id $t.id -io $io -range 0

    lorf_in $t.infile [keylget gap_defs CHECK_ASSEMBLY.INFILE] \
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

    frame $t.cf -bd 2 -relief groove
    entrybox $t.cf.min_len \
	-title   "[keylget l MINLEN.NAME]" \
	-default "[keylget l MINLEN.VALUE]" \
	-width 5 \
	-type "CheckIntMin 0"

    keylset ws WINSIZE [keylget gap_defs CHECK_ASSEMBLY.HIDDEN.WINSIZE]
    scalebox $t.cf.win_size \
            -title [keylget ws WINSIZE.NAME] \
            -orient horizontal \
            -to [keylget ws WINSIZE.MAX] \
            -from [keylget ws WINSIZE.MIN] \
            -default [keylget ws WINSIZE.VALUE]\
	    -bd 2 -relief groove \
            -width 5 \
            -type CheckInt

    keylset md MAXDASH [keylget gap_defs CHECK_ASSEMBLY.HIDDEN.MAXDASH]
    scalebox $t.cf.max_dash \
            -title [keylget md MAXDASH.NAME] \
            -orient horizontal \
            -to [keylget md MAXDASH.MAX] \
            -from [keylget md MAXDASH.MIN]\
            -default [keylget md MAXDASH.VALUE]\
	    -bd 2 -relief groove \
            -width 5 \
            -type CheckInt

    frame $t.us -bd 2 -relief groove
    keylset ws WINSIZE [keylget gap_defs CHECK_ASSEMBLY.USED.WINSIZE]
    scalebox $t.us.win_size \
            -title [keylget ws WINSIZE.NAME] \
            -orient horizontal \
            -to [keylget ws WINSIZE.MAX] \
            -from [keylget ws WINSIZE.MIN] \
            -default [keylget ws WINSIZE.VALUE]\
	    -bd 2 -relief groove \
            -width 5 \
            -type CheckInt

    yes_no $t.cutoff \
	-title   "[keylget l USE_CUTOFF.NAME]" \
	-default "[keylget l USE_CUTOFF.VALUE]" \
	-bd 2 \
	-relief groove \
	-orient horizontal \
	-ycommand "entrybox_configure $t.cf.min_len  -state normal;
		   scalebox_configure $t.cf.win_size -state normal;
		   scalebox_configure $t.cf.max_dash -state normal;
		   scalebox_configure $t.us.win_size -state disabled" \
	-ncommand "entrybox_configure $t.cf.min_len  -state disabled;
		   scalebox_configure $t.cf.win_size -state disabled;
		   scalebox_configure $t.cf.max_dash -state disabled;
		   scalebox_configure $t.us.win_size -state normal"

    okcancelhelp $t.ok \
	-ok_command "CheckAssembly2 $io $t $t.infile $t.id \
			$t.cutoff $t.mp.max_perc $t.us.win_size \
			$t.cf.min_len $t.cf.win_size $t.cf.max_dash" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Check Assembly}" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.cutoff $t.mp $t.mp.max_perc $t.us $t.us.win_size \
	$t.cf $t.cf.min_len $t.cf.win_size $t.cf.max_dash $t.ok \
	-side top -fill x
}

proc CheckAssembly2 {io t lorf id cutoff maxperc u_winsize
		     minlen c_winsize maxdash} {
    global gap_defs
    global NGRec

    set args ""

    if {[set mperc   [scalebox_get $maxperc]] == ""} {bell; return}
    if {[set cutoff [yes_no_get   $cutoff]]} {
        if {[set minl    [entrybox_get $minlen]]    == ""} {bell; return}
        if {[set wsize   [scalebox_get $c_winsize]] == ""} {bell; return}
        if {[set mdash   [scalebox_get $maxdash]]   == ""} {bell; return}
	set args "-min_len $minl -win_size $wsize -max_dashes $mdash"
    } else {
        if {[set wsize   [scalebox_get $u_winsize]] == ""} {bell; return}
	set args "-win_size $wsize"
    }

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

    eval check_assembly \
	-io $io \
	-max_pmismatch $mperc \
	-cutoff $cutoff \
	-contigs {$list} \
	$args
    ClearBusy
}
