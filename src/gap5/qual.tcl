#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

proc ListConfidence {io} {
    global gap5_defs

    set t [keylget gap5_defs LIST_CONFIDENCE.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "List confidence"

    contig_id $t.id -io $io

    lorf_in $t.infile [keylget gap5_defs LIST_CONFIDENCE.INFILE] \
	"{contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state normal}
	" -bd 2 -relief groove

    set l [keylget gap5_defs LIST_CONFIDENCE.SUMMARY]
    yes_no $t.summary \
	-title "[keylget l NAME]" \
	-orient horizontal \
	-default "[keylget l VALUE]"

    okcancelhelp $t.ok \
	-ok_command "ListConfidence2 $io $t $t.infile $t.id $t.summary" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 {Con-Evaluation}" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.summary $t.ok -side top -fill x
}

proc ListConfidence2 {io t lorf id summary} {
    if {[lorf_in_get $lorf] == 4} {
        if {[set contign [contig_id_gel $id]] == ""} {bell; return}
       	if {[set lreg  [contig_id_lreg $id]] == ""} {bell; return}
	if {[set rreg  [contig_id_rreg $id]] == ""} {bell; return}
	set list "{$contign $lreg $rreg}"
	SetContigGlobals $io $contign
    } elseif {[lorf_in_get $lorf] == 3} {
	set list [CreateAllContigList $io]
    } else {
	if {[set list [lorf_get_list $lorf]] == ""} {bell; return}
    }
    set summary [yes_no_get $summary]

    destroy $t

    SetBusy
    list_confidence \
	-io $io \
	-contigs $list \
	-summary $summary
    ClearBusy
}

proc ListBaseConfidence {io} {
    global gap5_defs

    set t [keylget gap5_defs LIST_CONFIDENCE.WIN]2
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "List base confidence"

    contig_id $t.id -io $io -range 0

    lorf_in $t.infile [keylget gap5_defs LIST_CONFIDENCE.INFILE] \
	"{contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state normal}
	" -bd 2 -relief groove

    okcancelhelp $t.ok \
	-ok_command "ListBaseConfidence2 $io $t $t.infile $t.id" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 {Con-ListBaseConf}" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.ok -side top -fill x
}


proc ListBaseConfidence2 {io t lorf id} {
    if {[lorf_in_get $lorf] == 4} {
        if {[set contign [contig_id_gel $id]] == ""} {bell; return}
	set list "{$contign}"
	SetContigGlobals $io $contign
    } elseif {[lorf_in_get $lorf] == 3} {
	set list [CreateAllContigList $io]
    } else {
	if {[set list [lorf_get_list $lorf]] == ""} {bell; return}
    }
    destroy $t

    SetBusy
    list_base_confidence \
	-io $io \
	-contigs $list
    ClearBusy

}

