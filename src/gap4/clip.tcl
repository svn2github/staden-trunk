#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc QualityClip {io} {
    global gap_defs

    set l [keylget gap_defs QUALITY_CLIP]
    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Quality clip"

    contig_id $t.id -io $io

    lorf_in $t.infile [keylget l INFILE] \
	"{contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state normal}
	" -bd 2 -relief groove

    scalebox $t.minqual \
	-orient horizontal -width 5 \
	-title   "[keylget l MINQUAL.NAME]" \
	-from    "[keylget l MINQUAL.MIN]" \
	-to      "[keylget l MINQUAL.MAX]" \
	-default "[keylget l MINQUAL.VALUE]"

    okcancelhelp $t.ok \
	-ok_command "QualityClip2 $io $t $t.infile $t.id $t.minqual" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Clip-Quality}" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.minqual $t.ok \
	-side top -fill x
}

proc QualityClip2 {io t lorf id minqual} {
    if {[set minq [scalebox_get $minqual]]  == ""} {bell; return}

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

    destroy $t

    if {![quit_displays $io "quality clip"]} {
	return
    }

    SetBusy
    quality_clip \
	-io $io \
	-contigs $list \
	-quality $minq
    catch {ContigInitReg $io}

    ClearBusy
}


proc DifferenceClip {io} {
    global gap_defs

    set l [keylget gap_defs DIFFERENCE_CLIP]
    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Difference clip"

    contig_id $t.id -io $io

    lorf_in $t.infile [keylget l INFILE] \
	"{contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state normal}
	" -bd 2 -relief groove

    yes_no $t.tag \
	-title [keylget l TAG.TITLE] \
	-relief groove \
	-bd 2 \
	-orient horizontal \
	-default [keylget l TAG.DEFAULT]

    okcancelhelp $t.ok \
	-ok_command "DifferenceClip2 $io $t $t.infile $t.id $t.tag" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Clip-Difference}" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.tag $t.ok -side top -fill x
}

proc DifferenceClip2 {io t lorf id tag} {
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

    destroy $t

    if {![quit_displays $io "difference clip"]} {
	return
    }

    SetBusy
    difference_clip \
	-io $io \
	-contigs $list \
	-tag [yes_no_get $tag]
    catch {ContigInitReg $io}

    ClearBusy
}

proc NClip {io} {
    global gap_defs

    set l [keylget gap_defs N_CLIP]
    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "N-base clip"

    contig_id $t.id -io $io

    lorf_in $t.infile [keylget l INFILE] \
	"{contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state normal}
	" -bd 2 -relief groove

    okcancelhelp $t.ok \
	-ok_command "NClip2 $io $t $t.infile $t.id" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Clip-NBases}" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.ok -side top -fill x
}

proc NClip2 {io t lorf id} {
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

    destroy $t

    if {![quit_displays $io "N-base clip"]} {
	return
    }

    SetBusy
    N_clip \
	-io $io \
	-contigs $list
    catch {ContigInitReg $io}

    ClearBusy
}

proc QClipEnds {io} {
    global gap_defs

    set l [keylget gap_defs Q_CLIP_ENDS]
    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Quality clip ends"

    contig_id $t.id -io $io

    lorf_in $t.infile [keylget l INFILE] \
	"{contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state normal}
	" -bd 2 -relief groove

    okcancelhelp $t.ok \
	-ok_command "QClipEnds2 $io $t $t.infile $t.id" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Clip-QClipEnds}" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.ok -side top -fill x
}

proc QClipEnds2 {io t lorf id} {
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

    destroy $t

    if {![quit_displays $io "N-base clip"]} {
	return
    }

    SetBusy
    quality_clip_ends \
	-io $io \
	-contigs $list
    catch {ContigInitReg $io}

    ClearBusy
}


