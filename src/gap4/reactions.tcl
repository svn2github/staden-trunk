#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc LongGelsDialog {io f} {
    global gap_defs

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Suggest long readings"

    ###########################################################################
    contig_id $f.id \
	    -state normal\
	    -io $io

    lorf_in $f.infile [keylget gap_defs LONGGELS.INFILE] \
	    "{contig_id_configure $f.id -state disabled} \
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state normal}" -bd 2 -relief groove

    ###########################################################################
    keylset gl GLEN [keylget gap_defs LONGGELS.GLEN]
    entrybox $f.len \
	    -title [keylget gl GLEN.NAME]\
	    -default [keylget gl GLEN.VALUE] \
	    -width 5 \
	    -type "CheckIntMin 1"
	    
    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "OK_Pressed_lg $io $f $f.id $f.len $f.infile"\
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Suggest Long}" \
	    -bd 2 \
	    -relief groove

    ###########################################################################
    pack $f.infile -side top -fill both
    pack $f.id -fill x
    pack $f.len -fill x
    pack $f.ok_cancel -fill x
}

proc OK_Pressed_lg {io f id len infile} {
    if {[lorf_in_get $infile] == 4} {
	set gel_name [contig_id_gel $id]
	set lreg [contig_id_lreg $id]
	set rreg [contig_id_rreg $id]
	
	SetContigGlobals $io $gel_name $lreg $rreg
	set list "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $infile] == 3 } {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $infile]
    }

    if {[set avg_len [entrybox_get $len]] == ""} {return}

    SetBusy
    find_long_gels \
	-io $io \
	-contigs $list \
	-avg_len $avg_len
    ClearBusy

    destroy $f
}
proc TaqTerminatorDialog {io f} {
    global gap_defs

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Compressions and stops"

    ###########################################################################
    contig_id $f.id \
	    -state normal\
	    -io $io

    lorf_in $f.infile [keylget gap_defs TTERM.INFILE] \
	    "{contig_id_configure $f.id -state disabled} \
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state normal}" -bd 2 -relief groove

    ###########################################################################
    keylset tt TLEN [keylget gap_defs TTERM.TLEN]
    entrybox $f.len \
	    -title [keylget tt TLEN.NAME]\
	    -default [keylget tt TLEN.VALUE] \
	    -width 5 \
	    -type "CheckIntMin 1"
	    
    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "OK_Pressed_tt $io $f $f.id $f.len $f.infile"\
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Compressions}" \
	    -bd 2 \
	    -relief groove

    ###########################################################################
    pack $f.infile -side top -fill both
    pack $f.id -fill x
    pack $f.len -fill x
    pack $f.ok_cancel -fill x
}

proc OK_Pressed_tt {io f id len infile} {
    if {[lorf_in_get $infile] == 4} {
	set gel_name [contig_id_gel $id]
	set lreg [contig_id_lreg $id]
	set rreg [contig_id_rreg $id]
	
	SetContigGlobals $io $gel_name $lreg $rreg
	set list "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $infile] == 3 } {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $infile]
    }
    if {[set avg_len [entrybox_get $len]] == ""} {return}

    SetBusy
    find_taq_terminator \
	-io $io \
	-contigs $list \
	-avg_len $avg_len
    ClearBusy

    destroy $f
}

