#-----------------------------------------------------------------------------
# Shuffle Pads

proc ShufflePads {io} {
    global gap5_defs

    set t [keylget gap5_defs SHUFFLE_PADS.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Shuffle Pads"
    contig_id $t.id -io $io -range 1

    lorf_in $t.infile [keylget gap5_defs SHUFFLE_PADS.INFILE] \
	"{contig_id_configure $t.id -state disabled} \
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state normal}" -bd 2 -relief groove

    xentry $t.band_size \
	-label "Band size" \
	-default [keylget gap5_defs SHUFFLE_PADS.BAND_SIZE]

    okcancelhelp $t.ok \
	-ok_command "ShufflePads2 $io $t $t.infile $t.id $t.band_size" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 {Tidying up alignments}" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.band_size $t.ok -side top -fill x
}

;proc ShufflePads2 {io t infile id band_size} {
    if {[lorf_in_get $infile] == 4} {
	set list [list [contig_id_gel $id]]
	set lreg [contig_id_lreg $id]
	set rreg [contig_id_rreg $id]
	set list "{$list $lreg $rreg}"
    } elseif {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $infile]
    }
    if {[set band_size [$band_size get]] < 0} {
	bell
	return
    }

    if {![quit_displays -io $io -msg "shuffle_pads"]} {
	# Someone's too busy to shutdown?
	return
    }

    destroy $t

    SetBusy
    shuffle_pads -io $io -contigs $list -band $band_size
    ClearBusy

    ContigInitReg $io
    InitContigGlobals $io
}

#-----------------------------------------------------------------------------
# Remove Pad Columns

proc RemovePadColumns {io} {
    global gap5_defs

    set t [keylget gap5_defs REMOVE_PAD_COLUMNS.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Remove Pad Columns"
    contig_id $t.id -io $io -range 0

    lorf_in $t.infile [keylget gap5_defs REMOVE_PAD_COLUMNS.INFILE] \
	"{contig_id_configure $t.id -state disabled} \
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state normal}" -bd 2 -relief groove

    xentry $t.percent_pad \
	-label "Percentage pad needed" \
	-default [keylget gap5_defs REMOVE_PAD_COLUMNS.PERCENT_PAD]

    okcancelhelp $t.ok \
	-ok_command "RemovePadColumns2 $io $t $t.infile $t.id $t.percent_pad" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap5 {Remove Pad Columns}" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.percent_pad $t.ok -side top -fill x
}

;proc RemovePadColumns2 {io t infile id percent_pad} {
    if {[lorf_in_get $infile] == 4} {
	set list [list [contig_id_gel $id]]
    } elseif {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $infile]
    }
    if {[set percent_pad [$percent_pad get]] < 0} {
	bell
	return
    }

    if {![quit_displays -io $io -msg "shuffle_pads"]} {
	# Someone's too busy to shutdown?
	return
    }

    destroy $t

    SetBusy
    remove_pad_columns -io $io -contigs $list -percent_pad $percent_pad
    ClearBusy

    ContigInitReg $io
    InitContigGlobals $io
}
