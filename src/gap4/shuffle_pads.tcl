proc ShufflePads {io} {
    global gap_defs

    set t [keylget gap_defs SHUFFLE_PADS.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Shuffle Pads"
    contig_id $t.id -io $io -range 0

    lorf_in $t.infile [keylget gap_defs SHUFFLE_PADS.INFILE] \
	"{contig_id_configure $t.id -state disabled} \
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state disabled}\
	 {contig_id_configure $t.id -state normal}" -bd 2 -relief groove

    xentry $t.band_size \
	-label "Band size" \
	-default [keylget gap_defs SHUFFLE_PADS.BAND_SIZE]

    okcancelhelp $t.ok \
	-ok_command "ShufflePads2 $io $t $t.infile $t.id $t.band_size" \
	-cancel_command "destroy $t" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.band_size $t.ok -side top -fill x
}

;proc ShufflePads2 {io t infile id band_size} {
    if {[lorf_in_get $infile] == 4} {
	set list [list [contig_id_gel $id]]
    } elseif {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $infile]
    }
    if {[set band_size [$band_size get]] < 0} {
	bell
	return
    }

    if {![quit_displays $io "disassemble_readings"]} {
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