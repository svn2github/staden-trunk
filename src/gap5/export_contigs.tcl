proc ExportSequences {io} {
    global gap5_defs

    set f .export_sequences
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Export Sequences"

    #--- contig identifier widget
    contig_id $f.id -io $io

    lorf_in $f.infile [keylget gap5_defs CONSENSUS.INFILE] \
        "{contig_id_configure $f.id -state disabled}
             {contig_id_configure $f.id -state disabled}
         {contig_id_configure $f.id -state disabled}
         {contig_id_configure $f.id -state normal}
        " -bd 2 -relief groove

    #--- formats
    radiolist $f.format \
	-title "Select format" \
	-default 1 \
	-orient horizontal \
	-buttons "sam ace baf fastq fasta"

    #--- output filename
    getFname $f.outfile "Output filename" save {} aln.out

    #--- OK/cancel/help
    okcancelhelp $f.ok \
	-ok_command "ExportSequences2 $io $f" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 ExportSequences" \
	-bd 2 \
	-relief groove

    #--- Packing
    pack $f.infile $f.id $f.format $f.outfile $f.ok -side top -fill both
}

proc ExportSequences2 {io f} {
    set format [lindex "x sam ace baf fastq fasta" [radiolist_get $f.format]]
    if {[lorf_in_get $f.infile] == 4} {
	set gel_name [contig_id_gel $f.id]
	set lreg [contig_id_lreg $f.id]
	set rreg [contig_id_rreg $f.id]
	set list [list [list $gel_name $lreg $rreg]]
    } elseif {[lorf_in_get $f.infile] == 3} {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $f.infile]
    }

    set fn [getFname_in_name $f.outfile]

    if {$list == ""} {
	raise $f
	return
    }

    export_contigs -io $io -contigs $list -format $format -outfile $fn
    destroy $f
}