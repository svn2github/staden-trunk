#-----------------------------------------------------------------------------
# Export Sequences

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

    #--- output filename
    getFname $f.outfile "Output filename" save {} [$io name].sam

    #--- formats
    radiolist $f.format \
	-title "Select format" \
	-default 1 \
	-orient horizontal \
	-buttons [list \
	     [list sam   -command "ExportSequences_format $io $f"] \
	     [list ace   -command "ExportSequences_format $io $f"] \
	     [list caf   -command "ExportSequences_format $io $f"] \
	     [list baf   -command "ExportSequences_format $io $f"] \
	     [list fastq -command "ExportSequences_format $io $f"] \
	     [list fasta -command "ExportSequences_format $io $f"]]
	

    #--- Output options
    checkbutton $f.fixmates \
	-text "Fix mate-pair information (SAM only)" \
	-variable $f.FixMates \
	-anchor w
    global $f.FixMates
    set $f.FixMates 0

    #--- OK/cancel/help
    okcancelhelp $f.ok \
	-ok_command "ExportSequences2 $io $f" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 ExportSequences" \
	-bd 2 \
	-relief groove

    #--- Packing
    pack $f.infile $f.id $f.format $f.outfile $f.fixmates $f.ok \
	-side top -fill both
}

# Callback for when the output format is changed
proc ExportSequences_format {io f} {
    set format [lindex "x sam ace caf baf fastq fasta" [radiolist_get $f.format]]
    set entry [entrybox_path $f.outfile.entry]
    set fn [$entry get]
    if {[regexp {(.*)\.(sam|fna|fa|fasta|fastq|baf|ace|caf)$} $fn _ pfix sfix]} {
	set fn $pfix.$format
    } else {
	set fn $fn.$format
    }
    $entry delete 0 end
    $entry insert 0 $fn
}

proc ExportSequences2 {io f} {
    global $f.FixMates

    set format [lindex "x sam ace caf baf fastq fasta" [radiolist_get $f.format]]
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

    export_contigs -io $io -contigs $list -format $format -outfile $fn \
	-fixmates [set $f.FixMates]
    destroy $f
}


#-----------------------------------------------------------------------------
# Export Tags


proc ExportTags {io} {
    global gap5_defs

    set f .export_tags
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Export Tags"

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
	-default 2 \
	-orient horizontal \
	-buttons "baf gff"

    #--- Output options
    checkbutton $f.unpadded \
	-text "Unpadded coordinates" \
	-variable $f.Unpadded \
	-anchor w
    global $f.Unpadded
    set $f.Unpadded 1

    checkbutton $f.consensus \
	-text "Map sequence tags to consensus" \
	-variable $f.Consensus \
	-anchor w
    global $f.Consensus
    set $f.Consensus 1

    #--- output filename
    getFname $f.outfile "Output filename" save {} tags.gff

    #--- OK/cancel/help
    okcancelhelp $f.ok \
	-ok_command "ExportTags2 $io $f" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 ExportTags" \
	-bd 2 \
	-relief groove

    #--- Packing
    pack $f.infile $f.id $f.format $f.unpadded $f.consensus \
	$f.outfile $f.ok -anchor w -side top -fill both
}

proc ExportTags2 {io f} {
    global $f.Unpadded
    global $f.Consensus

    set format [lindex "x baf gff" [radiolist_get $f.format]]
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

    export_tags \
	-io $io \
	-contigs $list \
	-format $format \
	-outfile $fn \
	-unpadded  [set $f.Unpadded] \
	-consensus [set $f.Consensus]
    destroy $f
}