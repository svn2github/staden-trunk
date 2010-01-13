# Run a command and display the output.
# Returns 0 for success, or error code for failure
# code on failure
;proc MapReads_run {cmd} {
    vmessage ""
    vmessage_tagged $cmd CMD_NAME
    update idletasks

    # If the command has it's own redirect for stdout, then assume stderr
    # is the normal diagnostic output and only report that. Otherwise assume
    # stdout is normal output and stderr is erroneous.
    if {[regexp {[^2>]>} $cmd]} {
	catch {eval exec $cmd} err
	if {$err != ""} {
	    if {$::errorCode != "NONE"} {
		vmessage_tagged $err CMD_STDERR
	    } else {
		vmessage_tagged $err CMD_STDOUT
	    }
	}
    } else {
	set o_fn [tmpnam]
	catch {eval exec $cmd > $o_fn} err
	set fd [open $o_fn r]
	set m [read $fd]
	close $fd
	file delete $o_fn
	if {$m != ""} {
	    vmessage_tagged [read $fd] CMD_STDOUT
	}
	if {$err != ""} {
	    vmessage_tagged $err CMD_STDERR
	}
    }
	update idletasks

    if {$::errorCode != "NONE"} {
	return 1
    }
    
    return 0
}

proc MapReads_tidyup {prefix} {
    foreach fn [glob $prefix.*] {
	file delete $fn
    }

    return 1
}

#-----------------------------------------------------------------------------
# BWA using aln and samse/sampe mode.

proc MapReads_bwa_aln {io} {
    global gap5_defs

    set w .map_reads
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Mapped assembly"

    # Sequences to map against
    contig_id $w.id -io $io 
    lorf_in $w.contigs [keylget gap5_defs MAP_READS.INFILE] \
	"{contig_id_configure $w.id -state disabled} \
	 {contig_id_configure $w.id -state disabled}\
	 {contig_id_configure $w.id -state disabled}\
	 {contig_id_configure $w.id -state normal}" -bd 2 -relief groove
	

    # fasta input file name, or something to generate a sam file from
    getFname $w.fwd "Forward read fastq file" load
    getFname $w.rev "Reverse read fastq (optional)" load_optional

    #save failures as file or list
#    lorf_out $w.fails [keylget gap5_defs AUTO_ASSEMBLE.FAILS] \
	"" -bd 2 -relief groove

    entrybox $w.aln_opt   -title "'bam aln' options"
    entrybox $w.sampe_opt -title "'bam sampe' options"

    #OK and Cancel buttons
    okcancelhelp $w.ok_cancel \
	    -ok_command "MapReads_bwa_aln2 $io $w" \
	    -cancel_command "destroy $w" \
	    -help_command "show_help gap5 {Assembly-Map}" \
	    -bd 2 \
	    -relief groove

    pack $w.contigs $w.id $w.fwd $w.rev $w.aln_opt $w.sampe_opt $w.ok_cancel \
	-side top -fill both
}

proc MapReads_bwa_aln2 {io w} {
    # Contigs to map against
    if {[lorf_in_get $w.contigs] == 4} {
	set gel_name [contig_id_gel $w.id]
	set lreg [contig_id_lreg $w.id]
	set rreg [contig_id_rreg $w.id]
	
	SetContigGlobals $io $gel_name $lreg $rreg
	set contigs "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $w.contigs] == 3 } {
	set contigs [CreateAllContigList $io]
    } else {
	set contigs [lorf_get_list $w.contigs]
    }

    # and its input sequences
    set fwd [entrybox_get $w.fwd.entry]
    if {$fwd == ""} { 
	bell
	return
    }
    set rev [entrybox_get $w.rev.entry]

    # options for various stages
    set aln_opt   [entrybox_get $w.aln_opt]
    set sampe_opt [entrybox_get $w.sampe_opt]

    destroy $w
    SetBusy

    # Do the work
    vfuncheader "Map Reads - bwa"
    set prefix [tmpnam]

    # Calculate and save the consensus as fastq
    get_consensus -io $io \
	-contigs $contigs \
	-format 1 \
	-outfile $prefix.cons \
	-strip_pads 1

    # Build a bwa index from it
    if {[MapReads_run "bwa index $prefix.cons"]} {
	ClearBusy
	return [MapReads_tidyup $prefix]
    }

    # Align our input data.
    if {[MapReads_run "bwa aln $aln_opt $prefix.cons $fwd > $prefix.f.sai"]} {
	ClearBusy
	return [MapReads_tidyup $prefix]
    }
    if {$rev == ""} {
	if {[MapReads_run "bwa samse $prefix.cons $prefix.f.sai $fwd > $prefix.sam"]} {
	    ClearBusy
	    return [MapReads_tidyup $prefix]
	}
    } else {
	if {[MapReads_run "bwa aln $aln_opt $prefix.cons $rev > $prefix.r.sai"]} {
	    ClearBusy
	    return [MapReads_tidyup $prefix]
	}
	if {[MapReads_run "bwa sampe $sampe_opt $prefix.cons $prefix.f.sai $prefix.r.sai $fwd $rev > $prefix.sam"]} {
	    ClearBusy
	    return [MapReads_tidyup $prefix]
	}
	ClearBusy
	return [MapReads_tidyup $prefix]
    }

    # Convert to bam and sort
    if {[MapReads_run "samtools view -b -S $prefix.sam > $prefix.bam"]} {
	ClearBusy
	return [MapReads_tidyup $prefix]
    }
    if {[MapReads_run "samtools sort $prefix.bam $prefix.srt"]} {
	ClearBusy
	return [MapReads_tidyup $prefix]
    }
    
    # Import
    vmessage "Importing reads..."
    import_reads \
	-io $io \
	-append 1 \
	-merge_contigs 1 \
	-repad 1 \
	-file $prefix.srt.bam \
	-format bam

    vmessage "Flushing"
    $io flush

    vmessage "(done)"

    MapReads_tidyup $prefix
    ClearBusy

    return 0;
}


#-----------------------------------------------------------------------------
# BWA using dbwtsw option

proc MapReads_bwa_dbwtsw {io} {
    global gap5_defs

    set w .map_reads
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Mapped assembly"

    # Sequences to map against
    contig_id $w.id -io $io 
    lorf_in $w.contigs [keylget gap5_defs MAP_READS.INFILE] \
	"{contig_id_configure $w.id -state disabled} \
	 {contig_id_configure $w.id -state disabled}\
	 {contig_id_configure $w.id -state disabled}\
	 {contig_id_configure $w.id -state normal}" -bd 2 -relief groove
	

    # fasta input file name, or something to generate a sam file from
    getFname $w.fwd "Input fastq file" load

    #save failures as file or list
#    lorf_out $w.fails [keylget gap5_defs AUTO_ASSEMBLE.FAILS] \
	"" -bd 2 -relief groove

    entrybox $w.dbwtsw_opt   -title "'bam dbwtsw' options"

    #OK and Cancel buttons
    okcancelhelp $w.ok_cancel \
	    -ok_command "MapReads_bwa_dbwtsw2 $io $w" \
	    -cancel_command "destroy $w" \
	    -help_command "show_help gap5 {Assembly-Map}" \
	    -bd 2 \
	    -relief groove

    pack $w.contigs $w.id $w.fwd $w.dbwtsw_opt $w.ok_cancel \
	-side top -fill both
}

proc MapReads_bwa_dbwtsw2 {io w} {
    # Contigs to map against
    if {[lorf_in_get $w.contigs] == 4} {
	set gel_name [contig_id_gel $w.id]
	set lreg [contig_id_lreg $w.id]
	set rreg [contig_id_rreg $w.id]
	
	SetContigGlobals $io $gel_name $lreg $rreg
	set contigs "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $w.contigs] == 3 } {
	set contigs [CreateAllContigList $io]
    } else {
	set contigs [lorf_get_list $w.contigs]
    }

    # and its input sequences
    set fwd [entrybox_get $w.fwd.entry]
    if {$fwd == ""} { 
	bell
	return
    }

    # options for various stages
    set dbwtsw_opt   [entrybox_get $w.dbwtsw_opt]

    destroy $w
    SetBusy

    # Do the work
    vfuncheader "Map Reads - bwa"
    set prefix [tmpnam]

    # Calculate and save the consensus as fastq
    get_consensus -io $io \
	-contigs $contigs \
	-format 1 \
	-outfile $prefix.cons \
	-strip_pads 1

    # Build a bwa index from it
    if {[MapReads_run "bwa index $prefix.cons"]} {
	ClearBusy
	return [MapReads_tidyup $prefix]
    }

    # Align our input data.
    if {[MapReads_run "bwa dbwtsw $dbwtsw_opt $prefix.cons $fwd > $prefix.sam"]} {
	ClearBusy
	return [MapReads_tidyup $prefix]
    }

    # Convert to bam and sort
    if {[MapReads_run "samtools view -b -S $prefix.sam > $prefix.bam"]} {
	ClearBusy
	return [MapReads_tidyup $prefix]
    }
    if {[MapReads_run "samtools sort $prefix.bam $prefix.srt"]} {
	ClearBusy
	return [MapReads_tidyup $prefix]
    }
    
    # Import
    vmessage "Importing reads..."
    import_reads \
	-io $io \
	-append 1 \
	-merge_contigs 1 \
	-repad 1 \
	-file $prefix.srt.bam \
	-format bam

    vmessage "Flushing"
    $io flush

    vmessage "(done)"

    MapReads_tidyup $prefix
    ClearBusy

    return 0;
}
