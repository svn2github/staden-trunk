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
	    vmessage_tagged $m CMD_STDOUT
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
    return 1

    foreach fn [glob $prefix.*] {
	file delete $fn
    }

    return 1
}

# Creates a fastq file from a fofn
proc generate_fastq {f out} {
    set fd [open $f r]
    set fq [open $out w]
    while {[gets $fd line] != -1} {
	foreach {s q} [read_seq_trace $line] break
	puts $fq "@$line\n$s\n+\n$q"
    }
    close $fd
    close $fq
}

#-----------------------------------------------------------------------------
# BWA using aln and samse/sampe mode.

proc MapReads_bwa_aln {io} {
    global gap5_defs

    set w .map_reads
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Mapped assembly - bwa aln"

    # Sequences to map against
    contig_id $w.id -io $io 
    lorf_in $w.contigs [keylget gap5_defs MAP_READS.INFILE] \
	"{contig_id_configure $w.id -state disabled} \
	 {contig_id_configure $w.id -state disabled}\
	 {contig_id_configure $w.id -state disabled}\
	 {contig_id_configure $w.id -state normal}" -bd 2 -relief groove
	

    # fasta input file name, or something to generate a sam file from
    radiolist $w.format \
	-title "Input readings from" \
	-orient horizontal \
	-buttons {fofn fasta fastq} \
	-default 3

    xentry $w.fwd \
	-label "Forward read file" \
	-checkcommand "check_fileinput"
    xentry $w.rev \
	-label "Reverse read file (optional)" \
	-checkcommand "check_fileinput 1"

    frame $w.sep1 -bd 2 -relief groove -height 2

    radiolist $w.out_format \
	-title "Output failed readings as" \
	-orient horizontal \
	-buttons {fasta fastq} \
	-default 2

    xentry $w.out_fn \
	-label "Sequence output filename" \
	-checkcommand "check_fileoutput" \
	-default [keylget gap5_defs MAP_READS.OUTPUT_FN]

    frame $w.sep2 -bd 2 -relief groove -height 2

    entrybox $w.aln_opt   -title "'bam aln' options"
    entrybox $w.sampe_opt -title "'bam sampe' options"

    xyn $w.index_names \
	-label "Index sequence names" \
	-orient horiz \
	-default [keylget gap5_defs MAP_READS.INDEX_SEQUENCE_NAMES]

    #OK and Cancel buttons
    okcancelhelp $w.ok_cancel \
	    -ok_command "MapReads_bwa_aln2 $io $w" \
	    -cancel_command "destroy $w" \
	    -help_command "show_help gap5 {Assembly-Map-bwa-aln}" \
	    -bd 2 \
	    -relief groove

    pack $w.contigs $w.id $w.format $w.fwd $w.rev $w.sep1 \
	$w.out_format $w.out_fn $w.sep2 \
	$w.index_names $w.aln_opt $w.sampe_opt $w.ok_cancel \
	-side top -fill both
}

proc MapReads_bwa_aln2 {io w} {
    set fq_tmp ""

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

    set no_tree [$w.index_names get]

    set prefix [tmpnam]

    # and its input sequences
    if {[radiolist_get $w.format] == 1} {
	set fwd $prefix.fq_in
	generate_fastq [$w.fwd get] $fwd
    } else {
	set fwd [$w.fwd get]
	if {$fwd == ""} { 
	    bell
	    return
	}
    }
    set rev [$w.rev get]

    # Output name for failure
    set out_fn [$w.out_fn get]
    set out_fmt [radiolist_get $w.out_format]

    # options for various stages
    set aln_opt   [entrybox_get $w.aln_opt]
    set sampe_opt [entrybox_get $w.sampe_opt]

    destroy $w
    SetBusy

    # Do the work
    vfuncheader "Map Reads - bwa"

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
    log_call import_reads \
	-io $io \
	-append 1 \
	-merge_contigs 1 \
	-repad 1 \
	-file $prefix.srt.bam \
	-format bam \
	-index_names $no_tree

    # Process failures from SAM file. Convert back to fasta or fastq.
    # $w.out_format: 1=fasta, 2=fastq
    if {$out_fn != ""} {
	sam_failures $prefix.sam $out_fn [expr {$out_fmt-1}]
    }

    vmessage "Flushing"
    $io flush

    vmessage "(done)"

    MapReads_tidyup $prefix
    ClearBusy

    return 0;
}


#-----------------------------------------------------------------------------
# BWA using dbwtsw option (now bwasw)

proc MapReads_bwa_bwasw {io} {
    global gap5_defs

    set w .map_reads2
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Mapped assembly - bwa bwasw"

    # Sequences to map against
    contig_id $w.id -io $io 
    lorf_in $w.contigs [keylget gap5_defs MAP_READS.INFILE] \
	"{contig_id_configure $w.id -state disabled} \
	 {contig_id_configure $w.id -state disabled}\
	 {contig_id_configure $w.id -state disabled}\
	 {contig_id_configure $w.id -state normal}" -bd 2 -relief groove
	
    # fasta input file name, or something to generate a sam file from
    radiolist $w.format \
	-title "Input readings from" \
	-orient horizontal \
	-buttons {fofn fasta fastq} \
	-default 3

    xentry $w.fwd \
	-label "Sequence file or fofn" \
	-checkcommand "check_fileinput"

    frame $w.sep1 -bd 2 -relief groove -height 2

    radiolist $w.out_format \
	-title "Output failed readings as" \
	-orient horizontal \
	-buttons {fasta fastq} \
	-default 2

    xentry $w.out_fn \
	-label "Sequence output filename" \
	-checkcommand "check_fileoutput" \
	-default [keylget gap5_defs MAP_READS.OUTPUT_FN]

    frame $w.sep2 -bd 2 -relief groove -height 2

    #save failures as file or list
#    lorf_out $w.fails [keylget gap5_defs AUTO_ASSEMBLE.FAILS] \
	"" -bd 2 -relief groove

    entrybox $w.bwasw_opt   -title "'bam bwasw' options"

    xyn $w.index_names \
	-label "Index sequence names" \
	-orient horiz \
	-default [keylget gap5_defs MAP_READS.INDEX_SEQUENCE_NAMES]

    #OK and Cancel buttons
    okcancelhelp $w.ok_cancel \
	    -ok_command "MapReads_bwa_bwasw2 $io $w" \
	    -cancel_command "destroy $w" \
	    -help_command "show_help gap5 {Assembly-Map-bwa-bwasw}" \
	    -bd 2 \
	    -relief groove

    pack $w.contigs $w.id $w.format $w.fwd $w.sep1 \
	$w.out_format $w.out_fn $w.sep2 $w.index_names \
	$w.bwasw_opt $w.ok_cancel -side top -fill both
}

proc MapReads_bwa_bwasw2 {io w} {
    set prefix [tmpnam]

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

    set no_tree [$w.index_names get]

    # and its input sequences
    if {[radiolist_get $w.format] == 1} {
	set fwd $prefix.fq_in
	generate_fastq [$w.fwd get] $fwd
    } else {
	set fwd [$w.fwd get]
	if {$fwd == ""} { 
	    bell
	    return
	}
    }

    # Output name for failure
    set out_fn [$w.out_fn get]
    set out_fmt [radiolist_get $w.out_format]

    # options for various stages
    set bwasw_opt   [entrybox_get $w.bwasw_opt]

    destroy $w
    SetBusy

    # Do the work
    vfuncheader "Map Reads - bwa"

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
    if {[MapReads_run "bwa bwasw $bwasw_opt $prefix.cons $fwd > $prefix.sam"]} {
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
    log_call import_reads \
	-io $io \
	-append 1 \
	-merge_contigs 1 \
	-repad 1 \
	-file $prefix.srt.bam \
	-format bam \
	-index_names $no_tree

    # Process failures from SAM file. Convert back to fasta or fastq.
    # $w.out_format: 1=fasta, 2=fastq
    if {$out_fn != ""} {
	sam_failures $prefix.sam $out_fn [expr {$out_fmt-1}]
    }

    vmessage "Flushing"
    $io flush

    MapReads_tidyup $prefix
    ClearBusy

    return 0;
}

# Reads the sam file extracting unmapped reads, outputting them to $out in
# either fastq ($fastq=1) or fasta format.
#
# Returns the number of unmapped reads.
proc sam_failures {sam out {fastq 1}} {
    set unmapped 0
    set in_fd  [open $sam]
    set out_fd [open $out w]

    while {[gets $in_fd line] != -1} {
	if {[string match "@*" $line]} continue
	if {[lindex $line 1] & 4} {
	    # Unmapped
	    if {$fastq} {
		set s [lindex $line 9]
		set q [lindex $line 10]
		if {$q == "*"} {
		    set q [string repeat "!" [string length $s]]
		}
		puts $out_fd "@[lindex $line 0]\n$s\n+\n$q"
	    } else {
		puts $out_fd ">[lindex $line 0]\n[lindex $line 9]"
	    }
	    incr unmapped
	}
    }
    close $in_fd
    close $out_fd;

    return $unmapped
}
