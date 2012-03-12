array set _ioh_ {}
proc ioh {io {val {}}} {
    if {$val != {}} {
	set _ioh_($io) 1
    }
    return $_ioh_($io)
}

proc check_database {args} {
    #puts [info level [info level]]
    return 0
}

proc io_read_database {io} {
    keylset db num_contigs 1
    keylset db num_readings 1
    keylset db max_gel_len 32768
    keylset db actual_db_size 1
    keylset db Nannotations 1
    keylset db Ntemplates 1
    keylset db Nclones 1
    keylset db Nvectors 1
    return $db
}

proc io_read_contig {io cnum} {
    set c [$io get_contig $cnum]
    set start [$c get_start]
    set end   [$c get_end]
    set left  [$c seqs_in_range $start [expr {$start+1}]]
    set right [$c seqs_in_range [expr {$end-1}] $end]
    keylset contig \
	length [expr {$end-$start+1}] \
	left [lindex [lindex $left 0] 2] \
	right [lindex [lindex $right 0] 2]
    
    return $contig
}

proc cname2crec {io name} {
    # Check B-Tree indices if present
    if {[set crec [db_info get_contig_num $io $name]] != -1} {
	return $crec
    }

    # =contig_rec => rec
    if {[string match "=*" $name]} {
	return [string range $name 1 end]
    }

    # #sequence_rec => seq2contig conversion
    if {[string match "\#*" $name]} {
	set rec [string range $name 1 end]
	if {[$io rec_exists 18 $rec]} {
	    set s [$io get_seq $rec]
	    set cr [$s get_contig]
	    $s delete
	    return $cr
	}
    }

    # else "contig_name"
    set nc [$io num_contigs]
    for {set i 0} {$i < $nc} {incr i} {
	set crec [$io contig_order $i]
	set c [$io get_contig $crec]
	if {[string match $name [$c get_name]]} {
	    $c delete
	    return $crec
	}
	$c delete
    }

    return -1
}