load_package gap
#load_package finish
load libprefinish.so

# -----------------------------------------------------------------------------
# Sanger rules

proc finishing_rules_sanger {bits} {
    # Problem 1: contig left end
    set p1	[expr {($bits & 0x0100)>0}]

    # Problem 2: contig right end
    set p2	[expr {($bits & 0x0200)>0}]

    # Problem 3: low sequence or low template coverage
    set p3	[expr {!($bits & 0x0004) || !($bits & 0x0008)}]

    # Problem 4: no BigDye terminator sequences
    set p4	[expr {!($bits & 0x0020)}]

    # Problem 5: low consensus confidence
    set p5	[expr {!($bits & 0x0010)}]

    # Problem 6: single stranded
    set p6	[expr {($bits & 0x0003) != 3}]

    return [expr {$p1 | ($p2<<1) | ($p3<<2) | ($p4<<3) | ($p5<<4) | ($p6<<5)}]
}

#
# Chemistry (from gap-dbstruct.h) - pick ONE of
# x/y => Primer/terminator
# 0       unknown         primer
# 1       unknown         terminator
# 2       ABI rhodamine   primer
# 3       ABI rhodamine   terminator
# 4       ABI dRhodamine  primer
# 5       ABI dRhodamine  terminator
# 6       ABI BigDye v2   primer
# 7       ABI BigDye v2   terminator
# 8       Energy transfer primer
# 9       Energy transfer terminator
# 10      Licor           primer
# 11      Licor           terminator
# 12      MegaBACE ET     primer
# 13      MegaBACE ET     terminator
# 14      ABI BigDye v1   primer
# 15      ABI BigDye v1   terminator
# 16      ABI BigDye v3   primer
# 17      ABI BigDye v3   terminator
#
# Strand - pick ONE of
# 0 Any
# 1 +
# 2 -
#
# Solution type - Bit pattern of one or more
# 0x01 None (skip)
# 0x02 Re-sequence
# 0x04 Vector primer walk
# 0x08 Long gel
# 0x10 PCR
# 0x20 Chromosomal primer walk
# 0x40 Reverse sequences
#

# All solutions are base around primer walking
# Nasty problems which need addressing first
proc find_solutions_sanger1 {base_bits problem_bits} {
    # Optimisation - no problems => no solutions required
    if {$problem_bits == 0} {return 0}

    set type 0;	  # None
    set chem 17;  # BDv3 term
    set strand 0; # Any

    # Problem 1: left contig end
    if {[expr {$problem_bits&0x01}]} {
	# Primer walk, bot strand.
	set type [expr {$type|0x04}]
	set strand 2
    }

    # Problem 2: right contig end
    if {[expr {$problem_bits&0x02}]} {
	# Primer walk, top strand.
	set type [expr {$type|0x04}]
	set strand 1
    }

    # Problem 3: low template/sequence coverage
    if {[expr {$problem_bits&0x04}]} {
	# Primer walk, any strand
	set type [expr {$type|0x04}]
    }

    # Problem 4: no BigDye terminator sequences
    if {[expr {$problem_bits&0x08}]} {
	# Primer walk, any strand
	set type [expr {$type|0x04}]
    }

    # Problem 5: do nothing
    if {[expr {$problem_bits&0x10}]} {
	# Do nothing
	set type [expr {$type|0x01}]
    }

    return [expr $type | ($strand << 16) | ($chem << 24)]
}

# All solutions are base around primer walking
# Just minor unresolved problems
proc find_solutions_sanger2 {base_bits problem_bits} {
    # Optimisation - no problems => no solutions required
    if {$problem_bits == 0} {return 0}

    set type 0;	  # None
    set strand 0; # Any
    set chem 7;	  # Big Dye terminator

    # Problem 5: low consensus confidence
    if {[expr {$problem_bits&0x10}]} {
	set type [expr {$type|0x04}]
    }

# The following optimisation sometimes prunes the best experiment.
#
#    # Problem 5: low consensus confidence
#    if {[expr {$problem_bits&0x10}]} {
#	if {[expr {!($base_bits&0x0001)}]} {
#	    # No top strand
#	    # Primer walk, top strand
#	    set type [expr {$type|0x04}]
#	    set strand 1
#	} elseif {[expr {!($base_bits&0x0002)}]} {
#	    # No bot strand
#	    # Primer walk, bot strand
#	    set type [expr {$type|0x04}]
#	    set strand 2
#	} else {
#	    # Primer walk, any strand
#	    set type [expr {$type|0x04}]
#	    set strand 0
#	}
#    }

    return [expr $type | ($strand << 16) | ($chem << 24)]
}

# -----------------------------------------------------------------------------
# Dumps the problem arrays.
proc dump_problem {io fin fd cnum} {
    set pos 1
    set plist [$fin dump_problems]
    set rname [left_gel $io $cnum]
    set pstart 1
    set plast [lindex $plist 0]
    foreach prob $plist {
	if {$prob != 0} {
	    if {$prob != $plast} {
		puts $fd "$rname $pstart $pstart..[expr {$pos-1}] problem flags=$plast"
		set pstart $pos
		set plast $prob
	    }
	}
	incr pos
    }
    if {$prob != 0} {
	puts $fd "$rname $pstart $pstart..[expr {$pos-1}] problem flags=$plast"
    }
}

# -----------------------------------------------------------------------------
# Main entry point

set add_tags ""
set dump_problems 0
set skip_solutions 0
set from 0
set to 0
set contig ""
set all_contigs 1

while {$argc > 0 && "[string index [lindex $argv 0] 0]" == "-"} {
    set arg [lindex $argv 0];
    set argv [lrange $argv 1 $argc]
    incr argc -1;

    if {$arg == "--"} {
	break;

    } elseif {$arg == "-add_tags"} {
	set add_tags [lindex $argv 0]
	set argv [lrange $argv 1 $argc]
	incr argc -1;

    } elseif {$arg == "-dump_problems"} {
	set dump_problems 1
	set prob_fd [open [lindex $argv 0] w]
	set argv [lrange $argv 1 $argc]
	incr argc -1

    } elseif {$arg == "-skip_solutions"} {
	set skip_solutions 1

    } elseif {$arg == "-contig" || $arg == "-contigs"} {
	set contig [lindex $argv 0]
	set argv [lrange $argv 1 $argc]
	incr argc -1;
	set all_contigs 0

    } elseif {$arg == "-from"} {
	set from [lindex $argv 0]
	set argv [lrange $argv 1 $argc]
	incr argc -1;

    } elseif {$arg == "-to"} {
	set to [lindex $argv 0]
	set argv [lrange $argv 1 $argc]
	incr argc -1;

    } else {
	break
    }
}

puts "*** DBname [lindex $argv 0]"
puts "*** Date	 [clock format [clock seconds]]"
puts ""

# Load vector file
set vector ""
if {![catch {set fd [open vector]}]} {
    set vector [read $fd]
    regsub -all "\[#>\]\[^\n\]*\n" $vector {} vector
    regsub -all {[^ACGTNacgtn]} $vector {} vector
    close $fd
    puts "Read [string length $vector] bytes of vector"
    puts ""
}

# Open database, choose contig, select consensus mode
foreach {dbname dbvers} [split [lindex $argv 0] .] {}
if {$add_tags != ""} {
    set io [open_db -name $dbname -version $dbvers -access rw]
} else {
    set io [open_db -name $dbname -version $dbvers -access r]
}

set consensus_cutoff 0.02
set quality_cutoff 1
set consensus_mode 2
set min_contig_len 2000

#
# Base classification bits used for this example:
#
# 0	Has data on top strand
# 1	Has data on bottom strand
# 2	Covered by 2 or more sequences (no effect in this example due to bit 3)
# 3	Covered by 2 or more templates
# 4	Consensus confidence is 15 or more
# 5	Is a BigDye or ET terminator sequence
# 6	Is a non-BigDye/ET terminator sequence
# 7	Is a dye-primer sequence
# 8	At contig left end
# 9	At contig right end


set class_bits {
    {0 strand_top}
    {1 strand_bottom}
    {2 sequence_depth_gt 1}
    {3 template_depth_gt 1}
    {4 confidence_gt 29}
    {5 chemistry 7}
    {5 chemistry 13}
    {5 chemistry 15}
    {5 chemistry 17}
    {6 chemistry 1}
    {6 chemistry 3}
    {6 chemistry 5}
    {6 chemistry 9}
    {6 chemistry 11}
    {7 chemistry 0}
    {7 chemistry 2}
    {7 chemistry 4}
    {7 chemistry 6}
    {7 chemistry 8}
    {7 chemistry 10}
    {7 chemistry 12}
    {7 chemistry 14}
    {7 chemistry 16}
    {8 contig_left_end}
    {9 contig_right_end}
}

set clist [CreateAllContigList $io]
set db [io_read_database $io]
set num_contigs [keylget db num_contigs]
set tcl [db_info t_contig_length $io]
set maxseq [expr {round(($tcl + 20*$num_contigs)*1.1)}]

# Allocate a 'finish' Tcl_Obj object. (Note that this can grow quite big.)
# This contains consensus, confidence values, virtual sequences, etc.
# It's the main data block passed between the various finishing functions.
finish .f \
    -io $io \
    -check_contigs $clist \
    -external_seq $vector \
    -use_avg_insert 0 \
    -reseq_cost 1.2 \
    -reseq_length 500 \
    -reseq_nsolutions 4 \
    -long_cost 2.0 \
    -long_length 700 \
    -long_nsolutions 4 \
    -vpwalk_cost 2.0 \
    -cpwalk_cost 8.0 \
    -pwalk_search_dist 100 \
    -pwalk_max_match 10 \
    -pwalk_osp_score 28 \
    -pwalk_noligos 8 \
    -pwalk_ntemplates 4 \
    -pwalk_offset1 130 \
    -pwalk_offset2 55 \
    -pwalk_length 400 \
    -pwalk_nsolutions 1 \
    -pwalk_seq_gap 55 \
    -pwalk_consistent_only 0 \
    -pwalk_max_err 0.02 \
    -pwalk_min_qual 30 \
    -pwalk_max_err2 0.1 \
    -pwalk_min_qual2 15 \
    -pwalk_end_dist 700 \
    -pwalk_use_template 1 \
    -pwalk_use_template_score 0 \
    -pwalk_prob_mask 7 \
    -mandatory_ratio 0.4 \
    -prob_mandatory 7 \
    -max_score_drop 0.2 \
    -min_template_score 0.4 \
    -min_score 6.0 \
    -avail_template_file avail_templates \
    -skip_template_file skip_templates \
    -pscores {0.5 0.5 0.8 0.8 1.5 0.01} \
    -mscores {10 10 0.8 0.8 1.5 0.01} \
    -dust_level 14 \
    -min_extension 50 \
    -primer_min_tm 50 \
    -primer_max_tm 56 \
    -primer_opt_tm 53 \
    -primer_min_len 17 \
    -primer_max_len 23 \
    -primer_opt_len 18 \
    -primer_min_gc 30 \
    -primer_max_gc 70 \
    -primer_self_any 6 \
    -primer_self_end 3 \
    -primer_gc_clamp 1 \
    -primer_max_poly_x 4 \
    -primer_max_end_stability 9 \
    -output_file $dbname.$dbvers.experiments \
    -debug0 2 \
    -debug1 2 \
    -debug2 2 \
    -debug3 2 \
    -debug4 2 \
    -debug5 2 \
    -debug6 2 \
    -debug7 2 \
    -debug8 2 \
    -debug9 2

if {$add_tags != ""} {
    .f configure -pwalk_tag_type $add_tags
}

# Produce a list of contigs to process. "contigs" is a list of contig numbers
# The -min_contig_len option only applies when -contig(s) is not explicitly
# used.
set contigs ""
if {$contig != ""} {
    foreach c $contig {
	set cnum [db_info get_contig_num $io $c]
	if {$cnum == -1} {
	    puts "Unknown contig $c"
	    exit
	}
	lappend contigs $cnum
    }
} else {
    for {set cnum 1} {$cnum <= $num_contigs} {incr cnum} {
	set c [io_read_contig 1 $cnum]
	if {[keylget c length] < $min_contig_len} {
	    continue
	}
	lappend contigs $cnum
    }
}

# Loop through selected contigs
foreach cnum $contigs {
    set c [io_read_contig 1 $cnum]

    set start [expr {$from ? $from : 1}]
    set end   [expr {$to ? $to : [keylget c length]}]

    set contig "#[keylget c left]"

    puts "### CONTIG ID $contig  (=$cnum [left_gel $io $cnum]) ###"

    # Reconfigure the finish object to work on this specific contig
    .f configure \
	-io $io \
	-contig "{$contig $start $end}"

    # Classify the bases to produce a bit-pattern
    .f classify \
	-bits $class_bits

    # -- Round 1

    # Identify problems and solutions from bit-classifications
    .f find_problems \
	-problem_command finishing_rules_sanger \
	-solution_command find_solutions_sanger1 \
	-tag_types {MASK CVEC}

    if {$dump_problems} {
	dump_problem $io .f $prob_fd $cnum
    }

    if {$skip_solutions} {
	continue
    }

    # Produce solutions
    set tags [.f implement_solutions \
		  -tag_types {AMBG OLIG PRIM MASK CVEC}]
    
    # -- Round 2

    # Identify problems and solutions from bit-classifications
    .f find_problems \
	-problem_command finishing_rules_sanger \
	-solution_command find_solutions_sanger2 \
	-tag_types {MASK CVEC}

    # Produce solutions
    append tags " [.f implement_solutions \
		      -tag_types {OLIG PRIM MASK CVEC}]"

    if {$add_tags != ""} {
	add_tags -io $io -tags $tags
    }
    set fd [open tags a]
    puts $fd $tags
    close $fd
    
    flush stdout
}

if {$dump_problems} {
    close $prob_fd
}

# Reclaim memory
.f delete

close_db -io $io

exit

