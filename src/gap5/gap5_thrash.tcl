# This simulates user edits with lots of complementing, joining, breaking,
# insertions and deletions.

# The purpose is simply to thrash the code and get contigs in a messy state
# and to stress test the code to make sure it all works fine.

proc random_contig {io} {
    global db
    set nc [$db get_num_contigs]

    while {[set r [expr int(rand()*$nc)]] >= $nc} {
    }

    return [$io contig_order $r]
}

proc test_complement {io} {
    set crec [random_contig $io]
    puts "/// complement contig $crec ///"
    complement_contig -io $io -contigs =$crec
}

proc test_break {io} {
    set crec [random_contig $io]
    set c [$io get_contig $crec]
    set e [$c get_end]
    set s [$c get_start]
    set p [expr {int(rand()*($e-$s+1)+$s)}]
    $c delete
    
    puts "/// break contig $crec ($s..$e) at $p ///"
    catch {break_contig -io $io -contig $crec -pos $p -break_holes 1}
}

proc test_join {io} {
    set crec1 [random_contig $io]
    set crec2 [random_contig $io]
    if {$crec1 == $crec2} return

    set c [$io get_contig $crec1]
    set e [$c get_end]
    set s [$c get_start]
    set l [$c get_length]
    set p [expr {rand()}]
    set pos [expr {int($p*$l + $s)}]
    $c delete

    puts "/// join contig $crec1 to $crec2 at offset $pos ///"
    join_contigs -io $io -contig1 $crec1 -contig2 $crec2 -pos1 $pos
    return

#    # defunct version, it crashed anyway!
#    package require Tk 
#    frame .e
#    editor .e.e1
#    set io1 [$io child]
#    ednames .e.n1
#    .e.e1 init $io $crec1 $crec1 0 .e.n1
#
#    editor .e.e2
#    set io2 [$io child]
#    ednames .e.n2
#    .e.e2 init $io $crec2 $crec2 0 .e.n2
#
#    sheet .e.d
#    .e.e1 link_to .e.e2 .e.d
#
#    .e.e1 xview moveto 0
#    .e.e2 xview moveto $p
#
#    grid .e.n1 .e.e1 -sticky nsew
#    grid .e.d        -sticky nsew -columnspan 2
#    grid .e.n2 .e.e2 -sticky nsew
#    grid .e -sticky nsew
#    update idletasks
#    update
#
#    #.e.e1 decr_contig
#    .e.e2 incr_contig
#    puts joining
#    .e.e1 join
#    puts joined
#
#    for {set i 0} {$i < 10} {incr i} {
#	update idletasks
#	update
#    }
#    puts e2
#    destroy .e.e2
#    puts e2-gone
#
#    destroy .e
#    puts io1-close
#
#    $io1 close
#    $io2 close
}

proc test_insertions {io {cycle 0}} {
    set crec [random_contig $io]
    puts "/// inserting to contig $crec ///"
    set c [$io get_contig $crec]
    set l [$c get_visible_length]
    set s [$c get_visible_start]
    set e [$c get_visible_end]
    set end 100
    #if {$cycle == 67} {set end 50}
    for {set i 0} {$i < $end} {incr i; incr l 1} {
	set p [expr {int(rand()*$l)+$s}]
	#puts "   ///Ins $p/$s..$e"
	$c insert_base $p * 11
	#$io flush
	#if {[$io check 0 2] != 0} exit
    }
    $c delete
}

proc test_deletions {io {cycle 0}} {
    set crec [random_contig $io]
    puts "/// deleting from contig $crec ///"
    set c [$io get_contig $crec]
    set l [$c get_visible_length]
    set s [$c get_visible_start]
    set e [$c get_visible_end]
    set end 100
    #if {$cycle == 40} {set end 21}
    for {set i 0} {$i < $end && $l > 10} {incr i; incr l -1} {
	set p [expr {int(rand()*$l)+$s}]
	puts "   ///Del $p/$s..$e"
	$c delete_base $p
	#$io flush
	#if {[$io check 0 1] != 0} {exit}
    }
    $c delete
}

proc test_tag_creation {io} {
    set crec [random_contig $io]
    puts "/// creating tags in contig $crec ///"

    set c [$io get_contig $crec]
    set cstart [$c get_visible_start]
    set cend   [$c get_visible_end]
    set crec   [$c get_rec]
    
    # Get sequence list
    set seqs ""

    foreach l [$c seqs_in_range $cstart $cend] {
	lappend seqs [lindex $l 2]
    }
    $c delete

    if {$seqs == ""} return

    # Create tags
    for {set i 0} {$i < 50} {incr i} {
	if {rand() > 0.5} {
	    set msg ""
	} else {
	    set msg [string repeat a [expr {int(rand()*100)}]]
	}
	if {rand() > 0.5} {
	    # Consensus tag
	    set otype 17
	    set orec $crec
	    set st   [expr {int(rand()*($cend-$cstart))+$cstart}]
	    if {rand() > 0.5} {
		set en $st
	    } else {
		set en   [expr {int(rand()*($cend-$st))+$st}]
	    }
	} else {
	    # Sequence tag
	    set otype 18
	    set orec [lindex $seqs [expr {int(rand()*[llength $seqs])}]]
	    set s [$io get_seq $orec]
	    set len [expr {abs([$s get_len])}]
	    $s delete
	    set st   [expr {int(rand()*$len)}]
	    if {rand() > 0.5} {
		set en $st
	    } else {
		set en   [expr {int(rand()*($len-$st))+$st}]
	    }
	}
	set rec [$io new_anno_ele $otype $orec $st $en]
	# puts "Creating tag $rec on $otype/$orec at $st..$en"
	set t [$io get_anno_ele $rec]
	$t set_comment $msg
	$t set_type COMM
	$t delete
    }
}

proc test_tag_deletion {io} {
    # Randomly remove half the tags in any contig

    set crec [random_contig $io]
    puts "/// removing tags from contig $crec ///"

    set c [$io get_contig $crec]
    set cstart [$c get_start]
    set cend   [$c get_end]
    set crec   [$c get_rec]

    # Pick tags to go
    set to_del {}
    foreach l [$c anno_in_range $cstart $cend] {
	if {rand() > 0.5} {
	    lappend to_del [lindex $l 2]
	}
    }
    $c delete

    # Delete them
    foreach rec $to_del {
	# puts "Removing tag $rec"
	set tag [$io get_anno_ele $rec]
	$tag remove
    }
}

proc test_clipping {io} {
    # Skip for now
    return;

    set crec [random_contig $io]
    puts "/// Adjusting seq clips in contig $crec ///"

    set c [$io get_contig $crec]
    set cstart [$c get_start]
    set cend   [$c get_end]
    set crec   [$c get_rec]
    
    # Get sequence list
    set seqs ""
    foreach l [$c seqs_in_range $cstart $cend] {
	lappend seqs [lindex $l 2]
    }
    $c delete

    set nseq [llength $seqs]
    for {set i 0} {$i < 100 && $i < $nseq} {incr i} {
	set srec [lindex $seqs [expr {int(rand()*$nseq)}]]

	set s [$io get_seq $srec]
	set len [expr {abs([$s get_length])}]
	set l [expr {int(rand()*$len)+1}]
	set r [expr {int(rand()*$len)+1}]
	if {$l == $r} continue
	if {$l > $r} {
	    set t $l
	    set l $r
	    set r $t
	}

	#puts "SEQ$srec set_clips $l $r"
	$s set_clips $l $r
	$s delete
    }
}

proc test_disassembly {io} {
    set crec [random_contig $io]
    puts "/// Disassembly from contig $crec ///"

    set c [$io get_contig $crec]
    set cstart [$c get_start]
    set cend   [$c get_end]
    set crec   [$c get_rec]
    
    # Get sequence list
    set seqs ""
    foreach l [$c seqs_in_range $cstart $cend] {
	lappend seqs [lindex $l 2]
    }
    $c delete

    set nseq [llength $seqs]
    if {$nseq/2 > 100} {
	set ndel 100
    } else {
	set ndel [expr {$nseq/2}]
    }

    if {$ndel == 0} {
	return
    }

    # Pick ndel reads from seq list
    set srec {}
    for {set i 0} {$i < $ndel} {incr i} {
	set s [expr {int(rand()*$nseq)}]
	lappend srec #[lindex $seqs $s]
	set seqs [lreplace $seqs $s $s]
	incr nseq -1
    }

    set opt1 [expr {rand()>0.5}]
    set opt2 [expr {rand()>0.5}]
    puts "Disassembling opt $opt1, $opt2, $srec"

    # Disassemble them
    set r [disassemble_readings \
	       -io $io \
	       -readings $srec \
	       -move 2 \
	       -remove_holes $opt1 \
	       -duplicate_tags $opt2]
    puts "Disassemble returned $r"
}

proc test_consensus {io} {
    set cl {}

    for {set i 0} {$i < 10} {incr i} {
	set crec [random_contig $io]
	set c [$io get_contig $crec]
	set cstart [$c get_start]
	set cend   [$c get_end]
	$c delete

	set len [expr {$cend-$cstart+1}]
	set l [expr {int(rand()*$len)+$cstart}]
	set r [expr {int(rand()*$len)+$cstart}]

	if {$l > $r} {
	    set t $l
	    set l $r
	    set r $t
	}

	puts "/// Computing consensus for contig $crec at $l to $r ///"
	calc_consensus -io $io -contigs [list =$crec $l $r]
    }
}

#-----------------------------------------------------------------------------
# MAIN

# Startup code
source $env(STADTABL)/shlib.conf
load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}
load_package tk_utils
tk_utils_init
load_package gap5

# Ensure reproducability
if {[lindex $argv 1] != ""} {
    expr srand([lindex $argv 1])
} else {
    expr srand(100)
}

# Open DB
exec cp "[lindex $argv 0].g5d" _tmp.g5d
exec cp "[lindex $argv 0].g5x" _tmp.g5x
set db [lindex $argv 0]
if {[catch {set io [g5::open_database -name _tmp -access rw]} err]} {
    puts stderr "Couldn't open database '$db': $err"
    exit 1
}
set db [$io get_database]
#$io check 0 2
$io debug_level 1

if {[llength $argv] > 2} {
    set ncycles [lindex $argv 2]
} else {
    set ncycles 100
}

# Perform N edits and keep checking.
for {set cycle 0} {$cycle < $ncycles} {incr cycle} {
    set r [expr int(rand()*10)]
    #if {$r != 3} continue

    puts "///$cycle r=$r"

    # Other tests to do:
    # - Compute & cache consensus
    switch $r {
	0 { test_complement $io }
	1 { test_break $io }
	2 { test_join $io }
	3 { test_insertions $io $cycle }
	4 { test_deletions $io $cycle }
	5 { test_tag_creation $io }
	6 { test_tag_deletion $io }
	7 { test_clipping $io }
	8 { test_consensus $io }
	9 { test_disassembly $io }
    }

    $io flush

    set err [$io check 0 1]
    if {$err != 0} {
	$io close
	puts stderr "ERROR: corrupted database\n"
	exit 1
    }
}

$io flush

puts FINAL_CHECK:[$io check 0 2]

$io close
exit 0
