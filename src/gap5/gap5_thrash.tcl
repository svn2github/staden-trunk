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
    break_contig -io $io -contig $crec -pos $p
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
    set pos [expr {$p*$l + $s}]
    $c delete

    puts "/// join contig $crec1 to $crec2 at offset $pos ///"
    join_contigs -io $io -contig1 $crec1 -contig2 $crec2 -pos1 $pos
    return

    # defunct version, it crashed anyway!
    package require Tk 
    frame .e
    editor .e.e1
    set io1 [$io child]
    ednames .e.n1
    .e.e1 init $io $crec1 $crec1 0 .e.n1

    editor .e.e2
    set io2 [$io child]
    ednames .e.n2
    .e.e2 init $io $crec2 $crec2 0 .e.n2

    sheet .e.d
    .e.e1 link_to .e.e2 .e.d

    .e.e1 xview moveto 0
    .e.e2 xview moveto $p

    grid .e.n1 .e.e1 -sticky nsew
    grid .e.d        -sticky nsew -columnspan 2
    grid .e.n2 .e.e2 -sticky nsew
    grid .e -sticky nsew
    update idletasks
    update

    #.e.e1 decr_contig
    .e.e2 incr_contig
    puts joining
    .e.e1 join
    puts joined

    for {set i 0} {$i < 10} {incr i} {
	update idletasks
	update
    }
    puts e2
    destroy .e.e2
    puts e2-gone

    destroy .e
    puts io1-close

    $io1 close
    $io2 close
}

proc test_insertions {io} {
    set crec [random_contig $io]
    puts "/// inserting to contig $crec ///"
    set c [$io get_contig $crec]
    set l [$c get_length]
    set s [$c get_start]
    set e [$c get_end]
    for {set i 0} {$i < 100} {incr i; incr l 1} {
	set p [expr {int(rand()*$l)+$s}]
	#puts "   ///Ins $p/$s..$e"
	$c insert_base $p * 11
	# set e [$c check]
	# if {$e} {
	#     puts "check_after=$e"
	#     puts "\n\n\n ...  ooo  OOO  ooo  ...  ooo  OOO  ooo  ...  ooo  OOO  ooo  ...   \n\n\n"
	#     break
	# }
    }
    $c delete
}

proc test_deletions {io} {
    set crec [random_contig $io]
    puts "/// deleting from contig $crec ///"
    set c [$io get_contig $crec]
    set l [$c get_length]
    set s [$c get_start]
    for {set i 0} {$i < 100 && $l > 10} {incr i; incr l -1} {
	set p [expr {int(rand()*$l)+$s}]
	# puts "   ///Del $p/$s+$l"
	$c delete_base $p
	set e [$c check]
	# if {$e} {
	#     puts "check_after=$e"
	#     puts "\n\n\n ...  ooo  OOO  ooo  ...  ooo  OOO  ooo  ...  ooo  OOO  ooo  ..    \n\n\n"
	#     break
	# }
    }
    $c delete
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
$io check 0 2

# Perform N edits and keep checking.
for {set cycle 0} {$cycle < 100} {incr cycle} {
    set r [expr int(rand()*5)]
    puts "///$cycle r=$r"

    # Other tests to do:
    # - Add tags
    # - Delete tags
    # - Adjust seq clips
    # - Compute & cache consensus
    switch $r {
	0 { test_complement $io }
	1 { test_break $io }
	2 { test_join $io }
	3 { test_insertions $io }
	4 { test_deletions $io }
    }

    $io flush
    set err [$io check 0 2]
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
