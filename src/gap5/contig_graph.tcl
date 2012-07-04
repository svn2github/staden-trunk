#!/usr/bin/wish
# Draws a contig graph and allows manipulation of links

# 
# The main data structure we hold is:
#
# contig_list:  contig ...                    A list of "contig" items
# contig:       rec name start end link ...   A single contig and a
#                                             list of links
# link:         rec2 pos1 pos2 end1 end2 orientation size type score
#                                             As produced by [$c get_links]
#
# We also have an associative array indexed by contig record ID holding the
# layout information. Each member consists of
#
#     pos y len style
#
# Where pos is in base pairs and y is just a unitless value used to prevent
# clashes on the line.

# Some debug data for pure wish development
if {1} {
proc debug_contig_list {} {
    return [list \
		[list 1 c1 1 100 \
		     [list 2 90 0 1 0 0 50 0 0] \
		     [list 6 90 0 1 0 0 610 0 0] \
		    ] \
		[list 2 c2 1 100 \
		     [list 1 0 90 0 1 0 50 0 0] \
		     [list 3 90 0 1 0 0 50 0 0] \
		     [list 7 90 0 1 0 0 50 0 0] \
		    ] \
		[list 3 c3 1 100 \
		     [list 2 0 90 0 1 0 50 0 0] \
		     [list 4 90 0 1 0 0 50 0 0] \
		    ] \
		[list 4 c4 1 100 \
		     [list 3 0 90 0 1 0 50 0 0] \
		     [list 5 90 0 1 0 0 50 0 0] \
		    ] \
		[list 5 c5 1 100 \
		     [list 4 0 90 0 1 0 50 0 0] \
		     [list 6 90 0 1 0 0 50 0 0] \
		     [list 8 0 90 0 1 0 50 0 0] \
		    ]  \
		[list 6 c6 1 100 \
		     [list 5 0 90 0 1 0 50 0 0] \
		     [list 1 0 90 0 1 0 610 0 0] \
		    ] \
		[list 7 c7 1 100 \
		     [list 2 0 90 0 1 0 50 0 0] \
		     [list 8 90 0 1 0 0 50 0 0] \
		     [list 9 90 0 1 0 0 50 0 0] \
		    ] \
		[list 8 c8 1 100 \
		     [list 7 0 90 0 1 0 50 0 0] \
		     [list 5 90 0 1 0 0 50 0 0] \
		    ] \
		[list 9 c9 1 100 \
		     [list 7 0 90 0 1 0 50 0 0] \
		    ] \
	   ]
}

} else {

proc debug_contig_list {} {
    return [list \
		[list 1 c1 1 100 \
		     [list 2 90 5 1 0 0 300 0 0] \
		    ] \
		[list 2 c2 1 100 \
		     [list 1 5 90 0 1 0 300 0 0] \
		     [list 3 90 5 1 0 0 300 0 0] \
		    ] \
		[list 3 c3 1 100 \
		     [list 2 5 90 0 1 0 300 0 0] \
		    ] \
	       ]
}
}

# Updates $lv() to hold layout information for the contig_list graph
#
# It is indexed by record number. The contents are a list containing
# start, end and y-line.
proc layout {contig_list lv} {
    puts [info level [info level]]
    upvar $lv pos

    # Initialise all positions to 1
    foreach contig $contig_list {
	set crec [lindex $contig 0]
	set pos($crec) 1
    }

    # Update positions based on links.
    # These could be inconsistent with respect to each other
    # (ctg2 is 500 away from ctg1, but ctg1 is 200 away from ctg2), but
    # we'll just end up with one or the other after this.
    foreach contig $contig_list {
	puts c=$contig
	set crec [lindex $contig 0]
	set size [expr {[lindex $contig 3]-[lindex $contig 2]+1}]
	foreach lnk [lrange $contig 4 end] {
	    foreach {c2 p1 p2 e1 e2 o sz} $lnk break;
	    if {$e1 == 1} {
		set pos($c2) [expr {$pos($crec)+$p1+$sz-$p2}]
	    } else {
		set pos($c2) [expr {$pos($crec)+$p1-$sz-$p2}]
	    }
	    puts "Setting pos($c2) to $pos($c2)"
	}
    }

    # Now allocate x1 and x2 instead of just pos.
    set rec_list {}
    foreach contig $contig_list {
	foreach {crec name start end} $contig break
	lappend rec_list $crec
	set size [expr {$end-$start+1}]
	lappend pos($crec) [expr {$pos($crec)+$size-1}]
    }

    # Sort by X
    set rec_list [lsort -command rec_x_sort $rec_list]
    puts $rec_list

    # Allocate non-clashing Y coords.
    set nl 0

    foreach crec $rec_list {
	foreach {x1 x2} $pos($crec) break

	# Try to keep the same line as the input link
	# FIXME:

	set blocked 1
	for {set l 0} {$l < $nl} {incr l} {
	    if {$line($l) < $x1} {
		set blocked 0
		break;
	    }
	}

	if {$blocked} {
	    lappend pos($crec) $nl
	    set line($nl) $x2
	    incr nl 
	} else {
	    lappend pos($crec) $l
	    set line($l) $x2
	}
    }
}

# Callback from lsort -command
proc rec_x_sort {r1 r2} {
    upvar pos p
    set p1 [lindex $p($r1) 0]
    set p2 [lindex $p($r2) 0]

    if {$p1 < $p2} {
	return -1
    } elseif {$p1 > $p2} {
	return +1
    } else {
	return 0
    }
}

# Draw the results
proc draw_contig_links {contig_list coord_var} {
    upvar $coord_var pos

    pack [canvas .c -width 900] -fill both -expand 1

    # Nodes
    foreach r [array names pos] {
	foreach {x1 x2 y} $pos($r) break;
	set y1 [expr {$y * 40 + 50}]
	.c create line $x1 $y1 $x2 $y1 -fill blue -width 10
	set px($r) $x1
	set py($r) $y1
    }

    # Edges
    foreach contig $contig_list {
	set crec [lindex $contig 0]
	set lnk($crec) [lrange $contig 4 end]
    }
    foreach c1 [array names pos] {
	foreach l $lnk($c1) {
	    foreach {c2 p1 p2 e1 e2 o sz} $l break;
	    puts $r=>$l,$px($c1),$px($c2)
	    set x1 [expr {$px($c1)+$p1}]
	    set y1 [expr {$py($c1)}]
	    set x2 [expr {$px($c2)+$p2}]
	    set y2 [expr {$py($c2)}]
	    set d [expr {int(2*log(1+abs($x2-$x1)))}]
	    if {$e1} {
		#.c create line $x1 $y1 [expr {($x1+$x2)/2}] [expr {($y1+$y2)/2+$d}] $x2 $y2 -fill black -width 2 -smooth true
		if {$y1 == $y2} {
		    .c create line $x1 $y1 $x1 [expr {($y1+$y2)/2-$d}] $x2 [expr {($y1+$y2)/2-$d}] $x2 $y2 -fill black -width 2 -smooth true
		} else {
		    .c create line $x1 $y1 $x1 [expr {($y1+$y2)/2+$d}] $x2 [expr {($y1+$y2)/2+$d}] $x2 $y2 -fill black -width 2 -smooth true
		}
	    } else {
		#.c create line $x1 $y1 [expr {($x1+$x2)/2}] [expr {($y1+$y2)/2-$d}] $x2 $y2 -fill orange4 -width 2 -smooth true
	    }
	    #.c create arc [expr {$x1-($x2-$x1)}] $y1 $x2 [expr {$y2+$y2-$y1}] -style arc -outline grey
	}
    }
}

# Draw layout in gap5
proc gap5_layout {} {
    global io

    if {[winfo exists .c]} {
	destroy .c
    }

    set nc [$io num_contigs]

    set l {}
    for {set i 0} {$i < $nc} {incr i} {
	set crec [$io contig_order $i]
	set c [$io get_contig $crec]
	set cc [list $crec [$c get_name] [$c get_start] [$c get_end]]
	foreach x [$c get_links] {puts x=$x;lappend cc $x}
	lappend l $cc
	$c delete
    }

    layout $l coords

    puts "\n\n========"
    puts l=$l
    puts "==Coords=="
    parray coords

    draw_contig_links $l coords

    update idletasks
    foreach {x1 y1 x2 y2} [.c bbox all] break
    set width [winfo width .c]
    set height [winfo height .c]
    puts ".c scale all $x1 $y1 [expr {double($width)/($x2-$x1)}] [expr {double($height)/($y2-$y1)}]"
    .c scale all $x1 $y1 [expr {double($width)/($x2-$x1)}] [expr {double($height)/($y2-$y1)}]
    .c scale all 0 0 .7 .7
    .c move all 20 20
}

if {[info commands g5::open_database] == ""} {
    set contig_list [debug_contig_list]
    layout $contig_list coords
    parray coords

    draw_contig_links $contig_list coords
}