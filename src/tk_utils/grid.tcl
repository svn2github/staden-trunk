# insert rows or columns into a grid
#  grid:  the geometry master
#  what:  row or column
#  index: where to insert
#  count: how many rows/cols to insert
      
proc grid_insert {grid what index {count 1}} {
    #puts "grid_insert $grid $what $index"

    foreach slave [grid slaves $grid] {
	array set info [grid info $slave]

	#puts "$info(-$what) $index"
	if {$info(-$what) >= $index} {
                 incr info(-$what) $count
	    eval {grid $slave} [array get info]
	} elseif {$info(-$what) + $info(-${what}span) > $index} {
	    incr info(-${what}span) $count
	    eval {grid $slave} [array get info]
	}
    }
}

proc grid_delete {grid what index {count 1}} {
    #puts "grid_delete $grid $what $index"

   foreach slave [grid slaves $grid] {
	array set info [grid info $slave]
	if {$info(-$what) >= [expr $index + $count]} {
	    incr info(-$what) -$count
	    eval {grid $slave} [array get info]
	} elseif {$info(-$what)+$info(-${what}span) > [expr $index + $count + 1]} {
	    incr info(-$what) -$count
	    eval {grid $slave} [array get info]
	}
    }
}

proc grid_size {grid what} {

    set num 0

    foreach slave [grid slaves $grid] {
	array set info [grid info $slave]
	#puts "slave $slave $info(-$what)"
	if {$info(-$what) > $num} {
	    set num $info(-$what)
	}
    }
    return $num
}

#return list of all windows in row or column
proc grid_find {grid what index} {
    #puts "****************grid_find $grid $what"

    set list ""
    foreach slave [grid slaves $grid] {
	array set info [grid info $slave]
	#puts "$slave $info(-$what)"
	if {$info(-$what) == $index} {
	    lappend list $slave
	}
    }
    #puts "FOUND list $list"
    return $list
}

#insert single element in row column
proc grid_insert_element {grid row column} {
    #puts "grid_insert $grid $row $column"

    foreach slave [grid slaves $grid] {
	array set info [grid info $slave]

	#puts "$info(-$what) $index"
	if {$info(-$what) >= $index} {
                 incr info(-$what) $count
	    eval {grid $slave} [array get info]
	} elseif {$info(-$what) + $info(-${what}span) > $index} {
	    incr info(-${what}span) $count
	    eval {grid $slave} [array get info]
	}
    }
}

