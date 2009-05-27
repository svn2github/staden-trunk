proc ListLibraries {io} {
    set t [xtoplevel .list_libraries]
    if {$t == ""} return
    wm title $t "List Libraries"
    wm geometry $t 550x200

    # Create our tablelist
    tablelist $t.list \
	-columns {7 Index 8 "Name" 10 "Count" 8 "Type" 10 "Insert size" 6 "s.d." 15 "Orientation"} \
        -labelcommand tablelist::sortByColumn \
	-exportselection 0 \
	-stretch 0 \
        -yscrollcommand [list $t.yscroll set]

    $t.list columnconfigure 0 -sortmode integer
    $t.list columnconfigure 2 -sortmode integer
    $t.list columnconfigure 4 -sortmode integer
    $t.list columnconfigure 5 -sortmode real
    
    # Add scrollbars
    scrollbar $t.yscroll -command "$t.list yview"

    # Layout
    grid columnconfigure $t 0 -weight 1
    grid rowconfigure $t 0 -weight 1
    grid $t.list $t.yscroll -sticky nsew

    # Add some bindings
    bind [$t.list bodypath] <<use>> "+ListLibrariesDetailed $io %W $t.list %x %y"

    # Populate the list
    ListLibrariesPopulate $io $t.list
}

;proc ListLibrariesPopulate {io w} {
    # Clear any existing data
    $w selection clear 0 end
    $w delete 0 end

    set db [$io get_database]
    set nc [$db get_num_libraries]
    for {set i 0} {$i < $nc} {incr i} {
	set lib [$io get_library [$db get_library_rec $i]]
	set name   "lib\#$i"
	set orient  [$lib get_orient]
	set type    [$lib get_machine]
	set mean    [$lib get_insert_size]
	set sd      [$lib get_insert_sd]
	set count   [$lib get_count]
	set type [lindex [list unknown sanger solexa solid 454] $type]
	
	# Find most likely library orientation
	set max 0
	set o 0
	set tot 0
	for {set j 0} {$j < 3} {incr j} {
	    if {$max < [lindex $count $j]} {
		set max [lindex $count $j]
		set o $j
	    }
	    incr tot [lindex $count $j]
	}
	set mean [lindex $mean $o]
	set sd   [lindex $sd   $o]
	set or   [lindex [list "-> <-" "<- ->" "-> -> / <- <-"] $o]

	$w insert end [list $i $name $tot $type $mean $sd $or]
	$lib delete
    }

    if {[$w sortcolumn] != -1} {
        $w sortbycolumn [$w sortcolumn] -[$w sortorder]
    }

    #$db delete
}

# Double-click callback from the library list.
;proc ListLibrariesDetailed {io w t x y} {
    # Move x,y to parent tablelist coordinates
    incr y [lindex [split [winfo geometry $w] +] 2]
    incr x [lindex [split [winfo geometry $w] +] 1]

    set row [$t containing $y]
    if {$row == -1} return

    incr row
    set idx [lindex [$w get $row.0 $row.end] 0]
    set db [$io get_database]
    set rec [$db get_library_rec $idx]

    draw_lib $io $t.lib_$rec $rec
}

# Detailed data about a library including the distribution graph
;proc draw_lib {io w rec} {
    if {[xtoplevel $w] == ""} return

    set lib [$io get_library $rec]
    set dist [$lib get_dist]

    set c $w.c
    scrollbar $w.xs -orient horiz -command "$c xview"
    scrollbar $w.ys -orient vert  -command "$c yview"
    canvas $c \
	-xscrollcommand "$w.xs set" \
	-yscrollcommand "$w.ys set" \
	-height 400 \
	-width 600
    grid rowconfigure    $w 0 -weight 1 -minsize 0
    grid columnconfigure $w 0 -weight 1 -minsize 0
    grid $c    -row 0 -column 0 -sticky nsew
    grid $w.xs -row 1 -column 0 -sticky nsew
    grid $w.ys -row 0 -column 1 -sticky nsew

    # Create the plots
    set maxy 1
    set maxx 1
    for {set j 0} {$j < 3} {incr j} {
	set line [list 0 0]
	set col [lindex {black red blue} $j]
	foreach {x y} [lindex $dist $j] {
	    lappend line $x $y
	    if {$maxy < $y} {set maxy $y}
	    if {$maxx < $x} {set maxx $x}
	}
	lappend line [lindex $line end-1] 0
	eval $c create line $line -tags type_$j -fill $col
    }

    $c scale all 0 0 1 -1
    $c move all 0 $maxy
    $c scale all 0 0 1 [expr {400.0/$maxy}]

    # Create the ruler
    puts $maxx,$maxy
    for {set x 0} {$x < $maxx+100} {incr x 10} {
	if {[expr {$x%1000}] == 0} {
	    set h 40
	} elseif {[expr {$x%100}] == 0} {
	    set h 20
	} else {
	    set h 10
	}
	$c create line $x 0 $x $h -tags ruler
    }

    $c configure -scrollregion [$c bbox all]

    $lib delete
}