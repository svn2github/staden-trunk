#
# Creates an keybox with an associated label.
# Command line arguments are as per the frame widget, with the addition of:
#	-state	 normal/disabled
#       -title   title
#
proc keybox {path args} {
    global key_

    # Create the frame, label and keybox
    frame $path -class KeyBox 
    frame $path.h

    #a terribly hack to work around the problem of the key frames not packing
    #to the north
    grid $path.h -row 0 -column 0 -sticky nw
    grid rowconfigure $path 0 -weight 1
    grid columnconfigure $path 0 -weight 1

    if {$args != ""} {
	eval $path $args
    }
    return $path
}

#HACK - to tidy up
proc key_enter {path args} {

    $path configure -cursor hand1
    eval $args
}

proc key_leave {path args} {
    
    #puts "KEYLEAVE $path"
    eval $args
}

proc key_drop {path X Y args} {

    eval $args $X $Y
}

proc key_motion {path X Y args } {

    eval $args $X $Y
}

proc key_menu {path X Y args} {
    
    if [winfo exists $path.m] {destroy $path.m}

    set m [create_popup $path.m Commands]
    
    eval $args $m

    tk_popup $m $X $Y

}

proc keybox_grid {path} {
    global $path.keynum $path.keyarray

    set win_ht [winfo height $path]
    if {![winfo exists $path.h.f0]} {
	return
    }

    set key_ht [winfo height $path.h.f0]
    set num_rows [expr ($win_ht / $key_ht) - 1]

    #puts "win_ht $win_ht key_ht $key_ht num_rows $num_rows"
    #always have at least 1 row
    if {$num_rows < 1} {
	set num_rows 1
    }

    for {set i 0} {$i < [set $path.keynum]} {incr i} {
	set column [expr ($i / $num_rows) * 2]
	set row [expr $i % $num_rows]
	set cnt [lindex [set $path.keyarray($i)] 0]
	grid forget $path.h.f$cnt
	grid forget $path.h.l$cnt

	grid $path.h.f$cnt -row $row -column $column
	grid $path.h.l$cnt -row $row -column [expr $column + 1] -sticky w
		    
    }
}

proc keybox_entryconf { args } {
    global KeyArray

    set in_arg 0
    set arglist ""
    array set KeyArray {text "" background "" enter "" motion "" leave "" \
			    drop "" menu "" width "-width 20" height "-height 20"}

    # Process command line args
    foreach i $args {
	if {$in_arg} {
	     if {$option == "-text"} {
		 array set KeyArray "text \"{$i}\""
	     } elseif {$option == "-background"} {
		array set KeyArray "background \"-background {$i}\""
	     } elseif {$option == "-width"} {
		array set KeyArray "width \"-width {$i}\""
	     } elseif {$option == "-height"} {
		array set KeyArray "height \"-height {$i}\""
	     } elseif {$option == "-enter"} {
		array set KeyArray "enter \"$i\""
	     } elseif {$option == "-leave"} {
		array set KeyArray "leave \"$i\""
	     } elseif {$option == "-drop"} {
		array set KeyArray "drop \"$i\""
	     } elseif {$option == "-motion"} {
		array set KeyArray "motion \"$i\""
	     } elseif {$option == "-menu"} {
		array set KeyArray "menu \"$i\""
	     } else {
		lappend arglist $option $i
	     }
	     set in_arg 0
	} else {
	    set option $i
	    set in_arg 1
	}
    }
}

proc keybox_entryconfigure {path index args} {
    global $path.keynum $path.keyarray KeyArray

    if {[regexp {^[0-9]+$} "$index"]} {
	set cnt $index
    } else {
	for {set i 0} {$i < [set $path.keynum]} {incr i} {
	    #puts "index {$index} [lindex [set $path.keyarray($i)] 1]"

	    if {[string compare [lindex [set $path.keyarray($i)] 1] "$index"] == 0} {
		set cnt $i
		break
	    }
	}
    }

    set cnt [lindex [set $path.keyarray($cnt)] 0]

    #configure keybox
    eval keybox_entryconf $args
    
    if {"$KeyArray(background)" != ""} {
	eval $path.h.f$cnt config $KeyArray(background)
    }

    if {"$KeyArray(text)" != ""} {
	eval $path.h.l$cnt config -text $KeyArray(text)
    }
    

}


#keynum - number of keys 
#keycnt - unique key identifier - always increasing
#keyarray - array containing the unique identifier and the unique text
proc keybox_add {path args} {
    global key_
    global $path.keynum $path.keyarray $path.keycnt
    global KeyArray

    #configure keybox
    eval keybox_entryconf $args

    if {![info exists $path.keynum]} {
	set $path.keynum 0
	set $path.keycnt 0
	set $path.keyarray(0) ""
    }

    #add new element to end of array
    set cnt [set $path.keycnt]

    #add unique text to array to allow for searching later on
    #cnt is the name given to the frame and label which remains constant even
    #when the array gets shuffled around due to deletions

    set $path.keyarray([set $path.keynum]) "$cnt $KeyArray(text)"

    eval frame $path.h.f$cnt $KeyArray(height) $KeyArray(width)

    if {"$KeyArray(background)" != ""} {
	eval $path.h.f$cnt config $KeyArray(background)
    }

    label $path.h.l$cnt 
    if {"$KeyArray(text)" != ""} {
	eval $path.h.l$cnt config -text $KeyArray(text)
    }

    #puts "!!!!!!!!!!!!!!!$path keynum [set $path.keynum] cnt $cnt"

   # set m [menu $path.h.l$cnt.opts -tearoff 0]

    grid $path.h.f$cnt -row [set $path.keynum] -column 0
    grid $path.h.l$cnt -row [set $path.keynum] -column 1 -sticky w

    bind $path.h.f$cnt <Any-Enter> "key_enter %W $KeyArray(enter)"
    bind $path.h.f$cnt <B2-Leave> "key_leave %W $KeyArray(leave)"
    bind $path.h.f$cnt <Alt-Button-1> "key_leave %W $KeyArray(leave)"
    bind $path.h.f$cnt <<move-release>> "key_drop %W %X %Y $KeyArray(drop)"
    bind $path.h.f$cnt <<move-drag>> "key_motion %W %X %Y $KeyArray(motion)"

    #add leave binding to left button to allow drag and drop
    bind $path.h.f$cnt <B1-Leave> "key_leave %W $KeyArray(leave)"

    bind $path.h.f$cnt <<menu>> "key_menu %W %X %Y $KeyArray(menu)"
    
    incr $path.keynum
    incr $path.keycnt

    #HACK - don't like this    set r_win [winfo parent $raster]
    update idletasks
    keybox_grid $path

    #for {set i 0} {$i < [set $path.keynum]} {incr i} {
	#puts "index $i  [lindex [set $path.keyarray($i)] 1]"
	#puts "grid [grid info $path.h.f[lindex [set $path.keyarray($i)] 0]]"
    #}
}

#HACK to finish - need to be able to delete from to
proc keybox_delete {path index1 args} {
    global $path.keyarray $path.keynum

    set cnt -1
    
    if {[regexp {^[0-9]+$} "$index1"]} {
	set cnt $index1
    } else {
	for {set i 0} {$i < [set $path.keynum]} {incr i} {
	    #puts "index {$index1} [lindex [set $path.keyarray($i)] 1]"

	    if {[string compare [lindex [set $path.keyarray($i)] 1] "$index1"] == 0} {
		set cnt $i
		break
	    }
	}
    }

    if {$cnt == -1} {
	return
    }
    #puts "?????????keybox_delete $path.h.f[lindex [set $path.keyarray($cnt)] 0] cnt $cnt"
    bind $path.h.f[lindex [set $path.keyarray($cnt)] 0] <Any-Enter> {}
    bind $path.h.f[lindex [set $path.keyarray($cnt)] 0] <B2-Leave> {}
    bind $path.h.f[lindex [set $path.keyarray($cnt)] 0] <Alt-Button-1> {}
    bind $path.h.f[lindex [set $path.keyarray($cnt)] 0] <<move-release>> {}
    bind $path.h.f[lindex [set $path.keyarray($cnt)] 0] <<move-drag>> {}

    destroy $path.h.f[lindex [set $path.keyarray($cnt)] 0]
    destroy $path.h.l[lindex [set $path.keyarray($cnt)] 0]

    #shift all cells up
    for {set i [expr $cnt + 1]} {$i < [set $path.keynum]} {incr i} {
	set $path.keyarray($cnt) [set $path.keyarray($i)]
	incr cnt
    }

    incr $path.keynum -1

    if {[set $path.keynum] > 0 } {
	keybox_grid $path
    } else {
	set $path.keycnt 0
    }
}

proc keybox_entrycget {path index option} {
    global $path.keyarray $path.keynum 

    if {[regexp {^[0-9]+$} "$index"]} {
	set cnt $index
    } else {
	for {set i 0} {$i < [set $path.keynum]} {incr i} {
	    if {[string compare [set $path.keyarray($i)] "$index"] == 0} {
		set cnt $i
		break
	    }
	}
    }

    set cnt [lindex [set $path.keyarray($cnt)] 0]

    if {$option == "-text"} {
	return [eval $path.h.l$cnt cget -text]
    } elseif {$option == "-background"} {
        return [eval $path.h.f$cnt cget -background]
    } elseif {$option == "-enter"} {
	set enter "[bind $path.h.f$cnt <Any-Enter>]"
	return [lrange $enter 2 end] 
    } elseif {$option == "-leave"} {
	set leave "[bind $path.h.f$cnt <B2-Leave>]"
	return [lrange $leave 2 end]
    } elseif {$option == "-drop"} {
	set drop "[bind $path.h.f$cnt <<move-release>>]"
	return [lrange $drop 4 end]
    } elseif {$option == "-motion"} {
	set motion "[bind $path.h.f$cnt <<move-motion>>]"
	return [lrange $motion 4 end]
    } else {
	lappend arglist $option $i
    }
    #puts "invalid option"
}

#return the number of keys 
proc keybox_size {path} {
    global $path.keynum

    return [set $path.keynum]
}
