# UNUSED

#set skip_list "Test Utils"
set skip_list ""

proc read_menu_file {file} {
    set mlist ""
    set depth 0

    set fd [open $file r]
    while {[gets $fd line] != -1} {
	if {$line == ""} {
	    continue
	}
	regexp {^(.)	*(.*)} $line dummy mode prog
	switch -- $mode {
	    + {
		append mlist "\{ [list "M {$prog}"] \{ "
	    }
	    - {
		append mlist "\} \} "
	    }
	    = {
		append mlist "[list [list C $prog]] "
	    }
	    default {
		puts "Badly formatted line '$line'"
	    }
	}
    }
    close $fd
    return $mlist
}

proc list_menu {mlist {depth ""}} {
    foreach i $mlist {
	if {[lindex [lindex $i 0] 0] == "M"} {
	    puts $depth+[lindex [lindex $i 0] 1]
	    list_menu [lindex $i 1] "$depth	"
	} else {
	    puts $depth[lindex $i 1]
	}
    }
}

proc create_menu {mlist w {skip_list {}}} {
    set count 0
    set mlist [lsort -dictionary $mlist]
    foreach i $mlist {
	if {[lindex [lindex $i 0] 0] == "M"} {
	    if {[lsearch -exact $skip_list [lindex [lindex $i 0] 1]] != -1} {
		continue
	    }
	    menu [set m $w.m$count]
	    incr count
	    $w add cascade -menu $m -label [lindex [lindex $i 0] 1]
	    create_menu [lindex $i 1] $m $skip_list
	} else {
	    $w add command -label [lindex $i 1]
	}
    }
}

if {$argc != 1} {
    puts stderr "Usage: wish disp_menu.tcl menu_file_name"
    exit 1
}

set mlist [read_menu_file [lindex $argv 0]]
catch {tkinit}

# Horizontal menu
pack [text .t -width 90] -side bottom -fill both -expand 1
menu .menu
. configure -menu .menu
create_menu $mlist .menu $skip_list

# A single EMBOSS menu with cascading items
toplevel .t2 -menu .t2.m
menu .t2.m
menu .t2.m.w
.t2.m add cascade -label Emboss -menu .t2.m.w
create_menu $mlist .t2.m.w $skip_list

#exit

