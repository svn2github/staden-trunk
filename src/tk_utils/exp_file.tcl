#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1998. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

# Reads an experiment file into an array structure. Returns the array
# in [array get] format
proc read_exp_file {file} {
    # Init
    array set lines {}
    set line_count 1

    # Open and read file
    set fd [open $file r]
    set last_type ""
    while {[gets $fd line] != -1} {
	regsub {^(..).*} $line {\1} line_type
	if {[regsub {^..        (.*)} $line {\1} line_val] &&
	    $line_type == $last_type} {
	    set append 1
	} else {
	    regsub {^..   (.*)} $line {\1} line_val
	    set append 0
	}
	if {$line_type == "SQ"} {
	    set line_val ""
	    while {[gets $fd line] != -1 && $line != "//"} {
		regsub -all "\[ \n\t\]*" $line {} line
		append line_val $line
	    }
	    set append 0
	}

	if {$append} {
	    set last [lindex $lines($line_type) end]
	    set last2 [lindex $last 1]
	    set lines($line_type) [lreplace $lines($line_type) \
		end end [list $line_count "$last2\n$line_val"]]
	} else {
	    if {[info exists lines($line_type)]} {
		lappend lines($line_type) [list $line_count $line_val]
	    } else {
		set lines($line_type) [list [list $line_count $line_val]]
	    }
	}
	incr line_count

	set last_type $line_type
    }
    set lines(__) $line_count

    close $fd
    array get lines
}

# Write an array (passed by reference) containing an experiment file to
# the file named 'fname'
proc write_exp_file {arr_name fname} {
    upvar $arr_name lines

    set fd [open $fname w]

    array set lcount {}
    foreach l [array names lines] {
	foreach i $lines($l) {
	    set lcount([lindex $i 0]) [list $l [lindex $i 1]]
	}
    }

    #set highest $lines(__)

    set highest [count_lines lines]

    for {set i 1} {$i <= $highest} {incr i} {
	if {![info exists lcount($i)]} { continue }
	set type [lindex $lcount($i) 0]
	if {$type == "SQ"} {
	    puts $fd "$type   "
	    set seq [lindex $lcount($i) 1]
	    set len [string length $seq]
	    set count 0
	    for {set j 0} {$j < $len} {incr j 10} {
		if {$count == 0} {
		    puts -nonewline $fd "     "
		} else {
		    puts -nonewline $fd " "
		}
		puts -nonewline $fd [string range $seq $j [expr $j+9]]
		if {[incr count] == 6} {
		    set count 0
		    puts $fd ""
		}
	    }
	    if {$count != 0} {
		puts $fd ""
	    }
	    puts $fd "//"
	} else {
	    set first 1
	    foreach item [split [lindex $lcount($i) 1] \n] {
		if {$first} {
		    puts $fd "$type   $item"
		    set first 0
		} else {
		    puts $fd "$type        $item"
	 	}
	    }
	}
    }

    close $fd
}

proc count_lines {arr_name} {
    upvar $arr_name lines

    set line_count 0
    if {[info exists lines(__)]} {
	set line_count $lines(__)
    } else {
	foreach l [array names lines] {
	    foreach i $lines($l) {
		incr line_count
	    }
	}
    }
    return $line_count
}

# Adds a line to an experiment file
proc add_to_exp_file {arr_name line_type line_val} {
    upvar $arr_name lines

    set line_count [count_lines lines]

    if {[info exists lines($line_type)]} {
	lappend lines($line_type) [list $line_count $line_val]
    } else {
	set lines($line_type) [list [list $line_count $line_val]]
    }

    incr line_count
    set lines(__) $line_count
}

# Queries a line from an experiment file.
# In general it returns the last value. For the case of TG, TC, CC, EX and PS
# lines it returns a list
proc query_exp_file {arr_name line_type} {
    upvar $arr_name lines

    if {![info exists lines($line_type)]} {
	return ""
    }

    set l $lines($line_type)
    if {$line_type == "TG" ||
	$line_type == "TC" ||
	$line_type == "CC" ||
	$line_type == "EX" ||
	$line_type == "PS"} {
	set items {}
	foreach i $l {
	    lappend items [lindex $i 1]
	}
	return $items
    } else {
	return [lindex [lindex $l end] 1]
    }
}

#array set l [read_exp_file a.exp]
#query_exp_file l SQ
#query_exp_file l ID
#query_exp_file l TG
#write_exp_file l z
