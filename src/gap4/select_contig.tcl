#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
##############################################################################
#check if reading is a name or a number
#return 1 for name
#return 0 for number
proc CheckReadingID { reading } {

    set letter [string index $reading 0]
    if { ([string compare $letter # ] == 0)} {
	#number
	return 0
    } else {
	#name
	return 1
    }

}

##############################################################################
#check the gel name exists in the database
#return 0 for success
#return 1 for failure
proc CheckGelName { io name {p .} } {

    set num [db_info get_contig_num $io $name]
    if { $num == -1 } {
	tk_messageBox \
		-icon error \
		-title "Bad name" \
		-message "Reading name does not exist in the database" \
		-type ok \
	        -parent $p
	set press_ok 1
    } else {
	return 0
    }
    tkwait variable press_ok
    return 1

}

##############################################################################
#check for valid database filenames.
#return 0 for success
#return 1 for failure
proc CheckDBFilename {filename {p .}} {
    set fname [file tail $filename]
    if {[regexp "\[ \t\n\]" $fname] != 0} {
	tk_messageBox \
		-icon error \
		-title "Bad filename" \
		-message "The use of whitespace in the database name is not permitted" \
		-type ok \
	        -parent $p
        set press_ok 1
    } elseif {[string length $fname] > 256-9} {
	tk_messageBox \
		-icon error \
		-title "Bad filename" \
		-message "The database name must be less than or equal to 247 characters long" \
		-type ok \
	        -parent $p
	set press_ok 1
    } else {
        return 1
    }
    tkwait variable press_ok
    return 0

}

##############################################################################
#check the reading name exists in the database
#return 0 for success
#return 1 for failure
proc CheckEntryName { io entry name c_len c_num } {

    upvar $c_len len $c_num num
    set press_ok 0

    set num [db_info get_contig_num $io $name]
 
    if { $num == -1 } {
	update
	tk_messageBox \
		-icon error \
		-title "Bad entry" \
		-message "Reading name does not exist in the database" \
		-type ok \
		-parent $entry
	# The focus below has been removed as it causes oddities with
	# multiple error messages. It's an interaction between the odd
	# behaviour of tk_dialog (it uses widthdraw, update, etc) and
	# the leave/focus-out bindings on the contig_id. Be warned that
	# this current solution works, but it's tricky to understand.
	# focus $entry
	$entry icursor end
	return -1
    } else {
	set c [io_read_contig $io $num]
	set len [keylget c length]
	return 0
    }
    #tkwait variable never_return
}

##############################################################################
#if the user alters the contig name then update the start & end of contig
proc UpdateContigLimits { io start end entry} {
    set name [$entry get]
    set contig_len 0
    set contig_num 0

    if { [CheckEntryName $io $entry $name contig_len contig_num] == 0 } {
	if {[$end cget -to] == [$end get]} {
	    set atend 1
	} else {
	    set atend 0
	}
	$start configure -to $contig_len
	$end configure -to $contig_len
	if {[$end get] > $contig_len || $atend} {
	    $end set $contig_len
	}
    }
}

##############################################################################
#if the user has not pressed enter after changing the contig identifier, then
#must check that the current selected values are correct
proc CheckEntry { io path name rreg } {

    set contig_len 0
    set contig_num 0

    set checked_ok 0

    if {[CheckEntryName $io $path.entry $name contig_len contig_num] == 0 } { 
	if { $rreg > $contig_len } {
	    tk_messageBox \
		    -icon error \
		    -title "Bad contig limits" \
		    -message "You have chosen an end \
		    point greater than the length of the contig. Try \
		    pressing return after entering the contig identifier \
		    to update the contig limits" \
		    -type ok \
		    -parent $path
	    set checked_ok 1
	} else {
	    #input is all OK
	    return 0
	}
	tkwait variable checked_ok
	return 1
    }
    return 1
}
