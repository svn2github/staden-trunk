#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
##############################################################################
#returns the length of a contig
proc c_length { io contig_num } {
    if {$contig_num <= 0} {
       	stack_dump
       	return 0
    }

    set c [io_read_contig $io $contig_num]
    return [keylget c length]

}
##############################################################################
#returns the length of a reading
proc r_length { io r_num } {

    set r [io_read_reading $io $r_num]
    return [keylget r sequence_length]

}
##############################################################################
#returns the name of the left most gel of a contig given a contig number
proc left_gel { io contig_num } {
    #set c [io_read_contig $io $contig_num]
    #set r_num [keylget c left]

    # Use chain_left now as this doesn't do any disk IO.
    set r_num [db_info chain_left $io =$contig_num]
    return [io_read_reading_name $io $r_num]
}

##############################################################################
#returns the reading name from the reading number
proc r_name {io r_num} {
    return [io_read_reading_name $io $r_num]
}

##############################################################################
#returns the template name from the template number
proc t_name {io t_num} {
    set t [io_read_template $io $t_num]
    return [io_read_text $io [keylget t name]]
}

##############################################################################
#initialises the global contig variables for a new database
proc InitContigGlobals { io } {
    global CurContig
    global LREG
    global RREG
    global NGELS
 
    #set LREG and RREG and CurContig

    set long_contig [db_info longest_contig $io]
    if {$long_contig == 0} {
	# No data
	return
    }
    set name [left_gel $io $long_contig]    
    set length [c_length $io $long_contig]

    set LREG 1
    set RREG $length
    set CurContig "$name"
    set db [io_read_database $io]
    set NGELS [keylget db num_readings]
}

##############################################################################
#set the contig global variables 
proc SetContigGlobals { io gel_id args} {
    global CurContig
    global LREG
    global RREG

    set gel_num [db_info get_read_num $io $gel_id]

    #convert gel_number to gel_name
    set gel_name "[r_name $io $gel_num]"

    set CurContig $gel_name
    if {"$args" != ""} {	
        set LREG [lindex $args 0]
        set RREG [lindex $args 1]
    } else {
	# Make sure we reset LREG and RREG even when in simple dialogues
	ContigParams $io
    }
}


# 
# Parses a tag list (in a suitable format for the add_tags command) and 
# complements the tag positions only for the tags on complemented sequences. 
# (The reason is that tags on complemented sequences are stored in their 
# original positions, but in some applications it's not easy to take this 
# into account until the tag list has already been constructed). 
# 
proc complement_tag_list {io tags} { 
    set new_tags "" 
    foreach tag $tags { 
	# Parse tag 
	scan $tag "%d %s %s %d..%d %n" num type dir start end pos 
	set com [string range $tag $pos end] 
 
	if {$num > 0} { 
	    set r [io_read_reading $io $num] 
	    if {[keylget r sense] == 0} { 
		lappend new_tags $tag 
	    } else { 
		# Complement start/end postions 
		set len [keylget r length] 
		set new_start [expr {$len-$end+1}] 
		set new_end   [expr {$len-$start+1}] 
		lappend new_tags \
		    "$num $type $dir $new_start..$new_end $com" 
	    } 
	} else { 
	    lappend new_tags $tag 
	} 
    } 
 
    return $new_tags 
} 
