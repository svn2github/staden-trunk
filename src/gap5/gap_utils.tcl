#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

##############################################################################
#Arbitrary one-time counter
set .counter 0
proc counter {} {
    global .counter
    incr .counter
    return [set .counter]
}

##############################################################################
#returns the length of a contig
proc c_length { io contig_num } {
    if {$contig_num <= 0} {
       	stack_dump
       	return 0
    }

    set c [$io get_contig $contig_num]
    set len [$c get_length]
    $c delete
    return $len
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
    return [r_name $io $r_num]
}

##############################################################################
#returns the reading name from the reading number
proc r_name {io r_num} {
    set s [$io get_seq $r_num]
    set name [$s get_name]
    $s delete
    return $name
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
    set c [$io get_contig $long_contig]
    set name [$c get_name]
    set length [$c get_length]

    set LREG [$c get_visible_start]
    set RREG [$c get_visible_end]

    set CurContig "$name"
    set NGELS 1
}

##############################################################################
#set the contig global variables 
proc SetContigGlobals { io gel_id args} {
    global CurContig
    global LREG
    global RREG

    set CurContig $gel_id
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


# Renames a contig
#
# If $w is specified it will assume this is being called from within a
# dialogue box and will generate appropriate parented errors. It will also
# destroy the parent on success.
#
# If you do not wish the parent to be destroyed, skip the $w parameter.
#
# The final "auto" parameter is used to determine whether contig names should
# be disambiguated automatically.
#
# Returns new name on success
#         "" on failure
proc contig_rename {io crec name {w {}} {auto 0}} {
    if {[$io contig_name2rec $name] > 0} {
	if {$auto} {
	    set cid 1
	    while {[$io contig_name2rec $name#$cid] > 0} {
		incr cid
	    }
	    append name #$cid
	} else {
	    bell
	    if {$w != ""} {
		tk_messageBox -type ok -icon error -parent $w \
		    -message "Contig name already exists"
	    } else {
		verror ERR_WARN "rename_contig" "Contig name already exists"
	    }
	    return ""
	}
    }

    set c [$io get_contig $crec]
    $c set_name $name
    $c delete

    if {$w != ""} {
	destroy $w
    }

    contig_notify \
	-io $io \
	-type RENAME \
	-cnum $crec \
	-args [list name $name]

    if {[$io base] != $io} {
	return $name
    }

    $io flush

    return $name
}

# Adds a contig to a named scaffold
proc contig_add_to_scaffold {io crec scaf_name {w {}}} {
    set c [$io get_contig $crec]
    if {[$c get_scaffold]} {
	# Remove from old scaffold
	$c remove_from_scaffold
    }
    if {$scaf_name != ""} {
	# Add to new scaffold
	$c add_to_scaffold $scaf_name 0 0 0
    }
    $c delete

    $io flush

    if {$w != ""} {
	destroy $w
    }
}