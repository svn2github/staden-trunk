#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1998. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#
# Query module of pregap.
# We'll probably need a different version of this file for Unix fs, Sybase and
# Oracle access.
#

#
# Queries a specific field for sequence 'name'.
# Depending on the data source, we may wish to look up file_id(name) to get
# the unique sequence identifier for this name.
#
# UNIX file system version
proc pg_open {name {mode a+}} {
    global file_fd

    set file_fd($name) [open $name $mode]
}

proc pg_query {name field} {
    global file_id
    global lines

    if {[info exists lines]} {
	unset lines
    }
    set id $file_id($name)
    array set lines [array get ::SEQ_${id}::lines]
    global $field
    if {[info exists $field] && [set $field] != ""} {
	return [set $field]
    } elseif {[info commands ${field}_com] != ""} {
	if {[set c [${field}_com]] == ""} {
	    verror ERR_WARN "SEQ $id" \
		"Warning, ${field} command returned empty string"
	}
	return $c
    }
}

proc pg_update {name field value} {
    global file_id
    global file_type
    global file_fd

    set ::SEQ_$file_id($name)::lines($field) $value

    if {$file_type($name) == "EXP"} {
	set value [string trimright $value \n] 
	set lines [split $value \n] 
	set nlines [llength $lines] 
	set str "$field   [lindex $lines 0]" 
	for {set i 1} {$i < $nlines} {incr i} { 
	    append str "\n$field        [lindex $lines $i]" 
	}
	puts $file_fd($name) "$str"
    }
}

proc pg_close {name} {
    global file_id
    global file_fd

    close $file_fd($name)
}
