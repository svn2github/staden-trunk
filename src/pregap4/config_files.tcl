#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1998. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#
# Config creation part of pregap.
#
# This handles writing to UNIX files, for now it just prints to stdout.
#

#
# Reads a config file into memory. It is sourced, so that the values and
# procedures are added to this interpreter.
# It is also stored in the conf_section array and conf_sections list.
#
proc read_conf_file {file} {
    global conf_file conf_section conf_sections conf_format
    global ns_regexp fofn_dir

    set conf_file $file
    set section "user_start"
    set conf_section($section) ""
    set conf_sections user_start
    set conf_format 1.0

    if {![file exists $conf_file]} {
	return
    }

    
    catch {unset ns_regexp}

    set fd [open $conf_file r]
    set all_file [read $fd]
    foreach line [split $all_file \n] {
	if {[regsub {^# Config file format version (.*)$} $line {\1} \
		c]} {
	    if {$c != $conf_format} {
		puts "Unknown config file format \"$c\""
	    }
	} elseif {[regsub {^\[(.*)\]$} $line {\1} ns]} {
	    if {$ns == "user_start"} {
		continue
	    }	
	    if {[string match "::*" $section]} {
		namespace eval $section $conf_section($section)
	    } else {
	       	uplevel #0 eval [list $conf_section($section)]
	    }
	    set section $ns
	    lappend conf_sections $section
	    set conf_section($section) ""
        } else {
	    append conf_section($section) $line\n
       	}
    }
    if {$section != "user_end"} {
	set conf_section(user_end) ""
	lappend conf_sections user_end
    }
    if {[string match "::*" $section]} {
	namespace eval $section $conf_section($section)
    } else {
       	uplevel #0 eval [list $conf_section($section)]
    }

    set_name_scheme

    close $fd

    # puts_conf_file
}


#
# Debugging function
#
proc puts_conf_file {} {
    global conf_file conf_section conf_sections

    foreach section $conf_sections {
	if {$conf_section($section) != ""} {
	    puts "********** Section '$section' **********"
	    puts $conf_section($section)
	}
    }
}

#
# Writes a config file. With no arguments it'll be written back to the same
# file loaded with read_conf_file. Otherwise the argument gives a new filename
# to write the file to.
#
proc write_conf_file {} {
    global conf_file conf_section conf_sections conf_format fofn_dir

    set cf $conf_file

    set fd [open $cf w]
    puts $fd "# Config file format version $conf_format"
    foreach section $conf_sections {
	if {$conf_section($section) != ""} {
	    puts $fd "\[$section\]"
	    puts $fd $conf_section($section)
	}
    }
    close $fd

    # puts_conf_file
}

#
# Returns an named section of the loaded configuration file
#
proc get_conf_data_section {section} {
    global conf_section

    if {[info exists conf_section($section)]} {
        return $conf_section($section)
    } else {
	return ""
    }
}

#
# Sets an named section of the loaded configuration file
#
proc set_conf_data_section {section data} {
    global conf_section conf_sections

    if {[lsearch $conf_sections $section] == -1} {
	lappend conf_sections $section
    }
    set conf_section($section) $data
}

#
# Sets an individual parameter of a module listed in the loaded configuration.
# The section name for a module is "::$module". This keeps track of which
# parameters have already been set, and rewrites them if appropriate.
#
proc set_conf_data_module_param {module var val} {
    # Get old configuration
    set c [get_conf_data_section ::$module]

    # Produce new config, replacing 'var val' if found, appending if not
    set new_conf ""
    set done 0
    foreach line [split $c \n] {
	if {$line == ""} {
	    continue
	}
	if {[string match "set $var *" $line]} {
	     append new_conf "set $var [list $val]\n"
	     set done 1
	} else {
	    append new_conf $line\n
	}
    }
    if {!$done} {
	append new_conf "set $var [list $val]\n"
    }

    # Rewrite new config
    set_conf_data_section ::$module $new_conf
}

#
# Sets a global variable
#
proc set_conf_data_globals {var val} {
    # Get old configuration
    set c [get_conf_data_section global_variables]

    # Produce new config, replacing 'var val' if found, appending if not
    set new_conf ""
    set done 0
    foreach line [split $c \n] {
	if {$line == ""} {
	    continue
	}
	if {[string match "set $var *" $line]} {
	     append new_conf "set $var [list $val]\n"
	     set done 1
	} else {
	    append new_conf $line\n
	}
    }
    if {!$done} {
	append new_conf "set $var [list $val]\n"
    }

    # Rewrite new config
    set_conf_data_section global_variables $new_conf
}
