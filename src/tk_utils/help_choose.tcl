# The guts of the show_help command, common to all implementations.
# This does the topic/index searching and returns the filename
#
# There are two styles of help system. The old method uses double indirection
# via searching for $topic in the $file.topic file. This returns a heading
# name which is then searched for in the $file.index file.
#
# The new mode works from texinfo node names. It doens't have double
# indirection. It searches only in ${file}.index.
#
# For backwards compatibility with third party modules (eg Pride), we
# look for the presence of a .topic file. If we find one then we'll
# use the old method.
#

proc show_help_common {file {topic Contents}} {
    global help_path package_dir

    if {$file == ""} {
	return ""
    }

    if {$file == "url"} {
	return $topic
    }

    set dir $help_path
	
    if {[string range $file 0 0] == "%"} {
	set file [string range $file 1 end]
	set dir $package_dir($file)
    }
    set file $dir/$file

    # This bit's inefficient on large indexes (eg ~10000 entries). (2.8 secs)
    # Could try building an awk script and executing it.
    # eg:
    #
    # awk '$1 == "{contig_register_ini}" {print $2}'
    # (1.9 for awk, 0.4 for gawk)
    #
    # or:
    #
    # sed -n 's/^{contig_register_init}[    ]*\(.*\)/\1/p'
    # (0.9 secs)
    #
    # or:
    # egrep '^\{contig_register_ini\}' < registration.index
    # (1.5 secs)
    #
    # However, we're typically going to have indexes well under a tenth of
    # those tested above.

    if {[file exists $file.topic]} {
	# Old method

	# Convert from topic to section
	set fd [open ${file}.topic]
	while {[gets $fd l] != -1} {
	    if {[lindex $l 0] == $topic} {
		close $fd
		
		set section [lindex $l 1]
		# Convert from section to URL
		set fd [open ${file}.index]
		while {[gets $fd l] != -1} {
		    if {[lindex $l 0] == $section} {
			close $fd
			return file:$dir/[lindex $l 1]
		    }
		}
		break
	    }
	}
    } else {
	# New method

	# Convert from topic to section
	set fd [open ${file}.index]
	while {[gets $fd l] != -1} {
	    if {[string compare [lindex $l 0] $topic] == 0} {
		close $fd
		return file:$dir/[lindex $l 1]
	    }
	}
    }

    bell
    puts "Couldn't find help for subject '$topic'"
    close $fd
    
    return ""
}

proc load_help_system {dir} {
    global tk_utils_defs

    set version [keylget tk_utils_defs HELP_PROGRAM]

    if {$version == "tcl-external"} {
        uplevel #0 "source {$dir/help_ext.tcl}"
    } elseif {$version == "netscape" || \
	      $version == "mozilla" || \
	      $version == "konqueror"} {
        uplevel #0 "source {$dir/help_netscape.tcl}"
    } elseif {$version == "galeon"} {
	uplevel #0 "source {$dir/help_galeon.tcl}"
    } elseif {$version == "windows"} {
	# 7/1/99 johnt - added windows support
        uplevel #0 "source {$dir/help_windows.tcl}"
    } elseif {$version == "macosx"} {
        uplevel #0 "source {$dir/help_macosx.tcl}"
    } elseif {$version == "tcl-internal"} {
        uplevel #0 "source {$dir/help_int.tcl}"
    } else {
	# Auto
	global errorCode
	catch {exec mozilla -version}
	if {$errorCode == "NONE"} {
	    keylset tk_utils_defs HELP_PROGRAM mozilla
	} else {
	    catch {exec netscape -version}
	    if {$errorCode == "NONE"} {
		keylset tk_utils_defs HELP_PROGRAM netscape
	    } else {
		keylset tk_utils_defs HELP_PROGRAM tcl-internal
	    }
	}
	load_help_system $dir
    }
}

load_help_system $env(STADTCL)/tk_utils
