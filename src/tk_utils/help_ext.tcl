#!/usr/local/bin/wish4.0
#
# Copyright (c) 1995 Medical Research Council, Laboratory of Molecular Biology.
# All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#
# Library to communicate with the help viewer. This should be the only file
# with internal knowledge of the help code; and it should be minimal knowledge
# at that.
#
# The "show_help" function requests for help to be viewed from a given file
# and topic. The help is viewed using a standalone help viewer (help.tcl).
#

source $env(STADTABL)/help_config

proc find_help_interp {} {
    foreach i [winfo interps] {
	if [string match help.tcl* $i] {
	    return $i
	}
    }

    return ""
}

proc show_help {file {topic Contents}} {
    if {[set i [find_help_interp]]==""} {
	global auto_index

	regsub "source " $auto_index(help_init) "" z
	exec stash -f $z -init &
	exec sleep 1
	return [show_help $file $topic]
    }

    set url [show_help_common $file $topic]
    if {$url != ""} {
	send -async $i "help_display $url"
    }
}
