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
# and topic. The help is viewed using a window in the current tcl application.
#

source $env(STADTABL)/help_config

proc show_url {url} {
    if {![winfo exists .help]} {
	help_init ".help"
    }

    if {$url != ""} {
	help_display $url
    }
}

proc show_help {file {topic Contents}} {
    set url [show_help_common $file $topic]
    show_url $url
}
