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
# Library to communicate with netscape as a help viewer.
#
# The "show_help" function requests for help to be viewed from a given file
# and topic.
#

source $env(STADTABL)/help_config

proc show_url {url} {
    global tk_utils_defs
    set prog [keylget tk_utils_defs HELP_PROGRAM]

    set cmd [format {[ "`%s -remote 'OpenUrl(%s,new-tab)' 2>&1`" != "" ] \
	    && %s '%s'} $prog $url $prog $url]

    exec sh -c "$cmd" &
}

proc show_help {file {topic Contents}} {
    set url [show_help_common $file $topic]
    if {$url != ""} {
	show_url $url
    }
}

