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
    global env
    global env tcl_platform


    # Windows file associations do not like urls such as
    # "document.html#section" - we need to remove the #section bit.
    regsub {#[^#]*$} $url {} url



    # The start programs are different on windows 9x and windows NT/2000
    # On NT/2000 we must insert quotes as the first argument (window title).
    # On win9x there is no such option.
    if { $tcl_platform(os) == "Windows NT" } {
	if [ regexp -nocase -- {COMMAND\.COM$} $env(COMSPEC) ] {
	    exec start "" $url &
	} else {
	    exec $env(COMSPEC) /c start "" $url &
        }
    } else {
	if [ regexp -nocase -- {COMMAND\.COM$} $env(COMSPEC) ] {
	    exec start $url &
	} else {
	    exec $env(COMSPEC) /c start $url &
        }
    }
}

proc show_help {file {topic Contents}} {
    set url [show_help_common $file $topic]
    if {$url != ""} {
	show_url $url
    }
}

