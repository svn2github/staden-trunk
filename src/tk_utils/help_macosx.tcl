#
# Copyright (c) 1995 Medical Research Council, Laboratory of Molecular Biology.
# All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#
# The "show_help" function requests for help to be viewed from a given file
# and topic.
#

source $env(STADTABL)/help_config

proc show_url {url} {
    global env

    # MacOS file associations do not like urls such as
    # "document.html#section" - we need to remove the #section bit.
    regsub {#[^#]*$} $url {} url
    regsub {^file:/+} $url {/} url

    # The mac open command can fail if we are in fullscreen X windows mode
    # and we try to bring up a help browser that cannot deal with X (eg
    # internet explorer). 
    if {[catch {exec open $url}]} {
        source $env(STADTCL)/tk_utils/help_int.tcl
	return [show_url $url]
    }
}

proc show_help {file {topic Contents}} {
    set url [show_help_common $file $topic]
    if {$url != ""} {
	show_url $url
    }
}

