#!/bin/sh
#\
exec stash -f "$0" ${@+"$@"} || exit 1

catch {console hide}

#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1998. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

# Make pregap4 take priority over tk_utils. This is temporary until we fold
# back in some of our local changes to tk_utils.
set auto_path "$env(STADTCL)/pregap4 $auto_path"
if [info exists env(STADEN_AUTO_PATH)] {
    set auto_path "$env(STADEN_AUTO_PATH) $auto_path"
}

#-----------------------------------------------------------------------------
# Parse command line arguments
set interactive 1
if {$tcl_platform(platform) == "windows"} {
    set conf_file config.pg4
} else {
    set conf_file pregap4.config
}
set conf_data ""
set fofn "pregap"
set fofn_dir [pwd]
set orig_dir [pwd]
set files ""

while {[string match {-*} [lindex $argv 0]]} {
    if {[lindex $argv 0] == "-no_win" ||
	[lindex $argv 0] == "-nowin"} {
	set interactive 0
	set argv [lrange $argv 1 end]
    } elseif {[lindex $argv 0] == "-config"} {
	set conf_file [lindex $argv 1]
	set argv [lrange $argv 2 end]
    } elseif {[lindex $argv 0] == "-fofn"} {
	set fofn [lindex $argv 1]
        set fd [open $fofn r]
	if {[set dir "[file dirname $fofn]/"] == "./"} {set dir ""}
	if {$dir == "[pwd]/"} {set dir ""}
        while {[gets $fd line] != -1} {
	    set line [string trim $line]
	    if {$line != ""} {
	        lappend files $dir$line
	    }
        }
        close $fd
	set argv [lrange $argv 2 end]
    } elseif {[lindex $argv 0] == "-out_dir"} {
	set fofn_dir [lindex $argv 1]
	if {[file pathtype $fofn_dir] != "absolute"} {
	    set fofn_dir [file join [pwd] $fofn_dir]
	}
	set argv [lrange $argv 2 end]
    } elseif {[lindex $argv 0] == "-win_compact"} {
    	keylset pregap4_defs WINDOW_STYLE compact
	set argv [lrange $argv 1 end]
    } elseif {[lindex $argv 0] == "-win_separate"} {
    	keylset pregap4_defs WINDOW_STYLE separate
	set argv [lrange $argv 1 end]
    } elseif {[lindex $argv 0] == "--"} {
	set argv [lrange $argv 1 end]
	break
    } else {
	puts "ERROR: unknown parameter '[lindex $argv 0]'"
	exit
    }
}

#-----------------------------------------------------------------------------
# Initialise variables
set modules ""

source $env(STADTABL)/shlib.conf
load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}

load_package tk_utils
load_package pregap4

if {$interactive} {
    if {[catch tkinit]} {
	package require Tk
    }
    tk_utils_init
    init_images
}

append files " $argv"
read_conf_file $conf_file

#-----------------------------------------------------------------------------
# Find and load pregap modules
load_modules

#-----------------------------------------------------------------------------
# Do it.
# Either go through the module steps and then exit (non interactive), or
# build the GUI and wait for the GUI commands.
if {!$interactive} {
    init_modules
    if {[check_modules 1] > 0} {
	verror ERR_WARN check_modules \
	    "One or more modules have not been configured correctly. Exiting"
	exit 1
    }
    set files [run_modules $files]
    set files [shutdown_modules $files]
    exit
} else {
    set tk_strictMotif 0
    #package require Ttk
    wm withdraw .
    wm protocol . WM_DELETE_WINDOW exit
    if {$licence(type) != "f"} {
	LicenceSplash Pregap4 2
    }
    build_gui {}
    update
    wm deiconify .
}

raise .
focus -force .
