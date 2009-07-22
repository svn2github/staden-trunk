#!/bin/sh
#\
exec stash "$0" ${@+"$@"} || exit 1

#-----------------------------------------------------------------------------
# Main entry point

# Package loading
if {[info exists env(STADEN_AUTO_PATH)]} {
    set auto_path "$env(STADEN_AUTO_PATH) $auto_path"
}

if {[catch tkinit]} {
    package require Tk
}
wm title . "Prefinish 1.0b1"
source $env(STADTABL)/shlib.conf
load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}

load_package tk_utils
tk_utils_init
package require Iwidgets
namespace import itcl::*
load_package prefinish

# Do the actual business
::prefinish::init
set h [prefinish::main_gui ""]

while {$argc > 0} {
    set arg [lindex $argv 0];
    set argv [lrange $argv 1 $argc]
    incr argc -1;

    switch -exact -- $arg {
	"-config" {
	    set arg [lindex $argv 0]
	    set argv [lrange $argv 1 $argc]
	    incr argc -1;

	    ::prefinish::prefin_load $h $arg
	}
    }
}
#-----------------------------------------------------------------------------

