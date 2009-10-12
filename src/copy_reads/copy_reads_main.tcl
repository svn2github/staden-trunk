#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#The exit command
proc Copy_Reads_Exit { } {
    exit
}

##############################################################################
#set up main tcl menu bar
#NB: cs isn't obviously used here, but create_menus uses it (via uplevel)
#when adding the menu items.
proc CopyReadsCreateMainWin { } {
    global copy_reads_menu

    wm protocol . WM_DELETE_WINDOW {Copy_Reads_Exit $io_from $io_to}

    ##########################################################################
    # Main Menu Bar and initialise busy code
    pack [frame .mainwin] -expand 1 -fill both
    menu .mainwin.menus
    . configure -menu .mainwin.menus
    create_menus $copy_reads_menu .mainwin.menus
    InitBusy .mainwin .mainwin.menus copy_reads_menu

    # Create the windows
    tout_create_wins ".mainwin"
    # Create tags and tag bindings
    #FIXME
    #gap4_text_init .mainwin.stdout.t

}
#end CreateMain


##############################################################################
#                                  main program                              #
##############################################################################
if {[info exists env(STADEN_AUTO_PATH)]} {
    set auto_path "$env(STADEN_AUTO_PATH) $auto_path"
}

if {$argc >= 2 && [lindex $argv 0] == "-menu_file"} {
    keylset gap_defs MENU_FILE [lindex $argv 1]
    set argv [lrange $argv 2 $argc]
    incr argc -2
}


if {$tcl_platform(os) == "Darwin"} {
    load $env(STADLIB)/libtk_utils.dylib
}

source $env(STADTABL)/shlib.conf
load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}
load_package tk_utils
tk_utils_init

load_package gap
load_package copy_reads

if {$licence(type) != "f"} {
    verror ERR_FATAL "copy reads" "ERROR: Unable to copy reads in Demonstation or Viewer mode"
    exit
}

set COPY_READS_VERSION 1.0

set win 0
set list_from ""
set list_to ""
set source_trace_dir ""
set args ""

proc copy_reads_usage {} {
    verror ERR_FATAL "Copy reads" {Usage: 
    copy_reads [-win] 
               [-source_trace_dir ("")]
               [-contigs_from <file> (all contigs)] 
               [-min_contig_len (2000)] 
               [-min_average_qual (30.0)] 
               [-contigs_to <file> (all contigs)] 
               [-mask <none mask> (none)] 
               [-tag_types <list> ("")] 
               [-word_length (8)] 
               [-min_overlap (20)] 
               [-max_pmismatch (30.0)] 
               [-min_match (20)] 
               [-band (1)] 
               [-display_cons] 
               [-align_max_mism (10.0)] 
               [-display_seq] 
               [source database] 
               [destination database]}
    exit
}

while {$argc > 0 && "[string index [lindex $argv 0] 0]" == "-"} {
    set arg [lindex $argv 0];

    if {$arg == "--"} {
        set argv [lrange $argv 1 $argc]
        incr argc -1
	break;

    } elseif {$arg == "-w" || $arg == "-win"} {
	set win 1
	package require Tk
	tk_utils_init
	wm title . "COPY_READS v$COPY_READS_VERSION"
	wm iconname . "COPY_READS v$COPY_READS_VERSION"
	tk appname Copy_reads

    } elseif {$arg == "-contigs_from"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    set list_from [lindex $argv 0]
	} else {
	    copy_reads_usage
	}
    } elseif {$arg == "-contigs_to"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    set list_to [lindex $argv 0]
	} else {
	    copy_reads_usage
	}
     } elseif {$arg == "-source_trace_dir"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    set source_trace_dir [lindex $argv 0]
	} else {
	    copy_reads_usage
	}
     } elseif {$arg == "-mask"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    lappend args $arg [lindex $argv 0]
	} else {
	    copy_reads_usage
	}
     } elseif {$arg == "-min_overlap"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    lappend args $arg [lindex $argv 0]
	} else {
	    copy_reads_usage
	}
     } elseif {$arg == "-max_pmismatch"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    lappend args $arg [lindex $argv 0]
	} else {
	    copy_reads_usage
	}
     } elseif {$arg == "-word_length"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    lappend args $arg [lindex $argv 0]
	} else {
	    copy_reads_usage
	}
     } elseif {$arg == "-min_match"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    lappend args $arg [lindex $argv 0]
	} else {
	    copy_reads_usage
	}
     } elseif {$arg == "-band"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    set band [lindex $argv 0]
	} else {
	    copy_reads_usage
	}
     } elseif {$arg == "-tag_types"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    lappend args $arg [lindex $argv 0]
	} else {
	    copy_reads_usage
	}
     } elseif {$arg == "-align_max_mism"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    lappend args $arg [lindex $argv 0]
	} else {
	    copy_reads_usage
	}
     } elseif {$arg == "-min_contig_len"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    lappend args $arg [lindex $argv 0]
	} else {
	    copy_reads_usage
	}

     } elseif {$arg == "-min_average_qual"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    lappend args $arg [lindex $argv 0]
	} else {
	    copy_reads_usage
	}
     } elseif {$arg == "-display_cons"} {
	 lappend args $arg 1
	 set display_cons 1
     } elseif {$arg == "-display_seq"} {
	 lappend args $arg 1
     } else {
	 verror ERR_WARN "Copy reads" "ERROR: Invalid argument \"$arg\""
	 copy_reads_usage
    }

    set argv [lrange $argv 1 $argc]
    incr argc -1
}

#allow the user to input 2 databases on the command line
if {$argc > 2} {
    copy_reads_usage
}

InitTagArray

if {$win} {
    #display top level menu
    CopyReadsCreateMainWin 
} else {
    if {$argc != 2} {
	verror ERR_FATAL "copy reads" "ERROR: need 2 databases"
	copy_reads_usage
    } else { 
	do_copy_reads [lindex $argv 0] [lindex $argv 1] $list_from $list_to $source_trace_dir $args
    }
    exit
}
