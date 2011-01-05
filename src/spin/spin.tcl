#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

if {[catch tkinit]} {
    package require Tk
}
catch {console hide}

##############################################################################
#sets the menu item states after Get sequences
proc SeqActivateMenu_Open { } {
    menu_state_set seq_menu 4 .mainwin.menus
}

##############################################################################
#deactivate menus when no sequences
proc SeqActivateMenu_Initial { } {
    menu_state_set seq_menu -12 .mainwin.menus
    menu_state_set seq_menu 1 .mainwin.menus
}

##############################################################################
#sequence manager displays list of available sequences and which ones are
#currently horizontal and vertical
proc SequenceManager { } {
    
    sequence_list_create
    sequence_list_update
}

##############################################################################
#simple sequence browser interface
proc CreateSimpleBrowser {} {
    global spin_defs
    
    set f [keylget spin_defs SIMPLE.WIN]

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Load sequence"

    xget_fname $f.entry -type load_multiple

    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $f.button -bd 2 -relief groove \
	    -ok_command "SimpleBrowser2 $f $f.entry"\
	    -cancel_command "destroy $f" \
	    -help_command "show_help spin {SPIN-Read Sequences}"

    pack $f.button -side bottom -fill x
    pack $f.entry
}

#
# Returns a temporary filename
#

##############################################################################
proc SimpleBrowser2 { f entry } {
    global seq_menu 

    if {[set sname_list [$entry get]] == ""} {
	bell
	return
    }

    foreach sname $sname_list {
	set seq_id [read_sequence -file $sname -name $sname]
	if {$seq_id == -1} {
	    bell
	}
    }
    
    sequence_list_update
    SeqActivateMenu_Open

    destroy $f
}

proc ChangeDir {} {
    set dir [tk_chooseDirectory -initialdir [pwd]]
    if {$dir == ""} {
	#user pressed cancel
	return
    }
    if {[catch {cd $dir}]} {
	verror ERR_WARN ChangeDir "Couldn't change directory to '$dir'"
    } else {
	vmessage "Changed directory to '$dir'"
    }
}

##############################################################################

# Not used any more
proc CreateLibBrowser { } {
    global seqlib_defs

    set f [keylget seqlib_defs LIBRARY.WIN]
    global $f.load_lib

    if {[xtoplevel $f] == ""} return
    wm title $f "Library browser"
    SeqLibraries $f
    trace variable $f.load_lib w "LoadSequences"
}

##############################################################################
#callback from seq_lib browser - invoked via trace variable $f.load_lib
proc LoadSequences {name element op } {
    upvar $name result
    set seq_id [read_sequence -library [lindex $result 0] \
		    -entry_mode [lindex $result 1] \
		    -entry [lindex $result 2]\
		    -direction 0]

    if {($seq_id != -1)} {
	SeqActivateMenu_Open
	sequence_list_update
	return
    }
}

##############################################################################
proc CreatePersonalBrowser {} {
    global spin_defs

    set f [keylget spin_defs PERSONAL.WIN]
    if {[xtoplevel $f] == ""} return
    wm title $f "Personal file browser"

    frame $f.top
    entrybox $f.fn -title Filename -width 15 -command "GetPersonalFiles $f.fn $f.lists"

    button $f.browse -text Browse -command "InvokeFileBrowser $f.fn open; GetPersonalFiles $f.fn $f.lists \[entrybox_get $f.fn\]"

    listbox $f.lists -xscrollcommand "$f.scrollx set" -yscrollcommand "$f.scrolly set" -selectmode extended
    scrollbar $f.scrolly -command "$f.lists yview" -orient vertical
    scrollbar $f.scrollx -command "$f.lists xview" -orient horizontal

    okcancelhelp $f.button -bd 2 -relief groove \
	-ok_command "LoadPersonalSequences $f.fn $f.lists; destroy $f"\
	-cancel_command "destroy $f" \
	-help_command "show_help spin {SPIN-Personal search}"
    
    bind [entrybox_path $f.fn] <Any-Leave> \
	"GetPersonalFiles $f.fn $f.lists \[entrybox_get $f.fn\]"

    grid rowconfig $f 1 -weight 1
    grid columnconfig $f 0 -weight 1

    grid $f.top -row 0 -column 0 -sticky w
    grid $f.fn -in $f.top -row 0 -column 0 -sticky w
    grid $f.browse -in $f.top -row 0 -column 1 -sticky w
    grid $f.lists -row 1 -column 0 -sticky nesw -columnspan 1
    grid $f.scrolly -row 1 -column 1 -sticky ns
    grid $f.scrollx -row 2 -column 0 -sticky ew
    grid $f.button -row 3 -columnspan 2 -sticky ew

}   

proc GetPersonalFiles {fn lb filename} {

    set list [get_archive_list -file $filename]
    $lb delete 0 end
    foreach i $list {
	$lb insert end $i
    }
}

##############################################################################
proc LoadPersonalSequences {fn lb} {

    set filename [entrybox_get $fn]
    set files ""
    set file [$lb curselection]
    foreach f $file {
	lappend files [$lb get $f]
    }

    set seq_id [read_sequence -entry $files -file $filename  -direction 0]

    if {($seq_id != -1)} {
	SeqActivateMenu_Open
	sequence_list_update
	return
    }
}

##############################################################################
proc SeqCreateMainWin {f } {
    global seq_menu sip_defs rem_dup
    
    wm protocol . WM_DELETE_WINDOW {StartSeqShutdown}

    menu .mainwin.menus
    . configure -menu .mainwin.menus
    create_menus $seq_menu .mainwin.menus
    InitBusy .mainwin .mainwin.menus seq_menu

    #set up global rem_dup (remove duplicates) checkbutton variable
    set rem_dup [keylget sip_defs SIP.DUP.VALUE]
    set_remove_dup $rem_dup

    tout_create_wins ".mainwin" 90
    tout_config_popup {"Show input parameters" Print "Output to disk" "Output to command"}
}

##############################################################################
proc StartSeqShutdown { } {
    seq_quit_displays
    exit
}

proc seq_shutdown {seq_id} {
    global $seq_id.start $seq_id.end

    if {[info exists $seq_id.start]} {
	unset $seq_id.start
    }    
    if {[info exists $seq_id.end]} {
	unset $seq_id.end
    }
}

proc seq_exit {} {
    exit
}

##############################################################################
#                                  main program                              #
##############################################################################
if {[info exists env(STADEN_AUTO_PATH)]} {
    set auto_path "$env(STADEN_AUTO_PATH) $auto_path"
}

source $env(STADTABL)/shlib.conf
load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}

load_package tk_utils
tk_utils_init
load_package spin


if {$licence(type) == "f"} {
   set SPIN_VERSION 1.3
} else {
   set SPIN_VERSION "1.3 DEMO"
}
wm title . "SPIN v$SPIN_VERSION"
wm iconname . "SPIN v$SPIN_VERSION"
wm withdraw .

if {$licence(type) != "f"} {
    LicenceSplash Spin 1
}

set DNA 1
set PROTEIN 2
set HORIZONTAL 0
set VERTICAL 1


#set up sequence library array in C
#if {[keylget seqlib_defs SEQLIB.LIB_TYPE] != ""}  {
#    find_seq_libs
#}

#set replot temporary results to be false
set_replot_temp 0

#necessary for busy cursor
set w .mainwin
pack [frame $w] -expand yes -fill both
SeqCreateMainWin $w[keylget spin_defs SEQ.WIN]

#set up score matrices
set_score_matrix -file [keylget sip_defs SIP.PROT_MAT] -type $PROTEIN
set_score_matrix -file "<identity>" -type $DNA
set_max_matches [keylget sip_defs SIP.MAX_MATCHES]
set_def_matches [keylget sip_defs SIP.NUM_MATCHES]

#set up Raster class bindings
set x_format %d
set y_format %6f
RasterBindings $x_format $y_format
RasterGlobals

# Autoload EventHandler
catch EventHandler
tk appname Spin
update
wm deiconify .

set direction 0
foreach file $argv {
    #only allow personal filenames to be read in on the command line
    set seq_id [read_sequence -file $file -direction $direction]
    if {$seq_id != -1} {
	set direction [expr ($direction+1)%2]
	sequence_list_update
	SeqActivateMenu_Open
    }
}

raise .
focus -force .

if {$tcl_platform(os) == "Darwin"} {
    # EMBOSS needs to write files, so make sure we start in a temporary
    # directory somewhere if we can't write to our current directory.
    if {[catch {set fd [open .___ w]}]} {
	catch {cd $env(HOME)}
	vmessage "Changed directory to '$env(HOME)'"
    } else {
	close $fd
	file delete .___
    }
}
