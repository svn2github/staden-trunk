#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#


# tkinit strips off -h from argv and produces help messages. We'd rather handle
# this ourselves before it goes that far.

proc gap4_usage {} {
    verror ERR_FATAL "Gap4" {Usage: gap4 [-ro] [-maxseq <number>] [-maxdb <number>] [-skip_notes] [NAME.V]}
    verror ERR_FATAL "Gap4" {
Usage: gap4 [options] [DATABASE_NAME.V]
Options:
    -ro                Read-only
    -maxseq N          Initial maximum total consensus length
    -maxdb N           Initial maximum number of contigs + readings
    -(no_)check        (Do not) run Check Database upon DB open
    -(no_)exec_notes   (Do not) execute OPEN/CLOS notes upon database open
    -(no_)rawdata_note (Do not) use RAWD note instead of RAWDATA env. variable
    -(no_)csel         (Do not) start up the contig selector upon db open
    -display           X Display to use
    -bitsize N         Specifies the aux file bitsize for new databases (32/64)
    -sync              Use synchronous mode for display server (debugging)
    --                 End of argument list
    }
    exit
}

if {[lindex $argv 0] == "-h" || [lindex $argv 1] == "-help" || [lindex $argv 1] == "--help"} {
    gap4_usage
}

if {[catch tkinit err]} {
    if {[string match "*connect to display*" $err] == 1} {
	puts stderr "Can't open display: $env(DISPLAY)"
	exit 1
    }
    if {[catch {package require Tk} err]} {
	if {[string match "*connect to display*" $err] == 1} {
	    puts stderr "Can't open display: $env(DISPLAY)"
	} else {
	    puts stderr $err
	}
	exit 1
    }
}
catch {console hide}
wm withdraw .

proc splash_screen {} {
    toplevel [set w .splash]
    image create photo splash -file /nfs/team71/psg/jkb/logo.gif
    wm transient .splash
    wm overrideredirect .splash 1
    label $w.splash -image splash
    pack $w.splash -fill both
    set x [expr {([winfo screenwidth $w]-[image width splash])/2}]
    set y [expr {([winfo screenheight $w]-[image height splash])/2}]
    wm geometry $w +$x+$y
}

#splash_screen
#update idletasks

# Override the current X defaults - these are not obvious, and ORDER COUNTS!
# Probably only necessary for CDE/Solaris. Best to override with
# tk_setPalette anyway.
#option add *Background #b03060
#option add *foreground #c3c3c3
#option add *Foreground black
#option add *background #d9d9d9
#option add *troughColor #c3c3c3
#option add *activeForeground black
#option add *activeBackground #ececec
#option add *indicator #b03060
#option add *disabled #a3a3a3
#option add *selectBackground #c3c3c3
#option add *selectForeground black

#
# Returns a string describing the database
# Also (as this is a common place) bumps at maxseq as needed.
#
proc DatabaseInfo {io} {
    global maxseq

    set db [io_read_database $io]
    set nc [keylget db num_contigs]
    set nr [keylget db num_readings]
    set tcl [db_info t_contig_length $io]
    set trl [db_info t_read_length $io]
    set maxgel [keylget db max_gel_len]

    set i ""
    append i [format "Database size       %10d       Max reading length %10d\n" \
		 [keylget db actual_db_size] $maxgel]
    append i [format "No. Readings        %10d       No. Contigs        %10d\n" \
	         $nr $nc]
    append i [format "No. Annotations     %10d       No. Templates      %10d\n" \
		 [keylget db Nannotations] [keylget db Ntemplates]]
    append i [format "No. Clones          %10d       No. Vectors        %10d\n" \
		 [keylget db Nclones]      [keylget db Nvectors]]
    append i [format "Total contig length %10ld       Average length      %11.1f\n"\
		 $tcl [expr double($tcl)/$nc]]
    append i [format "Total characters in readings                       %15ld\n" \
		 $trl]
    append i [format "Average reading characters per consensus character    %15.2f\n" \
	        [expr double($trl)/$tcl]]
    append i [format "Average used length of reading                        %15.2f\n" \
		[expr double($trl) / $nr]]

    set maxgel 1024; # minor space saving
    if {[expr {($tcl + (2 * $maxgel + 20)*$nc)*1.1}] > $maxseq} {
	set maxseq [expr {round(($tcl + (2 * $maxgel + 20)*$nc)*1.1)}]
	verror ERR_WARN gap "increasing maxseq parameter to $maxseq"
    }
    
    append i "Current maximum consensus length is $maxseq\n"

    return $i
}

proc PrintDatabaseInfo {io} {
    vfuncheader "Database information"
    vmessage [DatabaseInfo $io]
}

##############################################################################
# Copies a database
proc CopyDatabase {io} {
    global gap_defs
    global licence

    set t [keylget gap_defs COPY_DATABASE.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Copy database"

    global $t.GCollect

    entrybox $t.entry -title "New version character" -width 4 \
	-command "global $t.GCollect; CopyDatabase2 $t \[set $t.GCollect\] $io"

    if {$licence(type) == "f"} {
	checkbutton $t.collect -text "Garbage collect" \
		-variable $t.GCollect
	set $t.GCollect [keylget gap_defs COPY_DATABASE.COLLECT]
    } else {
	set $t.GCollect 0
    }

    okcancelhelp $t.ok -bd 2 -relief groove \
	-ok_command "CopyDatabase2 $t \[set $t.GCollect\] $io \
		     \[entrybox_get $t.entry\]" \
	-cancel_command "unset $t.GCollect;destroy $t" \
	-help_command "show_help gap4 {GapDB-CopyDatabase}"

    if {$licence(type) == "f"} {
	pack $t.entry $t.collect $t.ok -side top -fill both
    } else {
	pack $t.entry $t.ok -side top -fill both
    }
}

##############################################################################
proc CopyDatabase2 {t collect io version} {
    global $t.GCollect

    set version [string index $version 0]
    if {[regexp {^[a-zA-Z0-9~+-]$} $version] == 0} {
	bell
	tk_messageBox -icon error -type ok -title "Error" \
		-message "Version character should be alpha-numeric or one of\
		~, + or -"
	return
    }

    destroy $t
    unset $t.GCollect
    
    if {[set dot [string last . [set dbn [db_info db_name $io]]]] == -1} {
	bell
	return
    }

    if {[CheckSaveFile "[string range $dbn 0 [expr $dot]]$version"] == 1} {
	copy_db -io $io -version $version -collect $collect
    }
}

##############################################################################
proc ActivateMenu_Open { } {
    menu_state_set gap_menu 12 .mainwin.menus
}

##############################################################################
#sets the menu item states after New has been chsen
proc ActivateMenu_New { } {
    menu_state_set gap_menu 4 .mainwin.menus
    menu_state_set gap_menu -8 .mainwin.menus
}

##############################################################################
#checks for read-only mode and disasbles menus accordingly
proc Menu_Check_RO { } {
    global read_only

    if {$read_only} {
        menu_state_set gap_menu -16 .mainwin.menus
    }
}

##############################################################################
#sets the menus to the state where we've closed the database
proc DisableMenu_Open { } {
    menu_state_set gap_menu -12 .mainwin.menus
}

##############################################################################
#gap examine quality option
proc PlotQuality { io } {
    global gap_defs

    #plot_quality -io $io
    QualityPlotDialog [keylget gap_defs QUALITY_PLOT.WIN] $io

}

##############################################################################
#gap find internal joins option
proc FindInternalJoins { io } {
    global gap_defs

    #puts "start tcl FindInternalJoins"
    FIJDialog [keylget gap_defs FIJ.WIN] $io
}

##############################################################################
#gap find read pairs option
proc FindReadPairs { io } {
    global NGRec
    global gap_defs

    ReadPairDialog $io [keylget gap_defs READPAIR.WIN]
}
##############################################################################
#gap extract readings option
proc ExtractReadings { io} {
    global gap_defs

    #puts "start tcl ExtractReadings"
    ExtractReadingsDialog $io [keylget gap_defs EXTRACT.WIN]
}
##############################################################################
#gap auto assemble option
proc PreAssemble { io } {
    global gap_defs

    PreAssembleDialog $io [keylget gap_defs PREASS.WIN]

    # Not needed now (hopefully)
    # SaveContigOrder $io

}


##############################################################################
#gap auto assemble option
proc AutoAssemble { io option } {
    global gap_defs

    if {$option == 1} {
	NormalShotgun [keylget gap_defs AUTO_ASSEMBLE.1.WIN] $io
    } elseif {$option == 2} {
	Screen [keylget gap_defs AUTO_ASSEMBLE.2.WIN] $io
    } elseif {$option == 3} {
	OneContig [keylget gap_defs AUTO_ASSEMBLE.3.WIN] $io
    } elseif {$option == 4} {
	NewContig [keylget gap_defs AUTO_ASSEMBLE.4.WIN] $io
    } elseif {$option == 5} {
	SSAssemble [keylget gap_defs AUTO_ASSEMBLE.5.WIN] $io
    } elseif {$option == 6} {
	IgnorePrev [keylget gap_defs AUTO_ASSEMBLE.6.WIN] $io
    }
}

##############################################################################
#gap disassemble readings option
proc DisReadings { io } {
    global gap_defs

    DisReadingsDialog $io [keylget gap_defs DIS_READINGS.WIN] [keylget gap_defs CONTIG_SEL.WIN]
}

##############################################################################
#gap find repeats option
proc FindRepeats { io} {
    global gap_defs

    FindRepeatsDialog $io [keylget gap_defs FINDREP.WIN]
}

##############################################################################
#gap show relationships
proc ShowRelationships { io } {
    global gap_defs

    ShowRelationshipsDialog $io

}
##############################################################################
#gap suggest long gels
proc LongGels { io } {
    global gap_defs

    LongGelsDialog $io [keylget gap_defs LONGGELS.WIN]
}
##############################################################################
#gap suggest taq terminators
proc TaqTerminator { io } {
    global gap_defs 
    
    TaqTerminatorDialog $io [keylget gap_defs TTERM.WIN]
}

##############################################################################
#write the current order of contigs in the contig selector to the database
proc SaveContigOrder { io } {

    flush_contig_order -io $io

}

##############################################################################
#automatic contig ordering based on readpair information
proc OrderContigs { io } {

    order_contigs -io $io

}

##############################################################################
proc ChangeDirectory { } {
    set t .change_dir
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Change directory"

    entrybox $t.entry -title "Change directory" -width 25 \
	-default [pwd] \
	-command "destroy $t; ChangeDirectory2"

    okcancelhelp $t.ok -bd 2 -relief groove \
	-ok_command "ChangeDirectory2 \[expandpath \[entrybox_get $t.entry\]\]; destroy $t" \
	-cancel_command "destroy $t"\
	-help_command "show_help gap4 {GapDB-Directories}"

    pack $t.entry $t.ok -side top -fill both
}


proc ChangeDirectory2 { dir} {

    catch {cd $dir} e
    if {$e != ""} {
	tk_messageBox -icon error -type ok -title "Change directory" \
		-message $e
    } else {
	vfuncheader "Change Directory"
	vmessage "Changing directory to [pwd]"
	catch {unset ::tk::dialog::file::__tk_filedialog(selectPath)}
    }
}

##############################################################################
#The exit command
proc Gap_Exit { io } {
    if {$io > 0} {
        if {![quit_displays $io "exit"]} {
	    bell
	    tk_messageBox -icon warning -type ok -title "Database in use" \
		    -message "Please shut down any editing displays before quitting"
	    return
        }
	close_db -io $io
    }
    exit
}

##############################################################################
proc Gap_Open { io } {
    #check if any display is busy
    if {[InitOpenAnotherDB] } {
	set file_type {
	    {"database"		    *.aux		    }
	}
	set file [tk_getOpenFile -filetypes $file_type]
	if {$file != ""} {
	    DBClose
	    DB_Load $file
	}
    } else {
	Gap_Exit $io
    }
}
##############################################################################
#set up main tcl menu bar
#NB: cs isn't obviously used here, but create_menus uses it (via uplevel)
#when adding the menu items.
proc CreateMain { } {
    global gap_defs gap_menu
    global GAP_VERSION

    wm protocol . WM_DELETE_WINDOW {Gap_Exit $io}

    ##########################################################################
    # Main Menu Bar and initialise busy code
    pack [frame .mainwin] -expand 1 -fill both
    menu .mainwin.menus
    . configure -menu .mainwin.menus
    create_menus $gap_menu .mainwin.menus [keylget gap_defs MENU_LEVEL]
    InitBusy .mainwin .mainwin.menus gap_menu

    # Create the windows
    tout_create_wins ".mainwin"

    # Create tags and tag bindings
    gap4_text_init .mainwin.stdout.t

    set mode "???"
    foreach m [keylget gap_defs MENU_LEVELS] {
	if {[lindex $m 1] == [keylget gap_defs MENU_LEVEL]} {
	    set mode [lindex $m 0]
	}
    }
    vmessage "Gap4 has started up in '$mode' mode. To select another menu level," 
    vmessage "please use the 'Configure menus' command in the 'Options' menu."
}
#end CreateMain

# Moves window "w" to be under window "u" (where possible).
# If force is true then it will make sure that w is under u even if this
# means moving u up a bit to make some room.
#
# Returns whether or not there is room without resorting to force.
proc MoveWinUnder {w u {force 0}} {
    set fit 1

    set w [winfo toplevel $w]
    set u [winfo toplevel $u]

    # Geometry specifies where the top-left corner of the window is
    regexp {(.*)x(.*)[-+](.*)[-+](.*)} [wm geometry $u] _ xs ys x1 y1
    
    # "winfo y" returns the coordinate of the top-left corner of the window
    # as available to the application (ie minus the wm title bar and menus)
    set dy [expr {[winfo y $u]-$y1}]
    set dx [expr {[winfo x $u]-$x1}]
    
    set x2 [expr {$x1+$xs+$dx}]
    set y2 [expr {$y1+$ys+$dy}]

    if {[expr {[winfo height $w]+$y2}] > [winfo vrootheight $w]} {
	set fit 0
    }

    if {$fit || $force} {
	wm geometry $w +$x1+$y2
    } elseif {$force} {
	# calculate y pos of $w if it is put at the bottom of the screen
	set ypos [expr {[winfo vrootheight $w]-[winfo height $w]}]
	# Now work out the new ypos of $u
	set y1 [expr {$ypos-$ys-$dy}]
	wm geometry $u +$x1+$y1
	return MoveWinUnder $w $u 1
    }

    return $fit
}


##############################################################################
#                                  main program                              #
##############################################################################
#lappend auto_path $env(TCLUTILS_LIBRARY)
if {[info exists env(STADEN_AUTO_PATH)]} {
    set auto_path "$env(STADEN_AUTO_PATH) $auto_path"
}

source $env(STADTABL)/shlib.conf
load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}
load_package tk_utils
tk_utils_init

if {$argc >= 2 && [lindex $argv 0] == "-menu_file"} {
    keylset gap_defs MENU_FILE [lindex $argv 1]
    set argv [lrange $argv 2 $argc]
    incr argc -2
}

load_package gap

if {[keylget gap_defs BACKGROUNDS] != ""} {
    if {[tk appname Gap4] != "Gap4"} {
        set app [tk appname]
        # Multiple gap4s - recolour!
        set gnum 2
        foreach col [keylget gap_defs BACKGROUNDS] {
	    set "colours(Gap4 #$gnum)" $col
	    incr gnum
        }
        keylset tk_utils_defs FOREGROUND black
        if {[info exists colours($app)]} {
	    keylset tk_utils_defs BACKGROUND $colours($app)
        } else {
	    keylset tk_utils_defs BACKGROUND $col
        }
	unset col
	unset gnum
	unset colours
    }
}

if {$licence(type) == "d"} {
    LicenceSplash Gap4 3
}

global maxseq
global maxdb
global db_namelen
set db_namelen 40; #Also see IO.h and legacy.f
set io -1
set create 0
set access "rw"
set read_only 0
set do_check_db [keylget gap_defs CHECKDB_AT_STARTUP]
set do_csel     [keylget gap_defs CONTIG_SEL.DISPLAY_AT_STARTUP]
set logging [keylget gap_defs LOGGING]
set exec_notes 0
set rawdata_note 1

set GAP_VERSION "4.11.2$svn_version"

switch $licence(type) {
    f		{}
    v 		{append GAP_VERSION " (VIEWER)"; set read_only 1; set access r}
    default	{append GAP_VERSION " (DEMO)"}
}

wm title . "[tk appname] v$GAP_VERSION"
wm iconname . "[tk appname] v$GAP_VERSION"

while {$argc > 0 && "[string index [lindex $argv 0] 0]" == "-"} {
    set arg [lindex $argv 0];

    if {$arg == "--"} {
        set argv [lrange $argv 1 $argc]
        incr argc -1
	break;

    } elseif {$arg == "-ro" || $arg == "-read_only"} {
	set access "r"

    } elseif {$arg == "-maxseq"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    set maxseq [lindex $argv 0]
	} else {
	    gap4_usage
	}

    } elseif {$arg == "-maxdb"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    set maxdb [lindex $argv 0]
	} else {
	    gap4_usage
	}

    } elseif {$arg == "-nocheck" || $arg == "-no_check"} {
	set do_check_db 0

    } elseif {$arg == "-check"} {
	set do_check_db 1

    } elseif {$arg == "-no_exec_notes"} {
	set exec_notes 0

    } elseif {$arg == "-exec_notes"} {
	set exec_notes 1

    } elseif {$arg == "-no_csel"} {
	set do_csel 0

    } elseif {$arg == "-csel"} {
	set do_csel 1

    } elseif {$arg == "-no_rawdata_note"} {
	set rawdata_note 0

    } elseif {$arg == "-rawdata_note"} {
	set rawdata_note 1

    } elseif {$arg == "-bitsize"} {
	if {$argc > 1} {
            set argv [lrange $argv 1 $argc]
            incr argc -1
	    set bitsize [lindex $argv 0]
	    set_db_bitsize $bitsize
	} else {
	    gap4_usage
	}

    } elseif {$arg == "-h" || $arg == "-help" || $arg == "--help"} {
	gap4_usage
    } else {
	verror ERR_WARN "Gap4" "ERROR: Invalid argument \"$arg\""
	gap4_usage
    }

    set argv [lrange $argv 1 $argc]
    incr argc -1
}

if {$argc == 1} {
#    scan [lindex $argv 0] "%s.%d" db_name version_num
    set a [lindex $argv 0]
    if {[regexp {(.*)\.(.)(\.aux)?$} $a tmp db_name version_num] == 0} {
	verror ERR_FATAL "Gap4" "ERROR: Invalid database name '$a'"
	exit
    }

    cd [file dirname $db_name]
    set db_name [file tail $db_name]

    #open database 
    set origtype $licence(type)
    set io [open_db -name $db_name -version $version_num -access $access \
	    -create $create] 

    if {"$io" == ""} {
#	puts "ERROR: Couldn't open the database '$db_name.$version_num' - exiting."
#	6/1/99 johnt - use verror so message is shown with Windows NT
	if {$licence(type) == "d"} {
	    verror ERR_FATAL "Gap4" "ERROR: Demonstration mode could not open this database. Try \"gap4viewer\" instead."
	} else {	
	    verror ERR_FATAL "Gap4" "ERROR: Couldn't open the database '$db_name.$version_num' - exiting."
	}
	
	exit
    }

    if {$origtype == "d" && $licence(type) == "v"} {
	viewer_mode
    }

    if {$read_only == 1 && $access == "rw" && $licence(type) != "v"} {
	set ret [tk_messageBox -icon warning -type okcancel \
		-title "Database open" \
		-message "Database is busy. Opened in read-only mode"]
	if {$ret == "cancel"} {
	    exit
	}
    }

} elseif {$argc > 1} {
    gap4_usage
}

set consensus_mode          [keylget gap_defs CONSENSUS_MODE]
set consensus_cutoff        [keylget gap_defs CONSENSUS_CUTOFF]
set quality_cutoff          [keylget gap_defs QUALITY_CUTOFF]
set chem_as_double          [keylget gap_defs CHEM_AS_DOUBLE]
set consensus_iub           [keylget gap_defs CONSENSUS_IUB]
set template_size_tolerance [keylget gap_defs TEMPLATE_TOLERANCE]
set min_vector_len          [keylget gap_defs MIN_VECTOR_LENGTH]
set align_open_cost         [keylget gap_defs ALIGNMENT.OPEN.COST]
set align_extend_cost       [keylget gap_defs ALIGNMENT.EXTEND.COST]
load_alignment_matrix       [keylget gap_defs ALIGNMENT.MATRIX_FILE]
set ignore_all_ptype        [keylget gap_defs IGNORE_ALL_PTYPE]
set ignore_custom_ptype     [keylget gap_defs IGNORE_CUSTOM_PTYPE]
if {[catch {load_genetic_code -filename [keylget gap_defs GENETIC_CODE]} err]} {
    verror ERR_WARNING load_genetic_code $err
}
set gap_fatal_errors 1

#initialise C and tcl tag id and tag type arrays
InitTagArray
InitCanvasConstants
InitLists
InitListTrace

#display top level menu
CreateMain

#display contig selector and menus if opened database on command line
if {$io > 0} {
    if {$read_only} {
	set extras "     *READ-ONLY*"
    } else {
	set extras ""
    }
    wm title . "GAP v$GAP_VERSION: [db_info db_name $io]$extras"
    wm iconname . "GAP: [db_info db_name $io]"
    if {[db_info num_contigs $io] > 0} {
	if {$do_check_db == 1 || ($do_check_db == -1 && !$read_only)} {
	    if {[check_database -io $io] == 0} {
		if {$do_csel} {
		    ContigSelector $io
		    MoveWinUnder . [keylget gap_defs CONTIG_SEL.WIN]
		}
	    }
	} else {
	    if {$do_csel} {
		ContigSelector $io
		MoveWinUnder . [keylget gap_defs CONTIG_SEL.WIN]
	    }
	}
	ActivateMenu_Open
	InitContigGlobals $io
	update idletasks
	PrintDatabaseInfo $io
    } else {
	ActivateMenu_New
    }
}

Menu_Check_RO
update
wm deiconify .

# Some handy bindings

#tkEntryBind dummy_event
#bind Entry <Delete> [bind Entry <BackSpace>]

#tkTextBind dummy_event
#bind Text <Delete> [bind Text <BackSpace>]

# A sort of tk_focusFollowsMouse, but only within a window (not between them)
bind allfocus <Any-Enter> {
    catch {
	if {[winfo toplevel [focus]] == [winfo toplevel %W]} {
	    focus %W
	}
    }
}

# Force auto-loading of the EventHandler
catch EventHandler

#destroy .splash
raise .
focus -force .
 
