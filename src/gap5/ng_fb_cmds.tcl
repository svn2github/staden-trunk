#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.


# Clean up old windows that are hanging around when closing the database.
# This prevents accidental re-use of old $io handles, causing gap5 to crash

;proc CleanUpOldWindows { } {
    set kids [winfo children .]
    array set toplevels {}
    foreach win $kids {
	set toplevels([winfo toplevel $win]) 1
    }
    foreach top [array names toplevels] {
	if { $top ne "." } { destroy $top }
    }
}

proc DBClose {} {
    global io

    CleanUpOldWindows

    if {$io != ""} {
	$io close
	set io ""
    }
    
    DisableMenu_Open
}

proc InitOpenAnotherDB {} {
    global gap5_defs io

    if {$io != "" && ![quit_displays -io $io -msg "open_database"]} {
	# Someone's too busy to shutdown?
	return 0
    }
    CleanUpOldWindows
    update idletasks

    return 1
}

##############################################################################
#open a new database
proc NewFile {} {
    global gap5_defs io

    # Shut down any existing displays.
    if {$io != ""} {
	if {![quit_displays -io $io -msg "create new database"]} {
	    # Someone's too busy to shutdown?
	    return
	}
	CleanUpOldWindows
    }

    set f .new
    global $f.browsed_path
    set $f.browsed_path ""

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "New database"
    frame $f.files
    xentry $f.files.entry \
	-label "Enter new database name" \
	-width 15

    pack $f.files.entry -side left -expand 1 -fill both

    button $f.files.browse \
	-text "Browse" \
	-command "NewFileBrowser_invoke $f.files.entry $f.browsed_path"

    pack $f.files.browse -side right -fill both
    

    ###########################################################################
    #OK and Cancel buttons
     okcancelhelp $f.ok_cancel \
	-ok_command "New_OK_Pressed $f \[$f.files.entry get\] $f.browsed_path" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 {GapDB-New}" \
	-bd 2 \
	-relief groove
   

    pack $f.files
    pack $f.ok_cancel -fill x
}

proc NewFileBrowser_invoke { entry browsed_path } {
    global $browsed_path

    set types {
	{{Gap5 databases} {.g5d}}
    }
    set file [tk_getSaveFile -filetypes $types \
		  -defaultextension .g5d -parent $entry]
        set name $file
    if { $file != "" } {
	if [string match [pwd]* $file] {
	    set name [string range $file [expr [string length [pwd]]+1] end]
	}
	if [file exists $file] {
	    # If file exists, tk_getSaveFile will have asked if the user
	    # wants to replace it, so we don't want to again.
	    set $browsed_path $name
	} else {
	    set $browsed_path ""
	}
    }
    #if user has pressed cancel don't want to remove name from xentry
    if {$name != ""} {
	$entry delete 0 end
	$entry insert 0 $name
	$entry xview [string last / $name]
    }
}

##############################################################################
proc New_OK_Pressed { f entry_content browsed_path } {
    if {[New_Open_Database $f $entry_content $browsed_path] == ""} {
	return
    }

    PostLoadSetup
    destroy $f
}

proc New_Open_Database { f entry_content browsed_path } {
    global gap5_defs
    global io
    global GAP_VERSION
    global licence
    global $browsed_path

    if {$licence(type) == "v"} {
        set access "r"
    } else {
        set access "rw"
    }
    set create 1
    set file [expandpath $entry_content]
    set ask_overwrite 1
    if { $entry_content eq [set $browsed_path] } {
	set ask_overwrite 0
    }

    if {$file == ""} {
	bell
	return
    }

    regsub {\.(g5d|g5x)$} $file {} file
    # Ask if the user wants to overwrite an existing database, but
    # only if the save dialogue hasn't already.
    if { $ask_overwrite } {
	set overwrite [Overwrite "$file.g5d" $f]
	if { 0 == $overwrite || 2 == $overwrite } {
	    return
	}
    }
    
    # Parse $file to extract directory, db_name
    set dir [file dirname $file]
    set file [file tail $file]

    # if {[regexp {(.*)\.(.)$} $file dummy db_name version] == 0} {
    # 	set db_name $file
    # 	set version 0
    # }
    # if {[regexp {\.} $db_name]} {
    # 	tk_messageBox \
    # 		-icon error \
    # 		-title "Bad filename" \
    # 		-message "The use of \".\" in the database name is not permitted" \
    # 		-type ok \
    # 	        -parent $f
    # 	return
    # }

    #if a database if already open then close it
    InitOpenAnotherDB
    DBClose

    #change directory
    set cur_dir [pwd]
    cd $dir
    if {$dir != "." && $dir != $cur_dir} {
	vmessage "Changing directory to [pwd]"
    }

    # Do it...
    catch {file delete -force $file.g5d}
    catch {file delete -force $file.g5x}
    catch {file delete -force $file.g5d.BUSY}
    #catch {file delete -force $db_name.$version.log}
    #catch {file delete -force $db_name.$version.BUSY}
    set io [g5::open_database -name $file -access $access -create $create]

    return $io
}

##############################################################################
#open file and load in database
proc DB_Load { file } {
    global gap5_defs
    global io
    global licence

    if {$licence(type) == "v"} {
        set access "r"
    } else {
        set access "rw"
    }

    set response [CheckOpenFile $file]
    if {$response == 0} {
	    bell
	    #wait until user has re-entered a value
	    tkwait variable re_enter
	    return ""
    }

    #strip off .aux if is exists
    regsub {\.(aux|g5d|g5x)$} $file "" filename

    set db_name $filename

    #remove any remaining displays
    if {![InitOpenAnotherDB] } {
	Gap_Exit $io
    }
    DBClose

    set orig_type $licence(type)
    set new_io [g5::open_database -name "$db_name" -access $access]
    global debug_level
    $new_io debug_level $debug_level

    if {$new_io == ""} {
	tk_messageBox \
	    -icon error \
	    -title "Invalid database" \
	    -message "This database could not be opened." \
	    -type ok
	return
    }

    if {[$new_io read_only] == 1 && $access == "rw"} {
	set ret [tk_messageBox \
		-icon warning \
		-title "Database open" \
		-message "Database is busy. Opened in read-only mode." \
		-type okcancel]
	if {$ret == "cancel"} {
	    $new_io close
	    return
	}
    }

    set io $new_io

    #change directory
    set dir [file dirname $file]
    set cur_dir [pwd]
    cd $dir

    if {$dir != "." && $dir != $cur_dir} {
	vmessage "Changing directory to [pwd]"
    }
    
    PostLoadSetup
}

proc PostLoadSetup { } {
    # Post-load set-up of menus, lists, contig globals and the contig selector.
    global io
    global gap5_defs
    global GAP_VERSION

    # Set window title
    if {[$io read_only]} {
	set extras "     *READ-ONLY*"
    } else {
	set extras ""
    }
    wm title . "GAP v$GAP_VERSION: [db_info db_name $io]$extras"
    wm iconname . "GAP: [db_info db_name $io]"

    #check to see if the database is empty. If it is, open as if used the
    #'New' option and remove contig selector if it exists
    #else set up contig globals, contig_selector and all menus
    if {[db_info num_contigs $io] == 0} {
	DisableMenu_Open
	if {[winfo exists [keylget gap5_defs CONTIG_SEL.WIN]]} {
	    destroy [keylget gap5_defs CONTIG_SEL.WIN]
	}
	ActivateMenu_New
    } else {
	InitContigGlobals $io

	#draw the contig selector if it does not already exist
	if {[check_database -io $io] == 0} {
	    set cs_win [keylget gap5_defs CONTIG_SEL.WIN]
	    global do_csel
	    if {![winfo exists $cs_win] && $do_csel} {
		if {$do_csel == 2 || [db_info num_contigs $io] <= 1000} {
		    ContigSelector $io
		} else {
		    vmessage "\nSkipping contig selector due to large number of contigs."
		}
	    } else {
		ContigInitReg $io
		catch {raise $cs_win}
	    }
	}

	#PrintDatabaseInfo $io

	ActivateMenu_Open 
    }
    Menu_Check_RO $io
    InitLists
    #UpdateListContigs $io [keylget gap5_defs CONTIG_SEL.WIN]    
}
