#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc DBClose {} {
    global io

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
    }

    set f .new
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "New database"
    frame $f.files
    entrybox $f.files.entry \
	-title "Enter new filename" \
	-width 15 \
	-type CheckDBOutput
    pack $f.files.entry -side left -expand 1 -fill both
    button $f.files.browse \
	-text "Browse" \
	-command "NewFile_browse $f.files.entry \
                    \[tk_getSaveFile \
                                 -filetypes {{database *.g5d}} \
                  		 -title {New Database} \
				 -parent $f \
		    \]"
    pack $f.files.browse -side right -fill both
    

    ###########################################################################
    #OK and Cancel buttons
     okcancelhelp $f.ok_cancel \
	-ok_command "New_OK_Pressed $f \[expandpath \[entrybox_get $f.files.entry\]\]" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 {GapDB-New}" \
	-bd 2 \
	-relief groove
   

    pack $f.files
    pack $f.ok_cancel -fill x
}

proc NewFile_browse {f fname} {
    if {$fname != ""} {
        entrybox_delete $f 0 end
	entrybox_insert $f 0 $fname
	[entrybox_path $f] xview end
    }
}


##############################################################################
proc New_OK_Pressed { f file } {
    global gap5_defs
    global io
    global GAP_VERSION
    global licence

    if {$licence(type) == "v"} {
        set access "r"
    } else {
        set access "rw"
    }
    set create 1

    if {$file == ""} {
	bell
	return
    }

    # Parse $file to extract directory, db_name and version
    set dir [file dirname $file]
    set file [file tail $file]
    regsub {\.g5d$} $file {} file
    if {[regexp {(.*)\.(.)$} $file dummy db_name version] == 0} {
	set db_name $file
	set version 0
    }
    if {[regexp {\.} $db_name]} {
	tk_messageBox \
		-icon error \
		-title "Bad filename" \
		-message "The use of \".\" in the database name is not permitted" \
		-type ok
	return
    }

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
    catch {file delete -force $db_name.$version.g5d}
    catch {file delete -force $db_name.$version.g5x}
    #catch {file delete -force $db_name.$version.log}
    #catch {file delete -force $db_name.$version.BUSY}
    set io [open_db -name $db_name -version $version -access $access \
	    -create $create]

    wm title . "GAP v$GAP_VERSION: [db_info db_name $io]"
    wm iconname . "GAP: [db_info db_name $io]"
    DisableMenu_Open
    ActivateMenu_New
    InitLists
    destroy $f
}

##############################################################################
#open file and load in database
proc DB_Load { file } {
    global gap5_defs
    global io
    global GAP_VERSION
    global licence

    if {$licence(type) == "v"} {
        set access "r"
    } else {
        set access "rw"
    }

    #strip off .aux if is exists
    regsub {\.(aux|g5d|g5x)$} $file "" filename

    set response [CheckOpenFile $filename]
    if {$response == 0} {
	    bell
	    #wait until user has re-entered a value
	    tkwait variable re_enter
	    return ""
    }

    set db_name $filename

    #remove any remaining displays
    if {![InitOpenAnotherDB] } {
	Gap_Exit $io
    }
    DBClose

    set orig_type $licence(type)
    set new_io [g5::open_database -name "$db_name" -access $access]

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

    if {[$io read_only]} {
	set extras "     *READ-ONLY*"
    } else {
	set extras ""
    }
    wm title . "GAP v$GAP_VERSION: [db_info db_name $io]$extras"
    wm iconname . "GAP: [db_info db_name $io]"

    #change directory
    set dir [file dirname $file]
    set cur_dir [pwd]
    cd $dir

    if {$dir != "." && $dir != $cur_dir} {
	vmessage "Changing directory to [pwd]"
    }
    
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
