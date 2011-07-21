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

    if {$io > 0} {
	close_db -io $io
	set io 0
    }

    DisableMenu_Open
}

proc InitOpenAnotherDB {} {
    global gap_defs io

    if {$io > 0 && ![quit_displays $io "open_database"]} {
	# Someone's too busy to shutdown?
	return 0
    }
    update idletasks

    return 1
}

##############################################################################
#open a new database
proc NewFile {} {
    global gap_defs io
    
    # Shut down any existing displays.
    if {$io > 0} {
	if {![quit_displays $io "create new database"]} {
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
                                 -filetypes {{database *.aux}} \
                  		 -title {New Database} \
				 -parent $f \
		    \]"
    pack $f.files.browse -side right -fill both
    

    ###########################################################################
    #OK and Cancel buttons
     okcancelhelp $f.ok_cancel \
	-ok_command "New_OK_Pressed $f \[expandpath \[entrybox_get $f.files.entry\]\]" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap4 {GapDB-New}" \
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
    global gap_defs
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
    regsub {\.aux$} $file {} file
    if {[regexp {(.*)\.(.)$} $file dummy db_name version] == 0} {
	set db_name $file
	set version 0
    }
    if {[regexp {\.} $db_name]} {
	tk_messageBox \
		-icon error \
		-title "Bad filename" \
		-message "The use of \".\" in the database name is not permitted" \
		-type ok \
	        -parent $f
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
    catch {file delete -force $db_name.$version}
    catch {file delete -force $db_name.$version.aux}
    catch {file delete -force $db_name.$version.log}
    catch {file delete -force $db_name.$version.BUSY}
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
    global gap_defs
    global io
    global GAP_VERSION
    global licence

    if {$licence(type) == "v"} {
        set access "r"
    } else {
        set access "rw"
    }

    set create 0

    #strip off .aux if is exists
    regsub \.aux$ $file "" filename

    set response [CheckOpenFile $filename]
    if {$response == 0} {
	    bell
	    #wait until user has re-entered a value
	    tkwait variable re_enter
	    return ""
    }

    set version_num [file extension $filename]
    set version_num [string range $version_num 1 end]
    set db_name [file rootname $filename]

    #FIXME - how to reset Busy if there is an error in opening the db
    global read_only

    #remove any remaining displays
    if {![InitOpenAnotherDB] } {
	Gap_Exit $io
    }
    DBClose

    set orig_type $licence(type)
    set new_io [open_db -name $db_name -version $version_num -access $access \
	    -create $create] 

    if {$new_io == ""} {
	if {$licence(type) == "d"} {
	    tk_messageBox \
		    -icon error \
		    -title "Invalid database" \
		    -message "Demonstration mode may not open this database.\
			      Try \"gap4viewer\" for a read-only viewer instead." \
		    -type ok
	} else {
	    tk_messageBox \
		    -icon error \
		    -title "Invalid database" \
		    -message "This database could not be opened." \
		    -type ok
	}
	return
    }

    # If licence type has changed from demo to viewer, rebuild the menus
    if {$orig_type == "d" && $licence(type) == "v"} {
	viewer_mode
    }

    if {$read_only == 1 && $access == "rw" && $licence(type) != "v"} {
	set ret [tk_messageBox \
		-icon warning \
		-title "Database open" \
		-message "Database is busy. Opened in read-only mode." \
		-type okcancel]
	if {$ret == "cancel"} {
	    close_db -io $new_io
	    return
	}
    }

    set io $new_io

    if {$read_only} {
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
	if {[winfo exists [keylget gap_defs CONTIG_SEL.WIN]]} {
	    destroy [keylget gap_defs CONTIG_SEL.WIN]
	}
	ActivateMenu_New
    } else {
	InitContigGlobals $io

	#draw the contig selector if it does not already exist
	if {[check_database -io $io] == 0} {
	    set cs_win [keylget gap_defs CONTIG_SEL.WIN]
	    global do_csel
	    if {![winfo exists $cs_win] && $do_csel} {
		ContigSelector $io
	    } else {
		ContigInitReg $io
		catch {raise $cs_win}
	    }
	}

	PrintDatabaseInfo $io

	ActivateMenu_Open 
    }
    Menu_Check_RO
    InitLists
    #UpdateListContigs $io [keylget gap_defs CONTIG_SEL.WIN]
}

# Enables viewer mode from demo mode. This regenerates the menus and changes
# the title bar.
proc viewer_mode {} {
    global gap_defs env GAP_VERSION read_only

    tk_messageBox \
	-icon info \
	-title "Gap4 Viewer" \
	-message "This database may not be used in demonstration mode.\
		  Gap4 will now switch to a \"view-only\" mode with\
		  restricted functionality." \
	-type ok

    # Force read-only mode
    set read_only 1
    
    # Udpate title bar
    regsub {\(DEMO\)} $GAP_VERSION {(VIEWER)} GAP_VERSION
    wm title . "GAP v$GAP_VERSION"
    wm iconname . "GAP v$GAP_VERSION"

    # Clear old menus
    foreach mname {gap select_notes edit_note contig_editor_editmodes \
		       contig_editor_commands contig_editor_settings \
		       contig_editor_help template selector consistency} {
	global ${mname}_menu
	set ${mname}_menu ""
    }

    # Reload the menus specs
    set_defs gap
    load_menus

    if {![winfo exists .mainwin.menus]} {
	return
    }

    # Destroy current menu widgets
    set mpath .mainwin.menus
    $mpath delete 0 end
    foreach child [winfo children $mpath] {
	destroy $child
    }

    # Recreate the menus
    create_menus $gap_menu $mpath [keylget gap_defs MENU_LEVEL]
    reset_menu_state gap_menu .mainwin.menus

    foreach child [winfo children .] {
	if {[winfo toplevel $child] == $child} {
	    destroy $child
	}
    }
}