#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1998. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

proc set_module_state {nb mod} {
    if {[namespace eval ::$mod {info exists mandatory}]} {
	set img "   "
    } elseif {[set ${mod}::enabled]} {
	set img {[x]}
    } else {
	set img {[ ]}
    }
    
    $nb update $mod 1 $img
    $nb row_state $mod [set ${mod}::enabled]
}

proc enable_module {nb mod_num} {
    set mod [$nb number_to_name $mod_num]

    if {[namespace eval ::$mod {info exists mandatory}]} {
	return
    }

    if {[set ${mod}::enabled]} {
	set ${mod}::enabled 0
    } else {
	set ${mod}::enabled 1
    }

    set_module_state $nb $mod

    $nb raise $mod

    config_panel_col1 $nb $mod [winfo parent $nb]
}

proc enable_changed {nb mod args} {
    if {![winfo exists $nb]} { return }

    set_module_state $nb $mod
}

proc save_all_pressed {but_win mod_win mod} {
    ${mod}::configure_dialogue $mod_win save_all
    namespace eval ::$mod {mod_save enabled [set enabled]}
    write_conf_file
}

proc menu_panel {w} {
    global pregap4_menu
    global pregap4_defs
    keylset pregap4_defs MENU.WIN $w.menus

    if {$w == ""} {
	set top .
    } else {
	set top [winfo toplevel $w]
    }

    menu $w.menus
    $top configure -menu $w.menus
    create_menus $pregap4_menu $w.menus
}

proc module_panel {} {
    global pregap4_defs
    select_modules_create [keylget pregap4_defs MODULE.WIN]
}

proc text_panel {} {
    global pregap4_defs
    set w [keylget pregap4_defs OUTPUT.WIN]
    tout_create_wins $w
    tout_config_popup {Print "Output to disk" "Output to command"}
}

proc config_panel {{raise 1}} {
    global env
    global modules
    global pregap4_defs

    set w [keylget pregap4_defs CONFIG.WIN]

    if {![winfo exists $w]} {
	xtoplevel $w
	wm title $w "Configure modules"
	wm protocol $w WM_DELETE_WINDOW \
	    "if {\[ok_gui 0\]} {destroy $w}"
    }

    if {[winfo exists $w.nb]} {
	if {!$raise} { return }
	if {[keylget pregap4_defs WINDOW_STYLE] == "compact"} {
	    [keylget pregap4_defs NB.WIN] select 1
	} else {
	    raise $w
	}

	return
    }

    xlistnotebook $w.nb \
	-exportselection 0 \
	-height 10 \
	-bd 1 -relief raised \
	-column_pad 0

    $w.nb bindcol 1 <ButtonRelease-1> "enable_module $w.nb"

#Give up on this for the time being as it's hard to find working dingbats
#fonts. Some linux systems have a font named "dingbats" which is completely
#blank - hence giving "" and "" as the two tick/cross symbols!
#Instead we now use +/-.
#
#    # Find zapfdingbats - not an easy task. On windows it's called "monotype
#    # sorts", otherwise it's under various guises (itc zapf dingbats, zapf
#    # dingbats, zapfdingbats, itc zapfdingbats, etc).
#    set fam [font families]
#    if {[set ind [lsearch -glob $fam "Pregap"]] == -1} {
#	if {[set ind [lsearch -glob $fam "*Monotype Sorts*"]] == -1} {
#	    if {[set ind [lsearch -glob $fam "*dingbats*"]] == -1} {
#		set ind 0; # Not correct, but we have to have a font
#	    }
#	}
#    }
#    set font [lindex $fam $ind]
#    $w.nb columnadd -font [list $font 12 normal]

    $w.nb columnadd -font {-*-courier-bold-r-normal-*-14-*-*-*-*-*-*-*}
    $w.nb columnadd -font {-*-helvetica-bold-r-normal-*-12-*-*-*-*-*-*-*}
    $w.nb columnadd -font {-*-courier-bold-r-normal-*-12-*-*-*-*-*-*-*}

    # Get the name for each module, creating a tab for each, and then call
    # it's dialogue module
    foreach mod $modules {
	if {[catch {${mod}::name} err]} {
	    verror ERR_WARN init "Could not load module $mod: $err"
	    continue
	}
	if {![info exists ${mod}::hidden]} {
	    trace variable ${mod}::enabled w "enable_changed $w.nb $mod"
	    set f [$w.nb add $mod [list "   " "[${mod}::name]    " "      "] \
		      -raisecommand "config_notebook_raise $w.nb $mod"]
	    config_panel_col1 $w.nb $mod
	    frame $f.indent -bd 0 -relief flat
	    pack $f.indent -padx 5 -pady 5 -fill both -expand 1
	    set f $f.indent

	    set_module_state $w.nb $mod
	}

	if {[info commands ${mod}::configure_dialogue] != ""} {
	    if {![winfo exists $f.top_frame]} {
	        frame $f.top_frame -bd 0 -relief flat
		pack $f.top_frame -side top -fill both
	    }
	    button $f.top_frame.save \
		-text "Save these parameters" \
		-padx 4 \
		-command "save_all_pressed $f.top_frame $f $mod"
            frame $f.top_separator -bd 2 -relief raised -height 2
	    pack $f.top_frame.save -side left
	    pack $f.top_separator -side top -fill x -padx 10 -pady 5
	}

	if {[info commands ${mod}::create_dialogue] != ""} {
	    ${mod}::create_dialogue $f
	} elseif {![info exists ${mod}::hidden]} {
	    label $f.msg -text "(No configurable parameters)"
	    pack $f.msg -padx 10 -pady 10
	}
    }

    frame $w.button_pad -height 5 -relief flat
    frame $w.buttons -relief flat
    button $w.buttons.run -text "Run" -command run_gui
    label $w.buttons.status_pos -text "" -bd 0
    label $w.buttons.status -text ""
    button $w.buttons.help -text "Help on module"
    pack $w.buttons.run -side left -padx 5 -fill both
    

    if {[keylget pregap4_defs WINDOW_STYLE] != "compact"} {
        button $w.buttons.ok -text "OK" \
		-command "if {\[ok_gui 0\]} {destroy $w}"
        button $w.buttons.cancel -text "Cancel" \
		-command "destroy $w"
	pack $w.buttons.ok $w.buttons.cancel -side left -padx 5 -fill both
    }

    pack $w.buttons.status_pos -side left -padx 5
    place $w.buttons.status -in $w.buttons.status_pos

    pack $w.nb -expand 1 -fill both -side top
    pack $w.button_pad -side top -fill both
    pack $w.buttons -fill both -side top

    $w.nb raise [lindex $modules 0]
}

proc config_panel_col1 {w mod {wmain {}}} {
    set has_dialog 0
    set is_ok 1

    if {![set ${mod}::enabled]} {
	$w update $mod 3 ""
	if {$wmain != ""} {
	    set_status $wmain
	}
	return
    }

    if {[info commands ${mod}::create_dialogue] != ""} {
	set has_dialog 1
    } else {
	set has_dialog 0
    }
    
    if {[info commands ::${mod}::check_params] != ""} {
	if {[::${mod}::check_params] != ""} {
	    set is_ok 0
	} else {
	    set is_ok 1
	}
    }

    if {$has_dialog && $is_ok} {
	set msg "ok    "
	if {$wmain != ""} {
	    set_status $wmain
	}
    } elseif {$has_dialog} {
	set msg "edit  "
	if {$wmain != ""} {
	    set_status $wmain "Module $mod needs configuring"
	}
    } else {
	set msg "-     "
	if {$wmain != ""} {
	    set_status $wmain
	}
    }

    $w update $mod 3 $msg
}

proc config_notebook_raise {w mod} {
    global last_tab

    check_module_status [winfo parent $w] $mod

    if {[info exists last_tab]} {
	if {[info commands ${last_tab}::process_dialogue] != ""} {
	    ${last_tab}::process_dialogue [$w subpanel $last_tab].indent
	}
	config_panel_col1 $w $last_tab
    }
    set last_tab $mod

    set w [winfo parent $w]
    if {[info exists ${mod}::help]} {
	catch {pack $w.buttons.help -side right -fill both} var
	$w.buttons.help configure \
	    -command "eval show_help [list [set ${mod}::help]]"
    } else {
	catch {pack forget $w.buttons.help}
    }
}

proc file_panel {} {
    global pregap4_defs files fofn fofn_dir

    set w [keylget pregap4_defs FILES.WIN]
    if {![winfo exists $w]} {
	xtoplevel $w
	wm title $w "Files to process"
    }

    if {[winfo exists $w.files]} {
	if {[keylget pregap4_defs WINDOW_STYLE] == "compact"} {
	    [keylget pregap4_defs NB.WIN] select 0
	} else {
	    raise $w
	}

	return
    }

    frame $w.files -bd 1 -relief raised
    frame $w.files.text
    label $w.files.text.label -text "\tList of files to process" -bd 5
    text $w.files.text.list -width 40 -height 20 \
	-yscrollcommand "$w.files.text.yscroll set" \
	-xscrollcommand "$w.files.text.xscroll set" \
	-wrap none
    scrollbar $w.files.text.yscroll \
	-orient vert \
	-command "$w.files.text.list yview"
    scrollbar $w.files.text.xscroll \
	-orient hori \
	-command "$w.files.text.list xview"
    foreach f $files {
	$w.files.text.list insert end "$f\n"
    }

    grid rowconfigure $w.files.text 1 -weight 1
    grid columnconfigure $w.files.text 0 -weight 1
    grid $w.files.text.label   -row 0 -column 0 -sticky nsew
    grid $w.files.text.list    -row 1 -column 0 -sticky nsew
    grid $w.files.text.yscroll -row 1 -column 1 -sticky ns
    grid $w.files.text.xscroll -row 2 -column 0 -sticky ew
    pack $w.files.text -side left -fill both -anchor w -expand 1
    pack $w.files -side top -fill both -expand 1

    set right [frame $w.files.right]
    pack $right -side right -fill both -expand 0 -padx 10

    label $right.dummy
    xentry $right.output_fofn \
	-label "Output filename prefix" \
	-width 15 \
	-default $fofn \
	-textvariable fofn
    xentry $right.output_dir \
	-label "Output directory prefix" \
	-width 15 \
	-default $fofn_dir \
	-textvariable fofn_dir
    $right.output_dir xview moveto 1
    button $right.select \
	-text "Add files" \
	-command "add_files $w.files.text.list $right.output_dir"
    button $right.biolims \
	-text "Add BioLIMS files" \
	-command "add_files $w.files.text.list $right.output_dir"
    button $right.add \
	-text "Add file of filenames" \
	-command "add_fofn $w.files.text.list $right.output_fofn \
		  $right.output_dir"
    button $right.clear \
	-text "Clear current list" \
	-command "$w.files.text.list delete 1.0 end"
    button $right.save \
	-text "Save current list to..." \
	-command "save_fofn \[$w.files.text.list get 1.0 end\]"
    pack $right.dummy \
	$right.output_fofn \
	    $right.output_dir \
	    $right.select \
            $right.biolims \
	    $right.add \
	    $right.clear \
	    $right.save \
	    -side top -pady 5 -fill both

    frame $w.button_pad -height 5 -relief flat
    pack $w.button_pad -side top -fill both

    frame $w.buttons -bd 0
    button $w.buttons.run -text "Run" -command run_gui
    button $w.buttons.help -text "Help" \
	-command "show_help pregap4 {Pregap4-Files}"
    pack $w.buttons.run -side left -padx 5 -fill both
    pack $w.buttons.help -side right -padx 5 -fill both
    pack $w.buttons -fill both -side top

    if {[keylget pregap4_defs WINDOW_STYLE] != "compact"} {
	button $w.buttons.ok -text "OK" -command \
		"update_file_list; destroy $w"
	pack $w.buttons.ok -side left -padx 5 -fill both
    }
}

proc config_panel_restart {} {
    global pregap4_defs last_tab

    set w [keylget pregap4_defs CONFIG.WIN]
    set raised [$w.nb raise]
    catch {unset last_tab}
    $w.nb raise init
    foreach win [winfo children $w] {
	destroy $win
    }
    config_panel
    catch {$w.nb raise $raised}
}

proc build_gui {w} {
    global pregap4_defs svn_version

    init_modules

    wm title . "Pregap4 version 1.6$svn_version"
    wm protocol . WM_DELETE_WINDOW exit

    menu_panel $w

    if {[keylget pregap4_defs WINDOW_STYLE] == "compact"} {
	keylset pregap4_defs NB.WIN $w.nb
	ttk::notebook $w.nb \
	    -height 450 \
	    -width 735
        pack $w.nb -side top -fill both -expand 1

	set f [frame $w.nb.files]
	$w.nb add $f -text "Files to Process"
	$f configure -bd 5
	keylset pregap4_defs FILES.WIN $f
	file_panel

	keylset pregap4_defs MODULE.WIN $w.module

	set f [frame $w.nb.configure]
        $w.nb add $f -text "Configure Modules"
	$f configure -bd 5
        keylset pregap4_defs CONFIG.WIN $f
        config_panel

	set f [frame $w.nb.output]
        $w.nb add $f -text "Textual Output"
	$f configure -bd 5
        keylset pregap4_defs OUTPUT.WIN $f
        text_panel

	$w.nb select 0
    } else {
        keylset pregap4_defs CONFIG.WIN $w.config
        keylset pregap4_defs MODULE.WIN $w.module
        keylset pregap4_defs OUTPUT.WIN $w.mainwin
	keylset pregap4_defs FILES.WIN  $w.files

        pack [frame $w.mainwin] -expand 1 -fill both -side bottom
        text_panel

        config_panel
    }

    InitBusy \
	[keylget pregap4_defs OUTPUT.WIN] \
	[keylget pregap4_defs MENU.WIN] \
	pregap4_menu

    # Hack to stop resize flickering
    after idle {after 1000 {wm geometry . [wm geometry .]}}
}

proc build_gui_compact {w} {
    init_modules
    menu_panel $w
}

# Switches to the module panel and raises module "$mod". Used for error
# conditions. Message "$message" is displayed in the status line.
proc raise_mod {mod message} {
    global pregap4_defs

    set w [keylget pregap4_defs CONFIG.WIN]

    if {[keylget pregap4_defs WINDOW_STYLE] == "compact"} {
	[keylget pregap4_defs NB.WIN] select 1
    } else {
	catch {raise [keylget pregap4_defs CONFIG.WIN]}
    }
    $w.nb raise $mod
    bell
    set_status $w $message
}


# Checks that the GUI parameters are OK. If they are it returns 1.
# Otherwise it raises that GUI component and returns 0.
proc ok_gui {check_files} {
    global modules pregap4_defs

    set w [keylget pregap4_defs CONFIG.WIN]

    # If the configure gui dialog does not exist then just call check_modules
    # instead. If that fails, bring up the configure dialogue and recheck
    # to get to the first problematic pane.
    if {![winfo exists $w]} {
	if {[check_modules $check_files] == 0} {
	    return 1
	}
	config_panel
    }

    foreach mod $modules {
	if {[set ${mod}::enabled] == 0} {
	    continue
	}
	if {[info commands ${mod}::process_dialogue] != ""} {
	    if {[${mod}::process_dialogue [$w.nb subpanel $mod].indent] != 1} {
		verror ERR_WARN Run "Module $mod needs configuring"
		raise_mod $mod "Module $mod needs configuring"
		return 0
	    }
	}
    }

    set_status $w

    return 1
}

proc run_gui {} {
    global files interactive modules pregap4_defs fofn_dir

    set w [keylget pregap4_defs CONFIG.WIN]

    if {[keylget pregap4_defs WINDOW_STYLE] == "compact"} {
	[keylget pregap4_defs NB.WIN] select 2
    } else {
	raise [keylget pregap4_defs OUTPUT.WIN]
    }
    update

    if {![ok_gui 1]} {
	return
    }
    update_file_list
    if {[regexp "^\[ \n\t\r\]*$" $files]} {
	verror ERR_WARN check_modules "No input files specified"
	return	
    }

    if {[keylget pregap4_defs WINDOW_STYLE] != "compact"} {
	destroy $w
    }

    SetBusy
    set cur_pwd [pwd]

    catch {
	if {[set modname [init_modules]] == ""} {
	    set orig_files $files
	    set files [run_modules $files]
	    shutdown_modules $files
	    set files $orig_files
	} else {
	    global errorInfo
	    raise_mod $modname "${modname}::init: $errorInfo"
	}
    }
    cd $cur_pwd
    
    ClearBusy
}

proc save_all_params {} {
    global modules pregap4_defs

    set w [keylget pregap4_defs CONFIG.WIN]

    vfuncheader "Saving parameters"
    store_module_list

    config_panel 0
    foreach mod $modules {
	namespace eval ::$mod {mod_save enabled [set enabled]}

	if {[set ${mod}::enabled] == 0} {
	    continue
	}

	if {[info commands ${mod}::configure_dialogue] != ""} {
	    set f [$w.nb subpanel $mod].indent
	    ${mod}::configure_dialogue $f save_all
	    vmessage "Wrote parameters for module $mod."
	}
    }

    write_conf_file
}

proc save_all_params_to {} {
    global conf_file

    set fname [tk_getSaveFile -initialdir [pwd]]

    if {$fname == ""} {
	return
    }

    if {[file exists $fname]} {
	file delete $fname
    }

    set conf_file $fname
    save_all_params
}

proc save_module_list {} {
    vfuncheader "Saving parameters"

    store_module_list
    store_module_states
    write_conf_file
}

proc add_biolims_files {w dir_w} {

    set flist [spGetOpenBiolims \
		-multiple true \
		-parent [winfo toplevel $w]]

    if {$flist == ""} {
	return
    }

    lappend files $flist
    foreach file $flist {
	$w insert end "$file\n"
    }

    $dir_w delete 0 end
    $dir_w insert end [file dirname [lindex $flist 0]]
    $dir_w xview moveto 1

}

proc add_files {w dir_w} {
    set types {
	{"ABI"			*.ab1		}
	{"ALF"			*.alf		}
	{"EXP"			*.exp		}
	{"SCF"			*.scf		}
	{"ZTR"			*.ztr		}
	{"PLN"			*		}
	{"Any"			*		}
    }

    set flist [tk_getOpenFile \
		   -initialdir [pwd] \
		   -multiple 65000 \
		   -filetypes $types \
		   -parent [winfo toplevel $w]]
    
    if {$flist == ""} {
	return
    }

    lappend files $flist
    foreach file $flist {
	$w insert end "$file\n"
    }

    $dir_w delete 0 end
    $dir_w insert end [file dirname [lindex $flist 0]]
    $dir_w xview moveto 1

    catch {cd [file dirname [lindex $flist 0]]}
}

proc add_fofn {w prefix_w dir_w} {
    global files fofn tcl_platform

    set fname [tk_getOpenFile -initialdir [pwd] \
		-parent [winfo toplevel $w]]
    if {$fname == ""} {
	return
    }

    set dir [file dirname $fname]; # Always absolute

    if {$tcl_platform(platform) == "windows"} {
	regsub {^(.):.*} $dir {\1:} fofn_vol
    }

    set fd [open $fname r]
    while {[gets $fd line] != -1} {
	set line [string trim $line]
	if {$line != ""} {
	    switch [file pathtype $line] {
		absolute {
		    set n $line
		}
		relative {
		    set n [file join $dir $line]
		}
		volumerelative {
		    if {[string match /* $line]} {
			# Absolute, but no volume ID
			set n $fofn_vol$line
		    } else {
			set n $line
		    }
		}
	    }
	    lappend files $n
	    $w insert end "$n\n"
	}
    }
    close $fd

    $prefix_w delete 0 end
    $prefix_w insert end [file tail $fname]
    $prefix_w xview moveto 1

    $dir_w delete 0 end
    $dir_w insert end [file dirname $fname]
    $dir_w xview moveto 1

    catch {cd [file dirname $fname]}
}

proc save_fofn {files} {
    set fname [tk_getSaveFile -initialdir [pwd]]

    if {$fname == ""} {
	return
    }

    vfuncheader "Save file of filenames"
    if {[catch {set fd [open $fname w]}]} {
        tk_messageBox \
	    -icon error \
	    -title Error \
	    -message "Could not write to \"$fname\"."
	return
    }

    puts $fd $files
    close $fd
}

proc update_file_list {} {
    global files pregap4_defs

    set w [keylget pregap4_defs FILES.WIN].files.text.list

    if {![winfo exists $w]} {
	return
    }

    set files [$w get 1.0 end]
    regsub -all "(\n+$)|((\n)\n+)" $files {\3} files
    set files [split $files \n]
}

proc load_new_config {} {
    global files pregap4_defs

    set fname [tk_getOpenFile -initialdir [pwd]]
    if {$fname == ""} {
	return
    }

    if {[catch {read_conf_file $fname} err]} {
	tk_messageBox \
	    -icon error \
	    -title "Error" \
	    -message "Could not load configuration file.\
		     Received error \"$err\"." \
	    -type ok
    }

    load_modules
    init_modules
    if {[winfo exists [keylget pregap4_defs CONFIG.WIN]]} {
	config_panel_restart
    }
}

# Updates the user's $HOME/.pregap4rc file with the WINDOW_STYLE
proc set_window_style {style} {
    global env
    if {[catch {set fd [open $env(HOME)/.pregap4rc a+]}]} {
	verror ERR_WARN set_window_style \
	    "Could not create $env(HOME)/.pregapr4rc"
    }

    puts $fd "set_def WINDOW_STYLE $style"
    close $fd

    tk_messageBox -message "The change of window style will take affect the next time that Pregap4 is started." -type ok
}

proc set_status {w {message {}}} {
    $w.buttons.status configure -text $message
}

proc check_module_status {w mod} {
    if {[set ::${mod}::enabled] && \
	    [info commands ::${mod}::create_dialogue] != "" && \
	    [info commands ::${mod}::check_params] != "" && \
	    [::${mod}::check_params] != ""} {
	set_status $w "Module $mod needs configuring"
    } else {
	set_status $w
    }
}

# Select (state=1) or deselect (state=0) all modules. 
proc select_all_modules {state} { 
    global modules 
    foreach mod $modules { 
        if {![info exists ::${mod}::mandatory]} { 
            set ::${mod}::enabled $state 
        } 
    } 
} 
