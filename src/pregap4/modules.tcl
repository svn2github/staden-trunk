#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1998. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#-----------------------------------------------------------------------------
# Source .p4m files found in MODULE_PATH
#
# This creates a name space for each module, sources the .p4m file to create
# functions within that module, and then also adds a mod_preset function in
# that module (for handling module parameters).
#
# init and shutdown are special modules that are always first and last.
# They do not need to be specified in MODULES
proc load_modules {} {
    global modules MODULE_PATH MODULES tcl_platform env

    set modules "init $MODULES shutdown"
    set new_modules ""
    array set mod_arr1 {}
    foreach mod $modules {
	if {$tcl_platform(platform) == "windows" && $mod == "compress_trace"} {
	    continue
	}
	if {[info exists mod_arr1($mod)]} {
	    incr mod_arr1($mod)
	    lappend new_modules $mod#$mod_arr1($mod)
	} else {
	    set mod_arr1($mod) 1
	    lappend new_modules $mod
	}
    }
    set modules $new_modules

    # Scan through MODULE_PATH looking for modules
    foreach path $MODULE_PATH {
	set path [subst $path]
	if {[catch {set mods [glob $path/*.p4m]}]} {
	     continue
	}

	foreach mod $mods {
	    set mod [file tail $mod]
	    regsub {(.*)\.p4m$} $mod {\1} mod_name
	    regsub -all { } $mod_name {_} mod_name

	    # If it's already loaded - skip to next one
	    if {![info exists mod_arr1($mod_name)]} {
		continue
	    }

	    # Source and setup the module
	    for {set i 1} {$i <= $mod_arr1($mod_name)} {incr i} {
		if {$i > 1} {
		    set mname $mod_name#$i
		} else {
		    set mname $mod_name
		}
		namespace eval $mname [list source $path/$mod]
		if {![info exists ${mname}::enabled]} {
		    set ${mname}::enabled 1
		}

		namespace eval $mname {
		    proc mod_preset {var val} {
			variable $var
			if {![info exists $var]} {
			    set $var $val
			}
		    }

		    proc mod_save {var {val {}}} {
			set_conf_data_module_param \
			    [string trimleft [namespace current] ::] \
			    $var $val
		    }

		    proc glob_save {var {val {}}} {
			set_conf_data_globals $var $val
		    }
		}
	    }
	    unset mod_arr1($mod_name)
	}
    }
}

# Call init function
#
# Returns "" if all modules successfully initialise.
# Returns "module_name" if one or more modules fail to initialise. (Only the
#         first module if several fail.)
# errorInfo will also be left set to whatever the first module errors with.
#
proc init_modules {} {
    global modules
    global errorInfo

    set ret ""
    set einfo ""
    foreach mod $modules {
	set ${mod}::report ""
	if {[info commands ${mod}::init] != ""} {
	    update idletasks
	    set errorInfo ""
	    if {[catch {${mod}::init} var]} {
		regsub "\n.*" $errorInfo {} errorInfo
		verror ERR_WARN ${mod}::init $errorInfo
		if {$ret == ""} {
		    set ret $mod
		    set einfo $errorInfo
		}
	    }
	}
    }

    set errorInfo $einfo
    return $ret
}

# call check_param function. Returns the number of problematic modules
# 'input_files' determines whether to check for no input files.
proc check_modules {input_files} {
    global modules files interactive

    vfuncheader "Checking modules"
    if {!$interactive} {
	puts "=== Checking Module Parameters ==="
    }
    set problem 0
    foreach mod $modules {
	if {[set ${mod}::enabled] == 0} {
	    continue
	}
	if {[info commands ${mod}::check_params] != ""} {
	    vmessage "- [${mod}::name] -"
	    if {[set inv [${mod}::check_params]] != ""} {
		vmessage "Invalid value for parameter '$inv'."
		incr problem
	    }
	    update idletasks
	}
    }

    if {$input_files && [regexp "^\[ \n\t\r\]*$" $files]} {
	verror ERR_WARN check_modules "No input files specified"
	incr problem
    }
    return $problem
}

# call run function
proc run_modules {files} {
    global modules interactive
    vfuncheader "Running modules"
    if {!$interactive} {
	puts "\n=== Running Modules ==="
    }
    foreach mod $modules {
	if {[set ${mod}::enabled] == 0} {
	    continue
	}
	if {[info commands ${mod}::run] != ""} {
	    vmessage "- [${mod}::name] -"
	    update idletasks
	    if {[catch {set files [${mod}::run $files]} var]} {
		verror ERR_WARN ${mod}::run $var
	    }
	}
    }

    return $files
}

# call shutdown function
proc shutdown_modules {files} {
    global modules interactive
    vmessage ""
    vfuncheader "Terminating modules"
    if {!$interactive} {
	puts "=== Module Shutdown ==="
    }
    foreach mod $modules {
	if {[info commands ${mod}::shutdown] != ""} {
	    vmessage "- [${mod}::name] -"
	    update idletasks
	    if {[catch {set files [${mod}::shutdown $files]} var]} {
		verror ERR_WARN ${mod}::shutdown $var
	    }
	}
    }

    return $files
}

# The module command can be used to execute a script within a module.
# It's nothing more than "namespace eval", but hides the implementation away
# from the user.
proc module {modname script} {
    namespace eval $modname $script
}


#
# -----------------------------------------------------------------------------
# Produces a window for editing the MODULES and MODULE_PATH variables
#
proc select_modules {} {
    global pregap4_defs

    set win [keylget pregap4_defs MODULE.WIN]
    if {[keylget pregap4_defs WINDOW_STYLE] == "compact" &&
	[winfo exists $win]} {
	[keylget pregap4_defs NB.WIN] select 1
    } else {
	if {[modal $win -resizable 1] == ""} {
	    return
	}
	wm title $win "Add/Remove Modules"
	
	select_modules_create $win
	wm withdraw $win
	set x [expr [winfo screenwidth $win]/2-[winfo reqwidth $win]/2 \
		- [winfo vrootx [winfo parent $win]]]
	set y [expr [winfo screenheight $win]/2-[winfo reqheight $win]/2 \
		- [winfo vrooty [winfo parent $win]]]
	wm geometry $win +$x+$y
	wm deiconify $win
    }
}

proc select_modules_create {w} {
    global MODULE_PATH MODULES

    label $w.label \
	-text "(IMPORTANT: Please read the online help for this dialogue)" \
	-font {Helvetica -15 bold} \
	-fg "#800000"
    pack $w.label -side top -pady 3

    frame $w.separator0 -bd 2 -relief raised -height 2
    pack $w.separator0 -side top -fill x -padx 10 -pady 5

    # The search path
    frame $w.path
    xentry $w.path.e \
	-label "Module search path" \
	-default $MODULE_PATH \
	-width 40

    bind $w.path.e <Any-Leave> \
	[namespace code "update_modules $w \[%W get\]"]
    bind $w.path.e.entry <KeyPress-Return> \
	[namespace code "update_modules $w \[%W get\]"]

    pack $w.path.e -side left -fill both -expand 1
    pack $w.path -side top -fill x

    frame $w.separator1 -bd 2 -relief raised -height 2
    pack $w.separator1 -side top -fill x -padx 10 -pady 5

    # Create and pack the two listboxes
    frame $w.l2
    frame $w.l2.to
    label $w.l2.to.t -text "Modules to use:"
    listbox $w.l2.to.l -yscroll "$w.l2.to.s set" -width 30 -height 22
    scrollbar $w.l2.to.s -orient vert -command "$w.l2.to.l yview"

    frame $w.l2.from
    label $w.l2.from.t -text "Modules available:"
    listbox $w.l2.from.l -yscroll "$w.l2.from.s set" -width 30 -height 22
    scrollbar $w.l2.from.s -orient vert -command "$w.l2.from.l yview"

    frame $w.l2.inst
    label $w.l2.inst.l1 -text "Use the left mouse button to select and drag the modules to and from the left list." -wraplength 100
    label $w.l2.inst.l2 -text "Note that the order of modules is important and that Pregap4 will not check this." -wraplength 100

    grid $w.l2.to.t -sticky ew
    grid $w.l2.to.l -sticky nsew
    grid $w.l2.to.s -row 1 -column 1 -sticky ns
    grid $w.l2.from.t -sticky ew
    grid $w.l2.from.l -sticky nsew
    grid $w.l2.from.s -row 1 -column 1 -sticky ns
    grid rowconfigure $w.l2.to	 1 -weight 1
    grid rowconfigure $w.l2.from 1 -weight 1
    grid columnconfigure $w.l2.to   0 -weight 1
    grid columnconfigure $w.l2.from 0 -weight 1

    pack $w.l2 -side top -fill both -expand 1
    pack $w.l2.to $w.l2.from -side left -fill both -expand 1
    pack $w.l2.inst -side left -fill both -anchor c
    pack $w.l2.inst.l1 $w.l2.inst.l2 -side top -fill both -expand 1

    update_modules $w $MODULE_PATH 1

    bind $w.l2.from.l <B1-Motion>	  "dnd_motion  %X %Y; break"
    bind $w.l2.from.l <B2-Motion>	  "dnd_motion  %X %Y; break"
    bind $w.l2.from.l <B3-Motion>	  "dnd_motion  %X %Y; break"
    bind $w.l2.from.l <Any-ButtonPress>	  "dnd_grab    $w.l2.to.l \
						       $w.l2.from.l \
						       %x %y 1 0 1 0"
    bind $w.l2.to.l   <B1-Motion>	  "dnd_motion  %X %Y; break"
    bind $w.l2.to.l   <B2-Motion>	  "dnd_motion  %X %Y; break"
    bind $w.l2.to.l   <B3-Motion>	  "dnd_motion  %X %Y; break"

    bind $w.l2.from.l <Any-ButtonRelease> "dnd_release"
    bind $w.l2.to.l   <Any-ButtonPress>	  "dnd_grab    $w.l2.from.l \
						       $w.l2.to.l \
						       %x %y 0 1 0 1"
    bind $w.l2.to.l   <Any-ButtonRelease> "dnd_release"

    global grabbed_data
    set grabbed_data ""

    # Create the buttons
    frame $w.buttons
    button $w.buttons.done \
	-text "OK" \
	-command [namespace code "select_modules_done $w"]
    button $w.buttons.cancel \
	-text "Cancel" \
	-command [namespace code "select_modules_cancel $w 0"]
    button $w.buttons.help \
	-text "Help" \
	-command "show_help pregap4 {Pregap4-ModAdd}"
    pack $w.buttons.done $w.buttons.cancel -side left -padx 5
    pack $w.buttons.help -side right -padx 5
    pack $w.buttons -side bottom -fill both
}

proc select_modules_done {w} {
    global MODULE_PATH MODULES modules

#    global pregap4_defs
#    set nbw [keylget pregap4_defs CONFIG.WIN].nb
#    $nbw raise init

    set MODULE_PATH [$w.path.e get]
    set MODULES {}
    foreach mod [$w.l2.to.l get 0 end] {
	regsub "\.p4m$" $mod {} mod
	lappend MODULES $mod
    }

    select_modules_cancel $w 1
}

proc select_modules_cancel {w recreate} {
    global pregap4_defs

    if {$recreate} {
	load_modules
	init_modules
	if {[winfo exists [keylget pregap4_defs CONFIG.WIN]]} {
	    config_panel_restart
	}
    }

    if {[keylget pregap4_defs WINDOW_STYLE] == "compact"} {
	[keylget pregap4_defs NB.WIN] select 1
    }

    destroy $w
}

# Initialise/update data in the listboxes, eg upon creation of the display
# or when the MODULE_PATH changes
proc update_modules {w path {first_time 0}} {
    global MODULES env

    # Get complete list of known modules
    array set avail_mods {}
    foreach dir $path {
	set dir [subst $dir]
	if {[catch {set mods [glob $dir/*.p4m]}]} {
	     continue
	}
	foreach mod $mods {
	    set mod [file tail $mod]
	    set avail_mods($mod) 1
	}
    }
    catch {unset avail_mods(init.p4m)}
    catch {unset avail_mods(shutdown.p4m)}

    #$w.l2.to.l delete 0 end
    $w.l2.from.l delete 0 end

    # Check current settings (if any) and remove unknown modules
    foreach lb [list $w.l2.to.l $w.l2.from.l] {
	for {set s [expr [$lb size]-1]} {$s >= 0} {incr s -1} {
	    set mod [$lb get $s]
	    if {![info exists avail_mods($mod)]} {
		$lb delete $s
	    } else {
		# unset avail_mods($mod)
	    }
	}
    }

    # Add modules listed in MODULES (if we can find them) to the 'to' list
    if {$first_time} {
        foreach mod $MODULES {
	    if {[info exists avail_mods($mod.p4m)]} {
	        set mod_to($mod.p4m) 1
	        $w.l2.to.l insert end $mod.p4m
	        # unset avail_mods($mod.p4m)
	    }
	}
    }

    # Add all others to the 'from' list
    foreach mod [lsort -dictionary [array names avail_mods]] {
	if {[info exists avail_mods($mod)]} {
	    $w.l2.from.l insert end $mod
	}
    }
}

proc store_module_states {} {
    global modules

    foreach mod $modules {
	namespace eval ::$mod {mod_save enabled [set enabled]}
    }
    vmessage "Wrote module list on/off statuses."
}

proc store_module_list {} {
    global modules MODULE_PATH env

    set SRlen [string length $env(STADENROOT)]
    set new_mod_path ""
    foreach path $MODULE_PATH {
	if {[string first $env(STADENROOT) $path] == 0} {
	    set path "\$env(STADENROOT)[string range $path $SRlen end]"
	}
	append new_mod_path "$path "
    }
    set str "set MODULE_PATH [list $new_mod_path]\n"
    append str "set MODULES {\n"
    foreach mod $modules {
	if {$mod == "init" || $mod == "shutdown"} {
	    continue
	}
	regsub #.* $mod {} mod
	append str "	$mod\n"
    }
    append str "}\n"

    set_conf_data_section module_list $str    
    vmessage "Wrote selected module list."
}
