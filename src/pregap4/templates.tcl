#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1998. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#
# This component adds standard template components to the memory or to
# a configuration file.
#
# Templates consist of standard chunks of configuration which may be useful
# to load. Examples are:
#	A naming scheme parse
#	Fetching fields from a database
#

proc add_template_gui {is_name_scheme} {
    global env

    set w .add_template
    if {[modal $w] == ""} {
	return
    }
    if {$is_name_scheme} {
	wm title $w "Load naming scheme"
    } else {
	wm title $w "Include config component"
    }

    # Remember existing selectPath for filebrowser as we do not want it
    # switching to the naming_schemes or templates directory.
    if {[info exists ::tk::dialog::file::__tk_filedialog(selectPath)]} {
	set old_path $::tk::dialog::file::__tk_filedialog(selectPath)
    } else {
	set old_path ""
    }

    if {$is_name_scheme} {
	get_fname $w.file \
		-text "Naming scheme file name" \
		-type load \
		-initialdir $env(STADTCL)/pregap4/naming_schemes \
		-filetypes {{"template"	*.p4t}}
	set help "Pregap4-Naming"
    } else {
	get_fname $w.file \
		-text "Component file name" \
		-type load \
		-initialdir $env(STADTCL)/pregap4/templates \
		-filetypes {{"template"	*.p4t}}
	set help "Pregap4-Components"
    }

    xyn $w.save \
	-label "Save to config file" \
	-orient horiz \
	-default 1

    okcancelhelp $w.ok \
	-ok_command "add_template_gui2 $w $w.file \[$w.save get\] [list $old_path]" \
	-cancel_command "add_template_cancel $w [list $old_path]" \
	-help_command "show_help pregap4 $help" \
	-bd 2 -relief groove

    pack $w.file $w.save $w.ok -side top -fill both
}

proc add_template_cancel {w old_path} {
    if {$old_path == ""} {
	catch {unset ::tk::dialog::file::__tk_filedialog(selectPath)}
    } else {
	set ::tk::dialog::file::__tk_filedialog(selectPath) $old_path
    }
    
    destroy $w
}

proc add_template_gui2 {w file save old_path} {
    if {[set file [$w.file get]] == ""} { bell; return }

    if {$old_path == ""} {
	catch {unset ::tk::dialog::file::__tk_filedialog(selectPath)}
    } else {
	set ::tk::dialog::file::__tk_filedialog(selectPath) $old_path
    }

    vfuncheader "Insert Configuration Component"
    add_template $file $save

    if {$save} {
	vmessage "Copied [file tail $file] into pregap4 configuration file."
    } else {
	vmessage "Loaded [file tail $file] into memory."
    }

    destroy $w
}

proc add_template {fname save} {
    set fd [open $fname r]
    set data [read $fd]
    close $fd
    set lines [split $data \n]
    if {[regsub {^\[(.*)\]$} [lindex $lines 0] {\1} ns] == 0} {
	set ns "template $fname"
    } else {
	set data [join [lrange $lines 1 end] \n]
    }

    if {[catch {uplevel #0 eval [list $data]} var]} {
	vmessage $var
    }

    if {$save} {
	set_conf_data_section $ns $data\n
	write_conf_file
    }
}