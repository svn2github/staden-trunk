if {[catch {package require Iwidgets}]} {
    tk_dialog .no_iwidgets "Iwidgets error" \
"No working copy of the Tcl Iwidgets package was found. This may be due to a version conflict with Tcl/Tk or the lack of IncrTcl, IncrTk and/or Iwidgets packages installed.

Without these existing, this GUI dialogue will be absent." error 0 Ok

    proc prefinish {args} {}
} else {

namespace import itcl::*

# -----------------------------------------------------------------------------
# A Gap4 interface to the prefinish library and GUI.
#
proc prefinish {io} {
    global gap_defs
    load_package prefinish
    ::prefinish::init

    set w .prefinish
    if {[xtoplevel $w] == ""} return
    wm title $w "Prefinish"

    # Contig(s) to process
    contig_id $w.id -io $io

    global defs_c_in
    lorf_in $w.infile $defs_c_in \
	"{contig_id_configure $w.id -state disabled} \
	    {contig_id_configure $w.id -state disabled}\
	    {contig_id_configure $w.id -state disabled}\
	    {contig_id_configure $w.id -state normal}" -bd 2 -relief groove

    pack $w.infile $w.id -side top -fill both

    # Prefinish configuration file
    foreach {dbname vers} [split [db_info db_name $io] .] {}
    xget_fname $w.output \
	-text "Prefinish output file" \
	-type save \
	-default $dbname.$vers.prefinish
    xget_fname $w.config \
	-text "Prefinish configuration file" \
	-type load \
	-default [keylget gap_defs PREFINISH.CONFIG_FILE]
    button $w.edit_config \
	-text "Edit prefinish configuration" \
	-command "prefinish_config $io $w"
    pack $w.output $w.config -side top -fill both
    pack $w.edit_config -side top -anchor w

    okcancelhelp $w.ok \
	-ok_command "prefinish_run $io $w" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap4 {Prefinish}" \
	-bd 2 -relief groove

    pack $w.ok -side top -fill both
}

# -----------------------------------------------------------------------------
# PRIVATE functions below

# Invokes the prefinish main dialogue, with a little tweaking to remove the
# "Exit" command and to add a "Run" command.
proc prefinish_config {io wtop} {
    # Create the config window
    set w $wtop.edit
    if {[xtoplevel $w] == ""} return
    wm title $w "Prefinish Configuration"
    set h [prefinish::main_gui $w]

    if {[set fname [$wtop.config get2]] != ""} {
	::prefinish::prefin_load $h $fname
    }

    # Replace the File menu
    destroy $w.m.file
    menu $w.m.file
    $w.m.file add command \
	-label Load \
	-command "::prefinish::prefin_load $h"
    $w.m.file add command \
	-label Save \
	-command "prefinish_save $h $wtop"
    $w.m.file add separator
    $w.m.file add command \
	-label Close \
	-command "prefinish_close $h $w $wtop"

    wm protocol $w WM_DELETE_WINDOW "prefinish_close $h $w $wtop"
}

# Overloads the normal prefinish dialogue save command to also update the
# filename in the dialogue window that invoked this window.
proc prefinish_save {h w} {
    set fname [::prefinish::prefin_save $h]
    if {$fname == ""} {
	return ""
    }
    
    $w.config delete 0 end
    $w.config insert end $fname

    return $fname
}

# Overloads the normal prefinish dialogue exit command to also update the
# filename in the dialogue window that invoked this window.
proc prefinish_close {h w wtop} {
    if {[prefinish::check_saved $h]} {
	destroy $w
    }

    global $h.LastFilename
    if {[info exists $h.LastFilename]} {
	$wtop.config delete 0 end
	$wtop.config insert end [set $h.LastFilename]
    }
}

#
# The prefinish "Run" command.
# We generate the prefinish procedures using prefinish_lib.tcl and then append
# to this 'script' a series of process_contig calls.
# This script is written to a temporary file and we open a new stash process
# via a pipe to run it.
#
# Output from the pipe is sent to the text output and error windows.
# Any additional input on the pipe is used to indicate the controlling process
# (gap4) wishes to terminate the child process. (We use this method as windows
# does not appear to support signals in the same manner as Unix.)
#
proc prefinish_run {io w} {
    global tcl_platform

    # Find the config file name.
    set fname [$w.config get]
    if {$fname == ""} {
	bell
	return
    }

    # Find the output file name
    set output_file [$w.output get]
    if {$output_file == ""} {
	bell
	return
    }

    vfuncheader "Prefinish"

    # Get the list of contigs
    if {[lorf_in_get $w.infile] == 4} {
	set gel_name [contig_id_gel $w.id]
	set lreg [contig_id_lreg $w.id]
	set rreg [contig_id_rreg $w.id]
	
	SetContigGlobals $io $gel_name $lreg $rreg
	set clist "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $w.infile] == 3 } {
	set clist [CreateAllContigList $io]
    } else {
	set clist [lorf_get_list $w.infile]
    }

    # Create some temporary filenames for output and the script to run
    foreach {dbname vers} [split [db_info db_name $io] .] {}
    set tmpbase [tmpnam prefinish]
    set stderr_fn $tmpbase.stderr
    set script_fn $tmpbase.script
    set script_fd [open $script_fn w]

    # Load the config file and generate the code for our child process
    xtoplevel $w.gui
    wm withdraw $w.gui
    set h [::prefinish::main_gui $w.gui]
    ::prefinish::prefin_load $h $fname
    set appcode [::prefinish::generate_app2 $h]
    destroy $w.gui

    # Load the necessary libraries and packages
    if { $tcl_platform(platform) == "windows" } {
	puts $script_fd {catch {wm withdraw .}}
    }
    puts $script_fd {source $env(STADTABL)/shlib.conf}
    puts $script_fd {load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}}
    puts $script_fd {load_package gap}
    puts $script_fd {load_package prefinish}
    puts $script_fd {set consensus_cutoff 0.02}
    puts $script_fd {set quality_cutoff 1}
    puts $script_fd {set consensus_mode 2}

    # Send over the generated code
    puts $script_fd $appcode
    puts $script_fd {
	array set opt {
	    -add_tags 0
	    -dump_problems ""
	    -skip_solutions 0
	}
    }

    # Open the database
    puts $script_fd "set io \[open_db -name $dbname -version $vers -access ro\]"
    puts $script_fd "finish .f -io \$io -output_file $output_file"

    set contigs ""
    foreach contig $clist {
	lappend contigs [lindex $contig 0]
    }

    destroy $w

    set clist2 ""
    foreach contig $clist {
	foreach {cname start end} $contig {}
	set cnum [db_info get_contig_num $io $cname]
	if {$start == ""} {
	    set start 1
	}
	if {$end == ""} {
	    set c [io_read_contig $io $cnum]
	    set end [keylget c length]
	}
	lappend clist2 $cnum $start $end
    }

    puts $script_fd "fileevent stdin readable {
        if {\[read stdin\] == \"STOP\\n\"} {puts stderr Terminated; exit}
    }"
    puts $script_fd "foreach {cnum start end} [list $clist2] {
    process_contig \$io .f opt \$cnum \$start \$end \$class_bits
    update
}\n"
    puts $script_fd {puts "--- FINISHED ---"}
    close $script_fd

    # We've created the script, so now run it.
    if { $tcl_platform(platform) != "windows" } {
	set fd [open "|stash $script_fn 2> $stderr_fn" w+]
    } else {
	set fd [open "|wish $script_fn 2> $stderr_fn" w+]
    }
    fconfigure $fd -blocking 0

    # Create stop/animation window
    set t [xtoplevel .mainwin.stdout.stop]
    wm title $t "Working..."
    frame $t.dummy -width 300 -height 0 -bd 0
    pack $t.dummy -side top
    frame $t.labels -bd 0
    pack $t.labels -side top -fill both -expand 1
    label $t.labels.l -text "Starting"
    label $t.labels.r -text "" -font fixed
    pack $t.labels.l $t.labels.r -side left
    button $t.stop \
	-text Cancel \
	-command "prefinish_interrupt $fd $w.Done"
    pack $t.stop -side top -pady 3
    after 1000 "prefinish_animate $t.labels.r"

    global $w.Done
    set $w.Done 0
    fileevent $fd readable \
	"prefinish_run_output $fd $w.Done $t.labels.l [list $contigs]"

    # And await finishing
    SetBusy
    vwait $w.Done
    ClearBusy
    destroy $t

    # Load stderr and display
    set fd [open $stderr_fn r]
    verror ERR_WARN prefinish [read $fd]
    close $fd
    catch {file delete -force $stderr_fn}
    catch {file delete -force $script_fn}

    return
}

#
# An 'after <time>' callback procedure; used to animate the window containing
# the cancel button so that we can see it hasn't crashed.
#
proc prefinish_animate {label} {
    if {![winfo exists $label]} {
	return
    }

#    set states {
#	"|    "
#	">>   "
#	" >>  "
#	"  >> "
#	"   >>"
#	"    |"
#	"   <<"
#	"  << "
#	" <<  "
#	"<<   "
#    }

    set states {
	"   "
	".  "
	".. "
	"..."
    }
    set ind [expr {[lsearch $states [$label cget -text]]+1}]
    if {$ind >= [llength $states]} {
	set ind 0
    }
    $label configure -text [lindex $states $ind]
    
    after 707 "prefinish_animate $label"
}

proc prefinish_run_output {fd done label contigs} {
    global $done

    flush $fd
    set data [read $fd]

    set ncontigs [llength $contigs]

    foreach {_ cnum cname} [regexp -all -line -inline \
			   {^\#\#\# CONTIG .* \(=([0-9]+) (.*)\)} $data] {
        set count 1
        foreach contig $contigs {
	    if {$cname == [lindex $contig 0]} {
		$label configure \
		    -text "Processing contig $cname ($count of $ncontigs)"
	    }
	    incr count
	}
    }

    vmessage $data
    if {[eof $fd]} {
	close $fd
	set $done 1
    }
}

#
# Called by the cancel button. We just close the file descriptor which will
# cause an error on the child process the next time it attempts to output.
# Not ideal, but sufficient (I think).
#
proc prefinish_interrupt {fd done} {
    puts $fd STOP
    close $fd
    vmessage "--- MANUALLY INTERRUPTED ---"
    global $done
    set $done 1
}

}; #Iwidgets load