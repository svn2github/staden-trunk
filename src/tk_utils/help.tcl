#
# Copyright (c) 1995 Medical Research Council, Laboratory of Molecular Biology.
# All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#
# Standalone HTML viewer (v 0.1) to use Stephen Uhler's html-library.
# The HTML library is Copyright (c) 1995 by Sun Microsystems.
#

#
# External functions:
#	help_display {name}
#	help_init {{top ""}}
#
# Internal functions:
#	help_* {}
#	HM* {}
#
# internal global variables
#	help_history
#	help_future
#	help_win
#	help_toplevel
#	help_file
#	help_recurse

# Maximum number of links to remember for the "Back" command.
global MAX_LINK_HISTORY
set MAX_LINK_HISTORY 10

source $env(STADTABL)/help_config

#-----------------------------------------------------------------------------
# Externally callable interfaces
proc help_display {name} {
    global help_path help_win help_toplevel

    if {$help_toplevel != ""} {
       	wm deiconify $help_toplevel
    	raise $help_toplevel
    } else {
	wm deiconify .
	raise .
    }
    HMlink_callback $help_win $name
}

proc help_init {{top ""} {external 0}} {
    global help_file help_recurse help_history help_future help_win
    global help_toplevel help_path

    set help_file ""
    set help_history ""
    set help_future ""
    set help_recurse 0
    set help_toplevel $top

    if {$top != "" && ![winfo exists $top]} {
	xtoplevel $top
	wm title $top "Help"
    } else {
	set top .
    }

    # The main window. There's quite a bit to do to disable editing features
    # for this text display whilst retaining editing ability from the code.
    set t [text $help_toplevel.text -yscroll \
	"$help_toplevel.scroll set" -insertontime -1]
    set help_win $t
    bind $t <Any-Key> {break}
    bind $t <ButtonRelease-2> {break}
    scrollbar $help_toplevel.scroll \
	-command "$help_toplevel.text yview"
    wm protocol $top WM_DELETE_WINDOW "help_quit $top"

    # The menu bar
    set w [frame $help_toplevel.top -bd 2 -relief raised]
    menubutton $w.file -text "File" -menu $w.file.m
    set m [menu $w.file.m]

    if {$external} {
        $m add command -label "New window" -command "help_new_win"
        $m add separator
    }
    foreach i [glob $help_path/*_toc.html] {
	regsub "${help_path}/(.*)_toc.html" $i {\1} x
	$m add command -label "File '$x'" -command "HMlink_callback $t $i"
    }
    $m add separator
    $m add command -label "Quit" -command "help_quit $top"

    menubutton $w.font -text "Fonts" -menu $w.font.m
    set m [menu $w.font.m]
    $m add command -label "small" \
	-command "HMset_state $t -size -5; help_redraw_link $t"
    $m add command -label "medium" \
	-command "HMset_state $t -size -2; help_redraw_link $t"
    $m add command -label "large" \
	-command "HMset_state $t -size 3; help_redraw_link $t"

    button $w.back -text "Back" -command "help_back_link $t" -state disabled
    button $w.forw -text "Forward" -command "help_forward_link $t" \
	-state disabled
    button $w.stop -text "Stop" -command "HMstop_parse $t"

    # The status line
    entry $w.status -text ""
    bind $w.status <Return> \
	"HMlink_callback $help_toplevel.text \[%W get\]"

    pack $w.file $w.font -side left
    pack $w.stop $w.forw $w.back -side right
    pack $w.status -padx 0.5c -side left -fill both -expand 1
    pack $w -side top -fill both
    pack $help_toplevel.scroll -side right -fill both
    pack $help_toplevel.text -side left -fill both -expand 1

    #tk_focusFollowsMouse
    #bind Entry <Delete> [bind Entry <BackSpace>]

    trace variable help_history w help_history_trace
    trace variable help_future w help_future_trace
	
    HMinit_win $t
    HM_our_init
    HMset_state $t -size -2
    HMset_indent $t 1.5

    return $t
}

#-----------------------------------------------------------------------------
# Internal help display functions
proc help_new_win {} {
    exec stash -f ./help.tcl -init &
}

proc help_quit {top} {
    global help_win help_history help_future
    HMstop_parse $help_win
    catch {unset help_history help_future}; #remove traces
    destroy $top
}

proc help_back_link {win} {
    global help_history help_future
    set l [expr [llength $help_history]-2]
    if {$l < 0} {
	bell
	return
    }

    lappend help_future [lrange $help_history end end]

    set url [lrange $help_history $l $l]
    incr l -1
    set help_history [lrange $help_history 0 $l]

    HMlink_callback $win $url keep
}

proc help_forward_link {win} {
    global help_future

    set l [expr [llength $help_future]-1]
    if {$l < 0} {
	bell
	return
    }

    set url [lrange $help_future $l $l]
    incr l -1
    set help_future [lrange $help_future 0 $l]
    HMlink_callback $win $url keep
}

proc help_redraw_link {win} {
    global help_history help_file

    set l [expr [llength $help_history]-1]
    set url [lrange $help_history $l $l]
    incr l -1
    set help_history [lrange $help_history 0 $l]

    set help_file ""
    HMlink_callback $win $url keep
}

proc help_set_status {text} {
    global help_toplevel
    $help_toplevel.top.status delete 0 end
    $help_toplevel.top.status insert end $text
}

proc help_set_stop {mode} {
    global help_toplevel
    if [winfo exist $help_toplevel] {
        $help_toplevel.top.stop configure -state $mode
    }
}

proc help_history_trace {name element op} {
    global $name help_toplevel
    if {[llength [set $name]] < 2} {
	$help_toplevel.top.back configure -state disabled
    } else {
	$help_toplevel.top.back configure -state normal
    }
}

proc help_future_trace {name element op} {
    global $name help_toplevel
    if {[set $name] == ""} {
	$help_toplevel.top.forw configure -state disabled
    } else {
	$help_toplevel.top.forw configure -state normal
    }
}

#-----------------------------------------------------------------------------
# html-library callbacks

# Wrap these up in a procedure so we can redefine them again
# regardless of which order autoloading takes place.
proc HM_our_init {} {

# Stop and start parsing of html. These should be done using
# "HMset_state $win -stop 1" (or -stop 0), but this doesn't appear to work.
# For the time being we're hacking directly into the library.
proc HMstop_parse {win} {
    upvar #0 HM$win var
    set var(stop) 1
}

proc HMstart_parse {win} {
    upvar #0 HM$win var
    set var(stop) 0
}

# Callback for when an image is loaded.
proc HMset_image {win handle src} {
    update idletasks
    global TRANSPARENT_GIF_COLOR help_path
    set TRANSPARENT_GIF_COLOR [$win cget -bg]

    image create photo $src -file $help_path/$src -palette 7/7/4
    HMgot_image $handle $src
}

# Callback for when a link is followed. Complicated greatly by the fact that
# it needs to be reentrant for the case where people follow a link before
# the displaying of the current link has finished.
proc HMlink_callback {win link {future clear}} {
    global help_file help_recurse help_history help_future
    global MAX_LINK_HISTORY help_toplevel

    if [string match *:* $link] {
	set x [string first : $link]
	set service [string range $link 0 [expr $x-1]]
	set link [string range $link [expr $x+1] end]

	if {$service != "file"} {
#	    tk_dialog $help_toplevel.error "Error" \
#		"Couldn't follow link '$link'. \
#		Only file: based URLs are supported" \
#		error 0 Ok
	    bell
	    return
	}
    }

    if {[string index $link 0] != "/"} {
	global help_path
	set link "$help_path/$link"
    }

    if [string match *#* $link] {
	set x [string first # $link]
	set file [string range $link 0 [expr $x-1]]
	set mark [string range $link [expr $x+1] end]
    } else {
	set file $link
	set mark ""
    }

    update idletasks

    # Add to history. Chop off oldest element if we're getting too long.
    set old_hist $help_history
    set old_future $help_future
    lappend help_history $link
    if {$future == "clear"} {
	set help_future ""
    }
    if {[llength $help_history] > $MAX_LINK_HISTORY+1} {
	set help_history [lrange $help_history 1 end]
    }

    # If we've got a new file, load and render it
    if {$file != $help_file} {
	# Read the file into $html
        if [catch {set fd [open $file]}] {
	    set help_history $old_hist
	    set help_future $old_future
	    tk_dialog $help_toplevel.error "Error" \
		"Couldn't follow link '$link'. The file doesn't exist" \
		error 0 Ok
	    return
	}
        set html [read $fd]
        close $fd

        HMreset_win $win

	# Flag our loop
	incr help_recurse
	set help_file $file

	# Goto the correct position (which is allowed to be done before
	# it's visible) and then render.
	if {$mark != ""} {
	    HMgoto $win $mark
	}
	help_set_status $link

	help_set_stop normal
       	HMparse_html $html "HMrender $win"

	# This bit's horrid. We need to check whether we've been told to
	# redraw before finishing the current rendering. In this case this
	# routine will have been called again, the HMparse_html will have
	# finished, and recurse will be > 1. Hence we then disable any
	# future rendering thus causing the HMparse_html call of higher stacks
	# (which are still present) to exit. When the original HMparse_html has
	# exited we must reenable the rendering once more. Sorry if it seems
	# complicated, but it is.
	if {$help_recurse > 1} {
	    # stopping rendering for higher stacks
	    HMstop_parse $win
	}
	if {[incr help_recurse -1] == 0} {
	    # reenabling rendering when original call of this routine ends.
	    HMstart_parse $win
	    help_set_stop disabled
	}

    # Otherwise, (ie when we've ONLY changed mark) perform a goto request.
    } elseif {$mark != ""} {
        help_set_status $link
    	HMgoto $win $mark
    } else {
	$win yview 1.0
    }
}

# from "proc HM_our_init {} ..."
}

# Called after a goto; normally flashes orange, but we've disabled it.
proc HMwent_to {win where} {
}

#-----------------------------------------------------------------------------
# Main startup code

# When running as a standalone program, we bring up the display automatically
if {[info exists argv] && $argv == "-init"} {
    set auto_path "$env(STADLIB)/tkutils $auto_path"
    tkinit
    help_init "" 1
    cd $help_path
}
#HMlink_callback [help_init] registration_toc.html
#HMlink_callback [help_init] test.html
