catch {font create tooltip_font    -family Helvetica -size -12}

namespace eval ::tooltip {
    variable after_id ""
    variable tooltip_path ""
    variable pointerx -1
    variable pointery -1
    variable enabled 1

    # Attaches tooltip message "msg" to windor "win"
    proc attach {win msg} {
	variable old_enter [bind $win <Any-Enter>]
	variable old_leave [bind $win <Any-Leave>]
	# FIX to speed up "font measure" command
	if {![winfo exists .::tooltip::label]} {
	    label .::tooltip::label -font tooltip_font
	}
	set msg [wrap $msg 300 tooltip_font]
        bind $win <Any-Enter> +[namespace code [list enter $win $msg]]
        bind $win <Any-Leave> +[namespace code [list leave $win]]

	# Edit/'use' events will cancel the tooltip by simulating a leave
	set wtop [winfo toplevel $win]
	variable old_button [bind $wtop <Any-ButtonPress>]
	variable old_key [bind $wtop <Any-KeyPress>]
	bind $wtop <Any-Button> +[namespace code [list leave $win now]]
	bind $wtop <Any-KeyPress> +[namespace code [list leave $win now]]
    }

    # Detaches any tooltips from "win"
    # NOTE: If window or toplevel bindings have changed since the attach was
    # called then this will overwrite them with old copies.
    proc detach {win} {
	variable old_enter
	variable old_leave
	bind $win <Any-Enter> $old_enter
	bind $win <Any-Leave> $old_leave

	variable old_button
	variable old_key
	set wtop [winfo toplevel $win]
	bind $win <Any-Button> $old_button
	bind $win <Any-KeyPress> $old_key
    }

    # Enables tooltips
    proc enable {} {
	variable enabled
	set enabled 1
    }

    # Disables tooltips
    proc disable {} {
	variable enabled
	set enabled 0
    }

    # Entry binding to the dialogue component; brings up the tooltip after a
    # certain length of time
    proc enter {win msg} {
	variable enabled
	if {!$enabled} return

	# Don't redisplay the tip unless we have moved.
	variable pointerx
	variable pointery
	if {[winfo pointerx $win] == $pointerx &&
	    [winfo pointery $win] == $pointery} {
	    return
	}

	# If a tip is already present then bring up the new tip immediately
	variable after_id
	variable tooltip_path
	if {$tooltip_path != "" && [winfo exists $tooltip_path]} {
	    after cancel "catch {destroy $tooltip_path}"
	    if {$tooltip_path != "$win.tooltip"} {
		destroy $tooltip_path
		popup $win $msg
	    }
	} else {
	    set after_id [after 1500 [namespace code [list popup $win $msg]]]
	}
    }

    # Leave binding from the dialogue component. Cancel the request to bring
    # up a tooltip and schedule a request to destroy any existing tooltip.
    # (This scheduled event will be subsequently removed if the leave event
    # is due to an entry invent to the tooltip window itself.)
    proc leave {win {when 100}} {
	variable after_id
	if {$after_id != ""} {
	    after cancel $after_id
	    if {$when == "now"} {
		catch {destroy $win.tooltip}
	    } else {
		after $when "catch {destroy $win.tooltip}"
	    }
	}
    }

    # Creates the tooltip window itself. Called from the Entry binding.
    proc popup {win msg} {
	variable tooltip_path

	# Create a toplevel window
	set w $win.tooltip
	if {[winfo exists $w]} {destroy $w}
	toplevel $w
	bind $w <1> "
	    set [namespace current]::pointerx \[winfo pointerx %W\]
	    set [namespace current]::pointery \[winfo pointery %W\]
	    destroy $w"
	set tooltip_path $w
	wm withdraw $w
	wm overrideredirect $w 1

	set padx 5
	set bd 1

	# Work out coordinates, line-wrapping by words first
	set x [expr {[winfo pointerx $win]+5}]
	set y [expr {[winfo pointery $win]+15}]
	set xdim 0
	#set msg [wrap $msg 300 tooltip_font]
	set lines 0
	foreach line [split $msg "\n"] {
	    set dim [font measure tooltip_font $line]
	    if {$dim > $xdim} {
		set xdim $dim
	    }
	    incr lines
	}
	incr xdim [expr {2*($padx+$bd)}]
	array set metrics [font metrics tooltip_font]
	set ydim [expr {$lines*$metrics(-linespace) + 2*$bd+2}]

	# Ensure the window is entirely on the screen
	if {[expr {$x+$xdim}] > [winfo screenwidth $w]} {
	    set x [expr {[winfo screenwidth $w]-$xdim}]
	}
	if {[expr {$y+$ydim}] > [winfo screenheight $w]} {
	    set y [expr {[winfo pointery $win]-$ydim-10}]
	}
	wm geometry $w +$x+$y

	label $w.l \
	    -bg lightyellow \
	    -justify left \
	    -padx $padx \
	    -bd $bd -relief solid \
	    -font tooltip_font \
	    -text $msg
	pack $w.l -fill both -expand 1

	# Make it visible
	wm deiconify $w

	# Bindings for this tooltip. Cancel any scheduled deletion of this
	# window on entry, and schedule a destroy on exit.
	bind $w <Any-Enter> {after cancel {catch {destroy %W}}}
	bind $w <Any-Leave> "
 	    if {\"%W\" != \"$w\"} continue
	    after 100 {catch {destroy %W}}"
    }

    # Takes a string 'str' and word wraps it to a line of 'wrap_len' pixels
    # long using 'font'.
    # The new wrapped string is returned.
    proc wrap {str wrap_len font} {
	set xwid 0
	set new_str ""
	set first ""
	
	# Split first on space and tabs
	foreach w2 [split $str " \t"] {
	    set nl 0
	    # And then split on newlines, as we need to reset counters for
	    # these.
	    foreach word [split $w2 "\n"] {
		# newline => new width
		if {$nl} {
		    append newstr "\n"
		    set first ""
		    set xwid 0
		}
		
		# Increment width, and wrap if necessary
		incr xwid [font measure $font "$first$word"]
		if {$xwid > $wrap_len} {
		    append newstr "\n"
		    set xwid [font measure $font "$word"]
		    set first ""
		}
		append newstr "$first$word"
		if {$first == ""} {
		    set first " "
		}
		incr nl
	    }
	}

	return $newstr
    }
}
