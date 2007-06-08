# -----------------------------------------------------------------------------
#
# Creates a modal top-level window.
#
# The arguments are as per the toplevel command, except:
#	?-resizable boolean?
#		Specifies whether this window may be resized. Defaults to 0
#
# At present, this option needs to be at the start of the argument list.
#
# Success: returns the pathname of the window 
# Failure: returns "" (if the window already exists)
#
namespace eval modal {}
proc modal {w args} {
    global tk_utils_defs

    # Check if the window already exists. Return "" if so.
    if {[winfo exists $w]} {
	wm deiconify $w
	raise $w
	return ""
    }

    # Parse "-resizable" boolean option. In the current implementation this
    # has to be the first argument.
    set resizable 0
    if {[lindex $args 0] == "-resizable"} {
	set resizable [lindex $args 1]
	set args [lreplace $args 0 1]
    }

    # Create the window
    eval toplevel $w $args
    wm group $w [winfo parent $w]

    if {$resizable == 0} {
	wm resizable $w 0 0
    }

    # Take focus and grab, but remember old ones first.
    catch {unset ::modal::$w}
    upvar #0 ::modal::$w data
    set data(oldGrab) [grab current $w]
    if {$data(oldGrab) != ""} {
	set data(grabStatus) [grab status $data(oldGrab)]
    }
    set data(oldFocus) [focus]
    if {[catch {grab $w}]} {
	after idle "catch {grab $w}"
    }
    focus $w

    # Add Modal bindings, to deal with destroy and tidyup.
    bindtags $w [list $w Modal Toplevel all]

    # Move window to offset from parent.
    if {[keylget tk_utils_defs PLACE_WINDOWS] && \
	    [set p [winfo parent $w]] != ""} {
	set rootx [winfo rootx [winfo toplevel $p]]
	set rooty [winfo rooty [winfo toplevel $p]]
	incr rootx 20
	incr rooty 20
	
	if {$rootx < 0} {set rootx 0}
	if {$rooty < 0} {set rooty 0}

	wm geometry $w +$rootx+$rooty
    }

    return $w
}

# Modal Destroy binding.
# This restores the saved grab and focus.
bind Modal <Destroy> {
    upvar #0 ::modal::%W data
    grab release %W
    if {[string compare "$data(oldGrab)" ""] != 0} {
	if {[string compare "$data(grabStatus)" "global"] == 0} {
	    grab -global $data(oldGrab)
	} else {
	    grab $data(oldGrab)
	}
    }
    catch {focus $data(oldFocus)}
    unset data
}


# -----------------------------------------------------------------------------
#
# Creates a top-level window with the focus automatically set.
# The window is not modal (ie no grab) and is resizable.
#
# The arguments are as per the toplevel command, except:
#	?-resizable boolean?
#		Specifies whether this window may be resized. Defaults to 0
#       ?-focus boolean?
#		Specifies whether to take the focus. Defaults to 1
#
# Success: returns the pathname of the window 
# Failure: returns "" (if the window already exists)
#
namespace eval xtoplevel {}
proc xtoplevel {w args} {
    global tk_utils_defs
    set take_focus 1

    # Check if the window already exists. Return "" if so.
    if {[winfo exists $w]} {
	wm deiconify $w
	raise $w
	return ""
    }

    # Parse "-resizable" boolean option. In the current implementation this
    # has to be the first argument.
    set resizable 1
    if {[lindex $args 0] == "-resizable"} {
	set resizable [lindex $args 1]
	set args [lreplace $args 0 1]
    }
    if {[lindex $args 0] == "-focus"} {
	set take_focus [lindex $args 1]
	set args [lreplace $args 0 1]
    }

    # Create the window
    eval toplevel $w $args
    wm group $w [winfo parent $w]

    if {$resizable == 0} {
	wm resizable $w 0 0
    }
    # Take focus
    catch {unset ::xtoplevel::$w}
    upvar #0 ::xtoplevel::$w data
    if {$take_focus} {
	set data(oldFocus) [focus]
	focus $w
    }

    # Add bindings, to deal with destroy and tidyup.
    bindtags $w [list $w XToplevel Toplevel all]

    # Move window to offset from parent.
    if {[keylget tk_utils_defs PLACE_WINDOWS] && \
	    [set p [winfo parent $w]] != ""} {
	set rootx [winfo rootx [winfo toplevel $p]]
	set rooty [winfo rooty [winfo toplevel $p]]
	incr rootx 20
	incr rooty 20

	if {$rootx < 0} {set rootx 0}
	if {$rooty < 0} {set rooty 0}

	wm geometry $w +$rootx+$rooty
    }

    return $w
}

# XToplevel Destroy binding.
# This restores the saved focus.
bind XToplevel <Destroy> {
    upvar #0 ::xtoplevel::%W data
    catch {focus $data(oldFocus)}
    catch {unset data}

#    # Simulates the windows focus bugstash
#    if {[wm state %W] != "withdrawn"} {
#	raise [winfo parent %W]
#    }
}


# -----------------------------------------------------------------------------
#
# Test code
# 
# button .b -text foo -command {puts [info globals];exit}
# pack .b
# modal .w
# button .w.b -text foo -command "com1"
# button .w.b2 -text foo2 -command "destroy .w"
# pack .w.b .w.b2
# 
# proc com1 {} {
#     modal .w.w -resizable 1
#     button .w.w.b -text bar -command "destroy .w.w"
#     pack .w.w.b
# }
