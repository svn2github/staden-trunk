# search.tcl --
#
# This demonstration script creates a collection of widgets that
# allow you to load a file into a text widget, then perform searches
# on that file.
#
# @(#) search.tcl 1.2 95/06/21 09:17:11

# textSearchPos --
# Search for all instances of a given string in a text widget and
# apply a given tag to each instance found and move the cursor to the
# position
#
# Arguments:
# w -		The window in which to search.  Must be a text widget.
# string -	The string to search for.  The search is done using
#		exact matching only;  no special characters.
# tag -		Tag to apply to each instance of a matching string.

proc textSearchPos {w string tag direction} {
    global prev_cur cur 

    #set position of search to the current position of the cursor
    if {![info exists cur]} {
	set prev_cur [$w index insert]
    }
    set cur [$w index insert]

    if {$string == ""} {
	return
    }

    if {$direction == 1} {
	set cur [$w index "$cur + 1 char"]
	set cur [$w search -forwards -nocase -count length -- $string $cur end]
    } else { 
	set cur [$w index "$cur - 1 char"]
	set cur [$w search -backwards -nocase -count length -- $string $cur 1.0]
    }
    if {$cur == ""} {
	bell
	set cur $prev_cur
	return
    } else {
	set prev_cur $cur
    }
    $w see $cur
    $w tag remove search 0.0 end
    $w tag add $tag $cur "$cur + $length char"
    $w mark set insert $cur
}

# textSearch --
# Search for all instances of a given string in a text widget and
# apply a given tag to each instance found.
#
# Arguments:
# w -		The window in which to search.  Must be a text widget.
# string -	The string to search for.  The search is done using
#		exact matching only;  no special characters.
# tag -		Tag to apply to each instance of a matching string.

proc textSearch {w string tag} {
    $w tag remove search 0.0 end
    if {$string == ""} {
	return
    }
    set cur 1.0
    while 1 {
	set cur [$w search -count length -- $string $cur end]

	if {$cur == ""} {
	    break
	}
	$w tag add $tag $cur "$cur + $length char"
	set cur [$w index "$cur + $length char"]
    }
}

# textToggle --
# This procedure is invoked repeatedly to invoke two commands at
# periodic intervals.  It normally reschedules itself after each
# execution but if an error occurs (e.g. because the window was
# deleted) then it doesn't reschedule itself.
#
# Arguments:
# cmd1 -	Command to execute when procedure is called.
# sleep1 -	Ms to sleep after executing cmd1 before executing cmd2.
# cmd2 -	Command to execute in the *next* invocation of this
#		procedure.
# sleep2 -	Ms to sleep after executing cmd2 before executing cmd1 again.

proc textToggle {cmd1 sleep1 cmd2 sleep2} {
    catch {
	eval $cmd1
	after $sleep1 [list textToggle $cmd2 $sleep2 $cmd1 $sleep1]
    }
}
