#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
##############################################################################
# displays a stack dump for tcl
proc stack_dump {} {
    puts "ERROR!!! - Tcl stackframe follows"
    for {set i [info level]} {$i > 0} {incr i -1} {
        puts "Level $i: [info level $i]"
    }
}

##############################################################################
# checks if an input is an integer
proc isinteger { value } {
  return [regexp {^[+-]?[0-9]+$} $value]
}

##############################################################################
# checks if an input is an float
proc isfloat { value } {
    return [regexp {^[+-]?[0-9]*(\.[0-9]*)?([Ee][+-]?[0-9]+)?$} $value]
}

##############################################################################
#Set busy mode
proc InitBusy {main menu name} {
    global busy_main_path busy_menu_path busy_menu_name
    set busy_main_path $main
    set busy_menu_path $menu
    set busy_menu_name $name
}

proc SetBusy {} {
    global busy_main_path busy_menu_path busy_menu_name busy_menu_state

    set busy_menu_state [menu_state_save $busy_menu_path $busy_menu_name]
    menu_state_set $busy_menu_name -2 $busy_menu_path

    foreach win "[winfo children .]" {
	if {$win != "$busy_main_path"} {
            catch {$win configure -cursor watch}
	}
    }

    grab $busy_main_path
}

proc ClearBusy {} {
    global busy_main_path busy_menu_path busy_menu_name busy_menu_state

    menu_state_restore $busy_menu_path $busy_menu_name $busy_menu_state

    foreach win "[winfo children .]" {
    	catch {$win configure -cursor top_left_arrow}
    }

    grab release $busy_main_path
}

##############################################################################
#Creates a popup menu
proc create_popup {w title} {
    if {[winfo exists $w]} {destroy $w}
    menu $w -tearoff 0 -disabledforeground blue
    set bg [lindex [$w configure -bg] 4]
    $w add command -state disabled -label "$title" \
        -background [tk::Darken $bg 80] \
        -font menu_title_font

    return $w
}

##############################################################################
#two functions to set and get the "current frame" - useful for instances
#when you have two frames and you wish to differentiate between them by,
#say clicking in one
proc SetCurFrame {s frame} {
    global $s.frame $s.frame_index
    set $s.frame_index 0
    set $s.frame $frame
}

proc GetCurFrame {s} {
    global $s.frame $s.frame_index

    set f [lindex [set $s.frame] [set $s.frame_index]]
    incr $s.frame_index
    if {[set $s.frame_index] >= [llength [set $s.frame]]} {
	set $s.frame_index 0
    }
    return $f
}

##############################################################################
#deletes a file with error checking
proc DeleteFile { file } {

    catch {file delete $file} e

    if {$e != ""} {
	tk_messageBox -icon error -type ok -title "Delete file" \
		-message $e
    }
}

##############################################################################
# Fixes the maximum size of a toplevel window to take into account screen
# borders, such as the Windows task bar or a CDE desktop.
proc fix_maxsize {w} {
    global tk_utils_defs

    set border_x [keylget tk_utils_defs X_BORDER_SIZE]
    set border_y [keylget tk_utils_defs Y_BORDER_SIZE]

    foreach {width height} [wm maxsize $w] {}

    if {$width > [winfo screenwidth .]} {
	set width [winfo screenwidth .]
    }

    if {$height > [winfo screenheight .]} {
	set height [winfo screenheight .]
    }
    
    incr width -$border_x
    incr height -$border_y

    wm maxsize $w $width $height
}

##############################################################################
# Fixes the maximum size of a toplevel window that contains a gridded text
# window (so needs character coords) to take into account screen borders,
# such as the Windows task bar or a CDE desktop.
proc fix_maxsize_text {w font_width font_height extra_width extra_height} {
        global tk_utils_defs

        set border_x [keylget tk_utils_defs X_BORDER_SIZE]
        set border_y [keylget tk_utils_defs Y_BORDER_SIZE]

        set width [winfo screenwidth .]
        set height [winfo screenheight .]

        incr width -$border_x
        incr height -$border_y

        set width [expr ($width - $extra_width) / $font_width]
        set height [expr ($height - $extra_height) / $font_height]

        wm maxsize $w $width $height
}

#
# Force window size using wm geometry. This is needed in addition to
# fix_maxsize as on some window managers (AfterStep, MacOS X, etc) the
# wm maxsize command is ignored.
#
proc fit_on_screen2 {w} {
    global tk_utils_defs

    puts fit_on_screen2

    #10.10.02 (added but commented out - see fit_on_screen comment below)
    #wm geometry $w {}

    set border_x [keylget tk_utils_defs X_BORDER_SIZE]
    set border_y [keylget tk_utils_defs Y_BORDER_SIZE]

    foreach {width height} [wm maxsize $w] {}

    if {$width > [winfo screenwidth .]} {
        set width [winfo screenwidth .]
    }

    if {$height > [winfo screenheight .]} {
        set height [winfo screenheight .]
    }
    
    incr width -$border_x
    incr height -$border_y

    update idletasks

    set wid [lindex [split [wm geometry $w] x+] 0]
    set hei [lindex [split [wm geometry $w] x+] 1]

    if {$wid > $width} {
        set wid $width
    }

    if {$hei > $height} {
        set hei $height
    }

    wm geometry $w ${wid}x$hei
}

proc fit_on_screen {w} {
    # FIXME: MacOS X hack to deal with ignoring wm maxsize. This fixed
    # delay may still cause problems on slow macs, but this code will
    # be rewritten once the container class has been implemented.

    #after 1000 "catch {fit_on_screen2 $w}"

    #kfs/jkb 10.10.02 
    #fit_on_screen currently causes problems - especially the 1 second delay
    #which makes bringing up lots of plots at the same time (eg codon pref)
    #very slow. Also, if you bring up 2 comparison plots, separate them out 
    #and then superimpose them again, the new plot does not shrink in size as
    #it should. Tried adding a wm geometry $w {} (see above) which solved this
    #problem but it also lost manual resizing information. 
    #We think fit_on_screen's only purpose  was to solve resizing issues on
    #on the mac, specifically, using wm geometry to force resizing windows 
    #when they grow too large for the screen. Ideally wm maxsize will solve 
    #this, but apparently not on all window managers.
    #Since we are about to upgrade the mac, we need to see if this is still 
    #necessary.
    return
}

#
# Implements a "do <script> ??until|while? <expression>?" loop
#
# It is as fast as builtin "while" command for loops with
# more than just a few iterations.
#
# From http://mini.net/tcl/917.html
#
proc do {script {arg2 {}} {arg3 {}}} {
    if {![string length $arg2$arg3]} {set arg2 0}

    if {[string compare $arg3 {}]} {
        switch -- $arg2 {
	    until   {set bool "!"}
	    while   {set bool ""}
	    default {return -code 1 "usage: do script ??until|while? expr?"}
        }
    }

    set ret [catch { uplevel $script } result]
    switch $ret {
        0 -
        4 {}
        3 {return}
        default {
            return -code $ret $result
        }
    }

    set ret [catch {uplevel [list while ${bool}($arg3) $script]} result]
    return -code $ret $result
}

#
# Implements a tmpnam function. Prefix is optional, but if set then it
# is the start of the temporary filename (excluding the directory portion).
#
proc tmpnam {{prefix tmp}} {
    global tcl_platform env

    if { "$tcl_platform(platform)" != "windows" } {
	set tdir "/tmp/"
    } else {
	if {[info exists env(TMP)]} {
	    set tdir $env(TMP)/
	} elseif {[info exists env(TEMP)]} {
	    set tdir $env(TEMP)/
	} else {
	    set tdir "C:/"
	}
	regsub -all {\\} $tdir / tdir
    }

    set pid [pid]
    set count -1
    do {
	incr count
	set fname "${tdir}${prefix}${pid}_$count"
    } while {[file exists $fname]}

    return $fname
}

