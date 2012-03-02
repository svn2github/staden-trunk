#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#
# Adds text output windows - called from CreateMain
#
proc XXConsole {} {
    set w .xxconsole
    toplevel $w
    wm title $w "Tcl console window"
    pack [text $w.t] -fill both -expand 1
    $w.t insert end "% "
    bind $w.t <Return> {
        set __line [%W get "current linestart" "current lineend"]
        regsub "^%%\[ \n\]*" $__line {} __line
        %W insert end "\n"
        catch {uplevel #0 eval $__line} __result
        %W insert end "$__result\n"
        %W insert end "%% "
    }
}

proc tout_create_wins {f {width 81}} {
    global tk_utils_defs
    global $f.stdout.Scr $f.stderr.Scr $f.stdout.Redir $f.stderr.Redir \
	   $f.stderr.Bell

    #
    # The output window section
    #
    frame $f.stdout -bd 2 -relief raised
    frame $f.stderr -bd 2 -relief raised

    frame $f.stdout.bar
    frame $f.stderr.bar
    frame $f.stdout.x
    frame $f.stderr.x
    frame $f.stdout.x.padding
    frame $f.stderr.x.padding

    # The button bars
    label $f.stdout.bar.label -text "Output window:"
    label $f.stdout.bar.redirl -textvariable $f.stdout.Redir \
	-fg blue
    label $f.stderr.bar.label -text "Error window:"
    label $f.stderr.bar.redirl -textvariable $f.stderr.Redir \
	-fg blue

    xmenubutton $f.stdout.bar.redir -text "Redirect >>" \
	-menu $f.stdout.bar.redir.m
    xmenubutton $f.stderr.bar.redir -text "Redirect >>" \
	-menu $f.stderr.bar.redir.m
    menu $f.stdout.bar.redir.m -tearoff 0
    $f.stdout.bar.redir.m add command -label "Open"\
       -command "tout_redir_open $f.stdout.Redir stdout $f.stdout.bar.redir.m"
    $f.stdout.bar.redir.m add command -label "Close" -state disabled \
       -command "tout_redir_close $f.stdout.Redir stdout $f.stdout.bar.redir.m"
    menu $f.stderr.bar.redir.m -tearoff 0
    $f.stderr.bar.redir.m add command -label "Open" \
       -command "tout_redir_open $f.stderr.Redir stderr $f.stderr.bar.redir.m"
    $f.stderr.bar.redir.m add command -label "Close" -state disabled  \
       -command "tout_redir_close $f.stderr.Redir stderr $f.stderr.bar.redir.m"

    button $f.stdout.bar.clear -text "Clear" \
	-command "tout_clear $f.stdout.t"
    button $f.stderr.bar.clear -text "Clear" \
	-command "tout_clear $f.stderr.t"
    checkbutton $f.stdout.bar.scroll -text "Scroll on output " \
	-variable $f.stdout.Scr \
	-command "tout_scroll stdout \[set $f.stdout.Scr\]" \
	-bd 2 -relief raised -padx 4
    checkbutton $f.stderr.bar.scroll -text "Scroll on output" \
	-variable $f.stderr.Scr \
	-command "tout_scroll stderr \[set $f.stderr.Scr\]" \
	-bd 2 -relief raised -padx 4

    button $f.stdout.bar.search -text "Search" \
	-command "tout_search $f.stdout.t"
    button $f.stderr.bar.search -text "Search" \
	-command "tout_search $f.stderr.t"

    checkbutton $f.stderr.bar.bell -text "Bell" \
	-variable $f.stderr.Bell \
	-command "error_bell \[set $f.stderr.Bell\]" \
	-bd 2 -relief raised -padx 4

    set $f.stdout.Scr [keylget tk_utils_defs OUTPUT_SCROLL]
    set $f.stderr.Scr [keylget tk_utils_defs OUTPUT_SCROLL]
    set $f.stderr.Bell [keylget tk_utils_defs ERROR_BELL]
    error_bell [set $f.stderr.Bell]
    tout_scroll stdout [set $f.stdout.Scr]
    tout_scroll stderr [set $f.stderr.Scr]

    # The text displays
    text $f.stdout.t -width $width -height 16 -wrap none \
	-yscrollcommand "$f.stdout.y set" \
	-xscrollcommand "$f.stdout.x.sb set"
    text $f.stderr.t -width $width -height 4  -wrap none \
	-yscrollcommand "$f.stderr.y set" \
	-xscrollcommand "$f.stderr.x.sb set"

    if {[winfo screenwidth $f] < 1024 || \
            [winfo screenheight $f] < 768} {
	$f.stderr.t tag configure error -font title_font
        $f.stderr.t insert end \
            "WARNING: This program operates best at resolutions of 1024x768\n\
            and higher.\n" error
    }

    # The scrollbars
    scrollbar $f.stdout.y -orient vert  -command "$f.stdout.t yview"
    scrollbar $f.stderr.y -orient vert  -command "$f.stderr.t yview"
    scrollbar $f.stdout.x.sb -orient horiz -command "$f.stdout.t xview"
    scrollbar $f.stderr.x.sb -orient horiz -command "$f.stderr.t xview"

    # Packing
    pack $f.stdout -fill both -expand 1
    pack $f.stderr -fill both
    pack $f.stdout.bar $f.stderr.bar -side top -fill x
    pack $f.stdout.bar.label $f.stdout.bar.redirl -side left -fill both
    pack $f.stderr.bar.label $f.stderr.bar.redirl -side left -fill both
    pack $f.stdout.bar.redir $f.stdout.bar.clear $f.stdout.bar.scroll $f.stdout.bar.search\
	-side right -fill both
    pack $f.stderr.bar.redir $f.stderr.bar.clear $f.stderr.bar.scroll $f.stderr.bar.bell $f.stderr.bar.search\
	-side right -fill both
    pack $f.stdout.x $f.stderr.x -side bottom -fill x
    pack $f.stdout.y $f.stderr.y -side right -fill y
    pack $f.stdout.t $f.stderr.t -fill both -expand 1
    pack $f.stdout.x.sb $f.stderr.x.sb -side left -fill x -expand 1
    pack $f.stdout.x.padding $f.stderr.x.padding  -side right

    # Here come those padding hacks again...
    pack propagate $f.stdout.x.padding 0
    pack propagate $f.stderr.x.padding 0
    $f.stdout.x.padding configure -width 22
    $f.stderr.x.padding configure -width 22
#    update idletasks
#    $f.stdout.x.padding configure \
#	-width  [winfo width $f.stdout.y] \
#	-height [winfo width $f.stdout.y]
#    $f.stderr.x.padding configure \
#	-width  [winfo width $f.stderr.y] \
#	-height [winfo width $f.stderr.y]

    tout_init $f.stdout.t $f.stderr.t

    bind $f.stderr.t <c><o><n> "console show"
}


proc tout_tag_params {w tag params} {
    global param_$tag

    set param_$tag $params

}

#
# Adds a new header to the output window
#
proc tout_new_header {w tag name} {

    $w tag bind ${tag}_h <<menu>> "tout_popup $w $tag %X %Y"
    $w tag bind ${tag}_t <<menu>> "tout_popup $w $tag %X %Y"
    $w tag bind ${tag}_p <<menu>> "tout_popup $w $tag %X %Y"
    set bg [lindex [$w configure -bg] 4]
    $w tag configure ${tag}_h -background [tk::Darken $bg 80]
    $w tag configure ${tag}_p -background [tk::Darken $bg 60] -wrap word

}


#
# Pops up a menu for the text output display
#
proc tout_popup {w tag X Y} {
    global tcl_platform .tout_popup
    if {[winfo exists $w.m]} {destroy $w.m}
     
    create_popup $w.m Commands
    if {![info exists .tout_popup]} {
	set ".tout_popup(Show input parameters)" 1
	set ".tout_popup(Print)" 1
	set ".tout_popup(Output to disk)" 1
	set ".tout_popup(Output to list)" 1
	set ".tout_popup(Output to command)" 1
    }

    if {[info exists ".tout_popup(Show input parameters)"]} {
	$w.m add command -label "Show input parameters" \
	    -command "destroy $w.m; tout_in_params $w $tag"
    }
    if {$tcl_platform(platform) == "unix" && \
	    [info exists ".tout_popup(Print)"]} {
	$w.m add command -label "Print" \
	    -command "destroy $w.m; tout_print $w $tag"
    }
    if {[info exists ".tout_popup(Output to disk)"]} {
	$w.m add command -label "Output to disk" \
	    -command "destroy $w.m; tout_output_disk $w $tag"
    }
    if {[info exists ".tout_popup(Output to list)"]} {
	$w.m add command -label "Output to list" \
	    -command "destroy $w.m; tout_output_list $w $tag"
    }
    if {$tcl_platform(platform) == "unix" && \
	    [info exists ".tout_popup(Output to command)"]} {
	$w.m add command -label "Output to command" \
	    -command "destroy $w.m; tout_output_command $w $tag"
    }
    $w.m add separator
    $w.m add command -label "Remove" \
	-command "destroy $w.m; tout_remove $w $tag"

    tk_popup $w.m [expr $X-20] [expr $Y-10]
}


#
# Inserts input parameters associated with header
#
proc tout_in_params {w tag} {
    global param_$tag

    if {![info exists param_$tag]} {
	return
    }
    set end [lindex [$w tag ranges ${tag}_h] 1]

    #param_$tag is the the input parameter text with tag ${tag}_p
    $w insert $end [set param_$tag] ${tag}_p 
}

#
# Removes some text from the output display
#
proc tout_remove {w tag} {
    global param_$tag

    if {[set r [$w tag ranges ${tag}_h]] != ""} {
        eval $w delete $r
    }
    if {[set r [$w tag ranges ${tag}_t]] != ""} {
        eval $w delete $r
    }
    if {[set r [$w tag ranges ${tag}_p]] != ""} {
        eval $w delete $r
    }
    $w tag delete ${tag}_h
    $w tag delete ${tag}_t
    $w tag delete ${tag}_p

    if {[info exists param_$tag]} {
	unset param_$tag
    }
}


proc tout_output_checks { t } {
    global tk_utils_defs $t.b1 $t.b2 $t.b3

    frame $t
    set b1 [keylget tk_utils_defs TEXT_OUTPUT.OUTPUT.1.NAME]
    set b2 [keylget tk_utils_defs TEXT_OUTPUT.OUTPUT.2.NAME]
    set b3 [keylget tk_utils_defs TEXT_OUTPUT.OUTPUT.3.NAME]

    set $t.b1 [keylget tk_utils_defs TEXT_OUTPUT.OUTPUT.1.VALUE]
    set $t.b2 [keylget tk_utils_defs TEXT_OUTPUT.OUTPUT.2.VALUE]
    set $t.b3 [keylget tk_utils_defs TEXT_OUTPUT.OUTPUT.3.VALUE]

    checkbutton $t.b1 -text $b1 -variable $t.b1
    checkbutton $t.b2 -text $b2 -variable $t.b2
    checkbutton $t.b3 -text $b3 -variable $t.b3
	
    pack $t.b1 $t.b2 $t.b3 -anchor w
}

proc tout_output_get_checks {w t tag output } {
    global $t.Header $t.Params $t.Text
    global $output.b1 $output.b2 $output.b3

    set $t.Header ""
    set $t.Params ""
    set $t.Text ""

    if {[set $output.b1]} {
	if {[set r [$w tag ranges ${tag}_h]] != ""} {
	    set $t.Header [eval $w get $r]
	}
    }
    if {[set $output.b2]} {
	if {[set r [$w tag ranges ${tag}_p]] != ""} {
	    set $t.Params [eval $w get $r]
	}
    }
    if {[set $output.b3]} {
	if {[set r [$w tag ranges ${tag}_t]] != ""} {
	    set $t.Text [eval $w get $r]
	}
    }
}

#
# Prints a section of output
#
proc tout_print {w tag} {
    global tk_utils_defs 
    set t [keylget tk_utils_defs TEXT_OUTPUT.WIN]
    xtoplevel $t -resizable 0
    wm title $t "Print"

    tout_output_checks $t.output

    okcancelhelp $t.ok \
        -ok_command "tout_print2 $w $t $tag $t.output" \
        -cancel_command "destroy $t" \
	-help_command "show_help interface {UI-Output}"\
        -bd 2 -relief groove

    pack $t.output $t.ok -side top -fill x
}

proc tout_print2 {w t tag output} {
    global $t.Header $t.Params $t.Text tk_utils_defs

    tout_output_get_checks $w $t $tag $output
 
    if {[catch {set f [open "|[keylget tk_utils_defs PRINT_COMMAND]" w]}]} {
	destroy $t
	bell
	return
    }

    puts $f "[set $t.Header][set $t.Params][set $t.Text]"
    close $f

    destroy $t
}

#
# Outputs a section to disk
#
proc tout_output_disk {w tag} {
    global tk_utils_defs
    set t [keylget tk_utils_defs TEXT_OUTPUT.WIN]
    xtoplevel $t -resizable 0
    wm title $t "Output to disk"

    getFname $t.list [keylget tk_utils_defs TEXT_OUTPUT.FILENAME] save
    tout_output_checks $t.output

    okcancelhelp $t.ok \
        -ok_command "tout_output_disk2 $t $w $tag $t.output $t.list.entry" \
        -cancel_command "destroy $t" \
	-help_command "show_help interface {UI-Output}"\
        -bd 2 -relief groove

    pack $t.list $t.output $t.ok -side top -fill x
}


proc tout_output_disk2 {t w tag output entryw} {
    global $t.Header $t.Params $t.Text

    if {[set name [entrybox_get $entryw]] == ""} {
	return
    }

    tout_output_get_checks $w $t $tag $output
    if {![string match file* [set f [open $name w]]]} {
	bell
	return
    }

    destroy $t
    puts $f "[set $t.Header][set $t.Params][set $t.Text]"
    close $f
    unset $t.Header $t.Params $t.Text
}

#
# Outputs a section to list
#
proc tout_output_list {w tag} {
    global tk_utils_defs 
    set t [keylget tk_utils_defs TEXT_OUTPUT.WIN]
    xtoplevel $t -resizable 0
    wm title $t "Output to list"

    getLname $t.list [keylget tk_utils_defs TEXT_OUTPUT.LISTNAME] create
    tout_output_checks $t.output

    okcancelhelp $t.ok \
        -ok_command "tout_output_list2 $w $t $tag $t.list.entry	$t.output" \
        -cancel_command "destroy $t" \
	-help_command "show_help interface {UI-Output}"\
        -bd 2 -relief groove

    pack $t.list $t.output $t.ok -side top -fill x
}

proc tout_output_list2 {w t tag namew output} {
    global $t.Header $t.Params $t.Text

    if {[set name [entrybox_get $namew]] == ""} {
	tk_messageBox -icon error -type ok -title "No list" \
	    -message "No list name has been entered"
	return
    }

    tout_output_get_checks $w $t $tag $output
 
    if [ListCreate $name [StringToList "[set $t.Header][set $t.Params][set $t.Text]"]] {
	destroy $t
	unset $t.Header $t.Params $t.Text
	ListEdit $name
    }
}

proc tout_search_OK_Pressed { text entry direction} {

    set dir [radiolist_get $direction]
    set string [entrybox_get $entry]
    textSearchPos $text $string search $dir
}

#
# Searches for a string text display
#
proc tout_search {text} {

    set s $text.search
    if {[xtoplevel $s -resizable 0] == ""} {
	return
    }
    wm title $s "Search"
    entrybox $s.entry \
	    -title "Search string"\
	    -width 20 \
	    -type CheckString

    set b1 "forward"
    set b2 "backward"

    radiolist $s.direction \
	    -orient horizontal \
	    -default 1\
	    -buttons [format { \
	    { %s } \
	    { %s } } \
	    [list $b1] \
	    [list $b2] ]
	    

    okcancelhelp $s.ok_cancel \
	    -ok_command "tout_search_OK_Pressed $text $s.entry $s.direction"\
	    -cancel_command "destroy $s"\
	    -help_command "show_help interface {UI-Output}"\
	    -bd 2 \
	    -relief groove

    pack $s.entry -fill both -expand 1
    pack $s.direction -fill both -expand 1
    pack $s.ok_cancel -fill both -expand 1

    $text tag configure search -background #ffff00 
}
#
# Clears the text display
#
proc tout_clear {text} {
    $text delete 1.0 end
}


#
# Sets the "scroll to bottom upon text output" mode for this 
#
proc tout_scroll {win_name value} {
    tout_set_scroll $win_name $value
}


#
# Redirect all output to a specific file
#
proc tout_redir_open {label_var stream menu} {
    global tk_utils_defs

    set t [keylget tk_utils_defs TEXT_REDIR.WIN]

    xtoplevel $t -resizable 0
    wm title $t "Open redirection file"

    getFname $t.file [keylget tk_utils_defs TEXT_REDIR.FILENAME] save
    
    okcancelhelp $t.ok \
        -ok_command "tout_redir_open2 $t \[entrybox_get $t.file.entry\] \
		     $label_var $stream $menu" \
        -cancel_command "destroy $t" \
	-help_command "show_help interface {UI-Output}"\
        -bd 2 -relief groove

    pack $t.file $t.ok -side top -fill x
}

proc tout_redir_open2 {t name label_var stream menu} {
    global $label_var

    if {[tout_set_redir $stream $name]} {
	set $label_var [file tail $name]
        $menu entryconfigure 1 -state normal
    } else {
	bell
    }

    destroy $t
}

proc tout_redir_close {label_var stream menu} {
    global $label_var

    if {[tout_set_redir $stream ""]} {
	set $label_var ""
	$menu entryconfigure 1 -state disabled 
    } else {
	bell
    }
}


#
# Outputs to a command
#
proc tout_output_command {w tag} {
    global tk_utils_defs 

    set t [keylget tk_utils_defs TEXT_OUTPUT.WIN]
    xtoplevel $t -resizable 0
    wm title $t "Output to command"

    entrybox $t.command -title "Command"
    tout_output_checks $t.output
    
    okcancelhelp $t.ok \
        -ok_command "tout_output_command2 $t $w $tag $t.output \[entrybox_get $t.command\]" \
        -cancel_command "destroy $t" \
	-help_command "show_help interface {UI-Output}"\
        -bd 2 -relief groove

    pack $t.command $t.output $t.ok -side top -fill x    
}

proc tout_output_command2 {t w tag output comm} {
    global $t.Header $t.Params $t.Text

    destroy $t
    tout_output_get_checks $w $t $tag $output

    tout_pipe $comm "[set $t.Header][set $t.Params][set $t.Text]" 0
    unset $t.Header $t.Params $t.Text
}

# Configures which items appear in text-output popup menu.
#
# Commands is a Tcl list matching the command names listed in tout_popup.
proc tout_config_popup {commands} {
    global .tout_popup
    catch {unset .tout_popup}
    foreach com $commands {
	set .tout_popup($com) 1
    }
}
