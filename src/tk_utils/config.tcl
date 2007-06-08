#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#-----------------------------------------------------------------------------
# Changing the fonts
proc update_font_metrics {} {
    global env tk_utils_defs _font_metrics
    
    array set _font_metrics [keylget tk_utils_defs FONT_METRICS]

    set new_fonts 0
    foreach font [lsort [font families]] {
	if {![info exists _font_metrics($font)]} {
	    array set d [font metrics [list $font]]
	    set _font_metrics($font) $d(-fixed)
	    incr new_fonts
	}
    }

    if {$new_fonts} {
	keylset tk_utils_defs FONT_METRICS [array get _font_metrics]
	update_defs tk_utils_defs $env(HOME)/.tk_utilsrc FONT_METRICS
    }
}

proc font_name_change {w n} {
    $w.l.l1 configure -font $n
    $w.l.l2 configure -font $n
    $w.l.l3 configure -font $n
    array set fn [font configure $n]

    do {
	set ind [lsearch -exact [string tolower [$w.f2.family get 0 end]] \
		     [string tolower $fn(-family)]]
	if {$ind == -1} {
	    $w.f2.family insert end $fn(-family)
	}
    } while {$ind == -1}

    $w.f2.family activate $ind
    $w.f2.family selection clear 0 end
    $w.f2.family selection set $ind
    $w.f2.family see $ind

    global $w.f3.Bold $w.f3.Italic $w.f3.Fixed $w.f3.Underline
    set $w.f3.Bold [font configure $n -weight]
    set $w.f3.Italic [font configure $n -slant]
    set $w.f3.Underline [font configure $n -underline]
    font_size_change $w [expr {abs([font configure $n -size])}]
}

proc font_family_change {w f} {
    do {
	set ind [lsearch -exact [string tolower [$w.f2.family get 0 end]] \
		     [string tolower $f]]
	if {$ind == -1} {
	    $w.f2.family insert end $f
	}
    } while {$ind == -1}

    $w.f2.family activate $ind
    $w.f2.family selection clear 0 end
    $w.f2.family selection set $ind
    $w.f2.family see $ind

    set fn [$w.f2.name get [$w.f2.name curselection]]
    font configure $fn -family $f
}

proc font_size_change {w s} {
    do {
	set ind [lsearch -exact -integer [$w.f2.size get 0 end] $s]
	if {$ind == -1} {
	    $w.f2.size insert end $s
	}
    } while {$ind == -1}

    $w.f2.size activate $ind
    $w.f2.size selection clear 0 end
    $w.f2.size selection set $ind
    $w.f2.size see $ind

    set fn [$w.f2.name get [$w.f2.name curselection]]
    font configure $fn -size -[expr {abs($s)}]
}

proc font_style_change {w name var} {
    set fn [$w.f2.name get [$w.f2.name curselection]]

    global $var
    eval font configure $fn $name [set $var]
}

proc font_filter {w} {
    global _font_metrics $w.f3.Fixed
    update_font_metrics

    $w.f2.family delete 0 end
    foreach i [lsort [array names _font_metrics]] {
	if {[set $w.f3.Fixed] && $_font_metrics($i) == 0} {
	    continue
	}
	$w.f2.family insert end $i
    }
}

proc SetFonts {{w .font_config}} {
    if {$w != ""} {
	if {[xtoplevel $w -resizable 0] == ""} return
	wm title $w "Set fonts"
    }

    set old ""
    foreach i [font names] {
	lappend old "$i [font config $i]"
    }

    pack [frame $w.f -bd 2 -relief raised]
    okcancelhelp $w.ok -bd 2 -relief groove \
	-ok_command "destroy $w" \
	-perm_command "destroy $w; SetFontsPerm" \
	-cancel_command "destroy $w;SetFontsReset [list $old]" \
	-help_command "show_help interface {UI-Fonts}"
    pack $w.ok -side bottom -fill both
    set w $w.f

    frame $w.l -height 100 -width 400 -bd 1 -relief groove
    pack propagate $w.l 0
    label $w.l.l1 -text "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    label $w.l.l2 -text "abcdefghijklmnopqrstuvwxyz"
    label $w.l.l3 -text "01234567890!@#$%^&*()_-=_+"
    pack $w.l.l1 $w.l.l2 $w.l.l3
    set m [frame $w.f2]

    frame $w.def
    button $w.def.small -text "Default: small" -command "font_small $w"
    button $w.def.large -text "Default: large" -command "font_large $w"
    pack $w.def.small $w.def.large -side left -expand 1

    listbox $m.name -yscrollcommand "$m.namey set" -exportselection 0
    scrollbar $m.namey -orient vert -command "$m.name yview"
    bind $m.name <<select>> "+%W activate \[%W index @%x,%y\]"
    bind $m.name <<ListboxSelect>> \
	"+font_name_change $w \[%W get active\];"
    bind $m.name <Any-Enter> {focus %W}
    foreach i [font names] {
	$m.name insert end $i
    }

    listbox $m.family -yscrollcommand "$m.familyy set" -exportselection 0 \
	-selectmode browse
    scrollbar $m.familyy -orient vert -command "$m.family yview"
    bind $m.family <<select>> "+%W activate \[%W index @%x,%y\]"
    bind $m.family <<ListboxSelect>> \
	"+font_family_change $w \[%W get active\];"
    bind $m.family <Any-Enter> {focus %W}

    listbox $m.size -yscrollcommand "$m.sizey set" -exportselection 0 -width 10
    scrollbar $m.sizey -orient vert -command "$m.size yview"
    bind $m.size <<select>> "+%W activate \[%W index @%x,%y\]"
    bind $m.size <<ListboxSelect>> \
	"+font_size_change $w \[%W get active\];"
    bind $m.size <Any-Enter> {focus %W}
    for {set i 5} {$i <= 24} {incr i} {
	$m.size insert end "$i"
    }

    frame $w.f3
    checkbutton $w.f3.bold -text "Bold" \
	-command [list font_style_change $w -weight $w.f3.Bold] \
	-variable $w.f3.Bold \
	-onvalue bold -offvalue normal
    checkbutton $w.f3.italic -text "Italic" \
	-command [list font_style_change $w -slant $w.f3.Italic] \
	-variable $w.f3.Italic \
	-onvalue italic -offvalue roman
    checkbutton $w.f3.fixed -text "Fixed" \
	-command [list font_filter $w] \
	-variable $w.f3.Fixed
    checkbutton $w.f3.underline -text "Underline" \
	-command [list font_style_change $w -underline $w.f3.Underline] \
	-variable $w.f3.Underline

    pack $m.name $m.namey $m.family $m.familyy $m.size $m.sizey \
	-side left -fill both
    pack $w.f3.bold $w.f3.italic $w.f3.fixed $w.f3.underline -side left
    pack $w.def $w.l $w.f2 $w.f3  -side bottom -fill both
    
    $w.f2.name activate 0
    $w.f2.name selection set 0
    font_filter $w
    font_name_change $w [lindex [font names] 0]
}

proc font_small {w} {
    font configure menu_font         -family Helvetica -weight bold -size -12
    font configure menu_title_font   -family Helvetica -weight bold -slant italic -size -12
    font configure text_font         -family Courier -size -12
    font configure button_font       -family Helvetica -weight bold -size -12
    font configure listbox_font      -family Helvetica -weight bold -size -12
    font configure title_font        -family Helvetica -size -16 -weight bold
    font configure sheet_font        -family Fixed -size -20
    font configure trace_font        -family Helvetica -size -12
    font configure enzyme_font       -family Helvetica -size -11
    font configure trace_conf_font   -family Helvetica -size -9

    font_name_change $w [lindex [font names] 0]
}

proc font_large {w} {
    font configure menu_font         -family Helvetica -weight bold -size -14
    font configure menu_title_font   -family Helvetica -weight bold -slant italic -size -16
    font configure text_font         -family Courier -size -16
    font configure button_font       -family Helvetica -weight bold -size -15
    font configure listbox_font      -family Helvetica -weight bold -size -16
    font configure title_font        -family Helvetica -size -16 -weight bold
    font configure sheet_font        -family Fixed -size -20
    font configure trace_font        -family Helvetica -size -12
    font configure enzyme_font       -family Helvetica -size -11
    font configure trace_conf_font   -family Helvetica -size -9

    font_name_change $w [lindex [font names] 0]
}

proc SetFontsReset {old} {
    foreach i $old {
	eval font configure $i
    }
}

proc SetFontsPerm {} {
    global env tk_utils_defs
    foreach i [font names] {
	regsub {_font$} $i {} fn
        set fn FONT.[string toupper $fn]
	keylset tk_utils_defs $fn [font configure $i]
        update_defs tk_utils_defs $env(HOME)/.tk_utilsrc $fn
    }
}


#-----------------------------------------------------------------------------
# Sets the general colour scheme - foreground and background.
proc ConfigureColours {} {
    set t .colour
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Set colours"

    label $t.eg -text "Example text" -bd 2 -relief raised -height 2
    frame $t.fg -bd 2 -relief groove
    frame $t.bg -bd 2 -relief groove
    label $t.fg.label -text Foreground
    label $t.bg.label -text Background
    frame $t.fg.col
    frame $t.bg.col

    ColourSel $t.fg.col foreground "ColourSet $t.eg foreground"
    ColourSel $t.bg.col background "ColourSet $t.eg background"

    okcancelhelp $t.ok -bd 2 -relief groove \
	-ok_command "ConfigureColours2 0 $t" \
	-perm_command "ConfigureColours2 1 $t" \
	-cancel_command "destroy $t" \
	-help_command "show_help interface {UI-Colour}"

    pack $t.ok -side bottom -fill both
    pack $t.eg -side bottom -fill both -pady 2m -padx 2m
    pack $t.fg.label $t.bg.label -side top -fill both
    pack $t.fg.col $t.bg.col -side top -fill both -expand 1
    pack $t.fg $t.bg -side left -expand 1 -fill both
}

proc RecolourWidgets {fg bg} {

  #Dont support colour changing under windows
  #othwerwise things get just too complicated
  #users can still change colours using the standard windows properties
  global tcl_platform
  if { "$tcl_platform(platform)" != "windows" } {

    foreach {fga(r) fga(g) fga(b)} [winfo rgb . $fg] {}
    foreach {bga(r) bga(g) bga(b)} [winfo rgb . $bg] {}

    set bg_dark [format #%04x%04x%04x \
	[expr int(.9*$bga(r))] [expr int(.9*$bga(g))] [expr int(.9*$bga(b))]]

    set fg_dis [format #%04x%04x%04x \
	    [expr {int(.75*$bga(r)+.25*$fga(r))}] \
	    [expr {int(.75*$bga(g)+.25*$fga(g))}] \
	    [expr {int(.75*$bga(b)+.25*$fga(b))}]]

    foreach i {r g b} {
	set inc1 [expr ($bga($i)*15)/100]
	set inc2 [expr (65535-$bga($i))/3]
	set bgl($i) [expr $bga($i) + ($inc1>$inc2?$inc1:$inc2)]
	if {$bgl($i) > 65535} { set bgl($i) 65535 }
    }
    set bg_light [format #%04x%04x%04x $bgl(r) $bgl(g) $bgl(b)]

    tk_setPalette \
	background $bg \
	backPageColor $bg \
	highlightBackground $bg \
	sliderForeground $bg \
	inactiveBackground $bg_dark \
	selectBackground $bg_dark \
	troughColor $bg_dark \
	activeBackground $bg_light \
	sliderBackground $bg_light \
	foreground $fg \
	activeForeground $fg \
	insertBackground $fg \
	selectForeground $fg \
	highlightColor $fg \
	disabledForeground $fg_dis

    catch {destroy .tmp}
    button .tmp
    option add *Entry.background [tk::Darken [.tmp cget -background] 115] 100
    option add *Text.background [tk::Darken [.tmp cget -background] 115] 100
    option add *Listbox.background [tk::Darken [.tmp cget -background] 115] 100
    option add *EdNames.background [tk::Darken [.tmp cget -background] 115] 100
    option add *Editor.background [tk::Darken [.tmp cget -background] 115] 100
    destroy .tmp

  }
}

proc ConfigureColours2 {perm t} {
    global tk_utils_defs env

    set fg [format #%02x%02x%02x \
	[$t.fg.col.r get] [$t.fg.col.g get] [$t.fg.col.b get]]
    set bg [format #%02x%02x%02x \
	[$t.bg.col.r get] [$t.bg.col.g get] [$t.bg.col.b get]]

    RecolourWidgets $fg $bg

    keylset tk_utils_defs FOREGROUND $fg
    keylset tk_utils_defs BACKGROUND $bg

    if {$perm} {
	update_defs tk_utils_defs $env(HOME)/.tk_utilsrc FOREGROUND BACKGROUND
    }

    destroy $t
}

proc ColourSet {eg part w args} {
    $eg configure -$part [format #%02x%02x%02x \
	[$w.r get] [$w.g get] [$w.b get]]
}

proc ColourSel {w part cmd} {
    scale $w.r -label red -from 0 -to 255 -orient horiz -length 5c \
	-command "$cmd $w"
    scale $w.g -label green -from 0 -to 255 -orient horiz -length 5c \
	-command "$cmd $w"
    scale $w.b -label blue -from 0 -to 255 -orient horiz -length 5c \
	-command "$cmd $w"
    set rgb [winfo rgb . [$w.r cget -$part]]
    $w.r set [expr [lindex $rgb 0]/256]
    $w.g set [expr [lindex $rgb 1]/256]
    $w.b set [expr [lindex $rgb 2]/256]
    pack $w.r $w.g $w.b -side top -fill both
}

#-----------------------------------------------------------------------------
# Handles storing data in the .<prog>rc files.

# Processes a *rc config file to look for the automatic positioning lines.
# Anything within this region is considered fair game for auto-
# configuration.
#
# The data returned is the contents of the config file in a structure used only
# by update_config and write_config.
proc read_config {file} {
    if {[catch {set fd [open $file r]}]} {
	return {{} {} {}}
    }

    set part1 {}
    set part2 {}
    set part3 {}
    set part part1
    while {[gets $fd line] != -1} {
	if {$line == "# DO NOT EDIT: start of auto-config"} {
	    set part part2
	    continue
        }
	if {$line == "# DO NOT EDIT: end of auto-config"} {
	    set part part3
	    continue
        }
	lappend $part $line
    }
    close $fd

    return [list $part1 $part2 $part3]
}

proc display_config {cf} {
    puts part1:[lindex $cf 0]
    puts part2:[lindex $cf 1]
    puts part3:[lindex $cf 2]
}

proc write_config {cf file} {
    catch {file rename -force $file $file.bak}
    if {[catch {set fd [open $file w 0644]}]} {
	puts "Couldn't save config to $file"
	return
    }

    foreach i [lindex $cf 0] {
	puts $fd $i
    }
    puts $fd "# DO NOT EDIT: start of auto-config"
    foreach i [lindex $cf 1] {
	puts $fd $i
    }
    puts $fd "# DO NOT EDIT: end of auto-config"
    foreach i [lindex $cf 2] {
	puts $fd $i
    }

    close $fd
    return
}

proc update_config {cfp name value} {
    upvar $cfp cf
    set p1 [lindex $cf 0]
    set p2 [lindex $cf 1]
    set p3 [lindex $cf 2]
    set p2new {}
    set done 0
    
    foreach i $p2 {
	if {[string match "set_def $name\[ \t\]*" $i] == 1} {
	    lappend p2new "set_def $name [list $value]"
	    set done 1
	} else {
	    lappend p2new $i
	}
    }
    if {!$done} {
        lappend p2new "set_def $name [list $value]"
    }
    set cf [list $p1 $p2new $p3]
}

proc update_config_def {cfp def name} {
    global $def
    update_config $cfg $name [keylget $def $name]
}

# Updates a series of defaults
proc update_defs {def_name file args} {
    global $def_name
    set cf [read_config $file]
    foreach i $args {
	update_config cf $i [keylget $def_name $i]
    }
    write_config $cf $file
}
