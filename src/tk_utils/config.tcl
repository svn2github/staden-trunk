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
proc font_name_change {w n} {
    $w.l1 configure -font $n
    $w.l2 configure -font $n
    $w.l3 configure -font $n
    $w.f2.name configure -text $n
    $w.f2.family configure -text [font configure $n -family]
    global $w.f3.Bold $w.f3.Italic $w.f3.Overstrike $w.f3.Underline
    set $w.f3.Bold [font configure $n -weight]
    set $w.f3.Italic [font configure $n -slant]
    set $w.f3.Overstrike [font configure $n -overstrike]
    set $w.f3.Underline [font configure $n -underline]
    font_size_change $w [font configure $n -size]
}

proc font_family_change {w f} {
    set fn [$w.f2.name cget -text]
    $w.f2.family configure -text $f
    font configure $fn -family $f
}

proc font_size_change {w s} {
    set fn [$w.f2.name cget -text]
    if {$s > 0} {
	$w.f2.size configure -text "$s points"
    } else {
	$w.f2.size configure -text "[expr -$s] pixels"
    }
    font configure $fn -size $s
}

proc font_style_change {w name var} {
    set fn [$w.f2.name cget -text]

    global $var
    eval font configure $fn $name [set $var]
}

proc SetFonts {{w .font_config}} {
    if {$w != ""} {
	if {[xtoplevel $w -resizable 0] == ""} return
	wm title $w "Set fonts"
    }

    set font_names [lsort [font families]]

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

    label $w.l1 -text "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    label $w.l2 -text "abcdefghijklmnopqrstuvwxyz"
    label $w.l3 -text "01234567890!@#$%^&*()_-=_+"
    set m [frame $w.f2]

    xmenubutton $m.name -indicatoron 1 -text name -menu $m.name.m
    xmenubutton $m.family -indicatoron 1 -text family -menu $m.family.m
    xmenubutton $m.size -indicatoron 1 -text size -menu $m.size.m

    menu $m.name.m -tearoff 0
    foreach i [font names] {
	$m.name.m add command -label $i \
		-command [list font_name_change $w $i]
    }

    menu $m.family.m -tearoff 0
    foreach i $font_names {
	$m.family.m add command -label $i \
		-command [list font_family_change $w $i]
    }

    menu $m.size.m -tearoff 0
    for {set i 5} {$i < 16} {incr i} {
	$m.size.m add command -label "$i pixels" \
		-command [list font_size_change $w -$i]
    }
    $m.size.m add separator
    for {set i 8} {$i < 20} {incr i} {
	$m.size.m add command -label "$i points" \
		-command [list font_size_change $w $i]
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
    checkbutton $w.f3.overstrike -text "Overstrike" \
	-command [list font_style_change $w -overstrike $w.f3.Overstrike] \
	-variable $w.f3.Overstrike
    checkbutton $w.f3.underline -text "Underline" \
	-command [list font_style_change $w -underline $w.f3.Underline] \
	-variable $w.f3.Underline

    pack $w.l1 $w.l2 $w.l3 $w.f2 $w.f3 -side top
    pack $m.name $m.family $m.size -side left
    pack $w.f3.bold $w.f3.italic $w.f3.overstrike $w.f3.underline -side left

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
    if {[catch {set fd [open $file w]}]} {
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
