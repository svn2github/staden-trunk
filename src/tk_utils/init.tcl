#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

lappend auto_path $env(STADLIB)

if {[info command button] != ""} {
    catch {::tk::Darken}
    # FIXME!
    # On Aqua-Tk the minimum valid size for button's -padx is
    # 11. Trying lower than that yields clipped text.
    if {$tcl_platform(os) == "Darwin"} {
	rename button _button
	proc button {args} {
	    set ind [lsearch $args "-padx"]
	    if {$ind != -1} {
		incr ind
		set padx [lindex $args $ind]
		if {$padx < 11} {set padx 11}
	    } else {
		set padx 11
	    }
	    return [uplevel _button $args -padx $padx]
	}
    }

#-----------------------------------------------------------------------------
# Virtual bindings
#-----------------------------------------------------------------------------
if {$tcl_platform(os) == "Darwin"} {
    # Selecting items
    event add <<select>>		<1> <Command-Key-1>
    event add <<select-menu>>		<1> <Command-Key-1>
    event add <<not-select>>		<2> <3> <Command-Key-2> <Command-Key-3>
    event add <<use>>			<Double-1> <Double-Command-Key-1>
    event add <<select-drag>>		<B1-Motion>
    event add <<select-release>>	<B1-ButtonRelease>
    event add <<select-to>>		<Shift-1> <Shift-Command-1>
    
    # Moving items
    event add <<move>>	\
	<2> <Command-Key-2> <Alt-Button-1> <Option-Button-1>
    event add <<move-create>> \
	<Double-2> <Double-Command-Key-2> <Alt-Double-Button-1> <Option-Double-Button-1>
    event add <<move-drag>>		<B2-Motion> <Any-B1-Motion>
    event add <<move-release>> \
	<ButtonRelease-2> <Any-B1-ButtonRelease>
    event add <<move-autoscroll>>	<B2-Leave> <Any-B1-Leave>
    event add <<stop-autoscroll>>	<Any-Enter>
    
    # Zooming
    event add <<zoom>> \
	<Control-Button-3> <Control-Command-Button-1> <Control-Command-Key-3>
    event add <<zoom-drag>>		<Any-B3-Motion>
    event add <<zoom-release>>	\
	<Control-B3-ButtonRelease> <B3-ButtonRelease>
    
    # Popup menus
    event add <<menu>>			<3> <Command-Key-3> <Command-Button-1>
    
    # Editor bits
    event add <<trace>>	\
	<Double-2> <Double-1> <Double-Command-Key-1>
    event add <<save>>	\
	<Control-Key-x>s <Control-Key-x><Control-Key-s>
    event add <<toggle-read-only>>	<Control-Key-x><Control-Key-q>
    event add <<search>>		<Control-Key-s>
    event add <<rsearch>>		<Control-Key-r> <Escape><Control-Key-s>
    
    # Cut and Pasting
    event add <<copy>>	\
	<Control-Key-c> <Control-Insert> <Command-Key-c>
    event add <<paste>>	\
	<Control-Key-v> <Shift-Insert> <Command-Key-v>
    event add <<cut>>			<Control-Key-x> <Command-Key-x>
} else {
    # Selecting items
    event add <<select>>		<1>
    event add <<select-menu>>		<1>
    event add <<not-select>>		<2> <3>
    event add <<use>>			<Double-1>
    event add <<select-drag>>		<B1-Motion>
    event add <<select-release>>	<B1-ButtonRelease>
    event add <<select-to>>		<Shift-1>
    
    # Moving items
    event add <<move>>			<2> <Alt-Button-1>
    event add <<move-create>>		<Double-2>
    event add <<move-drag>>		<B2-Motion> <Any-B1-Motion>
    event add <<move-release>>		<ButtonRelease-2> <Any-B1-ButtonRelease>
    event add <<move-autoscroll>>	<B2-Leave> <Any-B1-Leave>
    event add <<stop-autoscroll>>	<Any-Enter>
    
    # Zooming
    event add <<zoom>>			<Control-Button-3>
    event add <<zoom-drag>>		<Any-B3-Motion>
    event add <<zoom-release>>		<Control-B3-ButtonRelease> <B3-ButtonRelease>
    
    # Popup menus
    event add <<menu>>			<3>
    
    # Editor bits
    event add <<trace>>			<Double-2> <Double-1>
    event add <<save>>			<Control-Key-x>s <Control-Key-x><Control-Key-s>
    event add <<toggle-read-only>>	<Control-Key-x><Control-Key-q>
    event add <<search>>		<Control-Key-s>
    event add <<rsearch>>		<Control-Key-r> <Escape><Control-Key-s>
    
    # Cut and Pasting
    event add <<copy>>			<Control-Key-c> <Control-Insert>
    event add <<paste>>			<Control-Key-v> <Shift-Insert>
    event add <<cut>>			<Control-Key-x>
}

#-----------------------------------------------------------------------------
# Font specifications
#-----------------------------------------------------------------------------
global tk_utils_defs
eval font create trace_font  		[keylget tk_utils_defs FONT.TRACE]
eval font create trace_conf_font  	[keylget tk_utils_defs FONT.TRACE_CONF]
eval font create title_font  		[keylget tk_utils_defs FONT.TITLE]
eval font create text_font  		[keylget tk_utils_defs FONT.TEXT]
eval font create sheet_font  		[keylget tk_utils_defs FONT.SHEET]
eval font create sheet_bold_font 	[keylget tk_utils_defs FONT.SHEET] -weight bold
eval font create menu_title_font  	[keylget tk_utils_defs FONT.MENU_TITLE]
eval font create menu_font 		[keylget tk_utils_defs FONT.MENU]
eval font create listbox_font  		[keylget tk_utils_defs FONT.LISTBOX]
eval font create button_font  		[keylget tk_utils_defs FONT.BUTTON]
eval font create enzyme_font		[keylget tk_utils_defs FONT.ENZYME]

# Sometimes asking for a bold font will give us a different shaped font as
# the exact match wasn't available. Double check this and revert to normal
# if necessary.
if {[font metrics sheet_font] != [font metrics sheet_bold_font] || \
    [font measure sheet_font lllM] != [font measure sheet_bold_font lllM]} {
    puts stderr "Warning: Normal and bold sheet fonts differ in size."
    font delete sheet_bold_font
    eval font create sheet_bold_font 	[keylget tk_utils_defs FONT.SHEET]
}

option add *Menu.font		menu_font
option add *Text.font		text_font
option add *Entry.font		text_font
option add *Button.font		button_font
option add *Menubutton.font	button_font
option add *Checkbutton.font	button_font
option add *Radiobutton.font	button_font
option add *Label.font		button_font
option add *labelFont		button_font
option add *Listbox.font	listbox_font

option add *trace*colour1	[keylget tk_utils_defs TRACE.COLOUR_A]
option add *trace*colour2 	[keylget tk_utils_defs TRACE.COLOUR_C]
option add *trace*colour3 	[keylget tk_utils_defs TRACE.COLOUR_G]
option add *trace*colour4 	[keylget tk_utils_defs TRACE.COLOUR_T]
option add *trace*cursorColour	[keylget tk_utils_defs TRACE.COLOUR_CURSOR]
option add *trace*qualColour	[keylget tk_utils_defs TRACE.COLOUR_QUALITY]
option add *trace*vectorColour	[keylget tk_utils_defs TRACE.COLOUR_VECTOR]
option add *trace*vectorColour	[keylget tk_utils_defs TRACE.COLOUR_VECTOR]
option add *trace*lineWidth	[keylget tk_utils_defs TRACE.LINE_WIDTH]

foreach binding [bind Text] {
    bind TextRO $binding [bind Text $binding]
}
bind TextRO <Tab> {}
bind TextRO <Control-i> {}
bind TextRO <Return> {}
bind TextRO <Insert> {}
bind TextRO <KeyPress> {}
bind TextRO <<Cut>> {}
bind TextRO <<Paste>> {}
bind TextRO <<PasteSelection>> {}
bind TextRO <Control-d> {}
bind TextRO <Control-k> {}
bind TextRO <Control-o> {}
bind TextRO <Meta-d> {}
bind TextRO <Meta-BackSpace> {}
bind TextRO <Meta-Delete> {}
bind TextRO <Delete> {}
bind TextRO <BackSpace> {}
bind TextRO <Control-h> {}
bind TextRO <Control-t> {}
}

proc fix_raise {} {
    bind Toplevel <Visibility> {set visibilityState(%W) %s}
,    bind Toplevel <Destroy> {catch {unset visibilityState(%W)}}
  
    rename raise raise_orig
  
    proc raise {toplevel {aboveThis ""}} {
	if {![info exists ::visibilityState($toplevel)] ||
	    ![string equal $::visibilityState($toplevel) "VisibilityUnobscured"]} {
	    eval raise_orig $toplevel $aboveThis
	}
    }
}

proc tk_utils_init {} {
    global tk_utils_defs tcl_platform
    global tcl_pkgPath
    global env

    # Add STADLIB to package path
    lappend tcl_pkgPath $env(STADLIB)

    # Abort now if we're being used in a non-graphical environment.
    if {![namespace exists ::tk]} {
	return
    }

    global ::tk::MotifFilebrowser
    set ::tk::MotifFilebrowser 1

    entry .e
    set bgg [.e cget -background]
    set fgg [.e cget -foreground]
    button .tmp
    set bg [.tmp cget -background]
    set fg [.tmp cget -foreground]
    if {$tcl_platform(platform) == "unix"} {
	#set bg [keylget tk_utils_defs BACKGROUND] 
	#if {$bg == ""} {set bg "#d9d9d9"} 
	option add *background $bg userDefault 

	#button .tmp
	#option add *disabledForeground [.tmp cget -disabledforeground] widgetDefault
	#set bg [.tmp cget -background]
	#set fg [.tmp cget -foreground]
	option add *background $bg userDefault
	option add *highlightBackground $bg userDefault
	option add *tabBackground $bg userDefault
	option add *foreground $fg userDefault
	option add *tabForeground $fg userDefault
	option add *selectBackground lightblue userDefault
	option add *selectBorderWidth 0 userDefault
	set bgl [tk::Darken $bg 115]
	option add *backdrop $bg userDefault
	option add *Entry.background $bgl userDefault
	option add *Text.background $bgl userDefault
	option add *Listbox.background $bgl userDefault
	option add *Spinbox.background $bgl userDefault
	# option add *Xentry.background $bgl userDefault
	option add *TixHList.background $bgl userDefault
	option add *Tablelist.background $bgl userDefault
	#option add *Comborange*entry.disabledForeground $fg userDefault
	#option add *Comborange*entry.disabledBackground $bgl userDefault
	option add *textBackground $bgl widgetDefault 
	destroy .tmp
	if {[keylget tk_utils_defs FOREGROUND] != {} || \
	    [keylget tk_utils_defs BACKGROUND] != {}} {
	        RecolourWidgets [keylget tk_utils_defs FOREGROUND] \
		    [keylget tk_utils_defs BACKGROUND]
	}

	#Fix raise for fvwm.
	#fix_raise
    } else {
	# Fix bug in Windows event propagation mechanisms
	# MJ: See CHANGES\TK803\CHANGE1.TXT - Fix no longer required here.
	# global env
	# source $env(STADENROOT)/lib/tk_utils/winfix.tcl
	option add *Tablelist.background SystemWindow userDefault
    }
    option add *Combobox*entry.disabledForeground $fgg userDefault
    option add *Combobox*entry.disabledBackground $bgg userDefault
    option add *Tabnotebook.background $bg userDefault
    option add *Tabnotebook.backdrop $bg userDefault
    option add *Tabnotebook.tabBackground $bg userDefault
    
    catch {
	if {[wm protocol . WM_DELETE_WINDOW] == ""} {
	    wm protocol . WM_DELETE_WINDOW exit
	}
    }
    destroy .e
    destroy .tmp
}

#-----------------------------------------------------------------------------
# Map new 8.5 widgets to 8.4 coded alternatives if not available
#-----------------------------------------------------------------------------
if {[info command ttk::notebook] == {}} {
    # Emulate ttk::notebook via Jeffrey Hobbs' widget package
    namespace eval ttk {
	proc notebook {args} {
	    set nb [eval tabnotebook $args -linewidth 1]
	    if {[info command ::_$nb] != {}} {
		rename ::_$nb {}
	    }
	    rename ::$nb ::_$nb
	    global ::$nb
	    set ::${nb}(count) 0

	    proc ::$nb {cmd args} {
		set nb [lindex [info level [info level]] 0]
		upvar \#0 ::$nb data
		switch $cmd {
		    "add" {
			foreach {k v} [lrange $args 1 end] {
			    set a($k) $v
			}
			set data(map_$data(count)) $a(-text)
			incr data(count)
			return [eval [list _$nb add $a(-text) \
					  -window [lindex $args 0]]]
		    }
		    "select" {
			if {[llength $args] == 0} {
			    return [eval [list _$nb select]]
			} else {
			    return [eval [list _$nb activate \
					      $data(map_[lindex $args 0])]]
			}
		    }
		    default {
			return [eval [list _$nb $cmd $args]]
		    }
		}
	    }

	    return $nb
	}
    }
}
