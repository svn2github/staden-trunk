# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

proc usage {} {
    puts stderr "Usage: trev \[options\] file ...\n"
    puts stderr "Valid options are:"
    puts stderr "    -abi/-alf/-ctf         Specify input format (def. any)"
    puts stderr "    -exp/-pln/-ztr/-any    ... \" ..."
    puts stderr "    -write_scf             Allow writing back to SCF, ZTR, CTF"
    puts stderr "    -editscf               Synonym for -write_scf"
    puts stderr "    -xmag scale            X magnification (def. 150)"
    puts stderr "    -ymag scale            Y magnification (def. 10)"
    puts stderr "    -trace_scale           Set max-trace-value (def. 0 => auto)"
    puts stderr "    -fofn fn               Load files from \"fn\""
    puts stderr "    -pregap_mode           Used from old pregap"
    puts stderr "    -restrict              Used from pregap4"
    puts stderr "    --                     Force end of argument list"

    exit 1
}

if {[lindex $argv 0] == "-h" || [lindex $argv 1] == "-help" || [lindex $argv 1]
    == "--help"} {
    usage
}

if {[catch tkinit]} {
    package require Tk
}
catch {console hide}

# Override the current X defaults - these are not obvious, and ORDER COUNTS!
#option add *Background #b03060
#option add *foreground #c3c3c3
#option add *Foreground black
#option add *background #d9d9d9
#option add *troughColor #c3c3c3
#option add *activeForeground black
#option add *activeBackground #ececec
#option add *indicator #b03060
#option add *disabled #a3a3a3
#option add *selectBackground #c3c3c3
#option add *selectForeground black

##############################################################################
proc about {} {
    global VERSION

    set w [modal .about]
    wm title $w "About Trev"
    wm transient $w [winfo toplevel [winfo parent $w]]

    text $w.top -bd 0 -width 40 -height 15 -wrap word
    $w.top tag configure title \
	    -font "-*-Times-Medium-R-Normal--*-180-*-*-*-*-*-*" \
	    -foreground blue \
	    -justify c
    $w.top insert end "Trev version $VERSION\n\n" title
    $w.top insert end "Trev is part of the Staden Package. Copyright (c) Medical Research Council, Laboratory of Molecular Biology, 1995 - 2001. All rights reserved.\n\n"
    $w.top insert end "Written by James Bonfield, Kathryn Beal and Rodger Staden.\n\n"
    $w.top insert end "For further details, see http://www.mrc-lmb.cam.ac.uk/pubseq/"

    button $w.bottom -text "Ok" -command "destroy $w"
    pack $w.top -side top -fill both -expand 1 -padx 5 -pady 5
    pack $w.bottom -side bottom
    wm withdraw $w
    update idletasks
    set x [expr [winfo screenwidth $w]/2 - [winfo reqwidth $w]/2 \
	    - [winfo vrootx [winfo parent $w]]]
    set y [expr [winfo screenheight $w]/2 - [winfo reqheight $w]/2 \
	    - [winfo vrooty [winfo parent $w]]]
    wm geom $w +$x+$y
    wm deiconify $w

    set g [grab current $w]
    set f [focus]
    grab $w
    focus $w
}

##############################################################################
#check to see if the current in_format is acceptable as an out_format ie
#can the file be saved in its current format
proc CheckSaveFormat {informat} {
    global out_types

    foreach i $out_types {
	if {[string compare $informat [lindex $i 0]] == 0} {
	    return 1
	}
    }
    return 0
}

##############################################################################
#set current trace filename
proc SetFilename { file } {
    global filename

    set filename $file
}

##############################################################################
#get current trace filename
proc GetFilename { } {
    global filename

    return $filename
}

##############################################################################
#set status of edit_flag: 1 = trace has been editted; 0 = trace has not
#been editted
proc SetEditFlag { state } {
    global edit_flag

    set edit_flag $state
}

##############################################################################
#get status of edit_flag
proc GetEditFlag { } {
    global edit_flag

    return $edit_flag
}

##############################################################################
proc LoadFile { filename informat pregap_mode} {
    global trev_defs trev_menu tk_utils_defs history

    set trace [keylget trev_defs TRACE.WIN].t
    set load_status [catch {$trace load $filename $informat} out]

    if {$load_status == 1} {
	tk_messageBox -icon error -type ok -title "Load failed" \
		-message "Unable to load $filename with format $informat"
	menu_state_set trev_menu -2 [keylget trev_defs MENU.FRAME]
    } else {
	SetRejectStatus
	SetFilename $filename
	SetEditFlag 0
	set history {}
	menu_state_set trev_menu 2 [keylget trev_defs MENU.FRAME]
	menu_state_set trev_menu -64 [keylget trev_defs MENU.FRAME]
	#enable save menus
	if {[CheckSaveFormat [$trace orig_format]]} {
	    #Can save - eg an experiment file
	    menu_state_set trev_menu 4 [keylget trev_defs MENU.FRAME]
	} else {
	    menu_state_set trev_menu -4 [keylget trev_defs MENU.FRAME]
	}

	if {$pregap_mode} {
	    menu_state_set trev_menu -8 [keylget trev_defs MENU.FRAME]
	} else {
	    menu_state_set trev_menu 8 [keylget trev_defs MENU.FRAME]
	}

	#bit of a hack to deal with setting -edits on the command line
	SetTraceDisplay tdisp

	# Configure the options needed by the print command
	set status [catch {$trace ps_configure}]
	set status [expr $status | [catch {$trace ps_configure_trace \
                -title [file tail $filename]}]]
	set status [expr $status | [catch {$trace ps_configure_trace_line A \
		-line_width [keylget tk_utils_defs PS.A.LW] \
		-colour [keylget tk_utils_defs PS.A.COLOUR]}]]
	set status [expr $status | [catch {$trace ps_configure_trace_line C \
		-line_width [keylget tk_utils_defs PS.C.LW] \
		-colour [keylget tk_utils_defs PS.C.COLOUR]}]]
	set status [expr $status | [catch {$trace ps_configure_trace_line G \
		-line_width [keylget tk_utils_defs PS.G.LW] \
		-colour [keylget tk_utils_defs PS.G.COLOUR]}]]
	set status [expr $status | [catch {$trace ps_configure_trace_line T \
		-line_width [keylget tk_utils_defs PS.T.LW] \
		-colour [keylget tk_utils_defs PS.T.COLOUR]}]]
	set status [expr $status | [catch {$trace ps_configure_trace_line N \
		-colour [keylget tk_utils_defs PS.N.COLOUR]}]]

	if {$status != 0} {
	    tk_messageBox -icon error -type ok -title "Open $filename failed" \
		    -message "Unable to configure print options"
	}
    }
}

##############################################################################
proc fileTypes { block_scf } {
    global in_types out_types print_types

    #   Type names		Extension(s)	Mac File Type(s)
    #
    #---------------------------------------------------------
    set in_types {
	{"ABI"			*.ab1		}
	{"ALF"			*.alf		}
	{"CTF"                  *.ctf           }
	{"EXP"			*.exp		}
	{"SCF"			*.scf		}
	{"ZTR"			*.ztr		}
	{"PLN"			*			}
	{"Any"			*			}
    }
    if {$block_scf} {
	set out_types {
	    {"EXP"		    *.exp		    }
	    {"PLN"		    *			    }
	}
    } else {
	set out_types {
	    {"CTF"                  *.ctf                   }
	    {"EXP"		    *.exp		    }
	    {"SCF"		    *.scf		    }
	    {"ZTR"		    *.ztr                   }
	    {"PLN"		    *			    }
	}
    }
    set print_types {
	{"PS"			*.ps			}
    }
}

##############################################################################
proc fileDialog {operation} {
    global in_types out_types print_types
    global trev_defs

    set file_type ""
    if {$operation == "open"} {
	lappend file_type [tk_getOpenFile -filetypes $in_types \
			       -multiple 65000]
	lappend file_type "any"

    } elseif {$operation == "print"} {
	set file_type [tk_getSaveFile -filetypes $print_types]

    } else {
	set file_type [list [tk_getSaveFile]]
	if {[catch {lappend file_type [tk_getFileType]}]} {
	    # Guess type
	    if {[regexp {\.([^.]*$)} [lindex $file_type 0] dummy type] == 1} {
		lappend file_type $type
	    } else {
		lappend file_type EXP
	    }
	}
    }

    return $file_type
}

##############################################################################
proc biolimsDialog {operation} {
    global trev_defs

    if {$operation == "open"} {
      lappend file_type [spGetOpenBiolimsFile -multiple true]
    } else {
      set file_type [spGetSaveBiolimsFile]
    }

    lappend file_type "Any"
    return $file_type
}


##############################################################################

proc OpenFile { } {

    if {![CheckSaved]} { return }
    OpenFileImplementation [fileDialog open]
}

##############################################################################

proc OpenBiolims { } {

    if {![CheckSaved]} { return }
    OpenFileImplementation [biolimsDialog open]
}

##############################################################################

proc OpenDroppedFiles { filelist } {

    if {![CheckSaved]} { return }
    OpenFileImplementation [list $filelist Any]
}

##############################################################################

proc OpenFileImplementation { file_type } {
    global file_list file_ind file_count file_status

    #check that have 2 values: filename list and format
    assert { [llength $file_type] == 2 }

    set filelist [lindex $file_type 0]
    set informat [lindex $file_type 1]

    if {[string compare $filelist ""]} {
	set file_list $filelist
	set file_ind 0
	set file_count [llength $file_list]
	foreach f $file_list {
	    set file_status($f) 1
	}
	catch {destroy .buttons}

	CreatePrevNext
	NextFile $informat
    }
}


##############################################################################

proc SaveFile {interactive format} {
    global trev_defs

    set trace [keylget trev_defs TRACE.WIN].t
    if {$interactive} {
	#Save As - bring up file browser
	set file_type [fileDialog save]

	set filename [lindex $file_type 0]
	# set format [lindex $file_type 1]
    } else {
	#Save - save in current format
	# set format [$trace format]
	set filename [GetFilename]
    }

    if {[string compare $filename ""] == 0} {
	# Cancel pressed
	return
    }

    set status [catch {$trace save $filename $format}]

    if {$status == 0} {
	SetEditFlag 0
    } else {
	tk_messageBox -icon error -type ok -title "Save failed" \
		-message "Unable to save $filename with format $format"

    }
}

##############################################################################
proc PrintFile {interactive} {
    global trev_defs

    set trace [keylget trev_defs TRACE.WIN].t
    if {$interactive} {
	#Save As - bring up file browser
	set file_type [fileDialog print]

	set filename [lindex $file_type 0]

    } else {
	#Save - save in current format
	set filename [GetFilename]
    }

    if {[string compare "" $filename] == 0} {
	return
    }

    set status [catch {$trace print $filename}]

    if {$status == 0} {
	SetEditFlag 0
    } else {
	tk_messageBox -icon error -type ok -title "Print failed" \
		-message "Unable to print $filename"
    }
}

##############################################################################
proc ps_pg_size_change {w s} {
    global tk_utils_defs $w.PS

    keylset $w.PS PAGE [lindex [keylget tk_utils_defs PS.PAGE_SIZES.$s] 0]
    keylset $w.PS PAGE_HEIGHT [lindex [keylget tk_utils_defs PS.PAGE_SIZES.$s] 1]
    keylset $w.PS PAGE_WIDTH [lindex [keylget tk_utils_defs PS.PAGE_SIZES.$s] 2]

    if {[string compare "Landscape" [keylget $w.PS ORIENTATION]] == 0} {
	set temp [keylget $w.PS PAGE_HEIGHT]
	keylset $w.PS PAGE_HEIGHT [keylget $w.PS PAGE_WIDTH]
	keylset $w.PS PAGE_WIDTH $temp
    }

    $w.page.size configure -text [lindex [keylget tk_utils_defs PS.PAGE_SIZES.$s] 0]
}

proc ps_orient_change {w s} {
    global $w.PS

    if {[string compare [keylget $w.PS ORIENTATION] $s] != 0} {
	set temp [keylget $w.PS PAGE_HEIGHT]
	keylset $w.PS PAGE_HEIGHT [keylget $w.PS PAGE_WIDTH]
	keylset $w.PS PAGE_WIDTH $temp
    }

    keylset $w.PS ORIENTATION $s
    $w.page.orient configure -text $s
}

proc ps_font_name_change {w s} {
    global $w.PS

    keylset $w.PS FONT $s
    $w.font.name configure -text $s
}

proc ps_font_size_change {w s} {
    global $w.PS

    keylset $w.PS FONT_SIZE $s
    $w.font.size configure -text "$s point"
}

proc PrintSetup {} {
    global tk_utils_defs

    set w [keylget tk_utils_defs PS.WIN]
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Page Print Options"

    global $w.PS
    set $w.PS [keylget tk_utils_defs PS]

    # Page size and orientation.
    frame $w.page
    menubutton $w.page.size -width 9 -indicatoron 1 -text [keylget tk_utils_defs PS.PAGE] \
	-relief raised -menu $w.page.size.m -anchor w
    menu $w.page.size.m -tearoff 0
    foreach i [keylkeys tk_utils_defs PS.PAGE_SIZES] {
	$w.page.size.m add command -label [lindex [keylget tk_utils_defs PS.PAGE_SIZES.$i] 0] \
	    -command [list ps_pg_size_change $w $i]
    }

    menubutton $w.page.orient -width 9 -indicatoron 1 -text [keylget $w.PS ORIENTATION] \
	-relief raised -menu $w.page.orient.m -anchor w
    menu $w.page.orient.m -tearoff 0
    $w.page.orient.m add command -label "Portrait" -command [list ps_orient_change $w "Portrait"]
    $w.page.orient.m add command -label "Landscape" -command [list ps_orient_change $w "Landscape"]

    pack $w.page.size $w.page.orient -side left -expand 1 -fill x
    pack $w.page -side top -padx 2 -pady 2

    # Margin and panel sizes.
    frame $w.mp

    frame $w.mp.margin
    label $w.mp.margin.l -text "Margins"
    entrybox $w.mp.margin.top -title "Top:" \
	-type CheckInt -width 3 -default [keylget $w.PS TOP_MARGIN]
    entrybox $w.mp.margin.bottom -title "Bottom:" \
	-type CheckInt -width 3 -default [keylget $w.PS BOTTOM_MARGIN]
    entrybox $w.mp.margin.left -title "Left:" \
	-type CheckInt -width 3 -default [keylget $w.PS LEFT_MARGIN]
    entrybox $w.mp.margin.right -title "Right:" \
	-type CheckInt -width 3 -default [keylget $w.PS RIGHT_MARGIN]

    pack $w.mp.margin.l $w.mp.margin.top $w.mp.margin.bottom $w.mp.margin.left $w.mp.margin.right \
	-side top -expand 1 -fill x

    frame $w.mp.panel
    label $w.mp.panel.l -text "Panels"
    entrybox $w.mp.panel.height -title "Height:" \
	-type "CheckIntRange 0 [keylget $w.PS PAGE_HEIGHT]" \
 	-width 3 -default [keylget $w.PS PANEL_HEIGHT]
    entrybox $w.mp.panel.separation -title "Separation:" \
	-type CheckInt -width 3 -default [keylget $w.PS PANEL_SEPARATION]

    pack $w.mp.panel.l $w.mp.panel.height $w.mp.panel.separation \
	-side top -expand 1 -fill x

    pack $w.mp.margin $w.mp.panel -side left -expand 1 -fill x -padx 2 -pady 2
    pack $w.mp -side top -padx 2 -pady 2

    # Font types and sizes.
    frame $w.font

    menubutton $w.font.name -indicatoron 1 -text [keylget $w.PS FONT] \
	-relief raised -menu $w.font.name.m -anchor w
    menu $w.font.name.m -tearoff 0
    foreach i [keylget tk_utils_defs PS.FONTS] {
	$w.font.name.m add command -label $i \
	    -command [list ps_font_name_change $w $i]
    }

    menubutton $w.font.size -indicatoron 1 -text "[keylget $w.PS FONT_SIZE] point"\
	-relief raised -menu $w.font.size.m -anchor w
    menu $w.font.size.m -tearoff 0
    for {set i [keylget tk_utils_defs PS.FONT_MIN]} \
	{$i <= [keylget tk_utils_defs PS.FONT_MAX]} \
	{incr i} {
	$w.font.size.m add command -label $i \
	    -command [list ps_font_size_change $w $i]
    }

    pack $w.font.name $w.font.size -side left -expand 1 -fill x
    pack $w.font -side top -padx 2 -pady 2

    # Add OK/cancel/help buttons
    okcancelhelp $w.ok \
	-ok_command "keylset tk_utils_defs PS \[set $w.PS\]; PrintSetup2 $w; destroy $w" \
	-cancel_command "destroy $w" \
        -help_command "show_help trev {Trev-Print-PageOptions}"
    pack $w.ok -side top -expand 1 -fill x -in $w -padx 2 -pady 2
}

proc PrintSetup2 {w} {
    global tk_utils_defs trev_defs

    keylset tk_utils_defs PS.TOP_MARGIN [entrybox_get $w.mp.margin.top]
    keylset tk_utils_defs PS.BOTTOM_MARGIN [entrybox_get $w.mp.margin.bottom]
    keylset tk_utils_defs PS.LEFT_MARGIN [entrybox_get $w.mp.margin.left]
    keylset tk_utils_defs PS.RIGHT_MARGIN [entrybox_get $w.mp.margin.right]
    keylset tk_utils_defs PS.PANEL_HEIGHT [entrybox_get $w.mp.panel.height]
    keylset tk_utils_defs PS.PANEL_SEPARATION [entrybox_get $w.mp.panel.separation]
    keylset tk_utils_defs PS.PANEL_WIDTH [expr [keylget tk_utils_defs PS.PAGE_WIDTH] \
					  - [keylget tk_utils_defs PS.LEFT_MARGIN] \
					  - [keylget tk_utils_defs PS.RIGHT_MARGIN]]


    # Limit margins to sensible values, 2-inches=144 points
    set max_margin 144
    if {[keylget tk_utils_defs PS.TOP_MARGIN] > $max_margin} {
	keylset tk_utils_defs PS.TOP_MARGIN $max_margin
    }
    if {[keylget tk_utils_defs PS.BOTTOM_MARGIN] > $max_margin} {
	keylset tk_utils_defs PS.BOTTOM_MARGIN $max_margin
    }
    if {[keylget tk_utils_defs PS.LEFT_MARGIN] > $max_margin} {
	keylset tk_utils_defs PS.LEFT_MARGIN $max_margin
    }
    if {[keylget tk_utils_defs PS.RIGHT_MARGIN] > $max_margin} {
	keylset tk_utils_defs PS.RIGHT_MARGIN $max_margin
    }


    # Avoid divide by zero below
    if {[keylget tk_utils_defs PS.PANEL_HEIGHT]==0} {
	keylset tk_utils_defs PS.PANEL_HEIGHT 1
    }
    keylset tk_utils_defs PS.N_PANEL [expr ([keylget tk_utils_defs PS.PAGE_HEIGHT] \
					    - [keylget tk_utils_defs PS.TOP_MARGIN] \
					    - [keylget tk_utils_defs PS.BOTTOM_MARGIN]) \
				      / ([keylget tk_utils_defs PS.PANEL_HEIGHT] \
					     + [keylget tk_utils_defs PS.PANEL_SEPARATION])]


    set trace [keylget trev_defs TRACE.WIN].t
    set status [catch {$trace ps_configure \
			   -page_height		[keylget tk_utils_defs PS.PAGE_HEIGHT] \
			   -page_width		[keylget tk_utils_defs PS.PAGE_WIDTH] \
			   -orientation		[keylget tk_utils_defs PS.ORIENTATION] \
			   -top_margin		[keylget tk_utils_defs PS.TOP_MARGIN] \
			   -bottom_margin	[keylget tk_utils_defs PS.BOTTOM_MARGIN] \
			   -left_margin		[keylget tk_utils_defs PS.LEFT_MARGIN] \
			   -right_margin	[keylget tk_utils_defs PS.RIGHT_MARGIN] \
			   -panel_height	[keylget tk_utils_defs PS.PANEL_HEIGHT] \
			   -panel_width		[keylget tk_utils_defs PS.PANEL_WIDTH] \
			   -panel_separation	[keylget tk_utils_defs PS.PANEL_SEPARATION] \
			   -n_panel		[keylget tk_utils_defs PS.N_PANEL] \
			   -font		[keylget tk_utils_defs PS.FONT] \
			   -font_size		[keylget tk_utils_defs PS.FONT_SIZE] \
		       } err]
#    set status [catch {$trace ps_configure}]

    if {$status != 0} {
	tk_messageBox -icon error -type ok -title "Print Setup Failed" \
	    -message "Unable to configure print options: $err"
    }
}

##############################################################################
proc ps_line_width_change {line s} {
    global $line.OPTIONS

    keylset $line.OPTIONS LW $s
    $line.eg itemconfigure line -width $s
}

proc ps_line_dash_change {line s} {
    global $line.OPTIONS

    keylset $line.OPTIONS DASH $s
    $line.dash configure -text "$s"
}

proc ps_line_configure {line name} {
    global tk_utils_defs $line.OPTIONS

    canvas $line.eg -width 3.0c -height 0.6c
    $line.eg create line 0.5c 0.3c 2.5c 0.3c \
	-width [keylget $line.OPTIONS LW] \
	-fill [keylget $line.OPTIONS COLOUR] \
	-tags line

    label $line.l -text "$name"
    set ok_cmd "global $line.OPTIONS; keylset $line.OPTIONS COLOUR"
    set update_cmd "$line.eg itemconfigure line -fill"
    set cancel_cmd ""
    button $line.colour -text "Colour" \
	-command "if {[winfo exists $line.cbox]} {raise $line.cbox; return}; \
                  xtoplevel $line.cbox; \
                  ColourBox $line.cbox \[keylget $line.OPTIONS COLOUR\] \
                  {$ok_cmd} {$update_cmd} {$cancel_cmd}"
    label $line.width_l -text "Line width:"
    scale $line.width -from [keylget tk_utils_defs PS.LW_MIN] -to [keylget tk_utils_defs PS.LW_MAX] \
	-length 3.0c -orient horizontal \
	-command "ps_line_width_change $line"
    $line.width set [keylget $line.OPTIONS LW]
    label $line.dash_l -text "Dash pattern:"
    menubutton $line.dash -indicatoron 1 -text "[keylget $line.OPTIONS DASH]" \
	-relief raised -menu $line.dash.m
    menu $line.dash.m -tearoff 0
    foreach i [keylget tk_utils_defs PS.DASHES] {
	$line.dash.m add command -label $i -command "ps_line_dash_change $line \"$i\""
    }

    pack $line.l $line.colour $line.width_l $line.width $line.dash_l $line.dash $line.eg \
	-side left -padx 2 -pady 2
}

proc PrintTraceSetup {} {
    global tk_utils_defs trev_defs
    global file_list file_ind file_status

    set filename [lindex $file_list [expr $file_ind - 1]]
    set status $file_status($filename)
    if {$status} {
        set status_l ""
	set label_prefix "Accepting"
    } else {
        set status_l " (rejected)"
	set label_prefix "Rejecting"
    }

    # Create a dialogue window
    set w [keylget tk_utils_defs PS.TRACE_WIN]
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Trace Print Options"

    entrybox $w.title -title "Title: " \
	-width 40 -default "[file tail $filename]$status_l"

    frame $w.a
    global $w.a.OPTIONS
    set $w.a.OPTIONS [keylget tk_utils_defs PS.A]
    ps_line_configure $w.a "A"

    frame $w.c
    global $w.c.OPTIONS
    set $w.c.OPTIONS [keylget tk_utils_defs PS.C]
    ps_line_configure $w.c "C"

    frame $w.g
    global $w.g.OPTIONS
    set $w.g.OPTIONS [keylget tk_utils_defs PS.G]
    ps_line_configure $w.g "G"

    frame $w.t
    global $w.t.OPTIONS
    set $w.t.OPTIONS [keylget tk_utils_defs PS.T]
    ps_line_configure $w.t "T"

    set trace [keylget trev_defs TRACE.WIN].t

    frame $w.range
    entrybox $w.range.first -title "Print bases: " \
	-type "CheckIntRange [lindex [$trace visible_region] 0] [lindex [$trace visible_region] 1]" \
 	-width 4 -default [lindex [$trace visible_region] 0]
    entrybox $w.range.last -title "to" \
	-type "CheckIntRange [lindex [$trace visible_region] 0] [lindex [$trace visible_region] 1]" \
 	-width 4 -default [lindex [$trace visible_region] 1]
    button $w.range.all -text "All" -command \
	"entrybox_configure $w.range.first -default \[lindex \[$trace visible_region\] 0\]; \
         entrybox_configure $w.range.last -default \[lindex \[$trace visible_region\] 1\];"
    button $w.range.visible -text "Visible" -command \
	"entrybox_configure $w.range.first -default \[lindex \[$trace visible_region\] 2\]; \
         entrybox_configure $w.range.last -default \[lindex \[$trace visible_region\] 3\]"
    pack $w.range.first $w.range.last -side left
    pack $w.range.all $w.range.visible -side left -padx 4

    pack $w.title $w.a $w.c $w.g $w.t $w.range -side top -anchor w

    okcancelhelp $w.ok \
	-ok_command "keylset tk_utils_defs PS.A \[set $w.a.OPTIONS\]; \
                     keylset tk_utils_defs PS.C \[set $w.c.OPTIONS\]; \
                     keylset tk_utils_defs PS.G \[set $w.g.OPTIONS\]; \
                     keylset tk_utils_defs PS.T \[set $w.t.OPTIONS\]; \
                     PrintTraceSetup2 $w; destroy $w" \
	-cancel_command "destroy $w" \
        -help_command "show_help trev {Trev-Print-TraceOptions}"
    pack $w.ok -side bottom -fill x -in $w
}

proc PrintTraceSetup2 {w} {
    global trev_defs tk_utils_defs

    set trace [keylget trev_defs TRACE.WIN].t
    set status [catch {$trace ps_configure_trace \
			   -title      [entrybox_get $w.title] \
			   -first_base [entrybox_get $w.range.first] \
			   -last_base  [entrybox_get $w.range.last] \
		       } err]
    set status [expr $status | [catch {$trace ps_configure_trace_line A \
					   -line_width	[keylget tk_utils_defs PS.A.LW] \
					   -colour	[keylget tk_utils_defs PS.A.COLOUR] \
					   -dash	[keylget tk_utils_defs PS.A.DASH] \
				       } err]]
    set status [expr $status | [catch {$trace ps_configure_trace_line C \
					   -line_width	[keylget tk_utils_defs PS.C.LW] \
					   -colour	[keylget tk_utils_defs PS.C.COLOUR] \
					   -dash	[keylget tk_utils_defs PS.C.DASH] \
				       } err]]
    set status [expr $status | [catch {$trace ps_configure_trace_line G \
					   -line_width	[keylget tk_utils_defs PS.G.LW] \
					   -colour	[keylget tk_utils_defs PS.G.COLOUR] \
					   -dash	[keylget tk_utils_defs PS.G.DASH] \
				       } err]]
    set status [expr $status | [catch {$trace ps_configure_trace_line T \
					   -line_width	[keylget tk_utils_defs PS.T.LW] \
					   -colour	[keylget tk_utils_defs PS.T.COLOUR] \
					   -dash	[keylget tk_utils_defs PS.T.DASH] \
				       } err]]

    if {$status != 0} {
	tk_messageBox -icon error -type ok -title "Print Trace Setup Failed" \
	    -message "Unable to configure print trace options: $err"
    }
}

##############################################################################
#
proc ymag_proc {w v} {
    $w ymag [expr double($v) / 10]
}

##############################################################################
proc CreateTrevWin {xmag ymag pregap_mode block_scf edit trace_scale} {
    global trev_menu trev_defs VERSION
    keylset trev_defs MENU.FRAME .menu
    set mf [keylget trev_defs MENU.FRAME]

    #main window parameters
    wm minsize . 600 180
    wm title . "Trev $VERSION"
    set f [frame [keylget trev_defs TREV.WIN]]

  # pack [frame [keylget trev_defs MENU.FRAME] -relief raised -bd 2] -fill both
  # create_menus $trev_menu [keylget trev_defs MENU.FRAME]
    . configure -menu $mf
    menu $mf
    create_menus $trev_menu $mf

    frame .place_holder -bd 0 -width 0 -height 0
    pack .place_holder -side top

    CreatePrevNext $pregap_mode

    if {$pregap_mode} {
	# Disable "save as" and "load" when called from pregap
        menu_state_set trev_menu -16 $mf

	# Add Accept and Reject buttons
	button .buttons.reject -text "Reject file" -command RejectFile

        pack .buttons.reject -side left -fill y -padx 15
    }

    fileTypes $block_scf
    SetEditFlag 0

    frame $f.trace
    CreateTraceWin [keylget trev_defs TRACE.WIN] $xmag $ymag $trace_scale
    pack $f.trace -side top -fill both -expand yes
    pack $f -fill both -expand yes

    SetBinding $edit

    catch {dnd bindtarget $f.trace.t text/uri-list <Drop> {OpenDroppedFiles %D}}

    bind . <Control-KeyPress-s> "
	if {\[lindex \[$f.trace.t configure -showedits\] 4\] == 1} {
	    SetBinding Seq
	}
    "
    bind . <Control-KeyPress-l> {SetBinding LCutOff}
    bind . <Control-KeyPress-r> {SetBinding RCutOff}

    # Hidden complement function. Disabled from normal view as it can cause
    # corruptions or crashes if editing is performed in this mode. You have
    # been warned...
    bind .trev.trace.t <Control-c><Control-o><Control-m><Control-p> {
	 .trev.trace.t complement
    }
}

proc CreatePrevNext {{pregap_mode 0}} {
    global file_count

    frame .buttons -bd 2 -relief raised
    pack .buttons -side top -fill both -after .place_holder

    if {$file_count > 1} {
	button .buttons.next -text "Next file >>"     -command NextFile \
		-state disabled
	button .buttons.prev -text "<< Previous file" -command PrevFile \
		-state disabled
	button .buttons.goto -text "Goto file ..."    -command GotoFile
	pack .buttons.prev .buttons.next .buttons.goto -side left -fill y

    }

    if {$file_count > 1 || $pregap_mode} {
	label .buttons.label
        pack .buttons.label -side right -fill both -padx 5
    }

    if {$file_count > 1} {
	catch {
	    label .buttons.label
	    pack .buttons.label -side right -fill both -padx 5
	}
    }
}

###########################################################################
proc CreateTraceWin { trace xmag ymag trace_scale} {
    global trev_defs tdisp

    dnatrace $trace.t \
	-relief ridge -bd 3 \
	-xscrollcommand "$trace.s set" \
	-shownumbers $tdisp(n) \
	-showsequence $tdisp(s) \
	-showedits $tdisp(e) \
	-showtrace $tdisp(t) \
	-showconf $tdisp(c) \
	-showends 0 \
	-trace_scale $trace_scale \
	-style $tdisp(style)
    focus $trace.t

    scrollbar $trace.s -command "$trace.t xview" -orient horizontal

    frame $trace.xmag
    label $trace.xmag.l -text "X"
    scale $trace.xmag.s -from 0 -to 3000 -command "$trace.t xmag" \
	    -showvalue 0
    $trace.xmag.s set $xmag
    after idle "$trace.t xview 0"

    frame $trace.ymag
    label $trace.ymag.l -text "Y"
    scale $trace.ymag.s -from 10 -to 500 -showvalue 0 -command \
	    "ymag_proc $trace.t"
    $trace.ymag.s set $ymag

    pack $trace.xmag.l $trace.ymag.l -side top -fill both
    pack $trace.xmag.s $trace.ymag.s -side bottom -fill both -expand 1
    pack $trace.xmag $trace.ymag -side left -fill both
    pack $trace.s -fill x
    pack $trace.t -fill both -expand yes
}

##############################################################################
# Display options
proc SetTraceDisplay {aname} {
    global trev_defs trev_menu
    upvar #0 $aname tdisp

    set trace [keylget trev_defs TRACE.WIN].t
    $trace configure -shownumbers $tdisp(n)
    $trace configure -showsequence $tdisp(s)
    $trace configure -showedits $tdisp(e)
    $trace configure -showtrace $tdisp(t)
    $trace configure -showconf $tdisp(c)

    #configure the edit menu depending on whether the editted sequence is
    #displayed
    if {$tdisp(e) == 0} {
	menu_state_set trev_menu -32 [keylget trev_defs MENU.FRAME]
	if {[keylget trev_defs DEFAULT_EDIT] == "Seq"} {
	  SetBinding RCutOff
	} else {
	  SetBinding [keylget trev_defs DEFAULT_EDIT]
	}
    } else {
	menu_state_set trev_menu 32 [keylget trev_defs MENU.FRAME]
    }
}
#end SetTraceDisplay


##############################################################################
#return 0 if cursor not on screen
#return 1 if cursor is on screen
proc SeeCursor { w } {
    if {[expr [lindex [$w position e] 0] + [lindex [$w position e] 1] + 1< \
	    [$w icursor]] || [lindex [$w position e] 0] > [$w icursor]} {
	return 0
    } else {
	return 1
    }
}
#end SeeCursor

##############################################################################
#redisplay the trace so that the cursor is in the centre of the display
proc FindCursor { w } {
    set width [lindex [$w position e] 1]

    if {[expr [lindex [$w position e] 0] + [lindex [$w position e] 1] < \
	    [$w icursor]] || [lindex [$w position e] 0] >= [$w icursor]} {
	$w xview %[expr [$w icursor] - ($width/2) ]
    }
}
#end FindCursor

##############################################################################
# search procedure
proc SearchProc { } {

    #check have not already got a search window open
    #if {[info commands .search] != ""} {return}
    set s .search
    if {[xtoplevel $s -resizable 0] == ""} return
    wm title $s "Search"

    label $s.label -text "Search for: "
    entry $s.entry -width 30 -relief sunken -bd 2 -textvariable tmp_search

    $s.entry delete 0 end

    button $s.next -text "Next" -command  \
	    {FindNext [string toupper $tmp_search]}
    button $s.previous -text "Previous" -command \
	    {FindPrevious [string toupper $tmp_search]}
    button $s.cancel -text "Cancel" -command "destroy $s"

    bind $s.entry <Return> {FindNext [string toupper $tmp_search]}
    focus $s.entry

    pack $s.label $s.entry $s.next $s.previous $s.cancel \
	    -side left -padx 1m -pady 2m

}
#end SearchProc

##############################################################################
#find the previous occurence of sub_str in str
proc FindPrevious {sub_str } {
    global trev_defs

    set trace [keylget trev_defs TRACE.WIN].t
    set str [string toupper [$trace sequence]]
    set str_pos [$trace icursor]
    set prev_str [string range $str 0 [expr $str_pos - 1]]

    #find the last position of "sub_str" in "prev_str"
    set prev_pos [string last $sub_str $prev_str]

    #if there is an occurence of sub_str in prev_str
    if { $prev_pos >= 0} {

	#position the cursor in the original string
	set str_pos $prev_pos
	$trace icursor $str_pos

	#position the window around the cursor
	FindCursor $trace
    } else {
	bell
    }
}
#end FindPrevious

##############################################################################
#find the next occurence of sub_str in str
proc FindNext { sub_str } {
    global trev_defs

    set trace [keylget trev_defs TRACE.WIN].t
    set str [string toupper [$trace sequence]]
    set str_len [string length $str]
    set str_pos [$trace icursor]
    set remaining_str [string range $str [expr $str_pos + 1] $str_len]

    #find the first position of "sub_str" in "remaining_str"
    set rem_pos [string first $sub_str $remaining_str]

    #if sub_str is in the remaining_str
    if { $rem_pos >= 0 } {

	#find the length of the "remaining_str"
	set rem_len [string length $remaining_str]

	#find the position of the sub_str in the original string
	set str_pos [expr $str_len - $rem_len + $rem_pos]

	#position the cursor in the original string
	$trace icursor $str_pos
	FindCursor $trace

    } else {
	bell
    }
}
#end FindNext


##############################################################################
# Assign edit menu options to mouse button 1
proc SetBinding { ledit } {
    global edit history trev_menu trev_defs
    set edit $ledit

    if {$edit == "Seq"} {
	bind DnaTrace <<select>> {
	    %W icursor @%x
	}
	bind DnaTrace <Delete> {
	    if {[%W loaded]} {
		if {![SeeCursor %W]} {
		    FindCursor %W
		} else {
		    %W delete
		    SetEditFlag 1
		}
	    }
	}
	bind DnaTrace <Any-KeyPress> {
	    if {[%W loaded] && [lsearch {a A c C g G t T -} %A] != -1} {
		if {![SeeCursor %W]} {
		    FindCursor %W
		} else {
		    %W insert %A
		    SetEditFlag 1
		}
	    }
	}
	bind DnaTrace <Right> {
	    %W icursor [expr [%W icursor] + 1]
	    FindCursor %W
	}
	bind DnaTrace <Left> {
	    %W icursor [expr [%W icursor] - 1]
	    FindCursor %W
	}

    } else {

	bind DnaTrace <Delete> {}
	bind DnaTrace <Any-KeyPress> {}
	bind DnaTrace <Right> {}
	bind DnaTrace <Left> {}

	if {$edit == "LCutOff"} {
	    bind DnaTrace <<select>> {
		if {[%W loaded]} {
		    lappend history "%W left_cutoff [%W left_cutoff]"
		    menu_state_set trev_menu 64 [keylget trev_defs MENU.FRAME]
		    %W left_cutoff @%x
		    SetEditFlag 1
		}
	    }
	} elseif {$edit == "RCutOff"} {
	    bind DnaTrace <<select>> {
		if {[%W loaded]} {
  		    lappend history "%W right_cutoff [%W right_cutoff]"
		    menu_state_set trev_menu 64 [keylget trev_defs MENU.FRAME]
		    %W right_cutoff @%x
		    SetEditFlag 1
		}
	    }
	} elseif {$edit == "LVector"} {
	    bind DnaTrace <<select>> {
		if {[%W loaded]} {
		    lappend history "%W left_vector [%W left_vector]"
		    menu_state_set trev_menu 64 [keylget trev_defs MENU.FRAME]
		    %W left_vector @%x
		    SetEditFlag 1
		}
	    }
	} elseif {$edit == "RVector"} {
	    bind DnaTrace <<select>> {
		if {[%W loaded]} {
		    lappend history "%W right_vector [%W right_vector]"
		    menu_state_set trev_menu 64 [keylget trev_defs MENU.FRAME]
		    %W right_vector @%x
		    SetEditFlag 1
		}
	    }
	}
    }
}
#end SetBinding

##############################################################################
# display trace information
proc DisplayInfo { } {
    global trev_defs

    set trace [keylget trev_defs TRACE.WIN].t
    set string [$trace info]
    set i .info

    #if an info box is already displayed then update with info on current seq
    if {[winfo exists $i]} {
	raise $i
	wm deiconify $i
	$i.top.info configure -text $string
	return
    }

    xtoplevel $i
    wm title $i Information

    frame $i.top -relief raised -bd 1
    pack $i.top -side top -fill both -expand 1
    frame $i.bot -relief raised -bd 1
    pack $i.bot -side bottom -fill both

    message $i.top.info -text $string
    pack $i.top.info -fill both -expand 1

    button $i.bot.button -text "OK" -command "destroy $i"
    pack $i.bot.button
}
#end DisplayInfo


##############################################################################
# Undo - for cutoff adjustment only at present
proc UndoClipping {} {
    global history trev_menu trev_defs

    set h [lindex $history end]
    if {$h != ""} {
	eval $h
	SetEditFlag 1
	set history [lreplace $history end end]
	if {$history == ""} {
	    menu_state_set trev_menu -64 [keylget trev_defs MENU.FRAME]
	}
    }
}

##############################################################################
# Multiple file handling

# Skip to next and previous trace files (specified on command line)
proc NextFile {{format any}} {
    global file_list file_ind file_count trev_defs

    set f [lindex $file_list $file_ind]
    incr file_ind

    if {[SwitchFile $f $format] == 0} {
	incr file_ind -1
	return; # cancel pressed
    }

    SwitchFile2
}

proc PrevFile {{format Any}} {
    global file_list file_ind file_count trev_defs

    incr file_ind -1
    set f [lindex $file_list [expr $file_ind-1]]

    if {[SwitchFile $f $format] == 0} {
	incr file_ind
	return; # cancel pressed
    }

    SwitchFile2
}

proc SwitchFile {f {format Any}} {
    global pregap_mode

    if {[CheckSaved]} {
	LoadFile $f $format $pregap_mode
	return 1
    } else {
	return 0
    }
}

proc SwitchFile2 {} {
    global file_ind file_count

    # Update the "Goto File" list
    UpdateFileList

    # Update Next and Prev buttons
    if {$file_count <= 1} {
	return
    }

    if {$file_ind == 1} {
	.buttons.prev configure -state disabled
    }
    if {$file_ind == $file_count} {
	.buttons.next configure -state disabled
    }
    if {$file_ind > 1} {
	.buttons.prev configure -state normal
    }
    if {$file_ind < $file_count} {
	.buttons.next configure -state normal
    }
}

proc GotoFile {{format any}} {
    global file_list file_status file_ind file_count

    set w .file_list
    if {[xtoplevel $w] == ""} return

    # X & Y scrollbars
    scrollbar $w.xs -orient hori -command "$w.t xview"
    scrollbar $w.ys -orient vert -command "$w.t yview"

    # Text widget with tags and tag bindings
    text $w.t -width 40 -height 20 -xscrollcommand "$w.xs set" \
	-yscrollcommand "$w.ys set" -wrap none
    $w.t tag configure current -foreground blue
    $w.t tag configure selectable
    $w.t tag configure highlight -relief raised -borderwidth 2

    $w.t tag bind selectable <<select>> "
       	set ind \[%W index @%x,%y\]
	set ind2 \[%W tag nextrange selectable \
		\"\$ind linestart\" \"\$ind lineend\"\]
	set f \[%W get \[lindex \$ind2 0\] \[lindex \$ind2 1\]\]
	global file_ind
	regsub {\\..*} \$ind {} file_ind
        if {\[SwitchFile \$f $format\] == 0} {
	    incr file_ind -1
	    return; # cancel pressed
        }
	SwitchFile2
    "

    $w.t tag bind selectable <Any-Motion> {
       	set ind [%W index @%x,%y]
	%W tag remove highlight 1.0 end
	%W tag add highlight "$ind linestart" "$ind lineend"
    }

    # Remove most default Text bindings
    bindtags $w.t "$w.t . all"
    bind $w.t <2> [bind Text <2>]
    bind $w.t <B2-Motion> [bind Text <B2-Motion>]
    bind $w.t <ButtonRelease-2> [bind Text <ButtonRelease-2>]

    button $w.ok -text OK -command "destroy $w"

    grid rowconfigure $w 0 -weight 1
    grid columnconfigure $w 0 -weight 1
    grid $w.t $w.ys -sticky nsew
    grid $w.xs -sticky nsew
    grid $w.ok

    UpdateFileList
}

proc UpdateFileList {} {
    global file_list file_ind file_status file_count file_info

    set t .file_list.t
    if {![winfo exists $t]} {
	return
    }

    $t delete 1.0 end
    set count 0
    foreach f $file_list {
	incr count
	set tags "selectable"
	set st ""
	set en ""
	if {$count == $file_ind} {
	    lappend tags "current"
	}
	if {$file_status($f) == 0} {
	    set en " (rejected)"
	}
	# FIXME: set file_info($f) to be "LANE XX".
	if {[info exists file_info($f)]} {
	    set st "\[$file_info($f)\]\t"
	}
	$t insert end $st {} $f $tags "$en\n"
    }
}

##############################################################################
# Pregap 'reject' and 'accept' commands

proc SetRejectStatus {} {
    global file_list file_ind file_status trev_defs VERSION
    global pregap_mode

    set filename [lindex $file_list [expr $file_ind-1]]
    set status $file_status($filename)
    if {$status} {
        set status_l ""
	set label_prefix "Accepting"
    } else {
        set status_l " (rejected)"
	set label_prefix "Rejecting"
    }
    if {!$pregap_mode} { set label_prefix "" }
    wm title . "Trev $VERSION: [file tail $filename]$status_l"

    if {[winfo exists .buttons.label]} {
	global file_ind file_count
	.buttons.label configure \
		-text "$label_prefix file $file_ind of $file_count"
    }

    if {!$pregap_mode} {
	return
    }

    if {$status} {
	# Accepted
	.buttons.reject configure -text "Reject file" -command RejectFile
    } else {
	# Rejected
	.buttons.reject configure -text "Accept file" -command AcceptFile
    }
}

proc RejectFile {} {
    global file_list file_ind file_count file_status

    set file_status([lindex $file_list [expr $file_ind-1]]) 0
    SetRejectStatus

    if {$file_ind != $file_count} {
        NextFile
    }
}

proc AcceptFile {} {
    global file_list file_ind file_count file_status

    set file_status([lindex $file_list [expr $file_ind-1]]) 1
    SetRejectStatus

    if {$file_ind != $file_count} {
        NextFile
    }
}

##############################################################################
# Exit handling

proc CheckSaved {} {
    global trev_defs

    if {[GetEditFlag] == 1} {
        set trace [keylget trev_defs TRACE.WIN].t

	#editted previous sequence
	set response [tk_messageBox -icon question -type yesnocancel \
		-title "File changed" \
		-message "Do you wish to save your file before loading a new \
		file?"]

	if {[string compare $response "yes"] ==0 } {
	    if {[CheckSaveFormat [$trace orig_format]]} {
		#Save
		SaveFile 0 [$trace orig_format]
	    } else {
		#Save As
		SaveFile 1 ZTR
	    }
	}
	if {[string compare $response "cancel"] == 0} {
	    return 0
	}
        if {[string compare $response "no"] == 0} {
            SetEditFlag 0
        }
    }

    return 1
}

proc ExitTrev {} {
    global file_list file_status pregap_mode

    if {[CheckSaved]} {
	if {$pregap_mode == 1} {
	  catch {
	    puts -nonewline STATUS
	    foreach f $file_list {
		puts -nonewline " [list [list $f $file_status($f)]]"
	    }
	    puts ""
	  }
	} elseif {$pregap_mode == 2} {
	    exit [expr 1-$file_status([lindex $file_list 0])]
	}
        exit 0
    }
}

##############################################################################
#set up global variables
##############################################################################
proc InitGlobals { } {
    global edit
    global dir
    global pregap_mode
    global block_scf
    global history
    global tdisp
    global trev_defs

    set edit [keylget trev_defs DEFAULT_EDIT]
    set tdisp(n) [keylget trev_defs SHOW_NUMBERS]
    set tdisp(s) [keylget trev_defs SHOW_SEQUENCE]
    set tdisp(e) [keylget trev_defs SHOW_EDITS]
    set tdisp(t) [keylget trev_defs SHOW_TRACE]
    set tdisp(c) [keylget trev_defs SHOW_CONFIDENCE]
    set tdisp(style) 0
    set pregap_mode 0
    set block_scf 1
    set history {}

    #default directory to current directory
    set dir [pwd]

}
#end InitGlobals

##############################################################################
#                                  main program                              #
##############################################################################

source $env(STADTABL)/shlib.conf
load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}

load_package tk_utils
tk_utils_init
load_package trev

set VERSION "1.9$svn_version"

InitGlobals
set xmag 150
set xmag_changed 0
set ymag 10
set filename ""
set informat Any
set file_list ""
set trace_scale 0

# Parse command line arguments
while {$argc >= 1} {
    set arg [string tolower [lindex $argv 0]]
    if {$arg == "-pregap_mode"} {
	set pregap_mode 1
	set tdisp(e) 0
	set block_scf 0
    } elseif {$arg == "-restrict"} {
	# Old pregap mode
	set pregap_mode 2
	set tdisp(e) 0
	set block_scf 0
    } elseif {$arg == "-editscf"} {
	set block_scf 0
    } elseif {$arg == "-xmag"} {
	set xmag [lindex $argv 1]
	set xmag_changed 1
        incr argc -1
	set argv [lrange $argv 1 end]
    } elseif {$arg == "-ymag"} {
	set ymag [lindex $argv 1]
        incr argc -1
	set argv [lrange $argv 1 end]
    } elseif {$arg == "-trace_scale"} {
	set trace_scale [lindex $argv 1]
        incr argc -1
	set argv [lrange $argv 1 end]
    } elseif {$arg == "-style"} {
	set tdisp(style) [lsearch {chroma filled pyro stick} [lindex $argv 1]]
        incr argc -1
	set argv [lrange $argv 1 end]
    } elseif {$arg == "-write_scf"} {
	set block_scf 0
    } elseif {$arg == "-abi" ||
	      $arg == "-alf" ||
	      $arg == "-ctf" ||
	      $arg == "-exp" ||
	      $arg == "-scf" ||
	      $arg == "-pln" ||
	      $arg == "-ztr" ||
	      $arg == "-any"} {
	regsub -- - $arg {} informat
    } elseif {$arg == "-fofn"} {
	incr argc -1
	set argv [lrange $argv 1 end]
	set fd [open [lindex $argv 0]]
	while {[gets $fd line] != -1} {
	    lappend file_list $line
	}
	close $fd
    } elseif {$arg == "--"} {
	incr argc -1
	set argv [lrange $argv 1 end]
	break;
    } elseif {[string match "-*" $arg]} {
	usage
    } else {
        break
    }

    incr argc -1
    set argv [lrange $argv 1 end]
}

append file_list " $argv"
set file_ind 0
set file_count [llength $file_list]
foreach f $file_list {
    set file_status($f) 1
}

CreateTrevWin $xmag $ymag $pregap_mode $block_scf $edit $trace_scale

if {$file_count > 0} {
    NextFile $informat
}

# Enable the Delete key to work like backspace - essential for some keyboards
bind Entry <Delete> [bind Entry <BackSpace>]
bind Text <Delete> [bind Text <BackSpace>]
wm protocol . WM_DELETE_WINDOW {ExitTrev}

raise .
focus -force [keylget trev_defs TRACE.WIN].t
