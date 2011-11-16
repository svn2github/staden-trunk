#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc GapRestrictionEnzymeDialogue {io } {
    global tk_utils_defs gap5_defs

    set w [keylget tk_utils_defs R_ENZ.WIN]
    if {[xtoplevel $w] == ""} return
    wm title $w "Select restriction enzymes"

    global $w.list

    #contig id
    if {![contig_renzbox $w $io]} {
	#pressed cancel
	return
    }

    set contig_name [renzbox_contig_name $w]
    set list [set $w.list]
    set filename [renzbox_filename $w]
    set num_items [llength $list]

    GapCreateREnzDisplay $io $contig_name $list $filename [renzbox_from $w] [renzbox_to $w]
}

proc gap_resize {io id} {
    if {$id != -1} {
	resize_canvas -io $io -id $id
    }
}

proc GapCreateREnzDisplay {io contig_name list filename from to} {
    global tk_utils_defs 
    global CurContig

    set w .gap_r_enzyme_map
    #generate a new display
    set num_display [next_renz_display]
    set w $w$num_display
   
    global $w.renz_id
    set $w.renz_id -1

    set x_scroll "gc_scroll_x $io \[set $w.renz_id\]"
    set y_scroll "gc_scroll_y $io \[set $w.renz_id\]"
    set zoomback_command "gap_zoom $io \[set $w.renz_id\] b"
    set zoom_command "gap_zoom $io \[set $w.renz_id\] b"
    set cursor_command "draw_canvas_cursor_x -io $io -id \[set $w.renz_id\]"

    set invoke_command "canvas_cursor_editor $io \
                [db_info get_contig_num $io $CurContig] \[set $w.renz_id\]"

    set renz_name_command "get_r_enz_name -id \[set $w.renz_id\] -io $io -cnum \[db_info get_contig_num $io $CurContig\]"
    set renz_info_command "get_r_enz_info -id \[set $w.renz_id\] -io $io -cnum \[db_info get_contig_num $io $CurContig\]"

    set config_command "gap_resize $io \[set $w.renz_id\]"

    set min_win_ht [keylget tk_utils_defs R_ENZ.PLOT_HEIGHT]

    set num_enz [llength $list]
    set tick_ht [keylget tk_utils_defs R_ENZ.TICK_HEIGHT]

    #HACK - want a better way of deciding on the text offset
    set text_offset [expr $tick_ht * 1.5]

    #set optimal window height for names and plot
    set height [expr ($num_enz + 2) * $tick_ht]

    #HACK to ensure window does not get too big!
    if { $height > $min_win_ht} {
	set height $min_win_ht
    }

    renz_map $w -zoom_command $zoom_command \
	-zoomback_command $zoomback_command\
	-scrollbar_x_cmd $x_scroll -scrollbar_y_cmd $y_scroll \
	-cursor_cmd $cursor_command \
	-renz_name_cmd $renz_name_command \
	-renz_info_cmd $renz_info_command \
	-width [keylget tk_utils_defs R_ENZ.PLOT_WIDTH]\
	-height $height \
	-ruler_height [keylget tk_utils_defs R_ENZ.RULER.PLOT_HEIGHT]\
	-names_width [keylget tk_utils_defs R_ENZ.NAME_WIDTH] \
	-selectbackground [keylget tk_utils_defs R_ENZ.SELECT_COLOUR] \
	-config_cmd $config_command \
	-invoke_cmd $invoke_command \
	-tick_ht [keylget tk_utils_defs R_ENZ.TICK_HEIGHT] \
	-text_offset $text_offset \
	-text_fill [keylget tk_utils_defs R_ENZ.TEXT_COLOUR]
 
    #HACK - to add popup menu to gap restriction enzyme ruler
    bind $w.ruler <<menu>> "PopUpSingleContigMenu $io %W $contig_name %X %Y"

    # Main Menu Bar
    GapCreateREnzMenu $w $io $CurContig $renz_name_command

    wm protocol $w WM_DELETE_WINDOW "GapREnzStartShutdown $w $io"

#    update idletasks
#    6/1/99 johnt - need just update for windows
    update
    set $w.renz_id [renz_map_plot $w "plot_renz -enzymes {$list} -num_enzymes $num_enz -file {$filename} -contigs \"{$contig_name $from $to}\" -io $io]"]

    wm title $w "Restriction enzyme map: $CurContig #[db_info get_read_num $io $CurContig]"
    #update result list
    #seq_result_list_update [keylget tk_utils_defs RASTER.RESULTS.WIN]

    wm minsize $w [winfo width $w] [winfo height $w]
    wm maxsize $w [winfo screenwidth $w] [winfo height $w]
    
}

##############################################################################
#stand alone restriction enzyme display
#menu display
proc GapCreateREnzMenu {w io contig_name name_cmd} {
    global $w.renz_id read_only

    #set up a File menu
    menubutton $w.menubar.file -text "File" -menu $w.menubar.file.opts
    menu $w.menubar.file.opts
    $w.menubar.file.opts add command -label "Exit" \
	    -command "GapREnzStartShutdown $w $io"

#    menubutton $w.menubar.edit -text "Edit" -menu $w.menubar.edit.opts
#    menu $w.menubar.edit.opts
#    $w.menubar.edit.opts add command -label "Create tags" \
#	-command "CreateREnzTags $w $io $contig_name $w.renz $w.brief \[set $w.renz_id\] {$name_cmd}"

    #set up a View menu
    menubutton $w.menubar.view -text "View" -menu $w.menubar.view.opts
    menu $w.menubar.view.opts
    $w.menubar.view.opts add command -label "Results manager" \
	    -command "ListResults $io"

    #set up a Results menu
    menubutton $w.menubar.results -text "Results" \
	-menu $w.menubar.results.opts
    menu $w.menubar.results.opts -postcommand "GetREnzResults $io $w"
    
    #Help menu
    menubutton $w.menubar.help -text "Help" -menu [set m $w.menubar.help.opts]
    menu $m

    $m add command -label "Introduction" \
	-command "show_help gap4 {Restrict}"
    $m add command -label "Selecting Enzymes" \
	-command "show_help gap4 {Restrict-Selecting}"
    $m add command -label "Examining the Plot" \
	-command "show_help gap4 {Restrict-Examining}"
    $m add command -label "Reconfiguring the Plot" \
	-command "show_help gap4 {Restrict-Reconfig}"
#    $m add command -label "Creating Tags for Cut Sites" \
#	-command "show_help gap4 {Restrict-Tags}"

#    if {$read_only} {
#	$w.menubar.edit configure -state disabled
#    }

    #do the packing
    pack $w.menubar.file -side left
#    pack $w.menubar.edit -side left
    pack $w.menubar.view -side left
    pack $w.menubar.results -side left
    pack $w.menubar.help -side right
}

##############################################################################
proc GetREnzResults {io w } {
    global $w.renz_id
    
    $w.menubar.results.opts delete 1 end
    
    result_list_popup_single $io [set $w.renz_id]\
	[reg_get_ops -io $io -id [set $w.renz_id]] \
	$w.menubar.results.opts

}


##############################################################################
#executed when the exit command is chosen from the File menu of stand alone
#restriction enzyme plot
proc GapREnzStartShutdown {w io} {
    global $w.renz_id

    if {[info exists $w.renz_id]} {
	result_delete -io $io -id [set $w.renz_id]
    }
}
##############################################################################
#stand alone restriction enzyme display
#create tags from the currently selected enzymes
proc CreateREnzTags { f io contig_name plot status_line id name_cmd} {
    global NGTag gap5_defs
    global $f.sel_enz

    if [winfo exists $f.create_restriction_tags] {
	raise $f.create_restriction_tags
	return
    }

    set enz_list ""
    set enz_names ""

    if {![info exists $f.sel_enz] || [llength [set $f.sel_enz]] == 0} {
	bell
	tk_messageBox -icon error -type ok -title "Error" \
		-message "No selection has been made" \
		-parent $plot
	raise $f
	tkwait variable press_ok
    }

    #NB only enzymes with cuts can have tags created!
    foreach item [set $f.sel_enz] {
	set index [GetREnzIndex $plot $item]
	if {$index > -1} {
	    lappend enz_list $index
	    lappend enz_names [GetREnzName $f $plot $item $name_cmd]
	}
    }
    #check that there are enzymes in the list
    if {[string compare $enz_list ""] == 0} {
	bell
	return
    }

    set t [xtoplevel $f.create_restriction_tags -resizable 0]
    wm title $f.create_restriction_tags "Create tags"
    frame $t.left
    frame $t.right

    set count 0

    foreach i $enz_list {
	set def NONE
	set enz_name [lindex $enz_names $count]

	#set up enz name label and menu button
	frame $t.f$i
	label $t.f$i.l -text $enz_name
	for {set j 0} {$j < $NGTag(num_tags)} {incr j} {
	    if {[string match [string toupper $enz_name] [string toupper $NGTag($j,tagtype)]] == 1} {
		set def $j
		break
	    }
	}

	if {[string compare $def NONE] == 0} {
	    menubutton $t.f$i.mb -text $def -indicatoron 1 -menu $t.f$i.mb.opts -bd 2 -relief raised -width 23
	} else {
	    menubutton $t.f$i.mb -text "$NGTag($def,tagid): $NGTag($def,tagtype)" -indicatoron 1 -menu $t.f$i.mb.opts -bd 2 -relief raised -width 23
	}

	#set up menu of available tags
	set m [menu $t.f$i.mb.opts -tearoff 0]
	for {set j 0; set k 1} {$j < $NGTag(num_tags)} {incr j; incr k} {
	    $m add command \
		    -label "$NGTag($j,tagid): $NGTag($j,tagtype)"\
		    -command "$t.f$i.mb configure -text \"$NGTag($j,tagid): $NGTag($j,tagtype)\""
	    if {$k >= [keylget gap5_defs MAX_MENU_ITEMS]} {
		$m add cascade -label "More..." -menu $m.more
		set m [menu $m.more -tearoff 0]
		set k 0
	    }
	}

	pack $t.f$i.l -side left
	pack $t.f$i.mb -side right
	pack $t.f$i -fill x
	incr count
    }
    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $t.ok_cancel \
	    -ok_command "CreateREnzTags_OK_Pressed $t $f $io $contig_name {$enz_list} $status_line; destroy $t" \
	    -cancel_command "destroy $t" \
	    -help_command "show_help gap4 {Restrict-Tags}" \
	    -bd 2 \
	    -relief groove
    pack $t.ok_cancel -side bottom -fill both
}

##############################################################################
#stand alone restriction enzyme display
proc CreateREnzTags_OK_Pressed {t f io contig_name enz_list status_line} {
    global $f.renz_id

    set enz_id ""
    foreach i $enz_list {
	set name [lindex [$t.f$i.mb configure -text] 4]
	if {[string compare $name NONE] == 0} {
	    bell
	    tkwait variable done
	} else {
	    lappend enz_id $name
	}
    }

    set num_tags [create_renz_tags -io $io -contigs $contig_name \
	    -id [set $f.renz_id] -enum $enz_list -tag_types $enz_id]

    if {$num_tags == 1} {
	$status_line configure -text "$num_tags tag created "	
    } else {
	$status_line configure -text "$num_tags tags created "
    }
}

##############################################################################
#                    template restriction enzyme display                     #
##############################################################################

##############################################################################
#template display restriction enzyme plot
proc templateRenzPopup {f canvas X Y x io renz_id renz_info_cmd} {
    global $f.template_id

    #set contig [template_contig -io $io -id [set $f.template_id] -x $x]
    #PopUpREnzMenu $f $canvas $X $Y $io $renz_id $contig
    PopUpREnzMenu $f $canvas $X $Y $renz_info_cmd
}

##############################################################################
#bindings specific to the restriction enzyme plot on the template display
proc SetTemplateREnzBindings { f canvas label dist io renz_id contig renz_name_cmd renz_info_cmd} {
    global $f.template_id

    $canvas bind S <<select>> "FindDistance $f %W $label $dist cut"

    #any mouse motion - highlight nearest cut line
    $canvas bind S <Any-Motion> \
	"HighlightREnz $f %W $label $renz_name_cmd"
    $canvas bind S <Shift-Motion> {;}

    $canvas bind S <<menu>> "templateRenzPopup $f $canvas %X %Y %x $io $renz_id $renz_info_cmd"
	
    bind $canvas <Any-Leave> "delete_canvas_cursor -io $io -id [set $f.template_id]"
    bind $canvas <Any-Motion> "AddTemplateCrossHair $io [set $f.template_id] $f  $canvas %x"

    bind $canvas <<move-create>> "
	canvas_cursor_editor $io [db_info get_contig_num $io $contig] \
		[set $f.template_id] %x $canvas
    "
    bind $canvas <<use>> "
	canvas_cursor_editor $io [db_info get_contig_num $io $contig] \
		[set $f.template_id] %x $canvas
    "
}

##############################################################################
#shutdown procedure of template restriction enzyme plot
#called from C
proc DeleteTemplateREnzPlot {f re_win } {
    global gap5_defs
    global $f.num_r_win
    global $f.pos
    global $f.plot

    incr $f.num_r_win -1

    set r $f[keylget gap5_defs TEMPLATE.RENZ.WIN]
    set num [string range $re_win [string length $r] end]

    global config$f.renz$num $f.renz_t_id$num

    unset $f.renz_t_id$num
    unset config$f.renz$num
    if {[info exists $f.pos]} {
	unset $f.pos
    }
    if {[info exists $f.plot]} {
	unset $f.plot
    }

    #if deleting exitting gap4, $f may not exist
    if {![winfo exists $f]} {
	return
    }

    set new_height [expr [winfo height $f] - [winfo height $f.re$num]]
    $re_win delete all
    destroy $re_win
    destroy $f.re$num
    if {[set $f.num_r_win] == 0} {
	destroy $f.renz
	unset $f.num_r_win
    }
    if {[winfo exists $f]} {
	wm geometry $f [winfo width $f]x$new_height
    }
}

##############################################################################
#single line in template display
proc CreateTemplateREnzPlot {io f plot label dist num} {
    global gap5_defs
    global $f.template_id

     global $f.num_r_win

    if {![info exists $f.num_r_win]} {
	set $f.num_r_win 0
    }
    incr $f.num_r_win

   #puts "Start CREATERENZPLOT"
    #only scroll in x direction
    set scroll x
    set borderwidth [keylget gap5_defs TEMPLATE.BORDERWIDTH]
    set height [keylget gap5_defs TEMPLATE.RENZ.PLOT_HEIGHT]
    set width [winfo width $f[keylget gap5_defs TEMPLATE.TEMPLATE.WIN]]

    #allow -height and -width to have affect
    wm geometry $f {}
    
    if {![winfo exists $f.renz]} {
	frame $f.renz
    }
    frame $f.re$num -bd $borderwidth -relief groove
    set zoom_cmd [list "gap_zoom $io \[set $f.template_id\] $scroll"]
    canvasbox $plot -width $width -height $height \
	-borderwidth 0 \
	-highlightthickness 0 \
	-xscrollcommand "$f.hscroll set" \
	-zoom_command $zoom_cmd 

    #set toplevel geometry when resized window
    set new_height [expr [winfo height $f] + [winfo reqheight $plot] + \
	   2 * $borderwidth]

    grid $f.renz -row 5 -column 0 -sticky ew
    pack $f.re$num -in $f.renz -fill x -expand yes
    pack $plot -in $f.re$num -padx 5 -fill x -expand yes
    update idletasks

    #set toplevel geometry when resized window
    if {[winfo height $f] != 1} {
	wm geometry $f [winfo width $f]x$new_height
	#wm geometry $parent [winfo width $parent]x[winfo height $parent]
    }
    #SetCanvasBindings $io [set $f.template_id] $plot $scroll
    SetCanvasBindings $plot $zoom_cmd
}

##############################################################################
proc CreateTemplateREnzDisplay { io f renz_win template_id label dist \
	contig_list num contig} {
    global tk_utils_defs 
    global $f.renz_t_id$num
    global $renz_win.prev_item
    global $renz_win.item_num
    global gap5_defs

    set l [keylget tk_utils_defs R_ENZ.WIN]
    if {[xtoplevel $l] == ""} return
    wm title $l "Restriction enzyme map"
    global $l.list

    if {![renzbox $l]} {
	#pressed cancel
	return 0
    }

    set offset [keylget gap5_defs TEMPLATE.RENZ.TICK_HEIGHT]
    set $renz_win.prev_item 0
    set $renz_win.item_num 0

    #HACK to work round the horrible problem in zooming, turning ruler
    #off, scrolling and turning ruler on again.
    set scroll [lindex [$f.hscroll get] 0]
   
    if {![winfo exists $renz_win]} {
	CreateTemplateREnzPlot $io $f $renz_win $label $dist $num 
    }

    set list [set $l.list]

    #set up registration of contig
    set $f.renz_t_id$num [plot_template_renz -window $renz_win -enzymes $list \
            -num_enzymes [llength $list] -file [renzbox_filename $l] \
	    -yoffset $offset \
            -contigs $contig_list -io $io -frame $f -id $template_id]
    if {[set $f.renz_t_id$num] == -1} {
	DeleteTemplateREnzPlot $f $renz_win
    }

    uplevel #0 eval [$f.hscroll cget -command] moveto $scroll

    set renz_name_command [list "get_r_enz_name -id [set $f.renz_t_id$num] -io $io -cnum \[db_info get_contig_num $io $contig\]"]
    set renz_info_command [list "get_r_enz_info -id [set $f.renz_t_id$num] -io $io -cnum \[db_info get_contig_num $io $contig\]"]

    SetTemplateREnzBindings $f $renz_win $label $dist $io [set $f.renz_t_id$num] $contig $renz_name_command $renz_info_command
    return 1
}

##############################################################################
#template quality display
proc next_template_renz_display { } {
    global next_t_r_display

    if {![info exists next_t_r_display]} {
	set next_t_r_display 0
	return $next_t_r_display
    }
    incr next_t_r_display
    return $next_t_r_display
}

proc unset_next_template_renz_display { } {
    global next_t_r_display

    unset next_t_r_display
}

##############################################################################
#template menu option to display restriction enzyme plot
proc TemplateREnzPlot { io f id label dist contig num} {
    global gap5_defs
    global $f.contig_list 
    global $f.renz_t_id$num
    global config$f.renz$num

    set renz_win $f[keylget gap5_defs TEMPLATE.RENZ.WIN]$num

    if {[set config$f.renz$num]} {
	
	set config$f.renz$num 0

	if {![template_win_free -io $io -id $id]} {
	    verror ERR_WARN template_display "no more windows available"
	    bell
	    set config$f.renz$num 0
	    return
	}
	if {[CreateTemplateREnzDisplay $io $f $renz_win $id $label $dist \
		$contig $num $contig] == 1} {
	    set config$f.renz$num 1
	}
    } else {
	delete_window -io $io -id $id -window $renz_win
	result_delete -io $io -id [set $f.renz_t_id$num]
    }
}
