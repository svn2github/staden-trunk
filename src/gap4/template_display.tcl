#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
##############################################################################
proc TempDisp_OkPressed { io f id infile l} {
    global CurContig

    set db [io_read_database $io]
    if {[keylget db Ntemplates] == 0} {
	bell
	verror ERR_WARN template_display "This database has no templates"
	destroy $f
	return
    }

    set contig_list ""
    if {[lorf_in_get $infile] == 4} {
	#single contig
	set contig_list [contig_id_gel $id]
	SetContigGlobals $io $contig_list
    } elseif {[lorf_in_get $infile] == 3 } {
	#all contigs
	set contig_list [CreateAllContigList $io]
    } else {
	#list or file
	set contig_list [lorf_get_list $infile]
	set contig_list [remove_contig_duplicates -io $io -contigs $contig_list]
    }

    if {$contig_list == ""} {
	raise $f
	return
    }

    # stop windows from hiding the plot
    destroy $f

    CreateTemplateDisplay $io $contig_list 
}

##############################################################################
#gaprc menu item
proc TempDispDialog { io } {
    global gap_defs InitConfig

    set f [keylget gap_defs TEMPLATE.WIN]
    set l [keylget gap_defs R_ENZ_LIST.WIN]
    xtoplevel $f -resizable 0
    wm title $f "Show templates"

    contig_id $f.id\
	    -io $io \
	    -range 0

    lorf_in $f.infile [keylget gap_defs TEMPLATE.INFILE] \
	    "{contig_id_configure $f.id -state disabled} \
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state normal}" -bd 2 -relief groove


    #default configuration
    InitTemplateConfigs 

    frame $f.cb
    checkbutton $f.cb.t -text "Templates" \
	    -variable InitConfig(template) -offvalue 0 -onvalue 1

    checkbutton $f.cb.r -text "Readings" \
	    -variable InitConfig(reading) -offvalue 0 -onvalue 1

    checkbutton $f.cb.mt -text "Ignore 'single' templates" \
	    -variable InitConfig(multi_template) -offvalue 0 -onvalue 1

    checkbutton $f.cb.rp -text "Show only read pairs"\
	    -variable InitConfig(read_pairs) -offvalue 0 -onvalue 1

    checkbutton $f.cb.srp -text "Show only spanning read pairs"\
	    -variable InitConfig(span_read_pairs) -offvalue 0 -onvalue 1

    checkbutton $f.cb.ccp -text "Calculate contig positions"\
	    -variable InitConfig(calc_contig_pos) -offvalue 0 -onvalue 1

    pack $f.cb.t $f.cb.r $f.cb.mt $f.cb.rp $f.cb.srp $f.cb.ccp -anchor w
    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "TempDisp_OkPressed $io $f $f.id $f.infile $l" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Template-Display}" \
	    -bd 2 \
	    -relief groove
    
    pack  $f.infile $f.id $f.cb $f.ok_cancel -side top -fill both
}

##############################################################################
proc InitTemplateConfigs { } {
    global InitConfig gap_defs

    set InitConfig(template) [keylget gap_defs TEMPLATE.CONFIG.TEMPLATE]
    set InitConfig(reading) [keylget gap_defs TEMPLATE.CONFIG.READING]
    set InitConfig(multi_template) [keylget gap_defs TEMPLATE.CONFIG.MULTI_TEMPLATE]
    set InitConfig(read_pairs) [keylget gap_defs TEMPLATE.CONFIG.READ_PAIR]
    set InitConfig(span_read_pairs) [keylget gap_defs TEMPLATE.CONFIG.SPAN_READ_PAIR]
    set InitConfig(calc_contig_pos) [keylget gap_defs TEMPLATE.CONFIG.CALC_CONTIG_POS]
}

##############################################################################
proc TagCheckList { parent path io id t_win r_win} {
    TagDialog TEMPLATE.TAGS $path \
	"DisplayReadingTags $parent $io $id $t_win $r_win" $parent
}

##############################################################################
proc DisplayReadingTags {f io id t_win r_win} {

    if {[winfo exists $t_win]} {
	$t_win delete tag
    } 
    if {[winfo exists $r_win]} {
	$r_win delete tag
    }
    display_reading_tags -io $io -id $id

}

##############################################################################
proc UpdateTemplateDisplay { io f t_win} {
    global config$f.template config$f.reading
    global gap_defs
    global $f.template_id

    if {![set config$f.template] && ![set config$f.reading]} {
	#delete template plot
	#DeleteTemplatePlot $io [set $f.template_id] $f $t_win
    } else {
	
	if {[grid info $f.t] == ""} {
	    #add new template plot
	    set height [keylget gap_defs TEMPLATE.PLOT_HEIGHT]
	    set borderwidth [keylget gap_defs TEMPLATE.BORDERWIDTH]
	    set new_height [expr [winfo height $f] + $height + \
		    2 * $borderwidth]

	    grid $f.t       -row 2 -column 0 -sticky nsew
	    grid $f.vscroll -row 2 -column 1 -sticky ns

	    add_template_plot -io $io -id [set $f.template_id]
	    wm geometry $f [winfo width $f]x$new_height
	    return
	}
    }
    update_template_display -io $io -id [set $f.template_id]
    
}

##############################################################################
proc RefreshContigOrder { io f } {
    global $f.template_id

    refresh_contig_order -io $io -id [set $f.template_id]
}
##############################################################################
#menu for template display
proc CreateTemplateMenu {io f t_win r_win hscroll scale_label brief} {
    global gap_defs read_only 
    global $f.contig_list
    global $f.template_id
    global template_menu

    #set the active tags for the template display to be the same as the
    #default values
    SetDefaultTags TEMPLATE.TAGS $f

    $f configure -menu $f.menubar
    menu $f.menubar
    create_menus $template_menu $f.menubar

    foreach i [set $f.contig_list] {
	set num [db_info get_contig_num $io $i]
	set path $f.menubar.[menu_path {View.Quality Plot}]

	$path add checkbutton \
		-label "contig $i" \
		-variable config$f.quality$num -offvalue 0 -onvalue 1 \
		-command "TemplateQualityPlot $io $f \[set $f.template_id\] \
		$hscroll $scale_label $i $num"
    }

    foreach i [set $f.contig_list] {
	set num [db_info get_contig_num $io $i]
	set path $f.menubar.[menu_path {View.Restriction Enzyme Plot}]

	$path add checkbutton \
		-label "contig $i" \
		-variable config$f.renz$num -offvalue 0 -onvalue 1 \
		-command "TemplateREnzPlot $io $f \[set $f.template_id\] \
		$f.brief $f.buttons.dist $i $num"
    }

    if {$read_only} {
        menu_state_set template_menu -16 $f.menubar
    } else {
	menu_state_set template_menu 12 $f.menubar
    }
}

##############################################################################
#update titlebar if changed contig
proc SetTitleBar { io f CurContig } {

    #set title bar of window to display the current contig identifier
    wm title $f "Template display: $CurContig #[db_info get_read_num $io $CurContig]"
}

##############################################################################
#initialise configuration settings for the template display
proc SetTemplateConfig {f} {
    global InitConfig
    global config$f.template config$f.reading 
    global config$f.multi_template config$f.read_pairs config$f.span_read_pairs
    global config$f.calc_contig_pos
    global config$f.ruler config$f.renz config$f.ticks 

    set config$f.reading $InitConfig(reading)
    set config$f.template $InitConfig(template)
    set config$f.multi_template $InitConfig(multi_template)
    set config$f.read_pairs $InitConfig(read_pairs)
    set config$f.span_read_pairs $InitConfig(span_read_pairs)
    set config$f.calc_contig_pos $InitConfig(calc_contig_pos)
    set config$f.ruler 1
    set config$f.ticks 1
    set config$f.renz 0

}

##############################################################################
proc UpdateRuler {io f id hscroll label } {
    global $f.contig_list $f.template_id
    global gap_defs
    global config$f.ruler
    global restoreCmd

    set t_win $f[keylget gap_defs TEMPLATE.TEMPLATE.WIN]
    set r_win $f[keylget gap_defs TEMPLATE.RULER.WIN]

    if { [set config$f.ruler] } {
	#plot ruler

	#HACK to work round the horrible problem in zooming, turning ruler
	#off, scrolling and turning ruler on again.
	set scroll [lindex [$hscroll get] 0]
	if { ![winfo exists $r_win] } {
	    set zoom_cmdr [list "gap_zoom $io \[set $f.template_id\] x"]
	    CreateTemplateRuler $io $f $r_win $hscroll [set $f.contig_list] \
		    [winfo width $t_win] $zoom_cmdr
	    SetRulerBindings $io $f $r_win $t_win $label
	}
	display_ruler -io $io -id $id -win_ruler $r_win
	DisplayReadingTags $f $io $id $t_win $r_win
	uplevel #0 eval [$hscroll cget -command] moveto $scroll
    } else {
	#delete ruler

	set new_height [expr [winfo height $f] - [winfo height $f.r]]
	
	delete_window -io $io -id $id -window $r_win
	$r_win delete all
	destroy $r_win 
	destroy $f.r
	wm geometry $f [winfo width $f]x$new_height
	if {[info exists restoreCmd($f,contig)]} {
	    unset restoreCmd($f,contig)
	}
    }
    
}

##############################################################################
proc DisplayRulerTicks {io f id hscroll } {
    global config$f.ticks
    global gap_defs

    set r_win $f[keylget gap_defs TEMPLATE.RULER.WIN]

    if {[set config$f.ticks]} {
	display_ruler_ticks -io $io -id $id -ticks 1
    } else {
	if {[winfo exists $r_win]} {
	    $r_win delete tick
	    RulerWindowSize 0 $f $r_win
	}
    }
}

proc RulerWindowSize {disp_ticks f r_win } {
    global gap_defs

    if {$disp_ticks} {
	set bbox [$r_win bbox tick]
    } else {
	set bbox [$r_win bbox contig]
    }

    if {[lindex $bbox 3] <= [winfo height $r_win]} {
	
	if {[lindex $bbox 3] > [keylget gap_defs TEMPLATE.RULER.PLOT_HEIGHT]} {
	    set height [lindex $bbox 3]
	} else {
	    set height [keylget gap_defs TEMPLATE.RULER.PLOT_HEIGHT]
	}
	set new_height [expr abs([winfo height $r_win] - $height)]

	$r_win config -height $height
	set new_height [expr [winfo height $f] - $new_height]
	wm geometry $f [winfo width $f]x$new_height
	#update idletasks
    }

    if {[lindex $bbox 3] > [winfo height $r_win]} {

	set height [lindex $bbox 3]
	set new_height [expr abs([winfo height $r_win] - $height)]
	$r_win config -height $height
	set new_height [expr [winfo height $f] + $new_height]
	wm geometry $f [winfo width $f]x$new_height
	#update idletasks
    }
}

##############################################################################
proc CreateTemplateRuler {io f r_win hscroll contig_list width zoom_cmd} {
    global gap_defs

    set NGContig(min_x) 100
    set scroll x
    set borderwidth [keylget gap_defs TEMPLATE.BORDERWIDTH]
    set height [keylget gap_defs TEMPLATE.RULER.PLOT_HEIGHT]

    #allow -height and -width to have affect
    wm geometry $f {}
    
    frame $f.r -bd $borderwidth -relief groove
    canvasbox $r_win -width $width \
	    -height $height \
	    -xscrollcommand "$hscroll set" \
	    -highlightthickness 0 \
	    -bd 0 -zoom_command $zoom_cmd
    
    #set toplevel geometry when resized window
    set new_height [expr [winfo height $f] + [winfo reqheight $r_win] + \
	    2 * $borderwidth]
    
    #pack $f.r -in $f.mid.left -fill x -side top
    grid rowconfig $f 3
    grid $f.r -row 3 -column 0 -sticky nsew
    pack $r_win -in $f.r -padx 5 -fill both -expand yes
    update idletasks

    #set toplevel geometry when resized window
    #check that things have been packed!
    if {[winfo height $f] != 1} {
	wm geometry $f [winfo width $f]x$new_height
    }
}

##############################################################################
proc GetItemInfo {io plot nearest} {

    foreach tag [$plot gettags $nearest] {
	if {[string compare [string range $tag 0 1] r_] == 0} {
	    set r_num [string trim $tag r_]
	    return "Reading: [r_name $io $r_num] (#$r_num)   \
		    Length: [r_length $io $r_num]"
	} elseif {[string compare [string range $tag 0 2] te_] == 0} {
	    set t_num [string trim $tag te_]
	    return "Template: [t_name $io $t_num]"
	} elseif {[string compare [string range $tag 0 3] num_] == 0} {
	    set c_num [string trim $tag num_]
	    set r_num [db_info get_read_num $io =$c_num]
	    set r [io_read_reading $io $r_num]
	    set dir [lindex {+ -} [keylget r sense]]
	    return "Contig: [left_gel $io $c_num]\
		    ($dir#$r_num) \
		    Length: [c_length $io $c_num]\
                    Num readings: [num_readings_in_contig $io $c_num]"
	}
    }
}

##############################################################################
#highlight readings/templates as move mouse over
proc HighlightTemplateDisplay {io f plot t_win r_win brief} {
    global restoreCmd
    global gap_defs
    global $f.prev_item

    if {![info exists $f.prev_item]} {
	set $f.prev_item ""
    }
    set nearest [$plot find withtag current]

    #only do this code if the nearest item is different from the previous one
    if {$nearest == [set $f.prev_item]} {
	return 
    }
    #unhighlight object
    if {[info exists restoreCmd($f,reading)]} {
	eval $restoreCmd($f,reading)
	
	#keep readings always raised
	set item [lindex $restoreCmd($f,reading) 2]
	if {[string compare [lindex [GetItemInfo $io $plot $item] 0] \
		Reading:] == 0} {
	    $plot raise $item

	}
    
    }
    if {[info exists restoreCmd($f,contig)]} {
	    eval $restoreCmd($f,contig)
    }

    if {$nearest != 0} {

	if {$plot == $t_win} {
	    InitialSettings $f $plot $nearest reading
	    set tags [$t_win gettags $nearest]
	    $plot itemconfig $nearest \
		    -fill [keylget gap_defs TEMPLATE.HIGHLIGHT_COLOUR]
	    $plot raise $nearest
	}
	if {[winfo exists $r_win]} {
	    foreach tag [$plot gettags $nearest] {
		if {[string compare [string range $tag 0 3] num_] == 0} {
		    set c_num [string trim $tag num_]
		    InitialSettings $f $r_win $tag contig
		    $r_win itemconfig num_$c_num \
			    -fill [keylget gap_defs TEMPLATE.HIGHLIGHT_COLOUR]

		}
	    }
	}
	$brief configure -text [GetItemInfo $io $plot $nearest]
    }

    #set previous item
    set $f.prev_item $nearest
}

#need to raise readings when get an leave event - which includes bringing up
#a popup menu!
proc UnHighlightTemplateDisplay {f t_win} {
    global $f.prev_item

    if {[info exists $f.prev_item]} {
	set prev [set $f.prev_item]
    } else {
	UnHighlightItems $f
	return
    }
    set tags [$t_win gettags $prev]
    
    UnHighlightItems $f
    if {[lsearch reading $tags] != 0} {
	$t_win raise $prev
    }
}

##############################################################################
#print contig selector tag info
proc PrintTagDetails { canvas io t_num } {

    set a [io_read_annotation $io $t_num]
    
    set str ""
    append str  "position [keylget a position] \n"
    append str "length [keylget a length] \n"
    append str "type [keylget a type] \n"
    if {[keylget a annotation] != 0} {
	append str "comment [io_read_text $io [keylget a annotation]] \n"
    }

    vfuncgroup 2 "Template display"
    vmessage $str
    messagebox $canvas $str
}

##############################################################################
proc PopUpTagMenu { io canvas obj X Y} {

    if {[winfo exists $canvas.m]} {destroy $canvas.m}
    
    create_popup $canvas.m "Tag Commands"
    $canvas.m add command -label information \
	    -command "destroy $canvas.m; \
	    PrintTagDetails $canvas $io $obj"
    tk_popup $canvas.m [expr $X-20] [expr $Y-10]
}

##############################################################################
#print contig selector tag info
proc PrintTemplateDetails {canvas io template_id t_num} {

    set a [io_read_template $io $t_num]

    set v [io_read_vector $io [keylget a vector]]
    set vector [io_read_text $io [keylget v name]] 
    set c [io_read_clone $io [keylget a clone]]
    set clone [io_read_text $io [keylget c name]] 

    set str ""
    append str "Template name\t\t[t_name $io $t_num] \n"
    append str "number\t\t\t$t_num \n"
    append str "strands\t\t\t[keylget a strands] \n"
    append str "vector\t\t\t$vector ([keylget a vector]) \n"
    append str "clone\t\t\t$clone ([keylget a clone]) \n"
    append str "min insert length\t[keylget a insert_length_min] \n"
    append str "max insert length\t[keylget a insert_length_max] \n"
    append str [print_template_readings -io $io -id $template_id -tnum $t_num]
    
    vfuncgroup 2 "Template display"
    vmessage $str
    messagebox $canvas $str
    #tk_messageBox -message $str -type ok -icon info
}

##############################################################################
#print contig selector tag info
proc PrintReadingDetails {canvas io r_num} {

    set a [io_read_reading $io $r_num]
    set s [keylget a sense]
    if {$s == 1} {
	set sense "-"
    } else {
	set sense "+"
    }

    set str ""
    append str "Reading name\t[r_name $io $r_num] \n"
    append str "number\t\t$r_num \n"
    append str "left neighbour\t[keylget a left] \n"
    append str "right neighbour\t[keylget a right] \n"
    append str "position\t[keylget a position] \n"
    append str "total length\t[keylget a length] \n"
    append str "sequence length\t[keylget a sequence_length] \n"
    append str "used start\t[keylget a start] \n"
    append str "used end\t[keylget a end] \n"
    append str "template\t[t_name $io [keylget a template]] ([keylget a template]) \n"
    append str "sense\t\t$sense \n"

    #if primer is custom but strand is 1, then primer is custom reverse (4)

    if {[keylget a strand] == 1 && [keylget a primer] == 3} {
	append str "primer\t\t[expr [keylget a primer] + 1] \n"
    } else {
	append str "primer\t\t[keylget a primer] \n"
    }

    append str "\n"

    vfuncgroup 2 "Template display"
    vmessage $str
    messagebox $canvas $str 
}

##############################################################################
proc PrintContigDetails {canvas io c_num args} {
    set str ""

    set c [io_read_contig $io $c_num]

    append str "Contig\t [left_gel $io $c_num]\n"
    append str "left\t [keylget c left]\n"
    append str "right\t [keylget c right]\n"
    append str "length\t [keylget c length]\n"

    vfuncgroup 2 "Template display"
    vmessage $str
    messagebox $canvas $str
}

##############################################################################
proc PrintTemplateDisplayDetails {io template_id canvas obj} {

    foreach item [$canvas gettags $obj] {
	if {[string compare [string range $item 0 1] r_] == 0} {
	    set r_num [string trim $item r_]
	    PrintReadingDetails $canvas $io $r_num
	} elseif {[string compare [string range $item 0 2] te_] == 0} {
	    set te_num [string trim $item te_]
	    PrintTemplateDetails $canvas $io $template_id $te_num
	} elseif {[string compare [string range $item 0 1] t_] == 0} {
	    set t_num [string trim $item t_]
	    PrintTagDetails $canvas $io $t_num
	} elseif {[string compare [string range $item 0 2] hl_] == 0} {
	    set c_num [string trim $item hl_]
	    PrintContigDetails $canvas $io $c_num
	}
    }
}

##############################################################################
proc td_config {io f t_win c_num} {
    global $f.template_id

    if {[set t [xtoplevel $f.lbox$c_num -resizable 0]] == ""} {
	return
    }
    wm title $f.lbox$c_num "Highlight templates"

    set width [$t_win itemcget c_$c_num -width]
    
    # Line width
    frame $t.lw -bd 3 -relief raised
    scale $t.lw.scale \
	    -label "Line width" \
	    -from 1 \
	    -to 10 \
	    -orient horiz \
	    -command "$t_win itemconfigure c_$c_num -width "
    $t.lw.scale set $width

    ###########################################################################
    #OK and Cancel buttons
    frame $t.ok_cancel -bd 2 -relief groove
    button $t.ok_cancel.ok -text OK \
	    -command "$t_win itemconfigure c_$c_num -width \[$t.lw.scale get\]; destroy $t"
    button $t.ok_cancel.cancel -text Cancel\
	    -command "destroy $t; $t_win itemconfigure c_$c_num -width $width"

    pack $t.ok_cancel.ok $t.ok_cancel.cancel -side left
    
    pack $t.lw -side top -fill both
    pack $t.lw.scale -side top -fill x
    pack $t.ok_cancel -side bottom -fill x
}

##############################################################################
proc HighlightContigTemplates {io f r_win t_win obj} {
    global config$f.template

    if {![set config$f.template]} {
	return
    }
    foreach item [$r_win gettags $obj] {
	if {[string compare [string range $item 0 3] num_] == 0} {
	    set c_num [string trim $item num_]
	}
    }
    if {![winfo exists $t_win]} {
	return
    }
    td_config $io $f $t_win $c_num
}

##############################################################################
proc PopUpTemplateMenu { io template_id canvas x y X Y} {

    if {[winfo exists $canvas.m]} {destroy $canvas.m}
    
    set nearest [$canvas find withtag current]
    set obj [itemNearest $canvas $x $y]

    foreach item [$canvas gettags $obj] {
	if {[string compare $item template]==0} {
	    set type Template
	}
	if {[string compare $item reading]==0} {
	    set type Reading
	}
	if {[string compare $item tag]==0} {
	    set type Tag
	}
    }
    if {![info exists type]} {
	return 
    }
    create_popup $canvas.m "$type Commands"

    $canvas.m add command -label information \
	    -command "destroy $canvas.m; \
	    PrintTemplateDisplayDetails $io $template_id $canvas $obj"
    if {$type == "Reading"} {
        $canvas.m add command -label "List notes" \
	    -command "destroy $canvas.m; \
		      ListReadNotes $io $template_id $canvas $obj"
    }
    
    tk_popup $canvas.m [expr $X-20] [expr $Y-10]
}

proc ListReadNotes {io read_id canvas obj} {
    foreach item [$canvas gettags $obj] {
	if {[string compare [string range $item 0 1] r_] == 0} {
	    set r_num [string trim $item r_]
	    NoteSelector $io reading "#$r_num"
	}
    }
}

##############################################################################
proc PopUpTemplateContigMenu {io f r_win t_win x y X Y } {
    global read_only
    global $f.template_id
    
    if {[winfo exists $r_win.m]} {destroy $r_win.m}
    
    set obj [itemNearest $r_win $x $y]

    #determine type of object selected
    #if its a tag, call PopUpTagMenu and return
    foreach item [$r_win gettags $obj] {
	if {[string compare $item tag] == 0} {
	    PopUpTagMenu $io $r_win $obj $X $Y
	    return
	}
    }
    
    #otherwise, popup the contig menu
    create_popup $r_win.m "Contig Commands"

    # Work out position clicked within the contig
    set current [$r_win find withtag current]
    set cnum [GetContigNum $r_win $current]
    set c [io_read_contig $io $cnum]
    set clen [keylget c length]
    foreach {x1 y1 x2 y2} [$r_win coords $current] {}
    set x [$r_win canvasx $x]
    set cpos [expr {int(((double($x-$x1))/($x2-$x1)) * $clen)}]
    if {$cpos < 1} {set cpos 1}
    if {$cpos > $clen} {set cpos $clen}

    $r_win.m add command -label Information \
	    -command "destroy $r_win.m; \
	    PrintTemplateDisplayDetails $io [set $f.template_id] $r_win $obj"
    $r_win.m add command -label "Edit contig" \
	    -command "popup_cs_contig_1 $io $r_win $current $cpos"
    if {!$read_only} {
	$r_win.m add command -label "Complement contig" \
		-command "popup_cs_contig_3 $io $r_win [$r_win find withtag current]"
    }
    $r_win.m add command -label "Highlight templates" \
	    -command "HighlightContigTemplates $io $f $r_win $t_win [$r_win find withtag current]"

    $r_win.m add command -label "List notes" \
	    -command "popup_cs_contig_cnotes $io $r_win [$r_win find withtag current]"

    tk_popup $r_win.m [expr $X-20] [expr $Y-10]
    
}

##############################################################################
#display the template crosshair
proc AddTemplateCrossHair {io id f canvas x} {
    global $f.template_id
    global $f.cursor
    
    if {[set $f.cursor]} {
	
	draw_canvas_cursor_x -io $io -id [set $f.template_id] \
		-x [$canvas canvasx $x]
    } else {
	delete_canvas_cursor -io $io -id [set $f.template_id]
    }
}

##############################################################################
#called from C
proc DrawTemplateCursor {f canvas cx} {
    global gap_defs

    set height [$canvas canvasy [winfo height $canvas]]
    if {[$canvas find withtag cursor] == ""} {
	$canvas create line $cx [$canvas canvasy 0] $cx $height -tag cursor\
		-fill [keylget gap_defs TEMPLATE.CURSOR_COLOUR]
	$canvas lower cursor
    } else {
	$canvas coords cursor $cx [$canvas canvasy 0] $cx $height
    }
}

##############################################################################
#Creates and/or moves an editor cursor from within the template display
proc template_cursor_editor {io win id x} {
    set cnum [template_contig -io $io -id $id -x [$win canvasx $x]]
    canvas_cursor_editor $io $cnum $id $x $win
}

##############################################################################
#bindings specific to item display
proc SetTemplateBindings { io f t_win r_win  label} {
    global $f.template_id
    global gap_defs
    
    #$canvas bind S <Any-Motion> "HighlightItems $io $f %W $label"
    $t_win bind S <Any-Motion> "HighlightTemplateDisplay $io $f %W $t_win $r_win $label"
    $t_win bind S <Shift-Motion> {;}
    bind $t_win <Any-Leave> "+UnHighlightTemplateDisplay $f $t_win"
	
    bind $t_win <<menu>> "PopUpTemplateMenu $io [set $f.template_id] %W %x %y %X %Y"

    bind $t_win <Any-Leave> "+delete_canvas_cursor -io $io -id [set $f.template_id]"
    bind $t_win <Any-Motion> "+AddTemplateCrossHair $io [set $f.template_id] $f  $t_win %x"

    # Double button 2 to move or create an editor
    bind $t_win <<move-create>> "
	template_cursor_editor $io $t_win [set $f.template_id] %x
    "

    # Double button 1 to move or create an editor
    bind $t_win <<use>> "
	template_cursor_editor $io $t_win [set $f.template_id] %x
    "

    #reading selection which adds to a reading list and highlights the
    #contig editor if it is running
    $t_win bind reading <<select>> "SelectSingleRead $io $t_win"
    bind $t_win <<select>> "itemMark $t_win %x %y"
    bind $t_win <<select-drag>> "itemStroke $t_win %x %y"
    bind $t_win <<select-release>> "SelectReadings $f $t_win $io"
    
}

##############################################################################
proc InputSelectReadingList { io } {
    global gap_defs

    set name [keylget gap_defs READING_LIST.WIN]
    if {[xtoplevel $name -resizable 0] == ""} {
	return
    }
    wm title $name "Highlight reading list"

    lorf_in $name.infile [keylget gap_defs READING_LIST.INFILE] \
	"" -bd 2 -relief groove
 
    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $name.ok_cancel \
	    -ok_command  "CreateReadingList $io \[lorf_get_list $name.infile\]; destroy $name"\
	    -cancel_command "destroy $name" \
	    -help_command "show_help gap4 {Template-Templates-Operations}" \
	    -bd 2 \
	    -relief groove
    
    pack $name.infile -fill both
    pack $name.ok_cancel -fill both

}

##############################################################################
proc CreateReadingList { io r_list } {

    SetReadingList $r_list
    foreach r_name $r_list {
	reg_notify_highlight -io $io -reading $r_name -highlight 1
    }

}

##############################################################################
proc SelectSingleRead {io canvas } {
    set item [$canvas gettags [$canvas find withtag current]]
    foreach tag $item {
	if {[string compare [string range $tag 0 1] r_] == 0} {
	    set r_list [r_name $io [string trim $tag r_]]
	}
    }
    set h_list [ToggleReadingList $r_list]
    foreach r $h_list {
	reg_notify_highlight -io $io -reading [lindex $r 0] \
		-highlight [lindex $r 1]
    }
}

##############################################################################
#when drag out region in template display, want to toggle between highlighted
#and non-highlighted readings
#update "Readings" list 
#notify changes to all other template displays and contig editors
proc SelectReadings {f canvas io } {
    
    set r_list ""
    set list [itemsUnderArea $f $canvas reading]

    foreach i $list {
	set item [$canvas gettags $i]
	foreach tag $item {
	    if {[string compare [string range $tag 0 1] r_] == 0} {
		set r_num [string trim $tag r_]
		lappend r_list [r_name $io $r_num]
	    }
	}
    }
    set h_list [ToggleReadingList $r_list]
    foreach r $h_list {
	reg_notify_highlight -io $io -reading [lindex $r 0] \
		-highlight [lindex $r 1]
    }

}

##############################################################################
#highlights the readings in the current "Readings" list
#call when create new template display
proc SelectReadingList { io } {
    
    set r_list [GetReadingList]
    if {$r_list == ""} {
	return
    }
    foreach r_name $r_list {
	reg_notify_highlight -io $io -reading $r_name -highlight 1
    }
}

proc templateItemBind { f r_win x y} {
    global $f.sel_contigs gap_defs

    set cnum [GetContigNum $r_win current]

    $r_win itemconfig hl_$cnum -width [keylget gap_defs TEMPLATE.LINE_BOLD]
   
    itemMark $r_win $x $y
    lappend $f.sel_contigs "$cnum"
}

proc templateItemMove {f r_win x y} {
    global $r_win.areaX1 $r_win.areaY1
    global $f.sel_contigs

    set min_x 1
    set max_x [winfo width $r_win]
    set min_y 1
    set max_y [winfo height $r_win]

    if {($x > $min_x) && ($x < $max_x) && ($y > $min_y) && ($y < $max_y)} {
	set x [$r_win canvasx $x]
	set dx [expr $x - [set $r_win.areaX1]]

	foreach cnum [set $f.sel_contigs] {
	    $r_win move num_$cnum $dx 0
	}
	set $r_win.areaX1 $x
	set $r_win.areaY1 $y
    }
}


#list of contigs in template in left to right order - used in quality plot
#and restriction enzyme plots
proc UpdateTemplateContigList {io f c_list} {
    global $f.contig_list $f.template_id $f.contig_order

    if {[llength [set $f.contig_list]] > 1} {
	set $f.contig_order 1
    }

    set hscroll $f.hscroll
    set scale_label $f.scale_label

    #complementing sends over a new contig list, moving a contig has already
    #set the $f.contig_list 
    if {$c_list != ""} {
	set $f.contig_list $c_list
    }
    
    set path $f.menubar.[menu_path {View.Quality Plot}]
    
    $path delete 1 [llength [set $f.contig_list]]
    foreach i [set $f.contig_list] {
	set num [db_info get_contig_num $io $i]
	$path add checkbutton \
		-label "contig $i" \
		-variable config$f.quality$num -offvalue 0 -onvalue 1 \
		-command "TemplateQualityPlot $io $f \[set $f.template_id\] \
		$hscroll $scale_label $i $num"
    }

    set path $f.menubar.[menu_path {View.Restriction Enzyme Plot}]

    $path delete 1 [llength [set $f.contig_list]]
    foreach i [set $f.contig_list] {
	set num [db_info get_contig_num $io $i]
	$path add checkbutton \
		-label "contig $i" \
		-variable config$f.renz$num -offvalue 0 -onvalue 1 \
		-command "TemplateREnzPlot $io $f \[set $f.template_id\] \
		$f.brief $f.buttons.dist $i $num"
    }
}

proc templateItemDrop {io f r_win x} {
    global $f.template_id $f.contig_list $f.sel_contigs gap_defs

    set x [$r_win canvasx $x]
    
    set $f.contig_list [update_template_contig_order -io $io -id [set $f.template_id] \
	    -contigs "=[set $f.sel_contigs]" -x $x]

    #update the list of templates to maintain the same left to right order
    UpdateTemplateContigList $io $f ""

    foreach cnum $f.sel_contigs {
	$r_win itemconfig hl_$cnum -width [keylget gap_defs TEMPLATE.LINE_WIDTH]
    }
    set $f.sel_contigs ""
}

##############################################################################
proc RulerEditorCursor {io f r_win x y} {
    global $f.template_id

    set obj [$r_win find closest $x $y] 
    foreach tag [$r_win gettags $obj] {
	if {[string compare [string range $tag 0 3] num_] == 0} {
	    set contig [string trim $tag num_]
	    canvas_cursor_editor $io $contig [set $f.template_id] $x $r_win
	    return
	}
    }
}

##############################################################################
proc SetRulerBindings { io f r_win t_win label} {
    global $f.template_id $r_win.Move

    $r_win bind S <Any-Motion> "HighlightTemplateDisplay $io $f %W $t_win $r_win $label"
    $r_win bind S <Shift-Motion> {;}
    bind $r_win <Any-Leave> "+UnHighlightItems $f"
    
    bind $r_win <Any-Leave> "+delete_canvas_cursor -io $io -id [set $f.template_id]"
    bind $r_win <Any-Motion> "+AddTemplateCrossHair $io [set $f.template_id] $f  $r_win %x"
    
    $r_win bind contig <<menu>> "PopUpTemplateContigMenu $io $f %W $t_win %x %y %X %Y"
    $r_win bind tag <<menu>> "+PopUpTemplateMenu $io [set $f.template_id] %W %x %y %X %Y"

    # Double button 1 and 2 to move or create an editor
    bind $r_win <<move-create>> "RulerEditorCursor $io $f $r_win %x %y"
    bind $r_win <<use>> "RulerEditorCursor $io $f $r_win %x %y"

    # Disable panning
    bind $r_win <<move-drag>> {break}

    set $r_win.Move 0
    #button 2 for moving - select and move
    $r_win bind contig <<move>> \
	    "set $r_win.Move 1; templateItemBind $f %W %x %y"
	
    #moving a contig
    $r_win bind contig <<move-drag>> \
	    "if {\[set $r_win.Move\] == 1} {
	        templateItemMove $f %W %x %y
	    }"

    $r_win bind contig <<move-release>> \
	    "if {\[set $r_win.Move\] == 1} {
		templateItemDrop $io $f %W %x
		set $r_win.Move 0
	    }"
}

##############################################################################
#display the tag menu if readings are displayed
#HACK - not work yet
proc DisplayTagMenu {f menubar} {
    global config$f.reading
    
    #if reading are displayed allow choosing of tags
    # if reading are not displayed, disallow tags
    if {[set config$f.reading]} {
	$menubar.view.opts entryconfigure "select tags" -state normal
    } else {
	$menubar.view.opts entryconfigure "select tags" -state disabled
    }
}

##############################################################################
proc next_template_display { } {
    global next_t_display
    
    if {![info exists next_t_display]} {
	set next_t_display 0
	return $next_t_display
    }
    incr next_t_display
    return $next_t_display
    
}


##############################################################################
#display templates and readings
proc CreateTemplateDisplay {io contig_list} {
    global gap_defs
    global CurContig
    global restoreCmd
    global InitConfig
    
    set db [io_read_database $io]
    if {[keylget db Ntemplates] == 0} {
	bell
	verror ERR_WARN template_display "This database has no templates"
	return
    }

    set st [keylget gap_defs TEMPLATE.WIN]
    set num_display [next_template_display]

    set f $st$num_display
    set t_win $f[keylget gap_defs TEMPLATE.TEMPLATE.WIN]
    set r_win $f[keylget gap_defs TEMPLATE.RULER.WIN]
    
    #set the contig_list for each template display
    global $f.contig_list
    global $f.template_id
    
    set $f.contig_list $contig_list
    
    xtoplevel $f
    fix_maxsize $f
    
    #needed when TemplateDisplay is invoked using readpairs on contig selector
    set contig [lindex $contig_list 0] 

    #set title bar of window to display the current contig identifier
    SetTitleBar $io $f $contig

    #highlighting
    global $f.prev_item
    set $f.prev_item 0
    set restoreCmd($f,item_num) 0
    
    #set display both x and y scrollbars
    set scroll b
    set borderwidth [keylget gap_defs TEMPLATE.BORDERWIDTH]
    set width [keylget gap_defs TEMPLATE.PLOT_WIDTH]
    set height [keylget gap_defs TEMPLATE.PLOT_HEIGHT]

    #want to set config if not called routine from dialogue box
    if {![info exists InitConfig(template)]} {
	InitTemplateConfigs
    }
    SetTemplateConfig $f

    ##########################################################################
    #template display
    frame $f.t -bd $borderwidth -relief groove
    set zoom_cmdt [list "gap_zoom $io \[set $f.template_id\] $scroll"]
    canvasbox $t_win -width $width -height $height \
	    -bd 0 -highlightthickness 0\
	    -xscrollcommand "$f.hscroll set" \
	    -yscrollcommand "$f.vscroll set" \
	    -closeenough 4 -zoom_command $zoom_cmdt
	    
    scrollbar $f.hscroll -orient horizontal -relief sunken -command \
	    "gc_scroll_x $io \[set $f.template_id\]"
    scrollbar $f.vscroll -relief sunken -command "gc_scroll_y $io \[set $f.template_id\]"

    ##########################################################################
    # Main Menu Bar
    #frame $f.menubar -relief raised -borderwidth 2
    CreateTemplateMenu $io $f $t_win $r_win $f.hscroll $f.scale_label $f.brief
    #DisplayTagMenu $f $f.menubar
    
    ##########################################################################
    #button bar
    frame $f.buttons

    #zoom back button. use ruler window since this can no longer be removed
    button $f.buttons.back -text "zoom out" -command "ZoomBackCanvas $io \[set $f.template_id\]"
    button $f.buttons.zoomin10 -text "+10%" \
	-command "ZoomInCanvas $r_win 0.05" 
    button $f.buttons.zoomin50 -text "+50%" \
	-command "ZoomInCanvas $r_win 0.1666" 
   
    #cursor checkbutton
    global $f.cursor
    checkbutton $f.buttons.cursor -text crosshairs -variable $f.cursor

    #cursor local and total position labels
    set cursor_t [keylget gap_defs TEMPLATE.CURSOR1]
    set cursor_l [keylget gap_defs TEMPLATE.CURSOR2]
    label $f$cursor_t -bd 2 -relief sunken -width 6
    label $f$cursor_l -bd 2 -relief sunken -width 6

    #restriction enzyme distance label
    label $f.buttons.dist -relief sunken -bd 2 -width 6

    pack $f.buttons.zoomin10 $f.buttons.zoomin50 $f.buttons.back \
	 $f.buttons.cursor -expand yes -side left
    pack $f$cursor_l $f$cursor_t -in $f.buttons -side left -expand yes 
    pack $f.buttons.dist -side left -expand yes 
    
    ##########################################################################
    label $f.brief
    
    grid columnconfig $f 0 -weight 1
    grid rowconfig $f 2 -weight 2

    #grid $f.menubar -row 0 -column 0 -sticky ew -columnspan 2
    grid $f.buttons -row 1 -column 0 -sticky ew -columnspan 2
    grid $f.t       -row 2 -column 0 -sticky nsew
    pack $t_win -in $f.t -padx 5 -pady 5 -fill both -expand yes
    grid $f.vscroll -row 2 -column 1 -sticky ns
    grid $f.hscroll -row 6 -column 0 -sticky ew
    grid $f.brief   -row 7 -column 0 -sticky ew

    #must ensure that the packing is complete before calling template_reg
    #which interrogates the canvas width and height
    tkwait visibility $t_win

    #set zoom command for ruler
    set zoom_cmdr [list "gap_zoom $io \[set $f.template_id\] x"]
    
    CreateTemplateRuler $io $f $r_win $f.hscroll $contig_list $width $zoom_cmdr

    #register template display and do first display
    set $f.template_id [display_templates -io $io \
			    -contigs $contig_list \
			    -frame $f \
			    -win_template $t_win \
			    -win_ruler $r_win ]

    #if user tries to destroy window 
    wm protocol $f WM_DELETE_WINDOW "TemplateStartShutdown $io $f \[set $f.template_id\]"

    #bind the configure actions to the toplevel
    bind $f <Any-Configure> "
    if {\[winfo toplevel %W\] == \"%W\"} {
	update idletasks
	resize_canvas -io $io -id [set $f.template_id]
    }
    "
    #SetCanvasBindings $io [set $f.template_id] $t_win $scroll 
    SetCanvasBindings $t_win $zoom_cmdt
    SetTemplateBindings $io $f $t_win $r_win $f.brief
    
    #SetCanvasBindings $io [set $f.template_id] $r_win x
    SetCanvasBindings $r_win $zoom_cmdr
    SetRulerBindings $io $f $r_win $t_win $f.brief
    SelectReadingList $io

    # Quality plot
    if {[llength $contig_list] == 1 && \
	[keylget gap_defs TEMPLATE.QUALITY.AT_STARTUP] == 1} {
	set i [lindex $contig_list 0]
	set num [db_info get_contig_num $io $i]
       	global config$f.quality$num 1
	set config$f.quality$num 1
	TemplateQualityPlot $io $f [set $f.template_id] \
		$f.hscroll $f.scale_label $i $num
    }
    
    if {[info exists InitConfig(template)]} {
	unset InitConfig(template)
    }

    #HACK - don't know why this doesn't work
    #parray InitConfig
    #puts "INIT [info exists InitConfig]"
    #unset InitConfig
    #parray InitConfig

}

##############################################################################
proc TemplateStartShutdown { io f template_id } {
    global $f.contig_order read_only
    global $f.ShuttingDown

    if {[info exists $f.ShuttingDown]} {
	return
    }
    set $f.ShuttingDown 1
    
    #check if want to save updated contig order
    if {[info exists $f.contig_order]} {
	if {$read_only} {
	    tk_messageBox -icon warning -type ok \
		    -title "Contig order has changed" \
		    -message "Database is read only. Unable to save contig order" \
		    -parent $f
	} else {
	    set ret [tk_messageBox -icon warning -type yesnocancel \
		    -title "Contig order has changed" \
		    -message "Do you wish to save the contig order?" \
		    -parent $f]
	    if {$ret == "yes"} {
		RefreshContigOrder $io $f
	    }
	}
	unset $f.contig_order
    }

    if {[info exists template_id]} {
	clear_template -io $io -id $template_id
	result_quit -io $io -id $template_id
    }
}

##############################################################################
#HACK - to tidy up
#delete just the template and reading plot
proc DeleteTemplatePlot { io id f t_win} {
    global restoreCmd

    $t_win delete all

    if {[info exists restoreCmd($f,reading)]} {
	unset restoreCmd($f,reading)
    }

    set new_height [expr [winfo height $f] - [winfo height $f.t]]
    #delete_window -io $io -id $id -window $t_win
    grid forget $f.t
    grid forget $f.vscroll

    wm geometry $f [winfo width $f]x$new_height
    update idletasks
}

##############################################################################
#HACK - to tidy up
#delete entire template display and unset global variables
proc DeleteTemplateDisplay { f t_win id  } {
    global config$f.template config$f.reading config$f.multi_template
    global config$f.read_pairs config$f.span_read_pairs
    global $f.contig_list
    global restoreCmd
    
    if {[info exists restoreCmd($f,reading)]} {
	unset restoreCmd($f,reading)
    }
    if {[info exists restoreCmd($f,contig)]} {
	unset restoreCmd($f,contig)
    }

    #unset all global variables
    unset config$f.template       
    unset config$f.reading        
    unset config$f.multi_template 
    unset config$f.read_pairs    
    unset config$f.span_read_pairs    

    unset $f.contig_list

    destroy $f
}

##############################################################################
proc FindReadingYCoord { canvas r_num} {
    
    set item [$canvas find withtag r_$r_num]

}

##############################################################################
proc SetTSizeTolerance {io f t} {
    global template_size_tolerance

    set w $f.tss
    if {[xtoplevel $w -resizable 0] == ""} {
	return
    }
    wm title $w "Set Template Size Tolerance"

    xentry $w.size \
	-label "Limits scale factor" \
	-default $template_size_tolerance
    
    okcancelhelp $w.ok \
	-ok_command "SetTSizeTolerance2 $io $f $t $w" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap4 FIXME"
    
    pack $w.size $w.ok -side top -fill both
}

proc SetTSizeTolerance2 {io f t w} {
    global template_size_tolerance
    
    set new [$w.size get]
    if {$new <= 0} {
	bell
	return
    }

    destroy $w
    set template_size_tolerance $new
    UpdateTemplateDisplay $io $f $t
}
