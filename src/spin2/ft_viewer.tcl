proc feature_key_array { } {
    global feat_names

    set list ""
    set cnt 0
    set key_names [get_feature_keys]
    foreach key $key_names {
	lappend list $key $cnt
	incr cnt
    }
    array set feat_names $list
}

proc feature_name2id {items} {
    global feat_names

    set list ""
    foreach i $items {
	set name [keylget i type]
	set id $feat_names($name)
	keylset i index $id
	lappend list $i
    }
    return $list
}

proc next_ft_viewer { } {
    global num_ft

    if {![info exists num_ft]} {
	set num_ft 0
    } else {
	incr num_ft
    }
    return $num_ft
}

proc get_ft_name { } {

    set t .ft_viewers
    set num [next_ft_viewer]
    return $t$num

}

proc save_default_feats {list name {version {}}} {
    global feats_$name$version

    set feats_$name$version $list
}

proc get_default_feats {name {version {}}} {
    global feat_names feats_$name$version env

    #set default to all
    if {![info exists feats_$name$version]} {
	read_feature_db -file $env(STADTABL)/FEATUREDB
	set key_info [get_feature_db]
	set list_all [feature_name2id $key_info]
	save_default_feats $list_all $name $version	
    }
    return [set feats_$name$version]
}

proc ft_save_selection {e_id fs name {version {}}} {

    set items [$fs get]
    set list [feature_name2id $items]
    save_default_feats $list $name $version

    ft_viewer_update -id $e_id
}

proc ft_selector {name parent t e_id} {
    global spin_defs

    #first check if toplevel is hiding and deiconify it
    if {[winfo exists $t]} {
	if {[string compare [wm state $t] withdrawn] == 0} {
	    wm deiconify $t
	    return
	}
    }
    if {[xtoplevel $t -resizable 0] == ""} {return}

    wm protocol $t WM_DELETE_WINDOW "wm withdraw $t"
    wm title $t "Feature selector"

    feature_selector $t.fs

    okcancelhelp $t.button -bd 2 -relief groove\
	    -ok_command " ft_save_selection $e_id $t.fs $name $e_id"\
	    -cancel_command "wm withdraw $t"\
	    -help_command ""

    pack $t.fs
    pack $t.button -side bottom -fill x
}

proc ft_create_menubar {w} {
    global spin_defs
    global ft_viewer_menu

    set menubar $w[keylget spin_defs FT.MENUBAR]

    $w configure -menu $menubar
    menu $menubar
    create_menus $ft_viewer_menu $menubar
}

proc ft_create_buttonbar {w } {
    global spin_defs

    set buttons $w[keylget spin_defs FT.BUTTONBAR]
    set crosshair [keylget spin_defs FT.CROSSHAIR]
    set cursor [keylget spin_defs FT.CURSOR]
    set window [keylget spin_defs FT.SINGLE.F]

    frame $buttons
    #zoom back button
    button $buttons.back -text "zoom out"
    button $buttons.zoomin10 -text "+10%" \
	-command "ZoomInCanvas $w$window 0.05" 
    button $buttons.zoomin50 -text "+50%" \
	-command "ZoomInCanvas $w$window 0.1666" 
  
    #cursor checkbutton
    global $w.crosshair
    checkbutton $buttons$crosshair -text crosshairs -variable $w.crosshair

    #cursor position label
    label $w$cursor -bd 2 -relief sunken -width 6
    pack $buttons.zoomin10 $buttons.zoomin50 $buttons.back \
	 $buttons$crosshair -expand yes -side left
    pack $w$cursor -in $buttons -side left -expand yes 

    grid $buttons -row 1 -column 0 -sticky ew -columnspan 3
}

proc create_ft_viewer_horizontal {ft_width win_height ruler_height seq_id} {
    global spin_defs

    set w [get_ft_name]
    toplevel $w
    #fix_maxsize $w

    global $w.ft_viewer_id
    set $w.ft_viewer_id -1

    #set up names
    set forward $w[keylget spin_defs FT.SINGLE.F]
    set reverse $w[keylget spin_defs FT.SINGLE.R]
    set ruler $w[keylget spin_defs FT.RULER.WIN]
    set hscroll $w[keylget spin_defs FT.HSCROLL]
    set brief $w[keylget spin_defs FT.BRIEF] 
    set vscroll_f $w[keylget spin_defs FT.VSCROLL.F]
    set vscroll_r $w[keylget spin_defs FT.VSCROLL.R]

    set cursor_command "nip_canvas_cursor_x -id \[set $w.ft_viewer_id\]"
    set zoomback_command "nip_zoomback $seq_id \[set $w.ft_viewer_id\]"
    set zoom_command "nip_zoom $seq_id \[set $w.ft_viewer_id\] -direction x" 
    set invoke_command "nip_canvas_cursor_editor $seq_id \[set $w.ft_viewer_id\]"
    set config_command "nip_resize \[set $w.ft_viewer_id\]"

    ft_create_menubar $w
    ft_create_buttonbar $w 

    canvasbox $forward -bd 0 -highlightthickness 0 -width $ft_width \
	    -height $win_height -xscrollcommand "$hscroll set " \
	    -bd 2 -relief groove -zoom_command $zoom_command\
	    -yscrollcommand "$vscroll_f set"

    $forward bind all <Any-Motion> "highlight_ft_item $forward $brief" 
    
    ##########################################################################
    #create ruler canvas
    canvasbox $ruler -height $ruler_height -bd 2 -relief groove -highlightthickness 0\
	    -width $ft_width -xscrollcommand "$hscroll set " -zoom_command $zoom_command

    canvasbox $reverse -bd 0 -highlightthickness 0 -width $ft_width \
	    -height $win_height -xscrollcommand "$hscroll set " \
	    -bd 2 -relief groove -zoom_command $zoom_command\
	    -yscrollcommand "$vscroll_r set"

    $reverse bind all <Any-Motion> "highlight_ft_item $reverse $brief" 

    scrollbar $hscroll -orient horizontal 
    scrollbar $vscroll_f -orient vertical
    scrollbar $vscroll_r -orient vertical

    label $brief 

    grid columnconfig $w 1 -weight 1
    grid rowconfig $w 2 -weight 1

    grid $forward -row 2 -column 1 -sticky nsew
    grid $ruler -row 3 -column 1 -sticky ew
    grid $reverse -row 4 -column 1 -sticky nsew
    grid $hscroll -row 5 -column 1 -sticky ew
    grid $brief -row 6 -column 1 -sticky ew

    grid $vscroll_f -row 2 -column 2 -sticky ns
    grid $vscroll_r -row 4 -column 2 -sticky ns

    fit_on_screen $w
    return $w
}

proc create_ft_viewer_vertical {ft_width win_height ruler_height seq_id} {
    global spin_defs

    set w [get_ft_name]
    toplevel $w
    fix_maxsize $w

    global $w.ft_viewer_id
    set $w.ft_viewer_id -1

    #set up names
    set forward $w[keylget spin_defs FT.SINGLE.F]
    set reverse $w[keylget spin_defs FT.SINGLE.R]
    set ruler $w[keylget spin_defs FT.RULER.WIN]
    set vscroll $w[keylget spin_defs FT.HSCROLL]
    set brief $w[keylget spin_defs FT.BRIEF] 

    set cursor_command "nip_canvas_cursor_x -id \[set $w.ft_viewer_id\]"
    set zoomback_command "nip_zoomback $seq_id \[set $w.ft_viewer_id\]"
    set zoom_command "nip_zoom $seq_id \[set $w.ft_viewer_id\] -direction x" 
    set invoke_command "nip_canvas_cursor_editor $seq_id \[set $w.ft_viewer_id\]"
    set config_command "nip_resize \[set $w.ft_viewer_id\]"

    ft_create_menubar $w
    ft_create_buttonbar $w 

    canvasbox $forward -bd 0 -highlightthickness 0 -width $win_height \
	    -height $ft_width -yscrollcommand "$vscroll set " \
	    -bd 2 -relief groove -zoom_command $zoom_command

    $forward bind all <Any-Motion> "highlight_ft_item $forward $brief" 
    
    ##########################################################################
    #create ruler canvas
    canvasbox $ruler -height $ft_width -bd 2 -relief groove -highlightthickness 0\
	    -width $ruler_height -xscrollcommand "$vscroll set " -zoom_command $zoom_command

    canvasbox $reverse -bd 0 -highlightthickness 0 -width $win_height \
	    -height $ft_width -yscrollcommand "$vscroll set " \
	    -bd 2 -relief groove -zoom_command $zoom_command

    $reverse bind all <Any-Motion> "highlight_ft_item $reverse $brief" 

    scrollbar $vscroll -orient vertical 

    label $brief 

    grid columnconfig $w 1 -weight 1
    grid rowconfig $w 2 -weight 1

    grid $forward -row 2 -column 1 -sticky nsew
    grid $ruler -row 2 -column 2 -sticky ns
    grid $reverse -row 2 -column 3 -sticky nsew
    grid $vscroll -row 2 -column 4 -sticky ns
    grid $brief -row 3 -column 1 -sticky ew
    return $w
}

proc create_ft_circle {ft_width forward_height ruler_height seq_id} { 
    global spin_defs

    set w [get_ft_name]
    toplevel $w
    fix_maxsize $w

    global $w.ft_viewer_id
    set $w.ft_viewer_id -1

    #set up names
    set circle $w[keylget spin_defs FT.CIRCLE.WIN]
    set hscroll $w[keylget spin_defs FT.HSCROLL]
    set vscroll $w[keylget spin_defs FT.VSCROLL.C]
    set brief $w[keylget spin_defs FT.BRIEF] 

    ft_create_menubar $w
    ft_create_buttonbar $w 

    canvasbox $circle -bd 0 -highlightthickness 0 -width $ft_width \
	    -height $ft_width -xscrollcommand "$hscroll set " \
	    -yscrollcommand "$vscroll set " \
	    -bd 2 -relief groove

    scrollbar $hscroll -orient horizontal
    scrollbar $vscroll -orient vertical
    label $brief 

    grid columnconfig $w 1 -weight 1
    grid rowconfig $w 2 -weight 1

    grid $circle -row 2 -column 1 -sticky nsew
    grid $hscroll -row 3 -column 1 -sticky ew
    grid $vscroll -row 2 -column 2 -sticky ns
    grid $brief -row 4 -column 1 -sticky ew

    return $w
}

proc display_ft_viewer {seq_id result_id results orientation strand mode offset} {
    global spin_defs nip_defs tk_utils_defs TOP_S BOTTOM_S CIRCLE HORIZONTAL 
    global ft_viewer_menu

    set type [keylget nip_defs FT_VIEWER]

    set frame 0
    if {$orientation == $CIRCLE} {
	set element_info [create_seq_element -1 -1 $type $strand $frame \
		$orientation XY "CANVAS" [keylget spin_defs CONTAINER.TITLE] \
		[keylget spin_defs FT.CIRCLE.PLOT_WIDTH] \
		[keylget spin_defs FT.CIRCLE.PLOT_WIDTH]]

    } else {
	set element_info [create_seq_element $seq_id -1 $type $strand $frame \
		$orientation X "CANVAS" [keylget spin_defs CONTAINER.TITLE] \
		[keylget spin_defs FT.PLOT_WIDTH] \
		[keylget spin_defs FT.SINGLE.PLOT_HEIGHT]]
    }
    set c_win [keylget element_info container_win]
    set c_id [keylget element_info container_id]
    set e_win [keylget element_info element_win]
    set e_id [keylget element_info element_id]

    #set the default orientation depending on orientation of seq_id except
    #with a circle!
    if {$orientation != $CIRCLE} {
	set orientation [keylget element_info orientation]
    }

    container_menubar $c_win $c_id $ft_viewer_menu
    set_element_specs $c_win$e_win "menu [list $ft_viewer_menu]"

    update idletasks

    ft_viewer plot -element $c_win$e_win\
	    -container $c_win\
	    -seq_id $seq_id \
	    -result_id $result_id\
	    -results $results\
	    -container_id [keylget element_info container_id]\
	    -element_id [keylget element_info element_id]\
	    -element_type "CANVAS"\
	    -orientation $orientation\
	    -display_mode $mode\
	    -offset $offset\
	    -origin [keylget spin_defs FT.CIRCLE.ORIGIN]

    seqed_element_bindings $c_id $c_win$e_win $e_id

    #update result list
    result_list_update $c_win     

    #set brief $c_win[keylget tk_utils_defs CONTAINER.BRIEF.WIN]
    #$c_win$e_win bind all <Any-Motion> "+highlight_ft_item $c_win$e_win $brief" 
}

proc feature_table { } {
    global spin_defs

    set t .feature_table

    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Feature table"
  
    #select plot range
    set seq_id [get_active_seq_id] 
    global $seq_id.start $seq_id.end

    set seq_length [seq_info $seq_id length] 
    set seq_start [seq_info $seq_id start] 
    set seq_end [seq_info $seq_id end] 
    if {[info exists $seq_id.start]} {
	set seq_start [set $seq_id.start]
    }    
    if {[info exists $seq_id.end]} {
	set seq_end [set $seq_id.end]
    }

    keylset us RANGE [keylget spin_defs FT.RANGE]
    seq_id $t.range -range 1 -browse 1 -from 1 -to $seq_length \
	-start_value $seq_start -end_value $seq_end -min_value 1 \
	-default [seq_info $seq_id name]\
 	-update_cmd [list [list seq_range_updates $t.range]]\
	-browse_cmd nip_seq_browser

    keylset ft TYPE [keylget spin_defs FT.TYPE]
    set b1 [keylget ft TYPE.BUTTON.1]
    global $t.type_var
    frame $t.t
    checkbutton $t.t.type -text [keylget ft TYPE.BUTTON.1] -variable $t.type_var

    #set b1 [keylget ft TYPE.BUTTON.1]
    #set b2 [keylget ft TYPE.BUTTON.2]
    #set b3 [keylget ft TYPE.BUTTON.3]
    #radiolist $t.type \
	    -title [keylget ft TYPE.NAME]\
	    -default [keylget ft TYPE.VALUE]\
	    -orient horizontal\
	    -buttons [format {{%s} {%s} {%s}} \
		      [list $b1] [list $b2] [list $b3]]
	    
    #strand selection
    strand_both $t.strand    

    #display mode
    keylset ft MODE [keylget spin_defs FT.MODE.DMODE]
    set b1 [keylget ft MODE.BUTTON.1]
    set b2 [keylget ft MODE.BUTTON.2]
    set b3 [keylget ft MODE.BUTTON.3]
    radiolist $t.mode \
	    -title [keylget ft MODE.NAME]\
	    -default [keylget ft MODE.VALUE]\
	    -orient vertical\
	    -buttons [format {{%s} {%s} {%s}} \
		      [list $b1] [list $b2] [list $b3]]

    keylset ft OFFSET [keylget spin_defs FT.MODE.OFFSET]
    entrybox $t.offset \
	    -title [keylget ft OFFSET.NAME]\
	    -default [keylget ft OFFSET.SINGLE.VALUE]\
	    -width 5\
	    -type "CheckFloat"

    #########################################################################
    #OK and Cancel buttons
    okcancelhelp $t.ok_cancel \
	    -ok_command "feature_table2 $t $t.range \[strand_get $t.strand\] $t.mode $t.offset; destroy $t" \
	    -cancel_command "seq_id_destroy $t.range; destroy $t" \
	    -help_command "show_help spin {SPIN-????}" \
	    -bd 2 \
	    -relief groove

    #final packing
    pack $t.range -fill both
    pack $t.strand -fill both
    pack $t.t.type -anchor w 
    pack $t.t -fill x 
    pack $t.mode -fill both
    pack $t.offset -fill both
    pack $t.ok_cancel -fill x
}


proc feature_table2 {t range strand mode offset} {
    global CIRCLE HORIZONTAL TOP_S BOTTOM_S $t.type_var

    set seq_id [name_to_seq_id [seq_id_name $range]]
    set start [seq_id_from $range]
    set end [seq_id_to $range]
    #set orientation [radiolist_get $type]

    set plot_mode [radiolist_get $mode]
    set plot_offset [entrybox_get $offset]
    
    if {[set $t.type_var]} {
	set orientation $CIRCLE
    } else {
	set orientation $HORIZONTAL
    }

    SetBusy


    if {$orientation == $CIRCLE} {

	puts "orientation $orientation"
	set res [ft_viewer create -seq_id $seq_id \
		-start [seq_id_from $range] \
		-end [seq_id_to $range]\
		-strand $strand\
		-orientation $orientation]

	#remove id from start of list
	set result_id [lindex $res 0]
	set results [lindex $res 1]    

	display_ft_viewer $seq_id $result_id $results $orientation $strand $plot_mode $plot_offset
    } else {

	if {[expr $strand & $TOP_S]} {
	    set res [ft_viewer create -seq_id $seq_id \
		    -start [seq_id_from $range] \
		    -end [seq_id_to $range]\
		    -strand $TOP_S\
		    -orientation $orientation]

	    #remove id from start of list
	    set result_id [lindex $res 0]
	    set results [lindex $res 1]    
	    
	    display_ft_viewer $seq_id $result_id $results $orientation $TOP_S $plot_mode $plot_offset
	} 
	if {[expr $strand & $BOTTOM_S]} {
	    set res [ft_viewer create -seq_id $seq_id \
		    -start [seq_id_from $range] \
		    -end [seq_id_to $range]\
		    -strand $BOTTOM_S\
		    -orientation $orientation]

	    #remove id from start of list
	    set result_id [lindex $res 0]
	    set results [lindex $res 1]    

	    display_ft_viewer $seq_id $result_id $results $orientation $BOTTOM_S $plot_mode $plot_offset
	} 
    }

    ClearBusy
    

}

proc highlight_ft_item {window} {

    set nearest [$window find withtag current]
    set tags [$window gettags $nearest] 

    set i [lsearch $tags current]
    set tags [lreplace $tags $i $i]

    return $tags
    #$brief configure -text $tags

}

proc highlight_ft_itemOLD {window brief} {

    set nearest [$window find withtag current]
    set tags [$window gettags $nearest] 

    set i [lsearch $tags current]
    set tags [lreplace $tags $i $i]

    $brief configure -text $tags

}

proc update_ft_viewer {w mode} {

    global $w.ft_viewer_id

    puts "ft_viewer_update $w $mode"
    ft_viewer_update -id [set $w.ft_viewer_id] -display_mode $mode

}

proc ft_viewer_display_mode {e_win} {
    global spin_defs

    set t .display_mode

    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Display mode"
  
    keylset ft MODE [keylget spin_defs FT.MODE.DMODE]
    set b1 [keylget ft MODE.BUTTON.1]
    set b2 [keylget ft MODE.BUTTON.2]
    set b3 [keylget ft MODE.BUTTON.3]
    radiolist $t.mode \
	    -title [keylget ft MODE.NAME]\
	    -default [keylget ft MODE.VALUE]\
	    -orient vertical\
	    -buttons [format {{%s} {%s} {%s}} \
		      [list $b1] [list $b2] [list $b3]]
	    
    strand_both $t.strand

    keylset ft OFFSET [keylget spin_defs FT.MODE.OFFSET]
    entrybox $t.offset \
	    -title [keylget ft OFFSET.NAME]\
	    -default [keylget ft OFFSET.SINGLE.VALUE]\
	    -width 5\
	    -type "CheckFloat"


    #########################################################################
    #OK and Cancel buttons
    okcancelhelp $t.ok_cancel \
	    -ok_command "ft_viewer_display_mode2 $e_win $t $t.mode \[strand_get $t.strand] $t.offset; destroy $t" \
	    -cancel_command "seq_id_destroy $t.range; destroy $t" \
	    -help_command "show_help spin {SPIN-????}" \
	    -bd 2 \
	    -relief groove
    
    pack $t.mode -fill x
    pack $t.strand -fill both
    pack $t.offset -fill x
    pack $t.ok_cancel -fill x

}

proc ft_viewer_display_mode2 {e_win t mode strand offset} {
    global $w.ft_viewer_id

    set mode [radiolist_get $mode]
    set strand [strand_get $strand]
    set offset [entrybox_get $offset]
    
    set e_id [get_element_id $e_win]

    ft_viewer_update -id $e_id -display_mode $mode -offset $offset -strand $strand
}
