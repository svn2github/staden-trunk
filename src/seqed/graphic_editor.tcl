###########################
proc ed_scroll_x {id command x args} {
    ed_nip_scroll_canvas -id $id -xscrollcommand "$command $x $args"
}
proc ed_scroll_y {id command y args} {
    ed_nip_scroll_canvas -id $id -yscrollcommand "$command $y $args"
}
proc ed_zoomback {seq_id id} {    
    ed_nip_zoom_canvas -seq_id $seq_id -id $id
}
proc ed_zoom {seq_id id args} {
    eval ed_nip_zoom_canvas -seq_id $seq_id -id $id $args
}
proc ed_resize {id} {
    if {$id != -1} {
	ed_nip_resize_canvas -id $id
    }
}

proc REnzymeMapDia { } {

    set w .renzyme_map

    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "R.Enzyme Search"        
    
    renzymebox $w.tl
    pack $w.tl -fill x -padx 4 -pady 6

    iwidgets::combobox $w.sl \
	-labeltext "Seq identifier:" \
	-labelpos w \
	-completion false \
	-grab global
    pack $w.sl -fill x -padx 4 -pady 10

    frame $w.se -bd 2 -relief sunken -height 2
    pack $w.se -fill x

    iwidgets::buttonbox $w.ok  
    pack $w.ok -expand yes -fill both -padx 4 -pady 1
    $w.ok add ok -text "OK" \
	    -command "renzyme_search_preview_dia $w $w.tl $w.sl"
    $w.ok add cancel -text "Cancel" -command "destroy $w"
    $w.ok add help -text "Help" -command "puts Help"

    set iden [get_sequences_iden]
    $w.sl delete list 0 end
    $w.sl delete entry 0 end
    $w.sl insert entry end [lindex $iden 0]
    foreach id $iden {
	$w.sl insert list end $id
    }
    $w.sl configure -editable false
}

proc renzyme_search_preview_dia {w tl sl} {

    set selRenz [$tl get_name]
    if {$selRenz == ""} {
	tk_messageBox -icon error -type ok \
		-title "restriction enzyme search" \
		-message "No renzyme has been selected"
	raise $w
	return
    }

    set selSeq [$sl get]
    if {$selSeq == ""} {
	tk_messageBox -icon error -type ok \
		-title "restriction enzyme search" \
		-message "No sequence identifier has been selected"
	raise $w
	return
    }
    destroy $w
    set rsp .renzyme_search_preview
    if {[xtoplevel $rsp -resizable 0] == ""} return
    wm title $rsp "R.enzyme search result preview"

    set matches ""
    set matches [renzyme_search_preview $selSeq $selRenz]
   
    set i 0
    set columns { 6 "Enzyme"   
                 22 "Recognition sequence/cut site"
                 10 "Prototype"   
	          6 "Supplier_codes"                  
                  6 "EFS"}
    foreach item [lindex [lindex $matches 0] 5] {
	set n [lindex $selSeq $i]
	lappend columns 5 "$n"
	incr i
    }    
    tablelist::tablelist $rsp.tl \
	    -columns $columns\
            -labelcommand tablelist::sortByColumn \
	    -xscrollcommand [list $rsp.hsb set] \
	    -yscrollcommand [list $rsp.vsb set] \
	    -selectbackground navy -selectforeground white \
	    -height 10 -width 95 -stretch all \
	    -selectmode extended \
	    -exportselection 0

    set num_col [llength $columns]
    set num_col [expr $num_col/2]
    for {set i 1} {$i < $num_col} {incr i 2} {
	$rsp.tl columnconfigure $i -background beige
    }
    for {set i 4} {$i < $num_col} {incr i} {
	$rsp.tl columnconfigure $i -sortmode real
    }

    
    scrollbar $rsp.vsb -orient vertical -command [list $rsp.tl yview]
    scrollbar $rsp.hsb -orient horizontal -command [list $rsp.tl xview]
    
    grid rowconfigure    $rsp.tl 0 -weight 1
    grid columnconfigure $rsp.tl 0 -weight 1

    frame $rsp.se -bd 2 -relief sunken -height 2
  
    iwidgets::buttonbox $rsp.ok  
    $rsp.ok add ok -text "OK" \
	    -command "renzyme_search_draw_map $rsp [list $selSeq]"
    $rsp.ok add cancel -text "Cancel" -command "destroy $rsp"
    $rsp.ok add help -text "Help" -command "puts Help"
  
    grid $rsp.tl -row 0 -column 0 -sticky news
    grid $rsp.vsb -row 0 -column 1 -sticky ns
    grid $rsp.hsb -row 1 -column 0 -sticky ew 
    grid $rsp.se -row 2 -column 0 -sticky w -columnspan 2
    grid $rsp.ok -row 3 -column 0 -sticky ew -columnspan 2

    $rsp.tl delete 0 end
    foreach items $matches {
	set item [lrange $items 0 4]
	foreach i [lindex $items 5] {
	    lappend item $i
	}
	$rsp.tl insert end $item
    }
}

proc renzyme_search_draw_map {path seq_ids} {

    set num_ids [llength $seq_ids]     
    set cur_sel [$path.tl curselection]
    if {$cur_sel == ""} {
	bell
	return 
    }
  
    set num_sel [llength $cur_sel]
    set sel_ren ""
    for {set i 0} {$i < $num_sel} {incr i} {
	set row_sel [lindex $cur_sel $i]
	set row_single [$path.tl get $row_sel]
	set name [lindex $row_single 0]
	set sel_ren [concat $sel_ren $name]
    }
    set seq_name [lindex $seq_ids 0];# FIXME
    set seq_id [get_id_from_name $seq_name]

    destroy $path
    renzyme_draw_map $sel_ren $seq_id $seq_name
}

proc renzyme_draw_map {sel_ren seq_id seq_name} {

    global tk_utils_defs

    set ed_id [add_editor $seq_name] ;# does ed_id needed? 
    set w .renzyme_map

    set num_display [next_renz_display]
    set w $w$num_display
    global $w.renz_id
    set $w.renz_id -1

    set x_scroll "ed_scroll_x \[set $w.renz_id\]"
    set y_scroll "ed_scroll_y \[set $w.renz_id\]"
    set zoomback_command "ed_zoomback $seq_id \[set $w.renz_id\]"
    set zoom_command "ed_zoom $seq_id \[set $w.renz_id\]"
    set cursor_command "ed_nip_canvas_cursor_x -id \[set $w.renz_id\]"
    set invoke_command "ed_nip_canvas_cursor_editor $seq_id \[set $w.renz_id\]"
    set renz_name_command "ed_nip_get_renz_name -id \[set $w.renz_id\]"
    set renz_info_command "ed_nip_get_renz_info -id \[set $w.renz_id\]"
    set config_command "ed_resize \[set $w.renz_id\]"
    set min_win_ht [keylget tk_utils_defs R_ENZ.PLOT_HEIGHT]
    set num_enz [llength $sel_ren]
    set tick_ht [keylget tk_utils_defs R_ENZ.TICK_HEIGHT]
    set text_offset [expr $tick_ht * 1.5]
    set height [expr ($num_enz + 2) * $tick_ht]

    #HACK to ensure window does not get too big!
    if { $height > $min_win_ht} {
	set height $min_win_ht
    }

    #create restriction enzyme canvas
    ed_renz_map $w -zoom_command $zoom_command \
	    -zoomback_command $zoomback_command \
	    -scrollbar_x_cmd $x_scroll \
	    -scrollbar_y_cmd $y_scroll \
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

    EdNipCreateREnzMenu $w

    wm protocol $w WM_DELETE_WINDOW "EdNipREnzStartShutdown $w"
    update
 
    #Do search and plot for selected sequence and R.enzyme

    set $w.renz_id [renz_map_plot $w "create_graphic_editor \
	    -enzymes [list $sel_ren]\
	    -num_enzymes [llength $sel_ren] -seq_id $seq_id"]

    if {[set $w.renz_id] == -1} {
	return
    }

    wm title $w "Restriction Enzyme Editor: $seq_name"

    wm minsize $w [winfo width $w] [winfo height $w]
    wm maxsize $w [winfo screenwidth $w] [winfo height $w]
   
}

##############################################################################

proc EdNipCreateREnzMenu {w } {
    
    global $w.renz_id

    #set up a File menu
    menubutton $w.menubar.file -text "File" -menu $w.menubar.file.opts
    menu $w.menubar.file.opts
    $w.menubar.file.opts add command -label "Exit" \
	    -command "EdNipREnzStartShutdown $w"

    menubutton $w.menubar.help -text "Help" -menu [set m $w.menubar.help.opts]
    menu $m
    #$m add command -label "Introduction" 
	#-command "show_help spin {SPIN-Restrict-Introduction}"
    
    #do the packing
    pack $w.menubar.file -side left
    pack $w.menubar.help -side right
}

##############################################################################
#executed when the exit command is chosen from the File menu of stand alone
#restriction enzyme plot
proc EdNipREnzStartShutdown {w } {

    global $w.renz_id

    if {[info exists $w.renz_id]} {
	editor_result_update -index [set $w.renz_id] -job QUIT
    }
}




