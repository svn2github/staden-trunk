# This file uses mega-widget comborange and tabnotebook(one page for using 
# "Selected range" and another one for using "Feature table").
#
#########################################################################
proc  NipTranslate { } {

    global tk_utils_defs nip_defs 

    set c .translation
    global $c.seq_name

    if {[xtoplevel $c -resizable 0] == ""} return
    wm title $c "Translation"

    iwidgets::tabnotebook $c.tnb -tabpos n \
	    -height 220 \
	    -width 350 
    pack $c.tnb -padx 2 -pady 2 -fill both 
    
# Page1
    set page1 [$c.tnb add -label "Selected Range"]
    frame $page1.f
    pack $page1.f -fill both -expand yes
    comborange $page1.f.seq -textvariable $c.seq_name
    pack $page1.f.seq -fill both
    
    radiobutton .r
    set se [.r cget -selectcolor]
    orientedradiobox $page1.f.dis -labeltext "Display as" \
	    -orient horizontal \
	    -selectcolor $se
    pack $page1.f.dis -fill x -expand yes
    $page1.f.dis add lt -text "3 letter"
    $page1.f.dis add lo -text "1 letter"
    $page1.f.dis select lo
    destroy .r

# Page2
    set page2 [$c.tnb add -label "Feature Table"] 

    frame $page2.f
    pack $page2.f -fill both -expand yes
    comborange $page2.f.iden -feature yes -trange no -textvariable $c.seq_name 
    pack $page2.f.iden -fill both -side top

    $c.tnb view "Selected Range"

    #$page1.f configure -background [$page1.f.seq cget -background]
    #$page2.f configure -background [$page2.f.iden cget -background]
    #$c.tnb configure -background [$page2.f.iden cget -background]
    #$c.tnb configure -backdrop [$page2.f.iden cget -background]
    #$c.tnb configure -tabbackground [$page2.f.iden cget -background]

#############################################################    
    dataentry $c.line \
              -labeltext "Line length" \
              -sticky e \
	      -validate integer \
	      -default 60 \
	      -width 6 
 
    pack $c.line -fill both -padx 6 -pady 6
#############################################################
    frame $c.separate -bd 2 -relief sunken -height 2 -width 10
    pack $c.separate -fill x -expand yes
#############################################################
    iwidgets::buttonbox $c.ok 
    $c.ok add OK -text OK \
                 -command "NipTranslate2 $c $page1.f.seq $page1.f.dis $page2.f.iden " 
    $c.ok add Cancel -text Cancel -command "destroy $c"
    $c.ok add Help -text Help -command "show_help spin {SPIN-Translation-General}"
    pack $c.ok -expand yes -fill both -padx 1 -pady 1
}

proc set_tab { path } {

    set ind  [$path.tnb index select]

    if {$ind == 1} {
	$path.dis buttonconfigure 0 -state disable
	$path.dis select 1
    } else {
	$path.dis buttonconfigure 0 -state normal 
    }
}

proc  NipTranslate2 {c seq dis iden} {
       
    global PROTEIN

    #set size [$c.dis get]
    set size [$dis get]
    set selcds [$iden get_curselection] 
    set length [$c.line get]
    set name [$seq get_seqname]

    if {$length > 600000} {
	set length 600000
    }

    set seq_id [name_to_seq_id $name]
    if {$size == "lt"} {
	set size 3
    } else {set size 1}

    set feat [$c.tnb view]

    if {$feat != 1 } {
	set feat 2
    } 
    
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "translate" "unable to process protein sequences"
	destroy $c
	return
    }
    if {$length <= 0} {
	bell
	verror ERR_WARN "Translation" "line length must be positive"
	return
    }
    
    if {$feat == 1} {
	if {[llength $selcds] == 0} {
	    tk_messageBox -title "Error" -message "No selection has been made" \
                                         -icon error -type ok
	raise $c
	return
	}
    }
    SetBusy
    nip_translate -seq_id $seq_id \
	    -start [$seq get_s] \
	    -end  [$seq get_e]\
	    -line_length $length \
	    -size $size \
	    -feat $feat \
	    -selcds $selcds
    ClearBusy
    destroy $c
}

                          
