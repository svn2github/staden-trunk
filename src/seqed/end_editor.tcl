proc EditHang {} {

    set w .edithang

    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Fragment Inseration" 

    label $w.l -text "The fragment ends are incompatible, do you want:"
    pack $w.l -fill both -expand yes -padx 10 -pady 20

    frame $w.se -bd 2 -relief sunken -height 2
    pack $w.se -fill x

    iwidgets::buttonbox $w.b 
    pack $w.b -fill x -expand yes 
    
    foreach b {edit cancel help} {
	$w.b add $b -text [string totitle $b] \
		-command "set $w.Pressed $b; destroy $w"
    }   
    global $w.Pressed
    set $w.Pressed ""
    vwait $w.Pressed
    return [set $w.Pressed]
}

#called from c
proc EndEditor {seq_id} {

    set w .hangeditor
    if {[xtoplevel $w -resizable 0] == ""} return
 
    wm title $w "Ends Editor"

    global $w.ee_id
    set $w.ee_id -1
                   
    set $w.ee_id [end_editor_register $w $w.ee $seq_id]
    if {[set $w.ee_id] == -1} {
	return
    }
    create_end_editor $w $seq_id
}

#called from c
proc DeleteEndEditor {w} {

    global $w.ee_id
    
    unset $w.ee_id
    destroy $w
}

proc create_end_editor {w seq_id} {
    
    ##########
    entry .eq
    set bgg [.eq cget -background]
    set fgg [.eq cget -foreground]
    destroy .eq    
    ##########
    #FIXME: use default configuretion
    #option add *Hangeditor*canvas.background $bgg userDefault
    hangeditor $w.ee -seqid $seq_id
    pack $w.ee -expand yes -fill both -padx 6 -pady 6
    [$w.ee component canvas] configure -background $bgg

    frame $w.se -bd 2 -relief sunken -height 2
    pack $w.se -fill x

    iwidgets::buttonbox $w.b 
    pack $w.b -fill x -padx 2 -pady 2

    foreach b {ok cancel help} {
	$w.b add $b -text [string totitle $b] \
		-command "set $w.buttonP $b; destroy $w"
    }
    global $w.buttonP
    set $w.buttonP ""
    eval vwait $w.buttonP
    return [set $w.buttonP]
}

proc EndEditorRedisplay {win seq_id} {

    $win init_end_editor $seq_id

}

proc HangEditor {seq_id} {

    set w .hangeditor
    if {[xtoplevel $w -resizable 0] == ""} return
 
    wm title $w "Ends Editor"
    ##########
    entry .eq
    set bgg [.eq cget -background]
    set fgg [.eq cget -foreground]
    destroy .eq    
    ##########
    #FIXME: use default configuretion
    #option add *Hangeditor*canvas.background $bgg userDefault

    hangeditor $w.c -seqid $seq_id
    pack $w.c -expand yes -fill both -padx 6 -pady 6
    [$w.c component canvas] configure -background $bgg

    frame $w.se -bd 2 -relief sunken -height 2
    pack $w.se -fill x
  
    iwidgets::buttonbox $w.b 
    pack $w.b -fill x -padx 2 -pady 2
    
    foreach b {ok cancel help} {
	$w.b add $b -text [string totitle $b] \
		-command "set $w.buttonP $b; destroy $w "
    }

    global $w.buttonP
    set $w.buttonP ""
    vwait $w.buttonP
    return [set $w.buttonP]
}


