proc messagebox_append { f str } {
    
    $f.text config -state normal
    $f.text insert end $str
    $f.text config -state disabled
}

#message box for information
proc messagebox {w str} {

    set f $w.mbox

    if {[winfo exists $f]} {messagebox_append $f $str; return}
    toplevel $f

    wm title $f Information

    wm transient $f [winfo toplevel [winfo parent $f]]

    text $f.text -wrap none -xscrollcommand "$f.sb_h set" \
	    -yscrollcommand "$f.sb_v set" \
	    -width 60 -height 15 -bd 2 -relief groove 
    scrollbar $f.sb_h -orient horizontal -relief sunken \
	    -command "$f.text xview"
    scrollbar $f.sb_v -orient vertical -relief sunken \
	    -command "$f.text yview"

    button $f.ok -text OK -command "grab release $f; destroy $f"

    grid columnconfig $f 0 -weight 1
    grid rowconfig $f 0 -weight 1

    grid $f.text -row 0 -column 0 -sticky nsew
    grid $f.sb_v -row 0 -column 1 -sticky ns
    grid $f.sb_h -row 1 -column 0 -sticky ew
    grid $f.ok -row 2 -column 0 

    # 6. Withdraw the window, then update all the geometry information
    # so we know how big it wants to be, then center the window in the
    # display and de-iconify it.

    wm withdraw $f
    update idletasks
    
    set x [expr [winfo screenwidth $f]/2 - [winfo reqwidth $f]/2 \
	    - [winfo vrootx [winfo parent $f]]]
    set y [expr [winfo screenheight $f]/2 - [winfo reqheight $f]/2 \
	    - [winfo vrooty [winfo parent $f]]]
    wm geom $f +$x+$y
    wm deiconify $f

    # 7. Set a grab and claim the focus too.
    #grab $f 
    focus $f.ok

    $f.text insert end $str
    $f.text config -state disabled
}


