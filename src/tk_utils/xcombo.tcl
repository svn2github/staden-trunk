image create bitmap xcombo_bitmap \
-data "#define down_arrow_width 15
\#define down_arrow_height 14
static unsigned char down_arrow_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0xfe, 0x3f, 0xfc, 0x1f, 0xf8, 0x0f, 0xf0, 0x07,
   0xe0, 0x03, 0xc0, 0x01, 0x80, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00};"

# Poor mans combobox until we switch to tk8.5
proc xcombo {w args} {
    global $w

    set ${w}(-text) "label"
    set ${w}(-textvariable) ${w}(value)
    set ${w}(-valuesvariable) ${w}(values)
    set ${w}(-values) ""
    set ${w}(-command) ""
    set ${w}(-validate) ""
    set ${w}(-validatecommand) ""
    foreach {k v} $args {
	set ${w}($k) $v
    }
    
    uplevel \#0 [list set [set ${w}(-valuesvariable)] [set ${w}(-values)]]

    frame $w
    label $w.label \
	-text [set ${w}(-text)]
    entry $w.entry \
	-textvariable [set ${w}(-textvariable)]
    #button $w.dot -text "..." -padx 1m -command "xcombo_dots $w"
    button $w.dot \
	-image xcombo_bitmap \
	-padx 1m \
	-command "xcombo_dots $w"

    pack $w.label -side left
    pack $w.dot $w.entry -side right

    return $w
}

;proc xcombo_set {w val} {
    $w.entry delete 0 end
    $w.entry insert 0 $val
}

;proc xcombo_configure {w args} {
    global $w

    set r [eval $w.entry configure $args]
    set ${w}(-textvariable) [$w.entry cget -textvariable]

    return $r
}

;proc xcombo_set_values {w args} {
    global $w
    set v [set ${w}(-valuesvariable)]
    eval set ::[list $v] $args
}

;proc xcombo_dots {w} {
    global $w

    if {[winfo exists $w.list]} {
	raise $w.list
	wm deiconify $w.list
	return
    }

    # Create a position a top-level window
    set x [winfo rootx $w.dot]
    set y [winfo rooty $w.dot]
    incr x 5
    incr y 5

    set W [winfo screenwidth $w]
    set H [winfo screenheight $w]
    if {$x+150 > $W} {
	set x [expr {$W-150}]
    }
    if {$y+200 > $H} {
	if {$y-200 > 0} {
	    incr y -200
	} else {
	    set y [expr {$H-200}]
	}
    }

    toplevel $w.list -width 150 -height 200 -class ComboList
    wm geometry $w.list 150x200+$x+$y
    wm transient $w.list $w
    wm overrideredirect $w.list 1
    
    set l [listbox $w.list.l \
	       -exportselection 0 \
	       -yscrollcommand "$w.list.ys set" \
	       -width 10]
    scrollbar $w.list.ys -command "$l yview" -orient vertical
    set v [uplevel \#0 [list set [set ${w}(-valuesvariable)]]]
    eval $l insert end $v
    if {"$v" == ""} {
	$l insert end ""
    }
    pack $w.list.l -fill both -expand 1 -side left
    pack $w.list.ys -fill both -side right

    update idletasks
    grab -global $w.list

    bind ComboList <ButtonRelease-1> "destroy $w.list"
    focus $w.list
    bind ComboList <Key-Escape> "destroy $w.list"

    bind $l <<ListboxSelect>> "
        set \[set ${w}(-textvariable)\] \[%W get \[%W curselection\]\]
        if {\[set ${w}(-command)\] != {}} {
            eval \[set ${w}(-command)\] \[list $w \[%W get \[%W curselection\]\]\]
        }
    "
    bind $l <ButtonRelease-1> "+after idle {destroy $w.list}"
    
}


# xcombo .test \
#     -text "Combobox test" \
#     -textvariable "fu bar" \
#     -values [list a b "c d" e]
# 
# pack .test -side bottom -fill both
# 
# pack [button .b -text Press\ me -command {puts [set "fu bar"]}]