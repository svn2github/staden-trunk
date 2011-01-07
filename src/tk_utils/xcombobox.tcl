## #TEST
## source $env(STADTABL)/shlib.conf
## load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}
## load_package tk_utils
## tk_utils_init

package provide Xcombobox 1.0

# Create this to make sure there are registered in auto_mkindex
# these must come before the [widget create ...]
proc Xcombobox args {}
proc xcombobox args {}

if {[info command ttk::combobox] == ""} {

#-----------------------------------------------------------------------------
# Our own combobox construction

widget create Xcombobox -type frame -base entry -components {
    xlabel
    {button button button {-image ::Widget::Xcombobox:bitmap \
			   -command "[namespace current]::dots $w"}}
} -options {
    {-label		label		Label		{}}
    {-default		defaultText	DefaultText	0}
    {-state		state		State		normal}
    {-valuesvariable	valuesVariable	ValuesVariable	{}}
    {-bg		-background}
    {-background	background	Background	{#d9d9d9}}
    {-entrybg		ALIAS entry -background}
    {-fg		-foreground}
    {-foreground	ALIAS entry -foreground}
    {-values 		values		Values		{}}
    {-command 		command		Command		{}}
}

namespace eval ::Widget::Xcombobox {;

;proc construct {w} {
    global tcl_platform
    variable $w
    upvar 0 $w data

    # Some defaults
    set data(-valuesvariable) ::${w}(values)
    uplevel \#0 [list set $data(-valuesvariable) ""]

    pack $data(xlabel) -side left -fill both
    pack $data(button) $data(entry) -side right

    bind $data(entry) <Return> "$w get"

    # Hack for windows
    if {$tcl_platform(platform) == "windows"} {
	catch {[winfo parent $data(entry)] configure -bg SystemButtonFace}
    }
}

;proc configure {w args} {
    variable $w
    upvar 0 $w data

    set truth {^(1|yes|true|on)$}

    foreach {key val} $args {
	switch -- $key {
	    -label {
		$data(xlabel) configure -text $val
	    }
	    -default {
		$data(entry) delete 0 end
		$data(entry) insert 0 $val
	    }
	    -state {
		$data(xlabel) configure -state $val
		$data(entry)  configure -state $val
		$data(button) configure -state $val
	    }
	    -background {
		$data(xlabel)    configure -bg $val
		$data(button)    configure -bg $val
		$data(container) configure -bg $val
	    }
	    -entrybg {
		$data(entry)     configure -bg $val
	    }
	    -foreground {
		$data(xlabel) configure -fg $val
		$data(button) configure -fg $val
		$data(entry)  configure -fg $val
	    }
	    -values {
		uplevel \#0 [list set $data(-valuesvariable) $val]
	    }
	}
	set data($key) $val
    }
}

;proc dots {w} {
    variable $w
    upvar 0 $w data

    if {[winfo exists $w.list]} {
	raise $w.list
	wm deiconify $w.list
	return
    }

    # Create a position a top-level window
    set x [winfo rootx $w.entry]
    set y [winfo rooty $w.entry]
    incr y [winfo height $w.entry]
    incr y -2

    set W [expr {[winfo width $w.entry] + [winfo width $w.button]}]
    set H [winfo screenheight $w]

    if {$y+200 > $H} {
	if {$y-200 > 0} {
	    incr y -200
	} else {
	    set y [expr {$H-200}]
	}
    }

    toplevel $w.list -width $W -height 200 -class ComboList
    wm geometry $w.list ${W}x200+$x+$y
    wm transient $w.list $w
    wm overrideredirect $w.list 1
    
    set l [listbox $w.list.l \
	       -exportselection 0 \
	       -yscrollcommand "$w.list.ys set" \
	       -width 10]
    scrollbar $w.list.ys -command "$l yview" -orient vertical
    if {[uplevel \#0 [list info exists [set ${w}(-valuesvariable)]]]} {
	set v [uplevel \#0 [list set [set ${w}(-valuesvariable)]]]
	eval $l insert end $v
    } else {
	set v ""
    }
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
	upvar 0 [namespace current]::$w data
        \$data(entry) delete 0 end
        \$data(entry) insert 0 \[%W get \[%W curselection\]\]
        if {\$data(-command) != {}} {
            eval \$data(-command) \[list $w \[%W get \[%W curselection\]\]\]
        }
    "
    bind $l <ButtonRelease-1> "+after idle {destroy $w.list}"
}

;proc _set {w val} {
    variable $w
    upvar 0 $w data

    $data(entry) delete 0 end
    $data(entry) insert 0 $val
}

# Returns the index into the listbox; for compatibility with iwidgets version
;proc _curselection {w} {
    variable $w
    upvar 0 $w data

    set ind [lsearch -exact [$w cget -values] [$data(entry) get]]
    if {$ind >= 0} {
	return $ind
    } else {
	return ""
    }
}

}; # end namespace eval ::Widget::Xcombobox

# Arrow bitmap for invoking the listbox
image create bitmap ::Widget::Xcombobox:bitmap \
-data "#define down_arrow_width 15
\#define down_arrow_height 14
static unsigned char down_arrow_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0xfe, 0x3f, 0xfc, 0x1f, 0xf8, 0x0f, 0xf0, 0x07,
   0xe0, 0x03, 0xc0, 0x01, 0x80, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00};"


#-----------------------------------------------------------------------------
# Otherwise wrap up the ttk::combobox widget

} else {

proc ttkcombobox {args} {
    return [eval ttk::combobox $args]
}

widget create Xcombobox -type frame -base ttkcombobox -components {
    xlabel
} -options {
    {-label		label		Label		{}}
    {-default		defaultText	DefaultText	0}
    {-values		values		Values		{}}
    {-valuesvariable	valuesVariable	ValuesVariable	{}}
    {-command 		command		Command		{}}
}

namespace eval ::Widget::Xcombobox {;

;proc construct {w} {
    global tcl_platform
    variable $w
    upvar 0 $w data

    pack $data(xlabel) -side left -fill both
    pack $data(ttkcombobox) -side right -fill both
    bind $data(ttkcombobox) <Return> "$w get"
    bind $data(ttkcombobox) <<ComboboxSelected>> "[namespace current]::selected $w"

    # Hack for windows
    if {$tcl_platform(platform) == "windows"} {
	catch {[winfo parent $data(ttkcombobox)] configure -bg SystemButtonFace}
    }
}

;proc destruct {w} {
    variable $w
    upvar 0 $w data

    if {$data(-valuesvariable) != ""} {
	uplevel \#0 [list trace remove variable $data(-valuesvariable) write "[namespace current]::values_changed $w"]
    }
}

;proc configure {w args} {
    variable $w
    upvar 0 $w data

    set truth {^(1|yes|true|on)$}

    foreach {key val} $args {
	switch -- $key {
	    -label {
		$data(xlabel) configure -text $val
	    }
	    -default {
		$data(ttkcombobox) set $val
	    }
	    -values {
		$data(ttkcombobox) configure -values $val
	    }
	    -state {
		$data(xlabel) configure -state $val
		$data(ttkcombobox)  configure -state $val
	    }
	    -valuesvariable {
		if {$data(-valuesvariable) != ""} {
		    uplevel \#0 [list trace remove variable $data(-valuesvariable) write "[namespace current]::values_changed $w"]
		}
		uplevel \#0 [list trace add variable $val write "[namespace current]::values_changed $w"]
	    }
	}
	set data($key) $val
    }
}

;proc values_changed {w varname args} {
    variable $w
    upvar 0 $w data

    global $varname
    after idle [list $data(ttkcombobox) configure -values [set $varname]]
}

;proc selected {w} {
    variable $w
    upvar 0 $w data

    if {$data(-command) != {}} {
	eval $data(-command) [list $w [$w get]]
    }
}
}; # end namespace eval ::Widget::Xcombobox

}; # end if ttk::combobox check

## #TEST
## xcombobox .e -label "Foo bar"
## .e configure -valuesvariable lll -default "c c"
## 
## pack .e
## pack [button .f -text "Press me" -command {puts :[.e get]:; set lll {1 2 3}}]
## 
## set lll [list a b c d e f g h {a a} {b b} {c c} {d d} {e e} {f f}] 
