#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

proc set_menu {name} {
    global menu_name $name
    set menu_name $name
    if ![info exists $name] {
	global $name
	set $name ""
    }
}

# 'pos' is now ignored.
proc add_menu {name onval offval {pos {}} {class {}}} {
    global menu_name
    global $menu_name
    regsub -all { } $name _ n
    set n [string tolower $n]
    regsub {.*\.} $name {} name
    keylset $menu_name $n.M [list $onval $offval $name $class]
}

proc add_cascade {name onval offval {class {}}} {
    global menu_name
    global $menu_name
    regsub -all { } $name _ n
    set n [string tolower $n]
    regsub {.*\.} $name {} name
    keylset $menu_name $n.m [list $onval $offval $name $class]
}

proc add_command {name onval offval comm {class {}}} {
    global menu_name
    global $menu_name
    regsub -all { } $name _ n
    set n [string tolower $n]
    regsub {.*\.} $name {} name
    keylset $menu_name $n.C [list $onval $offval $name $comm $class]
}

proc add_separator {name {class {}}} {
    global menu_name
    global $menu_name
    regsub -all { } $name _ n
    set n [string tolower $n]
    keylset $menu_name $n.S [list $class]
}

proc add_radio {name onval offval variable value comm {class {}}} {
    global menu_name
    global $menu_name
    regsub -all { } $name _ n
    set n [string tolower $n]
    regsub {.*\.} $name {} name
    keylset $menu_name $n.r [list $onval $offval $name $variable $value $comm \
	    $class]
}

proc add_check {name onval offval variable comm {class {}}} {
    global menu_name
    global $menu_name
    regsub -all { } $name _ n
    set n [string tolower $n]
    regsub {.*\.} $name {} name
    keylset $menu_name $n.c [list $onval $offval $name $variable $comm $class]
}

proc add_label {name} {
    global menu_name
    global $menu_name
    regsub -all { } $name _ n
    set n [string tolower $n]
    regsub {.*\.} $name {} name
    keylset $menu_name $n.L [list $name]
}

proc del_menu_item {name} {
    global menu_name
    global $menu_name

    regsub -all { } $name _ n
    set n [string tolower $n]
    regsub {.*\.} $name {} name

    keyldel $menu_name $n
}

proc menu_path {name} {
    regsub -all { } $name _ n
    set n [string tolower $n]
    return $n
}

#proc menu_path {name} {
#    regsub -all {[ 	]} $name _ name
#    regsub -all {\.} $name ._ name
#    return $name.m
#}

# Create a menu from the 'menu' specification. All the menu buttons are created
# and packed in 'name'. This is probably the pathname of a frame - eg .menubar.
# The mn and level parameters are used internally for recursion and shouldn't
# be specified.

# Create a menu from the 'menu' specification. All the top level menus are
# created as children of 'name'. This menu name needs to be given as the
# "-menu" option to a toplevel window.
# The mn and level parameters are used internally for recursion and shouldn't
# be specified.
proc create_menus {menu {mn {}} {class {}} {level {2}}} {
    global tcl_platform

    if {$tcl_platform(os) == "Darwin"} {
	set tearoffs 0
    } else {
	set tearoffs 1
    }
    set menu_state 0

    if {$class == ""} {
	set use_class 0
    } else {
	set use_class 1
    }

    foreach k [keylkeys menu] {
	switch [lindex [set j [lindex [set i [keylget menu $k]] 0]] 0] {
	    M {
		set j [lindex $j 1]
		if {!($use_class && [lindex $j 3] != "" && \
			[lsearch -exact [lindex $j 3] $class] == -1)} {
		    $mn add cascade -label [lindex $j 2] -menu $mn.$k
		    menu $mn.$k -title [lindex $j 2] -tearoff $tearoffs
		    set s [create_menus "[lrange $i 1 end]" $mn.$k $class]
		    set s [expr ([lindex $j 0]&1)|$s]
		    $mn entryconfigure [lindex $j 2] \
			    -state [lindex {disabled normal} $s]
		    global $mn.MS
		    set $mn.MS([lindex $j 2]) $s
		}
	    }
	    m {
		set j [lindex $j 1]
		if {!($use_class && [lindex $j 3] != "" && \
			[lsearch -exact [lindex $j 3] $class] == -1)} {
		    set s [expr [lindex $j 0]&1]
		    global $mn.MS
		    set $mn.MS([lindex $j 2]) $s
		    $mn add cascade \
			    -label [lindex $j 2] \
			    -menu $mn.$k \
			    -state [lindex {disabled normal} $s]
		    menu $mn.$k -tearoff $tearoffs
		    create_menus "[lrange $i 1 end]" $mn.$k $class \
			    [expr $level+1]
		}
	    }
	    C {
		set j [lindex $j 1]
		if {!($use_class && [lindex $j 4] != "" && \
			[lsearch -exact [lindex $j 4] $class] == -1)} {
		    regsub -all {"} [lindex $j 3] {\"} com
		    set s [expr [lindex $j 0]&1]
		    global $mn.MS
		    set $mn.MS([lindex $j 2]) $s
		    set c "$mn add command \
			    -label \"[lindex $j 2]\" \
			    -command \"eval [list $com]\" \
			    -state [lindex {disabled normal} $s]"
		    uplevel $level "$c"
		}
	    }
	    L {
		set j [lindex $j 1]
    		regsub -all {"} [lindex $j 3] {\"} com
		set c "$mn add command \
		    -label \"[lindex $j 0]\" \
		    -state disabled"
		uplevel $level "$c"
	    }
	    r {
		set j [lindex $j 1]
		if {!($use_class && [lindex $j 6] != "" && \
			[lsearch -exact [lindex $j 6] $class] == -1)} {
		    regsub -all {"} [lindex $j 5] {\"} com
		    set s [expr [lindex $j 0]&1]
		    global $mn.MS
		    set $mn.MS([lindex $j 2]) $s
		    set c "$mn add radiobutton \
			    -label  \"[lindex $j 2]\" \
			    -command \"eval [list $com]\" \
			    -state [lindex {disabled normal} $s] \
			    -variable [lindex $j 3] \
			    -value [lindex $j 4]"
		    uplevel $level "$c"
		}
	    }
	    c {
		set j [lindex $j 1]
		if {!($use_class && [lindex $j 5] != "" && \
			[lsearch -exact [lindex $j 5] $class] == -1)} {
		    regsub -all {"} [lindex $j 4] {\"} com
		    set s [expr [lindex $j 0]&1]
		    global $mn.MS
		    set $mn.MS([lindex $j 2]) $s
		    set c "$mn add checkbutton \
			    -label  \"[lindex $j 2]\" \
			    -command \"eval [list $com]\" \
			    -state [lindex {disabled normal} $s] \
			    -variable [lindex $j 3]"
		    uplevel $level "$c"
		}
	    }
	    S {
		set j [lindex $j 1]
		if {!($use_class && [lindex $j 0] != "" && \
			[lsearch -exact [lindex $j 0] $class] == -1)} {
		    $mn add separator
		}
	    }
	}
    }

    return $menu_state
}

proc merge_menus {menu {mn {}} {class {}} {level {3}}} {
    global tcl_platform menulist

    if {$tcl_platform(os) == "Darwin"} {
	set tearoffs 0
    } else {
	set tearoffs 1
    }
    set menu_state 0

    if {$class == ""} {
	set use_class 0
    } else {
	set use_class 1
    }

    lappend menulist $mn.MS

    foreach k [keylkeys menu] {
	switch [lindex [set j [lindex [set i [keylget menu $k]] 0]] 0] {
	    M {
		set j [lindex $j 1]
		if {!($use_class && [lindex $j 3] != "" && \
			[lsearch -exact [lindex $j 3] $class] == -1)} {

		    if {![winfo exists $mn.$k]} {
			$mn add cascade -label [lindex $j 2] -menu $mn.$k
			menu $mn.$k -title [lindex $j 2] -tearoff $tearoffs
		    }
		    set s [merge_menus "[lrange $i 1 end]" $mn.$k $class]
		    set s [expr ([lindex $j 0]&1)|$s]
		    $mn entryconfigure [lindex $j 2] \
			    -state [lindex {disabled normal} $s]
		    global $mn.MS
		    set $mn.MS([lindex $j 2]) $s
		}
	    }
	    m {
		set j [lindex $j 1]
		if {!($use_class && [lindex $j 3] != "" && \
			[lsearch -exact [lindex $j 3] $class] == -1)} {
		    set s [expr [lindex $j 0]&1]
		    global $mn.MS
		    if {![info exists $mn.MS([lindex $j 2])]} {
			set $mn.MS([lindex $j 2]) $s
			$mn add cascade \
			    -label [lindex $j 2] \
			    -menu $mn.$k \
			    -state [lindex {disabled normal} $s]
			menu $mn.$k -tearoff $tearoffs
			merge_menus "[lrange $i 1 end]" $mn.$k $class \
				[expr $level+1]
		    }
		}
	    }
	    C {
		set j [lindex $j 1]
		if {!($use_class && [lindex $j 4] != "" && \
			[lsearch -exact [lindex $j 4] $class] == -1)} {
		    regsub -all {"} [lindex $j 3] {\"} com
		    set s [expr [lindex $j 0]&1]
		    global $mn.MS
		    if {![info exists $mn.MS([lindex $j 2])]} {
			set $mn.MS([lindex $j 2]) $s
			set c "$mn add command \
				-label \"[lindex $j 2]\" \
				-command \"eval [list $com]\" \
				-state [lindex {disabled normal} $s]"
			uplevel $level "$c"
		    }
		}
	    }
	    L {
		set j [lindex $j 1]
    		regsub -all {"} [lindex $j 3] {\"} com
		set c "$mn add command \
		    -label \"[lindex $j 0]\" \
		    -state disabled"
		uplevel $level "$c"
	    }
	    r {
		set j [lindex $j 1]
		if {!($use_class && [lindex $j 6] != "" && \
			[lsearch -exact [lindex $j 6] $class] == -1)} {
		    regsub -all {"} [lindex $j 5] {\"} com
		    set s [expr [lindex $j 0]&1]
		    global $mn.MS
		    if {![info exists $mn.MS([lindex $j 2])]} {
			set $mn.MS([lindex $j 2]) $s
			set c "$mn add radiobutton \
				-label  \"[lindex $j 2]\" \
				-command \"eval [list $com]\" \
				-state [lindex {disabled normal} $s] \
				-variable [lindex $j 3] \
				-value [lindex $j 4]"
			uplevel $level "$c"
		    }
		}
	    }
	    c {
		set j [lindex $j 1]
		if {!($use_class && [lindex $j 5] != "" && \
			[lsearch -exact [lindex $j 5] $class] == -1)} {
		    regsub -all {"} [lindex $j 4] {\"} com
		    set s [expr [lindex $j 0]&1]
		    global $mn.MS
		    if {![info exists $mn.MS([lindex $j 2])]} {
			set $mn.MS([lindex $j 2]) $s
			set c "$mn add checkbutton \
				-label  \"[lindex $j 2]\" \
				-command \"eval [list $com]\" \
				-state [lindex {disabled normal} $s] \
				-variable [lindex $j 3]"
			uplevel $level "$c"
		    }
		}
	    }
	    S {
		set j [lindex $j 1]
		if {!($use_class && [lindex $j 0] != "" && \
			[lsearch -exact [lindex $j 0] $class] == -1)} {
		    $mn add separator
		}
	    }
	}
    }

    return $menu_state
}

proc menu_state_on {menu b {mn {}}} {
    global busy_mode
    set menu_state 0

    foreach k [keylkeys menu] {
	switch -glob [lindex [set j [lindex [set i [keylget menu $k]] 0]] 0] {
	    M {
		set j [lindex $j 1]
		set s [menu_state_on "[lrange $i 1 end]" $b $mn.$k]
		if {[expr (([lindex $j 0]&$b)!=0)|$s]} {
		    catch {$mn entryconfigure [lindex $j 2] -state normal}
		    global $mn.MS
		    set $mn.MS([lindex $j 2]) 1
		}
	    }
	    m {
		set j [lindex $j 1]
		if {[expr ([lindex $j 0]&$b)!=0]} {
		    set menu_state 1
		    catch {$mn entryconfigure [lindex $j 2] -state normal}
		    global $mn.MS
		    set $mn.MS([lindex $j 2]) 1
		    menu_state_on "[lrange $i 1 end]" $b $mn.$k
		}
	    }
	    {[Ccr]} {
		set j [lindex $j 1]
		if {[expr ([lindex $j 0]&$b)!=0]} {
		    set menu_state 1
		    catch {$mn entryconfigure [lindex $j 2] -state normal}
		    global $mn.MS
		    set $mn.MS([lindex $j 2]) 1
		}
	    }
	}
    }
    return $menu_state
}

proc menu_state_off {menu b {mn {}}} {
    set menu_state 0

    foreach k [keylkeys menu] {
	switch -glob [lindex [set j [lindex [set i [keylget menu $k]] 0]] 0] {
	    M {
		set j [lindex $j 1]
		set s [menu_state_off "[lrange $i 1 end]" $b $mn.$k]
		if {$s==0 && [expr (([lindex $j 1]&$b)!=0)|$s]} {
		    catch {$mn entryconfigure [lindex $j 2] -state disabled}
		    global $mn.MS
		    set $mn.MS([lindex $j 2]) 0
		}
	    }
	    m {
		set j [lindex $j 1]
		if {[expr ([lindex $j 1]&$b)!=0]} {
		    catch {$mn entryconfigure [lindex $j 2] -state disabled}
		    global $mn.MS
		    set $mn.MS([lindex $j 2]) 0
		    menu_state_off "[lrange $i 1 end]" $b $mn.$k
		} else {
		    set menu_state 1
		}
	    }
	    {[Ccr]} {
		set j [lindex $j 1]
		if {[expr ([lindex $j 1]&$b)!=0]} {
		    catch {$mn entryconfigure [lindex $j 2] -state disabled}
		    global $mn.MS
		    set $mn.MS([lindex $j 2]) 0
		} else {
		    set menu_state 1
		}
	    }
	}
    }
    return $menu_state
}

proc menu_state_save {menu menu_name {state {}}} {
    global $menu.MS $menu.LIST ${menu_name}_list
    if {$menu_name != "" && [info exists ${menu_name}_list]} {
	global ${menu_name}_list
	set $menu.LIST [set ${menu_name}_list]
    }

    if {[info exists $menu.MS]} {
        lappend state [array get $menu.MS]
    }

    foreach w [winfo children $menu] {
	if {[winfo children $w] != {}} {
	    set state [menu_state_save $w {} $state]
	} else {
	    global $w.MS
	    if {[info exists $w.MS]} {
		lappend state [array get $w.MS]
	    }
        }
    }

    return $state
}

proc menu_state_restore {menu menu_name state {first 1}} {
    global $menu.MS $menu.LIST ${menu_name}_list

    if {$menu_name != "" && [info exists ${menu_name}_list]} {
	global ${menu_name}_list
	catch {set ${menu_name}_list [set $menu.LIST]}
    }


    if {$first == 1 && [info exists $menu.MS]} {
	array set s [lindex $state 0]
	foreach i [array names s] {
	    catch {$menu entryconfigure $i \
		    -state [lindex {disabled normal} $s($i)]}
	    global $menu.MS
	    set $menu.MS($i) $s($i)
	}
	set state [lrange $state 1 end]
    }

    foreach w [winfo children $menu] {
	if [info exists s] {unset s}
	array set s [lindex $state 0]
	foreach i [array names s] {
	    catch {$w entryconfigure $i \
		    -state [lindex {disabled normal} $s($i)]}
	    global $w.MS
	    set $w.MS($i) $s($i)
	}
	set state [lrange $state 1 end]
	if {[winfo children $w] != {}} {
	    set state [menu_state_restore $w {} $state 0]
	}
    }

    return $state
}

# Set the menu enable/disable information to the state specified in 'bit'.
# Note that 'menu_name' here is the _name_ of the global menu variable rather
# than the contents. This is so that we can store extra information in the
# name pertaining to the current state.
#
# The 'bit' (state) value should be a power of 2 rather than the bit number.
#
# The 'prefix' parameter is a prefix to add to the pathname (eg it's the frame
# pathname that the menus are packed in).
proc menu_state_set {menu_name bit {prefix {}}} {
    global $menu_name
    global ${menu_name}_list

    lappend ${menu_name}_list $bit

    if {$bit > 0} {
        menu_state_on [set $menu_name] $bit $prefix
    } else {
        menu_state_off [set $menu_name] [expr -$bit] $prefix
    }
}

#
# Sets the menu states to their current settings.
# This may sound useless, but it's needed when the menus change so that we
# can add any new menu items with the correct menu state.
#
proc reset_menu_state {menu_name {prefix {}}} {
    global $menu_name
    global ${menu_name}_list

    if {[info exists ${menu_name}_list]} {
	set tmp_list [set ${menu_name}_list]
	foreach bit [set ${menu_name}_list] {
	    menu_state_set $menu_name $bit $prefix
	}
	set ${menu_name}_list $tmp_list
    }
}