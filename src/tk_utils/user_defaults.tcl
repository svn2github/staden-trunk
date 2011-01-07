#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
source $env(STADTABL)/shlib.conf

set defn ""
set packages() ""
set package_dir() ""

# Load a package
# tcldir is the directory containing the tcl sources to add to auto_path
# libdir is the location of a dynamically loadable module to link with
# name is the stub for the initialisation routine of libdir
# The presence of urther 'args' implies to ignore errors when loading.
proc load_package {tcldir {libdir {}} {name {}} {init {}}} {
    global packages package_dir env

    if {$libdir == {} && $name == {} && $init == {}} {
	simple_load_package $tcldir
	return
    }

    global auto_path lib_prefix lib_suffix packages env

    #puts "Loading package $name"
    if {[info exists packages($name)]} {
	if {$packages($name)==3} {
	    return
	}
    } else {
	set packages($name) 0
    }

    # Replace %L by $env(STADLIB), %S with src, %% with %
    regsub -all {%L} $tcldir $env(STADLIB) tcldir
    regsub -all {%L} $libdir $env(STADLIB) libdir
    regsub -all {%T} $tcldir $env(STADTCL) tcldir
    regsub -all {%T} $libdir $env(STADTCL) libdir
    regsub -all {%S} $tcldir $env(STADENROOT)/src tcldir
    regsub -all {%S} $libdir $env(STADENROOT)/src libdir
    regsub -all {%%} $tcldir %% tcldir
    regsub -all {%%} $libdir %% libdir

    set package_dir($name) $tcldir

    if {$tcldir != ""} {
        lappend auto_path $tcldir
	if {!($packages($name) & 1)} {
	    set packages($name) [expr $packages($name) | 1]
	    SetDefaults $env(STADTABL) ${name}rc ${name}_defs $tcldir/${name}rc
#	    catch {eval ${name}_init}
	}
    }

    if {!($packages($name) & 2) && $libdir != "-"} {
	if {$libdir != ""} {
	    set libdir "$libdir/"
	}
	if {[info exists env(STADEN_DEBUG)]} {
	    puts -nonewline "load ${libdir}${lib_prefix}${name}${lib_suffix}"
	    flush stdout
	}
	set err ""
        if {$init == 1} {
	    if {[file exists ${libdir}${lib_prefix}${name}${lib_suffix}]} {
		load ${libdir}${lib_prefix}${name}${lib_suffix} $name
	    }
        } else {
	    if {[file exists ${libdir}${lib_prefix}${name}${lib_suffix}]} {
		catch {load ${libdir}${lib_prefix}${name}${lib_suffix}} err
		set err "($err)"
	    }
        }
	if {[info exists env(STADEN_DEBUG)]} {
	    puts " => $err"
	}

	global ${name}_defs
        if {[info exists ${name}_defs] && [set ${name}_defs] == "NULL"} {
	    set ${name}_defs ""
	}

	set packages($name) [expr $packages($name) | 2]
    }

    # puts "Finished loading $name"
}

#
# A simpler form of load_package; called when load_package is only given
# one argument (the package name).
#
# We always load the rc file, which we assume to be called ${name}rc and
# to be held within the $STADTABL directory.
#
# Tcl code is assumed to be within $STADLIB/name, unless a ${name}_tcldir
# variable has previously been set.
#
# A dynamic library of (eg) $STADLIB/${lib_prefix}${name}${lib_suffix} is
# loaded, if it exists.
#
proc simple_load_package {name} {
    global packages env ${name}_tcldir lib_prefix lib_suffix auto_path
    global package_dir tcl_platform tcl_version
    # puts "start simple_load_package $name"

    regsub -all {%L} $name $env(STADLIB) name
    regsub -all {%T} $name $env(STADTCL) name
    regsub -all {%S} $name $env(STADENROOT)/src name
    regsub -all {%%} $name %% name

    if {[file pathtype $name] != "relative"} {
	set dir [file dirname $name]
	set name [file tail $name]
	set tdir $dir
	set ldir $dir
	set adir $dir
	set bdir $ldir
    } else {
	set tdir $env(STADTABL)
	set ldir $env(STADLIB)
	set adir $env(STADTCL)/$name
	if {$tcl_platform(os) == "Darwin" || $tcl_version < "8.3"} {
	    set bdir "$ldir"
	} else {
	    set bdir ""
	}
    }

    set package_dir($name) $ldir
    lappend auto_path $adir
    if {![info exists packages($name)]} {
	set packages($name) 0
    }
    if {!($packages($name) & 1)} {
	set packages($name) [expr $packages($name) | 1]
	SetDefaults $tdir ${name}rc ${name}_defs
    }

    if {$bdir != ""} {
        set bdir "$bdir/"
    }
    set fname $bdir${lib_prefix}${name}${lib_suffix}

    #read local libraries
    if {[info exists env(STADEN_LIB_PATH)]} {
	foreach path "$env(STADEN_LIB_PATH)" {
	    set tmp_name $path/[file tail $fname]
	    if {[file exists $tmp_name]} {
		set fname $tmp_name
	    }
	}
    }

    if {!($packages($name) & 2)} {
	if {[info exists env(STADEN_DEBUG)]} {
	    puts -nonewline "load $fname"
	    flush stdout
	}
	set err ""
	if {($tcl_platform(os) == "Darwin" && [file exists $fname]) || \
	    $tcl_platform(os) != "Darwin"} {
	    catch {load $fname} err
	}
	if {[info exists env(STADEN_DEBUG)]} {
	    puts " => $err"
	} elseif {[regexp {.*[Uu]n(resolv|defin).*} $err]} {
	    puts "load $fname => $err"
	}
	set packages($name) [expr $packages($name) | 2]

	global ${name}_defs
        if {[info exists ${name}_defs] && [set ${name}_defs] == "NULL"} {
	    set ${name}_defs ""
	}
    }

    # puts "end simple_load_package $name"
}

proc set_defs {name} {
    global defn
    set defn $name
}

proc set_def {a1 a2} {
    global defn
    global $defn
    keylset $defn $a1 $a2
}

proc set_defx {a1 a2 a3} {
    uplevel 1 global $a1\; keylset $a1 $a2 \"$a3\"
}

proc SetDefaults {dirname filename defname args} {
    # puts "setdefaults $defname start"
    global packages env $defname defn
    global normal_bg root_dir
    set normal_bg "#d9d9d9"

    set root_dir $dirname
    set tmp $defn
    set defn $defname
    global $defname
    if {![info exists $defname] || [set $defname] == "NULL"} {
	set $defname ""
    }
    if ![info exists packages] {
	set packages ""
    }

    # 11/1/99 johnt - convert WINNT \ to proper /
    regsub -all {\\} $env(HOME) / env(HOME)

    # We need to do this out of order. The sourcing of the rc file in STADTABL
    # is dependent on some user's values, so we parse their rc file first
    # to look for MENU_FILE lines.
    if {$filename != ""} {
	foreach i [list $env(HOME)/.$filename .$filename] {
	    if {[file exists $i]} {
		set fd [open $i r]
		foreach line [split [read $fd] "\n"] {
		    if {[regexp "^\[ \t\]*set_def\[ \t\]+MENU_FILE\[ \t\]+(\[^ \t\]+)" \
			$line dummy level]} {
			keylset $defname MENU_FILE $level
		    }
		}
		close $fd
	    }
	}
    }

    # Load RC file
    if {$filename != ""} {
        foreach i [list $args $dirname/$filename $env(HOME)/.$filename .$filename] {
	    if {$i != "" && [file exists $i]} {
	        uplevel #0 [list source $i]
            }
	    if {$i != "" && [file exists $i.local]} {
	        uplevel #0 [list source $i.local]
            }
        }
    }

    set defn $tmp
}

proc LoadLibs {name} {
    global packages env

    set packages($name) 1
    SetDefaults $env(STADTABL) ${name}rc ${name}_defs
}

proc load_menus {} {
    global defn root_dir licence
    global $defn

    regsub {_defs$} $defn {} tmp

    if {[catch {keylget $defn MENU_FILE}]} {
	switch $licence(type) {
	    f {set_def MENU_FILE full}
	    v {set_def MENU_FILE viewer}
	    default {set_def MENU_FILE demo}
	}
    }
    set f $root_dir/${tmp}rc_menu_[keylget $defn MENU_FILE]
    if {[file exists $f]} {
	uplevel #0 [list source $f]
    }

#    set body [info body add_command]
#    append body {
#	global ${menu_name}_extras
#	keylset ${menu_name}_extras $n.C [list $onval $offval $name $comm]
#    }
#    proc add_command [info args add_command] $body
}
