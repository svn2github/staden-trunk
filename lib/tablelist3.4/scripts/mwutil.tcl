#==============================================================================
# Contains utility procedures for mega-widgets.
#
# Structure of the module:
#   - Namespace initialization
#   - Public utility procedures
#
# Copyright (c) 2000-2004  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

package require Tcl 8
package require Tk  8

#
# Namespace initialization
# ========================
#

namespace eval mwutil {
    #
    # Public variables:
    #
    variable version	1.7
    variable library	[file dirname [info script]]

    #
    # Public procedures:
    #
    namespace export	wrongNumArgs defineKeyNav generateEvent \
			configure fullConfigOpt fullOpt enumOpts \
			setConfigVals configSubCmd attribSubCmd 
}

#
# Public utility procedures
# =========================
#

#------------------------------------------------------------------------------
# mwutil::wrongNumArgs
#
# Generates a "wrong # args" error message.
#------------------------------------------------------------------------------
proc mwutil::wrongNumArgs args {
    set optList {}
    foreach arg $args {
	lappend optList \"$arg\"
    }
    return -code error "wrong # args: should be [enumOpts $optList]"
}

#------------------------------------------------------------------------------
# mwutil::defineKeyNav
#
# For a given mega-widget class, the procedure defines the binding tag
# ${class}KeyNav as a partial replacement for "all", by substituting the
# scripts bound to the events <Key-Tab>, <Shift-Key-Tab>, and <<PrevWindow>>
# with new ones which propagate these events to the mega-widget of the given
# class containing the widget to which the event was reported.  (The event
# <Shift-Key-Tab> was replaced with <<PrevWindow>> in Tk 8.3.0.)  This tag is
# designed to be inserted before "all" in the list of binding tags of a
# descendant of a mega-widget of the specified class.
#------------------------------------------------------------------------------
proc mwutil::defineKeyNav class {
    foreach event {<Key-Tab> <Shift-Key-Tab> <<PrevWindow>>} {
	bind ${class}KeyNav $event \
	     [list mwutil::generateEvent %W $class $event]
    }
}

#------------------------------------------------------------------------------
# mwutil::generateEvent
#
# This procedure generates the given event for the mega-widget of the specified
# class containing the widget w if that mega-widget is not the only widget
# receiving the focus during keyboard traversal within its top-level widget.
#------------------------------------------------------------------------------
proc mwutil::generateEvent {w class event} {
    while {[string compare [winfo class $w] $class] != 0} {
	set w [winfo parent $w]
    }

    if {[string compare [tk_focusNext $w] $w] != 0} {
	focus $w				;# necessary on Windows
	event generate $w $event
    }

    return -code break ""
}

#------------------------------------------------------------------------------
# mwutil::configure
#
# Configures the widget win by processing the command-line arguments specified
# in optValPairs and, if the value of initialize is true, also those database
# options that don't match any command-line arguments.
#------------------------------------------------------------------------------
proc mwutil::configure {win configSpecsName configValsName \
			configCmd optValPairs initialize} {
    upvar $configSpecsName configSpecs
    upvar $configValsName configVals

    #
    # Process the command-line arguments
    #
    set cmdLineOpts {}
    set savedVals {}
    set failed 0
    set count [llength $optValPairs]
    foreach {opt val} $optValPairs {
	if {[catch {fullConfigOpt $opt configSpecs} result] != 0} {
	    set failed 1
	    break
	}
	if {$count == 1} {
	    set result "value for \"$opt\" missing"
	    set failed 1
	    break
	}
	set opt $result
	lappend cmdLineOpts $opt
	lappend savedVals $configVals($opt)
	if {[catch {eval $configCmd [list $win $opt $val]} result] != 0} {
	    set failed 1
	    break
	}
	incr count -2
    }

    if {$failed} {
	#
	# Restore the saved values
	#
	foreach opt $cmdLineOpts val $savedVals {
	    eval $configCmd [list $win $opt $val]
	}

	return -code error $result
    }

    if {$initialize} {
	#
	# Process those configuration options that were not
	# given as command-line arguments; use the corresponding
	# values from the option database if available
	#
	foreach opt [lsort [array names configSpecs]] {
	    if {[llength $configSpecs($opt)] == 1 ||
		[lsearch -exact $cmdLineOpts $opt] >= 0} {
		continue
	    }
	    set dbName [lindex $configSpecs($opt) 0]
	    set dbClass [lindex $configSpecs($opt) 1]
	    set dbValue [option get $win $dbName $dbClass]
	    if {[string compare $dbValue ""] != 0} {
		if {[catch {
		    eval $configCmd [list $win $opt $dbValue]
		} result] != 0} {
		    return -code error $result
		}
	    } else {
		set default [lindex $configSpecs($opt) 3]
		eval $configCmd [list $win $opt $default]
	    }
	}
    }

    return ""
}

#------------------------------------------------------------------------------
# mwutil::fullConfigOpt
#
# Returns the full configuration option corresponding to the possibly
# abbreviated option opt.
#------------------------------------------------------------------------------
proc mwutil::fullConfigOpt {opt configSpecsName} {
    upvar $configSpecsName configSpecs

    if {[info exists configSpecs($opt)]} {
	if {[llength $configSpecs($opt)] == 1} {
	    return $configSpecs($opt)
	} else {
	    return $opt
	}
    }

    set optList [lsort [array names configSpecs]]
    set count 0
    foreach elem $optList {
	if {[string first $opt $elem] == 0} {
	    incr count
	    if {$count == 1} {
		set option $elem
	    } else {
		break
	    }
	}
    }

    switch $count {
	0 {
	    ### return -code error "unknown option \"$opt\""
	    return -code error \
		   "bad option \"$opt\": must be [enumOpts $optList]"
	}

	1 {
	    if {[llength $configSpecs($option)] == 1} {
		return $configSpecs($option)
	    } else {
		return $option
	    }
	}

	default {
	    ### return -code error "unknown option \"$opt\""
	    return -code error \
		   "ambiguous option \"$opt\": must be [enumOpts $optList]"
	}
    }
}

#------------------------------------------------------------------------------
# mwutil::fullOpt
#
# Returns the full option corresponding to the possibly abbreviated option opt.
#------------------------------------------------------------------------------
proc mwutil::fullOpt {kind opt optList} {
    if {[lsearch -exact $optList $opt] >= 0} {
	return $opt
    }

    set count 0
    foreach elem $optList {
	if {[string first $opt $elem] == 0} {
	    incr count
	    if {$count == 1} {
		set option $elem
	    } else {
		break
	    }
	}
    }

    switch $count {
	0 {
	    return -code error \
		   "bad $kind \"$opt\": must be [enumOpts $optList]"
	}

	1 {
	    return $option
	}

	default {
	    return -code error \
		   "ambiguous $kind \"$opt\": must be [enumOpts $optList]"
	}
    }
}

#------------------------------------------------------------------------------
# mwutil::enumOpts
#
# Returns a string consisting of the elements of the given list, separated by
# commas and spaces.
#------------------------------------------------------------------------------
proc mwutil::enumOpts optList {
    set optCount [llength $optList]
    set n 1
    foreach opt $optList {
	if {$n == 1} {
	    set str $opt
	} elseif {$n < $optCount} {
	    append str ", $opt"
	} else {
	    if {$optCount > 2} {
		append str ","
	    }
	    append str " or $opt"
	}

	incr n
    }

    return $str
}

#------------------------------------------------------------------------------
# mwutil::setConfigVals
#
# Sets the elements of the array specified by configValsName to the values
# returned by passing the widget name win and the relevant options to the
# command given by cgetCmd.
#------------------------------------------------------------------------------
proc mwutil::setConfigVals {win configSpecsName configValsName
			    cgetCmd argList} {
    upvar $configSpecsName configSpecs
    upvar $configValsName configVals

    set optList {}
    if {[llength $argList] == 0} {
	foreach opt [array names configSpecs] {
	    if {[llength $configSpecs($opt)] > 1} {
		lappend optList $opt
	    }
	}
    } else {
	foreach {opt val} $argList {
	    lappend optList [fullConfigOpt $opt configSpecs]
	}
    }

    foreach opt $optList {
	set configVals($opt) [eval $cgetCmd [list $win $opt]]
    }
}

#------------------------------------------------------------------------------
# mwutil::configSubCmd
#
# This procedure is invoked to process configuration subcommands.
#------------------------------------------------------------------------------
proc mwutil::configSubCmd {win configSpecsName configValsName
			   configCmd argList} {
    upvar $configSpecsName configSpecs
    upvar $configValsName configVals

    switch [llength $argList] {
	0 {
	    #
	    # Return a list describing all available configuration options
	    #
	    foreach opt [lsort [array names configSpecs]] {
		if {[llength $configSpecs($opt)] == 1} {
		    set alias $configSpecs($opt)
		    if {$::tk_version < 8.1} {
			set dbName [lindex $configSpecs($alias) 0]
			lappend result [list $opt $dbName]
		    } else {
			lappend result [list $opt $alias]
		    }
		} else {
		    set dbName [lindex $configSpecs($opt) 0]
		    set dbClass [lindex $configSpecs($opt) 1]
		    set default [lindex $configSpecs($opt) 3]
		    lappend result [list $opt $dbName $dbClass $default \
				    $configVals($opt)]
		}
	    }
	    return $result
	}

	1 {
	    #
	    # Return the description of the specified configuration option
	    #
	    set opt [fullConfigOpt [lindex $argList 0] configSpecs]
	    set dbName [lindex $configSpecs($opt) 0]
	    set dbClass [lindex $configSpecs($opt) 1]
	    set default [lindex $configSpecs($opt) 3]
	    return [list $opt $dbName $dbClass $default $configVals($opt)]
	}

	default {
	    #
	    # Set the specified configuration options to the given values
	    #
	    return [configure $win configSpecs configVals $configCmd $argList 0]
	}
    }
}

#------------------------------------------------------------------------------
# mwutil::attribSubCmd
#
# This procedure is invoked to process the attrib subcommand.
#------------------------------------------------------------------------------
proc mwutil::attribSubCmd {win argList} {
    set classNs [string tolower [winfo class $win]]
    upvar ::${classNs}::ns${win}::attribVals attribVals

    set argCount [llength $argList]
    switch $argCount {
	0 {
	    #
	    # Return the current list of attribute names and values
	    #
	    set result {}
	    foreach attr [lsort [array names attribVals]] {
		lappend result [list $attr $attribVals($attr)]
	    }
	    return $result
	}

	1 {
	    #
	    # Return the value of the specified attribute
	    #
	    set attr [lindex $argList 0]
	    if {[info exists attribVals($attr)]} {
		return $attribVals($attr)
	    } else {
		return ""
	    }
	}

	default {
	    #
	    # Set the specified attributes to the given values
	    #
	    if {$argCount % 2 != 0} {
		return -code error "value for \"[lindex $argList end]\" missing"
	    }
	    array set attribVals $argList
	    return ""
	}
    }
}
