#==============================================================================
# Contains utility procedures for mega-widgets.
#
# Structure of the module:
#   - Namespace initialization
#   - Public utility procedures
#   - Private procedures used in bindings
#
# Copyright (c) 2000-2001  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
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
    variable version	1.3
    variable library	[file dirname [info script]]

    #
    # Public procedures:
    #
    namespace export	wrongNumArgs configure fullConfigOpt fullOpt enumOpts \
			attribSubCmd configSubCmd

    #
    # Define the binding tag allModif to be a duplicate of all (having the
    # same events and scripts), but replace the script corresponding to the
    # events <Shift-Key-Tab> and <<PrevWindow>> with a new one which skips
    # the widget's parent when searching for the previous window in "focus
    # order" (<Shift-Key-Tab> was replaced with <<PrevWindow>> in Tk 8.3.0)
    # 
    foreach event [bind all] {
	if {[string compare $event <Shift-Key-Tab>] == 0 ||
	    [string compare $event <<PrevWindow>>] == 0} {
	    bind allModif $event [list mwutil::skipParent %W $event]
	} else {
	    bind allModif $event [bind all $event]
	}
    }
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
proc mwutil::wrongNumArgs msg {
    return -code error "wrong # args: should be \"$msg\""
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
    set failed no
    set count [llength $optValPairs]
    foreach {opt val} $optValPairs {
	if {[catch {fullConfigOpt $opt configSpecs} result] != 0} {
	    set failed yes
	    break
	}
	if {$count == 1} {
	    set result "value for \"$opt\" missing"
	    set failed yes
	    break
	}
	set opt $result
	lappend cmdLineOpts $opt
	lappend savedVals $configVals($opt)
	if {[catch {eval $configCmd [list $win $opt $val]} result] != 0} {
	    set failed yes
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
		eval $configCmd [list $win $opt $dbValue]
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
proc mwutil::fullOpt {name opt optList} {
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
		   "bad $name \"$opt\": must be [enumOpts $optList]"
	}

	1 {
	    return $option
	}

	default {
	    return -code error \
		   "ambiguous $name \"$opt\": must be [enumOpts $optList]"
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
	    return [configure $win configSpecs configVals $configCmd $argList no]
	}
    }
}

#
# Private procedures used in bindings
# ===================================
#

#------------------------------------------------------------------------------
# mwutil::skipParent
#
# This procedure handles <Shift-Key-Tab> and <<PrevWindow>> events in the child
# w of a mega-widget.  It generates the given event for the parent of the
# widget.
#------------------------------------------------------------------------------
proc mwutil::skipParent {w event} {
    set parent [winfo parent $w]
    focus $parent			;# necessary on Windows
    event generate $parent $event
}
