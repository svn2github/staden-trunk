#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1998. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#
# Code for handling a naming scheme
# The naming scheme is specified as follows:
#
# set ns_regexp <regular_expression>
# set ns_lt(<line_type>) <subst expression>
# (etc)
# The format of subst expression is either a literal string (eg $1 or _${1}_)
# or "subst" followed by a string to match, a tcl list of "string match"
# patterns and replacement strings, and a default value for when no "string
# match" patterns are matched.

# Eg:
#	set ns_regexp {([^.]*)\.(..)(.).*}
#	set ns_lt(TN) {$1}
#	set ns_lt(PR) {subst {$2 {[sp]1 1} {[qr]1 2} {[sp]* 3} {[qr]* 4} 0}}
#	set ns_lt(CH) {subst {$3 {[tdc] 1} 0}}
#
# This expands PR_com into:
#
#	proc PR_com {} {
#	        global lines
#	        regexp {([^.]*)\.(..)(.).*} $lines(ID) matched 1 2 3
#	        if {[string match {[sp]1} $2]} {return 1}
#	        if {[string match {[qr]1} $2]} {return 2}
#	        if {[string match {[sp]*} $2]} {return 3}
#	        if {[string match {[qr]*} $2]} {return 4}
#	        return 0
#	}

proc set_name_scheme {} {
    global ns_regexp ns_lt
    if {![info exists ns_regexp]} {
	return
    }

    # Count number of open brackets
    regsub -all {[^(]*} $ns_regexp {} tmp
    set brackets [string length $tmp]

    set body_all "\tglobal lines\n\tif {\[regexp [list $ns_regexp] \$lines(ID) matched"
    for {set i 1} {$i <= $brackets} {incr i} {
	append body_all " $i"
    }
    append body_all "\] == 0} {\n\t\treturn \"\"\n\t}\n"

    foreach lt [array names ns_lt] {
	set body "$body_all"
	if {[regexp {^subst (.*)} $ns_lt($lt) dummy sub]} {
	    set var [lindex [lindex $sub 0] 0]
	    foreach pat [lrange [lindex $sub 0] 1 end] {
		if {[llength $pat] == 1} {
		    append body "\treturn $pat\n"
		} else {
		    set pat_l [lindex $pat 0]
		    set pat_r [lindex $pat 1]
		    append body "\tif {\[string match [list $pat_l] $var\]} {return $pat_r}\n"
		}
	    }
	} else {
	    append body "\treturn $ns_lt($lt)\n"
	}

	eval proc ${lt}_com {{}} [list $body]
    }
}
