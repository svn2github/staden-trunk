#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#==============================================================================
# Creates a new version of the file "tablelistWidget.tcl" by defining the
# procedure "arrElemExists" and replacing all invocations of "[info exists
# <array>(<name>)]" with "[arrElemExists <array> <name>]".  This works around a
# bug in Tcl versions 8.2, 8.3.0 - 8.3.2, and 8.4a1 (fixed in Tcl 8.3.3 and
# 8.4a2), which causes excessive memory use when calling "info exists" on non-
# existent array elements.
#
# Copyright (c) 2001  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

file copy tablelistWidget.tcl tablelistWidget.tcl.BAK

set fi [open tablelistWidget.tcl.BAK r]
set fo [open tablelistWidget.tcl     w]

set procDef {
    #
    # The following procedure returns 1 if arrName($name) exists and
    # 0 otherwise.  It is a (partial) replacement for [info exists
    # arrName($name)], which -- due to a bug in Tcl versions 8.2,
    # 8.3.0 - 8.3.2, and 8.4a1 (fixed in Tcl 8.3.3 and 8.4a2) --
    # causes excessive memory use if arrName($name) doesn't exist.
    # The first version of the procedure assumes that the second
    # argument doesn't contain glob-style special characters.
    #
    if {[regexp {^8\.(2\.[0-3]|3\.[0-2]|4a1)$} $tk_patchLevel]} {
	proc arrElemExists {arrName name} {
	    upvar $arrName arr
	    return [llength [array names arr $name]]
	}
    } else {
	proc arrElemExists {arrName name} {
	    upvar $arrName arr
	    return [info exists arr($name)]		;# this is much faster
	}
    }
}

for {set n 1} {[gets $fi line] >= 0} {incr n} {
    if {$n == 21} {
	puts -nonewline $fo $line
	puts $fo $procDef
    } else {
	regsub -all {\[info exists ([^\(]+)\(([^\]]+)\)\]} $line \
		    {[arrElemExists \1 \2]} line
	puts $fo $line
    }
}

puts "Made backup copy \"tablelistWidget.tcl.BAK\"."
puts "Created new version of the file \"tablelistWidget.tcl\"."
