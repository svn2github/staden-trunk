#!/bin/sh
#\
exec gap4sh -f "$0" ${@+"$@"} || exit 1

load_package gap

proc open_database {argv} {
    foreach {dummy dbname dbvers} \
	[regexp -inline {(.*)\.(.*)} [lindex $argv 0]] break
    return [open_db -name $dbname -version $dbvers -access rw]
}

load_package gap
set io [open_database $argv]
set db [io_read_database $io]
set nr [keylget db num_readings]

set chem_array(p) 2
set chem_array(t) 3
set chem_array(e) 8
set chem_array(f) 9
set chem_array(m) 12
set chem_array(n) 13
set chem_array(d) 5
set chem_array(b) 6
set chem_array(c) 7
set chem_array(l) 11
set chem_array(k) 17

for {set rnum 1} {$rnum <= $nr} {incr rnum} {
    set r [io_read_reading $io $rnum]
    set name [io_read_text $io [keylget r name]]
    set chem -1

    if {[regexp {^[0-9A-Z]{14}(_(left|right))??(\.(to|fm|pr|[0-9]+-)[0-9]+)*?$} $name] == 1} {
	# 454
	set chem 23
    } elseif {[regexp {^IL[0-9]+_[0-9]+:} $name] != 0} {
	# Solexa
	set chem 19
	keylset r chemistry 19
    } elseif {[regexp {\...(.).*$} $name _dummy chem] != 0} {
	# Capillary
	if {[info exists chem_array($chem)]} {
	    # puts $name=$chem_array($chem)
	    set chem $chem_array($chem)
	}
    }
    
    if {$chem != -1} {
	keylset r chemistry $chem
	io_write_reading $io $rnum $r
    } else {
        puts "ERROR: malformed name '$name'"
        continue
    }
}
close_db -io $io
exit
