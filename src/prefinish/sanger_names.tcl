#!/bin/sh
#\
exec gap4sh -f "$0" ${@+"$@"} || exit 1

proc open_database {argv} {
    foreach {dbname dbvers} [split [lindex $argv 0] .] {}
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
    if {[regexp {\...(.).*$} $name _dummy chem] == 0} {
        puts "ERROR: malformed name '$name'"
        continue
    }
    if {[info exists chem_array($chem)]} {
	# puts $name=$chem_array($chem)
	keylset r chemistry $chem_array($chem)
	io_write_reading $io $rnum $r
    }
}
close_db -io $io
exit
