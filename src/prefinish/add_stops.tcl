#!/bin/sh
#\
exec stash -f "$0" ${@+"$@"} || exit 1

set stops stops
set opts "-w 100 -t 4"

proc open_database {argv} {
    foreach {dbname dbvers} [split [lindex $argv 0] .] {}
    return [open_db -name $dbname -version $dbvers]
}

load_package gap
set io [open_database $argv]
set db [io_read_database $io]
set nr [keylget db num_readings]

for {set rnum 1} {$rnum <= $nr} {incr rnum} {
    set r [io_read_reading $io $rnum]
    set tname [io_read_text $io [keylget r trace_name]]

    puts $rnum:$tname...
    if {[catch {set out [eval exec [list $stops] $opts [list $tname]]} err]} {
	puts stderr "$stops: $err"
	continue
    }

    if {$out == ""} {
	continue
    }

    set len [keylget r length]
    set opos [keylget r orig_positions]
    if {$opos == 0} {
	puts "No original position data; skipping"
    }
    set opos [io_read_data $io $opos [expr {$len*2}] 2]
    binary scan $opos s* opos

    foreach line [split $out "\n"] {
	scan $line "%s Peak %d at %d / %d height %f" _ _ _ pos score
	set lpos [lsearch $opos $pos]
	if {$lpos == 0} {
	    puts "Couldn't find position in opos array"
	    continue
	}
	if {[keylget r sense] == 1} {
	    set tpos [expr {[llength $opos]-$lpos}]
	} else {
	    set tpos [expr {$lpos+1}]
	}
	# add_tags -io $io -tags "{$rnum STOP = $tpos..$tpos score = $score}"
	puts $pos,$score,$tpos
    }
}

close_db -io $io
exit

