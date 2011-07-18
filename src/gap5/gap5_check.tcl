#!/bin/sh

# A command line tool to run the export consensus sequences code of Gap5.
#
# See the usage function below or run the program with no output for
# usage information.

#\
exec tclsh $0 ${@+"$@"}

#-----------------------------------------------------------------------------
# Startup code
source $env(STADTABL)/shlib.conf
load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}
load_package tk_utils
tk_utils_init
load_package gap5

set fix 0
set acc ro
if {[lindex $argv 0] == "-f"} {
    set fix 1
    set acc rw
    set argv [lrange $argv 1 end]
}

set db [lindex $argv 0]

if {[catch {set io [g5::open_database -name $db -access $acc]} err]} {
    puts stderr "Couldn't open database '$db': $err"
    exit 1
}

set err 0

if {[llength $argv] == 1} {
    puts "=== checking entire DB ==="
    $io check $fix 2
} else {
    foreach crec [lrange $argv 1 end] {
	set c [$io get_contig $crec]
	puts "=== checking contig #$crec ==="
	incr err [$c check $fix]
	$c delete
    }
}

puts "Total errors: $err"

# Avoid the closing output spam if possible
if {$acc == "rw"} {
    puts "Saving changes"
    $io flush
    $io close
}

exit [expr {$err == 0 ? 0 : 1}]
