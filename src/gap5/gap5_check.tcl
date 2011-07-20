#!/bin/sh

# A command line tool to run the export consensus sequences code of Gap5.
#
# See the usage function below or run the program with no output for
# usage information.

#\
exec tclsh $0 ${@+"$@"}

#-----------------------------------------------------------------------------
proc main {args} {
    set fmt {
	h 0 0
	help 0 0
	contigs 1 ""
	level 1 2
	f 0 0
	fix 0 0}

    set args [eval parse_args [list $fmt] opt $args]

    if {[llength $args] != 1} {
	usage 1
    }

    if {$opt(h) || $opt(help)} {
	usage 0
    }

    if {$opt(f) || $opt(fix)} {
	set opt(fix) 1
	set acc rw
    } else {
	set acc ro
    }

    set io [open_db [lindex $args 0] $acc]

    set err 0

    if {$opt(contigs) == ""} {
	puts "=== checking entire DB ==="
	set err [$io check $opt(fix) $opt(level)]
    } else {
	foreach crec $opt(contigs) {
	    set c [$io get_contig $crec]
	    puts "=== checking contig #$crec ==="
	    incr err [$c check $opt(fix) $opt(level)]
	    $c delete
	}
    }

    if {$acc == "rw"} {
	puts "Saving changes"
	$io flush
	$io close
    }

    exit [expr {$err == 0 ? 0 : 1}]
}

proc usage {e} {
    puts ""
    puts {Usage:  gap5_command check [options] DBNAME.VERS}
    puts ""
    puts {Where [options] are any of the following:}
    puts "    -h|-help             Shows this help."
    puts "    -f|-fix              Attempts to fix the database"
    puts "                         *PLEASE BACK UP THE DB FIRST*"
    puts "    -contigs 'list'      Output check specific contig record numbers"
    puts "    -level 'num'         Set check level to '1' or '2'"
    puts ""
    exit $e
}

# Opens a gap5 db and returns the database handle
proc open_db {name {acc ro}} {
    if {[catch {set io [g5::open_database -name $name -access $acc]} err]} {
	puts stderr "Couldn't open database '$name': $err"
	exit 1
    }

    return $io
}

# Parse -minus args according to $v_opts
proc parse_args {opts v_opts args} {
    upvar $v_opts opt

    # Set default values
    foreach {tag val def} $opts {
	set opt($tag) $def
    }

    while {[string match "-*" [lindex $args]]} {
	set opcode [string range [lindex $args 0] 1 end]

	set skip 0
	foreach {tag val def} $opts {
	    if {$opcode != $tag} continue
	    
	    if {$val == 0} {
		set opt($tag) 1
		set skip 1
		break
	    } else {
		if {[llength $args] < 2} {
		    puts stderr "Option '$opcode' requires an argument"
		    exit 1
		}
		set opt($tag) [lindex $args 1]
		set skip 2
		break
	    }
	}

	if {$skip == 0} {
	    puts stderr "Unknown option '$opcode'"
	    usage 1
	}

	set args [lrange $args $skip end]
    }

    return $args
}

#-----------------------------------------------------------------------------
# Startup code
source $env(STADTABL)/shlib.conf
load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}
load_package tk_utils
tk_utils_init
load_package gap5

eval main $argv

