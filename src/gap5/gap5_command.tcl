#!/bin/sh

# A command line tool to run the export consensus sequences code of Gap5.
#
# See the usage function below or run the program with no output for
# usage information.

#\
exec tclsh $0 ${@+"$@"}

namespace eval subcmd {}

#-----------------------------------------------------------------------------
# "consensus" saving

namespace eval subcmd::consensus {

proc main {args} {
    set fmt {format 1 fastq
	h 0 0
	help 0 0
	contigs 1 "*"
	strip_pads 0 0
	out 1 "cons.*"}

    set args [eval parse_args [list $fmt] opt $args]

    if {[llength $args] != 1} {
	usage 1
    }

    if {$opt(h) || $opt(help)} {
	usage 0
    }

    set io [open_db [lindex $args 0]]

    if {$opt(contigs) == "*"} {
	set opt(contigs) [CreateAllContigList $io]
    }

    set opt(out) [regsub {\*} $opt(out) $opt(format)]
    set format [lsearch {- fastq fasta exp} $opt(format)]

    if {[catch {get_consensus \
		    -io $io \
		    -contigs $opt(contigs) \
		    -format $format \
		    -strip_pads $opt(strip_pads) \
		    -outfile $opt(out)} err]} {
	puts stderr "Failed in get_consensus call: $err"
	exit 1
    }
    
    #$io close
    
    return 0
}

proc usage {e} {
    puts ""
    puts {Usage:  gap5_command consensus [options] DBNAME.VERS}
    puts ""
    puts {Where [options] are any of the following:}
    puts "    -h                   Shows this help."
    puts "    -help                Shows this help."
    puts "    -format 'fmt'        Controls the output format. 'fmt' should be one"
    puts "                           of fasta or fastq (default)."
    puts "    -contigs 'list'      Output only specific contigs. 'list' is a space"
    puts "                           separated list of contig names."
    puts "    -out 'filename'      Where to write the ouput. Defaults to out.'fmt'."
    puts "    -strip_pads          Removes padding characters."
    puts ""
    exit $e
}

}; #namespace eval

#-----------------------------------------------------------------------------
# "export" sequences
namespace eval subcmd::export {

proc main {args} {
    set fmt {format 1 sam
	h 0 0
	help 0 0
	contigs 1 "*"
	fixmates 0 0
	out 1 "-"}

    set args [eval parse_args [list $fmt] opt $args]

    if {[llength $args] != 1} {
	usage 1
    }

    if {$opt(h) || $opt(help)} {
	usage 0
    }

    set io [open_db [lindex $args 0]]

    if {$opt(contigs) == "*"} {
	set opt(contigs) [CreateAllContigList $io]
    }

    set opt(out) [regsub {\*} $opt(out) $opt(format)]

    if {[catch {export_contigs \
		    -io $io \
		    -contigs $opt(contigs) \
		    -format $opt(format) \
		    -fixmates $opt(fixmates) \
		    -outfile $opt(out)} err]} {
	puts stderr "Failed in export_contigs call: $err"
	#$io close
	exit 1
    }

    #$io close

    return 0
}

proc usage {e} {
    puts ""
    puts {Usage:  gap5_command export [options] DBNAME.VERS}
    puts ""
    puts {Where [options] are any of the following:}
    puts "    -h                   Shows this help."
    puts "    -help                Shows this help."
    puts "    -format 'fmt'        Controls the output format. 'fmt' should be one"
    puts "                           of sam, ace, baf, fasta or fastq."
    puts "    -contigs 'list'      Output only specific contigs. 'list' is a space"
    puts "                           separated list of contig names."
    puts "    -out 'filename'      Where to write the ouput. Defaults to out.'fmt'."
    puts ""
    exit $e
}

}; #namespace eval

#-----------------------------------------------------------------------------
# "check" database
namespace eval subcmd::check {

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

}; #namespace eval

#-----------------------------------------------------------------------------
# Utility functions

# Error reporting, override the tk route
catch {rename tk_messageBox {}}
proc tk_messageBox {args} {
    foreach {a b} $args {
	set opt($a) $b
    }

    if {[info exists opt(-icon)] && $opt(-icon) == "error"} {
	global errorCode errorInfo
	puts stderr "ERROR: $opt(-message)"

	puts "\nError code: $errorCode"
	puts "\nError message:\n$errorInfo"
    } else {
	puts $opt(-message)
    }
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
# Startup code. Main()
source $env(STADTABL)/shlib.conf
load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}
load_package tk_utils
tk_utils_init
load_package gap5

if {[llength $argv] < 1} {
    puts stderr {gap5_command [sub-command] [arguments]... }
    exit 1
}

set cmd [lindex $argv 0]
if {[info command subcmd::${cmd}::main] == ""} {
    foreach _ [namespace children ::subcmd] {
	lappend cmds [lindex [regexp -inline {::subcmd::(.*)} $_] 1]
    }
    set l [join $cmds ", "]
    puts stderr "Unknown sub-command \"$cmd\"- should be one of: $l."
    exit 1
}

exit [eval subcmd::${cmd}::main [lrange $argv 1 end]]
