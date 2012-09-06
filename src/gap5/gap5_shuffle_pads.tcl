#!/bin/sh

# A command line tool to run the shuffle pads command of Gap5.
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


#-----------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------
# Argument parsing
proc usage {e} {
    puts ""
    puts {Usage:  gap5_shuffle_pads [options] DBNAME.VERS}
    puts ""
    puts {Where [options] are any of the following:}
    puts "    -h                   Shows this help."
    puts "    -help                Shows this help."
    puts "    -contigs 'list'      Output only specific contigs. 'list' is a space"
    puts "                           separated list of contig names."
    puts "    -band_size num       Set the alignment band size."
    puts ""
    exit $e
}

set opts [list format 1 fastq \
          h 0 0 \
          help 0 0 \
          contigs 1 "*" \
          band_size 0 [keylget gap5_defs SHUFFLE_PADS.BAND_SIZE]]
array set arg {}

foreach {tag val def} $opts {
    set opt($tag) $def
}

while {[string match "-*" [lindex $argv]]} {
    set opcode [string range [lindex $argv 0] 1 end]

    set skip 0
    foreach {tag val def} $opts {
	if {$opcode != $tag} continue

	if {$val == 0} {
	    set opt($tag) 1
	    set skip 1
	    break
	} else {
	    if {[llength $argv] < 2} {
		puts stderr "Option '$opcode' requires an argument"
		exit 1
	    }
	    set opt($tag) [lindex $argv 1]
	    set skip 2
	    break
	}
    }

    if {$skip == 0} {
	puts stderr "Unknown option '$opcode'"
	usage 1
    }

    set argv [lrange $argv $skip end]
}

#-----------------------------------------------------------------------------
# Main entry point
if {$opt(h) || $opt(help)} {
    usage 0
}

if {[llength $argv] != 1} {
    usage 1
}
set db [lindex $argv 0]

if {[catch {set io [g5::open_database -name $db -access rw]} err]} {
    puts stderr "Couldn't open database '$db': $err"
    exit 1
}

if {$opt(contigs) == "*"} {
    set opt(contigs) [CreateAllContigList $io]
}

if {[catch {shuffle_pads \
		-io $io \
		-contigs $opt(contigs) \
		-band $opt(band_size)} err]} {
    puts stderr "Failed in shuffle_pads call: $err"
    #$io close
    exit 1
}

$io flush
$io close

exit
