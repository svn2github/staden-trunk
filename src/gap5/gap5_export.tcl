#!/bin/sh

# A command line tool to run the export sequences code of Gap5.
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
load $env(STADLIB)/${lib_prefix}tgap${lib_suffix} g5


#-----------------------------------------------------------------------------
# Argument parsing
proc usage {e} {
    puts ""
    puts {Usage:  gap5_export [options] DBNAME.VERS}
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

set opts {format 1 sam
          h 0 0
          help 0 0
          contigs 1 "*"
          fixmates 0 0
          out 1 "out.*"}
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

if {[catch {set io [g5::open_database -name $db -access ro]} err]} {
    puts stderr "Couldn't open database '$db': $err"
    exit 1
}

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

exit
