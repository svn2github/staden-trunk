#!/bin/sh
#\
exec stash "$0" ${@+"$@"} || exit 1

#
# This converts SVEC tags into CVEC tags.
# The purpose for this is that the prefinishing code stops trying to extend
# a contig when it hits the cloning vector, but attempts to hunt
# around for other templates. This also prevents the first pass from
# finding a 'problem' (non existant) which isnt solved, but picking
# perhaps another experiment anyway in order to double-strand.
#
load_package gap

proc open_database {argv} {
    foreach {dbname dbvers} [split [lindex $argv 0] .] {}
    return [open_db -name $dbname -version $dbvers -access rw]
}

set io [open_database $argv]
set db [io_read_database $io]
set na [keylget db Nannotations]
for {set anno 1} {$anno <= $na} {incr anno} {
    set a [io_read_annotation $io $anno]
    if {[keylget a type] == "SVEC"} {
	puts "Updating tag $anno"
        keylset a type CVEC
	io_write_annotation $io $anno $a
    }
}
    
close_db -io $io
exit
