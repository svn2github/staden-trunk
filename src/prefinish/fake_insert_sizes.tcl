#!/bin/sh
#\
STADENROOT=/usr/local/badger/gap4_test; export STADENROOT; . $STADENROOT/staden.profile; exec stash "$0" ${@+"$@"} || exit 1

#
# This sets the insert size fields to min 1 max 100000 for all
# templates. Currently the fake database fields for these EST
# means that this info gets missed out of the experiment files
# and hence we end up with 0-0 as the min-max size in the database.
load_package gap

proc open_database {argv} {
    foreach {dbname dbvers} [split [lindex $argv 0] .] {}
    return [open_db -name $dbname -version $dbvers -access rw]
}

set io [open_database $argv]
set db [io_read_database $io]
set nt [keylget db Ntemplates]
for {set t 1} {$t <= $nt} {incr t} {
    set temp [io_read_template $io $t]
    keylset temp insert_length_min 1
    keylset temp insert_length_max 100000
    io_write_template $io $t $temp
}
    
io_flush $io
close_db -io $io
exit
