#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc ReadPairDialog { io f} {
    global gap5_defs

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Find read pairs"

    ###########################################################################
    #input 
    lorf_in $f.infile [keylget gap5_defs READPAIR.INFILE] \
	"" -bd 2 -relief groove

    # Comparison mode
    radiolist $f.mode \
	-title "Comparison mode" \
	-bd 2 -relief groove -orient horizontal \
	-default [keylget gap5_defs READPAIR.MODE] \
	-buttons { \
	    {{End vs End   }}
	    {{End vs All   }}
	    {{All vs All   }}
	}

    xentry $f.end_size \
	-label "Size of contig 'ends'" \
	-default [keylget gap5_defs READPAIR.END_SIZE]

    scalebox $f.min_mq \
	-title "Minimum mapping quality" \
	-from 0 -to 100 \
	-default [keylget gap5_defs READPAIR.MIN_MAP_QUAL] \
	-width 5 \
	-type CheckInt \
	-orient horiz

    scalebox $f.min_freq \
	-title "Minimum Spanning Frequency" \
	-from 0 -to 100 \
	-default [keylget gap5_defs READPAIR.MIN_FREQ] \
	-width 5 \
	-type CheckInt \
	-orient horiz

    # Produce a listbox of library names
    label $f.spacer -text ""
    labelframe $f.libs -text ""
    label $f.libs.label -text "By default templates in all libraries are used to spot read pairs. To restrict to specific libraries, select them from the list below." -wrap 400 -justify left
    tablelist $f.libs.tl \
	-height 5 \
	-columns {15 Name 10 "Pair count" 10 "Insert size"} \
	-selectmode extended \
	-exportselection 0 \
	-stretch 0 \
	-yscrollcommand [list $f.libs.ys set]
    scrollbar $f.libs.ys -command "$f.libs.tl yview" -orient vertical
    pack $f.libs.label -side top -fill both
    pack $f.libs.tl -side left -expand 1 -fill both
    pack $f.libs.ys -side right -fill both
    set db [$io get_database]
    set nl [$db get_num_libraries]
    for {set i 0} {$i < $nl} {incr i} {
	set rec [$db get_library_rec $i]
	set lib [$io get_library $rec]
	$lib update_stats

	set count [$lib get_count]
	set size  [$lib get_insert_size]
	for {set max 0; set k 0; set j 0} {$j < 3} {incr j} {
	    if {$max < [lindex $count $j]} {
		set max [lindex $count $j]
		set k $j
	    }
	}
	set count [lindex $count $k]
	set size  [lindex $size $k]

	$f.libs.tl insert end [list [$lib get_name] $count $size]

	$lib delete
	upvar \#0 $f.libs.tl l_rec
	set l_rec($i) $rec
    }

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "ReadPairs_OK_Pressed $io $f $f.infile $f.mode $f.end_size $f.min_mq $f.min_freq $f.libs.tl"\
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap5 {Read Pairs}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################
    #final packing

    pack $f.infile $f.mode $f.end_size $f.min_mq $f.min_freq $f.spacer \
	$f.libs $f.ok_cancel -fill x

}

proc ReadPairs_OK_Pressed {io f infile mode end_size min_mq min_freq lib_w} {
    global gap5_defs
    upvar \#0 $lib_w l_rec

    if {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $infile]
    }

    set end_size [$end_size get]
    set min_mq [scalebox_get $min_mq]
    set min_freq [scalebox_get $min_freq]
    set mode [lindex {end_end end_all all_all} [expr {[radiolist_get $mode]-1}]]
    set libs {}
    foreach idx [$lib_w curselection] {
	lappend libs $l_rec($idx)
    }

    set db [$io get_database]
    if {[llength $libs] == [$db get_num_libraries]} {
	puts "All libs"
	set libs {}
    }

    destroy $f

    ContigComparator $io

    SetBusy
    find_read_pairs \
	-io           $io \
	-contigs      $list \
	-mode         $mode \
	-end_size     $end_size \
	-min_map_qual $min_mq \
	-min_freq     $min_freq \
	-libraries    $libs
    ClearBusy
}
