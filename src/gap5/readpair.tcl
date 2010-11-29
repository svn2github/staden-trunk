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

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "ReadPairs_OK_Pressed $io $f $f.infile $f.mode $f.end_size $f.min_mq"\
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Read Pairs}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################
    #final packing

    pack $f.infile $f.mode $f.end_size $f.min_mq $f.ok_cancel -fill x

}

proc ReadPairs_OK_Pressed {io f infile mode end_size min_mq} {
    global gap5_defs

    if {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $infile]
    }

    set end_size [$end_size get]
    set min_mq [scalebox_get $min_mq]
    set mode [lindex {end_end end_all all_all} [expr {[radiolist_get $mode]-1}]]

    destroy $f

    ContigComparator $io

    SetBusy
    find_read_pairs \
	-io           $io \
	-contigs      $list \
	-mode         $mode \
	-end_size     $end_size \
	-min_map_qual $min_mq
    ClearBusy
}
