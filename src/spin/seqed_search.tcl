#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc SeqedSearchDialog { t search_cmd cancel_cmd help_cmd} {
    global nip_defs

    if {[xtoplevel $t -resizable 0] == ""} {
	return -1
    }
    wm title $t "String search"

    #strand
    strand $t.strand

    #search algorithm - either use iub codes or literal search
    keylset ll USE_IUB [keylget nip_defs SEQED.SEARCH.USE_IUB]
    set b1 [keylget ll USE_IUB.BUTTON1.NAME]
    set b2 [keylget ll USE_IUB.BUTTON2.NAME]

    radiolist $t.use_iub \
	-title [keylget ll USE_IUB.NAME]  \
	-default [keylget ll USE_IUB.VALUE]\
	-orient horizontal \
	-buttons [format {{%s -value %s} {%s -value %s}} \
		      [list $b1] [list [keylget ll USE_IUB.BUTTON1.VALUE]] [list $b2] [list [keylget ll USE_IUB.BUTTON2.VALUE]]]

    #percentage match
    keylset ll MATCH [keylget nip_defs SEQED.SEARCH.MATCH]
    entrybox $t.match \
	-title "[keylget ll MATCH.NAME]" \
	-default [keylget ll MATCH.VALUE]\
	-width 5 \
	-type "CheckFloatRange 1.0 100.0"

    #search string
    keylset ll STRING [keylget nip_defs SEQED.SEARCH.STRING]
    entrybox $t.string \
	-title "[keylget ll STRING.NAME]" \
	-default [keylget ll STRING.VALUE]\
	-width 20 \
	-type "CheckIUBCString"

    #########################################################################
    #search cancel help buttons 
    frame $t.bu -bd 2 -relief groove
    button $t.bu.search -text "Search" \
	-command $search_cmd
    button $t.bu.quit -text "Cancel" \
	-command $cancel_cmd
    button $t.bu.help -text "Help" \
	-command "$help_cmd"
    pack $t.bu.search $t.bu.quit $t.bu.help -side left -expand 1
    pack $t.strand $t.use_iub $t.match $t.string $t.bu -fill x

    return 0
}

proc SeqedSearch_string {t } {

    return [entrybox_get $t.string]
}

proc SeqedSearch_direction {t } {

    if {[radiolist_get $t.direction] == 2} {
	set direct backward
    } else {
	set direct forward
    }

    return $direct
}

proc SeqedSearch_strand {t} {

    if {[radiolist_get $t.strand] == 1} {
	set strand +
    } else {
	set strand -
    }

    return [strand_get $t.strand]
}

proc SeqedSearch_match {t} {

    return [entrybox_get $t.match]
    
}
proc SeqedSearch_use_iub {t } {

    return [radiolist_get $t.use_iub]
}
