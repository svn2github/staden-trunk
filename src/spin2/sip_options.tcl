#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

##############################################################################
#change the default protein substitution matrix
proc ChangeScoreMatrix { } {
    global sip_defs PROTEIN DNA

    if {[winfo exists .score_matrix]} {raise .score_matrix; return}

    set f .score_matrix
    xtoplevel $f -resizable 0
    wm title $f "Score matrix"

    getFname $f.file "matrix name" load_optional_default {} "<identity>"

    set name [get_score_matrix -type $PROTEIN]
    entrybox_configure $f.file.entry -default $name -width 25
    $f.file.entry.entry xview [string last / $name]

    #score matrix type
    keylset sm TYPE [keylget sip_defs SIP.SCORE.TYPE]
    set b1 [keylget sm TYPE.BUTTON.1]
    set b2 [keylget sm TYPE.BUTTON.2]

    #radiolist $f.type \
	    #-title [keylget sm TYPE.NAME]\
	    #-orient horizontal \
	    #-default [keylget sm TYPE.VALUE] \
	    #-buttons [format { \
	    #{ %s -command "entrybox_configure %s \
	                   -default \[get_score_matrix -type %s\]" } \
	    #{ %s -command "entrybox_configure %s \
	                  -default \[get_score_matrix -type %s\]" } } \
	    #[list $b1] [list $f.file.entry] [list $DNA] \
	    #[list $b2] [list $f.file.entry] [list $PROTEIN] ]

    #radiolist $f.type \
	    -title [keylget sm TYPE.NAME]\
	    -orient horizontal \
	    -default [keylget sm TYPE.VALUE] \
	    -buttons [format { \
	    { %s -command "entrybox_configure %s \
	                   -default \[get_score_matrix -type %s\]" } } \
	    [list $b2] [list $f.file.entry] [list $PROTEIN] ]

    #pack $f.type -fill x
    pack $f.file

    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $f.button -bd 2 -relief groove \
	    -ok_command "set_score_matrix\
	                 -file \[getFname_in_name $f.file\] \
	                 -type $PROTEIN; \
	                 destroy $f"\
	    -cancel_command "destroy $f" \
	    -help_command "show_help spin {SPIN-Changing the score matrix}"
    pack $f.button -fill x
}

proc ConfigMatches { } {

    set f .config_matches
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Configure matches"

    entrybox $f.mm \
	    -title "maximum matches" \
	    -default [get_max_matches] \
	    -width 10 \
	    -type "CheckIntRange 1 [INT_MAX]" 
    
    entrybox $f.dm \
	    -title "default matches" \
	    -default [get_def_matches] \
	    -width 10 \
	    -type "CheckIntRange 1 [INT_MAX]" 

    okcancelhelp $f.button -bd 2 -relief groove \
	    -ok_command "set_max_matches \[entrybox_get $f.mm]; \
	    set_def_matches \[entrybox_get $f.dm]; \
	    destroy $f"\
	    -cancel_command "destroy $f" \
	    -help_command "show_help spin {SPIN-Changing Max Match Number}"

    pack $f.mm $f.dm $f.button -fill x
}

