#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc FindOligos { io } {
    global gap_defs

    set f [keylget gap_defs FINDOLIGO.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Sequence search"

    ###########################################################################
    #input 
    #contig identifier widget
    contig_id $f.id -io $io 

    lorf_in $f.infile [keylget gap_defs FINDOLIGO.INFILE] \
	    "{contig_id_configure $f.id -state disabled} \
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state normal}" -bd 2 -relief groove

    ###########################################################################
    #mis-match scale
    keylset mm MAXMIS [keylget gap_defs FINDOLIGO.MAXMIS]
    scalebox $f.mis_match \
	    -title [keylget mm MAXMIS.NAME]\
	    -orient horizontal \
	    -to [keylget mm MAXMIS.MAX]\
	    -from [keylget mm MAXMIS.MIN] \
	    -resolution [keylget mm MAXMIS.RES] \
	    -default [keylget mm MAXMIS.VALUE] \
	    -width 5 \
	    -type CheckFloat

    ###########################################################################
    #select mode
    SetDefaultTags FINDOLIGO.TAGS

    keylset sm SELMODE [keylget gap_defs FINDOLIGO.SELMODE]
    keylset ms MATCHSEQ [keylget gap_defs FINDOLIGO.MATCHSEQ]
    set b1 [keylget sm SELMODE.BUTTON.1]
    set b2 [keylget sm SELMODE.BUTTON.2]
    set b3 [keylget sm SELMODE.BUTTON.3]
    frame $f.sel_mode -bd 2 -relief groove
    frame $f.sel_mode.l
    frame $f.sel_mode.r
    button $f.sel_mode.r.but \
	    -text "Select consensus tags" \
	    -command "TagDialog FINDOLIGO.TAGS \
			$f[keylget gap_defs SELECT_TAGS.WIN] {}"
    xentry $f.sel_mode.r.seq \
	-width 30 \
	-default [keylget ms MATCHSEQ.VALUE]

    xget_fname $f.sel_mode.r.fasta \
	-type load \
	-default [keylget ms MATCHSEQ.FASTA]

    radiolist $f.cons_or_seq \
	-orient horizontal \
	-title "Where to search" \
	-buttons {Consensus Sequences} \
	-default [keylget gap_defs FINDOLIGO.CONSONLY.VALUE]

    radiolist $f.hidden_data \
	-orient horizontal \
	-title "Use hidden data" \
	-buttons {No Yes} \
	-default [keylget gap_defs FINDOLIGO.CUTOFFS.VALUE]

    radiolist $f.sel_mode.l.rl \
	    -title [keylget sm SELMODE.NAME] \
	    -orient vertical\
	    -default [keylget sm SELMODE.VALUE]\
	    -buttons [format { \
	    { %s -command { \
		SetDefaultTags FINDOLIGO.TAGS; \
	        %s configure -state normal; \
		%s configure -state disabled; \
		%s configure -state disabled }\
	    } \
	    { %s -command { \
		%s configure -state disabled;\
		%s configure -state normal;\
		%s configure -state disabled }\
	    } \
	    { %s -command { \
		%s configure -state disabled;\
		%s configure -state disabled;\
		%s configure -state normal }\
	    }} \
	    [list $b1] \
	        [list $f.sel_mode.r.but] \
		[list $f.sel_mode.r.seq] \
		[list $f.sel_mode.r.fasta] \
	    [list $b2] \
	        [list $f.sel_mode.r.but] \
		[list $f.sel_mode.r.seq] \
		[list $f.sel_mode.r.fasta] \
	    [list $b3] \
	        [list $f.sel_mode.r.but] \
		[list $f.sel_mode.r.seq] \
		[list $f.sel_mode.r.fasta]]
    pack $f.sel_mode.l.rl
    pack $f.sel_mode.r.fasta -side bottom
    pack $f.sel_mode.r.seq -side bottom
    pack $f.sel_mode.r.but -side bottom -anchor e
    pack $f.sel_mode.l -side left
    pack $f.sel_mode.r -side right -fill y

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "FindOligo_OK_Pressed $io $f $f.infile $f.id \
	    		 $f.sel_mode.l.rl $f.mis_match $f.sel_mode.r.seq \
			 $f.cons_or_seq $f.hidden_data $f.sel_mode.r.fasta" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Find Oligos}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################
    #final packing
    
    pack $f.infile -side top -fill both
    pack $f.id -side top -fill both
    pack $f.cons_or_seq -side top -fill both
    pack $f.hidden_data -side top -fill both
    pack $f.mis_match -side top -fill both
    pack $f.sel_mode -side top -fill both
    pack $f.ok_cancel -side top -fill both
}

proc FindOligo_OK_Pressed {io f infile id sel_mode mis_match seq cons_or_seq hidden fasta} {
    global CurContig
    global NGRec
    global LREG
    global RREG
    global gap_defs

    set active_tags {}

    if {[lorf_in_get $infile] == 4} {
	set gel_name [contig_id_gel $id]
	set lreg [contig_id_lreg $id]
	set rreg [contig_id_rreg $id]
	
	SetContigGlobals $io $gel_name $lreg $rreg
	set list "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $infile] == 3 } {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $infile]
    }
	
    set tags [radiolist_get $sel_mode]
    set sequence {}
    set fastafn {}
    #if tags mode is 1 (use tags)
    #if tags mode is 2 (use sequence)
    #if tags mode is 3 (use fasta file)
    if {($tags == 1)} {
	set active_tags [GetDefaultTags FINDOLIGO.TAGS]
    } elseif {$tags == 2} {
	set sequence [$seq get]
	set sequence [string map {{*} {}} $sequence]
    } elseif {$tags == 3} {
	set fastafn [$fasta get]
	if {$fastafn == ""} {
	    bell
	    return
	}
    }

    if {$tags == 1 && $active_tags == ""} {
	set re_enter 0
	tk_messageBox -icon error -type ok -title "No active tags" \
		-message "No tags have been selected" \
	        -parent $f
	bell
	#wait forever...
	tkwait variable re_enter

    } elseif {$tags == 2 && $sequence == ""} {
	set re_enter 0
	tk_messageBox -icon error -type ok -title "No search sequence" \
		-message "No sequence has been entered" \
	        -parent $f
	bell
	#wait forever...
	tkwait variable re_enter
    }


    set mis_match [scalebox_get $mis_match]
    set use_hidden [expr {[radiolist_get $hidden]-1}]
    set cons_only [expr {[radiolist_get $cons_or_seq]==1?1:0}]

    destroy $f

    ContigComparator $io

    # If repeats are found, this also sets the tag_list variable
    SetBusy
    find_oligo -io $io \
            -min_pmatch $mis_match \
	    -contigs $list \
	    -seq $sequence \
	    -tag_types $active_tags \
	    -file $fastafn \
    	    -consensus_only $cons_only \
    	    -cutoffs $use_hidden
    ClearBusy
}

