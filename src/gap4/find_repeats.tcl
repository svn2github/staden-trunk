#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc FindRepeatsDialog { io f} {
    global gap_defs

    #puts "start FINDREPEATSDIALOG"
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Find repeats"

    ###########################################################################
    #input 
    #contig identifier widget
    contig_id $f.id -io $io 

    lorf_in $f.infile [keylget gap_defs FINDREP.INFILE] \
	    "{contig_id_configure $f.id -state disabled} \
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state disabled}\
	    {contig_id_configure $f.id -state normal}" -bd 2 -relief groove

    ###########################################################################

    keylset mr MINREP [keylget gap_defs FINDREP.MINREP]
    entrybox $f.min_rpt \
	    -title "[keylget mr MINREP.NAME] ([keylget mr MINREP.MIN], \
	    [keylget mr MINREP.MAX])"\
	    -width 5 \
	    -default [keylget mr MINREP.VALUE] \
	    -type "CheckIntRange [keylget mr MINREP.MIN] [keylget mr MINREP.MAX]"

    ###########################################################################
    #select task

    keylset st SELTASK [keylget gap_defs FINDREP.SELTASK]
    set b1 [keylget st SELTASK.BUTTON.1]
    set b2 [keylget st SELTASK.BUTTON.2]
    set b3 [keylget st SELTASK.BUTTON.3]
    radiolist $f.sel_task \
	    -title [keylget st SELTASK.NAME] \
	    -bd 2 \
	    -relief groove \
	    -default [keylget st SELTASK.VALUE] \
	    -buttons [format { {%s} {%s} {%s} } \
	    [list $b1] [list $b2] [list $b3] ]

    ###########################################################################
    #select mode
    SetDefaultTags FINDREP.TAGS

    keylset sm SELMODE [keylget gap_defs FINDREP.SELMODE]
    set b1 [keylget sm SELMODE.BUTTON.1]
    set b2 [keylget sm SELMODE.BUTTON.2]
    frame $f.sel_mode -bd 2 -relief groove
    button $f.sel_mode.but \
	    -text "Select tags" \
	    -command "TagDialog FINDREP.TAGS \
			 $f[keylget gap_defs SELECT_TAGS.WIN] {}"

    radiolist $f.sel_mode.rl \
	    -title [keylget sm SELMODE.NAME] \
	    -orient horizontal\
	    -default [keylget sm SELMODE.VALUE]\
	    -buttons [format { \
	    { %s -command { %s configure -state normal; \
	    SetDefaultTags %s }} \
	    {%s -command { %s configure -state disabled}}} \
	    [list $b1] [list $f.sel_mode.but] FINDREP.TAGS \
	    [list $b2] [list $f.sel_mode.but] ]
    pack $f.sel_mode.rl -side left
    pack $f.sel_mode.but -side right

    ###########################################################################
    #add tags to db

    keylset at ADDTAGS [keylget gap_defs FINDREP.ADDTAGS]
    frame $f.add_tags -relief raised -bd 2
    yes_no $f.add_tags.yn \
	    -title [keylget at ADDTAGS.NAME] \
	    -default [keylget at ADDTAGS.VALUE]
    #pack $f.add_tags.yn -side top

    ###########################################################################
    #save tags to file

    frame $f.save_tags -bd 2 -relief groove
    entrybox $f.save_tags.name \
	-title   [keylget gap_defs FINDREP.SAVETAGS.SAVENAME.NAME] \
	-default [keylget gap_defs FINDREP.SAVETAGS.SAVENAME.VALUE] \
	-type CheckOutput

    yes_no $f.save_tags.yn \
	-title   [keylget gap_defs FINDREP.SAVETAGS.SAVEYN.NAME] \
	-default [keylget gap_defs FINDREP.SAVETAGS.SAVEYN.VALUE] \
	-ycommand "entrybox_configure $f.save_tags.name -state normal" \
	-ncommand "entrybox_configure $f.save_tags.name -state disabled" \
	-orient horizontal

    pack $f.save_tags.yn $f.save_tags.name -side top -fill both

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "Find_Rep_OK_Pressed $io $f $f.infile $f.id \
	    		 $f.sel_task $f.sel_mode.rl $f.add_tags.yn \
			 \[entrybox_get $f.min_rpt\] $f.save_tags.yn \
			 $f.save_tags.name" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Repeats}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################
    #final packing
    
    pack $f.infile -side top -fill both
    #pack $f.proc_mode -side top -fill both
    pack $f.id -side top -fill both
    pack $f.min_rpt -side top -fill both
    pack $f.sel_task -side top -fill both
    pack $f.sel_mode -side top -fill both
    pack $f.add_tags -side top -fill both
    pack $f.save_tags -side top -fill both
    pack $f.ok_cancel -side top -fill both
}

proc Find_Rep_OK_Pressed {io f infile id sel_task sel_mode add_tags min_rpt\
	save_tags_yn save_tags_name} {
    global gap_defs

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
	
    set masking [radiolist_get $sel_mode]
    #if masking mode is 1 (mask active tags)
    #if masking mode is 2 (no masking) set masking variable to be 0
    if {$masking == 1} {
	set active_tags [GetDefaultTags FINDREP.TAGS]
    } else {
	set active_tags {}
    }

    if {[yes_no_get $save_tags_yn]} {
	set outfile [entrybox_get $save_tags_name]
    } else {
	set outfile ""
    }


    set sel_task [radiolist_get $sel_task]

    destroy $f

    ContigComparator $io


    # If repeats are found, this also sets the tag_list variable
    SetBusy
    find_repeats \
	-io $io \
	-direction $sel_task \
        -min_match $min_rpt\
	-contigs $list \
	-outfile $outfile \
	-tag_types $active_tags
    ClearBusy
}

