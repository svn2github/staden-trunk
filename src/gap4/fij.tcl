#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc FIJDialog { f io } {
    global gap_defs
    global LREG

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Find internal joins"

    ###########################################################################
    #input 

    lorf_in $f.infile [keylget gap_defs FIJ.INFILE] "" -bd 2 -relief groove

    ###########################################################################
    #match scales

    frame $f.match -relief groove -bd 2
    frame $f.match.s -relief groove -bd 2
    frame $f.match.f -relief groove -bd 2
    label $f.match.s.label -text "Sensitive algorithm:"
    label $f.match.f.label -text "Quick algorithm:"


    ###########################################################################
    set mm [keylget gap_defs FIJ.MINOVERLAP]
    scalebox $f.match.min_overlap \
	    -title [keylget mm NAME]\
	    -orient horizontal \
	    -to [keylget mm MAX] \
	    -from [keylget mm MIN]\
	    -default [keylget mm VALUE] \
	    -width 5 \
	    -type CheckInt
    
    set mm [keylget gap_defs FIJ.MINMATCH]
    scalebox $f.match.f.min_match \
	    -title [keylget mm NAME]\
	    -orient horizontal \
	    -to [keylget mm MAX]\
	    -from [keylget mm MIN] \
	    -default [keylget mm VALUE] \
	    -width 5 \
	    -type CheckInt

    set mm [keylget gap_defs FIJ.MAXDIAG]
    entrybox $f.match.s.max_diag \
	-title "[keylget mm NAME] ([keylget mm MIN] to [keylget mm MAX])"\
	-default [keylget mm VALUE]\
	-type "CheckFloatRange [keylget mm MIN] [keylget mm MAX]"

    set mm [keylget gap_defs FIJ.BANDSIZE]
    scalebox $f.match.s.band_size \
	    -title [keylget mm NAME]\
	    -orient horizontal \
	    -to [keylget mm MAX] \
	    -from [keylget mm MIN]\
	    -default [keylget mm VALUE] \
	    -width 5 \
	    -type CheckInt

    set mm [keylget gap_defs FIJ.USEBAND] 
    yes_no $f.match.f.use_band \
    	    -title [keylget mm NAME] \
	    -orient horizontal \
	    -bd 0 \
	    -default [keylget mm VALUE]

    set mm [keylget gap_defs FIJ.MAXMIS]
    scalebox $f.match.max_mis \
	    -title [keylget mm NAME]\
	    -orient horizontal \
	    -to [keylget mm MAX]\
	    -from [keylget mm MIN] \
	    -resolution [keylget mm RES] \
	    -default [keylget mm VALUE] \
	    -width 5 \
	    -type CheckFloat


    ###########################################################################
    #select word length

    set st [keylget gap_defs FIJ.WORDLENGTH]
    set b1 [keylget st BUTTON.1]
    set b2 [keylget st BUTTON.2]

    radiolist $f.match.word_length \
	    -title [keylget st NAME]\
	    -bd 0 \
	    -relief groove \
	    -orient horizontal \
	    -default [keylget st VALUE] \
	    -buttons [format { {%s} {%s} } \
	    [list $b1] [list $b2] ]

    frame $f.match.padding -relief groove -bd 2 -height 2


    set ff $f.match
    radiolist $f.match.blocks \
	-title "Alignment algorithm" \
	-bd 0 \
	-relief groove \
	-orient horizontal \
	-default [keylget gap_defs FIJ.ALIGN_MODE]\
	-buttons [format { \
            {quick     -command {scalebox_configure %s -state disabled; \
				 entrybox_configure %s -state disabled;
				 scalebox_configure %s -state normal;
				 yes_no_configure %s -state normal}} \
            {sensitive -command {scalebox_configure %s -state normal; \
				 entrybox_configure %s -state normal;
				 scalebox_configure %s -state disabled;
				 yes_no_configure %s -state disabled}}} \
            $ff.s.band_size $ff.s.max_diag $ff.f.min_match $ff.f.use_band \
            $ff.s.band_size $ff.s.max_diag $ff.f.min_match $ff.f.use_band]

    pack $f.match.word_length -fill x
    pack $f.match.min_overlap -fill x
    pack $f.match.max_mis -fill x
    pack $f.match.padding -fill x -pady 5
    pack $f.match.blocks -fill x

    pack $f.match.s -fill both -expand 1
    pack $f.match.s.label -anchor w -padx 50
    pack $f.match.s.max_diag -fill x
    pack $f.match.s.band_size -fill x

    pack $f.match.f -fill both -expand 1
    pack $f.match.f.label -anchor w -padx 50
    pack $f.match.f.min_match -fill x
    pack $f.match.f.use_band -fill x

    ###########################################################################
    #contig identifier widget
    set start $LREG
    contig_id $f.sc \
	    -io $io \
	    -range 0

    ###########################################################################
    #select task

    set st [keylget gap_defs FIJ.SELTASK]
    set b1 [keylget st BUTTON.1]
    set b2 [keylget st BUTTON.2]

    radiolist $f.sel_task \
	    -title [keylget st NAME]\
	    -bd 2 \
	    -relief groove \
	    -orient horizontal \
	    -default [keylget st VALUE] \
	    -buttons [format { \
	    {%s -command {contig_id_configure %s -state disabled}} \
	    {%s -command {contig_id_configure %s -state normal}}} \
	    [list $b1] [list $f.sc] \
	    [list $b2] [list $f.sc] ]
    
    ###########################################################################
    #hidden data

    frame $f.hidden -bd 2 -relief groove

    HiddenParametersDialogInit $f.ops
    button $f.hidden.options \
	-text "Edit hidden data paramaters" \
	-command "HiddenParametersDialog $f $f.ops"

    set hd [keylget gap_defs FIJ.HIDDEN]
    yes_no $f.hidden.yn \
	    -title [keylget hd NAME] \
	    -bd 0 \
	    -orient horizontal \
	    -ycommand "$f.hidden.options configure -state normal" \
	    -ncommand "$f.hidden.options configure -state disabled" \
	    -default [keylget hd VALUE]

    pack $f.hidden.yn -side top -fill x
    pack $f.hidden.options -side top -anchor w
     
    ###########################################################################
    # alignment displays
    entrybox $f.align_length \
	-title "Maximum alignment length to list (bp)" \
	-default [keylget gap_defs FIJ.MAX_ALIGNMENT] \
	-type "CheckIntRange 0 10000000"

    ###########################################################################
    #select mode
    SetDefaultTags FIJ.TAGS

    set sm [keylget gap_defs FIJ.SELMODE]
    set b1 [keylget sm BUTTON.1]
    set b2 [keylget sm BUTTON.2]
    set b3 [keylget sm BUTTON.3]
    frame $f.sel_mode -bd 2 -relief groove
    button $f.sel_mode.but \
	    -text "Select tags" \
	    -command "TagDialog FIJ.TAGS $f[keylget gap_defs SELECT_TAGS.WIN] \
			{}"

    radiolist $f.sel_mode.rl \
	    -title [keylget sm NAME] \
	    -default [keylget sm VALUE]\
	    -buttons [format { \
	    {%s -command { %s configure -state disabled}} \
	    { %s -command { %s configure -state normal; \
	    SetDefaultTags %s }} \
	    { %s -command { %s configure -state normal; \
	    SetDefaultTags %s } } } \
	    [list $b1] [list $f.sel_mode.but] \
	    [list $b2] [list $f.sel_mode.but] FIJ.TAGS \
	    [list $b3] [list $f.sel_mode.but] FIJ.TAGS ]
    pack $f.sel_mode.rl -side left
    pack $f.sel_mode.but -side right

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "FIJ_OK_Pressed $f $io $f.infile $f.sel_task \
	    $f.match.blocks $f.match.min_overlap $f.match.word_length\
	    $f.match.f.min_match $f.match.f.use_band $f.match.s.max_diag\
	    $f.match.s.band_size $f.match.max_mis \
	    $f.sc $f.hidden.yn $f.sel_mode.rl $f.ops $f.align_length" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {FIJ-Dialogue}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################

    pack $f.infile -fill x
    pack $f.sel_task -fill x
    pack $f.sc -fill x
    pack $f.sel_mode -fill x
    pack $f.align_length -fill x
    pack $f.hidden -fill x
    pack $f.match -fill x
    pack $f.ok_cancel -fill x

}

###########################################################################
proc FIJ_OK_Pressed { f io infile sel_task blocks min_overlap word_length \
		    min_match use_band max_diag band_size max_mis sc \
		    yn sel_mode hidden_ops align_length} {
    
    global CurContig
    global LREG
    global RREG
    global NGRec
    global gap_defs

    set segment ""
    set active_tags {}

    if {[lorf_in_get $infile] == 3} {
	set list [CreateAllContigList $io]
    } else {
	set list [lorf_get_list $infile]
    }

    if {[radiolist_get $sel_task] == 1} {
	# All against all
	set mode "all_all"
    } else {
	# Single against all
	set mode "segment"
	set list [linsert $list 0 [contig_id_gel $sc]]
	SetContigGlobals $io [contig_id_gel $sc]
    }

    if {[yes_no_get $yn]} {
        upvar #0 $hidden_ops data
	set win_size $data($data(which_mode)_win_size)
	set max_dash $data(base_max_dash)
	set min_conf $data(conf_min_conf)
	set use_conf [lsearch {base conf} $data(which_mode)]
	set use_hidden 1
    } else {
	set win_size 0
	set max_dash 0
	set min_conf 0
	set use_conf 0
	set use_hidden 0
    }
    set masking [radiolist_get $sel_mode]

    #if masking mode is 2 or 3 (mark or mask active tags)
    if {($masking == 2) || ($masking == 3)} {
        set active_tags [GetDefaultTags FIJ.TAGS]
    }


    set word_length [lindex {? 8 4} [radiolist_get $word_length]]
    
    set max_prob [entrybox_get $max_diag]
    set max_align_length [entrybox_get $align_length]

    if {[radiolist_get $blocks] == 1} {
        # Quick method
	set min_match [scalebox_get $min_match]
	set band_size [yes_no_get $use_band]
    } else {
	# Sensitive method
	set min_match 0
	set band_size [scalebox_get $band_size]
    }


    set min_overlap [scalebox_get $min_overlap]
    set max_mis [scalebox_get $max_mis]


    # Destroy dialog before showing plot to avoid window activation problems
    destroy $f


    #draw the contig comparator and dot plot
    ContigComparator $io 


    SetBusy
    find_internal_joins -io $io \
	    -mask [lindex {"" none mark mask} $masking] \
	    -mode $mode \
	    -min_overlap $min_overlap \
	    -max_pmismatch $max_mis \
	    -word_length $word_length \
	    -max_prob $max_prob\
	    -min_match $min_match \
	    -band $band_size \
	    -win_size $win_size \
	    -max_dashes $max_dash \
	    -min_conf $min_conf \
	    -use_conf $use_conf \
	    -use_hidden $use_hidden \
	    -tag_types $active_tags \
	    -max_display $max_align_length \
	    -contigs $list
    ClearBusy
}

# Hidde data parameters box, used by Find Internal Joins and Calculate
# Consensus
# $w is the parent window
# $aname is a global array name in which to fill out the results
proc HiddenParametersDialog {w aname} {
    global gap_defs
    upvar #0 $aname data

    set f $w.hiddenparam
    if {![winfo exists $f]} {
	xtoplevel $f -resizable 0
	wm title $f "Hidden data parameters"
    } else {
	raise $f
	wm deiconify $f
    }

    ###################################################
    # extend by confidence values

    set w [frame $f.conf -bd 2 -relief groove]

    label $w.label -text [keylget gap_defs HIDDEN.CONF.LABEL]

    set ws [keylget gap_defs HIDDEN.CONF.WINSIZE]
    scalebox $w.win_size \
	    -title [keylget ws NAME] \
	    -orient horizontal \
	    -to [keylget ws MAX] \
	    -from [keylget ws MIN] \
	    -default $data(conf_win_size) \
	    -width 5 \
	    -type CheckInt

    set md [keylget gap_defs HIDDEN.CONF.MINCONF]
    scalebox $w.min_conf \
	    -title [keylget md NAME] \
	    -orient horizontal \
	    -to [keylget md MAX] \
	    -from [keylget md MIN]\
	    -default $data(conf_min_conf) \
	    -width 5 \
	    -type CheckInt

    pack $w.label -side top -anchor c
    pack $w.win_size -side top -fill x
    pack $w.min_conf -side top -fill x

    ###################################################
    # extend by base calls

    set w [frame $f.unc -bd 2 -relief groove]

    label $w.label -text [keylget gap_defs HIDDEN.UNC.LABEL]

    set ws [keylget gap_defs HIDDEN.UNC.WINSIZE]
    scalebox $w.win_size \
	    -title [keylget ws NAME] \
	    -orient horizontal \
	    -to [keylget ws MAX] \
	    -from [keylget ws MIN] \
	    -default $data(base_win_size) \
	    -width 5 \
	    -type CheckInt

    set md [keylget gap_defs HIDDEN.UNC.MAXDASH]
    scalebox $w.max_dash \
	    -title [keylget md NAME] \
	    -orient horizontal \
	    -to [keylget md MAX] \
	    -from [keylget md MIN]\
	    -default $data(base_max_dash) \
	    -width 5 \
	    -type CheckInt

    pack $w.label -side top -anchor c
    pack $w.win_size -side top -fill x
    pack $w.max_dash -side top -fill x

    ###################################################
    # Main frame and question
    set yn [keylget gap_defs HIDDEN.MODE]
    yes_no $f.which \
	-title [keylget yn NAME] \
	-bd 2 \
	-relief groove \
	-orient horizontal \
	-ycommand "scalebox_configure $f.conf.win_size -state normal;\
                   scalebox_configure $f.conf.min_conf -state normal;\
                   scalebox_configure $f.unc.win_size  -state disabled;\
                   scalebox_configure $f.unc.max_dash  -state disabled" \
	-ncommand "scalebox_configure $f.conf.win_size -state disabled;\
                   scalebox_configure $f.conf.min_conf -state disabled;\
                   scalebox_configure $f.unc.win_size  -state normal;\
                   scalebox_configure $f.unc.max_dash  -state normal" \
	-default [lsearch {base conf} $data(which_mode)]

    ###################################################
    # Final bits

    okcancelhelp $f.ok_cancel \
	-ok_command "HiddenParametersDialog_OK $f $aname" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap4 {FIJ-Dialogue}" \

    pack $f.which $f.conf $f.unc $f.ok_cancel \
	-side top -fill both -expand 1
}

proc HiddenParametersDialogInit {aname} {
    global gap_defs
    upvar #0 $aname data

    set data(which_mode)    [lindex {base conf} [keylget gap_defs HIDDEN.MODE.VALUE]]
    set data(base_win_size) [keylget gap_defs HIDDEN.UNC.WINSIZE.VALUE]
    set data(base_max_dash) [keylget gap_defs HIDDEN.UNC.MAXDASH.VALUE]
    set data(conf_win_size) [keylget gap_defs HIDDEN.CONF.WINSIZE.VALUE]
    set data(conf_min_conf) [keylget gap_defs HIDDEN.CONF.MINCONF.VALUE]
}

proc HiddenParametersDialog_OK {w aname} {
    upvar #0 $aname data

    set data(which_mode)    [lindex {base conf} [yes_no_get   $w.which]]
    set data(base_win_size) [scalebox_get $w.unc.win_size]
    set data(base_max_dash) [scalebox_get $w.unc.max_dash]
    set data(conf_win_size) [scalebox_get $w.conf.win_size]
    set data(conf_min_conf) [scalebox_get $w.conf.min_conf]

    destroy $w
}
