#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc FIJDialog { f io } {
    global gap5_defs
    global LREG

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Find internal joins"

    ###########################################################################
    # contig_id boxes for selecting single contigs
    # need both first as FIJ_config_contig_ids gets called by lorf_in

    contig_id $f.id1 -io $io -range 0 -trace 2
    contig_id $f.id2 -io $io -range 0 -trace 0

    ###########################################################################
    #input set 1

    lorf_in $f.infile1 [keylget gap5_defs FIJ.INFILE1] \
	"{FIJ_config_contig_ids $f id1 disabled} \
	 {FIJ_config_contig_ids $f id1 disabled} \
	 {FIJ_config_contig_ids $f id1 disabled} \
	 {FIJ_config_contig_ids $f id1 normal}" \
	-bd 2 -relief groove

    ###########################################################################
    #input set 2

    lorf_in $f.infile2 [keylget gap5_defs FIJ.INFILE2] \
	"{FIJ_config_contig_ids $f id2 disabled} \
	 {FIJ_config_contig_ids $f id2 disabled} \
	 {FIJ_config_contig_ids $f id2 disabled} \
	 {FIJ_config_contig_ids $f id2 normal}" \
	-bd 2 -relief groove

    # selecting a contig_id box makes it the next to be updated
    bind [entrybox_path $f.id1.ent] <<select>> "FIJ_config_contig_ids $f id1"
    bind [entrybox_path $f.id2.ent] <<select>> "FIJ_config_contig_ids $f id2"
    
    FIJ_config_contig_ids $f id1

    ###########################################################################
    #match scales

    frame $f.match -relief groove -bd 2
    frame $f.match.s -relief groove -bd 2
    frame $f.match.f -relief groove -bd 2
    label $f.match.s.label -text "Sensitive algorithm:"
    label $f.match.f.label -text "Quick algorithm:"


    ###########################################################################
    set mm [keylget gap5_defs FIJ.MINOVERLAP]
    scalebox $f.match.min_overlap \
	    -title [keylget mm NAME]\
	    -orient horizontal \
	    -to [keylget mm MAX] \
	    -from [keylget mm MIN]\
	    -default [keylget mm VALUE] \
	    -width 5 \
	    -type CheckInt
    
    set mm [keylget gap5_defs FIJ.MINMATCH]
    scalebox $f.match.f.min_match \
	    -title [keylget mm NAME]\
	    -orient horizontal \
	    -to [keylget mm MAX]\
	    -from [keylget mm MIN] \
	    -default [keylget mm VALUE] \
	    -width 5 \
	    -type CheckInt \
	    -variable $f.MinMatch

    set mm [keylget gap5_defs FIJ.MAXDIAG]
    entrybox $f.match.s.max_diag \
	-title "[keylget mm NAME] ([keylget mm MIN] to [keylget mm MAX])"\
	-default [keylget mm VALUE]\
	-type "CheckFloatRange [keylget mm MIN] [keylget mm MAX]"

    set mm [keylget gap5_defs FIJ.BANDSIZE]
    scalebox $f.match.s.band_size \
	    -title [keylget mm NAME]\
	    -orient horizontal \
	    -to [keylget mm MAX] \
	    -from [keylget mm MIN]\
	    -default [keylget mm VALUE] \
	    -width 5 \
	    -type CheckInt

    set mm [keylget gap5_defs FIJ.USEBAND] 
    yes_no $f.match.f.use_band \
    	    -title [keylget mm NAME] \
	    -orient horizontal \
	    -bd 0 \
	    -default [keylget mm VALUE]

    set mm [keylget gap5_defs FIJ.MAXMIS]
    scalebox $f.match.max_mis \
	    -title [keylget mm NAME]\
	    -orient horizontal \
	    -to [keylget mm MAX]\
	    -from [keylget mm MIN] \
	    -resolution [keylget mm RES] \
	    -default [keylget mm VALUE] \
	    -width 5 \
	    -type CheckFloat

    set mm [keylget gap5_defs FIJ.USEFILTERWORDS]
    yes_no $f.match.f.use_filter \
	-title [keylget mm NAME] \
	-orient horizontal \
	-bd 0 \
	-default [keylget mm VALUE]

    set mm [keylget gap5_defs FIJ.FILTERWORDS]
    xentry $f.match.f.filter_cutoff \
	-label [keylget mm NAME] \
	-default [keylget mm VALUE]


    ###########################################################################
    #select word length

    set st [keylget gap5_defs FIJ.WORDLENGTH]
    set b1 [keylget st BUTTON.1]
    set b2 [keylget st BUTTON.2]
    set b3 [keylget st BUTTON.3]

    radiolist $f.match.word_length \
	    -title [keylget st NAME]\
	    -bd 0 \
	    -relief groove \
	    -orient horizontal \
	    -default [keylget st VALUE] \
	    -buttons [format { {%s} {%s} {%s} } \
			  [list $b1] [list $b2] [list $b3]]

    frame $f.match.padding -relief groove -bd 2 -height 2


    set ff $f.match
    radiolist $f.match.blocks \
	-title "Alignment algorithm" \
	-bd 0 \
	-relief groove \
	-orient horizontal \
	-default [keylget gap5_defs FIJ.ALIGN_MODE]\
	-buttons [format { \
            {fastest   -command {scalebox_configure %s -state disabled; \
				 entrybox_configure %s -state disabled;
				 scalebox_configure %s -state normal;
		                 yes_no_configure %s -state normal;
		                 yes_no_configure %s -state normal;
		                 %s configure -state normal;
		                 set %s.MinMatch 25}} \
            {quick     -command {scalebox_configure %s -state disabled; \
				 entrybox_configure %s -state disabled;
				 scalebox_configure %s -state normal;
				 yes_no_configure %s -state normal;
		                 yes_no_configure %s -state normal;
		                 %s configure -state normal;
	                         set %s.MinMatch 20}} \
            {sensitive -command {scalebox_configure %s -state normal; \
				 entrybox_configure %s -state normal;
				 scalebox_configure %s -state disabled;
		                 yes_no_configure %s -state disabled;
		                 yes_no_configure %s -state disabled;
		                 %s configure -state disabled}}} \
            $ff.s.band_size $ff.s.max_diag $ff.f.min_match $ff.f.use_band \
		      $ff.f.use_filter $ff.f.filter_cutoff $f \
            $ff.s.band_size $ff.s.max_diag $ff.f.min_match $ff.f.use_band \
		      $ff.f.use_filter $ff.f.filter_cutoff $f \
            $ff.s.band_size $ff.s.max_diag $ff.f.min_match $ff.f.use_band \
		      $ff.f.use_filter $ff.f.filter_cutoff]

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
    pack $f.match.f.use_band -fill x
    pack $f.match.f.min_match -fill x
    pack $f.match.f.use_filter -fill x
    pack $f.match.f.filter_cutoff -fill x
    
    ###########################################################################
    #hidden data

    frame $f.hidden -bd 2 -relief groove

    HiddenParametersDialogInit $f.ops
    button $f.hidden.options \
	-text "Edit hidden data paramaters" \
	-command "HiddenParametersDialog $f $f.ops"

    set hd [keylget gap5_defs FIJ.HIDDEN]
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
	-default [keylget gap5_defs FIJ.MAX_ALIGNMENT] \
	-type "CheckIntRange 0 10000000"

#    ###########################################################################
#    #select mode
#    SetDefaultTags FIJ.TAGS
#
#    set sm [keylget gap5_defs FIJ.SELMODE]
#    set b1 [keylget sm BUTTON.1]
#    set b2 [keylget sm BUTTON.2]
#    set b3 [keylget sm BUTTON.3]
#    frame $f.sel_mode -bd 2 -relief groove
#    button $f.sel_mode.but \
#	    -text "Select tags" \
#	    -command "TagDialog FIJ.TAGS $f[keylget gap5_defs SELECT_TAGS.WIN] \
#			{}"
#
#    radiolist $f.sel_mode.rl \
#	    -title [keylget sm NAME] \
#	    -default [keylget sm VALUE]\
#	    -buttons [format { \
#	    {%s -command { %s configure -state disabled}} \
#	    { %s -command { %s configure -state normal; \
#	    SetDefaultTags %s }} \
#	    { %s -command { %s configure -state normal; \
#	    SetDefaultTags %s } } } \
#	    [list $b1] [list $f.sel_mode.but] \
#	    [list $b2] [list $f.sel_mode.but] FIJ.TAGS \
#	    [list $b3] [list $f.sel_mode.but] FIJ.TAGS ]
#    pack $f.sel_mode.rl -side left
#    pack $f.sel_mode.but -side right

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "FIJ_OK_Pressed $f $io $f.infile1 $f.id1 \
            $f.infile2 $f.id2 \
	    $f.match.blocks $f.match.min_overlap $f.match.word_length\
	    $f.match.f.min_match $f.match.f.use_band $f.match.s.max_diag\
	    $f.match.s.band_size $f.match.max_mis \
	    $f.hidden.yn $f.sel_mode.rl $f.ops $f.align_length \
            $f.match.f.use_filter $f.match.f.filter_cutoff" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap5 {FIJ-Dialogue}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################

    pack $f.infile1 -fill x
    pack $f.id1 -fill x
    pack $f.infile2 -fill x
    pack $f.id2 -fill x
#    pack $f.sel_task -fill x
#    pack $f.sc -fill x
    pack $f.align_length -fill x
    pack $f.hidden -fill x
    pack $f.match -fill x
    pack $f.ok_cancel -fill x

}

proc FIJ_config_contig_ids { f id { state unchanged } } {
    # Set which of the contig_id boxes get updated when clicking in the
    # contig list or contig selector, and in which order.
    if { $state != "unchanged" } {
	contig_id_configure "$f.$id" -state $state
    }
    set id1_state [eval [entrybox_path $f.id1.ent] cget -state]
    set id2_state [eval [entrybox_path $f.id2.ent] cget -state]

    set boxen []
    if { "$id" == "id1" } {
	if { $id1_state == "normal" } { lappend boxen "$f.id1" }
	if { $id2_state == "normal" } { lappend boxen "$f.id2" }
    } else {
	if { $id2_state == "normal" } { lappend boxen "$f.id2" }
	if { $id1_state == "normal" } { lappend boxen "$f.id1" }	
    }
    # Need at least one item or we get errors...
    if { [llength $boxen] == 0 } { lappend boxen "$f.id1" }
    SetCurFrame $f $boxen
}

###########################################################################
proc FIJ_OK_Pressed { f io infile1 id1 infile2 id2 blocks min_overlap
		      word_length \
		      min_match use_band max_diag band_size max_mis \
		      yn sel_mode hidden_ops align_length \
		      use_filter filter_cutoff} {
    
    global CurContig
    global LREG
    global RREG
    global NGRec
    global gap5_defs

    set segment ""
    set active_tags {}

    if {[lorf_in_get $infile1] == 3} {
	set list1 [CreateAllContigList $io]
    } elseif {[lorf_in_get $infile1] == 4} {
	set gel_name [contig_id_gel $id1]
	set list1 "{$gel_name}"
    } else {
	set list1 [lorf_get_list $infile1]
    }

    if {[lorf_in_get $infile2] == 3} {
	set list2 [CreateAllContigList $io]
    } elseif {[lorf_in_get $infile2] == 4} {
	set gel_name [contig_id_gel $id2]
	set list2 "{$gel_name}"
    } else {
	set list2 [lorf_get_list $infile2]
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

    set word_length [lindex {? 12 8 4} [radiolist_get $word_length]]
    
    set max_prob [entrybox_get $max_diag]
    set max_align_length [entrybox_get $align_length]

    set fast_mode 0
    if {[radiolist_get $blocks] <= 2} {
        # Quick method
	set min_match [scalebox_get $min_match]
	set band_size [yes_no_get $use_band]

	if {[radiolist_get $blocks] == 1} {
	    set fast_mode 1
	}
    } else {
	# Sensitive method
	set min_match 0
	set band_size [scalebox_get $band_size]
    }

    set filter_words [yes_no_get $use_filter]
    if {$filter_words != 0} {
	set filter_words [$filter_cutoff get]
    }

    set min_overlap [scalebox_get $min_overlap]
    set max_mis [scalebox_get $max_mis]


    # Destroy dialog before showing plot to avoid window activation problems
    destroy $f


    #draw the contig comparator and dot plot
    ContigComparator $io 


    SetBusy
    find_internal_joins -io $io \
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
	    -max_display $max_align_length \
	    -fast_mode $fast_mode \
	    -filter_words $filter_words \
	    -contigs1 $list1 \
	    -contigs2 $list2

    ClearBusy
}

# Hidde data parameters box, used by Find Internal Joins and Calculate
# Consensus
# $w is the parent window
# $aname is a global array name in which to fill out the results
proc HiddenParametersDialog {w aname} {
    global gap5_defs
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

    label $w.label -text [keylget gap5_defs HIDDEN.CONF.LABEL]

    set ws [keylget gap5_defs HIDDEN.CONF.WINSIZE]
    scalebox $w.win_size \
	    -title [keylget ws NAME] \
	    -orient horizontal \
	    -to [keylget ws MAX] \
	    -from [keylget ws MIN] \
	    -default $data(conf_win_size) \
	    -width 5 \
	    -type CheckInt

    set md [keylget gap5_defs HIDDEN.CONF.MINCONF]
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

    label $w.label -text [keylget gap5_defs HIDDEN.UNC.LABEL]

    set ws [keylget gap5_defs HIDDEN.UNC.WINSIZE]
    scalebox $w.win_size \
	    -title [keylget ws NAME] \
	    -orient horizontal \
	    -to [keylget ws MAX] \
	    -from [keylget ws MIN] \
	    -default $data(base_win_size) \
	    -width 5 \
	    -type CheckInt

    set md [keylget gap5_defs HIDDEN.UNC.MAXDASH]
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
    set yn [keylget gap5_defs HIDDEN.MODE]
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
	-help_command "show_help gap5 {FIJ-Dialogue}" \

    pack $f.which $f.conf $f.unc $f.ok_cancel \
	-side top -fill both -expand 1
}

proc HiddenParametersDialogInit {aname} {
    global gap5_defs
    upvar #0 $aname data

    set data(which_mode)    [lindex {base conf} [keylget gap5_defs HIDDEN.MODE.VALUE]]
    set data(base_win_size) [keylget gap5_defs HIDDEN.UNC.WINSIZE.VALUE]
    set data(base_max_dash) [keylget gap5_defs HIDDEN.UNC.MAXDASH.VALUE]
    set data(conf_win_size) [keylget gap5_defs HIDDEN.CONF.WINSIZE.VALUE]
    set data(conf_min_conf) [keylget gap5_defs HIDDEN.CONF.MINCONF.VALUE]
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
