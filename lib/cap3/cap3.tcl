proc HuangAssembly_3 { io option} {
    global cap3_defs

    #assemble only
    if {$option == 1} {
	set t [keylget cap3_defs HUANG_ASSEMBLY1.WIN]
    } else {
	#assembly and import
	set t [keylget cap3_defs HUANG_ASSEMBLY3.WIN]
    }
    if {[xtoplevel $t -resizable 0] == ""} return

    if {$option == 1} {
	label $t.label -text "Perform CAP3 assembly" -bd 2 -relief groove
	wm title $t "Perform CAP3 assembly"
    } else {
	label $t.label -text "Perform and import CAP3 assembly" -bd 2 -relief groove	
	wm title $t "Perform CAP3 and import assembly"
    }
    lorf_in  $t.infile  [keylget cap3_defs HUANG_ASSEMBLY.INFILE] \
	"" -bd 2 -relief groove

#    keylset fm FORMAT [keylget cap3_defs HUANG_ASSEMBLY.FORMAT]
#    set b1 [keylget fm FORMAT.BUTTON.1]
#    set b2 [keylget fm FORMAT.BUTTON.2]
#    radiolist $t.format \
	    -title [keylget fm FORMAT.NAME] \
	    -relief groove \
	    -bd 2 \
	    -default [keylget fm FORMAT.VALUE] \
	    -orient horizontal \
	    -buttons [format { { %s } { %s } } [list $b1] [list $b2] ]

    #Constraint options
    #create constraint file from exp file yes/no dialogue
    global $t.constraint
    keylset ex CONS [keylget cap3_defs HUANG_ASSEMBLY.CONS]
    yes_no $t.constraint \
	-title   [keylget ex CONS.NAME] \
	-default [keylget ex CONS.VALUE] \
	-ycommand "set $t.constraint 1" \
	-ncommand "set $t.constraint 0" \
	-orient horizontal \
	-relief groove \
	-bd 2

    keylset ex DEST [keylget cap3_defs HUANG_ASSEMBLY.DEST]
    entrybox $t.dest_dir \
	    -title [keylget ex DEST.NAME ] \
	    -default [keylget ex DEST.VALUE] \
	    -width 15 \
	    -type CheckString \
	    -relief groove \
	    -bd 2

    #option values
    keylset be BAND_EXP [keylget cap3_defs HUANG_ASSEMBLY.BAND_EXP]
    entrybox $t.band_exp \
	-title "[keylget be BAND_EXP.NAME] ( > [keylget be BAND_EXP.MIN])"\
	    -default [keylget be BAND_EXP.VALUE] \
	    -width 5 \
	    -type "CheckIntMin [keylget be BAND_EXP.MIN]"


    keylset be QUAL_DIFF [keylget cap3_defs HUANG_ASSEMBLY.QUAL_DIFF]
    entrybox $t.qual_diff \
	-title "[keylget be QUAL_DIFF.NAME] ([keylget be QUAL_DIFF.MIN] to [keylget be QUAL_DIFF.MAX])"\
	    -default [keylget be QUAL_DIFF.VALUE] \
	    -width 5 \
	    -type "CheckIntRange [keylget be QUAL_DIFF.MIN] [keylget be QUAL_DIFF.MAX]"

    keylset be QUAL_CLIP [keylget cap3_defs HUANG_ASSEMBLY.QUAL_CLIP]
    entrybox $t.qual_clip \
	-title "[keylget be QUAL_CLIP.NAME] ([keylget be QUAL_CLIP.MIN] to [keylget be QUAL_DIFF.MAX])"\
	    -default [keylget be QUAL_CLIP.VALUE] \
	    -width 5 \
	    -type "CheckIntRange [keylget be QUAL_CLIP.MIN] [keylget be QUAL_CLIP.MAX]"

    keylset be QSCORE [keylget cap3_defs HUANG_ASSEMBLY.QSCORE]
    entrybox $t.qscore \
	-title "[keylget be QSCORE.NAME] ( > [keylget be QSCORE.MIN])"\
	    -default [keylget be QSCORE.VALUE] \
	    -width 5 \
	    -type "CheckIntMin [keylget be QSCORE.MIN]"

    keylset be CLEARANCE [keylget cap3_defs HUANG_ASSEMBLY.CLEARANCE]
    entrybox $t.clearance \
	-title "[keylget be CLEARANCE.NAME] ( > [keylget be CLEARANCE.MIN])"\
	    -default [keylget be CLEARANCE.VALUE] \
	    -width 5 \
	    -type "CheckIntMin [keylget be CLEARANCE.MIN]"

    keylset be GAP_PENALTY [keylget cap3_defs HUANG_ASSEMBLY.GAP_PENALTY]
    entrybox $t.gap_penalty \
	-title "[keylget be GAP_PENALTY.NAME] ( > [keylget be GAP_PENALTY.MIN])"\
	    -default [keylget be GAP_PENALTY.VALUE] \
	    -width 5 \
	    -type "CheckIntMin [keylget be GAP_PENALTY.MIN]"

    keylset be MATCH_SCORE [keylget cap3_defs HUANG_ASSEMBLY.MATCH_SCORE]
    entrybox $t.match_score \
	-title "[keylget be MATCH_SCORE.NAME] ( > [keylget be MATCH_SCORE.MIN])"\
	    -default [keylget be MATCH_SCORE.VALUE] \
	    -width 5 \
	    -type "CheckIntMin [keylget be MATCH_SCORE.MIN]"

    keylset be MISMATCH_SCORE [keylget cap3_defs HUANG_ASSEMBLY.MISMATCH_SCORE]
    entrybox $t.mismatch_score \
	-title "[keylget be MISMATCH_SCORE.NAME] ( < [keylget be MISMATCH_SCORE.MAX])"\
	-default [keylget be MISMATCH_SCORE.VALUE] \
	-width 5 \
	-type "CheckIntMax [keylget be MISMATCH_SCORE.MAX]"

    keylset be OVERLAP_LENGTH [keylget cap3_defs HUANG_ASSEMBLY.OVERLAP_LENGTH]
    entrybox $t.overlap_length \
	-title "[keylget be OVERLAP_LENGTH.NAME] ( > [keylget be OVERLAP_LENGTH.MIN])"\
	-default [keylget be OVERLAP_LENGTH.VALUE] \
	-width 5 \
	-type "CheckIntMin [keylget be OVERLAP_LENGTH.MIN]"

    keylset be OVERLAP_IDENTITY [keylget cap3_defs HUANG_ASSEMBLY.OVERLAP_IDENTITY]
    entrybox $t.overlap_identity \
	-title "[keylget be OVERLAP_IDENTITY.NAME] ( > [keylget be OVERLAP_IDENTITY.MIN])"\
	-default [keylget be OVERLAP_IDENTITY.VALUE] \
	-width 5 \
	-type "CheckIntMin [keylget be OVERLAP_IDENTITY.MIN]"

    keylset be OVERLAP_SIMILARITY [keylget cap3_defs HUANG_ASSEMBLY.OVERLAP_SIMILARITY]
    entrybox $t.overlap_similarity \
	-title "[keylget be OVERLAP_SIMILARITY.NAME] ( > [keylget be OVERLAP_SIMILARITY.MIN])"\
	-default [keylget be OVERLAP_SIMILARITY.VALUE] \
	-width 5 \
	-type "CheckIntMin [keylget be OVERLAP_SIMILARITY.MIN]"
    
    keylset be MIN_CORRECTION [keylget cap3_defs HUANG_ASSEMBLY.MIN_CORRECTION]
    entrybox $t.min_correction \
	-title "[keylget be MIN_CORRECTION.NAME] ( > [keylget be MIN_CORRECTION.MIN])"\
	-default [keylget be MIN_CORRECTION.VALUE] \
	-width 5 \
	-type "CheckIntMin [keylget be MIN_CORRECTION.MIN]"

    keylset be MIN_LINKING [keylget cap3_defs HUANG_ASSEMBLY.MIN_LINKING]
    entrybox $t.min_linking \
	-title "[keylget be MIN_LINKING.NAME] ( > [keylget be MIN_LINKING.MIN])"\
	-default [keylget be MIN_LINKING.VALUE] \
	-width 5 \
	-type "CheckIntMin [keylget be MIN_LINKING.MIN]"

    keylset be PREFIX [keylget cap3_defs HUANG_ASSEMBLY.PREFIX]
    entrybox $t.prefix \
	-title [keylget be PREFIX.NAME] \
	-default [keylget be PREFIX.VALUE] \
	-width 5 \
	-type CheckString

    #assemble and import
    if {$option == 2} {
	lorf_out $t.outfile [keylget cap3_defs HUANG_ASSEMBLY.OUTFILE] \
	"" -bd 2 -relief groove
	okcancelhelp $t.ok \
		-ok_command "HuangAssemblyImport_3 $io $t $t.infile $t.outfile \
		$t.dest_dir \[set $t.constraint\] \[entrybox_get $t.band_exp\] \[entrybox_get $t.qual_diff\] \[entrybox_get $t.qual_clip\] \[entrybox_get $t.qscore\] \[entrybox_get $t.clearance\] \[entrybox_get $t.gap_penalty\] \[entrybox_get $t.match_score\] \[entrybox_get $t.mismatch_score\] \[entrybox_get $t.overlap_length\] \[entrybox_get $t.overlap_identity\] \[entrybox_get $t.overlap_similarity\] \[entrybox_get $t.min_correction\] \[entrybox_get $t.min_linking\] \[entrybox_get $t.prefix\]; destroy $t" \
		-cancel_command "destroy $t" \
		-help_command "show_help gap4 {Assembly-Perform and import CAP3 assembly}" \
		-bd 2 -relief groove
    } else {

	okcancelhelp $t.ok \
		-ok_command "HuangAssembly2_3 $io $t $t.infile \
		$t.dest_dir \[set $t.constraint\] \[entrybox_get $t.band_exp\] \[entrybox_get $t.qual_diff\] \[entrybox_get $t.qual_clip\] \[entrybox_get $t.qscore\] \[entrybox_get $t.clearance\] \[entrybox_get $t.gap_penalty\] \[entrybox_get $t.match_score\] \[entrybox_get $t.mismatch_score\] \[entrybox_get $t.overlap_length\] \[entrybox_get $t.overlap_identity\] \[entrybox_get $t.overlap_similarity\] \[entrybox_get $t.min_correction\] \[entrybox_get $t.min_linking\] \[entrybox_get $t.prefix\]; destroy $t" \
		-cancel_command "destroy $t" \
		-help_command "show_help gap4 {Assembly-Perform CAP3 assembly}" \
		-bd 2 -relief groove
    }

    pack $t.label -side top -fill both
    pack $t.infile $t.constraint $t.dest_dir $t.band_exp $t.qual_diff $t.qual_clip $t.qscore $t.clearance $t.gap_penalty $t.match_score $t.mismatch_score $t.overlap_length $t.overlap_identity $t.overlap_similarity $t.min_correction $t.min_linking $t.prefix -side top -fill x    

    if {$option == 2} {
	pack $t.outfile -side top -fill x
    }
    pack $t.ok -side top -fill x

}

proc HuangAssembly2_3 { io t infile dest_dir constraint band_exp qual_diff qual_clip qscore clearance gap_penalty match_score mismatch_score overlap_length overlap_identity overlap_similarity min_correction min_linking prefix} {

    if {[set lin  [lorf_in_name $infile]]  == ""} {bell; return}
    set dir [entrybox_get $dest_dir]

#    if {[radiolist_get $informat] == 1} {
#	set format "exp"
#    } elseif {[radiolist_get $informat] == 2} {
#	set format "fasta"
#    } else {
#	puts "Invalid format"
#	return
#    }
    set format "exp"
    SetBusy

#    set fs [open "|/home7/kfs/assembly/work/cap2/sgi-binaries/cap2 -$format $lin -out $dir -r" r]
#
#    while {![eof $fs]} {
#	vmessage [gets $fs]
#	update idletasks
#    }
    if {$constraint} {
	tout_pipe "cap3_create_exp_constraints $lin" "" 1
    } else {
	#need to remove .con file to prevent doing any constraints when I
	#don't want to.
	file delete $lin.con
    }

    tout_pipe "cap3_s $lin -a $band_exp -b $qual_diff -c $qual_clip -d $qscore -e $clearance -g $gap_penalty -m $match_score -n $mismatch_score -o $overlap_length -p $overlap_identity -s $overlap_similarity -u $min_correction -v $min_linking -x $prefix -f exp -y $dir" "" 1

    ClearBusy
}

proc HuangAssemblyImport_3 { io t infile outfile dest_dir constraint band_exp qual_diff qual_clip qscore clearance gap_penalty match_score mismatch_score overlap_length overlap_identity overlap_similarity min_correction min_linking prefix} {

    HuangAssembly2_3 $io $t $infile $dest_dir $constraint $band_exp $qual_diff $qual_clip $qscore $clearance $gap_penalty $match_score $mismatch_score $overlap_length $overlap_identity $overlap_similarity $min_correction $min_linking $prefix
    HuangImport2_3 $io $t $dest_dir $outfile
}

proc HuangImport_3 { io } {
    global cap3_defs

    set t [keylget cap3_defs HUANG_ASSEMBLY2.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Import CAP3 assembly data"

    label $t.label -text "Import CAP3 assembly data" -bd 2 -relief groove

    pack $t.label -side top -fill both

    keylset ex INDIR [keylget cap3_defs HUANG_ASSEMBLY.INDIR]
    entrybox $t.indir \
	    -title [keylget ex INDIR.NAME ] \
	    -default [keylget cap3_defs HUANG_ASSEMBLY.DEST.VALUE] \
	    -width 15 \
	    -type CheckString \
	    -relief groove \
	    -bd 2

    lorf_out $t.outfile [keylget cap3_defs HUANG_ASSEMBLY.OUTFILE] \
	"" -bd 2 -relief groove

    okcancelhelp $t.ok \
	-ok_command "HuangImport2_3 $io $t $t.indir $t.outfile; destroy $t" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Assembly-Import CAP3 assembly}" \
	-bd 2 -relief groove

    pack $t.label -side top -fill both
    pack $t.indir $t.outfile $t.ok -side top -fill x
    

}

proc HuangImport2_3 {io t indir outfile} {
    global cap3_defs gap_defs

    set mism 100
    set display 0

    if {[set lout [lorf_out_name $outfile]] == ""} {bell; return}
    set dir [entrybox_get $indir]
    set infile $dir/[keylget cap3_defs HUANG_ASSEMBLY.FOFN]

    if {[catch {set f [open "$infile"]}]} {
	bell
	return
    }
    set lin ""
    while {[gets $f line] >= 0} {
	lappend lin $line
    }

    if {$lin == ""} {
	tk_messageBox \
		-icon error \
		-title "Empty list" \
		-message "'$infile' contains zero names. Aborting." \
		-type ok \
		-parent $t

	bell; return 
    }
    set orig_dir [exec pwd]
    cd $dir

    if {![quit_displays $io "auto_assemble"]} {
        # Someone's too busy to shutdown?
        return
    }
    SetBusy

#    update idletasks
    ListCreate2 $lout [assemble_direct \
	-io $io \
	-files $lin \
	-max_pmismatch $mism \
	-output_mode $display\
	-align 0]
    ClearBusy

    cd $orig_dir
    if {[lorf_out_get $outfile] == 2} {
	lorf_out_save $lout
    }

    ActivateMenu_Open
    InitContigGlobals $io

    #draw the contig selector if it does not already exist
    set cs_win [keylget gap_defs CONTIG_SEL.WIN]
    if {![winfo exists $cs_win]} {
        ContigSelector $io
    } else {
        ContigInitReg $io
        raise $cs_win
    }
}
