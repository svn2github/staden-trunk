proc HuangAssembly { io option} {
    global cap2_defs

    #assemble only
    if {$option == 1} {
	set t [keylget cap2_defs HUANG_ASSEMBLY1.WIN]
    } else {
	#assembly and import
	set t [keylget cap2_defs HUANG_ASSEMBLY3.WIN]
    }
    if {[xtoplevel $t -resizable 0] == ""} return

    if {$option == 1} {
	label $t.label -text "Perform CAP2 assembly" -bd 2 -relief groove
	wm title $t "Perform CAP2 assembly"
    } else {
	label $t.label -text "Perform and import CAP2 assembly" -bd 2 -relief groove	
	wm title $t "Perform ans import CAP2 assembly"
    }
    lorf_in  $t.infile  [keylget cap2_defs HUANG_ASSEMBLY.INFILE] \
	"" -bd 2 -relief groove

#    keylset fm FORMAT [keylget cap2_defs HUANG_ASSEMBLY.FORMAT]
#    set b1 [keylget fm FORMAT.BUTTON.1]
#    set b2 [keylget fm FORMAT.BUTTON.2]
#    radiolist $t.format \
	    -title [keylget fm FORMAT.NAME] \
	    -relief groove \
	    -bd 2 \
	    -default [keylget fm FORMAT.VALUE] \
	    -orient horizontal \
	    -buttons [format { { %s } { %s } } [list $b1] [list $b2] ]
    
    keylset ex DEST [keylget cap2_defs HUANG_ASSEMBLY.DEST]
    entrybox $t.dest_dir \
	    -title [keylget ex DEST.NAME ] \
	    -default [keylget ex DEST.VALUE] \
	    -width 15 \
	    -type CheckString \
	    -relief groove \
	    -bd 2

    keylset ex REP [keylget cap2_defs HUANG_ASSEMBLY.REPEAT]
    yes_no $t.repeats \
	-title   [keylget ex REP.NAME] \
	-default [keylget ex REP.VALUE] \
	-ycommand "set $t.repeat 1" \
	-ncommand "set $t.repeat 0" \
	-orient horizontal \
	-relief groove \
	-bd 2

    #assemble and import
    if {$option == 2} {
	lorf_out $t.outfile [keylget cap2_defs HUANG_ASSEMBLY.OUTFILE] \
	"" -bd 2 -relief groove
	okcancelhelp $t.ok \
		-ok_command "HuangAssemblyImport $io $t $t.infile $t.outfile \
		$t.dest_dir \[set $t.repeat\]; destroy $t" \
		-cancel_command "destroy $t" \
		-help_command "show_help gap4 {Assembly-Perform and import CAP2 assembly}" \
		-bd 2 -relief groove
    } else {

	okcancelhelp $t.ok \
		-ok_command "HuangAssembly2 $io $t $t.infile \
		$t.dest_dir \[set $t.repeat\]; destroy $t" \
		-cancel_command "destroy $t" \
		-help_command "show_help gap4 {Assembly-Perform CAP2 assembly}" \
		-bd 2 -relief groove
    }

    pack $t.label -side top -fill both
    pack $t.infile $t.dest_dir $t.repeats -side top -fill x
    if {$option == 2} {
	pack $t.outfile -side top -fill x
    }
    pack $t.ok -side top -fill x

}

proc HuangAssembly2 { io t infile dest_dir repeat} {

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
    tout_pipe "cap2_s -$format $lin -out $dir -r" "" 1

    ClearBusy
}

proc HuangAssemblyImport { io t infile outfile dest_dir repeat} {

    HuangAssembly2 $io $t $infile $dest_dir $repeat
    HuangImport2 $io $t $dest_dir $outfile
}

proc HuangImport { io } {
    global cap2_defs

    set t [keylget cap2_defs HUANG_ASSEMBLY2.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Import CAP2 assembly data"

    label $t.label -text "Import CAP2 assembly data" -bd 2 -relief groove

    pack $t.label -side top -fill both

    keylset ex INDIR [keylget cap2_defs HUANG_ASSEMBLY.INDIR]
    entrybox $t.indir \
	    -title [keylget ex INDIR.NAME ] \
	    -default [keylget cap2_defs HUANG_ASSEMBLY.DEST.VALUE] \
	    -width 15 \
	    -type CheckString \
	    -relief groove \
	    -bd 2

    lorf_out $t.outfile [keylget cap2_defs HUANG_ASSEMBLY.OUTFILE] \
	"" -bd 2 -relief groove

    okcancelhelp $t.ok \
	-ok_command "HuangImport2 $io $t $t.indir $t.outfile; destroy $t" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Assembly-Import CAP2 assembly}" \
	-bd 2 -relief groove

    pack $t.label -side top -fill both
    pack $t.indir $t.outfile $t.ok -side top -fill x
    

}

proc HuangImport2 {io t indir outfile} {
    global cap2_defs gap_defs

    set mism 100
    set display 0

    if {[set lout [lorf_out_name $outfile]] == ""} {bell; return}
    set dir [entrybox_get $indir]
    set infile $dir/[keylget cap2_defs HUANG_ASSEMBLY.FOFN]

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
	-output_mode $display]
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
