global env fak2_defs
if {![info exists env(FAKII)]} {
    set env(FAKII) [keylget fak2_defs BINARIES]
}

proc Fak2AssemblyDialogue {io t label } {
    global fak2_defs

    lorf_in  $t.infile  [keylget fak2_defs FAK2_ASSEMBLY.INFILE] \
	"" -bd 2 -relief groove
   
   # keylset fm FORMAT [keylget fak2_defs FAK2_ASSEMBLY.FORMAT]
   # set b1 [keylget fm FORMAT.BUTTON.1]
   # set b2 [keylget fm FORMAT.BUTTON.2]
   # radiolist $t.format \
	    -title [keylget fm FORMAT.NAME] \
	    -relief groove \
	    -bd 2 \
	    -default [keylget fm FORMAT.VALUE] \
	    -orient horizontal \
	    -buttons [format { { %s } { %s } } [list $b1] [list $b2] ]

    #Graph options
    frame $t.graph -bd 2 -relief groove
    label $t.graph.label -text "Create graph"

    keylset el E_LMT [keylget fak2_defs FAK2_ASSEMBLY.GRAPH.E_LMT]
    entrybox $t.graph.e_lmt \
	    -title [keylget el E_LMT.NAME] \
	    -default [keylget el E_LMT.VALUE] \
	    -width 15 \
	    -type "CheckFloatRange [keylget el E_LMT.MIN] \
	    [keylget el E_LMT.MAX]" \
	    -type "CheckFloat" \
	    -relief sunken \
	    -bd 2

    keylset ot O_THR [keylget fak2_defs FAK2_ASSEMBLY.GRAPH.O_THR]
    entrybox $t.graph.o_thr \
	    -title [keylget ot O_THR.NAME] \
	    -default [keylget ot O_THR.VALUE] \
	    -width 15 \
	    -type CheckFloat \
	    -relief sunken \
	    -bd 2 

    keylset dl D_LMT [keylget fak2_defs FAK2_ASSEMBLY.GRAPH.D_LMT]
    entrybox $t.graph.d_lmt \
	    -title [keylget dl D_LMT.NAME] \
	    -default [keylget dl D_LMT.VALUE] \
	    -width 15 \
	    -type "CheckFloatRange [keylget dl D_LMT.MIN] \
	    [keylget dl D_LMT.MAX]" \
	    -relief sunken \
	    -bd 2
    
    #Constraint options
    #create constraint file from exp file yes/no dialogue
    global $t.constraint
    keylset ex CONS [keylget fak2_defs FAK2_ASSEMBLY.GRAPH.CONS]
    yes_no $t.graph.constraint \
	-title   [keylget ex CONS.NAME] \
	-default [keylget ex CONS.VALUE] \
	-ycommand "set $t.constraint 1" \
	-ncommand "set $t.constraint 0" \
	-orient horizontal \
	-relief groove \
	-bd 2

    #Assembly options
    frame $t.assem -bd 2 -relief groove
    label $t.assem.label -text "Assembly"

    keylset an A_NUM [keylget fak2_defs FAK2_ASSEMBLY.ASSEM.A_NUM]
    entrybox $t.assem.a_num \
	    -title [keylget an A_NUM.NAME] \
	    -default [keylget an A_NUM.VALUE] \
	    -width 15 \
	    -type CheckInt \
	    -relief sunken \
	    -bd 2
   
    keylset er E_RATE [keylget fak2_defs FAK2_ASSEMBLY.ASSEM.E_RATE]
    entrybox $t.assem.e_rate \
	    -title [keylget er E_RATE.NAME] \
	    -default [keylget er E_RATE.VALUE] \
	    -width 15 \
	    -type CheckFloat \
	    -relief sunken \
	    -bd 2

    keylset ot O_THR [keylget fak2_defs FAK2_ASSEMBLY.ASSEM.O_THR]
    entrybox $t.assem.o_thr \
	    -title [keylget ot O_THR.NAME] \
	    -default [keylget ot O_THR.VALUE] \
	    -width 15 \
	    -type CheckFloat \
	    -relief sunken \
	    -bd 2 

    keylset dt D_THR [keylget fak2_defs FAK2_ASSEMBLY.ASSEM.D_THR]
    entrybox $t.assem.d_thr \
	    -title [keylget dt D_THR.NAME] \
	    -default [keylget dt D_THR.VALUE] \
	    -width 15 \
	    -type "CheckFloatRange [keylget dt D_THR.MIN] \
	    [keylget dt D_THR.MAX]" \
	    -relief sunken \
	    -bd 2 

    #destination directory
    keylset ex DEST [keylget fak2_defs FAK2_ASSEMBLY.ASSEM.DEST]
    entrybox $t.assem.dest_dir \
	    -title [keylget ex DEST.NAME ] \
	    -default [keylget ex DEST.VALUE] \
	    -width 15 \
	    -type CheckString \
	    -relief sunken \
	    -bd 2

    #constraints file
    #keylset l CONS [keylget fak2_defs FAK2_ASSEMBLY.ASSEM.CONS]
    #eval getFname \{$t.assem.constraints\} \{[keylget l CONS.NAME]\} load_optional {} \
#	[keylget l CONS.VALUE]
     
    #view options
    global $t.cb1 $t.cb2

    keylset v VIEW [keylget fak2_defs FAK2_ASSEMBLY.ASSEM.VIEW]
    set b1 [keylget v VIEW.BUTTON.1]
    set b2 [keylget v VIEW.BUTTON.2]
    set $t.cb1 [keylget v VIEW.VALUE]

    frame $t.assem.view
    checkbutton $t.assem.view.cb1 -text $b1 -variable $t.cb1
    checkbutton $t.assem.view.cb2 -text $b2 -variable $t.cb2

    pack $label -side top -fill both
    pack $t.infile -side top -fill x
    pack $t.graph $t.graph.label $t.graph.e_lmt $t.graph.o_thr $t.graph.d_lmt \
	    $t.graph.constraint -side top -fill x
    pack $t.assem $t.assem.label $t.assem.a_num $t.assem.e_rate $t.assem.o_thr \
	    $t.assem.d_thr $t.assem.dest_dir -side top -fill x
    pack $t.assem.view -side top -fill x
    pack $t.assem.view.cb1 $t.assem.view.cb2 -side left -fill x


}

proc Fak2Assembly { io } {
    global fak2_defs

    set t [keylget fak2_defs FAK2_ASSEMBLY1.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Perform FAKII assembly"

    label $t.label -text "Perform FAKII assembly" -bd 2 -relief groove
    Fak2AssemblyDialogue $io $t $t.label

    okcancelhelp $t.ok \
	    -ok_command "Fak2Assembly2 $io $t $t.infile $t.graph.e_lmt $t.graph.o_thr $t.graph.d_lmt $t.assem; destroy $t" \
	    -cancel_command "destroy $t" \
	    -help_command "show_help gap4 {Assembly-Perform FAKII assembly}" \
	    -bd 2 -relief groove

    pack $t.ok -side top -fill x

    
}

proc Fak2Assembly2 {io t infile e_lmt o_thr d_lmt assem} {
    global fak2_defs $t.cb1 $t.cb2 $t.constraint

    set bindir [keylget fak2_defs BINARIES]

    if {[set lin  [lorf_in_name $infile]]  == ""} {bell; return}
    set e_limit [entrybox_get $e_lmt]
    set o_threshold [entrybox_get $o_thr]
    set d_limit [entrybox_get $d_lmt]
    set assem_e_rate [entrybox_get $assem.e_rate]
    set assem_o_threshold [entrybox_get $assem.o_thr]
    set assem_d_threshold [entrybox_get $assem.d_thr]
    set a_num [entrybox_get $assem.a_num]
    #set constraints [getFname_in_name $assem.constraints]

    set dir [expandpath [entrybox_get $assem.dest_dir]]
    if {$dir == ""} {
	return
    }
    mkdir $dir

    #   if {[radiolist_get $informat] == 1} {
    #	set format "exp"
    #    } elseif {[radiolist_get $informat] == 2} {
    #	set format "fasta"
    #    } else {
    #	puts "Invalid format"
    #	return
    #    }
    #only allow experiment file format
    set format "exp"

    set layout [set $t.cb1]
    set multi [set $t.cb2]
    set binary [keylget fak2_defs FAK2_ASSEMBLY.BINARY]
    set g_binary [keylget fak2_defs FAK2_ASSEMBLY.G_BINARY]
    set const_b [keylget fak2_defs FAK2_ASSEMBLY.CONSTRAINT_B]
    set const_a [keylget fak2_defs FAK2_ASSEMBLY.CONSTRAINT_A]

    SetBusy

    vfuncheader "FAKII assembly"

    set g_error ""
    set c_error ""
    set a_error ""
    if {[set $t.constraint]} {
	#tout_pipe "create_graph $e_limit $o_threshold $d_limit -$format $lin\
	#	> $dir/$g_binary" "" 1
	#tout_pipe "create_exp_constraints $dir/$const_b $dir/$const_a \
	#	< $dir/$g_binary" "" 1
	#tout_pipe "assemble -$a_num -c$dir/$const_b $assem_e_rate \
	#	$assem_o_threshold $assem_d_threshold < $dir/$g_binary > \
	#$dir/$binary" "" 1
	catch {exec $bindir/create_graph $e_limit $o_threshold $d_limit \
		-$format $lin 2> $dir/graph_stderr > $dir/$g_binary} g_error
	catch {exec $bindir/create_exp_constraints $dir/$const_b $dir/$const_a\
		< $dir/$g_binary 2> $dir/constraints_stderr} c_error
	catch {exec $bindir/assemble -$a_num -c$dir/$const_b $assem_e_rate \
		$assem_o_threshold $assem_d_threshold < $dir/$g_binary \
		> $dir/$binary 2>$dir/assemble_stderr} a_error

    } else {
	#tout_pipe "create_graph $e_limit $o_threshold $d_limit \
	#	-$format $lin > x; assemble -$a_num $assem_e_rate \
	#	$assem_o_threshold $assem_d_threshold < x > $dir/$binary" "" 1

	catch {exec $bindir/create_graph $e_limit $o_threshold $d_limit \
		-$format $lin 2> $dir/graph_stderr > $dir/$g_binary} g_error
        catch {exec $bindir/assemble -$a_num $assem_e_rate $assem_o_threshold \
		$assem_d_threshold \
		< $dir/$g_binary > $dir/$binary 2>$dir/assemble_stderr} a_error
	
    }
    vmessage "--- graph output ---"
    set fs [open $dir/graph_stderr]
    while {![eof $fs]} {
	vmessage [gets $fs]
    }
    close $fs
    
    if {[set $t.constraint]} {
	vmessage "--- constraint output ---"
	set fs [open $dir/constraints_stderr]
	while {![eof $fs]} {
	    vmessage [gets $fs]
	}
	close $fs
    }	

    vmessage "--- assemble output ---"
    set fs [open $dir/assemble_stderr]
    while {![eof $fs]} {
	vmessage [gets $fs]
    }
    close $fs
    
    if {$g_error != ""} {
	ClearBusy
	return
    }
    if {$c_error != ""} {
	ClearBusy
	return
    }
    if {$a_error != ""} {
	ClearBusy
	return
    }

    if {$layout} {
	#tout_pipe "show_layout < $dir/$binary" "" 1
	set fs [open "|$bindir/show_layout < $dir/$binary \
		2> $dir/show_layout_stderr" r]
	while {![eof $fs]} {
	    vmessage [gets $fs]
	}
	close $fs
	set fs [open $dir/show_layout_stderr]
	while {![eof $fs]} {
	    verror ERR_WARN [gets $fs]
	}
	close $fs
    } 
    if {$multi} {
	#tout_pipe "show_multi < $dir/$binary" "" 1
	set fs [open "|$bindir/show_multi < $dir/$binary \
		2> $dir/show_multi_stderr" r]
	while {![eof $fs]} {
	    vmessage [gets $fs]
	}
	close $fs
	set fs [open $dir/show_multi_stderr]
	while {![eof $fs]} {
	    verror ERR_WARN [gets $fs]
	}
	close $fs
    }

    ClearBusy

    #set fs [open "|$bindir/create_graph $e_limit $o_threshold $d_limit \
	-$format $lin | assemble $assem_e_rate $assem_o_threshold \
	$assem_d_threshold | show_multi" r]

    #while {![eof $fs]} {
#	vmessage [gets $fs]
#	update idletasks
#    }
}

proc Fak2ImportDialogue { t } {
    global fak2_defs

    #use existing experiment files yes/no dialogue
    #global $t.flag
    #keylset ex FLAG [keylget fak2_defs FAK2_ASSEMBLY.IMPORT.FLAG]
    #yes_no $t.flag \
	-title   [keylget ex FLAG.NAME] \
	-default [keylget ex FLAG.VALUE] \
	-ycommand "set $t.flag 1" \
	-ncommand "set $t.flag 0" \
	-orient horizontal \
	-relief groove \
	-bd 2
   
    lorf_out $t.outfile [keylget fak2_defs FAK2_ASSEMBLY.OUTFILE] \
	"" -bd 2 -relief groove

    pack $t.label -side top -fill both
    pack $t.outfile -side top -fill x
}


proc Fak2Import {io } {
    global fak2_defs

    set t [keylget fak2_defs FAK2_ASSEMBLY2.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Import FAKII assembly data"

    label $t.label -text "Import FAKII assembly data" -bd 2 -relief groove
    pack $t.label -side top -fill both

    keylset ex INDIR [keylget fak2_defs FAK2_ASSEMBLY.INDIR]
    entrybox $t.indir \
	    -title [keylget ex INDIR.NAME ] \
	    -default [keylget fak2_defs FAK2_ASSEMBLY.ASSEM.DEST.VALUE] \
	    -width 15 \
	    -type CheckString \
	    -relief sunken \
	    -bd 2

    pack $t.indir -side top -fill both

    #select contig number or all
    keylset ex C_NUM [keylget fak2_defs FAK2_ASSEMBLY.IMPORT.C_NUM]
    entrybox $t.c_num \
	    -title [keylget ex C_NUM.NAME ] \
	    -width 15 \
	    -type CheckInt \
	    -relief sunken \
	    -bd 2
   
    lorf_in $t.infile [keylget fak2_defs FAK2_ASSEMBLY.IMPORT.INFILE] \
	    "{entrybox_configure $t.c_num -state disabled}
             {entrybox_configure $t.c_num -state disabled}
	     {entrybox_configure $t.c_num -state disabled}
	     {entrybox_configure $t.c_num -state normal}
	" -bd 2 -relief groove
    pack $t.infile $t.c_num -side top -fill x


    Fak2ImportDialogue $t 

    okcancelhelp $t.ok \
	-ok_command "Fak2Import2 $io $t $t.indir $t.infile $t.c_num $t.outfile; \
	destroy $t" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Assembly-Import FAKII assembly}" \
	-bd 2 -relief groove

    pack $t.ok -side top -fill x
}


proc Fak2Import2 {io t indir lorf c_num outfile} {
    global fak2_defs gap_defs

    set bindir [keylget fak2_defs BINARIES]
    set mism 100
    set display 0
    set contig ""

    #lorf is NULL if do both assembly and import together
    if {[string compare $lorf "NULL"] == 0} {
	set contig 0
    } else {
	if {[lorf_in_get $lorf] == 4} {
	    #single contig
	    if {[set contig [entrybox_get $c_num]] == ""} {bell; return}
	} elseif {[lorf_in_get $lorf] == 3} {
	    #all contings
	    set contig 0
	} else {
	    if {[set l [lorf_get_list $lorf]] == ""} {bell; return}
	    set contig $l
	}
    }

    if {[set lout [lorf_out_name $outfile]] == ""} {bell; return}

    set dir [expandpath [entrybox_get $indir]]
    if {$dir == ""} {
	return
    }
    set infile $dir/[keylget fak2_defs FAK2_ASSEMBLY.FOFN]

    set orig_dir [exec pwd]

    if {![quit_displays $io "auto_assemble"]} {
        # Someone's too busy to shutdown?
        return
    }

    SetBusy

    vfuncheader "FAKII import assembly"

    #tout_pipe "write_exp_file -c \"$contig\" -o $dir -f $infile \
    #< $dir/[keylget fak2_defs FAK2_ASSEMBLY.BINARY]" "" 1

    set w_error ""
    if {[catch {set fs [open "|$bindir/write_exp_file -c \"$contig\" -o $dir \
		-f $infile < $dir/[keylget fak2_defs FAK2_ASSEMBLY.BINARY] \
		2>$dir/exp_stderr" r]} w_error]} {
	bell
	verror ERR_WARN "FAKII import" $w_error
	ClearBusy
	return
    }
    while {![eof $fs]} {
	vmessage [gets $fs]
    }
    close $fs
    set fs [open $dir/exp_stderr]
    while {![eof $fs]} {
	verror ERR_WARN "FAKII import" [gets $fs]
    }
    close $fs
    
    if {[catch {set f [open "$infile"]}]} {
	bell
	ClearBusy
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
	bell; ClearBusy; return 
    }
	
    #change directory to where experiment files are and assemble them into
    #gap
    cd $dir
    ListCreate2 $lout [assemble_direct \
			   -io $io \
			   -files $lin \
			   -max_pmismatch $mism \
			   -output_mode $display]

    ClearBusy

    #change back to original directory
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

proc Fak2AssemblyImport { io } {
    global fak2_defs

    set t [keylget fak2_defs FAK2_ASSEMBLY3.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Perform and import FAKII assembly"

    label $t.label -text "Perform and import FAKII assembly" -bd 2 -relief groove	

    Fak2AssemblyDialogue $io $t $t.label
    frame $t.import
    label $t.import.label -text "Import FAKII assembly data" 
    pack $t.import.label -side top -fill both
    Fak2ImportDialogue $t.import

    pack $t.import -side top -fill x
    okcancelhelp $t.ok \
	    -ok_command "Fak2AssemblyImport2 $io $t; destroy $t" \
	    -cancel_command "destroy $t" \
	    -help_command "show_help gap4 {Assembly-Perform and import FAKII assembly}"  \
	    -bd 2 -relief groove

    pack $t.ok -side top -fill x
}

proc Fak2AssemblyImport2 {io t} {

    Fak2Assembly2 $io $t $t.infile $t.graph.e_lmt $t.graph.o_thr $t.graph.d_lmt $t.assem
    Fak2Import2 $io $t $t.assem.dest_dir "NULL" "NULL" $t.import.outfile
}
