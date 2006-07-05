#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#
# Returns a string describing the database
# Also (as this is a common place) bumps at maxseq as needed.
#
proc DirectAssembly {io} {
    global gap_defs

    set t [keylget gap_defs DIRECT_ASSEMBLY.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Directed assembly"

    set d [keylget gap_defs DIRECT_ASSEMBLY.DISPLAY]
    frame $t.display -bd 2 -relief groove
    checkbutton $t.display.c \
	-text [keylget d NAME] \
    	-variable $t.Display_alignments
    pack $t.display.c -side left
    set $t.Display_alignments [keylget d VALUE]

    set d [keylget gap_defs DIRECT_ASSEMBLY.IGNORE_VEC]
    frame $t.ignore_vec -bd 2 -relief groove
    checkbutton $t.ignore_vec.c \
	-text [keylget d NAME] \
    	-variable $t.Ignore_vec
    pack $t.ignore_vec.c -side left
    set $t.Ignore_vec [keylget d VALUE]

    set mm [keylget gap_defs DIRECT_ASSEMBLY.MAXMIS]
    scalebox $t.mism \
	-orient horizontal \
	-width 5 \
	-title [keylget mm NAME] \
	-from [keylget mm MIN] \
	-to [keylget mm MAX] \
	-resolution [keylget mm RES] \
	-default [keylget mm VALUE] \
	-type CheckFloat

    lorf_in  $t.infile  [keylget gap_defs DIRECT_ASSEMBLY.INFILE] \
	"{#list} {#fofn} {#selection}" -bd 2 -relief groove
    lorf_out $t.outfile [keylget gap_defs DIRECT_ASSEMBLY.OUTFILE] \
	"" -bd 2 -relief groove

    #select alignment failure mode
    keylset fm FAILURE [keylget gap_defs DIRECT_ASSEMBLY.FAILURE]
    set b1 [keylget fm FAILURE.BUTTON.1]
    set b2 [keylget fm FAILURE.BUTTON.2]
    radiolist $t.failure_mode \
	    -title [keylget fm FAILURE.NAME] \
	    -relief groove \
	    -bd 2 \
	    -default [keylget fm FAILURE.VALUE] \
	    -buttons [format { { %s } { %s } } [list $b1] [list $b2] ]

    okcancelhelp $t.ok \
	-ok_command "DirectAssembly2 $io $t $t.infile $t.outfile \
		         $t.mism \[set $t.Display_alignments\] \
			 $t.failure_mode \[set $t.Ignore_vec\]" \
	-cancel_command "destroy $t; catch {unset $t.Display_alignments}" \
	-help_command "show_help gap4 {Assembly-Directed}" \
	-bd 2 -relief groove

    pack $t.display $t.ignore_vec -side top -fill both
    pack $t.mism $t.infile $t.outfile $t.failure_mode $t.ok -side top -fill x
}

proc DirectAssembly2 {io t infile outfile mism_w display failure_mode
		      ignore_vec} {
    global gap_defs

    #special case for a single file name
    if {[lorf_in_get $infile] == 3} {
	if {[set lin [lorf_in_name $infile]]  == ""} {bell; return}
    } else {
	if {[set lin [lorf_get_list $infile]]  == ""} {bell; return}
    }
    if {[set lout [lorf_out_name $outfile]] == ""} {bell; return}
    if {[set mism [scalebox_get $mism_w]]   == ""} {bell; return}

    if {![quit_displays $io "auto_assemble"]} {
        # Someone's too busy to shutdown?
        return
    }

    # Expand $inlist to include sub-files if input is a hashed archive.
    set lin [expand_hash_archives $lin]

    SetBusy

    destroy $t
    catch {unset $t.Display_alignments}

#    update idletasks
    set al [expr {$mism?1:0}]
    ListCreate2 $lout [assemble_direct \
        -align $al \
	-io $io \
	-files $lin \
	-max_pmismatch $mism \
	-output_mode $display \
	-enter_failures [expr [radiolist_get $failure_mode]-1] \
	-ignore_vec $ignore_vec]
    ClearBusy

    reset_exp_path

    if {[lorf_out_get $outfile] == 2} {
	lorf_out_save $lout
    }

    #check that the database has contigs in it!
    if {[db_info num_contigs $io] > 0} {
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
}
