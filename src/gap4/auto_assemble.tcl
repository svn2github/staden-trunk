#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
##############################################################################
#display mode options
proc AADispMode { f } {
    global gap_defs

    ###########################################################################
    #select display mode

    keylset dm DISPMODE [keylget gap_defs AUTO_ASSEMBLE.DISPMODE]
    set b1 [keylget dm DISPMODE.BUTTON.1]
    set b2 [keylget dm DISPMODE.BUTTON.2]
    set b3 [keylget dm DISPMODE.BUTTON.3]
    set b4 [keylget dm DISPMODE.BUTTON.4]

    eval radiolist $f.disp_mode \
	    -title {[keylget dm DISPMODE.NAME]} \
	    -bd 2 \
	    -relief groove \
	    -default [keylget dm DISPMODE.VALUE] \
	    -buttons  \{ \{\"$b1\"\} \{\"$b2\"\} \{\"$b3\"\} \{\"$b4\"\}\}
    
}

##############################################################################
#matching scaleboxes 
proc AAMatch { f } {
    global gap_defs

    ###########################################################################
    #match scales

    frame $f.match -relief groove -bd 2

    keylset mm MINMATCH [keylget gap_defs AUTO_ASSEMBLE.MINMATCH]

    scalebox $f.match.min_match \
	    -title [keylget mm MINMATCH.NAME] \
	    -orient horizontal \
	    -to [keylget mm MINMATCH.MAX] \
	    -from [keylget mm MINMATCH.MIN] \
	    -default [keylget mm MINMATCH.VALUE] \
	    -width 5 \
	    -type CheckInt
    
    keylset mp MAXPADS [keylget gap_defs AUTO_ASSEMBLE.MAXPADS]
    scalebox $f.match.max_pads \
	    -title [keylget mp MAXPADS.NAME]\
	    -orient horizontal \
	    -to [keylget mp MAXPADS.MAX] \
	    -from [keylget mp MAXPADS.MIN] \
	    -default [keylget mp MAXPADS.VALUE] \
	    -width 5 \
	    -type CheckInt

    keylset mat MISMATCH [keylget gap_defs AUTO_ASSEMBLE.MISMATCH]
    scalebox $f.match.max_mis \
	    -title [keylget mat MISMATCH.NAME]\
	    -orient horizontal \
	    -to [keylget mat MISMATCH.MAX]\
	    -from [keylget mat MISMATCH.MIN] \
	    -default [keylget mat MISMATCH.VALUE] \
	    -resolution [keylget mat MISMATCH.RES] \
	    -width 5 \
	    -type CheckFloat

    pack $f.match.min_match -fill x
    pack $f.match.max_pads -fill x
    pack $f.match.max_mis -fill x
}

proc AAMasking { parent f } {
    global gap_defs

    frame $f -bd 2 -relief groove
    keylset am MASKING [keylget gap_defs AUTO_ASSEMBLE.MASKING]

    button $f.but \
	    -text "Select tags" \
	    -command "TagDialog AUTO_ASSEMBLE.TAGS\
			$parent[keylget gap_defs SELECT_TAGS.WIN] {}"

    yes_no $f.yn \
	    -title [keylget am MASKING.NAME] \
	    -orient horizontal \
	    -default [keylget am MASKING.VALUE] \
	    -ycommand "$f.but configure -state normal; \
		       SetDefaultTags AUTO_ASSEMBLE.TAGS
" \
	    -ncommand "$f.but configure -state disabled"

    pack $f.yn -side left
    pack $f.but -side right

}

##############################################################################
#normal shotgun assembly
proc NormalShotgun { f io } {
    global gap_defs

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Normal shotgun assembly"

    ###########################################################################
    #apply masking
    AAMasking $f $f.masking

    ###########################################################################
    #display mode
    AADispMode $f

    #matching scaleboxes
    AAMatch $f

    #create file of filenames entry windows
    lorf_in $f.infile [keylget gap_defs AUTO_ASSEMBLE.INFILE] \
	"{#list} {#fofn} {#selection} {#biolims}" -bd 2 -relief groove

    #save failures as file or list
    lorf_out $f.fails [keylget gap_defs AUTO_ASSEMBLE.FAILS] \
	"" -bd 2 -relief groove

    ###########################################################################
    #select alignment failure mode
    keylset fm FAILURE [keylget gap_defs AUTO_ASSEMBLE.FAILURE]
    set b1 [keylget fm FAILURE.BUTTON.1]
    set b2 [keylget fm FAILURE.BUTTON.2]
    radiolist $f.failure_mode \
	    -title [keylget fm FAILURE.NAME] \
	    -relief groove \
	    -bd 2 \
	    -default [keylget fm FAILURE.VALUE] \
	    -buttons [format { { %s } { %s } } [list $b1] [list $b2] ]
    
    ###########################################################################
    #permit joins
    keylset pj PERMITJOINS [keylget gap_defs AUTO_ASSEMBLE.PERMITJOINS]
    yes_no $f.permit_joins \
	    -title [keylget pj PERMITJOINS.NAME] \
	    -orient horizontal \
	    -relief groove \
	    -bd 2 \
	    -default [keylget pj PERMITJOINS.VALUE]

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "OK_Pressed1 $io $f $f.masking.yn $f.disp_mode \
            $f.match $f.infile $f.fails $f.permit_joins $f.failure_mode" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Assembly-Shot}" \
	    -bd 2 \
	    -relief groove
    
    ###########################################################################
    pack $f.masking -fill x
    pack $f.disp_mode -fill x
    pack $f.match -fill x
    pack $f.infile -fill x
    pack $f.fails -fill x
    pack $f.permit_joins -fill x
    pack $f.failure_mode -fill x
    pack $f.ok_cancel -fill x
    
}

##############################################################################
proc OK_Pressed1 { io f masking disp_mode match infile fails \
	permit_joins failure_mode} {
    global gap_defs

    #masking determines fortran variable IOPT to be 1 or 2
    set masking [yes_no_get $masking]
    if {$masking} {
	set active_tags [GetDefaultTags AUTO_ASSEMBLE.TAGS]
    } else {
	set active_tags {}
    }

    if {[set fail_name [lorf_out_name $fails]] == ""} {bell; return}
    set format    [lorf_out_get  $fails]

    #special case for a single file name
    set old_dir [pwd]
    if {[lorf_in_get $infile] == 3} {
	set inlist [lorf_in_name $infile]
	cd [file dirname $inlist]
	set inlist [file tail $inlist]
    } else {
	set inlist [lorf_get_list $infile]
	if {[file pathtype [lindex $inlist 0]] == "absolute"} {
	    cd [file dirname [lindex $inlist 0]]
	    set newlist ""
	    foreach file $inlist {
		lappend newlist [file tail $file]
	    }
	    set inlist $newlist
	} else {
	    cd [file dirname [lorf_in_name $infile]]
	}
    }

    if {![quit_displays $io "auto_assemble"]} {
	# Someone's too busy to shutdown?
	return
    }

    # Expand $inlist to include sub-files if input is a hashed archive.
    set inlist [expand_hash_archives $inlist]

    SetBusy
    set result [assemble_shotgun \
	-io $io \
	-files $inlist \
	-output_mode [radiolist_get $disp_mode] \
	-min_match [scalebox_get $match.min_match] \
	-max_pads [scalebox_get $match.max_pads] \
	-max_pmismatch [scalebox_get $match.max_mis] \
	-joins [yes_no_get $permit_joins] \
	-enter_failures [expr [radiolist_get $failure_mode]-1] \
        -tag_types $active_tags]
    catch {cd $old_dir}
    ListCreate2 $fail_name $result

    destroy $f
    update idletasks
    ClearBusy

    reset_exp_path

    #write failures to list or file
    if {"$format" == 2} {
        lorf_out_save $fail_name
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


proc Screen { f io } {
    global gap_defs

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Screen only"

    ###########################################################################
    #apply masking
    AAMasking $f $f.masking

    ###########################################################################
    AADispMode $f
    AAMatch $f
    
    ###########################################################################
    #save alignment scores in file
    keylset sa SAVEALIGN [keylget gap_defs AUTO_ASSEMBLE.SAVEALIGN]
    lorf_out $f.save_align [keylget gap_defs AUTO_ASSEMBLE.SAVEALIGN] \
	"" -bd 2 -relief groove

    ###########################################################################
    #create file of filenames entry windows
    lorf_in $f.infile [keylget gap_defs AUTO_ASSEMBLE.INFILE] \
	"{#list} {#fofn} {#selection}" -bd 2 -relief groove
    
    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "OK_Pressed2 $io $f $f.masking.yn $f.disp_mode \
	    $f.match $f.save_align $f.infile" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Assembly-Screen}" \
	    -bd 2 \
	    -relief groove
    
    ###########################################################################
    pack $f.masking -fill x
    pack $f.disp_mode -fill x
    pack $f.match -fill x
    pack $f.infile -fill x
    pack $f.save_align -fill x
    pack $f.ok_cancel -fill x
    
}

proc OK_Pressed2 { io f masking disp_mode match save_align infile} {
    global gap_defs

    SetDefaultTags AUTO_ASSEMBLE.TAGS

    #masking determines fortran variable IOPT to be 1 or 2
    set masking [yes_no_get $masking]
    if {$masking} {
	set active_tags [GetDefaultTags AUTO_ASSEMBLE.TAGS]
    } else {
	set active_tags {}
    }

    if {[set align_name [lorf_out_name $save_align]] == ""} {bell; return}
    set format     [lorf_out_get  $save_align]

    #special case for a single file name
    if {[lorf_in_get $infile] == 3} {
	set inlist [lorf_in_name $infile]
    } else {
	set inlist     [lorf_get_list $infile]
    }

    #screening doesn't do any entering of data so don't need to shut down
    #any displays
    #if {![quit_displays $io "auto_assemble"]} {
	# Someone's too busy to shutdown?
	#return
    #}

    # Expand $inlist to include sub-files if input is a hashed archive.
    set inlist [expand_hash_archives $inlist]

    SetBusy
    set result [assemble_screen -io $io \
	-files $inlist \
	-output_mode [radiolist_get $disp_mode] \
	-min_match [scalebox_get $match.min_match] \
	-max_pads [scalebox_get $match.max_pads] \
	-max_pmismatch [scalebox_get $match.max_mis] \
        -save_align 1 \
        -tag_types $active_tags]
    ListCreate2 $align_name $result

    destroy $f
    update idletasks
    ClearBusy

    reset_exp_path

    #write alignments to list or file
    if {"$format" == 2} {
	lorf_out_save $align_name
    }
}

##############################################################################
proc OneContig { f io } {
    global gap_defs

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Stack readings"

    ###########################################################################
    #create file of filenames entry windows
    lorf_in $f.infile [keylget gap_defs AUTO_ASSEMBLE.INFILE] \
	"{#list} {#fofn} {#selection}" -bd 2 -relief groove

    #save failures as file or list
    lorf_out $f.fails [keylget gap_defs AUTO_ASSEMBLE.FAILS] \
	"" -bd 2 -relief groove

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "OK_Pressed3 $io $f $f.infile $f.fails 3 one_contig" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Assembly-One}" \
	    -bd 2 \
	    -relief groove
    
    ###########################################################################
    pack $f.infile -fill x
    pack $f.fails -fill x
    pack $f.ok_cancel -fill x
    
}

proc OK_Pressed3 { io f infile fails option mode} {
    global gap_defs

    if {[set fail_name [lorf_out_name $fails]] == ""} {bell; return}
    set format    [lorf_out_get  $fails]

    #special case for a single file name
    if {[lorf_in_get $infile] == 3} {
	set inlist [lorf_in_name $infile]
    } else {
	set inlist    [lorf_get_list $infile]
    }
    if {![quit_displays $io "auto_assemble"]} {
	# Someone's too busy to shutdown?
	return
    }

    # Expand $inlist to include sub-files if input is a hashed archive.
    set inlist [expand_hash_archives $inlist]

    SetBusy
    set result [assemble_$mode -io $io\
	-files $inlist]
    ListCreate2 $fail_name $result
    
    destroy $f
    update idletasks
    ClearBusy

    reset_exp_path

    #write alignments to list or file
    if {"$format" == 2} {
        lorf_out_save $fail_name
    }
    #check that the database has contigs in it!
    if {[db_info num_contigs $io] > 0} {
        ActivateMenu_Open
        InitContigGlobals $io
        
        set cs_win [keylget gap_defs CONTIG_SEL.WIN]
        if {![winfo exists $cs_win]} {
    	ContigSelector $io
        } else {
    	raise $cs_win
    	ContigInitReg $io
        }
    }    
}

proc NewContig { f io } {
    global gap_defs

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Put all readings in separate contigs"

    ###########################################################################
    #create file of filenames entry windows
    lorf_in $f.infile [keylget gap_defs AUTO_ASSEMBLE.INFILE] \
	"{#list} {#fofn} {#selection}" -bd 2 -relief groove

    #save failures as file or list
    lorf_out $f.fails [keylget gap_defs AUTO_ASSEMBLE.FAILS] \
	"" -bd 2 -relief groove

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "OK_Pressed3 $io $f $f.infile $f.fails 4 new_contigs" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Assembly-New}" \
	    -bd 2 \
	    -relief groove
    
    ###########################################################################
    pack $f.infile -fill x
    pack $f.fails -fill x
    pack $f.ok_cancel -fill x
    
}

##############################################################################
#assembly into single stranded region
proc SSAssemble { f io } {
    global gap_defs

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Assembly into single stranded regions"

    ###########################################################################
    #apply masking
    #AAMasking $f $f.masking

    ###########################################################################
    #display mode
    AADispMode $f

    #matching scaleboxes
    AAMatch $f

    #create file of filenames entry windows
    lorf_in $f.infile [keylget gap_defs AUTO_ASSEMBLE.INFILE] \
	"{#list} {#fofn} {#selection}" -bd 2 -relief groove

    #save failures as file or list
    lorf_out $f.fails [keylget gap_defs AUTO_ASSEMBLE.FAILS] \
	"" -bd 2 -relief groove

    ###########################################################################
    #permit joins
    keylset pj PERMITJOINS [keylget gap_defs AUTO_ASSEMBLE.PERMITJOINS]
    yes_no $f.permit_joins \
	    -title [keylget pj PERMITJOINS.NAME] \
	    -orient horizontal \
	    -relief groove \
	    -bd 2 \
	    -default [keylget pj PERMITJOINS.VALUE]

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "OK_Pressed5 $io $f $f.disp_mode \
            $f.match $f.infile $f.fails $f.permit_joins" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Assembly-Single}" \
	    -bd 2 \
	    -relief groove
    
    ###########################################################################
    pack $f.disp_mode -fill x
    pack $f.match -fill x
    pack $f.infile -fill x
    pack $f.fails -fill x
    pack $f.permit_joins -fill x
    pack $f.ok_cancel -fill x
    
}

##############################################################################
proc OK_Pressed5 { io f disp_mode match infile fails \
	permit_joins} {
    global gap_defs

    if {[set fail_name [lorf_out_name $fails]] == ""} {bell; return}
    set format    [lorf_out_get  $fails]

    #special case for a single file name
    if {[lorf_in_get $infile] == 3} {
	set inlist [lorf_in_name $infile]
    } else {
	set inlist    [lorf_get_list $infile]
    }
    set option 5
    set failure 0

    if {![quit_displays $io "auto_assemble"]} {
	# Someone's too busy to shutdown?
	return
    }
 
    # Expand $inlist to include sub-files if input is a hashed archive.
    set inlist [expand_hash_archives $inlist]

    SetBusy
    set result [assemble_single_strand -io $io \
	    -files $inlist \
	    -output_mode [radiolist_get $disp_mode] \
	    -min_match [scalebox_get $match.min_match] \
	    -max_pads [scalebox_get $match.max_pads] \
	    -max_pmismatch [scalebox_get $match.max_mis] \
	    -joins [yes_no_get $permit_joins] \
	    -enter_failures $failure]
    ListCreate2 $fail_name $result

    destroy $f
    update idletasks
    ClearBusy

    reset_exp_path

    #write failures to list or file
    if {"$format" == 2} {
        lorf_out_save $fail_name
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

##############################################################################
#normal shotgun assembly but ignore previously assembled data
proc IgnorePrev { f io } {
    global gap_defs

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Assemble independently"

    ###########################################################################
    #apply masking
    AAMasking $f $f.masking

    ###########################################################################
    #display mode
    AADispMode $f

    #matching scaleboxes
    AAMatch $f

    #create file of filenames entry windows
    lorf_in $f.infile [keylget gap_defs AUTO_ASSEMBLE.INFILE] \
	"{#list} {#fofn} {#selection}" -bd 2 -relief groove

    #save failures as file or list
    lorf_out $f.fails [keylget gap_defs AUTO_ASSEMBLE.FAILS] \
	"" -bd 2 -relief groove

    ###########################################################################
    #select alignment failure mode
    keylset fm FAILURE [keylget gap_defs AUTO_ASSEMBLE.FAILURE]
    set b1 [keylget fm FAILURE.BUTTON.1]
    set b2 [keylget fm FAILURE.BUTTON.2]
    radiolist $f.failure_mode \
	    -title [keylget fm FAILURE.NAME] \
	    -relief groove \
	    -bd 2 \
	    -default [keylget fm FAILURE.VALUE] \
	    -buttons [format { { %s } { %s } } [list $b1] [list $b2] ]
    
    ###########################################################################
    #permit joins
    keylset pj PERMITJOINS [keylget gap_defs AUTO_ASSEMBLE.PERMITJOINS]
    yes_no $f.permit_joins \
	    -title [keylget pj PERMITJOINS.NAME] \
	    -orient horizontal \
	    -relief groove \
	    -bd 2 \
	    -default [keylget pj PERMITJOINS.VALUE]

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "OK_Pressed6 $io $f $f.masking.yn $f.disp_mode \
            $f.match $f.infile $f.fails $f.permit_joins $f.failure_mode" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Assembly-Ind}" \
	    -bd 2 \
	    -relief groove
    
    ###########################################################################
    pack $f.masking -fill x
    pack $f.disp_mode -fill x
    pack $f.match -fill x
    pack $f.infile -fill x
    pack $f.fails -fill x
    pack $f.permit_joins -fill x
    pack $f.failure_mode -fill x
    pack $f.ok_cancel -fill x
    
}

##############################################################################
proc OK_Pressed6 { io f masking disp_mode match infile fails \
	permit_joins failure_mode} {
    global gap_defs

    #masking determines fortran variable IOPT to be 1 or 2
    set masking [yes_no_get $masking]
    if {$masking} {
	set active_tags [GetDefaultTags AUTO_ASSEMBLE.TAGS]
    } else {
	set active_tags {}
    }

    if {[set fail_name [lorf_out_name $fails]] == ""} {bell; return}
    set format    [lorf_out_get  $fails]
    #special case for a single file name
    if {[lorf_in_get $infile] == 3} {
	set inlist [lorf_in_name $infile]
    } else {
	set inlist    [lorf_get_list $infile]
    }
    if {![quit_displays $io "auto_assemble"]} {
	# Someone's too busy to shutdown?
	return
    }

    # Expand $inlist to include sub-files if input is a hashed archive.
    set inlist [expand_hash_archives $inlist]

    SetBusy
    set result [assemble_independent -io $io \
	-files $inlist \
	-output_mode [radiolist_get $disp_mode] \
	-min_match [scalebox_get $match.min_match] \
	-max_pads [scalebox_get $match.max_pads] \
	-max_pmismatch [scalebox_get $match.max_mis] \
	-joins [yes_no_get $permit_joins] \
	-enter_failures [expr [radiolist_get $failure_mode]-1] \
        -tag_types $active_tags]
    ListCreate2 $fail_name $result

    destroy $f
    update idletasks
    ClearBusy

    reset_exp_path

    #write failures to list or file
    if {"$format" == 2} {
        lorf_out_save $fail_name
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

#############################################################################
# Checks the files named in $inlist to see if any appear to be an
# hashed archive.
# If so it expands them up to contain the contents
proc expand_hash_archives {inlist} {
    global env old_exp_path
    global gap_defs

    if {![keylget gap_defs AUTO_ASSEMBLE.USE_EXP_ARCHIVES]} {
	return $inlist
    }

    set outlist {}
    if {[info exists env(EXP_PATH)]} {
	set old_exp_path $env(EXP_PATH)
    }

    foreach f $inlist {
	set fd [open [list "|hash_list" $f]]
	set hl [read $fd]
	catch {close $fd}
	set done 0
	foreach newfile $hl {
	    if {$newfile != {}} {
		lappend outlist $newfile
		set done 1
	    }
	}
	if {!$done} {
	    lappend outlist $f
	} else {
	    if {[info exists env(EXP_PATH)]} {
		set env(EXP_PATH) $env(EXP_PATH):HASH=$f
	    } else {
		set env(EXP_PATH) HASH=$f
	    }
	}
    }

    return $outlist
}

proc reset_exp_path {} {
    global env old_exp_path
    if {[info exists old_exp_path]} {
	set env(EXP_PATH) $old_exp_path
	unset old_exp_path
    }
}