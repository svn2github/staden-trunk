#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
##############################################################################
proc DatabaseInfo {io} {
    global maxseq

    set db [io_read_database $io]
    set nc [keylget db num_contigs]
    set nr [keylget db num_readings]
    set tcl [db_info t_contig_length $io]
    set trl [db_info t_read_length $io]
    set maxgel [keylget db max_gel_len]

    set i ""
    append i [format "Database size       %7d       Max reading length %7d\n" \
		 [keylget db actual_db_size] $maxgel]
    append i [format "No. Readings        %7d       No. Contigs        %7d\n" \
	         $nr $nc]
    append i [format "No. Annotations     %7d       No. Templates      %7d\n" \
		 [keylget db Nannotations] [keylget db Ntemplates]]
    append i [format "No. Clones          %7d       No. Vectors        %7d\n" \
		 [keylget db Nclones]      [keylget db Nvectors]]
    append i [format "Total contig length %7d       Average length    %8.1f\n"\
		 $tcl [expr double($tcl)/$nc]]
    append i [format "Total characters in readings                    %12d\n" \
		 $trl]
    append i [format "Average reading characters per consensus character \
	    %8.2f\n" \
	        [expr double($trl)/$tcl]]
    append i [format "Average used length of reading                     \
	    %8.2f\n" \
		[expr double($trl) / $nr]]

    set maxgel 1024; # minor space saving
    if {[expr {($tcl + (2 * $maxgel + 20)*$nc)*1.1}] > $maxseq} {
	set maxseq [expr {round(($tcl + (2 * $maxgel + 20)*$nc)*1.1)}]
	verror ERR_WARN gap "increasing maxseq parameter to $maxseq"
    }
    
    append i "Current maximum consensus length is $maxseq\n"

    return $i
}

proc PrintDatabaseInfo {io} {
    vfuncheader "Database information [db_info db_name $io]"
    vmessage [DatabaseInfo $io]
}

proc copy_reads_open_database {a} {
    global read_only 
    set create 0

    if {[regexp {(.*)\.(.)(\.aux)?$} $a tmp db_name version_num] == 0} {
	verror ERR_FATAL "Gap4" "ERROR: Invalid database name '$a'"
	return -1
    }

    cd [file dirname $db_name]
    set db_name [file tail $db_name]

    #open database 
    set io [open_db -name $db_name -version $version_num -access rw \
	    -create $create] 

    if {$read_only == 1} {
	tk_messageBox -icon warning -type ok \
		-title "Database open" \
		-message "Database $db_name.$version_num is busy. Database must be writable"
	return -1
    }

    if {"$io" == ""} {
#	6/1/99 johnt - use verror so message is shown with Windows NT
	verror ERR_FATAL "copy reads" "ERROR: Couldn't open the database '$db_name.$version_num'."
	return -1
    }

    return $io
}

proc InvokeDBBrowser {fn } {

    set file_type {
	{"database"		    *.aux		    }
    }

    set file [tk_getOpenFile -filetypes $file_type -parent $fn]

    raise [winfo toplevel $fn]

    #if user has pressed cancel don't want to remove name from entrybox
    if {$file != ""} {
	entrybox_delete $fn 0 end
	entrybox_insert $fn 0 $file
	[entrybox_path $fn] xview [string last / $file]
    }
}

proc getDBname {t title value} {
    global $t.text

    set $t.text $value

    entrybox $t.entry \
	-title $title \
	-width 15 \
	-textvariable $t.text

    button $t.browse -text Browse -command "InvokeDBBrowser $t.entry"
    pack $t.entry -side left -expand 1 -fill x
    pack $t.browse -side right
}

proc getDBname_focus {t} {
    entrybox_focus $t.entry
}

proc getDBname_in_name {t} {
    return [expandpath [entrybox_get $t.entry]]
}

proc CheckDirectory { path } {
    set dir [expandpath [$path.entry get]]
    if {$dir == ""} {
	return 1
    }
    set response [file isdirectory $dir]
    return $response
}

proc CopyReads {io_to} {
    global gap_defs copy_reads_defs
    global argc argv

    set f .copy_reads
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Copy reads"

    ###########################################################################
    #read in FROM and TO databases
    frame $f.db -bd 2 -relief groove
    frame $f.db.from 
    frame $f.db.to
    
    if {$io_to != -1} {
	#interface from gap4
	getDBname $f.db.from "Open source database" "" 
    } elseif {$argc == 2} {
	#read from command line if available
	getDBname $f.db.from "Open source database" [lindex $argv 0] 
	getDBname $f.db.to "Open destination database" [lindex $argv 1]	
    } else {
	getDBname $f.db.from "Open source database" "" 
	getDBname $f.db.to "Open destination database" ""
    }
    
    pack $f.db.from -expand 1 -fill x
    if {$io_to == -1} {
	pack $f.db.to -expand 1 -fill x
    }
    ###########################################################################
    #trace files
    #FIXME - need directory browser
    frame $f.traces -bd 2 -relief groove

    set tf [keylget copy_reads_defs COPY_READS.TRACES.FROM]
    entrybox $f.traces.from_dir \
	-title "[keylget tf NAME]"\
	-default [keylget tf VALUE]\
	-type CheckDirectory

    set tf [keylget copy_reads_defs COPY_READS.TRACES.RADIO]
    set b1 [keylget tf BUTTON.1]
    set b2 [keylget tf BUTTON.2]
    radiolist $f.traces.rl \
	    -title [keylget tf NAME] \
	    -default [keylget tf VALUE]\
	    -orient horizontal \
	    -buttons [format { \
	    {%s -command {entrybox_configure %s -state normal}} \
	    {%s -command {entrybox_configure %s -state disabled}} } \
	    [list $b1] [list $f.traces.from_dir] [list $b2] [list $f.traces.from_dir]]

    pack $f.traces.rl -anchor w -fill x
    pack $f.traces.from_dir -anchor w -fill x

    ###########################################################################
    #select contigs FROM and TO
    frame $f.c_from -bd 2 -relief groove
    lorf_in $f.c_from.list [keylget copy_reads_defs COPY_READS.INFILE1] "#fofn #" 

    #minimum contig length of FROM db
    set mc [keylget copy_reads_defs COPY_READS.MINCL]
    entrybox $f.c_from.min_cl \
	-title "[keylget mc NAME]"\
	-default [keylget mc VALUE]\
	-type "CheckIntMin [keylget mc MIN]"

    #minimum average quality
    set mq [keylget copy_reads_defs COPY_READS.MIN_QUAL]
    entrybox $f.c_from.min_qual \
	-title "[keylget mq NAME]"\
	-default [keylget mq VALUE]\
	-type "CheckFloatRange [keylget mq MIN] [keylget mq MAX]"
    
    pack $f.c_from.list -fill x
    pack $f.c_from.min_cl -fill x
    pack $f.c_from.min_qual -fill x
    
    lorf_in $f.c_to [keylget copy_reads_defs COPY_READS.INFILE2] "#fofn #" -bd 2 -relief groove

    ###########################################################################
    #select mode
    SetDefaultTags FIJ.TAGS

    set sm [keylget gap_defs FIJ.SELMODE]

    set b1 [keylget sm BUTTON.1]
    set b2 [keylget sm BUTTON.3]
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
	    SetDefaultTags %s } } } \
	    [list $b1] [list $f.sel_mode.but] \
	    [list $b2] [list $f.sel_mode.but] FIJ.TAGS ]
    pack $f.sel_mode.rl -side left
    pack $f.sel_mode.but -side right
    
    ###########################################################################
    #consensus alignment parameters
    frame $f.c_align -relief groove -bd 2
    label $f.c_align.label -text [keylget copy_reads_defs COPY_READS.ALIGN_CONS]

    ###########################################################################
    #select word length

    set st [keylget gap_defs FIJ.WORDLENGTH]
    set b1 [keylget st BUTTON.1]
    set b2 [keylget st BUTTON.2]

    radiolist $f.c_align.word_length \
	    -title [keylget st NAME]\
	    -bd 0 \
	    -relief groove \
	    -orient horizontal \
	    -default [keylget st VALUE] \
	    -buttons [format { {%s} {%s} } \
	    [list $b1] [list $b2] ]

    ###########################################################################
    set mm [keylget gap_defs FIJ.MINOVERLAP]
    scalebox $f.c_align.min_overlap \
	    -title [keylget mm NAME]\
	    -orient horizontal \
	    -to [keylget mm MAX] \
	    -from [keylget mm MIN]\
	    -default [keylget mm VALUE] \
	    -width 5 \
	    -type CheckInt
    
    set mm [keylget gap_defs FIJ.MINMATCH]
    scalebox $f.c_align.min_match \
	    -title [keylget mm NAME]\
	    -orient horizontal \
	    -to [keylget mm MAX]\
	    -from [keylget mm MIN] \
	    -default [keylget mm VALUE] \
	    -width 5 \
	    -type CheckInt

    set mm [keylget gap_defs FIJ.MAXMIS]
    scalebox $f.c_align.max_mis \
	    -title [keylget mm NAME]\
	    -orient horizontal \
	    -to [keylget mm MAX]\
	    -from [keylget mm MIN] \
	    -resolution [keylget mm RES] \
	    -default [keylget mm VALUE] \
	    -width 5 \
	    -type CheckFloat
    
    set mm [keylget gap_defs FIJ.USEBAND] 
    yes_no $f.c_align.use_band \
    	    -title [keylget mm NAME] \
	    -orient horizontal \
	    -bd 0 \
	    -default [keylget mm VALUE]

    set mm [keylget copy_reads_defs COPY_READS.DISPLAY_CONS]
    global $f.Display_cons_alignments
    checkbutton $f.c_align.display \
	-text [keylget mm NAME] \
    	-variable $f.Display_cons_alignments

    set $f.Display_cons_alignments [keylget mm VALUE]

    ###########################################################################

    pack $f.c_align.label
    pack $f.c_align.word_length -fill x
    pack $f.c_align.min_overlap -fill x
    pack $f.c_align.max_mis -fill x
    pack $f.c_align.min_match -fill x
    pack $f.c_align.use_band -fill x
    pack $f.c_align.display -fill x -side left 

    ###########################################################################
    #Reading alignment parameters

    frame $f.r_align -bd 2 -relief groove
    label $f.r_align.label -text [keylget copy_reads_defs COPY_READS.ALIGN_READS]

    #maximum alignment mismatch
    set mm [keylget gap_defs DIRECT_ASSEMBLY.MAXMIS]
    entrybox $f.r_align.align_mism \
	-title "[keylget mm NAME] ([keylget mm MIN] to [keylget mm MAX])"\
	-default [keylget mm VALUE]\
	-type "CheckFloatRange [keylget mm MIN] [keylget mm MAX]"

    set mm [keylget copy_reads_defs COPY_READS.DISPLAY_SEQ]
    global $f.Display_seq_alignments
    checkbutton $f.r_align.display \
	-text [keylget mm NAME] \
    	-variable $f.Display_seq_alignments

    set $f.Display_seq_alignments [keylget mm VALUE]

    ###########################################################################
    
    pack $f.r_align.label
    pack $f.r_align.align_mism -fill x
    pack $f.r_align.display -fill x -side left

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "CopyReads2 $f $f.db.from $f.db.to $io_to $f.traces.rl $f.traces.from_dir $f.c_from.list $f.c_to $f.sel_mode.rl $f.c_align.word_length $f.c_align.min_overlap $f.c_align.min_match $f.c_align.use_band $f.c_align.max_mis $f.r_align.align_mism $f.c_from.min_cl $f.c_from.min_qual \[set $f.Display_cons_alignments\] \[set $f.Display_seq_alignments\]"\
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Copy Reads-Dialogue}" \
	    -bd 2 \
	    -relief groove
    ###########################################################################

    pack $f.db -expand 1 -fill x
    pack $f.traces -fill x
    pack $f.c_from -fill x
    pack $f.c_to -fill x
    pack $f.sel_mode -fill x
    pack $f.c_align -fill x
    pack $f.r_align -fill x
    pack $f.ok_cancel -fill x

}

proc CopyReads2 {f db_from db_to io_to traces dir c_from c_to sel_mode word_length min_overlap min_match use_band max_mis align_mis min_cl min_qual display_cons display_seq} {
    global env

    # flag for whether standalone or via gap4 interface
    set standalone 0

    ###########################################################################
    set source_mode [radiolist_get $traces]
    set rawdata ""

    #user provides source dir
    if {$source_mode == 1} {
	set rawdata [expandpath [entrybox_get $dir]]
    }

    ###########################################################################
    #select mode 
    set active_tags {}

    set masking [radiolist_get $sel_mode]

    #if masking mode is 2 (mask active tags)
    if {($masking == 2)} {
        set active_tags [GetDefaultTags FIJ.TAGS]
    }

    ###########################################################################
    #match info
    set word_length [lindex {? 8 4} [radiolist_get $word_length]]
    set min_overlap [scalebox_get $min_overlap]
    set min_match [scalebox_get $min_match]
    set band_size [yes_no_get $use_band]
    set max_mis [scalebox_get $max_mis]
    
    set align_mis [entrybox_get $align_mis]
    set min_cl [entrybox_get $min_cl]
    set min_qual [entrybox_get $min_qual]

    #open the databases last so that any errors trapped before this doesn't
    #mean that the databases are already open and busy

    ###########################################################################
    #open databases
    set file_from [getDBname_in_name $db_from]
    puts file_from=$file_from
    if {$io_to == -1} {
	set file_to [getDBname_in_name $db_to]
    }

    #must have db name
    if {$file_from == ""} {
	getDBname_focus $db_from; return
    }
    if {$io_to == -1} {
	if {$file_to == ""} {
	    getDBname_focus $db_to; return
	}
    }
    
    #open db and check is writable
    set io_from [copy_reads_open_database $file_from]
    if {$io_from == -1} {
	getDBname_focus $db_from; return
    }

    if {$io_to == -1} {
	set io_to [copy_reads_open_database $file_to]
	set standalone 1
    }

    if {$io_to == -1} {
	#must close $io_from otherwise I will always get an error
	close_db -io $io_from
	getDBname_focus $db_to; return
    }

    PrintDatabaseInfo $io_from
    PrintDatabaseInfo $io_to

    ###########################################################################
    #read in contig lists
    if {[lorf_in_get $c_from] == 2} {
	set list_from [CreateAllContigList $io_from]
    } else {
	set list_from [lorf_get_list $c_from]
    }

    if {[lorf_in_get $c_to] == 2} {
	set list_to [CreateAllContigList $io_to]
    } else {
	set list_to [lorf_get_list $c_to]
    }

    if {!$standalone} {
	if {![quit_displays $io_to "copy_reads"]} {
	    # Someone's too busy to shutdown?
	    close_db -io $io_from
	    return
	}
    }
 
   # Destroy dialog before showing plot to avoid window activation problems
    destroy $f

    #puts "-io_from $io_from -io_to $io_to \
	    -contigs_from $list_from -contigs_to $list_to\
	    -mask [lindex {"" none mask} $masking] \
	    -min_overlap $min_overlap \
	    -max_pmismatch $max_mis \
	    -word_length $word_length \
	    -min_match $min_match \
	    -band $band_size \
	    -tag_types $active_tags \
	    -align_max_mism $align_mis\
	    -min_contig_len $min_cl\
	    -min_average_qual $min_qual\
	    -display_cons $display_cons\
	    -display_seq $display_seq"

    SetBusy

    set list [copy_reads -io_from $io_from -io_to $io_to \
	    -contigs_from $list_from -contigs_to $list_to\
	    -mask [lindex {"" none mask} $masking] \
	    -min_overlap $min_overlap \
	    -max_pmismatch $max_mis \
	    -word_length $word_length \
	    -min_match $min_match \
	    -band $band_size \
	    -tag_types $active_tags \
	    -align_max_mism $align_mis\
	    -min_contig_len $min_cl\
	    -min_average_qual $min_qual\
	    -display_cons $display_cons\
	    -display_seq $display_seq]

    ClearBusy

    ###########################################################################
    #only add rawdata note if readings were added
    if {[llength $list] > 0} { 
	#read rawdata note from source database
	if {$source_mode == 2} {
	    set nn 0
	    if {[set nn [FindRawDataNote $io_from]] != 0} {
		set n [io_read_note $io_from $nn]
		if {[keylget n annotation] != 0} {
		    set rawdata [io_read_text $io_from [keylget n annotation]]
		}
	    }
	}
	
	if {$rawdata != ""} {
	    #leave - add location of FROM traces to TO database
	    AppendRawDataNote $io_to $rawdata
	}
    }
    close_db -io $io_from
    if {$standalone} {
	close_db -io $io_to
    }

    if {!$standalone} {
	global gap_defs
	#check that the database has contigs in it!
	if {[db_info num_contigs $io_to] > 0} {
	    ActivateMenu_Open
	    InitContigGlobals $io_to
	    
	    #draw the contig selector if it does not already exist
	    set cs_win [keylget gap_defs CONTIG_SEL.WIN]
	    if {![winfo exists $cs_win]} {
		ContigSelector $io_to
	    } else {
		ContigInitReg $io_to
		raise $cs_win
	    }
	}
    }
}

proc do_copy_reads {file_from file_to list_from list_to source_trace_dir arg} {

    #open db and check is writable
    set io_from [copy_reads_open_database $file_from]
    if {$io_from == -1} {
	return
    }

    set io_to [copy_reads_open_database $file_to]
    if {$io_to == -1} {
	#must close $io_from otherwise I will always get an error
	close_db -io $io_from
	return
    }

    PrintDatabaseInfo $io_from
    PrintDatabaseInfo $io_to

    if {$list_from == ""} {
	set list_from [CreateAllContigList $io_from]
    } else {
	set list_from [ListLoad $list_from from]
    }

    if {$list_to == ""} {
	set list_to [CreateAllContigList $io_to]
    } else {
	set list_to [ListLoad $list_to to]
    }

    if {$source_trace_dir == ""} {
	set source_mode 2
    } else {
	set source_mode 1
    }
    set rawdata ""
    
    #user provides source dir
    if {$source_mode == 1} {
	set rawdata [expandpath $source_trace_dir]
    }

    set list [eval copy_reads -io_from $io_from -io_to $io_to \
	    -contigs_from [list $list_from] -contigs_to [list $list_to] $arg]

    ###########################################################################
    #only add rawdata note if readings were added
    if {[llength $list] > 0} { 
	#read rawdata note from source database
	if {$source_mode == 2} {
	    set nn 0
	    if {[set nn [FindRawDataNote $io_from]] != 0} {
		set n [io_read_note $io_from $nn]
		if {[keylget n annotation] != 0} {
		    set rawdata [io_read_text $io_from [keylget n annotation]]
		}
	    }
	}
	
	if {$rawdata != ""} {
	    #leave - add location of FROM traces to TO database
	    AppendRawDataNote $io_to $rawdata
	}
    }
    close_db -io $io_from
    close_db -io $io_to
}
