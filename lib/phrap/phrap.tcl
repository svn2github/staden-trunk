# -----------------------------------------------------------------------------
# Utility functions

proc phrap_out_handler {fd} {
    vmessage -nonewline [read $fd]
    if {[eof $fd]} {
	global phrap_done
	set phrap_done 1
    }
}

proc PhrapAssemble_phrap {dir params} {
    global tcl_platform
    vfuncheader "Phrap assembly - phrap"
    SetBusy
    file mkdir $dir
    if {$tcl_platform(platform) == "unix"} {
        set stdout_fd [open "|sh -c {gcphrap $params 2>&1 > $dir/stdout}"]
	fconfigure $stdout_fd -blocking 0
	fileevent $stdout_fd readable "phrap_out_handler $stdout_fd"
	vwait phrap_done
        close $stdout_fd
    } else {
        vmessage "Phrap running - please wait"
	update idletasks
        catch {eval exec gcphrap $params > [list $dir/stdout]} err
	vmessage $err
    }
    ClearBusy
}

proc PhrapAssemble_enter {io dest_dir fofn} {
    vfuncheader "Phrap assembly - entering"

    set pwd [pwd]
    cd $dest_dir
    if {[catch {set fd [open $fofn.ph]}] == 1} {
	verror ERR_WARN "Could not open $fofn.ph"
	cd $pwd
	return
    }
    set lin [read $fd]
    close $fd

    SetBusy
    catch {assemble_direct \
	-io $io \
	-files $lin \
	-max_pmismatch 100 \
	-output_mode 0 \
	-align 0}
    ClearBusy

    cd $pwd
}

# -----------------------------------------------------------------------------
# Main assembly driver
#
# Strategy is:
#    Phrap assembly
#    Directed assembly

proc PhrapAssemble2 {io t infile params dest_dir doqclip minqual dodclip} {
    global gap_defs

    if {[set do_qclipping [yes_no_get $doqclip]] == 1} {
	set quality [scalebox_get $minqual]
    }
    set do_dclipping [yes_no_get $dodclip]

    # Get and/or create file of filenames
    if {[set lin [lorf_in_name $infile]] == ""} {
	bell
	return
    }

    if {[lorf_in_get $infile] == 1} {
	set fd [open "@$lin" w]
	regsub -all "\[ \t\]+" [ListGet $lin] "\n" f
	puts $fd $f
	close $fd
	set fofn "@$lin"
    } else {
	set fofn $lin
    }

    if {![quit_displays $io "phrap_assemble"]} {
        # Someone's too busy to shutdown?
        return
    }

    # Get phrap parameters
    set params [entrybox_get $params]
    set dir [entrybox_get $dest_dir]

    set db [io_read_database $io]
    set old_num_contigs [keylget db num_contigs]

    # Run phrap
    PhrapAssemble_phrap $dir "$params -exp $dir $fofn"

    # Enter phrap assembly
    PhrapAssemble_enter $io $dir $fofn

    # Clip on quality or differences
    if {$do_qclipping == 1 || $do_dclipping == 1} {
        set db [io_read_database $io]
        set num_contigs [keylget db num_contigs]
	set list {}
	for {set i $old_num_contigs} {$i < $num_contigs} {incr i} {
	    lappend list =[expr $i+1]
	}

        if {$do_qclipping} {
	    quality_clip -io $io -contigs $list -quality $quality
        }
        if {$do_dclipping} {
	    difference_clip -io $io -contigs $list
        }
    }

    # Tidy up
    destroy $t

    if {[lorf_in_get $infile] == 1} {
	file delete "@$lin"
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

# -----------------------------------------------------------------------------
# Reassembly driver
#
# Strategy is:
#    Backup database
#    Extract readings
#    Disassemble readings
#    Phrap assembly
#    Directed assembly

proc PhrapReassemble2 {io t infile params dest_dir doqclip minqual dodclip} {
    global gap_defs

    if {[set do_qclipping [yes_no_get $doqclip]] == 1} {
	set quality [scalebox_get $minqual]
    }
    set do_dclipping [yes_no_get $dodclip]

    # Get and/or create file of filenames
    if {[set lin [lorf_in_name $infile]] == ""} {
	bell
	return
    }

    if {[lorf_in_get $infile] == 1} {
	set fd [open "@$lin" w]
	regsub -all {[ \t]+} [ListGet $lin] "\n" f
	puts $fd $f
	close $fd
	set fofn "@$lin"
    } else {
	set fofn $lin
    }

    set list [lorf_get_list $infile]

    if {![quit_displays $io "phrap_assemble"]} {
        # Someone's too busy to shutdown?
        return
    }

    # Get phrap parameters
    set params [entrybox_get $params]
    set dir [entrybox_get $dest_dir]

    SetBusy

    # Create a database backup
    if {[catch {copy_db \
	-io $io \
	-version ~ \
	-collect 0}]} {
	verror ERR_WARN "Could not backup database to version ~"
	ClearBusy
	return
    }

    # Extract readings
    if {[catch {extract_readings \
	-io $io \
	-readings $list \
	-directory $dir \
	-format 2}]} {
	verror ERR_WARN "Failed when extracting readings"
	ClearBusy
	return
    }

    # Disassembly
    if {[catch {disassemble_readings \
	-io $io \
	-readings $list \
	-move 0}]} {
	verror ERR_WARN "Failed when disassembling readings."
	verror ERR_WARN "Remember that there is a backup - version \"~\"."
	ClearBusy
	return
    }
    ClearBusy

    set db [io_read_database $io]
    set old_num_contigs [keylget db num_contigs]

    # Run phrap
    set pwd [pwd]
    cd $dir
    if {[catch {PhrapAssemble_phrap $dir "$params -exp . fofn"}]} {
	verror ERR_WARN "Phrap assembly failed."
	verror ERR_WARN "Remember that there is a backup - version \"~\"."
	cd $pwd
	ClearBusy
	return
    }
    cd $pwd

    # Enter phrap assembly
    if {[catch {PhrapAssemble_enter $io $dir fofn}]} {
	verror ERR_WARN "Directed assembly failed."
	verror ERR_WARN "Remember that there is a backup - version \"~\"."
	ClearBusy
	return
    }

    # Clip on quality or differences
    if {$do_qclipping == 1 || $do_dclipping == 1} {
        set db [io_read_database $io]
        set num_contigs [keylget db num_contigs]
	set list {}
	for {set i $old_num_contigs} {$i < $num_contigs} {incr i} {
	    lappend list =[expr $i+1]
	}

        if {$do_qclipping} {
	    quality_clip -io $io -contigs $list -quality $quality
        }
        if {$do_dclipping} {
	    difference_clip -io $io -contigs $list
        }
    }

    # Tidy up
    destroy $t

    if {[lorf_in_get $infile] == 1} {
	file delete "@$lin"
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

# -----------------------------------------------------------------------------
# Phrap Assemble and Reassemble GUI

proc PhrapAssemble {io} {
    global phrap_defs gap_defs

    set t [keylget phrap_defs PHRAP.ASS.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Phrap assembly"

    # Input file(s)
    lorf_in $t.infile [keylget phrap_defs PHRAP.INFILE] "" \
	-bd 2 -relief groove

    #destination directory
    frame $t.dest_dir -bd 2 -relief groove
    set l [keylget phrap_defs PHRAP.DEST_DIR]
    entrybox $t.dest_dir.e -bd 2 -relief sunken\
	    -title [keylget l NAME ] \
	    -default [keylget l VALUE] \
	    -width 15 \
	    -type CheckString

    # Other phrap parameters
    frame $t.params -bd 2 -relief groove
    set l [keylget phrap_defs PHRAP.PARAM]
    entrybox $t.params.e -bd 2 -relief sunken \
	-title [keylget l NAME] \
	-default [keylget l VALUE] \
	-width 15

    # Quality clipping
    frame $t.qclip -bd 2 -relief groove
    set l [keylget gap_defs QUALITY_CLIP]
    scalebox $t.qclip.minqual \
	-orient horizontal -width 5 \
	-title   "[keylget l MINQUAL.NAME]" \
	-from    "[keylget l MINQUAL.MIN]" \
	-to      "[keylget l MINQUAL.MAX]" \
	-default "[keylget l MINQUAL.VALUE]"
    set l [keylget phrap_defs PHRAP.QUALITY_CLIP]
    yes_no $t.qclip.doclip \
	-title "[keylget l NAME]" \
	-default "[keylget l VALUE]" \
	-orient horizontal \
	-ycommand "scalebox_configure $t.qclip.minqual -state normal" \
	-ncommand "scalebox_configure $t.qclip.minqual -state disabled"

    # Difference clipping
    frame $t.dclip -bd 2 -relief groove
    set l [keylget phrap_defs PHRAP.DIFFERENCE_CLIP]
    yes_no $t.dclip.doclip \
	-title "[keylget l NAME]" \
	-default "[keylget l VALUE]" \
	-orient horizontal

    okcancelhelp $t.ok -bd 2 -relief groove \
	-ok_command "PhrapAssemble2 $io $t $t.infile $t.params.e \
		$t.dest_dir.e $t.qclip.doclip $t.qclip.minqual \
		$t.dclip.doclip" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Assembly-Phrap Assemble}"

    pack $t.infile $t.dest_dir $t.dest_dir.e $t.params $t.params.e \
	$t.qclip $t.qclip.doclip $t.qclip.minqual $t.dclip $t.dclip.doclip \
	$t.ok -side top -fill both
}

proc PhrapReassemble {io} {
    global phrap_defs gap_defs

    set t [keylget phrap_defs PHRAP.REASS.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Phrap reassembly"

    # Input file(s)
    lorf_in $t.infile [keylget phrap_defs PHRAP.INFILE] "" \
	-bd 2 -relief groove

    #destination directory
    frame $t.dest_dir -bd 2 -relief groove
    set l [keylget phrap_defs PHRAP.DEST_DIR]
    entrybox $t.dest_dir.e -bd 2 -relief sunken\
	    -title [keylget l NAME ] \
	    -default [keylget l VALUE] \
	    -width 15 \
	    -type CheckString

    # Other phrap parameters
    frame $t.params -bd 2 -relief groove
    set l [keylget phrap_defs PHRAP.PARAM]
    entrybox $t.params.e -bd 2 -relief sunken \
	-title [keylget l NAME] \
	-default [keylget l VALUE] \
	-width 15

    # Quality clipping
    frame $t.qclip -bd 2 -relief groove
    set l [keylget gap_defs QUALITY_CLIP]
    scalebox $t.qclip.minqual \
	-orient horizontal -width 5 \
	-title   "[keylget l MINQUAL.NAME]" \
	-from    "[keylget l MINQUAL.MIN]" \
	-to      "[keylget l MINQUAL.MAX]" \
	-default "[keylget l MINQUAL.VALUE]"
    set l [keylget phrap_defs PHRAP.QUALITY_CLIP]
    yes_no $t.qclip.doclip \
	-title "[keylget l NAME]" \
	-default "[keylget l VALUE]" \
	-orient horizontal \
	-ycommand "scalebox_configure $t.qclip.minqual -state normal" \
	-ncommand "scalebox_configure $t.qclip.minqual -state disabled"

    # Difference clipping
    frame $t.dclip -bd 2 -relief groove
    set l [keylget phrap_defs PHRAP.DIFFERENCE_CLIP]
    yes_no $t.dclip.doclip \
	-title "[keylget l NAME]" \
	-default "[keylget l VALUE]" \
	-orient horizontal

    okcancelhelp $t.ok -bd 2 -relief groove \
	-ok_command "PhrapReassemble2 $io $t $t.infile $t.params.e \
		$t.dest_dir.e $t.qclip.doclip $t.qclip.minqual \
		$t.dclip.doclip" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Assembly-Phrap Reassemble}"

    pack $t.infile $t.dest_dir $t.dest_dir.e $t.params $t.params.e \
	$t.qclip $t.qclip.doclip $t.qclip.minqual $t.dclip $t.dclip.doclip \
	$t.ok -side top -fill both
}
