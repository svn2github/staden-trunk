#-----------------------------------------------------------------------------
# Import sequences

proc ImportSequences { io job } {
    global gap5_defs

    # Make ourselves a window
    set f .import_sequences
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Import Sequences"

    # Input file
    getFname $f.infile "Filename" load {}
   
    # Packing list.  Note $f.format is created at the end
    set to_pack [list $f.infile $f.format]

    if { $job == "merge" } {
	# Merge options
	checkbutton $f.new_contigs_btn \
	    -text "Always make new contigs" \
	    -variable $f.new_contigs \
	    -command "ImportSequences_configure $f $job" \
	    -anchor w
	global $f.new_contigs
	set $f.new_contigs 1

	checkbutton $f.repad_btn \
	    -text "Re-pad sequences when merging" \
	    -variable $f.repad \
	    -anchor w
	global $f.repad
	set $f.repad 1

	lappend to_pack $f.new_contigs_btn $f.repad_btn
    } else {
	global $f.browsed_name
	set $f.browsed_name ""

	frame $f.database
	xentry $f.database.entry -label "Database name" -width 15
	button $f.database.browse -text "Browse" \
	    -command "NewFileBrowser_invoke $f.database.entry $f.browsed_path"

	pack $f.database.entry -side left -expand 1 -fill x
	pack $f.database.browse -side right

	lappend to_pack $f.database
    }

    checkbutton $f.refpos_btn \
	-text "Store reference position data (SAM/BAM only)" \
	-variable $f.refpos \
	-anchor w
    global $f.refpos
    set $f.refpos 1

    checkbutton $f.rmdup_btn \
	-text "Ignore reads marked as duplicates (SAM/BAM only)" \
	-variable $f.rmdup \
	-anchor w
    global $f.rmdup
    set $f.rmdup 1
    
    lappend to_pack $f.refpos_btn $f.rmdup_btn

    okcancelhelp $f.ok \
	-ok_command "ImportSequences2 \$io $f $job" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 ImportSequences" \
	-bd 2 \
	-relief groove

    lappend to_pack $f.ok

    # Input format - has to be last as radiolist calls ImportSequences_configure
    radiolist $f.format \
	-title "Input file type" \
	-default 5 \
	-orient horizontal \
	-buttons [list \
		      [list ace -command "ImportSequences_configure $f $job"] \
		      [list caf -command "ImportSequences_configure $f $job"] \
		      [list sam -command "ImportSequences_configure $f $job"] \
		      [list bam -command "ImportSequences_configure $f $job"] \
		      [list auto -command "ImportSequences_configure $f $job"]]

    eval pack $to_pack -side top -fill both
}

proc ImportSequences_configure { f job } {
    # Enable / disable SAM/BAM specific options.  Eanbled for auto

    # Radio buttons: auto ace caf sam bam auto
    set sambam [lindex "1 0 0 1 1 1" [radiolist_get $f.format]]

    if { $sambam } {
	$f.refpos_btn configure -state normal
	$f.rmdup_btn configure -state normal
    } else {
	$f.refpos_btn configure -state disabled
	$f.rmdup_btn configure -state disabled
    }
    
    # Enable / disable merge-specific options
    if { $job eq "merge" } {
	global $f.new_contigs
	if {[set $f.new_contigs]} {
	    $f.repad_btn configure -state disabled
	} else {
	    $f.repad_btn configure -state normal
	}
    }
}

proc ImportSequences2 { old_io f job } {
    global io $f.refpos $f.rmdup

    # Get settings
    if { $job == "new" } {
	set merge 0
	set repad 0
    } else {
	global $f.new_contigs $f.repad
	set merge [ expr ! [set $f.new_contigs] ]
	set repad [set $f.repad]
    }
    set refpos [set $f.refpos]
    set rmdup  [set $f.rmdup]

    # Format codes:
    # ace=A caf=C sam=s bam=b auto=a
    set fmt [lindex "a A C s b a" [radiolist_get $f.format]]

    SetBusy
    if { $old_io != "" } {
	if {![quit_displays -io $io -msg "break contig"]} {
	    # Someone's too busy to shutdown?
	    ClearBusy
	    return
	}
    }

    set infile [getFname_in_name $f.infile]
    if { ![file readable $infile]} {
	tk_messageBox -icon error -type ok -title "Can't read file" \
		-message "Couldn't open $infile"
	ClearBusy
	return
    }
    
    if { $job == "new"} {
	set new_io [New_Open_Database $f [$f.database.entry get] \
			$f.browsed_path]
	if { $new_io == "" } {
	    ClearBusy
	    return
	}
	set io $new_io
    }

    vmessage "Importing reads..."

    log_call import_reads              \
	-io $io               \
	-file $infile         \
	-format $fmt          \
	-append 1             \
	-merge_contigs $merge \
	-repad $repad         \
	-store_refpos $refpos \
	-remove_dups $rmdup

    vmessage "Flushing..."
    $io flush

    vmessage "Import finished."

    ClearBusy
    PostLoadSetup

    destroy $f
}
