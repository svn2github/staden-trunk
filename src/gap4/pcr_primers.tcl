#
# FIXME:
# * PCR Primers within a single contig (amplify a region)
# 
# * Primers from non-genomic DNA. Eg from a subclone
#
# * Produce text output to feed into Karen's primer DB program
# 
#

package require tablelist
load_package prefinish

namespace eval pcr_primers {
    global env

    variable checkedImg
    variable uncheckedImg
    variable instance 0

    set checkedImg [image create photo -file $env(STADTABL)/checked.gif]
    set uncheckedImg [image create photo -file $env(STADTABL)/unchecked.gif]
}

# -----------------------------------------------------------------------------
# The main entry point and the only public function
proc ::pcr_primers::GUI {io} {
    global gap_defs
    array set pdefs [get_primer_defs]

    set t .pcr_primers

    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Pick PCR primers"

    # List or File
    lorf_in $t.infile [keylget gap_defs PCR_PRIMERS.INFILE] \
	"" -bd 2 -relief groove


    # Manually selected contig pair
    set f [frame $t.id -bd 2 -relief groove]
    contig_id $f.id1 -io $io -range 0 -default "" -trace 2 \
	-title "Left Contig Identifier" 
    contig_id $f.id2 -io $io -range 0 -default "" -trace 0 \
	-title "Right Contig Identifier"
    $f.id1 configure -bd 0
    $f.id2 configure -bd 0
    pack $f.id1 $f.id2 -side top -fill both


    # Radiolist to pick between pair or list/file
    radiolist $t.task \
	-title "Select task:" \
	-bd 2 \
	-relief groove \
	-orient horizontal \
	-default [keylget gap_defs PCR_PRIMERS.TASK] \
	-buttons [format { \
	    {{Specific pair} -command \
		{contig_id_configure %s -state normal; \
		 contig_id_configure %s -state normal;
		 lorf_in_configure %s -state disabled}} \
	    {{Neighbouring contigs} -command  \
		{contig_id_configure %s -state disabled; \
		 contig_id_configure %s -state disabled; \
		 lorf_in_configure %s -state normal}} \
	 } [list $t.id.id1] [list $t.id.id2] [list $t.infile] \
	   [list $t.id.id1] [list $t.id.id2] [list $t.infile] \
           ]


    # Ranges for where to find primers
    frame $t.search -bd 2 -relief groove
    label $t.search.label -text "Offsets from contig ends"
    scalebox $t.search.start \
	-orient horizontal -width 5 \
        -title   "[keylget gap_defs PCR_PRIMERS.SEARCH_START.NAME]" \
        -from    "[keylget gap_defs PCR_PRIMERS.SEARCH_START.MIN]" \
        -to      "[keylget gap_defs PCR_PRIMERS.SEARCH_START.MAX]" \
        -default "[keylget gap_defs PCR_PRIMERS.SEARCH_START.VALUE]" \
	-command "CheckStartLimits $t.search.start $t.search.end 0"

    scalebox $t.search.end \
	-orient horizontal -width 5 \
        -title   "[keylget gap_defs PCR_PRIMERS.SEARCH_END.NAME]" \
        -from    "[keylget gap_defs PCR_PRIMERS.SEARCH_END.MIN]" \
        -to      "[keylget gap_defs PCR_PRIMERS.SEARCH_END.MAX]" \
        -default "[keylget gap_defs PCR_PRIMERS.SEARCH_END.VALUE]" \
	-command "CheckEndLimits $t.search.start $t.search.end 0"
    pack $t.search.label $t.search.start $t.search.end -side top -fill both


    # General "oligo parameters" control panel
    frame $t.params -bd 2 -relief groove
    label $t.params.label -text "Oligo parameters"
    pack $t.params.label -side top -fill both
    select_oligo_params $t.params

    xentry $t.repeats -label "External repeats file"

    okcancelhelp $t.but \
	-ok_command "::pcr_primers::OK $io $t $t.params \
                     $t.task $t.infile $t.id.id1 $t.id.id2 \
                     $t.search.start $t.search.end $t.repeats" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 FIXME" \
	-bd 2 -relief groove

    SetCurFrame $t.id $t.id.id1
    bind [entrybox_path $t.id.id1.ent] <<select>> "SetCurFrame $t.id $t.id.id1"
    bind [entrybox_path $t.id.id2.ent] <<select>> "SetCurFrame $t.id $t.id.id2"
    pack $t.task $t.infile $t.id $t.search $t.params $t.repeats \
	-side top -fill both
    pack $t.but -side bottom -fill both
}

# -----------------------------------------------------------------------------
# Private functions below here

# Fishes out results from dialogue to launch the PCR picking code itself
proc ::pcr_primers::OK {io t p task infile id1 id2 start_w end_w repeats} {
    # Get contigs to work on
    if {[radiolist_get $task] == 1} {
	if {[set contig1 [contig_id_gel $id1]] == ""} return
	if {[set contig2 [contig_id_gel $id2]] == ""} return
	set list [list $contig1 $contig2]
	SetContigGlobals $io $contig1
    } else {
	if {[lorf_in_get $infile] == 3} {
	    set list [CreateAllContigList $io]
	} else {
	    set list [lorf_get_list $infile]
	}
    }
    
    # Set primer paramters from GUI
    array set pdefs [get_primer_defs]

    set pdefs(min_tm)  [entrybox_get $p.tm.min] 
    set pdefs(opt_tm)  [entrybox_get $p.tm.opt]
    set pdefs(max_tm)  [entrybox_get $p.tm.max]
    set pdefs(min_len) [entrybox_get $p.length.min]
    set pdefs(opt_len) [entrybox_get $p.length.opt]
    set pdefs(max_len) [entrybox_get $p.length.max]
    set pdefs(min_gc)  [entrybox_get $p.gc.min]
    set pdefs(opt_gc)  [entrybox_get $p.gc.opt]
    set pdefs(max_gc)  [entrybox_get $p.gc.max]
    set pdefs(gc_clamp) [yes_no_get $p.gc_clamp]

    set_primer_defs [array get pdefs]

    if {[set start [scalebox_get $start_w]] == ""} {return}
    if {[set end   [scalebox_get $end_w]  ] == ""} {return}

    # Generate a prefinish object and initialise it
    finish .f \
	-io $io \
	-check_contigs [CreateAllContigList $io] \
	-dust_level 14 \
	-pcr_offset1 $end \
	-pcr_offset2 $start \
	-pwalk_max_match 12 \
	-debug0 2 \
	-debug1 2 \
	-debug2 2 \
	-debug3 2 \
	-debug4 2 \
	-debug5 2 \
	-debug6 2 \
	-debug7 2 \
	-debug8 2 \
	-debug9 2
	
    # Hash external repeat file
    set repfile [$repeats get]
    if {$repfile != ""} {
	.f configure -external_seq_file $repfile
    }

    # Run it
    destroy $t
    update idletasks

    set x [.f pick_pcr_primers \
	       -contigs $list \
	       -p_args [get_primer_defs]]

    create_pcr_list $io $x
}


# Creates the dialogue showing the PCR results
proc ::pcr_primers::create_pcr_list {io pcr_details} {
    variable instance
    variable checkedImg
    variable uncheckedImg
    variable RowInfo
    
    # Create the GUI framework
    incr instance
    set t .pcr_primers$instance
    xtoplevel $t
    
    # Ok/Help buttons
    set f [frame $t.buttons]
    button $f.ok -text OK -command "destroy $t"
    button $f.help -text Help -command "show_help gap4 FIXME"
    pack $f.ok $f.help -side left -expand 1 -anchor c
    pack $f    -side bottom -fill both

    # Text info panel
    text $t.info -height 4 -state disabled
    pack $t.info -side bottom -fill both


    # Tablelist
    pack [set f [frame $t.list]] -fill both -expand 1
    scrollbar $f.y -orient vertical -command "$f.l yview"
    set tl $f.l

    tablelist::tablelist $tl \
	-width 0 -height 15 \
	-columns {
	    0 "Use" 0 "Score"
	    0 "Contig" 0 "Pos" 0 "Tm" 0 "Primer"
	    0 "Contig" 0 "Pos" 0 "Tm" 0 "Primer"} \
	-labelcommand "$tl finishediting; tablelist::sortByColumn" \
	-exportselection 0 \
	-stretch 0 \
	-selectmode single \
	-editendcommand ::pcr_primers::editEndCmd \
	-yscrollcommand "$f.y set"

    $tl columnconfigure 0 \
	-sortmode integer \
	-editable yes \
	-editwindow checkbutton \
	-formatcommand ::pcr_primers::emptystr

    for {set c 2} {$c <= 5} {incr c} {
	$tl columnconfigure $c -bg bisque
    }
    for {set c 6} {$c <= 9} {incr c} {
	$tl columnconfigure $c -bg cornsilk
    }

    pack $f.y -side right -fill both
    pack $f.l -fill both -expand 1 -side left

    # Populate the tableist
    foreach p $pcr_details {
	foreach {x l r} $p {}
	foreach {use score} $x {break}
	foreach {lPrim lcon lpos llen lTm lGC} $l {break}
	foreach {rPrim rcon rpos rlen rTm rGC} $r {break}
	set lTm [format "%4.1f" $lTm]
	set rTm [format "%4.1f" $rTm]
	set lcon [left_gel $io $lcon]
	set rcon [left_gel $io $rcon]
	set line [list $use $score \
		      $lcon $lpos $lTm $lPrim \
		      $rcon $rpos $rTm $rPrim]
	$tl insert end $line 
	$tl cellconfigure end,0 \
	    -image [expr {$use ? $checkedImg : $uncheckedImg}]
	set RowInfo([lrange $line 1 end]) $p
    }
    set RowInfo() ""

    # Bindings.
    bind [$tl bodypath] <<select>> "
      set y \[$tl nearest \[expr {%y + \[winfo y %W\]}\]\]
      set line \[lrange \[$tl get \$y\] 1 end\]
      ::pcr_primers::update_info $io \$::pcr_primers::RowInfo(\$line) $t.info
    "

    bind [$tl bodypath] <<use>> "
      set y \[$tl nearest \[expr {%y + \[winfo y %W\]}\]\]
      set line \[lrange \[$tl get \$y\] 1 end\]
      ::pcr_primers::invoke_editors $io \$::pcr_primers::RowInfo(\$line)
    "
}

# Callback from the tablelist
proc ::pcr_primers::editEndCmd {tbl row col text} {
    variable checkedImg
    variable uncheckedImg

    set img [expr {$text ? $checkedImg : $uncheckedImg}]
    $tbl cellconfigure $row,$col -image $img
    return $text
}

proc ::pcr_primers::emptystr {args} { return "" }


# Called by the binding on the tableist to update the text information box
proc ::pcr_primers::update_info {io text w} {
    if {$text == ""} {
	bell
	return
    }

    # Compute the information lines
    foreach {x l r} $text {}
    foreach {use qual comp difftm prodtm proddifftm} $x {}

    set comp [expr {$comp/100}]
    set info ""
    append info [format "PRODUCT   Quality:  %9.3f     Product Tm:  %5.2f'C\n" \
		     $qual $prodtm]
    append info [format "PRODUCT   Self-compl:%7.2f      L/R Tm diff: %5.2f'C     Prod-oligo Tm diff:    %5.2f'C\n" \
		     $comp $difftm $proddifftm]

    foreach {seq contig pos len tm gc} $l {break}
    set c [io_read_contig $io $contig]
    set dist [expr {[keylget c length]-($pos+$len-1)}]
    append info [format "LEFT      GC %%:        %5.2f%%     Tm:          %5.2f'C     End distance: %14d     %s\n" \
		     $gc $tm $dist [left_gel $io $contig]]

    foreach {seq contig pos len tm gc} $r {break}
    append info [format "RIGHT     GC %%:        %5.2f%%     Tm:          %5.2f'C     End distance: %14d     %s" \
		     $gc $tm $pos [left_gel $io $contig]]

    $w configure -state normal
    $w delete 1.0 end
    $w insert end $info
    $w configure -state disabled
}


# Called by double clicking on the tablelist.
# This brings up a pair of contig editors showing the location of the PCR
# primers.
proc ::pcr_primers::invoke_editors {io text} {
    if {$text == ""} {
	bell
	return
    }

    foreach {x l r} $text {break}
    foreach {- contig(0) pos(0) len(0) - -} $l {break}
    foreach {- contig(1) pos(1) len(1) - -} $r {break}

    # Create editors
    set ed0 [edit_contig \
		 -io $io \
		 -contig [left_gel $io $contig(0)] \
		 -pos $pos(0) \
		 -reuse 1 \
		 -nojoin 1
	     ]

    set ed1 [edit_contig \
		 -io $io \
		 -contig [left_gel $io $contig(1)] \
		 -pos $pos(1) \
		 -reuse 1 \
		 -nojoin 1
	     ]

    # Temporarily highlight the primers
    $ed0 create_tmp_anno 0 $pos(0) $len(0) OLIG left 0
    $ed1 create_tmp_anno 0 $pos(1) $len(1) OLIG right 1

    # Move editors to stack on top of each other
    wm geometry [winfo toplevel $ed0] +0+0
    update idletasks
    set geom [wm geometry [winfo toplevel $ed0]]
    regexp {[0-9]+x[0-9]+\+([0-9]+)\+([0-9]+)} $geom _ xpos ypos
    set height [winfo height [winfo toplevel $ed0]]
    set bd [winfo rooty [winfo toplevel $ed0]]
    puts "wm geometry [winfo toplevel $ed1] +0+[expr {$ypos+$height+$bd}]"
    wm geometry [winfo toplevel $ed1] +0+[expr {$ypos+$height+$bd}]
    update idletasks
    puts done
}
