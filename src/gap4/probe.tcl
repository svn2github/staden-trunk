#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc SuggestProbes {io} {
    global gap_defs

    set l [keylget gap_defs FIND_PROBES]
    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Suggest probes"

    contig_id $t.id -io $io -range 0

    lorf_in $t.infile [keylget gap_defs FIND_PROBES.INFILE] \
        "{contig_id_configure $t.id -state disabled}
         {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
         {contig_id_configure $t.id -state normal}
	" -bd 2 -relief groove

    scalebox $t.maxperc \
	-orient horizontal -width 5 \
	-title   "[keylget l MAXPERC.NAME]" \
	-from    "[keylget l MAXPERC.MIN]" \
	-to      "[keylget l MAXPERC.MAX]" \
	-default "[keylget l MAXPERC.VALUE]" \
	-bd 2 -relief groove

    scale_range $t.size -bd 2 -relief groove \
	-start_value [keylget l OLIGOSIZE.MIN_DEF] \
	-end_value   [keylget l OLIGOSIZE.MAX_DEF] \
	-start_name  [keylget l OLIGOSIZE.MIN_NAME] \
	-end_name    [keylget l OLIGOSIZE.MAX_NAME] \
	-from        [keylget l OLIGOSIZE.MIN] \
	-to          [keylget l OLIGOSIZE.MAX]

    scale_range $t.range -bd 2 -relief groove \
	-start_value [keylget l OLIGOPOS.MIN_DEF] \
	-end_value   [keylget l OLIGOPOS.MAX_DEF] \
	-start_name  [keylget l OLIGOPOS.MIN_NAME] \
	-end_name    [keylget l OLIGOPOS.MAX_NAME] \
	-from        [keylget l OLIGOPOS.MIN] \
	-to          [keylget l OLIGOPOS.MAX]

    eval getFname \{$t.vectors\} \{[keylget l VECTOR.NAME]\} load_optional {} \
	[keylget l VECTOR.VALUE]

    okcancelhelp $t.ok \
	-ok_command "SuggestProbes2 $io $t $t.infile $t.id $t.maxperc \
		     $t.size $t.range $t.vectors" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Suggest Probes}" \
	-bd 2 -relief groove

    pack $t.infile $t.id $t.maxperc $t.size $t.range $t.vectors $t.ok \
	-side top -fill x
}

proc SuggestProbes2 {io t infile id maxperc oligosize position vectors} {
    set args ""
    if {[set maxp [scalebox_get $maxperc]]       == ""} {bell; return}
    if {[set mins [scale_range_from $oligosize]] == ""} {bell; return}
    if {[set maxs [scale_range_to $oligosize]] == ""} {bell; return}
    if {[set from [scale_range_from $position]]  == ""} {bell; return}
    if {[set to   [scale_range_to $position]]  == ""} {bell; return}

    # Single contig
    if {[lorf_in_get $infile] == 4} {
	if {[set id [contig_id_gel $id]] == ""} {
	    bell; contig_id_focus; return
       	}

	lappend args -contigs $id

    # All contigs
    } elseif {[lorf_in_get $infile] == 3} {
	lappend args -contigs [CreateAllContigList $io]
    
    # List or File of contigs
    } else {
	if {[set in_l [lorf_get_list $infile]] == ""} {
	    bell; lorf_focus $infile; return
	}
	lappend args -contigs $in_l
    }

    set vectors [getFname_in_name $vectors]
    if {$vectors != ""} {
	lappend args -vectors $vectors
    }

    SetBusy
    set ret [eval find_probes \
            -io $io \
            -max_pmatch $maxp \
            -min_size $mins \
            -max_size $maxs \
            -from     $from \
            -to       $to \
		 -primer_arg [list [get_primer_defs]] \
	    $args]
    ClearBusy
    destroy $t

    probes_display $ret $t $io
}

proc probes_display {results toplevel io} {
    set text [probe_create_win $toplevel $io]
    
    foreach cend $results {
        set end [lindex $cend 2]
        set cnum [lindex $cend 1]
        $text insert end \
            "\nContig [lindex $cend 0]($cnum): $end\n" title
        set selected "selected"
	if {[llength [lindex $cend 3]] == 0} {
	    $text insert end "    No oligos found\n"
	}
        foreach oligo [lindex $cend 3] {
            $text insert end [format \
		"    Pos %6d, Dist %3d, primer=%2.0f, Tm=%2.0f, match=%02.0f%%, %s" \
		[lindex $oligo 0] \
		[lindex $oligo 1] \
		[lindex $oligo 2] \
		[lindex $oligo 3] \
		[lindex $oligo 4] \
		[lindex $oligo 5]] "selectable c_$cnum e_$end $selected"
            $text insert end "\n"
            set selected ""
        }
        if {[lindex $cend 4] > 0} {
            $text insert end \
                 "    Rejected [lindex $cend 4] oligos due to non uniqueness\n"
        }
    }
}
    

#
# Scan through selected items outputting them to disk and adding tags
#
proc probe_use_selections {t io text outfile add_tags} {
    global gap_defs

    set file [getFname_in_name $outfile]
    set tag_list ""

    if {$file != ""} {
	if {[set fd [open $file w+ 0666]] == -1} {
	    bell
	    return
	}
    }

    set tag_type [keylget gap_defs FIND_PROBES.TAG_TYPE]

    set toggle 0
    foreach range [$text tag ranges selected] {
	if {$toggle == 0} {
	    set toggle $range
	} else {
	    set tags [$text tag names $toggle]
	    if {[set tc [lsearch -regexp $tags {^c_}]] != -1 &&
		[set te [lsearch -regexp $tags {^e_}]] != -1} {
		set c [string range [lindex $tags $tc] 2 end]
		scan [$text get $toggle $range] \
		    "    Pos %d, Dist %d, primer=%d, Tm=%d, match=%d%%, %s" \
			p d o tm m s
		if {$file != ""} {
		    puts $fd "Contig [r_name $io $c]   position $p   Tm $tm   sequence $s"
		}
		if {$add_tags} {
		    if {[string range [lindex $tags $te] 2 end] == "Start"} {
	        	set e "+"
			set st $p
			set en [expr $p+[string length $s]-1]
		    } else {
			set e "-"
			set st [expr $p-[string length $s]+1]
			set en $p
		    }
		    set tt "-[db_info get_contig_num $io #$c]\
			$tag_type\
			$e\
			$st..$en\
			\nscore=$o\nbest_match=$m%"
		    lappend tag_list $tt
		}
	    } else {
		puts "Failed to find tag for selected item"
	    }
	    set toggle 0
	}
    }

    if {$file != ""} {
	close $fd
    }
    if {$add_tags && $tag_list != ""} {
	add_tags -io $io -tags $tag_list
    }

    destroy $t
}


#
# Create the probe windows for selecting the oligos to use
#
proc probe_create_win {toplevel io} {
    global gap_defs read_only

    set l [keylget gap_defs FIND_PROBES]

    if {"$toplevel" != ""} {
	if {[winfo exists $toplevel]} {destroy $toplevel}
	xtoplevel $toplevel
	wm title $toplevel "Probes found"
    }

    # Command buttons
    frame [set b $toplevel.f] -bd 2 -relief groove
    global $b.Value
    frame $toplevel.f2
    set g [grab current]
    if {$g != ""} {
	set gcmd "grab $g"
    } else {
	set gcmd "grab release $toplevel"
    }
    okcancelhelp $toplevel.f2.ok \
	-ok_command "$gcmd; \
		     probe_use_selections $toplevel $io \
			$toplevel.t $b.file \[set $b.Value\]" \
	-cancel_command "destroy $toplevel" \
	-help_command "show_help gap4 {Suggest Probes}" \
	-bd 2 -relief groove
    eval getFname $b.file \{[keylget l OUTFILE.NAME]\} save_optional {} \
	[keylget l OUTFILE.VALUE]
    checkbutton $b.tags \
	-text [keylget l TAG.NAME] \
	-variable $b.Value \
	-bd 2 -relief raised -padx 5
    set $b.Value [keylget l TAG.VALUE]

    if {$read_only} {
	set $b.Value 0
	$b.tags configure -state disabled
    }

    pack $toplevel.f2 -side bottom -fill both
    pack $toplevel.f2.ok -side left -fill both -expand 1
    pack $b -side bottom -fill both 
    pack $b.file -side left -fill both
    pack $b.tags -side right -fill both

    # Text widget and scrollbar
    text $toplevel.t -width 80 \
            -yscrollcommand "$toplevel.s set"
    scrollbar $toplevel.s -orient vert -command "$toplevel.t yview"
    pack $toplevel.s -side right -fill both
    pack $toplevel.t -side bottom -fill both -expand 1

    # Add our text tags
    $toplevel.t tag configure selected -foreground blue
    $toplevel.t tag configure selectable
    $toplevel.t tag configure highlight -relief raised -borderwidth 2
    $toplevel.t tag configure title \
	-font title_font
    
    $toplevel.t tag bind selectable <<select>> {
        set lastLine [%W index {@%x,%y linestart}]
        if {[lsearch -exact [%W tag names [%W index @%x,%y]] selected] == -1} {
            %W tag add selected $lastLine "$lastLine lineend"
        } else {
            %W tag remove selected $lastLine "$lastLine lineend"
        }
    }
    
    $toplevel.t tag bind selectable <Any-Enter> {
        set lastLine [%W index {@%x,%y linestart}]
        %W tag add highlight $lastLine "$lastLine lineend"
    }
    
    $toplevel.t tag bind selectable <Any-Motion> {
        set lastLine [%W index {@%x,%y linestart}]
        %W tag remove highlight 1.0 end
        %W tag add highlight $lastLine "$lastLine lineend"
    }
    
    $toplevel.t tag bind selectable <Any-Leave> {
        %W tag remove highlight 1.0 end
    }
    
    # Remove most default Text bindings
    bindtags $toplevel.t "$toplevel.t . all"
    bind $toplevel.t <2> [bind Text <2>]
    bind $toplevel.t <B2-Motion> [bind Text <B2-Motion>]
    bind $toplevel.t <ButtonRelease-2> [bind Text <ButtonRelease-2>]

    grab $toplevel
    return $toplevel.t
}
