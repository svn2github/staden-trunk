#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#-----------------------------------------------------------------------------
# The parameter adjustments for the contig editor "Select Oligo" command.
#-----------------------------------------------------------------------------

proc oligo_param {f search_ahead search_back read_length} {
    global $search_ahead $search_back $read_length

    set t $f.oligo_param
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Oligo parameters"

    entrybox $t.e1 \
	-title "Search window bases ahead" \
	-default [set $search_ahead] \
	-width 5 \
	-type {CheckIntRange 0 1000}

    entrybox $t.e2 \
	-title "Search window bases back" \
	-default [set $search_back] \
	-width 5 \
	-type {CheckIntRange 0 1000}

    entrybox $t.e4 \
	-title "Average read length" \
	-default [set $read_length] \
	-width 5 \
	-type {CheckIntRange 1 5000}

    frame $t.p
    select_oligo_params $t.p

    okcancelhelp $t.but \
	    -bd 2 -relief groove \
	    -ok_command "oligo_paramok $t $search_ahead $search_back \
	                 $read_length; \
			 select_oligo_params_OK_Pressed $t.p; \
			 destroy $t" \
	    -cancel_command "destroy $t" \
	    -help_command "show_help gap4 {Editor-Primer Selection}"
    pack $t.but $t.p $t.e4 $t.e2 $t.e1 -side bottom -fill x
}


# Obtains and checks the most generally used (ie not constraints) parameters
# from the oligo code.
#
proc oligo_paramok {t search_ahead search_back read_length} {
    global $search_ahead $search_back $read_length gap_defs

    # Get and check for validity our parameters
    if {[set search_ahead_val [entrybox_get $t.e1]] == ""} {
	entrybox_focus $t.e1; return
    }

    if {[set search_back_val [entrybox_get $t.e2]] == ""} {
	entrybox_focus $t.e2; return
    }

    if {[set read_length_val [entrybox_get $t.e4]] == ""} {
	entrybox_focus $t.e4; return
    }

    # Set our values
    set $search_ahead $search_ahead_val
    set $search_back  $search_back_val
    set $read_length  $read_length_val

    keylset gap_defs SELECT_OLIGOS.SEARCH_AHEAD $search_ahead_val
    keylset gap_defs SELECT_OLIGOS.SEARCH_BACK  $search_back_val
    keylset gap_defs SELECT_OLIGOS.READ_LENGTH  $read_length_val
}

# Returns the primer_defs Tcl global from the PRIMER keyed list.
# The returned list is in 'array get' syntax
proc get_primer_defs {} {
    global gap_defs

    set primer_defs ""
    foreach i [keylget gap_defs PRIMER] {
	lappend primer_defs [lindex $i 0] [lindex $i 1]
    }
    
    return $primer_defs
}

# Takes a list of primer picking parameters in 'array get' syntax and stores
# this back in the gap_defs keyedlist syntax.
proc set_primer_defs {defs} {
    global gap_defs

    eval keylset primer_defs $defs
    keylset gap_defs PRIMER $primer_defs
}

#-----------------------------------------------------------------------------
# The "Select oligos" function within the contig editor
#-----------------------------------------------------------------------------
proc select_oligos {ed} {
    global oligo_arr gap_defs

    set t $ed.oligo
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Select primers"

    # Initialise defaults
    global $t.Direction $t.SearchAhead $t.SearchBack $t.ReadLength
    set $t.Direction   [keylget gap_defs SELECT_OLIGOS.DIRECTION]
    set $t.SearchAhead [keylget gap_defs SELECT_OLIGOS.SEARCH_AHEAD]
    set $t.SearchBack  [keylget gap_defs SELECT_OLIGOS.SEARCH_BACK]
    set $t.ReadLength  [keylget gap_defs SELECT_OLIGOS.READ_LENGTH]

    # The direction arrow
    frame $t.top
    label $t.top.name -text "Direction:"
    button $t.top.value -command "seloli_set_dir $t.Direction $t.top.value"
    seloli_set_dir $t.Direction $t.top.value 0

    # Pops up the parameter adjustment window(s)
    button $t.top.edit_params -text "Edit parameters" \
	-command "oligo_param $t $t.SearchAhead $t.SearchBack $t.ReadLength"

    # A 'doit' command
    button $t.top.find_oligos -text "Find oligos" -command "seloli_doit $ed $t"

    # An entrybox for the template name together with a menu of suggested
    # templates to choose from.
    frame $t.tname -bd 2 -relief groove
    entrybox $t.tname.name \
	-title "Template name" \
	-width 20
    xmenubutton $t.tname.choice \
	-text "Choose from >>" \
	-menu $t.tname.choice.menu
    menu $t.tname.choice.menu -tearoff 0

    # The bottom buttons "prev,next, accept, cancel"
    frame $t.bot -bd 2 -relief groove
    button $t.bot.prev   -text "Prev"   -command "seloli_next   $ed $t prev" \
       	-state disabled
    button $t.bot.next   -text "Next"   -command "seloli_next   $ed $t" \
       	-state disabled
    button $t.bot.accept -text "Accept" -command "seloli_accept $ed $t" \
	-state disabled
    button $t.bot.cancel -text "Cancel"   -command "seloli_quit   $ed $t"
    button $t.bot.help -text "Help" \
	-command "show_help gap4 {Editor-Comm-Primer Selection}"
    
	
    pack $t.top.name $t.top.value $t.top.edit_params $t.top.find_oligos \
	-side left
    pack $t.tname.name -side left -fill both
    pack $t.tname.choice -side right
    pack $t.bot.prev $t.bot.next $t.bot.accept $t.bot.cancel -side left
    pack $t.bot.help -side right

    pack $t.top $t.tname $t.bot -side top -fill both
}

proc seloli_set_dir {var button args} {
    global $var

    if {"$args" != ""} {
	set $var $args
    } else {
	set $var [expr [set $var]==0 ? 1 : 0]
    }

    $button configure -text [lindex "----> <----" [set $var]]
}

proc seloli_doit {ed t} {
    global $t.Direction $t.SearchAhead $t.SearchBack $t.ReadLength

    if {[$ed select_oligo generate "[set $t.Direction]" \
				   "[set $t.SearchAhead]" \
				   "[set $t.SearchBack]" \
				   "[set $t.ReadLength]"\
	                           [get_primer_defs]] > 0} {
        seloli_next $ed $t

       	$t.bot.prev configure -state normal
       	$t.bot.next configure -state normal
       	$t.bot.accept configure -state normal
	# editor_set_status_line $ed ""
    } else {
	bell
	verror ERR_WARN select_primers "No suitable primers found"
	editor_set_status_line $ed "No suitable primers found" 3000
    }
}

proc seloli_set_tname {t tname} {
    set e $t.tname.name

    entrybox_delete $e 0 end
    entrybox_insert $e 0 $tname
}

proc seloli_set_menu {t list} {
    set m $t.tname.choice.menu

    $m delete 0 last
    foreach i $list {
        $m add command -label $i -command "seloli_set_tname $t $i"
    }
    $m add command -label "<None>" -command "seloli_set_tname $t {}"
}

proc seloli_next {ed t {mode next}} {
    set l [$ed select_oligo $mode]

    if {"$l" != ""} {
        seloli_set_menu $t [lrange $l 1 end]
        seloli_set_tname $t [lindex $l 0]
    } else {
	bell
    }
}

proc seloli_accept {ed t} {
    set e $t.tname.name
    
    if {![ListExists2 oligos]} {
	ListCreate2 oligos ""
	ListEdit oligos
    }

    set tname [entrybox_get $e]
    ListAppend oligos [list [$ed select_oligo accept "$tname"]]

    seloli_quit $ed $t
}

proc seloli_quit {ed t} {
    global $t.Direction $t.SearchAhead $t.SearchBack $t.ReadLength

    $ed select_oligo quit
    destroy $t

    unset $t.Direction $t.SearchAhead $t.SearchBack $t.ReadLength
}

#-----------------------------------------------------------------------------
# The "Suggest Primers" functions on the main menu
#-----------------------------------------------------------------------------
proc SuggestPrimers {io} {
    global gap_defs

    #initialise defaults
    set l [keylget gap_defs SUGGEST_PRIMERS]
    set t [keylget l WIN]

    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Suggest primers"
    
    contig_id $t.id -io $io -range 0
    
    lorf_in $t.infile [keylget gap_defs SUGGEST_PRIMERS.INFILE] \
	"{contig_id_configure $t.id -state disabled}
         {contig_id_configure $t.id -state disabled}
	 {contig_id_configure $t.id -state disabled}
         {contig_id_configure $t.id -state normal}
	" -bd 2 -relief groove

    lorf_out $t.outfile [keylget gap_defs SUGGEST_PRIMERS.OUTFILE] \
	"" -bd 2 -relief groove

    frame $t.search -bd 2 -relief groove
    label $t.search.label -text "Offset from problem to analyse"
    scalebox $t.search.start \
	-orient horizontal -width 5 \
        -title   "[keylget l SEARCH_START.NAME]" \
        -from    "[keylget l SEARCH_START.MIN]" \
        -to      "[keylget l SEARCH_START.MAX]" \
        -default "[keylget l SEARCH_START.VALUE]" \
	-command "CheckStartLimits $t.search.start $t.search.end 0"

    scalebox $t.search.end \
	-orient horizontal -width 5 \
        -title   "[keylget l SEARCH_END.NAME]" \
        -from    "[keylget l SEARCH_END.MIN]" \
        -to      "[keylget l SEARCH_END.MAX]" \
        -default "[keylget l SEARCH_END.VALUE]" \
	-command "CheckEndLimits $t.search.start $t.search.end 0"

    entrybox $t.start_num \
	-title   "[keylget l SEARCH_NUM.NAME]" \
	-default "[keylget l SEARCH_NUM.VALUE]" \
	-width 5 \
	-type "CheckIntMin 1"

    entrybox $t.num_primers \
	-title   "[keylget l NUM_PRIMERS.NAME]" \
	-default "[keylget l NUM_PRIMERS.VALUE]" \
	-width 5 \
	-type "CheckIntMin 1"

    frame $t.p -bd 2 -relief groove
    button $t.p.but -text "Edit parameters" \
	-command "select_oligo_params $t.params" 

    okcancelhelp $t.ok \
        -ok_command "SuggestPrimers2 $io $t $t.infile $t.id $t.outfile \
	           $t.search.start $t.search.end $t.start_num $t.num_primers" \
        -cancel_command "destroy $t" \
	-help_command "show_help gap4 {Suggest Primers}" \
        -bd 2 -relief groove

    pack $t.search.label $t.search.start $t.search.end -side top -fill both
    pack $t.p.but -side left
    pack $t.infile $t.id $t.outfile $t.search $t.num_primers $t.start_num \
	$t.p $t.ok -side top -fill both
}

proc SuggestPrimers2 {io t infile_w id_w outfile_w start_w end_w \
	stnum_w numol_w} {
    if {[set start [scalebox_get $start_w]] == ""} {return}
    if {[set end [scalebox_get $end_w]] == ""} {return}

    if {[set stnum [entrybox_get $stnum_w]] == ""} {
	entrybox_focus $stnum_w; return
    }

    if {[set numol [entrybox_get $numol_w]] == ""} {
	entrybox_focus $numol_w; return
    }

    if {[lorf_in_get $infile_w] == 4} {
	if {[set id [contig_id_gel $id_w]] == ""} {
	    contig_id_focus; return
       	}

#      	if {[set lreg  [contig_id_lreg $id_w]] == ""} {return}
#	if {[set rreg  [contig_id_rreg $id_w]] == ""} {return}
	SetContigGlobals $io $id
	
#	set list "{$id $lreg $rreg}"
	set list "{$id}"
    } elseif {[lorf_in_get $infile_w] == 3} {
	set list [CreateAllContigList $io]
    } else {
	if {[set list [lorf_get_list $infile_w]] == ""} {
	    lorf_focus $infile_w; return
	}
    }

    if {[set out_l  [lorf_out_name $outfile_w]] == ""} {return}
    if {[set format [lorf_out_get  $outfile_w]] == ""} {return}

    destroy $t
    update idletasks

    SetBusy
    ListCreate2 $out_l [find_primers \
			    -io $io \
			    -search_from $start \
			    -search_to $end \
			    -num_primers $numol \
			    -primer_start $stnum \
			    -params [get_primer_defs] \
			    -contigs $list]
    ClearBusy

    if {"$format" == 2} {
	lorf_out_save $out_l
    }
}

proc select_oligo_params { f } {
    global gap_defs

    array set pdefs [get_primer_defs]

    if {![winfo exists $f]} {
	xtoplevel $f -resizable 0
	wm title $f "Primer parameters"
    }

    frame $f.tm -bd 0 -relief groove
    label $f.tm.label -text [keylget gap_defs SUGGEST_PRIMERS.PARAMS.TM.NAME]
    entrybox $f.tm.min \
	-default $pdefs(min_tm) \
	-title "Min" \
	-width 5 \
	-type {CheckFloatRange 0 100}

    entrybox $f.tm.opt \
	-default $pdefs(opt_tm) \
	-title "Opt" \
	-width 5 \
	-type {CheckFloatRange 0 100}

    entrybox $f.tm.max \
	-default $pdefs(max_tm) \
	-title "Max" \
	-width 5 \
	-type {CheckFloatRange 0 100}

    pack $f.tm.label -side left 
    pack $f.tm.max $f.tm.opt $f.tm.min -side right 

    frame $f.length -bd 0 -relief groove
    label $f.length.label -text [keylget gap_defs SUGGEST_PRIMERS.PARAMS.LEN.NAME]
    entrybox $f.length.min \
	-default $pdefs(min_len) \
	-title "Min"\
	-width 5 \
	-type {CheckIntRange 1 100}

    entrybox $f.length.opt \
	-default $pdefs(opt_len) \
	-title "Opt"\
	-width 5 \
	-type {CheckIntRange 1 100}

    entrybox $f.length.max \
	-default $pdefs(max_len) \
	-title "Max" \
	-width 5 \
	-type {CheckIntRange 1 100}

    pack $f.length.label -side left
    pack $f.length.max $f.length.opt $f.length.min -side right 

    frame $f.gc -bd 0 -relief groove
    label $f.gc.label -text "GC content (%)"
    entrybox $f.gc.min \
	-default $pdefs(min_gc) \
	-title "Min"\
	-width 5 \
	-type {CheckIntRange 1 100}

    entrybox $f.gc.opt \
	-default $pdefs(opt_gc) \
	-title "Opt"\
	-width 5 \
	-type {CheckIntRange 1 100}

    entrybox $f.gc.max \
	-default $pdefs(max_gc) \
	-title "Max" \
	-width 5 \
	-type {CheckIntRange 1 100}

    yes_no $f.gc_clamp \
    	    -title "GC Clamp" \
	    -orient horizontal \
	    -bd 0 \
	    -default $pdefs(gc_clamp)

    pack $f.gc.label -side left 
    pack $f.gc.max $f.gc.opt $f.gc.min -side right 

    pack $f.tm $f.length $f.gc $f.gc_clamp -fill x

    if {[winfo class $f] == "Toplevel"} {
	frame $f.but -bd 2 -relief raised
	button $f.but.ok -text OK \
		-command "select_oligo_params_OK_Pressed $f; destroy $f" 
	button $f.but.cancel -text Cancel -command "destroy $f"
	pack $f.but.ok $f.but.cancel -side left -expand 1 -anchor c
	pack $f.but -side bottom -side left -fill x -expand 1
    }

}

proc select_oligo_params_OK_Pressed {p} {
    global gap_defs

    array set pdefs [get_primer_defs]

    set pdefs(min_tm) [entrybox_get $p.tm.min] 
    set pdefs(opt_tm) [entrybox_get $p.tm.opt]
    set pdefs(max_tm) [entrybox_get $p.tm.max]
    
    set pdefs(min_len) [entrybox_get $p.length.min]
    set pdefs(opt_len) [entrybox_get $p.length.opt]
    set pdefs(max_len) [entrybox_get $p.length.max]
    
    set pdefs(min_gc) [entrybox_get $p.gc.min]
    set pdefs(opt_gc) [entrybox_get $p.gc.opt]
    set pdefs(max_gc) [entrybox_get $p.gc.max]

    set pdefs(gc_clamp) [yes_no_get $p.gc_clamp]

    set_primer_defs [array get pdefs]    
}
