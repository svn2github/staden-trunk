#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc SeqedDisplay_d { } {
    
    set t .seqed_display_d
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Sequence display"

    #set default as first sequence in list
    #set active_id [get_seq_id -seq_num 0]
############################################################
    #set default as horizontal sequence in list

    set sequences [sequence_names]

    foreach i $sequences {
	set dir [lindex $i 0]
	if {[llength $dir] == 1} {
	    if {$dir == "H"} {
		set hor_name [lindex $i 1]
	    }
	}
    }
    set active_id [name_to_seq_id $hor_name]
############################################################
    seq_id $t.seq_id -range 0 -default [seq_info $active_id name]\
	    -browse_cmd seq_browser

    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "SeqedDisplay_d2 $t $t.seq_id"\
	-cancel_command "seq_id_destroy $t.seq_id; destroy $t" \
	-help_command "show_help spin {SPIN-Sequence-Display}"

    pack $t.seq_id -fill x
    pack $t.button -fill x
}

proc SeqedDisplay_d2 {t seq_id} {
    global HORIZONTAL

    set id [name_to_seq_id [seq_id_name $seq_id]]
    set pos [seq_info $id start]

    SeqedDisplay $pos $id
    seq_id_destroy $seq_id
    destroy $t
}

proc SeqedDisplay { pos seq_id} {
    global se_id tk_utils_defs

    if {![info exists se_id]} {
	set se_id 0
    } else {
	incr se_id
    }

    set f .graph_sequence_display$se_id
    if {[xtoplevel $f] == ""} return
 #   if {[xtoplevel $f -resizable 0] == ""} return

    set font_height [font metrics sheet_font -linespace]
    set font_width [font measure sheet_font "A"]

#    fix_maxsize_text $f $font_width $font_height 40 80

    create_seqed $f $f.se $pos $seq_id

    # Trap window quit so that clean-up is done properly
    wm protocol $f WM_DELETE_WINDOW "seqed_quit $f $f.se $f.REnz"

    #fill in widget info and register
    set seqed_id [seqed_display -window $f.se -seq_id $seq_id]

    set name [seq_info $seq_id name]
    wm title $f "Sequence Display (\#$seqed_id) $name"
    #update result list
    seq_result_list_update [keylget tk_utils_defs RASTER.RESULTS.WIN]
}

proc create_seqed { w se pos seq_id} {
    global $w.Ruler $w.Complement DNA

    set $w.Ruler 1
    set $w.Complement 0
    seqed $se -xscrollcommand "$w.sb_h set" -yscrollcommand "$w.sb_v set" \
	-height 5 -foreground black -bd 2 -relief raised -pos $pos -font sheet_font

    #seqed $se -xscrollcommand "$w.sb_h set" -yscrollcommand "$w.sb_v set" \
	    -height 5 -ruler 1 -complement 1 -foreground black \
		-bd 5 -relief raised -heightmin 5 -heightmax 20 -pos $pos

    $se ruler [set $w.Ruler]

    #only complement if sequence is dna
    if {[seq_info $seq_id type] == $DNA} {
	set $w.Complement 1
    }
    $se complement [set $w.Complement]
    scrollbar $w.sb_h -orient horizontal -command "$se xview"
    scrollbar $w.sb_v -orient vertical -command "$se yview"

    #seqednames $w.names -width 23 -height 1 -bd 2 -relief raised
 
    frame $w.buttons
    create_seqed_menus $w $se $w.buttons $seq_id
   # pack $f.names -side left -expand yes

    grid columnconfig $w 0 -weight 1
    grid rowconfig $w 1 -weight 1

    grid $w.buttons -row 0 -column 0 -sticky ew
    grid $se -row 1 -column 0 -sticky nsew
    grid $w.sb_h -row 2 -column 0 -sticky ew
    grid $w.sb_v -row 1 -column 1 -sticky ns
}

proc create_seqed_menus { w se f seq_id} {
    global spin_defs $w.AminoMode $w.REnz $w.Ruler $w.Complement DNA

    set $w.Status0	  0
    set $w.Status1	  0
    set $w.Status2	  0
    set $w.Status3	  0
    set $w.Status4	  0
    set $w.Status5	  0
    set $w.Status6	  0
    set $w.AminoMode	  [keylget spin_defs SEQED.AMINO_ACID_MODE]
    set $w.REnz           0

    xmenubutton $f.settings -relief raised -text "Settings >>"\
	    -menu $f.settings.menu
#	    -menu $w.buttons.settings.menu
    button $f.search -relief raised -text "Search"\
	    -command "SeqedSearch $w $se"
    button $f.save -relief raised -text "Save"\
	    -command "SeqedSave $w $se $seq_id"
    button $f.quit -relief raised -text "Quit"\
	    -command "seqed_quit $w $se $w.REnz"
    xmenubutton $f.help -relief raised -text "Help >>" -menu $f.help.menu -padx 2

    menu [set m1 $f.settings.menu]
    menu [set m2 $f.help.menu]

    pack $f.settings $f.search $f.save $f.quit -side left -fill both
    pack $f.help -side right

####### yc #######
    if {[seq_info $seq_id type] != $DNA} {
	$f.search configure -state disabled
    }
###### yc #######

    # Settings menu

    #display ruler
    $m1 add checkbutton -label "ruler" \
	-command "$se ruler \[set $w.Ruler\]" \
	-variable $w.Ruler

     if {[seq_info $seq_id type] == $DNA} {
	#display reverse strand
	$m1 add checkbutton -label "reverse strand" \
		-command "$se complement \[set $w.Complement\]" \
		-variable $w.Complement
   
	 #select 3 or 1 letter amino acid names
	 $m1 add checkbutton -label "3 Character Amino Acids" \
		 -command "$se translation_mode \[lindex {1 3}\
		 \[set $w.AminoMode\]\]" \
		 -variable $w.AminoMode
	 
	 $m1 add checkbutton -label "Restriction enzyme map" \
	    -command "seqed_r_enz_dialogue $w $se \[set $w.REnz\]" \
	    -variable $w.REnz
	 
	 $m1 add cascade -label "Translate" -menu $m1.translate
    
	 set c [menu $m1.translate -tearoff 0]
	 #$c add checkbutton -label "Automatic Translation" -variable $w.Trans0 \
		 -command "seqed_set_translate $se $w \[set $w.Trans0\] 0"
	 $c add checkbutton -label "Translate frame 1+" -variable $w.Trans1 \
		 -command "seqed_set_translate $se $w \[set $w.Trans1\] 1"
	 $c add checkbutton -label "Translate frame 2+" -variable $w.Trans2 \
		 -command "seqed_set_translate $se $w \[set $w.Trans2\] 2"
	 $c add checkbutton -label "Translate frame 3+" -variable $w.Trans3 \
		 -command "seqed_set_translate $se $w \[set $w.Trans3\] 3"
	 $c add checkbutton -label "Translate frame 1-" -variable $w.Trans4 \
		 -command "seqed_set_translate $se $w \[set $w.Trans4\] 4"
	 $c add checkbutton -label "Translate frame 2-" -variable $w.Trans5 \
		 -command "seqed_set_translate $se $w \[set $w.Trans5\] 5"
	 $c add checkbutton -label "Translate frame 3-" -variable $w.Trans6 \
		 -command "seqed_set_translate $se $w \[set $w.Trans6\] 6"
	 
	 $c add separator
	 $c add command -label "Translate + frames" \
		 -command "seqed_set_translate $se $w 1 1 2 3"
	 $c add command -label "Translate - frames" \
		 -command "seqed_set_translate $se $w 1 4 5 6"
	 $c add command -label "Translate all frames" \
		 -command "seqed_set_translate $se $w 1 1 2 3 4 5 6"
	 $c add separator
	 $c add command -label "Remove all" \
		 -command "seqed_set_translate $se $w 0 0 1 2 3 4 5 6"
     }
    # Help menu
    $m2 add command -label "Introduction" \
	-command "show_help spin {SPIN-Sequence-Display}"
    $m2 add separator
    $m2 add command -label "String search" \
	-command "show_help spin {SPIN-Sequence-Display-Search}"
    $m2 add command -label "Save" \
	-command "show_help spin {SPIN-Sequence-Display-Save}"

}

#
# quit the sequence editor
#
proc seqed_quit {w se renz } {
    $se quit
    destroy $w
}

proc SeqedSearch { w se } {
    global spin_defs

    set t $w.seqed_search
    if {[SeqedSearchDialog $t "SeqedSearch2 $se $t $t.direction" "destroy $t" "show_help spin {SPIN-Sequence-Display-Search}"] == -1} {
	return
    }

   #direction
    keylset di DIRECTION [keylget spin_defs SEQED.SEARCH.DIRECTION]
    set b1 [keylget di DIRECTION.BUTTON.1]
    set b2 [keylget di DIRECTION.BUTTON.2]
    radiolist $t.direction \
	-title [keylget di DIRECTION.NAME] \
	-default [keylget di DIRECTION.VALUE] \
	-orient horizontal \
	-buttons [format {{%s} {%s}} \
		      [list $b1] [list $b2]]

    pack $t.direction -before $t.strand -fill x
} 

proc SeqedSearch2 {se t direction} {
    global new_search search_string search_match search_strand search_iub

    if {[radiolist_get $direction] == 1} {
	set direct +
    } else {
	set direct -
    }

    set strand [SeqedSearch_strand $t]
    set match [SeqedSearch_match $t]
    set string [SeqedSearch_string $t]
    set use_iub [SeqedSearch_use_iub $t]

    if {![info exists new_search]} {
	set new_search 1
	set search_string $string
	set search_match $match
	set search_strand $strand
	set search_iub [SeqedSearch_use_iub $t]
    }

    #check to see if user has altered the search string
    if {[string compare $search_string $string] != 0} {
	set new_search 1
    }

    if {$search_match != $match} {
	set new_search 1
    }

   if {$search_strand != $strand} {
	set new_search 1
    }

   if {$search_iub != $use_iub} {
	set new_search 1
    }

    $se search $string $direct $strand $match $new_search $use_iub

    set new_search 0
    set search_string $string
    set search_match $match
    set search_strand $strand
    set search_iub $use_iub
}

proc SeqedSearchShutdown {se } {
    global new_search
    
    if {[info exists new_search]} {
	unset new_search 
	$se search_destroy
    }
}

#save the contents of the sequence editor to a file
proc SeqedSave {w se seq_id} {
    global spin_defs

#    set t .seqed_save
    set t $w.seqed_save

    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t Save

    set seq_length [s_length -seq_id $seq_id]

    keylset us RANGE [keylget spin_defs SEQED.SAVE.RANGE]
    scale_range $t.range -title [keylget us RANGE.NAME] \
	-from 1 -to $seq_length -start_value 1 -end_value $seq_length\
	-min_value 1

    keylset ll LINE_LEN [keylget spin_defs SEQED.SAVE.LINE_LEN]
    entrybox $t.line_len \
	-title "[keylget ll LINE_LEN.NAME]" \
	-default [keylget ll LINE_LEN.VALUE]\
	-width 5 \
	-type "CheckIntMin 1"

    keylset fn FILENAME [keylget spin_defs SEQED.SAVE.FILENAME]
    entrybox $t.name \
	-title   [keylget fn FILENAME.NAME] \
	-default [keylget fn FILENAME.VALUE] \
	-type CheckOutput


    #########################################################################
    #ok cancel help buttons 
    okcancelhelp $t.button -bd 2 -relief groove \
	    -ok_command "SeqedSave2 $se $t.range $t.line_len $t.name;\
	    destroy $t"\
	    -cancel_command "destroy $t" \
	    -help_command "show_help spin {SPIN-Sequence-Display-Save}"

    pack $t.range 
    pack $t.line_len -fill x -anchor w
    pack $t.name -fill x
    pack $t.button -side bottom -fill x
}

proc SeqedSave2 {se range line_len filename} {

    #need to do this because entrybox_get checks if the file exists, and if
    #if does, a dialogue box is displayed which does a global grab. For 
    #large amounts of data, everything hangs
    set file [entrybox_get $filename]
    update
    $se save $file [scale_range_from $range] \
	[scale_range_to $range] [entrybox_get $line_len]

}

#
# Turns on and off the translation lines
#

proc seqed_set_translate {se w value args} {
    global $w.Trans0 $w.Trans1 $w.Trans2 $w.Trans3 $w.Trans4 $w.Trans5 $w.Trans6

    if {$value == 0} {
	foreach i $args {
    	    $se translate delete $i
	    set $w.Trans$i 0
	}
    } else {
	foreach i $args {
    	    $se translate add $i
	    set $w.Trans$i 1
	}
    }
}


proc seqed_r_enz_dialogue {w se mode } {
     global tk_utils_defs 
  
    #remove restriction enzyme plot
    if {$mode == 0} {
	$se restriction_enzyme delete
	return
    }
  
#    set l [keylget tk_utils_defs R_ENZ.WIN]
    set l $w.enzymes

    if {[xtoplevel $l] == ""} return
    wm title $l  "Select restriction enzymes"
    global $l.list

######### yc ###########

    set oldFocus [focus]
    set aa [winfo exists $l]
    grab set $l
    focus $l
#    tkwait variable button
#    grab release $l
#    destroy $l
    focus $oldFocus
######## yc ###########

    if {![renzbox $l]} {
	#pressed cancel
	return
    }

    #tkwait variable $l.list

    set list [set $l.list]
    set filename [renzbox_filename $l]
    set num_items [llength $list]

    #if the user has not selected any enzymes
    #if {[llength $list] == 0} {
	#tk_dialog .error "Error" "No selection has been made" error 0 OK
	#raise $l
	#tkwait variable press_ok
    #}

    $se restriction_enzyme add $filename $list $num_items
    raise $w
}


proc seqed_r_enz_dialogueOLD { se mode } {


    #remove restriction enzyme plot
    if {$mode == 0} {
	$se restriction_enzyme delete
	return
    }

    set l ".list_of_restriction_enzymes"
    set cancel_cmd "{}"

    if {[xtoplevel $l] == ""} return
    wm protocol $l WM_DELETE_WINDOW "eval $cancel_cmd; destroy $l"

    #########################################################################
    #OK and Cancel buttons
    okcancelhelp $l.ok_cancel \
	    -ok_command "seqed_renzbox_OK_Pressed $se $l; destroy $l" \
	    -cancel_command "eval $cancel_cmd;  destroy $l" \
	    -help_command "puts help"\
	    -bd 2 \
	    -relief groove

    pack $l.ok_cancel -side bottom -fill x
    renzbox $l

    

}

proc seqed_renzbox_OK_Pressed {se l} {

    set list_name [renzbox_path $l]
    set filename [renzbox_filename $l]

    set list [GetSelection $list_name]
    set num_items [llength $list]

    #if the user has not selected any enzymes
    if {[llength $list] == 0} {
	tk_messageBox \
		-icon error \
		-title "Error" \
		-message "No selection has been made" \
		-type ok \
		-parent $l
	raise $l
	tkwait variable press_ok
    }
    
    $se restriction_enzyme add $filename $list $num_items
}

proc display_seqed_info {w x y} {

    set info [$w seq_info $x $y]

}

#
# Add bindings
#
bind Seqed <Any-Enter>		{focus %W}

bind Seqed <1>	{
    %W cursor_set @%x @%y
    #if {[%W cursor_set @%x @%y] == 0} {
	#%W select from @%x
    #} else {
	#%W select clear
    #}
}
bind Seqed <3>		        {display_seqed_info %W @%x @%y}
bind Seqed <Key-Left>		{%W cursor_left}
bind Seqed <Control-Key-b>	{%W cursor_left}

bind Seqed <Key-Right>		{%W cursor_right}
bind Seqed <Control-Key-f>	{%W cursor_right}

bind Seqed <Key-Up>		{%W cursor_up}
bind Seqed <Control-Key-p>	{%W cursor_up}

bind Seqed <Key-Down>		{%W cursor_down}
bind Seqed <Control-Key-n>	{%W cursor_down}
