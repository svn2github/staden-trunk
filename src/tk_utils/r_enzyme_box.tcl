#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

##############################################################################
proc ReadREnzFile { input file list l } {
    global tk_utils_defs $l.filename

    set ok_pressed 0
    keylset f E_FILE [keylget tk_utils_defs R_ENZ.E_FILE]
    set filetype [radiolist_get $input]
    
    #read standard filenames defined in user-defaults
    set $l.filename [keylget f E_FILE.$filetype]

    #empty the list box
    $list delete 0 end

    #user has entered own filename
    if {[set $l.filename] == ""} {
        entrybox_delete $l.infile.entry 0 end
	return
    }
    #fill in list box with contents of the file
    set renz_list [read_enz_file -file [set $l.filename]]
    eval $list insert end $renz_list
}

##############################################################################
#execute this bit of code if user has entered personal r enzyme file
proc ReadUserFileCheck {entry list l} {
    set e [entrybox_path $entry]
    if {[$e cget -state] != "normal" || [$e get] == ""} {
	return 0
    }

    ReadUserFile $entry $list $l [$e get]
}

proc ReadUserFile {entry list l file_name} {
    global $l.filename

    set e [entrybox_path $entry]

    if {[$e get] == "" || ![file exists [$e get]]} {
	bell
	return 0
    }

    set $l.filename $file_name

    #fill in list box with contents of the file
    #read_enz_file -file [set $l.filename] -window $list

    set renz_list [read_enz_file -file [set $l.filename]]
    if {[string compare $renz_list [$list get 0 end]]} {
	#empty the list box
	$list delete 0 end
	
	eval $list insert end $renz_list
    }

    return 1
}

##############################################################################
proc GetSelection {list} {
    return [$list curselection]
}

#############################################################################
#create a list box for restriction enzymes
proc renzBox { l contig io seq_id scale from to start_value end_value} {
    global tk_utils_defs $l.retval

    #if [winfo exists $l] {raise $l; return}
    #toplevel $l 
    wm protocol $l WM_DELETE_WINDOW "destroy $l"
    wm resizable $l 0 0

    #########################################################################
    #create listbox
    frame $l.t -borderwidth 0
    listbox $l.t.lists -selectmode extended -yscrollcommand "$l.t.scrolly set"
    #canvas $l.col -yscrollcommand "$l.scrolly set" -width 50 -height 218 
    bind $l.t.lists <1> {+focus %W}
    bind $l.t.lists <Control-a> {%W selection set 0 end}
    focus $l.t.lists

    scrollbar $l.t.scrolly -command "$l.t.lists yview" -orient vertical

    #########################################################################
    #file name selection for entering personal file
    getFname $l.infile [keylget tk_utils_defs R_ENZ.INFILE.NAME] load

    #need to redefine the browse button command to update the list box with
    #the renz names in the personal file
    $l.infile.browse configure -command "InvokeFileBrowser $l.infile.entry open;\
	    ReadUserFileCheck $l.infile.entry $l.t.lists $l"

    
    #add binding for entrybox
    entrybox_configure $l.infile.entry \
	-command "ReadUserFile $l.infile.entry $l.t.lists $l"
    bind [entrybox_path $l.infile.entry] <Any-Leave> \
	"ReadUserFileCheck $l.infile.entry $l.t.lists $l"

    #remove toplevel Return key binding which invokes the ok button
    bindtags [entrybox_path $l.infile.entry] "[entrybox_path $l.infile.entry] Entry all"

    #########################################################################
    #keyboard entry of enzyme name and target sequence
    keylset n E_NAME [keylget tk_utils_defs R_ENZ.E_NAME]
    entrybox $l.r_name \
	-title   "[keylget n E_NAME.NAME]" \
	-width 20 \
	-type "CheckString"

    entrybox $l.r_seq \
	-title   "[keylget n E_NAME.SEQ]" \
	-width 20 \
	-type "CheckString"

    #########################################################################
    #radiobuttons to select file of enzyme names
    keylset st SELFILE [keylget tk_utils_defs R_ENZ.SELFILE]

    set b1 [keylget st SELFILE.BUTTON.1]
    set b2 [keylget st SELFILE.BUTTON.2]
    set b3 [keylget st SELFILE.BUTTON.3]
    set b4 [keylget st SELFILE.BUTTON.4]

    frame $l.but
    radiolist $l.sel_file \
	    -title [keylget st SELFILE.NAME]\
	    -bd 2 \
	    -relief groove \
	    -orient vertical \
	    -default [keylget st SELFILE.VALUE] \
	    -buttons [format { \
	    {%s -command {getFname_configure %s -state disabled; \
                 ReadREnzFile %s %s %s %s} } \
	    {%s -command {getFname_configure %s -state disabled; \
                 ReadREnzFile %s %s %s %s} } \
	    {%s -command {getFname_configure %s -state disabled; \
                 ReadREnzFile %s %s %s %s} } \
	    {%s -command {getFname_configure %s -state normal;\
	         ReadREnzFile %s %s %s %s} } } \
	    [list $b1] [list $l.infile] [list $l.sel_file] [list $l.infile]\
	    [list $l.t.lists] [list $l]\
	    [list $b2] [list $l.infile] [list $l.sel_file] [list $l.infile]\
	    [list $l.t.lists] [list $l]\
	    [list $b3] [list $l.infile] [list $l.sel_file] [list $l.infile]\
	    [list $l.t.lists] [list $l]\
	    [list $b4] [list $l.infile] [list $l.sel_file] [list $l.infile]\
	    [list $l.t.lists] [list $l] ]

    #########################################################################
    #if contig_id required
    if {$contig} {
	contig_id $l.sc -io $io 
    }

    if {$seq_id > -1} {
	#FIXME - if I ever want to do restriction enzyme search in anything
	#other than nip
	seq_id $l.range -range 1 -browse 1 -from $from -to $to \
	    -start_value $start_value -end_value $end_value -min_value 1 \
	    -default [seq_info $seq_id name] \
	    -update_cmd [list [list seq_range_updates $l.range]] \
	    -browse_cmd nip_seq_browser
    }

    #########################################################################
    #if scale range is required
    keylset us RANGE [keylget tk_utils_defs R_ENZ.RANGE]

    if {$scale} {
	scale_range $l.range -title [keylget us RANGE.NAME] \
	    -from $from -to $to -start_value $start_value \
	    -end_value $end_value -bd 2 -relief groove\
	    -min_value 1
    }

    if {$contig} {
	set help "show_help gap4 {Restrict-Selecting}"
    } else {
	set help "show_help spin {SPIN-Restrict-Selecting}"
    }

    #########################################################################
    #OK and Cancel buttons
    okcancelhelp $l.ok_cancel \
	    -ok_command "renzbox_OK_Pressed $l $seq_id" \
	    -cancel_command "renzbox_Cancel_Pressed $l $seq_id" \
	    -help_command $help \
	    -bd 2 \
	    -relief groove

    #final packing
    #pack $l.ok_cancel -side bottom -fill x
    #pack $l.infile -side bottom -fill x
    #pack $l.sel_file -side bottom -fill x
    #pack $l.scrolly -side right -fill y
    #pack $l.lists -fill both

    pack $l.t -fill both 
    pack $l.t.scrolly -side right -fill y
    pack $l.t.lists -fill both
    pack $l.sel_file -fill x
    pack $l.infile -fill x
    if {$contig} {
	pack $l.sc -fill x
    }
    if {$seq_id > -1} {
	pack $l.range -fill x
    }
    if {$scale} {
	pack $l.range -fill x
    }
    pack $l.ok_cancel -fill x

    #number of items in the list box
    set num_items [$l.t.lists size]

    tkwait variable $l.retval
    destroy $l
    return [set $l.retval]
}

proc renzbox_OK_Pressed { l seq_id} {
    global $l.list $l.retval

    set list_name [renzbox_path $l]
    set filename [renzbox_filename $l]

    set $l.list [GetSelection $list_name]

    if {[winfo exists $l.sc]} {
	global $l.contig_name $l.from $l.to
	set $l.contig_name [contig_id_gel $l.sc]
	set $l.from [contig_id_lreg $l.sc]
	set $l.to [contig_id_rreg $l.sc]
    }

    if {$seq_id > -1} {
	global $l.seqid
	set $l.seqid "{[seq_id_name $l.range]} [seq_id_from $l.range] [seq_id_to $l.range]"
	seq_id_destroy $l.range
    } elseif {[winfo exists $l.range]} {
	global $l.from $l.to
	set $l.from [scale_range_from $l.range]
	set $l.to [scale_range_to $l.range]
    }

    set num_items [llength [set $l.list]]

    #if the user has not selected any enzymes
    if {[llength [set $l.list]] == 0} {
	tk_messageBox -title "Error" -message "No selection has been made" -icon error -type ok
	raise $l
	return
    }

    set $l.retval 1
}

proc renzbox_Cancel_Pressed { l seq_id} {
    global $l.retval
    set $l.retval 0

    if {$seq_id > -1} {
	seq_id_destroy $l.range
    } 
}

proc renzbox_path { l } {

    return $l.t.lists
}

proc renzbox_filename { l } {
    global $l.filename

    return [set $l.filename]
}

proc renzbox_contig_name { l } {
    global $l.contig_name

    return [set $l.contig_name]
}

proc renzbox_from {l} {
    global $l.from

    return [set $l.from]
}

proc renzbox_to {l} {
    global $l.to

    return [set $l.to]
}

proc renzbox_get_seqid {l} {
    global $l.seqid

    return [set $l.seqid]
}

#restriction enzyme selection box only
proc renzbox { l } {

    return [renzBox $l 0 0 -1 0 0 0 0 0]

}

#restriction enzyme selection box plus contig_id selection
proc contig_renzbox {l io } {
    return [renzBox $l 1 $io -1 0 0 0 0 0]
}

#restriction enzyme selection box with ranges only
proc scale_renzbox { l from to start_value end_value} {
    return [renzBox $l 0 0 -1 1 $from $to $start_value $end_value]
}

#restriction enzyme selection box plus seq_id selection
proc seq_renzbox { l seq_id from to start_value end_value} {
    return [renzBox $l 0 0 $seq_id 0 $from $to $start_value $end_value]
}
