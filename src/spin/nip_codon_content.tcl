#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

proc CodonUsage { } {
    global nip_defs

    set t .codon_usage

    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Calculate and write codon table to disk"
    
    ##################################################################
    #reading in existing table either for adding data to the table or
    #concatenating with a new table

    frame $t.s

    ##################################################################
    #use existing codon table 
    keylset ct C_TABLE [keylget nip_defs NIP.CODON_USAGE.C_TABLE]
    frame $t.s.c_table
    getFname $t.s.c_table.name [keylget ct C_TABLE.NAME] load 

    keylset ct C_TABLE_YN [keylget nip_defs NIP.CODON_USAGE.C_TABLE_YN]
    yes_no $t.s.c_table.yn \
	-title   [keylget ct C_TABLE_YN.NAME] \
	-default [keylget ct C_TABLE_YN.VALUE] \
	-ycommand "getFname_configure $t.s.c_table.name -state normal" \
	-ncommand "getFname_configure $t.s.c_table.name -state disabled" \
	-orient horizontal
    pack $t.s.c_table.yn $t.s.c_table.name -side top -fill both


    keylset op OUTPUT [keylget nip_defs NIP.CODON_USAGE.OUTPUT]
    set b1 [keylget op OUTPUT.BUTTON.1]
    set b2 [keylget op OUTPUT.BUTTON.2]
#    radiolist $t.s.output \
#	-title   [keylget op OUTPUT.NAME] \
#	-default [keylget op OUTPUT.VALUE] \
#	-orient vertical \
#	-buttons [format {\
#	{%s -command {getFname_configure %s -state disabled}} \
#	{%s -command {getFname_configure %s -state normal}}} \
#		      [list $b1] [list $t.s.c_table.name]\
#		      [list $b2] [list $t.s.c_table.name]]
    
    radiolist $t.s.output \
	-title   [keylget op OUTPUT.NAME] \
	-default [keylget op OUTPUT.VALUE] \
	-orient vertical \
	-buttons [format {\
	{%s -command {yes_no_configure %s -default 0;\
	getFname_configure %s -state disabled}} \
	{%s -command {yes_no_configure %s -default 1;\
	getFname_configure %s -state normal}}} \
		      [list $b1] [list $t.s.c_table.yn]\
		      [list $t.s.c_table.name]\
		      [list $b2] [list $t.s.c_table.yn]\
		      [list $t.s.c_table.name]]
    #output filename
    keylset fn FILENAME [keylget nip_defs NIP.CODON_USAGE.FILENAME]
    entrybox $t.s.filename \
	-title   [keylget fn FILENAME.NAME] \
	-default [keylget fn FILENAME.VALUE] \
	-type CheckOutputOptional

    pack $t.s.output $t.s.c_table $t.s.filename -fill both -expand 1
    ##################################################################
    #calculate codon total as sum or percentage
    keylset ct CODON_TOTAL [keylget nip_defs NIP.CODON_USAGE.CODON_TOTAL]
    set b1 [keylget ct CODON_TOTAL.BUTTON.1]
    set b2 [keylget ct CODON_TOTAL.BUTTON.2]
    radiolist $t.format \
	-title   [keylget ct CODON_TOTAL.NAME] \
	-default [keylget ct CODON_TOTAL.VALUE] \
	-orient horizontal \
	-buttons [format {{%s} {%s}} \
		      [list $b1] [list $b2]]

    ##################################################################
    #which strand
    strand $t.strand

    ##################################################################
    #range 
    frame $t.range
    keylset ra RANGE [keylget nip_defs NIP.CODON_USAGE.RANGE]
    set seq_id [get_active_seq_id] 
    global $seq_id.start $seq_id.end

    set seq_length [seq_info $seq_id length] 
    set seq_start [seq_info $seq_id start] 
    set seq_end [seq_info $seq_id end] 
    if {[info exists $seq_id.start]} {
	set seq_start [set $seq_id.start]
    }    
    if {[info exists $seq_id.end]} {
	set seq_end [set $seq_id.end]
    }
  
    seq_id $t.range.single -range 1 -browse 1 -from 1 -to $seq_length \
	-start_value $seq_start -end_value $seq_end -min_value 2 \
	-default [seq_info $seq_id name] \
	-update_cmd [list [list seq_range_updates $t.range.single]]\
	-browse_cmd nip_seq_browser

    #keylset ra INFILE [keylget nip_defs NIP.CODON_USAGE.INFILE]
    #lorf_in $t.range.infile [keylget ra INFILE] \
	"{scale_range_configure $t.range.single -state disabled}\
         {scale_range_configure $t.range.single -state disabled}\
         {scale_range_configure $t.range.single -state normal}" 
    #pack $t.range.infile -fill x -expand 1
    pack $t.range.single -fill x -expand 1

    frame $t.separator -bd 2 -relief raised -height 2

    ##################################################################
    #ok cancel help buttons 
    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "CodonUsage2 $t $t.s.output $t.s.filename $t.s.c_table $t.format $t.strand $t.range"\
	-cancel_command "destroy $t" \
	-help_command "show_help spin {SPIN-Codon-Usage-Tables}"
    
    pack $t.range -fill x -expand 1
    pack $t.strand -fill x -expand 1
    pack $t.format -fill x -expand 1
    pack $t.separator -side top -fill x -padx 10 -pady 5
    pack $t.s -fill x -expand 1
    pack $t.button -side bottom -fill x
}


proc CodonUsage2 {t output filename c_table format strand range} {
    global PROTEIN

    set concat 0
    set out [radiolist_get $output]
    if {$out == 2} {
	set concat 1
    }

    set ranges ""

    if {[yes_no_get $c_table.yn] == 1} {

	set c_table [getFname_in_name $c_table.name]
	if {$c_table == ""} {
	    verror ERR_WARN "codon usage" "no filename entered"
	}
    } else {
	set c_table ""
    }

    set seq_id [name_to_seq_id [seq_id_name $range.single]]
    if {[seq_info $seq_id type] == $PROTEIN} {
	verror ERR_WARN "codon content" "unable to process protein sequences"
	seq_id_destroy $range.single
	return
    }

    #set range_type [lorf_in_get $range.infile]
    set range_type 3
    if {$range_type == 3} {
	#single 
	set ranges [list "[seq_id_from $range.single] [seq_id_to $range.single]"]
    } else {
	#list or file
	set ranges [lorf_get_list $range.infile 0]
    } 

    if {$ranges == ""} {
	raise $t
	return
    }

    set fn [entrybox_get $filename]

    SetBusy

    codon_usage -seq_id $seq_id -outfile $fn \
	    -concatenate $concat \
	    -format [radiolist_get $format] \
	    -strand [strand_get $strand] -table $c_table -range $ranges

    ClearBusy

    seq_id_destroy $range.single
    global $seq_id.start $seq_id.end
    set $seq_id.start [seq_id_from $range.single]
    set $seq_id.end [seq_id_to $range.single]
    destroy $t
}
