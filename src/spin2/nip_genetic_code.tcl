# -----------------------------------------------------------------------------
# Selection of genetic code

proc SetGeneticCode {} {
    global nip_defs
    set w .genetic_code

    if {[xtoplevel $w] == ""} return
    wm title $w "Set genetic code"

    frame $w.codes
    pack $w.codes -side top -fill both -expand 1

    listbox $w.codes.l \
	    -xscrollcommand "$w.codes.xs set" \
	    -yscrollcommand "$w.codes.ys set" \
	    -width 25 \
	    -height 15
    scrollbar $w.codes.xs \
	    -orient horizontal \
	    -command "$w.codes.l xview"
    scrollbar $w.codes.ys \
	    -orient vertical \
	    -command "$w.codes.l yview"
    grid rowconfigure $w.codes 0 -weight 1
    grid columnconfigure $w.codes 0 -weight 1
    grid $w.codes.l $w.codes.ys -sticky nsew
    grid $w.codes.xs -sticky nsew

    set fd [open [keylget nip_defs GENETIC_CODE_DIR]/code_index r]
    set i 0
    set index ""
    while {[gets $fd line] != -1} {
	regexp {([^ ]*) *(.*)} $line tmp file name
	$w.codes.l insert end "$name"
	lappend index $file
	incr i
    }
    close $fd

    okcancelhelp $w.ok -bd 2 -relief groove \
	-ok_command "SetGeneticCode2 0 $w $w.codes.l [list $index]" \
	-perm_command "SetGeneticCode2 1 $w $w.codes.l [list $index]" \
	-cancel_command "destroy $w" \
	-help_command "show_help spin {SPIN-Set-Genetic-Code}"
    pack $w.ok -side top -fill both
}

proc SetGeneticCode2 {perm w listbox index} {
    global nip_defs env

    if {[set cur [$listbox curselection]] == ""} {
	bell
	return
    }
    set file [lindex $index [lindex $cur 0]]

    if {[load_genetic_code \
	    -filename [keylget nip_defs GENETIC_CODE_DIR]/$file] != -1} {
	vfuncheader "Set genetic code"
	vmessage "[$listbox get [lindex $cur 0]] genetic code:\n"
	set fd [open [keylget nip_defs GENETIC_CODE_DIR]/$file r]
	vmessage [read $fd]
	close $fd
    }

    keylset nip_defs GENETIC_CODE [keylget nip_defs GENETIC_CODE_DIR]/$file

    if {$perm} {
        update_defs nip_defs $env(HOME)/.niprc GENETIC_CODE
    }

    destroy $w
}

proc SetGeneticCodeOLD { } {
    global nip_defs

    set t .genetic_code

    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Set genetic code"

    ##################################################################
    #enter name of codon table
    keylset ct C_TABLE [keylget nip_defs NIP.GENETIC_CODE.C_TABLE]
    getFname $t.c_table [keylget ct C_TABLE.NAME] load 
   
    ##################################################################
    #select genetic code
    keylset ct TYPE [keylget nip_defs NIP.GENETIC_CODE.TYPE]
    set b1 [keylget ct TYPE.BUTTON.1]
    set b2 [keylget ct TYPE.BUTTON.2]
    set b3 [keylget ct TYPE.BUTTON.3]
    set b4 [keylget ct TYPE.BUTTON.4]
    radiolist $t.format \
	-title   [keylget ct TYPE.NAME] \
	-default [keylget ct TYPE.VALUE] \
	-orient vertical \
	-bd 2 -relief groove \
	-buttons [format {\
		 {%s -command {getFname_configure %s -state disabled} }\
		 {%s -command {getFname_configure %s -state disabled} }\
		 {%s -command {getFname_configure %s -state disabled} }\
		 {%s -command {getFname_configure %s -state normal} } } \
		      [list $b1] [list $t.c_table]\
		      [list $b2] [list $t.c_table]\
		      [list $b3] [list $t.c_table]\
		      [list $b4] [list $t.c_table]]

    ##################################################################
    #ok cancel help buttons 
    okcancelhelp $t.button -bd 2 -relief groove \
	-ok_command "SetGeneticCode2 $t $t.format $t.c_table"\
	-cancel_command "destroy $t" \
	-help_command "show_help spin {SPIN-Set-Genetic-Code}"

    pack $t.format -fill x
    pack $t.c_table -fill x
    pack $t.button -fill x
}


proc SetGeneticCode2OLD {t format c_table} {

    set type [radiolist_get $format]

    if {$type == 4} {
	set codon_table [getFname_in_name $c_table]
	if {$codon_table == ""} {
	verror ERR_WARN "set genetic code" "no filename entered"
	}
	SetBusy
	
	set_genetic_code -code_type $type -table $codon_table
	
	ClearBusy
	
    } else {
	SetBusy
	
	set_genetic_code -code_type $type
	
	ClearBusy
    }
    destroy $t
}