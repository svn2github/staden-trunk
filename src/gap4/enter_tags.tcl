#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc EnterTags {io} {
    global gap_defs

    set t [keylget gap_defs ENTER_TAGS.WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Enter tags"

    getFname $t.file [keylget gap_defs ENTER_TAGS.FILENAME] load

    xyn $t.unpadded \
	-label "Unpadded tag positions?" \
	-orient horiz \
	-default 0

    okcancelhelp $t.ok \
	-ok_command "EnterTags2 $io $t $t.file $t.unpadded" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Enter Tags}" \
	-bd 2 -relief groove

    pack $t.file $t.unpadded $t.ok -side top -fill both
}

proc EnterTags2 {io t file unpadded} {
    if {[set name [entrybox_get $t.file.entry]] == ""} {bell; return}
    set unpadded [$unpadded get]
    destroy $t

    if {![quit_displays $io "enter_tags"]} {
	# Someone's too busy to shutdown?
	return
    }

    SetBusy
    set tags [reformat_tag_file $io $name]
    add_tags -io $io -tags $tags -unpadded $unpadded
    ContigInitReg $io
    ClearBusy

#    SetBusy
#    enter_tags -io $io -file $name
#
#    ContigInitReg $io
#    ClearBusy
}






# Converts an experiment file format file of tags into a Tcl list ready to
# be passed to the add_tags function.
proc reformat_tag_file {io file} {
    if {[catch {set fd [open $file r]}]} {
	return ""
    }
    set rnum() ""
    set cnum() ""
    while {[gets $fd line] != -1} {
	regexp {^([^ ][^ ])   (     )?(.*)$} $line dummy id cont data

	if {![info exists id]} {
	    continue
	}
	switch $id {
	    ID {
		set index $data
	    }
	    TC -
	    TG {
		set type [expr {$id == "TC" ? "contig" : "seq"}]
		if {$type == "contig"} {
		    if {![info exists cnum($index)]} {
			set cnum($index) [db_info get_contig_num $io $index]
		    }
		    set inum -$cnum($index)
		} else {
		    if {![info exists rnum($index)]} {
			set rnum($index) [db_info get_read_num $io $index]
		    }
		    set inum $rnum($index)
		}
		if {$cont == ""} {
		    if {![info exists arr($index.$type)]} {
			set arr($index.$type) ""
		    }
		    lappend arr($index.$type) "$inum $data"
		} else {
		    set last [lindex $arr($index.$type) end]
		    append last "\n$data"
		    set arr($index.$type) \
			[lreplace $arr($index.$type) end end $last]
		}
	    }
	}
    }
    close $fd

    set result ""
    foreach i [array names arr] {
	append result "$arr($i)\n"
    }
    return $result
}
