#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc ExtractReadingsDialog { io f } {
    global gap_defs

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Extract readings"

    contig_id $f.id -title "Reading identifier" -io $io -range 0

    ###########################################################################
    #file of reading names
    lorf_in $f.infile [keylget gap_defs EXTRACT.INFILE] \
	"{contig_id_configure $f.id -state disabled} \
	 {contig_id_configure $f.id -state disabled}\
	 {contig_id_configure $f.id -state normal}" -bd 2 -relief groove

    set ex [keylget gap_defs EXTRACT.DEST]
    entrybox $f.dest_dir \
	-title [keylget ex NAME ] \
	-default [keylget ex VALUE] \
	-width 15 \
	-type CheckString

    set fm [keylget gap_defs EXTRACT.FORMAT]
    radiolist $f.format \
	-title [keylget fm NAME] \
	-bd 2 -relief groove \
	-default [keylget fm VALUE] \
	-buttons [format {{%s} {%s} {%s}} \
	[list [keylget fm B1]] \
	[list [keylget fm B2]] \
	[list [keylget fm B3]]]

   ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "Extract_OkPressed $io $f $f.infile $f.id $f.dest_dir \
	    $f.format"\
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Extract Readings}" \
	    -bd 2 \
	    -relief groove

    pack $f.infile -fill x
    pack $f.id -fill x
    pack $f.dest_dir -fill x
    pack $f.format -fill x
    pack $f.ok_cancel -fill both

}

proc Extract_OkPressed { io f infile id dest_dir format} {
    global FBRec
    global NGList
    global gap_defs

    #special case for a single reading
    if {[lorf_in_get $infile] == 3} {
	set list [contig_id_gel $id]
    } else {
	set list [lorf_get_list $infile]
    }
    set dir [entrybox_get $dest_dir]

    if {$list == ""} {
	tk_messageBox -icon error -type ok -title "Error" \
		-message "No reading names specified" \
		-parent $f
	return
    }

    set format [lindex {0 0 2 3} [radiolist_get $format]]
    SetBusy
    extract_readings \
	-io $io \
	-readings $list \
	-directory $dir \
	-format $format
    ClearBusy

    destroy $f
}
