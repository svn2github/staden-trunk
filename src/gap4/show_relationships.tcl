#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
##############################################################################
#gap show_relationships option
proc ShowRelationshipsDialog { io } {
    global gap_defs

    set f [keylget gap_defs SHOWREL.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Show relationships"

    contig_id $f.id \
	    -range 1 \
	    -io $io

    keylset p POS [keylget gap_defs SHOWREL.POS]
    yes_no $f.yn \
	    -title [keylget p POS.NAME] \
	    -default [keylget p POS.VALUE] \
	    -bd 2 \
	    -relief groove \
	    -state normal \
	    -orient horizontal

     lorf_in $f.infile [keylget gap_defs SHOWREL.INFILE] \
	    "{contig_id_configure $f.id -state disabled; \
	      yes_no_configure $f.yn -state normal} \
	    {contig_id_configure $f.id -state disabled; \
	      yes_no_configure $f.yn -state normal}\
	    {contig_id_configure $f.id -state disabled ; \
	      yes_no_configure $f.yn -state normal}\
	    {contig_id_configure $f.id -state normal; \
              yes_no_configure $f.yn -state disabled}" -bd 2 -relief groove

    ###########################################################################
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "SR_OK_Pressed $f $io $f.id $f.infile $f.yn" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Show Relationships}" \
	    -bd 2 \
	    -relief groove
    pack $f.infile $f.id $f.yn $f.ok_cancel -side top -fill both
} 
#end ShowRelationships 

proc SR_OK_Pressed { f io id infile yn} {

    set ordered [yes_no_get $yn]

    set contig_list ""
    if {[lorf_in_get $infile] == 4} {
	#single contig
	set gel_name [contig_id_gel $id]
	set lreg [contig_id_lreg $id]
	set rreg [contig_id_rreg $id]
	
	SetContigGlobals $io $gel_name $lreg $rreg
	set contig_list "{$gel_name $lreg $rreg}"
    } elseif {[lorf_in_get $infile] == 3 } {
	#all contigs
	set contig_list [CreateAllContigList $io]
    } else {
	#list or file
	set contig_list [lorf_get_list $infile]
    }

    if {$contig_list == ""} {
	raise $f
	return
    }

    show_relationships \
	-io $io \
	-contigs $contig_list \
	-order $ordered

    destroy $f
}


