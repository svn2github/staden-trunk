#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc ComplementContig {io} {
    global gap5_defs

    set l [keylget gap5_defs COMPLEMENT_CONTIG]
    set t [keylget l WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Complement contig / scaffold"

    contig_id $t.id \
	-command "ComplementContig2 $io $t $t.id" \
	-io $io \
	-range 0 \
	-scaffold 1

    okcancelhelp $t.ok \
	-ok_command "ComplementContig2 $io $t $t.id" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Complement}" \
	-bd 2 -relief groove

    pack $t.id $t.ok -side top -fill x
}

proc ComplementContig2 {io t id} {
    if {[set c [contig_id_gel $id]] == ""} {bell; return}

    destroy $t
    update idletasks

    if {[db_info get_scaffold_num $io $c] > 0} {
	log_call complement_scaffold -io $io -scaffolds "$c"
    } else {
	set cnum [db_info get_contig_num $io $c]
	log_call complement_contig -io $io -contigs "$c"
	SetContigGlobals $io [left_gel $io $cnum]
    }
}


proc ContigRenameBulk {io} {
    global gap5_defs

    set w .rename_contigs
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Bulk contig rename"

    frame $w.f -bd 2 -relief groove -highlightthickness 2
    xentry $w.f.pattern \
	-label "Match pattern" \
	-default [keylget gap5_defs CONTIG_BULK_RENAME.SEARCH]

    lorf_in $w.f.infile [keylget gap5_defs CONTIG_BULK_RENAME.INFILE] \
        "{$w.f.pattern configure -state disabled}
         {$w.f.pattern configure -state disabled}
         {$w.f.pattern configure -state normal}" \
	-bd 0

    pack $w.f.infile $w.f.pattern -fill both -expand 1

    frame $w.g -bd 2 -relief groove -highlightthickness 2
    xentry $w.g.replace \
	-label "Replace pattern" \
	-default [keylget gap5_defs CONTIG_BULK_RENAME.REPLACE]

    xentry $w.g.start \
	-label "Auto-increment starting value" \
	-default [keylget gap5_defs CONTIG_BULK_RENAME.INDEX] \
	-type int

    pack $w.g.replace $w.g.start -fill both -expand 1


    okcancelhelp $w.ok \
	-ok_command "ContigRenameBulk2 $io $w" \
	-cancel_command "destroy $w" \
	-help_command "show_help gap5 {Contig Bulk Rename}" \
	-bd 2 -relief groove -highlightthickness 2

    pack $w.f $w.g $w.ok -fill both -expand 1
}

proc ContigRenameBulk2 {io w} {
    if {[lorf_in_get $w.f.infile] != 3} {
	foreach n [lorf_get_list $w.f.infile] {
	    set c_arr($n) 1
	}
	set pattern "*"
    } else {
	set pattern [$w.f.pattern get]
    }

    set replace [$w.g.replace get]
    set start   [$w.g.start get]

    if {$pattern == "" || $replace == "" || $start == ""} {
	bell
	return
    }
    if {[regexp {\s+} $replace]} {
	tk_messageBox -icon warning -type ok -parent $w \
	    -title "Bulk Rename Contig" \
	    -message "Sorry, the replacement pattern may not contain spaces"
	return
    }

    # Convert pattern from C-shell style filename glob to regexp
    regsub -all {\.} $pattern {\\.} pattern
    regsub -all {\?} $pattern {(.)}   pattern
    regsub -all {\*} $pattern {(.*)}  pattern
    set pattern "^$pattern\$"

    # Convert ? and * in replace pattern with numeric expansions
    regsub -all {\?} $replace {%s} replace
    regsub -all {\*} $replace {%s} replace

    vfuncheader "Bulk rename contig"
    
    # Iterate through contigs checking them.
    set db [$io get_database]
    set nc [$db get_num_contigs]
    for {set i 0} {$i < $nc} {incr i} {
        set cnum [$io contig_order $i]
        set c [$io get_contig $cnum]
	set name [$c get_name]
	$c delete

	if {[info exists c_arr]} {
	    if {![info exists c_arr($name)]} continue
	}

	if {[regsub $pattern $name [format $replace $start] name2]} {
	    # Matched, so rename it
	    if {$name == $name2} {
		vmessage "Skipping renaming of contig $name to itself"
		incr start
		continue
	    }
	    set name2 [log_call contig_rename $io $cnum $name2 {} 1]
	    vmessage -nonewline "Renaming contig $name to "
	    vmessage_tagged "$name2" SEQID
	    incr start
	}
    }

    $io flush

    destroy $w
}
