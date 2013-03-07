proc ContigExtend {io} {
    global gap5_defs

    set f .contig_extend
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Contig Extend"

    #--- contig identifier widget
    contig_id $f.id -io $io

    lorf_in $f.infile [keylget gap5_defs EXTEND_CONTIGS.INFILE] \
        "{contig_id_configure $f.id -state disabled}
         {contig_id_configure $f.id -state disabled}
         {contig_id_configure $f.id -state disabled}
         {contig_id_configure $f.id -state normal}
        " -bd 2 -relief groove

    #--- Scoring extension parameters
    xentry $f.min_depth \
	-label "Minimum depth" \
	-default [keylget gap5_defs EXTEND_CONTIGS.MIN_DEPTH] \
	-type int

    xentry $f.match_score \
	-label "Score for a match" \
	-default [keylget gap5_defs EXTEND_CONTIGS.MATCH_SCORE] \
	-type int

    xentry $f.mismatch_score \
	-label "Score for a mismatch" \
	-default [keylget gap5_defs EXTEND_CONTIGS.MISMATCH_SCORE] \
	-type int

    #--- OK/cancel/help
    okcancelhelp $f.ok \
	-ok_command "ContigExtend2 $io $f" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 ContigExtend" \
	-bd 2 \
	-relief groove

    #--- Packing
    pack $f.infile $f.id $f.min_depth $f.match_score $f.mismatch_score $f.ok \
	-side top -fill both
}

proc ContigExtend2 {io f} {
    if {[lorf_in_get $f.infile] == 4} {
	set gel_name [contig_id_gel $f.id]
	set lreg [contig_id_lreg $f.id]
	set rreg [contig_id_rreg $f.id]
	set list [list [list $gel_name $lreg $rreg]]
    } elseif {[lorf_in_get $f.infile] == 3} {
	set list [CreateAllContigList=Numbers $io]
    } else {
	set list [lorf_get_list $f.infile]
    }

    set min_depth [$f.min_depth get]
    set match_score [$f.match_score get]
    set mismatch_score [$f.mismatch_score get]

    log_call contig_extend \
	-io $io \
	-contigs $list \
	-min_depth $min_depth \
	-match_score $match_score \
	-mismatch_score $mismatch_score

    destroy $f
}
