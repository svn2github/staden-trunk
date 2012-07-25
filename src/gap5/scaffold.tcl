#-----------------------------------------------------------------------------
# AGP support

proc ImportAGP {io} {
    set f .import_scaffold

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Import Scaffold from AGP"

    #--- Input filename
    getFname $f.infile "AGP file" load {} "scaffolds.agp"

    #--- OK/cancel/help
    okcancelhelp $f.ok \
	-ok_command "ImportAGP2 $io $f" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 ImportAGP" \
	-bd 2 \
	-relief groove

    pack $f.infile $f.ok -side top -fill both -expand 1
}

;proc ImportAGP2 {io f} {
    set fn [getFname_in_name $f.infile]
    if {$fn == ""} {
	raise $f
	return
    }

    set r [scaffold_from_agp -io $io -filename $fn]
    if {$r != 0} {
	tk_messageBox \
	    -icon error \
	    -title "Import AGP" \
	    -message "Failed to import AGP file" \
	    -type ok
	return
    }

    scaffold_from_agp -io $io -filename $fn
    $io flush

    destroy $f
}