#-----------------------------------------------------------------------------
# AGP support

proc ImportAGP {io} {
    set f .import_scaffold

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Import Scaffold from AGP"

    #--- Input filename
    getFname $f.infile "AGP file" load {} scaffolds.agp

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

#-----------------------------------------------------------------------------
proc ExportAGP {io} {
    set f .export_scaffold

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Export Scaffold to AGP"

    #--- Output filename
    getFname $f.outfile "AGP file" save {} scaffolds.agp

    #--- OK/cancel/help
    okcancelhelp $f.ok \
	-ok_command "ExportAGP2 $io $f" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 ExportAGP" \
	-bd 2 \
	-relief groove

    pack $f.outfile $f.ok -side top -fill both -expand 1
}

;proc ExportAGP2 {io f} {
    set fn [getFname_in_name $f.outfile]
    if {$fn == ""} {
	raise $f
	return
    }

    set r [scaffold_to_agp -io $io -filename $fn]
    if {$r != 0} {
	tk_messageBox \
	    -icon error \
	    -title "Export AGP" \
	    -message "Failed to export AGP file" \
	    -type ok
	return
    }

    scaffold_to_agp -io $io -filename $fn
    $io flush

    destroy $f
}