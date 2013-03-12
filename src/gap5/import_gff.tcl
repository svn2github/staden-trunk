#-----------------------------------------------------------------------------
# Export Tags in GFF format

proc ImportGFF {io} {
    global gap5_defs

    set f .import_gff
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Import GFF Annotations"

    #--- Input filename
    getFname $f.infile "GFF filename" load {} "tags.gff"

    #--- Padding options
    checkbutton $f.padded \
	-text "GFF coordinates are already padded" \
	-variable $f.Padded \
	-anchor w
    set $f.Padded 0

    #--- OK/cancel/help
    okcancelhelp $f.ok \
	-ok_command "ImportGFF2 $io $f" \
	-cancel_command "destroy $f" \
	-help_command "show_help gap5 ImportGFF" \
	-bd 2 \
	-relief groove

    #--- Packing
    pack $f.infile $f.padded $f.ok -side top -fill both
}

proc ImportGFF2 {io f} {
    global $f.Padded

    set fn [getFname_in_name $f.infile]
    if {$fn == ""} {
	raise $f
	return
    }

    log_call import_gff -io $io -infile $fn -padded [set $f.Padded]
    destroy $f
}
