proc SetAlignmentSymbols { } {
        
        global cutoff_align1
        global cutoff_align2
        global cutoff_align3
        global symbol_align0
        global symbol_align1
        global symbol_align2
        global symbol_align3
       
	global sip_defs PROTEIN DNA 
        set s .setalignmentsymbols
        if {[xtoplevel $s -resizable 0] == ""} return

        wm title $s "Set Protein Alignment Symbols"

        foreach col { 1 2  } {
            frame $s.col$col
	    pack $s.col$col -side top -fill both  
       }

	label $s.col1.1 -text "Identity"
        entrybox $s.col2.1 -width 4  -default $symbol_align0

	label $s.col1.2 -text "\>"
        label $s.col2.2 -text "    "

        entrybox $s.col1.3 -width 4  -default $cutoff_align1 -type "CheckIntRange -8 17"
        entrybox $s.col2.3 -width 4  -default $symbol_align1 
        
        label $s.col1.4 -text "\>"
        label $s.col2.4 -text "  "

        entrybox $s.col1.5 -width 4  -default $cutoff_align2  -type "CheckIntRange -8 17"
        entrybox $s.col2.5 -width 4  -default $symbol_align2 

        label $s.col1.6 -text "\>"
        label $s.col2.6 -text "  "

        entrybox $s.col1.7 -width 4  -default $cutoff_align3 -type "CheckIntRange -8 17"
        entrybox $s.col2.7 -width 4  -default $symbol_align3 
       
        label $s.col1.8 -text "  Others  "
	label $s.col2.8 -text "  Space  "

        pack $s.col1.1 $s.col1.2 $s.col1.3 $s.col1.4  $s.col1.5 \
             $s.col1.6 $s.col1.7 $s.col1.8 -side left -fill x
        pack $s.col2.1 $s.col2.2 $s.col2.3 $s.col2.4 $s.col2.5 \
             $s.col2.6 $s.col2.7 $s.col2.8 -side left -fill x -anchor center

        #ok cancel help buttons 
        okcancelhelp $s.button -bd 2 -relief groove \
	    -ok_command "ChangeAlignmentSymbols $s $s.col1.3 $s.col2.3 $s.col1.5 $s.col2.5 $s.col1.7 $s.col2.7" \
	    -cancel_command " destroy $s" \
             -help_command "show_help spin {SPIN-Set protein alignment symbols}"
        pack $s.button -fill x -side bottom
#        tkwait 
    }
proc ChangeAlignmentSymbols {t col1.3 col2.3 col1.5 col2.5 col1.7 col2.7 } {     
     global symbol_align0
     global cutoff_align1 
     global symbol_align1
     global cutoff_align2 
     global symbol_align2
     global cutoff_align3 
     global symbol_align3
     
     
     set symbol_align0 [entrybox_get $t.col2.1]
     set cutoff_align1 [entrybox_get $t.col1.3]
     set symbol_align1 [entrybox_get $t.col2.3]
     set cutoff_align2 [entrybox_get $t.col1.5]
     set symbol_align2 [entrybox_get $t.col2.5]
     set cutoff_align3 [entrybox_get $t.col1.7]
     set symbol_align3 [entrybox_get $t.col2.7]
    
    if { $cutoff_align1 <= $cutoff_align2 || $cutoff_align1 <= $cutoff_align3 || $cutoff_align2 <= $cutoff_align3} {
#	raise [winfo toplevel $t.col1.3]
#	focus $t.col1.3
        bell
        return 0
    }
        destroy $t
}









