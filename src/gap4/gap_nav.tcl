# Copyright Genome Research Limited (GRL). All rights reserved

proc ViewNavigationData {io} {
    
    set w .nav
    global gap_defs $w.Count

    if {[xtoplevel $w] == ""} return
    wm title $w "Contig Navigation File"

    if {![info exists $w.Count]} {
        set $w.Count 0
    }

    xget_fname $w.navfile \
        -text "Navigation file" \
        -type text \

    pack $w.navfile -side top -fill both

    okcancelhelp $w.ok \
        -ok_command "OpenNavFile $io $w" \
        -cancel_command "destroy $w" \
        -help_command "show_help gap4 {Contig Navigation}" \
        -bd 2 -relief groove
    
    pack $w.ok -side top -fill both
}


proc OpenNavFile { io w } {

    global $w.Count
    set t .navtable[set $w.Count]
    incr $w.Count
    global $t.File $t.Auto $t.Edwin $t.FileName

    set fname [$w.navfile get2]
    if {$fname == ""} {
        bell
        return
    }
    if {[catch {set fd [open $fname r]}]} {
        bell
        tk_messageBox -message "Couldn't open $fname" \
            -type ok \
            -parent $w
        return
    }
    
    set $t.File [split [read $fd] "\n"]
    set $t.FileName $fname
    
    close $fd

    CreateTableList $io $t
    destroy $w

    vfuncheader "Contig Navigation"
}

# Call this from FocusOut events.
# When window 'w' loses focus it checks if the mouse window is still in the
# window. If it is then it means the window manager stole our focus by moving
# it another window (probably a new one - gnome wm's love to grant focus to
# any newly created window), and in this case we simply take focus back again.
proc keep_focus {w} {
    set x [winfo pointerx .]
    set y [winfo pointery .]
    if {$w == [winfo containing $x $y]} {
	after idle "focus $w"
    }
}

proc CreateTableList { io t } {
    package require Tablelist
    
    global gap_defs $t.File $t.Auto $t.Edwin $t.Traces

    set $t.Auto   [keylget gap_defs CONTIG_NAVIGATOR.AUTO_CLOSE] 
    set $t.Traces [keylget gap_defs CONTIG_NAVIGATOR.SHOW_TRACES] 

    if {[xtoplevel $t] == ""} return
    wm title $t "Navigate Regions"
    
    frame $t.f
    
    set tbl $t.f.tbl
    
    tablelist::tablelist $tbl \
        -columns {12 "ContigID"  12 "Position"  0 "Problem Type"} \
        -showseparators yes \
        -background gray96 -selectbackground navy -selectforeground white \
        -width 100 -stretch 2 \
        -labelcommand tablelist::sortByColumn -sortcommand compareAsSet \
        -yscrollcommand [list $t.f.yscroll set] \
        -selectmode single \
        -exportselection 0

    $tbl columnconfigure 1 -sortmode integer
    
    
    foreach row [set $t.File] {        
        if {[regexp {^([^ ]*) +([^ ]*) *([^\x0]*)(\x0)?} $row dummy cid pos comm hit] == 0} {
            continue
        }
        $tbl insert end \
            [list "$cid" "$pos" "$comm"]        
        if {$hit == "\0"} {
            set bg [keylget gap_defs NAVIGATION.SELECT_COLOUR]
            $tbl rowconfigure end -bg $bg
        }
    }

    set body [$tbl bodypath]
    bind $body <Double-1> [list RaiseAndMoveContig $io $tbl $t]

    scrollbar $t.f.yscroll -command "$tbl yview"
    pack $t.f.yscroll -side right -fill both
    pack $t.f.tbl -side top -fill both -expand 1
    
    frame $t.b
    pack $t.b -side bottom -expand yes -fill both -padx 5

    frame $t.s -relief groove -borderwidth 2
    pack $t.s -side bottom -expand yes -fill both -padx 5 -pady 5
    
    label $t.s.head -text "Contig Length: "
    pack $t.s.head -side left    

    label $t.s.info
    pack $t.s.info -side left -padx 10

    canvas $t.s.can -width 400 -height 20
    pack $t.s.can -side right -padx 5
    
    label $t.s.cpos -text "Region Indicator: "
    pack $t.s.cpos -side right  


    $t.s.can create line 2 10 2 20 -width 2 -fill black
    $t.s.can create line 400 10 400 20 -width 2 -fill black
    $t.s.can create line 0 15 400 15 -width 2 -fill black 

    button $t.b.save -text "Save"    -command "SaveTable $t $tbl"
    button $t.b.close -text "Close"  -command "destroy $t"
    button $t.b.help -text "Help"    -command "show_help gap4 {Contig Navigation}"
    button $t.b.reset -text "Reset"  -command "ResetTable $tbl $t"
    button $t.b.prev -text "  <<-  " -command "PrevProblem $io $tbl $t" 
    button $t.b.next -text "  ->>  " -command "NextProblem $io $tbl $t" 
    bind $body <Key-p> "PrevProblem $io $tbl $t" 
    bind $body <Key-n> "NextProblem $io $tbl $t" 
    bind $body <Key-Page_Up> "PrevProblem $io $tbl $t;break" 
    bind $body <Key-Page_Down> "NextProblem $io $tbl $t;break" 
    focus $body
    bind $body <Any-FocusOut> {keep_focus %W}
    checkbutton $t.b.auto -text "Auto-close editors" \
        -variable $t.Auto -offvalue 0 -onvalue 1
    checkbutton $t.b.traces -text "Show traces" \
	-variable $t.Traces -offvalue 0 -onvalue 1
        
    pack $t.b.reset  -side left  -padx 10
    pack $t.b.prev   -side left
    pack $t.b.next   -side left
    pack $t.b.auto   -side left  -padx 10
    pack $t.b.traces -side left  -padx 10
    pack $t.b.help   -side right -padx 10
    pack $t.b.close  -side right -padx 10
    pack $t.b.save   -side right -padx 10
    
    pack $t.f -side top -expand yes -fill both -padx 5 -pady 5
}


proc SaveTable { w tbl } {

    global gap_defs $w.FileName
    
    set fname [set $w.FileName]
    if {$fname == ""} {
        bell
        return
    }
    if {[catch {set fd [open $fname w]}]} {
        bell
        tk_messageBox -message "Couldn't open $fname for writing." 
        -type ok 
        -parent $w
        return
    }
    
    set rowCount [$tbl size]
    for {set row 0} {$row < $rowCount} {incr row} {
	set ct [$tbl cellcget $row,0 -text]
        set po [$tbl cellcget $row,1 -text]
        set cm [$tbl cellcget $row,2 -text]        
        set col [$tbl rowcget $row -bg] 

        if {[$tbl rowcget $row -bg] == [keylget gap_defs NAVIGATION.SELECT_COLOUR]} {
            puts $fd "$ct $po $cm \0"
        } else {
            puts $fd "$ct $po $cm" 
        }
    }       
    close $fd  
}

proc PrevProblem { io tbl w } {
    global $w.Edwin
    
    if {[$tbl curselection] == ""} {
        $tbl selection set 0        
        RaiseAndMoveContig $io $tbl $w
        return
    }    
    set sel [$tbl curselection] 
    set pn [expr {$sel-1}]
    if {[expr {$pn < 0}]} {
        bell
        return
    }
    $tbl selection clear $sel
    $tbl selection set $pn
    $tbl see $pn
    RaiseAndMoveContig $io $tbl $w
}

proc NextProblem { io tbl w } {
    global $w.Edwin 
    
    
    if {[$tbl curselection] == ""} {
        $tbl selection set 0         
        RaiseAndMoveContig $io $tbl $w
        return
    }
    set sel [$tbl curselection]
    set pn [expr {$sel+1}]
    if {[expr {$pn == [$tbl size]}]} {
        bell
        return
    }
    $tbl selection clear $sel
    $tbl selection set $pn
    $tbl see $pn
    RaiseAndMoveContig $io $tbl $w
}

proc ResetTable { tbl w } {
    global $w.File 

    $tbl delete 0 end

    foreach row [set $w.File] {        
        if {[regexp {^([^ ]*) +([^ ]*) *([^\x0]*)} $row dummy cid pos comm] == 0} {
            continue
        }
        $tbl insert end \
            [list "$cid" "$pos" "$comm"]        
    }     
}

proc RaiseAndMoveContig {io tbl w} {
    
    global gap_defs $w.Edwin $w.Auto $w.Traces
    
    set sel [$tbl curselection]
    set ctg [$tbl cellcget [$tbl curselection],0 -text]
    set pos [$tbl cellcget [$tbl curselection],1 -text]
    set com [$tbl cellcget [$tbl curselection],2 -text]

    if {[db_info get_read_num $io $ctg] <= 0} {
        bell
        tk_messageBox -message "Contig identifier not found! Is this the correct file for this database?" \
            -type ok \
            -parent $w   
        return
    }    
    set ed [edit_contig -io $io -contig $ctg -pos $pos -reuse 1 -nojoin 1]
    
    if {[expr {[set $w.Auto] == 1}]} {
        if {[info exists $w.Edwin] && [expr {[set $w.Edwin] != $ed}]} {
	    raise [winfo toplevel [set $w.Edwin]]
	    editor_quit [winfo toplevel [set $w.Edwin]] [winfo parent [set $w.Edwin]] [set $w.Edwin] $w
        }
        set $w.Edwin $ed
    }
    set bg [keylget gap_defs NAVIGATION.SELECT_COLOUR]
    $tbl rowconfigure $sel -bg $bg

    set ctgnum [db_info get_contig_num $io $ctg]
    set ctg_len [c_length $io $ctgnum]
    $w.s.info configure -text "$ctg_len bp"
   
    if {[regexp {To:([0-9]+)$} $com _dummy_ end] != 1} {
	set end $pos
    }

    $w.s.can delete POS
    set pstart [expr {int((400.0 * $pos)/$ctg_len)}]
    set pend   [expr {int((400.0 * $end)/$ctg_len)}]
    $w.s.can create rectangle $pstart 10 $pend 20 -fill red -width 0 -tags POS

    if {[set $w.Traces] == 1} {
	$ed show_problem_traces $pos
    }
}

proc compareAsSet {item1 item2} {
    foreach {opt1 dbName1 dbClass1 default1 current1} $item1 \
        {opt2 dbName2 dbClass2 default2 current2} $item2 {
            set changed1 [expr {[string compare $default1 $current1] != 0}]
            set changed2 [expr {[string compare $default2 $current2] != 0}]
            if {$changed1 == $changed2} {
                return [string compare $opt1 $opt2]
            } elseif {$changed1} {
                return -1
            } else {
                return 1
            }
        }
}
