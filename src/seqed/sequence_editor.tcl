#create sequence display interface

proc get_num_editor {ed_id} {

    global num_exist_editor

    if {![info exists num_exist_editor] && $ed_id != -1} {
	set num_exist_editor 0
	return $num_exist_editor
    }
    if {$ed_id != -1} {
	incr num_exist_editor
    }
    return $num_exist_editor
}

proc SequenceEditorDia { } {

    set w .seqed_d
    if {[xtoplevel $w -resizable 1] == ""} return
    wm title $w "Sequence Editor"

    iwidgets::scrolledlistbox $w.sl \
	    -labeltext "Sequence identifier:" \
	    -labelpos nw \
	    -vscrollmode dynamic \
	    -hscrollmode dynamic \
	    -selectmode extended \
	    -exportselection 0 \
	    -visibleitems 20x6 \
	    -borderwidth 2 \
	    -relief sunken     
    pack $w.sl -fill both -padx 4 -pady 4

    frame $w.se -bd 2 -relief sunken -height 2
    pack $w.se -fill x

    iwidgets::buttonbox $w.ok  
    pack $w.ok -expand yes -fill both -padx 4 -pady 1
    $w.ok add ok -text "OK" \
	    -command "SequenceEditorDia_1 $w"
    $w.ok add cancel -text "Cancel" -command "destroy $w"
    $w.ok add help -text "Help" -command "puts Help" 

    set iden [get_sequences_iden]
    $w.sl delete 0 end
    foreach id $iden {
	$w.sl insert end $id
    }
}

proc SequenceEditorDia_1 {w} {
    
    set seq_ids [$w.sl getcurselection]
    if {$seq_ids == ""} {
	tk_messageBox -icon error -type ok \
		-title "sequence editor" \
		-message "No sequence identifier has been selected"
	raise $w
	return
    }
  
    CreateSequenceEditor $seq_ids 1
    destroy $w
}

proc CreateSequenceEditor {seq_ids pos} {

    set ed_id [add_editor $seq_ids]
    set editor_num [get_num_editor $ed_id]
    
    set w .seqeditor$editor_num
    if {[xtoplevel $w -resizable 1] == ""} return
    set e_name [get_editor_name $ed_id]
    wm title $w "Sequence Editor: $e_name"
    wm protocol $w WM_DELETE_WINDOW [list quit $w]

    global $w.ted_id $w.se
    set $w.ted_id -1
    create_editor $w $ed_id $pos
    set $w.ted_id [seqed_register $w $w.se $ed_id $pos] 
    
    if {[set $w.ted_id] == -1} {
	return
    }
}

proc create_editor {w ed_id pos} {

    global $w.seq_info $w.ed_id $w.se
    #
    # create menus 
    #
    frame $w.bar -borderwidth 2 -relief flat 
    pack $w.bar -fill x 
    frame $w.bar.menu -borderwidth 2 -relief raised
    pack $w.bar.menu -fill x -side left -expand yes
    button $w.bar.undo -borderwidth 2 -relief raised \
	-text "Undo" -state disabled \
	-command [list Undo $w] \
	-takefocus 0
    pack $w.bar.undo -fill both -side right
    create_menus $w $w.bar.menu
    
    #
    # create sequence editor
    #
    set num_seq [get_num_seqs_in_editor $ed_id]

    #set height 40 ;#for ruler and scrollbar
    #           19 ;#for consensus
    #           95 ;#for every sequence and its child window
    #           10 ;#for border 
    set height [expr 40 + 19 + 114*$num_seq]
    
    texteditor $w.se \
	    -editorid $ed_id \
	    -cursorpos $pos \
	    -height $height 
  
    pack $w.se -fill both -expand yes
    set $w.seq_info [$w.se get_sequence_info $ed_id]
    set $w.ed_id $ed_id
    
    #Using this option to get signal feedback
    $w.se configure -signaltextB1Move [list sigtextB1Move $w.bar.menu]
    $w.se configure -signaltextUndo [list reset_undo $w]
    reset_undo $w
}

# ------------------------------------------------------------------
# PROCEDURE: Cut
#
# This procedure moves a selection from the document to the clipboard.
# ------------------------------------------------------------------
proc cut {w} {

    $w.se cut 

    $w.bar.menu.edit.menu entryconfigure 3 -state disabled
    $w.bar.menu.edit.menu entryconfigure 4 -state disabled
    $w.bar.menu.edit.menu entryconfigure 5 -state normal
    $w.bar.menu.edit.menu entryconfigure 6 -state disabled
}

# ------------------------------------------------------------------
# PROCEDURE: Copy
#
# This procedure makes a copy to the clipboard. (no Undo)
# ------------------------------------------------------------------
proc copy {w} {

    $w.se copy

    $w.bar.menu.edit.menu entryconfigure 3 -state disabled
    $w.bar.menu.edit.menu entryconfigure 4 -state disabled
    $w.bar.menu.edit.menu entryconfigure 5 -state normal
    $w.bar.menu.edit.menu entryconfigure 6 -state disabled
}

# ------------------------------------------------------------------
# PROCEDURE: Paste
#
# This procedure copies from the clipboard to the active sequence.
# ------------------------------------------------------------------
proc paste {w} {

    $w.se paste

    $w.bar.menu.edit.menu entryconfigure 3 -state disabled
    $w.bar.menu.edit.menu entryconfigure 4 -state disabled
    $w.bar.menu.edit.menu entryconfigure 5 -state normal
    $w.bar.menu.edit.menu entryconfigure 6 -state disabled
    $w.bar.menu.edit.menu entryconfigure 1 -state normal
    $w.bar.undo configure -state normal
}

# ------------------------------------------------------------------
# PROCEDURE:Clear
#
# This procedure removes a selection without copying to the clipboard.
# ------------------------------------------------------------------
proc clear {w} {

    $w.se clear

    $w.bar.menu.edit.menu entryconfigure 3 -state disabled
    $w.bar.menu.edit.menu entryconfigure 4 -state disabled
    $w.bar.menu.edit.menu entryconfigure 6 -state disabled
    $w.bar.undo configure -state normal
}

# ------------------------------------------------------------------
# PROCEDURE reset_undo
#
# This procedure is invoked when create new editor.
# ------------------------------------------------------------------
proc reset_undo {w} {

    global $w.ed_id

    set undo_info [$w.se get_undo_info [set $w.ed_id]]
    if {$undo_info == 1} {
	$w.bar.menu.edit.menu entryconfigure 1 -state disabled
	$w.bar.undo configure -state disabled
    }
    if {$undo_info == 0} {
	$w.bar.menu.edit.menu entryconfigure 1 -state normal
	$w.bar.undo configure -state normal 
    }   
}
# ------------------------------------------------------------------
# PROCEDURE:Undo
#
# This procedure is invoked by Undo.
# ------------------------------------------------------------------
proc Undo {w} {

    global $w.ed_id

    $w.se Undo
}

# ------------------------------------------------------------------
# PROCEDURE:set_undo
#
# This procedure is to reset undo button after editing.
# ------------------------------------------------------------------
proc set_undo {m} {

    $m.menu.edit.menu entryconfigure 1 -state normal
    $m.undo configure -state normal 
}

# ------------------------------------------------------------------
# PROCEDURE: save
#
# This procedure save the content of the editor.
# ------------------------------------------------------------------
proc editor_save {w} {

    global $w.ed_id
    $w.se save [set $w.ed_id]
    
}

# ------------------------------------------------------------------
# PROCEDURE: exit
#
# This procedure quit the editor.
# ------------------------------------------------------------------
proc quit {w} {

    global $w.ed_id
    $w.se exit [set $w.ed_id]
    
}
# ------------------------------------------------------------------
# PROCEDURE:signaltextB1Move
# This procedure is used to communicate with seqdisplay.
# To re_configure EDIT menu button after doing selection.
# ------------------------------------------------------------------
proc sigtextB1Move {m} {

    $m.edit.menu entryconfigure 3 -state normal
    $m.edit.menu entryconfigure 4 -state normal	
    $m.edit.menu entryconfigure 6 -state normal
}

# ------------------------------------------------------------------
# PROCEDURE:replace_dia
#
# This procedure to do string search and replace.
# ------------------------------------------------------------------
proc replace_dia {w} {

    global $w.ed_id

    set sr $w.string_replace	
    if {[xtoplevel $sr] == ""} return
    wm title $sr "String Replace"
    
    dataentry $sr.r\
              -labeltext "Replace" \
              -sticky e \
	      -validate alphabetic \
	      -width 20 	      
    pack $sr.r -fill both -padx 4 -pady 4
    
    dataentry $sr.w\
              -labeltext "With" \
              -sticky e \
	      -validate alphabetic \
	      -width 20 	      
    pack $sr.w -fill both -padx 4 -pady 4

    frame $sr.se -bd 2 -relief sunken -height 2
    pack $sr.se -fill x

    iwidgets::buttonbox $sr.ok  
    pack $sr.ok -expand yes -fill both -padx 4 -pady 1
    $sr.ok add yes -text "Replace" \
	    -command "replace_yes $w $sr"
    $sr.ok add no -text "No" \
	    -command "$w.se replace_no"
    $sr.ok add cancel -text "Cancel" -command "destroy $sr"
}

proc replace_yes {w path} {
    
    set_undo $w.bar
    set replace [$path.r get]	
    set with [$path.w get]
    
    $w.se replace_yes $replace $with
}

# ------------------------------------------------------------------
# PROCEDURE:string_dia
#
# This procedure to do string search.
# ------------------------------------------------------------------
proc string_dia {w} {

    set ss $w.string_search   
    if {[xtoplevel $ss] == ""} return
    wm title $ss "String Search"
   
    orientedradiobox $ss.d \
	    -labelpos nw \
	    -orient horizontal \
	    -labeltext Direction
    # orientedradiobox $ss.s \
	    -labelpos nw \
	    -orient horizontal \
	    -labeltext Strand       
    orientedradiobox $ss.a \
	    -labelpos nw \
	    -orient horizontal \
	    -labeltext "Search Algorithm"

   pack $ss.d $ss.a -side top -fill both -padx 4 -pady 2

   $ss.d add bw -text "backward"
   $ss.d add fw -text "forward"	
   $ss.d select fw   
   #$ss.s add rvs -text "reverse"
   #$ss.s add fwd -text "forward"
   #$ss.s select fwd
   $ss.a add lt -text "literal"
   $ss.a add ic -text "iub codes"
   $ss.a select ic
  
   dataentry $ss.m\
              -labeltext "Minimum percent match" \
              -sticky e \
	      -validate real \
	      -range {0 100.0} \
	      -default 75.0 \
	      -width 6 	      
    pack $ss.m -fill both -padx 4 -pady 2
   
    dataentry $ss.ss\
              -labeltext "Search string" \
              -sticky e \
	      -validate alphabetic \
	      -textvariable search_string \
	      -width 20 	      
    pack $ss.ss -fill both -padx 4 -pady 2

    frame $ss.se -bd 2 -relief sunken -height 2
    pack $ss.se -fill x

    iwidgets::buttonbox $ss.ok  
    pack $ss.ok -expand yes -fill both -padx 4 -pady 1

    $ss.ok add search -text "Search" \
	    -command "string_search $w.se $ss"
    $ss.ok add cancel -text "Cancel" -command "destroy $ss"
    $ss.ok add help -text "Help" -command "puts Help"
}

proc string_search {editor path} {

    set direction [$path.d get] ;#fw bw
    expr {($direction == "fw") ? [set dir 0 ] : [set dir 1]}
    set strand 0
    set sa  [$path.a get]        ;#lt ic
    expr {($sa == "ic") ? [set sa 1] : [set sa 0]}
    set per [$path.m get]	
    set string [$path.ss get]
    if {$string == " "} {
	bell
	return
    }
    # FIXME: should be option here to get argument
    $editor string_search $dir $strand $sa $per $string
}

# ------------------------------------------------------------------
# PROCEDURE:name_search_dia
#
# This procedure to invoke a dialogue to do sequence identifier search.
# ------------------------------------------------------------------
proc name_search_dia {w} {

    global $w.seq_info

    set ns $w.name_search
    if {[xtoplevel $ns -resizable 0] == ""} return
    wm title $ns "Name Search"
        
    iwidgets::combobox $ns.id -labeltext "Seq identifier" \
	    -labelpos w \
	    -editable false \
 	    -completion false \
	    -grab global 
    pack $ns.id -fill x -padx 4 -pady 6

    set info [set $w.seq_info]
    foreach n $info {
	$ns.id insert list end [lindex $n 0]
    }

    frame $ns.se -bd 2 -relief sunken -height 2
    pack $ns.se -fill x

    iwidgets::buttonbox $ns.ok  
    pack $ns.ok -expand yes -fill both -padx 4 -pady 1
    $ns.ok add ok -text "OK" \
	    -command "name_search $w.se $ns"
	   
    $ns.ok add cancel -text "Cancel" -command "puts Cancel"
    $ns.ok add help -text "Help" -command "puts Help"
}
proc name_search {editor path} {
    
    set iden [$path.id get]
    $editor name_search $iden
}

# ------------------------------------------------------------------
# PROCEDURE:key_search_dia
#
# This procedure to invoke a dialogue to do feature key word search.
# ------------------------------------------------------------------
proc key_search_dia {w} {

    set ks $w.key_search
    if {[xtoplevel $ks -resizable 0] == ""} return
    wm title $ks "Feature Key Search"
    
    iwidgets::combobox $ks.lb \
	    -labeltext "Keyword:" \
	    -labelpos w \
	    -completion false \
	    -grab global
    pack $ks.lb -side top -fill both -padx 4 -pady 2 
    set key [get_keyword]
    $ks.lb insert entry end [lindex $key 8]
    foreach k $key {
	$ks.lb insert list end $k
    }
    $ks.lb configure -editable false
    frame $ks.se -bd 2 -relief sunken -height 2
    pack $ks.se -fill x

    iwidgets::buttonbox $ks.ok  
    pack $ks.ok -expand yes -fill both -padx 4 -pady 1
    $ks.ok add ok -text "Search" \
	    -command "key_search $w.se $ks"
	   
    $ks.ok add cancel -text "Cancel" -command "destroy $ks"
    $ks.ok add help -text "Help" -command "puts Help"
}

proc key_search {editor path} {
    
    set key [$path.lb get]  
    $editor keyword_search $key
}

# ------------------------------------------------------------------
# PROCEDURE:qual_search_dia
#
# This procedure to invoke a dialogue to do  feature qualifier search.
# ------------------------------------------------------------------
proc qual_search_dia {w} {

    set qs $w.qual_search
    if {[xtoplevel $qs -resizable 0] == ""} return
    wm title $qs "Feature Qualifier Search"
        
    iwidgets::combobox $qs.lb \
	    -labeltext "Qualifier:" \
	    -labelpos w \
	    -completion false \
	    -grab global
    pack $qs.lb -side top -fill both -padx 4 -pady 2 
    #set key [get_keyword]
    set qual [get_qualifier]
    $qs.lb insert entry end [lindex $qual 0]
    foreach q $qual {
	$qs.lb insert list end $q
    }
    $qs.lb configure -editable false
    
    frame $qs.se -bd 2 -relief sunken -height 2
    pack $qs.se -fill x

    iwidgets::buttonbox $qs.ok  
    pack $qs.ok -expand yes -fill both -padx 4 -pady 1
    $qs.ok add ok -text "Search" \
	    -command "qual_search $w.se $qs"
	   
    $qs.ok add cancel -text "Cancel" -command "destroy $qs"
    $qs.ok add help -text "Help" -command "puts Help"
}

proc qual_search {editor path} {
    set qual [$path.lb get]
    $editor qualifier_search $qual
}

# ------------------------------------------------------------------
# PROCEDURE:show_renzyme_map
#
# This procedure to .FIXME
# ------------------------------------------------------------------


proc show_renzyme_map {w} {

    global $w.renzyme

    if {[set $w.renzyme] == 1} {
	renzyme_search_dia $w
    } 
    if {[set $w.renzyme] == 0} {
	$w.se configure -renzymemap " "
	$w.se renzyme_search 
    } 
}

# ------------------------------------------------------------------
# PROCEDURE:renzyme_search_dia
#
# This procedure to invoke a dialogue for renzyme search.
# ------------------------------------------------------------------
proc renzyme_search_dia {w} {

    set rens $w.renzyme_search

    if {[xtoplevel $rens -resizable 0] == ""} return
    wm title $rens "R.Enzyme Search"
        
    renzymebox $rens.ld
    pack $rens.ld -fill x -padx 4 -pady 6

    frame $rens.se -bd 2 -relief sunken -height 2
    pack $rens.se -fill x

    iwidgets::buttonbox $rens.ok  
    pack $rens.ok -expand yes -fill both -padx 4 -pady 1
    $rens.ok add ok -text "OK" \
	    -command "search_result_preview $w $rens"
    $rens.ok add cancel -text "Cancel" -command "destroy $rens"
    $rens.ok add help -text "Help" -command "puts Help"
}

# ------------------------------------------------------------------
# PROCEDURE:search_result_preview
#
# This procedure to pop up a window for displaying of selected 
# renzyme search result.
# ------------------------------------------------------------------
proc search_result_preview {w rens} {

    global $w.seq_info $w.ed_id

    set rsr $w.renzyme_search_result
    if {[xtoplevel $rsr -resizable 0] == ""} return
    wm title $rsr "Selected R.Enzyme Search Result"

    set sel_ren_search [$rens.ld get_name]
    destroy $rens

    set mat_all [renzyme_search_select [set $w.ed_id] $sel_ren_search]
    set i 0
    set columns { 6 "Enzyme"   
                 22 "Recognition sequence/cut site"
                 10 "Prototype"   
	          6 "Supplier_codes"                  
                  6 "EFS"}
    foreach item [lindex [lindex $mat_all 0] 5] {
	incr i
	set n [lindex [lindex [set $w.seq_info] $i] 0]
	lappend columns 5 "$n"
    }    
    tablelist::tablelist $rsr.tl \
	    -columns $columns\
            -labelcommand tablelist::sortByColumn \
	    -xscrollcommand [list $rsr.hsb set] \
	    -yscrollcommand [list $rsr.vsb set] \
	    -selectbackground navy -selectforeground white \
	    -height 10 -width 95 -stretch all \
	    -selectmode extended \
	    -exportselection 0

    set num_col [llength $columns]
    set num_col [expr $num_col/2]
    for {set i 1} {$i < $num_col} {incr i 2} {
	$rsr.tl columnconfigure $i -background beige
    }
    for {set i 4} {$i < $num_col} {incr i} {
	$rsr.tl columnconfigure $i -sortmode real
    }

    scrollbar $rsr.vsb -orient vertical -command "$rsr.tl yview"
    scrollbar $rsr.hsb -orient horizontal -command "$rsr.tl xview"
    
    grid rowconfigure    $rsr.tl 0 -weight 1
    grid columnconfigure $rsr.tl 0 -weight 1

    frame $rsr.se -bd 2 -relief sunken -height 2
  
    iwidgets::buttonbox $rsr.ok  
    $rsr.ok add ok -text "OK" \
	    -command "renzyme_search $w $rsr"

    $rsr.ok add cancel -text "Cancel" -command "destroy $rsr"
    $rsr.ok add help -text "Help" -command "puts Help"
  
    grid $rsr.tl -row 0 -column 0 -sticky news
    grid $rsr.vsb -row 0 -column 1 -sticky ns
    grid $rsr.hsb -row 1 -column 0 -sticky ew 
    grid $rsr.se -row 2 -column 0 -sticky w -columnspan 2
    grid $rsr.ok -row 3 -column 0 -sticky ew -columnspan 2

    $rsr.tl delete 0 end
    foreach items $mat_all {
	set item [lrange $items 0 4]
	foreach i [lindex $items 5] {
	    lappend item $i
	}
	$rsr.tl insert end $item
    }
    #destroy $rens
}

proc renzyme_search {w path} {

    set sel_ren ""
    set cur_sel [$path.tl curselection]
    set num_sel [llength $cur_sel]

    if {$num_sel == 0} {
	bell
	return ""
    }
    for {set i 0} {$i < $num_sel} {incr i} {
	set row_sel [lindex $cur_sel $i]
	set row_single [$path.tl get $row_sel]
	set na [lindex $row_single 0]
	set sel_ren [concat $sel_ren $na]
    }
   
    #lower $path
    destroy $path
    $w.se configure -renzymemap $sel_ren
    $w.se renzyme_search
}

# ------------------------------------------------------------------
# PROCEDURE:show_feature
#
# This procedure to .
# ------------------------------------------------------------------
proc show_feature {w} {

    global $w.features
    $w.se configure -showfeature [set $w.features]
    $w.se show_feature
}

proc SetInsertionMode {w path} {
    
    global $w.ed_id
    set sel [$path get]
    set_ft_insertion_mode [set $w.ed_id] $sel
}

proc set_insertion_mode {w} {

    set sim $w.set_insertion_mode
    
    if {[xtoplevel $sim -resizable 0] == ""} return
    wm title $sim "Set Insertion Mode"

    radiobutton .r
    set sc [.r cget -selectcolor]
    destroy .r

    iwidgets::radiobox $sim.r \
	    -labeltext "Select a Mode" \
            -selectcolor $sc
    pack $sim.r -fill both -expand yes
   
    foreach m {"extend_feature" "break_feature" "delete_feature"} {
	$sim.r add $m -text [string totitle $m]
    }
    $sim.r select "extend_feature"

    frame $sim.se -bd 2 -relief sunken -height 2
    pack $sim.se -fill x
  
    iwidgets::buttonbox $sim.ok  
    pack $sim.ok -expand yes -fill both
    $sim.ok add ok -text "OK"\
	    -command "SetInsertionMode $w $sim.r"
    $sim.ok add cancel -text "Cancel" -command "destroy $sim"
    $sim.ok add help -text "Help" -command "puts Help"
}
# ------------------------------------------------------------------
# PROCEDURE:translate
#
# This procedure to .
# ------------------------------------------------------------------
proc translate {w frameVar value args} {

    upvar #0 $frameVar frame
    set translation_frame ""

    #2:single frame
    if {$value == 2} {
	set value $frame($args)
    }
    #0:frame+ frame- all remove. 1:frame+ frame- all add.
    if {$value == 0} {
	foreach f $args {
	    set frame($f) 0 
	}
    } else {
	foreach f $args {    
	    set frame($f) 1
	}
    }
    foreach i {1 2 3 4 5 6} {
	set translation_frame [concat $translation_frame $frame($i)]
    }

    $w.se configure -translation $translation_frame
    $w.se translate 
}

proc DeleteTextEditor {w} {

    global $w.ted_id
    unset $w.ted_id 
    destroy $w
}

proc tEditorShutdown {w} {
    
    global $w.ted_id
   
    if {[info exists $w.ted_id]} {
	text_editor_result_update -index [set $w.ted_id] -job QUIT
    }
}

proc create_menus {w m} {

    global $w.features
    global $w.renzyme

    menubutton $m.file -text "File" -menu $m.file.menu
    pack $m.file -side left
    menu $m.file.menu
    $m.file.menu add command -label "Save" \
	-command [list editor_save $w]
    $m.file.menu add separator
    $m.file.menu add command -label "Exit" \
	-command [list quit $w]
    
    ###################################
    menubutton $m.edit -text "Edit" -menu $m.edit.menu
    menu $m.edit.menu
    pack $m.edit -side left
    $m.edit.menu add command -label "Undo" -state disabled \
	    -command [list Undo $w]
    $m.edit.menu add separator
    $m.edit.menu add command -label "Cut" -state disabled \
	    -command [list cut $w]
    $m.edit.menu add command -label "Copy" -state disabled \
	    -command [list copy $w]
    $m.edit.menu add command -label "Paste" -state disabled \
	    -command [list paste $w]
    $m.edit.menu add command -label "Clear" -state disabled \
	    -command [list clear $w]
    $m.edit.menu add separator
    $m.edit.menu add command -label "Replace" \
	    -command "replace_dia $w"

    #######################################
    menubutton $m.search -text "Search" -menu $m.search.menu
    menu $m.search.menu
    pack $m.search -side left
    $m.search.menu add command -label "String" \
	    -command "string_dia $w"
    ######### TO DO ################
    #$m.search.menu add command -label "Name" \
	    \-command "name_search_dia $w"
    $m.search.menu add cascade -label "Feature" \
	    -menu $m.search.menu.feat
    menu $m.search.menu.feat
    $m.search.menu.feat add command -label "Feature Key" \
	    -command "key_search_dia $w"
    $m.search.menu.feat add command -label "Qualifier" \
	    -command "qual_search_dia $w"
    $m.search.menu add check -label "Restriction enzyme map" \
	    -variable $w.renzyme \
	    -command "show_renzyme_map $w"
	    #-command "renzyme_search_changed $w"

    #########################################
    menubutton $m.feature -text "Feature" -menu $m.feature.menu
    menu $m.feature.menu
    pack $m.feature -side left
    $m.feature.menu add check -label "Show Features" \
	    -variable $w.features \
	    -command "show_feature $w"
    
    ############# TO DO ####################
    #$m.feature.menu add command -label "Add Feature" 
	    #-command "add_feat_dia"
    #$m.feature.menu add command -label "Move Feature" 
	    #-command "add_feat"    
    $m.feature.menu add command -label "Set Insertion Mode" \
	    -command "set_insertion_mode $w"
    #########################################
    menubutton $m.translate -text "Translation" -menu $m.translate.menu
    menu $m.translate.menu
    pack $m.translate -side left
    $m.translate.menu add cascade -label "Translate" \
	    -menu $m.translate.menu.tran
    menu $m.translate.menu.tran -tearoff 0
    $m.translate.menu.tran add checkbutton -label "Translate frame 1+" \
	    -variable ::frame(1) -command "translate $w frame 2 1"
    $m.translate.menu.tran add checkbutton -label "Translate frame 2+" \
	    -variable ::frame(2) -command "translate $w frame 2 2"
    $m.translate.menu.tran add checkbutton -label "Translate frame 3+" \
	    -variable ::frame(3) -command "translate $w frame 2 3"
    $m.translate.menu.tran add checkbutton -label "Translate frame 1-" \
	    -variable ::frame(4) -command "translate $w frame 2 4"
    $m.translate.menu.tran add checkbutton -label "Translate frame 2-" \
	    -variable ::frame(5) -command "translate $w frame 2 5"
    $m.translate.menu.tran add checkbutton -label "Translate frame 3-" \
	    -variable ::frame(6) -command "translate $w frame 2 6"
    $m.translate.menu.tran add separator
    $m.translate.menu.tran add command -label "Translate + frame" \
	    -command "translate $w frame 1 1 2 3"
    $m.translate.menu.tran add command -label "Translate - frame" \
	    -command "translate $w frame 1 4 5 6"
    $m.translate.menu.tran add command -label "Translat all frames" \
	    -command "translate $w frame 1 1 2 3 4 5 6"
    $m.translate.menu.tran add separator
    $m.translate.menu.tran add command -label "Remove All" \
	    -command "translate $w frame 0 1 2 3 4 5 6"

    ################## TO DO ###################
    #menubutton $m.sequence -text "Sequence" -menu $m.sequence.menu
    #menu $m.sequence.menu
    #pack $m.sequence -side left
    #$m.sequence.menu add command -label "Add Sequence" -command "puts add"
    #$m.sequence.menu add command -label "Move Sequence" -command "puts move"

}

