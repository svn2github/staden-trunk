
####################
if {[catch tkinit]} {
    package require Tk
}

if {[info exists env(STADEN_AUTO_PATH)]} {
    set auto_path "$env(STADEN_AUTO_PATH) $auto_path"
}

source $env(STADTABL)/shlib.conf
load $env(STADLIB)/$env(MACHINE)-binaries/${lib_prefix}tk_utils${lib_suffix}
load_package tk_utils
tk_utils_init
load_package seq_utils
load_package tk_utils
load_package iwidgets
load_package spin
load_package seqed
package require tablelist

wm title . "Sequence Editor v1.0b1"
wm geometry . 400x30

frame .mbar -borderwidth 2 -relief raised
pack .mbar -fill x -side top

set editor(file) [menubutton .mbar.file -text "File" -menu .mbar.file.menu]
pack .mbar.file -side left

menu .mbar.file.menu
.mbar.file.menu add command -label "Open" -command "OpenFileBrowser"
.mbar.file.menu add command -label "Save" -state disable \
    -command "SaveFileBrowser"

.mbar.file.menu add separator
.mbar.file.menu add command -label "Exit" -command "ExitSpinEditor"

set editor(view) [menubutton .mbar.view -text "Edit" \
	-menu .mbar.view.menu -state disabled]
pack .mbar.view -side left
menu .mbar.view.menu
.mbar.view.menu add command -label "Sequence" \
	-command "SequenceEditorDia"
.mbar.view.menu add command -label "Restriction enzyme" \
	-command "REnzymeMapDia"
.mbar.view.menu add command -label "Feature table" \
	-command "OpenFeatureEditorDia"

#set editor(search) [menubutton .mbar.search -text "Graphic" \
	\-menu .mbar.search.menu -state disabled]
#pack .mbar.search -side left
#menu .mbar.search.menu


#set editor(feature) [menubutton .mbar.feature -text "Feature" \
	-menu .mbar.feature.menu -state disabled]
#pack .mbar.feature -side left
#menu .mbar.feature.menu


set editor(help) [menubutton .mbar.help -text "Help" \
	-menu .mbar.help.menu]
pack .mbar.help -side right
menu .mbar.help.menu
.mbar.help.menu add command -label "Introduction" 

###############################

#create sequence browser interface
proc OpenFileBrowser { } {
   
    set f .file

    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "Load multi_sequences"

    #getfname $f.s -labeltext "Filename"
    checkedgetfname $f.s -labeltext "Filename" -type load -width 20
    pack $f.s -fill both -side top -pady 10 -padx 5

    frame $f.separator -bd 2 -relief raised -height 2
    pack $f.separator -fill x -padx 10 -pady 2

    set entry [ [$f.s component entry] component entry]

    #OkCancelHelp button
    iwidgets::buttonbox $f.b
    pack $f.b -side bottom -padx 2 -pady 2
  
    $f.b add OK -text OK -command "GetSequence $f $entry"
    $f.b add Cancel -text Cancel -command "destroy $f"
    $f.b add Help -text Help -command "puts Help"

    bind $entry <Return> "GetSequence $f $entry" 
}

proc GetSequence {f entry} {

    set fname ""
    set fname [$f.s get]    
    if {$fname == ""} {
	raise $f
	focus $entry
	return
    }

    #read sequence file using C function     
    set err [$f.s read_file $fname] 
    destroy $f
    if {$err != -1} {
	.mbar.view configure -state normal
	.mbar.file.menu entryconfigure 2 -state normal
	#.mbar.search configure -state normal
	#.mbar.feature configure -state normal
    }
}

proc SaveFileBrowser { } {

    set w .file_save
    
    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w Save

    iwidgets::combobox $w.cb \
	    -labeltext "Seq identifier:  " \
	    -labelpos w \
	    -completion false \
	    -grab global 
	  
    pack $w.cb -fill x -padx 2 -pady 10
    InitIdentifier $w.cb
    
    set entry [getFname $w.file "Output filename" save]
    entrybox_configure $entry -width 20
    pack $w.file -fill x

    frame $w.separator -bd 2 -relief raised -height 2
    pack $w.separator -fill x -padx 10 -pady 2

    iwidgets::buttonbox $w.b
    pack $w.b -side bottom -padx 2 -pady 2  
    $w.b add OK -text OK -command "Save" \
	-command "EditorSaveFile $w.cb $w.file; destroy $w"
    
    $w.b add Cancel -text Cancel -command "destroy $w"
    $w.b add Help -text Help -command "puts Help"
}

proc EditorSaveFile {in out} {

    set file_in [$in get]
    set file_out [getFname_in_name $out]
    
    editor_save_file $file_in $file_out

}

proc InitIdentifier {f} {

    set iden [get_sequences_iden]
    
    $f delete list 0 end
    $f delete entry 0 end
    $f insert entry end [lindex $iden 0]
    foreach id $iden {
	$f insert list end $id
    }

}

proc ExitSpinEditor {} {

    spin_editor_exit
    exit
}



