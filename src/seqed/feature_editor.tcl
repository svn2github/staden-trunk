

proc OpenFeatureEditorDia  { } {

    set w .ft_editor

    if {[xtoplevel $w -resizable 0] == ""} return
    wm title $w "Feature Editor" 
   
    FtEditorMenu $w $w.m
    pack $w.m -fill x 

    featureeditor $w.f    
    pack $w.f -fill both 

}

proc FtEditorMenu {w mbar} {

    frame $mbar -borderwidth 2 -relief raised
    pack $mbar -fill x

    menubutton $mbar.file -text "File" -menu $mbar.file.f
    pack $mbar.file -side left
    menu $mbar.file.f
    $mbar.file.f add command -label "Save" -command "[list SaveFT $w]"
    $mbar.file.f add separator
    $mbar.file.f add command -label "Exit" -command "destroy $w"

    menubutton $mbar.edit -text "Edit" -menu $mbar.edit.e
    pack $mbar.edit -side left
    menu $mbar.edit.e
    $mbar.edit.e add command -label "Create" -command [list CreateFT $w]
    $mbar.edit.e add command -label "Delete" -command [list DeleteFT $w]
}

proc DeleteFT {w} {
    $w.f delete
}

proc CreateFT {w} {
    $w.f create
}

proc SaveFT {w} {
    $w.f save
}



