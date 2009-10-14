# biolimsbox.tcl --
#
#	Implements the "Motif" style file selection dialog 
#       to allow selection from a biolims database
#       created from SCCS: @(#) xmfbox.tcl 1.6 97/10/01 15:06:07
#

# spBiolimsDialog --
#
#	Implements a file dialog similar to the standard Motif file
#	selection box.
#
# Return value:
#
#	A list of two members. The first member is the absolute
#	pathname of the selected file or "" if user hits cancel. The
#	second member is the name of the selected file type, or ""
#	which stands for "default file type"
#
#       if -multiple is set to true, allows multiple selection, and returns
#       a list of selected files, and no filetype
#

load_package biolimstcl

proc spGetOpenBiolimsFile {args} {
    return [eval spCommonBiolimsDialog open $args]
}

proc spGetSaveBiolimsFile {args} {
    return [eval spCommonBiolimsDialog save $args]
}

proc spCommonBiolimsDialog {type args} {
    global tkPriv
    set w __sp_biolimsdialog
    upvar #0 $w data

    set oldGrab [grab current]
    grab release $oldGrab	

    spBiolimsDialog_Config $w $type $args

    if ![string compare $data(-parent) .] {
        set w .$w
    } else {
        set w $data(-parent).$w
    }

    # (re)create the dialog box if necessary
    #
    if ![winfo exists $w] {
	spBiolimsDialog_Create $w
    } else {
	destroy $w
	spBiolimsDialog_Create $w
    }     

    # Connect to Biolims, and move to the Correct Collection
    blbrowse start
    blbrowse set $data(selectPath)

    wm transient $w $data(-parent)

    if {[spBiolimsDialog_FileTypes $w] == 0} {
        spBiolimsDialog_Update $w
    }

    # 5. Withdraw the window, then update all the geometry information
    # so we know how big it wants to be, then center the window in the
    # display and de-iconify it.

    wm withdraw $w
    update idletasks
    set x [expr {[winfo screenwidth $w]/2 - [winfo reqwidth $w]/2 \
	    - [winfo vrootx [winfo parent $w]]}]
    set y [expr {[winfo screenheight $w]/2 - [winfo reqheight $w]/2 \
	    - [winfo vrooty [winfo parent $w]]}]
    wm geom $w +$x+$y
    wm deiconify $w
    wm title $w $data(-title)

    # 6. Set a grab and claim the focus too.

    catch {grab $oldGrab}
    set oldFocus [focus]
    set oldGrab [grab current $w]
    if {$oldGrab != ""} {
	set grabStatus [grab status $oldGrab]
    }
    grab $w
    focus $data(sEnt)
    $data(sEnt) select from 0
    $data(sEnt) select to   end

    # 7. Wait for the user to respond, then restore the focus and
    # return the index of the selected button.  Restore the focus
    # before deleting the window, since otherwise the window manager
    # may take the focus away so we can't redirect it.  Finally,
    # restore any grab that was in effect.

    tkwait variable tkPriv(selectFilePath)
    catch {focus $oldFocus}
    grab release $w
    wm withdraw $w
    if {$oldGrab != ""} {
	if {$grabStatus == "global"} {
	    grab -global $oldGrab
	} else {
	    grab $oldGrab
	}
    }

    return $tkPriv(selectFilePath)
}

proc spBiolimsDialog_Config {w type argList} {
    upvar #0 $w data

    set data(type) $type

    # 1: the configuration specs
    #
    set specs {
	{-defaultextension "" "" ""}
	{-filetypes "" "" ""}
	{-initialdir "" "" ""}
	{-initialfile "" "" ""}
	{-parent "" "" "."}
	{-title "" "" ""}
	{-defaulttype "" "" "0"}
	{-multiple "" "" ""}
        {-assembly "" "" ""}
    }

    # 2: default values depending on the type of the dialog
    #
    if {![info exists data(selectPath)]} {
	# first time the dialog has been popped up
	set data(selectPath) "/"
	set data(selectFile) ""
    }

    # 3: parse the arguments
    #
    tclParseConfigSpec $w $specs "" $argList

    if {![string compare $data(-title) ""]} {
	if {![string compare $type "open"]} {
            if ![string compare $data(-multiple) "true"] {
		set data(-title) "Open Multiple BioLIMS Files"
            } else {
	        set data(-title) "Open BioLIMS File"
	    }  
        } else {
	    set data(-title) "Save to BioLIMS"
	}
    }

    # 4: set the default directory and selection according to the -initial
    #    settings
    #
    set data(selectFile) $data(-initialfile)

    # 5. Parse the -filetypes option. It is not used by the motif
    #    file dialog, but we check for validity of the value to make sure
    #    the application code also runs fine with the TK file dialog.
    #
    set data(-filetypes) [tkFDGetFileTypes $data(-filetypes)]

    if {![info exists data(filter)]} {
	set data(filter) *
    }
    if {![winfo exists $data(-parent)]} {
	error "bad window path name \"$data(-parent)\""
    }
}

proc spBiolimsDialog_Create {w} {
    set dataName [lindex [split $w .] end]
    upvar #0 $dataName data

    # 1: Create the dialog ...
    #
    toplevel $w -class spBiolimsDialog
    set top [frame $w.top -relief raised -bd 1]
    set bot [frame $w.bot -relief raised -bd 1]

    pack $w.bot -side bottom -fill x
    pack $w.top -side top -expand yes -fill both

    set f1 [frame $top.f1]
    set f2 [frame $top.f2]
    set f3 [frame $top.f3]

    pack $f1 -side top    -fill x
    pack $f3 -side bottom -fill x
    pack $f2 -expand yes -fill both

    set f2a [frame $f2.a]
    set f2b [frame $f2.b]

    grid $f2a -row 0 -column 0 -rowspan 1 -columnspan 1 -padx 4 -pady 4 \
	-sticky news
    grid $f2b -row 0 -column 1 -rowspan 1 -columnspan 1 -padx 4 -pady 4 \
	-sticky news
    grid rowconfig $f2 0    -minsize 0   -weight 1
    grid columnconfig $f2 0 -minsize 0   -weight 1
    grid columnconfig $f2 1 -minsize 150 -weight 2

    # The Filter box
    #
    label $f1.lab -text "Filter:" -under 3 -anchor w
    entry $f1.ent
    pack $f1.lab -side top -fill x -padx 6 -pady 4
    pack $f1.ent -side top -fill x -padx 4 -pady 0
    set data(fEnt) $f1.ent

    # The file and directory lists
    #
    set data(dList) [spBiolimsDialog_MakeSList $w $f2a Directory: 0 DList browse]
    if ![string compare $data(-multiple) "true"] {
      set data(fList) [spBiolimsDialog_MakeSList $w $f2b Files:     2 FList extended]
    } else {
      set data(fList) [spBiolimsDialog_MakeSList $w $f2b Files:     2 FList browse]
    }
    # The Selection box
    #
    label $f3.lab -text "Selection:" -under 0 -anchor w
    entry $f3.ent
    pack $f3.lab -side top -fill x -padx 6 -pady 0
    pack $f3.ent -side top -fill x -padx 4 -pady 4
    set data(sEnt) $f3.ent

    # The buttons
    #
    set data(okBtn) [button $bot.ok     -text OK     -width 6 -under 0 \
	-command "spBiolimsDialog_OkCmd $w"]
    set data(filterBtn) [button $bot.filter -text Filter -width 6 -under 0 \
	-command "spBiolimsDialog_FilterCmd $w"]
    set data(cancelBtn) [button $bot.cancel -text Cancel -width 6 -under 0 \
	-command "spBiolimsDialog_CancelCmd $w"]

    pack $bot.ok $bot.filter $bot.cancel -padx 10 -pady 10 -expand yes \
	-side left

    # Create the bindings:
    #
    bind $w <Alt-t> "focus $data(fEnt)"
    bind $w <Alt-d> "focus $data(dList)"
    bind $w <Alt-l> "focus $data(fList)"
    bind $w <Alt-s> "focus $data(sEnt)"

    bind $w <Alt-o> "tkButtonInvoke $bot.ok    "
    bind $w <Alt-f> "tkButtonInvoke $bot.filter"
    bind $w <Alt-c> "tkButtonInvoke $bot.cancel"

    bind $data(fEnt) <Return> "spBiolimsDialog_ActivateFEnt $w"
    bind $data(sEnt) <Return> "spBiolimsDialog_ActivateSEnt $w"

    wm protocol $w WM_DELETE_WINDOW "spBiolimsDialog_CancelCmd $w"
}

# Set the file types. Returns 1 if this has already udated the file list.
# Otherwise returns 0
proc spBiolimsDialog_FileTypes {w} {
    upvar #0 [winfo name $w] data

    if {$data(-filetypes) == ""} {
	set data(filter) *
	return 0
    }

    set f $w.top.f3.types
    if [winfo exists $f] {
	destroy $f
    }
    frame $f

    # The filetypes radiobuttons
    global fileType
    set fileType $data(-defaulttype)

    spBiolimsDialog_SetFilter $w [lindex $data(-filetypes) $fileType]
    set cnt 0

    #don't produce radiobuttons for only one filetype
    if {[llength $data(-filetypes)] == 1} {
	return 1
    }
    if {$data(-filetypes) != {}} {
	foreach type $data(-filetypes) {
	    set title  [lindex [lindex $type 0] 0]
	    set filter [lindex $type 1]
	    radiobutton $f.b$cnt -text $title -variable fileType \
		    -value $cnt \
		    -command "[list spBiolimsDialog_SetFilter $w $type]"
	    pack $f.b$cnt -side left
	    incr cnt
	}
    }

    pack $f -side bottom -fill both
    return 1
}


# This proc gets called whenever data(filter) is set
#
proc spBiolimsDialog_SetFilter {w type} {
    upvar #0 [winfo name $w] data
    global tkpriv

    set data(filter) [lindex $type 1]
    set tkpriv(selectFileType) [lindex [lindex $type 0] 0]

    spBiolimsDialog_Update $w
}

proc spBiolimsDialog_MakeSList {w f label under cmd selectmode} {
    label $f.lab -text $label -under $under -anchor w
    listbox $f.l -width 12 -height 5 -selectmode $selectmode -exportselection 0\
	-xscrollcommand "$f.h set" \
	-yscrollcommand "$f.v set" 
    scrollbar $f.v -orient vertical   -takefocus 0 \
	-command "$f.l yview"
    scrollbar $f.h -orient horizontal -takefocus 0 \
	-command "$f.l xview"
    grid $f.lab -row 0 -column 0 -sticky news -rowspan 1 -columnspan 2 \
	-padx 2 -pady 2
    grid $f.l -row 1 -column 0 -rowspan 1 -columnspan 1 -sticky news
    grid $f.v -row 1 -column 1 -rowspan 1 -columnspan 1 -sticky news
    grid $f.h -row 2 -column 0 -rowspan 1 -columnspan 1 -sticky news

    grid rowconfig    $f 0 -weight 0 -minsize 0
    grid rowconfig    $f 1 -weight 1 -minsize 0
    grid columnconfig $f 0 -weight 1 -minsize 0

    # bindings for the listboxes
    #
    set list $f.l
    bind $list <Up>        "spBiolimsDialog_Browse$cmd $w"
    bind $list <Down>      "spBiolimsDialog_Browse$cmd $w"
    bind $list <space>     "spBiolimsDialog_Browse$cmd $w"
    bind $list <1>         "spBiolimsDialog_Browse$cmd $w"
    bind $list <B1-Motion> "spBiolimsDialog_Browse$cmd $w"
    bind $list <Double-1>  "spBiolimsDialog_Activate$cmd $w"
    bind $list <Return>    "spBiolimsDialog_Browse$cmd $w; spBiolimsDialog_Activate$cmd $w"

    bindtags $list "Listbox $list [winfo toplevel $list] all"
    tkListBoxKeyAccel_Set $list

    return $f.l
}


proc spBiolimsDialog_BrowseDList {w} {
    upvar #0 [winfo name $w] data

    focus $data(dList)
    if {![string compare [$data(dList) curselection] ""]} {
	return
    }
    set subdir [$data(dList) get [$data(dList) curselection]]
    if {![string compare $subdir ""]} {
	return
    }

    $data(fList) selection clear 0 end
    $data(sEnt) delete 0 end

    set list [spBiolimsDialog_InterpFilter $w]
    set data(filter) [lindex $list 1]

    case $subdir {
	. {
	    set newSpec "$data(selectPath)$data(filter)"
	}
	.. {
	    set newSpec \
		    "[string range $data(selectPath) 0 \
		      [string last / [string range $data(selectPath) 0 \
	                [expr [string length $data(selectPath)]-2]]]]$data(filter)"
	}
	default {
	    set newSpec "$data(selectPath)$subdir/$data(filter)"
	}
    }

    $data(fEnt) delete 0 end
    $data(fEnt) insert 0 $newSpec
}

proc spBiolimsDialog_ActivateDList {w} {
    upvar #0 [winfo name $w] data

    if {![string compare [$data(dList) curselection] ""]} {
	return
    }
    set subdir [$data(dList) get [$data(dList) curselection]]
    if {![string compare $subdir ""]} {
	return
    }
    if {![string compare $subdir "."]} {
	return
    }

    $data(fList) selection clear 0 end

    case $subdir {
	.. {
	    blbrowse set ..
	    set newDir [string range $data(selectPath) 0 \
		    [string last / [string range $data(selectPath) 0 \
	               [expr [string length $data(selectPath)]-2]]]]
	}
	default {
	    blbrowse set $subdir
	    set newDir "$data(selectPath)$subdir/"
	}
    }

    set data(selectPath) "$newDir"
    set data(selectFile) ""
    spBiolimsDialog_Update $w

    if {[string compare $subdir ..]} {
	$data(dList) selection set 0
	$data(dList) activate 0
    } else {
	$data(dList) selection set 1
	$data(dList) activate 1
    }
}

proc spBiolimsDialog_BrowseFList {w} {
    upvar #0 [winfo name $w] data

    focus $data(fList)
    if {![string compare [$data(fList) curselection] ""]} {
	return
    }
    set sep ""
    set data(selectFile) ""
    foreach { sel } [$data(fList) curselection] {
       set data(selectFile)  [concat $data(selectFile)$sep[$data(fList) get $sel]]
       set sep ":"
    }
    if {![string compare $data(selectFile) ""]} {
	return
    }

    $data(dList) selection clear 0 end

    $data(fEnt) delete 0 end
    $data(fEnt) insert 0 [file join $data(selectPath) $data(filter)]
    $data(fEnt) xview end

    $data(sEnt) delete 0 end
    # if it's a multiple selection box, just put in the filenames
    # otherwise put in the full path as usual
    if { ![ string compare $data(-multiple) "true" ] } {
        $data(sEnt) insert 0 $data(selectFile)
    } else {
        $data(sEnt) insert 0 [file join $data(selectPath) $data(selectFile)]
    }
    $data(sEnt) xview end
}

proc spBiolimsDialog_ActivateFList {w} {
    upvar #0 [winfo name $w] data

    if {![string compare [$data(fList) curselection] ""]} {
	return
    }
    set data(selectFile) [$data(fList) get [$data(fList) curselection]]
    if {![string compare $data(selectFile) ""]} {
	return
    } else {
	spBiolimsDialog_ActivateSEnt $w
    }
}

proc spBiolimsDialog_ActivateFEnt {w} {
    upvar #0 [winfo name $w] data

    set list [spBiolimsDialog_InterpFilter $w]
    set data(selectPath) [lindex $list 0]
    blbrowse set $data(selectPath)
    set data(filter)    [lindex $list 1]

    spBiolimsDialog_Update $w
}

proc spBiolimsDialog_InterpFilter {w} {
    upvar #0 [winfo name $w] data

    set text [string trim [$data(fEnt) get]]
    # Perform tilde substitution
    #
    if {![string compare [string index $text 0] ~]} {
	set list [file split $text]
	set tilde [lindex $list 0]
	catch {
	    set tilde [glob $tilde]
	}
	set text [eval file join [concat $tilde [lrange $list 1 end]]]
    }

    set resolved [file join [file dirname $text] [file tail $text]]

    if {[file isdirectory $resolved]} {
	set dir $resolved
	set fil $data(filter)
    } else {
	set dir [file dirname $resolved]
	set fil [file tail    $resolved]
    }

    return [list $dir $fil]
}

proc spBiolimsDialog_ActivateSEnt {w} {
    global tkPriv
    upvar #0 [winfo name $w] data

    #if not multiple file, behave as usual
    #if multiple files have been entered throw out any dirnames in filelist
    #and use the path selected
    set selectFilePath [string trim [$data(sEnt) get]]
    set selectFile     [file tail    $selectFilePath]
    set selectPath     [file dirname $selectFilePath]
    if ![string compare $data(-assembly) "true"] {
	set fileType "assemblies"
    } else {
	set fileType "samples"
    }

    #If multiple filenames are allowed
    if ![string compare $data(-multiple) "true"] {

        #If multiple files are specified
        if {[string first : $selectFilePath] != -1} {
            #makesure there are no pathnames
            if {[string first / $selectFilePath] != -1 } {
		tk_messageBox -icon warning -type ok \
			-message "Pathnames not allowed in Selection with multiple files"
		return
	    }
	}
        #use stored selectPath, as there shouldn't be one in the select box
        set selectPath  $data(selectPath)
    } else {
        #use path in the select box
	append selectPath "/"
    }

    #nothing specified
    if {![string compare $selectFilePath ""]} {
	spBiolimsDialog_FilterCmd $w
	return
    }

    #see if it is a directory
    if ![catch {blbrowse set $selectFilePath}] {
        #Add the ful path if required
        if { [string first / $selectFilePath] != 0 } {
	    set selectFilePath "$selectPath$selectFilePath"
	}
	set data(selectPath) "$selectFilePath/"
	set data(selectFile) ""
	spBiolimsDialog_Update $w
	return
    }

    #check the path
    if { [string first / $selectPath] != 0 } {
	tk_messageBox -icon warning -type ok \
	    -message "\"$selectPath\" must be an absolute pathname"
	return
    }

    #cd to the required directory
    if [catch {blbrowse set $selectPath}] {
	tk_messageBox -icon warning -type ok \
	    -message "Collection \"$selectPath\" does not exist."
	return
    }
    set validSamples [blbrowse list $fileType]

    #check each file
    if ![string compare $data(type) open] {
        set selectFilePath ""

	foreach {filename} [split $selectFile : ] {
	    set ffound 0
	    foreach f $validSamples {
		if ![string compare $f $filename] {
		    set ffound 1
		    break;
		}
	    }
            if {$ffound != 1} {
	        tk_messageBox -icon warning -type ok \
		    -message "\"$filename\" does not exist."
	        return
	    }
            lappend selectFilePath "Biolims=$selectPath$filename"
        }
    } else {
        #can't save to multiple files, so don't need to loop here
	set ffound 0
	foreach f $validSamples {
	    if  ![string compare $f $selectFile] {
		set ffound 1
		break;
	    }
	}
	if {$ffound == 1} {
	    set message [format %s%s \
		"\"$selectFilePath\" already exists.\n\n" \
		"Replace existing ?"]
	    set answer [tk_messageBox -icon warning -type yesno \
		-message $message]
	    if ![string compare $answer "no"] {
		return
	    }
	}
	set selectFilePath "Biolims=$selectFilePath"
    }

    set tkPriv(selectFilePath) $selectFilePath
    set tkPriv(selectFile)     $selectFile
    set tkPriv(selectPath)     $selectPath
    blbrowse stop
}


proc spBiolimsDialog_OkCmd {w} {
    upvar #0 [winfo name $w] data

    spBiolimsDialog_ActivateSEnt $w
}

proc spBiolimsDialog_FilterCmd {w} {
    upvar #0 [winfo name $w] data

   spBiolimsDialog_ActivateFEnt $w
}

proc spBiolimsDialog_CancelCmd {w} {
    global tkPriv

    set tkPriv(selectFilePath) ""
    set tkPriv(selectFile)     ""
    set tkPriv(selectPath)     ""
    blbrowse stop
}

# spBiolimsDialog_Update
#
#	Load the files and synchronize the "filter" and "selection" fields
#	boxes.
#
# popup:
#	If this is true, then update the selection field according to the
#	"-selection" flag
#
proc spBiolimsDialog_Update {w} {
    upvar #0 [winfo name $w] data

    $data(fEnt) delete 0 end
    $data(fEnt) insert 0 [file join $data(selectPath) $data(filter)]
    $data(sEnt) delete 0 end
    $data(sEnt) insert 0 [file join $data(selectPath) $data(selectFile)]
 
    spBiolimsDialog_LoadFiles $w
}

# JKB: Major rewrite to use our own dir_or_file function. This is MUCH faster
# than the original Tk version.
proc spBiolimsDialog_LoadFiles {w} {
    upvar #0 [winfo name $w] data

    if ![string compare $data(-assembly) "true"] {
	set fileType "assemblies"
    } else {
	set fileType "samples"
    }

    $data(dList) delete 0 end
    $data(fList) delete 0 end

    # Make the dir list
    $data(dList) insert end "."
    if {[string compare $data(selectPath) /]} {
	$data(dList) insert end ".."
    }
    eval $data(dList) insert end [blbrowse list collections]

    # Make the file list
    if {![string compare $data(filter) *]} {
	set files [blbrowse list $fileType]
    } else {
	set files {}
	foreach f [blbrowse list $fileType] {
	    if {[string match $data(filter) $f]} {
		lappend files $f
	    }
	}
    }
    set tlist [lsort -ascii $files]

    set top [lsearch -regexp $tlist {^[^.]}]
    if {$top == -1} {
	set top [llength $tlist]
    }
    eval $data(fList) insert end $tlist

    # The user probably doesn't want to see the . files. We adjust the view
    # so that the listbox displays all the non-dot files
    $data(fList) yview $top

}

proc tkListBoxKeyAccel_Set {w} {
    bind Listbox <Any-KeyPress> ""
    bind $w <Destroy> "tkListBoxKeyAccel_Unset $w"
    bind $w <Any-KeyPress> "tkListBoxKeyAccel_Key $w %A"
}

proc tkListBoxKeyAccel_Unset {w} {
    global tkPriv

    catch {after cancel $tkPriv(lbAccel,$w,afterId)}
    catch {unset tkPriv(lbAccel,$w)}
    catch {unset tkPriv(lbAccel,$w,afterId)}
}

proc tkListBoxKeyAccel_Key {w key} {
    global tkPriv

    append tkPriv(lbAccel,$w) $key
    tkListBoxKeyAccel_Goto $w $tkPriv(lbAccel,$w)
    catch {
	after cancel $tkPriv(lbAccel,$w,afterId)
    }
    set tkPriv(lbAccel,$w,afterId) [after 500 tkListBoxKeyAccel_Reset $w]
}

#johnt
#fixed so that control/shift dont reset selections
proc tkListBoxKeyAccel_Goto {w string} {
    global tkPriv

    if { [string compare $string ""] } {
      set string [string tolower $string]
      set end [$w index end]
      set theIndex -1

      for {set i 0} {$i < $end} {incr i} {
	set item [string tolower [$w get $i]]
	if {[string compare $string $item] >= 0} {
	    set theIndex $i
	}
	if {[string compare $string $item] <= 0} {
	    set theIndex $i
	    break
	}
      }

      if {$theIndex >= 0} {
	$w selection clear 0 end
	$w selection set $theIndex $theIndex
	$w activate $theIndex
	$w see $theIndex
      }
    }
}

proc tkListBoxKeyAccel_Reset {w} {
    global tkPriv

    catch {unset tkPriv(lbAccel,$w)}
}

proc tk_getFileType {} {
    global tkpriv

    return $tkpriv(selectFileType)
}




