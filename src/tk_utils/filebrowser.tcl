#invoke a filebrowser from a browse button attached to entrybox fn
#displays ALL files, no filetype has been implemented 
proc InvokeFileBrowser {fn type} {
    set filelist ""
    if {$type == "openmulti"} {
	set filelist [tk_getOpenFile -multiple 65000 -parent $fn]

	set fcount [llength $filelist]
	if { $fcount == 1 } {
	    set file [lindex $filelist 0]
	} else {
	    set file ""
	}
    } elseif {$type == "open"} {
	set file [tk_getOpenFile -parent $fn]
    } else {
	set file [tk_getSaveFile -parent $fn]
    }

    if { $file != "" || $filelist == ""} {
      if [string match [pwd]* $file] {
          set name [string range $file \
	      [expr [string length [pwd]]+1] end]
      } else {
          set name $file
      }
    } else {
      #Create a list, and keep the files there ...
      global $fn.Radio

      set baselist "#List_"
      for { set count 0 } { [ListExists2 "$baselist$count"] } { incr count } {
	  #tk_messageBox -icon error -type ok -title "DEBUG" \
		  # -message "$baselist$count"
      }
      set name "$baselist$count"
      
      ListCreate $name $filelist
      #eval set \$\{$fn.Radio\} 1
    }

    raise [winfo toplevel $fn]

    #if user has pressed cancel don't want to remove name from entrybox
    if {$name != ""} {
	entrybox_delete $fn 0 end
	entrybox_insert $fn 0 $name
	[entrybox_path $fn] xview [string last / $name]
    }

    if {$type == "save" && $file != ""} {
	DeleteFile $file
    }

}

#invoke a biolimsbrowser from a browse button attached to entrybox fn
proc InvokeBiolimsBrowser {fn type} {

    set filelist ""
    if {$type == "openmulti"} {
	set filelist [spGetOpenBiolimsFile -multiple true -parent $fn]
	set fcount [llength $filelist]
	if { $fcount == 1 } {
	    set file [lindex $filelist 0]
	} else {
	    set file ""
	}
    } elseif {$type == "open"} {
	set file [spGetOpenBiolimsFile -parent $fn]
    } else {
	set file [spGetSaveBiolimsFile -parent $fn]
    }

    if { $file != "" || $filelist == ""} {
	set name $file
    } else {
      #Create a list, and keep the files there ...
      global $fn.Radio

      set baselist "#List_"
      for { set count 0 } { [ListExists2 "$baselist$count"] } { incr count } {
      }
      set name "$baselist$count"
      
      ListCreate $name $filelist
      #eval set \$\{$fn.Radio\} 1
    }

    raise [winfo toplevel $fn]

    #if user has pressed cancel don't want to remove name from entrybox
    if {$name != ""} {
	entrybox_delete $fn 0 end
	entrybox_insert $fn 0 $name
	[entrybox_path $fn] xview [string last / $name]
    }

}

##############################################################################
#error handling procedures
##############################################################################

##############################################################################
#check the filename to be loaded already exists and is readable by user
#return 0 for failure
#return 1 for success
proc CheckOpenFile { filename {p .}} { 

    #check file exists and have read permissions
    if { [FileExists $filename $p] && [FileReadable $filename $p] } {
	    return 1
    }
    #unable to load file
    return 0
}
#end CheckOpenFile

##############################################################################
#check the filename to be saved already exists and is writable by user
#return 0 for failure
#return 1 for success
proc CheckSaveFile { filename {p .} } { 
    #check to see if file already exists and if the user wishes to
    #overwrite the existing file
    #result = no if file exists and user does not wish to overwrite
    #result = yes if file exists and user does wish to overwrite
    #result = cancel if cancel

    if {[file isdir $filename]} {
	tk_messageBox -icon error -type ok -title "Permission denied" \
		-message "$filename is a directory" \
	        -parent $p
	return 0
    }
    set result [Overwrite "$filename" $p]

    case $result {
	0 {return 0}
	1 {if {![FileWritable $filename $p]} {return 0}}
	2 {return 0}
	3 {
	    if {[catch {set fd [open $filename w]} err]} {
		tk_messageBox \
		    -icon error \
		    -type ok \
		    -title "error" \
		    -message "$err" \
		    -parent $p
		return 0
	    } else {
		close $fd
		file delete $filename
	    }
	}
    }

    return 1
}
#end CheckSaveFile

###############################################################################
#check input file exists
proc FileExists { filename {p .} } {

    set stem [file tail $filename]

    if {$stem == ""} {
	tk_messageBox -icon error -type ok -title "File does not exist" \
		-message "No filename has been entered" \
	        -parent $p
	return 0
    }

    if {![file exists $filename]} {
	tk_messageBox -icon error -type ok -title "File does not exist" \
	        -parent $p \
		-message "$filename cannot be opened. \
	Please check the filename and directory"
	return 0
    }
    return 1
}
#end FileExists


##############################################################################
#check the user has permission to read the file. 
proc FileReadable { filename {p .} } {

    if {[file readable $filename] == 0} {
	tk_messageBox -icon error -type ok -title "Permission denied" \
		-message "You do not have permission to read $filename" \
	        -parent $p
	return 0
    }
    #do have permission
    return 1
    
}
#end FileReadable

##############################################################################
#check the user has permission to write the file. 
proc FileWritable { filename {p .} } {

    if {[file writable $filename] == 0} {
	tk_messageBox -icon error -type ok -title "Permission denied" \
		-message "You do not have permission to write to $filename" \
	        -parent $p
	return 0
    }
    #do have permission
    return 1
    
}
#end FileWritable

##############################################################################
#check to see if file already exists
#return 0 if file exists and user does not wish to overwrite
#return 1 if file exists and user does wish to overwrite
#return 0 if Cancel
#return 3 if file does not exist
proc Overwrite { filename {p .} } {

    set stem [file tail $filename]

    if {$stem == ""} {
	tk_messageBox -icon error -type ok -title "File does not exist" \
		-message "No filename has been entered" \
	        -parent $p
	return 0
    }

    if {[file exists $filename]} {
	return [lsearch {no yes cancel} [tk_messageBox -icon warning -type yesnocancel \
		-default no -title "File Exists" \
		-message "Do you wish to overwrite $filename" \
		-parent $p]]
	#return [tk_dialog .fileexists "File Exists" "Do you wish to overwrite $filename" warning 0 No Yes Cancel]
    }
    
    #file does not already exist
    return 3
}
#end Overwrite
