#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1998. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#invoke a filebrowser from a browse button attached to xentry fn
#displays ALL files, no filetype has been implemented 
proc XInvokeFileBrowser {fn type args} {
    if {[set initdir [lsearch -exact $args -initialdir]] != -1} {
        set cwd [pwd]
	catch {cd [lindex $args [expr $initdir+1]]}
    }

    switch $type {
	"open" {
	    set file [eval tk_getOpenFile -parent $fn $args]
	}
	"open_multiple" {
	    set file [eval tk_getOpenFile -parent $fn $args -multiple 65000]
	}
	"save" {
	    set file [eval tk_getSaveFile -parent $fn $args]
	}
	"default" {
	    return -code error -errorinfo "Unknown action type $type"
	}
    }

    if {$initdir != -1} {
	catch {cd $cwd}
    }

    if {$file == ""} {
	return
    }

#    if {[string match [pwd]* $file]} {
#        set name [string range $file \
#	    [expr [string length [pwd]]+1] end]
#    } else {
#        set name $file
#    }
    set name $file

    raise [winfo toplevel $fn]

    $fn delete 0 end
    $fn insert 0 $name
    $fn xview [string last / $name]

    if {$type == "save" && $file != ""} {
	DeleteFile $file
    }

}


##############################################################################
#error handling procedures
##############################################################################

##############################################################################
#check the filename to be loaded already exists and is readable by user
#return 0 for failure
#return 1 for success
proc XCheckOpenFile { filename {p .} } { 

    #check file exists and have read permissions
    if { [XFileExists $filename $p] && [XFileReadable $filename $p] } {
	    return 1
    }
    #unable to load file
    return 0
}
#end XCheckOpenFile

##############################################################################
#check the filename to be saved already exists and is writable by user
#return 0 for failure
#return 1 for success
proc XCheckSaveFile { filename {p .} } { 
    set result -1

    #check to see if file already exists and if the user wishes to
    #overwrite the existing file
    #result = no if file exists and user does not wish to overwrite
    #result = yes if file exists and user does wish to overwrite
    #result = cancel if cancel
    set result [XOverwrite "$filename" $p]

    case $result {
	yes {if {![XFileWritable $filename $p]} {return 0}}
	no {return 0}
	cancel {return 0}
    }

    return 1
}
#end XCheckSaveFile

###############################################################################
#check input file exists
proc XFileExists { filename {p .} } {

    set stem [file tail $filename]

    if {$stem == ""} {
	#tk_messageBox -icon error -type ok -title "File does not exist" \
	#	-message "No filename has been entered"
	return 0
    }

    if {![file exists $filename]} {
	#tk_messageBox -icon error -type ok -title "File does not exist" \
	#	-message "$filename cannot be opened. \
	#Please check the filename and directory"
	return 0
    }
    return 1
}
#end XFileExists


##############################################################################
#check the user has permission to read the file. 
proc XFileReadable { filename {p .} } {

    if {[file readable $filename] == 0} {
	tk_messageBox -icon error -type ok -title "Permission denied" \
		-message "You do not have permission to read $filename" \
	        -parent $p
	return 0
    }
    #do have permission
    return 1
    
}
#end XFileReadable

##############################################################################
#check the user has permission to write the file. 
proc XFileWritable { filename {p .} } {

    if {[file writable $filename] == 0} {
	tk_messageBox -icon error -type ok -title "Permission denied" \
		-message "You do not have permission to write to $filename" \
	        -parent $p
	return 0
    }
    #do have permission
    return 1
    
}
#end XFileWritable

##############################################################################
#check to see if file already exists
#return no if file exists and user does not wish to overwrite
#return yes if file exists and user does wish to overwrite
#return cancel if Cancel
#return 3 if file does not exist
proc XOverwrite { filename {p .} } {

    set stem [file tail $filename]

    if {$stem == ""} {
	#tk_messageBox -icon error -type ok -title "File does not exist" \
	#	-message "No filename has been entered"
	return 0
    }

    if {[file exists $filename]} {
	return [tk_messageBox -icon warning -type yesnocancel \
		-default no -title "File Exists" \
		-message "Do you wish to overwrite $filename" \
		-parent $p]
	#return [tk_dialog .fileexists "File Exists" "Do you wish to overwrite $filename" warning 0 No Yes Cancel]
    }
    
    #file does not already exist
    return 3

}
#end XOverwrite
