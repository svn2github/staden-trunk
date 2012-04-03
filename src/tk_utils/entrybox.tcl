#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc CheckInput { path } {

    set filename [expandpath [$path.entry get]]
    set response [CheckOpenFile $filename $path]

    return $response
}

proc CheckInputOptional { path } {
    set filename [expandpath [$path.entry get]]
    if {$filename != ""} {
        set response [CheckOpenFile $filename $path]
    } else {
	set response 1
    }

    return $response
}

proc CheckInputOptionalDefault { default path } {
    set filename [expandpath [$path.entry get]]
    if {$filename != $default} {
        set response [CheckOpenFile $filename $path]
    } else {
	set response 1
    }

    return $response
}

proc CheckContig {io path } {
    set stem [file rootname $path]
    set rreg [scalebox_get $stem.rreg]
    set name [$path.entry get]
	
    eval set response \{[CheckEntry $io $path $name $rreg]\}
    #response = 1 if not OK, 0 if OK
    if { $response == 1} {
	    return 0
    }
    return 1
}

proc CheckContigName {io path } {
    set name [$path.entry get]
    eval set response \{[CheckGelName $io $name $path]\}
    if { $response == 1} {
	    return 0
    }
    return 1
}

proc CheckOutput { path } {
    set filename [expandpath [$path.entry get]]
    set response [CheckSaveFile "$filename" $path]
    if {$response && [file exists "$filename"]} {
	DeleteFile $filename
    }
    return $response
}

proc CheckOutputOptional { path } {
    set filename [expandpath [$path.entry get]]
    if {$filename == ""} {
	return 1;
    }
    set response [CheckSaveFile "$filename" $path]
    if {$response && [file exists "$filename"]} {
	DeleteFile $filename
    }
    return $response
}

proc CheckDBOutput { path } {
    #HACK set file extension for new database to be ALWAYS 0
    set ext 0

    set filename [expandpath [$path.entry get]]
    set response1 [CheckSaveFile "$filename.$ext" $path]
    set response2 [CheckDBFilename $filename $path]

    #user wishes to overwrite file, must remove existing file
    if {$response1} {
	if {[file exists $filename.$ext] } {
	    DeleteFile "$filename.$ext"
	}
	if {[file exists $filename.$ext.aux] } {
	    DeleteFile "$filename.$ext.aux"
	}
	if {[file exists $filename.$ext.BUSY] } {
	    DeleteFile "$filename.$ext.BUSY"
	}
    }

    #only return true if response1 and response2 are true
    if {$response1 && $response2} {
	return 1
    } else { 
	return 0
    }

}

proc CheckStringExists { path } {
    set current [$path.entry get]
    if {[string length $current] == 0} {
	return 0
    }
    return 1
}

proc CheckString { path } {
    return 1
}

proc CheckDNAString {path} {
    set current [$path.entry get]

    if {[regexp {^[ACGTacgt]+$} $current]} {
	return 1
    }
    return 0
}

proc CheckIUBCString {path} {
    set current [$path.entry get]

    if {[regexp {^[ACGTRYMKSWBDHVNacgtrymkswbdhvn?-]+$} $current]} {
	return 1
    }
    return 0
}

proc CheckInt { path } {
    global entry_

    set current [$path.entry get]
    if {![isinteger $current]} {
	return 0
    }
    return 1
}
proc CheckIntRange { min max path } {
    global entry_

    set current [$path.entry get]
    if {![isinteger $current]} {
	return 0
    }

    if {$current < $min || $current > $max} {
	return 0
    }
	
    return 1
}

proc CheckIntMin { min path } {
    global entry_

    set current [$path.entry get]
    if {![isinteger $current]} {
	return 0
    }

    if {$current < $min} {
	return 0
    }
	
    return 1
}
proc CheckIntMax { max path } {
    global entry_

    set current [$path.entry get]
    if {![isinteger $current]} {
	return 0
    }

    if {$current > $max} {
	return 0
    }
	
    return 1
}

proc CheckFloat { path } {
    global entry_

    set current [$path.entry get]
    if {![isfloat $current]} {
	return 0
    }
    return 1
}

proc CheckFloatRange {  min max path} {
    global entry_

    set current [$path.entry get]
    if {![isfloat $current]} {
	return 0
    }

    if {$current < $min || $current > $max} {
	return 0
    }
    return 1
}

proc CheckList {path} {
    return [ListExists [$path.entry get]]
}

#
# Configures an entrybox
#
proc entrybox_configure {path args} {
    global entry_
    global CurContig

    set in_arg 0
    set arglist ""
    set title ""
    set foreground ""
    set state ""
    set default ""
    set default_set 0
    set width ""
    set font ""
    set $path.Type ""

    # Process command line args
    foreach i $args {
	if {$in_arg} {

	     if {$option == "-title"} {
		set title "-text {$i}"
	     } elseif {$option == "-state"} {
		set state "-state $i"
	     } elseif {$option == "-type"} {
		 if {"$i" != ""} {
		     set entry_($path,type) "$i"
		 } else {
		     catch {unset entry_($path,type)}
		 }
	     } elseif {$option == "-default"} {
		 set default $i
		 incr default_set
	     } elseif {$option == "-width"} {
		set width "-width $i"
	     } elseif {$option == "-font"} {
		set font "-font {$i}"
	     } elseif {$option == "-command"} {
		set entry_($path,command) "$i"
	     } else {
		lappend arglist $option $i
	     }
	     set in_arg 0
	} else {
	     set option $i
	     set in_arg 1
	}
    }
    eval $path.entry configure $arglist

    if {"$title" != "" || "$foreground" != "" || "$font" != ""} {
        eval $path.label configure $title $foreground $font
    }

    if {"$state" != "" || "$foreground" != "" || "$width" != ""} {
        eval $path.entry configure $state $foreground $width
	eval $path.label configure $state $foreground
    }

    if {$default_set != 0} {
	$path.entry delete 0 end
	$path.entry insert 0 $default
	#eval $path.entry insert 0 $default
    }

    if {[info exists entry_($path,command)]} {
        bind $path.entry <Return> "entrybox_command $path"
        catch {bind $path.entry <KP_Enter> "entrybox_command $path"}
    } else {
	bind $path.entry <Return> "entrybox_get $path"
	catch {bind $path.entry <KP_Enter> "entrybox_get $path"}
    }
}

#
# Creates an entrybox with an associated label.
# Command line arguments are as per the frame widget, with the addition of:
#	-state	 normal/disabled
#       -title   title
#       -type    check function - eg "CheckInt", "CheckIntRange 1 10".
#	-width	 width
#	-font	 font for label only
#       -command for binding entrybox to Return
#	-default default value
#
bind EntryBox <Destroy> {
    catch {unset entry_(%W,type)};
    catch {unset entry_(%W,command)};
}

proc entrybox {path args} {
    global entry_

    # Create the frame, label and entrybox
    frame $path -class EntryBox
    xlabel $path.label -anchor w
    entry $path.entry

#    bind $path.entry <Delete> {tkEntryBackspace %W}

    # Configure	
    eval entrybox_configure $path $args

    # Pack and return
    pack $path.label -side left -fill x
    pack $path.entry -side right -fill x

    return $path
}

#
# entrybox_get:
#
# Given a entrybox path we return the current value
#
proc entrybox_get {path} {
    global entry_
    set re_enter 0

    if {[info exists entry_($path,type)]} {
	if {[eval $entry_($path,type) $path] == 0} {
	    raise [winfo toplevel $path]
	    focus $path.entry
	    $path.entry icursor end
	    bell
#	    #wait forever...
#	    tkwait variable re_enter
	    return ""
	}
    }
    return [$path.entry get]
}

#
# entrybox_insert:
#
# Given a entrybox path and text we insert the text in the entry box
#
proc entrybox_insert {path args} {

    eval $path.entry insert $args
}

#
# entrybox_delete:
#
# Given a entrybox path and text we delete the text in the entry box
#
proc entrybox_delete {path args} {

    eval $path.entry delete $args
}

#
# entrybox_focus:
#
# Given a entrybox path, put focus on entry box
#
proc entrybox_focus {path} {

    focus $path.entry
}

#
# entrybox_destroy:
#
# Destroys a entrybox path
#
proc entrybox_destroy {path} {
    destroy $path
}

#
# entrybox_command:
#
# Called when Return is pressed in the entry
#
proc entrybox_command {path} {
    global entry_

    # Check for entry_(path,command) and execute
    if {[info exists entry_($path,command)]} {
	eval $entry_($path,command) \{[$path.entry get]\}
	#eval $entry_($path,command) \{[entrybox_get $path]\}
    }
}

#
#returns the pathname for the widget
#
proc entrybox_path {path} {
    return $path.entry
}
