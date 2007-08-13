# This program parses the EMBOSS acd files to produce a 'spin' Tcl/Tk
# dialogue for each program along with the necessary menu definitions.
#
# These will be written to
#     $STADENROOT/lib/spin_emboss/acdtcl/*
#     $STADENROOT/tables/emboss_menu
#
# To run it, simply run this script. NOTE that it requires DISPLAY to be set
# correctly, even though no graphics are obvious. This is because certain
# widgets are temporarily created in order to query their pathnames.

set acdtcl_dir acdtcl
set menu_file  ../../tables/emboss_menu

tkinit
wm withdraw .

# Returns the directory holding the ACD files.
#
# Uses embossdata to get the EMBOSS_DATA directory. Assumes the acd dir
# is a sibling of this.
proc acd_dir {} {
    set acd_dir ""
    set fd [open "|embossdata . 2>/dev/null"]
    set out [read $fd]
    close $fd
    foreach line [split $out "\n"] {
	if {[regexp {^([^\#].*?)/data/? *Exists$} $line dummy dir]} {
	    puts dir=$dir
	    set acd [file dirname $dir]/../share/EMBOSS/acd
	    if {![file exists $acd/wossname.acd]} {
		continue
	    }
	    set acd_dir $acd
	}
    }
    return $acd_dir
}

proc process_files {acd_dir out_dir} {
    global e_prog e_menu errs errorCode errorInfo
    foreach file [glob $acd_dir/*.acd] {
	set name [file tail $file]
	set pname [file root $name]
#	if {[catch {exec acdc $pname -acdpretty}]} {
#	    puts "Skipping $pname"
#	    continue
#	}
#	set file $pname.acdpretty
	puts "processing $name"
	catch {file delete $tnam}
	set errorCode NONE
	catch {exec stash acd2tcl.tcl $file > $out_dir/$name} menu_line
	if {$errorCode != "NONE"} {
	    puts "ERROR: could not parse"
	    puts $errorInfo
	    continue	   
	}
	if {[string match "lappend *e_menu(*" $menu_line] || \
	    [string match "set *e_prog(*" $menu_line]} {
	    # (Emboss fix) Replace PROTEIN with Protein
	    regsub {e_menu\(PROTEIN} $menu_line {e_menu(Protein} menu_line
	    eval $menu_line
	} else {
	    set errs($name) $menu_line
	    puts "ERROR: $menu_line"
	}
#	catch {file delete $pname.acdpretty}
    }
}

proc _do_emboss_menu {{out {}} {prefix {}}} {
    global e_menu e_prog
    global e_casc

    #sort menus into alphabetical order
    set menu_list [array names e_menu $prefix*]
    set menu_list [lsort -dictionary $menu_list]

    foreach i $menu_list {
	if {[regexp {^(Emboss\..*)\.[^.]*$} Emboss.$i dummy parent]} {
	    if {![info exists e_casc($parent)]} {
		append out [list add_cascade $parent] { $::EMBOSS::init 4} \n
		set e_casc($parent) 1
            }
        }
        append out [list add_cascade "Emboss.$i"] { $::EMBOSS::init 4} \n
	set e_casc(Emboss.$i) 1
	
	if {[array names e_menu $prefix$i.*] != ""} {
	    set out [_do_emboss_menu $out $prefix$i.]
	}
        foreach item [lsort -dictionary $e_menu($i)] {
	    append out [list add_command "Emboss.$i.$item"] \
		{ $::EMBOSS::init 4 } \
		[list "::EMBOSS::$e_prog($item)::create_dialogue"]\n
        }
    }

    return $out
}

# Find EMBOSS ACD source
set acd_dir [acd_dir]

# Clear things ready for new files
catch {file delete -force $acdtcl_dir}
catch {file mkdir $acdtcl_dir}

# Convert ACD to Tcl
process_files $acd_dir $acdtcl_dir
set curdir [pwd]
cd $acdtcl_dir
auto_mkindex . *.acd
cd $curdir

# Delete Test and Demo menu entries
catch {unset e_menu(Test)}
catch {unset e_menu(test)}
catch {unset e_menu(demo)}
catch {unset e_menu(Demo)}

# Build the emboss_menu file
set menu_def [_do_emboss_menu "add_menu Emboss 1 0 left\n"]
set fd [open $menu_file w]
puts $fd $menu_def
close $fd

exit
