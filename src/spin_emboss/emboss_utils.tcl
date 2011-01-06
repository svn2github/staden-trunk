namespace eval ::EMBOSS {

# -----------------------------------------------------------------------------
# Dialogue generation functions

#
# Calculates the appropriate window size for a notebook, so as to avoid
# needing to scroll the tabs.
#
# Assumes that all packing within the pages is -side top or -side bottom.
# If this is not the case then enclose the children within a single frame.
proc resizebook {book} {
    return; #disabled for now to see if ttk::notebook has issues.

    update idletasks
    update
    set bd [expr {2*[$book cget -borderwidth]}]
    set maxheight 0
    set maxwidth [expr {int(ceil([lindex [[$book component tabset] bbox] 2]))}]

    set ts [$book component tabset]
    set ntabs [expr {[$ts index end]+1}]
    # Calculate the width due to the sloping tab. For N tabs we have N+1 of
    # these width portions.
    # Assume the same tab angle and height is used throughout
    set angle [$ts tabcget 0 -angle]
    set labelHeight [expr {[$ts tabcget 0 -height] + 2*[$ts tabcget 0 -pady]}]
    set angleOffset [expr {int(ceil($labelHeight * tan($angle/57.2957795133)))}]
    set total_wid $angleOffset
    for {set tab 0} {$tab < $ntabs} {incr tab} {
	set font [$ts tabcget $tab -font]
	set padx [$ts tabcget $tab -padx]
	set text [$ts tabcget $tab -label]
	set labelWidth [font measure $font $text]
	# +7? not sure what this is, but we need a bit more per tab...
	incr total_wid [expr {$angleOffset + $labelWidth+7 + 2*$padx}]
    }
    incr total_wid [$ts cget -start]
    set maxwidth $total_wid

    # As an alternative to using bbox:
    # set ts [$book component tabset]
    # set ntabs [expr {[$ts index end]+1}]
    # set _labelWidth  [expr {[$ts tabcget 0 -width] +2*[$ts tabcget 0 -padx]}]
    # set _labelHeight [expr {[$ts tabcget 0 -height]+2*[$ts tabcget 0 -pady]}]
    # set angle [$ts tabcget 0 -angle]
    # set angleOffset [expr {$_labelHeight * tan($angle/57.2957795133)}]
    # set maxwidth [expr {int(ceil($ntabs * ($_labelWidth + $angleOffset) + $angleOffset + 4))}]

    foreach page [$book childsite] {
        set height [expr {$bd+[winfo reqheight [$book component tabset]]}]
        foreach w [winfo children $page] {
	    if {[winfo class $w] == "Tabnotebook"} {
	        resizebook $w
	    }
            set height [expr {$height+[winfo reqheight $w]}]
	    #puts $w:[winfo reqwidth $w]:[winfo class $w]
	    if {$maxwidth < [winfo reqwidth $w]} {
	        set maxwidth [winfo reqwidth $w]
	    }
        }
	if {$maxheight < $height} {
	    set maxheight $height
        }
    }
    $book configure -height $maxheight
    $book configure -width $maxwidth
    update idletasks
}

#
# Fetches a list of all the DNA or Protein matrices. Type is either "protein"
# or "dna".
#
proc list_matrices {type} {
    global env
    if {[string match p* $type]} {
	set l [glob $env(EMBOSS_DATA)/EPAM* $env(EMBOSS_DATA)/EBLOSUM*]
    } else {
	set l [glob $env(EMBOSS_DATA)/ENUC*]
    }
    set r ""
    foreach f $l {
	if {[string match *~ $f]} {
	    continue
	}
	lappend r [file tail $f]
    }
    return $r
}

#
# Generates a list of all the Codon Usage tables.
#
proc list_codon_tables {} {
    global env
    set l [lsort -dictionary [glob $env(EMBOSS_DATA)/CODONS/*.cut]]
    set r ""
    foreach f $l {
	lappend r [file tail $f]
    }
    return $r
}

#
# Generates a list of all the emboss graphics output formats.
#
proc list_graph_types {} {
    return [list colourps data hp7470 hp7580 hpgl meta none png postscript \
	    tektronics tek4107t text xterm xwindows]
}

#
# Generates a list of all the emboss file formats
#
proc list_file_formats {} {
    return [list acedb asn1 clustal codata embl fasta fitch gcg genbank ig \
	    msf ncbi nbrf phylip staden strider swiss text]
}

# -----------------------------------------------------------------------------
# Trace callbacks, for tracing dialogue component modifications.

#
# Modification of a variable that vars($vname) depends on in namespace $nspace.
#
proc reset_value {nspace vname args} {
    upvar ${nspace}::vars vars
    if {$vars($vname) == $vars($vname.orig)} {
	eval set vars($vname) $vars($vname.expr)
	set vars($vname.orig) $vars($vname)
    }
}

# Modification of a variable used by a single-element list
proc reset_list {nspace vname args} {
    upvar ${nspace}::vars vars
    if {$vars($vname) != $vars($vname.orig)} {
	return
    }

    eval set newval $vars($vname.expr)
    array set a $vars($vname.mapping2)
    if {[info exists a($newval)]} {
	set vars($vname) $newval
	set vars($vname.orig) $vars($vname)
	set vars($vname.name) $a($newval)
    }
    unset a
}

# Modification of a variable used by a selection list
proc reset_select {nspace vname args} {
    upvar ${nspace}::vars vars
    if {$vars($vname) == $vars($vname.orig)} {
	eval set t $vars($vname.expr)
	$vars($vname.path) delete entry 0 end
	$vars($vname.path) insert entry end $t
	set vars($vname.orig) $vars($vname)
    }
}

# A different sequence name has been picked
proc sequence_changed {nspace name args} {
    upvar ${nspace}::vars vars

    set vars($name) [name_to_seq_id $vars($name.name)]
    if {$vars($name) != -1} {
	set vars($name.begin)   [seq_info $vars($name) start]
	set vars($name.end)     [seq_info $vars($name) end]
	set vars($name.length)  [expr {$vars($name.end)-$vars($name.begin)+1}]
	set vars($name.protein) [expr {[seq_info $vars($name) type]==1?0:1}]
	set vars(acdprotein)    $vars($name.protein)
	set vars($name.nucleotide) [expr {[seq_info $vars($name) type]==1?1:0}]
    }
}

# The 'needed' flag has changed
proc reset_needed {nspace vname args} {
    upvar ${nspace}::vars vars

    # puts exp=$vars($vname.needed_expr)
    if {[subst $vars($vname.needed_expr)] != 0} {
	catch {$vars($vname.path) configure -state normal} err
    } else {
	catch {$vars($vname.path) configure -state disabled} err
    }
}

# Modification of a variable that the 'info' attribute depends on
proc reset_name {nspace vname args} {
    upvar ${nspace}::vars vars
    set t [subst $vars($vname.info.expr)]
    if {[catch {$vars($vname.path) configure -labeltext $t} err1]} {
	if {[catch {$vars($vname.path) configure -label $t} err2]} {
	    if {[catch {$vars($vname.path) configure -text $t} err3]} {
		puts "Failed to configure $vname.path -labeltext $t"
		puts "    1. $err1"
		puts "    2. $err2"
		puts "    3. $err3"
	    }
	}
    }
}


# A 'list' selection has been modified
proc list_changed {nspace vname args} {
    upvar ${nspace}::vars vars
    set newsel $vars($vname.name)
    array set a $vars($vname.mapping1)
    if {[info exists a($newsel)]} {
	set vars($vname) $a($newsel)
    }
}

# A multiple list box has been modified
proc list_multi_changed {nspace vname args} {
    upvar ${nspace}::vars vars
    $vars($vname.path) selection clear 0 end
    set list [$vars($vname.path) get 0 end]
    array set a $vars($vname.mapping2)
    foreach item [split $vars($vname) ", "] {
	if {$item == ""} {
	    continue
	}
	catch {$vars($vname.path) selection set [lsearch $list $a($item)]} err
    }
    unset a
}

# A multiple selection box has been modified
proc select_multi_changed {nspace vname args} {
    upvar ${nspace}::vars vars
    $vars($vname.path) selection clear 0 end
    set list [$vars($vname.path) get 0 end]
    foreach item [split $vars($vname) ", "] {
	if {$item == ""} {
	    continue
	}
	catch {$vars($vname.path) selection set [lsearch -glob $list $item*]} err
    }
}

# A scrolledlistbox -selectioncommand callback - updates the variable
proc listbox_selected {nspace vname args} {
    upvar ${nspace}::vars vars
    set sel [$vars($vname.path) curselection]
    array set a $vars($vname.mapping1)
    foreach i $sel {
	catch {lappend l $a([$vars($vname.path) get $i])}
    }
    unset a
    set vars($vname) $l
}

# A scrolledlistbox -selectioncommand callback - updates the variable
proc selection_selected {nspace vname args} {
    upvar ${nspace}::vars vars
    set sel [$vars($vname.path) curselection]
    foreach i $sel {
	lappend l [$vars($vname.path) get $i]
    }
    set vars($vname) $l
}

# -----------------------------------------------------------------------------
# Other callback functions

#
# Callback from xentry (ACD "outfile" type).
#
proc check_outfile {path value} {
    global $path.overwrite

    set $path.overwrite 1

    if {$value == "" || ([string compare $value "stdout"] == 0)} {
	#puts FAILED
	return 0
    }
    #puts "check_outfile $value [file exists $value]"
    set response [CheckSaveFile "$value"]

    if {$response == 0} {
	set $path.overwrite 0
    } 

    if {$response && [file exists "$value"]} {
	DeleteFile $value
    }
    return $response
}

proc check_infile {path value} {
    global $path.overwrite

    set $path.overwrite 1

    if {$value == "" || ([string compare $value "stdout"] == 0)} {
	#puts FAILED
	return 0
    }
    #puts "check_outfile $value [file exists $value]"
    return [FileExists "$value"]
}

#
# Callback from ACD dirlist type
#
proc check_directory {path value} {
    global $path.overwrite

    set $path.overwrite 1

    if {$value == ""} {
	#puts FAILED
	return 0
    }

    # FIXME: No checking for now
    return 1
}

# Called whenever an update to the sequence is generated. For efficiencies
# sake this only sets the variables when they are really different (as this
# can be called many many times and may have variable traces set on the
# elements of var). We do not handle the sequence name change here as that's
# done separately by a variable trace of vars($name.name).
proc seq_updates {nspace name w} {
    upvar ${nspace}::vars vars
     
    set begin  [seq_id_from $w]
    set end    [seq_id_to   $w]
    set length [expr {$end-$begin+1}]
    if {$vars($name.length) != $length} {set vars($name.length) $length}
    if {$vars($name.begin)  != $begin}  {set vars($name.begin)  $begin}
    if {$vars($name.end)    != $end}    {set vars($name.end)    $end}
}

# -----------------------------------------------------------------------------
# Plotting functions

# Plots graphical output.
# There doesn't appear to be an easy way of identifying the graphical output
# produced except to search for $progname[0-9].dat files. This is clumsy as
# they could be old files, so we need to remove them before running the
# program.
#
# Arguments:
#	vname		Name and namespace of vars array 
#
# Returns
#       A list of output files that were not found to contain graphical data
proc plot_emboss1 {vname} {
    upvar $vname vars
    set outfiles ""

    set program $vars(application)

    # Find sequence and outfile types
    set nseqs 0
    set files ""
    foreach n [array names vars *._type] {
	if {[lsearch {sequence seqall seqset} $vars($n)] != -1} {
	    regsub {\._type$} $n {} n
	    set seq_id($nseqs)    [name_to_seq_id $vars($n.name)]
	    set seq_start($nseqs) $vars($n.begin)
	    set seq_end($nseqs)   $vars($n.end)
	    incr nseqs
	} elseif {[lsearch -exact "outfile featout report align" $vars($n)] != -1} {
	    regsub {\._type$} $n {} n
	    lappend files $vars($n)
	}
    }

    #find any data files
    if {![catch {set gfiles [glob $vars(application)*.dat]}]} {
	set files [concat $files $gfiles]
    }
    if {$files == ""} {
	# No output files
	return
    }

    #set up Raster class bindings
    set x_format %d
    set y_format %6f
    RasterBindings $x_format $y_format
    RasterGlobals

    # Iterate around all data files looking for ones that contain graphics.
    # Plot those that do, and return a list of those that do not.
    set id_count 0
    foreach f [lsort -ascii $files] {
	if {$f == "stdout"} continue

	# What type of data file is it? Parse first line.
	set fd [open $f r]
	set line [string tolower [gets $fd]]
	#puts "File=$f,line='$line'"

	set graphical 0
	if {[string match -nocase "##graphic" $line]} {
	    ##graphic file
	    plot_emboss_graphic $program $fd
	    close $fd
	    set graphical 1
	} elseif {[string match -nocase "##*2d plot*" $line]} {
	    ##2d plot / multi 2d plot / overlay 2d plot
	    close $fd
	    set graphical 1
	    switch $nseqs {
		1 {
		    set result_id($id_count) \
			[emboss create \
			     -seq_id_h $seq_id(0) \
			     -start_h $seq_start(0) \
			     -end_h $seq_end(0) \
			     -graph 0 \
			     -data $f]
		    if {$result_id($id_count) != -1} {
			plot_emboss $seq_id(0) $result_id($id_count) $program
		    }
		}
		2 {
		    set result_id($id_count) \
			[emboss create \
			     -seq_id_h $seq_id(0) \
			     -start_h $seq_start(0) \
			     -end_h $seq_end(0) \
			     -seq_id_v $seq_id(1) \
			     -start_v $seq_start(1) \
			     -end_v $seq_end(1) \
			     -graph 1 \
			     -data $f]
		    if {$result_id($id_count) != -1} {
			plot_emboss_dot $seq_id(0) $seq_id(1) \
				$result_id($id_count) $program
		    }
		}
		0 {
		    # Do nothing
		}
		default {
		    puts "FIXME: Unsupported number of sequences $nseqs"
		}
	    }

	    incr id_count
	}

	if {!$graphical} {
	    lappend outfiles $f
	}
    }

    return $outfiles
}

# -----------------------------------------------------------------------------
# Functions for handling the OK, Cancel, Help buttons.
#
# These are all running in the ::EMBOSS namespace, but are given the namespace
# of the dialogue to operate on.

# Tidy up function; a callback when the window is destroyed.
# This removes any namespace-global variables. This will also have the effect
# of removing the (possibly many) variable traces.
proc destroy_dialogue {nspace} {
    namespace eval $nspace {     
	variable vars
	variable arguments
	unset vars
	unset arguments
    }
}

#
# Called when a dialogue is created. Here we setup global variables expected
# by the ACD code.
# Currently the only one known is $(acdprotein).
proc init_dialogue {nspace} {
    upvar ${nspace}::vars vars
    set vars(acdprotein) 0
}

#
# Called when the OK button is pressed
# This actually runs the associated EMBOSS program, assuming that it manages
# to query all the values correctly.
#
proc run_dialogue {nspace w} {
    upvar ${nspace}::vars vars
    upvar ${nspace}::arguments arguments

    set outfile [tmpnam]-o
    set stderrfile [tmpnam]-1
    set stdoutfile [tmpnam]-2
    set tmpfiles [list $outfile $stderrfile $stdoutfile]
    set seqfiles ""

    foreach arg $arguments {
	# Disabled arguments can be skipped
	set d 0
	catch {set d [expr {[$vars($arg.path) cget -state] == "disabled"}]}
	if {$d} {
	    foreach i [array names vars $arg.*] {
		unset vars($i)
	    }
	    continue
	}

	# Otherwise parse and check
	switch $vars($arg._type) {
	    seqset -
	    seqall -
	    sequence {
		set type [string tolower $vars($arg.type)]
		#puts type=$type
		#puts $arg.protein=$vars($arg.protein)
	        if {($vars($arg.protein) == 0 && \
		    [regexp {protein$} $type] == 1) || \
		    ($vars($arg.protein) == 1 && \
		    [regexp {(protein$)|(.*any$)} $type] == 0)} {
		    tk_messageBox \
			-parent $w \
			-icon warning \
			-message "Invalid sequence type; should be $type"
		    return ""
		}

		# Save sequence to a temporary file
		set name $vars($arg.name)
		if {[CheckSequenceName $name] == 1} {
		    return ""
		}
		set id [name_to_seq_id $name]
		regsub -all {[/\\:]} $name _ name
		set name [tmpnam]
		lappend tmpfiles $name
		seq_file_save -seq_id $id -file $name -format 2
		set val $name
		lappend args -$arg $val -sbegin $vars($arg.begin) -send $vars($arg.end)
	    }
	    int -
	    integer {
		set val $vars($arg)
		if {([info exists vars($arg.minimum)] && $val < $vars($arg.minimum))||\
		    ([info exists vars($arg.maximum)] && $val > $vars($arg.maximum))} {
		    bell
		    tk_messageBox \
			-parent $w \
			-icon warning \
			-message "Value '$val' outside of legal range for \
				  setting '$arg'"
		    return ""
		}
		lappend args -$arg $val
	    }
	    seqout -
	    seqoutseq -
	    seqoutall {
		set format $vars($arg.format)
		set val $vars($arg)
		if {$format != ""} {
		    set val ${format}:$val
		}
		lappend args -$arg $val
		lappend seqfiles $val
	    }
	    featout -
	    report -
	    align -
	    outfile {
		if {$vars($arg) == ""} {
		    set vars($arg) $outfile
		}
		set val $vars($arg)
		lappend args -$arg $val
	    }
	    bool {
		set val $vars($arg)
		if {$val == "Y" || $val == "1"} {
		    lappend args -$arg
		} else {
		    lappend args -no$arg
		}
	    }
	    list {
		lappend args -$arg [$vars($arg.path) get]
	    }
	    list_multi {
		set val ""
		array set a $vars($arg.mapping1)
		foreach index [$vars($arg.path) curselection] {
		    lappend val $a([$vars($arg.path) get $index])
		}
		unset a
		lappend args -$arg [join $val $vars($arg.delimiter)]
	    }
	    selection_multi {
		set val ""
		foreach index [$vars($arg.path) curselection] {
		    lappend val [$vars($arg.path) get $index]
		}
		lappend args -$arg $val
	    }
	    default {
		set val $vars($arg)
		lappend args -$arg $val
	    }
	}
	if {$val == "" && $vars($arg.required) == 1} {
	    bell
	    tk_messageBox \
		-parent $w\
		-icon warning \
		-message "Invalid setting for '$arg'"
	    return ""
	}
    }

    # Remove any old output files
    set f ""
    catch {set f [glob $vars(application)*.dat]}
    foreach file $f {
	catch {file delete $file} err
    }

    # Run the program
    #puts "$vars(application) $args"
    if {[catch {eval exec $vars(application) $args >$stdoutfile 2>$stderrfile} err]} {
	verror ERR_WARN "EMBOSS:$vars(application) \"$err\""
	set fd [open $stderrfile]
	set err [read $fd]
	verror ERR_WARN $err
	close $fd
	tk_messageBox -parent $w -icon error \
	    -message "ERROR PRODUCED BY EMBOSS:\n\n$err"
	return ""
    }
    
    # Plot graphics and list text output
    vfuncheader "Output from EMBOSS '$vars(application)' tool"
    set fd [open $stdoutfile]
    vmessage [read $fd]
    close $fd
    foreach file [plot_emboss1 ${nspace}::vars] {
	#puts "display file $file"
	set fd [open $file]
	vmessage [read $fd]
	close $fd
    }

    # Load any newly generated sequences
    #puts seqfiles=$seqfiles
    foreach file $seqfiles {
	#puts file=$file
	if {[regexp {(?:([^:]*):)?(.*)$} $file dummy format name] == 0} {
	    continue
	}
	if {[file size $name] > 0} {
	    catch {read_sequence -file $name} err
	} else {
	    verror ERR_WARN $vars(application) \
	       "Program produced zero length sequence '$name'"
	}
    }

    foreach file $tmpfiles {
	file delete $file
    }

    destroy $w
}

# Returns the directory holding the data files.
# Uses embossdata to get the EMBOSS_DATA directory.
proc data_dir {} {
    global env
    if {[info exists env(EMBOSS_DATA)]} {
	return $env(EMBOSS_DATA)
    }

    set data_dir ""
    if {[catch {
	set fd [open "|embossdata . 2>/dev/null </dev/null"]
	set out [read $fd]
	close $fd
    }]} {
    	# verror ERR_WARN EMBOSS data directory not found
	return ""; #Unknown
    }
    foreach line [split $out "\n"] {
	if {[regexp {^([^\#].*?) *Exists$} $line dummy dir]} {
	    if {![file exists $dir/EPAM250]} {
		continue
	    }
	    set data_dir $dir
	}
    }
    if {$data_dir == ""} {
	verror ERR_WARN EMBOSS data directory not found
    }
    return $data_dir
}

# Returns the directory holding the ACD files.
#
# Uses embossdata to get the EMBOSS_DATA directory. Assumes the acd dir
# is a sibling of this.
proc acd_dir {} {
    set acd_dir ""
    set fd [open "|embossdata . 2>/dev/null </dev/null"]
    set out [read $fd]
    close $fd
    foreach line [split $out "\n"] {
	if {[regexp {^([^\#].*?) *Exists$} $line dummy dir]} {
	    set acd [file dirname $dir]/acd
	    if {![file exists $acd/wossname.acd]} {
		continue
	    }
	    set acd_dir $acd
	}
    }
    return $acd_dir
}

# Converts y/n, t/f into 1/0
proc convert_bool {text} {
    if {[string match {[YyTt]*} $text]} {
	return 1
    } elseif {[string match {[NnFf]*} $text]} {
	return 0
    } else {
	return $text
    }
}


namespace export *

}; # end namespace ::EMBOSS eval

#
# Initialises the emboss package
#
proc emboss_init {} {
    global tcl_platform auto_path env

    font create tabset_font -family Helvetica -weight bold -size -12
    font create label_font -family Helvetica -weight normal -size -12

    option add *Tabnotebook.font		tabset_font
    option add *Labeledframe.labelFont	label_font
    option add *Labeledframe.labelPos	nw

    if {$tcl_platform(platform) == "unix"} {
	button .tmp
	set bgl [tk::Darken [.tmp cget -background] 115]
	option add *textBackground $bgl widgetDefault
	destroy .tmp
    }

    lappend auto_path $env(STADTCL)/spin_emboss/acdtcl

    if {![info exists env(EMBOSS_DATA)]} {
	set dir [EMBOSS::data_dir]
	if {$dir != ""} {
	    set env(EMBOSS_DATA) $dir
	}
    }
  
    if {![info exists env(EMBOSS_DATA)] || $env(EMBOSS_DATA) == ""} {
        set ::EMBOSS::init 0
    } else {
        set ::EMBOSS::init 4
	regsub -all {\\} $env(EMBOSS_DATA) / env(EMBOSS_DATA)
    }
}

