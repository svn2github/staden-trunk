# ---------------------------------------------------------------------------
# A rewrite of the old incrtcl comborange.itk widget to use Jeffrey Hobbs'
# pure tcl widget.tcl package instead.


# Create this to make sure there are registered in auto_mkindex
# these must come before the [widget create ...]
proc Xcomborange args {}
proc xcomborange args {}

widget create Xcomborange -type frame -base entry -components {
    {labelframe labelframe f {-text {Seq ID} -padx 5 -pady 5}}
    {xcombobox xcombo f.xcombo {\
        -command [list ::Widget::Xcomborange::entry_changed $w] \
	-label {Seq identifier} \
        -textvariable ${w}(tvar) \
    }}
    {xtwinspin xtwinspin f.xtwinspin {}}
    {frame feature f.feature {}}
    {listbox listbox f.feature.lbox { \
        -yscrollcommand [list $data(scrollbar) set] \
	-selectmode extended \
        -exportselection 0 \
    }}
    {scrollbar scrollbar f.feature.sy {-orient vertical \
				      -command [list $data(listbox) yview]}}
} -options {
    {-textvariable	textvariable	Textvariable	{}}
    {-label		label		Label		{Seq ID}}
    {-label_from	label		Label		{}}
    {-label_to		label		Label		{}}
    {-label_single	label		Label		{}}
    {-from		from		From		1}
    {-to                to              To              1000}
    {-state		state		State           normal}
    {-single            single          Single          0}
    {-trange		trange		Trange		1}
    {-feature		feature		Feature		0}
    {-default		default		Default		{}}
    {-labeltext		-label}
    {-start_value	-from}
    {-end_value		-to}
}

namespace eval ::Widget::Xcomborange {

;proc construct {w} {
    variable $w
    upvar 0 $w data

    pack $data(labelframe) -side top -fill both -expand 1 -padx 5 -pady 5

    if {$data(-textvariable) == ""} {
	set data(-textvariable) ${w}(tvar)
    }

    pack $data(xcombo) -side top -fill both
    if {$data(-trange)} {
	pack $data(xtwinspin) -side top -fill both -pady 5
    }
    pack $data(listbox) -side left -fill both -expand 1
    pack $data(scrollbar) -side right -fill both
    if {$data(-feature)} {
	pack $data(feature) -side top -fill both
    }

    # Populate the comobobox with a list of sequence IDs
    set name_list {}
    set data(sequences) [sequence_names]
    foreach i $data(sequences) {
	lappend name_list [lindex $i 1]
    }

    $data(xcombo) configure -values $name_list

    # Set the default starting point
    foreach i [sequence_names] {
	if {[lindex $i 0] == "H"} {
	    set_seq $w $i
	    set seq [lindex $i 1]
	    set pos [lindex $i 2]
	    foreach {_ start end} [regexp -inline {(\d+)..(\d+)} $pos] break

	    $data(xcombo) set $seq
	    $data(xtwinspin) set $start $end
	    $data(xtwinspin) configure -from 1 -to [lindex $i 3]
	}
    }
}

;proc configure {w args} {
    variable $w
    upvar 0 $w data

    foreach {key val} $args {
	set data($key) $val
	switch -- $key {
	    -label {
		$data(labelframe) configure -text $val
	    }
	    -textvariable {
		$data(xcombo) configure -textvariable $val
		if {[info exists data(trace)]} {
		    eval trace remove variable $data(trace)
		}

		set data(trace) [list ::$val write \
				    "[namespace current]::entry_changed $w"]
		eval trace add variable $data(trace)
	    }
	    -trange {
		if {$val} {
		    catch {pack $data(xtwinspin) -side top -fill both}
		} else {
		    catch {pack forget $data(xtwinspin)}
		}
	    }
	    -feature {
		if {$val} {
		    catch {pack $data(feature) -side top -fill both}
		} else {
		    catch {pack forget $data(feature)}
		}
	    }
	    -from {
		$data(xtwinspin) configure -from $val
	    }
	    -to {
		$data(xtwinspin) configure -to $val
	    }
	    -single {
		$data(xtwinspin) configure -single $val
	    }
	    -label_single -
	    -label_from {
		$data(xtwinspin) configure -label1 $val
	    }
	    -label_to {
		$data(xtwinspin) configure -label2 $val
	    }
	    -default {
		$data(xcombo) set $val
	    }
	}
    }
}

;proc entry_changed {w args} {
    variable $w
    upvar 0 $w data

    set text [set ::$data(-textvariable)]

    foreach i $data(sequences) {
	if {[string match [lindex $i 1] $text]} {
	    set_seq $w $i
	}
    }
}

;proc set_seq {w i} {
    variable $w
    upvar 0 $w data

    # Entry window
    set range [lindex $i 2]
    foreach {_ start end} [regexp -inline {(\d+)..(\d+)} $range] break

    # twin range
    $data(xtwinspin) set $start $end
    $data(xtwinspin) configure -from 1 -to [lindex $i 3]

    # Feature list
    set seq_id [name_to_seq_id [lindex $i 1]]
    set seq_ncds [seq_info $seq_id numbercds]

    $data(listbox) delete 0 end
    for {set i 1} {$i <= $seq_ncds} {incr i} {
	$data(listbox) insert end [seq_info $seq_id key_index_cds $i]
    }
}

;proc _get_seqname {w} {
    variable $w
    upvar 0 $w data

    set seq_name [set ::$data(-textvariable)]
    set seq_id [name_to_seq_id $seq_name]

    if {$seq_id == -1} {
	return ""
    } 
    return $seq_name
}

;proc _get_s {w} {
    variable $w
    upvar 0 $w data

    return [$data(xtwinspin) get_start]
}

;proc _get_e {w} {
    variable $w
    upvar 0 $w data

    return [$data(xtwinspin) get_end]
}

;proc _get_single {w} {
    variable $w
    upvar 0 $w data

    return [$data(xtwinspin) get_start]
}

# Returns any selection from the listbox element
;proc _get_curselection {w} {
    variable $w
    upvar 0 $w data

    set l {}
    foreach pos [$data(listbox) curselection] {
	lappend l [$data(listbox) get $pos]
    }

    return $l
}


}; # matches namespace eval
