package provide Xmclistbox 1.0

# Create this to make sure there are registered in auto_mkindex
# these must come before the [widget create ...]
proc Xmclistbox args {}
proc xmclistbox args {}

widget create Xmclistbox -type frame -base text -components {
    {base	text	text	{-font {Courier 12}}}
    {scrollbar  yscroll yscroll {}}
} -options {
    {-column_pad	columnPad	ColumnPad	{0}}
}

namespace eval ::Widget::Xmclistbox {;

;proc construct {w} {
    variable $w
    upvar 0 $w data

    # Default configurations
    set data(num_columns) 0
    set data(nItems) 0
    set data(pressed) 0

    # Display and configure text window
    grid rowconfigure $w 0 -weight 1
    grid columnconfigure $w 0 -weight 1
    grid $data(text) -row 0 -column 0 -sticky nsew

    $data(text) configure -yscrollcommand "$data(yscroll) set"
    $data(text) configure -wrap none

    bindtags $data(text) [list $data(text) Xmclistbox . all]

    $data(text) tag configure rowbox \
	-border 1 \
	-relief solid

    $data(text) tag configure enabled \
	-foreground [$data(text) cget -foreground]

    button $w._tmp
    $data(text) tag configure disabled \
	-foreground [$w._tmp cget -disabledforeground]
    destroy $w._tmp

#    $data(text) tag configure selected \
#	-background blue

    bind Xmclistbox <ButtonPress-1> \
	"[namespace current]::select_event %W \[%W index @%x,%y\] press"

    bind Xmclistbox <Any-Motion> \
	"[namespace current]::select_event %W \[%W index @%x,%y\] motion"

    bind Xmclistbox <ButtonRelease-1> \
	"[namespace current]::select_event %W \[%W index @%x,%y\] release"

    bind Xmclistbox <Any-Configure> \
	"if {\"%W\" == \"$w\"} {after idle {[namespace current]::config_event %W}}"

    # Configure vertical scrollbar
    $data(yscroll) configure -command "$data(text) yview"
}

;proc configure {w args} {
    variable $w
    upvar 0 $w data

    set truth {^(1|yes|true|on)$}

    foreach {key val} $args {
	set data($key) $val
    }
}

;proc _insert {w pos args} {
    variable $w
    upvar 0 $w data

    if {$pos != "end"} {
	set pos $pos.0
    }

    set new_tabs 0

    # When computing the width we append an "M". The reason for this is
    # that there is a minimum size for a tab in the text window, and using
    # font measure on \t (or a physical tab) actually computes the size of
    # backslash and "t". Also see _columnconfigure.
    foreach {item tags} $args {
	set targs ""
	for {set col 0} {$col < $data(num_columns)} {incr col} {
	    set i [lindex $item $col]
	    lappend targs "$i\t" "[list col_$col items] $tags"
	    set width [font measure $data(col_${col}_font) "${i}M"]
	    if {$data(col_${col}_width) < $width} {
		set data(col_${col}_width) $width
		set new_tabs 1
	    }
	}
	incr data(nItems)
	regsub "\t$" $targs {} targs
	eval $data(text) insert $pos $targs {\n}
    }

    if {$new_tabs} {
	recompute_tabs $w
    }
}

;proc _delete {w first {last ""}} {
    variable $w
    upvar 0 $w data

    if {$first == ""} {
	return
    } elseif {$first == "end"} {
	set first $data(nItems)
    }

    if {$last == ""} {
	set last $first
    } elseif {$last == "end"} {
	set last $data(nItems)
    }
    
    incr data(nItems) [expr {-($last-$first+1)}]

    return [eval $data(text) delete $first.0 $last.0+1l]
}

;proc _columnconfigure {w col args} {
    variable $w
    upvar 0 $w data

    eval $data(text) tag configure col_$col $args

    set fn [$data(text) tag cget col_$col -font]
    if {$fn == ""} {
	set fn [$data(text) cget -font]
    }
    set data(col_${col}_font) $fn  


    # Check whether this configure resizes the font
    set data(col_${col}_width) 1
    set width 0
    for {set i 0} {$i < $data(nItems)} {incr i} {
	set t [lindex [split [$data(text) get $i.0 $i.end] \t] $col]
	set width [font measure $data(col_${col}_font) "${t}M"]
	if {$data(col_${col}_width) < $width} {
	    set data(col_${col}_width) $width
	}
    }
    recompute_tabs $w
}

;proc _columnadd {w args} {
    variable $w
    upvar 0 $w data

    set c $data(num_columns)
    eval $data(text) tag configure col_$c $args
    set data(col_${c}_width) 1
    set fn [$data(text) tag cget col_$c -font]
    if {$fn == ""} {
	set fn [$data(text) cget -font]
    }
    set data(col_${c}_font) $fn  

    incr data(num_columns)

    recompute_tabs $w
}

;proc _curselection {w args} {
    variable $w
    upvar 0 $w data

    set ranges [$data(text) tag ranges sel]
    set lines ""
    foreach {start stop} $ranges {
	set start [lindex [split $start .] 0]
	set stop [lindex [split $start .] 0]
	for {set line $start} {$line <= $stop} {incr line} {
	    lappend lines $line
	}
    }
    return $lines
}

;proc _get {w first {last ""}} {
    variable $w
    upvar 0 $w data

    if {$first == ""} {
	return
    }

    if {$last == ""} {
	set last $first
    }

    return [$data(text) get $first.0 $last.end]
}

;proc _row_state {w row state} {
    variable $w
    upvar 0 $w data
    if {$state == 1} {
	set state enabled
    } elseif {$state == 0} {
	set state disabled
    }

    $data(text) tag remove enabled $row.0 $row.end 
    $data(text) tag remove disabled $row.0 $row.end 
    $data(text) tag add $state $row.0 $row.end 
}

;proc _bindcol {w column event script} {
    variable $w
    upvar 0 $w data

    incr column -1
    $data(text) tag bind col_$column $event \
	"after idle {uplevel #0 [list $script \[$data(text) curselection\]]}"
}

;proc _update {w row column text} {
    variable $w
    upvar 0 $w data

    # Get old text line
    set old_text [_get $w $row]
    
    # Edit it
    incr column -1
    set new_text [lreplace [split $old_text \t] $column $column $text]

    # Find tag names for this row
    set tags [$data(text) tag names $row.0]
    regsub -all {col_0|items} $tags {} tags
    
    # Replace with edited line
    _delete $w $row
    _insert $w $row $new_text $tags
}

;proc _select {w row} {
    variable $w
    upvar 0 $w data

    $data(text) tag remove sel 1.0 end
    $data(text) tag add sel $row.0 $row.end
    $data(text) tag remove rowbox 1.0 end
    $data(text) tag add rowbox $row.0 $row.end
}

;proc select_event {w el mode} {
    set w [winfo parent $w]

    variable $w
    upvar 0 $w data

    set row [lindex [split $el .] 0]

    switch $mode {
	"press" {
	    $data(text) tag remove rowbox 1.0 end
	    $data(text) tag add rowbox $row.0 $row.end

	    set data(pressed) 1
	}
	"motion" {
	    if {[info exists data(pressed)] && $data(pressed) != 0} {
		$data(text) tag remove rowbox 1.0 end
		$data(text) tag add rowbox $row.0 $row.end
	    }
	}
	"release" {
	    _select $w $row

	    set data(pressed) 0
	}
    }
}

;proc config_event {w} {
    variable $w
    upvar 0 $w data

    if {[llength [$data(text) bbox 1.0]] == 4 && \
	    [llength [$data(text) bbox [expr {$data(nItems)+1}].end]] == 4} {
	grid forget $data(yscroll)
    } else {
	grid $data(yscroll) -row 0 -column 1 -sticky nsew
    }
}

;proc recompute_tabs {w} {
    variable $w
    upvar 0 $w data

    # Reset the tab settings within the text window
    for {set tot 0; set col 0} {$col < $data(num_columns)} {incr col} {
	if {$col < [expr {$data(num_columns)-1}]} {
	    set tab [expr {$tot+$data(col_${col}_width)+$data(-column_pad)}]
	} else {
	    set tab [expr {$tot+$data(col_${col}_width)}]
	}
	lappend tabs $tab
	set tot $tab
    }
    incr tot -$data(-column_pad)
    $data(text) configure -tabs $tabs

#    incr tot [$data(yscroll) cget -width]

    # Compute new size of text window - tricky with varying fonts.
    # Use the base text font, which we have set to be fixed width.
    set fwid [font measure [$data(text) cget -font] a]
    set len [expr {int(ceil(double($tot) / $fwid))}]

    $data(text) configure -width $len
}

}; # end namespace eval ::Widget::Xmclistbox

