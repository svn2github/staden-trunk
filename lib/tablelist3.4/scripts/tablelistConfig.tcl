#==============================================================================
# Contains private configuration procedures for tablelist widgets.
#
# Copyright (c) 2000-2004  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

#------------------------------------------------------------------------------
# tablelist::extendConfigSpecs
#
# Extends the elements of the array configSpecs.
#------------------------------------------------------------------------------
proc tablelist::extendConfigSpecs {} {
    variable helpLabel
    variable configSpecs
    variable winSys

    #
    # Append the default values of the configuration options
    # of a temporary, invisible listbox widget to the values
    # of the corresponding elements of the array configSpecs
    #
    set helpListbox .__helpListbox
    for {set n 0} {[winfo exists $helpListbox]} {incr n} {
	set helpListbox .__helpListbox$n
    }
    listbox $helpListbox
    foreach configSet [$helpListbox configure] {
	if {[llength $configSet] != 2} {
	    set opt [lindex $configSet 0]
	    if {[info exists configSpecs($opt)]} {
		lappend configSpecs($opt) [lindex $configSet 3]
	    }
	}
    }
    destroy $helpListbox

    #
    # Append the default values of some configuration options
    # of an invisible label widget to the values of the
    # corresponding -label* elements of the array configSpecs
    #
    set helpLabel .__helpLabel
    for {set n 0} {[winfo exists $helpLabel]} {incr n} {
	set helpLabel .__helpLabel$n
    }
    label $helpLabel
    foreach optTail {font height} {
	set configSet [$helpLabel config -$optTail]
	lappend configSpecs(-label$optTail) [lindex $configSet 3]
    }
    if {[catch {$helpLabel config -state disabled}] == 0 &&
	[catch {$helpLabel config -state normal}] == 0 &&
	[catch {$helpLabel config -disabledforeground} configSet] == 0} {
	lappend configSpecs(-labeldisabledforeground) [lindex $configSet 3]
    } else {
	unset configSpecs(-labeldisabledforeground)
    }
    if {[string compare $winSys win32] == 0} {
	lappend configSpecs(-labelpady) 0
    } else {
	set configSet [$helpLabel config -pady]
	lappend configSpecs(-labelpady) [lindex $configSet 3]
    }

    #
    # Steal the default values of some configuration
    # options from a temporary, invisible button widget
    #
    set helpButton .__helpButton
    for {set n 0} {[winfo exists $helpButton]} {incr n} {
	set helpButton .__helpButton$n
    }
    button $helpButton
    foreach opt {-disabledforeground -state} {
	set configSet [$helpButton config $opt]
	lappend configSpecs($opt) [lindex $configSet 3]
    }
    foreach optTail {background foreground} {
	set configSet [$helpButton config -$optTail]
	lappend configSpecs(-label$optTail) [lindex $configSet 3]
    }
    if {[string compare $winSys classic] == 0 ||
	[string compare $winSys aqua] == 0} {
	lappend configSpecs(-labelborderwidth) 1
    } else {
	set configSet [$helpButton config -borderwidth]
	lappend configSpecs(-labelborderwidth) [lindex $configSet 3]
    }
    destroy $helpButton

    #
    # Extend the remaining elements of the array configSpecs
    #
    lappend configSpecs(-activestyle)		underline
    lappend configSpecs(-arrowcolor)		{}
    lappend configSpecs(-arrowdisabledcolor)	{}
    lappend configSpecs(-columns)		{}
    lappend configSpecs(-editendcommand)	{}
    lappend configSpecs(-editstartcommand)	{}
    lappend configSpecs(-incrarrowtype)		up
    lappend configSpecs(-labelcommand)		{}
    lappend configSpecs(-labelrelief)		raised
    lappend configSpecs(-listvariable)		{}
    lappend configSpecs(-movablecolumns)	0
    lappend configSpecs(-movecolumncursor)	icon
    lappend configSpecs(-resizablecolumns)	1
    lappend configSpecs(-resizecursor)		sb_h_double_arrow
    lappend configSpecs(-showarrow)		1
    lappend configSpecs(-showlabels)		1
    lappend configSpecs(-showseparators)	0
    lappend configSpecs(-snipstring)		...
    lappend configSpecs(-sortcommand)		{}
    lappend configSpecs(-stretch)		{}
    lappend configSpecs(-stripebackground)	{}
    lappend configSpecs(-stripeforeground)	{}
    lappend configSpecs(-stripeheight)		1
    lappend configSpecs(-targetcolor)		black
}

#------------------------------------------------------------------------------
# tablelist::doConfig
#
# Applies the value val of the configuration option opt to the tablelist widget
# win.
#------------------------------------------------------------------------------
proc tablelist::doConfig {win opt val} {
    variable helpLabel
    variable configSpecs
    upvar ::tablelist::ns${win}::data data

    #
    # Apply the value to the widget(s) corresponding to the given option
    #
    switch [lindex $configSpecs($opt) 2] {
	c {
	    #
	    # Apply the value to all children and save the
	    # properly formatted value of val in data($opt)
	    #
	    foreach w [winfo children $win] {
		if {[regexp {^(body|hdr|sep[0-9]+)$} [winfo name $w]]} {
		    $w configure $opt $val
		}
	    }
	    $data(hdrTxt) configure $opt $val
	    $data(hdrLbl) configure $opt $val
	    foreach w [winfo children $data(hdrTxtFr)] {
		$w configure $opt $val
		foreach c [winfo children $w] {
		    $c configure $opt $val
		}
	    }
	    set data($opt) [$data(hdrLbl) cget $opt]
	}

	b {
	    #
	    # Apply the value to the body text widget and save
	    # the properly formatted value of val in data($opt)
	    #
	    set w $data(body)
	    $w configure $opt $val
	    set data($opt) [$w cget $opt]

	    switch -- $opt {
		-background {
		    #
		    # Apply the value to the frame (because of
		    # the shadow colors of its 3-D border), to
		    # the separators, and to the "disabled" tag
		    #
		    $win configure $opt $val
		    foreach c [winfo children $win] {
			if {[regexp {^sep[0-9]+$} [winfo name $c]]} {
			    $c configure $opt $val
			}
		    }
		    $w tag configure disabled $opt $val
		}
		-font {
		    #
		    # Apply the value to the listbox child, rebuild the lists
		    # of the column fonts and tag names, configure the edit
		    # window if present, set up and adjust the columns, and
		    # make sure the items will be redisplayed at idle time
		    #
		    $data(lb) configure $opt $val
		    set data(charWidth) [font measure $val -displayof $win 0]
		    makeColFontAndTagLists $win
		    if {$data(editRow) >= 0} {
			variable editWin
			set name $data($data(editCol)-editwindow)
			set item [lindex $data(itemList) $data(editRow)]
			set key [lindex $item end]
			catch {
			    $data(bodyFrEd) configure $editWin($name-fontOpt) \
				[cellFont $win $key $data(editCol)]
			}
		    }
		    for {set col 0} {$col < $data(colCount)} {incr col} {
			if {$data($col-maxwidth) > 0} {
			    set data($col-maxPixels) \
				[charsToPixels $win $val $data($col-maxwidth)]
			}
		    }
		    setupColumns $win $data(-columns) 0
		    adjustColumns $win allCols 1
		    redisplayWhenIdle $win
		}
		-foreground {
		    #
		    # Apply the value to the "disabled" tag if needed
		    #
		    if {[string compare $data(-disabledforeground) ""] == 0} {
			$w tag configure disabled $opt $val
		    }
		}
	    }
	}

	h {
	    #
	    # Apply the value to the header text widget and save
	    # the properly formatted value of val in data($opt)
	    #
	    set w $data(hdrTxt)
	    $w configure $opt $val
	    set data($opt) [$w cget $opt]
	}

	l {
	    #
	    # Apply the value to all not individually configured labels
	    # and save the properly formatted value of val in data($opt)
	    #
	    set optTail [string range $opt 6 end]	;# remove the -label
	    $data(hdrLbl) configure -$optTail $val
	    for {set col 0} {$col < $data(colCount)} {incr col} {
		set w $data(hdrTxtFrLbl)$col
		if {![info exists data($col$opt)]} {
		    configLabel $w -$optTail $val
		}
	    }
	    set data($opt) [$data(hdrLbl) cget -$optTail]

	    switch -- $opt {
		-labelbackground {
		    #
		    # Apply the value to the children of the header frame and
		    # conditionally both to the children of the labels (if
		    # any) and to the canvas displaying an up- or down-arrow
		    #
		    $data(hdrTxt) configure -$optTail $val
		    for {set col 0} {$col < $data(colCount)} {incr col} {
			set w $data(hdrTxtFrLbl)$col
			if {![info exists data($col$opt)]} {
			    foreach c [winfo children $w] {
				$c configure -$optTail $val
			    }
			}
		    }
		    if {$data(arrowCol) >= 0 &&
			![info exists data($data(arrowCol)$opt)]} {
			configCanvas $win
		    }
		}
		-labelborderwidth {
		    #
		    # Adjust the columns (including
		    # the height of the header frame)
		    #
		    adjustColumns $win allLabels 1
		}
		-labeldisabledforeground {
		    #
		    # Apply the value to the children of the labels (if any)
		    #
		    foreach w [winfo children $data(hdrTxtFr)] {
			foreach c [winfo children $w] {
			    $c configure $opt $val
			}
		    }
		}
		-labelfont {
		    #
		    # Conditionally apply the value to the children of
		    # the labels (if any), conditionally resize the canvas
		    # displaying an up- or down-arrow, and adjust the
		    # columns (including the height of the header frame)
		    #
		    for {set col 0} {$col < $data(colCount)} {incr col} {
			set w $data(hdrTxtFrLbl)$col
			if {![info exists data($col$opt)]} {
			    foreach c [winfo children $w] {
				$c configure -$optTail $val
			    }
			}
		    }
		    if {$data(arrowCol) >= 0 &&
			![info exists data($data(arrowCol)$opt)]} {
			configCanvas $win
			drawArrows $win
		    }
		    adjustColumns $win allLabels 1
		}
		-labelforeground {
		    #
		    # Conditionally apply the value to
		    # the children of the labels (if any)
		    #
		    for {set col 0} {$col < $data(colCount)} {incr col} {
			set w $data(hdrTxtFrLbl)$col
			if {![info exists data($col$opt)]} {
			    foreach c [winfo children $w] {
				$c configure -$optTail $val
			    }
			}
		    }
		}
		-labelheight -
		-labelpady {
		    #
		    # Adjust the height of the header frame
		    #
		    adjustHeaderHeight $win
		}
	    }
	}

	f {
	    #
	    # Apply the value to the frame and save the
	    # properly formatted value of val in data($opt)
	    #
	    $win configure $opt $val
	    set data($opt) [$win cget $opt]
	}

	w {
	    switch -- $opt {
		-activestyle {
		    #
		    # Configure the "active" tag and save the
		    # properly formatted value of val in data($opt)
		    #
		    variable activeStyles
		    set val [mwutil::fullOpt "active style" $val $activeStyles]
		    set w $data(body)
		    switch $val {
			frame {
			    $w tag configure active -relief solid -underline 0
			}
			none {
			    $w tag configure active -relief flat -underline 0
			}
			underline {
			    $w tag configure active -relief flat -underline 1
			}
		    }
		    set data($opt) $val
		}
		-arrowcolor {
		    #
		    # Set the color of the normal arrow and save the
		    # properly formatted value of val in data($opt)
		    #
		    set data($opt) [fillArrow $data(hdrTxtFrCanv) normal $val]
		}
		-arrowdisabledcolor {
		    #
		    # Set the color of the disabled arrow and save the
		    # properly formatted value of val in data($opt)
		    #
		    set data($opt) [fillArrow $data(hdrTxtFrCanv) disabled $val]
		}
		-columns {
		    #
		    # Set up and adjust the columns, rebuild
		    # the lists of the column fonts and tag
		    # names, and redisplay the items
		    #
		    setupColumns $win $val 1
		    adjustColumns $win allCols 1
		    makeColFontAndTagLists $win
		    redisplay $win
		}
		-disabledforeground {
		    #
		    # Configure the "disabled" tag in the body text widget and
		    # save the properly formatted value of val in data($opt)
		    #
		    set w $data(body)
		    if {[string compare $val ""] == 0} {
			$w tag configure disabled -fgstipple gray50 \
				-foreground $data(-foreground)
			set data($opt) ""
		    } else {
			$w tag configure disabled -fgstipple "" \
				-foreground $val
			set data($opt) [$w tag cget disabled -foreground]
		    }
		}
		-editendcommand -
		-editstartcommand -
		-labelcommand -
		-selectmode -
		-sortcommand {
		    set data($opt) $val
		}
		-exportselection {
		    #
		    # Save the boolean value specified by val in
		    # data($opt).  In addition, if the selection is
		    # exported and there are any selected rows in the
		    # widget then make win the new owner of the PRIMARY
		    # selection and register a callback to be invoked
		    # when it loses ownership of the PRIMARY selection
		    #
		    set data($opt) [expr {$val ? 1 : 0}]
		    if {$val &&
			[llength [$data(body) tag nextrange select 1.0]] != 0} {
			selection own -command \
				[list ::tablelist::lostSelection $win] $win
		    }
		}
		-height {
		    #
		    # Adjust the heights of the body text widget
		    # and of the listbox child, and save the
		    # properly formatted value of val in data($opt)
		    #
		    set val [format %d $val]	;# integer check with error msg
		    if {$val <= 0} {
			$data(body) configure $opt $data(itemCount)
			$data(lb) configure $opt $data(itemCount)
		    } else {
			$data(body) configure $opt $val
			$data(lb) configure $opt $val
		    }
		    set data($opt) $val
		}
		-incrarrowtype {
		    #
		    # Save the properly formatted value of val
		    # in data($opt) and draw the arrows if
		    # the canvas widget is presently mapped
		    #
		    variable arrowTypes
		    set data($opt) \
			[mwutil::fullOpt "arrow type" $val $arrowTypes]
		    if {$data(arrowCol) >= 0} {
			drawArrows $win
		    }
		}
		-listvariable {
		    #
		    # Associate val as list variable with the
		    # given widget and save it in data($opt)
		    #
		    makeListVar $win $val
		    set data($opt) $val
		    if {[string compare $val ""] == 0} {
			set data(hasListVar) 0
		    } else {
			set data(hasListVar) 1
		    }
		}
		-movablecolumns -
		-resizablecolumns {
		    #
		    # Save the boolean value specified by val in data($opt)
		    #
		    set data($opt) [expr {$val ? 1 : 0}]
		}
		-movecolumncursor -
		-resizecursor {
		    #
		    # Save the properly formatted value of val in data($opt)
		    #
		    $helpLabel configure -cursor $val
		    set data($opt) [$helpLabel cget -cursor]
		}
		-selectbackground -
		-selectforeground {
		    #
		    # Configure the "select" tag in the body text widget
		    # and save the properly formatted value of val in
		    # data($opt).  Don't use the built-in "sel" tag
		    # because on Windows the selection in a text widget only
		    # becomes visible when the window gets the input focus.
		    #
		    set w $data(body)
		    set optTail [string range $opt 7 end] ;# remove the -select
		    $w tag configure select -$optTail $val
		    set data($opt) [$w tag cget select -$optTail]
		}
		-selectborderwidth {
		    #
		    # Configure the "select" tag in the body text widget
		    # and save the properly formatted value of val in
		    # data($opt).  Don't use the built-in "sel" tag
		    # because on Windows the selection in a text widget only
		    # becomes visible when the window gets the input focus.
		    # In addition, adjust the line spacing accordingly and
		    # apply the value to the listbox child, too.
		    #
		    set w $data(body)
		    set optTail [string range $opt 7 end] ;# remove the -select
		    $w tag configure select -$optTail $val
		    set data($opt) [$w tag cget select -$optTail]
		    if {$val < 0} {
			set val 0
		    }
		    $w configure -spacing1 $val -spacing3 [expr {$val + 1}]
		    $data(lb) configure $opt $val
		}
		-setgrid {
		    #
		    # Apply the value to the listbox child and save
		    # the properly formatted value of val in data($opt)
		    #
		    $data(lb) configure $opt $val
		    set data($opt) [$data(lb) cget $opt]
		}
		-showarrow {
		    #
		    # Save the boolean value specified by val in
		    # data($opt) and conditionally unmanage the
		    # canvas displaying an up- or down-arrow
		    #
		    set data($opt) [expr {$val ? 1 : 0}]
		    if {!$data($opt) && $data(arrowCol) >= 0} {
			place forget $data(hdrTxtFrCanv)
			set oldArrowCol $data(arrowCol)
			set data(arrowCol) -1
			adjustColumns $win l$oldArrowCol 1
		    }
		}
		-showlabels {
		    #
		    # Save the boolean value specified by val in data($opt)
		    # and adjust the height of the header frame
		    #
		    set data($opt) [expr {$val ? 1 : 0}]
		    adjustHeaderHeight $win
		}
		-showseparators {
		    #
		    # Save the boolean value specified by val in data($opt),
		    # and create or destroy the separators if needed
		    #
		    set oldVal $data($opt)
		    set data($opt) [expr {$val ? 1 : 0}]
		    if {!$oldVal && $data($opt)} {
			createSeps $win
		    } elseif {$oldVal && !$data($opt)} {
			foreach w [winfo children $win] {
			    if {[regexp {^sep[0-9]+$} [winfo name $w]]} {
				destroy $w
			    }
			}
		    }
		}
		-snipstring {
		    #
		    # Save val in data($opt), adjust the columns, and make
		    # sure the items will be redisplayed at idle time
		    #
		    set data($opt) $val
		    adjustColumns $win {} 0
		    redisplayWhenIdle $win
		}
		-state {
		    #
		    # Apply the value to all labels and their children (if
		    # any), as well as to the edit window (if present),
		    # raise the corresponding arrow in the canvas, add/
		    # remove the "disabled" tag to/from the contents of
		    # the body text widget, configure the borderwidth
		    # of the "active" and "select" tags, and save the
		    # properly formatted value of val in data($opt)
		    #
		    variable states
		    set val [mwutil::fullOpt "state" $val $states]
		    catch {
			foreach w [winfo children $data(hdrTxtFr)] {
			    $w configure $opt $val
			    foreach c [winfo children $w] {
				$c configure $opt $val
			    }
			}
		    }
		    if {$data(editRow) >= 0} {
			catch {$data(bodyFrEd) configure $opt $val}
		    }
		    raiseArrow $data(hdrTxtFrCanv) $val
		    set w $data(body)
		    switch $val {
			disabled {
			    $w tag add disabled 1.0 end
			    $w tag configure active -borderwidth 0
			    $w tag configure select -borderwidth 0
			    set data(isDisabled) 1
			}
			normal {
			    $w tag remove disabled 1.0 end
			    $w tag configure active -borderwidth 1
			    $w tag configure select -borderwidth \
				    $data(-selectborderwidth)
			    set data(isDisabled) 0
			}
		    }
		    set data($opt) $val
		}
		-stretch {
		    #
		    # Save the properly formatted value of val in
		    # data($opt) and stretch the stretchable columns
		    #
		    if {[string first $val all] == 0} {
			set data($opt) all
		    } else {
			set data($opt) $val
			sortStretchableColList $win
		    }
		    set data(forceAdjust) 1
		    stretchColumnsWhenIdle $win
		}
		-stripebackground -
		-stripeforeground {
		    #
		    # Configure the "stripe" tag in the body text
		    # widget, save the properly formatted value of val
		    # in data($opt), and draw the stripes if necessary
		    #
		    set w $data(body)
		    set optTail [string range $opt 7 end] ;# remove the -stripe
		    $w tag configure stripe -$optTail $val
		    set data($opt) [$w tag cget stripe -$optTail]
		    makeStripesWhenIdle $win
		}
		-stripeheight {
		    #
		    # Save the properly formatted value of val in
		    # data($opt) and draw the stripes if necessary
		    #
		    set val [format %d $val]	;# integer check with error msg
		    set data($opt) $val
		    makeStripesWhenIdle $win
		}
		-targetcolor {
		    #
		    # Set the color of the header gap and save the
		    # properly formatted value of val in data($opt)
		    #
		    $data(hdrGap) configure -background $val
		    set data($opt) [$data(hdrGap) cget -background]
		}
		-width {
		    #
		    # Adjust the widths of the body text widget,
		    # header frame, and listbox child, and save the
		    # properly formatted value of val in data($opt)
		    #
		    set val [format %d $val]	;# integer check with error msg
		    $data(body) configure $opt $val
		    if {$val <= 0} {
			$data(hdr) configure $opt $data(hdrPixels)
			$data(lb) configure $opt \
				  [expr {$data(hdrPixels) / $data(charWidth)}]
		    } else {
			$data(hdr) configure $opt 0
			$data(lb) configure $opt $val
		    }
		    set data($opt) $val
		}
	    }
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::doColConfig
#
# Applies the value val of the column configuration option opt to the col'th
# column of the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::doColConfig {col win opt val} {
    upvar ::tablelist::ns${win}::data data

    switch -- $opt {
	-align {
	    #
	    # Set up and adjust the columns, and make sure
	    # the items will be redisplayed at idle time
	    #
	    set idx [expr {3*$col + 2}]
	    setupColumns $win [lreplace $data(-columns) $idx $idx $val] 0
	    adjustColumns $win {} 0
	    redisplayWhenIdle $win
	}

	-background -
	-foreground {
	    set w $data(body)
	    set tag $col$opt
	    if {[string compare $val ""] == 0} {
		if {[info exists data($tag)]} {
		    $w tag delete $tag
		    unset data($tag)
		}
	    } else {
		#
		# Configure the tag $tag in the body text widget
		#
		$w tag configure $tag $opt $val
		$w tag raise $tag stripe

		if {!$data($col-hide)} {
		    #
		    # Apply the tag to the elements of the given column
		    #
		    for {set idx 0; set line 1} {$idx < $data(itemCount)} \
			{incr idx; incr line} {
			if {[lsearch -exact [$w tag names $line.0] select]
			    < 0} {
			    findCellTabs $win $line $col tabIdx1 tabIdx2
			    $w tag add $tag $tabIdx1 $tabIdx2+1c
			}
		    }
		}

		#
		# Save the properly formatted value of val in data($tag)
		#
		set data($tag) [$w tag cget $tag $opt]
	    }

	    #
	    # Rebuild the lists of the column fonts and tag names
	    #
	    makeColFontAndTagLists $win
	}

	-editable -
	-resizable {
	    #
	    # Save the boolean value specified by val in data($col$opt)
	    #
	    set data($col$opt) [expr {$val ? 1 : 0}]
	}

	-editwindow {
	    variable editWin
	    if {[regexp {^(entry|checkbutton)$} $val] ||
		($::tk_version >= 8.4 && [string compare $val spinbox] == 0) ||
		[info exists editWin($val-creationCmd)]} {
		set data($col$opt) $val
	    } else {
		return -code error "name \"$val\" is not registered\
				    for interactive cell editing"
	    }
	}

	-font {
	    set w $data(body)
	    set tag $col$opt
	    if {[string compare $val ""] == 0} {
		if {[info exists data($tag)]} {
		    $w tag delete $tag
		    unset data($tag)
		}
	    } else {
		#
		# Configure the tag $tag in the body text widget
		#
		$w tag configure $tag $opt $val
		$w tag lower $tag

		if {!$data($col-hide)} {
		    #
		    # Apply the tag to the elements of the given column
		    #
		    for {set idx 0; set line 1} {$idx < $data(itemCount)} \
			{incr idx; incr line} {
			findCellTabs $win $line $col tabIdx1 tabIdx2
			$w tag add $tag $tabIdx1 $tabIdx2+1c
		    }
		}

		#
		# Save the properly formatted value of val in data($tag)
		#
		set data($tag) [$w tag cget $tag $opt]
	    }

	    #
	    # Rebuild the lists of the column fonts and tag names
	    #
	    makeColFontAndTagLists $win

	    #
	    # Adjust the columns, and make sure the items
	    # will be redisplayed at idle time if needed
	    #
	    adjustColumns $win $col 1
	    set pixels [lindex $data(colList) [expr {2*$col}]]
	    if {$pixels == 0} {			;# convention: dynamic width
		if {$data($col-maxPixels) > 0 &&
		    $data($col-reqPixels) > $data($col-maxPixels)} {
		    set pixels $data($col-maxPixels)
		}
	    }
	    if {$pixels != 0} {
		redisplayWhenIdle $win
	    }

	    if {$col == $data(editCol)} {
		#
		# Configure the edit window
		#
		variable editWin
		set name $data($col-editwindow)
		set item [lindex $data(itemList) $data(editRow)]
		set key [lindex $item end]
		catch {
		    $data(bodyFrEd) configure $editWin($name-fontOpt) \
			[cellFont $win $key $col]
		}
	    }
	}

	-formatcommand {
	    if {[string compare $val ""] == 0} {
		if {[info exists data($col$opt)]} {
		    unset data($col$opt)
		}
		set fmtCmdFlag 0
	    } else {
		set data($col$opt) $val
		set fmtCmdFlag 1
	    }

	    #
	    # Update the corresponding element of the list data(fmtCmdFlagList)
	    #
	    set data(fmtCmdFlagList) \
		[lreplace $data(fmtCmdFlagList) $col $col $fmtCmdFlag]

	    #
	    # Adjust the columns and make sure the
	    # items will be redisplayed at idle time
	    #
	    adjustColumns $win $col 1
	    redisplayWhenIdle $win
	}

	-hide {
	    #
	    # Save the boolean value specified by val in data($col$opt),
	    # adjust the columns, and redisplay the items
	    #
	    set oldVal $data($col$opt)
	    set data($col$opt) [expr {$val ? 1 : 0}]
	    adjustColumns $win $col 1
	    if {$data($col$opt) != $oldVal} {
		redisplay $win
	    }
	}

	-labelalign {
	    if {[string compare $val ""] == 0} {
		#
		# Unset data($col$opt)
		#
		set alignment [lindex $data(colList) [expr {2*$col + 1}]]
		if {[info exists data($col$opt)]} {
		    unset data($col$opt)
		}
	    } else {
		#
		# Save the properly formatted value of val in data($col$opt)
		#
		variable alignments
		set val [mwutil::fullOpt "label alignment" $val $alignments]
		set alignment $val
		set data($col$opt) $val
	    }

	    #
	    # Adjust the col'th label
	    #
	    set pixels [lindex $data(colList) [expr {2*$col}]]
	    if {$pixels == 0} {			;# convention: dynamic width
		if {$data($col-maxPixels) > 0 &&
		    $data($col-reqPixels) > $data($col-maxPixels)} {
		    set pixels $data($col-maxPixels)
		}
	    }
	    if {$pixels != 0} {	
		incr pixels $data($col-delta)
	    }
	    adjustLabel $win $col $pixels $alignment
	}

	-labelbackground {
	    set w $data(hdrTxtFrLbl)$col
	    set optTail [string range $opt 6 end]	;# remove the -label
	    if {[string compare $val ""] == 0} {
		#
		# Apply the value of the corresponding widget configuration
		# option to the col'th label and its children (if any)
		# and conditionally to the canvas displaying an up-
		# or down-arrow, and unset data($col$opt)
		#
		$w configure -$optTail $data($opt)
		foreach c [winfo children $w] {
		    $c configure -$optTail $data($opt)
		}
		if {$col == $data(arrowCol)} {
		    configCanvas $win
		}
		if {[info exists data($col$opt)]} {
		    unset data($col$opt)
		}
	    } else {
		#
		# Apply the given value to the col'th label and its
		# children (if any) and conditionally to the canvas
		# displaying an up- or down-arrow, and save the
		# properly formatted value of val in data($col$opt)
		#
		$w configure -$optTail $val
		foreach c [winfo children $w] {
		    $c configure -$optTail $val
		}
		if {$col == $data(arrowCol)} {
		    configCanvas $win
		}
		set data($col$opt) [$w cget -$optTail]
	    }
	}

	-labelborderwidth {
	    set w $data(hdrTxtFrLbl)$col
	    set optTail [string range $opt 6 end]	;# remove the -label
	    if {[string compare $val ""] == 0} {
		#
		# Apply the value of the corresponding widget configuration
		# option to the col'th label and unset data($col$opt)
		#
		configLabel $w -$optTail $data($opt)
		if {[info exists data($col$opt)]} {
		    unset data($col$opt)
		}
	    } else {
		#
		# Apply the given value to the col'th label and save the
		# properly formatted value of val in data($col$opt)
		#
		configLabel $w -$optTail $val
		set data($col$opt) [$w cget -$optTail]
	    }

	    #
	    # Adjust the columns (including the height of the header frame)
	    #
	    adjustColumns $win l$col 1
	}

	-labelcommand -
	-name -
	-sortcommand {
	    if {[string compare $val ""] == 0} {
		if {[info exists data($col$opt)]} {
		    unset data($col$opt)
		}
	    } else {
		set data($col$opt) $val
	    }
	}

	-labelfont {
	    set w $data(hdrTxtFrLbl)$col
	    set optTail [string range $opt 6 end]	;# remove the -label
	    if {[string compare $val ""] == 0} {
		#
		# Apply the value of the corresponding widget
		# configuration option to the col'th label and
		# its children (if any), and unset data($col$opt)
		#
		$w configure -$optTail $data($opt)
		foreach c [winfo children $w] {
		    $c configure -$optTail $data($opt)
		}
		if {[info exists data($col$opt)]} {
		    unset data($col$opt)
		}
	    } else {
		#
		# Apply the given value to the col'th label and
		# its children (if any), and save the properly
		# formatted value of val in data($col$opt)
		#
		$w configure -$optTail $val
		foreach c [winfo children $w] {
		    $c configure -$optTail $val
		}
		set data($col$opt) [$w cget -$optTail]
	    }

	    #
	    # Conditionally resize the canvas displaying an up- or down-arrow
	    # and adjust the columns (including the height of the header frame)
	    #
	    if {$col == $data(arrowCol)} {
		configCanvas $win
		drawArrows $win
	    }
	    adjustColumns $win l$col 1
	}

	-labelforeground {
	    set w $data(hdrTxtFrLbl)$col
	    set optTail [string range $opt 6 end]	;# remove the -label
	    if {[string compare $val ""] == 0} {
		#
		# Apply the value of the corresponding widget
		# configuration option to the col'th label and
		# its children (if any), and unset data($col$opt)
		#
		$w configure -$optTail $data($opt)
		foreach c [winfo children $w] {
		    $c configure -$optTail $data($opt)
		}
		if {[info exists data($col$opt)]} {
		    unset data($col$opt)
		}
	    } else {
		#
		# Apply the given value to the col'th label and
		# its children (if any), and save the properly
		# formatted value of val in data($col$opt)
		#
		$w configure -$optTail $val
		foreach c [winfo children $w] {
		    $c configure -$optTail $val
		}
		set data($col$opt) [$w cget -$optTail]
	    }
	}

	-labelheight -
	-labelpady {
	    set w $data(hdrTxtFrLbl)$col
	    set optTail [string range $opt 6 end]	;# remove the -label
	    if {[string compare $val ""] == 0} {
		#
		# Apply the value of the corresponding widget configuration
		# option to the col'th label and unset data($col$opt)
		#
		configLabel $w -$optTail $data($opt)
		if {[info exists data($col$opt)]} {
		    unset data($col$opt)
		}
	    } else {
		#
		# Apply the given value to the col'th label and save the
		# properly formatted value of val in data($col$opt)
		#
		configLabel $w -$optTail $val
		set data($col$opt) [$w cget -$optTail]
	    }

	    #
	    # Adjust the height of the header frame
	    #
	    adjustHeaderHeight $win
	}

	-labelimage {
	    set w $data(hdrTxtFrLbl)$col
	    if {[string compare $val ""] == 0} {
		foreach c [winfo children $w] {
		    destroy $c
		}
		if {[info exists data($col$opt)]} {
		    unset data($col$opt)
		}
	    } else {
		if {![winfo exists $w.il]} {
		    variable configSpecs
		    variable configOpts
		    foreach c [list $w.il $w.tl] {	;# image and text labels
			#
			# Create the label $c
			#
			label $c -borderwidth 0 -height 0 \
				 -highlightthickness 0 -padx 0 \
				 -pady 0 -takefocus 0 -width 0

			#
			# Apply to it the current configuration options
			#
			foreach opt2 $configOpts {
			    if {[string compare \
				 [lindex $configSpecs($opt2) 2] c] == 0} {
				$c configure $opt2 $data($opt2)
			    }
			}
			foreach opt2 {-background -font -foreground} {
			    $c configure $opt2 [$w cget $opt2]
			}
			foreach opt2 {-disabledforeground -state} {
			    catch {$c configure $opt2 [$w cget $opt2]}
			}

			#
			# Replace the binding tag Label with TablelistSubLabel
			# in the list of binding tags of the label $c
			#
			bindtags $c \
				 [lreplace [bindtags $c] 1 1 TablelistSubLabel]
		    }
		}

		#
		# Display the specified image in the label
		# $w.il and save val in data($col$opt)
		#
		$w.il configure -image $val
		set data($col$opt) $val
	    }

	    #
	    # Adjust the columns (including the height of the header frame)
	    #
	    adjustColumns $win l$col 1
	}

	-labelrelief {
	    set w $data(hdrTxtFrLbl)$col
	    set optTail [string range $opt 6 end]	;# remove the -label
	    if {[string compare $val ""] == 0} {
		#
		# Apply the value of the corresponding widget configuration
		# option to the col'th label and unset data($col$opt)
		#
		configLabel $w -$optTail $data($opt)
		if {[info exists data($col$opt)]} {
		    unset data($col$opt)
		}
	    } else {
		#
		# Apply the given value to the col'th label and save the
		# properly formatted value of val in data($col$opt)
		#
		configLabel $w -$optTail $val
		set data($col$opt) [$w cget -$optTail]
	    }
	}

	-maxwidth {
	    #
	    # Save the properly formatted value of val in
	    # data($col$opt), adjust the columns, and make
	    # sure the items will be redisplayed at idle time
	    #
	    set val [format %d $val]	;# integer check with error message
	    set data($col$opt) $val
	    if {$val > 0} {		;# convention: max. width in characters
		set pixels [charsToPixels $win $data(-font) $val]
	    } elseif {$val < 0} {	;# convention: max. width in pixels
		set pixels [expr {(-1)*$val}]
	    } else {			;# convention: no max. width
		set pixels 0
	    }
	    set data($col-maxPixels) $pixels
	    adjustColumns $win $col 1
	    redisplayWhenIdle $win
	}

	-selectbackground -
	-selectforeground {
	    set w $data(body)
	    set tag $col$opt
	    if {[string compare $val ""] == 0} {
		if {[info exists data($tag)]} {
		    $w tag delete $tag
		    unset data($tag)
		}
	    } else {
		#
		# Configure the tag $tag in the body text widget
		#
		set optTail [string range $opt 7 end]	;# remove the -select
		$w tag configure $tag -$optTail $val
		$w tag raise $tag select

		if {!$data($col-hide)} {
		    #
		    # Apply the tag to the elements of the given
		    # column in all currently selected rows
		    #
		    set selRange [$w tag nextrange select 1.0]
		    while {[llength $selRange] != 0} {
			set selStart [lindex $selRange 0]
			set selEnd [lindex $selRange 1]

			findCellTabs $win [expr {int($selStart)}] $col \
				     tabIdx1 tabIdx2
			$w tag add $tag $tabIdx1 $tabIdx2+1c

			set selRange [$w tag nextrange select $selEnd]
		    }
		}

		#
		# Save the properly formatted value of val in data($tag)
		#
		set data($tag) [$w tag cget $tag -$optTail]
	    }
	}

	-showarrow {
	    #
	    # Save the boolean value specified by val in data($col$opt) and
	    # conditionally unmanage the canvas displaying an up- or down-arrow
	    #
	    set data($col$opt) [expr {$val ? 1 : 0}]
	    if {!$data($col$opt) && $col == $data(arrowCol)} {
		place forget $data(hdrTxtFrCanv)
		set data(arrowCol) -1
		adjustColumns $win l$col 1
	    }
	}

	-sortmode {
	    #
	    # Save the properly formatted value of val in data($col$opt)
	    #
	    variable sortModes
	    set data($col$opt) [mwutil::fullOpt "sort mode" $val $sortModes]
	}

	-text {
	    if {$data(isDisabled)} {
		return ""
	    }

	    #
	    # Replace the column's contents in the internal list
	    #
	    set newItemList {}
	    set row 0
	    foreach item $data(itemList) text [lrange $val 0 $data(itemCount)] {
		set item [lreplace $item $col $col $text]
		lappend newItemList $item
	    }
	    set data(itemList) $newItemList

	    #
	    # Update the list variable if present
	    #
	    condUpdateListVar $win

	    #
	    # Adjust the columns and make sure the
	    # items will be redisplayed at idle time
	    #
	    adjustColumns $win $col 1
	    redisplayWhenIdle $win
	}

	-title {
	    #
	    # Save the given value in the corresponding
	    # element of data(-columns) and adjust the columns
	    #
	    set idx [expr {3*$col + 1}]
	    set data(-columns) [lreplace $data(-columns) $idx $idx $val]
	    adjustColumns $win l$col 1
	}

	-width {
	    #
	    # Set up and adjust the columns, and make sure
	    # the items will be redisplayed at idle time
	    #
	    set idx [expr {3*$col}]
	    if {$val != [lindex $data(-columns) $idx]} {
		setupColumns $win [lreplace $data(-columns) $idx $idx $val] 0
		adjustColumns $win $col 1
		redisplayWhenIdle $win
	    }
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::doColCget
#
# Returns the value of the column configuration option opt for the col'th
# column of the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::doColCget {col win opt} {
    upvar ::tablelist::ns${win}::data data

    switch -- $opt {
	-align {
	    return [lindex $data(-columns) [expr {3*$col + 2}]]
	}

	-text {
	    set result {}
	    foreach item $data(itemList) {
		lappend result [lindex $item $col]
	    }
	    return $result
	}

	-title {
	    return [lindex $data(-columns) [expr {3*$col + 1}]]
	}

	-width {
	    return [lindex $data(-columns) [expr {3*$col}]]
	}

	default {
	    if {[info exists data($col$opt)]} {
		return $data($col$opt)
	    } else {
		return ""
	    }
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::doRowConfig
#
# Applies the value val of the row configuration option opt to the row'th row
# of the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::doRowConfig {row win opt val} {
    upvar ::tablelist::ns${win}::data data

    set w $data(body)

    switch -- $opt {
	-background -
	-foreground {
	    set item [lindex $data(itemList) $row]
	    set key [lindex $item end]
	    set tag $key$opt

	    if {[string compare $val ""] == 0} {
		if {[info exists data($tag)]} {
		    $w tag delete $tag
		    unset data($tag)
		    incr data(tagCount) -1
		}
	    } else {
		#
		# Configure the tag $tag in the body text widget and
		# apply it to the given row if it is not selected
		#
		$w tag configure $tag $opt $val
		$w tag lower $tag disabled
		for {set col 0} {$col < $data(colCount)} {incr col} {
		    if {[info exists data($key-$col$opt)]} {
			$w tag lower $key-$col$opt disabled
		    }
		}
		set line [expr {$row + 1}]
		if {[lsearch -exact [$w tag names $line.0] select] < 0} {
		    $w tag add $tag $line.0 $line.end
		}

		#
		# Save the properly formatted value of val in data($tag)
		#
		if {![info exists data($tag)]} {
		    incr data(tagCount)
		}
		set data($tag) [$w tag cget $tag $opt]
	    }
	}

	-font {
	    #
	    # Save the current cell fonts in a temporary array
	    #
	    set item [lindex $data(itemList) $row]
	    set key [lindex $item end]
	    for {set col 0} {$col < $data(colCount)} {incr col} {
		set oldCellFonts($col) [cellFont $win $key $col]
	    }

	    set tag $key$opt
	    if {[string compare $val ""] == 0} {
		if {[info exists data($tag)]} {
		    $w tag delete $tag
		    unset data($tag)
		    incr data(tagCount) -1
		}
	    } else {
		#
		# Configure the tag $tag in the body text
		# widget and apply it to the given row
		#
		$w tag configure $tag $opt $val
		for {set col 0} {$col < $data(colCount)} {incr col} {
		    if {[info exists data($key-$col$opt)]} {
			$w tag raise $key-$col$opt
		    }
		}
		set line [expr {$row + 1}]
		$w tag add $tag $line.0 $line.end

		#
		# Save the properly formatted value of val in data($tag)
		#
		if {![info exists data($tag)]} {
		    incr data(tagCount)
		}
		set data($tag) [$w tag cget $tag $opt]
	    }

	    set dispItem [strToDispStr $item]
	    set colWidthsChanged 0
	    set colList {}
	    set line [expr {$row + 1}]
	    set textIdx1 $line.1
	    set col 0
	    foreach text [lrange $dispItem 0 $data(lastCol)] \
		    fmtCmdFlag $data(fmtCmdFlagList) \
		    {pixels alignment} $data(colList) {
		if {$data($col-hide)} {
		    incr col
		    continue
		}

		#
		# Adjust the cell text and the image width
		#
		if {$fmtCmdFlag} {
		    set text [uplevel #0 $data($col-formatcommand) \
			      [list [lindex $item $col]]]
		    set text [strToDispStr $text]
		}
		if {[info exists data($key-$col-image)]} {
		    set image $data($key-$col-image)
		    set imageWidth [image width $image]
		} else {
		    set image ""
		    set imageWidth 0
		}
		set cellFont [cellFont $win $key $col]
		set workPixels $pixels
		if {$pixels == 0} {		;# convention: dynamic width
		    if {$data($col-maxPixels) > 0 &&
			$data($col-reqPixels) > $data($col-maxPixels)} {
			set workPixels $data($col-maxPixels)
		    }
		}
		if {$workPixels != 0} {
		    incr workPixels $data($col-delta)
		}
		adjustElem $win text imageWidth $cellFont \
			   $workPixels $alignment $data(-snipstring)

		if {$row == $data(editRow) && $col == $data(editCol)} {
		    #
		    # Configure the edit window
		    #
		    variable editWin
		    set name $data($col-editwindow)
		    catch {
			$data(bodyFrEd) configure $editWin($name-fontOpt) \
			    [cellFont $win $key $col]
		    }
		} else {
		    #
		    # Delete the old cell contents between the
		    # two tabs, and insert the text and the image
		    #
		    set textIdx2 [$w search \t $textIdx1 $line.end]
		    $w delete $textIdx1 $textIdx2
		    insertElem $w $textIdx1 $text $image $imageWidth $alignment
		}

		if {$pixels == 0} {		;# convention: dynamic width
		    #
		    # Check whether the width of the current column has changed
		    #
		    if {$workPixels > 0} {
			adjustElem $win text imageWidth $cellFont \
				   $pixels $alignment $data(-snipstring)
		    }
		    set textWidth [font measure $cellFont -displayof $win $text]
		    set newElemWidth [expr {$imageWidth + $textWidth}]
		    if {$newElemWidth > $data($col-elemWidth)} {
			set data($col-elemWidth) $newElemWidth
			set data($col-widestCount) 1
			if {$newElemWidth > $data($col-reqPixels)} {
			    set data($col-reqPixels) $newElemWidth
			    set colWidthsChanged 1
			}
		    } else {
			set oldTextWidth [font measure $oldCellFonts($col) \
					  -displayof $win $text]
			set oldElemWidth [expr {$imageWidth + $oldTextWidth}]
			if {$oldElemWidth < $data($col-elemWidth) &&
			    $newElemWidth == $data($col-elemWidth)} {
			    incr data($col-widestCount)
			} elseif {$oldElemWidth == $data($col-elemWidth) &&
				  $newElemWidth < $oldElemWidth &&
				  [incr data($col-widestCount) -1] == 0} {
			    set colWidthsChanged 1
			    lappend colList $col
			}
		    }
		}

		set textIdx1 [$w search \t $textIdx1 $line.end]+2c
		incr col
	    }

	    #
	    # Adjust the columns if necessary
	    #
	    if {$colWidthsChanged} {
		adjustColumns $win $colList 1
	    }
	}

	-selectable {
	    set val [expr {$val ? 1 : 0}]
	    set item [lindex $data(itemList) $row]
	    set key [lindex $item end]

	    if {$val} {
		if {[info exists data($key$opt)]} {
		    unset data($key$opt)
		}
	    } else {
		#
		# Set data($key$opt) to 0 and deselect the row
		#
		set data($key$opt) 0
		selectionSubCmd $win clear $row $row
	    }
	}

	-selectbackground -
	-selectforeground {
	    set item [lindex $data(itemList) $row]
	    set key [lindex $item end]
	    set tag $key$opt

	    if {[string compare $val ""] == 0} {
		if {[info exists data($tag)]} {
		    $w tag delete $tag
		    unset data($tag)
		}
	    } else {
		#
		# Configure the tag $tag in the body text widget
		# and apply it to the given row if it is selected
		#
		set optTail [string range $opt 7 end]	;# remove the -select
		$w tag configure $tag -$optTail $val
		$w tag lower $tag disabled
		for {set col 0} {$col < $data(colCount)} {incr col} {
		    if {[info exists data($key-$col$opt)]} {
			$w tag lower $key-$col$opt disabled
		    }
		}
		set line [expr {$row + 1}]
		if {[lsearch -exact [$w tag names $line.0] select] >= 0} {
		    $w tag add $tag $line.0 $line.end
		}

		#
		# Save the properly formatted value of val in data($tag)
		#
		set data($tag) [$w tag cget $tag -$optTail]
	    }
	}

	-text {
	    if {$data(isDisabled)} {
		return ""
	    }

	    set colWidthsChanged 0
	    set colList {}
	    set oldItem [lindex $data(itemList) $row]
	    set key [lindex $oldItem end]
	    set newItem [adjustItem $val $data(colCount)]
	    set line [expr {$row + 1}]
	    set textIdx1 $line.1
	    set col 0
	    foreach text [strToDispStr $newItem] \
		    fmtCmdFlag $data(fmtCmdFlagList) \
		    {pixels alignment} $data(colList) {
		if {$data($col-hide)} {
		    incr col
		    continue
		}

		#
		# Adjust the cell text and the image width
		#
		if {$fmtCmdFlag} {
		    set text [uplevel #0 $data($col-formatcommand) \
			      [list [lindex $newItem $col]]]
		    set text [strToDispStr $text]
		}
		if {[info exists data($key-$col-image)]} {
		    set image $data($key-$col-image)
		    set imageWidth [image width $image]
		} else {
		    set image ""
		    set imageWidth 0
		}
		set cellFont [cellFont $win $key $col]
		set workPixels $pixels
		if {$pixels == 0} {		;# convention: dynamic width
		    if {$data($col-maxPixels) > 0 &&
			$data($col-reqPixels) > $data($col-maxPixels)} {
			set workPixels $data($col-maxPixels)
		    }
		}
		if {$workPixels != 0} {
		    incr workPixels $data($col-delta)
		}
		adjustElem $win text imageWidth $cellFont \
			   $workPixels $alignment $data(-snipstring)

		if {$row != $data(editRow) || $col != $data(editCol)} {
		    #
		    # Delete the old cell contents between the
		    # two tabs, and insert the text and the image
		    #
		    set textIdx2 [$w search \t $textIdx1 $line.end]
		    $w delete $textIdx1 $textIdx2
		    insertElem $w $textIdx1 $text $image $imageWidth $alignment
		}

		if {$pixels == 0} {		;# convention: dynamic width
		    #
		    # Check whether the width of the current column has changed
		    #
		    if {$workPixels > 0} {
			adjustElem $win text imageWidth $cellFont \
				   $pixels $alignment $data(-snipstring)
		    }
		    set textWidth [font measure $cellFont -displayof $win $text]
		    set newElemWidth [expr {$imageWidth + $textWidth}]
		    if {$newElemWidth > $data($col-elemWidth)} {
			set data($col-elemWidth) $newElemWidth
			set data($col-widestCount) 1
			if {$newElemWidth > $data($col-reqPixels)} {
			    set data($col-reqPixels) $newElemWidth
			    set colWidthsChanged 1
			}
		    } else {
			set oldText [lindex $oldItem $col]
			if {$fmtCmdFlag} {
			    set oldText [uplevel #0 $data($col-formatcommand) \
					 [list $oldText]]
			}
			set oldText [strToDispStr $oldText]
			adjustElem $win oldText imageWidth $cellFont \
				   $pixels $alignment $data(-snipstring)
			set oldTextWidth \
			    [font measure $cellFont -displayof $win $oldText]
			set oldElemWidth [expr {$imageWidth + $oldTextWidth}]
			if {$oldElemWidth < $data($col-elemWidth) &&
			    $newElemWidth == $data($col-elemWidth)} {
			    incr data($col-widestCount)
			} elseif {$oldElemWidth == $data($col-elemWidth) &&
				  $newElemWidth < $oldElemWidth &&
				  [incr data($col-widestCount) -1] == 0} {
			    set colWidthsChanged 1
			    lappend colList $col
			}
		    }
		}

		set textIdx1 [$w search \t $textIdx1 $line.end]+2c
		incr col
	    }

	    #
	    # Replace the row contents in the list variable if present
	    #
	    if {$data(hasListVar)} {
		trace vdelete ::$data(-listvariable) wu $data(listVarTraceCmd)
		upvar #0 $data(-listvariable) var
		set var [lreplace $var $row $row $newItem]
		trace variable ::$data(-listvariable) wu $data(listVarTraceCmd)
	    }

	    #
	    # Replace the row contents in the internal list
	    #
	    lappend newItem [lindex $oldItem end]
	    set data(itemList) [lreplace $data(itemList) $row $row $newItem]

	    #
	    # Adjust the columns if necessary
	    #
	    if {$colWidthsChanged} {
		adjustColumns $win $colList 1
	    }
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::doRowCget
#
# Returns the value of the row configuration option opt for the row'th row of
# the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::doRowCget {row win opt} {
    upvar ::tablelist::ns${win}::data data

    #
    # Return the value of the specified row configuration option
    #
    switch -- $opt {
	-text {
	    set item [lindex $data(itemList) $row]
	    return [lrange $item 0 $data(lastCol)]
	}

	-selectable {
	    set item [lindex $data(itemList) $row]
	    set key [lindex $item end]
	    if {[info exists data($key$opt)]} {
		return $data($key$opt)
	    } else {
		return 1
	    }
	}

	default {
	    set item [lindex $data(itemList) $row]
	    set key [lindex $item end]
	    if {[info exists data($key$opt)]} {
		return $data($key$opt)
	    } else {
		return ""
	    }
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::doCellConfig
#
# Applies the value val of the cell configuration option opt to the cell
# row,col of the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::doCellConfig {row col win opt val} {
    upvar ::tablelist::ns${win}::data data

    set w $data(body)

    switch -- $opt {
	-background -
	-foreground {
	    set item [lindex $data(itemList) $row]
	    set key [lindex $item end]
	    set tag $key-$col$opt

	    if {[string compare $val ""] == 0} {
		if {[info exists data($tag)]} {
		    $w tag delete $tag
		    unset data($tag)
		    incr data(tagCount) -1
		}
	    } else {
		#
		# Configure the tag $tag in the body text widget
		#
		$w tag configure $tag $opt $val
		$w tag lower $tag disabled

		if {!$data($col-hide)} {
		    #
		    # Apply the tag to the given cell if it is not selected
		    #
		    set line [expr {$row + 1}]
		    if {[lsearch -exact [$w tag names $line.0] select] < 0} {
			findCellTabs $win $line $col tabIdx1 tabIdx2
			$w tag add $tag $tabIdx1 $tabIdx2+1c
		    }
		}

		#
		# Save the properly formatted value of val in data($tag)
		#
		if {![info exists data($tag)]} {
		    incr data(tagCount)
		}
		set data($tag) [$w tag cget $tag $opt]
	    }
	}

	-editable {
	    #
	    # Save the boolean value specified by val in data($key-$col$opt)
	    #
	    set item [lindex $data(itemList) $row]
	    set key [lindex $item end]
	    set data($key-$col$opt) [expr {$val ? 1 : 0}]
	}

	-font {
	    #
	    # Save the current cell font
	    #
	    set item [lindex $data(itemList) $row]
	    set key [lindex $item end]
	    set oldCellFont [cellFont $win $key $col]

	    set tag $key-$col$opt
	    if {[string compare $val ""] == 0} {
		if {[info exists data($tag)]} {
		    $w tag delete $tag
		    unset data($tag)
		    incr data(tagCount) -1
		}
	    } else {
		#
		# Configure the tag $tag in the body text widget
		#
		$w tag configure $tag $opt $val
		$w tag raise $tag

		if {!$data($col-hide)} {
		    #
		    # Apply the tag to the given cell
		    #
		    findCellTabs $win [expr {$row + 1}] $col tabIdx1 tabIdx2
		    $w tag add $tag $tabIdx1 $tabIdx2+1c
		}

		#
		# Save the properly formatted value of val in data($tag)
		#
		if {![info exists data($tag)]} {
		    incr data(tagCount)
		}
		set data($tag) [$w tag cget $tag $opt]
	    }

	    set pixels [lindex $data(colList) [expr {2*$col}]]
	    set workPixels $pixels
	    if {$pixels == 0} {			;# convention: dynamic width
		if {$data($col-maxPixels) > 0 &&
		    $data($col-reqPixels) > $data($col-maxPixels)} {
		    set workPixels $data($col-maxPixels)
		}
	    }
	    if {$workPixels != 0} {
		incr workPixels $data($col-delta)
	    }
	    set alignment [lindex $data(colList) [expr {2*$col + 1}]]

	    #
	    # Adjust the cell text and the image width
	    #
	    set text [lindex $item $col]
	    if {[info exists data($col-formatcommand)]} {
		set text [uplevel #0 $data($col-formatcommand) [list $text]]
	    }
	    set text [strToDispStr $text]
	    if {[info exists data($key-$col-image)]} {
		set image $data($key-$col-image)
		set imageWidth [image width $image]
	    } else {
		set image ""
		set imageWidth 0
	    }
	    set cellFont [cellFont $win $key $col]
	    adjustElem $win text imageWidth $cellFont \
		       $workPixels $alignment $data(-snipstring)

	    if {!$data($col-hide)} {
		if {$row == $data(editRow) && $col == $data(editCol)} {
		    #
		    # Configure the edit window
		    #
		    variable editWin
		    set name $data($col-editwindow)
		    catch {
			$data(bodyFrEd) configure $editWin($name-fontOpt) \
			    [cellFont $win $key $col]
		    }
		} else {
		    #
		    # Delete the old cell contents between the
		    # two tabs, and insert the text and the image
		    #
		    findCellTabs $win [expr {$row + 1}] $col tabIdx1 tabIdx2
		    $w delete $tabIdx1+1c $tabIdx2
		    insertElem $w $tabIdx1+1c $text $image $imageWidth \
			       $alignment
		}
	    }

	    #
	    # Adjust the columns if necessary
	    #
	    if {$pixels == 0} {			;# convention: dynamic width
		if {$workPixels > 0} {
		    adjustElem $win text imageWidth $cellFont \
			       $pixels $alignment $data(-snipstring)
		}
		set textWidth [font measure $cellFont -displayof $win $text]
		set newElemWidth [expr {$imageWidth + $textWidth}]
		if {$newElemWidth > $data($col-elemWidth)} {
		    set data($col-elemWidth) $newElemWidth
		    set data($col-widestCount) 1
		    if {$newElemWidth > $data($col-reqPixels)} {
			set data($col-reqPixels) $newElemWidth
			adjustColumns $win {} 1
		    }
		} else {
		    set oldTextWidth \
			[font measure $oldCellFont -displayof $win $text]
		    set oldElemWidth [expr {$imageWidth + $oldTextWidth}]
		    if {$oldElemWidth < $data($col-elemWidth) &&
			$newElemWidth == $data($col-elemWidth)} {
			incr data($col-widestCount)
		    } elseif {$oldElemWidth == $data($col-elemWidth) &&
			      $newElemWidth < $oldElemWidth &&
			      [incr data($col-widestCount) -1] == 0} {
			adjustColumns $win $col 1
		    }
		}
	    }
	}

	-image {
	    if {$data(isDisabled)} {
		return ""
	    }

	    set pixels [lindex $data(colList) [expr {2*$col}]]
	    set workPixels $pixels
	    if {$pixels == 0} {			;# convention: dynamic width
		if {$data($col-maxPixels) > 0 &&
		    $data($col-reqPixels) > $data($col-maxPixels)} {
		    set workPixels $data($col-maxPixels)
		}
	    }
	    if {$workPixels != 0} {
		incr workPixels $data($col-delta)
	    }
	    set alignment [lindex $data(colList) [expr {2*$col + 1}]]

	    #
	    # Adjust the cell text and the image width
	    #
	    set item [lindex $data(itemList) $row]
	    set text [lindex $item $col]
	    if {[info exists data($col-formatcommand)]} {
		set text [uplevel #0 $data($col-formatcommand) [list $text]]
	    }
	    set text [strToDispStr $text]
	    set oldText $text			;# will be needed later
	    set image $val
	    if {[string compare $image ""] == 0} {
		set imageWidth 0
	    } else {
		set imageWidth [image width $image]
	    }
	    set key [lindex $item end]
	    set cellFont [cellFont $win $key $col]
	    adjustElem $win text imageWidth $cellFont \
		       $workPixels $alignment $data(-snipstring)

	    if {!$data($col-hide) &&
		!($row == $data(editRow) && $col == $data(editCol))} {
		#
		# Delete the old cell contents between the
		# two tabs, and insert the text and the image
		#
		findCellTabs $win [expr {$row + 1}] $col tabIdx1 tabIdx2
		$w delete $tabIdx1+1c $tabIdx2
		insertElem $w $tabIdx1+1c $text $image $imageWidth $alignment
	    }

	    #
	    # Save the old image width
	    #
	    set name $key-$col$opt
	    if {[info exists data($name)]} {
		set oldImageWidth [image width $data($name)]
	    } else {
		set oldImageWidth 0
	    }

	    #
	    # Delete data($name) or save the specified value in it
	    #
	    if {[string compare $val ""] == 0} {
		if {[info exists data($name)]} {
		    unset data($name)
		    incr data(imgCount) -1
		}
	    } else {
		if {![info exists data($name)]} {
		    incr data(imgCount)
		}
		set data($name) $val
	    }

	    #
	    # Adjust the columns if necessary
	    #
	    if {$pixels == 0} {			;# convention: dynamic width
		if {$workPixels > 0} {
		    adjustElem $win text imageWidth $cellFont \
			       $pixels $alignment $data(-snipstring)
		}
		set textWidth [font measure $cellFont -displayof $win $text]
		set newElemWidth [expr {$imageWidth + $textWidth}]
		if {$newElemWidth > $data($col-elemWidth)} {
		    set data($col-elemWidth) $newElemWidth
		    set data($col-widestCount) 1
		    if {$newElemWidth > $data($col-reqPixels)} {
			set data($col-reqPixels) $newElemWidth
			adjustColumns $win {} 1
		    }
		} else {
		    adjustElem $win oldText oldImageWidth $cellFont \
			       $pixels $alignment $data(-snipstring)
		    set oldTextWidth \
			[font measure $cellFont -displayof $win $oldText]
		    set oldElemWidth [expr {$oldImageWidth + $oldTextWidth}]
		    if {$oldElemWidth < $data($col-elemWidth) &&
			$newElemWidth == $data($col-elemWidth)} {
			incr data($col-widestCount)
		    } elseif {$oldElemWidth == $data($col-elemWidth) &&
			      $newElemWidth < $oldElemWidth &&
			      [incr data($col-widestCount) -1] == 0} {
			adjustColumns $win $col 1
		    }
		}
	    }
	}

	-selectbackground -
	-selectforeground {
	    set item [lindex $data(itemList) $row]
	    set key [lindex $item end]
	    set tag $key-$col$opt

	    if {[string compare $val ""] == 0} {
		if {[info exists data($tag)]} {
		    $w tag delete $tag
		    unset data($tag)
		}
	    } else {
		#
		# Configure the tag $tag in the body text widget
		#
		set optTail [string range $opt 7 end]	;# remove the -select
		$w tag configure $tag -$optTail $val
		$w tag lower $tag disabled

		if {!$data($col-hide)} {
		    #
		    # Apply the tag to the given cell if it is selected
		    #
		    set line [expr {$row + 1}]
		    if {[lsearch -exact [$w tag names $line.0] select] >= 0} {
			findCellTabs $win $line $col tabIdx1 tabIdx2
			$w tag add $tag $tabIdx1 $tabIdx2+1c
		    }
		}

		#
		# Save the properly formatted value of val in data($tag)
		#
		set data($tag) [$w tag cget $tag -$optTail]
	    }
	}

	-text {
	    if {$data(isDisabled)} {
		return ""
	    }

	    set pixels [lindex $data(colList) [expr {2*$col}]]
	    set workPixels $pixels
	    if {$pixels == 0} {			;# convention: dynamic width
		if {$data($col-maxPixels) > 0 &&
		    $data($col-reqPixels) > $data($col-maxPixels)} {
		    set workPixels $data($col-maxPixels)
		}
	    }
	    if {$workPixels != 0} {
		incr workPixels $data($col-delta)
	    }
	    set alignment [lindex $data(colList) [expr {2*$col + 1}]]

	    #
	    # Adjust the cell text and the image width
	    #
	    set text $val
	    set fmtCmdFlag [info exists data($col-formatcommand)]
	    if {$fmtCmdFlag} {
		set text [uplevel #0 $data($col-formatcommand) [list $text]]
	    }
	    set text [strToDispStr $text]
	    set oldItem [lindex $data(itemList) $row]
	    set key [lindex $oldItem end]
	    if {[info exists data($key-$col-image)]} {
		set image $data($key-$col-image)
		set imageWidth [image width $image]
	    } else {
		set image ""
		set imageWidth 0
	    }
	    set cellFont [cellFont $win $key $col]
	    adjustElem $win text imageWidth $cellFont \
		       $workPixels $alignment $data(-snipstring)

	    if {!$data($col-hide) &&
		!($row == $data(editRow) && $col == $data(editCol))} {
		#
		# Delete the old cell contents between the
		# two tabs, and insert the text and the image
		#
		findCellTabs $win [expr {$row + 1}] $col tabIdx1 tabIdx2
		$w delete $tabIdx1+1c $tabIdx2
		insertElem $w $tabIdx1+1c $text $image $imageWidth $alignment
	    }

	    #
	    # Replace the cell contents in the internal list
	    #
	    set newItem [lreplace $oldItem $col $col $val]
	    set data(itemList) [lreplace $data(itemList) $row $row $newItem]

	    #
	    # Replace the cell contents in the list variable if present
	    #
	    if {$data(hasListVar)} {
		trace vdelete ::$data(-listvariable) wu $data(listVarTraceCmd)
		upvar #0 $data(-listvariable) var
		set var [lreplace $var $row $row \
			 [lrange $newItem 0 $data(lastCol)]]
		trace variable ::$data(-listvariable) wu $data(listVarTraceCmd)
	    }

	    #
	    # Adjust the columns if necessary
	    #
	    if {$pixels == 0} {			;# convention: dynamic width
		if {$workPixels > 0} {
		    adjustElem $win text imageWidth $cellFont \
			       $pixels $alignment $data(-snipstring)
		}
		set textWidth [font measure $cellFont -displayof $win $text]
		set newElemWidth [expr {$imageWidth + $textWidth}]
		if {$newElemWidth > $data($col-elemWidth)} {
		    set data($col-elemWidth) $newElemWidth
		    set data($col-widestCount) 1
		    if {$newElemWidth > $data($col-reqPixels)} {
			set data($col-reqPixels) $newElemWidth
			adjustColumns $win {} 1
		    }
		} else {
		    set oldText [lindex $oldItem $col]
		    if {$fmtCmdFlag} {
			set oldText [uplevel #0 $data($col-formatcommand) \
				     [list $oldText]]
		    }
		    set oldText [strToDispStr $oldText]
		    adjustElem $win oldText imageWidth $cellFont \
			       $pixels $alignment $data(-snipstring)
		    set oldTextWidth \
			[font measure $cellFont -displayof $win $oldText]
		    set oldElemWidth [expr {$imageWidth + $oldTextWidth}]
		    if {$oldElemWidth < $data($col-elemWidth) &&
			$newElemWidth == $data($col-elemWidth)} {
			incr data($col-widestCount)
		    } elseif {$oldElemWidth == $data($col-elemWidth) &&
			      $newElemWidth < $oldElemWidth &&
			      [incr data($col-widestCount) -1] == 0} {
			adjustColumns $win $col 1
		    }
		}
	    }
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::doCellCget
#
# Returns the value of the cell configuration option opt for the cell row,col
# of the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::doCellCget {row col win opt} {
    upvar ::tablelist::ns${win}::data data

    #
    # Return the value of the specified cell configuration option
    #
    switch -- $opt {
	-editable {
	    return [isCellEditable $win $row $col]
	}

	-text {
	    set item [lindex $data(itemList) $row]
	    return [lindex $item $col]
	}

	default {
	    set item [lindex $data(itemList) $row]
	    set key [lindex $item end]
	    if {[info exists data($key-$col$opt)]} {
		return $data($key-$col$opt)
	    } else {
		return ""
	    }
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::makeListVar
#
# Arranges for the global variable specified by varName to become the list
# variable associated with the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::makeListVar {win varName} {
    upvar ::tablelist::ns${win}::data data

    if {[string compare $varName ""] == 0} {
	#
	# If there is an old list variable associated with the
	# widget then remove the trace set on this variable
	#
	if {$data(hasListVar)} {
	    synchronize $win
	    trace vdelete ::$data(-listvariable) wu $data(listVarTraceCmd)
	}
	return ""
    }

    #
    # The list variable may be an array element but must not be an array
    #
    if {![regexp {^(.*)\((.*)\)$} $varName dummy name1 name2]} {
	if {[array exists ::$varName]} {
	    return -code error "variable \"$varName\" is array"
	}
	set name1 $varName
	set name2 ""
    }

    if {[info exists ::$varName]} {
	#
	# Invoke the trace procedure associated with the new list variable
	#
	listVarTrace $win $name1 $name2 w
    } else {
	#
	# Set ::$varName according to the value of data(itemList)
	#
	set ::$varName {}
	foreach item $data(itemList) {
	    lappend ::$varName [lrange $item 0 $data(lastCol)]
	}
    }

    #
    # If there is an old list variable associated with the
    # widget then remove the trace set on this variable,
    # and in any case set a trace on the new list variable
    #
    if {$data(hasListVar)} {
	trace vdelete ::$data(-listvariable) wu $data(listVarTraceCmd)
    }
    trace variable ::$varName wu $data(listVarTraceCmd)
}
