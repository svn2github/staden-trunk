#==============================================================================
# Contains private procedures used in tablelist bindings.
#
# Copyright (c) 2000-2004  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

#
# Binding tag Tablelist
# =====================
#

#------------------------------------------------------------------------------
# tablelist::addActiveTag
#
# This procedure is invoked when the tablelist widget win gains the keyboard
# focus.  It adds the "active" tag to the line that displays the active element
# of the widget in its body text child.
#------------------------------------------------------------------------------
proc tablelist::addActiveTag win {
    upvar ::tablelist::ns${win}::data data

    set line [expr {$data(activeIdx) + 1}]
    $data(body) tag add active $line.0 $line.end

    set data(ownsFocus) 1
}

#------------------------------------------------------------------------------
# tablelist::removeActiveTag
#
# This procedure is invoked when the tablelist widget win loses the keyboard
# focus.  It removes the "active" tag from the line that displays the active
# element of the widget in its body text child.
#------------------------------------------------------------------------------
proc tablelist::removeActiveTag win {
    upvar ::tablelist::ns${win}::data data

    set line [expr {$data(activeIdx) + 1}]
    $data(body) tag remove active $line.0 $line.end

    set data(ownsFocus) 0
}

#------------------------------------------------------------------------------
# tablelist::cleanup
#
# This procedure is invoked when the tablelist widget win is destroyed.  It
# executes some cleanup operations.
#------------------------------------------------------------------------------
proc tablelist::cleanup win {
    upvar ::tablelist::ns${win}::data data

    #
    # Cancel the execution of all delayed adjustSeps, makeStripes,
    # stretchColumns, synchronize, and redisplay commands
    #
    foreach afterId {sepsId stripesId stretchId syncId redispId} {
	if {[info exists data($afterId)]} {
	    after cancel $data($afterId)
	}
    }

    #
    # If there is a list variable associated with the
    # widget then remove the trace set on this variable
    #
    if {$data(hasListVar)} {
	trace vdelete ::$data(-listvariable) wu $data(listVarTraceCmd)
    }

    namespace delete ::tablelist::ns$win
    catch {rename ::$win ""}
}

#
# Binding tag TablelistBody
# =========================
#

#------------------------------------------------------------------------------
# tablelist::defineTablelistBody
#
# Defines the binding tag TablelistBody to have the same events as Listbox and
# the binding scripts obtained from those of Listbox by replacing the widget %W
# with its parent as well as the %x and %y fields with the corresponding
# coordinates relative to the parent.
#------------------------------------------------------------------------------
proc tablelist::defineTablelistBody {} {
    bind TablelistBody <Button-1> {
	if {[winfo exists %W]} {
	    tablelist::condEditContainingCell %W %x %y
	}
    }

    bind TablelistBody <ButtonRelease-1> {
	if {[winfo exists %W]} {
	    tablelist::condEvalInvokeCmd %W %X %Y
	}
    }

    foreach event [bind Listbox] {
	set script [strMap {
	    %W $tablelist::W  %x $tablelist::x  %y $tablelist::y
	    tkListboxAutoScan   tablelist::tablelistAutoScan
	    tk::ListboxAutoScan tablelist::tablelistAutoScan
	} [bind Listbox $event]]

	bind TablelistBody $event +[format {
	    if {[winfo exists %%W]} {
		set tablelist::W [winfo parent %%W]
		set tablelist::x [expr {%%x + [winfo x %%W]}]
		set tablelist::y [expr {%%y + [winfo y %%W]}]
		%s
	    }
	} $script]
    }
}

#------------------------------------------------------------------------------
# tablelist::tablelistAutoScan
#
# This is a modified version of the Tk library procedure tk(::)ListboxAutoScan.
# It is invoked when the mouse leaves the body text child of a tablelist
# widget.  It scrolls the child and reschedules itself as an after command so
# that the child continues to scroll until the mouse moves back into the window
# or the mouse button is released.
#------------------------------------------------------------------------------
proc tablelist::tablelistAutoScan win {
    if {![winfo exists $win] || [string compare [::$win editwinpath] ""] != 0} {
	return ""
    }

    if {[array exists ::tk::Priv]} {
	set x $::tk::Priv(x)
	set y $::tk::Priv(y)
    } else {
	set x $::tkPriv(x)
	set y $::tkPriv(y)
    }

    set w [::$win bodypath]
    set _x [expr {$x - [winfo x $w]}]
    set _y [expr {$y - [winfo y $w]}]

    if {$_y >= [winfo height $w]} {
	yviewSubCmd $win {scroll 1 units}
    } elseif {$_y < 0} {
	yviewSubCmd $win {scroll -1 units}
    } elseif {$_x >= [winfo width $w]} {
	xviewSubCmd $win {scroll 2 units}
    } elseif {$_x < 0} {
	xviewSubCmd $win {scroll -2 units}
    } else {
	return ""
    }

    if {[array exists ::tk::Priv]} {
	tk::ListboxMotion $win [rowIndex $win @$x,$y 1]
	set ::tk::Priv(afterId) \
	    [after 50 [list tablelist::tablelistAutoScan $win]]
    } else {
	tkListboxMotion $win [rowIndex $win @$x,$y 1]
	set ::tkPriv(afterId) \
	    [after 50 [list tablelist::tablelistAutoScan $win]]
    }
}

#------------------------------------------------------------------------------
# tablelist::condEditContainingCell
#
# This procedure is invoked when mouse button 1 is pressed in the body w of a
# tablelist widget or in one of its separator frames.  If the mouse click
# occurred inside an editable cell and the latter is not already being edited,
# then the procedure starts the interactive editing in that cell.  Otherwise it
# finishes a possibly active cell editing.
#------------------------------------------------------------------------------
proc tablelist::condEditContainingCell {w x y} {
    set win [winfo parent $w]
    upvar ::tablelist::ns${win}::data data
    synchronize $win

    #
    # Get the containing cell from the coordinates relative to the parent
    #
    incr x [winfo x $w]
    incr y [winfo y $w]
    set row [containingSubCmd $win $y]
    set col [containingcolumnSubCmd $win $x]

    if {$row >= 0 && $col >= 0 && [isCellEditable $win $row $col]} {
	#
	# Get the coordinates relative to the
	# tablelist body and invoke editcellSubCmd
	#
	set w $data(body)
	incr x -[winfo x $w]
	incr y -[winfo y $w]
	scan [$w index @$x,$y] %d.%d line charPos
	editcellSubCmd $win $row $col 0 $charPos
    } else {
	#
	# Finish a possibly active cell editing
	#
	if {$data(editRow) >= 0} {
	    finisheditingSubCmd $win
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::condEvalInvokeCmd
#
# This procedure is invoked when mouse button 1 is released in the body w of a
# tablelist widget or in one of its separator frames.  If interactive cell
# editing is in progress in a column whose associated edit window has an invoke
# command that hasn't yet been called in the current edit session, then the
# procedure evaluates that command.
#------------------------------------------------------------------------------
proc tablelist::condEvalInvokeCmd {w X Y} {
    set win [winfo parent $w]
    upvar ::tablelist::ns${win}::data data

    if {$data(editCol) < 0} {
	return ""
    }

    variable editWin
    set name $data($data(editCol)-editwindow)
    if {[string compare $editWin($name-invokeCmd) ""] == 0 || $data(invoked)} {
	return ""
    }

    #
    # If the edit window is a checkbutton then ignore this event outside the
    # edit window or if the checkbutton's selection state has already been
    # toggled (the latter is the case if the windowing system equals "x11")
    #
    if {[string compare [winfo class $data(bodyFrEd)] Checkbutton] == 0} {
	set c [winfo containing -displayof $w $X $Y]
	if {[string compare $c $data(bodyFrEd)] != 0 ||
	    [string compare $data(editText) $data(origEditText)] != 0} {
	    return ""
	}
    }

    eval [strMap {%W $data(bodyFrEd)} $editWin($name-invokeCmd)]
    set data(invoked) 1
}

#
# Binding tags TablelistLabel, TablelistSubLabel, and TablelistArrow
# ==================================================================
#

#------------------------------------------------------------------------------
# tablelist::defineTablelistSubLabel
#
# Defines the binding tag TablelistSubLabel (for children of tablelist labels)
# to have the same events as TablelistLabel and the binding scripts obtained
# from those of TablelistLabel by replacing the widget %W with its parent as
# well as the %x and %y fields with the corresponding coordinates relative to
# the parent.
#------------------------------------------------------------------------------
proc tablelist::defineTablelistSubLabel {} {
    foreach event [bind TablelistLabel] {
	set script [strMap {
	    %W $tablelist::W  %x $tablelist::x  %y $tablelist::y
	} [bind TablelistLabel $event]]

	bind TablelistSubLabel $event [format {
	    set tablelist::W [winfo parent %%W]
	    set tablelist::x [expr {%%x + [winfo x %%W]}]
	    set tablelist::y [expr {%%y + [winfo y %%W]}]
	    %s
	} $script]
    }
}

#------------------------------------------------------------------------------
# tablelist::defineTablelistArrow
#
# Defines the binding tag TablelistArrow to have the same events as
# TablelistLabel and the binding scripts obtained from those of TablelistLabel
# by replacing the widget %W with the containing label as well as the %x and %y
# fields with the corresponding coordinates relative to the label
#------------------------------------------------------------------------------
proc tablelist::defineTablelistArrow {} {
    foreach event [bind TablelistLabel] {
	set script [strMap {
	    %W $tablelist::W  %x $tablelist::x  %y $tablelist::y
	} [bind TablelistLabel $event]]

	bind TablelistArrow $event [format {
	    if {$::tk_version < 8.4} {
		regexp {^.+ -in (.+)$} [place info %%W] \
		       tablelist::dummy tablelist::W
	    } else {
		set tablelist::W [lindex [place configure %%W -in] end]
	    }
	    set tablelist::x \
		[expr {%%x + [winfo x %%W] - [winfo x $tablelist::W]}]
	    set tablelist::y \
		[expr {%%y + [winfo y %%W] - [winfo y $tablelist::W]}]
	    %s
	} $script]
    }
}

#------------------------------------------------------------------------------
# tablelist::labelEnter
#
# This procedure is invoked when the mouse pointer enters the header label w of
# a tablelist widget, or is moving within that label.  It updates the cursor,
# depending on whether the pointer is on the right border of the label or not.
#------------------------------------------------------------------------------
proc tablelist::labelEnter {w x} {
    parseLabelPath $w win col
    upvar ::tablelist::ns${win}::data data

    configLabel $w -cursor $data(-cursor)
    if {$data(isDisabled)} {
	return ""
    }

    if {$data(-resizablecolumns) && $data($col-resizable) &&
	$x >= [winfo width $w] - [$w cget -borderwidth] - 4} {
	configLabel $w -cursor $data(-resizecursor)
    }
}

#------------------------------------------------------------------------------
# tablelist::labelB1Down
#
# This procedure is invoked when mouse button 1 is pressed in the header label
# w of a tablelist widget.  If the pointer is on the right border of the label
# then the procedure records its x-coordinate relative to the label, the width
# of the column, and some other data needed later.  Otherwise it saves the
# label's relief so it can be restored later, and changes the relief to sunken.
#------------------------------------------------------------------------------
proc tablelist::labelB1Down {w x} {
    parseLabelPath $w win col
    upvar ::tablelist::ns${win}::data data

    if {$data(isDisabled) ||
	[info exists data(x)]} {		;# resize operation in progress
	return ""
    }

    set data(labelClicked) 1
    set labelWidth [winfo width $w]

    if {$data(-resizablecolumns) && $data($col-resizable) &&
	$x >= $labelWidth - [$w cget -borderwidth] - 4} {
	set data(x) $x

	set data(oldStretchedColWidth) [expr {$labelWidth - 2*$data(charWidth)}]
	set data(oldColDelta) $data($col-delta)
	set data(configColWidth) [lindex $data(-columns) [expr {3*$col}]]

	if {$col == $data(arrowCol)} {
	    set data(minColWidth) $data(arrowWidth)
	} else {
	    set data(minColWidth) 1
	}

	set topWin [winfo toplevel $win]
	set data(topBindEsc) [bind $topWin <Escape>]
	bind $topWin <Escape> [list tablelist::escape [strMap {% %%} $win] $col]
    } else {
	set data(inClickedLabel) 1
	set data(relief) [$w cget -relief]

	if {[info exists data($col-labelcommand)] ||
	    [string compare $data(-labelcommand) ""] != 0} {
	    set data(changeRelief) 1
	    configLabel $w -relief sunken
	} else {
	    set data(changeRelief) 0
	}

	if {$data(-movablecolumns)} {
	    set topWin [winfo toplevel $win]
	    set data(topBindEsc) [bind $topWin <Escape>]
	    bind $topWin <Escape> \
		 [list tablelist::escape [strMap {% %%} $win] $col]
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::labelB1Motion
#
# This procedure is invoked to process mouse motion events in the header label
# w of a tablelist widget while button 1 is down.  If this event occured during
# a column resize operation then the procedure computes the difference between
# the pointer's new x-coordinate relative to that label and the one recorded by
# the last invocation of labelB1Down, and adjusts the width of the
# corresponding column accordingly.  Otherwise a horizontal scrolling is
# performed if needed, and the would-be target position of the clicked label is
# visualized if the columns are movable.
#------------------------------------------------------------------------------
proc tablelist::labelB1Motion {w x y} {
    parseLabelPath $w win col
    upvar ::tablelist::ns${win}::data data

    if {!$data(labelClicked)} {
	return ""
    }

    if {[info exists data(x)]} {		;# resize operation in progress
	set width [expr {$data(oldStretchedColWidth) + $x - $data(x)}]
	if {$width >= $data(minColWidth)} {
	    set idx [expr {3*$col}]
	    set data(-columns) [lreplace $data(-columns) $idx $idx -$width]
	    set idx [expr {2*$col}]
	    set data(colList) [lreplace $data(colList) $idx $idx $width]
	    set data($col-lastStaticWidth) $width
	    set data($col-delta) 0
	    adjustColumns $win {} 0
	    redisplayCol $win $col [rowIndex $win @0,0 0] \
				   [rowIndex $win @0,[winfo height $win] 0]
	}
    } else {
	#
	# Scroll the window horizontally if needed
	#
	set scroll 0
	set X [expr {[winfo rootx $w] + $x}]
	set hdrX [winfo rootx $data(hdr)]
	set rightX [expr {$hdrX + [winfo width $data(hdr)]}]
	if {($X >= $rightX) &&
	    (![info exists data(X)] || $data(X) < $rightX)} {
	    set scroll 1
	} elseif {($X < $hdrX) &&
		  (![info exists data(X)] || $data(X) >= $hdrX)} {
	    set scroll 1
	}
	set data(X) $X
	if ($scroll) {
	    horizAutoScan $win
	}

	if {$x >= 0 && $x < [winfo width $w] &&
	    $y >= 0 && $y < [winfo height $w]} {
	    #
	    # The following code is needed because the event can also
	    # occur in the canvas displaying an up- or down-arrow
	    #
	    set data(inClickedLabel) 1
	    $data(hdrTxtFrCanv) configure -cursor $data(-cursor)
	    configLabel $w -cursor $data(-cursor)
	    if {$data(changeRelief)} {
		configLabel $w -relief sunken
	    }

	    place forget $data(hdrGap)
	} else {
	    #
	    # The following code is needed because the event can also
	    # occur in the canvas displaying an up- or down-arrow
	    #
	    set data(inClickedLabel) 0
	    configLabel $w -relief $data(relief)

	    if {$data(-movablecolumns)} {
		$data(hdrTxtFrCanv) configure -cursor $data(-movecolumncursor)
		configLabel $w -cursor $data(-movecolumncursor)

		#
		# Get the target column index and visualize the
		# would-be target position of the clicked label
		#
		set contW [winfo containing -displayof $w $X [winfo rooty $w]]
		parseLabelPath $contW dummy targetCol
		if {[info exists targetCol]} {
		    set master $contW
		    if {$X < [winfo rootx $contW] + [winfo width $contW]/2} {
			set relx 0.0
		    } else {
			incr targetCol
			set relx 1.0
		    }
		} elseif {[string compare $contW $data(hdrGap)] == 0} {
		    set targetCol $data(targetCol)
		    set master $data(master)
		    set relx $data(relx)
		} elseif {$X < [winfo rootx $w]} {
		    for {set targetCol 0} {$targetCol < $data(colCount)} \
			{incr targetCol} {
			if {!$data($targetCol-hide)} {
			    break
			}
		    }
		    set master $data(hdrTxtFr)
		    set relx 0.0
		} else {
		    for {set targetCol $data(lastCol)} {$targetCol >= 0} \
			{incr targetCol -1} {
			if {!$data($targetCol-hide)} {
			    break
			}
		    }
		    incr targetCol
		    set master $data(hdrTxtFr)
		    set relx 1.0
		}
		set data(targetCol) $targetCol
		set data(master) $master
		set data(relx) $relx
		$data(hdrTxtFrCanv) configure -cursor $data(-movecolumncursor)
		configLabel $w -cursor $data(-movecolumncursor)
		place $data(hdrGap) -in $master -anchor n -bordermode outside \
				    -relheight 1.0 -relx $relx
	    }
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::labelB1Enter
#
# This procedure is invoked when the mouse pointer enters the header label w of
# a tablelist widget while mouse button 1 is down.  If the label was not
# previously clicked then nothing happens.  Otherwise, if this event occured
# during a column resize operation then the procedure updates the mouse cursor
# accordingly.  Otherwise it changes the label's relief to sunken.
#------------------------------------------------------------------------------
proc tablelist::labelB1Enter w {
    parseLabelPath $w win col
    upvar ::tablelist::ns${win}::data data

    if {!$data(labelClicked)} {
	return ""
    }

    configLabel $w -cursor $data(-cursor)

    if {[info exists data(x)]} {		;# resize operation in progress
	configLabel $w -cursor $data(-resizecursor)
    } else {
	set data(inClickedLabel) 1
	if {$data(changeRelief)} {
	    configLabel $w -relief sunken
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::labelB1Leave
#
# This procedure is invoked when the mouse pointer leaves the header label w of
# a tablelist widget while mouse button 1 is down.  If the label was not
# previously clicked then nothing happens.  Otherwise, if no column resize
# operation is in progress then the procedure restores the label's relief, and,
# if the columns are movable, then it changes the mouse cursor, too.
#------------------------------------------------------------------------------
proc tablelist::labelB1Leave {w x y} {
    parseLabelPath $w win col
    upvar ::tablelist::ns${win}::data data

    if {!$data(labelClicked) ||
	[info exists data(x)]} {		;# resize operation in progress
	return ""
    }

    #
    # The following code is needed because the event can also
    # occur in the canvas displaying an up- or down-arrow
    #
    if {$x >= 0 && $x < [winfo width $w] &&
	$y >= 0 && $y < [winfo height $w]} {
	return ""
    }

    set data(inClickedLabel) 0
    configLabel $w -relief $data(relief)
}

#------------------------------------------------------------------------------
# tablelist::labelB1Up
#
# This procedure is invoked when mouse button 1 is released, if it was
# previously clicked in a label of the tablelist widget win.  If this event
# occured during a column resize operation then the procedure redisplays the
# columns and stretches the stretchable ones.  Otherwise, if the mouse button
# was released in the previously clicked label then the procedure restores the
# label's relief and invokes the command specified with the -labelcommand
# configuration option, passing to it the widget name and the column number as
# arguments.  Otherwise the column of the previously clicked label is moved
# before the column containing the mouse cursor or to its right, if the columns
# are movable.
#------------------------------------------------------------------------------
proc tablelist::labelB1Up {w X} {
    parseLabelPath $w win col
    upvar ::tablelist::ns${win}::data data

    if {!$data(labelClicked)} {
	return ""
    }

    if {[info exists data(x)]} {		;# resize operation in progress
	configLabel $w -cursor $data(-cursor)
	bind [winfo toplevel $win] <Escape> $data(topBindEsc)
	redisplayWhenIdle $win
	if {$data(-width) <= 0} {
	    $data(hdr) configure -width $data(hdrPixels)
	} elseif {[info exists data(stretchableCols)] &&
		  [lsearch -exact $data(stretchableCols) $col] >= 0} {
	    set oldColWidth \
		[expr {$data(oldStretchedColWidth) - $data(oldColDelta)}]
	    set stretchedColWidth \
		[expr {[winfo width $w] - 2*$data(charWidth)}]
	    if {$oldColWidth < $data(stretchablePixels) &&
		$stretchedColWidth < $oldColWidth + $data(delta)} {
		#
		# Compute the new column width, using the following equations:
		#
		# $stretchedColWidth = $colWidth + $colDelta
		# $colDelta =
		#    ($data(delta) - $colWidth + $oldColWidth) * $colWidth /
		#    ($data(stretchablePixels) + $colWidth - $oldColWidth)
		#
		set colWidth [expr {
		    $stretchedColWidth *
		    ($data(stretchablePixels) - $oldColWidth) /
		    ($data(stretchablePixels) + $data(delta) -
		     $stretchedColWidth)
		}]
		if {$colWidth < 1} {
		    set colWidth 1
		}
		set idx [expr {3*$col}]
		set data(-columns) \
		    [lreplace $data(-columns) $idx $idx -$colWidth]
		set idx [expr {2*$col}]
		set data(colList) [lreplace $data(colList) $idx $idx $colWidth]
		set data($col-delta) [expr {$stretchedColWidth - $colWidth}]
	    }
	}
	stretchColumns $win $col
	unset data(x)
    } else {
	if {[info exists data(X)]} {
	    unset data(X)
	}
    	if {$data(-movablecolumns)} {
	    bind [winfo toplevel $win] <Escape> $data(topBindEsc)
	    place forget $data(hdrGap)
	}
	if {$data(inClickedLabel)} {
	    configLabel $w -relief $data(relief)
	    if {[info exists data($col-labelcommand)]} {
		uplevel #0 $data($col-labelcommand) [list $win $col]
	    } elseif {[string compare $data(-labelcommand) ""] != 0} {
		uplevel #0 $data(-labelcommand) [list $win $col]
	    }
	} elseif {$data(-movablecolumns)} {
	    $data(hdrTxtFrCanv) configure -cursor $data(-cursor)
	    if {$data(targetCol) != $col && $data(targetCol) != $col + 1} {
		movecolumnSubCmd $win $col $data(targetCol)
	    }
	}
    }

    set data(labelClicked) 0
}

#------------------------------------------------------------------------------
# tablelist::labelB3Down
#
# This procedure is invoked when mouse button 3 is pressed in the header label
# w of a tablelist widget.  It configures the width of the given column to be
# just large enough to hold all the elements (including the label).
#------------------------------------------------------------------------------
proc tablelist::labelB3Down w {
    parseLabelPath $w win col
    upvar ::tablelist::ns${win}::data data

    if {!$data(isDisabled) &&
	$data(-resizablecolumns) && $data($col-resizable)} {
	doColConfig $col $win -width 0
    }
}

#------------------------------------------------------------------------------
# tablelist::labelShiftB3Down
#
# This procedure is invoked when mouse button 3 together with the Shift key is
# pressed in the header label w of a tablelist widget.  It restores the last
# static width of the given column.
#------------------------------------------------------------------------------
proc tablelist::labelShiftB3Down w {
    parseLabelPath $w win col
    upvar ::tablelist::ns${win}::data data

    if {!$data(isDisabled) &&
	$data(-resizablecolumns) && $data($col-resizable)} {
	doColConfig $col $win -width -$data($col-lastStaticWidth)
    }
}

#------------------------------------------------------------------------------
# tablelist::escape
#
# This procedure is invoked to process <Escape> events in the top-level window
# containing the tablelist widget win during a column resize or move operation.
# The procedure cancels the action in progress and, in case of column resizing,
# it restores the initial width of the respective column.
#------------------------------------------------------------------------------
proc tablelist::escape {win col} {
    upvar ::tablelist::ns${win}::data data

    set w $data(hdrTxtFrLbl)$col
    if {[info exists data(x)]} {		;# resize operation in progress
	configLabel $w -cursor $data(-cursor)
	update idletasks
	bind [winfo toplevel $win] <Escape> $data(topBindEsc)
	set idx [expr {3*$col}]
	setupColumns $win [lreplace $data(-columns) $idx $idx \
				    $data(configColWidth)] 0
	adjustColumns $win $col 1
	redisplayCol $win $col [rowIndex $win @0,0 0] \
			       [rowIndex $win @0,[winfo height $win] 0]
	unset data(x)
	set data(labelClicked) 0
    } elseif {!$data(inClickedLabel)} {
	configLabel $w -cursor $data(-cursor)
	$data(hdrTxtFrCanv) configure -cursor $data(-cursor)
	bind [winfo toplevel $win] <Escape> $data(topBindEsc)
	place forget $data(hdrGap)
	if {[info exists data(X)]} {
	    unset data(X)
	}
	set data(labelClicked) 0
    }
}

#------------------------------------------------------------------------------
# tablelist::parseLabelPath
#
# Extracts the path name of the tablelist widget as well as the column number
# from the path name w of a header label.
#------------------------------------------------------------------------------
proc tablelist::parseLabelPath {w winName colName} {
    upvar $winName win $colName col

    regexp {^(.+)\.hdr\.t\.f\.l([0-9]+)$} $w dummy win col
}

#------------------------------------------------------------------------------
# tablelist::redisplayCol
#
# Redisplays the elements of the col'th column of the tablelist widget win, in
# the range specified by first and last.
#------------------------------------------------------------------------------
proc tablelist::redisplayCol {win col first last} {
    upvar ::tablelist::ns${win}::data data

    if {$data($col-hide) || $first < 0} {
	return ""
    }

    set snipStr $data(-snipstring)
    set fmtCmdFlag [info exists data($col-formatcommand)]
    set colFont [lindex $data(colFontList) $col]

    set w $data(body)
    set pixels [lindex $data(colList) [expr {2*$col}]]
    if {$pixels == 0} {				;# convention: dynamic width
	if {$data($col-maxPixels) > 0 &&
	    $data($col-reqPixels) > $data($col-maxPixels)} {
	    set pixels $data($col-maxPixels)
	}
    }
    if {$pixels != 0} {
	incr pixels $data($col-delta)
    }
    set alignment [lindex $data(colList) [expr {2*$col + 1}]]

    for {set idx $first; set line [expr {$first + 1}]} {$idx <= $last} \
	{incr idx; incr line} {
	if {$idx == $data(editRow) && $col == $data(editCol)} {
	    continue
	}

	#
	# Adjust the cell text and the image width
	#
	set item [lindex $data(itemList) $idx]
	set text [lindex $item $col]
	if {$fmtCmdFlag} {
	    set text [uplevel #0 $data($col-formatcommand) [list $text]]
	}
	set text [strToDispStr $text]
	set key [lindex $item end]
	if {[info exists data($key-$col-image)]} {
	    set image $data($key-$col-image)
	    set imageWidth [image width $image]
	} else {
	    set image ""
	    set imageWidth 0
	}
	if {[info exists data($key-$col-font)]} {
	    set cellFont $data($key-$col-font)
	} elseif {[info exists data($key-font)]} {
	    set cellFont $data($key-font)
	} else {
	    set cellFont $colFont
	}
	adjustElem $win text imageWidth $cellFont \
		   $pixels $alignment $snipStr

	#
	# Delete the old cell contents between the
	# two tabs, and insert the text and the image
	#
	findCellTabs $win $line $col tabIdx1 tabIdx2
	$w delete $tabIdx1+1c $tabIdx2
	insertElem $w $tabIdx1+1c $text $image $imageWidth $alignment
    }
}

#------------------------------------------------------------------------------
# tablelist::horizAutoScan
#
# This procedure is invoked when the mouse leaves the header frame of a
# tablelist widget.  It scrolls the child and reschedules itself as an after
# command so that the child continues to scroll until the mouse moves back into
# the window or the mouse button is released.
#------------------------------------------------------------------------------
proc tablelist::horizAutoScan win {
    if {![winfo exists $win]} {
	return ""
    }

    upvar ::tablelist::ns${win}::data data
    if {![info exists data(X)]} {
	return ""
    }

    set X $data(X)
    set hdrX [winfo rootx $data(hdr)]
    set rightX [expr {$hdrX + [winfo width $data(hdr)]}]

    if {$X >= $rightX} {
	xviewSubCmd $win {scroll 2 units}
    } elseif {$X < $hdrX} {
	xviewSubCmd $win {scroll -2 units}
    } else {
	return ""
    }

    after 50 [list tablelist::horizAutoScan $win]
}
