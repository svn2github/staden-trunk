#==============================================================================
# Contains the implementation of the tablelist::sortByColumn command as well as
# of the tablelist sort and sortbycolumn subcommands.
#
# Copyright (c) 2000-2004  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

#------------------------------------------------------------------------------
# tablelist::sortByColumn
#
# Sorts the contents of the tablelist widget win by its col'th column.  Returns
# the sorting order (increasing or decreasing).
#------------------------------------------------------------------------------
proc tablelist::sortByColumn {win col} {
    #
    # Check the arguments
    #
    if {![winfo exists $win]} {
	return -code error "bad window path name \"$win\""
    }
    if {[string compare [winfo class $win] Tablelist] != 0} {
	return -code error "window \"$win\" is not a tablelist widget"
    }
    if {[catch {::$win columnindex $col} result] != 0} {
	return -code error $result
    }
    if {$result < 0 || $result >= [::$win columncount]} {
	return -code error "column index \"$col\" out of range"
    }

    #
    # Determine the sorting order
    #
    set col $result
    if {$col == [::$win sortcolumn] &&
	[string compare [::$win sortorder] increasing] == 0} {
	set sortOrder decreasing
    } else {
	set sortOrder increasing
    }

    #
    # Sort the widget's contents based on the given column
    #
    if {[catch {::$win sortbycolumn $col -$sortOrder} result] == 0} {
	return $sortOrder
    } else {
	return -code error $result
    }
}

#------------------------------------------------------------------------------
# tablelist::sortSubCmd
#
# This procedure is invoked to process the tablelist sort and sortbycolumn
# subcommands.
#------------------------------------------------------------------------------
proc tablelist::sortSubCmd {win col order} {
    upvar ::tablelist::ns${win}::data data

    #
    # Save the keys corresponding to anchorIdx and activeIdx 
    #
    foreach type {anchor active} {
	set item [lindex $data(itemList) $data(${type}Idx)]
	set ${type}Key [lindex $item end]
    }

    #
    # Save the keys of the selected items
    #
    set selKeys {}
    set w $data(body)
    set selRange [$w tag nextrange select 1.0]
    while {[llength $selRange] != 0} {
	set selStart [lindex $selRange 0]
	set selEnd [lindex $selRange 1]
	set item [lindex $data(itemList) [expr {int($selStart) - 1}]]
	lappend selKeys [lindex $item end]

	set selRange [$w tag nextrange select $selEnd]
    }

    #
    # Save some data of the edit window if present
    #
    if {[set editCol $data(editCol)] >= 0} {
	set editKey $data(editKey)
	saveEditData $win
    }

    #
    # Sort the item list and update the sort info
    #
    if {$col < 0} {				;# not sorting by a column
	if {[string compare $data(-sortcommand) ""] == 0} {
	    return -code error \
		   "value of the -sortcommand option is empty"
	}

	set data(itemList) \
	    [lsort $order -command $data(-sortcommand) $data(itemList)]
    } else {					;# sorting by a column
	if {[string compare $data($col-sortmode) command] == 0} {
	    if {[info exists data($col-sortcommand)]} {
		set data(itemList) \
		    [lsort $order -index $col \
		     -command $data($col-sortcommand) $data(itemList)]
	    } else {
		return -code error \
		       "value of the -sortcommand option for\
			column $col is missing or empty"
	    }
	} else {
	    set data(itemList) \
		[lsort $order -index $col \
		 -$data($col-sortmode) $data(itemList)]
	}
    }
    set data(sortCol) $col
    set data(sortOrder) [string range $order 1 end]

    #
    # Replace the contents of the list variable if present
    #
    if {$data(hasListVar)} {
	trace vdelete ::$data(-listvariable) wu $data(listVarTraceCmd)
	upvar #0 $data(-listvariable) var
	set var {}
	foreach item $data(itemList) {
	    lappend var [lrange $item 0 $data(lastCol)]
	}
	trace variable ::$data(-listvariable) wu $data(listVarTraceCmd)
    }

    #
    # Update anchorIdx and activeIdx
    #
    foreach type {anchor active} {
	upvar 0 ${type}Key key
	if {[string compare $key ""] != 0} {
	    set data(${type}Idx) [lsearch $data(itemList) "* $key"]
	}
    }

    #
    # Unmanage the canvas and adjust the columns
    #
    set canvas $data(hdrTxtFrCanv)
    place forget $canvas
    set oldArrowCol $data(arrowCol)
    set data(arrowCol) -1
    adjustColumns $win l$oldArrowCol 1

    #
    # Check whether an up- or down-arrow is to be displayed
    #
    if {$col >= 0 && $data(-showarrow) && $data($col-showarrow)} {
	#
	# Configure the canvas and draw the arrows
	#
	set data(arrowCol) $col
	configCanvas $win
	drawArrows $win

	#
	# Make sure the arrow will fit into the column
	#
	set idx [expr {2*$col}]
	set pixels [lindex $data(colList) $idx]
	if {$pixels == 0 && $data($col-maxPixels) > 0 &&
	    $data($col-reqPixels) > $data($col-maxPixels) &&
	    $data($col-maxPixels) < $data(arrowWidth)} {
	    set data($col-maxPixels) $data(arrowWidth)
	    set data($col-maxwidth) -$data(arrowWidth)
	}
	if {$pixels != 0 && $pixels < $data(arrowWidth)} {
	    set data(colList) \
		[lreplace $data(colList) $idx $idx $data(arrowWidth)]
	    set idx [expr {3*$col}]
	    set data(-columns) \
		[lreplace $data(-columns) $idx $idx -$data(arrowWidth)]
	}

	#
	# Adjust the columns; this will also place the canvas into the label
	#
	adjustColumns $win l$col 1
    }

    #
    # Delete the items from the body text widget and insert the sorted ones.
    # Interestingly, for a large number of items it is much more efficient
    # to empty each line individually than to invoke a global delete command.
    #
    set widgetFont $data(-font)
    set snipStr $data(-snipstring)
    set isSimple [expr {$data(tagCount) == 0 && $data(imgCount) == 0 &&
			!$data(hasColTags)}]
    set line 1
    foreach item $data(itemList) {
	#
	# Empty the line, clip the elements if necessary,
	# and insert them with the corresponding tags
	#
	$w delete $line.0 $line.end
	set dispItem [strToDispStr $item]
	set col 0
	if {$isSimple} {
	    set insertStr ""
	    foreach text [lrange $dispItem 0 $data(lastCol)] \
		    fmtCmdFlag $data(fmtCmdFlagList) \
		    {pixels alignment} $data(colList) {
		if {$data($col-hide)} {
		    incr col
		    continue
		}

		#
		# Clip the element if necessary
		#
		if {$fmtCmdFlag} {
		    set text [uplevel #0 $data($col-formatcommand) \
			      [list [lindex $item $col]]]
		    set text [strToDispStr $text]
		}
		if {$pixels == 0} {		;# convention: dynamic width
		    if {$data($col-maxPixels) > 0 &&
			$data($col-reqPixels) > $data($col-maxPixels)} {
			set pixels $data($col-maxPixels)
		    }
		}
		if {$pixels != 0} {
		    incr pixels $data($col-delta)
		    set text [strRangeExt $win $text $widgetFont \
			      $pixels $alignment $snipStr]
		}

		append insertStr \t$text\t
		incr col
	    }

	    #
	    # Insert the item into the body text widget
	    #
	    $w insert $line.0 $insertStr

	} else {
	    set key [lindex $item end]
	    array set tagData [array get data $key*-\[bf\]*]	;# for speed

	    set rowTags [array names tagData $key-\[bf\]*]
	    foreach text [lrange $dispItem 0 $data(lastCol)] \
		    colFont $data(colFontList) \
		    colTags $data(colTagsList) \
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
		if {[info exists data($key-$col-font)]} {
		    set cellFont $data($key-$col-font)
		} elseif {[info exists data($key-font)]} {
		    set cellFont $data($key-font)
		} else {
		    set cellFont $colFont
		}
		if {$pixels == 0} {		;# convention: dynamic width
		    if {$data($col-maxPixels) > 0 &&
			$data($col-reqPixels) > $data($col-maxPixels)} {
			set pixels $data($col-maxPixels)
		    }
		}
		if {$pixels != 0} {
		    incr pixels $data($col-delta)
		}
		adjustElem $win text imageWidth $cellFont \
			   $pixels $alignment $snipStr

		#
		# Insert the text and the image
		#
		set cellTags [array names tagData $key-$col-\[bf\]*]
		set tagNames [concat $colTags $rowTags $cellTags]
		if {$imageWidth == 0} {
		    $w insert $line.end \t$text\t $tagNames
		} else {
		    $w insert $line.end \t\t $tagNames
		    insertElem $w $line.end-1c $text $image $imageWidth \
			       $alignment
		}

		incr col
	    }

	    unset tagData
	}

	incr line
    }

    #
    # Restore the stripes in the body text widget
    #
    makeStripes $win

    #
    # Select the items that were selected before
    #
    foreach key $selKeys {
	set idx [lsearch $data(itemList) "* $key"]
	selectionSubCmd $win set $idx $idx
    }

    #
    # Restore the edit window if it was present before
    #
    if {$editCol >= 0} {
	set editRow [lsearch $data(itemList) "* $editKey"]
	editcellSubCmd $win $editRow $editCol 1
    }

    #
    # Disable the body text widget if it was disabled before
    #
    if {$data(isDisabled)} {
	$w tag add disabled 1.0 end
	$w tag configure select -borderwidth 0
    }

    return ""
}
