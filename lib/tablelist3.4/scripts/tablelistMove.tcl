#==============================================================================
# Contains the implementation of the tablelist move and movecolumn subcommands.
#
# Copyright (c) 2003-2004  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

#------------------------------------------------------------------------------
# tablelist::moveSubCmd
#
# This procedure is invoked to process the tablelist move subcommand.
#------------------------------------------------------------------------------
proc tablelist::moveSubCmd {win source target} {
    upvar ::tablelist::ns${win}::data data

    if {$data(isDisabled) || $data(itemCount) == 0} {
	return ""
    }

    #
    # Adjust the indices to fit within the existing items and check them
    #
    if {$source > $data(lastRow)} {
	set source $data(lastRow)
    } elseif {$source < 0} {
	set source 0
    }
    if {$target > $data(itemCount)} {
	set target $data(itemCount)
    } elseif {$target < 0} {
	set target 0
    }
    if {$target == $source} {
	return -code error \
	       "cannot move item with index \"$source\" before itself"
    } elseif {$target == $source + 1} {
	return ""
    }

    #
    # Save some data of the edit window if present
    #
    if {[set editCol $data(editCol)] >= 0} {
	set editRow $data(editRow)
	set editKey $data(editKey)
	saveEditData $win
    }

    #
    # Check whether the source line is selected and then delete it
    #
    set w $data(body)
    set textIdx1 [expr {double($source + 1)}]
    set textIdx2 [expr {double($source + 2)}]
    if {[lsearch -exact [$w tag names $textIdx1] select] >= 0} {
	set selected 1
    } else {
	set selected 0
    }
    $w delete $textIdx1 $textIdx2

    #
    # Insert the source item before the target one
    #
    set target1 $target
    if {$source < $target} {
	incr target1 -1
    }
    set targetLine [expr {$target1 + 1}]
    $w insert $targetLine.0 \n
    set snipStr $data(-snipstring)
    set sourceItem [lindex $data(itemList) $source]
    set dispItem [strToDispStr $sourceItem]
    set key [lindex $sourceItem end]
    array set tagData [array get data $key*-\[bf\]*]		;# for speed
    set rowTags [array names tagData $key-\[bf\]*]
    set col 0
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
		      [list [lindex $sourceItem $col]]]
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
	if {$pixels == 0} {			;# convention: dynamic width
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
	    $w insert $targetLine.end \t$text\t $tagNames
	} else {
	    $w insert $targetLine.end \t\t $tagNames
	    insertElem $w $targetLine.end-1c $text $image $imageWidth $alignment
	}

	incr col
    }

    #
    # Update the item list
    #
    set data(itemList) [lreplace $data(itemList) $source $source]
    if {$target == $data(itemCount)} {
	lappend data(itemList) $sourceItem	;# this works much faster
    } else {
	set data(itemList) [linsert $data(itemList) $target1 $sourceItem]
    }

    #
    # Update the list variable if present
    #
    if {$data(hasListVar)} {
	trace vdelete ::$data(-listvariable) wu $data(listVarTraceCmd)
	upvar #0 $data(-listvariable) var
	set var [lreplace $var $source $source]
	set pureSourceItem [lrange $sourceItem 0 $data(lastCol)]
	if {$target == $data(itemCount)} {
	    lappend var $pureSourceItem		;# this works much faster
	} else {
	    set var [linsert $var $target1 $pureSourceItem]
	}
	trace variable ::$data(-listvariable) wu $data(listVarTraceCmd)
    }

    #
    # Update anchorIdx and activeIdx if needed
    #
    if {$data(anchorIdx) == $source} {
	set data(anchorIdx) $target1
    }
    if {$data(activeIdx) == $source} {
	set data(activeIdx) $target1
    }

    #
    # Restore the stripes in the body text widget
    #
    makeStripesWhenIdle $win

    #
    # Select the source item if it was selected before
    #
    if {$selected} {
	selectionSubCmd $win set $target1 $target1
    }

    #
    # Restore the edit window if it was present before
    #
    if {$editCol >= 0} {
	if {$editRow == $source} {
	    editcellSubCmd $win $target1 $editCol 1
	} else {
	    set data(editRow) [lsearch $data(itemList) "* $editKey"]
	}
    }

    return ""
}

#------------------------------------------------------------------------------
# tablelist::movecolumnSubCmd
#
# This procedure is invoked to process the tablelist movecolumn subcommand.
#------------------------------------------------------------------------------
proc tablelist::movecolumnSubCmd {win source target} {
    upvar ::tablelist::ns${win}::data data

    if {$data(isDisabled)} {
	return ""
    }

    #
    # Check the indices
    #
    if {$target == $source} {
	return -code error \
	       "cannot move column with index \"$source\" before itself"
    } elseif {$target == $source + 1} {
	return ""
    }

    #
    # Update the column list
    #
    set source3 [expr {3*$source}]
    set source3Plus2 [expr {$source3 + 2}]
    set target1 $target
    set target3 [expr {3*$target}]
    if {$source < $target} {
	incr target1 -1
	incr target3 -3
    }
    set sourceRange [lrange $data(-columns) $source3 $source3Plus2]
    set data(-columns) [lreplace $data(-columns) $source3 $source3Plus2]
    set data(-columns) [eval linsert {$data(-columns)} $target3 $sourceRange]

    #
    # Save some elements of data corresponding to source
    #
    array set tmp [array get data $source-*]
    array set tmp [array get data k*-$source-*]
    set tmp(editCol) $data(editCol)
    set tmp(arrowCol) $data(arrowCol)
    set tmp(sortCol) $data(sortCol)

    #
    # Remove source from the list of stretchable columns
    # if it was explicitly specified as stretchable
    #
    if {[string first $data(-stretch) all] != 0} {
	set sourceIsStretchable 0
	set stretchableCols {}
	foreach elem $data(-stretch) {
	    if {[string first $elem end] != 0 && $elem == $source} {
		set sourceIsStretchable 1
	    } else {
		lappend stretchableCols $elem
	    }
	}
	set data(-stretch) $stretchableCols
    }

    #
    # Build two lists of column numbers, neeeded
    # for shifting some elements of the data array
    #
    if {$source < $target} {
	for {set n $source} {$n < $target1} {incr n} {
	    lappend oldCols [expr {$n + 1}]
	    lappend newCols $n
	}
    } else {
	for {set n $source} {$n > $target} {incr n -1} {
	    lappend oldCols [expr {$n - 1}]
	    lappend newCols $n
	}
    }

    #
    # Move the elements of data corresponding to the columns in oldCols to the
    # elements corresponding to the columns with the same indices in newCols
    #
    foreach oldCol $oldCols newCol $newCols {
	moveColData $win data data imgs $oldCol $newCol
    }

    #
    # Move the elements of data corresponding to
    # source to the elements corresponding to target1
    #
    moveColData $win tmp data imgs $source $target1

    #
    # If the column given by source was explicitly specified as
    # stretchable then add target1 to the list of stretchable columns
    #
    if {[string first $data(-stretch) all] != 0 && $sourceIsStretchable} {
	lappend data(-stretch) $target1
	sortStretchableColList $win
    }

    #
    # Update the item list
    #
    set newItemList {}
    foreach item $data(itemList) {
	set sourceText [lindex $item $source]
	set item [lreplace $item $source $source]
	set item [linsert $item $target1 $sourceText]
	lappend newItemList $item
    }
    set data(itemList) $newItemList

    #
    # Update the list variable if present
    #
    condUpdateListVar $win

    #
    # Set up and adjust the columns, and rebuild
    # the lists of the column fonts and tag names
    #
    setupColumns $win $data(-columns) 0
    makeColFontAndTagLists $win
    adjustColumns $win {} 0

    #
    # Reconfigure the relevant column labels
    #
    foreach col [lappend newCols $target1] {
	reconfigColLabels $win imgs $col
    }

    #
    # Make sure the items will be redisplayed at idle time
    #
    redisplayWhenIdle $win
    return ""
}
