#==============================================================================
# Contains private utility procedures for tablelist widgets.
#
# Copyright (c) 2000-2004  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

#------------------------------------------------------------------------------
# tablelist::rowIndex
#
# Checks the row index idx and returns either its numerical value or an error.
# endIsSize must be a boolean value: if true, end refers to the number of items
# in the tablelist, i.e., to the element just after the last one; if false, end
# refers to 1 less than the number of items, i.e., to the last element in the
# tablelist.
#------------------------------------------------------------------------------
proc tablelist::rowIndex {win idx endIsSize} {
    upvar ::tablelist::ns${win}::data data

    set idxLen [string length $idx]
    if {[string first $idx active] == 0 && $idxLen >= 2} {
	return $data(activeIdx)
    } elseif {[string first $idx anchor] == 0 && $idxLen >= 2} {
	return $data(anchorIdx)
    } elseif {[string first $idx end] == 0} {
	if {$endIsSize} {
	    return $data(itemCount)
	} else {
	    return $data(lastRow)
	}
    } elseif {[string compare [string index $idx 0] @] == 0} {
	if {[catch {$data(body) index $idx}] == 0} {
	    if {$data(itemCount) == 0} {
		return -1
	    } else {
		scan $idx @%d,%d x y
		incr x -[winfo x $data(body)]
		incr y -[winfo y $data(body)]
		set textIdx [$data(body) index @$x,$y]
		return [expr {int($textIdx) - 1}]
	    }
	} else {
	    return -code error \
		   "bad row index \"$idx\": must be active, anchor,\
		    end, @x,y, a number, or a full key"
	}
    } elseif {[string compare [string index $idx 0] k] == 0} {
	if {[set index [lsearch $data(itemList) "* $idx"]] >= 0} {
	    return $index
	} else {
	    return -code error \
		   "bad row index \"$idx\": must be active, anchor,\
		    end, @x,y, a number, or a full key"
	}
    } elseif {[catch {format %d $idx} index] == 0} {
	return $index
    } else {
	return -code error \
	       "bad row index \"$idx\": must be active, anchor,\
	        end, @x,y, a number, or a full key"
    }
}

#------------------------------------------------------------------------------
# tablelist::colIndex
#
# Checks the column index idx and returns either its numerical value or an
# error.  checkRange must be a boolean value: if true, it is additionally
# checked whether the numerical value corresponding to idx is within the
# allowed range.
#------------------------------------------------------------------------------
proc tablelist::colIndex {win idx checkRange} {
    upvar ::tablelist::ns${win}::data data

    if {[string first $idx end] == 0} {
	return $data(lastCol)
    } elseif {[string compare [string index $idx 0] @] == 0 &&
	      [catch {$data(body) index $idx}] == 0} {
	scan $idx @%d x
	incr x -[winfo x $data(body)]
	set bodyWidth [winfo width $data(body)]
	if {$x >= $bodyWidth} {
	    set x [expr {$bodyWidth - 1}]
	} elseif {$x < 0} {
	    set x 0
	}
	set x [expr {$x + [winfo rootx $data(body)]}]

	set lastVisibleCol -1
	for {set col 0} {$col < $data(colCount)} {incr col} {
	    if {$data($col-hide)} {
		continue
	    }

	    set lastVisibleCol $col
	    set w $data(hdrTxtFrLbl)$col
	    set wX [winfo rootx $w]
	    if {$x >= $wX && $x < $wX + [winfo width $w]} {
		return $col
	    }
	}
	set index $lastVisibleCol
    } elseif {[catch {format %d $idx} index] != 0} {
	for {set col 0} {$col < $data(colCount)} {incr col} {
	    set hasName [info exists data($col-name)]
	    if {$hasName && [string compare $idx $data($col-name)] == 0 ||
		!$hasName && [string compare $idx ""] == 0} {
		set index $col
		break
	    }
	}
	if {$col == $data(colCount)} {
	    return -code error \
		   "bad column index \"$idx\": must be\
		    end, @x,y, a number, or a name"
	}
    }

    if {$checkRange && ($index < 0 || $index > $data(lastCol))} {
	return -code error "column index \"$idx\" out of range"
    } else {
	return $index
    }
}

#------------------------------------------------------------------------------
# tablelist::cellIndex
#
# Checks the cell index idx and returns either its value in the form row,col or
# an error.  checkRange must be a boolean value: if true, it is additionally
# checked whether the two numerical values corresponding to idx are within the
# respective allowed ranges.
#------------------------------------------------------------------------------
proc tablelist::cellIndex {win idx checkRange} {
    upvar ::tablelist::ns${win}::data data

    if {[string first $idx end] == 0} {
	set row [rowIndex $win $idx 0]
	set col [colIndex $win $idx 0]
    } elseif {[string compare [string index $idx 0] @] == 0} {
	if {[catch {rowIndex $win $idx 0} row] != 0 ||
	    [catch {colIndex $win $idx 0} col] != 0} {
	    return -code error \
		   "bad cell index \"$idx\": must be end, @x,y, or row,col,\
		    where row must be active, anchor, end, a number, or\
		    a full key, and col must be end, a number, or a name"
	}
    } else {
	set lst [split $idx ,]
	if {[llength $lst] != 2 ||
	    [catch {rowIndex $win [lindex $lst 0] 0} row] != 0 ||
	    [catch {colIndex $win [lindex $lst 1] 0} col] != 0} {
	    return -code error \
		   "bad cell index \"$idx\": must be end, @x,y, or row,col,\
		    where row must be active, anchor, end, a number, or\
		    a full key, and col must be end, a number, or a name"
	}
    }

    if {$checkRange && ($row < 0 || $row > $data(lastRow) || \
	$col < 0 || $col > $data(lastCol))} {
	return -code error "cell index \"$idx\" out of range"
    } else {
	return $row,$col
    }
}

#------------------------------------------------------------------------------
# tablelist::findCellTabs
#
# Searches for the tab characters within the col'th cell in the given line of
# the body text child of the tablelist widget win.  Assigns the index of the
# first tab to $idx1Name and the index of the second tab to $idx2Name.
#------------------------------------------------------------------------------
proc tablelist::findCellTabs {win line col idx1Name idx2Name} {
    upvar ::tablelist::ns${win}::data data
    upvar $idx1Name idx1 $idx2Name idx2

    set w $data(body)
    set idx1 $line.0
    set endIdx $line.end
    for {set n 0} {$n < $col} {incr n} {
	if {!$data($n-hide)} {
	    set idx1 [$w search \t $idx1+1c $endIdx]+1c
	}
    }
    set idx1 [$w index $idx1]
    set idx2 [$w search \t $idx1+1c $endIdx]

    return ""
}

#------------------------------------------------------------------------------
# tablelist::cellFont
#
# Returns the font to be used in the tablelist cell specified by win, key, and
# col.
#------------------------------------------------------------------------------
proc tablelist::cellFont {win key col} {
    upvar ::tablelist::ns${win}::data data

    if {[info exists data($key-$col-font)]} {
	return $data($key-$col-font)
    } elseif {[info exists data($key-font)]} {
	return $data($key-font)
    } else {
	return [lindex $data(colFontList) $col]
    }
}

#------------------------------------------------------------------------------
# tablelist::sortStretchableColList
#
# Replaces the column indices different from end in the list of the stretchable
# columns of the tablelist widget win with their numerical equivalents and
# sorts the resulting list.
#------------------------------------------------------------------------------
proc tablelist::sortStretchableColList win {
    upvar ::tablelist::ns${win}::data data

    if {[llength $data(-stretch)] == 0 ||
	[string first $data(-stretch) all] == 0} {
	return ""
    }

    set containsEnd 0
    foreach elem $data(-stretch) {
	if {[string first $elem end] == 0} {
	    set containsEnd 1
	} else {
	    set tmp([colIndex $win $elem 0]) ""
	}
    }

    set data(-stretch) [lsort -integer [array names tmp]]
    if {$containsEnd} {
	lappend data(-stretch) end
    }
}

#------------------------------------------------------------------------------
# tablelist::deleteColData
#
# Cleans up the data associated with the col'th column of the tablelist widget
# win.
#------------------------------------------------------------------------------
proc tablelist::deleteColData {win col} {
    upvar ::tablelist::ns${win}::data data

    if {$data(editCol) == $col} {
	set data(editCol) -1
	set data(editRow) -1
    }
    if {$data(arrowCol) == $col} {
	set data(arrowCol) -1
    }
    if {$data(sortCol) == $col} {
	set data(sortCol) -1
    }

    #
    # Remove the elements with names of the form $col-*
    #
    set w $data(body)
    foreach name [array names data $col-*] {
	unset data($name)
	$w tag delete $name
    }

    #
    # Remove the elements with names of the form k*-$col-*
    #
    foreach name [array names data k*-$col-*] {
	unset data($name)
	$w tag delete $name
	if {[string match k*-$col-\[bf\]* $name]} {
	    incr data(tagCount) -1
	} elseif {[string match k*-$col-image $name]} {
	    incr data(imgCount) -1
	}
    }

    #
    # Remove col from the list of stretchable columns if explicitly specified
    #
    if {[string first $data(-stretch) all] != 0} {
	set stretchableCols {}
	foreach elem $data(-stretch) {
	    if {[string first $elem end] == 0 || $elem != $col} {
		lappend stretchableCols $elem
	    }
	}
	set data(-stretch) $stretchableCols
    }
}

#------------------------------------------------------------------------------
# tablelist::moveColData
#
# Moves the elements of oldArrName corresponding to oldCol to those of
# newArrName corresponding to newCol.
#------------------------------------------------------------------------------
proc tablelist::moveColData {win oldArrName newArrName imgArrName
			     oldCol newCol} {
    upvar $oldArrName oldArr $newArrName newArr $imgArrName imgArr

    #
    # Move the elements of oldArr with names of the form $oldCol-*
    # to those of newArr with names of the form $newCol-*
    #
    set w $newArr(body)
    foreach newName [array names newArr $newCol-*] {
	unset newArr($newName)
	$w tag delete $newName
    }
    if {$newCol < $newArr(colCount)} {
	foreach c [winfo children $newArr(hdrTxtFrLbl)$newCol] {
	    destroy $c
	}
	set newArr(fmtCmdFlagList) \
	    [lreplace $newArr(fmtCmdFlagList) $newCol $newCol 0]
    }
    if {$oldArr(editCol) == $oldCol} {
	set newArr(editCol) $newCol
    }
    if {$oldArr(arrowCol) == $oldCol} {
	set newArr(arrowCol) $newCol
    }
    if {$oldArr(sortCol) == $oldCol} {
	set newArr(sortCol) $newCol
    }
    foreach oldName [array names oldArr $oldCol-*] {
	regsub $oldCol- $oldName $newCol- newName
	set newArr($newName) $oldArr($oldName)

	unset oldArr($oldName)
	$w tag delete $oldName

	set tail [lindex [split $newName -] 1]
	switch $tail {
	    background -
	    foreground {
		$w tag configure $newName -$tail $newArr($newName)
		$w tag raise $newName stripe
	    }
	    font {
		$w tag configure $newName -$tail $newArr($newName)
		$w tag lower $newName
	    }
	    formatcommand {
		if {$newCol < $newArr(colCount)} {
		    set newArr(fmtCmdFlagList) \
			[lreplace $newArr(fmtCmdFlagList) $newCol $newCol 1]
		}
	    }
	    labelimage {
		set imgArr($newCol-$tail) $newArr($newName)
		unset newArr($newName)
	    }
	    selectbackground -
	    selectforeground {
		set tail [string range $tail 6 end]	;# remove the select
		$w tag configure $newName -$tail $newArr($newName)
		$w tag raise $newName select
	    }
	}
    }

    #
    # Move the elements of oldArr with names of the form k*-$oldCol-*
    # to those of newArr with names of the form k*-$newCol-*
    #
    foreach newName [array names newArr k*-$newCol-*] {
	unset newArr($newName)
	$w tag delete $newName
    }
    foreach oldName [array names oldArr k*-$oldCol-*] {
	regsub -- -$oldCol- $oldName -$newCol- newName
	set newArr($newName) $oldArr($oldName)

	unset oldArr($oldName)
	$w tag delete $oldName

	set tail [lindex [split $newName -] 2]
	switch $tail {
	    background -
	    foreground {
		$w tag configure $newName -$tail $newArr($newName)
		$w tag lower $newName disabled
	    }
	    font {
		$w tag configure $newName -$tail $newArr($newName)
		$w tag raise $newName
	    }
	    selectbackground -
	    selectforeground {
		set tail [string range $tail 6 end]	;# remove the select
		$w tag configure $newName -$tail $newArr($newName)
		$w tag lower $newName disabled
	    }
	}
    }

    #
    # Replace oldCol with newCol in the list of
    # stretchable columns if explicitly specified
    #
    if {[info exists oldArr(-stretch)] &&
	[string first $oldArr(-stretch) all] != 0} {
	set stretchableCols {}
	foreach elem $oldArr(-stretch) {
	    if {[string first $elem end] != 0 && $elem == $oldCol} {
		lappend stretchableCols $newCol
	    } else {
		lappend stretchableCols $elem
	    }
	}
	set newArr(-stretch) $stretchableCols
    }
}

#------------------------------------------------------------------------------
# tablelist::condUpdateListVar
#
# Updates the list variable of the tablelist widget win if present.
#------------------------------------------------------------------------------
proc tablelist::condUpdateListVar win {
    upvar ::tablelist::ns${win}::data data

    if {$data(hasListVar)} {
	trace vdelete ::$data(-listvariable) wu $data(listVarTraceCmd)
	upvar #0 $data(-listvariable) var
	set var {}
	foreach item $data(itemList) {
	    lappend var [lrange $item 0 $data(lastCol)]
	}
	trace variable ::$data(-listvariable) wu $data(listVarTraceCmd)
    }
}

#------------------------------------------------------------------------------
# tablelist::reconfigColLabels
#
# Reconfigures the labels of the col'th column of the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::reconfigColLabels {win imgArrName col} {
    upvar ::tablelist::ns${win}::data data
    upvar $imgArrName imgArr

    foreach opt {-labelalign -labelbackground -labelborderwidth -labelfont
		 -labelforeground -labelheight -labelpady -labelrelief} {
	if {[info exists data($col$opt)]} {
	    doColConfig $col $win $opt $data($col$opt)
	} else {
	    doColConfig $col $win $opt ""
	}
    }

    if {[info exists imgArr($col-labelimage)]} {
	doColConfig $col $win -labelimage $imgArr($col-labelimage)
    }
}

#------------------------------------------------------------------------------
# tablelist::charsToPixels
#
# Returns the width in pixels of the string consisting of a given number of "0"
# characters.
#------------------------------------------------------------------------------
proc tablelist::charsToPixels {win font charCount} {
    ### set str [string repeat 0 $charCount]
    set str ""
    for {set n 0} {$n < $charCount} {incr n} {
	append str 0
    }
    return [font measure $font -displayof $win $str]
}

#------------------------------------------------------------------------------
# tablelist::strRange
#
# Returns the largest initial (for alignment = left or center) or final (for
# alignment = right) range of characters from str whose width, when displayed
# in the given font, is no greater than pixels.
#------------------------------------------------------------------------------
proc tablelist::strRange {win str font pixels alignment} {
    if {[font measure $font -displayof $win $str] <= $pixels} {
	return $str
    }

    set halfLen [expr {[string length $str] / 2}]
    if {$halfLen == 0} {
	return ""
    }

    if {[string compare $alignment right] == 0} {
	set rightStr [string range $str $halfLen end]
	set width [font measure $font -displayof $win $rightStr]
	if {$width == $pixels} {
	    return $rightStr
	} elseif {$width > $pixels} {
	    return [strRange $win $rightStr $font $pixels $alignment]
	} else {
	    set str [string range $str 0 [expr {$halfLen - 1}]]
	    return [strRange $win $str $font \
		    [expr {$pixels - $width}] $alignment]$rightStr
	}
    } else {
	set leftStr [string range $str 0 [expr {$halfLen - 1}]]
	set width [font measure $font -displayof $win $leftStr]
	if {$width == $pixels} {
	    return $leftStr
	} elseif {$width > $pixels} {
	    return [strRange $win $leftStr $font $pixels $alignment]
	} else {
	    set str [string range $str $halfLen end]
	    return $leftStr[strRange $win $str $font \
			    [expr {$pixels - $width}] $alignment]
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::strRangeExt
#
# Invokes strRange with the given arguments and returns a string obtained by
# appending (for alignment = left or center) or prepending (for alignment =
# right) (part of) the snip string to (part of) its result.
#------------------------------------------------------------------------------
proc tablelist::strRangeExt {win str font pixels alignment snipStr} {
    set subStr [strRange $win $str $font $pixels $alignment]
    set len [string length $subStr]
    if {$pixels < 0 || $len == [string length $str] ||
	[string compare $snipStr ""] == 0} {
	return $subStr
    }

    if {[string compare $alignment right] == 0} {
	set extSubStr $snipStr$subStr
	while {[font measure $font -displayof $win $extSubStr] > $pixels} {
	    if {$len > 0} {
		set subStr [string range $subStr 1 end]
		incr len -1
		set extSubStr $snipStr$subStr
	    } else {
		set extSubStr [string range $extSubStr 1 end]
	    }
	}
    } else {
	set last [expr {$len - 1}]
	set extSubStr $subStr$snipStr
	while {[font measure $font -displayof $win $extSubStr] > $pixels} {
	    if {$last >= 0} {
		incr last -1
		set subStr [string range $subStr 0 $last]
		set extSubStr $subStr$snipStr
	    } else {
		set extSubStr [string range $extSubStr 1 end]
	    }
	}
    }

    return $extSubStr
}

#------------------------------------------------------------------------------
# tablelist::adjustItem
#
# Returns the list obtained by adjusting the list specified by item to the
# length expLen.
#------------------------------------------------------------------------------
proc tablelist::adjustItem {item expLen} {
    set len [llength $item]
    if {$len < $expLen} {
	for {set n $len} {$n < $expLen} {incr n} {
	    lappend item ""
	}
	return $item
    } else {
	return [lrange $item 0 [expr {$expLen - 1}]]
    }
}

#------------------------------------------------------------------------------
# tablelist::adjustElem
#
# Prepares the text specified by $textName and the image width specified by
# $imageWidthName for insertion into a cell of the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::adjustElem {win textName imageWidthName font
			    pixels alignment snipStr} {
    upvar $textName text $imageWidthName imageWidth

    if {$pixels == 0} {				;# convention: dynamic width
	if {$imageWidth != 0 && [string compare $text ""] != 0} {
	    if {[string compare $alignment right] == 0} {
		set text "$text "
	    } else {
		set text " $text"
	    }
	}
    } elseif {$imageWidth == 0} {		;# no image
	set text [strRangeExt $win $text $font $pixels $alignment $snipStr]
    } elseif {[string compare $text ""] == 0} {	;# image w/o text
	if {$imageWidth > $pixels} {
	    set imageWidth 0			;# can't display the image
	}
    } else {					;# both image and text
	set gap [font measure $font -displayof $win " "]
	if {$imageWidth + $gap <= $pixels} {
	    incr pixels -[expr {$imageWidth + $gap}]
	    set text [strRangeExt $win $text $font $pixels $alignment $snipStr]
	    if {[string compare $alignment right] == 0} {
		set text "$text "
	    } else {
		set text " $text"
	    }
	} elseif {$imageWidth <= $pixels} {
	    set text ""				;# can't display the text
	} else {
	    set imageWidth 0			;# can't display the image
	    set text ""				;# can't display the text
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::insertElem
#
# Inserts the given text and image into the text widget w, just before the
# character position specified by index.  The image will follow the text if
# alignment is "right", and will precede it otherwise.
#------------------------------------------------------------------------------
proc tablelist::insertElem {w index text image imageWidth alignment} {
    set index [$w index $index]

    if {$imageWidth == 0} {
	$w insert $index $text
    } elseif {[string compare $alignment right] == 0} {
	$w image create $index -image $image
	$w insert $index $text
    } else {
	$w insert $index $text
	$w image create $index -image $image
    }
}

#------------------------------------------------------------------------------
# tablelist::isCharVisible
#
# Checks whether the given text widget character is visible.  It is assumed
# that the line containing the character is visible.
#------------------------------------------------------------------------------
proc tablelist::isCharVisible {w textIdx} {
    set lineEnd [$w index "$textIdx lineend"]
    if {[string compare $textIdx $lineEnd] == 0} {
	return [expr {[lindex [$w xview] 1] == 1.0}]
    } else {
	return [expr {[string compare [$w bbox $textIdx] ""] != 0}]
    }
}

#------------------------------------------------------------------------------
# tablelist::makeColFontAndTagLists
#
# Builds the lists data(colFontList) of the column fonts and data(colTagsList)
# of the column tag names.
#------------------------------------------------------------------------------
proc tablelist::makeColFontAndTagLists win {
    upvar ::tablelist::ns${win}::data data

    set widgetFont $data(-font)
    set data(colFontList) {}
    set data(colTagsList) {}
    set data(hasColTags) 0

    for {set col 0} {$col < $data(colCount)} {incr col} {
	set tagNames {}

	if {[info exists data($col-font)]} {
	    lappend data(colFontList) $data($col-font)
	    lappend tagNames $col-font
	    set data(hasColTags) 1
	} else {
	    lappend data(colFontList) $widgetFont
	}

	foreach opt {-background -foreground} {
	    if {[info exists data($col$opt)]} {
		lappend tagNames $col$opt
		set data(hasColTags) 1
	    }
	}

	lappend data(colTagsList) $tagNames
    }
}

#------------------------------------------------------------------------------
# tablelist::setupColumns
#
# Updates the value of the -colums configuration option for the tablelist
# widget win by using the width, title, and alignment specifications given in
# the columns argument, and creates the corresponding label (and separator)
# widgets if createLabels is true.
#------------------------------------------------------------------------------
proc tablelist::setupColumns {win columns createLabels} {
    variable configSpecs
    variable configOpts
    variable alignments
    upvar ::tablelist::ns${win}::data data

    set argCount [llength $columns]
    set colConfigVals {}

    #
    # Check the syntax of columns before performing any changes
    #
    for {set n 0} {$n < $argCount} {incr n} {
	#
	# Get the column width
	#
	set width [lindex $columns $n]
	set width [format %d $width]	;# integer check with error message

	#
	# Get the column title
	#
	if {[incr n] == $argCount} {
	    return -code error "column title missing"
	}
	set title [lindex $columns $n]

	#
	# Get the column alignment
	#
	set alignment left
	if {[incr n] < $argCount} {
	    set next [lindex $columns $n]
	    if {[catch {format %d $next}] == 0} {	;# integer check
		incr n -1
	    } else {
		set alignment [mwutil::fullOpt "alignment" $next $alignments]
	    }
	}

	#
	# Append the properly formatted values of width,
	# title, and alignment to the list colConfigVals
	#
	lappend colConfigVals $width $title $alignment
    }

    #
    # Save the value of colConfigVals in data(-columns)
    #
    set data(-columns) $colConfigVals

    #
    # Delete the labels and separators if requested
    #
    if {$createLabels} {
	set children [winfo children $data(hdrTxtFr)]
	foreach w [lrange [lsort $children] 1 end] {
	    destroy $w
	}
	foreach w [winfo children $win] {
	    if {[regexp {^sep[0-9]+$} [winfo name $w]]} {
		destroy $w
	    }
	}
	set data(fmtCmdFlagList) {}
    }

    #
    # Build the list data(colList), and create the labels if requested
    #
    set widgetFont $data(-font)
    set data(colList) {}
    set col 0
    foreach {width title alignment} $data(-columns) {
	#
	# Append the width in pixels and the
	# alignment to the list data(colList)
	#
	if {$width > 0} {		;# convention: width in characters
	    set pixels [charsToPixels $win $widgetFont $width]
	    set data($col-lastStaticWidth) $pixels
	} elseif {$width < 0} {		;# convention: width in pixels
	    set pixels [expr {(-1)*$width}]
	    set data($col-lastStaticWidth) $pixels
	} else {			;# convention: dynamic width
	    set pixels 0
	}
	lappend data(colList) $pixels $alignment

	if {$createLabels} {
	    foreach {name val} {delta 0  lastStaticWidth 0  maxPixels 0
				editable 0  editwindow entry  hide 0
				maxwidth 0  resizable 1  showarrow 1
				sortmode ascii} {
		if {![info exists data($col-$name)]} {
		    set data($col-$name) $val
		}
	    }
	    lappend data(fmtCmdFlagList) [info exists data($col-formatcommand)]

	    #
	    # Create the label
	    #
	    set w $data(hdrTxtFrLbl)$col
	    label $w -bitmap "" -highlightthickness 0 -image "" -takefocus 0 \
		     -text "" -textvariable "" -underline -1 -wraplength 0

	    #
	    # Apply to it the current configuration options
	    #
	    foreach opt $configOpts {
		set optGrp [lindex $configSpecs($opt) 2]
		if {[string compare $optGrp l] == 0} {
		    set optTail [string range $opt 6 end]
		    if {[info exists data($col$opt)]} {
			$w configure -$optTail $data($col$opt)
		    } else {
			$w configure -$optTail $data($opt)
		    }
		} elseif {[string compare $optGrp c] == 0} {
		    $w configure $opt $data($opt)
		}
	    }
	    catch {$w configure -state $data(-state)}

	    #
	    # Replace the binding tag Label with TablelistLabel
	    # in the list of binding tags of the label
	    #
	    bindtags $w [lreplace [bindtags $w] 1 1 TablelistLabel]

	    if {[info exists data($col-labelimage)]} {
		doColConfig $col $win -labelimage $data($col-labelimage)
	    }
	}

	#
	# Configure the edit window if present
	#
	if {$col == $data(editCol)} {
	    if {[string compare [winfo class $data(bodyFrEd)] Mentry] != 0} {
		catch {$data(bodyFrEd) configure -justify $alignment}
	    }
	}

	incr col
    }

    #
    # Save the number of columns in data(colCount)
    #
    set oldColCount $data(colCount)
    set data(colCount) $col
    set data(lastCol) [expr {$col - 1}]

    #
    # Clean up the data associated with the deleted columns
    #
    for {set col $data(colCount)} {$col < $oldColCount} {incr col} {
	deleteColData $win $col
    }

    #
    # Create the separators if needed
    #
    if {$createLabels && $data(-showseparators)} {
	createSeps $win
    }
}

#------------------------------------------------------------------------------
# tablelist::createSeps
#
# Creates and manages the separator frames in the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::createSeps win {
    upvar ::tablelist::ns${win}::data data

    for {set col 0} {$col < $data(colCount)} {incr col} {
	#
	# Create the col'th separator frame and attach it
	# to the right edge of the col'th header label
	#
	set w $data(sep)$col
	frame $w -background $data(-background) -borderwidth 1 -container 0 \
		 -highlightthickness 0 -relief raised -takefocus 0 -width 2
	place $w -in $data(hdrTxtFrLbl)$col -anchor ne -bordermode outside \
		 -relx 1.0

	#
	# Replace the binding tag Frame with TablelistBody
	# in the list of binding tags of the separator frame
	#
	bindtags $w [lreplace [bindtags $w] 1 1 TablelistBody]
    }
    
    adjustSepsWhenIdle $win
}

#------------------------------------------------------------------------------
# tablelist::adjustSepsWhenIdle
#
# Arranges for the height and vertical position of each separator frame in the
# tablelist widget win to be adjusted at idle time.
#------------------------------------------------------------------------------
proc tablelist::adjustSepsWhenIdle win {
    upvar ::tablelist::ns${win}::data data

    if {[info exists data(sepsId)]} {
	return ""
    }

    set data(sepsId) [after idle [list tablelist::adjustSeps $win]]
}

#------------------------------------------------------------------------------
# tablelist::adjustSeps
#
# Adjusts the height and vertical position of each separator frame in the
# tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::adjustSeps win {
    upvar ::tablelist::ns${win}::data data

    if {[info exists data(sepsId)]} {
	after cancel $data(sepsId)
	unset data(sepsId)
    }

    #
    # Get the height to be applied to the separator frames
    #
    set w $data(body)
    set textIdx [$w index @0,[winfo height $w]]
    set dlineinfo [$w dlineinfo $textIdx]
    if {$data(itemCount) == 0 || [string compare $dlineinfo ""] == 0} {
	set sepHeight 1
    } else {
	foreach {x y width height baselinePos} $dlineinfo {
	    set sepHeight [expr {$y + $height}]
	}
    }

    #
    # Set the height and vertical position of each separator frame
    #
    foreach w [winfo children $win] {
	if {[regexp {^sep[0-9]+$} [winfo name $w]]} {
	    $w configure -height $sepHeight
	    if {$data(-showlabels)} {
		place configure $w -rely 1.0 -y 0
	    } else {
		place configure $w -rely 0.0 -y 1
	    }
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::adjustColumns
#
# Applies some configuration options to the labels of the tablelist widget win,
# places them in the header frame, computes and sets the tab stops for the body
# text widget, and adjusts the width and height of the header frame.  The
# whichWidths argument specifies the dynamic-width columns or labels whose
# widths are to be computed when performing these operations.  The stretchCols
# argument specifies whether to stretch the stretchable columns.
#------------------------------------------------------------------------------
proc tablelist::adjustColumns {win whichWidths stretchCols} {
    upvar ::tablelist::ns${win}::data data

    set compAllColWidths [expr {[string compare $whichWidths allCols] == 0}]
    set compAllLabelWidths [expr {[string compare $whichWidths allLabels] == 0}]

    #
    # Configure the labels, place them in the header frame, and compute
    # the positions of the tab stops to be set in the body text widget
    #
    set data(hdrPixels) 0
    set tabs {}
    set col 0
    set x 0
    foreach {pixels alignment} $data(colList) {
	set w $data(hdrTxtFrLbl)$col
	if {$data($col-hide)} {
	    place forget $w
	    incr col
	    continue
	}

	#
	# Adjust the col'th label
	#
	if {$pixels != 0} {			;# convention: static width
	    incr pixels $data($col-delta)
	}
	if {[info exists data($col-labelalign)]} {
	    set labelAlignment $data($col-labelalign)
	} else {
	    set labelAlignment $alignment
	}
	adjustLabel $win $col $pixels $labelAlignment

	if {$pixels == 0} {			;# convention: dynamic width
	    #
	    # Compute the column or label width if requested
	    #
	    if {$compAllColWidths} {
		computeColWidth $win $col
	    } elseif {$compAllLabelWidths} {
		computeLabelWidth $win $col
	    } elseif {[lsearch -exact $whichWidths $col] >= 0} {
		computeColWidth $win $col
	    } elseif {[lsearch -exact $whichWidths l$col] >= 0} {
		computeLabelWidth $win $col
	    }

	    set pixels $data($col-reqPixels)
	    if {$data($col-maxPixels) > 0 && $pixels > $data($col-maxPixels)} {
		set pixels $data($col-maxPixels)
		incr pixels $data($col-delta)
		adjustLabel $win $col $pixels $labelAlignment
	    } else {
		incr pixels $data($col-delta)
	    }
	}

	if {$col == $data(editCol)} {
	    adjustEditWindow $win $pixels
	}

	if {$col == $data(arrowCol)} {
	    #
	    # Place the canvas to the left side of the label if the
	    # latter is right-justified and to its right side otherwise
	    #
	    set canvas $data(hdrTxtFrCanv)
	    if {[string compare $labelAlignment right] == 0} {
		place $canvas -in $w -anchor w -bordermode outside \
			      -relx 0.0 -x $data(charWidth) -rely 0.5
	    } else {
		place $canvas -in $w -anchor e -bordermode outside \
			      -relx 1.0 -x -$data(charWidth) -rely 0.5
	    }
	    raise $canvas
	}

	#
	# Place the label in the header frame
	#
	set labelPixels [expr {$pixels + 2*$data(charWidth)}]
	place $w -x $x -relheight 1.0 -width $labelPixels
	incr x $labelPixels

	#
	# Append a tab stop and the alignment to the tabs list
	#
	incr data(hdrPixels) $data(charWidth)
	switch $alignment {
	    left {
		lappend tabs $data(hdrPixels) left
		incr data(hdrPixels) $pixels
	    }
	    right {
		incr data(hdrPixels) $pixels
		lappend tabs $data(hdrPixels) right
	    }
	    center {
		lappend tabs [expr {$data(hdrPixels) + $pixels/2}] center
		incr data(hdrPixels) $pixels
	    }
	}
	incr data(hdrPixels) $data(charWidth)
	lappend tabs $data(hdrPixels) left

	incr col
    }
    place $data(hdrLbl) -x $data(hdrPixels)

    #
    # Apply the value of tabs to the body text widget
    #
    $data(body) configure -tabs $tabs

    #
    # Adjust the width and height of the frames data(hdrTxtFr) and data(hdr)
    #
    $data(hdrTxtFr) configure -width $data(hdrPixels)
    if {$data(-width) <= 0} {
	if {$stretchCols} {
	    $data(hdr) configure -width $data(hdrPixels)
	    $data(lb) configure -width \
		      [expr {$data(hdrPixels) / $data(charWidth)}]
	}
    } else {
	$data(hdr) configure -width 0
    }
    adjustHeaderHeight $win

    #
    # Stretch the stretchable columns if requested
    #
    if {$stretchCols} {
	stretchColumnsWhenIdle $win
    }
}

#------------------------------------------------------------------------------
# tablelist::adjustLabel
#
# Applies some configuration options to the col'th label of the tablelist
# widget win as well as to the label's children (if any), and places the
# children.
#------------------------------------------------------------------------------
proc tablelist::adjustLabel {win col pixels alignment} {
    upvar ::tablelist::ns${win}::data data

    #
    # Apply some configuration options to the label and its children (if any)
    #
    set w $data(hdrTxtFrLbl)$col
    switch $alignment {
	left	{ set anchor w }
	right	{ set anchor e }
	center	{ set anchor center }
    }
    set padX [expr {$data(charWidth) - [$w cget -borderwidth]}]
    $w configure -anchor $anchor -justify $alignment -padx $padX
    if {[info exists data($col-labelimage)]} {
	set imageWidth [image width $data($col-labelimage)]
	if {[string compare $alignment right] == 0} {
	    $w.il configure -anchor e -width 0
	} else {
	    $w.il configure -anchor w -width 0
	}
	$w.tl configure -anchor $anchor -justify $alignment
    } else {
	set imageWidth 0
    }

    #
    # Make room for the canvas displaying an an up- or down-arrow if needed
    #
    set title [lindex $data(-columns) [expr {3*$col + 1}]]
    set labelFont [$w cget -font]
    if {$col == $data(arrowCol)} {
	if {[font metrics $labelFont -displayof $w -fixed]} {
	    set spaces "   "				;# 3 spaces
	} else {
	    set spaces "     "				;# 5 spaces
	}
    } else {
	set spaces ""
    }
    set spacePixels [font measure $labelFont -displayof $w $spaces]

    if {$pixels == 0} {				;# convention: dynamic width
	#
	# Set the label text
	#
	if {$imageWidth == 0} {				;# no image
	    if {[string compare $title ""] == 0} {
		set text $spaces
	    } else {
		set lines {}
		foreach line [split $title \n] {
		    if {[string compare $alignment right] == 0} {
			lappend lines $spaces$line
		    } else {
			lappend lines $line$spaces
		    }
		}
		set text [join $lines \n]
	    }
	    $w configure -text $text
	} elseif {[string compare $title ""] == 0} {	;# image w/o text
	    $w configure -text ""
	    set text ""
	    $w.il configure -width [expr {$imageWidth + $spacePixels}]
	} else {					;# both image and text
	    $w configure -text ""
	    set lines {}
	    foreach line [split $title \n] {
		if {[string compare $alignment right] == 0} {
		    lappend lines $spaces$line
		} else {
		    lappend lines $line$spaces
		}
	    }
	    set text [join $lines \n]
	    $w.tl configure -text $text
	    set colFont [lindex $data(colFontList) $col]
	    set gap [font measure $colFont -displayof $win " "]
	    $w.il configure -width [expr {$imageWidth + $gap}]
	}
    } else {
	#
	# Clip each line of title according to pixels and alignment
	#
	set lessPixels [expr {$pixels - $spacePixels}]
	if {$imageWidth == 0} {				;# no image
	    if {[string compare $title ""] == 0} {
		set text $spaces
	    } else {
		set lines {}
		foreach line [split $title \n] {
		    set line [strRangeExt $win $line $labelFont \
			      $lessPixels $alignment $data(-snipstring)]
		    if {[string compare $alignment right] == 0} {
			lappend lines $spaces$line
		    } else {
			lappend lines $line$spaces
		    }
		}
		set text [join $lines \n]
	    }
	    $w configure -text $text
	} elseif {[string compare $title ""] == 0} {	;# image w/o text
	    $w configure -text ""
	    set text ""
	    if {$imageWidth <= $lessPixels} {
		$w.il configure -width [expr {$imageWidth + $spacePixels}]
	    } else {
		set imageWidth 0		;# can't display the image
	    }
	} else {					;# both image and text
	    $w configure -text ""
	    set colFont [lindex $data(colFontList) $col]
	    set gap [font measure $colFont -displayof $win " "]
	    if {$imageWidth + $gap <= $lessPixels} {
		incr lessPixels -[expr {$imageWidth + $gap}]
		set lines {}
		foreach line [split $title \n] {
		    set line [strRangeExt $win $line $labelFont \
			      $lessPixels $alignment $data(-snipstring)]
		    if {[string compare $alignment right] == 0} {
			lappend lines $spaces$line
		    } else {
			lappend lines $line$spaces
		    }
		}
		set text [join $lines \n]
		$w.tl configure -text $text
		$w.il configure -width [expr {$imageWidth + $gap}]
	    } elseif {$imageWidth <= $lessPixels} {	
		set text ""			;# can't display the text
		$w.il configure -width [expr {$imageWidth + $spacePixels}]
	    } else {
		set imageWidth 0		;# can't display the image
		set text ""			;# can't display the text
	    }
	}
    }

    #
    # Place the label's children (if any)
    #
    if {$imageWidth == 0} {
	if {[info exists data($col-labelimage)]} {
	    place forget $w.il
	    place forget $w.tl
	}
    } else {
	if {[string compare $text ""] == 0} {
	    place forget $w.tl
	}

	switch $alignment {
	    left {
		place $w.il -anchor nw -relx 0.0 -x $padX -relheight 1.0
		if {[string compare $text ""] != 0} {
		    set textX [expr {$padX + [winfo reqwidth $w.il]}]
		    place $w.tl -anchor nw -relx 0.0 -x $textX -relheight 1.0
		}
	    }

	    right {
		place $w.il -anchor ne -relx 1.0 -x -$padX -relheight 1.0
		if {[string compare $text ""] != 0} {
		    set textX [expr {-$padX - [winfo reqwidth $w.il]}]
		    place $w.tl -anchor ne -relx 1.0 -x $textX -relheight 1.0
		}
	    }

	    center {
		if {[string compare $text ""] == 0} {
		    place $w.il -anchor n -relx 0.5 -x 0 -relheight 1.0
		} else {
		    set halfWidth [expr {([winfo reqwidth $w.il] + \
					  [winfo reqwidth $w.tl]) / 2}]
		    place $w.il -anchor nw -relx 0.5 -x -$halfWidth \
				-relheight 1.0
		    place $w.tl -anchor ne -relx 0.5 -x $halfWidth \
				-relheight 1.0
		}
	    }
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::computeColWidth
#
# Computes the width of the col'th column of the tablelist widget win to be just
# large enough to hold all the elements of the column (including its label).
#------------------------------------------------------------------------------
proc tablelist::computeColWidth {win col} {
    upvar ::tablelist::ns${win}::data data

    set fmtCmdFlag [info exists data($col-formatcommand)]
    set colFont [lindex $data(colFontList) $col]

    set data($col-elemWidth) 0
    set data($col-widestCount) 0

    #
    # Column elements
    #
    foreach item $data(itemList) {
	if {$col >= [llength $item] - 1} {
	    continue
	}

	set text [lindex $item $col]
	if {$fmtCmdFlag} {
	    set text [uplevel #0 $data($col-formatcommand) [list $text]]
	}
	set text [strToDispStr $text]
	set key [lindex $item end]
	if {[info exists data($key-$col-image)]} {
	    set imageWidth [image width $data($key-$col-image)]
	} else {
	    set imageWidth 0
	}
	if {[info exists data($key-$col-font)]} {
	    set cellFont $data($key-$col-font)
	} elseif {[info exists data($key-font)]} {
	    set cellFont $data($key-font)
	} else {
	    set cellFont $colFont
	}
	adjustElem $win text imageWidth $cellFont 0 left ""
	set textWidth [font measure $cellFont -displayof $win $text]
	set elemWidth [expr {$imageWidth + $textWidth}]
	if {$elemWidth == $data($col-elemWidth)} {
	    incr data($col-widestCount)
	} elseif {$elemWidth > $data($col-elemWidth)} {
	    set data($col-elemWidth) $elemWidth
	    set data($col-widestCount) 1
	}
    }
    set data($col-reqPixels) $data($col-elemWidth)

    #
    # Column label
    #
    computeLabelWidth $win $col
}

#------------------------------------------------------------------------------
# tablelist::computeLabelWidth
#
# Computes the width of the col'th label of the tablelist widget win and
# adjusts the column's width accordingly.
#------------------------------------------------------------------------------
proc tablelist::computeLabelWidth {win col} {
    upvar ::tablelist::ns${win}::data data

    set w $data(hdrTxtFrLbl)$col
    if {[info exists data($col-labelimage)]} {
	set title [lindex $data(-columns) [expr {3*$col + 1}]]
	if {[string compare $title ""] == 0} {		;# image w/o text
	    set netLabelWidth [winfo reqwidth $w.il]
	} else {					;# both image and text
	    set netLabelWidth [expr {[winfo reqwidth $w.il] +
				     [winfo reqwidth $w.tl]}]
	}
    } else {						;# no image
	set netLabelWidth [expr {[winfo reqwidth $w] - 2*$data(charWidth)}]
    }

    if {$netLabelWidth < $data($col-elemWidth)} {
	set data($col-reqPixels) $data($col-elemWidth)
    } else {
	set data($col-reqPixels) $netLabelWidth
    }
}

#------------------------------------------------------------------------------
# tablelist::adjustHeaderHeight
#
# Sets the height of the header frame of the tablelist widget win to the max.
# height of its children.
#------------------------------------------------------------------------------
proc tablelist::adjustHeaderHeight win {
    upvar ::tablelist::ns${win}::data data

    #
    # Compute the max. label height
    #
    set maxLabelHeight [winfo reqheight $data(hdrLbl)]
    set children [winfo children $data(hdrTxtFr)]
    foreach w [lrange [lsort $children] 1 end] {
	if {[string compare [winfo manager $w] ""] == 0} {
	    continue
	}

	set reqHeight [winfo reqheight $w]
	if {$reqHeight > $maxLabelHeight} {
	    set maxLabelHeight $reqHeight
	}

	foreach c [winfo children $w] {
	    if {[string compare [winfo manager $c] ""] == 0} {
		continue
	    }

	    set reqHeight \
		[expr {[winfo reqheight $c] + 2*[$w cget -borderwidth]}]
	    if {$reqHeight > $maxLabelHeight} {
		set maxLabelHeight $reqHeight
	    }
	}
    }

    #
    # Set the height of the header frame and adjust the separators
    #
    $data(hdrTxtFr) configure -height $maxLabelHeight
    if {$data(-showlabels)} {
	$data(hdr) configure -height $maxLabelHeight
    } else {
	$data(hdr) configure -height 1
    }
    adjustSepsWhenIdle $win
}

#------------------------------------------------------------------------------
# tablelist::stretchColumnsWhenIdle
#
# Arranges for the stretchable columns of the tablelist widget win to be
# stretched at idle time.
#------------------------------------------------------------------------------
proc tablelist::stretchColumnsWhenIdle win {
    upvar ::tablelist::ns${win}::data data

    if {[info exists data(stretchId)]} {
	return ""
    }

    set data(stretchId) [after idle [list tablelist::stretchColumns $win -1]]
}

#------------------------------------------------------------------------------
# tablelist::stretchColumns
#
# Stretches the stretchable columns to fill the tablelist window win
# horizontally.  The colOfFixedDelta argument specifies the column for which
# the stretching is to be made using a precomputed amount of pixels.
#------------------------------------------------------------------------------
proc tablelist::stretchColumns {win colOfFixedDelta} {
    upvar ::tablelist::ns${win}::data data

    if {[info exists data(stretchId)]} {
	after cancel data(stretchId)
	unset data(stretchId)
    }

    set forceAdjust $data(forceAdjust)
    set data(forceAdjust) 0

    if {$data(hdrPixels) == 0 || $data(-width) <= 0} {
	return ""
    }

    #
    # Get the list data(stretchableCols) of the
    # numerical indices of the stretchable columns
    #
    set data(stretchableCols) {}
    if {[string first $data(-stretch) all] == 0} {
	for {set col 0} {$col < $data(colCount)} {incr col} {
	    lappend data(stretchableCols) $col
	}
    } else {
	foreach col $data(-stretch) {
	    lappend data(stretchableCols) [colIndex $win $col 0]
	}
    }

    #
    # Compute the total number data(delta) of pixels by which the
    # columns are to be stretched and the total amount
    # data(stretchablePixels) of stretchable column widths in pixels
    #
    set data(delta) [winfo width $data(hdr)]
    set data(stretchablePixels) 0
    set lastColToStretch -1
    set col 0
    foreach {pixels alignment} $data(colList) {
	if {$data($col-hide)} {
	    incr col
	    continue
	}

	if {$pixels == 0} {			;# convention: dynamic width
	    set pixels $data($col-reqPixels)
	    if {$data($col-maxPixels) > 0 && $pixels > $data($col-maxPixels)} {
		set pixels $data($col-maxPixels)
	    }
	}
	incr data(delta) -[expr {$pixels + 2*$data(charWidth)}]
	if {[lsearch -exact $data(stretchableCols) $col] >= 0} {
	    incr data(stretchablePixels) $pixels
	    set lastColToStretch $col
	}

	incr col
    }
    if {$data(delta) < 0} {
	set delta 0
    } else {
	set delta $data(delta)
    }
    if {$data(stretchablePixels) == 0 && !$forceAdjust} {
	return ""
    }

    #
    # Distribute the value of delta to the stretchable
    # columns, proportionally to their widths in pixels
    #
    set rest $delta
    set col 0
    foreach {pixels alignment} $data(colList) {
	if {$data($col-hide) ||
	    [lsearch -exact $data(stretchableCols) $col] < 0} {
	    set data($col-delta) 0
	} else {
	    set oldDelta $data($col-delta)
	    if {$pixels == 0} {			;# convention: dynamic width
		set pixels $data($col-reqPixels)
		if {$data($col-maxPixels) > 0 &&
		    $pixels > $data($col-maxPixels)} {
		    set pixels $data($col-maxPixels)
		    set dynamic 0
		} else {
		    set dynamic 1
		}
	    } else {
		set dynamic 0
	    }
	    if {$data(stretchablePixels) == 0} {
		set data($col-delta) 0
	    } else {
		if {$col != $colOfFixedDelta} {
		    set data($col-delta) \
			[expr {$delta*$pixels/$data(stretchablePixels)}]
		}
		incr rest -$data($col-delta)
	    }
	    if {$col == $lastColToStretch} {
		incr data($col-delta) $rest
	    }
	    if {!$dynamic && $data($col-delta) != $oldDelta} {
		redisplayWhenIdle $win
	    }
	}

	incr col
    }

    #
    # Adjust the columns
    #
    adjustColumns $win {} 0
    update idletasks
}

#------------------------------------------------------------------------------
# tablelist::redisplayWhenIdle
#
# Arranges for the items of the tablelist widget win to be redisplayed at idle
# time.
#------------------------------------------------------------------------------
proc tablelist::redisplayWhenIdle win {
    upvar ::tablelist::ns${win}::data data

    if {[info exists data(redispId)] || $data(itemCount) == 0} {
	return ""
    }

    set data(redispId) [after idle [list tablelist::redisplay $win]]
}

#------------------------------------------------------------------------------
# tablelist::redisplay
#
# Redisplays the items of the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::redisplay win {
    upvar ::tablelist::ns${win}::data data

    if {[info exists data(redispId)]} {
	after cancel $data(redispId)
	unset data(redispId)
    }

    #
    # Save some data of the edit window if present
    #
    if {[set editCol $data(editCol)] >= 0} {
	set editRow $data(editRow)
	saveEditData $win
    }

    set w $data(body)
    set widgetFont $data(-font)
    set snipStr $data(-snipstring)
    set isSimple [expr {$data(tagCount) == 0 && $data(imgCount) == 0 &&
			!$data(hasColTags)}]
    set newItemList {}
    set idx 0
    set line 1
    foreach item $data(itemList) {
	#
	# Check whether the line is selected
	#
	set tagNames [$w tag names $line.0]
	if {[lsearch -exact $tagNames select] >= 0} {
	    set selected 1
	} else {
	    set selected 0
	}

	#
	# Empty the line, clip the elements if necessary,
	# and insert them with the corresponding tags
	#
	$w delete $line.0 $line.end
	set keyIdx [expr {[llength $item] - 1}]
	set key [lindex $item end]
	set newItem {}
	set col 0
	if {$isSimple} {
	    set insertStr ""
	    foreach fmtCmdFlag $data(fmtCmdFlagList) \
		    {pixels alignment} $data(colList) {
		if {$col < $keyIdx} {
		    set text [lindex $item $col]
		} else {
		    set text ""
		}
		lappend newItem $text

		if {$data($col-hide)} {
		    incr col
		    continue
		}

		#
		# Clip the element if necessary
		#
		if {$fmtCmdFlag} {
		    set text [uplevel #0 $data($col-formatcommand) [list $text]]
		}
		set text [strToDispStr $text]
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

	}  else {
	    array set tagData [array get data $key*-\[bf\]*]	;# for speed

	    set rowTags [array names tagData $key-\[bf\]*]
	    foreach colFont $data(colFontList) \
		    colTags $data(colTagsList) \
		    fmtCmdFlag $data(fmtCmdFlagList) \
		    {pixels alignment} $data(colList) {
		if {$col < $keyIdx} {
		    set text [lindex $item $col]
		} else {
		    set text ""
		}
		lappend newItem $text

		if {$data($col-hide)} {
		    incr col
		    continue
		}

		#
		# Adjust the cell text and the image width
		#
		if {$fmtCmdFlag} {
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
	lappend newItem $key
	lappend newItemList $newItem

	#
	# Select the item if it was selected before
	#
	if {$selected} {
	    selectionSubCmd $win set $idx $idx
	}

	incr idx
	incr line
    }

    set data(itemList) $newItemList

    #
    # Restore the stripes in the body text widget
    #
    makeStripes $win

    #
    # Restore the edit window if it was present before
    #
    if {$editCol >= 0} {
	editcellSubCmd $win $editRow $editCol 1
    }
}

#------------------------------------------------------------------------------
# tablelist::makeStripesWhenIdle
#
# Arranges for the stripes in the body of the tablelist widget win to be
# redrawn at idle time.
#------------------------------------------------------------------------------
proc tablelist::makeStripesWhenIdle win {
    upvar ::tablelist::ns${win}::data data

    if {[info exists data(stripesId)] || $data(itemCount) == 0} {
	return ""
    }

    set data(stripesId) [after idle [list tablelist::makeStripes $win]]
}

#------------------------------------------------------------------------------
# tablelist::makeStripes
#
# Redraws the stripes in the body of the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::makeStripes win {
    upvar ::tablelist::ns${win}::data data

    if {[info exists data(stripesId)]} {
	after cancel $data(stripesId)
	unset data(stripesId)
    }

    set w $data(body)
    $w tag remove stripe 1.0 end
    if {[string compare $data(-stripebackground) ""] == 0 &&
	[string compare $data(-stripeforeground) ""] == 0} {
	return ""
    }

    set step [expr {2*$data(-stripeheight)}]
    for {set n [expr {$data(-stripeheight) + 1}]} {$n <= $step} {incr n} {
	for {set line $n} {$line <= $data(itemCount)} {incr line $step} {
	    $w tag add stripe $line.0 $line.end
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::synchronize
#
# This procedure is invoked either as an idle callback after the list variable
# associated with the tablelist widget win was written, or directly, upon
# execution of some widget commands.  It makes sure that the content of the
# widget is synchronized with the value of the list variable.
#------------------------------------------------------------------------------
proc tablelist::synchronize win {
    upvar ::tablelist::ns${win}::data data

    #
    # Nothing to do if the list variable was not written
    #
    if {![info exists data(syncId)]} {
	return ""
    }

    #
    # Here we are in the case that the procedure was scheduled for
    # execution at idle time.  However, it might have been invoked
    # directly, before the idle time occured; in this case we should
    # cancel the execution of the previously scheduled idle callback.
    #
    after cancel $data(syncId)	;# no harm if data(syncId) is no longer valid
    unset data(syncId)

    upvar #0 $data(-listvariable) var
    set newCount [llength $var]
    if {$newCount < $data(itemCount)} {
	#
	# Delete the items with indices >= newCount from the widget
	#
	set updateCount $newCount
	deleteRows $win $newCount $data(lastRow) 0
    } elseif {$newCount > $data(itemCount)} {
	#
	# Insert the items of var with indices
	# >= data(itemCount) into the widget
	#
	set updateCount $data(itemCount)
	insertSubCmd $win $data(itemCount) \
		     [lrange $var $data(itemCount) end] 0
    } else {
	set updateCount $newCount
    }

    #
    # Update the first updateCount items of the internal list
    #
    set itemsChanged 0
    for {set idx 0} {$idx < $updateCount} {incr idx} {
	set oldItem [lindex $data(itemList) $idx]
	set newItem [adjustItem [lindex $var $idx] $data(colCount)]
	lappend newItem [lindex $oldItem end]

	if {[string compare $oldItem $newItem] != 0} {
	    set data(itemList) [lreplace $data(itemList) $idx $idx $newItem]
	    set itemsChanged 1
	}
    }

    #
    # If necessary, adjust the columns and make sure
    # that the items will be redisplayed at idle time
    #
    if {$itemsChanged} {
	adjustColumns $win allCols 1
	redisplayWhenIdle $win
    }
}

#------------------------------------------------------------------------------
# tablelist::configLabel
#
# This procedure configures the label widget w according to the options and
# their values given in args.  It is needed for label widgets with children,
# managed by the place geometry manager, because - strangely enough - by just
# configuring the label causes its children to become invisible on Windows (but
# not on UNIX).  The procedure solves this problem by using a trick: after
# configuring the label, it applies a constant configuration value to its
# children, which makes them visible again.
#------------------------------------------------------------------------------
proc tablelist::configLabel {w args} {
    eval [list $w configure] $args

    foreach c [winfo children $w] {
	$c configure -borderwidth 0
    }
}

#------------------------------------------------------------------------------
# tablelist::create3DArrows
#
# Creates the items to be used later for drawing two up- or down-arrows with
# sunken relief and 3-D borders in the canvas w.
#------------------------------------------------------------------------------
proc tablelist::create3DArrows w {
    foreach state {normal disabled} {
	$w create polygon 0 0 0 0 0 0 -tags ${state}Triangle
	$w create line    0 0 0 0     -tags ${state}DarkLine
	$w create line    0 0 0 0     -tags ${state}LightLine
    }
}

#------------------------------------------------------------------------------
# tablelist::configCanvas
#
# Sets the background, width, and height of the canvas displaying an up- or
# down-arrow, fills the two arrows contained in the canvas, and saves its width
# in data(arrowWidth).
#------------------------------------------------------------------------------
proc tablelist::configCanvas win {
    upvar ::tablelist::ns${win}::data data

    set w $data(hdrTxtFrLbl)$data(arrowCol)
    set labelBg [$w cget -background]
    set labelFont [$w cget -font]
    if {[font metrics $labelFont -displayof $w -fixed]} {
	set spaces " "
    } else {
	set spaces "  "
    }

    set size [expr {[font measure $labelFont -displayof $w $spaces] + 2}]
    if {$size % 2 == 0} {
	incr size
    }

    set w $data(hdrTxtFrCanv)
    $w configure -background $labelBg -height $size -width $size
    fillArrow $w normal   $data(-arrowcolor)
    fillArrow $w disabled $data(-arrowdisabledcolor)

    set data(arrowWidth) $size
}

#------------------------------------------------------------------------------
# tablelist::drawArrows
#
# Draws the two arrows contained in the canvas associated with the tablelist
# widget win.
#------------------------------------------------------------------------------
proc tablelist::drawArrows win {
    upvar ::tablelist::ns${win}::data data

    switch $data(-incrarrowtype) {
	up {
	    switch $data(sortOrder) {
		increasing { set arrowType up }
		decreasing { set arrowType down }
	    }
	}

	down {
	    switch $data(sortOrder) {
		increasing { set arrowType down }
		decreasing { set arrowType up }
	    }
	}
    }

    set w $data(hdrTxtFrCanv)
    set maxX [expr {[$w cget -width] - 1}]
    set maxY [expr {[$w cget -height] - 1}]
    set midX [expr {$maxX / 2}]

    switch $arrowType {
	up {
	    foreach state {normal disabled} {
		$w coords ${state}Triangle  0 $maxY $maxX $maxY $midX 0
		$w coords ${state}DarkLine  $midX 0 0 $maxY
		$w coords ${state}LightLine 0 $maxY $maxX $maxY $midX 0
	    }
	}

	down {
	    foreach state {normal disabled} {
		$w coords ${state}Triangle  $maxX 0 0 0 $midX $maxY
		$w coords ${state}DarkLine  $maxX 0 0 0 $midX $maxY
		$w coords ${state}LightLine $midX $maxY $maxX 0
	    }
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::fillArrow
#
# Fills one of the two arrows contained in the canvas w with the given color,
# or with (a slightly darker color than) the background color of the canvas if
# color is an empty string.  Also fills the arrow's borders with the
# corresponding 3-D shadow colors.  The state argument specifies the arrow to
# be processed.  Returns the properly formatted value of color.
#------------------------------------------------------------------------------
proc tablelist::fillArrow {w state color} {
    if {[string compare $color ""] == 0} {
	set origColor $color
	set color [$w cget -background]

	#
	# To get a better contrast, make the color slightly
	# darker by cutting 5% from each of its components
	#
	set maxIntens [lindex [winfo rgb $w white] 0]
	set len [string length [format %x $maxIntens]]
	foreach comp [winfo rgb $w $color] {
	    lappend rgb [expr {95*$comp/100}]
	}
	set color [eval format "#%0${len}x%0${len}x%0${len}x" $rgb]
    }

    getShadows $w $color darkColor lightColor

    $w itemconfigure ${state}Triangle  -fill $color
    $w itemconfigure ${state}DarkLine  -fill $darkColor
    $w itemconfigure ${state}LightLine -fill $lightColor

    if {[info exists origColor]} {
	return $origColor
    } else {
	return [$w itemcget ${state}Triangle -fill]
    }
}

#------------------------------------------------------------------------------
# tablelist::getShadows
#
# Computes the shadow colors for a 3-D border from a given (background) color.
# This is a modified Tcl-counterpart of the function TkpGetShadows() in the
# Tk distribution file unix/tkUnix3d.c.
#------------------------------------------------------------------------------
proc tablelist::getShadows {w color darkColorName lightColorName} {
    upvar $darkColorName darkColor $lightColorName lightColor

    set maxIntens [lindex [winfo rgb $w white] 0]
    set len [string length [format %x $maxIntens]]

    set rgb [winfo rgb $w $color]
    foreach {r g b} $rgb {}

    #
    # Compute the dark shadow color
    #
    if {[string compare $::tk_patchLevel 8.3.1] >= 0 &&
	$r*0.5*$r + $g*1.0*$g + $b*0.28*$b < $maxIntens*0.05*$maxIntens} {
	#
	# The background is already very dark: make the dark
	# color a little lighter than the background by increasing
	# each color component 1/4th of the way to $maxIntens
	#
	foreach comp $rgb {
	    lappend darkRGB [expr {($maxIntens + 3*$comp)/4}]
	}
    } else {
	#
	# Compute the dark color by cutting 45% from
	# each of the background color components.
	#
	foreach comp $rgb {
	    lappend darkRGB [expr {55*$comp/100}]
	}
    }
    set darkColor [eval format "#%0${len}x%0${len}x%0${len}x" $darkRGB]

    #
    # Compute the light shadow color
    #
    if {[string compare $::tk_patchLevel 8.3.1] >= 0 && $g > $maxIntens*0.95} {
	#
	# The background is already very bright: make the
	# light color a little darker than the background
	# by reducing each color component by 10%
	#
	foreach comp $rgb {
	    lappend lightRGB [expr {9*$comp/10}]
	}
    } else {
	#
	# Compute the light color by boosting each background
	# color component by 45% or half-way to white, whichever
	# is greater (the first approach works better for
	# unsaturated colors, the second for saturated ones)
	#
	foreach comp $rgb {
	    set comp1 [expr {145*$comp/100}]
	    if {$comp1 > $maxIntens} {
		set comp1 $maxIntens
	    }
	    set comp2 [expr {($maxIntens + $comp)/2}]
	    lappend lightRGB [expr {($comp1 > $comp2) ? $comp1 : $comp2}]
	}
    }
    set lightColor [eval format "#%0${len}x%0${len}x%0${len}x" $lightRGB]
}

#------------------------------------------------------------------------------
# tablelist::raiseArrow
#
# Raises one of the two arrows contained in the canvas w, according to the
# state argument.
#------------------------------------------------------------------------------
proc tablelist::raiseArrow {w state} {
    $w raise ${state}Triangle
    $w raise ${state}DarkLine
    $w raise ${state}LightLine
}

#------------------------------------------------------------------------------
# tablelist::isCellEditable
#
# Checks whether the given cell of the tablelist widget win is editable.
#------------------------------------------------------------------------------
proc tablelist::isCellEditable {win row col} {
    upvar ::tablelist::ns${win}::data data

    set item [lindex $data(itemList) $row]
    set key [lindex $item end]
    if {[info exists data($key-$col-editable)]} {
	return $data($key-$col-editable)
    } else {
	return $data($col-editable)
    }
}
