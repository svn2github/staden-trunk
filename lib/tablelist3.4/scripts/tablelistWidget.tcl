#==============================================================================
# Contains the implementation of the tablelist widget.
#
# Structure of the module:
#   - Namespace initialization
#   - Public procedure
#   - Private configuration procedures
#   - Private procedures implementing the tablelist widget command
#   - Private callback procedures
#   - Private procedures used in bindings
#   - Private utility procedures
#
# Copyright (c) 2000-2004  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

#
# Namespace initialization
# ========================
#

namespace eval tablelist {
    #
    # The array configSpecs is used to handle configuration options.  The
    # names of its elements are the configuration options for the Tablelist
    # class.  The value of an array element is either an alias name or a list
    # containing the database name and class as well as an indicator specifying
    # the widget(s) to which the option applies: c stands for all children
    # (text widgets and labels), b for the body text widget, h for the header
    # text widget, l for the labels, f for the frame, and w for the widget
    # itself.
    #
    #	Command-Line Name	 {Database Name		  Database Class      W}
    #	------------------------------------------------------------------------
    #
    variable configSpecs
    array set configSpecs {
	-activestyle		 {activeStyle		  ActiveStyle	      w}
	-arrowcolor		 {arrowColor		  ArrowColor	      w}
	-arrowdisabledcolor	 {arrowDisabledColor	  ArrowDisabledColor  w}
	-background		 {background		  Background	      b}
	-bg			 -background
	-borderwidth		 {borderWidth		  BorderWidth	      f}
	-bd			 -borderwidth
	-columns		 {columns		  Columns	      w}
	-cursor			 {cursor		  Cursor	      c}
	-disabledforeground	 {disabledForeground	  DisabledForeground  w}
	-editendcommand		 {editEndCommand	  EditEndCommand      w}
	-editstartcommand	 {editStartCommand	  EditStartCommand    w}
	-exportselection	 {exportSelection	  ExportSelection     w}
	-font			 {font			  Font		      b}
	-foreground		 {foreground		  Foreground	      b}
	-fg			 -foreground
	-height			 {height		  Height	      w}
	-highlightbackground	 {highlightBackground	  HighlightBackground f}
	-highlightcolor		 {highlightColor	  HighlightColor      f}
	-highlightthickness	 {highlightThickness	  HighlightThickness  f}
	-incrarrowtype		 {incrArrowType		  IncrArrowType	      w}
	-labelbackground	 {labelBackground	  Background	      l}
	-labelbg		 -labelbackground
	-labelborderwidth	 {labelBorderWidth	  BorderWidth	      l}
	-labelbd		 -labelborderwidth
	-labelcommand		 {labelCommand		  LabelCommand	      w}
	-labeldisabledforeground {labelDisabledForeground DisabledForeground  l}
	-labelfont		 {labelFont		  Font		      l}
	-labelforeground	 {labelForeground	  Foreground	      l}
	-labelfg		 -labelforeground
	-labelheight		 {labelHeight		  Height	      l}
	-labelpady		 {labelPadY		  Pad		      l}
	-labelrelief		 {labelRelief		  Relief	      l}
	-listvariable		 {listVariable		  Variable	      w}
	-movablecolumns	 	 {movableColumns	  MovableColumns      w}
	-movecolumncursor	 {moveColumnCursor	  MoveColumnCursor    w}
	-relief			 {relief		  Relief	      f}
	-resizablecolumns	 {resizableColumns	  ResizableColumns    w}
	-resizecursor		 {resizeCursor		  ResizeCursor	      w}
	-selectbackground	 {selectBackground	  Foreground	      w}
	-selectborderwidth	 {selectBorderWidth	  BorderWidth	      w}
	-selectforeground	 {selectForeground	  Background	      w}
	-selectmode		 {selectMode		  SelectMode	      w}
	-setgrid		 {setGrid		  SetGrid	      w}
	-showarrow		 {showArrow		  ShowArrow	      w}
	-showlabels		 {showLabels		  ShowLabels	      w}
	-showseparators		 {showSeparators	  ShowSeparators      w}
	-snipstring		 {snipString		  SnipString	      w}
	-sortcommand		 {sortCommand		  SortCommand	      w}
	-state			 {state			  State		      w}
	-stretch		 {stretch		  Stretch	      w}
	-stripebackground	 {stripeBackground	  Background	      w}
	-stripebg		 -stripebackground
	-stripeforeground	 {stripeForeground	  Foreground	      w}
	-stripefg		 -stripeforeground
	-stripeheight		 {stripeHeight		  StripeHeight	      w}
	-takefocus		 {takeFocus		  TakeFocus	      f}
	-targetcolor		 {targetColor		  TargetColor	      w}
	-width			 {width			  Width		      w}
	-xscrollcommand		 {xScrollCommand	  ScrollCommand	      h}
	-yscrollcommand		 {yScrollCommand	  ScrollCommand	      b}
    }

    #
    # Get the current windowing system ("x11", "win32", "classic", or "aqua")
    #
    variable winSys
    if {[catch {tk windowingsystem} winSys] != 0} {
	switch $::tcl_platform(platform) {
	    unix	{ set winSys x11 }
	    windows	{ set winSys win32 }
	    macintosh	{ set winSys classic }
	}
    }

    #
    # Extend the elements of the array configSpecs
    #
    extendConfigSpecs 

    variable configOpts [lsort [array names configSpecs]]

    #
    # The array colConfigSpecs is used to handle column configuration options.
    # The names of its elements are the column configuration options for the
    # Tablelist widget class.  The value of an array element is either an alias
    # name or a list containing the database name and class.
    #
    #	Command-Line Name	{Database Name		Database Class	}
    #	-----------------------------------------------------------------
    #
    variable colConfigSpecs
    array set colConfigSpecs {
	-align			{align			Align		}
	-background		{background		Background	}
	-bg			-background
	-editable		{editable		Editable	}
	-editwindow		{editWindow		EditWindow	}
	-font			{font			Font		}
	-foreground		{foreground		Foreground	}
	-fg			-foreground
	-formatcommand		{formatCommand		FormatCommand	}
	-hide			{hide			Hide		}
	-labelalign		{labelAlign		Align		}
	-labelbackground	{labelBackground	Background	}
	-labelbg		-labelbackground
	-labelborderwidth	{labelBorderWidth	BorderWidth	}
	-labelbd		-labelborderwidth
	-labelcommand		{labelCommand		LabelCommand	}
	-labelfont		{labelFont		Font		}
	-labelforeground	{labelForeground	Foreground	}
	-labelfg		-labelforeground
	-labelheight		{labelHeight		Height		}
	-labelimage		{labelImage		Image		}
	-labelpady		{labelPadY		Pad		}
	-labelrelief		{labelRelief		Relief		}
	-maxwidth		{maxWidth		MaxWidth	}
	-name			{name			Name		}
	-resizable		{resizable		Resizable	}
	-selectbackground	{selectBackground	Foreground	}
	-selectforeground	{selectForeground	Background	}
	-showarrow		{showArrow		ShowArrow	}
	-sortcommand		{sortCommand		SortCommand	}
	-sortmode		{sortMode		SortMode	}
	-text			{text			Text		}
	-title			{title			Title		}
	-width			{width			Width		}
    }

    #
    # Extend some elements of the array colConfigSpecs
    #
    lappend colConfigSpecs(-align)	- left
    lappend colConfigSpecs(-editable)	- 0
    lappend colConfigSpecs(-editwindow)	- entry
    lappend colConfigSpecs(-hide)	- 0
    lappend colConfigSpecs(-maxwidth)	- 0
    lappend colConfigSpecs(-resizable)	- 1
    lappend colConfigSpecs(-showarrow)	- 1
    lappend colConfigSpecs(-sortmode)	- ascii
    lappend colConfigSpecs(-width)	- 0

    #
    # The array rowConfigSpecs is used to handle row configuration options.
    # The names of its elements are the row configuration options for the
    # Tablelist widget class.  The value of an array element is either an alias
    # name or a list containing the database name and class.
    #
    #	Command-Line Name	{Database Name		Database Class	}
    #	-----------------------------------------------------------------
    #
    variable rowConfigSpecs
    array set rowConfigSpecs {
	-background		{background		Background	}
	-bg			-background
	-font			{font			Font		}
	-foreground		{foreground		Foreground	}
	-fg			-foreground
	-selectable		{selectable		Selectable	}
	-selectbackground	{selectBackground	Foreground	}
	-selectforeground	{selectForeground	Background	}
	-text			{text			Text		}
    }

    #
    # Extend some elements of the array rowConfigSpecs
    #
    lappend rowConfigSpecs(-selectable) - 1

    #
    # The array cellConfigSpecs is used to handle cell configuration options.
    # The names of its elements are the cell configuration options for the
    # Tablelist widget class.  The value of an array element is either an alias
    # name or a list containing the database name and class.
    #
    #	Command-Line Name	{Database Name		Database Class	}
    #	-----------------------------------------------------------------
    #
    variable cellConfigSpecs
    array set cellConfigSpecs {
	-background		{background		Background	}
	-bg			-background
	-editable		{editable		Editable	}
	-font			{font			Font		}
	-foreground		{foreground		Foreground	}
	-fg			-foreground
	-image			{image			Image		}
	-selectbackground	{selectBackground	Foreground	}
	-selectforeground	{selectForeground	Background	}
	-text			{text			Text		}
    }

    #
    # Extend some elements of the array cellConfigSpecs
    #
    lappend cellConfigSpecs(-editable) - 1

    #
    # Use a list to facilitate the handling of the command options 
    #
    variable cmdOpts [list \
	activate attrib bbox bodypath cancelediting cellcget cellconfigure \
	cellindex cget columncget columnconfigure columncount columnindex \
	configure containing containingcell containingcolumn curselection \
	delete deletecolumns editcell editwinpath entrypath fillcolumn \
	finishediting get getcolumns getkeys index insert insertcolumnlist \
	insertcolumns insertlist labelpath labels move movecolumn nearest \
	nearestcell nearestcolumn rejectinput resetsortinfo rowcget \
	rowconfigure scan see seecell seecolumn selection separatorpath \
	separators size sort sortbycolumn sortcolumn sortorder xview yview]

    #
    # Use lists to facilitate the handling of miscellaneous options
    #
    variable activeStyles	[list frame none underline]
    variable alignments		[list left right center]
    variable arrowTypes		[list up down]
    variable states		[list disabled normal]
    variable sortModes		[list ascii command dictionary integer real]
    variable sortOrders		[list -increasing -decreasing]
    variable scanCmdOpts	[list mark dragto]
    variable selCmdOpts		[list anchor clear includes set]

    #
    # Define the procedure strToDispStr, which returns the string obtained
    # by replacing all \t and \n characters in its argument with \\t and
    # \\n, respectively, as well as the procedure strMap, needed because
    # the "string map" command is not available in Tcl 8.0 and 8.1.0.
    #
    if {[catch {string map {} ""}] == 0} {
	proc strToDispStr str {
	    if {[string first \t $str] >= 0 || [string first \n $str] >= 0} {
		return [string map {\t \\t  \n \\n} $str]
	    } else {
		return $str
	    }
	}

	proc strMap {charMap str} {
	    return [string map $charMap $str]
	}
    } else {
	proc strToDispStr str {
	    if {[string first \t $str] >= 0} {
		regsub -all \t $str \\t str
	    }
	    if {[string first \n $str] >= 0} {
		regsub -all \n $str \\n str
	    }

	    return $str
	}

	proc strMap {charMap str} {
	    foreach {key val} $charMap {
		#
		# We will only need this for noncritical key values
		#
		regsub -all $key $str $val str
	    }

	    return $str
	}
    }

    #
    # Define some Tablelist class bindings
    #
    bind Tablelist <KeyPress> continue
    bind Tablelist <FocusIn> {
	tablelist::addActiveTag %W

	if {[string compare [focus -lastfor %W] %W] == 0} {
	    if {[winfo exists [%W editwinpath]]} {
		focus [set tablelist::ns%W::data(editFocus)]
	    } else {
		focus [%W bodypath]
	    }
	}
    }
    bind Tablelist <FocusOut> {
	tablelist::removeActiveTag %W
    }
    bind Tablelist <Destroy> {
	tablelist::cleanup %W
    }

    #
    # Define the binding tags TablelistKeyNav and TablelistBody
    #
    mwutil::defineKeyNav Tablelist
    defineTablelistBody 

    #
    # Define the virtual events <<Button3>> and <<ShiftButton3>>
    #
    event add <<Button3>> <Button-3>
    event add <<ShiftButton3>> <Shift-Button-3>
    if {[string compare $winSys classic] == 0 ||
	[string compare $winSys aqua] == 0} {
	event add <<Button3>> <Control-Button-1>
	event add <<ShiftButton3>> <Shift-Control-Button-1>
    }

    #
    # Define some mouse bindings for the binding tag TablelistLabel
    #
    bind TablelistLabel <Enter>		  { tablelist::labelEnter    %W %x }
    bind TablelistLabel <Motion>	  { tablelist::labelEnter    %W %x }
    bind TablelistLabel <Button-1>	  { tablelist::labelB1Down   %W %x }
    bind TablelistLabel <B1-Motion>	  { tablelist::labelB1Motion %W %x %y }
    bind TablelistLabel <B1-Enter>	  { tablelist::labelB1Enter  %W }
    bind TablelistLabel <B1-Leave>	  { tablelist::labelB1Leave  %W %x %y }
    bind TablelistLabel <ButtonRelease-1> { tablelist::labelB1Up     %W %X}
    bind TablelistLabel <<Button3>>	  { tablelist::labelB3Down   %W }
    bind TablelistLabel <<ShiftButton3>>  { tablelist::labelShiftB3Down %W }

    #
    # Define the binding tags TablelistSubLabel and TablelistArrow
    #
    defineTablelistSubLabel 
    defineTablelistArrow 
}

#
# Public procedure
# ================
#

#------------------------------------------------------------------------------
# tablelist::tablelist
#
# Creates a new tablelist widget whose name is specified as the first command-
# line argument, and configures it according to the options and their values
# given on the command line.  Returns the name of the newly created widget.
#------------------------------------------------------------------------------
proc tablelist::tablelist args {
    variable configSpecs
    variable configOpts

    if {[llength $args] == 0} {
	mwutil::wrongNumArgs "tablelist pathName ?options?"
    }

    #
    # Create a frame of the class Tablelist
    #
    set win [lindex $args 0]
    if {[catch {
	frame $win -class Tablelist -container 0 -height 0 -width 0
    } result] != 0} {
	return -code error $result
    }

    #
    # Create a namespace within the current one to hold the data of the widget
    #
    namespace eval ns$win {
	#
	# The folowing array holds various data for this widget
	#
	variable data
	array set data {
	    hasListVar		 0
	    isDisabled		 0
	    ownsFocus		 0
	    charWidth		 1
	    hdrPixels		 0
	    oldActiveIdx	 0
	    activeIdx		 0
	    anchorIdx		 0
	    seqNum		-1
	    itemList		 {}
	    itemCount		 0
	    lastRow		-1
	    colList		 {}
	    colCount		 0
	    lastCol		-1
	    tagCount		 0
	    imgCount		 0
	    labelClicked	 0
	    arrowCol		-1
	    sortCol		-1
	    sortOrder		 {}
	    editRow		-1
	    editCol		-1
	    forceAdjust		 0
	    fmtCmdFlagList	 {}
	}

	#
	# The following array is used to hold arbitrary
	# attributes and their values for this widget
	#
	variable attribVals
    }

    #
    # Initialize some further components of data
    #
    upvar ::tablelist::ns${win}::data data
    foreach opt $configOpts {
	set data($opt) [lindex $configSpecs($opt) 3]
    }
    set data(colFontList)	[list $data(-font)]
    set data(listVarTraceCmd)	[list tablelist::listVarTrace $win]
    set data(body)		$win.body
    set data(bodyFr)		$data(body).f
    set data(bodyFrEd)		$data(bodyFr).e
    set data(hdr)		$win.hdr
    set data(hdrTxt)		$data(hdr).t
    set data(hdrTxtFr)		$data(hdrTxt).f
    set data(hdrTxtFrCanv)	$data(hdrTxtFr).c
    set data(hdrTxtFrLbl)	$data(hdrTxtFr).l
    set data(hdrLbl)		$data(hdr).l
    set data(hdrGap)		$data(hdr).g
    set data(lb)		$win.lb
    set data(sep)		$win.sep

    #
    # Create a child hierarchy used to hold the column labels.  The
    # labels will be created as children of the frame data(hdrTxtFr),
    # which is embedded into the text widget data(hdrTxt) (in order
    # to make it scrollable), which in turn fills the frame data(hdr)
    # (whose width and height can be set arbitrarily in pixels).
    #
    set w $data(hdr)			;# header frame
    frame $w -borderwidth 0 -container 0 -height 0 -highlightthickness 0 \
	     -relief flat -takefocus 0 -width 0
    bind $w <Configure> { tablelist::stretchColumnsWhenIdle [winfo parent %W] }
    pack $w -fill x
    set w $data(hdrTxt)			;# text widget within the header frame
    text $w -borderwidth 0 -highlightthickness 0 -insertwidth 0 \
	    -padx 0 -pady 0 -state normal -takefocus 0 -wrap none
    place $w -relheight 1.0 -relwidth 1.0
    bindtags $w [lreplace [bindtags $w] 1 1]
    frame $data(hdrTxtFr) -borderwidth 0 -container 0 -height 0 \
			  -highlightthickness 0 -relief flat -takefocus 0 \
			  -width 0
    $w window create 1.0 -window $data(hdrTxtFr)
    label $data(hdrTxtFrLbl)0 
    set w $data(hdrLbl)			;# filler label within the header frame
    label $w -bitmap "" -highlightthickness 0 -image "" -takefocus 0 \
	     -text "" -textvariable "" -underline -1 -wraplength 0
    place $w -relheight 1.0 -relwidth 1.0

    #
    # Create a canvas as a child of the frame data(hdrTxtFr),
    # needed for displaying an up- or down-arrow when
    # sorting the items by a column.   Set its width and
    # height to temporary values and create two 3-D arrows
    #
    set w $data(hdrTxtFrCanv)
    set size 9
    canvas $w -borderwidth 0 -height $size -highlightthickness 0 \
	      -relief flat -takefocus 0 -width $size
    create3DArrows $w

    #
    # Replace the binding tag Canvas with TablelistArrow
    # in the list of binding tags of the canvas
    #
    bindtags $w [lreplace [bindtags $w] 1 1 TablelistArrow]

    #
    # Create a frame used to display a gap between two
    # consecutive columns when moving a column interactively
    #
    frame $data(hdrGap) -borderwidth 1 -container 0 -highlightthickness 0 \
			-relief sunken -takefocus 0 -width 4

    #
    # Create the body text widget within the main frame
    #
    set w $data(body)
    text $w -borderwidth 0 -exportselection 0 -highlightthickness 0 \
	    -insertwidth 0 -padx 0 -pady 0 -state normal -takefocus 0 -wrap none
    bind $w <Configure> { tablelist::adjustSepsWhenIdle [winfo parent %W] }
    pack $w -expand 1 -fill both

    #
    # Modify the list of binding tags of the body text widget
    #
    bindtags $w [list $w TablelistBody [winfo toplevel $w] TablelistKeyNav all]

    #
    # Create the "active", "stripe", "select", and "disabled" tags
    # in the body text widget.  Don't use the built-in "sel" tag
    # because on Windows the selection in a text widget only
    # becomes visible when the window gets the input focus.
    #
    $w tag configure stripe -background "" -foreground ""
    $w tag configure active -borderwidth 1 -underline 1
    $w tag configure select -relief raised
    $w tag configure disabled -underline 0

    #
    # Create an unmanaged listbox child, used to handle the -setgrid option
    #
    listbox $data(lb)

    #
    # Configure the widget according to the command-line
    # arguments and to the available database options
    #
    if {[catch {
	mwutil::configure $win configSpecs data tablelist::doConfig \
			  [lrange $args 1 end] 1
    } result] != 0} {
	destroy $win
	return -code error $result
    }

    #
    # Move the original widget command into the current namespace
    # and build a new widget procedure in the global one
    #
    rename ::$win $win
    proc ::$win args [format {
	if {[catch {tablelist::tablelistWidgetCmd %s $args} result] == 0} {
	    return $result
	} else {
	    return -code error $result
	}
    } [list $win]]

    #
    # Register a callback to be invoked whenever the PRIMARY selection is
    # owned by the window win and someone attempts to retrieve it as a
    # UTF8_STRING or STRING (the type UTF8_STRING is only needed to work
    # around a bug in Tk 8.4.0 and 8.4.1 causing crashes under KDE 3.0)
    #
    selection handle -type UTF8_STRING $win \
	[list ::tablelist::fetchSelection $win]
    selection handle -type STRING $win \
	[list ::tablelist::fetchSelection $win]

    #
    # Set a trace on the array element data(activeIdx)
    #
    trace variable data(activeIdx) w [list tablelist::activeIdxTrace $win]

    return $win
}

#
# Private configuration procedures
# ================================
#
# See the module "tablelistConfig.tcl"
#

#
# Private procedures implementing the tablelist widget command
# ============================================================
#

#------------------------------------------------------------------------------
# tablelist::tablelistWidgetCmd
#
# This procedure is invoked to process the Tcl command corresponding to a
# tablelist widget.
#------------------------------------------------------------------------------
proc tablelist::tablelistWidgetCmd {win argList} {
    variable cmdOpts
    upvar ::tablelist::ns${win}::data data

    set argCount [llength $argList]
    if {$argCount == 0} {
	mwutil::wrongNumArgs "$win option ?arg arg ...?"
    }

    set cmd [mwutil::fullOpt "option" [lindex $argList 0] $cmdOpts]
    switch $cmd {
	activate -
	bbox -
	see {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd index"
	    }

	    synchronize $win
	    set index [rowIndex $win [lindex $argList 1] 0]
	    return [${cmd}SubCmd $win $index]
	}

	attrib {
	    return [mwutil::attribSubCmd $win [lrange $argList 1 end]]
	}

	bodypath {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    return $data(body)
	}

	cancelediting -
	curselection -
	finishediting {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    synchronize $win
	    return [${cmd}SubCmd $win]
	}

	cellcget {
	    if {$argCount != 3} {
		mwutil::wrongNumArgs "$win $cmd cellIndex option"
	    }

	    synchronize $win
	    scan [cellIndex $win [lindex $argList 1] 1] %d,%d row col
	    variable cellConfigSpecs
	    set opt [mwutil::fullConfigOpt [lindex $argList 2] cellConfigSpecs]
	    return [doCellCget $row $col $win $opt]
	}

	cellconfigure {
	    if {$argCount < 2} {
		mwutil::wrongNumArgs "$win $cmd cellIndex ?option? ?value?\
				      ?option value ...?"
	    }

	    synchronize $win
	    scan [cellIndex $win [lindex $argList 1] 1] %d,%d row col
	    variable cellConfigSpecs
	    set argList [lrange $argList 2 end]
	    mwutil::setConfigVals $win cellConfigSpecs cellConfigVals \
				  "tablelist::doCellCget $row $col" $argList
	    return [mwutil::configSubCmd $win cellConfigSpecs cellConfigVals \
		    "tablelist::doCellConfig $row $col" $argList]
	}

	cellindex {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd cellIndex"
	    }

	    synchronize $win
	    return [cellIndex $win [lindex $argList 1] 0]
	}

	cget {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd option"
	    }

	    #
	    # Return the value of the specified configuration option
	    #
	    variable configSpecs
	    set opt [mwutil::fullConfigOpt [lindex $argList 1] configSpecs]
	    return $data($opt)
	}

	columncget {
	    if {$argCount != 3} {
		mwutil::wrongNumArgs "$win $cmd columnIndex option"
	    }

	    synchronize $win
	    set col [colIndex $win [lindex $argList 1] 1]
	    variable colConfigSpecs
	    set opt [mwutil::fullConfigOpt [lindex $argList 2] colConfigSpecs]
	    return [doColCget $col $win $opt]
	}

	columnconfigure {
	    if {$argCount < 2} {
		mwutil::wrongNumArgs "$win $cmd columnIndex ?option? ?value?\
				      ?option value ...?"
	    }

	    synchronize $win
	    set col [colIndex $win [lindex $argList 1] 1]
	    variable colConfigSpecs
	    set argList [lrange $argList 2 end]
	    mwutil::setConfigVals $win colConfigSpecs colConfigVals \
				  "tablelist::doColCget $col" $argList
	    return [mwutil::configSubCmd $win colConfigSpecs colConfigVals \
		    "tablelist::doColConfig $col" $argList]
	}

	columncount {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    return $data(colCount)
	}

	columnindex {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd columnIndex"
	    }

	    synchronize $win
	    return [colIndex $win [lindex $argList 1] 0]
	}

	configure {
	    variable configSpecs
	    return [mwutil::configSubCmd $win configSpecs data \
		    tablelist::doConfig [lrange $argList 1 end]]
	}

	containing {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd y"
	    }

	    set y [lindex $argList 1]
	    format %d $y		;# integer check with error message
	    synchronize $win
	    return [containingSubCmd $win $y]
	}

	containingcell {
	    if {$argCount != 3} {
		mwutil::wrongNumArgs "$win $cmd x y"
	    }

	    set x [lindex $argList 1]
	    format %d $x		;# integer check with error message
	    set y [lindex $argList 2]
	    format %d $y		;# integer check with error message
	    synchronize $win
	    return [containingSubCmd $win $y],[containingcolumnSubCmd $win $x]
	}

	containingcolumn {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd x"
	    }

	    set x [lindex $argList 1]
	    format %d $x		;# integer check with error message
	    synchronize $win
	    return [containingcolumnSubCmd $win $x]
	}

	delete -
	get -
	getkeys {
	    if {$argCount < 2 || $argCount > 3} {
		mwutil::wrongNumArgs "$win $cmd firstIndex lastIndex" \
				     "$win $cmd indexList"
	    }

	    synchronize $win
	    set first [lindex $argList 1]
	    if {$argCount == 3} {
		set last [lindex $argList 2]
	    } else {
		set last $first
	    }
	    incr argCount -1
	    return [${cmd}SubCmd $win $first $last $argCount]
	}

	deletecolumns -
	getcolumns {
	    if {$argCount < 2 || $argCount > 3} {
		mwutil::wrongNumArgs "$win $cmd firstColumnIndex\
				      lastColumnIndex" \
				     "$win $cmd columnIndexList"
	    }

	    synchronize $win
	    set first [lindex $argList 1]
	    if {$argCount == 3} {
		set last [lindex $argList 2]
	    } else {
		set last $first
	    }
	    incr argCount -1
	    return [${cmd}SubCmd $win $first $last $argCount]
	}

	editcell {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd cellIndex"
	    }

	    synchronize $win
	    scan [cellIndex $win [lindex $argList 1] 1] %d,%d row col
	    return [${cmd}SubCmd $win $row $col 0]
	}

	editwinpath {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    if {[winfo exists $data(bodyFrEd)]} {
		return $data(bodyFrEd)
	    } else {
		return ""
	    }
	}

	entrypath {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    if {[winfo exists $data(bodyFrEd)]} {
		set class [winfo class $data(bodyFrEd)]
		if {[string compare $class Mentry] == 0 ||
		    [string compare $class Checkbutton] == 0} {
		    return ""
		} else {
		    return $data(editFocus)
		}
	    } else {
		return ""
	    }
	}

	fillcolumn {
	    if {$argCount != 3} {
		mwutil::wrongNumArgs "$win $cmd columnIndex text"
	    }

	    synchronize $win
	    set col [colIndex $win [lindex $argList 1] 1]
	    return [fillcolumnSubCmd $win $col [lindex $argList 2]]
	}

	index {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd index"
	    }

	    synchronize $win
	    return [rowIndex $win [lindex $argList 1] 1]
	}

	insert {
	    if {$argCount < 2} {
		mwutil::wrongNumArgs "$win $cmd index ?item item ...?"
	    }

	    synchronize $win
	    set index [rowIndex $win [lindex $argList 1] 1]
	    return [insertSubCmd $win $index [lrange $argList 2 end] \
		    $data(hasListVar)]
	}

	insertcolumnlist {
	    if {$argCount != 3} {
		mwutil::wrongNumArgs "$win $cmd columnIndex columnList"
	    }

	    synchronize $win
	    set arg1 [lindex $argList 1]
	    if {[string first $arg1 end] == 0 || $arg1 == $data(colCount)} {
		set col $data(colCount)
	    } else {
		set col [colIndex $win $arg1 1]
	    }
	    return [insertcolumnsSubCmd $win $col [lindex $argList 2]]
	}

	insertcolumns {
	    if {$argCount < 2} {
		mwutil::wrongNumArgs "$win $cmd columnIndex\
				      ?width title ?alignment?\
				       width title ?alignment? ...?"
	    }

	    synchronize $win
	    set arg1 [lindex $argList 1]
	    if {[string first $arg1 end] == 0 || $arg1 == $data(colCount)} {
		set col $data(colCount)
	    } else {
		set col [colIndex $win $arg1 1]
	    }
	    return [insertcolumnsSubCmd $win $col [lrange $argList 2 end]]
	}

	insertlist {
	    if {$argCount != 3} {
		mwutil::wrongNumArgs "$win $cmd index list"
	    }

	    synchronize $win
	    set index [rowIndex $win [lindex $argList 1] 1]
	    return [insertSubCmd $win $index [lindex $argList 2] \
		    $data(hasListVar)]
	}

	labelpath {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd columnIndex"
	    }

	    synchronize $win
	    set col [colIndex $win [lindex $argList 1] 1]
	    return $data(hdrTxtFrLbl)$col
	}

	labels {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    set children [winfo children $data(hdrTxtFr)]
	    return [lrange [lsort -dictionary $children] 1 end]
	}

	move {
	    if {$argCount != 3} {
		mwutil::wrongNumArgs "$win $cmd sourceIndex targetIndex"
	    }

	    synchronize $win
	    set source [rowIndex $win [lindex $argList 1] 0]
	    set target [rowIndex $win [lindex $argList 2] 1]
	    return [moveSubCmd $win $source $target]
	}

	movecolumn {
	    if {$argCount != 3} {
		mwutil::wrongNumArgs "$win $cmd sourceColumnIndex\
				      targetColumnIndex"
	    }

	    synchronize $win
	    set arg1 [lindex $argList 1]
	    set source [colIndex $win $arg1 1]
	    set arg2 [lindex $argList 2]
	    if {[string first $arg2 end] == 0 || $arg2 == $data(colCount)} {
		set target $data(colCount)
	    } else {
		set target [colIndex $win $arg2 1]
	    }
	    return [movecolumnSubCmd $win $source $target]
	}

	nearest {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd y"
	    }

	    set y [lindex $argList 1]
	    format %d $y		;# integer check with error message
	    synchronize $win
	    return [rowIndex $win @0,$y 0]
	}

	nearestcell {
	    if {$argCount != 3} {
		mwutil::wrongNumArgs "$win $cmd x y"
	    }

	    set x [lindex $argList 1]
	    format %d $x		;# integer check with error message
	    set y [lindex $argList 2]
	    format %d $y		;# integer check with error message
	    synchronize $win
	    return [cellIndex $win @$x,$y 0]
	}

	nearestcolumn {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd x"
	    }

	    set x [lindex $argList 1]
	    format %d $x		;# integer check with error message
	    synchronize $win
	    return [colIndex $win @$x,0 0]
	}

	rejectinput {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    set data(rejected) 1
	}

	resetsortinfo {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    set data(sortCol) -1
	    set data(sortOrder) {}

	    place forget $data(hdrTxtFrCanv)
	    set oldArrowCol $data(arrowCol)
	    set data(arrowCol) -1
	    synchronize $win
	    adjustColumns $win l$oldArrowCol 1
	    return ""
	}

	rowcget {
	    if {$argCount != 3} {
		mwutil::wrongNumArgs "$win $cmd index option"
	    }

	    #
	    # Check the row index
	    #
	    synchronize $win
	    set rowArg [lindex $argList 1]
	    set row [rowIndex $win $rowArg 0]
	    if {$row < 0 || $row > $data(lastRow)} {
		return -code error \
		       "row index \"$rowArg\" out of range"
	    }

	    variable rowConfigSpecs
	    set opt [mwutil::fullConfigOpt [lindex $argList 2] rowConfigSpecs]
	    return [doRowCget $row $win $opt]
	}

	rowconfigure {
	    if {$argCount < 2} {
		mwutil::wrongNumArgs "$win $cmd index ?option? ?value?\
				      ?option value ...?"
	    }

	    #
	    # Check the row index
	    #
	    synchronize $win
	    set rowArg [lindex $argList 1]
	    set row [rowIndex $win $rowArg 0]
	    if {$row < 0 || $row > $data(lastRow)} {
		return -code error \
		       "row index \"$rowArg\" out of range"
	    }

	    variable rowConfigSpecs
	    set argList [lrange $argList 2 end]
	    mwutil::setConfigVals $win rowConfigSpecs rowConfigVals \
				  "tablelist::doRowCget $row" $argList
	    return [mwutil::configSubCmd $win rowConfigSpecs rowConfigVals \
		    "tablelist::doRowConfig $row" $argList]
	}

	scan {
	    if {$argCount != 4} {
		mwutil::wrongNumArgs "$win $cmd mark|dragto x y"
	    }

	    set x [lindex $argList 2]
	    set y [lindex $argList 3]
	    format %d $x		;# integer check with error message
	    format %d $y		;# integer check with error message
	    variable scanCmdOpts
	    set opt [mwutil::fullOpt "option" [lindex $argList 1] $scanCmdOpts]
	    synchronize $win
	    return [scanSubCmd $win $opt $x $y]
	}

	seecell {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd cellIndex"
	    }

	    synchronize $win
	    scan [cellIndex $win [lindex $argList 1] 1] %d,%d row col
	    return [${cmd}SubCmd $win $row $col]
	}

	seecolumn {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd columnIndex"
	    }

	    synchronize $win
	    set col [colIndex $win [lindex $argList 1] 1]
	    return [seecellSubCmd $win [rowIndex $win @0,0 0] $col]
	}

	selection {
	    if {$argCount < 3 || $argCount > 4} {
		mwutil::wrongNumArgs "$win $cmd option firstIndex lastIndex" \
				     "$win $cmd indexList"
	    }

	    synchronize $win
	    variable selCmdOpts
	    set opt [mwutil::fullOpt "option" [lindex $argList 1] $selCmdOpts]
	    set first [lindex $argList 2]
	    switch $opt {
		anchor -
		includes {
		    if {$argCount != 3} {
			mwutil::wrongNumArgs "$win selection $opt index"
		    }
		    set index [rowIndex $win $first 0]
		    return [selectionSubCmd $win $opt $index $index]
		}
		clear -
		set {
		    if {$argCount == 3} {
			foreach elem $first {
			    set index [rowIndex $win $elem 0]
			    selectionSubCmd $win $opt $index $index
			}
			return ""
		    } else {
			set first [rowIndex $win $first 0]
			set last [rowIndex $win [lindex $argList 3] 0]
			return [selectionSubCmd $win $opt $first $last]
		    }
		}
	    }
	}

	separatorpath {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd columnIndex"
	    }

	    synchronize $win
	    set col [colIndex $win [lindex $argList 1] 1]
	    if {$data(-showseparators)} {
		return $data(sep)$col
	    } else {
		return ""
	    }
	}

	separators {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    set sepList {}
	    foreach w [winfo children $win] {
		if {[regexp {^sep[0-9]+$} [winfo name $w]]} {
		    lappend sepList $w
		}
	    }
	    return $sepList
	}

	size {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    synchronize $win
	    return $data(itemCount)
	}

	sort {
	    if {$argCount < 1 || $argCount > 2} {
		mwutil::wrongNumArgs "$win $cmd  ?-increasing|-decreasing?"
	    }

	    if {$argCount == 1} {
		set order -increasing
	    } else {
		variable sortOrders
		set order [mwutil::fullOpt "option" \
			   [lindex $argList 2] $sortOrders]
	    }
	    synchronize $win
	    return [sortSubCmd $win -1 $order]
	}

	sortbycolumn {
	    if {$argCount < 2 || $argCount > 3} {
		mwutil::wrongNumArgs "$win $cmd columnIndex\
				      ?-increasing|-decreasing?"
	    }

	    synchronize $win
	    set col [colIndex $win [lindex $argList 1] 1]
	    if {$argCount == 2} {
		set order -increasing
	    } else {
		variable sortOrders
		set order [mwutil::fullOpt "option" \
			   [lindex $argList 2] $sortOrders]
	    }
	    return [sortSubCmd $win $col $order]
	}

	sortcolumn {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    return $data(sortCol)
	}

	sortorder {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    return $data(sortOrder)
	}

	xview -
	yview {
	    synchronize $win
	    return [${cmd}SubCmd $win [lrange $argList 1 end]]
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::activateSubCmd
#
# This procedure is invoked to process the tablelist activate subcommand.
#------------------------------------------------------------------------------
proc tablelist::activateSubCmd {win index} {
    upvar ::tablelist::ns${win}::data data

    if {$data(isDisabled)} {
	return ""
    }

    #
    # Adjust the index to fit within the existing items
    #
    if {$index > $data(lastRow)} {
	set index $data(lastRow)
    }
    if {$index < 0} {
	set index 0
    }

    set data(activeIdx) $index
    return ""
}

#------------------------------------------------------------------------------
# tablelist::bboxSubCmd
#
# This procedure is invoked to process the tablelist bbox subcommand.
#------------------------------------------------------------------------------
proc tablelist::bboxSubCmd {win index} {
    upvar ::tablelist::ns${win}::data data

    set w $data(body)
    set dlineinfo [$w dlineinfo [expr {double($index + 1)}]]
    if {$data(itemCount) == 0 || [string compare $dlineinfo ""] == 0} {
	return {}
    }

    foreach {x y width height baselinePos} $dlineinfo {
	lappend bbox [expr {$x + [winfo x $w]}] \
		     [expr {$y + [winfo y $w] + $data(-selectborderwidth)}] \
		     $width [expr {$height - 2*$data(-selectborderwidth) - 1}]
    }
    return $bbox
}

#------------------------------------------------------------------------------
# tablelist::containingSubCmd
#
# This procedure is invoked to process the tablelist containing subcommand.
#------------------------------------------------------------------------------
proc tablelist::containingSubCmd {win y} {
    upvar ::tablelist::ns${win}::data data

    set row [rowIndex $win @0,$y 0]

    set w $data(body)
    incr y -[winfo y $w]
    set dlineinfo [$w dlineinfo [expr {double($row + 1)}]]
    if {$y < [lindex $dlineinfo 1] + [lindex $dlineinfo 3]} {
	return $row
    } else {
	return -1
    }
}

#------------------------------------------------------------------------------
# tablelist::containingcolumnSubCmd
#
# This procedure is invoked to process the tablelist containingcolumn
# subcommand.
#------------------------------------------------------------------------------
proc tablelist::containingcolumnSubCmd {win x} {
    upvar ::tablelist::ns${win}::data data

    set col [colIndex $win @$x,0 0]
    if {$col < 0} {
	return -1
    }

    set lbl $data(hdrTxtFrLbl)$col
    if {$x + [winfo rootx $win] < [winfo width $lbl] + [winfo rootx $lbl]} {
	return $col
    } else {
	return -1
    }
}

#------------------------------------------------------------------------------
# tablelist::curselectionSubCmd
#
# This procedure is invoked to process the tablelist curselection subcommand.
#------------------------------------------------------------------------------
proc tablelist::curselectionSubCmd win {
    upvar ::tablelist::ns${win}::data data

    #
    # Find the selected lines of the body text widget
    #
    set result {}
    set w $data(body)
    set selRange [$w tag nextrange select 1.0]
    while {[llength $selRange] != 0} {
	set selStart [lindex $selRange 0]
	set selEnd [lindex $selRange 1]
	lappend result [expr {int($selStart) - 1}]

	set selRange [$w tag nextrange select $selEnd]
    }
    return $result
}

#------------------------------------------------------------------------------
# tablelist::deleteSubCmd
#
# This procedure is invoked to process the tablelist delete subcommand.
#------------------------------------------------------------------------------
proc tablelist::deleteSubCmd {win first last argCount} {
    upvar ::tablelist::ns${win}::data data

    if {$data(isDisabled)} {
	return ""
    }

    if {$argCount == 1} {
	if {[llength $first] == 1} {			;# just to save time
	    set index [rowIndex $win [lindex $first 0] 0]
	    return [deleteRows $win $index $index $data(hasListVar)]
	} elseif {$data(itemCount) == 0} {		;# no items present
	    return ""
	} else {					;# a bit more work
	    #
	    # Sort the numerical equivalents of the
	    # specified indices in decreasing order
	    #
	    set indexList {}
	    foreach elem $first {
		set index [rowIndex $win $elem 0]
		if {$index < 0} {
		    set index 0
		} elseif {$index > $data(lastRow)} {
		    set index $data(lastRow)
		}
		lappend indexList $index
	    }
	    set indexList [lsort -integer -decreasing $indexList]

	    #
	    # Traverse the sorted index list and ignore any duplicates
	    #
	    set prevIndex -1
	    foreach index $indexList {
		if {$index != $prevIndex} {
		    deleteRows $win $index $index $data(hasListVar)
		    set prevIndex $index
		}
	    }
	    return ""
	}
    } else {
	set first [rowIndex $win $first 0]
	set last [rowIndex $win $last 0]
	return [deleteRows $win $first $last $data(hasListVar)]
    }
}

#------------------------------------------------------------------------------
# tablelist::deleteRows
#
# Deletes a given range of rows of a tablelist widget.
#------------------------------------------------------------------------------
proc tablelist::deleteRows {win first last updateListVar} {
    upvar ::tablelist::ns${win}::data data

    #
    # Adjust the range to fit within the existing items
    #
    if {$first < 0} {
	set first 0
    }
    if {$last > $data(lastRow)} {
	set last $data(lastRow)
    }
    set count [expr {$last - $first + 1}]
    if {$count <= 0} {
	return ""
    }

    #
    # Check whether the width of any dynamic-width
    # column might be affected by the deletion
    #
    set w $data(body)
    set itemListRange [lrange $data(itemList) $first $last]
    if {$count == $data(itemCount)} {
	set colWidthsChanged 1				;# just to save time
    } else {
	set colWidthsChanged 0
	set snipStr $data(-snipstring)
	foreach item $itemListRange {
	    set dispItem [strToDispStr $item]
	    set key [lindex $item end]
	    set col 0
	    foreach text [lrange $dispItem 0 $data(lastCol)] \
		    colFont $data(colFontList) \
		    fmtCmdFlag $data(fmtCmdFlagList) \
		    {pixels alignment} $data(colList) {
		if {$data($col-hide) || $pixels != 0} {
		    incr col
		    continue
		}

		if {$fmtCmdFlag} {
		    set text [uplevel #0 $data($col-formatcommand) \
			      [list [lindex $item $col]]]
		    set text [strToDispStr $text]
		}
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
		adjustElem $win text imageWidth $cellFont \
			   $pixels $alignment $snipStr
		set textWidth [font measure $cellFont -displayof $win $text]
		set elemWidth [expr {$imageWidth + $textWidth}]
		if {$elemWidth == $data($col-elemWidth) &&
		    [incr data($col-widestCount) -1] == 0} {
		    set colWidthsChanged 1
		    break
		}

		incr col
	    }

	    if {$colWidthsChanged} {
		break
	    }
	}
    }

    #
    # Delete the given items and their tags from the body text widget.
    # Interestingly, for a large number of items it is much more efficient to
    # delete each line individually than to invoke a global delete command.
    #
    set textIdx1 [expr {double($first + 1)}]
    set textIdx2 [expr {double($first + 2)}]
    foreach item $itemListRange {
	$w delete $textIdx1 $textIdx2

	set key [lindex $item end]
	array set tagData [array get data $key*-\[bfs\]*]	;# for speed

	foreach tag [array names tagData $key-\[bfs\]*] {
	    $w tag delete $tag
	    unset data($tag)
	    if {[string match $key-\[bf\]* $tag]} {
		incr data(tagCount) -1
	    }
	}

	for {set col 0} {$col < $data(colCount)} {incr col} {
	    foreach tag [array names tagData $key-$col-\[bfs\]*] {
		$w tag delete $tag
		unset data($tag)
		if {[string match $key-$col-\[bf\]* $tag]} {
		    incr data(tagCount) -1
		}
	    }
	    if {[info exists data($key-$col-image)]} {
		unset data($key-$col-image)
		incr data(imgCount) -1
	    }
	}

	unset tagData
    }

    #
    # Delete the given items from the internal list
    #
    set data(itemList) [lreplace $data(itemList) $first $last]
    incr data(itemCount) -$count
    incr data(lastRow) -$count

    #
    # Delete the given items from the list variable if needed
    #
    if {$updateListVar} {
	trace vdelete ::$data(-listvariable) wu $data(listVarTraceCmd)
	upvar #0 $data(-listvariable) var
	set var [lreplace $var $first $last]
	trace variable ::$data(-listvariable) wu $data(listVarTraceCmd)
    }

    #
    # Adjust the heights of the body text widget
    # and of the listbox child, if necessary
    #
    if {$data(-height) <= 0} {
	$w configure -height $data(itemCount)
	$data(lb) configure -height $data(itemCount)
    }

    #
    # Adjust the columns if necessary, adjust the separators,
    # and redraw the stripes in the body text widget
    #
    if {$colWidthsChanged} {
	adjustColumns $win allCols 1
    }
    adjustSepsWhenIdle $win
    makeStripesWhenIdle $win

    #
    # Update the indices anchorIdx and activeIdx
    #
    if {$first <= $data(anchorIdx)} {
	incr data(anchorIdx) -$count
	if {$data(anchorIdx) < $first} {
	    set data(anchorIdx) $first
	}
    }
    if {$last < $data(activeIdx)} {
	incr data(activeIdx) -$count
    } elseif {$first <= $data(activeIdx)} {
	set data(activeIdx) $first
	if {$data(activeIdx) > $data(lastRow) && $data(lastRow) >= 0} {
	    set data(activeIdx) $data(lastRow)
	}
    }

    #
    # Update data(editRow) if the edit window is present
    #
    if {$data(editRow) >= 0} {
	set data(editRow) [lsearch $data(itemList) "* $data(editKey)"]
    }

    return ""
}

#------------------------------------------------------------------------------
# tablelist::deletecolumnsSubCmd
#
# This procedure is invoked to process the tablelist deletecolumns subcommand.
#------------------------------------------------------------------------------
proc tablelist::deletecolumnsSubCmd {win first last argCount} {
    upvar ::tablelist::ns${win}::data data

    if {$data(isDisabled)} {
	return ""
    }

    if {$argCount == 1} {
	if {[llength $first] == 1} {			;# just to save time
	    set col [colIndex $win [lindex $first 0] 1]
	    deleteCols $win $col $col
	    redisplay $win
	} elseif {$data(colCount) == 0} {		;# no columns present
	    return ""
	} else {					;# a bit more work
	    #
	    # Sort the numerical equivalents of the
	    # specified column indices in decreasing order
	    #
	    set colList {}
	    foreach elem $first {
		lappend colList [colIndex $win $elem 1]
	    }
	    set colList [lsort -integer -decreasing $colList]

	    #
	    # Traverse the sorted column index
	    # list and ignore any duplicates
	    #
	    set deleted 0
	    set prevCol -1
	    foreach col $colList {
		if {$col != $prevCol} {
		    deleteCols $win $col $col
		    set deleted 1
		    set prevCol $col
		}
	    }
	    if {$deleted} {
		redisplay $win
	    }
	}
    } else {
	set first [colIndex $win $first 1]
	set last [colIndex $win $last 1]
	if {$first <= $last} {
	    deleteCols $win $first $last
	    redisplay $win
	}
    }

    return ""
}

#------------------------------------------------------------------------------
# tablelist::deleteCols
#
# Deletes a given range of columns of a tablelist widget.
#------------------------------------------------------------------------------
proc tablelist::deleteCols {win first last} {
    upvar ::tablelist::ns${win}::data data

    #
    # Delete the data corresponding to the given range
    #
    for {set col $first} {$col <= $last} {incr col} {
	deleteColData $win $col
    }

    #
    # Shift the elements of data corresponding to the column
    # indices > last to the left by last - first + 1 positions
    #
    for {set oldCol [expr {$last + 1}]; set newCol $first} \
	{$oldCol < $data(colCount)} {incr oldCol; incr newCol} {
	moveColData $win data data imgs $oldCol $newCol
    }

    #
    # Update the item list
    #
    set newItemList {}
    foreach item $data(itemList) {
	set item [lreplace $item $first $last]
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
    setupColumns $win \
	[lreplace $data(-columns) [expr {3*$first}] [expr {3*$last + 2}]] 1
    makeColFontAndTagLists $win
    adjustColumns $win {} 1

    #
    # Reconfigure the relevant column labels
    #
    for {set col $first} {$col < $data(colCount)} {incr col} {
	reconfigColLabels $win imgs $col
    }
}

#------------------------------------------------------------------------------
# tablelist::fillcolumnSubCmd
#
# This procedure is invoked to process the tablelist fillcolumn subcommand.
#------------------------------------------------------------------------------
proc tablelist::fillcolumnSubCmd {win colIdx text} {
    upvar ::tablelist::ns${win}::data data

    if {$data(isDisabled)} {
	return ""
    }

    #
    # Update the item list
    #
    set newItemList {}
    foreach item $data(itemList) {
	set item [lreplace $item $colIdx $colIdx $text]
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
    adjustColumns $win $colIdx 1
    redisplayWhenIdle $win
    return ""
}

#------------------------------------------------------------------------------
# tablelist::getSubCmd
#
# This procedure is invoked to process the tablelist get subcommand.
#------------------------------------------------------------------------------
proc tablelist::getSubCmd {win first last argCount} {
    upvar ::tablelist::ns${win}::data data

    #
    # Get the specified items from the internal list
    #
    set result {}
    if {$argCount == 1} {
	foreach elem $first {
	    set index [rowIndex $win $elem 0]
	    if {$index >= 0 && $index < $data(itemCount)} {
		set item [lindex $data(itemList) $index]
		lappend result [lrange $item 0 $data(lastCol)]
	    }
	}

	if {[llength $first] == 1} {
	    return [lindex $result 0]
	} else {
	    return $result
	}
    } else {
	set first [rowIndex $win $first 0]
	set last [rowIndex $win $last 0]

	#
	# Adjust the range to fit within the existing items
	#
	if {$first > $data(lastRow)} {
	    return {}
	}
	if {$last > $data(lastRow)} {
	    set last $data(lastRow)
	}
	if {$first < 0} {
	    set first 0
	}

	foreach item [lrange $data(itemList) $first $last] {
	    lappend result [lrange $item 0 $data(lastCol)]
	}
	return $result
    }
}

#------------------------------------------------------------------------------
# tablelist::getcolumnsSubCmd
#
# This procedure is invoked to process the tablelist getcolumns subcommand.
#------------------------------------------------------------------------------
proc tablelist::getcolumnsSubCmd {win first last argCount} {
    upvar ::tablelist::ns${win}::data data

    #
    # Get the specified columns from the internal list
    #
    set result {}
    if {$argCount == 1} {
	foreach elem $first {
	    set col [colIndex $win $elem 1]
	    set colResult {}
	    foreach item $data(itemList) {
		lappend colResult [lindex $item $col]
	    }
	    lappend result $colResult
	}

	if {[llength $first] == 1} {
	    return [lindex $result 0]
	} else {
	    return $result
	}
    } else {
	set first [colIndex $win $first 1]
	set last [colIndex $win $last 1]

	if {$first > $last} {
	    return {}
	}

	for {set col $first} {$col <= $last} {incr col} {
	    set colResult {}
	    foreach item $data(itemList) {
		lappend colResult [lindex $item $col]
	    }
	    lappend result $colResult
	}
	return $result
    }
}

#------------------------------------------------------------------------------
# tablelist::getkeysSubCmd
#
# This procedure is invoked to process the tablelist getkeys subcommand.
#------------------------------------------------------------------------------
proc tablelist::getkeysSubCmd {win first last argCount} {
    upvar ::tablelist::ns${win}::data data

    #
    # Get the specified keys from the internal list
    #
    set result {}
    if {$argCount == 1} {
	foreach elem $first {
	    set index [rowIndex $win $elem 0]
	    if {$index >= 0 && $index < $data(itemCount)} {
		set item [lindex $data(itemList) $index]
		lappend result [string range [lindex $item end] 1 end]
	    }
	}

	if {[llength $first] == 1} {
	    return [lindex $result 0]
	} else {
	    return $result
	}
    } else {
	set first [rowIndex $win $first 0]
	set last [rowIndex $win $last 0]

	#
	# Adjust the range to fit within the existing items
	#
	if {$first > $data(lastRow)} {
	    return {}
	}
	if {$last > $data(lastRow)} {
	    set last $data(lastRow)
	}
	if {$first < 0} {
	    set first 0
	}

	foreach item [lrange $data(itemList) $first $last] {
	    lappend result [string range [lindex $item end] 1 end]
	}
	return $result
    }
}

#------------------------------------------------------------------------------
# tablelist::insertSubCmd
#
# This procedure is invoked to process the tablelist insert and insertlist
# subcommands.
#------------------------------------------------------------------------------
proc tablelist::insertSubCmd {win index argList updateListVar} {
    upvar ::tablelist::ns${win}::data data

    if {$data(isDisabled)} {
	return ""
    }

    set argCount [llength $argList]
    if {$argCount == 0} {
	return ""
    }

    if {$index < 0} {
	set index 0
    }

    #
    # Insert the items into the body text widget and into the internal list
    #
    set w $data(body)
    set widgetFont $data(-font)
    set snipStr $data(-snipstring)
    set savedCount $data(itemCount)
    set colWidthsChanged 0
    set idx $index
    set line [expr {$index + 1}]
    foreach item $argList {
	set item [adjustItem $item $data(colCount)]
	if {$data(itemCount) != 0} {
	    $w insert $line.0 \n
	}
	set col 0

	if {$data(hasColTags)} {
	    set insertArgs {}
	    foreach text [strToDispStr $item] \
		    colFont $data(colFontList) \
		    colTags $data(colTagsList) \
		    fmtCmdFlag $data(fmtCmdFlagList) \
		    {pixels alignment} $data(colList) {
		if {$data($col-hide)} {
		    incr col
		    continue
		}

		#
		# Update the column width or clip the element if necessary
		#
		if {$fmtCmdFlag} {
		    set text [uplevel #0 $data($col-formatcommand) \
			      [list [lindex $item $col]]]
		    set text [strToDispStr $text]
		}
		if {$pixels == 0} {		;# convention: dynamic width
		    set textWidth \
			[font measure $colFont -displayof $win $text]
		    if {$data($col-maxPixels) > 0 &&
			$textWidth > $data($col-maxPixels)} {
			set pixels $data($col-maxPixels)
		    }
		    if {$textWidth == $data($col-elemWidth)} {
			incr data($col-widestCount)
		    } elseif {$textWidth > $data($col-elemWidth)} {
			set data($col-elemWidth) $textWidth
			set data($col-widestCount) 1
			if {$textWidth > $data($col-reqPixels)} {
			    set data($col-reqPixels) $textWidth
			    if {$pixels == 0} {
				set colWidthsChanged 1
			    }
			}
		    }
		}
		if {$pixels != 0} {
		    incr pixels $data($col-delta)
		    set text [strRangeExt $win $text $colFont \
			      $pixels $alignment $snipStr]
		}

		lappend insertArgs \t$text\t $colTags
		incr col
	    }

	    #
	    # Insert the item into the body text widget
	    #
	    if {[llength $insertArgs] != 0} {
		eval [list $w insert $line.0] $insertArgs
	    }

	} else {
	    set insertStr ""
	    foreach text [strToDispStr $item] \
		    fmtCmdFlag $data(fmtCmdFlagList) \
		    {pixels alignment} $data(colList) {
		if {$data($col-hide)} {
		    incr col
		    continue
		}

		#
		# Update the column width or clip the element if necessary
		#
		if {$fmtCmdFlag} {
		    set text [uplevel #0 $data($col-formatcommand) \
			      [list [lindex $item $col]]]
		    set text [strToDispStr $text]
		}
		if {$pixels == 0} {		;# convention: dynamic width
		    set textWidth \
			[font measure $widgetFont -displayof $win $text]
		    if {$data($col-maxPixels) > 0 &&
			$textWidth > $data($col-maxPixels)} {
			set pixels $data($col-maxPixels)
		    }
		    if {$textWidth == $data($col-elemWidth)} {
			incr data($col-widestCount)
		    } elseif {$textWidth > $data($col-elemWidth)} {
			set data($col-elemWidth) $textWidth
			set data($col-widestCount) 1
			if {$textWidth > $data($col-reqPixels)} {
			    set data($col-reqPixels) $textWidth
			    if {$pixels == 0} {
				set colWidthsChanged 1
			    }
			}
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
	}

	#
	# Insert the item into the list variable if needed
	#
	if {$updateListVar} {
	    trace vdelete ::$data(-listvariable) wu $data(listVarTraceCmd)
	    upvar #0 $data(-listvariable) var
	    if {$idx == $data(itemCount)} {
		lappend var $item		;# this works much faster
	    } else {
		set var [linsert $var $idx $item]
	    }
	    trace variable ::$data(-listvariable) wu $data(listVarTraceCmd)
	}

	#
	# Insert the item into the internal list
	#
	lappend item k[incr data(seqNum)]
	if {$idx == $data(itemCount)} {
	    lappend data(itemList) $item	;# this works much faster
	} else {
	    set data(itemList) [linsert $data(itemList) $idx $item]
	}

	incr idx
	incr line
	incr data(itemCount)
    }
    set data(lastRow) [expr {$data(itemCount) - 1}]

    #
    # Adjust the heights of the body text widget
    # and of the listbox child, if necessary
    #
    if {$data(-height) <= 0} {
	$w configure -height $data(itemCount)
	$data(lb) configure -height $data(itemCount)
    }

    #
    # Adjust the horizontal view in the body text 
    # widget if the tablelist was previously empty
    #
    if {$savedCount == 0} {
	$w xview moveto [lindex [$data(hdrTxt) xview] 0]
    }

    #
    # Adjust the columns if necessary, adjust the separators,
    # and redraw the stripes in the body text widget
    #
    if {$colWidthsChanged} {
	adjustColumns $win {} 1
    }
    adjustSepsWhenIdle $win
    makeStripesWhenIdle $win

    #
    # Update the indices anchorIdx and activeIdx
    #
    if {$index <= $data(anchorIdx)} {
	incr data(anchorIdx) $argCount
    }
    if {$index <= $data(activeIdx)} {
	incr data(activeIdx) $argCount
	if {$data(activeIdx) > $data(lastRow) && $data(lastRow) >= 0} {
	    set data(activeIdx) $data(lastRow)
	}
    }

    #
    # Update data(editRow) if the edit window is present
    #
    if {$data(editRow) >= 0} {
	set data(editRow) [lsearch $data(itemList) "* $data(editKey)"]
    }

    return ""
}

#------------------------------------------------------------------------------
# tablelist::insertcolumnsSubCmd
#
# This procedure is invoked to process the tablelist insertcolumns and
# insertcolumnlist subcommands.
#------------------------------------------------------------------------------
proc tablelist::insertcolumnsSubCmd {win colIdx argList} {
    variable alignments
    upvar ::tablelist::ns${win}::data data

    if {$data(isDisabled)} {
	return ""
    }

    #
    # Check the syntax of argList and get the number of columns to be inserted
    #
    set count 0
    set argCount [llength $argList]
    for {set n 0} {$n < $argCount} {incr n} {
	#
	# Check the column width
	#
	format %d [lindex $argList $n]	;# integer check with error message

	#
	# Check whether the column title is present
	#
	if {[incr n] == $argCount} {
	    return -code error "column title missing"
	}

	#
	# Check the column alignment
	#
	set alignment left
	if {[incr n] < $argCount} {
	    set next [lindex $argList $n]
	    if {[catch {format %d $next}] == 0} {	;# integer check
		incr n -1
	    } else {
		mwutil::fullOpt "alignment" $next $alignments
	    }
	}

	incr count
    }

    #
    # Shift the elements of data corresponding to the column
    # indices >= colIdx to the right by count positions
    #
    for {set oldCol $data(lastCol); set newCol [expr {$oldCol + $count}]} \
	{$oldCol >= $colIdx} {incr oldCol -1; incr newCol -1} {
	moveColData $win data data imgs $oldCol $newCol
    }

    #
    # Update the item list
    #
    set emptyStrs {}
    for {set n 0} {$n < $count} {incr n} {
	lappend emptyStrs ""
    }
    set newItemList {}
    foreach item $data(itemList) {
	set item [eval [list linsert $item $colIdx] $emptyStrs]
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
    setupColumns $win \
	[eval [list linsert $data(-columns) [expr {3*$colIdx}]] $argList] 1
    makeColFontAndTagLists $win
    set limit [expr {$colIdx + $count}]
    set cols {}
    for {set col $colIdx} {$col < $limit} {incr col} {
	lappend cols $col
    }
    adjustColumns $win $cols 1

    #
    # Reconfigure the relevant column labels
    #
    for {set col $limit} {$col < $data(colCount)} {incr col} {
	reconfigColLabels $win imgs $col
    }

    #
    # Redisplay the columns
    #
    redisplay $win
    return ""
}

#------------------------------------------------------------------------------
# tablelist::scanSubCmd
#
# This procedure is invoked to process the tablelist scan subcommand.
#------------------------------------------------------------------------------
proc tablelist::scanSubCmd {win opt x y} {
    upvar ::tablelist::ns${win}::data data

    set w $data(body)
    incr x -[winfo x $w]
    incr y -[winfo y $w]

    $w scan $opt $x $y
    $data(hdrTxt) scan $opt $x $y
    return ""
}

#------------------------------------------------------------------------------
# tablelist::seeSubCmd
#
# This procedure is invoked to process the tablelist see subcommand.
#------------------------------------------------------------------------------
proc tablelist::seeSubCmd {win index} {
    upvar ::tablelist::ns${win}::data data

    #
    # Adjust the view in the body text widget
    #
    set w $data(body)
    set fraction [lindex [$w xview] 0]
    $w see [expr {double($index + 1)}]
    $w xview moveto $fraction
    return ""
}

#------------------------------------------------------------------------------
# tablelist::seecellSubCmd
#
# This procedure is invoked to process the tablelist seecell subcommand.
#------------------------------------------------------------------------------
proc tablelist::seecellSubCmd {win row col} {
    upvar ::tablelist::ns${win}::data data

    if {$data($col-hide)} {
	return ""
    }

    set alignment [lindex $data(colList) [expr {2*$col + 1}]]
    set w $data(body)
    findCellTabs $win [expr {$row + 1}] $col tabIdx1 tabIdx2
    set nextIdx [$w index $tabIdx2+1c]

    if {[string compare $alignment right] == 0} {
	$w see $nextIdx

	#
	# Shift the view in the body text widget until the first tab
	# becomes visible but finish the scrolling before the character
	# (\t or \n) at the position nextIdx would become invisible
	#
	if {![isCharVisible $w $tabIdx1]} {
	    while 1 {
		$w xview scroll -1 units
		if {![isCharVisible $w $nextIdx]} {
		    $w xview scroll 1 units
		    break
		} elseif {[isCharVisible $w $tabIdx1]} {
		    break
		}
	    }
	}
    } else {
	$w see $tabIdx1

	#
	# Shift the view in the body text widget until the character
	# (\t or \n) at the position nextIdx becomes visible but finish
	# the scrolling before the first tab would become invisible
	#
	if {![isCharVisible $w $nextIdx]} {
	    while 1 {
		$w xview scroll 1 units
		if {![isCharVisible $w $tabIdx1]} {
		    $w xview scroll -1 units
		    break
		} elseif {[isCharVisible $w $nextIdx]} {
		    break
		}
	    }
	}
    }

    $data(hdrTxt) xview moveto [lindex [$w xview] 0]
    return ""
}

#------------------------------------------------------------------------------
# tablelist::selectionSubCmd
#
# This procedure is invoked to process the tablelist selection subcommand.
#------------------------------------------------------------------------------
proc tablelist::selectionSubCmd {win opt first last} {
    upvar ::tablelist::ns${win}::data data

    if {$data(isDisabled) && [string compare $opt includes] != 0} {
	return ""
    }

    switch $opt {
	anchor {
	    #
	    # Adjust the index to fit within the existing items
	    #
	    if {$first > $data(lastRow)} {
		set first $data(lastRow)
	    }
	    if {$first < 0} {
		set first 0
	    }

	    set data(anchorIdx) $first
	    return ""
	}

	clear {
	    #
	    # Swap the indices if necessary
	    #
	    if {$last < $first} {
		set tmp $first
		set first $last
		set last $tmp
	    }

	    #
	    # Find the selected lines of the body text widget
	    # in the text range specified by the two indices
	    #
	    set w $data(body)
	    set firstTextIdx [expr {$first + 1}].0
	    set lastTextIdx [expr {$last + 1}].end
	    set selRange [$w tag nextrange select $firstTextIdx $lastTextIdx]
	    while {[llength $selRange] != 0} {
		set selStart [lindex $selRange 0]
		set selEnd [lindex $selRange 1]

		$w tag remove select $selStart $selEnd

		#
		# Handle the -(select)background and -(select)foreground cell
		# and column configuration options for each element of the row
		#
		set item [lindex $data(itemList) [expr {int($selStart) - 1}]]
		set key [lindex $item end]
		set textIdx1 $selStart
		for {set col 0} {$col < $data(colCount)} {incr col} {
		    if {$data($col-hide)} {
			continue
		    }

		    set textIdx2 \
			[$w search \t $textIdx1+1c "$selStart lineend"]+1c
		    foreach optTail {background foreground} {
			foreach tag [list $col-select$optTail \
				     $key-select$optTail \
				     $key-$col-select$optTail] {
			    if {[info exists data($tag)]} {
				$w tag remove $tag $textIdx1 $textIdx2
			    }
			}
			foreach tag [list $col-$optTail $key-$optTail \
				     $key-$col-$optTail] {
			    if {[info exists data($tag)]} {
				$w tag add $tag $textIdx1 $textIdx2
			    }
			}
		    }
		    set textIdx1 $textIdx2
		}

		set selRange [$w tag nextrange select $selEnd $lastTextIdx]
	    }

	    return ""
	}

	includes {
	    set tagNames [$data(body) tag names [expr {double($first + 1)}]]
	    if {[lsearch -exact $tagNames select] >= 0} {
		return 1
	    } else {
		return 0
	    }
	}

	set {
	    #
	    # Swap the indices if necessary and adjust
	    # the range to fit within the existing items
	    #
	    if {$last < $first} {
		set tmp $first
		set first $last
		set last $tmp
	    }
	    if {$first < 0} {
		set first 0
	    }
	    if {$last > $data(lastRow)} {
		set last $data(lastRow)
	    }

	    set w $data(body)
	    for {set idx $first; set line [expr {$first + 1}]} \
		{$idx <= $last} {incr idx; incr line} {
		#
		# Nothing to do if the row is already selected
		#
		if {[lsearch -exact [$w tag names $line.0] select] >= 0} {
		    continue
		}

		#
		# Check whether the row is selectable
		#
		set item [lindex $data(itemList) $idx]
		set key [lindex $item end]
		if {[info exists data($key-selectable)]} {    ;# not selectable
		    continue
		}

		$w tag add select $line.0 $line.end

		#
		# Handle the -(select)background and -(select)foreground cell
		# and column configuration options for each element of the row
		#
		set textIdx1 $line.0
		for {set col 0} {$col < $data(colCount)} {incr col} {
		    if {$data($col-hide)} {
			continue
		    }

		    set textIdx2 [$w search \t $textIdx1+1c $line.end]+1c
		    foreach optTail {background foreground} {
			foreach tag [list $col-select$optTail \
				     $key-select$optTail \
				     $key-$col-select$optTail] {
			    if {[info exists data($tag)]} {
				$w tag add $tag $textIdx1 $textIdx2
			    }
			}
			foreach tag [list $col-$optTail $key-$optTail \
				     $key-$col-$optTail] {
			    if {[info exists data($tag)]} {
				$w tag remove $tag $textIdx1 $textIdx2
			    }
			}
		    }
		    set textIdx1 $textIdx2
		}
	    }

	    #
	    # If the selection is exported and there are any selected
	    # rows in the widget then make win the new owner of the
	    # PRIMARY selection and register a callback to be invoked
	    # when it loses ownership of the PRIMARY selection
	    #
	    if {$data(-exportselection) &&
		[llength [$w tag nextrange select 1.0]] != 0} {
		selection own -command \
			[list ::tablelist::lostSelection $win] $win
	    }

	    return ""
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::xviewSubCmd
#
# This procedure is invoked to process the tablelist xview subcommand.
#------------------------------------------------------------------------------
proc tablelist::xviewSubCmd {win argList} {
    variable winSys
    upvar ::tablelist::ns${win}::data data

    switch [llength $argList] {
	0 {
	    #
	    # Command: $win xview
	    #
	    return [$data(hdrTxt) xview]
	}

	1 {
	    #
	    # Command: $win xview units
	    #
	    set units [lindex $argList 0]
	    format %d $units		;# integer check with error message
	    foreach w [list $data(hdrTxt) $data(body)] {
		$w xview moveto 0
		$w xview scroll $units units
	    }
	    if {[string compare $winSys aqua] == 0} {
		update			;# because of a Tk bug on Mac OS X Aqua
	    }
	    return ""
	}

	default {
	    #
	    # Command: $win xview moveto fraction
	    #	       $win xview scroll number what
	    #
	    foreach w [list $data(hdrTxt) $data(body)] {
		eval [list $w xview] $argList
	    }
	    if {[string compare $winSys aqua] == 0} {
		update			;# because of a Tk bug on Mac OS X Aqua
	    }
	    return ""
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::yviewSubCmd
#
# This procedure is invoked to process the tablelist yview subcommand.
#------------------------------------------------------------------------------
proc tablelist::yviewSubCmd {win argList} {
    upvar ::tablelist::ns${win}::data data
    variable winSys

    set w $data(body)
    set argCount [llength $argList]
    switch $argCount {
	0 {
	    #
	    # Command: $win yview
	    #
	    return [$w yview]
	}

	1 {
	    #
	    # Command: $win yview index
	    #
	    set index [rowIndex $win [lindex $argList 0] 0]
	    $w yview $index
	    $w xview moveto [lindex [$w xview] 0]
	    if {[string compare $winSys aqua] == 0} {
		update			;# because of a Tk bug on Mac OS X Aqua
	    }
	    return ""
	}

	default {
	    #
	    # Command: $win yview moveto fraction
	    #	       $win yview scroll number what
	    #
	    set opt [lindex $argList 0]
	    if {[string first $opt moveto] == 0} {
		eval [list $w yview] $argList
	    } elseif {[string first $opt scroll] == 0} {
		if {$argCount != 3} {
		    #
		    # Let Tk report the error
		    #
		    return [eval [list $w yview] $argList]
		}

		set number [lindex $argList 1]
		set number [format %d $number]	;# integer check with error msg
		set what [lindex $argList 2]
		if {[string first $what units] == 0} {
		    $w yview scroll $number units
		} elseif {[string first $what pages] == 0} {
		    if {$number < 0} {
			$w yview scroll $number pages
		    } else {
			#
			# The following loop is needed because "$w yview scroll
			# $number pages" doesn't produce the expected effect.
			#
			for {set n 0} {$n < $number} {incr n} {
			    $w yview scroll 1 pages
			    if {[lindex [$w yview] 1] < 1.0} {
				$w yview scroll -1 units
			    }
			}
		    }
		} else {
		    #
		    # Let Tk report the error
		    #
		    return [eval [list $w yview] $argList]
		}
	    } else {
		return -code error \
		       "unknown option \"$opt\": must be moveto or scroll"
	    }
	    if {[string compare $winSys aqua] == 0} {
		update			;# because of a Tk bug on Mac OS X Aqua
	    }
	    return ""
	}
    }
}

#
# Private callback procedures
# ===========================
#

#------------------------------------------------------------------------------
# tablelist::fetchSelection
#
# This procedure is invoked when the PRIMARY selection is owned by the
# tablelist widget win and someone attempts to retrieve it as a STRING.  It
# returns part or all of the selection, as given by offset and maxChars.  The
# string which is to be (partially) returned is built by joining all of the
# visible elements of the selected rows together with tabs and the rows
# themselves with newlines.
#------------------------------------------------------------------------------
proc tablelist::fetchSelection {win offset maxChars} {
    upvar ::tablelist::ns${win}::data data

    if {!$data(-exportselection)} {
	return ""
    }

    set selection ""
    set gotItem 0
    foreach idx [curselectionSubCmd $win] {
	if {$gotItem} {
	    append selection \n
	}

	set item [lindex $data(itemList) $idx]
	set gotText 0
	for {set col 0} {$col < $data(colCount)} {incr col} {
	    if {$data($col-hide)} {
		continue
	    }

	    set text [lindex $item $col]
	    if {[info exists data($col-formatcommand)]} {
		set text [uplevel #0 $data($col-formatcommand) [list $text]]
	    }

	    if {$gotText} {
		append selection \t
	    }
	    append selection $text

	    set gotText 1
	}

	set gotItem 1
    }

    return [string range $selection $offset [expr {$offset + $maxChars - 1}]]
}

#------------------------------------------------------------------------------
# tablelist::lostSelection
#
# This procedure is invoked when the tablelist widget win loses ownership of
# the PRIMARY selection.  It deselects all items of the widget with the aid of
# the selectionSubCmd procedure if the selection is exported.
#------------------------------------------------------------------------------
proc tablelist::lostSelection win {
    upvar ::tablelist::ns${win}::data data

    if {$data(-exportselection)} {
	selectionSubCmd $win clear 0 $data(lastRow)
    }
}

#------------------------------------------------------------------------------
# tablelist::activeIdxTrace
#
# This procedure is executed whenever the array element data(activeIdx) is
# written.  It moves the "active" tag to the line that displays the active
# element of the widget in its body text child if the latter has the focus.
#------------------------------------------------------------------------------
proc tablelist::activeIdxTrace {win varName index op} {
    upvar ::tablelist::ns${win}::data data

    set w $data(body)
    if {$data(ownsFocus)} {
	set line [expr {$data(oldActiveIdx) + 1}]
	$w tag remove active $line.0 $line.end

	set line [expr {$data(activeIdx) + 1}]
	$w tag add active $line.0 $line.end
    }

    set data(oldActiveIdx) $data(activeIdx)
}

#------------------------------------------------------------------------------
# tablelist::listVarTrace
#
# This procedure is executed whenever the global variable specified by varName
# is written or unset.  It makes sure that the contents of the widget will be
# synchronized with the value of the variable at idle time, and that the
# variable is recreated if it was unset.
#------------------------------------------------------------------------------
proc tablelist::listVarTrace {win varName index op} {
    upvar ::tablelist::ns${win}::data data

    switch $op {
	w {
	    if {![info exists data(syncId)]} {
		#
		# Arrange for the contents of the widget to be synchronized
		# with the value of the variable ::$varName at idle time
		#
		set data(syncId) [after idle [list tablelist::synchronize $win]]

		#
		# Cancel the execution of all delayed redisplay
		# commands, to make sure that the synchronize command
		# will be invoked first; the latter will then schedule
		# a redisplay command for execution at idle time
		#
		if {[info exists data(redispId)]} {
		    after cancel $data(redispId)
		}
	    }
	}

	u {
	    #
	    # Recreate the variable ::$varName by setting it according to
	    # the value of data(itemList), and set the trace on it again
	    #
	    if {[string compare $index ""] != 0} {
		set varName ${varName}($index)
	    }
	    set ::$varName {}
	    foreach item $data(itemList) {
		lappend ::$varName [lrange $item 0 $data(lastCol)]
	    }
	    trace variable ::$varName wu $data(listVarTraceCmd)
	}
    }
}

#
# Private procedures used in bindings
# ===================================
#
# See the module "tablelistBind.tcl".
#

#
# Private utility procedures
# ==========================
#
# See the module "tablelistUtil.tcl"
#
