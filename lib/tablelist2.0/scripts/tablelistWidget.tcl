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
#   - Private helper procedures
#
# Copyright (c) 2000-2001  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

#
# Namespace initialization
# ========================
#

namespace eval tablelist {
    #
    # Some help variables needed below:
    #
    variable helpListbox	.__helpListbox
    variable helpLabel		.__helpLabel
    variable helpButton		.__helpButton
    variable n
    variable configSet
    variable opt
    variable optTail

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
	-arrowcolor		 {arrowColor		  ArrowColor	      w}
	-arrowdisabledcolor	 {arrowDisabledColor	  ArrowDisabledColor  w}
	-background		 {background		  Background	      b}
	-bg			 -background
	-borderwidth		 {borderWidth		  BorderWidth	      f}
	-bd			 -borderwidth
	-columns		 {columns		  Columns	      w}
	-cursor			 {cursor		  Cursor	      c}
	-disabledforeground	 {disabledForeground	  DisabledForeground  w}
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
	-relief			 {relief		  Relief	      f}
	-resizablecolumns	 {resizableColumns	  ResizableColumns    w}
	-resizecursor		 {resizeCursor		  ResizeCursor	      w}
	-selectbackground	 {selectBackground	  Foreground	      w}
	-selectborderwidth	 {selectBorderWidth	  BorderWidth	      w}
	-selectforeground	 {selectForeground	  Background	      w}
	-selectmode		 {selectMode		  SelectMode	      w}
	-showarrow		 {showArrow		  ShowArrow	      w}
	-showlabels		 {showLabels		  ShowLabels	      w}
	-setgrid		 {setGrid		  SetGrid	      w}
	-sortcommand		 {sortCommand		  SortCommand	      w}
	-state			 {state			  State		      w}
	-stretch		 {stretch		  Stretch	      w}
	-takefocus		 {takeFocus		  TakeFocus	      f}
	-width			 {width			  Width		      w}
	-xscrollcommand		 {xScrollCommand	  ScrollCommand	      h}
	-yscrollcommand		 {yScrollCommand	  ScrollCommand	      b}
    }

    #
    # Append the default values of the configuration options
    # of a temporary, invisible listbox widget to the values
    # of the corresponding elements of the array configSpecs
    #
    for {set n 0} {[winfo exists $helpListbox]} {incr n} {
	set helpListbox .__helpListbox$n
    }
    listbox $helpListbox
    foreach configSet [$helpListbox config] {
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
    for {set n 0} {[winfo exists $helpLabel]} {incr n} {
	set helpLabel .__helpLabel$n
    }
    label $helpLabel
    foreach optTail {background font foreground height pady} {
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

    #
    # Steal the default values of some configuration
    # options from a temporary, invisible button widget
    #
    for {set n 0} {[winfo exists $helpButton]} {incr n} {
	set helpButton .__helpButton$n
    }
    button $helpButton
    foreach opt {-disabledforeground -state} {
	set configSet [$helpButton config $opt]
	lappend configSpecs($opt) [lindex $configSet 3]
    }
    set configSet [$helpButton config -borderwidth]
    lappend configSpecs(-labelborderwidth) [lindex $configSet 3]
    destroy $helpButton

    #
    # Extend the remaining elements of the array configSpecs
    #
    lappend configSpecs(-arrowcolor)		gray45
    lappend configSpecs(-arrowdisabledcolor)	gray65
    lappend configSpecs(-columns)		{0 {} left}
    lappend configSpecs(-incrarrowtype)		up
    lappend configSpecs(-labelcommand)		{}
    lappend configSpecs(-labelrelief)		raised
    lappend configSpecs(-listvariable)		{}
    lappend configSpecs(-resizablecolumns)	1
    lappend configSpecs(-resizecursor)		sb_h_double_arrow
    lappend configSpecs(-showarrow)		1
    lappend configSpecs(-showlabels)		1
    lappend configSpecs(-sortcommand)		{}
    lappend configSpecs(-stretch)		{}

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
	-foreground		{foreground		Foreground	}
	-fg			-foreground
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
	-resizable		{resizable		Resizable	}
	-selectbackground	{selectBackground	Foreground	}
	-selectforeground	{selectForeground	Background	}
	-showarrow		{showArrow		ShowArrow	}
	-sortcommand		{sortCommand		SortCommand	}
	-sortmode		{sortMode		SortMode	}
	-title			{title			Title		}
	-width			{width			Width		}
    }

    #
    # Extend some elements of the array colConfigSpecs
    #
    lappend colConfigSpecs(-align)	- left
    lappend colConfigSpecs(-hide)	- 0
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
	-foreground		{foreground		Foreground	}
	-fg			-foreground
	-selectbackground	{selectBackground	Foreground	}
	-selectforeground	{selectForeground	Background	}
	-text			{text			Text		}
    }

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
	-foreground		{foreground		Foreground	}
	-fg			-foreground
	-image			{image			Image		}
	-selectbackground	{selectBackground	Foreground	}
	-selectforeground	{selectForeground	Background	}
	-text			{text			Text		}
    }

    #
    # Use a list to facilitate the handling of the command options 
    #
    variable cmdOpts [list \
	activate attrib bbox bodypath cellcget cellconfigure cellindex cget \
	columncget columnconfigure columncount columnindex configure \
	curselection delete get index insert labelpath labels nearest \
	nearestcell nearestcolumn resetsortinfo rowcget rowconfigure scan see \
	selection size sort sortbycolumn sortcolumn sortorder xview yview]

    #
    # Use lists to facilitate the handling of miscellaneous options
    #
    variable alignments		[list left right center]
    variable arrowTypes		[list up down]
    variable states		[list disabled normal]
    variable sortModes		[list ascii command dictionary integer real]
    variable sortOrders		[list -increasing -decreasing]
    variable scanCmdOpts	[list mark dragto]
    variable selCmdOpts		[list anchor clear includes set]

    #
    # Define some Tablelist class bindings
    #
    bind Tablelist <FocusIn> {
	tablelist::addActiveTag %W

	if {[string compare [focus -lastfor %W] %W] == 0} {
	    focus [%W bodypath]
	}
    }
    bind Tablelist <FocusOut> {
	tablelist::removeActiveTag %W
    }
    bind Tablelist <Destroy> {
	tablelist::cleanup %W
    }

    #
    # Define the binding tag TablelistBody to have the same events as Listbox
    # and the binding scripts obtained from those of Listbox by replacing the
    # widget %W with [winfo parent %W] as well as the %x and %y fields with
    # [expr {%x + [winfo x %W]}] and [expr {%y + [winfo y %W]}], respectively
    # 
    foreach event [bind Listbox] {
	set script [bind Listbox $event]
	regsub -all %W $script {$tablelist::win} script
	regsub -all %x $script {$tablelist::x} script
	regsub -all %y $script {$tablelist::y} script
	regsub -all tkListboxAutoScan $script \
		    tablelist::tablelistAutoScan script

	bind TablelistBody $event [format {
	    set tablelist::win [winfo parent %%W]
	    set tablelist::x [expr {%%x + [winfo x %%W]}]
	    set tablelist::y [expr {%%y + [winfo y %%W]}]
	    %s
	} $script]
    }

    #
    # Define the virtual event <<Button3>>
    #
    switch $::tcl_platform(platform) {
	unix -
	windows {
	    event add <<Button3>> <Button-3>
	}
	macintosh {
	    event add <<Button3>> <Control-Button-1>
	}
    }
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
	    frame $win -class Tablelist -colormap . -container 0 \
		       -height 0 -width 0
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
	    hdrPixels		 0
	    activeIdx		 0
	    anchorIdx		 0
	    seqNum		-1
	    itemList		 {}
	    itemCount		 0
	    colList		 {0 left}
	    colCount		 1
	    clickedLblCol	-1
	    arrowCol		-1
	    sortCol		-1
	    sortOrder		 {}
	    forceAdjust		 0
	    0-hide		 0
	    0-delta		 0
	}

	#
	# The following array holds the configuration
	# options and their values for this widget
	#
	variable configVals

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
    set data(body)		$win.body
    set data(hdr)		$win.hdr
    set data(hdrTxt)		$data(hdr).t
    set data(hdrTxtFr)		$data(hdrTxt).f
    set data(hdrTxtFrCanv)	$data(hdrTxtFr).c
    set data(hdrTxtFrLbl)	$data(hdrTxtFr).l
    set data(lb)		$win.lb

    #
    # Initialize the components of configVals
    #
    upvar ::tablelist::ns${win}::configVals configVals
    foreach opt $configOpts {
	set configVals($opt) [lindex $configSpecs($opt) 3]
    }

    #
    # Create a child hierarchy used to hold the column labels.  The
    # labels will be created as children of the frame data(hdrTxtFr),
    # which is embedded into the text widget data(hdrTxt) (in order
    # to make it scrollable), which in turn fills the frame data(hdr)
    # (whose width and height can be set arbitrarily in pixels).
    #
    frame $data(hdr) -borderwidth 0 -colormap . -container 0 -height 0 \
		     -highlightthickness 0 -relief flat -takefocus 0 \
		     -width 0
    bind $data(hdr) <Configure> [list tablelist::stretchColumnsWhenIdle $win]
    pack $data(hdr) -fill x
    text $data(hdrTxt) -borderwidth 0 -highlightthickness 0 -insertwidth 0 \
		       -padx 0 -pady 0 -state normal -takefocus 0 -wrap none
    frame $data(hdrTxtFr) -borderwidth 0 -colormap . -container 0 -height 0 \
			  -highlightthickness 0 -relief flat -takefocus 0 \
			  -width 0
    $data(hdrTxt) window create 1.0 -window $data(hdrTxtFr)
    label $data(hdrTxtFrLbl)0 

    #
    # Create a canvas as a child of the frame data(hdrTxtFr),
    # needed for displaying an up- or down-arrow when
    # sorting the items by a column.   Set its width and
    # height to temporary values and draw two up-arrows
    #
    set w $data(hdrTxtFrCanv)
    set size 5
    canvas $w -borderwidth 0 -height $size -highlightthickness 0 \
	      -relief flat -takefocus 0 -width $size
    set last [expr {$size - 1}]
    set half [expr {$size / 2}]
    $w create polygon 0 $last $half 0 $last $last -tags normalArrow
    $w create polygon 0 $last $half 0 $last $last -tags disabledArrow

    #
    # Remove the binding tag Text from the list
    # of binding tags of the header text widget
    #
    set w $data(hdrTxt)
    bindtags $w [list $w [winfo toplevel $w] all]

    #
    # Create the body text widget within the main frame
    #
    set w $data(body)
    text $w -borderwidth 0 -exportselection no -highlightthickness 0 \
	    -insertwidth 0 -padx 0 -pady 0 -state normal -takefocus 0 \
	    -takefocus tablelist::focusCtrl -wrap none
    pack $w -expand yes -fill both

    #
    # Replace the binding tags Text and all with TablelistBody and allModif,
    # respectively, in the list of binding tags of the body text widget
    # (see mwutil.tcl for the definition of the bindings for the tag allModif)
    #
    bindtags $w [list $w TablelistBody [winfo toplevel $w] allModif]

    #
    # Create the "active", "select", and "disabled" tags in
    # the body text widget.  Don't use the built-in "sel" tag
    # because on Windows the selection in a text widget only
    # becomes visible when the window gets the input focus.
    #
    $w tag configure active -underline yes
    $w tag configure select -relief raised
    $w tag configure disabled -underline no

    #
    # Create an unmanaged listbox child, used to handle the -setgrid option
    #
    listbox $data(lb)

    #
    # Configure the widget according to the command-line
    # arguments and to the available database options
    #
    if {[catch {
	    mwutil::configure $win configSpecs configVals tablelist::doConfig \
			      [lrange $args 1 end] yes
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
    # owned by the window win and someone attempts to retrieve it as a STRING
    #
    selection handle $win [list ::tablelist::fetchSelection $win]

    return $win
}

#
# Private configuration procedures
# ================================
#

#------------------------------------------------------------------------------
# tablelist::doConfig
#
# Applies the value val of the configuration option opt to the tablelist widget
# win.
#------------------------------------------------------------------------------
proc tablelist::doConfig {win opt val} {
    variable helpLabel
    variable configSpecs
    variable arrowTypes
    variable states
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    #
    # Apply the value to the widget(s) corresponding to the given option
    #
    switch [lindex $configSpecs($opt) 2] {
	c {
	    #
	    # Apply the value to all children and save the
	    # properly formatted value of val in configVals($opt)
	    #
	    foreach w [winfo children $win] {
		$w configure $opt $val
	    }
	    $data(hdrTxt) configure $opt $val
	    foreach w [winfo children $data(hdrTxtFr)] {
		$w configure $opt $val
		foreach c [winfo children $w] {
		    $c configure $opt $val
		}
	    }
	    $helpLabel configure $opt $val
	    set configVals($opt) [$helpLabel cget $opt]
	}

	b {
	    #
	    # Apply the value to the body text widget and save the
	    # properly formatted value of val in configVals($opt)
	    #
	    set w $data(body)
	    $w configure $opt $val
	    set configVals($opt) [$w cget $opt]

	    switch -- $opt {
		-background {
		    #
		    # Apply the value to the header frame (this will
		    # make it fully invisible if the labels are not
		    # being shown) and to the "disabled" tag
		    #
		    $data(hdr) configure $opt $val
		    $w tag configure disabled $opt $val
		}
		-font {
		    #
		    # Apply the value to the listbox child, adjust the columns,
		    # and make sure the items will be redisplayed at idle time
		    #
		    $data(lb) configure $opt $val
		    adjustColumns $win all yes
		    redisplayWhenIdle $win
		}
		-foreground {
		    #
		    # Apply the value to the "disabled" tag if needed
		    #
		    if {[string compare $configVals(-disabledforeground) ""]
			== 0} {
			$w tag configure disabled $opt $val
		    }
		}
	    }
	}

	h {
	    #
	    # Apply the value to the header text widget and save the
	    # properly formatted value of val in configVals($opt)
	    #
	    set w $data(hdrTxt)
	    $w configure $opt $val
	    set configVals($opt) [$w cget $opt]
	}

	l {
	    #
	    # Apply the value to all not individually configured labels and
	    # save the properly formatted value of val in configVals($opt)
	    #
	    set optTail [string range $opt 6 end]	;# remove the -label
	    for {set col 0} {$col < $data(colCount)} {incr col} {
		set w $data(hdrTxtFrLbl)$col
		if {![info exists data($col$opt)]} {
		    $w configure -$optTail $val
		}
	    }
	    $helpLabel configure -$optTail $val
	    set configVals($opt) [$helpLabel cget -$optTail]

	    switch -- $opt {
		-labelbackground {
		    #
		    # Apply the value to the text child of the header frame
		    # and conditionally both to the children of the labels (if
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
			$data(hdrTxtFrCanv) configure -$optTail $val
		    }
		}
		-labelborderwidth {
		    #
		    # Adjust the columns (including
		    # the height of the header frame)
		    #
		    adjustColumns $win all yes
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
			redrawArrows $win
		    }
		    adjustColumns $win all yes
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
	    # Apply the value to the frame and save the properly
	    # formatted value of val in configVals($opt)
	    #
	    $win configure $opt $val
	    set configVals($opt) [$win cget $opt]
	}

	w {
	    switch -- $opt {
		-arrowcolor {
		    #
		    # Set the color of the normal arrow and save the
		    # properly formatted value of val in configVals($opt)
		    #
		    set w $data(hdrTxtFrCanv)
		    $w itemconfigure normalArrow -fill $val -outline $val
		    set configVals($opt) [$w itemcget normalArrow -fill]
		}
		-arrowdisabledcolor {
		    #
		    # Set the color of the disabled arrow and save the
		    # properly formatted value of val in configVals($opt)
		    #
		    set w $data(hdrTxtFrCanv)
		    $w itemconfigure disabledArrow -fill $val -outline $val
		    set configVals($opt) [$w itemcget disabledArrow -fill]
		}
		-columns {
		    #
		    # Set up and adjust the columns, and make sure
		    # the items will be redisplayed at idle time
		    #
		    setupColumns $win $val yes
		    adjustColumns $win all yes
		    redisplayWhenIdle $win
		}
		-disabledforeground {
		    #
		    # Configure the "disabled" tag in the body
		    # text widget and save the properly formatted
		    # value of val in configVals($opt)
		    #
		    set w $data(body)
		    if {[string compare $val ""] == 0} {
			$w tag configure disabled -fgstipple gray50 \
			       -foreground $configVals(-foreground)
			set configVals($opt) ""
		    } else {
			$w tag configure disabled -fgstipple "" \
			       -foreground $val
			set configVals($opt) [$w tag cget disabled -foreground]
		    }
		}
		-exportselection {
		    #
		    # Save the boolean value specified by val in
		    # configVals($opt).  In addition, if the selection
		    # is exported and there are any selected rows in the
		    # widget then make win the new owner of the PRIMARY
		    # selection and register a callback to be invoked
		    # when it loses ownership of the PRIMARY selection
		    #
		    set configVals($opt) [expr {$val ? 1 : 0}]
		    if {$val &&
			[llength [$data(body) tag nextrange select 1.0]] != 0} {
			selection own -command \
				  [list ::tablelist::lostSelection $win] $win
		    }
		}
		-height {
		    #
		    # Adjust the height of the body text widget and save the
		    # properly formatted value of val in configVals($opt)
		    #
		    set val [format %d $val]	;# integer check with error msg
		    if {$val <= 0} {
			$data(body) configure $opt $data(itemCount)
		    } else {
			$data(body) configure $opt $val
		    }
		    set configVals($opt) $val
		}
		-incrarrowtype {
		    #
		    # Save the properly formatted value of val
		    # in configVals($opt) and redraw the arrows
		    # if the canvas widget is presently mapped
		    #
		    set configVals($opt) [mwutil::fullOpt "arrow type" \
					  $val $arrowTypes]
		    if {$data(arrowCol) >= 0} {
			redrawArrows $win
		    }
		}
		-labelcommand -
		-selectmode -
		-sortcommand {
		    #
		    # Save val in configVals($opt)
		    #
		    set configVals($opt) $val
		}
		-listvariable {
		    #
		    # Associate val as list variable with the
		    # given widget and save it in configVals($opt)
		    #
		    makeListVar $win $val
		    set configVals($opt) $val
		}
		-resizablecolumns {
		    #
		    # Save the boolean value specified
		    # by val in configVals($opt)
		    #
		    set configVals($opt) [expr {$val ? 1 : 0}]
		}
		-resizecursor {
		    #
		    # Save the properly formatted value
		    # of val in configVals($opt)
		    #
		    $helpLabel configure -cursor $val
		    set configVals($opt) [$helpLabel cget -cursor]
		}
		-selectbackground -
		-selectforeground {
		    #
		    # Configure the "select" tag in the body text widget and
		    # save the properly formatted value of val in
		    # configVals($opt).  Don't use the built-in "sel" tag
		    # because on Windows the selection in a text widget only
		    # becomes visible when the window gets the input focus.
		    #
		    set w $data(body)
		    set optTail [string range $opt 7 end] ;# remove the -select
		    $w tag configure select -$optTail $val
		    set configVals($opt) [$w tag cget select -$optTail]
		}
		-selectborderwidth {
		    #
		    # Configure the "select" tag in the body text widget
		    # and save the properly formatted value of val in
		    # configVals($opt).  Don't use the built-in "sel" tag
		    # because on Windows the selection in a text widget only
		    # becomes visible when the window gets the input focus.
		    # In addition, adjust the line spacing accordingly and
		    # apply the value to the listbox child, too.
		    #
		    set w $data(body)
		    set optTail [string range $opt 7 end] ;# remove the -select
		    $w tag configure select -$optTail $val
		    set configVals($opt) [$w tag cget select -$optTail]
		    if {$val < 0} {
			set val 0
		    }
		    $w configure -spacing1 $val -spacing3 [expr {$val + 1}]
		    $data(lb) configure $opt $val
		}
		-setgrid {
		    #
		    # Apply the value to the listbox child and save the
		    # properly formatted value of val in configVals($opt)
		    #
		    $data(lb) configure $opt $val
		    set configVals($opt) [$data(lb) cget $opt]
		}
		-showarrow {
		    #
		    # Save the boolean value specified by val in
		    # configVals($opt) and conditionally unmanage
		    # the canvas displaying an up- or down-arrow
		    #
		    set configVals($opt) [expr {$val ? 1 : 0}]
		    if {!$configVals($opt)} {
			place forget $data(hdrTxtFrCanv)
			set colOfUnknownWidth $data(arrowCol)
			set data(arrowCol) -1
			adjustColumns $win $colOfUnknownWidth yes
		    }
		}
		-showlabels {
		    #
		    # Save the boolean value specified by
		    # val in configVals($opt) and adjust
		    # the height of the header frame
		    #
		    set configVals($opt) [expr {$val ? 1 : 0}]
		    adjustHeaderHeight $win
		}
		-state {
		    #
		    # Apply the value to all labels and their children
		    # (if any), raise the corresponding arrow in the
		    # canvas, add/remove the "disabled" tag to/from the
		    # contents of the body text widget, and save the
		    # properly formatted value of val in configVals($opt)
		    #
		    set val [mwutil::fullOpt "state" $val $states]
		    catch {
			foreach w [winfo children $data(hdrTxtFr)] {
			    $w configure $opt $val
			    foreach c [winfo children $w] {
				$c configure $opt $val
			    }
			}
		    }
		    $data(hdrTxtFrCanv) raise ${val}Arrow
		    set w $data(body)
		    switch $val {
			disabled {
			    $w tag add disabled 1.0 end
			    $w tag configure select -borderwidth 0
			}
			normal {
			    $w tag remove disabled 1.0 end
			    $w tag configure select -borderwidth \
				   $configVals(-selectborderwidth)
			}
		    }
		    set configVals($opt) $val
		}
		-stretch {
		    #
		    # Save val in configVals($opt) and
		    # stretch the stretchable columns
		    #
		    if {[string first $val all] == 0} {
			set configVals($opt) all
		    } else {
			foreach col $val {
			    colIndex $win $col
			}
			set configVals($opt) $val
		    }
		    set data(forceAdjust) 1
		    stretchColumnsWhenIdle $win
		}
		-width {
		    #
		    # Adjust the widths of the body text widget and
		    # of the header frame, and save the properly
		    # formatted value of val in configVals($opt)
		    #
		    set val [format %d $val]	;# integer check with error msg
		    $data(body) configure $opt $val
		    if {$val <= 0} {
			$data(hdr) configure $opt $data(hdrPixels)
		    } else {
			$data(hdr) configure $opt 0
		    }
		    set configVals($opt) $val
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
    variable configSpecs
    variable configOpts
    variable alignments
    variable sortModes
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    switch -- $opt {
	-align {
	    #
	    # Set up and adjust the columns, and make sure the
	    # given column will be redisplayed at idle time
	    #
	    set idx [expr {3*$col + 2}]
	    setupColumns $win [lreplace $configVals(-columns) $idx $idx $val] no
	    adjustColumns $win -1 yes
	    redisplayColWhenIdle $win $col
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
		$w tag lower $tag

		if {!$data($col-hide)} {
		    #
		    # Apply the tag to the elements of the given column
		    #
		    for {set idx 0; set line 1} {$idx < $data(itemCount)} \
			{incr idx; incr line} {
			if {[lsearch -exact [$w tag names $line.0] select]
			    < 0} {
			    findCellTabs $win $line.0 $col tabIdx1 tabIdx2
			    $w tag add $tag $tabIdx1 $tabIdx2+1c
			}
		    }
		}

		#
		# Save the properly formatted value of val in data($tag)
		#
		set data($tag) [$w tag cget $tag $opt]
	    }
	}

	-hide {
	    #
	    # Save the boolean value specified by val in
	    # data($col$opt), adjust the columns, and make
	    # sure the items will be redisplayed at idle time
	    #
	    set data($col$opt) [expr {$val ? 1 : 0}]
	    adjustColumns $win $col yes
	    redisplayWhenIdle $win
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
		set val [mwutil::fullOpt "label alignment" $val $alignments]
		set alignment $val
		set data($col$opt) $val
	    }

	    #
	    # Adjust the col`th label
	    #
	    set charWidth [font measure $configVals(-font) -displayof $win 0]
	    set pixels [lindex $data(colList) [expr {2*$col}]]
	    if {$pixels != 0} {			;# convention: static width
		incr pixels $data($col-delta)
	    }
	    adjustLabel $win $col $charWidth $pixels $alignment
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
		$w configure -$optTail $configVals($opt)
		foreach c [winfo children $w] {
		    $c configure -$optTail $configVals($opt)
		}
		if {$col == $data(arrowCol)} {
		    $data(hdrTxtFrCanv) configure -$optTail $configVals($opt)
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
		    $data(hdrTxtFrCanv) configure -$optTail $val
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
		$w configure -$optTail $configVals($opt)
		if {[info exists data($col$opt)]} {
		    unset data($col$opt)
		}
	    } else {
		#
		# Apply the given value to the col'th label and save the
		# properly formatted value of val in data($col$opt)
		#
		$w configure -$optTail $val
		set data($col$opt) [$w cget -$optTail]
	    }

	    #
	    # Adjust the columns (including the height of the header frame)
	    #
	    adjustColumns $win $col yes
	}

	-labelcommand {
	    #
	    # Save val in data($col$opt)
	    #
	    set data($col$opt) $val
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
		$w configure -$optTail $configVals($opt)
		foreach c [winfo children $w] {
		    $c configure -$optTail $configVals($opt)
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
		redrawArrows $win
	    }
	    adjustColumns $win $col yes
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
		$w configure -$optTail $configVals($opt)
		foreach c [winfo children $w] {
		    $c configure -$optTail $configVals($opt)
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

	-labelrelief {
	    set w $data(hdrTxtFrLbl)$col
	    set optTail [string range $opt 6 end]	;# remove the -label
	    if {[string compare $val ""] == 0} {
		#
		# Apply the value of the corresponding widget configuration
		# option to the col'th label and unset data($col$opt)
		#
		$w configure -$optTail $configVals($opt)
		if {[info exists data($col$opt)]} {
		    unset data($col$opt)
		}
	    } else {
		#
		# Apply the given value to the col'th label and save the
		# properly formatted value of val in data($col$opt)
		#
		$w configure -$optTail $val
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
		$w configure -$optTail $configVals($opt)
		if {[info exists data($col$opt)]} {
		    unset data($col$opt)
		}
	    } else {
		#
		# Apply the given value to the col'th label and save the
		# properly formatted value of val in data($col$opt)
		#
		$w configure -$optTail $val
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
		    foreach c [list $w.il $w.tl] {	;# image and text labels
			#
			# Create the label $c
			#
			label $c -borderwidth 0 -height 0 \
				 -highlightthickness 0 \
				 -padx 0 -pady 0 -takefocus 0 -width 0

			#
			# Apply to it the current configuration options
			#
			foreach opt $configOpts {
			    if {[llength $configSpecs($opt)] != 1 &&
				[string compare [lindex $configSpecs($opt) 2] c]
				== 0} {
				$c configure $opt $configVals($opt)
			    }
			}
			foreach opt {-background -font -foreground} {
			    $c configure $opt [$w cget $opt]
			}
			foreach opt {-disabledforeground -state} {
			    catch {$c configure $opt [$w cget $opt]}
			}
			set opt -labelimage	;# restore the original value

			#
			# Define for the binding tag $c the binding scripts
			# obtained from those of $w by replacing the %x
			# field with [expr {%x + [winfo x $c]}] and the
			# %y field with [expr {%y + [winfo y $c]}]
			#
			foreach event [bind $w] {
			    set script [bind $w $event]
			    regsub -all %x $script {$tablelist::x} script
			    regsub -all %y $script {$tablelist::y} script

			    bind $c $event [format {
				set tablelist::x [expr {%%x + [winfo x %s]}]
				set tablelist::y [expr {%%y + [winfo y %s]}]
				%s
			    } $c $c $script]
			}
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
	    adjustColumns $win $col yes
	}

	-resizable {
	    #
	    # Save the boolean value specified by val in data($col$opt)
	    #
	    set data($col$opt) [expr {$val ? 1 : 0}]
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

			findCellTabs $win $selStart $col tabIdx1 tabIdx2
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
		adjustColumns $win $col yes
	    }
	}

	-sortcommand {
	    #
	    # Save val in data($col$opt)
	    #
	    set data($col$opt) $val
	}

	-sortmode {
	    #
	    # Save the properly formatted value of val in data($col$opt)
	    #
	    set data($col$opt) [mwutil::fullOpt "sort mode" $val $sortModes]
	}

	-title {
	    #
	    # Save the given value in the corresponding element
	    # of configVals(-columns) and adjust the columns
	    #
	    set idx [expr {3*$col + 1}]
	    set configVals(-columns) \
		[lreplace $configVals(-columns) $idx $idx $val]
	    adjustColumns $win $col yes
	}

	-width {
	    #
	    # Set up and adjust the columns, and make sure the
	    # given column will be redisplayed at idle time
	    #
	    set idx [expr {3*$col}]
	    setupColumns $win [lreplace $configVals(-columns) $idx $idx $val] no
	    adjustColumns $win $col yes
	    redisplayColWhenIdle $win $col
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
    upvar ::tablelist::ns${win}::configVals configVals

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
		set data($tag) [$w tag cget $tag $opt]
	    }
	}

	-text {
	    if {[string compare $configVals(-state) disabled] == 0} {
		return ""
	    }

	    set font $configVals(-font)
	    set colWidthsChanged no
	    set oldItem [lindex $data(itemList) $row]
	    set key [lindex $oldItem end]
	    set newItem {}
	    set line [expr {$row + 1}]

	    set textIdx1 $line.1
	    set col 0
	    foreach {pixels alignment} $data(colList) {
		set text [lindex $val $col]
		regsub -all "\t|\n" $text " " text
		lappend newItem $text

		if {$data($col-hide)} {
		    incr col
		    continue
		}

		#
		# Adjust the cell text and the image width
		#
		set name $key-$col-image
		if {[info exists data($name)]} {
		    set image $data($name)
		    set imageWidth [image width $image]
		} else {
		    set image ""
		    set imageWidth 0
		}
		set oldImageWidth $imageWidth
		if {$pixels != 0} {		;# convention: static width
		    incr pixels $data($col-delta)
		}
		adjustElem $win text imageWidth $font $pixels $alignment

		#
		# Delete the old cell contents between the
		# two tabs, and insert the text and the image
		#
		set textIdx2 [$w search \t $textIdx1 $line.end]
		$w delete $textIdx1 $textIdx2
		insertElem $w $textIdx1 $text $image $imageWidth $alignment

		if {$pixels == 0} {		;# convention: dynamic width
		    #
		    # Check whether the width of the current column has changed
		    #
		    set textWidth [font measure $font -displayof $win $text]
		    set newElemWidth [expr {$imageWidth + $textWidth}]
		    if {$newElemWidth > $data($col-width)} {
			set colWidthsChanged yes
		    } else {
			set oldText [lindex $oldItem $col]
			adjustElem $win oldText oldImageWidth \
				   $font $pixels $alignment
			set oldTextWidth [font measure $font \
					  -displayof $win $oldText]
			set oldElemWidth [expr {$oldImageWidth + $oldTextWidth}]
			if {$oldElemWidth == $data($col-width) &&
			    $newElemWidth < $oldElemWidth} {
			    set colWidthsChanged yes
			}
		    }
		}

		set textIdx1 [$w search \t $textIdx1 $line.end]+2c
		incr col
	    }

	    #
	    # Replace the row contents in the list variable if present
	    #
	    if {[string compare $configVals(-listvariable) ""] != 0} {
		set traceCmd [list tablelist::traceProc $win]
		trace vdelete ::$configVals(-listvariable) wu $traceCmd
		upvar #0 $configVals(-listvariable) var
		set var [lreplace $var $row $row $newItem]
		trace variable ::$configVals(-listvariable) wu $traceCmd
	    }

	    #
	    # Replace the row contents in the internal list
	    #
	    lappend newItem [lindex $oldItem $col]
	    set data(itemList) [lreplace $data(itemList) $row $row $newItem]

	    #
	    # Adjust the columns if necessary
	    #
	    if {$colWidthsChanged} {
		adjustColumns $win all yes
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
    upvar ::tablelist::ns${win}::configVals configVals

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
			findCellTabs $win $line.0 $col tabIdx1 tabIdx2
			$w tag add $tag $tabIdx1 $tabIdx2+1c
		    }
		}

		#
		# Save the properly formatted value of val in data($tag)
		#
		set data($tag) [$w tag cget $tag $opt]
	    }
	}

	-image {
	    if {[string compare $configVals(-state) disabled] == 0} {
		return ""
	    }

	    set font $configVals(-font)
	    set pixels [lindex $data(colList) [expr {2*$col}]]
	    if {$pixels != 0} {			;# convention: static width
		incr pixels $data($col-delta)
	    }
	    set alignment [lindex $data(colList) [expr {2*$col + 1}]]

	    #
	    # Adjust the cell text and the image width
	    #
	    set oldItem [lindex $data(itemList) $row]
	    set text [lindex $oldItem $col]
	    set image $val
	    if {[string compare $image ""] == 0} {
		set imageWidth 0
	    } else {
		set imageWidth [image width $image]
	    }
	    adjustElem $win text imageWidth $font $pixels $alignment

	    if {!$data($col-hide)} {
		#
		# Delete the old cell contents between the
		# two tabs, and insert the text and the image
		#
		findCellTabs $win [expr {$row + 1}].0 $col tabIdx1 tabIdx2
		$w delete $tabIdx1+1c $tabIdx2
		insertElem $w $tabIdx1+1c $text $image $imageWidth $alignment
	    }

	    #
	    # Save the old image width
	    #
	    set key [lindex $oldItem end]
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
		}
	    } else {
		set data($name) $val
	    }

	    #
	    # Adjust the columns if necessary
	    #
	    if {$pixels == 0} {			;# convention: dynamic width
		set textWidth [font measure $font -displayof $win $text]
		set newElemWidth [expr {$imageWidth + $textWidth}]
		if {$newElemWidth > $data($col-width)} {
		    set data($col-width) $newElemWidth
		    adjustColumns $win -1 yes
		} else {
		    set oldText [lindex $oldItem $col]
		    adjustElem $win oldText oldImageWidth \
			       $font $pixels $alignment
		    set oldTextWidth [font measure $font \
				      -displayof $win $oldText]
		    set oldElemWidth [expr {$oldImageWidth + $oldTextWidth}]
		    if {$oldElemWidth == $data($col-width) &&
			$newElemWidth < $oldElemWidth} {
			adjustColumns $win $col yes
		    }
		}
	    }
	}

	-text {
	    if {[string compare $configVals(-state) disabled] == 0} {
		return ""
	    }

	    set font $configVals(-font)
	    set pixels [lindex $data(colList) [expr {2*$col}]]
	    if {$pixels != 0} {			;# convention: static width
		incr pixels $data($col-delta)
	    }
	    set alignment [lindex $data(colList) [expr {2*$col + 1}]]

	    #
	    # Adjust the cell text and the image width
	    #
	    regsub -all "\t|\n" $val " " val
	    set text $val
	    set oldItem [lindex $data(itemList) $row]
	    set key [lindex $oldItem end]
	    set name $key-$col-image
	    if {[info exists data($name)]} {
		set image $data($name)
		set imageWidth [image width $image]
	    } else {
		set image ""
		set imageWidth 0
	    }
	    set oldImageWidth $imageWidth
	    adjustElem $win text imageWidth $font $pixels $alignment

	    if {!$data($col-hide)} {
		#
		# Delete the old cell contents between the
		# two tabs, and insert the text and the image
		#
		findCellTabs $win [expr {$row + 1}].0 $col tabIdx1 tabIdx2
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
	    if {[string compare $configVals(-listvariable) ""] != 0} {
		set traceCmd [list tablelist::traceProc $win]
		trace vdelete ::$configVals(-listvariable) wu $traceCmd
		upvar #0 $configVals(-listvariable) var
		set lastCol [expr {$data(colCount) - 1}]
		set var [lreplace $var $row $row [lrange $newItem 0 $lastCol]]
		trace variable ::$configVals(-listvariable) wu $traceCmd
	    }

	    #
	    # Adjust the columns if necessary
	    #
	    if {$pixels == 0} {			;# convention: dynamic width
		set textWidth [font measure $font -displayof $win $text]
		set newElemWidth [expr {$imageWidth + $textWidth}]
		if {$newElemWidth > $data($col-width)} {
		    set data($col-width) $newElemWidth
		    adjustColumns $win -1 yes
		} else {
		    set oldText [lindex $oldItem $col]
		    adjustElem $win oldText oldImageWidth \
			       $font $pixels $alignment
		    set oldTextWidth [font measure $font \
				      -displayof $win $oldText]
		    set oldElemWidth [expr {$oldImageWidth + $oldTextWidth}]
		    if {$oldElemWidth == $data($col-width) &&
			$newElemWidth < $oldElemWidth} {
			adjustColumns $win $col yes
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
			findCellTabs $win $line.0 $col tabIdx1 tabIdx2
			$w tag add $tag $tabIdx1 $tabIdx2+1c
		    }
		}

		#
		# Save the properly formatted value of val in data($tag)
		#
		set data($tag) [$w tag cget $tag -$optTail]
	    }
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::setupColumns
#
# Updates the value of the -colums configuration option for the tablelist
# widget win by using the width, title, and alignment specifications given in
# the columns argument, and creates the corresponding label widgets if
# createLabels is true.
#------------------------------------------------------------------------------
proc tablelist::setupColumns {win columns createLabels} {
    variable configSpecs
    variable configOpts
    variable alignments
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    set argCount [llength $columns]
    set colConfigVals {}

    #
    # Check the syntax of columns before performing any changes
    #
    if {$argCount == 0} {
	return -code error "expected at least one column specification"
    }
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
    # Save the value of colConfigVals in configVals(-columns)
    #
    set configVals(-columns) $colConfigVals

    #
    # Delete the labels if requested
    #
    if {$createLabels} {
	set children [winfo children $data(hdrTxtFr)]
	foreach w [lrange [lsort $children] 1 end] {
	    destroy $w
	}
    }

    #
    # Build the list data(colList), and create the labels if requested
    #
    set font $configVals(-font)
    set data(colList) {}
    set col 0
    foreach {width title alignment} $configVals(-columns) {
	#
	# Append the width in pixels and the
	# alignment to the list data(colList)
	#
	if {$width > 0} {		;# convention: width in characters
	    ### set zeroStr [string repeat 0 $count]
	    set zeroStr ""
	    for {set n 0} {$n < $width} {incr n} {
		append zeroStr 0
	    }
	    set pixels [font measure $font -displayof $win $zeroStr]
	} elseif {$width < 0} {		;# convention: width in pixels
	    set pixels [expr {(-1)*$width}]
	} else {			;# convention: dynamic width
	    set pixels 0
	}
	lappend data(colList) $pixels $alignment

	if {$createLabels} {
	    if {![info exists data($col-hide)]} {
		set data($col-hide) 0
	    }
	    if {![info exists data($col-resizable)]} {
		set data($col-resizable) 1
	    }
	    if {![info exists data($col-showarrow)]} {
		set data($col-showarrow) 1
	    }
	    if {![info exists data($col-sortmode)]} {
		set data($col-sortmode) ascii
	    }
	    if {![info exists data($col-delta)]} {
		set data($col-delta) 0
	    }

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
		if {[llength $configSpecs($opt)] != 1 &&
		    [regexp {[lc]} [lindex $configSpecs($opt) 2] optGrp]} {
		    if {[string compare $optGrp l] == 0} {
			set optTail [string range $opt 6 end]
			if {[info exists data($col$opt)]} {
			    $w configure -$optTail $data($col$opt)
			} else {
			    $w configure -$optTail $configVals($opt)
			}
		    } else {
			$w configure $opt $configVals($opt)
		    }
		}
	    }
	    if {[info exists data($col-labelimage)]} {
		doColConfig $col $win -labelimage $data($col-labelimage)
	    }
	    catch {$w configure -state $configVals(-state)}

	    #
	    # Define some mouse bindings for the label, needed for resizing
	    #
	    bind $w <Enter>	[list tablelist::labelEnter	 $win $col %x]
	    bind $w <Motion>	[list tablelist::labelEnter	 $win $col %x]
	    bind $w <Button-1>	[list tablelist::labelB1Down	 $win $col %x]
	    bind $w <B1-Motion>	[list tablelist::labelB1Motion	 $win %x %y]
	    bind $w <B1-Enter>	[list tablelist::labelB1Enter	 $win]
	    bind $w <B1-Leave>	[list tablelist::labelB1Leave	 $win %x %y]
	    bind $w <ButtonRelease-1> [list tablelist::labelB1Up $win]
	    bind $w <<Button3>>	[list tablelist::labelB3Down	 $win $col]
	}

	incr col
    }

    #
    # Save the number of columns in data(colCount)
    #
    set oldColCount $data(colCount)
    set data(colCount) $col

    if {$data(colCount) < $oldColCount} {
	#
	# Clean up the data associated with the deleted columns
	#
	for {set col $data(colCount)} {$col < $oldColCount} {incr col} {
	    if {[info exists data($col-redispId)]} {
		after cancel $data($col-redispId)
	    }
	    foreach name [array names data $col-*] {
		unset data($name)
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
# colOfUnknownWidth argument specifies the dynamic-width column(s) whose width
# is to be computed when performing these operations (a number or "all").  The
# stretchCols argument specifies whether to stretch the stretchable columns.
#------------------------------------------------------------------------------
proc tablelist::adjustColumns {win colOfUnknownWidth stretchCols} {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    #
    # Configure the labels, place them in the header frame, and compute
    # the positions of the tab stops to be set in the body text widget
    #
    set stretchCols [expr {$stretchCols ? 1 : 0}]	;# for versions < 8.3.3
    set font $configVals(-font)
    set charWidth [font measure $font -displayof $win 0]
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
	# Adjust the col`th label
	#
	if {$pixels != 0} {
	    incr pixels $data($col-delta)
	}
	if {[info exists data($col-labelalign)]} {
	    set labelAlignment $data($col-labelalign)
	} else {
	    set labelAlignment $alignment
	}
	adjustLabel $win $col $charWidth $pixels $labelAlignment

	if {$pixels == 0} {			;# convention: dynamic width
	    #
	    # Compute the column width if requested
	    #
	    if {[string compare $colOfUnknownWidth all] == 0 ||
		$colOfUnknownWidth == $col} {
		computeColWidth $win $col
	    }

	    set pixels $data($col-width)
	    incr pixels $data($col-delta)
	}

	if {$col == $data(arrowCol)} {
	    #
	    # Place the canvas to the left side of the label if the
	    # column is right-justified and to its right side otherwise
	    #
	    set canvas $data(hdrTxtFrCanv)
	    if {[string compare $alignment right] == 0} {
		place $canvas -in $w -anchor w -bordermode outside \
			      -relx 0.0 -x $charWidth -rely 0.5
	    } else {
		place $canvas -in $w -anchor e -bordermode outside \
			      -relx 1.0 -x -$charWidth -rely 0.5
	    }
	    raise $canvas
	}

	#
	# Place the label in the header frame
	#
	set labelPixels [expr {$pixels + 2*$charWidth}]
	place $w -x $x -relheight 1.0 -width $labelPixels
	incr x $labelPixels

	#
	# Append a tab stop and the alignment to the tabs list
	#
	incr data(hdrPixels) $charWidth
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
	incr data(hdrPixels) $charWidth
	lappend tabs $data(hdrPixels) left

	incr col
    }

    #
    # Apply the value of tabs to the body text widget
    #
    $data(body) configure -tabs $tabs

    #
    # Adjust the width and height of the frames data(hdrTxtFr) and data(hdr)
    #
    $data(hdrTxtFr) configure -width $data(hdrPixels)
    if {$configVals(-width) <= 0} {
	if {$stretchCols} {
	    $data(hdr) configure -width $data(hdrPixels)
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
proc tablelist::adjustLabel {win col charWidth pixels alignment} {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    #
    # Apply some configuration options to the label and its children (if any)
    #
    set w $data(hdrTxtFrLbl)$col
    switch $alignment {
	left	{ set anchor w }
	right	{ set anchor e }
	center	{ set anchor center }
    }
    $w configure -anchor $anchor -justify $alignment \
		  -padx [expr {$charWidth - [$w cget -borderwidth]}]
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
    set title [lindex $configVals(-columns) [expr {3*$col + 1}]]
    set labelFont [$w cget -font]
    if {$col == $data(arrowCol)} {
	if {[font metrics $labelFont -displayof $w -fixed]} {
	    set spaces "  "
	} else {
	    set spaces "    "
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
	    set gap [font measure $configVals(-font) -displayof $win " "]
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
		    set line [strRangeExt $win $line \
			      $labelFont $lessPixels $alignment]
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
	    set gap [font measure $configVals(-font) -displayof $win " "]
	    if {$imageWidth + $gap <= $lessPixels} {
		incr lessPixels -[expr {$imageWidth + $gap}]
		set lines {}
		foreach line [split $title \n] {
		    set line [strRangeExt $win $line \
			      $labelFont $lessPixels $alignment]
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
		place $w.il -anchor w -bordermode outside -relx 0.0 \
			    -x $charWidth -rely 0.5
		if {[string compare $text ""] != 0} {
		    set textX [expr {$charWidth + [winfo reqwidth $w.il]}]
		    place $w.tl -anchor w -bordermode outside -relx 0.0 \
				-x $textX -rely 0.5
		}
	    }

	    right {
		place $w.il -anchor e -bordermode outside -relx 1.0 \
			    -x -$charWidth -rely 0.5
		if {[string compare $text ""] != 0} {
		    set textX [expr {-$charWidth - [winfo reqwidth $w.il]}]
		    place $w.tl -anchor e -bordermode outside -relx 1.0 \
				-x $textX -rely 0.5
		}
	    }

	    center {
		if {[string compare $text ""] == 0} {
		    place $w.il -anchor center -bordermode outside \
				-relx 0.5 -x 0 -rely 0.5
		} else {
		    set halfWidth [expr {([winfo reqwidth $w.il] + \
					  [winfo reqwidth $w.tl]) / 2}]
		    place $w.il -anchor w -bordermode outside -relx 0.5 \
				-x -$halfWidth -rely 0.5
		    place $w.tl -anchor e -bordermode outside -relx 0.5 \
				-x $halfWidth -rely 0.5
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
    upvar ::tablelist::ns${win}::configVals configVals

    set data($col-width) 0
    set font $configVals(-font)

    #
    # Column elements
    #
    foreach item $data(itemList) {
	if {$col < [expr {[llength $item] - 1}]} {
	    set text [lindex $item $col]
	    set key [lindex $item end]
	    set name $key-$col-image
	    if {[info exists data($name)]} {
		set imageWidth [image width $data($name)]
	    } else {
		set imageWidth 0
	    }
	    adjustElem $win text imageWidth $font 0 left
	    set textWidth [font measure $font -displayof $win $text]
	    set elemWidth [expr {$imageWidth + $textWidth}]
	    if {$elemWidth > $data($col-width)} {
		set data($col-width) $elemWidth
	    }
	}
    }

    #
    # Column label
    #
    set w $data(hdrTxtFrLbl)$col
    if {[info exists data($col-labelimage)]} {
	set title [lindex $configVals(-columns) [expr {3*$col + 1}]]
	if {[string compare $title ""] == 0} {		;# image w/o text
	    set netLabelWidth [winfo reqwidth $w.il]
	} else {					;# both image and text
	    set netLabelWidth [expr {[winfo reqwidth $w.il] +
				     [winfo reqwidth $w.tl]}]
	}
    } else {						;# no image
	set netLabelWidth [expr {[winfo reqwidth $w] -
				 2*[font measure $font -displayof $win 0]}]
    }
    if {$netLabelWidth > $data($col-width)} {
	set data($col-width) $netLabelWidth
    }
}

#------------------------------------------------------------------------------
# tablelist::adjustHeaderHeight
#
# Sets the height of the header frame of the tablelist widget win to the max.
# height of its children and adjusts the y-coordinate of the canvas containing
# an up- or down-arrow.
#------------------------------------------------------------------------------
proc tablelist::adjustHeaderHeight win {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    #
    # Compute the max. label height
    #
    set maxLabelHeight 0
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

	    set reqHeight [expr {[winfo reqheight $c] +
				 2*[$w cget -borderwidth]}]
	    if {$reqHeight > $maxLabelHeight} {
		set maxLabelHeight $reqHeight
	    }
	}
    }

    #
    # Set the height of the header frame
    #
    $data(hdrTxtFr) configure -height $maxLabelHeight
    if {$configVals(-showlabels)} {
	$data(hdr) configure -height $maxLabelHeight
	place $data(hdrTxt) -relheight 1.0 -relwidth 1.0
    } else {
	$data(hdr) configure -height 1
	place forget $data(hdrTxt)
    }
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

    set data(stretchId) [after idle [list tablelist::stretchColumns $win]]
}

#------------------------------------------------------------------------------
# tablelist::stretchColumns
#
# Stretches the stretchable columns to fill the tablelist window win
# horizontally.
#------------------------------------------------------------------------------
proc tablelist::stretchColumns win {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    if {[info exists data(stretchId)]} {
	unset data(stretchId)
    }

    set forceAdjust $data(forceAdjust)
    set data(forceAdjust) 0

    if {$data(hdrPixels) == 0 || $configVals(-width) <= 0} {
	return ""
    }

    #
    # Compute the total number delta of pixels
    # by which the columns are to be stretched
    #
    update idletasks
    set font $configVals(-font)
    set charWidth [font measure $font -displayof $win 0]
    set delta [winfo width $data(hdr)]
    set col 0
    foreach {pixels alignment} $data(colList) {
	if {$data($col-hide)} {
	    incr col
	    continue
	}
	if {$pixels == 0} {			;# convention: dynamic width
	    set pixels $data($col-width)
	}
	incr delta -[expr {$pixels + 2*$charWidth}]

	incr col
    }
    if {$delta < 0} {
	set delta 0
    }

    #
    # Get the list of the numerical indices of the stretchable columns
    #
    set stretchableCols {}
    if {[string compare $configVals(-stretch) all] == 0} {
	for {set col 0} {$col < $data(colCount)} {incr col} {
	    lappend stretchableCols $col
	}
    } else {
	foreach col $configVals(-stretch) {
	    lappend stretchableCols [colIndex $win $col]
	}
    }

    #
    # Get the number of columns to be stretched and the mean value of delta
    #
    set count 0
    set lastColToStretch -1
    for {set col 0} {$col < $data(colCount)} {incr col} {
	if {[lsearch -exact $stretchableCols $col] >= 0 && !$data($col-hide)} {
	    incr count
	    set lastColToStretch $col
	}
    }
    if {$count == 0} {
	set meanVal 0
	set rest 0
    } else {
	set meanVal [expr {$delta / $count}]
	set rest [expr {$delta % $count}]
    }
    if {$count == 0 && !$forceAdjust} {
	return ""
    }

    #
    # Add the mean value of delta to the widths of the stretchable
    # columns and the rest to the width of the last such column
    #
    set col 0
    foreach {width title alignment} $configVals(-columns) \
	    {pixels alignment} $data(colList) {
	if {[lsearch -exact $stretchableCols $col] >= 0 && !$data($col-hide)} {
	    set oldDelta $data($col-delta)
	    set data($col-delta) $meanVal
	    if {$col == $lastColToStretch} {
		incr data($col-delta) $rest
	    }
	    if {$pixels != 0 && $data($col-delta) != $oldDelta} {
		redisplayColWhenIdle $win $col
	    }
	} else {
	    set data($col-delta) 0
	}

	incr col
    }

    #
    # Adjust the columns
    #
    adjustColumns $win -1 no
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

    set data(redispId) [after idle [list tablelist::redisplay $win 0 end]]

    #
    # Cancel the execution of all delayed redisplayCol commands
    #
    foreach name [array names data *-redispId] {
	after cancel $data($name)
	unset data($name)
    }
}

#------------------------------------------------------------------------------
# tablelist::redisplay
#
# Redisplays the items of the tablelist widget win, in the range specified by
# first and last.
#------------------------------------------------------------------------------
proc tablelist::redisplay {win first last} {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    if {$first == 0 && [string first $last end] == 0 &&
	[info exists data(redispId)]} {
	unset data(redispId)
    }

    if {$first < 0} {
	return ""
    }
    if {[string first $last end] == 0} {
	set last [expr {$data(itemCount) - 1}]
    }

    set w $data(body)
    set font $configVals(-font)
    for {set idx $first; set line [expr {$first + 1}]} {$idx <= $last} \
	{incr idx; incr line} {
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
	set item [lindex $data(itemList) $idx]
	set keyIdx [expr {[llength $item] - 1}]
	set key [lindex $item end]
	set newItem {}
	set col 0
	foreach {pixels alignment} $data(colList) {
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
	    set name $key-$col-image
	    if {[info exists data($name)]} {
		set image $data($name)
		set imageWidth [image width $image]
	    } else {
		set image ""
		set imageWidth 0
	    }
	    if {$pixels != 0} {			;# convention: static width
		incr pixels $data($col-delta)
	    }
	    adjustElem $win text imageWidth $font $pixels $alignment

	    #
	    # Insert the text and the image
	    #
	    set tagNames {}
	    foreach opt {-background -foreground} {
		foreach tag [list $col$opt $key$opt $key-$col$opt] {
		    if {[info exists data($tag)]} {
			lappend tagNames $tag
		    }
		}
	    }
	    if {$imageWidth == 0} {
		$w insert $line.end \t$text\t $tagNames
	    } else {
		$w insert $line.end \t\t $tagNames
		insertElem $w $line.end-1c $text $image $imageWidth $alignment
	    }

	    incr col
	}
	lappend newItem $key
	set data(itemList) [lreplace $data(itemList) $idx $idx $newItem]

	#
	# Select the item if it was selected before
	#
	if {$selected} {
	    selectionSubCmd $win set $idx $idx 2
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::redisplayColWhenIdle
#
# Arranges for the elements of the col'th column of the tablelist widget win to
# be redisplayed at idle time.
#------------------------------------------------------------------------------
proc tablelist::redisplayColWhenIdle {win col} {
    upvar ::tablelist::ns${win}::data data

    if {[info exists data($col-redispId)] || [info exists data(redispId)] ||
	$data(itemCount) == 0} {
	return ""
    }

    set data($col-redispId) \
	[after idle [list tablelist::redisplayCol $win $col 0 end]]
}

#------------------------------------------------------------------------------
# tablelist::redisplayCol
#
# Redisplays the elements of the col'th column of the tablelist widget win, in
# the range specified by first and last.
#------------------------------------------------------------------------------
proc tablelist::redisplayCol {win col first last} {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    if {$first == 0 && [string first $last end] == 0 &&
	[info exists data($col-redispId)]} {
	unset data($col-redispId)
    }

    if {$data($col-hide) || $first < 0} {
	return ""
    }
    if {[string first $last end] == 0} {
	set last [expr {$data(itemCount) - 1}]
    }

    set w $data(body)
    set font $configVals(-font)
    set pixels [lindex $data(colList) [expr {2*$col}]]
    if {$pixels != 0} {				;# convention: static width
	incr pixels $data($col-delta)
    }
    set alignment [lindex $data(colList) [expr {2*$col + 1}]]

    for {set idx $first; set line [expr {$first + 1}]} {$idx <= $last} \
	{incr idx; incr line} {
	#
	# Adjust the cell text and the image width
	#
	set item [lindex $data(itemList) $idx]
	set text [lindex $item $col]
	set key [lindex $item end]
	set name $key-$col-image
	if {[info exists data($name)]} {
	    set image $data($name)
	    set imageWidth [image width $image]
	} else {
	    set image ""
	    set imageWidth 0
	}
	adjustElem $win text imageWidth $font $pixels $alignment

	#
	# Delete the old cell contents between the
	# two tabs, and insert the text and the image
	#
	findCellTabs $win $line.0 $col tabIdx1 tabIdx2
	$w delete $tabIdx1+1c $tabIdx2
	insertElem $w $tabIdx1+1c $text $image $imageWidth $alignment
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
    upvar ::tablelist::ns${win}::configVals configVals

    if {[string compare $varName ""] == 0} {
	#
	# If there is an old list variable associated with the
	# widget then remove the trace set on this variable
	#
	if {[string compare $configVals(-listvariable) ""] != 0} {
	    trace vdelete ::$configVals(-listvariable) wu \
		  [list tablelist::traceProc $win]
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
	traceProc $win $name1 $name2 w
    } else {
	#
	# Set ::$varName according to the value of data(itemList)
	#
	set lastCol [expr {$data(colCount) - 1}]
	set ::$varName {}
	foreach item $data(itemList) {
	    lappend ::$varName [lrange $item 0 $lastCol]
	}
    }

    #
    # If there is an old list variable associated with the
    # widget then remove the trace set on this variable,
    # and in any case set a trace on the new list variable
    #
    set traceCmd [list tablelist::traceProc $win]
    if {[string compare $configVals(-listvariable) ""] != 0} {
	trace vdelete ::$configVals(-listvariable) wu $traceCmd
    }
    trace variable ::$varName wu $traceCmd
}

#------------------------------------------------------------------------------
# tablelist::traceProc
#
# This procedure is executed whenever the global variable specified by varName
# is written or unset.  It makes sure that the contents of the widget will be
# synchronized with the value of the variable at idle time, and that the
# variable is recreated if it was unset.
#------------------------------------------------------------------------------
proc tablelist::traceProc {win varName index op} {
    upvar ::tablelist::ns${win}::data data

    if {[string compare $index ""] != 0} {
	set varName ${varName}($index)
    }

    switch $op {
	w {
	    #
	    # Check whether the value of the variable
	    # ::$varName is a valid list of good lists
	    #
	    foreach item [set ::$varName] {
		llength $item
	    }

	    if {![info exists data(syncId)]} {
		#
		# Arrange for the contents of the widget to be synchronized
		# with the value of the variable ::$varName at idle time.
		#
		set data(syncId) [after idle [list tablelist::synchronize $win]]

		#
		# Cancel the execution of all delayed redisplay and
		# redisplayCol commands, to make sure that the synchronize
		# command will be invoked first; the latter will then
		# schedule a redisplay command for execution at idle time
		#
		foreach name [array names data *redispId] {
		    after cancel $data($name)
		    unset data($name)
		}
	    }
	}

	u {
	    #
	    # Recreate the variable ::$varName by setting it according to
	    # the value of data(itemList), and set the trace on it again
	    #
	    set lastCol [expr {$data(colCount) - 1}]
	    set ::$varName {}
	    foreach item $data(itemList) {
		lappend ::$varName [lrange $item 0 $lastCol]
	    }
	    trace variable ::$varName wu [list tablelist::traceProc $win]
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::synchronize
#
# This procedure is invoked either as an idle callback after the list variable
# associated with the tablelist widget win was written, or directly, upon
# execution of some widget commands.  It makes sure that the contents of the
# widget is synchronized with the value of the list variable.
#------------------------------------------------------------------------------
proc tablelist::synchronize win {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

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

    upvar #0 $configVals(-listvariable) var
    set newCount [llength $var]
    if {$newCount < $data(itemCount)} {
	#
	# Delete the items with indices >= newCount from the widget
	#
	set updateCount $newCount
	deleteSubCmd $win $newCount [expr {$data(itemCount) - 1}] no
    } elseif {$newCount > $data(itemCount)} {
	#
	# Insert the items of var with indices
	# >= data(itemCount) into the widget
	#
	set updateCount $data(itemCount)
	insertSubCmd $win $data(itemCount) \
		     [lrange $var $data(itemCount) end] no
    } else {
	set updateCount $newCount
    }

    #
    # Update the first updateCount items of the internal list
    #
    for {set idx 0} {$idx < $updateCount} {incr idx} {
	set oldItem [lindex $data(itemList) $idx]
	set varItem [lindex $var $idx]

	set newItem {}
	for {set col 0} {$col < $data(colCount)} {incr col} {
	    set elem [lindex $varItem $col]
	    regsub -all "\t|\n" $elem " " elem
	    lappend newItem $elem
	}
	lappend newItem [lindex $oldItem end]

	if {[string compare $oldItem $newItem] != 0} {
	    set data(itemList) \
		[lreplace $data(itemList) $idx $idx $newItem]
	}
    }

    #
    # Adjust the columns and make sure the
    # items will be redisplayed at idle time
    #
    adjustColumns $win all yes
    redisplayWhenIdle $win
}

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
    variable configSpecs
    variable labelConfigSpecs
    variable colConfigSpecs
    variable rowConfigSpecs
    variable cellConfigSpecs
    variable cmdOpts
    variable scanCmdOpts
    variable selCmdOpts
    variable sortOrders
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    set argCount [llength $argList]
    if {$argCount == 0} {
	mwutil::wrongNumArgs "$win option ?arg arg ...?"
    }

    set cmd [mwutil::fullOpt "option" [lindex $argList 0] $cmdOpts]
    switch $cmd {
	activate {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd index"
	    }

	    synchronize $win
	    set index [rowIndex $win [lindex $argList 1] no]
	    if {[string compare $configVals(-state) disabled] == 0} {
		return ""
	    } else {
		return [activateSubCmd $win $index]
	    }
	}

	attrib {
	    return [mwutil::attribSubCmd $win [lrange $argList 1 end]]
	}

	bbox {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd index"
	    }

	    synchronize $win
	    set index [rowIndex $win [lindex $argList 1] no]
	    return [bboxSubCmd $win $index]
	}

	bodypath {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    return $data(body)
	}

	cellcget {
	    if {$argCount != 3} {
		mwutil::wrongNumArgs "$win $cmd cellIndex option"
	    }

	    #
	    # Check the cell index
	    #
	    synchronize $win
	    set cell [lindex $argList 1]
	    set cellIdx [cellIndex $win $cell]
	    scan $cellIdx %d,%d row col
	    if {$row < 0 || $row >= $data(itemCount) ||
		$col < 0 || $col >= $data(colCount)} {
		return -code error \
		       "cell index \"$cell\" out of range"
	    }

	    #
	    # Return the value of the specified cell configuration option
	    #
	    set opt [mwutil::fullConfigOpt [lindex $argList 2] cellConfigSpecs]
	    switch -- $opt {
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

	cellconfigure {
	    if {$argCount < 2} {
		mwutil::wrongNumArgs "$win $cmd cellIndex ?option? ?value?\
				      ?option value ...?"
	    }

	    #
	    # Check the cell index
	    #
	    synchronize $win
	    set cell [lindex $argList 1]
	    set cellIdx [cellIndex $win $cell]
	    scan $cellIdx %d,%d row col
	    if {$row < 0 || $row >= $data(itemCount) ||
		$col < 0 || $col >= $data(colCount)} {
		return -code error \
		       "cell index \"$cell\" out of range"
	    }

	    #
	    # Build an array from the cell configuration options and their
	    # values, and pass it to the procedure mwutil::configSubCmd
	    #
	    set item [lindex $data(itemList) $row]
	    set cellConfigVals(-text) [lindex $item $col]
	    set key [lindex $item end]
	    foreach opt {-background -foreground -image
			 -selectbackground -selectforeground} {
		if {[info exists data($key-$col$opt)]} {
		    set cellConfigVals($opt) $data($key-$col$opt)
		} else {
		    set cellConfigVals($opt) ""
		}
	    }
	    return [mwutil::configSubCmd $win cellConfigSpecs cellConfigVals \
					 "tablelist::doCellConfig $row $col" \
					 [lrange $argList 2 end]]
	}

	cellindex {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd cellIndex"
	    }

	    synchronize $win
	    return [cellIndex $win [lindex $argList 1]]
	}

	cget {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd option"
	    }

	    #
	    # Return the value of the specified configuration option
	    #
	    set opt [mwutil::fullConfigOpt [lindex $argList 1] configSpecs]
	    return $configVals($opt)
	}

	columncget {
	    if {$argCount != 3} {
		mwutil::wrongNumArgs "$win $cmd columnIndex option"
	    }

	    #
	    # Check the column index
	    #
	    synchronize $win
	    set colArg [lindex $argList 1]
	    set col [colIndex $win $colArg]
	    if {$col < 0 || $col >= $data(colCount)} {
		return -code error \
		       "column index \"$colArg\" out of range"
	    }

	    #
	    # Return the value of the specified column configuration option
	    #
	    set opt [mwutil::fullConfigOpt [lindex $argList 2] colConfigSpecs]
	    switch -- $opt {
		-align {
		    return [lindex $configVals(-columns) [expr {3*$col + 2}]]
		}
		-title {
		    return [lindex $configVals(-columns) [expr {3*$col + 1}]]
		}
		-width {
		    return [lindex $configVals(-columns) [expr {3*$col}]]
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

	columnconfigure {
	    if {$argCount < 2} {
		mwutil::wrongNumArgs "$win $cmd columnIndex ?option? ?value?\
				      ?option value ...?"
	    }

	    #
	    # Check the column index
	    #
	    synchronize $win
	    set colArg [lindex $argList 1]
	    set col [colIndex $win $colArg]
	    if {$col < 0 || $col >= $data(colCount)} {
		return -code error \
		       "column index \"$colArg\" out of range"
	    }

	    #
	    # Build an array from the column configuration options and their
	    # values, and pass it to the procedure mwutil::configSubCmd
	    #
	    set n [expr {3*$col}]
	    foreach opt {-width -title -align} {
		set colConfigVals($opt) [lindex $configVals(-columns) $n]
		incr n
	    }
	    foreach opt {-background -foreground -hide -labelalign
			 -labelbackground -labelborderwidth -labelcommand
			 -labelfont -labelforeground -labelheight -labelimage
			 -labelpady -labelrelief -resizable -selectbackground
			 -selectforeground -showarrow -sortcommand -sortmode} {
		if {[info exists data($col$opt)]} {
		    set colConfigVals($opt) $data($col$opt)
		} else {
		    set colConfigVals($opt) ""
		}
	    }
	    return [mwutil::configSubCmd $win colConfigSpecs colConfigVals \
					 "tablelist::doColConfig $col" \
					 [lrange $argList 2 end]]
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
	    return [colIndex $win [lindex $argList 1]]
	}

	configure {
	    return [mwutil::configSubCmd $win configSpecs configVals \
					 tablelist::doConfig \
					 [lrange $argList 1 end]]
	}

	curselection {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    synchronize $win
	    return [curselectionSubCmd $win]
	}

	delete {
	    if {$argCount < 2 || $argCount > 3} {
		mwutil::wrongNumArgs "$win $cmd firstIndex ?lastIndex?"
	    }

	    synchronize $win
	    set first [rowIndex $win [lindex $argList 1] no]
	    if {$argCount == 3} {
		set last [rowIndex $win [lindex $argList 2] no]
	    } else {
		set last $first
	    }
	    if {[string compare $configVals(-state) disabled] == 0} {
		return ""
	    } else {
		return [deleteSubCmd $win $first $last yes]
	    }
	}

	get {
	    if {$argCount < 2 || $argCount > 3} {
		mwutil::wrongNumArgs "$win $cmd firstIndex ?lastIndex?"
	    }

	    synchronize $win
	    set first [rowIndex $win [lindex $argList 1] no]
	    if {$argCount == 3} {
		set last [rowIndex $win [lindex $argList 2] no]
		set idxCount 2
	    } else {
		set last $first
		set idxCount 1
	    }
	    return [getSubCmd $win $first $last $idxCount]
	}

	index {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd index"
	    }

	    synchronize $win
	    return [rowIndex $win [lindex $argList 1] yes]
	}

	insert {
	    if {$argCount < 2} {
		mwutil::wrongNumArgs "$win $cmd index ?item item ...?"
	    }

	    synchronize $win
	    set index [rowIndex $win [lindex $argList 1] yes]
	    if {[string compare $configVals(-state) disabled] == 0} {
		return ""
	    } else {
		return [insertSubCmd $win $index [lrange $argList 2 end] yes]
	    }
	}

	labelpath {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd columnIndex"
	    }

	    set col [colIndex $win [lindex $argList 1]]
	    return $data(hdrTxtFrLbl)$col
	}

	labels {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    set children [winfo children $data(hdrTxtFr)]
	    return [lrange [lsort $children] 1 end]
	}

	nearest {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd y"
	    }

	    set y [lindex $argList 1]
	    format %d $y		;# integer check with error message
	    synchronize $win
	    return [rowIndex $win @0,$y no]
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
	    return [cellIndex $win @$x,$y]
	}

	nearestcolumn {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd x"
	    }

	    set x [lindex $argList 1]
	    format %d $x		;# integer check with error message
	    synchronize $win
	    return [colIndex $win @$x,0]
	}

	resetsortinfo {
	    if {$argCount != 1} {
		mwutil::wrongNumArgs "$win $cmd"
	    }

	    set data(sortCol) -1
	    set data(sortOrder) {}

	    place forget $data(hdrTxtFrCanv)
	    set colOfUnknownWidth $data(arrowCol)
	    set data(arrowCol) -1
	    synchronize $win
	    adjustColumns $win $colOfUnknownWidth yes
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
	    set row [rowIndex $win $rowArg no]
	    if {$row < 0 || $row >= $data(itemCount)} {
		return -code error \
		       "row index \"$rowArg\" out of range"
	    }

	    #
	    # Return the value of the specified row configuration option
	    #
	    set opt [mwutil::fullConfigOpt [lindex $argList 2] rowConfigSpecs]
	    switch -- $opt {
		-text {
		    set item [lindex $data(itemList) $row]
		    return [lrange $item 0 [expr {$data(colCount) - 1}]]
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
	    set row [rowIndex $win $rowArg no]
	    if {$row < 0 || $row >= $data(itemCount)} {
		return -code error \
		       "row index \"$rowArg\" out of range"
	    }

	    #
	    # Build an array from the row configuration options and their
	    # values, and pass it to the procedure mwutil::configSubCmd
	    #
	    set item [lindex $data(itemList) $row]
	    set rowConfigVals(-text) \
		[lrange $item 0 [expr {$data(colCount) - 1}]]
	    set key [lindex $item end]
	    foreach opt {-background -foreground
			 -selectbackground -selectforeground} {
		if {[info exists data($key$opt)]} {
		    set rowConfigVals($opt) $data($key$opt)
		} else {
		    set rowConfigVals($opt) ""
		}
	    }
	    return [mwutil::configSubCmd $win rowConfigSpecs rowConfigVals \
					 "tablelist::doRowConfig $row" \
					 [lrange $argList 2 end]]
	}

	scan {
	    if {$argCount != 4} {
		mwutil::wrongNumArgs "$win $cmd mark|dragto x y"
	    }

	    set x [lindex $argList 2]
	    set y [lindex $argList 3]
	    format %d $x		;# integer check with error message
	    format %d $y		;# integer check with error message
	    set opt [mwutil::fullOpt "option" [lindex $argList 1] $scanCmdOpts]
	    synchronize $win
	    return [scanSubCmd $win $opt $x $y]
	}

	see {
	    if {$argCount != 2} {
		mwutil::wrongNumArgs "$win $cmd index"
	    }

	    synchronize $win
	    set index [rowIndex $win [lindex $argList 1] no]
	    return [seeSubCmd $win $index]
	}

	selection {
	    if {$argCount < 3 || $argCount > 4} {
		mwutil::wrongNumArgs "$win $cmd option index ?index?"
	    }

	    synchronize $win
	    set first [rowIndex $win [lindex $argList 2] no]
	    if {$argCount == 4} {
		set last [rowIndex $win [lindex $argList 3] no]
		set idxCount 2
	    } else {
		set last $first
		set idxCount 1
	    }
	    set opt [mwutil::fullOpt "option" [lindex $argList 1] $selCmdOpts]
	    if {[string compare $configVals(-state) disabled] == 0} {
		return ""
	    } else {
		return [selectionSubCmd $win $opt $first $last $idxCount]
	    }
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

	    #
	    # Check the column index
	    #
	    set colArg [lindex $argList 1]
	    set col [colIndex $win $colArg]
	    if {$col < 0 || $col >= $data(colCount)} {
		return -code error \
		       "column index \"$colArg\" out of range"
	    }

	    if {$argCount == 2} {
		set order -increasing
	    } else {
		set order [mwutil::fullOpt "option" \
			   [lindex $argList 2] $sortOrders]
	    }
	    synchronize $win
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

	xview {
	    synchronize $win
	    return [xviewSubCmd $win [lrange $argList 1 end]]
	}

	yview {
	    synchronize $win
	    return [yviewSubCmd $win [lrange $argList 1 end]]
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

    #
    # Adjust the index to fit within the existing items
    #
    if {$index >= $data(itemCount)} {
	set index [expr {$data(itemCount) - 1}]
    }
    if {$index < 0} {
	set index 0
    }

    #
    # If the focus is currently on the body text widget
    # then move the "active" tag to the corresponding line
    #
    set w $data(body)
    if {[string compare [focus -displayof $win] $w] == 0} {
	set line [expr {$data(activeIdx) + 1}]
	$w tag remove active $line.0 $line.end

	set line [expr {$index + 1}]
	$w tag add active $line.0 $line.end
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
    upvar ::tablelist::ns${win}::configVals configVals

    set w $data(body)
    set dlineinfo [$w dlineinfo [expr {$index + 1}].0]
    if {$data(itemCount) == 0 || [string compare $dlineinfo ""] == 0} {
	return {}
    }

    foreach {x y width height baselinePos} $dlineinfo {
	lappend bbox [expr {$x + [winfo x $w]}] \
		     [expr {$y + [winfo y $w]}] \
		     $width \
		     [expr {$height - 2*$configVals(-selectborderwidth) - 1}]
    }
    return $bbox
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
proc tablelist::deleteSubCmd {win first last updateListVar} {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    #
    # Adjust the range to fit within the existing items
    #
    if {$first < 0} {
	set first 0
    }
    if {$last >= $data(itemCount)} {
	set last [expr {$data(itemCount) - 1}]
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
    set font $configVals(-font)
    if {$count == $data(itemCount)} {
	set colWidthsChanged yes			;# just to save time
    } else {
	set colWidthsChanged no
	for {set idx $first} {$idx <= $last} {incr idx} {
	    set item [lindex $data(itemList) $idx]
	    set key [lindex $item end]
	    set col 0
	    foreach {pixels alignment} $data(colList) {
		if {$data($col-hide) || $pixels != 0} {
		    incr col
		    continue
		}

		set text [lindex $item $col]
		set name $key-$col-image
		if {[info exists data($name)]} {
		    set imageWidth [image width $data($name)]
		} else {
		    set imageWidth 0
		}
		adjustElem $win text imageWidth $font $pixels $alignment
		set textWidth [font measure $font -displayof $win $text]
		set elemWidth [expr {$imageWidth + $textWidth}]
		if {$elemWidth == $data($col-width)} {
		    set colWidthsChanged yes
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
    # Delete the given items and their tags from the body text widget
    #
    $w delete [expr {$first + 1}].0 [expr {$last + 2}].0
    for {set idx $first} {$idx <= $last} {incr idx} {
	set item [lindex $data(itemList) $idx]
	set key [lindex $item end]

	foreach opt {-background -foreground
		     -selectbackground -selectforeground} {
	    set tag $key$opt
	    if {[info exists data($tag)]} {
		$w tag delete $tag
		unset data($tag)
	    }
	}

	for {set col 0} {$col < $data(colCount)} {incr col} {
	    foreach opt {-background -foreground
			 -selectbackground -selectforeground} {
		set tag $key-$col$opt
		if {[info exists data($tag)]} {
		    $w tag delete $tag
		    unset data($tag)
		}
	    }
	    set name $key-$col-image
	    if {[info exists data($name)]} {
		unset data($name)
	    }
	}
    }

    #
    # Delete the given items from the internal list
    #
    set data(itemList) [lreplace $data(itemList) $first $last]

    #
    # Delete the given items from the list variable if needed
    #
    if {$updateListVar && [string compare $configVals(-listvariable) ""] != 0} {
	set traceCmd [list tablelist::traceProc $win]
	trace vdelete ::$configVals(-listvariable) wu $traceCmd
	upvar #0 $configVals(-listvariable) var
	set var [lreplace $var $first $last]
	trace variable ::$configVals(-listvariable) wu $traceCmd
    }
    incr data(itemCount) -$count

    #
    # Adjust the height of the body text widget if necessary
    #
    if {$configVals(-height) <= 0} {
	$w configure -height $data(itemCount)
    }

    #
    # Adjust the columns if necessary
    #
    if {$colWidthsChanged} {
	adjustColumns $win all yes
    }

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
	if {$data(activeIdx) >= $data(itemCount) && $data(itemCount) > 0} {
	    set data(activeIdx) [expr {$data(itemCount) - 1}]
	}
    }

    return ""
}

#------------------------------------------------------------------------------
# tablelist::getSubCmd
#
# This procedure is invoked to process the tablelist get subcommand.
#------------------------------------------------------------------------------
proc tablelist::getSubCmd {win first last idxCount} {
    upvar ::tablelist::ns${win}::data data

    #
    # Adjust the range to fit within the existing items
    #
    if {$first >= $data(itemCount)} {
	return {}
    }
    if {$last >= $data(itemCount)} {
	set last [expr {$data(itemCount) - 1}]
    }
    if {$first < 0} {
	set first 0
    }
    if {$first > $last} {
	return {}
    }

    #
    # Get the given items from the internal list
    #
    set lastCol [expr {$data(colCount) - 1}]
    if {$idxCount == 1} {
	set item [lindex $data(itemList) $first]
	return [lrange $item 0 $lastCol]
    } else {
	set result {}
	for {set idx $first} {$idx <= $last} {incr idx} {
	    set item [lindex $data(itemList) $idx]
	    lappend result [lrange $item 0 $lastCol]
	}
	return $result
    }
}

#------------------------------------------------------------------------------
# tablelist::insertSubCmd
#
# This procedure is invoked to process the tablelist insert subcommand.
#------------------------------------------------------------------------------
proc tablelist::insertSubCmd {win index argList updateListVar} {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    set argCount [llength $argList]
    if {$argCount == 0} {
	return ""
    }

    if {$index < 0} {
	set index 0
    }

    #
    # Check whether each item is a good list
    #
    foreach item $argList {
	llength $item
    }

    #
    # Insert the items into the body text widget and into the internal list
    #
    set w $data(body)
    set font $configVals(-font)
    set colWidthsChanged no
    set idx $index
    set line [expr {$index + 1}]
    foreach item $argList {
	if {$data(itemCount) != 0} {
	    $w insert $line.0 \n
	}

	set newItem {}
	set col 0
	foreach {pixels alignment} $data(colList) {
	    set elem [lindex $item $col]
	    regsub -all "\t|\n" $elem " " elem
	    lappend newItem $elem

	    if {$data($col-hide)} {
		incr col
		continue
	    }

	    #
	    # Update the column width or clip the element if necessary
	    #
	    if {$pixels == 0} {			;# convention: dynamic width
		set elemWidth [font measure $font -displayof $win $elem]
		if {$elemWidth > $data($col-width)} {
		    set data($col-width) $elemWidth
		    set colWidthsChanged yes
		}
	    } else {
		incr pixels $data($col-delta)
		set elem [strRangeExt $win $elem $font $pixels $alignment]
	    }

	    #
	    # Insert the element into the body text widget
	    # with the tags defined for the current column
	    #
	    set tagNames {}
	    foreach opt {-background -foreground} {
		if {[info exists data($col$opt)]} {
		    lappend tagNames $col$opt
		}
	    }
	    $w insert $line.end \t$elem\t $tagNames

	    incr col
	}

	#
	# Insert the item into the list variable if needed
	#
	if {$updateListVar &&
	    [string compare $configVals(-listvariable) ""] != 0} {
	    set traceCmd [list tablelist::traceProc $win]
	    trace vdelete ::$configVals(-listvariable) wu $traceCmd
	    upvar #0 $configVals(-listvariable) var
	    if {$idx == $data(itemCount)} {
		lappend var $newItem		;# this works much faster
	    } else {
		set var [linsert $var $idx $newItem]
	    }
	    trace variable ::$configVals(-listvariable) wu $traceCmd
	}

	#
	# Insert the item into the internal list
	#
	lappend newItem k[incr data(seqNum)]
	if {$idx == $data(itemCount)} {
	    lappend data(itemList) $newItem	;# this works much faster
	} else {
	    set data(itemList) [linsert $data(itemList) $idx $newItem]
	}

	incr idx
	incr line
	incr data(itemCount)
    }

    #
    # Adjust the height of the body text widget if necessary
    #
    if {$configVals(-height) <= 0} {
	$w configure -height $data(itemCount)
    }

    #
    # Adjust the columns if necessary
    #
    if {$colWidthsChanged} {
	adjustColumns $win -1 yes
    }

    #
    # The following is necessary if the tablelist was previously empty
    #
    set fraction [lindex [$data(hdrTxt) xview] 0]
    $w xview moveto $fraction

    #
    # Update the indices anchorIdx and activeIdx
    #
    if {$index <= $data(anchorIdx)} {
	incr data(anchorIdx) $argCount
    }
    if {$index <= $data(activeIdx)} {
	incr data(activeIdx) $argCount
	if {$data(activeIdx) >= $data(itemCount) && $data(itemCount) > 0} {
	    set data(activeIdx) [expr {$data(itemCount) - 1}]
	}
    }

    return ""
}

#------------------------------------------------------------------------------
# tablelist::scanSubCmd
#
# This procedure is invoked to process the tablelist scan subcommand.
#------------------------------------------------------------------------------
proc tablelist::scanSubCmd {win opt x y} {
    upvar ::tablelist::ns${win}::data data

    incr x -[winfo x $data(body)]
    incr y -[winfo y $data(body)]

    $data(body) scan $opt $x $y
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
    $w see [expr {$index + 1}].0
    $w xview moveto $fraction
    return ""
}

#------------------------------------------------------------------------------
# tablelist::selectionSubCmd
#
# This procedure is invoked to process the tablelist selection subcommand.
#------------------------------------------------------------------------------
proc tablelist::selectionSubCmd {win opt first last idxCount} {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    switch $opt {
	anchor {
	    if {$idxCount != 1} {
		mwutil::wrongNumArgs "$win selection $opt index"
	    }

	    #
	    # Adjust the index to fit within the existing items
	    #
	    if {$first >= $data(itemCount)} {
		set first [expr {$data(itemCount) - 1}]
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
	    if {$idxCount != 1} {
		mwutil::wrongNumArgs "$win selection $opt index"
	    }

	    set tagNames [$data(body) tag names [expr {$first + 1}].0]
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
	    if {$last >= $data(itemCount)} {
		set last [expr {$data(itemCount) - 1}]
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

		$w tag add select $line.0 $line.end

		#
		# Handle the -(select)background and -(select)foreground cell
		# and column configuration options for each element of the row
		#
		set item [lindex $data(itemList) $idx]
		set key [lindex $item end]
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
	    if {$configVals(-exportselection) &&
		[llength [$w tag nextrange select 1.0]] != 0} {
		selection own -command \
			  [list ::tablelist::lostSelection $win] $win
	    }

	    return ""
	}
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
    upvar ::tablelist::ns${win}::configVals configVals

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
    # Sort the item list and update the sort info
    #
    if {$col < 0} {				;# not sorting by a column
	if {[string compare $configVals(-sortcommand) ""] == 0} {
	    return -code error \
		   "value of the -sortcommand option is empty"
	}

	set data(itemList) [eval lsort $order \
	    -command $configVals(-sortcommand) {$data(itemList)}]
    } else {					;# sorting by a column
	if {[string compare $data($col-sortmode) command] == 0} {
	    if {[info exists data($col-sortcommand)] &&
		[string compare $data($col-sortcommand) ""] != 0} {
		set mode [list -command $data($col-sortcommand)]
	    } else {
		return -code error \
		       "value of the -sortcommand option for\
			column $col is missing or empty"
	    }
	} else {
	    set mode -$data($col-sortmode)
	}

	set data(itemList) \
	    [eval lsort -index $col $order $mode {$data(itemList)}]
    }
    set data(sortCol) $col
    set data(sortOrder) [string range $order 1 end]

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
    set colOfUnknownWidth $data(arrowCol)
    set data(arrowCol) -1
    adjustColumns $win $colOfUnknownWidth yes

    #
    # Check whether an up- or down-arrow is to be displayed
    #
    set font $configVals(-font)
    if {$col >= 0 && $configVals(-showarrow) && $data($col-showarrow)} {
	#
	# Configure the canvas and redraw the arrows
	#
	set data(arrowCol) $col
	configCanvas $win
	redrawArrows $win

	#
	# Make sure the arrow will fit into the static-width
	# column by enlarging the latter if necessary
	#
	set idx [expr {2*$col}]
	set pixels [lindex $data(colList) $idx]
	if {$pixels != 0 && $pixels < $data(arrowSize)} {
	    set data(colList) \
		[lreplace $data(colList) $idx $idx $data(arrowSize)]
	    set idx [expr {3*$col}]
	    set configVals(-columns) \
		[lreplace $configVals(-columns) $idx $idx -$data(arrowSize)]
	}

	#
	# Adjust the columns; this will also place the canvas into the label
	#
	adjustColumns $win $col yes

	#
	# Define for the binding tag $canvas the binding scripts
	# obtained from those of $label by replacing the %x field with
	# [expr {%x + [winfo x $canvas] - [winfo x $label]}] and the %y
	# field with [expr {%y + [winfo y $canvas] - [winfo y $label]}]
	# 
	set label $data(hdrTxtFrLbl)$col
	foreach event [bind $label] {
	    set script [bind $label $event]
	    regsub -all %x $script {$tablelist::x} script
	    regsub -all %y $script {$tablelist::y} script

	    bind $canvas $event [format {
		set tablelist::x [expr {%%x + [winfo x %s] - [winfo x %s]}]
		set tablelist::y [expr {%%y + [winfo y %s] - [winfo y %s]}]
		%s
	    } $canvas $label $canvas $label $script]
	}
    }

    #
    # Delete the items from the body text widget and insert the sorted ones
    #
    $w delete 1.0 end
    set line 1
    foreach item $data(itemList) {
	if {$line != 1} {
	    $w insert $line.0 \n
	}

	#
	# Insert the item; preserve any tags associated with its elements
	#
	set key [lindex $item end]
	set col 0
	foreach {pixels alignment} $data(colList) {
	    if {$data($col-hide)} {
		incr col
		continue
	    }

	    #
	    # Adjust the cell text and the image width
	    #
	    set text [lindex $item $col]
	    set name $key-$col-image
	    if {[info exists data($name)]} {
		set image $data($name)
		set imageWidth [image width $image]
	    } else {
		set image ""
		set imageWidth 0
	    }
	    if {$pixels != 0} {			;# convention: static width
		incr pixels $data($col-delta)
	    }
	    adjustElem $win text imageWidth $font $pixels $alignment

	    #
	    # Insert the text and the image
	    #
	    set tagNames {}
	    foreach opt {-background -foreground} {
		foreach tag [list $col$opt $key$opt $key-$col$opt] {
		    if {[info exists data($tag)]} {
			lappend tagNames $tag
		    }
		}
	    }
	    if {$imageWidth == 0} {
		$w insert $line.end \t$text\t $tagNames
	    } else {
		$w insert $line.end \t\t $tagNames
		insertElem $w $line.end-1c $text $image $imageWidth $alignment
	    }

	    incr col
	}

	incr line
    }

    #
    # Select the items that were selected before
    #
    foreach key $selKeys {
	set idx [lsearch $data(itemList) "* $key"]
	selectionSubCmd $win set $idx $idx 2
    }

    #
    # Disable the body text widget if it was disabled before
    #
    if {[string compare $configVals(-state) disabled] == 0} {
	$w tag add disabled 1.0 end
	$w tag configure select -borderwidth 0
    }

    #
    # Replace the contents in the list variable if present
    #
    if {[string compare $configVals(-listvariable) ""] != 0} {
	set traceCmd [list tablelist::traceProc $win]
	trace vdelete ::$configVals(-listvariable) wu $traceCmd
	set lastCol [expr {$data(colCount) - 1}]
	upvar #0 $configVals(-listvariable) var
	set var {}
	foreach item $data(itemList) {
	    lappend var [lrange $item 0 $lastCol]
	}
	trace variable ::$configVals(-listvariable) wu $traceCmd
    }

    return ""
}

#------------------------------------------------------------------------------
# tablelist::xviewSubCmd
#
# This procedure is invoked to process the tablelist xview subcommand.
#------------------------------------------------------------------------------
proc tablelist::xviewSubCmd {win argList} {
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

    set w $data(body)
    switch [llength $argList] {
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
	    set index [rowIndex $win [lindex $argList 0] no]
	    set fraction [lindex [$w xview] 0]
	    $w yview [expr {$index + 1}].0
	    $w xview moveto $fraction
	    return ""
	}

	default {
	    #
	    # Command: $win yview moveto fraction
	    #	       $win yview scroll number what
	    #
	    # We must check the word following "yview" before passing
	    # the arguments to the body text widget because of some
	    # additional options supported by its yview subcommand.
	    #
	    set opt [lindex $argList 0]
	    if {[string first $opt moveto] == 0 ||
		[string first $opt scroll] == 0} {
		return [eval [list $w yview] $argList]
	    } else {
		return -code error \
		       "unknown option \"$opt\": must be moveto or scroll"
	    }
	}
    }
}

#
# Private callback procedures
# ===========================
#

#------------------------------------------------------------------------------
# tablelist::focusCtrl
#
# Determines whether the body text child w of a tablelist widget will receive
# the input focus during keyboard traversal.
#------------------------------------------------------------------------------
proc tablelist::focusCtrl w {
    set win [winfo parent $w]
    upvar ::tablelist::ns${win}::configVals configVals
    set val $configVals(-takefocus)

    switch -- $val {
	0 - 1 - "" {
	    return $val
	}

	default {
	    return [uplevel #0 $val [list $win]]
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::fetchSelection
#
# This procedure is invoked when the PRIMARY selection is owned by the
# tablelist widget win and someone attempts to retrieve it as a STRING.  It
# returns part or all of the selection, as given by offset and maxChars.  The
# string which is to be (partially) returned is built by joining all of the
# elements of the selected rows together with tabs and the rows themselves with
# newlines.
#------------------------------------------------------------------------------
proc tablelist::fetchSelection {win offset maxChars} {
    upvar ::tablelist::ns${win}::configVals configVals

    if {!$configVals(-exportselection)} {
	return ""
    }

    set selection ""
    set gotItem false
    foreach idx [curselectionSubCmd $win] {
	if {$gotItem} {
	    append selection \n
	}

	append selection [join [getSubCmd $win $idx $idx 1] \t]
	set gotItem true
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
    upvar ::tablelist::ns${win}::configVals configVals

    if {$configVals(-exportselection)} {
	selectionSubCmd $win clear 0 [expr {$data(itemCount) - 1}] 2
    }
}

#
# Private procedures used in bindings
# ===================================
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
}

#------------------------------------------------------------------------------
# tablelist::cleanup
#
# This procedure is invoked when the tablelist widget win is destroyed.  It
# executes some cleanup operations.
#------------------------------------------------------------------------------
proc tablelist::cleanup win {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    #
    # Cancel the execution of all delayed stretchColumns,
    # synchronize, redisplay, and redisplayCol commands
    #
    if {[info exists data(stretchId)]} {
	after cancel $data(stretchId)
    }
    if {[info exists data(syncId)]} {
	after cancel $data(syncId)
    }
    foreach name [array names data *redispId] {
	after cancel $data($name)
    }

    #
    # If there is a list variable associated with the
    # widget then remove the trace set on this variable
    #
    if {[string compare $configVals(-listvariable) ""] != 0} {
	trace vdelete ::$configVals(-listvariable) wu \
	      [list tablelist::traceProc $win]
    }

    namespace delete ::tablelist::ns$win
    catch {rename ::$win ""}
}

#------------------------------------------------------------------------------
# tablelist::tablelistAutoScan
#
# This is a modified version of the Tk library procedure tkListboxAutoScan.  It
# is invoked when the mouse leaves the body text child of a tablelist widget.
# It scrolls the child and reschedules itself as an after command so that the
# child continues to scroll until the mouse moves back into the window or the
# mouse button is released.
#------------------------------------------------------------------------------
proc tablelist::tablelistAutoScan win {
    global tkPriv
    if {![winfo exists $win]} {
	return ""
    }
    set x $tkPriv(x)
    set y $tkPriv(y)

    set w [::$win bodypath]
    set _x [expr {$x - [winfo x $w]}]
    set _y [expr {$y - [winfo y $w]}]

    if {$_y >= [winfo height $w]} {
	::$win yview scroll 1 units
    } elseif {$_y < 0} {
	::$win yview scroll -1 units
    } elseif {$_x >= [winfo width $w]} {
	::$win xview scroll 2 units
    } elseif {$_x < 0} {
	::$win xview scroll -2 units
    } else {
	return ""
    }

    tkListboxMotion $win [::$win index @$x,$y]
    set tkPriv(afterId) [after 50 [list tablelist::tablelistAutoScan $win]]
}

#------------------------------------------------------------------------------
# tablelist::labelEnter
#
# This procedure is invoked when the mouse pointer enters the col'th label of
# the tablelist widget win, or is moving within that label.  It updates the
# cursor, depending on whether the pointer is on the right border of the label
# or not.
#------------------------------------------------------------------------------
proc tablelist::labelEnter {win col x} {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    set w $data(hdrTxtFrLbl)$col
    $w configure -cursor $configVals(-cursor)
    if {[string compare $configVals(-state) disabled] == 0} {
	return ""
    }

    if {$configVals(-resizablecolumns) && $data($col-resizable) &&
	$x >= [expr {[winfo width $w] - [$w cget -borderwidth] - 4}]} {
	$w configure -cursor $configVals(-resizecursor)
    }
}

#------------------------------------------------------------------------------
# tablelist::labelB1Down
#
# This procedure is invoked when mouse button 1 is pressed in the col'th label
# of the tablelist widget win.  If the pointer is on the right border of the
# label then the procedure records its x-coordinate relative to the label, the
# width of the column, and some other data needed later.  Otherwise it saves
# the label's relief so it can be restored later, and changes the relief to
# sunken.
#------------------------------------------------------------------------------
proc tablelist::labelB1Down {win col x} {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    if {[string compare $configVals(-state) disabled] == 0 ||
	[info exists data(x)]} {		;# resize operation in progress
	return ""
    }

    set w $data(hdrTxtFrLbl)$col
    set labelWidth [winfo width $w]

    if {$configVals(-resizablecolumns) && $data($col-resizable) &&
	$x >= [expr {$labelWidth - [$w cget -borderwidth] - 4}]} {
	set data(clickedLblCol) $col
	set data(x) $x

	set font $configVals(-font)
	set data(curColWidth) \
	    [expr {$labelWidth - 2*[font measure $font -displayof $win 0]}]
	set data(configColWidth) [lindex $configVals(-columns) [expr {3*$col}]]

	if {$col == $data(arrowCol)} {
	    set data(minColWidth) $data(arrowSize)
	} else {
	    set data(minColWidth) 1
	}

	set topWin [winfo toplevel $win]
	set data(topBindEsc) [bind $topWin <Escape>]
	bind $topWin <Escape> [list tablelist::escape $win]
    } elseif {([info exists data($col-labelcommand)] &&
	       [string compare $data($col-labelcommand) ""] != 0) ||
	      [string compare $configVals(-labelcommand) ""] != 0} {
	set data(clickedLblCol) $col
	set data(inClickedLabel) yes
	set data(relief) [$w cget -relief]
	$w configure -relief sunken
    }
}

#------------------------------------------------------------------------------
# tablelist::labelB1Motion
#
# This procedure is invoked to process mouse motion events in a label of the
# tablelist widget win while button 1 is down.  If this event occured during a
# column resize operation then the procedure computes the difference between
# the pointer's new x-coordinate relative to that label and the one recorded by
# the last invocation of labelB1Down, and adjusts the width of the
# corresponding column accordingly.  This column is only redisplayed in the
# currently visible rows of the widget.
#------------------------------------------------------------------------------
proc tablelist::labelB1Motion {win x y} {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    if {$data(clickedLblCol) < 0} {
	return ""
    }

    set col $data(clickedLblCol)

    if {[info exists data(x)]} {		;# resize operation in progress
	set width [expr {$data(curColWidth) + $x - $data(x)}]
	if {$width >= $data(minColWidth)} {
	    set idx [expr {3*$col}]
	    set configVals(-columns) \
		[lreplace $configVals(-columns) $idx $idx -$width]
	    set idx [expr {2*$col}]
	    set data(colList) [lreplace $data(colList) $idx $idx $width]
	    set data($col-delta) 0
	    adjustColumns $win -1 no
	    update idletasks
	    redisplayCol $win $col [rowIndex $win @0,0 no] \
				   [rowIndex $win @0,[winfo height $win] no]
	}
    } else {
	#
	# The following code is needed because the event can also
	# occur in the canvas displaying an up- or down-arrow
	#
	set w $data(hdrTxtFrLbl)$col
	if {$x >= 0 && $x < [winfo width $w] &&
	    $y >= 0 && $y < [winfo height $w]} {
	    $w configure -relief sunken
	    set data(inClickedLabel) yes
	} else {
	    $w configure -relief $data(relief)
	    set data(inClickedLabel) no
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::labelB1Enter
#
# This procedure is invoked when the mouse pointer enters a label of the
# tablelist widget win while mouse button 1 is down.  If no label was
# previously clicked then nothing happens.  Otherwise, if this event occured
# during a column resize operation then the procedure updates the mouse cursor
# accordingly.  Otherwise it changes the label's relief to sunken.
#------------------------------------------------------------------------------
proc tablelist::labelB1Enter win {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    if {$data(clickedLblCol) < 0} {
	return ""
    }

    set w $data(hdrTxtFrLbl)$data(clickedLblCol)
    $w configure -cursor $configVals(-cursor)

    if {[info exists data(x)]} {		;# resize operation in progress
	$w configure -cursor $configVals(-resizecursor)
	return ""
    }

    $w configure -relief sunken
    set data(inClickedLabel) yes
}

#------------------------------------------------------------------------------
# tablelist::labelB1Leave
#
# This procedure is invoked when the mouse pointer leaves a label of the
# tablelist widget win while mouse button 1 is down.  If no label was
# previously clicked then nothing happens.  Otherwise, if no column resize
# operation is in progress then the procedure restores the label's relief.
#------------------------------------------------------------------------------
proc tablelist::labelB1Leave {win x y} {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    if {$data(clickedLblCol) < 0 ||
	[info exists data(x)]} {		;# resize operation in progress
	return ""
    }

    set w $data(hdrTxtFrLbl)$data(clickedLblCol)

    #
    # The following code is needed because the event can also
    # occur in the canvas displaying an up- or down-arrow
    #
    if {$x >= 0 && $x < [winfo width $w] &&
	$y >= 0 && $y < [winfo height $w]} {
	return ""
    }

    $w configure -relief $data(relief)
    set data(inClickedLabel) no
}

#------------------------------------------------------------------------------
# tablelist::labelB1Up
#
# This procedure is invoked when mouse button 1 is released in a previously
# clicked label of the tablelist widget win.  If this event occured during a
# column resize operation then the procedure redisplays that column.  Otherwise
# it restores the label's relief and invokes the command specified with the
# -labelcommand configuration option, passing to it the widget name and the
# column number as arguments.
#------------------------------------------------------------------------------
proc tablelist::labelB1Up win {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    if {$data(clickedLblCol) < 0} {
	return ""
    }

    set col $data(clickedLblCol)
    set w $data(hdrTxtFrLbl)$col

    if {[info exists data(x)]} {		;# resize operation in progress
	$w configure -cursor $configVals(-cursor)
	bind [winfo toplevel $win] <Escape> $data(topBindEsc)
	redisplayColWhenIdle $win $col
	if {$configVals(-width) <= 0} {
	    $data(hdr) configure -width $data(hdrPixels)
	}
	stretchColumnsWhenIdle $win
	unset data(x)
    } elseif {$data(inClickedLabel)} {
	$w configure -relief $data(relief)
	if {[info exists data($col-labelcommand)] &&
	    [string compare $data($col-labelcommand) ""] != 0} {
	    uplevel #0 $data($col-labelcommand) [list $win $col]
	} elseif {[string compare $configVals(-labelcommand) ""] != 0} {
	    uplevel #0 $configVals(-labelcommand) [list $win $col]
	}
    }

    set data(clickedLblCol) -1
}

#------------------------------------------------------------------------------
# tablelist::labelB3Down
#
# This procedure is invoked when mouse button 3 is pressed in the col'th label
# of the tablelist widget win.  It configures the width of the given column to
# be just large enough to hold all the elements (including the label).
#------------------------------------------------------------------------------
proc tablelist::labelB3Down {win col} {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    if {[string compare $configVals(-state) normal] == 0 &&
	$configVals(-resizablecolumns) && $data($col-resizable)} {
	doColConfig $col $win -width 0
    }
}

#------------------------------------------------------------------------------
# tablelist::escape
#
# This procedure is invoked to process <Escape> events in the top-level window
# containing the tablelist widget win during a column resize operation.  The
# procedure cancels this operation and restores the initial width of the
# respective column.
#------------------------------------------------------------------------------
proc tablelist::escape win {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    set col $data(clickedLblCol)
    $data(hdrTxtFrLbl)$col configure -cursor $configVals(-cursor)
    update idletasks
    bind [winfo toplevel $win] <Escape> $data(topBindEsc)
    set idx [expr {3*$col}]
    setupColumns $win [lreplace $configVals(-columns) $idx $idx \
				$data(configColWidth)] no
    adjustColumns $win $col yes
    redisplayCol $win $col [rowIndex $win @0,0 no] \
			   [rowIndex $win @0,[winfo height $win] no]

    unset data(x)
    set data(clickedLblCol) -1
}

#
# Private helper procedures
# =========================
#

#------------------------------------------------------------------------------
# tablelist::configCanvas
#
# Sets the background, width, and height of the canvas displaying an up- or
# down-arrow, and saves its size in data(arrowSize).
#------------------------------------------------------------------------------
proc tablelist::configCanvas win {
    upvar ::tablelist::ns${win}::data data

    set w $data(hdrTxtFrLbl)$data(arrowCol)
    set labelFont [$w cget -font]
    if {[font metrics $labelFont -displayof $w -fixed]} {
	set spaces " "
    } else {
	set spaces "  "
    }

    set spacePixels [font measure $labelFont -displayof $w $spaces]
    if {[expr {$spacePixels % 2}] == 0} {
	incr spacePixels
    }

    $data(hdrTxtFrCanv) configure -background [$w cget -background] \
				  -height $spacePixels -width $spacePixels
    set data(arrowSize) $spacePixels
}

#------------------------------------------------------------------------------
# tablelist::redrawArrows
#
# Modifies the coordinates of the two filled polygons contained in the canvas.
#------------------------------------------------------------------------------
proc tablelist::redrawArrows win {
    upvar ::tablelist::ns${win}::data data
    upvar ::tablelist::ns${win}::configVals configVals

    switch $configVals(-incrarrowtype) {
	up {
	    switch $data(sortOrder) {
		increasing {
		    set arrowType up
		}
		decreasing {
		    set arrowType down
		}
	    }
	}

	down {
	    switch $data(sortOrder) {
		increasing {
		    set arrowType down
		}
		decreasing {
		    set arrowType up
		}
	    }
	}
    }

    set w $data(hdrTxtFrCanv)
    set last [expr {$data(arrowSize) - 1}]
    set half [expr {$data(arrowSize) / 2}]

    switch $arrowType {
	up {
	    $w coords normalArrow   0 $last $half 0 $last $last
	    $w coords disabledArrow 0 $last $half 0 $last $last
	}

	down {
	    $w coords normalArrow   0 0 $half $last $last 0
	    $w coords disabledArrow 0 0 $half $last $last 0
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::rowIndex
#
# Checks the row index idx and returns either its value or an error.  endIsSize
# must be a boolean value: if true, end refers to the number of items in the
# tablelist, i.e., to the element just after the last one; if false, end refers
# to 1 less than the number of items, i.e., to the last element in the
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
	    return [expr {$data(itemCount) - 1}]
	}
    } elseif {[string compare [string range $idx 0 0] @] == 0} {
	if {[catch {$data(body) index $idx}] == 0} {
	    if {$data(itemCount) == 0} {
		return -1
	    } else {
		scan $idx @%d,%d x y
		incr x -[winfo x $data(body)]
		incr y -[winfo y $data(body)]
		set bodyIdx [$data(body) index @$x,$y]
		return [expr {int($bodyIdx) - 1}]
	    }
	} else {
	    return -code error \
		   "bad row index \"$idx\": must be active, anchor,\
		    end, @x,y, or a number"
	}
    } elseif {[catch {format %d $idx} index] == 0} {
	return $index
    } else {
	return -code error \
	       "bad row index \"$idx\": must be active, anchor,\
	        end, @x,y, or a number"
    }
}

#------------------------------------------------------------------------------
# tablelist::colIndex
#
# Checks the column index idx and returns either its value or an error.
#------------------------------------------------------------------------------
proc tablelist::colIndex {win idx} {
    upvar ::tablelist::ns${win}::data data

    if {[string first $idx end] == 0} {
	return [expr {$data(colCount) - 1}]
    } elseif {[string compare [string range $idx 0 0] @] == 0} {
	if {[catch {$data(body) index $idx}] == 0} {
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
		if {$x >= $wX && $x < [expr {$wX + [winfo width $w]}]} {
		    return $col
		}
	    }
	    return $lastVisibleCol
	} else {
	    return -code error \
		   "bad column index \"$idx\": must be end, @x,y, or a number"
	}
    } elseif {[catch {format %d $idx} index] == 0} {
	return $index
    } else {
	return -code error \
	       "bad column index \"$idx\": must be end, @x,y, or a number"
    }
}

#------------------------------------------------------------------------------
# tablelist::cellIndex
#
# Checks the cell index idx and returns either its value (in the form row,col)
# or an error.
#------------------------------------------------------------------------------
proc tablelist::cellIndex {win idx} {
    if {[string first $idx end] == 0} {
	return [rowIndex $win $idx no],[colIndex $win $idx]
    } elseif {[string compare [string range $idx 0 0] @] == 0} {
	if {[catch {rowIndex $win $idx no} row] == 0 &&
	    [catch {colIndex $win $idx   } col] == 0} {
	    return $row,$col
	} else {
	    return -code error \
		   "bad cell index \"$idx\": must be end, @x,y, or\
		    row,col, where row must be active, anchor, end,\
		    or a number, and col must be end or a number"
	}
    } else {
	set lst [split $idx ,]
	if {[llength $lst] == 2 &&
	    [catch {rowIndex $win [lindex $lst 0] no} row] == 0 &&
	    [catch {colIndex $win [lindex $lst 1]   } col] == 0} {
	    return $row,$col
	} else {
	    return -code error \
		   "bad cell index \"$idx\": must be end, @x,y, or\
		    row,col, where row must be active, anchor, end,\
		    or a number, and col must be end or a number"
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::findCellTabs
#
# Searches for the tab characters within the col'th cell in the body text child
# of the tablelist widget win, starting at startIdx.  Assigns the index of the
# first tab to $idx1Name and the index of the second tab to $idx2Name.
#------------------------------------------------------------------------------
proc tablelist::findCellTabs {win startIdx col idx1Name idx2Name} {
    upvar ::tablelist::ns${win}::data data
    upvar $idx1Name idx1 $idx2Name idx2

    set w $data(body)
    set endIdx "$startIdx lineend"
    set idx1 $startIdx
    for {set n 0} {$n < $col} {incr n} {
	if {$data($n-hide)} {
	    continue
	}
	set idx1 [$w search \t $idx1+1c $endIdx]+1c
    }
    set idx2 [$w search \t $idx1+1c $endIdx]

    return ""
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
# Invokes strRange with the given arguments and returns a string with up to
# three trailing (for alignment = left or center) or leading (for alignment =
# right) dots built from its result.
#------------------------------------------------------------------------------
proc tablelist::strRangeExt {win str font pixels alignment} {
    set len [string length $str]
    set subStr [strRange $win $str $font $pixels $alignment]
    if {$pixels < 0 || [string length $subStr] == $len} {
	return $subStr
    }

    if {[string compare $alignment right] == 0} {
	set first 0
	set extSubStr ...$subStr
	while {[font measure $font -displayof $win $extSubStr] > $pixels} {
	    if {$first < $len} {
		incr first
		set subStr [string range $subStr $first end]
		set extSubStr ...$subStr
	    } else {
		set extSubStr [string range $extSubStr 1 end]
	    }
	}
    } else {
	set last [expr {$len - 1}]
	set extSubStr $subStr...
	while {[font measure $font -displayof $win $extSubStr] > $pixels} {
	    if {$last >= 0} {
		incr last -1
		set subStr [string range $subStr 0 $last]
		set extSubStr $subStr...
	    } else {
		set extSubStr [string range $extSubStr 1 end]
	    }
	}
    }

    return $extSubStr
}

#------------------------------------------------------------------------------
# tablelist::adjustElem
#
# Prepares the text specified by $textName and the image width specified by
# $imageWidthName for insertion into a cell of the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::adjustElem {win textName imageWidthName font pixels alignment} {
    upvar $textName text
    upvar $imageWidthName imageWidth

    if {$pixels == 0} {				;# convention: dynamic width
	if {$imageWidth != 0 && [string compare $text ""] != 0} {
	    if {[string compare $alignment right] == 0} {
		set text "$text "
	    } else {
		set text " $text"
	    }
	}
    } elseif {$imageWidth == 0} {		;# no image
	set text [strRangeExt $win $text $font $pixels $alignment]
    } elseif {[string compare $text ""] == 0} {	;# image w/o text
	if {$imageWidth > $pixels} {
	    set imageWidth 0			;# can't display the image
	}
    } else {					;# both image and text
	set gap [font measure $font -displayof $win " "]
	if {$imageWidth + $gap <= $pixels} {
	    incr pixels -[expr {$imageWidth + $gap}]
	    set text [strRangeExt $win $text $font $pixels $alignment]
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
