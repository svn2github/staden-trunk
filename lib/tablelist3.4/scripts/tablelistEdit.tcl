#==============================================================================
# Contains the implementation of interactive cell editing in tablelist widgets.
#
# Stucture of the module:
#   - Namespace initialization
#   - Public procedures related to interactive cell editing
#   - Private procedures implementing the interactive cell editing
#   - Private procedures used in bindings related to interactive cell editing
#
# Copyright (c) 2003-2004  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

#
# Namespace initialization
# ========================
#

namespace eval tablelist {
    #
    # Define some bindings for the binding tag TablelistEdit
    #
    bind TablelistEdit <Button-1>         { focus %W }
    bind TablelistEdit <Control-i>        { tablelist::insertChar    %W \t }
    bind TablelistEdit <Control-j>        { tablelist::insertChar    %W \n }
    bind TablelistEdit <Control-Return>   { tablelist::insertChar    %W \n }
    bind TablelistEdit <Control-KP_Enter> { tablelist::insertChar    %W \n }
    bind TablelistEdit <Escape>	          { tablelist::cancelEditing %W }
    bind TablelistEdit <Return>	          { tablelist::finishEditing %W }
    bind TablelistEdit <KP_Enter>         { tablelist::finishEditing %W }
    bind TablelistEdit <Tab>              { tablelist::goToNextCell  %W }
    bind TablelistEdit <Shift-Tab>        { tablelist::goToPrevCell  %W }
    bind TablelistEdit <<PrevWindow>>     { tablelist::goToPrevCell  %W }
    bind TablelistEdit <Alt-Left>         { tablelist::goLeft        %W }
    bind TablelistEdit <Meta-Left>        { tablelist::goLeft        %W }
    bind TablelistEdit <Alt-Right>        { tablelist::goRight       %W }
    bind TablelistEdit <Meta-Right>       { tablelist::goRight       %W }
    bind TablelistEdit <Alt-Up>           { tablelist::goUp          %W }
    bind TablelistEdit <Meta-Up>          { tablelist::goUp          %W }
    bind TablelistEdit <Alt-Down>         { tablelist::goDown        %W }
    bind TablelistEdit <Meta-Down>        { tablelist::goDown        %W }
    bind TablelistEdit <Alt-Prior>        { tablelist::goToPrevPage  %W }
    bind TablelistEdit <Meta-Prior>       { tablelist::goToPrevPage  %W }
    bind TablelistEdit <Alt-Next>         { tablelist::goToNextPage  %W }
    bind TablelistEdit <Meta-Next>        { tablelist::goToNextPage  %W }
    bind TablelistEdit <Control-Home>     { tablelist::goToNextCell  %W 0 -1 }
    bind TablelistEdit <Control-End>      { tablelist::goToPrevCell  %W 0  0 }
    bind TablelistEdit <Left> {
	if {![tablelist::isKeyReserved %W %K]} {
	    tablelist::goLeft %W
	}
    }
    bind TablelistEdit <Right> {
	if {![tablelist::isKeyReserved %W %K]} {
	    tablelist::goRight %W
	}
    }
    bind TablelistEdit <Up> {
	if {![tablelist::isKeyReserved %W %K]} {
	    tablelist::goUp %W
	}
    }
    bind TablelistEdit <Down> {
	if {![tablelist::isKeyReserved %W %K]} {
	    tablelist::goDown %W
	}
    }
    bind TablelistEdit <Prior> {
	if {![tablelist::isKeyReserved %W %K]} {
	    tablelist::goToPrevPage %W
	}
    }
    bind TablelistEdit <Next> {
	if {![tablelist::isKeyReserved %W %K]} {
	    tablelist::goToNextPage %W
	}
    }
    bind TablelistEdit <Control-Tab> {
	mwutil::generateEvent %W Tablelist <Tab>
    }
    bind TablelistEdit <Meta-Tab> {
	mwutil::generateEvent %W Tablelist <Tab>
    }
    bind TablelistEdit <Control-Shift-Tab> {
	mwutil::generateEvent %W Tablelist <Shift-Tab>
    }
    bind TablelistEdit <Meta-Shift-Tab> {
	mwutil::generateEvent %W Tablelist <Shift-Tab>
    }
    bind TablelistEdit <FocusIn> {
	set tablelist::ns[tablelist::parseEditWinPath %W]::data(editFocus) %W
    }

    #
    # Define some emacs-like key bindings for the binding tag TablelistEdit
    #
    bind TablelistEdit <Meta-b> {
	if {!$tk_strictMotif && ![tablelist::isKeyReserved %W Meta-b]} {
	    tablelist::goLeft %W
	}
    }
    bind TablelistEdit <Meta-f> {
	if {!$tk_strictMotif && ![tablelist::isKeyReserved %W Meta-f]} {
	    tablelist::goRight %W
	}
    }
    bind TablelistEdit <Control-p> {
	if {!$tk_strictMotif && ![tablelist::isKeyReserved %W Control-p]} {
	    tablelist::goUp %W
	}
    }
    bind TablelistEdit <Control-n> {
	if {!$tk_strictMotif && ![tablelist::isKeyReserved %W Control-n]} {
	    tablelist::goDown %W
	}
    }
    bind TablelistEdit <Meta-less> {
	if {!$tk_strictMotif && ![tablelist::isKeyReserved %W Meta-less]} {
	    tablelist::goToNextCell %W 0 -1
	}
    }
    bind TablelistEdit <Meta-greater> {
	if {!$tk_strictMotif && ![tablelist::isKeyReserved %W Meta-greater]} {
	    tablelist::goToPrevCell %W 0 0
	}
    }

    #
    # Define some bindings for the binding tag TablelistEdit that
    # propagate the mousewheel events to the tablelist's body
    #
    if {[string compare $::tk_patchLevel 8.0.4] >= 0} {
	bind TablelistEdit <MouseWheel> {
	    if {![tablelist::isComboTopMapped %W]} {
		tablelist::genMouseWheelEvent \
			[[tablelist::parseEditWinPath %W] bodypath] %D
	    }
	}
    }
    bind TablelistEdit <Button-4> {
	if {![tablelist::isComboTopMapped %W]} {
	    event generate [[tablelist::parseEditWinPath %W] bodypath] \
			<Button-4>
	}
    }
    bind TablelistEdit <Button-5> {
	if {![tablelist::isComboTopMapped %W]} {
	    event generate [[tablelist::parseEditWinPath %W] bodypath] \
			<Button-5>
	}
    }

    #
    # Create two images, needed for checkbuttons
    #
    proc createImages {} {
	variable winSys
	if {[catch {tk windowingsystem} winSys] != 0} {
	    switch $::tcl_platform(platform) {
		unix		{ set winSys x11 }
		windows		{ set winSys win32 }
		macintosh	{ set winSys classic }
	    }
	}

	variable library
	variable checkedImg [image create bitmap -file \
	    [file join $library images checked_$winSys.xbm]]
	variable uncheckedImg [image create bitmap -file \
	    [file join $library images unchecked_$winSys.xbm]]
    }
    createImages 

    #
    # Register the Tk core widgets entry, checkbutton,
    # and spinbox for interactive cell editing
    #
    proc addTkCoreWidgets {} {
	set name entry
	array set ::tablelist::editWin [list \
	    $name-creationCmd	"$name %W -width 0" \
	    $name-putTextCmd	"%W delete 0 end; %W insert 0 %T" \
	    $name-getTextCmd	"%W get" \
	    $name-putListCmd	"" \
	    $name-getListCmd	"" \
	    $name-selectCmd	"" \
	    $name-invokeCmd	"" \
	    $name-fontOpt	-font \
	    $name-useFormat	1 \
	    $name-useReqWidth	0 \
	    $name-usePadX	0 \
	    $name-focusWin	%W \
	    $name-reservedKeys	{Left Right} \
	]

	set name checkbutton
	array set ::tablelist::editWin [list \
	    $name-creationCmd	"createCheckbutton %W" \
	    $name-putTextCmd	{set [%W cget -variable] %T} \
	    $name-getTextCmd	{set [%W cget -variable]} \
	    $name-putListCmd	"" \
	    $name-getListCmd	"" \
	    $name-selectCmd	"" \
	    $name-invokeCmd	"%W invoke" \
	    $name-fontOpt	-font \
	    $name-useFormat	0 \
	    $name-useReqWidth	1 \
	    $name-usePadX	0 \
	    $name-focusWin	%W \
	    $name-reservedKeys	{} \
	]

	if {$::tk_version < 8.4} {
	    return ""
	}

	set name spinbox
	array set ::tablelist::editWin [list \
	    $name-creationCmd	"$name %W -width 0" \
	    $name-putTextCmd	"%W delete 0 end; %W insert 0 %T" \
	    $name-getTextCmd	"%W get" \
	    $name-putListCmd	"" \
	    $name-getListCmd	"" \
	    $name-selectCmd	"" \
	    $name-invokeCmd	"" \
	    $name-fontOpt	-font \
	    $name-useFormat	1 \
	    $name-useReqWidth	0 \
	    $name-usePadX	1 \
	    $name-focusWin	%W \
	    $name-reservedKeys	{Left Right Up Down} \
	]
    }
    addTkCoreWidgets 
}

#
# Public procedures related to interactive cell editing
# =====================================================
#

#------------------------------------------------------------------------------
# tablelist::addBWidgetEntry
#
# Registers the Entry widget from the BWidget package for interactive cell
# editing.
#------------------------------------------------------------------------------
proc tablelist::addBWidgetEntry {{name Entry}} {
    checkEditWinName $name

    array set ::tablelist::editWin [list \
	$name-creationCmd	"Entry %W -width 0" \
	$name-putTextCmd	"%W delete 0 end; %W insert 0 %T" \
	$name-getTextCmd	"%W get" \
	$name-putListCmd	"" \
	$name-getListCmd	"" \
	$name-selectCmd		"" \
	$name-invokeCmd		"" \
	$name-fontOpt		-font \
	$name-useFormat		1 \
	$name-useReqWidth	0 \
	$name-usePadX		0 \
	$name-focusWin		%W \
	$name-reservedKeys	{Left Right} \
    ]

    return $name
}

#------------------------------------------------------------------------------
# tablelist::addBWidgetSpinBox
#
# Registers the SpinBox widget from the BWidget package for interactive cell
# editing.
#------------------------------------------------------------------------------
proc tablelist::addBWidgetSpinBox {{name SpinBox}} {
    checkEditWinName $name

    array set ::tablelist::editWin [list \
	$name-creationCmd	"SpinBox %W -editable 1 -width 0" \
	$name-putTextCmd	"%W configure -text %T" \
	$name-getTextCmd	"%W cget -text" \
	$name-putListCmd	"" \
	$name-getListCmd	"" \
	$name-selectCmd		"" \
	$name-invokeCmd		"" \
	$name-fontOpt		-font \
	$name-useFormat		1 \
	$name-useReqWidth	0 \
	$name-usePadX		1 \
	$name-focusWin		%W.e \
	$name-reservedKeys	{Left Right Up Down Prior Next} \
    ]

    return $name
}

#------------------------------------------------------------------------------
# tablelist::addBWidgetComboBox
#
# Registers the ComboBox widget from the BWidget package for interactive cell
# editing.
#------------------------------------------------------------------------------
proc tablelist::addBWidgetComboBox {{name ComboBox}} {
    checkEditWinName $name

    array set ::tablelist::editWin [list \
	$name-creationCmd	"ComboBox %W -editable 1 -width 0" \
	$name-putTextCmd	"%W configure -text %T" \
	$name-getTextCmd	"%W cget -text" \
	$name-putListCmd	"" \
	$name-getListCmd	"" \
	$name-selectCmd		"" \
	$name-invokeCmd		"%W.a invoke" \
	$name-fontOpt		-font \
	$name-useFormat		1 \
	$name-useReqWidth	0 \
	$name-usePadX		1 \
	$name-focusWin		%W.e \
	$name-reservedKeys	{Left Right Up Down} \
    ]

    return $name
}

#------------------------------------------------------------------------------
# tablelist::addIncrEntryfield
#
# Registers the entryfield widget from the Iwidgets package for interactive
# cell editing.
#------------------------------------------------------------------------------
proc tablelist::addIncrEntryfield {{name entryfield}} {
    checkEditWinName $name

    array set ::tablelist::editWin [list \
	$name-creationCmd	"iwidgets::entryfield %W -width 0" \
	$name-putTextCmd	"%W clear; %W insert 0 %T" \
	$name-getTextCmd	"%W get" \
	$name-putListCmd	"" \
	$name-getListCmd	"" \
	$name-selectCmd		"" \
	$name-invokeCmd		"" \
	$name-fontOpt		-textfont \
	$name-useFormat		1 \
	$name-useReqWidth	0 \
	$name-usePadX		0 \
	$name-focusWin		{[%W component entry]} \
	$name-reservedKeys	{Left Right} \
    ]

    return $name
}

#------------------------------------------------------------------------------
# tablelist::addIncrDateTimeWidget
#
# Registers the datefield, dateentry, timefield, or timeentry widget from the
# Iwidgets package, with or without the -clicks option for its get subcommand,
# for interactive cell editing.
#------------------------------------------------------------------------------
proc tablelist::addIncrDateTimeWidget {widgetType args} {
    if {![regexp {^(datefield|dateentry|timefield|timeentry)$} $widgetType]} {
	return -code error \
	       "bad widget type \"$widgetType\": must be\
		datefield, dateentry, timefield, or timeentry"
    }

    switch [llength $args] {
	0 {
	    set useClicks 0
	    set name $widgetType
	}

	1 {
	    set arg [lindex $args 0]
	    if {[string compare $arg -seconds] == 0} {
		set useClicks 1
		set name $widgetType
	    } else {
		set useClicks 0
		set name $arg
	    }
	}

	2 {
	    set arg0 [lindex $args 0]
	    if {[string compare $arg0 -seconds] != 0} {
		return -code error "bad option \"$arg0\": must be -seconds"
	    }

	    set useClicks 1
	    set name [lindex $args 1]
	}

	default {
	    mwutil::wrongNumArgs "addIncrDateTimeWidget\
				  datefield|dateentry|timefield|timeentry\
				  ?-seconds? ?name?"
	}
    }
    checkEditWinName $name

    array set ::tablelist::editWin [list \
	$name-creationCmd	"iwidgets::$widgetType %W" \
	$name-putTextCmd	"%W show %T" \
	$name-getTextCmd	"%W get" \
	$name-putListCmd	"" \
	$name-getListCmd	"" \
	$name-selectCmd		"" \
	$name-invokeCmd		"" \
	$name-fontOpt		-textfont \
	$name-useReqWidth	1 \
	$name-usePadX		[string match *entry $widgetType] \
	$name-useFormat		1 \
	$name-reservedKeys	{Left Right Up Down} \
    ]
    if {$useClicks} {
	lappend ::tablelist::editWin($name-getTextCmd) -clicks
	set ::tablelist::editWin($name-useFormat) 0
    }
    if {[string match date* $widgetType]} {
	set ::tablelist::editWin($name-focusWin) {[%W component date]}
    } else {
	set ::tablelist::editWin($name-focusWin) {[%W component time]}
    }

    return $name
}

#------------------------------------------------------------------------------
# tablelist::addIncrSpinner
#
# Registers the spinner widget from the Iwidgets package for interactive cell
# editing.
#------------------------------------------------------------------------------
proc tablelist::addIncrSpinner {{name spinner}} {
    checkEditWinName $name

    array set ::tablelist::editWin [list \
	$name-creationCmd	"iwidgets::spinner %W -width 0" \
	$name-putTextCmd	"%W clear; %W insert 0 %T" \
	$name-getTextCmd	"%W get" \
	$name-putListCmd	"" \
	$name-getListCmd	"" \
	$name-selectCmd		"" \
	$name-invokeCmd		"" \
	$name-fontOpt		-textfont \
	$name-useFormat		1 \
	$name-useReqWidth	0 \
	$name-usePadX		1 \
	$name-focusWin		{[%W component entry]} \
	$name-reservedKeys	{Left Right} \
    ]

    return $name
}

#------------------------------------------------------------------------------
# tablelist::addIncrSpinint
#
# Registers the spinint widget from the Iwidgets package for interactive cell
# editing.
#------------------------------------------------------------------------------
proc tablelist::addIncrSpinint {{name spinint}} {
    checkEditWinName $name

    array set ::tablelist::editWin [list \
	$name-creationCmd	"iwidgets::spinint %W -width 0" \
	$name-putTextCmd	"%W clear; %W insert 0 %T" \
	$name-getTextCmd	"%W get" \
	$name-putListCmd	"" \
	$name-getListCmd	"" \
	$name-selectCmd		"" \
	$name-invokeCmd		"" \
	$name-fontOpt		-textfont \
	$name-useFormat		1 \
	$name-useReqWidth	0 \
	$name-usePadX		1 \
	$name-focusWin		{[%W component entry]} \
	$name-reservedKeys	{Left Right} \
    ]

    return $name
}

#------------------------------------------------------------------------------
# tablelist::addIncrCombobox
#
# Registers the combobox widget from the Iwidgets package for interactive cell
# editing.
#------------------------------------------------------------------------------
proc tablelist::addIncrCombobox {{name combobox}} {
    checkEditWinName $name

    array set ::tablelist::editWin [list \
	$name-creationCmd	"createIncrCombobox %W" \
	$name-putTextCmd	"%W clear entry; %W insert entry 0 %T" \
	$name-getTextCmd	"%W get" \
	$name-putListCmd	{eval [list %W insert list end] %L} \
	$name-getListCmd	"%W component list get 0 end" \
	$name-selectCmd		"%W selection set %I" \
	$name-invokeCmd		"%W invoke" \
	$name-fontOpt		-textfont \
	$name-useFormat		1 \
	$name-useReqWidth	0 \
	$name-usePadX		1 \
	$name-focusWin		{[%W component entry]} \
	$name-reservedKeys	{Left Right Up Down Control-p Control-n} \
    ]

    return $name
}

#------------------------------------------------------------------------------
# tablelist::addOakleyCombobox
#
# Registers Bryan Oakley's combobox widget for interactive cell editing.
#------------------------------------------------------------------------------
proc tablelist::addOakleyCombobox {{name combobox}} {
    checkEditWinName $name

    array set ::tablelist::editWin [list \
	$name-creationCmd	"createOakleyCombobox %W" \
	$name-putTextCmd	"%W delete 0 end; %W insert 0 %T" \
	$name-getTextCmd	"%W get" \
	$name-putListCmd	{eval [list %W list insert end] %L} \
	$name-getListCmd	"%W list get 0 end" \
	$name-selectCmd		"%W select %I" \
	$name-invokeCmd		"%W open" \
	$name-fontOpt		-font \
	$name-useFormat		1 \
	$name-useReqWidth	0 \
	$name-usePadX		1 \
	$name-focusWin		%W.entry \
	$name-reservedKeys	{Left Right Up Down Prior Next} \
    ]

    #
    # Patch the ::combobox::UpdateVisualAttributes procedure to make sure it
    # won't change the background and trough colors of the vertical scrollbar
    #
    catch {combobox::combobox}	;# enforces the evaluation of "combobox.tcl"
    if {[catch {rename ::combobox::UpdateVisualAttributes \
		::combobox::_UpdateVisualAttributes}] == 0} {
	proc ::combobox::UpdateVisualAttributes w {
	    set vsbBackground [$w.top.vsb cget -background]
	    set vsbTroughColor [$w.top.vsb cget -troughcolor]

	    ::combobox::_UpdateVisualAttributes $w

	    $w.top.vsb configure -background $vsbBackground
	    $w.top.vsb configure -troughcolor $vsbTroughColor
	}
    }

    return $name
}

#------------------------------------------------------------------------------
# tablelist::addDateMentry
#
# Registers the widget created by the mentry::dateMentry command from the
# Mentry package, with a given format and separator and with or without the
# "-gmt 1" option for the mentry::putClockVal and mentry::getClockVal commands,
# for interactive cell editing.
#------------------------------------------------------------------------------
proc tablelist::addDateMentry {fmt sep args} {
    #
    # Parse the fmt argument
    #
    if {![regexp {^([dmyY])([dmyY])([dmyY])$} $fmt dummy \
		 fields(0) fields(1) fields(2)]} {
	return -code error \
	       "bad format \"$fmt\": must be a string of length 3,\
		consisting of the letters d, m, and y or Y"
    }

    #
    # Check whether all the three date components are represented in fmt
    #
    for {set n 0} {$n < 3} {incr n} {
	set lfields($n) [string tolower $fields($n)]
    }
    if {[string compare $lfields(0) $lfields(1)] == 0 ||
	[string compare $lfields(0) $lfields(2)] == 0 ||
	[string compare $lfields(1) $lfields(2)] == 0} {
	return -code error \
	       "bad format \"$fmt\": must have unique components for the\
		day, month, and year"
    }

    #
    # Parse the remaining arguments (if any)
    #
    switch [llength $args] {
	0 {
	    set useGMT 0
	    set name dateMentry
	}

	1 {
	    set arg [lindex $args 0]
	    if {[string compare $arg -gmt] == 0} {
		set useGMT 1
		set name dateMentry
	    } else {
		set useGMT 0
		set name $arg
	    }
	}

	2 {
	    set arg0 [lindex $args 0]
	    if {[string compare $arg0 -gmt] != 0} {
		return -code error "bad option \"$arg0\": must be -gmt"
	    }

	    set useGMT 1
	    set name [lindex $args 1]
	}

	default {
	    mwutil::wrongNumArgs "addDateMentry format separator ?-gmt? ?name?"
	}
    }
    checkEditWinName $name

    array set ::tablelist::editWin [list \
	$name-creationCmd	[list mentry::dateMentry %W $fmt $sep] \
	$name-putTextCmd	"mentry::putClockVal %T %W -gmt $useGMT" \
	$name-getTextCmd	"mentry::getClockVal %W -gmt $useGMT" \
	$name-putListCmd	{eval [list %W put 0] %L} \
	$name-getListCmd	"%W getlist" \
	$name-selectCmd		"" \
	$name-invokeCmd		"" \
	$name-fontOpt		-font \
	$name-useFormat		0 \
	$name-useReqWidth	1 \
	$name-usePadX		1 \
	$name-focusWin		"" \
	$name-reservedKeys	{Left Right Up Down} \
    ]

    return $name
}

#------------------------------------------------------------------------------
# tablelist::addTimeMentry
#
# Registers the widget created by the mentry::timeMentry command from the
# Mentry package, with a given format and separator and with or without the
# "-gmt 1" option for the mentry::putClockVal and mentry::getClockVal commands,
# for interactive cell editing.
#------------------------------------------------------------------------------
proc tablelist::addTimeMentry {fmt sep args} {
    #
    # Parse the fmt argument
    #
    if {![regexp {^(H|I)(M)(S?)$} $fmt dummy fields(0) fields(1) fields(2)]} {
	return -code error \
	       "bad format \"$fmt\": must be a string of length 2 or 3\
		starting with H or I, followed by M and optionally by S"
    }

    #
    # Parse the remaining arguments (if any)
    #
    switch [llength $args] {
	0 {
	    set useGMT 0
	    set name timeMentry
	}

	1 {
	    set arg [lindex $args 0]
	    if {[string compare $arg -gmt] == 0} {
		set useGMT 1
		set name timeMentry
	    } else {
		set useGMT 0
		set name $arg
	    }
	}

	2 {
	    set arg0 [lindex $args 0]
	    if {[string compare $arg0 -gmt] != 0} {
		return -code error "bad option \"$arg0\": must be -gmt"
	    }

	    set useGMT 1
	    set name [lindex $args 1]
	}

	default {
	    mwutil::wrongNumArgs "addTimeMentry format separator ?-gmt? ?name?"
	}
    }
    checkEditWinName $name

    array set ::tablelist::editWin [list \
	$name-creationCmd	[list mentry::timeMentry %W $fmt $sep] \
	$name-putTextCmd	"mentry::putClockVal %T %W -gmt $useGMT" \
	$name-getTextCmd	"mentry::getClockVal %W -gmt $useGMT" \
	$name-putListCmd	{eval [list %W put 0] %L} \
	$name-getListCmd	"%W getlist" \
	$name-selectCmd		"" \
	$name-invokeCmd		"" \
	$name-fontOpt		-font \
	$name-useFormat		0 \
	$name-useReqWidth	1 \
	$name-usePadX		1 \
	$name-focusWin		"" \
	$name-reservedKeys	{Left Right Up Down} \
    ]

    return $name
}

#------------------------------------------------------------------------------
# tablelist::addFixedPointMentry
#
# Registers the widget created by the mentry::fixedPointMentry command from the
# Mentry package, with a given number of characters before and a given number
# of digits after the decimal point, with or without the -comma option, for
# interactive cell editing.
#------------------------------------------------------------------------------
proc tablelist::addFixedPointMentry {cnt1 cnt2 args} {
    #
    # Check the arguments cnt1 and cnt2
    #
    if {[catch {format %d $cnt1}] != 0 || $cnt1 <= 0} {
	return -code error "expected positive integer but got \"$cnt1\""
    }
    if {[catch {format %d $cnt2}] != 0 || $cnt2 <= 0} {
	return -code error "expected positive integer but got \"$cnt2\""
    }

    #
    # Parse the remaining arguments (if any)
    #
    switch [llength $args] {
	0 {
	    set useComma 0
	    set name fixedPointMentry_$cnt1.$cnt2
	}

	1 {
	    set arg [lindex $args 0]
	    if {[string compare $arg -comma] == 0} {
		set useComma 1
		set name fixedPointMentry_$cnt1,$cnt2
	    } else {
		set useComma 0
		set name $arg
	    }
	}

	2 {
	    set arg0 [lindex $args 0]
	    if {[string compare $arg0 -comma] != 0} {
		return -code error "bad option \"$arg0\": must be -comma"
	    }

	    set useComma 1
	    set name [lindex $args 1]
	}

	default {
	    mwutil::wrongNumArgs "addFixedPointMentry count1 count2\
				  ?-comma? ?name?"
	}
    }
    checkEditWinName $name

    array set ::tablelist::editWin [list \
	$name-creationCmd	[list mentry::fixedPointMentry %W $cnt1 $cnt2] \
	$name-putTextCmd	"mentry::putReal %T %W" \
	$name-getTextCmd	"mentry::getReal %W" \
	$name-putListCmd	{eval [list %W put 0] %L} \
	$name-getListCmd	"%W getlist" \
	$name-selectCmd		"" \
	$name-invokeCmd		"" \
	$name-fontOpt		-font \
	$name-useFormat		0 \
	$name-useReqWidth	1 \
	$name-usePadX		1 \
	$name-focusWin		"" \
	$name-reservedKeys	{Left Right} \
    ]
    if {$useComma} {
	lappend ::tablelist::editWin($name-creationCmd) -comma
    }

    return $name
}

#------------------------------------------------------------------------------
# tablelist::addIPAddrMentry
#
# Registers the widget created by the mentry::ipAddrMentry command from the
# Mentry package for interactive cell editing.
#------------------------------------------------------------------------------
proc tablelist::addIPAddrMentry {{name ipAddrMentry}} {
    checkEditWinName $name

    array set ::tablelist::editWin [list \
	$name-creationCmd	"mentry::ipAddrMentry %W" \
	$name-putTextCmd	"mentry::putIPAddr %T %W" \
	$name-getTextCmd	"mentry::getIPAddr %W" \
	$name-putListCmd	{eval [list %W put 0] %L} \
	$name-getListCmd	"%W getlist" \
	$name-selectCmd		"" \
	$name-invokeCmd		"" \
	$name-fontOpt		-font \
	$name-useFormat		0 \
	$name-useReqWidth	1 \
	$name-usePadX		1 \
	$name-focusWin		"" \
	$name-reservedKeys	{Left Right} \
    ]

    return $name
}

#
# Private procedures implementing the interactive cell editing
# ============================================================
#

#------------------------------------------------------------------------------
# tablelist::checkEditWinName
#
# Generates an error if the given edit window name is one of "entry",
# "spinbox", or "checkbutton".
#------------------------------------------------------------------------------
proc tablelist::checkEditWinName name {
    if {[regexp {^(entry|spinbox|checkbutton)$} $name]} {
	return -code error \
	       "edit window name \"$name\" is reserved for Tk $name widgets"
    }
}

#------------------------------------------------------------------------------
# tablelist::createCheckbutton
#
# Creates a checkbutton widget with the given path name for interactive cell
# editing in a tablelist widget.
#------------------------------------------------------------------------------
proc tablelist::createCheckbutton {w args} {
    variable checkedImg
    variable uncheckedImg
    set win [parseEditWinPath $w]
    upvar ::tablelist::ns${win}::data data

    eval [list checkbutton $w -indicatoron 0 -image $uncheckedImg \
	  -selectimage $checkedImg -selectcolor "" \
	  -variable ::tablelist::ns${win}::data(editText)] $args

    if {$::tk_version >= 8.4} {
	$w configure -offrelief sunken
    }
}

#------------------------------------------------------------------------------
# tablelist::createIncrCombobox
#
# Creates an [incr Widgets] combobox with the given path name for interactive
# cell editing in a tablelist widget.
#------------------------------------------------------------------------------
proc tablelist::createIncrCombobox {w args} {
    eval [list iwidgets::combobox $w -dropdown 1 -editable 1 -width 0] $args

    #
    # Make sure that the entry component will receive the input focus
    # whenever the list component (a scrolledlistbox widget) gets unmapped
    #
    bind [$w component list] <Unmap> +[list focus [$w component entry]]
}

#------------------------------------------------------------------------------
# tablelist::createOakleyCombobox
#
# Creates an Oakley combobox widget with the given path name for interactive
# cell editing in a tablelist widget.
#------------------------------------------------------------------------------
proc tablelist::createOakleyCombobox {w args} {
    eval [list combobox::combobox $w -editable 1 -width 0] $args

    #
    # Repack the widget's components, to make sure that the
    # button will remain visible when shrinking the combobox.
    # This patch is needed for combobox versions earlier than 2.3.
    #
    pack forget $w.entry $w.button
    pack $w.button -side right -fill y    -expand 0
    pack $w.entry  -side left  -fill both -expand 1
}

#------------------------------------------------------------------------------
# tablelist::editcellSubCmd
#
# This procedure is invoked to process the tablelist editcell subcommand.
# charPos stands for the character position component of the index in the body
# text widget of the character underneath the mouse cursor if this command was
# invoked by clicking mouse button 1 in the body of the tablelist widget.
#------------------------------------------------------------------------------
proc tablelist::editcellSubCmd {win row col restore {charPos -1}} {
    variable editWin
    upvar ::tablelist::ns${win}::data data

    if {$data(isDisabled) || $data($col-hide) ||
	![isCellEditable $win $row $col]} {
	return ""
    }
    if {$data(editRow) == $row && $data(editCol) == $col} {
	return ""
    }
    if {$data(editRow) >= 0 && ![finisheditingSubCmd $win]} {
	return ""
    }

    #
    # Create a frame to be embedded into the tablelist's
    # body, together with a child of column-specific type
    #
    set f $data(bodyFr)
    frame $f -borderwidth 0 -container 0 -highlightthickness 0 -relief flat \
	     -takefocus 0
    bind $f <Destroy> {
	array set tablelist::ns[winfo parent [winfo parent %W]]::data \
		  {editRow -1  editCol -1}
    }
    set name $data($col-editwindow)
    set creationCmd [strMap {%W $w} $editWin($name-creationCmd)]
    set item [lindex $data(itemList) $row]
    set key [lindex $item end]
    append creationCmd { $editWin($name-fontOpt) [cellFont $win $key $col] } \
	   {-borderwidth 2 -highlightthickness 0 -relief ridge -state normal}
    set w $data(bodyFrEd)
    if {[catch {eval $creationCmd} result] != 0} {
	destroy $f
	return -code error $result
    }
    set class [winfo class $w]
    set isMentry [expr {[string compare $class Mentry] == 0}]
    set isCheckbtn [expr {[string compare $class Checkbutton] == 0}]
    set alignment [lindex $data(colList) [expr {2*$col + 1}]]
    if {!$isMentry} {
	catch {$w configure -justify $alignment}
    }
    clearTakefocusOpt $w

    #
    # Replace the cell contents between the two tabs with the above frame
    #
    set b $data(body)
    set data(editKey) $key
    set data(editRow) $row
    set data(editCol) $col
    findCellTabs $win [expr {$row + 1}] $col tabIdx1 tabIdx2
    set editIdx [$b index $tabIdx1+1c]
    $b delete $editIdx $tabIdx2
    $b window create $editIdx -align bottom -pady -3 -stretch 1 -window $f
    $b mark set editMark $editIdx

    #
    # Insert the binding tag TablelistEdit in the list of binding tags
    # of some components of w, just before the respective path names
    #
    if {$isMentry} {
	set compList [$w entries]
    } else {
	set comp [subst [strMap {%W $w} $editWin($name-focusWin)]]
	set compList [list $comp]
	set data(editFocus) $comp
    }
    foreach comp $compList {
	set bindTags [bindtags $comp]
	set idx [lsearch -exact $bindTags $comp]
	bindtags $comp [linsert $bindTags $idx TablelistEdit]
    }

    #
    # Restore or initialize some of the edit window's data
    #
    if {$restore} {
	restoreEditData $win
    } else {
	#
	# Put the cell's contents to the edit window
	#
	set data(canceled) 0
	set text [lindex $item $col]
	if {$editWin($name-useFormat) &&
	    [info exists data($col-formatcommand)]} {
	    set text [uplevel #0 $data($col-formatcommand) [list $text]]
	}
	catch {eval [strMap {%W $w  %T $text} $editWin($name-putTextCmd)]}
	if {[string compare $data(-editstartcommand) ""] != 0} {
	    set text [uplevel #0 $data(-editstartcommand) \
		      [list $win $row $col $text]]
	    if {$data(canceled)} {
		return ""
	    }
	    catch {eval [strMap {%W $w  %T $text} $editWin($name-putTextCmd)]}
	}

	#
	# Save the edit window's text
	#
	if {$isMentry} {
	    set data(origEditText) [$w getstring]
	} elseif {$isCheckbtn} {
	    set data(origEditText) $text
	} else {
	    set data(origEditText) [$comp get]
	}
	set data(rejected) 0
	set data(invoked) 0

	if {[string compare $editWin($name-getListCmd) ""] != 0 &
	    [string compare $editWin($name-selectCmd) ""] != 0} {
	    #
	    # Select the edit window's item corresponding to text
	    #
	    set itemList [eval [strMap {%W $w} $editWin($name-getListCmd)]]
	    if {[set idx [lsearch -exact $itemList $text]] >= 0} {
		eval [strMap {%W $w  %I $idx} $editWin($name-selectCmd)]
	    }
	}

	#
	# Set the focus and the insertion cursor
	#
	seecellSubCmd $win $row $col
	if {$charPos >= 0} {
	    if {$isCheckbtn} {
		focus $w
	    } else {
		set image [doCellCget $row $col $win -image]
		if {[string compare $alignment right] == 0} {
		    scan $tabIdx2 %d.%d line tabCharIdx2
		    if {$isMentry} {
			set len [string length [$w getstring]]
		    } else {
			set len [$comp index end]
		    }
		    set number [expr {$len - $tabCharIdx2 + $charPos}]
		    if {[string compare $image ""] != 0} {
			incr number 2
		    }
		} else {
		    scan $tabIdx1 %d.%d line tabCharIdx1
		    set number [expr {$charPos - $tabCharIdx1 - 1}]
		    if {[string compare $image ""] != 0} {
			incr number -2
		    }
		}
		if {$isMentry} {
		    setMentryCursor $w $number
		} else {
		    focus $comp
		    $comp icursor $number
		}
	    }
	} else {
	    if {$isMentry || $isCheckbtn} {
		focus $w
	    } else {
		focus $comp
		$comp icursor end
		$comp selection range 0 end
	    }
	}
    }

    #
    # Adjust the frame's dimensions and paddings
    #
    update idletasks
    $f configure -height [winfo reqheight $w]
    set pixels [lindex $data(colList) [expr {2*$col}]]
    if {$pixels == 0} {				;# convention: dynamic width
	set pixels $data($col-reqPixels)
	if {$data($col-maxPixels) > 0 && $pixels > $data($col-maxPixels)} {
	    set pixels $data($col-maxPixels)
	}
    }
    incr pixels $data($col-delta)
    adjustEditWindow $win $pixels
    if {$isCheckbtn} {
	$b window configure $editIdx -align center -stretch 0
    }
    place $w -relheight 1.0 -relwidth 1.0

    #
    # If this command was invoked by clicking mouse button 1 in
    # the tablelist's body and the edit window is a checkbutton,
    # then generate a mouse click event in the checkbutton widget
    #
    if {$charPos >= 0 && $isCheckbtn} {
	after idle [list event generate $w <Button-1>]
    }

    return ""
}

#------------------------------------------------------------------------------
# tablelist::canceleditingSubCmd
#
# This procedure is invoked to process the tablelist cancelediting subcommand.
# Aborts the interactive cell editing and restores the cell's contents after
# destroying the edit window.
#------------------------------------------------------------------------------
proc tablelist::canceleditingSubCmd win {
    upvar ::tablelist::ns${win}::data data

    if {[set row $data(editRow)] < 0} {
	return ""
    }
    set col $data(editCol)

    destroy $data(bodyFr)
    set item [lindex $data(itemList) $row]
    doCellConfig $row $col $win -text [lindex $item $col]
    focus $data(body)
    set data(canceled) 1
    return ""
}

#------------------------------------------------------------------------------
# tablelist::finisheditingSubCmd
#
# This procedure is invoked to process the tablelist finishediting subcommand.
# Invokes the command specified by the -editendcommand option if needed, and
# updates the element just edited after destroying the edit window if the
# latter's content was not rejected.  Returns 1 on normal termination and 0
# otherwise.
#------------------------------------------------------------------------------
proc tablelist::finisheditingSubCmd win {
    variable editWin
    upvar ::tablelist::ns${win}::data data

    if {[set row $data(editRow)] < 0} {
	return 1
    }
    set col $data(editCol)

    #
    # Get the edit window's text, and invoke the command
    # specified by the -editendcommand option if needed
    #
    set w $data(bodyFrEd)
    set class [winfo class $w]
    set isMentry [expr {[string compare $class Mentry] == 0}]
    set isCheckbtn [expr {[string compare $class Checkbutton] == 0}]
    if {$isMentry} {
	set text [$w getstring]
    } elseif {$isCheckbtn} {
	set text $data(editText)
    } else {
	set text [$data(editFocus) get]
    }
    if {[string compare $text $data(origEditText)] == 0} {
	set item [lindex $data(itemList) $row]
	set text [lindex $item $col]
    } else {
	set name $data($col-editwindow)
	set getTextCmd [strMap {%W $w} $editWin($name-getTextCmd)]
	if {[catch {
	    eval [strMap {%W $w} $editWin($name-getTextCmd)]
	} text] != 0} {
	    set data(rejected) 1
	}
	if {[string compare $data(-editendcommand) ""] != 0} {
	    set text \
		[uplevel #0 $data(-editendcommand) [list $win $row $col $text]]
	}
    }

    #
    # Check whether the input was rejected (by the above "set data(rejected) 1"
    # statement or within the command specified by the -editendcommand option)
    #
    if {$data(rejected)} {
	seecellSubCmd $win $row $col
	if {!$isMentry} {
	    focus $data(editFocus)
	}
	set data(rejected) 0
	return 0
    } else {
	destroy $data(bodyFr)
	doCellConfig $row $col $win -text $text
	focus $data(body)
	return 1
    }
}

#------------------------------------------------------------------------------
# tablelist::clearTakefocusOpt
#
# Sets the -takefocus option of all members of the widget hierarchy starting
# with w to 0.
#------------------------------------------------------------------------------
proc tablelist::clearTakefocusOpt w {
    catch {$w configure -takefocus 0}
    foreach c [winfo children $w] {
	clearTakefocusOpt $c
    }
}

#------------------------------------------------------------------------------
# tablelist::adjustEditWindow
#
# Adjusts the width and the horizontal padding of the frame containing the edit
# window associated with the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::adjustEditWindow {win pixels} {
    variable editWin
    upvar ::tablelist::ns${win}::data data

    set class [winfo class $data(bodyFrEd)]
    set name $data($data(editCol)-editwindow)
    if {[string compare $class Checkbutton] == 0} {
	set width [winfo reqwidth $data(bodyFrEd)]
	set padX 0
    } elseif {$editWin($name-useReqWidth) &&
	      [set reqWidth [winfo reqwidth $data(bodyFrEd)]] <=
	      $pixels + 2*$data(charWidth)} {
	set width $reqWidth
	set padX [expr {$reqWidth <= $pixels ? -3 : ($pixels - $reqWidth) / 2}]
    } elseif {$editWin($name-usePadX)} {
	set width [expr {$pixels + 2*$data(charWidth)}]
	set padX -$data(charWidth)
    } else {
	set width [expr {$pixels + 6}]
	set padX -3
    }

    $data(bodyFr) configure -width $width
    $data(body) window configure editMark -padx $padX
}

#------------------------------------------------------------------------------
# tablelist::setMentryCursor
#
# Sets the focus to the entry child of the mentry widget w that contains the
# global character position specified by number, and sets the insertion cursor
# in that entry to the relative character position corresponding to number.  If
# that entry is not enabled then the procedure sets the focus to the last
# enabled entry child preceding the found one and sets the insertion cursor to
# its end.
#------------------------------------------------------------------------------
proc tablelist::setMentryCursor {w number} {
    #
    # Find the entry child containing the given character
    # position; if the latter is contained in a label child
    # then take the entry immediately preceding that label
    #
    set entryIdx -1
    set childIdx 0
    set childCount [llength [$w cget -body]]
    foreach c [winfo children $w] {
	set class [winfo class $c]
	switch $class {
	    Entry {
		set str [$c get]
		set entry $c
		incr entryIdx
	    }
	    Label {
		set str [$c cget -text]
	    }
	}
	set len [string length $str]

	if {$number < $len} {
	    break
	} elseif {$childIdx < $childCount - 1} {
	    incr number -$len
	}

	incr childIdx
    }

    #
    # If the entry's state is normal then set the focus to this entry and
    # the insertion cursor to the relative character position corresponding
    # to number; otherwise set the focus to the last enabled entry child
    # preceding the found one and set the insertion cursor to its end
    #
    switch $class {
	Entry { set relIdx $number }
	Label { set relIdx end }
    }
    if {[string compare [$entry cget -state] normal] == 0} {
	focus $entry
	$entry icursor $relIdx
    } else {
	for {incr entryIdx -1} {$entryIdx >= 0} {incr entryIdx -1} {
	    set entry [$w entrypath $entryIdx]
	    if {[string compare [$entry cget -state] normal] == 0} {
		focus $entry
		$entry icursor end
		return ""
	    }
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::saveEditData
#
# Saves some data of the edit window associated with the tablelist widget win.
#------------------------------------------------------------------------------
proc tablelist::saveEditData win {
    variable editWin
    upvar ::tablelist::ns${win}::data data

    set w $data(bodyFrEd)
    set entry $data(editFocus)
    set class [winfo class $w]
    set isMentry [expr {[string compare $class Mentry] == 0}]
    set isCheckbtn [expr {[string compare $class Checkbutton] == 0}]

    #
    # Miscellaneous data
    #
    if {!$isMentry && !$isCheckbtn} {
	set data(editText) [$entry get]
    }
    set name $data($data(editCol)-editwindow)
    if {[string compare $editWin($name-getListCmd) ""] != 0} {
	set data(editList) \
	    [eval [strMap {%W $w} $editWin($name-getListCmd)]]
    }
    if {!$isCheckbtn} {
	set data(entryPos) [$entry index insert]
	if {[set data(entryHadSel) [$entry selection present]]} {
	    set data(entrySelFrom) [$entry index sel.first]
	    set data(entrySelTo)   [$entry index sel.last]
	}
    }
    set data(entryHadFocus) \
	[expr {[string compare [focus -lastfor $entry] $entry] == 0}]

    #
    # Configuration options and widget callbacks
    #
    saveEditConfigOpts $w
    if {[info exists ::wcb::version] && !$isMentry && !$isCheckbtn} {
	foreach when {before after} {
	    foreach opt {insert delete motion} {
		set data(entryCb-$when-$opt) \
		    [::wcb::callback $entry $when $opt]
	    }
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::saveEditConfigOpts
#
# Saves the non-default values of the configuration options of the edit window
# w associated with a tablelist widget, as well as those of its descendants.
#------------------------------------------------------------------------------
proc tablelist::saveEditConfigOpts w {
    regexp {^(.+)\.body\.f\.(e.*)$} $w dummy win tail
    upvar ::tablelist::ns${win}::data data

    foreach configSet [$w configure] {
	if {[llength $configSet] != 2} {
	    set default [lindex $configSet 3]
	    set current [lindex $configSet 4]
	    if {[string compare $default $current] != 0} {
		set opt [lindex $configSet 0]
		set data($tail$opt) [lindex $configSet 4]
	    }
	}
    }

    foreach c [winfo children $w] {
	saveEditConfigOpts $c
    }
}

#------------------------------------------------------------------------------
# tablelist::restoreEditData
#
# Restores some data of the edit window associated with the tablelist widget
# win.
#------------------------------------------------------------------------------
proc tablelist::restoreEditData win {
    variable editWin
    upvar ::tablelist::ns${win}::data data

    set w $data(bodyFrEd)
    set entry $data(editFocus)
    set class [winfo class $w]
    set isMentry [expr {[string compare $class Mentry] == 0}]
    set isCheckbtn [expr {[string compare $class Checkbutton] == 0}]
    set isIncrDateTimeWidget [regexp {^(Date.+|Time.+)$} $class]

    #
    # Miscellaneous data
    #
    if {!$isMentry && !$isCheckbtn && !$isIncrDateTimeWidget} {
	$entry delete 0 end
	$entry insert 0 $data(editText)
    }
    set name $data($data(editCol)-editwindow)
    if {[string compare $editWin($name-putListCmd) ""] != 0 &&
	[string compare $data(editList) ""] != 0} {
	eval [strMap {%W $w  %L $data(editList)} $editWin($name-putListCmd)]
    }
    if {[string compare $editWin($name-selectCmd) ""] != 0 &&
	[set idx [lsearch -exact $data(editList) $data(editText)]] >= 0} {
	eval [strMap {%W $w  %I $idx} $editWin($name-selectCmd)]
    }
    if {!$isCheckbtn} {
	$entry icursor $data(entryPos)
	if {$data(entryHadSel)} {
	    $entry selection range $data(entrySelFrom) $data(entrySelTo)
	}
    }
    if {$data(entryHadFocus)} {
	focus $entry
    }

    #
    # Configuration options and widget callbacks
    #
    restoreEditConfigOpts $w
    if {[info exists ::wcb::version] && !$isMentry && !$isCheckbtn} {
	foreach when {before after} {
	    foreach opt {insert delete motion} {
		eval [list ::wcb::callback $entry $when $opt] \
		     $data(entryCb-$when-$opt)
	    }
	}
    }

    #
    # If the edit window is a datefield, dateentry, timefield, or timeentry
    # widget then restore its text here, because otherwise it would be
    # overridden when the above invocation of restoreEditConfigOpts sets
    # the widget's -format option.  Note that this is a special case; in
    # general we must restore the text BEFORE the configuration options.
    #
    if {$isIncrDateTimeWidget} {
	$entry delete 0 end
	$entry insert 0 $data(editText)
    }
}

#------------------------------------------------------------------------------
# tablelist::restoreEditConfigOpts
#
# Restores the non-default values of the configuration options of the edit
# window w associated with a tablelist widget, as well as those of its
# descendants.
#------------------------------------------------------------------------------
proc tablelist::restoreEditConfigOpts w {
    regexp {^(.+)\.body\.f\.(e.*)$} $w dummy win tail
    upvar ::tablelist::ns${win}::data data

    set isMentry [expr {[string compare [winfo class $w] Mentry] == 0}]

    foreach name [array names data $tail-*] {
	set opt [string range $name [string last - $name] end]
	if {!$isMentry || [string compare $opt -body] != 0} {
	    $w configure $opt $data($name)
	}
	unset data($name)
    }

    foreach c [winfo children $w] {
	restoreEditConfigOpts $c
    }
}

#
# Private procedures used in bindings related to interactive cell editing
# =======================================================================
#

#------------------------------------------------------------------------------
# tablelist::insertChar
#
# Inserts the string str into the entry or spinbox widget w at the point of the
# insertion cursor.
#------------------------------------------------------------------------------
proc tablelist::insertChar {w str} {
    set class [winfo class $w]
    if {[string compare $class Entry] != 0 &&
	[string compare $class Spinbox] != 0} {
	return ""
    }

    if {[string compare [info procs ::tkEntryInsert] ::tkEntryInsert] == 0} {
	tkEntryInsert $w $str
    } else {
	tk::EntryInsert $w $str
    }
}

#------------------------------------------------------------------------------
# tablelist::cancelEditing
#
# Invokes the canceleditingSubCmd procedure.
#------------------------------------------------------------------------------
proc tablelist::cancelEditing w {
    if {[isComboTopMapped $w]} {
	return ""
    }

    canceleditingSubCmd [parseEditWinPath $w]
    return -code break ""
}

#------------------------------------------------------------------------------
# tablelist::finishEditing
#
# Invokes the finisheditingSubCmd procedure.
#------------------------------------------------------------------------------
proc tablelist::finishEditing w {
    if {[isComboTopMapped $w]} {
	return ""
    }

    finisheditingSubCmd [parseEditWinPath $w]
    return -code break ""
}

#------------------------------------------------------------------------------
# tablelist::goToNextCell
#
# Moves the edit window into the next editable cell different from the one
# indicated by the given row and column, if there is such a cell.
#------------------------------------------------------------------------------
proc tablelist::goToNextCell {w args} {
    if {[isComboTopMapped $w]} {
	return ""
    }

    set win [parseEditWinPath $w]
    upvar ::tablelist::ns${win}::data data

    if {[llength $args] == 0} {
	set row $data(editRow)
	set col $data(editCol)
    } else {
	set row [lindex $args 0]
	set col [lindex $args 1]
    }

    set oldRow $row
    set oldCol $col

    while 1 {
	incr col
	if {$col > $data(lastCol)} {
	    incr row
	    if {$row > $data(lastRow)} {
		set row 0
	    }
	    set col 0
	}

	if {$row == $oldRow && $col == $oldCol} {
	    return -code break ""
	} elseif {!$data($col-hide) && [isCellEditable $win $row $col]} {
	    editcellSubCmd $win $row $col 0
	    return -code break ""
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::goToPrevCell
#
# Moves the edit window into the previous editable cell different from the one
# indicated by the given row and column, if there is such a cell.
#------------------------------------------------------------------------------
proc tablelist::goToPrevCell {w args} {
    if {[isComboTopMapped $w]} {
	return ""
    }

    set win [parseEditWinPath $w]
    upvar ::tablelist::ns${win}::data data

    if {[llength $args] == 0} {
	set row $data(editRow)
	set col $data(editCol)
    } else {
	set row [lindex $args 0]
	set col [lindex $args 1]
    }

    set oldRow $row
    set oldCol $col

    while 1 {
	incr col -1
	if {$col < 0} {
	    incr row -1
	    if {$row < 0} {
		set row $data(lastRow)
	    }
	    set col $data(lastCol)
	}

	if {$row == $oldRow && $col == $oldCol} {
	    return -code break ""
	} elseif {!$data($col-hide) && [isCellEditable $win $row $col]} {
	    editcellSubCmd $win $row $col 0
	    return -code break ""
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::goLeft
#
# Moves the edit window into the previous editable cell of the current row if
# the cell being edited is not the first editable one within that row.
# Otherwise sets the insertion cursor to the beginning of the entry and clears
# the selection in it.
#------------------------------------------------------------------------------
proc tablelist::goLeft w {
    if {[isComboTopMapped $w]} {
	return ""
    }

    set win [parseEditWinPath $w]
    upvar ::tablelist::ns${win}::data data

    set row $data(editRow)
    set col $data(editCol)

    while 1 {
	incr col -1
	if {$col < 0} {
	    set class [winfo class $w]
	    if {[string compare $class Entry] == 0 ||
		[string compare $class Spinbox] == 0} {
		#
		# On Windows the "event generate" command does not behave
		# as expected if a Tk version older than 8.2.2 is used.
		#
		if {[string compare $::tk_patchLevel 8.2.2] < 0} {
		    tkEntrySetCursor $w 0
		} else {
		    event generate $w <Home>
		}
	    }
	    return -code break ""
	} elseif {!$data($col-hide) && [isCellEditable $win $row $col]} {
	    editcellSubCmd $win $row $col 0
	    return -code break ""
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::goRight
#
# Moves the edit window into the next editable cell of the current row if the
# cell being edited is not the last editable one within that row.  Otherwise
# sets the insertion cursor to the end of the entry and clears the selection in
# it.
#------------------------------------------------------------------------------
proc tablelist::goRight w {
    if {[isComboTopMapped $w]} {
	return ""
    }

    set win [parseEditWinPath $w]
    upvar ::tablelist::ns${win}::data data

    set row $data(editRow)
    set col $data(editCol)

    while 1 {
	incr col
	if {$col > $data(lastCol)} {
	    set class [winfo class $w]
	    if {[string compare $class Entry] == 0 ||
		[string compare $class Spinbox] == 0} {
		#
		# On Windows the "event generate" command does not behave
		# as expected if a Tk version older than 8.2.2 is used.
		#
		if {[string compare $::tk_patchLevel 8.2.2] < 0} {
		    tkEntrySetCursor $w end
		} else {
		    event generate $w <End>
		}
	    }
	    return -code break ""
	} elseif {!$data($col-hide) && [isCellEditable $win $row $col]} {
	    editcellSubCmd $win $row $col 0
	    return -code break ""
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::goUp
#
# Invokes the goToPrevLine procedure.
#------------------------------------------------------------------------------
proc tablelist::goUp w {
    if {[isComboTopMapped $w]} {
	return ""
    }

    set win [parseEditWinPath $w]
    upvar ::tablelist::ns${win}::data data

    goToPrevLine $w $data(editRow) $data(editCol)
    return -code break ""
}

#------------------------------------------------------------------------------
# tablelist::goToPrevLine
#
# Moves the edit window into the last editable cell that is located in the
# specified column and has a row index less than the given one, if there is
# such a cell.
#------------------------------------------------------------------------------
proc tablelist::goToPrevLine {w row col} {
    set win [parseEditWinPath $w]

    while 1 {
	incr row -1
	if {$row < 0} {
	    return 0
	} elseif {[isCellEditable $win $row $col]} {
	    editcellSubCmd $win $row $col 0
	    return 1
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::goDown
#
# Invokes the goToNextLine procedure.
#------------------------------------------------------------------------------
proc tablelist::goDown w {
    if {[isComboTopMapped $w]} {
	return ""
    }

    set win [parseEditWinPath $w]
    upvar ::tablelist::ns${win}::data data

    goToNextLine $w $data(editRow) $data(editCol)
    return -code break ""
}

#------------------------------------------------------------------------------
# tablelist::goToNextLine
#
# Moves the edit window into the first editable cell that is located in the
# specified column and has a row index greater than the given one, if there is
# such a cell.
#------------------------------------------------------------------------------
proc tablelist::goToNextLine {w row col} {
    set win [parseEditWinPath $w]
    upvar ::tablelist::ns${win}::data data

    while 1 {
	incr row
	if {$row > $data(lastRow)} {
	    return 0
	} elseif {[isCellEditable $win $row $col]} {
	    editcellSubCmd $win $row $col 0
	    return 1
	}
    }
}

#------------------------------------------------------------------------------
# tablelist::goToPrevPage
#
# Moves the edit window up by one page within the current column if the cell
# being edited is not the first editable one within that column.
#------------------------------------------------------------------------------
proc tablelist::goToPrevPage w {
    if {[isComboTopMapped $w]} {
	return ""
    }

    set win [parseEditWinPath $w]
    upvar ::tablelist::ns${win}::data data

    #
    # Check whether there is any editable cell
    # above the current one, in the same column
    #
    set row $data(editRow)
    set col $data(editCol)
    while 1 {
	incr row -1
	if {$row < 0} {
	    return -code break ""
	} elseif {[isCellEditable $win $row $col]} {
	    break
	}
    }

    #
    # Scroll up the view by one page and get the corresponding row index
    #
    set row $data(editRow)
    seeSubCmd $win $row
    set bbox [bboxSubCmd $win $row]
    yviewSubCmd $win {scroll -1 pages}
    set newRow [rowIndex $win @0,[lindex $bbox 1] 0]

    if {$newRow < $row} {
	if {![goToPrevLine $w [expr {$newRow + 1}] $col]} {
	    goToNextLine $w $newRow $col
	}
    } else {
	goToNextLine $w -1 $col
    }

    return -code break ""
}

#------------------------------------------------------------------------------
# tablelist::goToNextPage
#
# Moves the edit window down by one page within the current column if the cell
# being edited is not the last editable one within that column.
#------------------------------------------------------------------------------
proc tablelist::goToNextPage w {
    if {[isComboTopMapped $w]} {
	return ""
    }

    set win [parseEditWinPath $w]
    upvar ::tablelist::ns${win}::data data

    #
    # Check whether there is any editable cell
    # below the current one, in the same column
    #
    set row $data(editRow)
    set col $data(editCol)
    while 1 {
	incr row
	if {$row > $data(lastRow)} {
	    return -code break ""
	} elseif {[isCellEditable $win $row $col]} {
	    break
	}
    }

    #
    # Scroll down the view by one page and get the corresponding row index
    #
    set row $data(editRow)
    seeSubCmd $win $row
    set bbox [bboxSubCmd $win $row]
    yviewSubCmd $win {scroll 1 pages}
    set newRow [rowIndex $win @0,[lindex $bbox 1] 0]

    if {$newRow > $row} {
	if {![goToNextLine $w [expr {$newRow - 1}] $col]} {
	    goToPrevLine $w $newRow $col
	}
    } else {
	goToPrevLine $w $data(itemCount) $col
    }

    return -code break ""
}

#------------------------------------------------------------------------------
# tablelist::genMouseWheelEvent
#
# Generates a <MouseWheel> event with the given delta on the widget w.
#------------------------------------------------------------------------------
proc tablelist::genMouseWheelEvent {w delta} {
    set focus [focus -displayof $w]
    focus $w
    event generate $w <MouseWheel> -delta $delta
    focus $focus
}

#------------------------------------------------------------------------------
# tablelist::parseEditWinPath
#
# Extracts the path name of the tablelist widget from the path name w of a
# widget that is assumed to be the edit window or one of its descendants.
#------------------------------------------------------------------------------
proc tablelist::parseEditWinPath w {
    regexp {^(.+)\.body\.f\.e.*$} $w dummy win
    return $win
}

#------------------------------------------------------------------------------
# tablelist::isKeyReserved
#
# Checks whether the given keysym is used in the standard binding scripts
# associated with the widget w, which is assumed to be the edit window or one
# of its descendants.
#------------------------------------------------------------------------------
proc tablelist::isKeyReserved {w keySym} {
    variable editWin
    set win [parseEditWinPath $w]
    upvar ::tablelist::ns${win}::data data

    set name $data($data(editCol)-editwindow)
    return [expr {[lsearch -exact $editWin($name-reservedKeys) $keySym] >= 0}]
}

#------------------------------------------------------------------------------
# tablelist::isComboTopMapped
#
# Checks whether the given widget is a component of an Oakley combobox having
# its toplevel child mapped.  This is needed in our binding scripts to make
# sure that the interactive cell editing won't be terminated prematurely,
# because Bryan Oakley's combobox keeps the focus on its entry child even if
# its toplevel component is mapped.
#------------------------------------------------------------------------------
proc tablelist::isComboTopMapped w {
    set par [winfo parent $w]
    if {[string compare [winfo class $par] Combobox] == 0 &&
	[winfo exists $par.top] && [winfo ismapped $par.top]} {
	return 1
    } else {
	return 0
    }
}
