#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

##############################################################################
proc get_colour { colframe } {

    set colour [lindex [$colframe.sample configure -background] 4]
    return $colour
}

##############################################################################
proc ConfigColour {colframe cmd dummy} {
    # Taken pretty much directly from the Tk book
    set colour [format #%02x%02x%02x \
	[$colframe.red get] \
	[$colframe.green get] \
	[$colframe.blue get]]
    $colframe.sample configure -background $colour

    #to update the item colour as move the scale
    eval $cmd $colour

}

##############################################################################
proc ColourBox { t colour ok_cmd update_cmd cancel_cmd } {

    # Colour
    frame $t.col -bd 3 -relief raised
    label $t.col.label -text "Colour"
    scale $t.col.red \
	-label "Red" \
	-from 0 \
	-to 255 \
	-orient horiz \
	-length 7c \
	-bd 1 \
	-relief groove \
	-command "ConfigColour $t.col {$update_cmd}"
    scale $t.col.green \
	-label "Green" \
	-from 0 \
	-to 255 \
	-orient horiz \
	-length 7c \
	-bd 1 \
	-relief groove \
	-command "ConfigColour $t.col {$update_cmd}"
    scale $t.col.blue \
	-label "Blue" \
	-from 0 \
	-to 255 \
	-orient horiz \
	-length 7c \
	-bd 1 \
	-relief groove \
	-command "ConfigColour $t.col {$update_cmd}"
    frame $t.col.sample -height 1.5c -width 6c

    set rgb [winfo rgb . $colour]
    set r [expr [lindex $rgb 0]/256]
    set g [expr [lindex $rgb 1]/256]
    set b [expr [lindex $rgb 2]/256]

    $t.col.red set $r
    $t.col.green set $g
    $t.col.blue set $b
    
    pack $t.col -side top -fill both
    pack $t.col.label -side top -fill x
    pack $t.col.red $t.col.green $t.col.blue -side top -fill both
    pack $t.col.sample -side bottom -padx 2m -pady 2m

    #########################################################################
    #OK and Cancel buttons
    if {$cancel_cmd != ""} {
	set cancel_cmd "eval $cancel_cmd"
    }
    okcancelhelp $t.ok_cancel \
	    -ok_command "ColSel_OK_Pressed $t.col {$ok_cmd}; destroy $t" \
	    -cancel_command "destroy $t; $cancel_cmd" \
	    -help_command "show_help interface {UI-Colour}" \
	    -bd 2 \
	    -relief groove
    pack $t.ok_cancel -side bottom -fill x

}

proc ColSel_OK_Pressed { colframe cmd} {

    set colour [get_colour $colframe]
    eval $cmd $colour
}
