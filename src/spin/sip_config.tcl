#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc UpdateSipPlot {index colour} {

} 

proc ConfigureSipPlot {index line_width colour} {
    raster_config -index $index -width $line_width -fill $colour
} 

proc ConfigureSipWidth {index line_width} {
    raster_config -index $index -width $line_width
} 

#called from C
proc SipConfig {index} {

    if {[xtoplevel [set t .cbox] -resizable 0] == ""} return

    #puts "SIPCONFIG index $index"

    set line_width [lindex [lindex [raster_getconfig -index $index] 1] 1]

    frame $f.lw
    # Line width
    scale $f.lw.scale \
	    -label "Line width" \
	    -from 0 \
	    -to 10 \
	    -orient horiz
    $f.lw.scale set $line_width

    pack $f.lw -side top -fill both
    pack $f.lw.scale -side top -fill x

    #starting colour
    set colour [lindex [lindex [raster_getconfig -index $index] 0] 1]

    #cmd to execute when ok button on colourbox pressed
    set ok_cmd "ConfigureSipPlot $index \[$f.lw.scale get\]"

    #cmd to execute when changing colours on colourbox
    set update_cmd "UpdateSipPlot $index"

    ColourBox $f $colour $ok_cmd $update_cmd {}
}
