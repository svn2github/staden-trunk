#==============================================================================
# Contains some Tk option database settings for the demo scripts "bwidget.tcl",
# "iwidgets.tcl", and "miscWidgets.tcl".
#
# Copyright (c) 2004  Csaba Nemethi (E-mail: csaba.nemethi@t-online.de)
#==============================================================================

#
# Get the current windowing system ("x11", "win32", "classic", or "aqua")
#
if {[catch {tk windowingsystem} winSys] != 0} {
    switch $::tcl_platform(platform) {
	unix      { set winSys x11 }
	windows   { set winSys win32 }
	macintosh { set winSys classic }
    }
}

#
# Add some entries to the Tk option database
#
switch $winSys {
    x11     { option add *Font		"Helvetica -12" }
    classic { option add *background	#dedede }
}
option add *Tablelist.activeStyle	frame
option add *Tablelist.setGrid		yes
option add *Tablelist.movableColumns	yes
option add *Tablelist.labelCommand	tablelist::sortByColumn
option add *Tablelist.background	gray96
option add *Tablelist.stripeBackground	#e0e8f0
option add *Tablelist*selectBackground	navy
option add *Tablelist*selectForeground	white
