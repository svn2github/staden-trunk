# A menubutton with a style which appears the same size and shape as the
# default for a normal button.
# The purpose of this is to make the packing look smart.
proc xmenubutton {path args} {
    global tcl_platform

    if {$tcl_platform(platform) == "unix"} {
	return [eval menubutton $path \
		    -bd 2 \
		    -relief raised \
		    -highlightthickness 1 \
		    -padx 2 \
		    $args]
    } else {
	return [eval menubutton $path \
		    -bd 2 \
		    -relief raised \
		    -highlightthickness 0 \
		    -padx 2 \
		    $args]
    }
}