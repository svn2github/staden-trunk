proc ShowLicence {} {
    global licence licence licence

    # Create the window
    set w .licence
    modal $w
    wm title $w "Licence"
    pack [frame $w.l -bd 2 -relief sunken] -side top -fill both

    set width 90
    # Show the Licence Type
    label $w.l.type_l -text "Licence Type"
    switch $licence(type) {
	f	{set type full}
	v	{set type viewer}
	default	{set type demo}
    }
    entry $w.l.type_v -width $width
    $w.l.type_v insert 0 $type
    $w.l.type_v configure -state disabled
    grid $w.l.type_l $w.l.type_v

    # Show the Licence OS
    label $w.l.os_l -text "Licence Type"
    switch $licence(os) {
	u	{set os Unix}
	w	{set os {MS Windows}}
	m	{set os {MacOS X}}
	default	{set os Unknown}
    }
    entry $w.l.os_v -width $width
    $w.l.os_v insert 0 $os
    $w.l.os_v configure -state disabled
    grid $w.l.os_l $w.l.os_v
    
    # Show the expiry date
    if {$licence(expire) == 0} {
	set expire "never"
    } else {
	set expire [clock format $licence(expire)]
    }
    label $w.l.expire_l -text "Expiry date"
    entry $w.l.expire_v -width $width
    $w.l.expire_v insert 0 $expire
    $w.l.expire_v configure -state disabled
    grid $w.l.expire_l $w.l.expire_v

    # Show the number of users
    if {$licence(users) == 0} {
	set users "unlimited"
    } else {
	set users $licence(users)
    }
    label $w.l.users_l -text "No. of users"
    entry $w.l.users_v -width $width
    $w.l.users_v insert 0 $users
    $w.l.users_v configure -state disabled
    grid $w.l.users_l $w.l.users_v

    # Show the Host ID
    label $w.l.hostid_l -text "Host ID"
    entry $w.l.hostid_v -width $width
    $w.l.hostid_v insert 0 $licence(hostid)
    grid $w.l.hostid_l $w.l.hostid_v -pady 5

    # Info. on how to obtain a valid licence
    if {$licence(type) == "d" || $expire != "never"} {
	label $w.message1 \
		-text "To obtain a full permanent licence please contact"
	entry $w.message2 -bd 0
	$w.message2 insert 0 "\tstaden-admin@mrc-lmb.cam.ac.uk"
	label $w.message3 -text "or read"
	entry $w.message4 -bd 0
	$w.message4 insert 0 "\thttp://www.mrc-lmb.cam.ac.uk/pubseq/licence.html"
	pack $w.message1 -side top -anchor w
	pack $w.message2 -side top -fill both
	pack $w.message3 -side top -anchor w
	pack $w.message4 -side top -fill both

	if {$licence(type) == "d"} {
	    label $w.message5 -text \
		    "Alternatively free time-limited full licences are also available." \
		    -justify left
	    pack $w.message5 -side top -anchor w
	}
    }

    button $w.button -text OK -command "destroy $w"
    pack $w.button -side bottom -anchor c
    focus $w.button
}


# Display licence information for demonstration versions.
#
# Prog is the program name
#
# demo_mode is how the program handles demonstration mode.
# 0 => does not work at all
# 1 => works fully with a limited data set
# 2 => limited functionality (any data)
# 3 => works fully with a limited data set and read-only with other data
proc LicenceSplash {prog demo_mode} {
    global licence tcl_platform

    xtoplevel [set w .licence_warning]
    wm resizable $w 0 0
    wm withdraw $w

    wm protocol $w WM_DELETE_WINDOW exit

    pack [text $w.t -width 100 -height 18 -wrap word] -fill both -expand 1
    $w.t tag configure title -font {Helvetica 10 bold} -justify center
    $w.t tag configure bold -font {Courier 10 bold} -justify center
    $w.t tag configure hid -font {Courier 10 bold} -justify center

    $w.t insert end "$prog Demonstration Version" title
    $w.t insert end "\n\nThis is a demonstration release of the Staden Package. "
    if {$demo_mode == 1} {
	$w.t insert end "$prog will only be able to load and edit the\
		supplied demonstration files.\n"
    } elseif {$demo_mode == 2} {
	$w.t insert end "$prog will run, but with limited functionality.\n"
    } elseif {$demo_mode == 3} {
	$w.t insert end "$prog will only be able to load and edit the\
		supplied demonstration files. \
                Other data may be loaded in a read-only mode.\n"
    } else {
	$w.t insert end "$prog is not available in demonstration mode.\n"
    }
    if {$tcl_platform(platform) == "unix"} {
        $w.t insert end "\nTo obtain a free time-limited licence enabling\
	    full use of the package, or to obtain (free to academics) an\
	    unlimited licence, please connect to the following URL:"
    } else {
        $w.t insert end "\nTo obtain a free time-limited licence enabling\
	    full use of the package, or to purchase an unlimited licence,\
	    please connect to the following URL:"
    }
    $w.t insert end "\n\nhttp://www.mrc-lmb.cam.ac.uk/pubseq/licence.html" bold
    $w.t insert end "\n\nTo obtain a licence you will need to know your\
	    Host ID. This is:\n\n"
    if {[llength [split $licence(hostid) /]] > 5} {
	regsub {(/[^/]*/[^/]*/[^/]*)} $licence(hostid) "\\1\n" h
    } else {
	set h $licence(hostid)
    }
    $w.t insert end $h hid

    frame $w.f
    pack $w.f -side bottom -fill both
    if {$tcl_platform(platform) == "unix"} {
	button $w.f.c -text "Copy HostID to cut-buffer" \
	    -command "LicenceSelectionOwn $w PRIMARY"
    } else {
	button $w.f.c -text "Copy HostID to clipboard" \
	    -command "LicenceSelectionOwn $w CLIPBOARD"
    }
    pack $w.f.c -side right -fill both -expand 1
    button $w.f.b -text "OK" -command "destroy $w; set .licence_acknowledged 1"
    pack $w.f.b -side left -fill both -expand 1
    focus $w.f.b

    update idletasks
    set x [expr {[winfo screenwidth $w]/2 - [winfo reqwidth $w]/2 \
	    - [winfo vrootx [winfo parent $w]]}]
    set y [expr {[winfo screenheight $w]/2 - [winfo reqheight $w]/2 \
	    - [winfo vrooty [winfo parent $w]]}]
    wm geom $w +$x+$y
    wm deiconify $w

    global .licence_acknowledged
    vwait .licence_acknowledged
    unset .licence_acknowledged

    if {$demo_mode == 0} {
	exit
    }
}

proc LicenceSelectionOwn {w where} {
    # Grab selection
    selection own -selection $where -command "LicenceSelectionCancel $w" .
    if {[selection own -selection $where] == "."} {
	$w.t tag configure hid -background #d0d0ff
    } 

   if {$where == "CLIPBOARD"} {
	global licence
	clipboard clear
	clipboard append $licence(hostid)
    } else {      
	selection handle \
	    -selection $where \
	    . LicenceSelectionHandler

	bind $w <Destroy> "selection handle -selection $where . {}"
    }
}

proc LicenceSelectionCancel {w} {
    $w.t tag configure hid -background [$w.t cget -background]
}

proc LicenceSelectionHandler {args} {
    global licence
    return $licence(hostid)
}

#tkinit
#wm withdraw .
#ShowLicence
