# A horribly hacky example of a tk canvas view of contig 0 from a tgap
# database. It's inefficient in caching data too, but despite that is a
# proof of concept that we're fast enough right now.
#
# Warning, this is full of hard coded Tk pathnames, Tcl global variables,
# assumptions of contig length, etc. All this ensures we'll not feel tempted
# to steal any of this hacked up demo code :-)

package require Tk

#
# Constructs a container window for drawing horizontal tracks of data.
# It consists of standard controls like menus, horizontal scrollbars and
# a scale control along with zero or more children plots, named "tracks".
#
# Each track has it's own ID value and function associated with it.
# The 1.5Dplot passes high level events (like scroll in X) down to all
# suitable child tracks.
#
proc 1.5Dplot {w io wid hei {cnum {}}} {
    global $w

    # Create the window
    if {[toplevel $w] == ""} return

    wm geometry $w ${wid}x${hei}

    if {$cnum == ""} {
	set cnum [$io contig_order 0]
    }

    set ${w}(io) $io
    set ${w}(cnum) $cnum
    set ${w}(contig) [$io get_contig $cnum]
    set ${w}(start) [[set ${w}(contig)] get_start]
    set ${w}(length) [[set ${w}(contig)] get_length]
    set ${w}(tracks) {}
    set ${w}(move_status) 0
    set ${w}(yzoom) 1
    set ${w}(xorigin) 0
    set ${w}(width) $wid
    set ${w}(pwidth) $wid; # 1st guess
    set ${w}(height) $hei
    set ${w}(border) 500
    set ${w}(x1)    0
    set ${w}(x2)    100
    set ${w}(last_x1) -9999999
    set ${w}(last_x2) -9999999
    set ${w}(ntracks) 0
    # something to store the common x range information
    set ${w}(grange) [g5::range -io $io -cnum $cnum]

    set ${w}(x1) 0
    set ${w}(x2) [expr {[set ${w}(x1)]+1000}]
    
    # starting tracks
    set ${w}(template) 1
    set ${w}(depth) 1

    wm title $w "Contig [[set ${w}(contig)] get_name]"
    
    # now to do the menu
    global new_template_menu
    menu $w.menubar
    $w configure -menu $w.menubar
    set m $w.menubar
    create_menus $new_template_menu $w.menubar
    

    # Bottom control panel
    set bc [frame $w.bcontrol -bd 0]
    grid columnconfigure $w 0 -weight 0
    grid columnconfigure $w 1 -weight 1
    grid $bc -column 1 -row 999 -sticky nsew
    grid rowconfigure $w 999 -weight 0

    scale $bc.xzoom -from 1 -to 250 -orient horiz -label "X Scale" \
	-resolution 0.1 -command "set_xzoom $w" -repeatinterval 20
    $bc.xzoom set 20
    pack $bc.xzoom -fill both -expand 1 -side left


    # X scrollbar
    scrollbar $w.xscroll \
	-command "scrollx1.5 $w" \
	-repeatinterval 5 \
	-orient horiz
    $w.xscroll set \
	[expr {[set ${w}(x1)]/double([set ${w}(length)])}] \
	[expr {[set ${w}(x2)]/double([set ${w}(length)])}]

    grid $w.xscroll -column 1 -row 998 -sticky nsew
    grid rowconfigure $w 998 -weight 0

    # An information label
    frame $w.label
    grid $w.label -column 1 -row 1000 -sticky nsew
    grid rowconfigure $w 1000 -weight 0
    label $w.label.l -textvariable ${w}(info)
    pack $w.label.l

    bind $w <5> "zoom1.5 $w 0 1 1.1"
    bind $w <4> "zoom1.5 $w 0 1 [expr {1/1.1}]"

    bind $w <Any-Configure> "after idle {resize1.5 $w}"

    set ${w}(pwidth) [expr {[set ${w}(width)]-21}]
    puts pwdth=[set ${w}(pwidth)]
    
    set ${w}(xorigin) [set ${w}(x1)]
    set ${w}(xzoom) [expr {([set ${w}(x2)]-[set ${w}(x1)]+1)/
			   double([set ${w}(pwidth)])}]


    # Contig registration
    set ${w}(reg) [contig_register \
		       -io $io \
		       -contig $cnum \
		       -command "1.5plot_contig_event $w" \
		       -flags [list QUERY_NAME CURSOR_NOTIFY]]
    wm protocol $w WM_DELETE_WINDOW "1.5plot_exit $w"

    redraw_plot $w
}

proc 1.5plot_exit {w} {
    global $w
    rename [set ${w}(grange)] ""
    contig_deregister -io [set ${w}(io)] -id [set ${w}(reg)]
    destroy $w
    unset $w
}

# Contig event handling
proc 1.5plot_contig_event {w type id cdata args} {
    global $w

#    puts [info level [info level]]
    switch $type {
	QUERY_NAME {
	    return "Template Display"
	}

	CURSOR_NOTIFY {
	    set c_col [list blue green orange red]
	    array set arg {}
	    foreach a $args {
		foreach {k v} $a break;
		set arg($k) $v
	    }
	    set cid $arg(id)

	    foreach id [set ${w}(tracks)] {
		set t $w.track$id
		global $t $t.Cursors $t.Id2Cid

		if {![info exists ${t}(Raster)]} continue;
		set r [set ${t}(Raster)]

		if {[lsearch $arg(job) DELETE] != -1} {
		    catch {unset $t.Cursors($cid)}
		    catch {destroy $t.cursor$cid}
		    continue
		}

		set $t.Id2Cid($arg(sent_by)) $cid

		if {![winfo exists $t.cursor$cid]} {
		    set col [lindex $c_col [expr {$cid % [llength $c_col]}]]
		    frame $t.cursor$cid -width 5 -height 10000 -bg $col
		    raise $t.cursor$cid
		    bind $t.cursor$cid <ButtonPress-1> \
			"1.5cursor_press $t $cid"
		    bind $t.cursor$cid <ButtonRelease-1> \
			"1.5cursor_release $t $cid"
		    bind $t.cursor$cid <Any-Motion> \
			"1.5cursor_motion $w $t $cid %X"
		}
		set $t.Cursors($cid) $arg(abspos)
		set x [lindex [$r topixmap $arg(abspos) 0] 0]
		incr x -2

		global $t.cursorx$cid
		set $t.cursorx$cid $x
		place $t.cursor$cid -x $x
	    }
	}

	default {
	    puts "depth.tcl: Event '$type $args' not handled"
	}
    }
}

# Redraws cursors
proc 1.5redraw_cursor {w t} {
    global $w $t $t.Cursors

    if {![info exists ${t}(Raster)]} return;
    set r [set ${t}(Raster)]

    foreach id [array names $t.Cursors] {
	place $t.cursor$id -x [lindex [$r topixmap [set $t.Cursors($id)] 0] 0]
    }
}

# Cursor drag events
proc 1.5cursor_press {t id} {
    global $t.cursor_selected_$id
    set $t.cursor_selected_$id 1
}

proc 1.5cursor_motion {w t id x} {
    global $w $t $t.cursor_selected_$id
    if {![info exists $t.cursor_selected_$id]} return

    if {![info exists ${t}(Raster)]} return;
    set r [set ${t}(Raster)]

    incr x -[winfo rootx $t]
    set bx [expr {round([lindex [$r toworld $x 0] 0])}]

    contig_notify \
	-io [set ${w}(io)] \
	-cnum [set ${w}(cnum)] \
	-type CURSOR_NOTIFY \
	-args [list id      $id \
	            job     MOVE \
		    seq     [set ${w}(cnum)] \
		    abspos  $bx \
		    pos     $bx \
		    sent_by 0]
}

proc 1.5cursor_release {t id} {
    global $t.cursor_selected_$id
    catch {unset $t.cursor_selected_$id}
}

# Resize events
proc resize1.5 {w} {
    global $w

    if {![info exists ${w}(xzoom)]} return

    # We get rogue events when highlighting. Skip these
    if {[winfo width  $w] == [set ${w}(width)] && \
	[winfo height $w] == [set ${w}(height)]} return

    set ${w}(width) [winfo width $w]
    set ${w}(height) [winfo height $w]
    set ${w}(pwidth) [winfo width $w.xscroll]

    foreach id [set ${w}(tracks)] {
	set t $w.track$id
	global $t
	set width  [winfo width  [set ${t}(canvas)]]
	set height [winfo height [set ${t}(canvas)]]
	set ${t}(width)  $width
	set ${t}(height) $height
    }

    set ${w}(x2) [expr {[set ${w}(xzoom)]*[set ${w}(pwidth)] +
			+[set ${w}(x1)] - 1}]
    redraw_plot $w
}

# Scroll in X
proc scrollx1.5 {w cmd args} {
    global $w

    puts [info level [info level]]

    set sbar $w.xscroll
    set clen [set ${w}(length)]

    puts ""
    #parray $w

    if {$cmd == "moveto"} {
	set xpos [expr {int([lindex $args 0]*$clen)}]
    } elseif {$cmd == "set_xpos"} {
	set xpos [lindex $args 0]
    } elseif {$cmd == "scroll"} {
	set xpos [expr {[lindex [$sbar get] 0]*$clen}]
	if {[lindex $args 1] == "pages"} {
	    set wid [expr {[c2x $w [set ${w}(pwidth)]] - [c2x $w 0]}]
	    set xpos [expr {$xpos + $wid/2*[lindex $args 0]}]
	} else {
	    set wid [expr {[c2x $w [lindex $args 0]] - [c2x $w 0]}]
	    set xpos [expr {$xpos+5*$wid}]
	}
    } else {
	set xpos $cmd
    }

    set wid [expr {[set ${w}(x2)] - [set ${w}(x1)]}]
    set ${w}(x1) $xpos
    set ${w}(x2) [expr {$xpos+$wid}]

    # Update scrollbar
    $sbar set \
	[expr {[set ${w}(x1)]/double([set ${w}(length)])}] \
	[expr {[set ${w}(x2)]/double([set ${w}(length)])}]

    redraw_plot $w
}

# X zoom only via scalebox
proc set_xzoom {w val} {
    global $w

    # Change the scale
    set wid [set ${w}(pwidth)]
    set scale [set ${w}(xzoom)]
    set mid [c2x $w [expr {$wid/2.0}]]
    set ${w}(xzoom) [expr {pow(($val+4)/10,4)}]

    # Adjust xorigin to ensure $mid stays the same
    set ${w}(xorigin) [expr {$mid - ($wid/2)*[set ${w}(xzoom)]}]

    # Reset x1/x2
    set ${w}(x1) [c2x $w 0]
    set ${w}(x2) [c2x $w $wid]

    # Simulate a scrollbar movement to update size, position and redraw.
    set sbar $w.xscroll

    $sbar set \
	[expr {[set ${w}(x1)]/double([set ${w}(length)])}] \
	[expr {[set ${w}(x2)]/double([set ${w}(length)])}]

    redraw_plot $w
}

# Zoom the plot. z>1 => magnify/zoom in, z<1 => zoom out
# x and y are booleans to govern whether we want to zoom in x, y or both.
proc zoom1.5 {w x y z} {
    global $w

    if {$y} {
	set ${w}(yzoom) [expr {[set ${w}(yzoom)]/$z}]
    }

    if {$x} {
	#set ${w}(xzoom) [expr {$z*[set ${w}(xzoom)]}]
	set x1 [set ${w}(x1)]
	set x2 [set ${w}(x2)]
	set mid [expr {($x1+$x2)/2.0}]
	set wid [expr {$x2-$x1+1}]
	set wid [expr {$wid * $z}];
    
	puts "old x=$x1..$x2"

	set x1 [expr {int($mid-$wid/2)}]
	set x2 [expr {int($mid+$wid/2)}]
	set ${w}(x1) $x1
	set ${w}(x2) $x2

	puts "new x=$x1..$x2"

	$w.xscroll set \
	    [expr {$x1/double([set ${w}(length)])}] \
	    [expr {$x2/double([set ${w}(length)])}]

	set ${w}(xorigin) [expr $x1]
	set ${w}(xzoom)   [expr {($x2-$x1+1)/double([set ${w}(pwidth)])}]
    }

    redraw_plot $w seq_depth
    return
}

#
# Functions to convert x/y base-level coordinates into canvas coordinates (c)
# and back again.
#
proc x2c {w pos} {
    global $w
    return [expr {($pos - [set ${w}(xorigin)]) / [set ${w}(xzoom)]}]
}

proc c2x {w pos} {
    global $w
    return [expr {$pos * [set ${w}(xzoom)] + [set ${w}(xorigin)]}]
}

# m* is the size of the margins.
# +ve = from that edge.
# -ve = from opposite edge
proc add_plot {w func height has_scroll has_scale args} {
    puts [info level [info level]]

    global $w
    incr ${w}(ntracks)
    set tnum [set ${w}(ntracks)]
    lappend ${w}(tracks) $tnum
    set t $w.track$tnum
    global $t

    set ${t}(w)      $w
    set ${t}(num)    $tnum
    set ${t}(func)   $func
    set ${t}(height) $height

    set weight 0    
    
    # wants scroll bars
    if {$has_scroll} { 
    	set weight 1
	
 	set ${t}(ys) $w.yscroll$tnum
	scrollbar $w.yscroll$tnum -orient vert -command "yscroll_plot $w $t"
	grid $w.yscroll$tnum -row $tnum -column 2 -sticky nsew
    }
    
    # add a Y scale to the left hand side
    if {$has_scale} {
    	set weight 1
	
	yscale_init $w $t $tnum $height
	grid $w.yscale$tnum -row $tnum -column 0 -sticky nsew
    }	    

    set ${t}(y1) 0
    set ${t}(y2) $height

    set ${t}(canvas) [eval canvas $t -height $height $args \
			  [list -yscrollcommand [list yscroll_plot_set $w $t]]]
    grid $t -row $tnum -column 1 -sticky nsew
    grid rowconfigure $w $tnum -weight $weight
}

#
# TEMP testing 
#
proc add_separator {w height} {
    puts [info level [info level]]
    
    global $w
    incr ${w}(ntracks)
    set tnum [set ${w}(ntracks)]
    
    grid rowconfigure $w $tnum -minsize $height -weight 0
}





proc yscale_init {w t tnum height} {
    global $w $t
    
    # create a scale for the y axis
    canvas $w.yscale$tnum -height $height -width 60
    set ${t}(yscale) $w.yscale$tnum
    set ys [set ${t}(yscale)]
    
    $ys create rectangle 0 0 60 $height -fill \#ffcccc
}

proc yscale_resize {ys height} {
    $ys delete all
    $ys create rectangle 0 0 60 $height -fill \#ffcccc
}    


proc yscale_seq {w t height} {
    global $w $t
    
    set ys [set ${t}(yscale)]
    yscale_resize $ys $height

    set w1 [expr {0.1 * 60}]
    set w2 [expr {0.2 * 60}]
    set w3 [expr {0.3 * 60}]
    set w4 [expr {0.4 * 60}]
    set w5 [expr {0.45 * 60}]
    set w6 [expr {0.7 * 60}]

    set td [set ${t}(track)]
    set zoom [$td yz]
    
    if {$zoom == 0} {
    	set zoom 1
    }

    set ymode [$td configure -ymode]
    set ylog  [set ${t}(YLog)]
    set strands [set ${t}(SeparateStrands)]
    
    # vertical base line
    $ys create line $w1 0 $w1 $height
    
    set y1 [set ${t}(y1)]
    
    if {$ylog && $ymode != 1} {
	set max $height
	set min 0

	set i    $min
	set val  $min
	incr max [expr {abs($y1)}]
	set jump 1
	set step 10
    	set yoffset [set ${t}(YOffset)]
	
	if {$strands} {
	    set midway [expr {($height / 2) - $y1}]
   	    $ys create line $w1 $midway $w6 $midway
	} else {
    	    set midway [expr {int(-$y1)}]
	}
	
	while {$i < $max} {
 
 	    set idown   [expr {$i + $midway}]
	    set iup [expr {$midway - $i}]
	    
    	    $ys create line $w1 $idown $w3 $idown
    	    $ys create line $w1 $iup $w3 $iup
	    
	    if {$val == $step} {
		set jump $step
		set step [expr {$step * 10}]

		foreach {fmt num} [number_format $val] break

   		$ys create line $w1 $idown $w4 $idown
		$ys create text $w5 $idown -anchor w -text [format $fmt $num]

    	    	$ys create line $w1 $iup $w4 $iup
		$ys create text $w5 $iup -anchor w -text [format $fmt $num]
	    }

	    incr val $jump

    	    if {$val < 0} {
		return
	    } 

	    set new [value_to_log $val $zoom $yoffset]

	    set i [expr {int($new)}]
	}
    } elseif {!$ylog && $ymode != 1} {
	# zoom   - zoom value (default 0.5)
	# y1     - offset (canvas coords)
	# height - window height (canvas)
	
	# min - lower value on screen (world)   = y1 * (1 / zoom)
	# max - highest value on screen (world) = (height + y1) * (1 /zoom)
	 
	# initially 0 point will be at height / 2 for separate strands
	
	set offset      0
	set down_strand 0
	set up_strand   0
	    
	if {$strands} {
            set offset [expr {($height / 2) - $y1}]
	    
	    if {$offset <= 0} {
		set down_strand 1
	    } elseif {$offset > $height} {
		set up_strand 1
	    } else {
		set down_strand 1
		set up_strand 1
	    }
	} else {
	    set down_strand 1
	}
	
	set zoom [expr {1 / $zoom}]
	
	if {$strands} {
	    set midway [expr {$height / 2}]
	
	    if {$down_strand && $up_strand} {
		set min 0
		set max [expr {$height * $zoom}]
	    } elseif {$down_strand} {
		set min [expr {($y1 - $midway) * $zoom}]
		set max [expr {($midway + $y1) * $zoom}]
	    } else {
	    	# up strand
		set max [expr {($midway - $y1) * $zoom}]
		set min [expr {($midway - $y1 - $height) * $zoom}]
	    }
	} else {
	    set offset [expr {$y1 * -1}]
	    set min [expr {$y1 * $zoom}]
	    set max [expr {($height + $y1) * $zoom}]
	}
	
        foreach {min step p1 pstep fmt} [nice_num $min [expr {$max - $min + 1}] 0] break
	set dist [expr {$step / $zoom}]
	
	if {$down_strand} {

	    for {set i $min; set j $p1} {$i < $max} {set i [expr {$i + $step}]; set j [expr {$j + $pstep}]} {
    		set yi [expr {($i / $zoom) + $offset}]

		$ys create line $w1 $yi $w4 $yi

		if {$dist < 200 && $j != 0} {
		    $ys create text $w5 $yi  -anchor w -text [format $fmt $j]
		}
	    }

	    set step2 [expr {$step / 2.0}]
	    set pstep2 [expr {$pstep / 2.0}]
	    for {set i $min; set j $p1} {$i < $max} {set i [expr {$i + $step2}]; set j [expr {$j + $pstep2}]} {
    		set yi [expr {($i / $zoom) + $offset}]

		$ys create line $w1 $yi $w3 $yi

		if {$dist >= 200} {
		    $ys create text $w5 $yi  -anchor w -text [format $fmt $j]
		}
	    }

	    set step3 [expr {$step / 10.0}]
	    for {set i $min} {$i < $max} {set i [expr {$i + $step3}]} {
    		set yi [expr {($i / $zoom) + $offset}]

		$ys create line $w1 $yi $w2 $yi
	    }
	}

	if {$up_strand} {

	    for {set i $min; set j $p1} {$i < $max} {set i [expr {$i + $step}]; set j [expr {$j + $pstep}]} {
    		set yi [expr {$offset - ($i / $zoom)}]

		$ys create line $w1 $yi $w4 $yi

		if {$dist < 200 && $j != 0} {
		    $ys create text $w5 $yi  -anchor w -text [format $fmt $j]
		}
	    }

	    set step2 [expr {$step / 2.0}]
	    set pstep2 [expr {$pstep / 2.0}]
	    for {set i $min; set j $p1} {$i < $max} {set i [expr {$i + $step2}]; set j [expr {$j + $pstep2}]} {
    		set yi [expr {$offset - ($i / $zoom)}]

		$ys create line $w1 $yi $w3 $yi

		if {$dist >= 200} {
		    $ys create text $w5 $yi  -anchor w -text [format $fmt $j]
		}
	    }

	    set step3 [expr {$step / 10.0}]
	    for {set i $min} {$i < $max} {set i [expr {$i + $step3}]} {
    		set yi [expr {$offset - ($i / $zoom)}]

		$ys create line $w1 $yi $w2 $yi
	    }
	}
	    	
    } 
}


proc number_format {val} {
    set level 0
    
    set e [expr {int(log10($val) + 1.0e-10)}]
    set i [expr {int(pow(10,$e))}]
    
    set val $i

    while {$val >= 1000} {
	incr level
	set val [expr {$val/1000.0}]
    }
    set suffix [lindex {{} K M G T P E Z} $level]
    set len [string length [regsub {(.*\.(0$)?)|(^[^.]*$)} [expr {$val/2}] {}]]
    set format "%-.${len}f$suffix $val"
    
    return $format
}


proc value_to_log {y zoom yoffset} {
    return [expr {(((50 * log(($y + 1))) + 50 - $yoffset) * $zoom)}]
}


proc log_to_value {y yz yoffset} {
    # from the template display code
    set tmp [expr {((($y / $yz) - 50 + $yoffset) / 50)}]
    set new_y [expr {exp($tmp) - 1}]

    return $new_y
}


# Picks a nice number for val. Ie value to only a few decimal places
proc nice_num {v r p} {
    set v [expr {int($v)}]
    set e [expr {int(log10($r) + 1.0e-10) - $p}]
    set i [expr {int(pow(10,$e))}]
    set n [expr {($v/$i)*$i}]

    # Construct sensible ways of printing the numbers too
    set n2 $n
    set i2 $i
    set level 0
    while {$i2 >= 1000} {
	incr level
	set i2 [expr {$i2/1000.0}]
	set n2 [expr {$n2/1000.0}]
    }
    set suffix [lindex {{} K M G T} $level]
    set len [string length [regsub {(.*\.(0$)?)|(^[^.]*$)} [expr {$i2/2}] {}]]
    set format "%-.${len}f$suffix"
    
    return [list $n $i $n2 $i2 $format]
}


proc yscale_depth {w t height} {
    global $w $t
    
    set ys [set ${t}(yscale)]
    yscale_resize $ys $height
    
    set w1 [expr {0.1 * 60}]
    set w2 [expr {0.2 * 60}]
    set w3 [expr {0.3 * 60}]
    set w4 [expr {0.45 * 60}]
    set w5 [expr {0.4 * 60}]

    set td [set ${t}(track)]

    set zoom [$td yz]
    
    if {$zoom == 0} {
    	set zoom 1
    }
    
    set zoom [expr {1 / $zoom}] 
    
    # This works okay for resizing depth, not for fixed/manual zoom
    set max [expr {$height * $zoom}]
    set min 0

    foreach {min step p1 pstep fmt} [nice_num $min [expr {$max - $min + 1}] 0] break

    # draw scale from the bottom up (0 bottom), clumsy though
    
    set dist [expr {$step / $zoom}]
    
    $ys create line $w1 0 $w1 $height
    
    for {set i $min; set j $p1} {$i < $max} {set i [expr {$i + $step}]; set j [expr {$j + $pstep}]} {
    	set yi [expr {$height - ($i / $zoom)}]
	$ys create line $w1 $yi $w5 $yi

	if {$dist < 200} {
	    $ys create text $w4 $yi -anchor w -text [format $fmt $j]
	}
    }

    set step2 [expr {$step / 2.0}]
    set pstep2 [expr {$pstep / 2.0}]
    
    for {set i $min; set j $p1} {$i < $max} {set i [expr {$i + $step2}]; set j [expr {$j + $pstep2}]} {
    	set yi [expr {$height - ($i / $zoom)}]
 	$ys create line $w1 $yi $w3 $yi
	
	if {$dist >= 200} {
	    $ys create text $w4 $yi -anchor w -text [format $fmt $j]
	}
    }
    
    set step3 [expr {$step / 10.0}]
    
    for {set i $min} {$i < $max} {set i [expr {$i + $step3}]} {
    	set yi [expr {$height - ($i / $zoom)}]
	$ys create line $w1 $yi $w2 $yi
    }
} 


proc yscroll_plot {w t cmd {opt1 {}} {opt2 {}}} {
    global $w $t
    puts [info level [info level]]
    
    set td   [set ${t}(track)]
    foreach {top bottom} [$td yrange] break
    
    if {[set ${t}(all_visible)] == 1} {
    	return; # no scroll needed
    }

    set y1 [set ${t}(y1)]
    set h [winfo height [set ${t}(canvas)]]
    set ymax [set ${t}(ymax)]
    
    switch $cmd {
	scroll {
	    if {$opt2 == "pages"} {
		return; # we don't do pages
	    } else {
	    	if {$top <= 0 && $opt1 <= 0} {
		    return
		}
		
		if {$bottom >= 1 && $opt1 >= 0} {
		    return
		}
		
	    	incr y1 [expr {int($opt1*1)}]
	    }
	}

	moveto {
	    if {$top == 0 && $opt1 < $top} {
	    	# scroll bar at top
		return
	    }
	    
	    if {$bottom == 1 && $opt1 > $top} {
	    	# scroll bar at bottom and nowhere else to go
		return
	    }
	
	    set y1 [expr {int($opt1 * [set ${t}(scroll_height)])}]
	}
    }
    
    set ${t}(y1) $y1
    set y2 [expr {$y1+$h}]

    puts "$y1 to $y2 out of [set ${t}(scroll_height)]"

    if {[info exists ${t}(Raster)]} {
	eval [set ${t}(Raster)] yview $cmd $opt1 $opt2
    }
    
    redraw_plot $w [set ${t}(func)]
}

proc yscroll_plot_set {w t args} {
    puts [info level [info level]]
}


proc show_track {w track_type height show} {
    global $w
    
    if {$show} {
    	#add_plot $w $track_type $height -bd 0 -relief raised
	show_plot $w $track_type
    } else {
    	remove_plot $w $track_type
    }
}


proc show_plot {w track_type} {
    global $w
    
    foreach id [set ${w}(tracks)] {
    	set t $w.track$id
	global $t
	
	set comp [string compare $track_type [set ${t}(func)]]
	
	if {$comp == 0} {
	    grid $w.yscale$id
	    grid $t
	    grid $w.yscroll$id
	    grid rowconfigure $w $id -weight 1
	}
    }
    
    redraw_plot $w 
}

proc remove_plot {w track_type} {
    global $w
    
    foreach id [set ${w}(tracks)] {
    	set t $w.track$id
	global $t
	
	set comp [string compare $track_type [set ${t}(func)]]

    	# width needs to be set so the raster can be resized
	# there is a width check in the seq_seqs proc (and equivalents)	
    	set ${t}(width)  1

	if {$comp == 0} {
	    grid remove $w.yscale$id
	    grid remove $t
	    grid remove $w.yscroll$id
	    grid rowconfigure $w $id -weight 0
	}
    }
    
    update idletasks
    redraw_plot $w
}


# The top-level redraw function. This works out what region is visible and
# fires off the individual redraw functions registered.
proc redraw_plot {w {track_types {}} args} {
    global $w

    if {$track_types != {}} {
	set tracks {}
	foreach id [set ${w}(tracks)] {
	    set t $w.track$id
	    global $t
	    if {[lsearch -exact $track_types [set ${t}(func)]] != -1} {
		lappend tracks $id
	    }
	}
    } else {
	set tracks [set ${w}(tracks)]
    }

    if {[info exists ${w}(RedrawPending)]} {
	append ${w}(RedrawPending) " $tracks"
	return
    }

    set ${w}(RedrawPending) $tracks
    after idle [list redraw_plot_doit $w]
}


# The actual redraw function, rather than a simple request to redraw
proc redraw_plot_doit {w} {
    global $w

    set bd [set ${w}(border)]

    set tracks [lsort -unique [set ${w}(RedrawPending)]]

    set x1 [set ${w}(x1)]
    set x2 [set ${w}(x2)]
    set ${w}(xorigin) [expr $x1]

#    puts ""
#    puts ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
#    puts "xorigin=[set ${w}(xorigin)] xzoom=[set ${w}(xzoom)] yzoom=[set ${w}(yzoom)]"
#    puts "base coord=$x1..$x2   canvas coord=[x2c $w $x1]..[x2c $w $x2]"
#    puts ""

    $w configure -cursor watch

    foreach id $tracks {
	set t $w.track$id
	global $t

	set d [set ${t}(canvas)]
	set y1 [set ${t}(y1)]
	set y2 [expr {[set ${t}(y1)]+[winfo height [set ${t}(canvas)]]}]
	puts [set ${t}(func)]:[time {[set ${t}(func)] $w $t $x1 $x2 $y1 $y2}]
    }
    
    $w configure -cursor {}
    

    set ${w}(last_x1) $x1
    set ${w}(last_x1) $x2

    catch {unset ${w}(RedrawPending)}
}


#
# range_track_init
#
# Common initilizations for range based
# tracks.
#

proc range_track_init {w t} {
    global $w $t
    
    set ${t}(Init) 1
    set ${t}(Accurate) 0
    set ${t}(YLog) 1
    set ${t}(YScale) 100
    set ${t}(OldYScale) [set ${t}(YScale)]
    set ${t}(YOffset) 50
    set ${t}(Simple) 0
    set ${t}(Y) "Template Size"
    set ${t}(Colour) "Combined mapping quality"
    set ${t}(Spread) 0
    set ${t}(PlotDepth) 0
    set ${t}(SeparateStrands) 0
    set ${t}(MinQual) 0
    set ${t}(MaxQual) 255
    set ${t}(MinYSize) 1024
    set ${t}(FilterPair) 0
    set ${t}(FilterConsistent) 0
    set ${t}(FilterSpanning) 0
    set ${t}(ReadsOnly) 0
    set ${t}(m_start) -1
    set ${t}(m_stop)  -1
    
    # Add raster component to canvas
    update idletasks

    set d    [set ${t}(canvas)]
    puts 1:[winfo width $d]x[winfo height $d]
    set r [raster $d.r \
	       -width [winfo width $d] \
	       -height [winfo height $d] \
	       -bg black \
	       -bd 0]

    set rid [$d create window 0 0 -anchor nw -window $r]
    set ${t}(Raster) $r
    set ${t}(R_id)   $rid
    $d raise [set ${t}(R_id)]

    $r world_scroll 1 1 [set ${w}(length)] 10000
    $r world 0 0 [expr {[set ${w}(pwidth)]*[set ${w}(xzoom)]}] \
	[winfo height $d]

    set ${t}(R_zoom) [set ${w}(xzoom)]
    
    bind $r <B1-Motion> "drag_x $w $t %x %y"
    bind $r <B1-ButtonRelease> "end_drag_x $w $t"
}
   
###############################################################################
# The tracks
###############################################################################

#
# Initialises the template display window - range track
#
proc seq_seqs_init {w t} {
    global $w $t

    puts START:[info level [info level]]

    range_track_init $w $t

    set r [set ${t}(Raster)]

    set ${t}(track) [g5::template_display \
			 -io 	 [set ${w}(io)] \
			 -cnum   [set ${w}(cnum)] \
			 -raster $r \
			 -range  [set ${w}(grange)]]
    set td   [set ${t}(track)]
	
    bind $r <Any-Motion> "
         foreach {x y} \[$td xhair %x %y\] break;
         set ${w}(info) \"X: \$x    Y: \$y)\"
    "
    bind $r <Any-Leave> "$td xhair; set ${w}(info) {}"

    # Add GUI elements
    top_controls    $w $t
    bottom_controls $w $t

    set d    [set ${t}(canvas)]
    $d bind all <Any-Enter> "seq_seqs_bind $t $w $d"

    bind $r <<use>> "invoke_editor $w $t %x"
    
    puts END:[info level [info level]]
}


proc bottom_controls {w t} {
    global $w $t
    
    set bc $w.bcontrol

    scale $bc.yzoom -from 1 -to 1000 -orient horiz -label "Y Magnification" \
	-variable ${t}(YScale) -command "scale_Y $w $t seq_seqs"

    scale $bc.yspread -from 0 -to 250 -orient horiz -label "Y Spread" \
	-variable ${t}(Spread) -command "redraw_plot $w seq_seqs"

    scale $bc.yoffset -from 0 -to 250 -orient horiz -label "Y Offset" \
	-variable ${t}(YOffset) -command "redraw_plot $w seq_seqs"

    pack $bc.yzoom $bc.yspread $bc.yoffset -fill both -expand 1 -side left

    scale $bc.minysize -from 100 -to 5000 -resolution 10 -orient horiz \
	-label "Stacking Y Size" -variable ${t}(MinYSize) \
	-command "redraw_plot $w seq_seqs"

    pack $bc.minysize -fill both -expand 1 -side left
}


proc top_controls {w t} {
    global $w $t
    
    set f [frame $w.controls]
    grid $f -column 1 -row 0 -sticky w

    button $f.template -text Template \
    	-command "template_dialog $w $t"
    button $f.filter -text Filter \
	-command "seq_seqs_filter $w $t"

    pack $f.template $f.filter -side left

#    tk_optionMenu $f1.y ${t}(Y) \
#	{Template Size} \
#	Stacking \
#	{Mapping quality}
#    tk_optionMenu $f1.col ${t}(Colour) \
#	{Combined mapping quality} \
#	{Minimum mapping quality} \
#	{Maximum mapping quality} \
#	{Reads}
#    trace add variable ${t}(Y) write "redraw_plot $w seq_seqs"
#    trace add variable ${t}(Colour) write "redraw_plot $w seq_seqs"
#    pack $f1.y $f1.col -side left
#    
#    checkbutton $f1.acc -text ">>Acc" -variable ${t}(Accurate) \
#	-command "redraw_plot $w seq_seqs"
#    checkbutton $f1.reads -text "Reads" -variable ${t}(ReadsOnly) \
#	-command "redraw_plot $w seq_seqs"
#    checkbutton $f1.log -text "Y-Log scale" -variable ${t}(YLog) \
#	-command "redraw_plot $w seq_seqs"
#    checkbutton $f1.sep_strands -text "Separate strands" \
#	-variable ${t}(SeparateStrands) \
#	-command "redraw_plot $w seq_seqs"
#    checkbutton $f1.plot_depth -text "Depth" \
#	-variable ${t}(PlotDepth) \
#	-command "redraw_plot $w seq_seqs"
#    pack $f1.y $f1.col $f1.template $f1.filter $f1.log $f1.acc $f1.reads $f1.sep_strands \
#	$f1.plot_depth -side left
}
    


#
# Plots the sequence read depth.
#
proc seq_seqs {w t x1 x2 y1 y2} {
    global $w $t

    if {![info exists ${t}(Init)]} {
	seq_seqs_init $w $t
    }

    set r [set ${t}(Raster)]
    set d [set ${t}(canvas)]

    if {[info exists ${t}(width)]} {
    	update idletasks
	
	$r configure -width [winfo width $d] -height [winfo height $d]
    }

    # Replacement for the scroll/zooming code - just reset world
    foreach {X1 Y1 X2 Y2} [$r world] break;
    # the raster coords can become fractions, which doesn't make much sense
    set Y1 [expr {int($Y1)}]
    # change the raster wy1 coord to match visible height
    set nheight [winfo height $d] 
    set Y2 [expr {$Y1 + $nheight}]
    
    $r world $x1 $Y1 $x2 $Y2
    
    set ${t}(y1) $Y1

    set cmode [lsearch {{Combined mapping quality} \
			    {Minimum mapping quality} \
			    {Maximum mapping quality} \
			    {Reads}} [set ${t}(Colour)]]
    set ymode [lsearch {{Template Size} Stacking {Mapping quality}} \
		   [set ${t}(Y)]]

    update idletasks

    # Forces a redraw too
    set flag [expr { [set ${t}(FilterPair)]
		    +[set ${t}(FilterConsistent)]
		    +[set ${t}(FilterSpanning)]}]
		     
    1.5redraw_cursor $w $t

    set td   [set ${t}(track)]

    $td configure \
	-accuracy    [set ${t}(Accurate)] \
	-logy        [set ${t}(YLog)] \
	-yoffset     [set ${t}(YOffset)] \
	-spread      [set ${t}(Spread)] \
	-ymode       $ymode \
	-cmode       $cmode \
	-reads_only  [set ${t}(ReadsOnly)] \
	-by_strand   [set ${t}(SeparateStrands)] \
	-filter      $flag \
	-yzoom       [set ${t}(YScale)] \
	-plot_depth  [set ${t}(PlotDepth)] \
 	-min_qual    [set ${t}(MinQual)] \
 	-max_qual    [set ${t}(MaxQual)] \
	-min_y_size  [set ${t}(MinYSize)]


    # Set Y scrollbar
    set ymax [set ${t}(ymax) [$td ymax]]
    set ymin [set ${t}(ymin) [$td ymin]]
    
    if {[string compare [$td yrange] "0.000000 1.000000"] == 0} {
    	set ${t}(all_visible) 1
    } else {
    	set ${t}(all_visible) 0
    }
 
    
    set ${t}(scroll_height) [expr {[$td ymax]-[$td ymin]}]
    eval [set ${t}(ys)] set [$td yrange]
    
    # change the yscale
    yscale_seq $w $t $nheight
    
    if {$ymax <= $Y1} {
    	puts "Out of range (max) $ymax (y1) $Y1 (y2) $Y2"

	# check for no data
	if {$ymax == 0 && $Y1 == 0} {
	    return
	}
	
	# bring the range back to normal range
	set diff [expr {$ymax - $Y2}]
	
	incr Y2 $diff
	incr Y1 $diff
	
	if {$Y1 < 0} {
	    set Y1 0
	}
	
	if {$Y1 == $Y2} {
	    incr Y2; # world needs a rectangle
	}
	
    	$r world $x1 $Y1 $x2 $Y2
	
    	# need to redraw and do it now, simple redraw_plot
	# stacks up redraw but does not.
	redraw_plot_doit $w
    }
}

#
# initialise the depth track - a range track
#
proc depth_track_init {w t} {
    global $w $t

    puts START:[info level [info level]]

    range_track_init $w $t

    set r [set ${t}(Raster)]
    set ${t}(track) [g5::depth_track \
			 -raster $r \
			 -range  [set ${w}(grange)]]

    set dt   [set ${t}(track)]
    bind $r <Any-Motion> "
         foreach {x y} \[$dt xhair %x %y\] break;
         set ${w}(info) \"X: \$x    Y: \$y)\"
    "
    bind $r <Any-Leave> "$dt xhair; set ${w}(info) {}"

    puts END:[info level [info level]]
}


#
# plot depth track - a range track
#
proc depth_track {w t x1 x2 y1 y2} {
    global $w $t

    if {![info exists ${t}(Init)]} {
	depth_track_init $w $t
    }

    set r    [set ${t}(Raster)]
    set d    [set ${t}(canvas)]

    if {[info exists ${t}(width)]} {
	$r configure -width [winfo width $d] -height [winfo height $d]
    }

    # Replacement for the scroll/zooming code - just reset world
    foreach {X1 Y1 X2 Y2} [$r world] break;

    set Y1 [expr {int($Y1)}]
    # change the raster wy1 coord to match visible height
    set nheight [winfo height $d] 
    set Y2 [expr {$Y1 + $nheight}]

    $r world $x1 $Y1 $x2 $Y2

    update idletasks

    # Forces a redraw too
    1.5redraw_cursor $w $t

    set td   [set ${t}(track)]
    $td replot

   # Set Y scrollbar
    if {[string compare [$td yrange] "0.000000 1.000000"] == 0} {
    	set ${t}(all_visible) 1
    } else {
    	set ${t}(all_visible) 0
    }

    set ${t}(ymax) [$td ymax]
    set ${t}(scroll_height) [expr {[$td ymax]-[$td ymin]}]
    
    eval [set ${t}(ys)] set [$td yrange]

    # change the yscale
    yscale_depth $w $t [winfo height $d]
}



#
# Plots the sequence ruler
#
proc seq_ruler {w t x1 x2 y1 y2} {
    global $w $t
    
    set d    [set ${t}(canvas)]

    set wx1 [x2c $w $x1]
    set wx2 [x2c $w $x2]

    set h  [expr {$y2-$y1}]
    set h1 [expr {$y1+0.9*($y2-$y1)}]
    set h2 [expr {$y1+0.1*($y2-$y1)}]
    set h3 [expr {$y1+0.7*($y2-$y1)}]
    set h4 [expr {$y1+0.3*($y2-$y1)}]
    set h5 [expr {$y1+0.5*($y2-$y1)}]

    $d delete all

    $d create rectangle $wx1 0 $wx2 $h -fill \#ffcccc

    $d create line $wx1 $h3 $wx2 $h3 -fill black

    foreach {x1 step p1 pstep fmt} [nice_num $x1 [expr {$x2-$x1+1}] 0] break
    set dist [expr {$step/[set ${w}(xzoom)]}]
    for {set i $x1; set j $p1} {$i < $x2} {set i [expr {$i+$step}]; set j [expr {$j+$pstep}]} {
    	set xi [x2c $w $i]
	$d create line $xi $h2 $xi $h3
	if {$dist < 200} {
	    $d create text $xi $h3 -anchor nw \
		-text [format $fmt $j]
	}
    }

    set step2 [expr {$step/2.0}]
    set pstep2 [expr {$pstep/2.0}]
    for {set i $x1; set j $p1} {$i < $x2} {set i [expr {$i+$step2}]; set j [expr {$j+$pstep2}]} {
    	set xi [x2c $w $i]
	$d create line $xi $h4 $xi $h3
	if {$dist >= 200} {
	    $d create text $xi $h3 -anchor nw \
		-text [format $fmt $j]
	}
    }

    set step3 [expr {$step/10.0}]
    for {set i $x1} {$i < $x2} {set i [expr {$i+$step3}]} {
    	set xi [x2c $w $i]
	$d create line $xi $h5 $xi $h3
    }
}


#
# functions used by tracks
#

proc drag_x {w t x y} {
    global $w $t
    
    set start [set ${t}(m_start)]
    
    if {$start == -1} {
    	set ${t}(m_start) $x
	$w configure -cursor fleur
    } else {
    	set ${t}(m_stop) $x
    }

    set td [set ${t}(track)]
    $td xhair $x $y
}

proc end_drag_x {w t} {
    global $w $t
    
    set start [set ${t}(m_start)]
    set end   [set ${t}(m_stop)]
    
    if {$start != -1 && $end != -1} {
	
	set c_end   [c2x $w $start]
	set c_start [c2x $w $end]
	
	set c_diff [expr {$c_end - $c_start}]
	
	set x1 [set ${w}(x1)]
	set x2 [set ${w}(x2)]
	
	set ${w}(x1) [expr {$x1 + $c_diff}] 	
	set ${w}(x2) [expr {$x2 + $c_diff}]
	
	# set the scroll bar
	set sbar $w.xscroll
    	$sbar set \
	[expr {[set ${w}(x1)]/double([set ${w}(length)])}] \
	[expr {[set ${w}(x2)]/double([set ${w}(length)])}]
    } else {
    	return
    }
    
    $w configure -cursor {}
    
    set ${t}(m_start) -1
    set ${t}(m_stop) -1
    
    redraw_plot $w
}


proc scale_Y {w t type mag} {
    global $w $t

    set d [set ${t}(canvas)]
    set r [set ${t}(Raster)]
    set height [winfo height $d]
    set new_zoom [set ${t}(YScale)]
    set old_zoom [set ${t}(OldYScale)] 
    
    # get the world coords
    foreach {X1 Y1 X2 Y2} [$r world] break;
    
    set wy_mid [expr {(($Y2 - $Y1) / 2) + $Y1}]
    set new_wy_mid [expr {($wy_mid * $new_zoom) / $old_zoom}]
    
    set diff [expr {int(($new_wy_mid - $wy_mid))}]
    incr Y1 $diff

    if {$Y1 < 0} {
	set Y1 0
    } 

    set Y2 [expr {$Y1 + $height}]
    $r world $X1 $Y1 $X2 $Y2
    
    set ${t}(OldYScale) $new_zoom
    redraw_plot $w $type
}


# Mouse-over event callback
proc seq_seqs_bind {t w canvas} {
    global $w $t

    # Find the current reading record number
    set tags [$canvas gettags current]
    $canvas raise current
    set rec [lindex [regexp -inline {(^| )rec_([0-9]+)} $tags] 2]
    if {$rec == ""} return

    # Produce some information about it
    set io [set ${w}(io)]
    set r1 [$io get_sequence $rec]
    set c1 [$r1 get_contig]
    set p1 [$r1 get_position]

    set info "[$r1 get_name] len [$r1 get_len] mq [$r1 get_mapping_qual]"

    set pair [$r1 get_pair]
    if {$pair} {
	set r2 [$io get_sequence $pair]
	set c2 [$r2 get_contig]
	set p2 [$r2 get_position]
	set l2 [$r2 get_length]
	append info ", pair [$r2 get_name] len $l2 mq [$r2 get_mapping_qual]"
	if {$c1 != $c2} {
	    append info " (contig\#$c2)"
	} else {
	    set size [expr {$p2+abs($l2)-$p1+1}]
	    append info " (insert size $size)"
	}
    } else {
	append info " single ended (rec \#$rec)"
    }

    set ${w}(info) $info
}

#
# Brings up a contig editor
#
proc invoke_editor {w t x} {
    global $w $t

    # Hunt for an existing editor for this contig
    set io [set ${w}(io)]
    set found 0
    foreach result [result_names -io $io] {
	foreach {crec id name} $result break;
	if {$name == "Contig Editor" && $crec == [set ${w}(cnum)]} {
	    set found $id
	    break
	}
    }

    set r    [set ${t}(Raster)]
    set x [lindex [$r toworld $x 0] 0]
    if {$found} {
	global $t.Id2Cid
	if {[info exists $t.Id2Cid($id)]} {
	    set cid [set $t.Id2Cid($id)]
	} else {
	    set cid 0
	}

	contig_notify -io $io -cnum [set ${w}(cnum)] -type CURSOR_NOTIFY \
	    -args [list id      $cid \
		        job     MOVE \
		        seq     [set ${w}(cnum)] \
		        abspos  $x \
		        pos     $x \
		        sent_by 0]
    } else {
	puts x=$x
	edit_contig \
	    -io $io \
	    -contig  [set ${w}(cnum)] \
	    -reading [set ${w}(cnum)] \
	    -pos $x
    }
}


#
# track dialog boxes
#


#
# Contols the filter dialogue
#
proc seq_seqs_filter {w t} {
    global $w $t

    set f [xtoplevel $w.filter_win -resizable 0]

    if {$f != ""} {
    	wm title $f "Filter"
    } else {
    	return
    }
	 
    set f2 [frame $f.sub -bd 2 -relief groove]

    # Initialise variable copies
    foreach v {FilterPair FilterConsistent FilterSpanning MinQual MaxQual} {
	set ${t}(_$v) [set ${t}($v)]
	set ${t}($v~) [set ${t}($v)]
    }

    # Auto-update?
    checkbutton $f2.auto \
	-text "Auto update" \
	-variable ${t}(FilterAutoUpdate) \
	-command "seq_seqs_filter_update $w $t $f"
    grid $f2.auto - - - -sticky w

    # Radio buttons
    set r 1
    foreach {label var list} {
	Pairs       FilterPair {
	    All 0
	    {Single only} 2
	    {Pair only} 1
	} \
	Consistency FilterConsistent {
	    All 0
	    {Consistent only} 4
	    {Inconsistent only} 8
	} \
        Spanning    FilterSpanning {
	    All 0
	    {Non-spanning only} 32
	    {Spanning only} 16
	}
    } {
	set c 1
	grid [label $f2.label$r -text "$label:"] -row $r -column 0 -sticky w

	foreach {name value} $list {
	    grid [radiobutton $f2.b$r$c \
		      -text $name \
		      -value $value \
		      -variable ${t}(_$var) \
		      -command "seq_seqs_filter_update $w $t $f" \
		     ] -row $r -column $c -sticky w
	    incr c
	}

	incr r
    }

#    foreach {name value} {All 0 {Single only} 2 {Pair only} 1} {
#	grid [radiobutton $f2.b$r$c \
#		  -text $name \
#		  -value $value \
#		  -variable ${t}(_FilterPair)] -row $r -column $c -sticky w
#	incr c
#    }
#
#    incr r
#    grid [label $f2.label$r -text "Consistency:"] -row $r -column 0 -sticky w
#    set c 1
#    foreach {name value} {All 0 {Consistent only} 4 {Inconsistent only} 8} {
#	grid [radiobutton $f2.b$r$c \
#		  -text $name \
#		  -value $value \
#		  -variable ${t}(_FilterConsistent)] \
#	    -row $r -column $c -sticky w
#	incr c
#    }
#
#    incr r
#    grid [label $f2.label$r -text "Spanning:"] -row $r -column 0 -sticky w
#    set c 1
#    foreach {name value} {All 0 {Non-spanning only} 32 {Spanning only} 16} {
#	grid [radiobutton $f2.b$r$c \
#		  -text $name \
#		  -value $value \
#		  -variable ${t}(_FilterSpanning)] -row $r -column $c -sticky w
#	incr c
#    }

    # Scale bars
    incr r
    label $f2.min_label -text "Min. Qual"
    scale $f2.min_qual \
	-from 0 \
	-to 255 \
	-orient horiz \
	-variable ${t}(_MinQual) \
	-command "seq_seqs_filter_update $w $t $f min"
    grid $f2.min_label $f2.min_qual - - -sticky ew
    label $f2.max_label -text "Max. Qual"
    scale $f2.max_qual \
	-from 0 \
	-to 255 \
	-orient horiz \
	-variable ${t}(_MaxQual) \
	-command "seq_seqs_filter_update $w $t $f max"
    grid $f2.max_label $f2.max_qual - - -sticky ew

    # Standard controls
    okcancelhelp $f.ok \
	-ok_command "seq_seqs_filter_ok $w $t $f" \
	-cancel_command "seq_seqs_filter_cancel $w $t $f" \
	-apply_command "seq_seqs_filter_apply $w $t $f 1" \
	-bd 2 \
	-relief groove

    pack $f2 $f.ok -side top -fill both -expand 1
}


# Filter::OK 
proc seq_seqs_filter_apply {w t f keep} {
    global $w $t
    foreach var [array names $t] {
	if {[regexp {^_(.*)} $var dummy v]} {
	    set ${t}($v)  [set ${t}($var)]
	    if {$keep} {
		set ${t}($v~) [set ${t}($var)]
	    }
	}
    }

    redraw_plot $w
}


proc seq_seqs_filter_ok {w t f} {
    seq_seqs_filter_apply $w $t $f 1
    destroy $f
}


proc seq_seqs_filter_cancel {w t f} {
    global $w $t
    foreach var [array names $t] {
	if {[regexp {^(.*)~} $var dummy v]} {
	    puts "set ${t}($v) [set ${t}($var)]"
	    set ${t}($v) [set ${t}($var)]
	}
    }

    destroy $f
    redraw_plot $w
}


proc seq_seqs_filter_update {w t f {cmd {}} args} {
    global $w $t

    if {$cmd == "min" && [set ${t}(_MinQual)] > [set ${t}(_MaxQual)]} {
	set ${t}(_MaxQual) [set ${t}(_MinQual)]
    }
    if {$cmd == "max" && [set ${t}(_MaxQual)] < [set ${t}(_MinQual)]} {
	set ${t}(_MinQual) [set ${t}(_MaxQual)]
    }

    if {[set ${t}(FilterAutoUpdate)]} {
	seq_seqs_filter_apply $w $t $f 0
    }
}


proc template_dialog {w t} {
    global $w $t
    
    set c [xtoplevel $w.coverage -resizable 0]
    
    if {$c != ""} {
    	wm title $c "Template Track"
    } else {
    	return
    }
    
    set cf [frame $c.sub -padx 5  -pady 3 -bd 2 -relief groove]

    # auto
    checkbutton $cf.auto \
	-text "Auto update" \
	-variable ${t}(FilterAutoUpdate) \
	-command "seq_seqs_filter_update $w $t $c"
    grid $cf.auto -row 0 -column 0 -columnspan 4 -sticky w
    
     # set up the Y and Colour label frames
    set col 0
    foreach {lflabel suf var lflist} {
    	{Y Postion} l Y {
	    {Template size} {Template Size}
	    Stacking Stacking
	    {Mapping quality} {Mapping Quality}
	} \
	Colour lc Colour {
	    {Combined mapping quality} {Combined mapping quality}
	    {Minimum mapping quality} {Minimum mapping quality}
    	    {Maximum mapping quality} {Maximum mapping quality}
    	    Reads Reads
	}
    } {
    	set r 0
	set ${t}(_$var) [set ${t}($var)]
	labelframe $cf.$suf -text $lflabel
	
	foreach {name value} $lflist {
	    grid [radiobutton $cf.$suf.r$r \
	    	    -text $name \
		    -value $value \
		    -variable ${t}(_$var) \
		    -command "seq_seqs_filter_update $w $t $c" \
		    ] -row $r -column 0 -stick w
	    incr r
	}
	
	# labelframe in main grid 
    	grid $cf.$suf  -row 1 -column $col -columnspan 2 -stick n
	incr col 2 
    }
    
    # view controls
    set col 0
    foreach {blabel var} {
    	{>>Acc} Accurate \
	Reads ReadsOnly \
	{Y-log scale} YLog \
	{Separate strands} SeparateStrands
    } {
    	set ${t}(_$var) [set ${t}($var)]
	
	grid [checkbutton $cf.c$col \
	    	-text $blabel \
		-variable ${t}(_$var) \
		-command "seq_seqs_filter_update $w $t $c" \
		] -row 2 -column $col -stick w -pady {5 0}
		
	incr col
    }
    
    # Ok for changes
    okcancelhelp $c.ok -ok_command "seq_seqs_filter_ok $w $t $c" \
    		       -apply_command "seq_seqs_filter_apply $w $t $c 1" \
    		       -cancel_command "seq_seqs_filter_cancel $w $t $c" \
		       -bd 2 \
		       -relief groove

    pack $cf $c.ok -side top -fill both -expand 1   
}










##
## Plots the sequence read depth.
##
#proc seq_depth {w t x1 x2 y1 y2} {
#    global $w $t
#
#    puts [info level [info level]]
#
#    set c    [set ${w}(contig)]
#    set d    [set ${t}(canvas)]
#    $d delete all
#    $d configure -scrollregion \
#	[list 0 0 \
#	     [expr [winfo width $d]] \
#	     [expr [winfo height $d]]]
#    set wid  [set ${w}(pwidth)]
#    set yz   [set ${w}(yzoom)]
#
#    set inc [expr {($x2-$x1+1)/double($wid)}]
#    if {$inc < 1} {set inc 1}
#
#    set l [$c read_depth [expr {int($x1)}] [expr {int($x2)}] $inc]
#    set len [llength $l]
#    puts ">>> Region $x1..$x2 y=$y1..$y2 in steps of $inc => $len items"
#
#    set h [expr {$y2-$y1+1}]
#
#    set wx1 [x2c $w $x1]
#    set wx2 [x2c $w $x2]
#    set wid [expr {$wx2-$wx1+1}]
#    if {$wid < 0} {set wid 1}
#    set bpv [expr {$len/double($wid)}]
#
#    set line {}
#    set max $h
#    if {$bpv < 1} {
#	for {set i 0} {$i < $len} {incr i} {
#	    set avg [expr {[lindex $l $i]*$yz}]
#	    if {$avg > $max} {set avg $max}
#	    lappend line [expr {$wx1+$i/$bpv}] [expr {$h-$avg}]
#	}
#    } else {
#	for {set i 0} {$i < $len} {incr i} {
#	    set avg [expr {[lindex $l $i]*$yz}]
#	    if {$avg > $max} {set avg $max}
#	    lappend line [expr {$wx1+$i}] [expr {$h-$avg}]
#	}
#    }
#
#    $d create line $line -fill black
#    puts BPV=$bpv/[lrange $line 0 1]..[lrange $line [expr {[llength $line]-4}] end]
#}

#proc piarray {n} {
#    upvar 1 $n a
#
#    foreach name [lsort -integer [array names a]] {
#	puts [format "${n}(%2d) %f" $name $a($name)]
#    }
#}

# Given a list of horizontal lines this distributes the elements around
# so they do not clash in Y
#proc compute_max_y {w lines} {
#    global $w
#
#    set ymax 5
#    set ymap [lindex [lindex $lines 0] 1]
#    for {set i 1} {$i < $ymax} {incr i} {
#	lappend ymap -999999
#    }
#    set ysize 1
#
#    set all_y {}
#
#    set n 1
#    foreach l $lines {
#	# First and last coords represent the range of lines
#	set x1 [lindex $l 1]
#	set x2 [lindex $l end-2]
#
#	while {$x1 < [lindex $ymap [expr {$n%$ymax}]]} {
#	    # Compute new and old modulus
#	    set old_mod [expr {$n%$ymax}]
#	    set new_mod [expr {$n%($ymax+1)}]
#	    set mod_delta [expr {$new_mod - $old_mod}]
#
#	    # Rotate if needed
#	    if {$mod_delta > 0} {
#		set list [concat \
#			    [lrange $ymap [expr {$ymax-$mod_delta}] end] \
#			    [lrange $ymap 0 [expr {$ymax-$mod_delta-1}]]]
#	    } elseif {$mod_delta < 0} {
#		set list [concat \
#			    [lrange $ymap [expr {-$mod_delta}] end] \
#			    [lrange $ymap 0 [expr {-$mod_delta-1}]]]
#	    } else {
#		set list $ymap
#	    }
#
#
#	    # Insert previous item
#	    if {[llength $all_y] > $ymax} {
#		set app [lindex $all_y end-$ymax]
#		set all_y [lreplace $all_y end-$ymax end-$ymax]
#	    } else {
#		set app -99999999
#	    }
#	    set ymap [linsert $list $new_mod $app]
#	    incr ymax
#	}
#
#	lappend all_y [expr {$x2+5}]
#	lset ymap [expr {$n%$ymax}] [expr {$x2+5}]
#	incr n
#    }
#
#    return $ymax
#}




#proc lsort_callback {a b} {
#    set base 4
#    set l1 [expr {int(log([lindex $a end-2]-[lindex $a 0])/log($base))}]
#    set l2 [expr {int(log([lindex $b end-2]-[lindex $b 0])/log($base))}]
#
#    if {$l2 > $l1} {
#	return 1;
#    } elseif {$l1 > $l2} {
#	return -1;
#    }
#
#    # Same approximate length, so sort by position still
#    set p1 [lindex $a 0]
#    set p2 [lindex $b 0]
#    
#    return [expr {$p2>$p1 ? -1 : ($p2<$p1 ? 1 : 0)}]
#}




##############################################################################
#user interface dialog box for reading coverage histogram
proc ReadingCoverage { io } {
    global gap5_defs

    set f [keylget gap5_defs READING_COVERAGE.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "reading coverage"

    contig_id $f.id \
	    -io $io \
	    -range 0
    
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "ReadingCoverage2 $io $f $f.id" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Consistency-ReadingCov}" \
	    -bd 2 \
	    -relief groove

    pack $f.id $f.ok_cancel -side top -fill both
}

proc ReadingCoverage2 { io f id} {
    
    if {[set contign [contig_id_gel $id]] == ""} {bell; return}
    SetContigGlobals $io $contign

    # stop windows from hiding the plot
    destroy $f


    set pwin .read_depth[counter]
    set cname [lindex $contign 0]
    1.5Dplot $pwin $io 800 300 [db_info get_contig_num $io $cname]
    add_plot $pwin seq_depth 150  -bd 2 -relief raised
    add_plot $pwin seq_ruler 50  -bd 2 -relief raised
}

##############################################################################
#user interface dialog box for reading coverage histogram
proc TemplateDisplay { io } {
    global gap5_defs

    set f [keylget gap5_defs READING_COVERAGE.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "reading coverage"

    contig_id $f.id \
	    -io $io \
	    -range 0
    
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "TemplateDisplay2 $io $f $f.id" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Consistency-ReadingCov}" \
	    -bd 2 \
	    -relief groove

    pack $f.id $f.ok_cancel -side top -fill both
}

proc TemplateDisplay2 { io f id} {
    
    if {[set contign [contig_id_gel $id]] == ""} {bell; return}
    SetContigGlobals $io $contign

    # stop windows from hiding the plot
    destroy $f

    set cname [lindex $contign 0]
    set cnum [db_info get_contig_num $io $cname]
    CreateTemplateDisplay $io $cnum
}

# THIS ONE
proc CreateTemplateDisplay {io cnum} {
    set pwin .read_depth[counter]
    1.5Dplot $pwin $io 900 600 $cnum
    add_plot $pwin seq_seqs    250  1 1 -bd 0 -relief raised
    add_separator $pwin 1
    add_plot $pwin depth_track 150  1 1 -bd 0 -relief raised
    add_plot $pwin seq_ruler    50  0 0 -bd 0 -relief raised
}


##############################################################################
# Test version when running as a script in its own right
if {[string match "*depth.tcl" $argv0]} {
    source $env(STADTABL)/shlib.conf
    load $env(STADLIB)/${lib_prefix}tk_utils${lib_suffix}
    load_package tk_utils
    tk_utils_init
    load_package gap5
    load libtgap.so g5

    #set dbname "/lustre/seq/scratch1/jkb/pe_rmdup.0"
    #set dbname "/lustre/seq/scratch1/jkb/1112_s_1"
    foreach dbname $argv {break;}
    set io [g5::open_database -name $dbname -access rw]

    wm withdraw .
    bind . <1> {puts "[llength [info command]] commands"}

    #add_plot $w bg_grid
    set pwin .plot
    1.5Dplot $pwin $io 800 600
    add_plot $pwin seq_depth  50 -bd 2 -relief raised
    add_plot $pwin seq_seqs -200 -bd 2 -relief raised
    add_plot $pwin seq_ruler  50 -bd 2 -relief raised
    #add_plot $pwin bg_grid
}




