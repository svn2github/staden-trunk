#
# This is the Tcl/Tk interface to the gap4 haplotype splitting module
#

# A poor mans class because I can't be bothered with the more heavyweight
# incrTcl solution.
#
# This basically just uses a single data matrix as instance-specific storage
# for the ::haplo functions. The "new" method takes "-arg val" parameter list
# and automatically sets them so that data(-arg)==val.

# The core data structure is passed by reference to the methods ($d typically,
# accessed via $data). It contains the following elements:
#
# -io 			GapIO handle
# -contig		Contig identifier
#  cnum			Contig number (computed from -contig)
# -lreg			Leftmost base (0 for all)
# -rreg			Rightmost base (0 for all)
# -fastmode		Fast mode param and also checkbox variable

# -discrep_cutoff 	Discrepancy cutoff
# -min_base_qual	Minimum base confidence when computing discrepancies
# -snp_cutoff		Minimum SNP score after adjustments (pad, polyX, etc)?
# -two_alleles		Boolean: whether to adjust snp scores for diploid orgs.

# -minscore		Minimum score acceptable during cluster merges
# -twopass		Boolean: whether to add zero-cost edges and recompute
#			during clustering algorithm

# snps			Tcl list of the form:
#			{SNP_pos SNP_score template template ...} ...
#			where template is {t_number t_score base confidence}
# snp_status		Tcl list with one item per SNP indicating SNP status
#			(unknown, low quality, high quality, etc)

# toplevel		Tk widget pathname for the toplevel window
# canvas		Tk widget pathname for the main canvas
# disp_bases		Boolean: whether to display bases
# disp_depth		Boolean: whether to display template depths
# disp_position		Boolean: whether to display SNP positions
# disp_ruler		Boolean: whether to display vertical ruler
# disp_score		Boolean: whether to display SNP score
# disp_selected		Boolean: whether to display SNP checkboxes
# disp_sets		Boolean: whether to display clustered sets
# disp_sets_check	Tk pathname of set checkboxes

# ruler_scale		Scale from bases to pixels for the Y ruler
# ruleritem,<X>		Y Coord for ruler items for SNP <X> (0 onwards)
# rulerx		Position of rightmost X coordinate for ruler plot
# selected,<X>		Boolean: whether SNP <X> is selected
# regid			io-reg ID for this plot
# cursor_id		ID for the cursor used by this plot
# cursoritem_<X>	Canvas item id for the cursor with ID <X>
# cursorapos_<X>	Absolute position of cursor <X>
# status_line		Textvariable for the status line at the window bottom

namespace eval haplo {
    variable counter 0

    proc new {args} {
	variable counter

	set data [namespace current]::instance$counter
	incr counter

	set ${data}(-lreg)           0
	set ${data}(-rreg)           0
	set ${data}(-discrep_cutoff) 40
	set ${data}(-min_base_qual)  15
	set ${data}(-two_alleles)    1
	set ${data}(-snp_cutoff)     10
	set ${data}(-minscore)       0
	set ${data}(-twopass)        0
	set ${data}(-fastmode)       0

	foreach {a b} $args {
	    set ${data}($a) $b
	}
	set ${data}(cnum) [db_info get_contig_num \
			       [set ${data}(-io)] \
			       [set ${data}(-contig)]]
	
	return $data
    }

    proc delete {d} {
	unset $d
    }
}

# ----------------------------------------------------------------------------
# Haplo methods

# An interface to the "haplo snps" method in C
proc haplo::snps {d {snp_cutoff {}} {discrep_cutoff {}}} {
    upvar $d data

    if {$snp_cutoff != ""} {
	set data(-snp_cutoff) $snp_cutoff
    }
    if {$discrep_cutoff != ""} {
	set data(-discrep_cutoff) $discrep_cutoff
    }

    set data(snps) [haplo snps \
			-io $data(-io) \
			-contigs "{$data(-contig) $data(-lreg) $data(-rreg)}"\
			-discrep_cutoff $data(-discrep_cutoff) \
			-snp_cutoff $data(-snp_cutoff) \
			-min_base_qual $data(-min_base_qual) \
			-two_alleles $data(-two_alleles)]
    set data(snps) [lindex $data(snps) 0]
    set status {}
    set num 0
    foreach snp $data(snps) {
	puts "SNP $num = $snp"
	lappend status unknown
	set data(selected,$num) 1
	incr num
    }
    set data(snp_status) $status

    return $data(snps)
}

# Creates and returns a new toplevel window
proc haplo::toplevel {} {
    global haplo_defs
    variable counter

    set t [keylget haplo_defs CANDIDATE_SNPS.WIN]$counter
    incr counter
    xtoplevel $t

    return $t
}

#-----------------------------------------------------------------------------
# Produces the vertical candidate SNP location window, but doesn't actually
# draw anything
proc haplo::create_display {d} {
    upvar $d data
    variable boldfont
    variable textfont
    variable xcoord
    variable ycoord

    set f [toplevel]
    set contig $data(cnum)

    set data(toplevel) $f

    grid columnconfigure $f 0 -weight 1
    grid rowconfigure $f 3 -weight 1

    # Candidate SNPs parameters
    set w [frame $f.control1 -bd 2 -relief raised]
    grid $w -sticky nsew -columnspan 2

    xentry $w.minsnp \
	-label "Min. SNP score" \
	-textvariable ${d}(-snp_cutoff) \
	-width 4 \
	-type "int 0"

    xentry $w.mindisc \
	-label "Min. discrepancy score" \
	-textvariable ${d}(-discrep_cutoff) \
	-width 4 \
	-type "int 0"

    xentry $w.minqual \
	-label "Min. base qual." \
	-textvariable ${d}(-min_base_qual) \
	-width 4 \
	-type "int 0"

    checkbutton $w.two_alleles \
	-text "2 alleles only" \
	-variable ${d}(-two_alleles)

    button $w.recalc \
	-text "Recalculate candidate SNPs" \
	-command "::haplo::snps $d; ::haplo::redisplay $d 0"
    pack $w.minsnp $w.mindisc $w.minqual $w.two_alleles -side left -padx 10
    pack $w.recalc -side right

    # Clustering parameters
    set w [frame $f.control2 -bd 2 -relief raised]
    grid $w -sticky nsew -columnspan 2

    xentry $w.minscore \
	-label "Min. merge score" \
	-textvariable ${d}(-minscore) \
	-width 4 \
	-type "int 0"

    checkbutton $w.twopass \
	-text "2nd inference pass" \
	-variable ${d}(-twopass)

    checkbutton $w.fast \
	-text "Fast mode" \
	-variable ${d}(-fastmode)

    button $w.cluster \
	-text "Cluster by SNPs" \
	-command [list [namespace current]::ClusterSNP $d]
    pack $w.minscore $w.twopass $w.fast -side left -padx 10
    pack $w.cluster -side right

    # Option buttons for data to display
    set w [frame $f.options -bd 2 -relief raised]
    grid $w -sticky nsew -columnspan 2

    checkbutton $w.depth    -variable ${d}(disp_depth)    -text Depth \
	-command "haplo::redisplay $d"
    checkbutton $w.ruler    -variable ${d}(disp_ruler)    -text Ruler \
	-command "haplo::redisplay $d"
    checkbutton $w.selected -variable ${d}(disp_selected) -text Use \
	-command "haplo::redisplay $d"
    checkbutton $w.position -variable ${d}(disp_position) -text Position \
	-command "haplo::redisplay $d"
    checkbutton $w.score    -variable ${d}(disp_score)    -text Score \
	-command "haplo::redisplay $d"
    checkbutton $w.bases    -variable ${d}(disp_bases)    -text Bases \
	-command "haplo::redisplay $d"
    checkbutton $w.sets     -variable ${d}(disp_sets)     -text Sets \
	-command "haplo::redisplay $d" -state disabled

    set data(disp_depth)      0
    set data(disp_ruler)      1
    set data(disp_selected)   1
    set data(disp_position)   1
    set data(disp_score)      1
    set data(disp_bases)      1
    set data(disp_sets)       0
    set data(disp_sets_check) $w.sets

    grid $w.depth $w.ruler $w.selected $w.position $w.score $w.bases $w.sets

    # Canvas
    set w [canvas $f.canvas \
	       -width 800 \
	       -height 500 \
	       -xscrollcommand "$f.xscroll set" \
	       -yscrollcommand "$f.yscroll set"]
    scrollbar $f.xscroll -command "$w xview" -orient horizontal
    scrollbar $f.yscroll -command "$w yview" -orient vertical
    set data(canvas) $w

    grid $f.canvas $f.yscroll -sticky nsew
    grid $f.xscroll -sticky nsew

    # Fonts
    set textfont {helvetica 8}
    set boldfont {helvetica 8 bold}

    # Coordinates/widths, xcoord(ruler) is an absolute X coordinate.
    # position and onwards are all relative widths for that display type.
    # Ycoords start with 0 being the first base in the contig and negative
    # values therefore being above the ruler graphics.
    set ycoord(heading)  -50
    set ycoord(setchk)   -11
    set xcoord(plot)      20
    set xcoord(depth)     80
    set xcoord(ruler)     80
    set xcoord(position)  80
    set xcoord(score)     70
    set xcoord(bases)    140
    set xcoord(bases+)    25
    set xcoord(selected)  50
    set xcoord(sets)     200
    set xcoord(sets+)     20
    set xcoord(width)    430

    # scroll-wheel
    bind $w <MouseWheel> {puts "%W yview scroll %D units";%W yview scroll %D units}

    # status display
    label $f.status -anchor w -textvariable ${d}(status_line)
    grid $f.status -sticky nsew -columnspan 2

    # Control buttons.
    frame $f.select
    grid $f.select -sticky nsew -columnspan 2
    label $f.select.l -text "Select:"
    button $f.select.all  -text All            \
	-command "haplo::select_snps $d all"
    button $f.select.none -text None           \
	-command "haplo::select_snps $d none"
    button $f.select.lq   -text "Low quality"  \
	-command "haplo::select_snps $d lq"
    button $f.select.hq   -text "High quality" \
	-command "haplo::select_snps $d hq"
    button $f.filter -text "Remove unselected" \
	-command "haplo::remove_snps $d"
    grid $f.select.l $f.select.all $f.select.none $f.select.lq $f.select.hq \
	$f.filter -sticky ew -padx 10

    frame $f.buttons
    button $f.buttons.move \
	-text "Make contigs from selected sets" \
	-command [list [namespace current]::MoveContigs $d]
    button $f.buttons.cancel \
	-text Cancel \
	-command "haplo::shutdown $d"

    grid $f.buttons -sticky nsew -columnspan 2
    grid $f.buttons.move $f.buttons.cancel -sticky ew -padx 10

    wm protocol $data(toplevel) WM_DELETE_WINDOW "haplo::shutdown $d"

    # Register with "io-reg" so that we can be notified of changes to this
    # contig.
    set data(regid) [contig_register \
			 -io $data(-io) \
			 -contig $contig \
			 -command "haplo::reg_callback $d" \
			 -flags {REQUIRED LENGTH JOIN_TO DELETE COMPLEMENT \
				 CURSOR_NOTIFY}]

    # Add a shared cursor
    set data(cursor_id) [create_cursor \
			     -io $data(-io) \
			     -cnum $contig \
			     -sent_by $data(regid)]

    return $f
}

proc haplo::reg_callback {d id type args} {
    if {[catch {

    upvar $d data
    switch $id {
	"QUERY_NAME" {
	    return "Candidate SNPs plot"
	}

	"QUERY_PARAMS" {
	    return ""
	}

	"DELETE" -
	"QUIT" {
	    shutdown $d
	}

	"LENGTH" -
	"JOIN_TO" -
	"COMPLEMENT" {
	    snps $d
	    set data(disp_sets) 0
	    unset data(rsets)
	    $data(disp_sets_check) configure -state disabled
	    redisplay $d 0
	}

	"CURSOR_NOTIFY" {
	    set pos  [keylget args pos]
	    set seq  [keylget args seq]
	    set apos [keylget args abspos]
	    set job  [keylget args job]
	    set cid  [keylget args id]
	    set refs [keylget args refs]
	    set sent_by [keylget args sent_by]
	    set cpos $apos

	    if {!$data(disp_ruler)} {
		break;
	    }

	    # DECREMENT is notified before refs is decremented, so we work
	    # out here the ref count after the event will have been processed.
	    if {[lsearch $job DECREMENT] != -1} {incr refs -1}

	    # It's possible we may be drawing cursors before the
	    # display has been created. In this case ruler_scale is
	    # undefined, but our cursor will be auto-scaled once the
	    # ruler is drawn.
	    if {[info exists data(ruler_scale)]} {
		set cpos [expr {$apos * $data(ruler_scale)}]
	    }

	    # We don't show the cursor when the reference count is 0 or
	    # the reference count is 1 and it is "our" cursor.
	    if {![info exists data(regid)] || \
		    ($data(regid) != $sent_by && $refs >= 1) || \
		    $refs >= 2} {
		set visible 1
	    } else {
		set visible 0
	    }

	    set data(cursorapos_$cid) $apos

	    # Create, display, move, or destroy the cursor line as
	    # appropriate.
	    if {[info exists data(cursoritem_$cid)]} {
		if {!$visible} {
		    $data(canvas) delete $data(cursoritem_$cid)
		    unset data(cursoritem_$cid)
		} else {
		    if {[lsearch $job MOVE] != -1} {
			$data(canvas) coords $data(cursoritem_$cid) \
			    0 $cpos $data(rulerx) $cpos
			$data(canvas) raise $data(cursoritem_$cid)
		    }
		}
	    } else {
		if {$visible} {
		    set data(cursoritem_$cid) \
		    [$data(canvas) create line 0 $cpos $data(rulerx) $cpos \
			 -width 2 -fill blue -tags cursor]
		    $data(canvas) bind $data(cursoritem_$cid) <<move-drag>> \
			"haplo::cursor_move $d $cid %y"
		}
	    }
	}
    }

    } err]} {
	puts err=$err
    }

    return ""
}

# Notifies the request to move a cursor. This does not actually update
# the displayed cursor. Instead it uses contig_notify to notify all users of
# the cursor (including us) which triggers a redraw.
proc haplo::cursor_move {d cid y} {
    upvar $d data

    set cpos [$data(canvas) canvasy $y]
    set apos [expr {int($cpos / $data(ruler_scale))}]
    keylset l id $cid seq -1 pos -1 abspos $apos sent_by $data(regid) job MOVE
    contig_notify \
	-io $data(-io) \
	-cnum $data(cnum) \
	-type CURSOR_NOTIFY \
	-args $l
}

# Destroys the window and tidies up any contig and cursor registrations
proc haplo::shutdown {d} {
    upvar $d data

    set w $data(toplevel)
    bind $w <Destroy> {}
    catch {contig_deregister -io $data(-io) -id $data(regid)}
    delete_cursor -io $data(-io) -cnum $data(cnum) -id $data(cursor_id)
    destroy $w
    unset data
}

proc haplo::select_snps {d filter} {
    upvar $d data

    set cutoff [lsearch {notfound unknown lq hq} $filter]
    if {$filter == "none"} {
	set cutoff 100
    }
    set num 0
    foreach stat $data(snp_status) {
	set status [lsearch {notfound unknown lq hq} $stat]
	set data(selected,$num) [expr {($status >= $cutoff) ? 1 : 0}]
	incr num
    }
}

# Takes SNP scores from the entry boxes and stores them in the data(snps)
# array
proc haplo::snps_from_display {d} {
    upvar $d data
    set w $data(canvas)
    
    set snps $data(snps)
    set new_snps {}
    set num 0
    foreach snp $data(snps) {
	if {[winfo exists $w.score_$num]} {
	    set snp [lreplace $snp 1 1 [$w.score_$num get]]
	    incr num
	}
	lappend new_snps $snp
    }
    set data(snps) $new_snps
}

proc haplo::redisplay {d {save_snp_scores 1}} {
    upvar $d data
    set w $data(canvas)

    if {$save_snp_scores} {
	snps_from_display $d
    }
    $w delete all

    display $d
}

proc haplo::display {d} {
    upvar $d data
    variable xcoord
    variable ycoord
    variable textfont

    set f $data(toplevel)
    set w $data(canvas)
    set contig $data(cnum)

    set io $data(-io)
    set snps $data(snps)

    set c [io_read_contig $io $contig]
    set clen [keylget c length]

    set nsnps [llength $snps]

    set headingfont {helvetica 12 bold}

    # Headings
    set xpos $xcoord(plot)
    foreach heading {Depth Ruler Selected Position Score Bases Sets} {
	set type [string tolower $heading]
	if {$data(disp_$type)} {
	    if {$heading == "Selected"} {
		set h Use
	    } else {
		set h $heading
	    }
	    $w create text $xpos $ycoord(heading) \
		-anchor nw -text $h -font $headingfont
	    incr xpos $xcoord($type)
	}
    }
   
    
    # Data types
    set xpos $xcoord(plot)
    foreach type {depth ruler selected position score bases sets} {
	if {$type == "ruler"} {set xrul1 $xpos}
	if {$data(disp_$type)} {
	    incr xpos [display_$type $d $xpos]
	}
	if {$type == "ruler"} {set xrul2 $xpos}
    }
    display_cursors $d


    # Heading underline. Do now that we know the plot width
    set x2 $xpos
    array set fnt [font metrics $headingfont]
    set ypos [expr {$ycoord(heading) + $fnt(-linespace)+1}]
    $w create line 0 $ypos $x2 $ypos -width 2

    # boundary box for text
    set num 0
    array set fnt [font metrics $textfont]
    foreach snp $snps {
 	set ypos1 [expr {$data(ypos,$num) - $fnt(-linespace)*0.75}]
 	set ypos2 [expr {$data(ypos,$num) + $fnt(-linespace)*0.75}]
	set status [lindex $data(snp_status) $num]
	switch $status {
	    "hq"       {set bg "#eeeeee"}
	    "lq"       {set bg "#cccccc"}
	    "notfound" {set bg "#aaaaaa"}
	    default    {set bg "#d9d9d9"}
	}
	set snpbg($num) $bg
 	set b [$w create rectangle $xrul2 $ypos1 $x2 $ypos2 \
		   -fill $bg -width 0 -tags "snp_$num snpbg_$num textcoord"]
	$w lower $b

 	incr num
    }


    # Rescale the ruler to match the scale dictated by the text elements
    set bbox [$w bbox textcoord]
    if {$bbox != ""} {
	set texth [expr {[lindex $bbox 3]-[lindex $bbox 1]}]
	set textw [expr {[lindex $bbox 2]-[lindex $bbox 0]}]
    } else {
	set texth 400
	set textw 0
    }
    set bbox [$w bbox ruler]
    if {$bbox != ""} {
	set rulerh [expr {[lindex $bbox 3]-[lindex $bbox 1]}]
	set rulery 0
	set scale [expr {double($texth)/$rulerh}]
	$w scale ruler 0 $rulery 1 $scale
	$w scale cursor 0 $rulery 1 $scale
	set data(ruler_scale) $scale

	# Draw the connecting lines from ruler to text plots
	if {$data(disp_ruler) && $textw > 20} {
	    set x1 [expr {$xrul1+10}]
	    set x3 $xrul2
	    set x2 [expr {int($x3)-10}]
	    for {set i 0} {$i < $nsnps} {incr i} {
		set y1 [lindex [$w coords $data(ruleritem,$i)] 1]
		set y3 $data(ypos,$i)
		set y2 $y3
		$w create line $x1 $y1 $x2 $y2 $x3 $y3 -tags "snp_$i snpx_$i"
	    }
	}
    }
    
    $w configure -scrollregion [$w bbox all]

    for {set i 0} {$i < $nsnps} {incr i} {
	$w bind snp_$i <Enter> "$w itemconfigure snpx_$i -width 3 -fill blue"
	$w bind snp_$i <Enter> "+$w itemconfigure snpbg_$i -fill \#b0f0ff"
	$w bind snp_$i <Leave> "$w itemconfigure snpx_$i -width 1 -fill black"
	$w bind snp_$i <Leave> "+$w itemconfigure snpbg_$i -fill $snpbg($i)"
	$w bind snp_$i <<use>> [list haplo::snp_invoke_editor $d $i]
    }
}

# -----------------------------------------------------------------------------
# Display components for each vertical track in our plot
proc haplo::display_depth {d xpos} {
    upvar $d data
    variable xcoord

    set depth [haplo tdepth -io $data(-io) -contig $data(-contig)]
    set max [lindex $depth 0]; puts depth=$max
    set max [expr {int(($max+9)/10) * 10}]
    set xscale [expr {($xcoord(depth)-20) / double($max)}]
    set w $data(canvas)
    set contig $data(cnum)
    set c [io_read_contig $data(-io) $contig]
    set clen [keylget c length]

    # Grey dividers in incremements of 10
    set y1 0
    set y2 $clen
    for {set d 0} {$d <= $max} {incr d 10} {
	set x [expr {$xpos+$xcoord(depth)-20-$d*$xscale}]
	$w create line $x $y1 $x $y2 -tags ruler -fill grey66
    }

    # Text for the axis
    set x [expr {$xpos+$xcoord(depth)-20-0*$xscale}]
    set y2 [expr {$y1-10}]
    $w create line $x $y1 $x $y2 -fill grey66
    $w create text $x $y2 -text 0 -anchor s
    set x [expr {$xpos+$xcoord(depth)-20-$max*$xscale}]
    set y2 [expr {$y1-10}]
    $w create line $x $y1 $x $y2 -fill grey66
    $w create text $x $y2 -text $max -anchor s

    # The plot itself
    set coords [list [expr {$xpos+$xcoord(depth)-20}] 0]
    foreach {st en d} [lrange $depth 1 end] {
	incr st 0
	incr en 0
	set x [expr {$xpos+$xcoord(depth)-20-$d*$xscale}]

	lappend coords $x $en
	#$w create line $x $st $x $en -tags ruler
    }
    $w create line $coords -tags ruler

    return $xcoord(depth)
}

proc haplo::display_ruler {d xpos} {
    upvar $d data
    variable xcoord
    variable textfont

    set snps $data(snps)
    set w $data(canvas)
    set contig $data(cnum)
    set c [io_read_contig $data(-io) $contig]
    set clen [keylget c length]

    # consensus ruler with SNP line
    $w create line $xpos 0 $xpos $clen -tags ruler
    set data(rulerx) [expr {$xpos+10}]

    set num 0
    array set fnt [font metrics $textfont]
    foreach snp $snps {
	set pos [lindex $snp 0]
	set data(ruleritem,$num) [$w create line \
				      [expr {$xpos-10}] $pos \
				      [expr {$xpos+10}] $pos \
	                              -tags "ruler snp_$num snpx_$num"]
	set data(ypos,$num) [expr {$fnt(-linespace)*1.5*$num}]

	# # fake text to keep scales working
	# $w create text 0 $data(ypos,$num) -tags textcoord

	incr num
    }

    return $xcoord(ruler)
}

proc haplo::display_position {d xpos} {
    upvar $d data
    variable textfont
    variable xcoord

    set snps $data(snps)
    set w $data(canvas)

    set num 0
    foreach snp $snps {
	$w create text $xpos $data(ypos,$num) \
	    -text "[lindex $snp 0]" \
	    -font $textfont \
	    -anchor w \
	    -tags "textcoord snp_$num"

	incr num
    }

    return $xcoord(position)
}

proc haplo::display_score {d xpos} {
    upvar $d data
    variable textfont
    variable xcoord

    set snps $data(snps)
    set w $data(canvas)

    set num 0
    foreach snp $snps {
# 	set i [$w create text $xpos $data(ypos,$num) \
# 		   -text [format "%5.2f" [lindex $snp 1]] \
# 		   -font $textfont \
# 		   -anchor w \
# 		   -tags "textcoord snp_$num"]
	if {[winfo exists $w.score_$num]} {
	    destroy $w.score_$num
	}
        entry $w.score_$num -width 5
	$w.score_$num insert end [format "%5.2f" [lindex $snp 1]]
	$w create window $xpos $data(ypos,$num) \
	    -window $w.score_$num \
	    -anchor w \
	    -tags "textcoord snp_$num"

	incr num
    }

    return $xcoord(score)
}

proc haplo::display_bases {d xpos} {
    upvar $d data
    variable textfont
    variable boldfont
    variable xcoord

    set snps $data(snps)
    set w $data(canvas)

    set num 0
    foreach snp $snps {
 	array set basea [list A 0 C 0 G 0 T 0 - 0 * 0]
 	foreach b [lrange $snp 2 end] {
 	    incr basea([lindex $b 2])
 	}
	set xpos2 $xpos
 	foreach b {A C G T *} {
 	    if {$basea($b) > 99} {
 		set nums X
 	    } else {
 		set nums $basea($b)
 	    }
 	    if {$basea($b) > 1} {
 		set fn $boldfont
	    } else {
 		set fn $textfont 
 	    }
	    $w create text $xpos2 $data(ypos,$num) \
 		-text "$nums$b" \
 		-font $fn \
 		-anchor w \
 		-tags "textcoord snp_$num"

	    incr xpos2 $xcoord(bases+)
 	}

	incr num
    }

    return $xcoord(bases)
}

proc haplo::display_selected {d xpos} {
    upvar $d data
    variable textfont
    variable xcoord

    set snps $data(snps)
    set w $data(canvas)

    set num 0
    foreach snp $snps {
	catch {checkbutton $w.check_$num -padx 0 -variable ${d}(selected,$num)}
	$w create window $xpos $data(ypos,$num) -window $w.check_$num \
	    -anchor w
	incr num
    }

    return $xcoord(selected)
}

proc haplo::display_sets {d xpos} {
    upvar $d data
    variable textfont
    variable boldfont
    variable xcoord

    puts [info level [info level]]

    set io $data(-io)
    set snps $data(snps)
    set contig $data(cnum)
    set sets $data(cons)
    set w $data(canvas)

    puts "Drawing [llength $sets] sets"

    # Create sets
    set snum 0
    array set fnt [font metrics $textfont]
    foreach snp $snps {
	set ypos $data(ypos,$snum)
	set setnum 0
	foreach g $sets {
	    set cons [lindex $g 1]
	    set qual [lindex $g 2]
	    set ind [expr {[lindex $snp 0]-1}]
	    set base [string index $cons $ind]
	    binary scan [string index $qual $ind] c q
	    set x [expr {$xpos + $setnum*$xcoord(sets+)}]
	    set y1 [expr {$ypos-$fnt(-linespace)*0.75}]
	    set y2 [expr {$ypos+$fnt(-linespace)*0.75}]
	    if {$base != "-"} {
		set x1 [expr {$x-3}]
		set x2 [expr {$x+10}]
		if {$q > 100} {set q 100}
		set q [expr {int(105+$q*1.5)}]
		set col [format "#%02x%02x%02x" $q $q $q]
		set item [$w create rectangle $x1 $y1 $x2 $y2 \
			      -fill $col -tags snp_$snum]
		$w bind $item <Enter> [list haplo::update_status $d $snum $setnum]
		set item [$w create text $x $ypos \
			      -text $base -anchor w \
			      -tags "snp_$snum set_$setnum"]
	    } else {
		set x1 [expr {$x+3}]
		$w create line $x1 $y1 $x1 $y2 -fill grey66
	    }
	    incr setnum
	}
	
	incr snum
    }

    # Set menus
    set setnum 0
    foreach g $sets {
	$w bind set_$setnum <<menu>> "haplo::set_menu $d $setnum %X %Y"
	$w bind set_$setnum <Enter> \
	    "set ${d}(status_line) \"Set $setnum: [llength [lindex $g 0]] templates\""
	incr setnum
    }

    # Set checkbuttons
    set setnum 0
    set y $ycoord(setchk)
    foreach g $sets {
	set x [expr {$xpos + $setnum*$xcoord(sets+)-4}]
	catch {checkbutton $w.set_$setnum -padx 0 -pady 0 \
		   -variable ${d}(set_selected,$setnum)}
	$w create window $x $y -window $w.set_$setnum -anchor sw
	bind $w.set_$setnum <<menu>> "haplo::set_menu $d $setnum %X %Y"

	incr setnum
    }

    return [expr {$setnum * $xcoord(sets+)}]
}

# Redraw cursors
proc haplo::display_cursors {d} {
    upvar $d data

    $data(canvas) delete cursor
    foreach x [array names data -glob cursoritem_*] {
	regexp {.*_(.*)} $x _ cid
	set cpos $data(cursorapos_$cid)
	set data(cursoritem_$cid) \
	    [$data(canvas) create line 0 $cpos $data(rulerx) $cpos \
		-tags cursor -width 2 -fill blue]
	$data(canvas) bind $data(cursoritem_$cid) <<move-drag>> \
	    "haplo::cursor_move $d $cid %y"
    }
}

# -----------------------------------------------------------------------------
# Display callbacks

# Controls the contig editor
proc haplo::snp_invoke_editor {d snum} {
    upvar $d data

    set snp [lindex $data(snps) $snum]

    if {[info exists data(rsets)]} {
	edit_contig \
	    -io $data(-io) \
	    -contig $data(-contig) \
	    -pos [lindex $snp 0] \
	    -reuse 1 \
	    -nojoin 1 \
	    -sets $data(rsets)
    } else {
	edit_contig \
	    -io $data(-io) \
	    -contig $data(-contig) \
	    -pos [lindex $snp 0] \
	    -reuse 1 \
	    -nojoin 1
    }
}


# Updates the status line with snp info
proc haplo::update_status {d snum setnum} {
    upvar $d data
    set f $data(toplevel).status
    set snp [lindex $data(snps) $snum]
    set set [lindex $data(cons) $setnum]
    
    set pos [expr {[lindex $snp 0]-1}]
    set base [string index [lindex $set 1] $pos]
    binary scan [string index [lindex $set 2] $pos] c qual

    $f configure -text "Pos:[expr {$pos+1}] Cons:$base Qual:$qual"
}


# Popup menu for sets
proc haplo::set_menu {d setnum x y} {
    upvar $d data
    set w $data(canvas)
    if {![winfo exists $w.menu]} {
	menu $w.menu -tearoff 0
    }
    $w.menu delete 0 end
    $w.menu add command \
	-label "Delete set" \
	-command "haplo::delete_set $d $setnum"
    $w.menu add command \
	-label "Merge selected sets" \
	-command "haplo::merge_sets $d"

    tk_popup $w.menu $x $y
}

# Removes a set
proc haplo::delete_set {d setnum} {
    upvar $d data
    set data(sets) [lreplace $data(sets) $setnum $setnum]
    set len [llength $data(sets)]
    for {} {$setnum < $len} {incr setnum} {
	set setnum2 [expr {$setnum+1}]
	set data(set_selected,$setnum) $data(set_selected,$setnum2)
    }
    compute_rsets $d

    compute_set_cons $d
    set_report $d
    redisplay $d
}

# Merges selected sets
proc haplo::merge_sets {d} {
    upvar $d data

    set merged {}
    set sets {}
    set setnum 0
    foreach g $data(sets) {
	if {$data(set_selected,$setnum)} {
	    foreach i $g {
		lappend merged $i
	    }
	} else {
	    lappend sets $g
	}
	incr setnum
    }
    lappend sets $merged

    set data(sets) $sets
    compute_rsets $d

    compute_set_cons $d
    set_report $d
    redisplay $d
}

# Creates sets of readings based on the existing sets of templates
# No check is done regarding whether either set happens to span contigs.
proc haplo::compute_rsets {d} {
    upvar $d data
    set rsets {}
    foreach g $data(sets) {
	lappend rsets [templates2readings $data(-io) $data(-contig) $g]
    }
    set data(rsets) $rsets
}

# Compute set consensus sequences
proc haplo::compute_set_cons {d} {
    upvar $d data

    set snps $data(snps)
    set contig $data(cnum)
    set gdata {}
    set tot [llength $data(sets)]
    set n 0
    foreach g $data(sets) {
	set data(status_line) "Computing consensus $n / $tot"
	update idletasks
	foreach {c q} [haplo consensus \
			   -io $data(-io) \
			   -templates $g \
			   -contig =$contig] {break;}
	lappend gdata [list $g $c $q]
	set data(cons,$n) $c
	set data(qual,$n) $q
	incr n
	foreach snp $snps {
	    set pos [lindex $snp 0]
	    append snp_data($pos) [string index $c [expr {$pos-1}]]
	}
    }
    set data(cons) $gdata
}

#-----------------------------------------------------------------------------
# Clusters a set of SNPs
proc haplo::sort_sets {sets} {
    set setnum 0
    set l {}
    foreach set $sets {
	lappend l [llength $set].$setnum
	set gdata($setnum) $set
	incr setnum
    }
    set gorder {}
    foreach g [lsort -real -decreasing $l] {
	lappend gorder [lindex [split $g .] 1]
    }

    set sorted {}
    foreach g $gorder {
	lappend sorted $gdata($g)
    }

    return $sorted
}

proc haplo::dump_sets {sets} {
    vfuncheader "Haplotypes"

    set setnum 0
    foreach set $sets {
	vmessage "=== Set $setnum, size [llength $set] ==="
	incr setnum
	foreach template $set {
	    vmessage "    $template"
	}
    }
}

# Returns a list of templates covering a given contig position along with the
# base calls and confidence values.
proc haplo::templates_at_pos {io contig pos} {
    set cnum [db_info get_contig_num $io $contig]

    # Get a list of base calls per template (maybe multiples reads on one
    # template) at this position
    set c [io_read_contig $io $cnum]
    for {set rn [keylget c left]} {$rn} {set rn [keylget r right]} {
	set r [io_read_reading $io $rn]
	set p [keylget r position]
	if {$p > $pos} {
	    break
	}
	if {$p+[keylget r sequence_length]-1 < $pos} {
	    continue
	}
	set i [expr {$pos-($p-[keylget r start])}]
	set base [string index [io_read_text $io [keylget r sequence]] $i]
	set conf [io_read_data $io [keylget r confidence] [keylget r length] 1]
	set conf [string index $conf $i]
	binary scan $conf c qval
	lappend tinfo([keylget r template]) $base $qval
    }

    # Compute the consensus and quality for the template base
    set tlist ""
    foreach t [array names tinfo] {
	if {[llength $tinfo($t)] == 2} {
	    foreach {bcall callq} $tinfo($t) break
	} else {
	    array set max {A 0 C 0 G 0 T 0 * 0 - 0}
	    foreach {base qual} $tinfo($t) {
		switch $base {
		    "A" - "a" { set b A }
		    "C" - "c" { set b C }
		    "G" - "g" { set b G }
		    "T" - "t" { set b T }
		    "*"       { set b * }
		    default   { set b - }
		}
		if {$max($b) < $qual} {
		    set max($b) $qual
		}
	    }
	    set bcall ""
	    set maxq 0
	    set totq 0
	    foreach {base} {A C G T * -} {
		if {$max($base) > $maxq} {
		    set maxq $max($base)
		    set bcall $base
		}
		incr totq $max($base)
	    }
	    set callq [expr {$maxq-($totq-$maxq)}]
	    if {$callq <= 0} {
		set bcall -
		set callq 0
	    }
	}

	# Build up our template info list
	lappend tlist [list $t 1 $bcall $callq]
    }

    return $tlist
}

proc haplo::set_report {d} {
    upvar $d data
    set qval 90

    set io $data(-io)
    set contig $data(-contig)
    set snps $data(snps)
    set sets $data(sets)

    # Convert snps list into an array indexed by position
    array set snpa {}
    set index 0
    foreach s $snps {
	set status [lindex $data(snp_status) $index]
	set selected $data(selected,$index)
	set snpa([lindex $s 0]) "$status $selected $s"
	incr index
    }

    # Identify high good quality differences
    set nsets [llength $sets]
    set clen [string length $data(cons,0)]
    set count 0
    array set confirmed {}
    for {set pos 0} {$pos < $clen} {incr pos} {
	catch {unset vector}
	array set vector {A 0 C 0 G 0 T 0 - 0 * 0}
	for {set i 0} {$i < $nsets} {incr i} {
	    incr vector([string index $data(cons,$i) $pos])
	}
	if {([expr {$vector(A)+$vector(-) != $nsets}]) &&
	    ([expr {$vector(C)+$vector(-) != $nsets}]) &&
	    ([expr {$vector(G)+$vector(-) != $nsets}]) &&
	    ([expr {$vector(T)+$vector(-) != $nsets}]) &&
	    ([expr {$vector(*)+$vector(-) != $nsets}])} {
	    # Possible conflict, so check quality scores
	    
	    set hq ""
	    for {set i 0} {$i < $nsets} {incr i} {
		binary scan [string index $data(qual,$i) $pos] c q
		if {$q >= $qval} {
		    append hq [string index $data(cons,$i) $pos]
		}
	    }
	    set p2 [expr {$pos+1}]
	    set confirmed($p2) 1
	    if {[string length $hq] > 1 && ![regexp {***:^(.)\1+$} $hq]} {
		#puts "=== [expr {$pos+1}] $hq ==="
		if {[info exists snpa($p2)]} {
		    set snpa($p2) [lreplace $snpa($p2) 0 0 hq]
		    #set snpa($p2) [lreplace $snpa($p2) 3 3 300]
		} else {
		    puts "+++ New snp at $p2: $q, $hq"
		    set tlist [templates_at_pos $io $contig $p2]
		    set snpa($p2) "hq 0 $p2 300 $tlist"
		}
		set snp($count) $pos
		incr count
	    } else {
		# low quality, but a difference still
		if {[info exists snpa($p2)]} {
		    set snpa($p2) [lreplace $snpa($p2) 0 0 lq]
		    #set snpa($p2) [lreplace $snpa($p2) 3 3 30]
		} else {
		    #puts "+++ New snp at $p2: $q, $hq"
		    #set tlist [templates_at_pos $io $contig $p2]
		    #set snpa($p2) "lq 0 $p2 100 $tlist"
		}
	    }
	}
    }
    puts "    $count HQ differences"

    # Turn snp array back into data(snps) list. Filter out existing SNPs
    # that were not confirmed as real.
    set snps ""
    set status ""
    set snum 0
    foreach p [lsort -integer [array names snpa]] {
	if {![info exists confirmed($p)]} {
	    puts "--- Unconfirmed SNP at $p"
	    lappend status notfound
	} else {
	    lappend status [lindex $snpa($p) 0]
	}
	set data(selected,$snum) [lindex $snpa($p) 1]
	lappend snps [lrange $snpa($p) 2 end]
	incr snum
    }
    set data(snps) $snps
    set data(snp_status) $status

    # Display contig plot
    puts "Contig overlap positions:"
    for {set units 1000000} {$units >= 1} {set units [expr {$units / 10}]} {
	puts -nonewline "     "
	for {set j 0} {$j < $count} {incr j} {
	    puts -nonewline [expr {(($snp($j)+1) / $units) % 10}]
	}
	puts ""
    }
    
    for {set i 0} {$i < $nsets} {incr i} {
	puts -nonewline [format "\#%3d " $i]
	for {set j 0} {$j < $count} {incr j} {
	    puts -nonewline "[string index $data(cons,$i) $snp($j)]"
	}
	puts ""
    }
}

# Filters snps to only those that are "selected" in data(selected,$num)
proc haplo::filter_snps {d} {
    upvar $d data

    set newsnps {}
    set num 0
    foreach snp $data(snps) {
	if {$data(selected,$num)} {
	    lappend newsnps $snp
	}
	incr num
    }

    return $newsnps
}

proc haplo::remove_snps {d} {
    upvar $d data

    snps_from_display $d
    
    set newsnps {}
    set newstatus {}
    set num 0
    foreach snp $data(snps) {
	if {$data(selected,$num)} {
	    lappend newsnps $snp
	    lappend newstatus [lindex $data(snp_status) $num]
	}
	incr num
    }

    set data(snps) $newsnps
    set data(snp_status) $newstatus

    set num 0
    foreach snp $data(snps) {
	set data(selected,$num) 1
	incr num
    }

    redisplay $d 0
}

proc haplo::ClusterSNP {d} {
    upvar $d data

    snps_from_display $d

    set io $data(-io)
    set snps $data(snps)
    set contig $data(cnum)

    set data(status_line) "Filtering SNPs"
    update idletasks
    set snps [filter_snps $d]

    if {[llength $snps] == 0} {
	return
    }

    SetBusy
    set data(status_line) "Producing splits"
    update idletasks
    puts "========" 
    puts "-minscore $data(-minscore)"
    puts "-twopass $data(-twopass)"
    puts "-fast $data(-fastmode)"
    set data(sets) [haplo split \
			  -io $io \
			  -snps $snps \
			  -verbosity 1 \
			  -minscore $data(-minscore) \
			  -twopass $data(-twopass) \
			  -fast $data(-fastmode)]
    ClearBusy
    set data(sets) [sort_sets $data(sets)]
    set data(status_line) "Obtaining reading lists"
    update idletasks
    compute_rsets $d
    dump_sets $data(sets)

    compute_set_cons $d

    set data(status_line) "Text report"
    update idletasks
    set_report $d
    # FIXME: need to now move data(snps) scores into the entry boxes.

    set data(disp_sets) 1
    $data(disp_sets_check) configure -state normal
    set data(status_line) "Redisplaying"
    update idletasks
    redisplay $d

    set data(status_line) "Done"
    update idletasks
}

# Converts a list of template names into a list of reading names.
# For efficiency all templates are assumed to be in a specific contig.
# Specify a zero contig and an empty templates list to free up any
# temporary memory used by this function.
set last_contig 0
proc haplo::templates2readings {io contig templates} {
    global tnames
    global last_contig

    if {$last_contig == 0} {
	# Free memory
	if {$templates == ""} {
	    catch {unset tnames}
	    return
	}

	# Initialise template to reading lookup table.
	set last_contig $contig
	set db [io_read_database $io]
	set nr [keylget db num_readings]
	for {set i 1} {$i <= $nr} {incr i} {
	    set r [io_read_reading $io $i]
	    set rname [io_read_reading_name $io $i]
	    if {[keylget r template]} {
		set t [io_read_template $io [keylget r template]]
		set tname [io_read_text $io [keylget t name]]
		lappend tnames($tname) $rname
	    }
	}
    }

    set reads {}
    foreach t $templates {
	foreach r $tnames($t) {
	    lappend reads $r
	}
    }
    return $reads
}

# Creates a fake reading and links it on to the start of contig $contig
proc haplo::add_fake {io cnum rname} {
    set c [io_read_contig $io $cnum]
    set clen [keylget c length]
    set seq [string repeat "-" $clen]

    # Create the reading structure
    set rnum [io_add_reading $io]
    set r [io_read_reading $io $rnum]
    io_write_reading_name $io $rnum $rname
    keylset r position 1
    keylset r length $clen
    keylset r sequence_length $clen
    keylset r start 0
    keylset r end [expr {$clen+1}]
    keylset r sequence [io_allocate $io text]
    io_write_text $io [keylget r sequence] $seq
    keylset r left 0
    keylset r right [keylget c left]
    io_write_reading $io $rnum $r

    # Link the previous left-most reading to our fake reading
    set r [io_read_reading $io [keylget c left]]
    keylset r left $rnum
    io_write_reading $io [keylget c left] $r
    
    # Update the old left-neighbour for the contig
    keylset c left $rnum
    io_write_contig $io $cnum $c
}

# Move selected contigs
proc haplo::MoveContigs {d} {
    upvar $d data

    set io $data(-io)
    set contig $data(cnum)
    set rsets $data(rsets)

#    if {![quit_displays $io "Split contigs"]} {
#	return
#    }
    puts "contig_lock_write -io $io -cnum $contig"
    if {[contig_lock_write -io $io -cnum $contig] != 0} {
	bell
	return
    }

    SetBusy
    # A fake sequence to ensure that the "remainder" is held together.
    set ng 0
    add_fake $io $contig fake-$ng

    foreach reads $rsets {
	if {!$data(set_selected,$ng)} {
	    incr ng
	    continue
	}
	incr ng
	add_fake $io $contig fake-$ng
	set reads [linsert $reads 0 fake-$ng]
	disassemble_readings -io $io -readings $reads -move 2 -duplicate_tags 0
    }
    templates2readings $io 0 {}
    ClearBusy

    puts "contig_notify -io $io -type LENGTH -cnum $contig -args {}"
    puts [contig_notify -io $io -type LENGTH -cnum $contig -args {}]
}


# -----------------------------------------------------------------------------
# DEBUG
proc haplo::dump_snps {snps} {
    foreach snp $snps {
	set pos [lindex $snp 0]
	set score [lindex $snp 1]
	puts -nonewline "SNP $pos\t$score\t"
	foreach t [lrange $snp 2 end] {
	    puts -nonewline [lindex $t 2]
	}
	puts ""
    }
}

# DEBUG
proc haplo::dump_sets {sets} {
    set setnum 0
    set l {}
    foreach set $sets {
	lappend l [llength $set].$setnum
	set gdata($setnum) $set
	incr setnum
    }
    set gorder {}
    foreach g [lsort -real -decreasing $l] {
	lappend gorder [lindex [split $g .] 1]
    }

    set setnum 0
    foreach g $gorder {
	set set $gdata($g)
	puts "=== Set $setnum, size [llength $set] ==="
	incr setnum
	foreach template $set {
	    puts "    $template"
	}
    }
}


#-----------------------------------------------------------------------------
# EXTERNAL:   Tcl dialogue for Candiate SNPS

# SNP candidates
proc CandidateSNPs {io} {
    global haplo_defs

    set defs [keylget haplo_defs CANDIDATE_SNPS]
    set t [keylget defs WIN]
    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "SNP Candidates"

    contig_id $t.id -io $io
    
    okcancelhelp $t.ok \
	-ok_command "CandidateSNPs2 $io $t $t.id" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {FIXME}" \
	-bd 2 -relief groove

    pack $t.id $t.ok -side top -fill x
}

set CandidateSNPs.win 0
proc CandidateSNPs2 {io t id} {
    global haplo_defs CandidateSNPs.win

    set defs [keylget haplo_defs CANDIDATE_SNPS]

    if {[set contign [contig_id_gel $id]] == ""} {bell; return}
    if {[set lreg  [contig_id_lreg $id]] == ""} {bell; return}
    if {[set rreg  [contig_id_rreg $id]] == ""} {bell; return}

    destroy $t

    SetBusy
    set h [::haplo::new \
	       -io $io \
	       -contig $contign \
	       -lreg $lreg \
	       -rreg $rreg \
	       -discrep_cutoff 40 \
	       -min_base_qual 15 \
	       -snp_cutoff 10]

    set snps [::haplo::snps $h]

    ClearBusy

    if {$snps == ""} {
	tk_messageBox \
	    -icon info \
	    -title "No SNPs found" \
	    -message "No candidate SNP sites were detected." \
	    -type ok
	return
    }

    ::haplo::create_display $h
    ::haplo::display $h
}
