#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#set CurContig "FIXME"
set AR_counter 0
set AR_active 0

proc NumReadings {io} {
    set db [io_read_database $io]
    return [keylget db num_readings]
}

proc NumContigs {io} {
    set db [io_read_database $io]
    return [keylget db num_contigs]
}

proc ContigLength {io contig} {
    set c [io_read_contig $io $contig]
    return [keylget c length]
}

#-----------------------------------------------------------------------------
# Routines for displaying a line containing particular information
#-----------------------------------------------------------------------------

#
# The simplest - a single integer
#
proc AR_draw_int {io f list item min max} {
    frame $f._$item
    if {$max > $min} {
	set type "CheckIntRange $min $max"
    } else {
	set type "CheckIntMin $min"
    }

    entrybox $f._$item.int \
	-title $item \
	-width 8 \
	-default [keylget list $item] \
	-type $type
    pack $f._$item.int -side left -fill x -expand 1

    pack $f._$item -side top -fill x
    return $f._$item
}

#
# Shows the string referenced by this record (uses io_read_text)
#
proc AR_draw_str {io f width list item args} {
    if {![winfo exists $f._$item]} {
        frame $f._$item

	if {[set x [keylget list $item]] != 0} {
	    set s [io_read_text $io $x]
	} else {
	    set s ""
	}

        entrybox $f._$item.str \
	    -title $item \
	    -width $width \
	    -default $s
         pack $f._$item.str -side left -fill x -expand 1
    
        entrybox $f._$item.int \
	    -width 8 \
	    -default [keylget list $item] \
	    -command "AR_draw_str $io $f $width {$list} $item"
        pack $f._$item.int -side left -fill x
    } else {
	if {$args != 0} {
	    set s [io_read_text $io $args]
	} else {
	    set s ""
	}

       	entrybox_configure $f._$item.str \
	    -default $s
    }

    pack $f._$item -side top -fill x
    return $f._$item
}

#
# Adds a 'step-to' button to display the structure referenced by this record
#
proc AR_draw_foll {io fr f list item name} {
    global env
    frame $f._$item
    label $f._$item.label -text $item
    pack $f._$item.label -side left -fill x

    entrybox $f._$item.int \
	-width 8 \
	-default [keylget list $item] \
	-type "CheckIntMin 0"
    button $f._$item.follow -bitmap @$env(STADTABL)/follow.bitmap \
	-command "if {[set x [entrybox_get $f._$item.int]] != \"\"} {\
	              AR_display $io $fr $name \"$x\"}"
    pack $f._$item.int $f._$item.follow -side right -fill x

    pack $f._$item -side top -fill x
    return $f._$item
}

#
# GReadings.sense
#
proc AR_draw_sense {io f list item args} {
    if {![winfo exists $f._$item]} {
        frame $f._$item
    
        label $f._$item.label -text $item
        pack $f._$item.label -side left -fill x
        
        set val [keylget list $item]
        entrybox $f._$item.int \
            -title [lindex {{original} {complemented}} $val] \
            -width 8 \
            -default $val\
            -type "CheckIntRange 0 1" \
	    -command "AR_draw_sense $io $f {$list} $item"
        pack $f._$item.int -side right -fill x
    } else {
	entrybox_configure $f._$item.int \
	    -title [lindex {{original} {complemented}} $args]
    }

    pack $f._$item -side top -fill x
    return $f._$item
}

#
# GReadings.strand
#
proc AR_draw_strand {io f list item args} {
    if {![winfo exists $f._$item]} {
        frame $f._$item
    
        label $f._$item.label -text $item
        pack $f._$item.label -side left -fill x
        
        set val [keylget list $item]
        entrybox $f._$item.int \
            -title [lindex {{forward} {reverse}} $val] \
            -width 8 \
            -default $val\
            -type "CheckIntRange 0 1" \
	    -command "AR_draw_strand $io $f {$list} $item"
        pack $f._$item.int -side right -fill x
    } else {
	entrybox_configure $f._$item.int \
	    -title [lindex {{original} {reverse}} $args]
    }

    pack $f._$item -side top -fill x
    return $f._$item
}

#
# GReadings.primer
#
proc AR_draw_primer {io f list item args} {
    if {![winfo exists $f._$item]} {
        frame $f._$item
    
        label $f._$item.label -text $item
        pack $f._$item.label -side left -fill x
        
        set val [keylget list $item]
        entrybox $f._$item.int \
            -title [lindex {{unknown} {forward} {reverse} {custom forward} {custom reverse}} $val] \
            -width 8 \
            -default $val\
            -type "CheckIntRange 0 4" \
	    -command "AR_draw_primer $io $f {$list} $item"
        pack $f._$item.int -side right -fill x
    } else {
	entrybox_configure $f._$item.int \
	    -title [lindex {{unknown} {forward} {reverse} {custom forward} {custom reverse}} $args] \
    }
    pack $f._$item -side top -fill x
    return $f._$item
}

#
# GAnnotations.type
#
proc AR_draw_type {io f list item} {
#    set itype [keylget list $item]
#    set type [format "%c%c%c%c" \
#	[expr $itype>>24] \
#	[expr ($itype>>16)&0xff] \
#	[expr ($itype>>8)&0xff] \
#	[expr $itype&0xff]]
    set type [keylget list $item]

    frame $f._$item
    entrybox $f._$item.int \
	-title $item \
	-width 8 \
	-default $type
    pack $f._$item.int -side left -fill x -expand 1

    pack $f._$item -side top -fill x
    return $f._$item
}

#
# GAnnotations.strand
#
proc AR_draw_astrand {io f list item args} {
    if {![winfo exists $f._$item]} {
        frame $f._$item
    
        label $f._$item.label -text $item
        pack $f._$item.label -side left -fill x
        
        set val [keylget list $item]
        entrybox $f._$item.int \
            -title [lindex {{forward} {reverse} {both}} $val] \
            -width 8 \
            -default $val\
            -type "CheckIntRange 0 2" \
	    -command "AR_draw_astrand $io $f {$list} $item"
        pack $f._$item.int -side right -fill x
    } else {
	entrybox_configure $f._$item.int \
	    -title [lindex {{original} {reverse}} $args]
    }

    pack $f._$item -side top -fill x
    return $f._$item
}


#
# Search for the annotation which links to this one.
# There may be more than one (if we have an error), so we don't go to that
# anno, we just list them.
#
proc AR_prev_anno {io anno_num} {
    set db [io_read_database $io]
    set nanno [keylget db Nannotations]

    set found 0
    for {set i 1} {$i <= $nanno} {incr i} {
	set anno [io_read_annotation $io $i]
	if {[keylget anno next] == $anno_num} {
	    vmessage "Linked to from annotation $i"
	    incr found
	}
    }
    if {$found == 0} {
	vmessage "Not linked to by any other annotation."
    } else {
	vmessage "Linked to by $found annotation."
    }
}

#
# Annotations referenced by:
#
proc AR_draw_ref {io fr f val} {
    global env
    set c 0

    set t [frame $f.previous_anno]
    button $t.but \
	    -text "List linking annotations" \
	    -command "AR_prev_anno $io $val"
    pack $t.but -side left
    pack $t -side top -fill x

    set r [annotation_address -io $io -annotation $val]
    foreach i $r {
        set l "Referenced by:"

        if {3 == [scan $i "%d %d %d" type num count]} {
	    set l "$l [lindex {{Reading} {Contig} {Free list}} $type] "
	    if {$count > 1} {
	        set l "${l}(loop) "
	    }
	    if {$type != 2} {
	        append l $num
	    }
        } else {
	    set l "$l Unknown!"
	    set type -1
        }

        set t [frame $f.referenced_by$c]
        pack [label $t.label$c -text $l] -side left
        pack $t -side top -fill x
    
        if {$type == 0} {
    	    pack [button $t.follow$c -bitmap @$env(STADTABL)/follow.bitmap \
    	        -command "AR_display $io $fr reading $num"] -side right
        } elseif {$type == 1} {
    	    pack [button $t.follow$c -bitmap @$env(STADTABL)/follow.bitmap \
    	        -command "AR_display $io $fr contig $num"] -side right
        }

       	incr c
    }    

    # Separator lookalike:
    pack [frame $f.separator -height 2 -bd 1 -relief groove] -side top -fill x
}

#-----------------------------------------------------------------------------
# Routines for displaying whole structures of information
#-----------------------------------------------------------------------------
proc AR_reading {io f o g} {
    AR_draw_str	  $io $o 16 $g	name
    AR_draw_str   $io $o 16 $g	trace_name
    AR_draw_str   $io $o  8 $g	trace_type
    AR_draw_foll  $io $f $o $g	left  reading
    AR_draw_foll  $io $f $o $g	right reading
    AR_draw_int   $io $o  $g	position        1 -1
    AR_draw_int   $io $o  $g	length          1 -1
    AR_draw_int   $io $o  $g	sequence_length 1 -1
    AR_draw_int   $io $o  $g	start           0 -1
    AR_draw_int   $io $o  $g	end             1 -1
    AR_draw_sense $io $o  $g	sense
    AR_draw_int   $io $o  $g	sequence        13 -1
    AR_draw_int   $io $o  $g	confidence      0 -1
    AR_draw_int   $io $o  $g	orig_positions  0 -1
    AR_draw_foll  $io $f $o $g	annotations annotation
    AR_draw_foll  $io $f $o $g  template template
    AR_draw_strand $io $o  $g	strand
    AR_draw_primer $io $o  $g	primer
    AR_draw_int   $io $o  $g    chemistry       0 -1
    AR_draw_foll  $io $f $o $g  notes note
}

proc AR_contig {io f o g} {
    AR_draw_foll $io $f $o $g left  reading
    AR_draw_foll $io $f $o $g right reading
    AR_draw_int  $io    $o $g length 1 -1
    AR_draw_foll $io $f $o $g annotations annotation
    AR_draw_foll $io $f $o $g  notes note
}

proc AR_annotation {io f o g} {
    AR_draw_type    $io    $o $g type
    AR_draw_int	    $io    $o $g position   1 -1
    AR_draw_int     $io    $o $g length     0 -1
    AR_draw_astrand $io    $o $g strand
    AR_draw_int     $io    $o $g annotation 0 -1
    AR_draw_foll    $io $f $o $g next annotation
}

proc AR_note {io f o g} {
    AR_draw_type    $io    $o $g type
    AR_draw_int     $io    $o $g ctime      0 -1
    AR_draw_int     $io    $o $g mtime      0 -1
    AR_draw_int     $io    $o $g annotation 0 -1
    AR_draw_foll    $io $f $o $g next note
    AR_draw_foll    $io $f $o $g prev note
    AR_draw_int     $io    $o $g prev_type  0 -1
}

proc AR_template {io f o g} {
    AR_draw_str  $io $o 16 $g name 
    AR_draw_int  $io    $o $g strands 1  2
    AR_draw_foll $io $f $o $g vector vector
    AR_draw_foll $io $f $o $g clone clone
    AR_draw_int  $io $o $g    insert_length_min 0 -1
    AR_draw_int  $io $o $g    insert_length_max 0 -1
}

proc AR_clone {io f o g} {
    AR_draw_str  $io $o 16 $g name
    AR_draw_foll $io $f $o $g vector vector
}

proc AR_vector {io f o g} {
    AR_draw_str $io $o 16 $g name
    AR_draw_int $io $o    $g level 0 -1
}

proc AR_database {io f o g} {
    AR_draw_int $io $o $g version	0 -1
    AR_draw_int $io $o $g maximum_db_size 0 -1
    AR_draw_int $io $o $g actual_db_size 0 -1
    AR_draw_int $io $o $g max_gel_len	0 -1
    AR_draw_int $io $o $g data_class	0 -1
    AR_draw_int $io $o $g num_contigs	0 -1
    AR_draw_int $io $o $g num_readings	0 -1
    AR_draw_int $io $o $g Nfreerecs	0 -1
    AR_draw_int $io $o $g freerecs	0 -1
    AR_draw_int $io $o $g Ncontigs	0 -1
    AR_draw_int $io $o $g contigs	0 -1
    AR_draw_int $io $o $g Nreadings	0 -1
    AR_draw_int $io $o $g readings	0 -1
    AR_draw_int $io $o $g Nannotations	0 -1
    AR_draw_int $io $o $g annotations	0 -1
    AR_draw_foll $io $f $o $g free_annotations annotation
    AR_draw_int $io $o $g Ntemplates	0 -1
    AR_draw_int $io $o $g templates	0 -1
    AR_draw_int $io $o $g Nclones	0 -1
    AR_draw_int $io $o $g clones	0 -1
    AR_draw_int $io $o $g Nvectors	0 -1
    AR_draw_int $io $o $g vectors	0 -1
    AR_draw_int $io $o $g contig_order	0 -1
    AR_draw_int $io $o $g Nnotes        0 -1
    AR_draw_int $io $o $g notes_a       0 -1
    AR_draw_foll $io $f $o $g notes note
    AR_draw_foll $io $f $o $g free_notes note
}

#-----------------------------------------------------------------------------
# Routines for accepting changes to a structure
#-----------------------------------------------------------------------------
proc AR_get_int {f l item} {
    upvar $l list
    keylset list $item "[entrybox_get $f._$item.int]"
}

proc AR_get_str {io f l item} {
    upvar $l list
    io_write_text $io [keylget list $item] "[entrybox_get $f._$item.str]"
}

proc AR_get_type {f l item} {
    upvar $l list
    
#    set a 32; #ASCII for space
#    set b 32
#    set c 32
#    set d 32
#
#    scan "[entrybox_get $f._$item.int]" "%c%c%c%c" a b c d
#    keylset list $item [expr ($a<<24)+($b<<16)+($c<<8)+$d]

     keylset list $item "[format %-4s [entrybox_get $f._$item.int]]"
}

proc AR_accept_reading {io f num} {
    # Fetch structure
    set l ""
    AR_get_int $f l name
    AR_get_int $f l trace_name
    AR_get_int $f l trace_type
    AR_get_int $f l left
    AR_get_int $f l right
    AR_get_int $f l position
    AR_get_int $f l length
    AR_get_int $f l sense
    AR_get_int $f l sequence
    AR_get_int $f l confidence
    AR_get_int $f l orig_positions
    AR_get_int $f l chemistry
    AR_get_int $f l annotations
    AR_get_int $f l sequence_length
    AR_get_int $f l start
    AR_get_int $f l end
    AR_get_int $f l template
    AR_get_int $f l strand
    AR_get_int $f l primer
    AR_get_int $f l notes

    # Write structure 
    io_write_reading $io $num $l

    # Write strings
    if {[keylget l name] >= 13} {
	io_write_reading_name $io $num "[entrybox_get $f._name.str]"
    }
    if {[keylget l trace_name] >= 13} {
        AR_get_str $io $f l trace_name
    }
    if {[keylget l trace_type] >= 13} {
        AR_get_str $io $f l trace_type
    }

    # Flush everything
    io_flush $io
}

proc AR_accept_contig {io f num} {
    # Fetch structure
    set l ""
    AR_get_int $f l left
    AR_get_int $f l right
    AR_get_int $f l length
    AR_get_int $f l annotations
    AR_get_int $f l notes

    # Write structure 
    io_write_contig $io $num $l

    # Flush everything
    io_flush $io
}

proc AR_accept_annotation {io f num} {
    # Fetch structure
    set l ""
    AR_get_type $f l type
    AR_get_int  $f l position
    AR_get_int  $f l length
    AR_get_int  $f l strand
    AR_get_int  $f l annotation
    AR_get_int  $f l next

    # Write structure 
    io_write_annotation $io $num $l

    # Flush everything
    io_flush $io
}

proc AR_accept_note {io f num} {
    # Fetch structure
    set l ""
    AR_get_type $f l type
    AR_get_int  $f l ctime
    AR_get_int  $f l mtime
    AR_get_int  $f l annotation
    AR_get_int  $f l next
    AR_get_int  $f l prev
    AR_get_int  $f l prev_type

    # Write structure 
    io_write_note $io $num $l

    # Flush everything
    io_flush $io
}

proc AR_accept_template {io f num} {
    # Fetch structure
    set l ""
    AR_get_int $f l name
    AR_get_int $f l strands
    AR_get_int $f l vector
    AR_get_int $f l clone
    AR_get_int $f l insert_length_min
    AR_get_int $f l insert_length_max

    # Write structure 
    io_write_template $io $num $l

    # Fetch and write strings
    if {[keylget l name] >= 13} {
        AR_get_str $io $f l name
    }

    # Flush everything
    io_flush $io
}

proc AR_accept_clone {io f num} {
    # Fetch structure
    set l ""
    AR_get_int $f l name
    AR_get_int $f l vector

    # Write structure 
    io_write_clone $io $num $l

    # Fetch and write strings
    if {[keylget l name] >= 13} {
        AR_get_str $io $f l name
    }

    # Flush everything
    io_flush $io
}

proc AR_accept_vector {io f num} {
    # Fetch structure
    set l ""
    AR_get_int $f l name
    AR_get_int $f l level

    # Write structure 
    io_write_vector $io $num $l

    # Fetch and write strings
    if {[keylget l name] >= 13} {
        AR_get_str $io $f l name
    }

    # Flush everything
    io_flush $io
}

proc AR_accept_database {io f num} {
    # Fetch structure
    set l ""
    AR_get_int $f l version
    AR_get_int $f l maximum_db_size
    AR_get_int $f l actual_db_size
    AR_get_int $f l max_gel_len
    AR_get_int $f l data_class
    AR_get_int $f l num_contigs
    AR_get_int $f l num_readings
    AR_get_int $f l Nfreerecs
    AR_get_int $f l freerecs
    AR_get_int $f l Ncontigs
    AR_get_int $f l contigs
    AR_get_int $f l Nreadings
    AR_get_int $f l readings
    AR_get_int $f l annotations
    AR_get_int $f l free_annotations
    AR_get_int $f l Ntemplates
    AR_get_int $f l templates
    AR_get_int $f l Nclones
    AR_get_int $f l clones
    AR_get_int $f l Nvectors
    AR_get_int $f l vectors
    AR_get_int $f l contig_order
    AR_get_int $f l Nnotes
    AR_get_int $f l notes_a
    AR_get_int $f l notes
    AR_get_int $f l free_notes

    # Write structure 
    io_write_database $io $l

    # Flush everything
    io_flush $io
}

#-----------------------------------------------------------------------------
# General user interface code
#-----------------------------------------------------------------------------
proc AR_display {io f name val} {
    global gap_fatal_errors
    set gap_fatal_errors 0
    global gap_defs

    if {[winfo exists $f.dia]} {destroy $f.dia}

    if {[winfo exists $f.tout]} {
        catch {set o [frame $f.tout2]}
        set old $f.tout
    } else {
        catch {set o [frame $f.tout]}
        set old $f.tout2
    }

    if {"$name" != "database"} {
        set l [io_read_${name} $io $val]
        label $o.title -text "${name} number $val" -bd 2 -relief raised
    } else {
        set l [io_read_${name} $io]
        label $o.title -text "${name} structure" -bd 2 -relief raised
    }

    pack $o.title -side top -fill x -anchor c

    if {"$l" != ""} {
	if {"$name" == "annotation"} {
	    AR_draw_ref $io $f $o $val
        }

        AR_${name} $io $f $o $l
        frame $o.but -bd 2 -relief raised
        button $o.but.accept -text "Accept" \
		-command "AR_accept $io $o $name $val"
        button $o.but.reload -text "Reload" \
		-command "AR_display $io $f $name $val"
        button $o.but.cancel -text "Cancel" -command "destroy $o"

        pack $o.but -side bottom -fill both
        pack $o.but.accept $o.but.reload $o.but.cancel -side left
    } else {
	bell
    }

    if {[winfo exists $old]} {destroy $old}
    pack $o -side bottom -fill both
}

proc AR_accept {io fr name num} {
    AR_accept_${name} $io $fr $num
#    destroy $fr
}

proc AR_which {io f name} {
    global gap_defs

    set db [io_read_database $io]
    if {"$name" == "contig" || "$name" == "reading"} {
        set n [keylget db num_${name}s]
    } else {
	set n [keylget db N${name}s]
    }

    if {[winfo exists $f.dia]} {destroy $f.dia}
    frame $f.dia
    entrybox $f.dia.e \
	-title "Which $name (1-$n)" \
	-command "AR_display $io $f $name" \
	-type "CheckIntRange 1 $n" \
	-default 1

    pack $f.dia -side top -after $f.mbar
    pack $f.dia.e -fill x
}

proc AR_quit {io t} {
    global AR_active
    global gap_fatal_errors
    global gap_defs

    if {[incr AR_active -1] == 0} {
        set gap_fatal_errors 1
       	grab release [grab current]
        destroy [keylget gap_defs DOCTOR_DATABASE.WIN]
	ContigInitReg $io
    } else {
        destroy $t
    }
}

proc AR_extend_many {io name count} {
    for {set i 0} {$i < $count} {incr i} {
	io_add_${name} $io
    }
    io_flush $io
}

proc AR_extend {io f name} {
    global gap_fatal_errors
    global gap_defs

    set gap_fatal_errors 0
    if [winfo exists $f.dia] {destroy $f.dia}
    frame $f.dia
    entrybox $f.dia.e \
	-title "Extend by how many?" \
	-command "entrybox_destroy $f.dia.e;AR_extend_many $io $name" \
	-type "CheckIntMin 0"\
	-default 1

    pack $f.dia -side top -after $f.mbar
    pack $f.dia.e -fill x
}

proc AR_delete_contig {io f} {
    global gap_fatal_errors
    global CurContig
    global gap_defs

    set gap_fatal_errors 0
    if [winfo exists $f.dia] {destroy $f.dia}
    frame $f.dia
    entrybox $f.dia.e \
	-title "Delete which contig?" \
	-type "CheckContigName $io" \
	-default "$CurContig"

    frame $f.dia.but
    button $f.dia.but.ok -text "Ok" -command "\
	delete_contig \
	    -io $io \
	    -contigs \[entrybox_get $f.dia.e\]; \
	destroy $f.dia"
    button $f.dia.but.cancel -text "Cancel" -command "destroy $f.dia"

    pack $f.dia -side top -after $f.mbar
    pack $f.dia.but.ok $f.dia.but.cancel -side left
    pack $f.dia.e $f.dia.but -side top -fill x
}

proc AR_fix_holes {io f} {
    global gap_fatal_errors
    global CurContig
    global gap_defs

    set gap_fatal_errors 0
    if [winfo exists $f.dia] {destroy $f.dia}
    frame $f.dia
    xyn $f.dia.e2 \
	-label "All contigs" \
	-orient horiz \
	-ycommand "entrybox_configure $f.dia.e -state disabled" \
	-ncommand "entrybox_configure $f.dia.e  -state normal"
    entrybox $f.dia.e \
	-title "Which contig?" \
	-type "CheckContigName $io" \
	-default "$CurContig"
    $f.dia.e2 set 1

    frame $f.dia.but
    button $f.dia.but.ok -text "Ok" -command "\
	if {\[$f.dia.e2 get\]} { \
            remove_contig_holes \
		-io $io \
		-contigs \[ListGet allcontigs\];\
        } else { \
            remove_contig_holes \
		-io $io \
		-contigs \[entrybox_get $f.dia.e\];\
        }; \
	destroy $f.dia"
    button $f.dia.but.cancel -text "Cancel" -command "destroy $f.dia"

    pack $f.dia -side top -after $f.mbar
    pack $f.dia.but.ok $f.dia.but.cancel -side left
    pack $f.dia.e $f.dia.e2 $f.dia.but -side top -fill x
}

proc AR_output_annos {io f} {
    global NGTag
    set type REPT
	
    if [winfo exists $f.dia] {destroy $f.dia}
    frame $f.dia
    frame $f.dia.x

    menubutton $f.dia.x.m -text "$type" -relief raised \
	-menu $f.dia.x.m.m -indicatoron 1
    menu $f.dia.x.m.m -tearoff 0
    for {set i 0} {$i < $NGTag(num_tags)} {incr i} {
        $f.dia.x.m.m add command \
            -label "$NGTag($i,tagid): $NGTag($i,tagtype)" \
            -command "$f.dia.x.m configure -text $NGTag($i,tagid)"
    }
    
    entrybox $f.dia.x.e \
	-title "Filename to save annotations to" \
	-type "CheckOutput" \
	-default anno_list

    frame $f.dia.but
    button $f.dia.but.ok -text "Ok" -command "\
	AR_output_annos2 $io \[entrybox_get $f.dia.x.e\] \
		\[$f.dia.x.m cget -text\]; destroy $f.dia"
    button $f.dia.but.cancel -text "Cancel" -command "destroy $f.dia"
    
    pack $f.dia -side top -after $f.mbar
    pack $f.dia.but.ok $f.dia.but.cancel -side left
    pack $f.dia.x.e $f.dia.x.m -side left -fill both
    pack $f.dia.x $f.dia.but -side top -fill both
}

proc AR_output_annos2 {io file type} {
    puts [set fd [open $file w]] [anno_list -io $io -type $type];
    close $fd
}

proc AR_delete_annos {io f} {
    if [winfo exists $f.dia] {destroy $f.dia}
    frame $f.dia

    entrybox $f.dia.e \
	-title "Filename of annotations to delete" \
	-type "CheckInput" \
	-default anno_list

    frame $f.dia.but
    button $f.dia.but.ok -text "Ok" \
	-command "AR_delete_annos2 $io \[entrybox_get $f.dia.e\]; \
		  destroy $f.dia"
    button $f.dia.but.cancel -text "Cancel" -command "destroy $f.dia"
    
    pack $f.dia -side top -after $f.mbar
    pack $f.dia.but.ok $f.dia.but.cancel -side left
    pack $f.dia.e $f.dia.but -side top -fill both
}

proc AR_delete_annos2 {io file} {
    delete_anno_list -io $io -annos [read [set fd [open $file r]]]
    close $fd
}

proc AlterRelationships {io} {
    global AR_counter
    global AR_active
    global gap_fatal_errors
    global gap_defs

    set tname [keylget gap_defs DOCTOR_DATABASE.WIN]

    vfuncheader "Doctor database"

    if {![quit_displays $io "doctor_database"]} {
	# Someone's too busy to shutdown?
	return
    }

    # Should this be blocking? I think so ...
    if {![winfo exists [keylget gap_defs DOCTOR_DATABASE.WIN]]} {
	frame $tname -width 0 -height 0 -bd 2
	pack $tname
    }

    set t $tname.alter_rel$AR_counter
    incr AR_counter
    if {[winfo exists $t]} {
	raise $t
	return
    }
    set gap_fatal_errors 0

    if {[xtoplevel $t -resizable 0] == ""} return
    wm title $t "Doctor Database"
    wm protocol $t WM_DELETE_WINDOW "AR_quit $io $t"
    incr AR_active

    # Create the menubar
    frame $t.mbar -bd 2 -relief raised
    set file $t.mbar.file
    menubutton $file -text "File" -menu $file.m
    menu $file.m
    $file.m add command -label "New" -command "AlterRelationships $io"
    $file.m add separator
    $file.m add command -label "Quit" -command "AR_quit $io $t"

    set inf $t.mbar.info
    menubutton $inf -text "Structures" -menu $inf.m
    menu $inf.m
    $inf.m add command -label "Database" \
	-command "AR_display $io $t database 1"
    $inf.m add command -label "Reading" \
	-command "AR_which $io $t reading"
    $inf.m add command -label "Contig" \
	-command "AR_which $io $t contig"
    $inf.m add command -label "Annotation" \
	-command "AR_which $io $t annotation"
    $inf.m add command -label "Template" \
	-command "AR_which $io $t template"
    $inf.m add command -label "Original clone" \
	-command "AR_which $io $t clone"
    $inf.m add command -label "Vector" \
	-command "AR_which $io $t vector"
    $inf.m add command -label "Note" \
	-command "AR_which $io $t note"

    set com $t.mbar.com
    menubutton $com -text "Commands" -menu $com.m
    menu $com.m
    $com.m add command -label "Check" \
	-command "check_database -io $io"
    $com.m add checkbutton -label "Ignore Check Database" \
	-variable ignore_checkdb
    $com.m add cascade -label "Extend structures" -menu $com.m.s
    $com.m add command -label "Delete contig" \
	-command "AR_delete_contig $io $t"
    $com.m add command -label "Fix contig holes" \
	-command "AR_fix_holes $io $t"
    $com.m add command -label "Reset contig order" \
	-command "reset_contig_order -io $io"
    $com.m add command -label "Output annotations to file" \
	-command "AR_output_annos $io $t"
    $com.m add separator
    $com.m add command -label "Delete annotations" \
	-command "AR_delete_annos $io $t"

    menu $com.m.s
    $com.m.s add command -label "Contig" \
	-command "AR_extend $io $t contig"
    $com.m.s add command -label "Reading" \
	-command "AR_extend $io $t reading"
    $com.m.s add command -label "Annotation" \
	-command "AR_extend $io $t annotation"
    $com.m.s add command -label "Note" \
	-command "AR_extend $io $t note"
    $com.m.s add command -label "Template" \
	-command "AR_extend $io $t template"
    $com.m.s add command -label "Clone" \
	-command "AR_extend $io $t clone"
    $com.m.s add command -label "Vector" \
	-command "AR_extend $io $t vector"

    set help $t.mbar.help
    menubutton $help -text "Help" -menu $help.m
    menu $help.m
    $help.m add command -label "Introduction" \
	-command "show_help gap4 {Doctor Database}"
    $help.m add command -label "Structure menu" \
	-command "show_help gap4 {Doctor-Structures}"
    $help.m add command -label "Extend structures" \
	-command "show_help gap4 {Doctor-Extend}"
    $help.m add command -label "Annotations" \
	-command "show_help gap4 {Doctor-Anno}"
    $help.m add command -label "Shift readings" \
	-command "show_help gap4 {Doctor-Shift}"
    $help.m add command -label "Delete contig" \
	-command "show_help gap4 {Doctor-Delete}"
    $help.m add command -label "Contig order" \
	-command "show_help gap4 {Doctor-Contig Order}"
	
    pack $file $inf $com -side left
    pack $help -side right
    pack $t.mbar -side top -fill both
    update idletasks
    grab $tname
}

#set auto_path "$env(GAP_LIBRARY) $auto_path"
#set io [open_db -name B0334 -version 0 -access rw]
#tk_focusFollowsMouse
#AlterRelationships $io
