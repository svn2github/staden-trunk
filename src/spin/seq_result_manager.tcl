#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
#
# Creates a results list window
# 
proc seq_result_list_create {t help_cmd} {

    listbox $t.l -yscrollcommand "$t.s set" -width 37
    scrollbar $t.s -orient vertical -command "$t.l yview"
    frame $t.bb
    button $t.bb.q -text "OK" -command "destroy $t"
    button $t.bb.h -text "Help" -command "$help_cmd"

    bind $t.l <Motion> {
	%W selection clear 0 end
	%W selection set [%W index @%x,%y]
	%W selection anchor [%W index @%x,%y]
    }

    pack $t.bb -side bottom -fill x
    pack $t.bb.q -side left -expand 1
    pack $t.bb.h -side right -expand 1
    pack $t.s -side right -fill y
    pack $t.l -side left -expand 1 -fill both
}

#destroy a results list window
proc seq_result_list_destroy { t } {
    if {[winfo exists $t]} {
	destroy $t
    }
}

#
# Updates the results list window
# 
proc seq_result_list_update { t } {

    #puts seq_result_list_update
    if {![winfo exists $t]} {
	return
    }
    #want to do list cos of Results menu even in results manager not yet 
    # created
    # Grab the registration list.
    set results [seq_result_names]

    # Clear and remove old binding
    $t.l delete 0 end
    bind $t.l <<menu>> {}

    foreach i $results {
	$t.l insert end $i
    }
    # Reenable our binding with our new list.
    bind $t.l <<menu>> {result_list_popup %W [%W index @%x,%y] %X %Y}
}

proc get_result_id {text} {

    regsub {.*#([^)]*).*} $text {\1} id
    return $id
}

#
# Displays a menu from a given result list item
#
proc result_list_popup_single {id ops m} {
    set count 0
    foreach i $ops {
	if {$i == "SEPARATOR"} {
	    $m add separator
	} else {
	    if {$i != "PLACEHOLDER"} {
	        $m add command \
	            -label $i \
	            -command "seq_invoke_op -index $id -job $count"
	    }
	    incr count
	}
    }
}

proc result_list_popup {lbox index x y} {
    if [winfo exists $lbox.m] {destroy $lbox.m}

    #get id number from list item
    #regsub {.*#([^)]*).*} [$lbox get $index] {\1} id
    set id [get_result_id [$lbox get $index]]
    set m [create_popup $lbox.m [$lbox get $index]]

    set ops [seq_get_ops -index $id]
    set count 0
    foreach i $ops {
	if {$i == "SEPARATOR"} {
	    $m add separator
	} else {
	    if {$i != "PLACEHOLDER"} {
	        $m add command \
	            -label $i \
	            -command "destroy $m; \
		              seq_invoke_op -index $id -job $count"
	    }
	    incr count
	}
    }

    $lbox selection clear 0 end
    $lbox selection set $index
    tk_popup $m $x $y
}

proc seq_result_menu_update { r_win results } {
    global nip_defs

    #set raster [keylget nip_defs RASTER.WIN]
    set r $r_win.menubar.results.opts
    #set r $r_win.menubar.[menu_path {Results}]

    #check the dot plot exists
    if {![winfo exists $r]} {
	return
    }

    $r delete 0 end

    set cnt 0

    foreach res $results {

	#puts "RESULTS $res"
	#get id number from current array item
	set id [get_result_id $res]
	#puts "index $id"
	set colour [lindex [lindex [raster_getconfig -index $id] 0] 1]

	#destroy cascade menu
	if {[winfo exists $r.ops$cnt]} {
	    destroy $r.ops$cnt
	}
	$r add cascade -label $res -menu $r.ops$cnt -foreground $colour \
		-activeforeground $colour
	menu $r.ops$cnt -tearoff 0

	#get id number from list item
	regsub {.*#([^)]*).*} $res {\1} id

	set ops [seq_get_ops -index $id]
	set count 0
	foreach i $ops {
	    if {$i == "SEPARATOR"} {
		$r.ops$cnt add separator
	    } else {
		if {$i != "PLACEHOLDER"} {
		    $r.ops$cnt add command \
			    -label $i \
			    -command "\
			    seq_invoke_op -index $id -job $count"
		}
		incr count
	    }
	}
	incr cnt
    }
}
   
proc seq_result_keybox_popup {r id } {
    
    set ops [seq_get_ops -index $id]
    set count 0
    foreach i $ops {
	if {$i == "SEPARATOR"} {
	    $r add separator
	} else {
	    if {$i != "PLACEHOLDER"} {
		$r add command \
		    -label $i \
		    -command "\
			    seq_invoke_op -index $id -job $count"
	    }
	    incr count
	}
    }
}


proc seq_result_keybox_update {r_win result_id results r} {
    global nip_defs

    #check the dot plot exists
    if {![winfo exists $r]} {
	return
    }

    $r delete 0 end

    set cnt 0

    foreach res $results {

	#get id number from list item
	regsub {.*#([^)]*).*} $res {\1} id
	
	#only show specific results. If result_id == -1 -> all
        if {$result_id > -1} {
	    
	    if {$id != $result_id} {
		continue
	    }
	    seq_result_keybox_popup $r $id
	    return
	}
	
	set colour [lindex [lindex [raster_getconfig -index $id] 0] 1]
	
	#destroy cascade menu
	if {[winfo exists $r.ops$cnt]} {
	    destroy $r.ops$cnt
	}
	$r add cascade -label $res -menu $r.ops$cnt -foreground $colour \
	    -activeforeground $colour
	menu $r.ops$cnt -tearoff 0
	
	set ops [seq_get_ops -index $id]
	set count 0
	foreach i $ops {
	    if {$i == "SEPARATOR"} {
		$r.ops$cnt add separator
	    } else {
		if {$i != "PLACEHOLDER"} {
		    $r.ops$cnt add command \
			-label $i \
			-command "\
			    seq_invoke_op -index $id -job $count"
		}
		incr count
	    }
	}
	incr cnt
    }
}