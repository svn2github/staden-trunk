proc init_contig_order_list { } {
    global c_order_list c_direction_list
    
    set c_order_list ""
    set c_direction_list ""
}

proc create_contig_order_list { io contig_list direction_list} {
    global c_order_list c_direction_list
    
    set dir_list ""

    #need to change direction list into absolute eg left most reading doesn't
    #need complementing, right most reading if does need complementing
    set num_contigs [llength $contig_list]
    for {set i 0} {$i < $num_contigs} {incr i} {
	set dir [lindex $direction_list $i]
	set c_name [lindex $contig_list $i]
	set c_num [db_info get_contig_num $io $c_name]
	set c [io_read_contig $io $c_num]

	if {$dir == 1} {
	    set d +[keylget c length]
	} else {
	    set d -[keylget c length]
	}

	vmessage -nonewline "$c_name ($d) "
	    
	#check for single read contig and save the sense
	if {[keylget c left] == [keylget c right]} {
	    set r [io_read_reading $io [keylget c left]]
	    set sense [keylget r sense]
	    if {$dir == -1} {
		#change sense if needs complementing
		if {$sense == 1} {
		    set sense 0
		} else {
		    set sense 1
		}
	    }
	    lappend dir_list $sense
	} else {
	    lappend dir_list -1
	}

	if {$dir == -1} {
	    set right [rightmost_read -io $io -contig $c_num]
	    set contig_list [lreplace $contig_list $i $i $right]
	}
    }

    vmessage "\n"

    if {![info exists c_order_list]} {
	set c_order_list ""
    }
    lappend c_order_list $contig_list

    if {![info exists c_direction_list]} {
	set c_direction_list ""
    }
    lappend c_direction_list $dir_list

    #puts "DIRECTION LIST $direction_list"
}

proc contig_order_OkPressed { io t lb } {
    global c_direction_list read_only

    #puts $c_direction_list

    if {[$lb curselection] == ""} {
	tk_messageBox -icon warning -type ok -title "Contig order" \
		-message "No selection has been made" \
		-parent $t
	raise $t
 	
    } else {
	#complement contigs if possible
	set c_list [$lb get [$lb curselection]]
	set dir_list [lindex $c_direction_list [$lb curselection]]
	set cnt 0
	set comp_list ""

	foreach c_name $c_list {

	    set c_num [db_info get_contig_num $io $c_name]
	    
	    #need check for contig of only one reading
	    set c [io_read_contig $io $c_num]
	    if {[keylget c left] == [keylget c right]} {
		set r [io_read_reading $io [keylget c left]]
		if {[lindex $dir_list $cnt] != [keylget r sense]} {
		    lappend comp_list $c_name
		}
	    } else {
		set left [left_gel $io $c_num]
		if {[string compare $left $c_name] != 0} {
		    lappend comp_list $c_name
		}
	    }
	    incr cnt
	}

	if {!$read_only} {
	    complement_contig -io $io -contigs $comp_list
	}
	CreateTemplateDisplay $io [$lb get [$lb curselection]]
    }
}

proc contig_order_listbox { io } {
    global c_order_list

    if {[llength $c_order_list] == 0} {
	tk_messageBox -icon warning -type ok -title "Contig order" \
		-message "No contig order was found"
	
	return
    }

    set t .c_order_listbox

    if {[xtoplevel $t] == ""} return
    wm title $t "Order contigs"

    listbox $t.lb -xscrollcommand "$t.sb_h set" -yscrollcommand "$t.sb_v set" \
	    -selectmode single -width 80
    scrollbar $t.sb_h -orient horizontal -command "$t.lb xview"
    scrollbar $t.sb_v -orient vertical -command "$t.lb yview"
    okcancelhelp $t.ok_cancel \
	-ok_command "contig_order_OkPressed $io $t $t.lb"\
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Order-Contigs}" \
	-bd 2 -relief groove

    grid columnconfig $t 0 -weight 1
    grid rowconfig $t 0 -weight 1
    grid $t.lb -row 0 -column 0 -sticky nsew
    grid $t.sb_h -row 1 -column 0 -sticky ew
    grid $t.sb_v -row 0 -column 1 -sticky ns
    grid $t.ok_cancel -row 2 -column 0 -sticky w

    $t.lb delete 0 end

    foreach c_list $c_order_list {
	$t.lb insert end $c_list
    }

    bind $t.lb <<use>> "contig_order_OkPressed $io $t $t.lb"
}

