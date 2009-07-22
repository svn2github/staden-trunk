namespace eval prefinish {
    array set rules {}

    variable opt_tooltips 0
}

#-----------------------------------------------------------------------------
# Load the rule database and store it in rules
proc prefinish::load_rules {&rules} {
    global prefinish_defs
    upvar ${&rules} rules

    set rule_dir [keylget prefinish_defs RULE_DIR]
    foreach f [glob $rule_dir/*.rule] {
	set name [lindex [regexp -inline {(.*)\.rule$} $f] 1]
	interp create tmp_interp
	array set vars ""
	foreach v [tmp_interp eval info globals] {
	    set vars($v) 1
	}
	tmp_interp eval source [list $f]
	catch {unset rule}
	set rule(name) $name
	foreach v [tmp_interp eval info globals] {
	    if {[info exists vars($v)]} {
		continue
	    }
	    if {[tmp_interp eval info exists $v]} {
		set rule($v) [tmp_interp eval set $v]
	    }
	}
	set rules($rule(name)) [array get rule]
	interp delete tmp_interp
    }
}

#-----------------------------------------------------------------------------
# Construct the basic screen layout
proc prefinish::load_tooltips {fname} {
    variable tooltips

    set fd [open $fname]
    set state OPTION

    while {[gets $fd line] != -1} {
	for {set looping 1} {$looping > 0} {incr looping -1} {
	    switch $state {
		OPTION {
		    set opt [lindex $line 0]
		    set args [lrange $line 1 end]
		    set tooltips($opt) ""
		    set state HELP
		}
		
		HELP {
		    if {[regexp "^\[^ \t\]" $line]} {
			set state OPTION
			incr looping
		    } else {
			if {[string trim $line] == ""} {
			    append tooltips($opt) "\n\n"
			} else {
			    append tooltips($opt) "[string trim $line] "
			}
		    }
		}
	    }
	}
    }

    foreach tip [array names tooltips] {
	set tooltips($tip) [string trimright $tooltips($tip)]
    }
}

#-----------------------------------------------------------------------------
# Construct the basic screen layout
proc prefinish::main_gui {w} {
    variable rules
    global env

    load_tooltips $env(STADTCL)/prefinish/help

    if {$w == ""} {
	set wt .
    } else {
	set wt $w
    }
    wm geometry $wt 700x400

    # Create the main hierarchybook
    set h [hierarchybook $w.h]
    pack $h -fill both -expand 1

    # Add->Rule menu
    set_menu prefinish_menu
    foreach rule [lsort [array name rules]] {
	array set tmp $rules($rule)
	add_command Add.Rule.$tmp(name) 1 0 \
	    [list prefinish::add_rule $h $tmp(name) $rules($rule)]
    }

    # Create the menus
    global prefinish_menu
    menu $w.m
    $wt configure -menu $w.m
    create_menus $prefinish_menu $w.m

    # Always have a set of parameters to configure.
    add_param $h
    global $h.StartState
    set $h.StartState [prefin_save_recurse [$h root]]

    wm protocol $wt WM_DELETE_WINDOW [code exit $h]

    # Options
    variable opt_tooltips
    global prefinish_defs
    set opt_tooltips [keylget prefinish_defs SHOW_TOOLTIPS]
    ::prefinish::show_tooltips

    return $h
}

#-----------------------------------------------------------------------------
# Functions for adding and removing rules, passes and parameters.

# Add a generic tree node 'obj' to the hbook root
proc prefinish::add_node {h obj pos {node {}}} {
    if {$node == {}} {
	set node [$h root]
    }
    $node add $obj $pos
    $h register $obj
    $h draw -eventually
    $h raise $obj
}

# Finds the currently selected pass. We find any pass and then query it to
# know the last raised pass.
# If no current pass exists, create one.
# Returns the current pass.
proc prefinish::current_pass {h} {
    foreach child [[$h root] contents] {
	if {[$child isa Pass]} {
	    return [$child current_pass]
	}
    }

    # None found, so create one.
    return [add_pass $h]
}

# Finds the currently selected rule. We find any rule and then query it to
# know the last raised rule.
# Returns the current rule, or "" if none.
proc prefinish::current_rule {h} {
    foreach pass [[$h root] contents] {
	foreach rule [$pass contents] {
	    if {[$rule isa Rule]} {
		return [$rule current_rule]
	    }
	}
    }

    # None found
    return ""
}

# Adds a 'pass' to the hbook root
proc prefinish::add_pass {h} {
    set pass [Pass ::\#auto]
    add_node $h $pass end
    return $pass
}

# Adds global 'param's to either the hbook root or the current pass.
proc prefinish::add_param {h {node {}}} {
    set inherit_from ""

    if {$node == ""} {
	# What node to add to
	if {![[$h root] contains Params]} {
	    set node [$h root]
	} elseif {![[current_pass $h] contains Params]} {
	    set inherit_from [lindex [[$h root] find_by_type Params] 0]
	    set node [current_pass $h]
	} else {
	    bell
	    return
	}
    }

    set params [Params ::\#auto]
    add_node $h $params start $node

    $params add [set pw [ParamWalk ::\#auto]]
    $params add [set p3 [ParamPrimer3 ::\#auto]]
    $params add [set pc [ParamCosts ::\#auto]]
}

# Adds a rule and solution to the last selected pass.
proc prefinish::add_rule {h name ruleconf {pass {}}} {
    # Create a pass
    if {$pass == ""} {
	set pass [current_pass $h]
    }
    
    # Obtain rule/solution arguments
    set args {}
    set sol_args {}
    foreach {param value} $ruleconf {
	if {$param == "solution_defaults"} {
	    set sol_args $value
	    continue
	}
	lappend args -$param $value
    }
    regsub -all "\n" $sol_args  " " sol_args

    # Create rule
    set r [eval Rule ::\#auto $args]
    $pass add $r
    
    # Create the solution
    set s [eval Solution ::\#auto "$sol_args"]
    $r add $s

    # Display
    $h raise $r
    $h register $r
    $h draw  -eventually
}

# Adds a solution to the last selected rule.
proc prefinish::add_solu {h {rule {}}} {
    if {$rule == ""} {
	# Find the last selected rule
	set rule [current_rule $h]
	if {$rule == ""} {
	    bell
	    return
	}
    }
    
    # Create the solution
    set s [Solution ::\#auto]
    $rule add $s

    # Display
    $h raise $s
    $h register $s
    $h draw -eventually
}

# Removes the currently selected node, if possible.
proc prefinish::remove_node {h} {
    set sel [$h selection get]
    if {$sel == ""} {
	return
    }

    set cl [lindex [split [$sel info class] ::] end]
    if {[lsearch {Rule Pass Params Solution} $cl] == -1} {
	return
    }

    $sel delete

#    delete object $sel
}

#-----------------------------------------------------------------------------
# Save and Load functions
#
# The list we build is a series of {type params childen} tuples, where
# children may be a list of the same type. Eg:
#
# {  config
#    {-blah 1 -foo 2}
#    {
#       {  primer3
#          {-...}
#       }
#      {  costs
#          {-...}
#       }
#     }
# }
# {  pass1
#    {..}
# }
#
proc prefinish::prefin_save_recurse {node} {
    set str ""
    foreach child [$node contents] {
	if {[catch {$child info function save} err]} {
	    set save ""
	} else {
	    set save [$child save]
	}
	lappend str \
	    [list [$child info class] \
		 $save \
		 [prefin_save_recurse $child]]
    }
    return $str
}

proc prefinish::dump_config {str {depth 0}} {
    foreach child $str {
	set ind ""
	for {set i 0} {$i < $depth} {incr i} {
	    append ind "    "
	}
	puts ${ind}Node=[lindex $child 0]
	regsub -all "\n" [lindex $child 1] "\n${ind}  +conf=" conf
	puts "${ind}  +conf=$conf"
	dump_config [lindex $child 2] [expr {$depth+1}]
    }
}

#
# Returns the filename saved, or "" for none
#
proc prefinish::prefin_save {h {fname {}}} {
    global $h.StartState
    global $h.LastFilename

    if {$fname == "" && [info exists $h.LastFilename]} {
	set fname [set $h.LastFilename]
    }
    if {$fname == "__DEFAULT__"} {
	set fname ""
    }

    puts $fname

    if {$fname == {}} {
	set fname [tk_getSaveFile]
    }

    if {$fname == ""} {
	return ""
    }

    #puts :[[$h root] copy_node]:
    set str [prefin_save_recurse [$h root]]

    if {[catch {
	set fd [open $fname w]
	puts $fd $str
	close $fd} err]} {
	tk_messageBox \
	    -type ok \
	    -icon error \
	    -parent $h \
	    -message "Failed to save to '$fname': $err"
	return ""
    }

    set $h.LastFilename $fname
    set $h.StartState $str

    return $fname
}

proc prefinish::prefin_load_recurse {parent str {topparams {}}} {
    foreach child $str {
	set node [eval [lindex $child 0] ::\#auto [lindex $child 1]]
	$parent add $node
	prefin_load_recurse $node [lindex $child 2] $topparams
    }
}

#
# Returns 1 for OK
#         0 for cancelled (so don't continue)
#
proc prefinish::check_saved {h} {
    global $h.StartState
    global $h.LastFilename

    set CurrentState [prefin_save_recurse [$h root]]
    if {$CurrentState != [set $h.StartState]} {
	set answer [tk_messageBox \
			-type yesnocancel \
			-icon question \
			-parent $h \
			-message "Do you wish to save changes to this\
				      configuration?"]
	if {$answer == "cancel"} {
	    return 0
	} elseif {$answer == "yes"} {
	    if {[info exists $h.LastFilename]} {
		set res [prefin_save $h [set $h.LastFilename]]
	    } else {
		set res [prefin_save $h]
	    }
	    if {$res == ""} {
		return 0
	    }
	}
    }

    return 1
}

proc prefinish::exit {h} {
    if {[check_saved $h]} {
	::exit
    }
}

proc prefinish::prefin_load {h {fname {}}} {
    global $h.StartState
    global $h.LastFilename

    if {![check_saved $h]} {
	return
    }

    if {$fname == {}} {
	set fname [tk_getOpenFile]
    }

    if {$fname == ""} {
	return
    }

    if {[catch {
	set fd [open $fname]
	set str [read $fd]
	close $fd} err]} {
	tk_messageBox \
	    -type ok \
	    -icon error \
	    -parent $h \
	    -message "Failed to open '$fname': $err"
	return
    }

    set $h.StartState ""
    set $h.LastFilename $fname

    $h configure -cursor watch
    frame $h.tmp
    grab $h.tmp
    $h clear
    if {[catch {prefin_load_recurse [$h root] $str} err]} {
	tk_messageBox \
	    -type ok \
	    -icon error \
	    -parent $h \
	    -message "Failed to read '$fname': $err"
	destroy $h.tmp
	$h configure -cursor ""
	return
    }
    $h register [$h root]
    $h draw -eventually
    update idletasks
    destroy $h.tmp
    $h configure -cursor ""

    set $h.StartState [prefin_save_recurse [$h root]]
}

#-----------------------------------------------------------------------------
# Generates a hard-coded CLI application to run prefinish on the specified
# configuration. The benefit of this is that the application created is
# independent of Tk (graphics) and may be manually adjusted to better suit
# local requirements.

# class_bits is an array indexed with dependency names and having the bit
# number as the array value.
proc prefinish::fetch_dependencies {h &class_bits} {
    upvar ${&class_bits} class_bits

    # Collate rule dependency lists
    set dependencies {}
    foreach pass [[$h root] contents] {
	if {![$pass isa Pass]} {
	    continue;
	}

	foreach rule [$pass contents] {
	    if {![$rule isa Rule]} {
		continue
	    }

	    $rule save
	    
	    set dependencies [concat $dependencies [$rule dependencies]]
	}
    }

    # Sort and uniq them.
    # Assign bit values.
    set text "set class_bits {\n"
    set dependencies [lsort -uniq $dependencies]
    set bit_num 0
    foreach dep $dependencies {
	if {[lindex $dep 0] == "chemistry"} {
	    foreach ch [lrange $dep 1 end] {
		append text "    {$bit_num chemistry $ch}\n"
	    }
	} else {
	    append text "    {$bit_num $dep}\n"
	}
	set class_bits($dep) $bit_num
	incr bit_num
    }
    append text "}\n\n"

    return $text
}

# Generates the finishing rules tcl function, based on the class bits
# NOTE: problem 0 and 1 MUST be extend left and extend right problems.
# If the user has not chosen these, then we allocate them, but do not
# generate code to set those bits.
proc prefinish::generate_rules {h &class_bits &problem_bits &problem_pass
				&pscores &mscores &pwalk_prob_mask
				&mandatory} {
    upvar ${&class_bits} class_bits
    upvar ${&problem_bits} problem_bits
    upvar ${&problem_pass} problem_pass
    upvar ${&pscores} pscores
    upvar ${&mscores} mscores
    upvar ${&pwalk_prob_mask} pwalk_prob_mask
    upvar ${&mandatory} mandatory

    set text "proc finishing_rules {bits} \{\n"
    set rulenum 0
    set return_statement ""
    set passnum 0

    set pscores ""
    set mscores ""
    set pwalk_prob_mask 0
    set mandatory 0

    # Generate a list of rules
    set rules [list "" ""]
    foreach pass [[$h root] contents] {
	if {![$pass isa Pass]} {
	    continue;
	}

	foreach rule [$pass contents] {
	    if {![$rule isa Rule]} {
		continue
	    }

	    if {[$rule cget -name] == "Extend contig leftwards"} {
		set rules [lreplace $rules 0 0 $rule]
	    } elseif {[$rule cget -name] == "Extend contig rightwards"} {
		set rules [lreplace $rules 1 1 $rule]
	    } else {
		lappend rules $rule
	    }

	    set problem_pass($rule) $passnum
	}

	incr passnum
    }

    foreach rule $rules {
	if {$rule == ""} {
	    set problem_bits($rulenum) ""
	    incr rulenum
	    lappend pscores 0
	    lappend mscores 0
	    continue
	}

	lappend pscores [$rule cget -pweight]
	lappend mscores [$rule cget -mweight]
	if {[$rule cget -avoid_primers] == 1} {
	    set pwalk_prob_mask [expr {$pwalk_prob_mask | (1<<$rulenum)}]
	}

	if {[$rule cget -mandatory] == 1} {
	    set mandatory [expr {$mandatory | (1<<$rulenum)}]
	}

	append text "    \# Problem $rulenum: [$rule cget -name]\n"

	set expr [$rule cget -expression]
	set regres [regexp -all -inline {([^<]*)(<([^>]*)>)?} $expr]
	set newexpr "\[expr {!("
	foreach {_ st _ ex} $regres {
	    if {$ex != ""} {
		set bitnum [format "0x%04x" [expr {1<<$class_bits($ex)}]]
		switch $ex {
		    contig_left_end -
		    contig_right_end {
			append newexpr "${st}((\$bits & $bitnum)==0)"
		    }
		    default {
			append newexpr "${st}(\$bits & $bitnum)"
		    }
		}
	    } else {
		append newexpr $st
	    }
	}
	append newexpr ")}\]"
	append text "    set p$rulenum $newexpr\n\n"

	if {$return_statement == ""} {
	    set return_statement "    return \[expr \{\$p$rulenum"
	} else {
	    append return_statement " | (\$p$rulenum<<$rulenum)"
	}
	
	set problem_bits($rulenum) $rule

	incr rulenum
    }

    append text "$return_statement\}\]\n\}\n\n"

    return $text
}

# Generates the solution picking functions
proc prefinish::generate_solutions {&problem_bits &problem_pass} {
    upvar ${&problem_bits} problem_bits
    upvar ${&problem_pass} problem_pass

    set probcount [llength [array names problem_bits]]
    set used 0
    set passcount 0
    set solu_str ""

    do {
	append solu_str "proc find_solutions_$passcount {base_bits problem_bits} \{\n"
	append solu_str "    \# Optimisation - no problems => no solutions required\n"
	append solu_str "    if {\$problem_bits == 0} {return 0}\n\n"
	append solu_str "    set type 0;   # None\n"
	append solu_str "    set chem 17;  # BDv3 term\n"
	append solu_str "    set strand 0; # Any\n\n"

	for {set pnum 0} {$pnum < $probcount} {incr pnum} {
	    set rule $problem_bits($pnum)

	    if {$rule == ""} {
		if {$passcount == 0} {
		    incr used
		}
		continue
	    }

	    if {$problem_pass($rule) != $passcount} {
		continue
	    }

	    set solu [lindex [$rule contents] 0]

	    if {$solu != ""} {
		append solu_str "    \#Problem $pnum: [$rule cget -name]\n"
		set bit [format "0x%02x" [expr {1<<$pnum}]]
		append solu_str "    if {\[expr {\$problem_bits&$bit}\]} \{\n"
		append solu_str "        set type \[expr {\$type|[$solu experiments]}\]\n"
		append solu_str "        set chem [$solu cget -chemistry]\n"
		if {[$solu cget -strand] != 0} {
		    append solu_str "        set strand [$solu cget -strand]\n"
		}
		append solu_str "    \}\n\n"
	    }

	    incr used
	}

	append solu_str {    return [expr $type | ($strand << 16) | ($chem << 24)]}
	append solu_str "\n\}\n\n"

	incr passcount
    } while {$used != $probcount}

    return $solu_str
}

#
# Creates a global_params procedure containing the top-level parameter
# definitions. This is called from within each pass before pass specific
# parameters are applied.
#
proc prefinish::global_params {h pscores mscores pwalk_prob_mask mandatory} {
    set master_conf ""

    foreach node [[$h root] contents] {
	if {[$node isa Params]} {
	    set master_conf "    \$f configure"
	    set switches "[$node save]"
	    foreach subnode [$node contents] {
		append switches " [$subnode save]"
	    }

	    foreach {var val} $switches {
		append master_conf " \\\n        [list $var] [list $val]"
	    }
	    append master_conf " \\\n        -pscores [list $pscores]"
	    append master_conf " \\\n        -mscores [list $mscores]"
	    append master_conf " \\\n        -pwalk_prob_mask $pwalk_prob_mask"
	    append master_conf " \\\n        -prob_mandatory $mandatory"
	    append master_conf "\n"
	}
    }

    return "proc global_params {f} {\n$master_conf}\n\n"
}

#
# Loops through the main tree performing the configurations and pass info
#
proc prefinish::generate_process_proc {h} {
    # Generate a list of rules
    set rules [list "" ""]
    set passnum 0

    set conf "\nproc process_contig {io f &opt cnum start end class_bits} \{\n"
    append conf "    upvar \$\{&opt\} opt\n\n"

    append conf "    set c \[io_read_contig \$io \$cnum\]\n"
    append conf "    set tags \"\"\n"
    append conf "    set contig \"#\[keylget c left\]\"\n\n"
    append conf "    puts \"### CONTIG ID \$contig  (=\$cnum \[left_gel \$io \$cnum\]) ###\"\n\n"

    append conf "    # Reconfigure the finish object to work on this specific contig\n"
    append conf "    \$f configure \\\n"
    append conf "	-io \$io \\\n"
    append conf "	-contig \"{\$contig \$start \$end}\"\n\n"

    append conf "    # Classify the bases to produce a bit-pattern\n"
    append conf "    \$f classify \\\n"
    append conf "	-bits \$class_bits\n\n"

    foreach node [[$h root] contents] {
	if {![$node isa Pass]} {
	    continue
	}

	append conf "    # -- Pass [expr {$passnum+1}]\n\n"
	append conf "    global_params \$f\n\n"

	foreach child [$node contents] {
	    if {[$child isa Params]} {
		set sub_conf ""
		set switches [$child save]
		foreach subnode [$child contents] {
		    append switches " [$subnode save]"
		}

		foreach {var val} $switches {
		    if {$val == "inherited"} {
			continue
		    }
		    append sub_conf " \\\n        [list $var] [list $val]"
		}

		if {$sub_conf != ""} {
		    append conf "    \$f configure$sub_conf\n\n"
		}
	    }
	}

	append conf "    # Identify problems and solutions from bit-classifications\n"
	append conf "    \$f find_problems \\\n"
	append conf "        -problem_command finishing_rules \\\n"
	append conf "        -solution_command find_solutions_$passnum \\\n"
	append conf "        -tag_types {MASK CVEC}\n\n"

	append conf "    if {\$opt(-dump_problems) != \"\"} \{\n"
	append conf "        dump_problem \$io \$f \$opt(-dump_problems) \$cnum\n"
	append conf "    \}\n\n"

	append conf "    if {\$opt(-skip_solutions)} \{\n"
	append conf "        return\n"
	append conf "    \}\n\n"
	
        append conf "    # Produce solutions\n"
        append conf "    append tags \" \[\$f implement_solutions \\\n"
        append conf "		    -tag_types {OLIG PRIM MASK CVEC}\]\"\n\n"

	incr passnum
    }

    append conf "    if {\$opt(-add_tags)} \{\n"
    append conf "	puts tags=\"\$tags\"\n"
    append conf "	add_tags -io \$io -tags \$tags\n"
    append conf "    \}\n"
    append conf "    set fd \[open tags a\]\n"
    append conf "    puts \$fd \$tags\n"
    append conf "    close \$fd\n\n"
    append conf "    flush stdout\n\}\n\n"

    return $conf
}

proc prefinish::generate_save {h {fname {}}} {
    if {$fname == ""} {
	set fname [tk_getSaveFile]
	if {$fname == ""} {
	    return
	}
    }

    set str [generate_app $h]

    if {[catch {
	set fd [open $fname w]
	puts $fd $str
	close $fd} err]} {
	tk_messageBox \
	    -type ok \
	    -icon error \
	    -parent $h \
	    -message "Failed to save to '$fname': $err"
    }

    catch {file attributes $fname -permission a+xr}
}

#
# IMPORTANT: Problem bit 0 and 1 are hard-coded to indicate extend-left and
# extend-right ends of contig. If the user has not selected these problem type
# then these bits should be left as zero.
#
# Returns a string containing the app.
#
proc prefinish::generate_app2 {h} {
    set str ""
    set class_str [fetch_dependencies $h class_bits]

    append str [generate_rules $h \
		    class_bits \
		    problem_bits \
		    problem_pass \
		    pscores \
		    mscores \
		    pwalk_prob_mask \
		    mandatory]

    append str [generate_solutions \
		    problem_bits \
		    problem_pass]

    append str $class_str

    append str [global_params $h \
		    $pscores \
		    $mscores \
		    $pwalk_prob_mask \
		    $mandatory]

    append str [generate_process_proc $h]

    return $str
}

proc prefinish::generate_app {h} {
    global env

    set str "\#!/bin/sh\n"
    append str "\#\\\nexec stash \"\$0\" \${@+\"\$@\"} || exit 1\n\n"

    append str [generate_app2 $h]
    
    set fd [open $env(STADTCL)/prefinish/args.template]
    append str "[read $fd]\n"
    close $fd

    return $str
}

# -----------------------------------------------------------------------------
# Checks the opt_tooltips variable to toggle tooltips on and off
proc ::prefinish::show_tooltips {} {
    variable opt_tooltips

    if {$opt_tooltips == 1} {
	::tooltip::enable
    } else {
	::tooltip::disable
    }
}

# -----------------------------------------------------------------------------
# Saves any user-adjustable options
proc ::prefinish::save_options {} {
    variable opt_tooltips
    global prefinish_defs env

    keylset prefinish_defs SHOW_TOOLTIPS $opt_tooltips
    update_defs prefinish_defs $env(HOME)/.prefinishrc SHOW_TOOLTIPS
}

# -----------------------------------------------------------------------------
#
# Main initialisation function. Call this to autoload the rest
#
proc ::prefinish::init {} {
    variable rules

    load_rules rules
}
