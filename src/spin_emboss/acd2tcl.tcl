#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 2001. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#
# Parser for ACD -> Tcl
#
# Reads ACD from argv1 and writes itkwish code to stdout.
# Also appends to the "menu_file" file.
#
# Recommended use:
#    /bin/rm menu_file
#    for i in acd/*.acd
#    do
#        a=`echo $i | sed 's:acd/::'`
#        echo $a
#        itkwish3.1 parseB.tcl acd/$a > acdtcl/$a
#    done
#    stash create_menu.tcl > menu

#source $env(STADTABL)/shlib.conf
#tkinit
#wm withdraw .

#-----------------------------------------------------------------------------
# The lexical analyser

# A lex style rule set.
#
# Each line consists of regexp token and a regsub list
# The regexps are searched in order, terminating at the first one that
# matches. Upon a match the token is stored along with the matching text.
# Specifying a blank token will simply throw away the matching data.
# The regsub list may be used to edit the matching text before storing it on
# the tokenised list. Specifying a blank regsub list will just store the
# unedited text with the token.
# Tokens of + and - are a mechanism for introducing blocks or hierarchies.
# + produces a token named "BLOCK" with a value of all tokenised text up to
# the next nested - token.

set tlist {
    {^.(.*).$} {\1}
    {\\[ \n\r]+} {\\n}
    {[ \n\r\t]+} { }
    {\\n} "\n"
    {\\(.)} {\1}
}

set rules [format {
    {#!}			{}		{}
    {#[^\n]*}			{}		{}
    {\[}			+		{}
    {\]}			-		{}
    {[$@][^] \t\n]+}		EXPR		{}
    {"(\\.|[^"\\])*"}		STRING		{%s}
    {'(\\.|[^'\\])*'}		STRING		{%s}
    {<(\\.|[^>\\])*>}		STRING		{%s}
    {\{(\\.|[^\}\\])*\}}	STRING		{%s}
    {var(iable)?:}		VAR		{{^(.*).$} {!\1}}
    {[^ \t\n:=]+[ \t\n]*[:=]}	ID		{{[ \t\n]*[:=]$} {}}
    {[^] \t\n:]+}		WORD		{}
    {[ \t\n]+}			{}		{}
    {:}				COLON		{}
    {.}				UNKNOWN		{}
} $tlist $tlist $tlist $tlist $tlist $tlist]

# Tokenises a string bases on a defined set of rules.
# Returns a list of {token line_num text} tuple.
proc tokenise_string {str rules} {
    set tokenised ""
    set line 1
    while {$str != ""} {
	foreach {r t p} $rules {
	    if {[regexp "^$r" $str match]} {
		incr line [expr [llength [split $match "\n"]]-1]
		set end [expr {[string length $match]-1}]
		set str [string replace $str 0 $end]
		if {$t != ""} {
		    foreach {l r} $p {
		        regsub -all $l $match $r match
		    }
		    if {[string match "+" $t]} {
		        append tokenised " \{BLOCK $line \{"
		    } elseif {[string match "-" $t]} {
		        append tokenised " \}\} "
		    } else {
		        append tokenised "[list [list $t $line $match]] "
		    }
		}
		break
	    }
	}
    }
    return $tokenised
}

proc print_tokens {toks {depth 0}} {
    foreach t $toks {
        foreach {l r} $t {}
        if {$l == "BLOCK"} {
	    incr depth
	    print_tokens $r $depth
	    incr depth -1
        } else {
            for {set i 0} {$i < $depth} {incr i} {
	        puts stderr -nonewline "    "
	    }
            puts stderr $t
	}
    }
}

#-----------------------------------------------------------------------------
# The parser, consisting of a parser for the tokenised acd and a conversion
# from acd expression syntax to tcl expr syntax.

# Converts acd strings, containing mixtures of real strings and acd
# expressions, into a tcl string that can be expanded using "subst $str".
proc convert_string {code dependvar} {
    upvar $dependvar depend
    set depend ""

    set brackets 0
    set newword 1
    set inexpr 0
    set text ""
    set char ""
    foreach next [split "$code " ""] {
        if {$char == ""} {
	    set char $next
	    continue
	}
	switch -- $char {
	    "\"" -
	    "\\" -
	    {[} -
	    {]} {
		if {!$inexpr} {
		    append text "\\$char"
		} else {
		    append expr "\\$char"
		}
	    }
	    "@" -
	    "$" {
		if {!$inexpr} {
		    if {$next == "("} {
			set brackets 0
			set inexpr 1
			set expr $char
		    } else {
			append text "\\$char"
		    }
		} else {
		    append expr "$char"
		}
	    }
	    "(" {
		if ($inexpr) {
		    incr brackets
		    append expr $char
		} else {
		    append text $char
		}
	    }
	    ")" {
		if ($inexpr) {
		    incr brackets -1
		    append expr $char
		    if {$brackets == 0} {
			set inexpr 0
			append text [convert_acd $expr dummy]
			set depend [concat $depend $dummy]
		    }
		} else {
		    append text $char
		}
	    }
	    default {
		if {!$inexpr} {
		    append text $char
		} else {
		    append expr $char
		}
	    }
	}
        set char $next
    }
    return $text
}

# Coverts acd statements to Tcl expressions
proc convert_acd {code varname} {
    upvar $varname vname

    if {![regexp {^[@$]} $code]} {
	return [list $code]
    }

    # Remove leading @
    regsub -all {@\(} $code {(} code
    
    set vars ""
    # Replace $(name) with $vars(name)
    set nul [binary format c 0]
    while {[regexp {\$\(([^)]*)\)} $code dummy name]} {
	regsub {\$\(([^)]*)\)} $code "\$vars($nul\\1$nul)" code
    }

    # Lowercase variable names in expressions. We've surrounded them with
    # nuls so that we can identify where they are.
    set newcode ""
    foreach {text var} [split $code $nul] {
        set var [string tolower $var]
        if {$var != ""} {
	     lappend vars $var
        }
        append newcode $text$var
    }
    set code $newcode

    if {$vars != ""} {
        set vname $vars
    }

    # Replace True and False
    regsub -all True  $code 1 code
    regsub -all False $code 0 code

    # Convert & to && and | to ||
    regsub -all {[&|]} $code {&&} code

    # Quote literals after ==, !=, ? :, < and > operators.
    regsub -all {([=?:<> ] *)([a-zA-Z_][a-zA-Z_0-9]*)} $code {\1"\2"} code

    # Switch statements - recognised by a single =
    # Replace with a series of ?: operators.
    if {[regexp {(\$vars\([^)]*\)) *=[^=]([^)]*)} $code dummy var sw]} {
        set newstr ""
	regsub -all {:} $sw { :} sw
        regsub -all { +} $sw { } sw
        set finished 0
        foreach {left colon right} [split $sw " "] {
            if {$left == {"else"}} {
		append newstr $right
		set finished 1
		break
	    } else {
		append newstr "$var == $left ? $right : "
	    }
        }
        if {!$finished} {
            append newstr {""}
        }
        regsub {(\$vars\([^)]*\)) *=[^=]([^)]*)} $code $newstr code
    }

    # is: statements (undocumented - assume to be "is set"
    regsub -all {is:} $code {"" != } code

    return "\[expr {$code}\]"
}

proc parse_error {msg} {
    puts stderr $msg
    return -code return
}

proc next_token {token} {
    upvar $token t
    set tok [lindex $t 0]
    set t [lrange $t 1 end]
    return $tok
}

proc peek_token {token} {
    upvar $token t
    return [lindex $t 0]
}

# Converts y/n into 1/0
proc convert_bool {text} {
    if {[string match {[YyTt]*} $text]} {
	return 1
    } elseif {[string match {[NnFf]*} $text]} {
	return 0
    } else {
	return $text
    }
}

# This parses the tokens list and generates an array of items to create.
# The array is index by megawidget name to generate using normal tk pathnames
# (starting at $prefix). The contents of array($pathname) is the widget type
# to create. The array is also populated by array($pathname.$variable) which
# holds information on the options supplied to $pathname. In contents of these
# array elements is the value of the associated variable. Finally an array
# pathname element starting with ! indicates a variable (eg .!var.options).
proc parse_acd {tokens aname {level 0} {prefix {}}} {
    upvar $aname data
    set next_prefix ""
    set item_list ""
    set item_list_ind ORDER:$prefix
    
    while {1} {
	set id_t ""
	foreach {id_t id_l id_v} [next_token tokens] break
	switch -- $id_t {
	    "" {
		break
	    }
	    VAR {
		foreach {name_t name_l name_v} [next_token tokens] break
		foreach {val_t val_l val_v} [next_token tokens] break
		set data($prefix$name_v) var
		set data($prefix$name_v.default) \
		    [convert_string $val_v data($prefix$name_v.depend)]
		lappend item_list $prefix$name_v
	    }
	    WORD -
	    ID {
		# Expand up shortenings, so "def" becomes "default" etc.
		foreach word {application information default required \
			      optional expected documentation outfile \
			      parameter needed delimiter codedelimiter \
			      values selection minimum maximum dirlist} {
		    if {[string match -nocase ${id_v}* $word]} {
			set id_v $word
			break
		    }
		}
		# Treat "parameter" as "required"
		if {$id_v == "parameter"} {
		    set id_v required
		}
		foreach {val_t val_l val_v} [next_token tokens] break

		# Lowercase the type name:
		# Eg "int: Foo [...]" to "int: foo [...]".
		if {$level == 0} {
		    set val_v [string tolower $val_v]
		}
		if {$id_v == "endsection"} {
		    set data(END$prefix$val_v) endsection
		    lappend item_list END$prefix$val_v
		}
		set next_t ""
		foreach {next_t next_l next_v} [peek_token tokens] break
		if {$next_t == "BLOCK"} {
		    set data($prefix$val_v) $id_v
		    lappend item_list $prefix$val_v
		} else {
		    set data(${prefix}$id_v) \
			[convert_string $val_v data(${prefix}$id_v.depend)]
		    set data(${prefix}$id_v.) ""; # Marks a non default value

		    if {$id_v == "prompt" && \
			    ![info exists data(${prefix}information.)]} {
			set data(${prefix}information) $data(${prefix}$id_v)
		    } elseif {$id_v == "needed"} {
			set data(${prefix}$id_v) \
			    [convert_bool $data(${prefix}$id_v)]
		    }

		    # FIXME: Hack for now
		    # required expressions => needed expression
		    if {$id_v == "required" && \
			    [string match {[@$](*} $val_v] && \
			    ![info exists data(${prefix}needed.)]} {
			set data(${prefix}needed) $data(${prefix}$id_v)
			set data(${prefix}needed.depend) $data(${prefix}$id_v.depend)
		    }

		    if {"$prefix" != "" && \
			    $data([string trimright $prefix .]) == "bool" && \
			    $id_v == "default"} {
			set data(${prefix}$id_v) \
			    [convert_bool $data(${prefix}$id_v)]
		    }
		}
		
		if {$level == 0 && ![info exists data($prefix$val_v._done)]} {
		    default_all data $prefix$val_v
		    set type [string tolower $id_v]
		    catch {default_$type data $prefix$val_v}
		}

		set next_prefix $prefix$val_v.
	    }
	    BLOCK {
		parse_acd $id_v data [expr $level+1] $next_prefix
	    }
	    default {
		if {[string length $id_v] > 20} {
		    set first [string range $id_v 0 20]...
		} else {
		    set first $id_v
		}
		puts "Unexpected token $id_t at line $id_l: $first"
	    }
	}
    }
    set data($item_list_ind) $item_list
}

proc print_parse_list {list {depth 0}} {
    foreach item $list {
	puts [format "%*s%s %s" \
		  $depth "" [lindex $item 0] [lindex $item 1]]
	set item [lrange $item 2 end]
	if {$item != ""} {
	    incr depth 4
	    print_parse_list $item $depth
	    incr depth -4
	}
    }
}


#-----------------------------------------------------------------------------
# Default values for ACD types - called from the parser.
#

proc set_def {dname index val} {
    upvar $dname data
    if {![info exists data($index)]} {
	set data($index) $val
    }
}

proc default_all {dname prefix} {
    upvar $dname data
    set data($prefix._done) 1
    set data($prefix.default)     ""
    set data($prefix.optional)	   0
    set data($prefix.required)    0
    set data($prefix.needed)      1
}

proc default_list {dname prefix} {
    upvar $dname data
    set data($prefix.default)     "A"
    set data($prefix.values)      "A:list"
    set data($prefix.delimiter)  ";"
    set data($prefix.codedelimiter)  ":"
    set data($prefix.minimum)	1
    set data($prefix.maximum)	1
    set data($prefix.header)   "Select several"
    set data($prefix.information) $prefix
}

proc default_selection {dname prefix} {
    upvar $dname data
    set data($prefix.default)     "list"
    set data($prefix.values)      "list"
    set data($prefix.delimiter)  ";"
    set data($prefix.minimum)	1
    set data($prefix.maximum)	1
    set data($prefix.header)   "Select several"
}

proc default_sequence {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Sequence"
    set data($prefix.type)        "any"
}

proc default_seqall {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Sequence"
    set data($prefix.type)        "any"
}

proc default_section {dname prefix} {
    upvar $dname data
    set data($prefix.information) ""
    set data($prefix.type)        "frame"
    set data($prefix.book)        ""
    set data($prefix.side)	  "top"
    set data($prefix.border)	  "1"
}

proc default_int {dname prefix} {
    upvar $dname data
    set data($prefix.information) $prefix
    set data($prefix.default)     "0"
}

proc default_integer {dname prefix} {
    upvar $dname data
    set data($prefix.information) $prefix
    set data($prefix.default)     "0"
}

proc default_float {dname prefix} {
    upvar $dname data
    set data($prefix.information) $prefix
    set data($prefix.default)     "0.0"
}

proc default_array {dname prefix} {
    upvar $dname data
    set data($prefix.information) $prefix
    set data($prefix.default)     "0.0"
}

proc default_dirlist {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Directory name"
    set data($prefix.default)     ""
}

proc default_outfile {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Output filename"
    set data($prefix.default)     "%N%C.out"
}

proc default_datafile {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Data filename"
}

proc default_infile {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Input filename"
}

proc default_directory {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Input directory"
}

proc default_outdir {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Output directory"
}

proc default_seqoutall {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Filename"
    set data($prefix.default)     "sequence%C.out"
}

proc default_features {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Features"
}

proc default_seqoutset {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Filename"
    set data($prefix.default)     "sequence%C.out"
}

proc default_seqset {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Filename"
    set data($prefix.totweight)	  0; #Unknown!
    set data($prefix.type)        "any"
}

proc default_featout {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Feature output filename"
    set data($prefix.default)     "feature%C.out"
}

proc default_outcodon {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Codon output filename"
    set data($prefix.default)     "codon%C.out"
}

proc default_report {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Report output filename"
    set data($prefix.default)     "report%C.out"
}

proc default_align {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Alignment output filename"
    set data($prefix.default)     "align%C.out"
}

proc default_seqout {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Sequence output filename"
    set data($prefix.default)     "sequence%C.out"
}

proc default_matrix {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Matrix filename"
    set data($prefix.protein)     1
    set data($prefix.default)     {[expr {$vars(%W.protein)?"EBLOSUM62":"EDNAMAT"}]}
    set data($prefix.default.depend)     {%W.protein}
}

proc default_matrixf {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Matrix filename"
    set data($prefix.protein)     1
    set data($prefix.default)     {[expr {$vars(%W.protein)?"EBLOSUM62":"EDNAMAT"}]}
    set data($prefix.default.depend)     {%W.protein}
}

proc default_codon {dname prefix} {
    upvar $dname data
    set data($prefix.information) "Codon usage table"
    set data($prefix.default)     Ehum.cut
}

proc default_string {dname prefix} {
    upvar $dname data
    set data($prefix.information) $prefix
    set data($prefix.default)     ""
}

proc default_xygraph {dname prefix} {
    upvar $dname data
    set data($prefix.default)     "data"
    set data($prefix.information) "Graphics output format"
}

proc default_graph {dname prefix} {
    upvar $dname data
    set data($prefix.default)     "data"
    set data($prefix.information) "Graphics output format"
}

proc default_bool {dname prefix} {
    upvar $dname data
    set data($prefix.default)     0
    set data($prefix.information) {boolean}
}

proc default_boolean {dname prefix} {
    upvar $dname data
    set data($prefix.default)     0
    set data($prefix.information) {boolean}
}

proc default_toggle {dname prefix} {
    upvar $dname data
    set data($prefix.default)     0
    set data($prefix.information) {toggle}
}

#-----------------------------------------------------------------------------
# Code generation for each specific ACD type. These have a one-to-one mapping
# with standard widgets or mega-widgets.

proc expand {dname name suffix {bool 0}} {
    upvar $dname data

    set text $data($name.$suffix)    
    regsub -all {%W} $text $name text
    regsub -all {%N} $text $data(appl) text
    set c [regsub -all {%C} $text $data(appl._%C_count) text]
    incr data(appl._%C_count) $c
    set data($name.$suffix) $text
    if {$bool && [string match {\[expr*} $text]} {
        return "\[::EMBOSS::convert_bool \[subst [list $text]\]\]"
    } else {
        return "\[subst [list $text]\]"
    }
}

proc generate_application {dname name args} {
    global cstr menu_file
    upvar $dname data
    set data(appl) $name
    set data(appl._%C_count) 0

    if {[info exists data($name.group)]} {
	set data($name.groups) $data($name.group)
    }
    foreach cat [split $data($name.groups) ,] {
	set cat [string trim $cat { }]
	if {$cat == ""} {
	    continue
	}
	regsub -all {\.} $data($name.documentation) {} desc
	regsub -all {: ?} $cat . cat
	set desc "$desc ($name)"
	if {[info exists menu_file]} {
	    puts $menu_file "lappend [list e_menu($cat)] [list $desc]"
	    puts $menu_file "set [list e_prog($desc)] [list $name]"
	}
    }

    append cstr "option add *e_$name*Xentry.entry.width 30\n"
    append cstr "option add *e_$name*Xcombobox.entry.width 27\n\n"

    append cstr "proc create_dialogue \{\} \{\n"
    append cstr "    variable vars\n"
    append cstr "    variable arguments\n"
    append cstr "    set vars(application) $name\n"
    append cstr "    set w \[xtoplevel .e_$name -resizable 0\]\n"
    append cstr "    if {\$w == {}} return\n"
    append cstr "    bind \$w <Destroy> \\
\"if {{\$w} == {%W}} {::EMBOSS::destroy_dialogue \[namespace current\]}\"\n"
    append cstr "    fix_maxsize \$w\n"
    append cstr "    wm title \$w {EMBOSS - $data(appl)}\n"
    append cstr "    label \$w._title -text [list $data($name.documentation)]\n"
    append cstr "    pack \$w._title -side top -fill both\n"
    append cstr "    ::EMBOSS::init_dialogue \[namespace current\]\n"
    return 1
}

proc generate_end {dname} {
    global cstr
    upvar $dname data
    set name $data(appl)

    append cstr "    okcancelhelp \$w._okcancelhelp -bd 2 -relief groove\\
\t-ok_command \"::EMBOSS::run_dialogue \[namespace current\] \$w\" \\
\t-cancel_command \"destroy \$w\" \\
\t-help_command {show_url \[::EMBOSS::data_dir\]/../doc/programs/html/$name.html}\n"
    append cstr "    pack \$w._okcancelhelp -side bottom -fill x\n"
    if {[info exists data(_tabnotebooks)]} {
        append cstr "\n"
#        append cstr "    wm withdraw \$w\n"
#        foreach book $data(_tabnotebooks) {
#	    append cstr "    ::EMBOSS::resizebook \$book($book).$book\n"
#        }
#        append cstr "    wm deiconify \$w\n"
    }
    append cstr "\}\n"
    return 1
}

proc generate_sequence {dname name args} {
    global cstr
    upvar $dname data
    if {![info exists data(_seqid_count)]} {
	set data(_seqid_count) 0
    } else {
	set data(_seqid_count) [expr {$data(_seqid_count) ? 0 : 1}]
    }
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    set vars($name)       \[get_active_seq_id $data(_seqid_count)\]\n"
    append cstr "    if {\$vars($name) == -1} {set vars($name) \[get_active_seq_id 0\]}\n"
    append cstr "    set vars($name.name)  \[seq_info \$vars($name) name\]\n"
    append cstr "    sequence_changed \[namespace current\] $name\n"
    append cstr "    set vars($name.type) [expand data $name type]\n"
    append cstr "    seq_id \$w.$name \\
\t-textvariable \[namespace current\]::vars($name.name)\\
\t-start_value \$vars($name.begin)\\
\t-end_value \$vars($name.end)\\
\t-to \[seq_info \$vars($name) length\]\\
\t-browse_cmd seq_browser\\
\t-update_cmd \"{::EMBOSS::seq_updates \[namespace current\] $name \$w.$name}\"\n"
    append cstr "    set vars($name) [expand data $name default]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    append cstr "    trace variable vars($name.name) w \
    	\"::EMBOSS::sequence_changed \[namespace current\] $name\"\n"
#    append cstr "    set w \[lindex \$wlist end\]\n"
#    append cstr "    set wlist \[lreplace \$wlist end end\]\n"

    generate_name_trace data $name .information
    generate_value_trace data $name.type
    generate_value_trace data $name .default
    return 1
}

# Treat as generate_sequence for now
proc generate_seqall {dname name args} {
    uplevel generate_sequence [list $dname] [list $name] $args
}

proc generate_value_trace {dname name {default ""}} {
    global cstr
    upvar $dname data

    if {![info exists data($name$default.depend)] || \
        $data($name$default.depend) == ""} {
        return
    }

    regsub -all {%W} $data($name$default.depend) $name data($name$default.depend)

    append cstr "    set vars($name.orig) \$vars($name)\n"

    append cstr "    set vars($name.expr) [list $data($name$default)]\n"
    foreach var $data($name$default.depend) {
        append cstr "    trace variable vars($var) w \
	       \"::EMBOSS::reset_value \[namespace current\] $name\"\n"
    }
    return 1
}

proc generate_name_trace {dname name {default ""}} {
    global cstr
    upvar $dname data

    if {![info exists data($name$default.depend)] || \
        $data($name$default.depend) == ""} {
        return
    }

    append cstr "    set vars($name.info.expr) [list $data($name$default)]\n"
    foreach var $data($name$default.depend) {
        append cstr "    trace variable vars($var) w \
	       \"::EMBOSS::reset_name \[namespace current\] $name\"\n"
    }
}

proc generate_var_trace {dname name} {
    global cstr
    upvar $dname data

    if {![info exists data($name.depend)] || \
        $data($name.depend) == ""} {
        return
    }

    regsub -all {%W} $data($name.depend) $name data($name.depend)

    append cstr "    set vars($name.orig) \$vars($name)\n"

    append cstr "    set vars($name.expr) [list $data($name.default)]\n"
    foreach var $data($name.depend) {
        append cstr "    trace variable vars($var) w \
	       \"::EMBOSS::reset_value \[namespace current\] $name\"\n"
    }
    return 1
}

proc generate_list_trace {dname name procname {default ""}} {
    global cstr
    upvar $dname data

    if {![info exists data($name$default.depend)] || \
        $data($name$default.depend) == ""} {
        return
    }

    regsub -all {%W} $data($name$default.depend) $name data($name$default.depend)

    append cstr "    set vars($name.orig) \$vars($name)\n"

    append cstr "    set vars($name.expr) [list $data($name$default)]\n"
    foreach var $data($name$default.depend) {
        append cstr "    trace variable vars($var) w \
	       \"$procname \[namespace current\] $name\"\n"
    }
    return 1
}

proc generate_integer {dname name args} {
    global cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    xentry \$w.$name \\
\t-type int \\
\t-textvariable \[namespace current\]::vars($name) \\
\t-label [expand data $name information]\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    set vars($name) [expand data $name default]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    if {[info exists data($name.minimum)]} {
	append cstr "    set vars($name.minimum) [expand data $name minimum]\n"
        generate_value_trace data $name.minimum
    }
    if {[info exists data($name.maximum)]} {
	append cstr "    set vars($name.maximum) [expand data $name maximum]\n"
	generate_value_trace data $name.maximum
    }
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"

    generate_name_trace data $name .information
    generate_value_trace data $name .default
    return 1
}

proc generate_int {dname name args} {
    uplevel generate_integer [list $dname] [list $name] $args
    return 1
}

proc generate_float {dname name args} {
    global cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    xentry \$w.$name \\
\t-type float \\
\t-textvariable \[namespace current\]::vars($name) \\
\t-label [expand data $name information]\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    set vars($name) [expand data $name default]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    if {[info exists data($name.minimum)]} {
	append cstr "    set vars($name.minimum) [expand data $name minimum]\n"
	generate_value_trace data $name.minimum
    }
    if {[info exists data($name.maximum)]} {
	append cstr "    set vars($name.maximum) [expand data $name maximum]\n"
	generate_value_trace data $name.maximum
    }
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"

    generate_name_trace data $name .information
    generate_value_trace data $name .default
    return 1
}

proc generate_string {dname name args} {
    global cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    xentry \$w.$name \\
\t-textvariable \[namespace current\]::vars($name) \\
\t-label [expand data $name information]\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    set vars($name) [expand data $name default]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    generate_value_trace data $name .default
    generate_name_trace data $name .information
    return 1
}

# Treat as generate_string for now
proc generate_array {dname name args} {
    uplevel generate_string [list $dname] [list $name] $args
    return 1
}

# Treat as generate_string for now
proc generate_regexp {dname name args} {
    uplevel generate_string [list $dname] [list $name] $args
    return 1
}

# Treat as generate_regexp for now
proc generate_pattern {dname name args} {
    uplevel generate_regexp [list $dname] [list $name] $args
    return 1
}

# Treat as generate_string for now
proc generate_range {dname name args} {
    uplevel generate_string [list $dname] [list $name] $args
    return 1
}

# Treat as generate_string for now
proc generate_features {dname name args} {
    uplevel generate_string [list $dname] [list $name] $args
    return 1
}

proc list_expand {dname name mname1 mname2} {
    upvar $dname data
    upvar $mname1 mapping1
    upvar $mname2 mapping2

    set list ""
    foreach item [split $data($name.values) $data($name.delimiter)] {
        set item [string trimleft $item]
	set delim [string first $data($name.codedelimiter) $item]
	set value [string range $item 0 [expr {$delim-1}]]
        set value [string trimright $value]
	set label [string range $item [expr {$delim+1}] end]
        set label [string trimleft $label]
        set mapping1($label) $value
        set mapping2($value) $label
	lappend list [list $label]
    }
    return $list
}

proc list_expand2 {dname name mname1 mname2} {
    upvar $dname data
    upvar $mname1 mapping1
    upvar $mname2 mapping2

    set list ""
    foreach item [split $data($name.values) $data($name.delimiter)] {
        set item [string trimleft $item]
	set delim [string first $data($name.codedelimiter) $item]
	set value [string range $item 0 [expr {$delim-1}]]
        set value [string trimright $value]
	set label [string range $item [expr {$delim+1}] end]
        set label [string trimleft $label]
        set mapping1($label) $value
        set mapping2($value) $label
	lappend list $label
    }
    return $list
}

proc generate_list_multi {dname name args} {
    global cstr
    upvar $dname data

    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    xscrolledlistbox \$w.$name \\
\t-exportselection 0\\
\t-label [expand data $name header] \\
\t-selectmode extended\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    bind \$w.$name <<ListboxSelect>> \"::EMBOSS::listbox_selected \[namespace current\] $name\"\n"
    append cstr "    pack \$w.$name -side top -fill both -expand 1\n"
    set l [list_expand data $name mapping1 mapping2]
    append cstr "    set vars($name.mapping1) [list [array get mapping1]]\n"
    append cstr "    set vars($name.mapping2) [list [array get mapping2]]\n"
    append cstr "    eval \$w.$name insert end $l\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    append cstr "    set vars($name.delimiter) [list $data($name.delimiter)]\n"

    append cstr "    trace variable vars($name) w \
    \"::EMBOSS::list_multi_changed \[namespace current\] $name\"\n"
    append cstr "    set vars($name) [expand data $name default]\n"

    generate_name_trace data $name .information
    append cstr "    set vars($name._type) list_multi\n"
    generate_value_trace data $name .default 
    return 0
}

proc generate_list {dname name args} {
    global cstr
    upvar $dname data

    if {$data($name.maximum) > 1} {
        return [generate_list_multi data $name]
    }

    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    set l [list_expand2 data $name mapping1 mapping2]
    append cstr "    xcombobox \$w.$name\\
\t-textvariable \[namespace current\]::vars($name.name)\\
\t-label [expand data $name information]\\
\t-values [list $l]\n"
    append cstr "    trace variable vars($name.name) w \
	   \"::EMBOSS::list_changed \[namespace current\] $name\"\n"
    append cstr "    set vars($name.mapping1) [list [array get mapping1]]\n"
    append cstr "    set vars($name.mapping2) [list [array get mapping2]]\n"
    append cstr "    array set tmpmap \$vars($name.mapping2)\n"
    append cstr "    set def [expand data $name default]\n"
    append cstr "    catch {set def \$tmpmap(\$def)}\n"
    append cstr "    set vars($name) \$def\n"
    append cstr "    \$w.$name set \$def\n"
    append cstr "    \$w.$name configure \\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    append cstr "    set vars($name.delimiter) [list $data($name.delimiter)]\n"
    generate_name_trace data $name .information
    generate_list_trace data $name ::EMBOSS::reset_list .default 
    return 1
}

proc selection_expand {dname name} {
    upvar $dname data

    set list ""
    foreach item [split $data($name.values) $data($name.delimiter)] {
        set item [string trimleft $item]
        lappend list [list $item]
    }
    return $list
}

proc generate_selection_multi {dname name args} {
    global cstr
    upvar $dname data

    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    xscrolledlistbox \$w.$name \\
\t-exportselection 0\\
\t-label [expand data $name header] \\
\t-selectmode extended\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    bind \$w.$name <<ListboxSelect>> \"::EMBOSS::selection_selected \[namespace current\] $name\"\n"
    append cstr "    pack \$w.$name -side top -fill both -expand 1\n"
    set l [selection_expand data $name]
    append cstr "    eval \$w.$name insert end $l\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"

    append cstr "    trace variable vars($name) w \
    \"::EMBOSS::select_multi_changed \[namespace current\] $name\"\n"
    append cstr "    set vars($name) [expand data $name default]\n"

    generate_name_trace data $name .information
    append cstr "    set vars($name._type) selection_multi\n"
    generate_value_trace data $name .default 
    return 0
}

proc generate_selection {dname name args} {
    global cstr
    upvar $dname data

    if {$data($name.maximum) > 1} {
        return [generate_selection_multi data $name]
    }

    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    set l [selection_expand data $name]
    append cstr "    xcombobox \$w.$name\\
\t-textvariable \[namespace current\]::vars($name)\\
\t-label [expand data $name information]\\
\t-values [list $l]\n"
    append cstr "    \$w.$name set [expand data $name default]\n"
    append cstr "    \$w.$name configure \\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    generate_name_trace data $name .information
    generate_list_trace data $name ::EMBOSS::reset_select .default 
    return 1
}

proc generate_xygraph {dname name args} {
    global cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    xcombobox \$w.$name\\
\t-textvariable \[namespace current\]::vars($name)\\
\t-label [expand data $name information]\\
\t-values \[list_graph_types\]\n"
    append cstr "    \$w.$name set [expand data $name default]\n"
    append cstr "    \$w.$name configure \\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    generate_name_trace data $name .information
    generate_value_trace data $name .default
    return 1
}

proc generate_graph {dname name args} {
    global cstr
    upvar $dname data

    # FIXME: shouldn't be needed really - fix EMBOSS
    # Specific check for bool: data [ info: "Display as data" ] query. If this
    # exists then set the default here to be xwindows instead of data, as it
    # probably means that this program cannot understand the data graph type.
    if {[info exists data(data)] && \
        $data(data) == "bool" && \
        $data($name.default) == "data"} {
        set data($name.default) xwindows
    }

    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    xcombobox \$w.$name\\
\t-textvariable \[namespace current\]::vars($name)\\
\t-label [expand data $name information]\\
\t-values \[list_graph_types\]\n"
    append cstr "    \$w.$name set [expand data $name default]\n"
    append cstr "    \$w.$name configure \\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    generate_name_trace data $name .information
    generate_value_trace data $name .default
    return 1
}

proc generate_outfile {dname name args} {
    global cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    xentry \$w.$name \\
\t-label [expand data $name information]\\
\t-textvariable \[namespace current\]::vars($name) \\
\t-checkcommand ::EMBOSS::check_outfile\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    \$w.$name delete 0 end\n"
    append cstr "    \$w.$name insert end [expand data $name default]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    generate_name_trace data $name .information
    generate_value_trace data $name .default
    return 1
}

# Treat as generate_outfile for now
proc generate_outdir {dname name args} {
    uplevel generate_outfile [list $dname] [list $name] $args
    return 1
}

# Treat as generate_outfile for now
proc generate_featout {dname name args} {
    uplevel generate_outfile [list $dname] [list $name] $args
    return 1
}

# Treat as generate_outfile for now
proc generate_report {dname name args} {
    uplevel generate_outfile [list $dname] [list $name] $args
    return 1
}

# Treat as generate_outfile for now
proc generate_align {dname name args} {
    uplevel generate_outfile [list $dname] [list $name] $args
    return 1
}

# Treat as generate_outfile for now
proc generate_outcodon {dname name args} {
    uplevel generate_outfile [list $dname] [list $name] $args
    return 1
}

proc generate_dirlist {dname name args} {
    global cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    xentry \$w.$name \\
\t-label [expand data $name information]\\
\t-textvariable \[namespace current\]::vars($name) \\
\t-checkcommand ::EMBOSS::check_directory\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    \$w.$name delete 0 end\n"
    append cstr "    \$w.$name insert end [expand data $name default]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    generate_name_trace data $name .information
    generate_value_trace data $name .default
    return 1
}

proc generate_seqoutall {dname name args} {
    global cstr
    upvar $dname data

    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    labelframe \$w.$name \\
\t-text [expand data $name information]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    append cstr "    lappend wlist \$w\n"
    append cstr "    set w \$w.$name\n"
    append cstr "    xcombobox \$w.format\\
\t-textvariable \[namespace current\]::vars($name.format)\\
\t-label {File format}\\
\t-values \[list_file_formats\]\n"
    append cstr "    \$w.format set fasta\n"
    append cstr "    \$w.format configure \\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    set vars($name.format.path) \$w.format\n"
    append cstr "    pack \$w.format -side top -fill both\n"
    append cstr "    xentry \$w.name \\
\t-label {Filename}\\
\t-textvariable \[namespace current\]::vars($name) \\
\t-checkcommand ::EMBOSS::check_outfile\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    set vars($name.name.path) \$w.name\n"
    append cstr "    pack \$w.name -side top -fill both\n"
    append cstr "    \$w.name delete 0 end\n"
    append cstr "    \$w.name insert end [expand data $name default]\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    append cstr "    set w \[lindex \$wlist end\]\n"
    append cstr "    set wlist \[lreplace \$wlist end end\]\n"
#    generate_name_trace data $name .information
    generate_value_trace data $name .default

    # FIX for now as seqoutall doesn't have a proper iwidget, so we have to
    # manually reconfigure the children.
    if {[info exists data($name.needed.depend)]} {
        foreach var $data($name.needed.depend) {
	    append cstr "    set vars($name.format.needed_expr) \
	      [list "\[::EMBOSS::convert_bool $data($name.needed)\]"]\n"
	    append cstr "    trace variable vars($var) w \
	      \"::EMBOSS::reset_needed \[namespace current\] $name.format\"\n"
	    append cstr "    set vars($name.name.needed_expr) \
	      [list "\[::EMBOSS::convert_bool $data($name.needed)\]"]\n"
	    append cstr "    trace variable vars($var) w \
	      \"::EMBOSS::reset_needed \[namespace current\] $name.name\"\n"
        }
    }

    return 1
}

# Treat as generate_seqoutall for now
proc generate_seqoutset {dname name args} {
    uplevel generate_seqoutall [list $dname] [list $name] $args
    return 1
}

# Treat as generate_outfile for now
proc generate_seqout {dname name args} {
    uplevel generate_seqoutall [list $dname] [list $name] $args
    return 1
}

# treat as generate_sequence for now
proc generate_seqset {dname name args} {
    global cstr
    upvar $dname data

    uplevel generate_sequence [list $dname] [list $name] $args
    append cstr "   set vars($name.totweight) [expand data $name totweight]\n"
    return 1
}

# treat as generate_infile for now
proc generate_filelist {dname name args} {
    global cstr
    upvar $dname data

    uplevel generate_infile [list $dname] [list $name] $args
    return 1
}

proc generate_directory {dname name args} {
    global cstr
    upvar $dname data

    uplevel generate_infile [list $dname] [list $name] $args
    return 1
}

proc generate_infile {dname name args} {
    global cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    xentry \$w.$name \\
\t-label [expand data $name information]\\
\t-textvariable \[namespace current\]::vars($name) \\
\t-checkcommand ::EMBOSS::check_infile\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
#\t-checkcommand ::EMBOSS::check_infile
    append cstr "    \$w.$name delete 0 end\n"
    append cstr "    \$w.$name insert end [expand data $name default]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    generate_name_trace data $name .information
    generate_value_trace data $name .default
    return 1
}

# Treat as generate_infile for now
proc generate_datafile {dname name args} {
    uplevel generate_infile [list $dname] [list $name] $args
    return 1
}

proc generate_matrix {dname name args} {
    global cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    set vars($name.protein) [expand data $name protein]\n"
    append cstr "    xcombobox \$w.$name\\
\t-textvariable \[namespace current\]::vars($name)\\
\t-label [expand data $name information]\\
\t-values \[list_matrices p\]\n"
    append cstr "    \$w.$name set [expand data $name default]\n"
    append cstr "    \$w.$name configure \\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n" 

    generate_value_trace data $name.protein
    generate_value_trace data $name .default
    generate_name_trace data $name .information
    return 1
}

# Treat as generate_matrix
proc generate_matrixf {dname name args} {
    uplevel generate_matrix [list $dname] [list $name] $args
    return 1
}

proc generate_codon {dname name args} {
    global cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    xcombobox \$w.$name\\
\t-textvariable \[namespace current\]::vars($name)\\
\t-label [expand data $name information]\\
\t-values \[list_codon_tables\]\n"
    append cstr "    \$w.$name set [expand data $name default]\n"
    append cstr "    \$w.$name configure \\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    generate_name_trace data $name .information
    generate_value_trace data $name .default
    return 1
}

proc generate_boolean {dname name args} {
    global cstr
    upvar $dname data

    # FIXME: override defaults. If this is a "Display as data" style question
    # Then set the default to be true.
    if {[string match "*isplay as data" $data($name.information)] && \
        $data($name.default) == 0} {
	set data($name.default) 1
    }

    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    checkbutton \$w.$name \\
\t-text [expand data $name information]\\
\t-variable \[namespace current\]::vars($name)\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    pack \$w.$name -side top -anchor w\n"
    append cstr "    set vars($name) [expand data $name default]\n"
    append cstr "    set vars($name.path) \$w.$name\n"
    append cstr "    set vars($name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    generate_name_trace data $name .information
    generate_value_trace data $name .default
    return 1
}

proc generate_bool {dname name args} {
    uplevel generate_boolean [list $dname] [list $name] $args
    return 1
}

proc generate_toggle {dname name args} {
    uplevel generate_boolean [list $dname] [list $name] $args
    return 1
}

proc generate_var {dname name args} {
    global cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    set vars($name) [expand data $name default]\n"
    generate_var_trace data $name
    return 1
}

proc generate_section {dname name args} {
    global page_depth
    global cstr
    upvar $dname data
    append cstr "\n"
    incr page_depth
    switch -- $data($name.type) {
	page {
	    set book $data($name.book)
	    if {$book == ""} {
	        set book book_$page_depth
	    }
	    if {![info exists data($book.created)]} {
		append cstr "    ttk::notebook \$w.$book\n"
		append cstr "    pack \$w.$book -side top -fill both\n"
		append cstr "    set book($book) \$w\n"
		set data($book.created) 1
		lappend data(_tabnotebooks) $book
	    } else {
		incr data($book.created)
	    }

	    append cstr "    set page \$w.$book.f_$data($book.created)\n"
	    append cstr "    frame \$page\n"
	    append cstr "    \$w.$book add \$page -text [expand data $name information]\n"

	    if {$data($book.created) == 1} {
		append cstr "    $\w.$book select 0\n"
	    }
	    append cstr "    lappend wlist \$w\n"
	    append cstr "    set w \$page\n"
	}
        frame {
            if {$data($name.border)} {
	        append cstr "    labelframe \$w.$name \\
\t-text [expand data $name information]\n"
	        append cstr "    pack \$w.$name -side $data($name.side) -fill both\n"
	        append cstr "    lappend wlist \$w\n"
	        append cstr "    append w .$name\n"
            } else {
	        append cstr "    frame \$w.$name -bd 0\n"
	        append cstr "    pack \$w.$name -side $data($name.side) -fill both\n"
	        append cstr "    lappend wlist \$w\n"
	        append cstr "    append w .$name\n"
            }
	}
	default {
	    puts stderr "Unknown section type '$data($name.type)'"
        }
    }
    return 0
}

proc generate_endsection {fname name args} {
    global cstr page_depth
    incr page_depth -1
    append cstr "\n"
    append cstr "    set w \[lindex \$wlist end\]\n"
    append cstr "    set wlist \[lreplace \$wlist end end\]\n"
    return 0
}

# Not needed at present
# Computes the maximum text length of prompts.
proc generate_pass1 {aname {prefix {}} {max_len 0}} {
    upvar $aname data
    set order $data(ORDER:$prefix)
    foreach i $order {
	set wargs [array names data $i.*]
	regsub -all "$i\." $wargs "" wargs
	if {[lsearch {string int matrix} $data($i)] != -1} {
	    set len [font measure button_font $data($i.information)]
	    if {$len > $max_len} {
		set max_len $len
	    }
	}
	if {[info exists data(ORDER:$prefix$i)]} {
	    set len [generate_pass1 data $prefix$i $max_len]
	    if {$len > $max_len} {
		set max_len $len
	    }
	}
    }
    return $max_len
}

proc generate_pass2 {aname {prefix {}} {level 0}} {
    global cstr page_depth
    upvar $aname data
    set order $data(ORDER:$prefix)
    set page_depth 0
    foreach i $order {
	set wargs [array names data $i.*]
	regsub -all "$i\." $wargs "" wargs
	set type [string tolower $data($i)]
        if {[eval generate_$type data $i $wargs] == 1} {
	    append cstr "    set vars($i._type) $data($i)\n"
	}
	# "needed" dependency on another variable
	if {[info exists data($i.needed.depend)]} {
	    foreach var $data($i.needed.depend) {
		append cstr "    set vars($i.needed_expr) \
			[list "\[::EMBOSS::convert_bool $data($i.needed)\]"]\n"
		append cstr "    trace variable vars($var) w \
			\"::EMBOSS::reset_needed \[namespace current\] $i\"\n"
	    }
	}
	if {[info exists data(ORDER:$prefix$i)]} {
	    generate_pass2 data $prefix$i [expr $level+1]
	}
    }
    if {$level == 0} {
	generate_end data
    }
}

# For on-the-fly usage: generate code for an EMBOSS program and return it
proc acd2tcl {prog} {
    variable rules
    variable cstr
    global menu_file

    # Just incase something else sets it.
    # generate_application uses this (when set).
    catch {unset menu_file}

    set fn [file join [::EMBOSS::acd_dir] $prog.acd]
    set fd [open $fn]
    set acd [read $fd]
    close $fd

    set toks [tokenise_string $acd $rules]
    #print_tokens $toks

    set cstr ""
    parse_acd $toks data
    #parray data

    set name [lindex $data(ORDER:) 0]

    append startstr "namespace eval ::EMBOSS::$name \{\n"
    append startstr "namespace import ::EMBOSS::*\n"
    append startstr "variable vars\n"
    append startstr "variable arguments\n"
    append startstr "\n"
    generate_pass2 data
    append endstr "\n"
    append endstr "\}; # namespace eval\n"

    return "$startstr $cstr $endstr"
}

#-----------------------------------------------------------------------------

# When being run on the command line as a program
if {[info level] == 0} {
    set f [open $argv]
    set acd [read $f]
    close $f

    set toks [tokenise_string $acd $rules]
    #print_tokens $toks

    parse_acd $toks data
    #parray data

    #measure maximum size of label
    #set max_len [generate_pass1 data]
    #eval font create button_font -family Helvetica -weight bold -size -12
    #set max_len [expr {1+int(ceil($max_len/[font measure button_font 0]))}]

    #set menu_file [open menu_file a]
    set menu_file stderr
    set name [lindex $data(ORDER:) 0]

    #set startstr "package require Iwidgets\n\n"
    append startstr "namespace eval ::EMBOSS::$name \{\n"
    append startstr "namespace import ::EMBOSS::*\n"
    append startstr "variable vars\n"
    append startstr "variable arguments\n"
    append startstr "\n"
    generate_pass2 data
    append endstr "\n"
    append endstr "\}; # namespace eval\n"

    set fp stdout
    puts $fp $startstr
    puts $fp $cstr
    puts $fp $endstr
    exit
}