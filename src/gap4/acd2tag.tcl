#
# Parser for ACD -> Tcl tag dialogues
# A ripoff from the Spin EMBOSS interface, but for simplicity of maintenance
# (given that Spin is no longer worked on) I've split this off from that
# source.
#

#-----------------------------------------------------------------------------
# The lexical analyser

namespace eval acd2tag {

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

variable tlist {
    {^.(.*).$} {\1}
    {\\[ \n\r]+} {\\n}
    {[ \n\r\t]+} { }
    {\\n} "\n"
    {\\(.)} {\1}
}

variable rules [format {
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

variable page_depth 0

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
	regsub {\$\(([^)]*)\)} $code "\$vars($nul\${namespace}::\\1$nul)" code
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
#    regsub -all {([^:][=?:<> ] *)([a-zA-Z_][a-zA-Z_0-9]*)} $code {\1"\2"} code
    regsub -all {([=?:<> ] *)([a-zA-Z_][a-zA-Z_0-9]*)} $code {\1"\2"} code
    # Unfortunately a::foo?x:y becomes a::"foo"?"x":"y", so remove quotes
    # again after ::
    regsub -all {::"([^"]*)"} $code {::\1} code

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
		foreach word {acdtag information default required \
			      parameter needed delimiter codedelimiter \
			      values selection minimum maximum} {
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
    set data($prefix.button)   0
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

proc default_string {dname prefix} {
    upvar $dname data
    set data($prefix.information) $prefix
    set data($prefix.default)     ""
}

proc default_label {dname prefix} {
    upvar $dname data
    set data($prefix.side)	  "top"
    set data($prefix.information) $prefix
    set data($prefix.default)     ""
    set data($prefix.width)       30
    set data($prefix.labelpos)    "w"
}

proc default_text {dname prefix} {
    upvar $dname data
    set data($prefix.information) ""
    set data($prefix.labelpos)    "n"
    set data($prefix.default)     ""
    set data($prefix.geometry)    40x5
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

proc default_radio {dname prefix} {
    upvar $dname data
    set data($prefix.side)	  "top"
    set data($prefix.default)     {}
    set data($prefix.information) {boolean}
    set data($prefix.store)	  radio
    set data($prefix.content)     {%W}
}

#-----------------------------------------------------------------------------
# Trace triggers called to keep various things up to date.

proc generate_value_trace {dname name {default ""}} {
    variable cstr
    upvar $dname data

    if {![info exists data($name$default.depend)] || \
        $data($name$default.depend) == ""} {
        return
    }

    regsub -all {%W} $data($name$default.depend) $name data($name$default.depend)

    append cstr "    set vars(\${namespace}::$name.orig) \$vars(\${namespace}::$name)\n"

    append cstr "    set vars(\${namespace}::$name.expr) \"[quote_str $data($name$default)]\"\n"
    foreach var [lsort -unique $data($name$default.depend)] {
        append cstr "    trace variable vars($var) w \
	       \"::acd2tag::reset_value \$varsp \${namespace}::$name\"\n"
    }
    return 1
}

proc generate_name_trace {dname name {default ""}} {
    variable cstr
    upvar $dname data

    if {![info exists data($name$default.depend)] || \
        $data($name$default.depend) == ""} {
        return
    }

    append cstr "    set vars(\${namespace}::$name.info.expr) \"[quote_str $data($name$default)]\"\n"
    foreach var $data($name$default.depend) {
        append cstr "    trace variable vars($var) w \
	       \"::acd2tag::reset_name \$varsp \${namespace}::$name\"\n"
    }
}

proc generate_var_trace {dname name} {
    variable cstr
    upvar $dname data

    if {![info exists data($name.depend)] || \
        $data($name.depend) == ""} {
        return
    }

    regsub -all {%W} $data($name.depend) $name data($name.depend)

    append cstr "    set vars(\${namespace}::$name.orig) \$vars(\${namespace}::$name)\n"

    append cstr "    set vars(\${namespace}::$name.expr) \"[quote_str $data($name.default)]\"\n"
    foreach var $data($name.depend) {
        append cstr "    trace variable vars($var) w \
	       \"::acd2tag::reset_value \$varsp $name\"\n"
    }
    return 1
}

proc generate_list_trace {dname name procname {default ""}} {
    variable cstr
    upvar $dname data

    if {![info exists data($name$default.depend)] || \
        $data($name$default.depend) == ""} {
        return
    }

    regsub -all {%W} $data($name$default.depend) $name data($name$default.depend)

    append cstr "    set vars(\${namespace}::$name.orig) \$vars(\${namespace}::$name)\n"

    append cstr "    set vars(\${namespace}::$name.expr) \"[quote_str $data($name$default)]\"\n"
    foreach var $data($name$default.depend) {
        append cstr "    trace variable vars($var) w \
	       \"$procname \$varsp \${namespace}::$name\"\n"
    }
    return 1
}


# Expands up any inline ACD expressions and % rules within a string
proc expand {dname name suffix {bool 0}} {
    upvar $dname data

    set text $data($name.$suffix)    
    regsub -all {%W} $text $name text
    regsub -all {%N} $text $data(tag_type) text
    set c [regsub -all {%C} $text $data(tag_type._%C_count) text]
    incr data(tag_type._%C_count) $c
    set data($name.$suffix) $text
    if {$bool && [string match {\[expr*} $text]} {
        return "\[::acd2tag::convert_bool \[subst [list $text]\]\]"
    } else {
        return "\[subst [list $text]\]"
    }
}

#-----------------------------------------------------------------------------
# Code generation for each specific ACD type. These have a one-to-one mapping
# with standard widgets or mega-widgets.

proc generate_acdtag {dname name args} {
    variable cstr
    upvar $dname data
    set data(tag_type) $name
    set data(tag_type._%C_count) 0

    append cstr "option add *e_$name*Xentry.entry.width 30\n"
    append cstr "option add *e_$name*Xcombobox.entry.width 27\n\n"

    append cstr "\n"
    append cstr "variable last_selection 0\n"
    append cstr "proc page_raise {book} {\n"
    append cstr "    variable last_selection\n"
    append cstr "    set last_selection \[\$book view\]\n"
    append cstr "}\n"
    append cstr "\n"

    append cstr "proc create_dialogue \{w varsp\} \{\n"
    append cstr "    upvar \$varsp vars\n"
    append cstr "    variable arguments\n"
    append cstr "    set vars(acdtag) $name\n"
    append cstr "    label \$w._title -text [list $data($name.information)]\n"
    append cstr "    pack \$w._title -side top -fill both\n"
    append cstr "    set namespace \"\"\n"
    append cstr "\n"
    append cstr "    if {!\[info exists vars(namespaces)\]} {\n"
    append cstr "        set vars(namespaces) {}\n"
    append cstr "    }\n"
    append cstr "    set ns \$vars(namespaces)\n"
    append cstr "    if {\$ns == {}} {\n"
    append cstr "	set ns \[list {}\]\n"
    append cstr "    } else {\n"
    append cstr "	ttk::notebook \$w.ns_book\n"
    append cstr "	pack \$w.ns_book -side top -fill both -expand 1\n"
    append cstr "	set book(ns_book) \$w\n"
    append cstr "    }\n"
    append cstr "\n"
    append cstr "    set tabnum 0\n"
    append cstr "    foreach namespace \$ns \{\n"
    append cstr "	if {\$namespace != {}} {\n"
    append cstr "	    lappend wlist \$w\n"
    append cstr "	    if {\[info exists vars(\${namespace}::TABNAME)\]} {\n"
    append cstr "		set tabname \$vars(\${namespace}::TABNAME)\n"
    append cstr "	    } else {\n"
    append cstr "		set tabname \$namespace\n"
    append cstr "	    }\n"
    append cstr "	    set page \$w.ns_book.p_\$tabnum\n"
    append cstr "	    incr tabnum\n"
    append cstr "	    frame \$page\n"
    append cstr "	    \$w.ns_book add \$page -text \$tabname -command \"\[namespace current\]::page_raise \$w.ns_book\"\n"
    append cstr "	    lappend wlist \$w\n"
    append cstr "	    set w \$page\n"
    append cstr "	}\n"
    append cstr "\n"

    return 1
}

proc generate_end {dname} {
    variable cstr
    append cstr "	if {\$namespace != {}} {\n"
    append cstr "	    set w \[lindex \$wlist end\]\n"
    append cstr "	    set wlist \[lreplace \$wlist end end\]\n"
    append cstr "	    set w \[lindex \$wlist end\]\n"
    append cstr "	    set wlist \[lreplace \$wlist end end\]\n"
    append cstr "	}\n"
    append cstr "    \}\n"
    append cstr "\n"
    append cstr "    if {\$vars(namespaces) != {}} {\n"
    append cstr "	variable last_selection\n"
    append cstr "\n"
    append cstr "	if {\$last_selection >= 0} {\n"
    append cstr "	    if {\[catch {\$w.ns_book select \$last_selection}\]} {\n"
    append cstr "               \$w.ns_book select 0\n"
    append cstr "           }\n"
    append cstr "	} else {\n"
    append cstr "	    \$w.ns_book select 0\n"
    append cstr "	}\n"
#    append cstr "	acd2tag::resizebook \$w.ns_book\n"
#    append cstr "    } else {\n"
#    append cstr "	if {\[info exists book\]} {\n"
#    append cstr "	    foreach b \[array names book\] {\n"
#    append cstr "	        acd2tag::resizebook \$w.\$b\n"
#    append cstr "           }\n"
#    append cstr "        }\n"
    append cstr "    }\n"
    append cstr "\}\n"
    append cstr "\n"
    append cstr "proc vars {} {\n"
    append cstr "    variable arguments\n"
    append cstr "    return \$arguments\n"
    append cstr "}\n"

    return 1
}

proc generate_integer {dname name args} {
    variable cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    if {!\[info exists vars(\${namespace}::$name)\]} {
        set vars(\${namespace}::$name) [expand data $name default]
    }\n"
    append cstr "    xentry \$w.$name \\
\t-type int \\
\t-textvariable \${varsp}(\${namespace}::$name) \\
\t-label [expand data $name information]\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    if {[info exists data($name.minimum)]} {
	append cstr "    set vars(\${namespace}::$name.minimum) [expand data $name minimum]\n"
        generate_value_trace data $name.minimum
    }
    if {[info exists data($name.maximum)]} {
	append cstr "    set vars(\${namespace}::$name.maximum) [expand data $name maximum]\n"
	generate_value_trace data $name.maximum
    }
    append cstr "    set vars(\${namespace}::$name.path) \$w.$name\n"
    append cstr "    set vars(\${namespace}::$name.required) \
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
    variable cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    if {!\[info exists vars(\${namespace}::$name)\]} {
        set vars(\${namespace}::$name) [expand data $name default]
    }\n"
    append cstr "    xentry \$w.$name \\
\t-type float \\
\t-textvariable \${varsp}(\${namespace}::$name) \\
\t-label [expand data $name information]\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    if {[info exists data($name.minimum)]} {
	append cstr "    set vars(\${namespace}::$name.minimum) [expand data $name minimum]\n"
	generate_value_trace data $name.minimum
    }
    if {[info exists data($name.maximum)]} {
	append cstr "    set vars(\${namespace}::$name.maximum) [expand data $name maximum]\n"
	generate_value_trace data $name.maximum
    }
    append cstr "    set vars(\${namespace}::$name.path) \$w.$name\n"
    append cstr "    set vars(\${namespace}::$name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"

    generate_name_trace data $name .information
    generate_value_trace data $name .default
    return 1
}

proc generate_string {dname name args} {
    variable cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    if {!\[info exists vars(\${namespace}::$name)\]} {
        set vars(\${namespace}::$name) [expand data $name default]
    }\n"
    append cstr "    xentry \$w.$name \\
\t-textvariable \${varsp}(\${namespace}::$name) \\
\t-label [expand data $name information]\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    append cstr "    set vars(\${namespace}::$name.path) \$w.$name\n"
    append cstr "    set vars(\${namespace}::$name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    generate_value_trace data $name .default
    generate_name_trace data $name .information
    return 1
}

proc generate_label {dname name args} {
    variable cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    if {!\[info exists vars(\${namespace}::$name)\]} {
        set vars(\${namespace}::$name) [expand data $name default]
    }\n"
    append cstr "    frame \$w.$name -bd 0 -class LabelPair\n"
    append cstr "    label \$w.$name.left\\
\t-text [expand data $name information]\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    label \$w.$name.right \\
\t-anchor [expand data $name labelpos] \\
\t-width [expand data $name width] \\
\t-textvariable \${varsp}(\${namespace}::$name) \\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    pack \$w.$name.left -side left\n"
    append cstr "    pack \$w.$name.right -side right -fill both\n"
    append cstr "    pack \$w.$name -side $data($name.side) -fill both\n"
    append cstr "    set vars(\${namespace}::$name.path) \$w.$name\n"
    generate_name_trace data $name .information
    generate_value_trace data $name .default
    return 1
}

proc generate_text {dname name args} {
    variable cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    if {!\[info exists vars(\${namespace}::$name)\]} {
        set vars(\${namespace}::$name) [expand data $name default]
    }\n"
    if {![regexp {^([0-9]+)x([0-9]+)} $data($name.geometry) _ W H]} {
	set W 80
	set H 24
    }
    append cstr "    xscrolledtext \$w.$name \\
\t-label [expand data $name information] \\
\t-wrap none \\
\t-width $W\\
\t-height $H\n"
    append cstr "    pack \$w.$name -side top -fill both -expand 1\n"
    append cstr "    set vars(\${namespace}::$name.path) \$w.$name\n"
    append cstr "    set vars(\${namespace}::$name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    append cstr "    \$w.$name insert end \$vars(\${namespace}::$name)\n"
    generate_value_trace data $name .default
    append cstr "    trace variable vars(\${namespace}::$name) w \"::acd2tag::reset_text \$varsp  \${namespace}::$name\"\n"
    generate_name_trace data $name .information
    append cstr "    \$w.$name configure -state \[lindex {disabled normal} [expand data $name needed 1]\]\n" 
    return 1
}

# Treat as generate_string for now
proc generate_regexp {dname name args} {
    uplevel generate_string [list $dname] [list $name] $args
    return 1
}

# Treat as generate_string for now
proc generate_range {dname name args} {
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

# FIXME: defaults do not work correctly for this yet
proc generate_list_multi {dname name args} {
    variable cstr
    upvar $dname data

    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    set vars(\${namespace}::$name) [expand data $name default]\n"
    append cstr "    xscrolledlistbox \$w.$name \\
\t-exportselection 0\\
\t-label [expand data $name header] \\
\t-selectmode extended\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    bind \$w.$name <<ListboxSelect>> \"::acd2tag::listbox_selected \$varsp $name\"\n"
    append cstr "    pack \$w.$name -side top -fill both -expand 1\n"
    set l [list_expand data $name mapping1 mapping2]
    append cstr "    set vars(\${namespace}::$name.mapping1) [list [array get mapping1]]\n"
    append cstr "    set vars(\${namespace}::$name.mapping2) [list [array get mapping2]]\n"
    append cstr "    eval \$w.$name insert end $l\n"
    append cstr "    set vars(\${namespace}::$name.path) \$w.$name\n"
    append cstr "    set vars(\${namespace}::$name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    append cstr "    set vars(\${namespace}::$name.delimiter) [list $data($name.delimiter)]\n"

    append cstr "    trace variable vars($name) w \
    \"::acd2tag::list_multi_changed \$varsp $name\"\n"
    append cstr "    set vars(\${namespace}::$name) [expand data $name default]\n"

    generate_name_trace data $name .information
    append cstr "    set vars(\${namespace}::$name._type) list_multi\n"
    generate_value_trace data $name .default 
    return 0
}

# FIXME: defaults do not work correctly for this yet
proc generate_list {dname name args} {
    variable cstr
    upvar $dname data

    if {$data($name.maximum) > 1} {
        return [generate_list_multi data $name]
    }

    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    set l [list_expand2 data $name mapping1 mapping2]
    append cstr "    xcombobox \$w.$name\\
\t-textvariable \${varsp}(\${namespace}::$name.name)\\
\t-label [expand data $name information]\\
\t-values [list $l]\n"
    append cstr "    trace variable vars($name.name) w \
	   \"::acd2tag::list_changed \$varsp $name\"\n"
    append cstr "    set vars(\${namespace}::$name.mapping1) [list [array get mapping1]]\n"
    append cstr "    set vars(\${namespace}::$name.mapping2) [list [array get mapping2]]\n"
    append cstr "    array set tmpmap \$vars(\${namespace}::$name.mapping2)\n"

    append cstr "    if {!\[info exists vars(\${namespace}::$name)\]} {\n"
    append cstr "        set vars(\${namespace}::$name) [expand data $name default]\n"
    append cstr "    }\n"

    append cstr "    set def [expand data $name default]\n"
    append cstr "    catch {set def \$tmpmap(\$def)}\n"
    append cstr "    set vars(\${namespace}::$name) \$def\n"
    append cstr "    \$w.$name set \$def\n"
    append cstr "    \$w.$name configure \\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    pack \$w.$name -side top -fill both\n"
    append cstr "    set vars(\${namespace}::$name.path) \$w.$name\n"
    append cstr "    set vars(\${namespace}::$name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    append cstr "    set vars(\${namespace}::$name.delimiter) [list $data($name.delimiter)]\n"
    generate_name_trace data $name .information
    generate_list_trace data $name ::acd2tag::reset_list .default 
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
    variable cstr
    upvar $dname data

    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    xscrolledlistbox \$w.$name \\
\t-exportselection 0\\
\t-label [expand data $name header] \\
\t-selectmode extended\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    bind \$w.$name <<ListboxSelect>> \"::acd2tag::selection_selected \$varsp $name\"\n"
    append cstr "    pack \$w.$name -side top -fill both -expand 1\n"
    set l [selection_expand data $name]
    append cstr "    eval \$w.$name insert end $l\n"
    append cstr "    set vars(\${namespace}::$name.path) \$w.$name\n"
    append cstr "    set vars(\${namespace}::$name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"

    append cstr "    trace variable vars(\${namespace}::$name) w \
    \"::acd2tag::select_multi_changed \$varsp $name\"\n"
    append cstr "    if {!\[info exists vars(\${namespace}::$name)\]} {\n"
    append cstr "        set vars(\${namespace}::$name) [expand data $name default]\n"
    append cstr "    }\n"
    append cstr "    set vars(\${namespace}::$name) \$vars(\${namespace}::$name)\n"

    generate_name_trace data $name .information
    append cstr "    set vars(\${namespace}::$name._type) selection_multi\n"
    generate_value_trace data $name .default 
    return 0
}

proc generate_selection {dname name args} {
    variable cstr
    upvar $dname data

    if {$data($name.maximum) > 1} {
        return [generate_selection_multi data $name]
    }

    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    if {$data($name.button) == 1} {
	set l [selection_expand data $name]
	set count 0
	append cstr "   label \$w.${name}_label -text [expand data $name information]\n"
	append cstr "   pack \$w.${name}_label -side top -fill both\n"
	foreach item $l {
	    append cstr "   radiobutton \$w.${name}_$count\\
\t-variable \${varsp}(\${namespace}::$name)\\
\t-text [list $item] -value [list $item]\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
            append cstr "    pack \$w.${name}_$count -side top -anchor w\n"
            incr count
	}
    } else {
        append cstr "    if {!\[info exists vars(\${namespace}::$name)\]} {
        set vars(\${namespace}::$name) [expand data $name default]
    }\n"
        set l [selection_expand data $name]
	append cstr "    xcombobox \$w.$name\\
\t-textvariable \${varsp}(\${namespace}::$name)\\
\t-label [expand data $name information]\\
\t-values [list $l]\n"
        append cstr "    set vars(\${namespace}::$name) \$vars(\${namespace}::$name)\n"
        append cstr "    \$w.$name configure \\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
        append cstr "    pack \$w.$name -side top -fill both\n"
    }
    append cstr "    set vars(\${namespace}::$name.path) \$w.$name\n"
    append cstr "    set vars(\${namespace}::$name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    generate_name_trace data $name .information
    generate_list_trace data $name ::acd2tag::reset_select .default 
    return 1
}

proc generate_boolean {dname name args} {
    variable cstr
    upvar $dname data

    # FIXME: override defaults. If this is a "Display as data" style question
    # Then set the default to be true.
    if {[string match "*isplay as data" $data($name.information)] && \
        $data($name.default) == 0} {
	set data($name.default) 1
    }

    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    if {!\[info exists vars(\${namespace}::$name)\]} {
        set vars(\${namespace}::$name) [expand data $name default]
    }\n"
    append cstr "    checkbutton \$w.$name \\
\t-text [expand data $name information]\\
\t-variable \${varsp}(\${namespace}::$name)\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    pack \$w.$name -side top -anchor w\n"
    append cstr "    set vars(\${namespace}::$name.path) \$w.$name\n"
    append cstr "    set vars(\${namespace}::$name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    generate_name_trace data $name .information
    generate_value_trace data $name .default
    return 1
}

proc generate_bool {dname name args} {
    uplevel generate_boolean [list $dname] [list $name] $args
    return 1
}

proc generate_radio {dname name args} {
    variable cstr
    upvar $dname data

    # FIXME: override defaults. If this is a "Display as data" style question
    # Then set the default to be true.
    if {[string match "*isplay as data" $data($name.information)] && \
        $data($name.default) == 0} {
	set data($name.default) 1
    }

    if {[catch {set val [expand data $name content]}]} {
	set val $name
    }
    append cstr "\n"
    append cstr "    lappend arguments [expand data $name store]\n"
    append cstr "    if {!\[info exists vars(\${namespace}::[expand data $name store])\] || \$vars(\${namespace}::[expand data $name store]) == {}} {\n"
    append cstr "        if {[expand data $name default] != {}} {\n"
    append cstr "            set vars(\${namespace}::[expand data $name store]) [expand data $name default]\n"
    append cstr "        }\n"
    append cstr "    }\n"
    append cstr "    radiobutton \$w.$name \\
\t-variable \${varsp}(\${namespace}::[expand data $name store])\\
\t-value $val\\
\t-text [expand data $name information]\\
\t-state \[lindex {disabled normal} [expand data $name needed 1]\]\n"
    append cstr "    pack \$w.$name -side $data($name.side) -anchor w\n"
    append cstr "    set vars(\${namespace}::$name.path) \$w.$name\n"
    append cstr "    set vars(\${namespace}::$name.required) \
                         [expr {$data($name.required)=="Y"?1:0}]\n"
    generate_name_trace data $name .information
    generate_value_trace data $name .default
    return 1
}

proc generate_var {dname name args} {
    variable cstr
    upvar $dname data
    append cstr "\n"
    append cstr "    lappend arguments $name\n"
    append cstr "    if {!\[info exists vars(\${namespace}::$name)\]} {
        set vars(\${namespace}::$name) [expand data $name default]
    }\n"
#    append cstr "    set vars(\${namespace}::$name) [expand data $name default]\n"
    generate_var_trace data $name
    return 1
}

proc generate_section {dname name args} {
    variable page_depth
    variable cstr
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
	        append cstr "    labelframe \$w.$name\\
\t-text [expand data $name information]\n"
	        append cstr "    pack \$w.$name -side $data($name.side) -fill both -expand 1\n"
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
    variable cstr
    variable page_depth
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

# Escapes all square brackets, quotes and dollars in a string except for
# ${namespace}
proc quote_str {str} {
    regsub -all {[]\"$[]} $str {\\&} str
    regsub -all {\\\${namespace}} $str {${namespace}} str
    return $str
}

proc generate_pass2 {aname {prefix {}} {level 0}} {
    variable cstr
    variable page_depth
    upvar $aname data
    set order $data(ORDER:$prefix)
    set page_depth 0
    foreach i $order {
	set wargs [array names data $i.*]
	regsub -all "$i\." $wargs "" wargs
	set type [string tolower $data($i)]
        if {[eval generate_$type data $i $wargs] == 1} {
	    append cstr "    set vars(\${namespace}::$i._type) $data($i)\n"
	}
	# "needed" dependency on another variable
	if {[info exists data($i.needed.depend)]} {
	    foreach var $data($i.needed.depend) {
		set qstr [quote_str "\[::acd2tag::convert_bool $data($i.needed)\]"]
		append cstr "    set vars(\${namespace}::$i.needed_expr) \"$qstr\"\n"
		append cstr "    trace variable vars($var) w \
			\"::acd2tag::reset_needed \$varsp \${namespace}::$i\"\n"
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

#-----------------------------------------------------------------------------
# Dialogue callback functions

#
# Calculates the appropriate window size for a notebook, so as to avoid
# needing to scroll the tabs.
#
# Assumes that all packing within the pages is -side top or -side bottom.
# If this is not the case then enclose the children within a single frame.
proc resizebook {book} {
    return; #disabled for now to see if ttk::notebook has issues.

    update idletasks
    set bd [expr {2*[$book cget -borderwidth]}]
    set maxheight 0
    set maxwidth [expr {int(ceil([lindex [[$book component tabset] bbox] 2]))}]

    set ts [$book component tabset]
    set ntabs [expr {[$ts index end]+1}]
    # Calculate the width due to the sloping tab. For N tabs we have N+1 of
    # these width portions.
    # Assume the same tab angle and height is used throughout
    set angle [$ts tabcget 0 -angle]
    set labelHeight [expr {[$ts tabcget 0 -height] + 2*[$ts tabcget 0 -pady]}]
    set angleOffset [expr {int(ceil($labelHeight * tan($angle/57.2957795133)))}]
    set total_wid $angleOffset
    for {set tab 0} {$tab < $ntabs} {incr tab} {
	set font [$ts tabcget $tab -font]
	set padx [$ts tabcget $tab -padx]
	set text [$ts tabcget $tab -label]
	set labelWidth [font measure $font $text]
	# +7? not sure what this is, but we need a bit more per tab...
	incr total_wid [expr {$angleOffset + $labelWidth+7 + 2*$padx}]
    }
    incr total_wid [$ts cget -start]
    set maxwidth $total_wid

    # compute tabheight
    array set metrics [font metrics $font]
    set th [expr {$metrics(-linespace)+2*[$ts cget -pady]+[$ts cget -margin]+9}]

    # As an alternative to using bbox:
    # set ts [$book component tabset]
    # set ntabs [expr {[$ts index end]+1}]
    # set _labelWidth  [expr {[$ts tabcget 0 -width] +2*[$ts tabcget 0 -padx]}]
    # set _labelHeight [expr {[$ts tabcget 0 -height]+2*[$ts tabcget 0 -pady]}]
    # set angle [$ts tabcget 0 -angle]
    # set angleOffset [expr {$_labelHeight * tan($angle/57.2957795133)}]
    # set maxwidth [expr {int(ceil($ntabs * ($_labelWidth + $angleOffset) + $angleOffset + 4))}]

    foreach page [$book childsite] {
        #set height [expr {$bd+[winfo reqheight [$book component tabset]]}]
        set height [expr {$bd+$th}]
        foreach w [winfo children $page] {
	    if {[winfo class $w] == "Tabnotebook"} {
	        resizebook $w
	    }
            set height [expr {$height+[winfo reqheight $w]}]
	    #puts $w:[winfo reqwidth $w]:[winfo class $w]:[winfo height $w]:[winfo reqheight $w]
	    if {$maxwidth < [winfo reqwidth $w]} {
	        set maxwidth [winfo reqwidth $w]
	    }
        }
	if {$maxheight < $height} {
	    set maxheight $height
        }
    }
    $book configure -height $maxheight
    $book configure -width $maxwidth
    update idletasks
}

# -----------------------------------------------------------------------------
# Trace callbacks, for tracing dialogue component modifications.

#
# Modification of a variable that vars($vname) depends on in namespace $nspace.
#
proc reset_value {varsp vname args} {
    upvar $varsp vars
    if {$vars($vname) == $vars($vname.orig)} {
	set vars($vname) [subst $vars($vname.expr)]
	set vars($vname.orig) $vars($vname)
    }
}

# Modification of a variable used by a single-element list
proc reset_list {varsp vname args} {
    upvar $varsp vars
    if {$vars($vname) != $vars($vname.orig)} {
	return
    }

    eval set newval $vars($vname.expr)
    array set a $vars($vname.mapping2)
    if {[info exists a($newval)]} {
	set vars($vname) $newval
	set vars($vname.orig) $vars($vname)
	set vars($vname.name) $a($newval)
    }
    unset a
}

# Modification of a variable used by a selection list
proc reset_select {varsp vname args} {
    upvar $varsp vars
    if {$vars($vname) == $vars($vname.orig)} {
	eval set t $vars($vname.expr)
	$vars($vname.path) delete entry 0 end
	$vars($vname.path) insert entry end $t
	set vars($vname.orig) $vars($vname)
    }
}


# The 'needed' flag has changed
proc reset_needed {varsp vname args} {
    upvar $varsp vars
    # puts exp=$vars($vname.needed_expr)
    if {[subst $vars($vname.needed_expr)] != 0} {
	catch {$vars($vname.path) configure -state normal} err
    } else {
	catch {$vars($vname.path) configure -state disabled} err
    }
}

# Modification of a variable that the 'info' attribute depends on
proc reset_name {varsp vname args} {
    upvar $varsp vars
    set t [subst $vars($vname.info.expr)]
    if {[winfo class $vars($vname.path)] == "LabelPair"} {
	catch {$vars($vname.path).left configure -text $t} err
    } else {
	if {[catch {$vars($vname.path) configure -labeltext $t} err1]} {
	    if {[catch {$vars($vname.path) configure -label $t} err2]} {
		if {[catch {$vars($vname.path) configure -text $t} err3]} {
		    puts "$vars, $vname"
		    puts "Failed to configure $vname.path -labeltext $t"
		    puts "    1. $err1"
		    puts "    2. $err2"
		    puts "    3. $err3"
		}
	    }
	}
    }
}


# A text window section has been changed
proc reset_text {varsp vname args} {
    upvar $varsp vars
    $vars($vname.path) delete 1.0 end
    $vars($vname.path) insert end $vars($vname)
}

# A text window has been modified
proc text_changed {varsp vname args} {
    upvar $varsp vars
    set vars($vname) [$vars($vname.path) get 1.0 end-1char]
}

# A 'list' selection has been modified
proc list_changed {varsp vname args} {
    upvar $varsp vars
    set newsel $vars($vname.name)
    array set a $vars($vname.mapping1)
    if {[info exists a($newsel)]} {
	set vars($vname) $a($newsel)
    }
}

# A multiple list box has been modified
proc list_multi_changed {varsp vname args} {
    upvar $varsp vars
    $vars($vname.path) selection clear 0 end
    set list [$vars($vname.path) get 0 end]
    array set a $vars($vname.mapping2)
    foreach item [split $vars($vname) ", "] {
	if {$item == ""} {
	    continue
	}
	catch {$vars($vname.path) selection set [lsearch $list $a($item)]} err
    }
    unset a
}

# A multiple selection box has been modified
proc select_multi_changed {varsp vname args} {
    upvar $varsp vars
    $vars($vname.path) selection clear 0 end
    set list [$vars($vname.path) get 0 end]
    foreach item [split $vars($vname) ", "] {
	if {$item == ""} {
	    continue
	}
	catch {$vars($vname.path) selection set [lsearch -glob $list $item*]} err
    }
}

# A scrolledlistbox -selectioncommand callback - updates the variable
proc listbox_selected {varsp vname args} {
    upvar $varsp vars
    set sel [$vars($vname.path) curselection]
    array set a $vars($vname.mapping1)
    foreach i $sel {
	catch {lappend l $a([$vars($vname.path) get $i])}
    }
    unset a
    set vars($vname) $l
}

# A scrolledlistbox -selectioncommand callback - updates the variable
proc selection_selected {varsp vname args} {
    upvar $varsp vars
    set sel [$vars($vname.path) curselection]
    foreach i $sel {
	lappend l [$vars($vname.path) get $i]
    }
    set vars($vname) $l
}

# Converts y/n, t/f into 1/0
proc convert_bool {text} {
    if {[string match {[YyTt]*} $text]} {
	return 1
    } elseif {[string match {[NnFf]*} $text]} {
	return 0
    } else {
	return $text
    }
}

# -----------------------------------------------------------------------------
# Functions for handling the OK, Cancel, Help buttons.
#
# These are all running in the ::acd2tag namespace, but are given the namespace
# of the dialogue to operate on.

# Tidy up function; a callback when the window is destroyed.
# This removes any namespace-global variables. This will also have the effect
# of removing the (possibly many) variable traces.
proc destroy_dialogue {nspace varsp} {
    namespace eval $nspace {     
	upvar $varsp vars
	unset vars

	variable arguments
	unset arguments
    }
}

#
# Reads elements from the dialogue and returns
# the associated tag contents string.
#
proc get_values {nspace varsp} {
    upvar $varsp vars
    upvar ${nspace}::arguments arguments

    set val ""
    set ns $vars(namespaces)
    if {$ns == ""} {set ns {""}}
    foreach n $ns {
	if {$n != ""} {
	    append val "NAMESPACE=$n\n"
	}
	catch {unset done}
	array set vals {}
	array set done ""
	if {[info exists vars(${n}::TABNAME)]} {
	    append val "TABNAME=$vars(${n}::TABNAME)\n"
	}
	foreach arg $arguments {
	    # Text windows do not automatically update the value upon edit,
	    # so do so now.
	    if {[info exists vars(${n}::$arg._type)] &&
		$vars(${n}::$arg._type) == "text"} {
		text_changed vars ${n}::$arg
	    }
	    if {[info exists done($arg)]} {
		continue
	    }
	    set vals($arg) $vars(${n}::$arg)
	    set done($arg) 1
	}
	append val [array2str [array get vals]]
    }

    # Strip off any trailing nl
    return [regsub {\n*$} $val {}]
}

#-----------------------------------------------------------------------------
# Converts tag=value syntax into array get/set syntax

# This replaces \n with newline and \\ with \.
proc unescape {str} {
    # not_backslash + (double backslash pairs) "\n" -> stuff_matched + nl
    regsub -all {((^|[^\\])(\\\\)*)\\n} $str "\\1\n" str
    regsub -all {\\\\} $str "\\" str
    return $str
}

# Escapes newline with \n, and also makes sure that previous \ are escaped
# with \\.
proc escape {str} {
    regsub -all {\\} $str {\\\\} str
    regsub -all "\n" $str {\\n} str
    return $str
}

# Replaces a "key=value" list (separated with newlines) with an array returned
# in [array get] syntax.
proc str2array {str} {
    array set a {namespaces ""}
    set namespace ""
    foreach line [split $str "\n"] {
	regexp {([^=]*)=(.*)} $line _ key value
	if {$key == "NAMESPACE"} {
	    set namespace $value
	    lappend a(namespaces) $value
	}
	set a(${namespace}::[::acd2tag::unescape $key]) [::acd2tag::unescape $value]
    }

    return [array get a]
}

# Converts a stringified array (from [array get]) into a tag=value list
# separated with newlines. Newlines in tag or value are escaped appropriately.
proc array2str {astr} {
    set str ""
    foreach {tag value} $astr {
	append str "[escape $tag]=[escape $value]\n"
    }

    return $str
}

#-----------------------------------------------------------------------------

proc parse {acd ns} {
    variable rules
    variable cstr

    set toks [tokenise_string $acd $rules]
    #print_tokens $toks

    set cstr ""
    parse_acd $toks data
    #parray data

    #measure maximum size of label
    #set max_len [generate_pass1 data]
    #eval font create button_font -family Helvetica -weight bold -size -12
    #set max_len [expr {1+int(ceil($max_len/[font measure button_font 0]))}]

    set name [lindex $data(ORDER:) 0]

    set startstr ""
    append startstr "namespace eval $ns \{\n"
    append startstr "variable arguments\n"
    append startstr "\n"
    generate_pass2 data
    append endstr "\n"
    append endstr "\}; # namespace eval\n"
    
    return "$startstr $cstr $endstr"
}

}; # namespace eval
