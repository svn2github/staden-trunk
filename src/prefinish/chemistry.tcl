#
# Sequencing chemistry functions.
# All functions and data is held in the ::chemistry namespace.
#
# Public functions provided:
#
# num_to_str {num}
#	Converts a gap4 chemistry code to a user-readable string
#
# str_to_num {num}
#	Converts a user-readable string to a gap4 chemistry code
#
# names {}
#	Returns a Tcl list of user-readable chemistry names
#


namespace eval ::chemistry {
    variable chemistries "\
        {Unknown primer} 0
	{Unknown terminator} 1
	{ABI rhodamine primer} 2
        {ABI rhodamine terminator} 3
	{ABI dRhodamine primer} 4
	{ABI dRhodamine terminator} 5
	{Energy transfer primer} 8
	{Energy transfer terminator} 9
	{ABI BigDye primer V1} 14
	{ABI BigDye terminator V1} 15
	{ABI BigDye primer V2} 6
	{ABI BigDye terminator V2} 7
	{ABI BigDye primer V3} 16
	{ABI BigDye terminator V3} 17
	{MegaBACE ET primer} 12
	{MegaBACE ET terminator} 13
	{LiCor primer} 10
	{LiCor terminator} 11
    "
}

proc ::chemistry::initialise {} {
    variable chemistries
    variable chnames
    variable chemistry2num
    variable num2chemistry

    set chnames {}
    foreach {name value} $chemistries {
	lappend chnames $name
	set chemistry2num($name) $value
	set num2chemistry($value) $name
    }
}

proc ::chemistry::num_to_str {num} {
    variable num2chemistry
    return $num2chemistry($num)
}

proc ::chemistry::str_to_num {str} {
    variable chemistry2num
    return $chemistry2num($str)
}

proc ::chemistry::names {} {
    variable chnames
    return $chnames
}

::chemistry::initialise
