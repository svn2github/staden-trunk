#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1998. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#----------
# Extracting the sample names from ABI, ALF or SCF files
proc ABI_get_sample_name {name} {
    # Extract the SMPL field from the ABI file - io_lib handles this for us
    set fd [open "|get_comment -c NAME < \"$name\""]
    set id [gets $fd]
    if {[catch {close $fd}]} {
	return ""
    }

    regsub -all / $id _ id
    return $id
}

proc ALF_get_sample_name {name} {
    # Currently no known method
    return ""
}

#31/05/00 johnt - added biolims support
# is this function required ?
proc BIO_get_sample_name {name} {
    # the leaf of the path is the lane name, which
    # seems to be the same as the sample name
    set fd [string range $name [expr [string last / $name]+1] end]

    return $id
}

proc SCF_get_sample_name {name} {
    # Extract the NAME field from the SCF file
    set fd [open "|get_comment -c NAME < \"$name\""]
    set id [gets $fd]
    if {[catch {close $fd}]} {
	return ""
    }

    # For old LiCor sequences, which had NAME as a full pathname:
    regsub -all {.*\\} $id {} id

    # For Mac sample names with slashes in them (eg in dates)
    regsub -all / $id _ id
    return $id
}

proc ZTR_get_sample_name {name} {
    set fd [open "|get_comment -c NAME < \"$name\""]
    set id [gets $fd]
    if {[catch {close $fd}]} {
	return ""
    }

    regsub -all / $id _ id
    return $id
}

proc CTF_get_sample_name {name} {
    set fd [open "|get_comment -c NAME < \"$name\""]
    set id [gets $fd]
    if {[catch {close $fd}]} {
	return ""
    }

    regsub -all / $id _ id
    return $id
}

#----------
# Getting an SCF file name
proc ABI_filename_to_scf_name {name} {
    global init::use_sample_name
    global fofn_dir

    if {$init::use_sample_name} {
	set new_name [ABI_get_sample_name $name]
	if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return [file join $fofn_dir [file tail $new_name.scf]]
	}
    }

    regexp {(.*/)?(.*)} $name _dummy d f
    regsub -all { } $f {_} f
    set name $d$f
    if {[regsub {\.?([Aa][Bb][Ii1])?(\.(gz|bz|bz2|Z|z))?$} $name {.scf} new_name] == 1} {
	return [file join $fofn_dir [file tail $new_name]]
    }
    return [file join $fofn_dir [file tail $name.scf]]
}

# 31/05/00 johnt added biolims support
# is this function required ?
proc BIO_filename_to_scf_name {name} {
    global fofn_dir

    #always use biolims lane name
    #get leaf
    set lane [string range $name [expr [string last / $name]+1] end]
    #remove any . extension
    if { [string last . $lane] != -1 } {
	set lane [string range $lane 0 [expr [string last . $lane]-1]]
    }

    return [file join $fofn_dir $lane.scf]
}

proc ALF_filename_to_scf_name {name} {
    global init::use_sample_name
    global fofn_dir

    if {$init::use_sample_name} {
	set new_name [ALF_get_sample_name $name]
	if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return [file join $fofn_dir [file tail $new_name.scf]]
	}
    }

    regexp {(.*/)?(.*)} $name _dummy d f
    regsub -all { } $f {_} f
    set name $d$f
    if {[regsub {\.?([Aa][Ll][Ff])?(\.(gz|bz|bz2|Z|z))?$} $name {.scf} new_name] == 1} {
	return [file join $fofn_dir [file tail $new_name]]
    }
    return [file join $fofn_dir [file tail $name.scf]]
}

proc SCF_filename_to_scf_name {name} {
    global fofn_dir

    regexp {(.*/)?(.*)} $name _dummy d f
    regsub -all { } $f {_} f
    set name $d$f
    if {[regsub {\.?([Ss][Cc][Ff])?(\.(gz|bz|bz2|Z|z))?$} $name {.scf} new_name] == 1} {
	if {$name == $new_name} {
	    regsub {\.scf$} $name {..scf} new_name
	}
	return [file join $fofn_dir [file tail $new_name]]
    }
    return [file join $fofn_dir [file tail $name.scf]]
}

proc ZTR_filename_to_scf_name {name} {
    global init::use_sample_name
    global fofn_dir

    if {$init::use_sample_name} {
	set new_name [ZTR_get_sample_name $name]
	if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return [file join $fofn_dir [file tail $new_name.scf]]
	}
    }

    regexp {(.*/)?(.*)} $name _dummy d f
    regsub -all { } $f {_} f
    set name $d$f
    if {[regsub {\.?([Zz][Tt][Rr])?(\.(gz|bz|bz2|Z|z))?$} $name {.scf} new_name] == 1} {
	return [file join $fofn_dir [file tail $new_name]]
    }
    return [file join $fofn_dir [file tail $name.scf]]
}

proc ZTR_filename_to_ctf_name {name} {
    global init::use_sample_name
    global fofn_dir

    if {$init::use_sample_name} {
	set new_name [ZTR_get_sample_name $name]
	if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return [file join $fofn_dir [file tail $new_name.ctf]]
	}
    }

    regexp {(.*/)?(.*)} $name _dummy d f
    regsub -all { } $f {_} f
    set name $d$f
    if {[regsub {\.?([Zz][Tt][Rr])?(\.(gz|bz|bz2|Z|z))?$} $name {.ctf} new_name] == 1} {
	return [file join $fofn_dir [file tail $new_name]]
    }
    return [file join $fofn_dir [file tail $name.ctf]]
}

proc CTF_filename_to_scf_name {name} {
    global init::use_sample_name
    global fofn_dir

    if {$init::use_sample_name} {
	set new_name [CTF_get_sample_name $name]
	if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return [file join $fofn_dir [file tail $new_name.scf]]
	}
    }

    regexp {(.*/)?(.*)} $name _dummy d f
    regsub -all { } $f {_} f
    set name $d$f
    if {[regsub {\.?([Cc][Tt][Ff])?(\.(gz|bz|bz2|Z|z))?$} $name {.scf} new_name] == 1} {
	return [file join $fofn_dir [file tail $new_name]]
    }
    return [file join $fofn_dir [file tail $name.scf]]
}

#----------
# Getting a ZTR file name
proc ABI_filename_to_ztr_name {name} {
    global init::use_sample_name
    global fofn_dir

    if {$init::use_sample_name} {
	set new_name [ABI_get_sample_name $name]
	if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return [file join $fofn_dir [file tail $new_name.ztr]]
	}
    }

    regexp {(.*/)?(.*)} $name _dummy d f
    regsub -all { } $f {_} f
    set name $d$f
    if {[regsub {\.?([Aa][Bb][Ii1])?(\.(gz|bz|bz2|Z|z))?$} $name {.ztr} new_name] == 1} {
	return [file join $fofn_dir [file tail $new_name]]
    }
    return [file join $fofn_dir [file tail $name.ztr]]
}

proc ALF_filename_to_ztr_name {name} {
    global init::use_sample_name
    global fofn_dir

    if {$init::use_sample_name} {
	set new_name [ALF_get_sample_name $name]
	if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return [file join $fofn_dir [file tail $new_name.ztr]]
	}
    }

    regexp {(.*/)?(.*)} $name _dummy d f
    regsub -all { } $f {_} f
    set name $d$f
    if {[regsub {\.?([Aa][Ll][Ff])?(\.(gz|bz|bz2|Z|z))?$} $name {.ztr} new_name] == 1} {
	return [file join $fofn_dir [file tail $new_name]]
    }
    return [file join $fofn_dir [file tail $name.ztr]]
}

proc SCF_filename_to_ztr_name {name} {
    global fofn_dir

    regexp {(.*/)?(.*)} $name _dummy d f
    regsub -all { } $f {_} f
    set name $d$f
    if {[regsub {\.?([Ss][Cc][Ff])?(\.(gz|bz|bz2|Z|z))?$} $name {.ztr} new_name] == 1} {
	if {$name == $new_name} {
	    regsub {\.ztr$} $name {..ztr} new_name
	}
	return [file join $fofn_dir [file tail $new_name]]
    }
    return [file join $fofn_dir [file tail $name.ztr]]
}

proc ZTR_filename_to_ztr_name {name} {
    global fofn_dir

    regexp {(.*/)?(.*)} $name _dummy d f
    regsub -all { } $f {_} f
    set name $d$f
    if {[regsub {\.?([Zz][Tt][Rr])?(\.(gz|bz|bz2|Z|z))?$} $name {.ztr} new_name] == 1} {
	if {$name == $new_name} {
	    regsub {\.ztr$} $name {..ztr} new_name
	}
	return [file join $fofn_dir [file tail $new_name]]
    }
    return [file join $fofn_dir [file tail $name.ztr]]
}

#----------
# Getting a CTF file name
proc ABI_filename_to_ctf_name {name} {
    global init::use_sample_name
    global fofn_dir

    if {$init::use_sample_name} {
	set new_name [ABI_get_sample_name $name]
	if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return [file join $fofn_dir [file tail $new_name.ctf]]
	}
    }

    regexp {(.*/)?(.*)} $name _dummy d f
    regsub -all { } $f {_} f
    set name $d$f
    if {[regsub {\.?([Aa][Bb][Ii1])?(\.(gz|bz|bz2|Z|z))?$} $name {.ctf} new_name] == 1} {
	return [file join $fofn_dir [file tail $new_name]]
    }
    return [file join $fofn_dir [file tail $name.ctf]]
}

proc ALF_filename_to_ctf_name {name} {
    global init::use_sample_name
    global fofn_dir

    if {$init::use_sample_name} {
	set new_name [ALF_get_sample_name $name]
	if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return [file join $fofn_dir [file tail $new_name.ctf]]
	}
    }

    regexp {(.*/)?(.*)} $name _dummy d f
    regsub -all { } $f {_} f
    set name $d$f
    if {[regsub {\.?([Aa][Ll][Ff])?(\.(gz|bz|bz2|Z|z))?$} $name {.ctf} new_name] == 1} {
	return [file join $fofn_dir [file tail $new_name]]
    }
    return [file join $fofn_dir [file tail $name.ctf]]
}

proc SCF_filename_to_ctf_name {name} {
    global fofn_dir

    regexp {(.*/)?(.*)} $name _dummy d f
    regsub -all { } $f {_} f
    set name $d$f
    if {[regsub {\.?([Ss][Cc][Ff])?(\.(gz|bz|bz2|Z|z))?$} $name {.ctf} new_name] == 1} {
	if {$name == $new_name} {
	    regsub {\.ctf$} $name {..ctf} new_name
	}
	return [file join $fofn_dir [file tail $new_name]]
    }
    return [file join $fofn_dir [file tail $name.ctf]]
}

#----------
# Getting an Experiment File filename
proc ABI_filename_to_exp_name {name} {
    global init::use_sample_name
    global fofn_dir

    if {$init::use_sample_name} {
	set new_name [ABI_get_sample_name $name]
	if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return [file join $fofn_dir [file tail $new_name.exp]]
	}
    }

    set n [file tail $name]
    regsub -all "\[ \t\n\]" $n {_} name
    if {[regsub {\.?([Aa][Bb][Ii1])?(\.(gz|bz|bz2|Z|z))?$} $name {.exp} new_name] == 1} {
	return [file join $fofn_dir $new_name]
    }
    return [file join $fofn_dir $name.exp]
}

# 31/05/00 johnt added biolims support
# is this function required ?
proc BIO_filename_to_exp_name {name} {
    global fofn_dir

    #always use biolims lane name
    #get leaf
    set lane [string range $name [expr [string last / $name]+1] end]
    #remove any . extension
    if { [string last . $lane] != -1 } {
	set lane [string range $lane 0 [expr [string last . $lane]-1]]
    }

    return [file join $fofn_dir $lane.exp]
}

proc ALF_filename_to_exp_name {name} {
    global init::use_sample_name
    global fofn_dir

    if {$init::use_sample_name} {
	set new_name [ALF_get_sample_name $name]
	if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return [file join $fofn_dir [file tail $new_name.exp]]
	}
    }

    set n [file tail $name]
    regsub -all "\[ \t\n\]" $n {_} name
    if {[regsub {\.?([Aa][Ll][Ff])?(\.(gz|bz|bz2|Z|z))?$} $name {.exp} new_name] == 1} {
	return [file join $fofn_dir $new_name]
    }
    return [file join $fofn_dir $name.exp]
}

proc SCF_filename_to_exp_name {name} {
    global init::use_sample_name
    global fofn_dir

    if {$init::use_sample_name} {
        set new_name [SCF_get_sample_name $name]
        if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return [file join $fofn_dir [file tail $new_name.exp]]
        }
    }

    set n [file tail $name]
    regsub -all "\[ \t\n\]" $n {_} name
    if {[regsub {\.*([Ss][Cc][Ff])?(\.(gz|bz|bz2|Z|z))?$} $name {.exp} new_name] == 1} {
	return [file join $fofn_dir $new_name]
    }
    return [file join $fofn_dir $name.exp]
}


proc ZTR_filename_to_exp_name {name} {
    global init::use_sample_name
    global fofn_dir

    if {$init::use_sample_name} {
        set new_name [ZTR_get_sample_name $name]
        if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return [file join $fofn_dir [file tail $new_name.exp]]
        }
    }

    set n [file tail $name]
    regsub -all "\[ \t\n\]" $n {_} name
    if {[regsub {\.*([Zz][Tt][Rr])?(\.(gz|bz|bz2|Z|z))?$} $name {.exp} new_name] == 1} {
	return [file join $fofn_dir $new_name]
    }
    return [file join $fofn_dir $name.exp]
}


proc CTF_filename_to_exp_name {name} {
    global init::use_sample_name
    global fofn_dir

    if {$init::use_sample_name} {
        set new_name [CTF_get_sample_name $name]
        if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return [file join $fofn_dir [file tail $new_name.exp]]
        }
    }

    set n [file tail $name]
    regsub -all "\[ \t\n\]" $n {_} name
    if {[regsub {\.*([Cc][Tt][Ff])?(\.(gz|bz|bz2|Z|z))?$} $name {.exp} new_name] == 1} {
	return [file join $fofn_dir $new_name]
    }
    return [file join $fofn_dir $name.exp]
}


#----------
# Getting a sequence entry (reading) name
proc ABI_filename_to_entry_name {name} {
    global init::use_sample_name

    if {$init::use_sample_name} {
	set new_name [ABI_get_sample_name $name]
	if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return $new_name
	}
    }

    regsub -all "\[ \t\n\]" [file tail $name] {_} name
    if {[regsub {\.?([Aa][Bb][Ii1])?(\.(gz|bz|bz2|Z|z))?$} $name {} new_name] == 1} {
	return $new_name
    }
    return $name
}

proc ALF_filename_to_entry_name {name} {
    global init::use_sample_name

    if {$init::use_sample_name} {
	set new_name [ALF_get_sample_name $name]
	if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
	    return $new_name
	}
    }

    regsub -all "\[ \t\n\]" [file tail $name] {_} name
    if {[regsub {\.?([Aa][Ll][Ff])?(\.(gz|bz|bz2|Z|z))?$} $name {} new_name] == 1} {
	return $new_name
    }
    return $name
}

proc SCF_filename_to_entry_name {name} {
    global init::use_sample_name

    if {$init::use_sample_name} {
        set new_name [SCF_get_sample_name $name]
        if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
            return $new_name
        }
    }

    regsub -all "\[ \t\n\]" [file tail $name] {_} name
    if {[regsub {\.*([Ss][Cc][Ff])?(\.(gz|bz|bz2|Z|z))?$} $name {} new_name] == 1} {
	return $new_name
    }
    return $name
}

proc ZTR_filename_to_entry_name {name} {
    global init::use_sample_name

    if {$init::use_sample_name} {
        set new_name [ZTR_get_sample_name $name]
        if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
            return $new_name
        }
    }

    regsub -all "\[ \t\n\]" [file tail $name] {_} name
    if {[regsub {\.*([Zz][Tt][Rr])?(\.(gz|bz|bz2|Z|z))?$} $name {} new_name] == 1} {
	return $new_name
    }
    return $name
}

proc CTF_filename_to_entry_name {name} {
    global init::use_sample_name

    if {$init::use_sample_name} {
        set new_name [CTF_get_sample_name $name]
        if {$new_name != ""} {
	    regsub -all "\[/ \t\n\]" $new_name {_} new_name
            return $new_name
        }
    }

    regsub -all "\[ \t\n\]" [file tail $name] {_} name
    if {[regsub {\.*([Cc][Tt][Ff])?(\.(gz|bz|bz2|Z|z))?$} $name {} new_name] == 1} {
	return $new_name
    }
    return $name
}

proc EXP_filename_to_entry_name {name} {
    regsub -all "\[ \t\n\]" [file tail $name] {_} name_
    if {[catch {set fd [open $name r]}]} {
	return $name_
    }

    set name $name_
    while {[gets $fd line] != -1} {
	set type [lindex $line 0]
	if {$type == "ID" || $type == "EN"} {
	    set name [lindex $line 1]
	    break;
        }
    }

    close $fd
    return $name
}

proc PLN_filename_to_entry_name {name} {
    regsub -all "\[ \t\n\]" [file tail $name] {_} name
    return $name
}

proc PLN_filename_to_exp_name {name} {
    regsub -all "\[ \t\n\]" [file tail $name] {_} name 
    return $name.exp
}

proc UNK_filename_to_entry_name {name} {
    regsub -all "\[ \t\n\]" [file tail $name] {_} name
    return $name
}

#31/05/00 johnt - added biolims support
proc BIO_filename_to_entry_name {name} {
    # the leaf of the path is the lane name, which
    # seems to be the same as the sample name
    set fd [string range $name [expr [string last / $name]+1] end]

    return $fd
}
