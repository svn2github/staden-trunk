# Tcl package index file, version 1.0

if {$tcl_platform(platform) == "windows"} {
    package ifneeded Itcl 3.3 [list load [file join $dir "itcl33.dll"] Itcl]
} else {
    package ifneeded Itcl 3.3 [list load [file join $dir "libitcl3.3.so"] Itcl]
}

