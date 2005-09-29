# Tcl package index file, version 1.0

if {$tcl_platform(platform) == "windows"} {
    package ifneeded Itk 3.3 [list load [file join $dir "itk33.dll"] Itk]
} else {
    package ifneeded Itk 3.3 [list load [file join $dir "libitk3.3.so"] Itk]
}
