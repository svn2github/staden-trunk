#####################
package require Iwidgets
source labeledwidget.itk
source dataentry.itk

dataentry .test0 -labeltext "Entrydata1" -default 20 -range {10 100} \
	-width 20 \
	-sticky e \
	-validate integer
pack .test0 -fill both


dataentry .test1 -labeltext "Entrydata2" -default 50 -range {50 500} \
	-width 20 \
	-sticky e \
	-validate real
pack .test1 -fill both

dataentry .test2 -labeltext "Depends on entrydata2 " -default 16 -range {10 100} \
	-width 20 \
	-sticky e \
	-validate integer
pack .test2 -fill both

dataentry .test3 -labeltext "Entrydata3" \
	-width 20 \
	-sticky e 
pack .test3 -fill both

#
# Records the default range setting, so that it can be re_used when needed.
#
set o_range [.test2 cget -range]
bind .test1 <Any-Leave> {test1_changed "$o_range"}
	                 
frame .separator -bd 2 -relief sunken -height 2
pack .separator -fill x 

iwidgets::buttonbox .ok     
.ok add ok -text OK -command "puts \[.test0 get\]; puts \[.test1 get\]; \
	                      puts \[.test2 get\]; puts \[.test3 get\]"

pack .ok -padx 2 -pady 2 -side bottom -fill x

proc test1_changed {default} {

    set current [.test1 get]
    set value [.test2 get]

    if {$current > 500} {
	.test2 configure -range {60 216}
	if {$value > 216 || $value < 60} {
	.test2.lwchildsite.entry configure -background mistyrose    
	}
    } else {
	.test2 configure -range $default
	set new_range [.test2 cget -range]
	set new_min [lindex $new_range 0]
	set new_max [lindex $new_range 1]
	set new_value [.test2 get]
	if {$new_value < $new_min || $new_value > $new_max} {	
	    .test2.lwchildsite.entry configure -background mistyrose
	} else {
	    .test2.lwchildsite.entry configure -background white 
	}
    } 
}
