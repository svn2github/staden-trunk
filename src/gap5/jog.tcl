;proc jog_callback {w call auto {pos {}}} {
    global $w.After $w.Down $w.Waiting $w.Rept

    if {$auto == 0 && [set $w.Waiting] == 1} return;
    catch {after cancel [set $w.After]}
    if {![set $w.Down]} return

    if {$pos == {}} {
	set pos [$w get]
    }
    set sign [expr {($pos-1000) >= 0 ? 1 : -1}]
    set val  [expr {abs(($pos-1000)/30.0)}]
    set dist [expr {$val >= 1 ? int($sign*pow(1.5, $val-1)) : 0}]
    # puts $dist
    
    eval $call $dist
    set $w.Waiting 1
    set $w.After [after [set $w.Rept] [list jog_callback $w $call 1]]
}

#
# The main job interface
#
proc jog {w args} {
    global $w.Down $w.Waiting $w.Rept

    # Extract our own options
    array set opts $args
    if {[info exists opts(-command)]} {
	set cmd $opts(-command)
    } else {
	bgerror "No -command option specified to jog widget"
	set cmd puts
    }

    if {[info exists opts(-repeatinterval)]} {
	set $w.Rept $opts(-repeatinterval)
    } else {
	set $w.Rept 100
    }
    
    # Scale options we reset
    set opts(-repeatinterval 5)
    set opts(-command) [list jog_callback $w $cmd 0]

    set $w.Down 0
    set $w.Waiting 0

    eval scale $w -showvalue 0 -orient horiz [array get opts] -from 0 -to 2000
    $w set 1000

    bind $w <ButtonPress-1> {
	global %W.Down
	set %W.Down 1
    }

    bind $w <ButtonRelease-1> {
	global %W.After %W.Down %w.Waiting
	catch {after cancel [set %W.After]}
	set %W.Down 0
	set %W.Waiting 0
	%W set 1000
    }
    
    catch {after cancel [set %W.After]}

    return $w
}

# proc foo {args} {
#     puts [info level [info level]]
# }
# 
# jog .foo \
#     -bd 2 \
#     -relief sunken \
#     -length 400 \
#     -command "foo fish" \
#     -repeatinterval 300
# pack .foo -expand 1 
