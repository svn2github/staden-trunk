#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
proc do_repeat {state script} {
    upvar #0 $state s
    if {$s > 0} {
        eval $script

        if {$s == 1} {
	    set s 2
	    after 500 "do_repeat $state \"$script\""
	} else {
	    after 50 "do_repeat $state \"$script\""
        }
    }
}


proc repeater {path callback args} {
    eval button $path $args

    bind $path <ButtonPress-1> "+
        set $path.State 1
        do_repeat $path.State \"$callback\"
    "
    
    bind $path <ButtonRelease-1> "+
        set $path.State 0
	after cancel \"do_repeat $path.State \\\"$callback\\\"\"
    "

    return $path
}
