#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
##############################################################################
#scroll all active windows in the x direction
proc scroll_x { f io id command x args} {
    global $f.template_id

    scroll_canvas -io $io -id $id -xscrollcommand "$command $x $args"
} 
#end scroll_x

##############################################################################
#scroll all active windows in the y direction
proc scroll_y { f io id command y args} {
    global NGWinType

    scroll_canvas -io $io -id $id -yscrollcommand "$command $y $args"
    return
} 
#end scroll_y

