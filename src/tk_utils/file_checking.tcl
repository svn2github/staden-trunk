#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1995. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#
##############################################################################
#put the current directory path in the Selection entry
proc ClearSel { path } {
    entrybox_delete $path 0 end
    entrybox_focus $path
}
#end ClearSel
