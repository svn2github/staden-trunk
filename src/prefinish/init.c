#include <tcl.h>

#include "tkEditor.h"
#include "tkEdNames.h"
#include "tkSheet.h"
#include "contigEditor.h"
#include "newgap_cmds.h"
#include "gap-tcl.h"
#include "active_tags.h"
#include "finish.h"

int Gap_Init(Tcl_Interp *interp) {
    Editor_Init(interp);
    EdNames_Init(interp);
    Sheet_Init(interp);
    Ced_Init(interp);
    NewGap_Init(interp);
    Db_Init(interp);
    Prefinish_Init(interp);

    get_tag_types();

    return Tcl_PkgProvide(interp, "gap4", "1.0");
}
