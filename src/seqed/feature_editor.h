#if !defined(FEATURE_EDITOR_H)
#define FEATURE_EDITOR_H

typedef struct {
    /* Must be first item */
    char *type;
    /* values taken from FQUALIFIERDB file */
    char *qual; /*search_id*/
    char *value;
    char *entry;
    char id[4];
} feat_qual_db;

int tcl_get_qual_array(ClientData clientData, 
		       Tcl_Interp *interp,
		       int argc, 
		       char **argv);
int FeatureEditor_Init(Tcl_Interp *interp);

#endif
