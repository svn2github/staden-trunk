#ifndef _NIP_SENDTO_H_
#define _NIP_SENDTO_H_
int sender(char *rid, Tcl_Interp *interp);
int nip_sender(Tcl_Interp *interp, char *rid, int seq_id);
#endif
