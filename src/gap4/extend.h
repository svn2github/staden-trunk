#ifndef _extend_h
#define _extend_h

#include "edUtils.h"

#define EXTEND_LEFT 1
#define EXTEND_RIGHT 2

extern int extend(EdStruct *xx, int seq, int dir);
extern int unextend(EdStruct *xx, int seq, int dir);
extern int undo_unextend(EdStruct *xx, int seq, int dir, int time);

/*
 * Zap (unextend) to right end
 */
void zap_Right (EdStruct *xx);

/*
 * Zap (unextend) to left end
 */
void zap_Left (EdStruct *xx);

/*
 * Handle cut-off adjust
 */
int meta_arrow (EdStruct *xx, int key);

/* Low level - avoid if possible */
int _adjust_ends(DBInfo *db, int seq, int start_bases, int end_bases,
		 int seq_flags);
#endif /* _extend_h */
