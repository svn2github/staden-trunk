#ifndef _edCommands_h
#define _edCommands_h

extern void dumpContig(EdStruct *xx, char *fn, int left, int right, int llen,
		       int nwidth);
extern void save_consensus_trace(EdStruct *xx, char *fn, int left, int right,
				 int strand, int matching);
extern int align_read(EdStruct *xx);
extern int strip_pads(EdStruct *xx, int con_mode, float con_cutoff);

#endif
