#ifndef _SEQEDTRANSLATE_H_
#define _SEQEDTRANSLATE_H_

#include "tkSeqedUtils.h"

void seqed_write_translation(char *sequence,
			     int frame,
			     int size,
			     int pos,
			     int line_length,
			     int overlap,
			     char *line);

void
seqed_translate_frame(tkSeqed *se,
		      char *sequence,
		      int pos,
		      int width,
		      int frame,
		      char *sline,
		      char *name,
		      int size,
		      XawSheetInk *splodge);

int
find_auto_lines(region **exons,
		int num_exons,
		int width,
		int pos,
		int complement);

void
seqed_auto_translate(tkSeqed *se,
		     char *sequence,
		     int pos,
		     int width,
		     char *sline,
		     char *name,
		     XawSheetInk *splodge,
		     int size,
		     region *exon,
		     int exon_num,
		     int start,
		     int end,
		     int line_num,
		     int complement);

void
seqed_auto_translateOLD(tkSeqed *se,
		     char *sequence,
		     int pos,
		     int width,
		     char *sline,
		     char *name,
		     XawSheetInk *splodge,
		     int size,
		     region *exon,
		     int num_exon);
#endif
