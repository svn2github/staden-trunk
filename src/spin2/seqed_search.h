#ifndef _SEQED_SEARCH_H_
#define _SEQED_SEARCH_H_

#include "tkSeqed.h"

int seqed_string_search_free(void);
int seqed_next_string(tkSeqed *se, int direction);
int seqed_string_search(char *sequence, int seq_len, char *seq_name, 
			char *string, int direction, int strand,
			double per_match, int cursorPos, int use_iub_code);
#endif
