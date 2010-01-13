#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>

#include "maqmap.h"

maqmap_t *maq_new_maqmap(void)
{
	maqmap_t *mm = (maqmap_t*)calloc(1, sizeof(maqmap_t));
	mm->format = MAQMAP_FORMAT_NEW;
	return mm;
}

void maq_delete_maqmap(maqmap_t *mm)
{
	int i;
	if (mm == 0) return;
	for (i = 0; i < mm->n_ref; ++i)
		free(mm->ref_name[i]);
	free(mm->ref_name);
	free(mm);
}

maqmap_t *maqmap_read_header(gzFile fp)
{
	maqmap_t *mm;
	int k, len;
	mm = maq_new_maqmap();
	gzread(fp, &mm->format, sizeof(int));
	if (mm->format != MAQMAP_FORMAT_NEW) {
		if (mm->format > 0) {
			fprintf(stderr, "** Obsolete map format is detected. Please use 'mapass2maq' command to convert the format.\n");
			exit(3);
		}
		assert(mm->format == MAQMAP_FORMAT_NEW);
	}
	gzread(fp, &mm->n_ref, sizeof(int));
	mm->ref_name = (char**)calloc(mm->n_ref+1, sizeof(char*));
	for (k = 0; k != mm->n_ref; ++k) {
		gzread(fp, &len, sizeof(int));
		mm->ref_name[k] = (char*)malloc(len * sizeof(char));
		gzread(fp, mm->ref_name[k], len);
	}
	/* read number of mapped reads */
	gzread(fp, &mm->n_mapped_reads, sizeof(bit64_t));
	return mm;
}

/*
 * Auto detects size to distinguish between short and long maq variants
 * Returns 64 for short
 *         128 for long
 *         -1 for error
 */
int maq_detect_size(gzFile fp) {
    z_off_t curr = gztell(fp);
    maqmap64_t m;
    int i;

    if (gzread(fp, &m, sizeof(maqmap64_t)) == -1)
	return -1;

    gzseek(fp, curr, SEEK_SET);

    if (m.size > 64)
	return 128;
    if (m.size <= 0)
	return 128;

    for (i = 0; i < MAX_NAMELEN; i++) {
	if (!m.name[i])
	    break;
	if (!isprint(m.name[i]))
	    return 128;
    }

    return 64;
}
