#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
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
	free(mm->mapped_reads);
	free(mm);
}
void maqmap_write_header(gzFile fp, const maqmap_t *mm)
{
	int i, len;
	gzwrite(fp, &mm->format, sizeof(int));
	gzwrite(fp, &mm->n_ref, sizeof(int));
	for (i = 0; i != mm->n_ref; ++i) {
		len = strlen(mm->ref_name[i]) + 1;
		gzwrite(fp, &len, sizeof(int));
		gzwrite(fp, mm->ref_name[i], len);
	}
	gzwrite(fp, &mm->n_mapped_reads, sizeof(bit64_t));
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
	mm->ref_name = (char**)calloc(mm->n_ref, sizeof(char*));
	for (k = 0; k != mm->n_ref; ++k) {
		gzread(fp, &len, sizeof(int));
		mm->ref_name[k] = (char*)malloc(len * sizeof(char));
		gzread(fp, mm->ref_name[k], len);
	}
	/* read number of mapped reads */
	gzread(fp, &mm->n_mapped_reads, sizeof(bit64_t));
	return mm;
}

/* mapview */

static void mapview_core(FILE *fpout, gzFile fpin)
{
	bit64_t i, j;
	maqmap_t *m = maqmap_read_header(fpin);
	for (i = 0; i != m->n_mapped_reads; ++i) {
		maqmap1_t *m1, mm1;
		m1 = &mm1;
		maqmap_read1(fpin, m1);
		fprintf(fpout, "%s\t%s\t%d\t%c\t%d\t%u\t%d\t%d\t%d\t%d\t%d\t%d\t",
				m1->name, m->ref_name[m1->seqid], (m1->pos>>1) + 1,
				(m1->pos&1)? '-' : '+', m1->dist, m1->flag, m1->map_qual, m1->seq[MAX_READLEN-1],
				m1->alt_qual, m1->i1>>5, m1->i2>>5, m1->size);
		for (j = 0; j != m1->size; ++j) {
			if (m1->seq[j] == 0) fputc('n', fpout);
			else if ((m1->seq[j]&0x3f) < 27) fputc("acgt"[m1->seq[j]>>6&3], fpout);
			else fputc("ACGT"[m1->seq[j]>>6&3], fpout);
		}
		fputc('\t', fpout);
		for (j = 0; j != m1->size; ++j)
			fputc((m1->seq[j]&0x3f) + 33, fpout);
		fputc('\n', fpout);
	}
	maq_delete_maqmap(m);
}

int ma_mapview(int argc, char *argv[])
{
	if (argc == 1) {
		fprintf(stderr, "Usage: maq mapview <in.map>\n");
		return 1;
	}
	gzFile fp = (strcmp(argv[optind], "-") == 0)? gzdopen(STDIN_FILENO, "r") : gzopen(argv[optind], "r");
	mapview_core(stdout, fp);
	gzclose(fp);
	return 0;
}
