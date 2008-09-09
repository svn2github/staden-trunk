/*
 * ======================================================================
 * This software has been created by Genome Research Limited (GRL).
 *
 * GRL hereby grants permission to use, copy, modify and distribute
 * this software and its documentation for non-commercial purposes
 * without fee at the user's own risk on the basis set out below.
 *
 * GRL neither undertakes nor accepts any duty whether contractual or
 * otherwise in connection with the software, its use or the use of
 * any derivative, and makes no representations or warranties, express
 * or implied, concerning the software, its suitability, fitness for
 * a particular purpose or non-infringement.
 *
 * In no event shall the authors of the software or GRL be responsible
 * or liable for any loss or damage whatsoever arising in any way
 * directly or indirectly out of the use of this software or its
 * derivatives, even if advised of the possibility of such damage.
 *
 * Our software can be freely distributed under the conditions set out
 * above, and must contain this copyright notice.
 * ======================================================================
 */

/*
 * This performs a linear (non-indexed) search for a trace in an SRF archive.
 *
 * It's not intended as a suitable production program or as a library of code
 * to use, but as a test and benchmark statistic.
 */

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <io_lib/Read.h>
#include <io_lib/misc.h>
#include <io_lib/ztr.h>
#include <io_lib/srf.h>
#include <io_lib/hash_table.h>

#define LEVEL_READ  (1 << 0)
#define LEVEL_CHUNK (1 << 1)
#define LEVEL_NAME  (1 << 2)
#define LEVEL_ALL   (LEVEL_READ | LEVEL_CHUNK | LEVEL_NAME);

/* only checks the first 10 traces */
#define LEVEL_CHECK 255 

#define READ_TOTAL 0
#define READ_GOOD  1
#define READ_BAD   2

#define NREADS 3

#define READ_TOTAL_STR  "TOTAL"
#define READ_GOOD_STR   "GOOD"
#define READ_BAD_STR    "BAD"

/* see ztr.h for a list of all possible ztr chunk types */

#define CHUNK_BASE 0
#define CHUNK_CNF1 1
#define CHUNK_CNF4 2
#define CHUNK_SAMP 3
#define CHUNK_SMP4 4
#define CHUNK_REGN 5

#define NCHUNKS 6

#define CHUNK_BASE_TYPE ZTR_TYPE_BASE
#define CHUNK_CNF1_TYPE ZTR_TYPE_CNF1
#define CHUNK_CNF4_TYPE ZTR_TYPE_CNF4
#define CHUNK_SAMP_TYPE ZTR_TYPE_SAMP
#define CHUNK_SMP4_TYPE ZTR_TYPE_SMP4
#define CHUNK_REGN_TYPE ZTR_TYPE_REGN

#define KEY_TYPE    0
#define KEY_VALTYPE 1
#define KEY_GROUP   2
#define KEY_OFFS    3
#define KEY_SCALE   4
#define KEY_COORD   5
#define KEY_NAME    6

#define NKEYS 7

#define KEY_TYPE_STR    "TYPE"
#define KEY_VALTYPE_STR "VALTYPE"
#define KEY_GROUP_STR   "GROUP"
#define KEY_OFFS_STR    "OFFS"
#define KEY_SCALE_STR   "SCALE"
#define KEY_COORD_STR   "COORD"
#define KEY_NAME_STR    "NAME"

#define TYPE_PROC 0
#define TYPE_SLXI 1
#define TYPE_SLXN 2
#define TYPE_0FAM 3
#define TYPE_1CY3 4
#define TYPE_2TXR 5
#define TYPE_3CY5 6

#define NTYPES 7

#define TYPE_PROC_STR "PROC"
#define TYPE_SLXI_STR "SLXI"
#define TYPE_SLXN_STR "SLXN"
#define TYPE_0FAM_STR "0FAM"
#define TYPE_1CY3_STR "1CY3"
#define TYPE_2TXR_STR "2TXR"
#define TYPE_3CY5_STR "3CY5"

/* regn chunk */
typedef struct {
    char coord;
    char *name;
    uint4 *bndy;
    int nbndy;
    int count;
} regn_t;


/* ------------------------------------------------------------------------ */

/*
 * Print usage message to stderr and exit with the given \"code\".
 */
void usage(int code) {
    fprintf(stderr, "Usage: srf_info [-level level_bitmap] input(s)\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -l level_bitmap \n");
    fprintf(stderr, "              1\tCount of good/bad reads.\n");
    fprintf(stderr, "              2\tCounts for selected chunk types.\n");
    fprintf(stderr, "              4\tTrace count and trace name prefix for each trace_header.\n");
    fprintf(stderr, "\n");

    exit(code);
}

/*
 * Parse the comma delimited list of chunk types and put them in the single character \"mode\".
 *
 * Returns 0 on success.
 */
int parse_regn(ztr_t *z, ztr_chunk_t *chunk, HashTable *regn_hash) {
    char key[1024];
    char *name;
    HashItem *hi;
    
    /* the hash key is a combination of the region names and boundaries */
    name = ztr_lookup_mdata_value(z, chunk, "NAME");
    sprintf(key, "%s+", name);
    if( chunk->dlength ){
	int i,j;
	for (i=1,j=strlen(key); i<chunk->dlength; i++, j++)
	    if(chunk->data[i])
		key[j] = chunk->data[i];
	key[j] = '\0';
    }
    
    if (NULL == (hi = (HashTableSearch(regn_hash, key, strlen(key))))) {
	int nbndy = 0;
	char *coord;
	HashData hd;
	regn_t *regn;

	if( chunk->dlength )
	    nbndy = (chunk->dlength-1)/4;
	coord = ztr_lookup_mdata_value(z, chunk, "COORD");

	if( NULL == (regn = (regn_t *)malloc(sizeof(regn_t)))) {
	    return 1;
	}
	regn->coord = (NULL == coord ? 'B' : *coord );
	regn->name = strdup(name);
	if( NULL == (regn->bndy = malloc(nbndy * sizeof(int4)))) {
	    free(regn->name);
	    free(regn);
	    return 1;
	}
	if( nbndy )
	    memcpy(regn->bndy, chunk->data+1, nbndy*4);
	regn->nbndy = nbndy;
	regn->count = 1;
            
	hd.p = regn;
	if (NULL == (hi = HashTableAdd(regn_hash, key, strlen(key), hd, NULL))) {
	    free(regn->bndy);
	    free(regn->name);
	    free(regn);
	    return 1;
	}
    } else {
	regn_t *regn = (regn_t *)(hi->data.p);
	regn->count++;
    }

    return 0;
}

/*
 * count the mdata keys
 *
 * Returns 0 on success.
 */
int count_mdata_keys(ztr_t *z, ztr_chunk_t *chunk, int ichunk, long key_count[NCHUNKS][NKEYS], long type_count[NCHUNKS][NTYPES]) {
    char *keys_str[] = {KEY_TYPE_STR, KEY_VALTYPE_STR, KEY_GROUP_STR, KEY_OFFS_STR, KEY_SCALE_STR, KEY_COORD_STR, KEY_NAME_STR};
    char *types_str[] = {TYPE_PROC_STR, TYPE_SLXI_STR, TYPE_SLXN_STR, TYPE_0FAM_STR, TYPE_1CY3_STR, TYPE_2TXR_STR, TYPE_3CY5_STR};
    int ikey, itype;

    if (z->header.version_major > 1 ||
	z->header.version_minor >= 2) {
	/* ZTR format 1.2 onwards */

	char *cp = chunk->mdata;
	int32_t dlen = chunk->mdlength;

	/*
	 * NB: we may wish to rewrite this using a dedicated state machine
	 * instead of strlen/strcmp as this currently assumes the meta-
	 * data is correctly formatted, which we cannot assume as the 
	 * metadata is external and outside of our control.
	 * Passing in non-nul terminated strings could crash this code.
	 */
	while (dlen > 0) {
	    size_t l;

	    /* key */
	    l = strlen(cp);
	    for (ikey=0; ikey<NKEYS; ikey++)
		if(0 == strcmp(cp, keys_str[ikey]))
		    break;

	    cp += l+1;
	    dlen -= l+1;

	    /* value */
	    if (ikey < NKEYS)
		key_count[ichunk][ikey]++;

	    /* for the type key check the value */
	    if (ikey == KEY_TYPE && (ichunk == CHUNK_SAMP || ichunk == CHUNK_SMP4)) {
		for (itype=0; itype<NTYPES; itype++)
		    if(0 == strcmp(cp, types_str[itype]))
			break;
		if(itype < NTYPES)
		    type_count[ichunk][itype]++;
	    }

	    l = strlen(cp);
	    cp += l+1;
	    dlen -= l+1;
	}

    } else {
	/* v1.1 and before only supported a few types, specifically coded
	 * per chunk type.
	 */

	switch (chunk->type) {
	case ZTR_TYPE_SAMP:
	case ZTR_TYPE_SMP4:
	    key_count[ichunk][KEY_TYPE]++;
	    for (itype=0; itype<NTYPES; itype++)
		if(0 == strcmp(chunk->mdata, types_str[itype])) {
		    type_count[ichunk][itype]++;
		}
	    break;

	default:
	    break;
	}
    }

    return 0;
}

/*
 * Given the archive name and the level_mode
 * generate information about the archive
 *
 * Note the generated srf file is NOT indexed
 *
 * Returns 0 on success.
 */
int srf_info(char *input, int level_mode, long *read_count, long *chunk_count, long key_count[NCHUNKS][NKEYS], long type_count[NCHUNKS][NTYPES], HashTable *regn_hash) {
    srf_t *srf;
    char name[1024];
    long trace_body_count = 0;

    int count = 0;
    
    if (NULL == (srf = srf_open(input, "rb"))) {
	perror(input);
	return 1;
    }

    do {
	int type;

	switch(type = srf_next_block_type(srf)) {
	case SRFB_CONTAINER:
	    if( trace_body_count ){
		if( level_mode & LEVEL_NAME )
		    printf( " ... %s x%ld\n", name+strlen(srf->th.id_prefix), trace_body_count);
		trace_body_count = 0;
	    }
	    if (0 != srf_read_cont_hdr(srf, &srf->ch)) {
		fprintf(stderr, "Error reading container header.\nExiting.\n");
		exit(1);
	    }
	    break;

        case SRFB_XML:
	    if( trace_body_count ){
		if( level_mode & LEVEL_NAME )
		    printf( " ... %s x%ld\n", name+strlen(srf->th.id_prefix), trace_body_count);
		trace_body_count = 0;
	    }
	    if (0 != srf_read_xml(srf, &srf->xml)) {
		fprintf(stderr, "Error reading XML.\nExiting.\n");
		exit(1);
	    }
	    break;

	case SRFB_TRACE_HEADER:
	    if( trace_body_count ){
		if( level_mode & LEVEL_NAME )
		    printf( " ... %s x%ld\n", name+strlen(srf->th.id_prefix), trace_body_count);
		trace_body_count = 0;
	    }
	    if (0 != srf_read_trace_hdr(srf, &srf->th)) {
		fprintf(stderr, "Error reading trace header.\nExiting.\n");
		exit(1);
	    }

	    if( 0 == (level_mode & LEVEL_CHUNK) )
		break;

	    /* Decode ZTR chunks in the header */
	    if (srf->mf)
		mfdestroy(srf->mf);

	    srf->mf = mfcreate(NULL, 0);
	    if (srf->th.trace_hdr_size)
		mfwrite(srf->th.trace_hdr, 1, srf->th.trace_hdr_size, srf->mf);
	    if (srf->ztr)
		delete_ztr(srf->ztr);
	    mrewind(srf->mf);

	    if (NULL != (srf->ztr = partial_decode_ztr(srf, srf->mf, NULL))) {
		srf->mf_pos = mftell(srf->mf);
		mfseek(srf->mf, 0, SEEK_END);
		srf->mf_end = mftell(srf->mf);
	    } else {
		/* Maybe not enough to decode or no headerBlob. */
		/* So delay until decoding the body. */
		srf->mf_pos = srf->mf_end = 0;
	    }

	    break;

	case SRFB_TRACE_BODY: {
	    srf_trace_body_t old_tb;
	    ztr_t *ztr_tmp;

	    if (0 != srf_read_trace_body(srf, &old_tb, 0)) {
		fprintf(stderr, "Error reading trace body.\nExiting.\n");
		exit(1);
	    }

	    if (-1 == construct_trace_name(srf->th.id_prefix,
					   (unsigned char *)old_tb.read_id,
					   old_tb.read_id_length,
					   name, 512)) {
		fprintf(stderr, "Error constructing trace name.\nExiting.\n");
		exit(1);
	    }

	    trace_body_count++;
	    if( 1 == trace_body_count ){
		if( level_mode & LEVEL_NAME )
		    printf( "trace_name: %s + %s", srf->th.id_prefix, name+strlen(srf->th.id_prefix));
	    }
          
	    read_count[READ_TOTAL]++;

	    if (old_tb.flags & SRF_READ_FLAG_BAD_MASK ){
		read_count[READ_BAD]++;
	    } else {
		read_count[READ_GOOD]++;
	    }
          
	    if( 0 == (level_mode & LEVEL_CHUNK) )
		break;

	    if (!srf->mf) {
		fprintf(stderr, "Error reading trace body.\nExiting.\n");
		exit(1);
	    }

	    mfseek(srf->mf, srf->mf_end, SEEK_SET);
	    if (old_tb.trace_size) {
		mfwrite(old_tb.trace, 1, old_tb.trace_size, srf->mf);
		free(old_tb.trace);
		old_tb.trace = NULL;
	    }
          
	    mftruncate(srf->mf, mftell(srf->mf));
	    mfseek(srf->mf, srf->mf_pos, SEEK_SET);

	    if (srf->ztr)
		ztr_tmp = ztr_dup(srf->ztr); /* inefficient, but simple */
	    else
		ztr_tmp = NULL;

	    if (NULL != partial_decode_ztr(srf, srf->mf, ztr_tmp)) {
		int i;
		for (i=0; i<ztr_tmp->nchunks; i++) {
		    int ichunk = -1;
		    switch (ztr_tmp->chunk[i].type) {
		    case ZTR_TYPE_BASE:
			ichunk = CHUNK_BASE;
			break;
		    case ZTR_TYPE_CNF1:
			ichunk = CHUNK_CNF1;
			break;
		    case ZTR_TYPE_CNF4:
			ichunk = CHUNK_CNF4;
			break;
		    case ZTR_TYPE_SAMP:
			ichunk = CHUNK_SAMP;
			break;
		    case ZTR_TYPE_SMP4:
			ichunk = CHUNK_SMP4;
			break;
		    case ZTR_TYPE_REGN:
			ichunk = CHUNK_REGN;
			if( parse_regn(srf->ztr, &ztr_tmp->chunk[i], regn_hash) ){
			    delete_ztr(ztr_tmp);
			    return 1;
			}
			break;
		    default:
			break;
		    }

		    if( ichunk > -1 ) {
			chunk_count[ichunk]++;
			count_mdata_keys(srf->ztr, &ztr_tmp->chunk[i], ichunk, key_count, type_count);
		    }
		}

	    }

	    if( ztr_tmp )
		delete_ztr(ztr_tmp);

	    count++;
	    if( (level_mode == LEVEL_CHECK) && (count == 10) ){
		printf( " ... %s x%ld\n", name+strlen(srf->th.id_prefix), trace_body_count);
		srf_destroy(srf, 1);
		return 0;
	    }
          
	    break;
        }

        default:
	    if( trace_body_count ){
		if( level_mode & LEVEL_NAME )
		    printf( " ... %s x%ld\n", name+strlen(srf->th.id_prefix), trace_body_count);
	    }
	    break;
	}

	if( type == -1 || type == SRFB_INDEX || type == SRFB_NULL_INDEX )
	    break;

    } while (1);

    srf_destroy(srf, 1);
    return 0;
}

/* ------------------------------------------------------------------------ */

/*
 * Main method.
 */
int main(int argc, char **argv) {
    int ifile, nfiles;
    char *input = NULL;

    int c;
    int errflg = 0;
    extern char *optarg;
    extern int optind, optopt;

    int level_mode = LEVEL_ALL;

    long read_count[NREADS];
    char *read_str[] = {READ_TOTAL_STR, READ_GOOD_STR, READ_BAD_STR};
    long chunk_count[NCHUNKS];
    uint4 chunk_type[] = {CHUNK_BASE_TYPE, CHUNK_CNF1_TYPE, CHUNK_CNF4_TYPE, CHUNK_SAMP_TYPE, CHUNK_SMP4_TYPE, CHUNK_REGN_TYPE};
    long key_count[NCHUNKS][NKEYS];
    char *keys_str[] = {KEY_TYPE_STR, KEY_VALTYPE_STR, KEY_GROUP_STR, KEY_OFFS_STR, KEY_SCALE_STR, KEY_COORD_STR, KEY_NAME_STR};
    long type_count[NCHUNKS][NTYPES];
    char *types_str[] = {TYPE_PROC_STR, TYPE_SLXI_STR, TYPE_SLXN_STR, TYPE_0FAM_STR, TYPE_1CY3_STR, TYPE_2TXR_STR, TYPE_3CY5_STR};
    int iread, ichunk, ikey, itype;

    while ((c = getopt(argc, argv, "l:")) != -1) {
        switch (c) {
        case 'l':
            if (1 != sscanf(optarg, "%d", &level_mode)) {
                fprintf(stderr,
                        "Otion -%c requires an operand\n", optopt);
                errflg++;
            }
	    break;
        case ':':       /* -? without operand */
            fprintf(stderr,
                    "Option -%c requires an operand\n", optopt);
            errflg++;
            break;
        case '?':
            fprintf(stderr,
                    "Unrecognised option: -%c\n", optopt);
            errflg++;
        }
    }

    if (errflg) {
	usage(1);
    }

    nfiles = (argc-optind);
    if( nfiles < 1 ){
        fprintf(stderr, "Please specify input archive name(s).\n");
        usage(1);
    }
    
    for (ifile=0; ifile<nfiles; ifile++) {
        HashTable *regn_hash;
        char type[5];

        input = argv[optind+ifile];
        printf("Reading archive %s.\n", input);

        for (iread=0; iread<NREADS; iread++)
	    read_count[iread] = 0;

        for (ichunk=0; ichunk<NCHUNKS; ichunk++)
	    chunk_count[ichunk] = 0;

        for (ichunk=0; ichunk<NCHUNKS; ichunk++)
            for (ikey=0; ikey<NKEYS; ikey++)
                key_count[ichunk][ikey] = 0;

        for (ichunk=0; ichunk<NCHUNKS; ichunk++)
            for (itype=0; itype<NTYPES; itype++)
                type_count[ichunk][itype] = 0;

        if (NULL == (regn_hash = HashTableCreate(0, HASH_DYNAMIC_SIZE|HASH_FUNC_JENKINS3))) {
	    return 1;
        }
    
        if( 0 == srf_info(input, level_mode, read_count, chunk_count, key_count, type_count, regn_hash) ){

            /* read counts */
            if( level_mode & LEVEL_READ ) {
                for (iread=0; iread<NREADS; iread++) {
                    if( read_count[iread] )
			printf("Reads: %s : %ld\n", read_str[iread], read_count[iread]);
                }
            }

            /* chunk, key and type counts */
            if( level_mode & LEVEL_CHUNK ) {
                for (ichunk=0; ichunk<NCHUNKS; ichunk++) {
                    if( chunk_count[ichunk] ) {
                        printf("Chunk: %s : %ld\n", ZTR_BE2STR(chunk_type[ichunk], type), chunk_count[ichunk]);
                        for (ikey=0; ikey<NKEYS; ikey++) {
                            if(key_count[ichunk][ikey]) {
                                printf("  Mdata key: %s : %ld\n", keys_str[ikey], key_count[ichunk][ikey]);
                                if (ikey == KEY_TYPE && (ichunk == CHUNK_SAMP || ichunk == CHUNK_SMP4)) {
                                    for (itype=0; itype<NTYPES; itype++)
                                        if(type_count[ichunk][itype])
                                            printf("    types: %s : %ld\n", types_str[itype], type_count[ichunk][itype]);
                                }
                                if (ikey == KEY_NAME && (ichunk == CHUNK_REGN)) {
                                    int ibucket;
                                    for (ibucket=0; ibucket<regn_hash->nbuckets; ibucket++) {
                                        HashItem *hi;
                                        for (hi = regn_hash->bucket[ibucket]; hi; hi = hi->next) {
                                            regn_t *regn = (regn_t *)hi->data.p;
                                            printf("    boundaries: coord=%c name=", regn->coord);
                                            if( regn->name ){
                                                int count = 0;
                                                char *type = strtok (regn->name,";");
                                                while(type) {
                                                    char *cp;
                                                    if(NULL == (cp = strchr(type,':'))) {
                                                        fprintf(stderr, "Invalid region name/code pair %s\n", type);
                                                        return 1;
                                                    }
                                                    *cp++ = '\0';
                                                    if( count ){
                                                        if( count > regn->nbndy ){
							    fprintf(stderr, "More name/code pairs than boundaries\n");
							    return 1;
                                                        }
                                                        printf(" %d ", be_int4(regn->bndy[count-1]));
                                                    }
                                                    printf("%s:%s", type, cp);
                                                    count++;
                                                    type = strtok (NULL, ";");
                                                }
                                            }
                                            printf(" x%d\n", regn->count);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}
