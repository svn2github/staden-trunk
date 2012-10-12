#include <staden_config.h>
#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tg_gio.h"
#include "tg_struct.h"
#include "tg_index_common.h"
#include "sam_index.h"

#include <staden_config.h>

#define _IOLIB 2
#include "depad_seq_tree.h"
#include "sam_pileup.h"
#include "consensus.h"

/*
 * Uncomment this if you want sam auxillary tags to be added as tag in
 * gap5.
 */
//#define SAM_AUX_AS_TAG

typedef struct bio_seq {
    struct bio_seq *next;
    char *seq;
    char *conf;
    int  *pad;
    int seq_len;
    int alloc_len;
    int mqual;
    int pos;
    int left;
    int padded_pos;
    tg_rec rec;
} bio_seq_t;

typedef struct rec_list {
    tg_rec rec;
    struct rec_list *next;
} rec_list_t;

typedef struct {
    GapIO *io;
    bam_file_t *fp;
    char *fn;
    bio_seq_t *seqs;
    bio_seq_t *free_seq;
    int max_seq;
    tg_pair_t *pair;
    HacheTable *libs;
    contig_t *c;
    int c_start;   /* First used based in an existing contig */
    int n_inserts; /* Insertions to seq? */
    int npads;     /* Cons inserts due to new data */
    int npads2;    /* Seq inserts due to old cons pads */
    int count;
    int total_count;
    int skip;
    //bam_header_t *header;
    tg_args *a;
    struct PAD_COUNT *tree; /* re-padding */
    int last_tid;
    void *rg2pl_hash;
} bam_io_t;


/*
 * Converts a string from the read-group PL field to a constant STECH_*
 * value for sequencing technology.
 * str may be null in which case STECH_UNKNOWN is returned.
 */
int stech_str2int(const char *str) {
    if (!str)
	return STECH_UNKNOWN;

    if (strcasecmp(str, "ILLUMINA") == 0)
	return STECH_SOLEXA;
    if (strcasecmp(str, "SOLEXA") == 0)
	return STECH_SOLEXA;

    if (strcasecmp(str, "ABI") == 0)
	return STECH_SANGER;
    if (strcasecmp(str, "CAPILLARY") == 0)
	return STECH_SANGER;
    if (strcasecmp(str, "SANGER") == 0)
	return STECH_SANGER;

    if (strcasecmp(str, "454") == 0)
	return STECH_454;
    if (strcasecmp(str, "LS454") == 0)
	return STECH_454;

    if (strcasecmp(str, "SOLID") == 0)
	return STECH_SOLID;

    return STECH_UNKNOWN;
}

/*
 * Attempts to guess a sequencing technology based on sequence name.
 * This is far from ideal and is also a bit Sanger Institute centric.
 */
int stech_guess_by_name(char *str) {
    size_t l;
    char *cp;
    int colons;

    if (!str || !*str)
	return STECH_UNKNOWN;

    l = strlen(str);

    /* 454 follow a rigid pattern, [0-9A-Z]{14}, first 7 being the plate */
    if (l == 14) {
	int i;
	for (i = 0; i < l; i++) {
	    if (!isalnum(str[i]))
		break;
	}
	if (i == l)
	    return STECH_454;
    }

    /* SOLID appear to start VAB_ (vendor = AB) */
    if (strncmp(str, "VAB_", 4) == 0)
	return STECH_SOLID;

    /*
     * Illumina tend to start with machine-name followed by lane, run,
     * tile, x, y and possibly #index. We look for 4 colons or for
     * machine name IL[0-9]+.
     */
    if (strncmp(str, "IL", 2) == 0 && isdigit(str[2]))
	return STECH_SOLEXA;

    colons = 0;
    cp = str;
    do {
	cp = strchr(cp, ':');
	if (cp) {
	    colons++;
	    cp++;
	}
    } while(cp);

    if (colons == 4) {
	return STECH_SOLEXA;
    }

    
    /*
     * Sanger capillary sequences tend to be template_name.[pq][0-9]k.
     * Very sanger specific, but there's just too much variation.
     */
    cp = strchr(str, '.');
    if (cp) {
	if (cp[1] == 'p' || cp[1] == 'q') {
	    if (isdigit(cp[2])) {
		if (cp[3] == 'k') {
		    return STECH_SANGER;
		}
	    }
	}
    }

    return STECH_UNKNOWN;
}


/*
 * Allocates a new bio_seq_t entry in bam_io_t struct.
 * We use a single indexed array counting from 0 to N-1 representing the
 * N active sequences in a bam_pileup1_t struct. This may mean we have
 * O(D^2) complexity for depth D, so if we find this becomes a bottleneck
 * then we can replace the bio_seq_t array with an ordered linked list
 * instead.
 *
 * Returns the bio_seq_t * on success
 *         NULL on failure
 */
bio_seq_t *bio_new_seq(bam_io_t *bio, pileup_t *p, int pos) {
    int i;
    bio_seq_t *s;
    int seq_offset;

    static struct char2 {
	char c1;
	char c2;
    } lc2[256];
    static int lookup_done = 0;

    if (!lookup_done) {
	int i;
	for (i = 0; i < 256; i++) {
	    lc2[i].c1 = "=ACMGRSVTWYHKDBN"[i >> 4];
	    lc2[i].c2 = "=ACMGRSVTWYHKDBN"[i & 0x0f];
	}

	lookup_done = 1;
    }

    if (bio->free_seq) {
	s = bio->free_seq;
	bio->free_seq = s->next;
    } else {
	if (!(s = malloc(sizeof(*s))))
	    return NULL;
	s->alloc_len = 0;
	s->seq = NULL;
	s->conf = NULL;
	s->pad = NULL;
    }

    s->next = NULL;

    /* Grow storage, first guess, but b->len can be shorter than cigar
     * padded sequence.
     */
    if (s->alloc_len < p->b->len+10)
	s->alloc_len = p->b->len+10;
    
    s->seq = (char *)realloc(s->seq, (int)(s->alloc_len * 1.2));
    s->conf = (char *)realloc(s->conf, (int)(s->alloc_len * 1.2));
    s->pad = (int *)realloc(s->pad, (int)(s->alloc_len * 1.2)*sizeof(int));
    s->padded_pos = 0;

    //    printf("New seq %p / %p => %p,%p\n", p, s, s->seq, s->conf);

    seq_offset = p->seq_offset;
    if (p->first_del)
	seq_offset++;


    /* left hand cutoff data */
    if (p->ref_skip) {
	/* 2nd or more portion of a sequence with 'N' cigar ops */
	i = 0;
	s->seq_len = 0;
    } else if (seq_offset >= 0) {
	unsigned char *seq  = (unsigned char *)bam_seq(p->b);
	unsigned char *qual = (unsigned char *)bam_qual(p->b);
	char *sp = s->seq;
	char *qp = s->conf;

        if (seq_offset >= s->alloc_len) {
            s->alloc_len = (seq_offset + 10);
            if (NULL == (s->seq  = (char *)realloc(s->seq,  s->alloc_len)))
                return NULL;
            if (NULL == (s->conf = (char *)realloc(s->conf, s->alloc_len)))
                return NULL;
            if (NULL == (s->pad = (int *)realloc(s->pad,
						 s->alloc_len*sizeof(int))))
                return NULL;

	    sp = s->seq;
	    qp = s->conf;
        }

	if ((bio->a->data_type & DATA_SEQ) && (p->b->len >= seq_offset)) {
	    for (i = 0; i < seq_offset; i+=2) {
		unsigned char coded = *seq++;
		*sp++ = lc2[coded].c1;
		*sp++ = lc2[coded].c2;
	    }
	} else {
	    for (i = 0; i < seq_offset; i++)
		*sp++ = 'N';
	}

	for (i = 0; i < seq_offset; i++)
	    s->pad[i] = i;
	s->padded_pos = seq_offset;

	if ((bio->a->data_type & DATA_QUAL) && (p->b->len >= seq_offset)) {
	    for (i = 0; i < seq_offset; i++)
		*qp++ = *qual++;
	} else {
	    for (i = 0; i < seq_offset; i++)
		*qp++ = 0;
	}

	i = seq_offset;
	s->seq_len = i;

	/* Equiv to:
	 * for (i = 0; i < p->seq_offset; i++) {
	 *    s->seq [i] = bam_nt16_rev_table[bam_seqi(seq, i)];
	 *    s->conf[i] = qual[i];
	 * }
	 */
    } else if (seq_offset == 0) {
	i = 0;
    } else {
	/*
	 * Either the input data is unsorted (caught elsewhere)
	 * or the CIGAR string starts with a deletion.
	 */
	i = 0;
    }

    s->seq_len = i;
    s->pos = pos - i;
    s->left = i;

    if (bio->a->data_type == DATA_BLANK) {
	static int fake_recno = 1;
	s->rec = fake_recno++;
    } else {
	s->rec = sequence_new_from(bio->io, NULL);
    }

    return s;
}


static char *append_int(char *cp, int i) {
    int j, k = 0;

    if (i < 0) {
	*cp++ = '-';
	i = -i;
    } else if (i == 0) {
	*cp++ = '0';
	return cp;
    }

    if (i < 1000)
	goto b1;
    if (i < 100000)
	goto b2;
    if (i < 100000000)
	goto b3;

    j = i / 1000000000;
    if (j || k) *cp++ = j + '0', k=1, i %= 1000000000;

    j = i / 100000000;
    if (j || k) *cp++ = j + '0', k=1, i %= 100000000;
    
 b3:
    j = i / 10000000;
    if (j || k) *cp++ = j + '0', k=1, i %= 10000000;
    
    j = i / 1000000;
    if (j || k) *cp++ = j + '0', k=1, i %= 1000000;
    
    j = i / 100000;
    if (j || k) *cp++ = j + '0', k=1, i %= 100000;
    
 b2:
    j = i / 10000;
    if (j || k) *cp++ = j + '0', k=1, i %= 10000;

    j = i / 1000;
    if (j || k) *cp++ = j + '0', k=1, i %= 1000;

 b1:
    j = i / 100;
    if (j || k) *cp++ = j + '0', k=1, i %= 100;

    j = i / 10;
    if (j || k) *cp++ = j + '0', k=1, i %= 10;

    if (i || k) *cp++ = i + '0';

    return cp;
}

#define APPEND_FMT(fmt) \
    *cp++ = s[0]; \
    *cp++ = s[1]; \
    *cp++ = ':'; \
    *cp++ = (fmt); \
    *cp++ = ':';

/* Experiments to speed up printf */
#if 0

/* A basic printf without a lot of the %-10 or %5.5 etc field width
 * specifiers.
 * It is, however, much faster.
 */
#define EVEN_SIMPLER
int fast_fprintf(FILE *fp, char *fmt, ...) {
    char *cp, c;
    long l;
    int i;
    double d; 
    char tmp[128], *ct;
    va_list ap;

    va_start(ap, fmt);

    for (cp = fmt; *cp; cp++) {
	switch(*cp) {

	/* A format specifier */
	case '%': {
	    char *endp;
	    long conv_len1=0, conv_len2=0, conv_len=0;
	    signed int arg_size;

	    cp++;
#ifndef EVEN_SIMPLER
	    /* Firstly, strip the modifier flags (+-#0 and [space]) */
	    for(; c=*cp;) {
		if ('#' == c)
		    ;
		else if ('-' == c || '+' == c || ' ' == c)
		    ;
		else
		    break;
	    }

	    /* Width specifier */
	    if (*cp >= '0' && *cp <= '9') {
		for (l = 0; *cp >= '0' && *cp <= '9'; cp++)
		    l = l*10 + *cp - '0';
		conv_len = conv_len1 = l;
	    } else if (*cp == '*') {
		conv_len = conv_len1 = (int)va_arg(ap, int);
		cp++;
	    }
#endif

	    /* Precision specifier */
	    if ('.' == *cp) {
		cp++;
		for (l = 0; *cp >= '0' && *cp <= '9'; cp++)
		    l = l*10 + *cp - '0';
		conv_len2 = l;
		if (*cp == '*') {
		    conv_len2 = (int)va_arg(ap, int);
		    cp++;
		}
		conv_len = MAX(conv_len1, conv_len2);
	    }

#ifndef EVEN_SIMPLER
	    /* Short/long identifier */
	    if ('h' == *cp) {
		arg_size = -1; /* short */
		cp++;
	    } else if ('l' == *cp) {
		arg_size = 1; /* long */
		cp++;
	    } else {
		arg_size = 0; /* int */
	    }
#endif

	    /* The actual type */
	    switch (*cp) {
	    case '%':
		/*
		 * Not real ANSI I suspect, but we'll allow for the
		 * completely daft "%10%" example.
		 */
		putc('%', fp);
		break;

	    case 'd':
	    case 'i':
	    case 'u':
		/* Remember: char and short are sent as int on the stack */
		if (arg_size == -1)
		    l = (long)va_arg(ap, int);
		else if (arg_size == 1)
		    l = va_arg(ap, long); 
		else 
		    l = (long)va_arg(ap, int);

		ct = append_int(tmp, (int)l);
		*ct++ = 0;
		fputs(tmp, fp);
		break;

	    case 'c':
		i = va_arg(ap, int);
		putc(i, fp);
		break;

	    case 'f':
		d = va_arg(ap, double);
		fprintf(fp, "%f", d);
		break;

	    case 'e':
	    case 'E':
	    case 'g':
	    case 'G':
		d = va_arg(ap, double);
		fprintf(fp, "%g", d);
		break;

	    case 'p':
	    case 'x':
	    case 'X':
		l = (long)va_arg(ap, void *);
		puts("TODO");
		break;

	    case 'n':
		/* produces no output */
		break;

	    case 's': {
		char *s = (char *)va_arg(ap, char *);
		if (conv_len2) {
		    fwrite(s, conv_len2, 1, fp);
		} else
		    fputs(s, fp);
		break;
	    }

	    default:
		/* wchar_t types of 'C' and 'S' aren't supported */
		puts("TODO");
	    }
	    
	}

	case '\0':
	    break;

	default:
	    putc(*cp, fp);
	}
    }

    va_end(ap);
    
    return 0;
}

/* Max 8k of output */
int faster_fprintf(FILE *fp, char *fmt, ...) {
    char *cp, c;
    long l;
    int i;
    double d; 
    char tmp[128], *ct;
    va_list ap;
    char max_line[8192], *out = max_line;

    va_start(ap, fmt);

    for (cp = fmt; *cp; cp++) {
	switch(*cp) {

	/* A format specifier */
	case '%': {
	    char *endp;
	    long conv_len1=0, conv_len2=0, conv_len=0;
	    signed int arg_size;

	    cp++;

#ifndef EVEN_SIMPLER
	    /* Firstly, strip the modifier flags (+-#0 and [space]) */
	    for(; c=*cp;) {
		if ('#' == c)
		    ;
		else if ('-' == c || '+' == c || ' ' == c)
		    ;
		else
		    break;
	    }

	    /* Width specifier */
	    if (*cp >= '0' && *cp <= '9') {
		for (l = 0; *cp >= '0' && *cp <= '9'; cp++)
		    l = l*10 + *cp - '0';
		conv_len = conv_len1 = l;
	    } else if (*cp == '*') {
		conv_len = conv_len1 = (int)va_arg(ap, int);
		cp++;
	    }
#endif

	    /* Precision specifier */
	    if ('.' == *cp) {
		cp++;
		for (l = 0; *cp >= '0' && *cp <= '9'; cp++)
		    l = l*10 + *cp - '0';
		conv_len2 = l;
		if (*cp == '*') {
		    conv_len2 = (int)va_arg(ap, int);
		    cp++;
		}
		conv_len = MAX(conv_len1, conv_len2);
	    }

#ifndef EVEN_SIMPLER
	    /* Short/long identifier */
	    if ('h' == *cp) {
		arg_size = -1; /* short */
		cp++;
	    } else if ('l' == *cp) {
		arg_size = 1; /* long */
		cp++;
	    } else {
		arg_size = 0; /* int */
	    }
#endif

	    /* The actual type */
	    switch (*cp) {
	    case '%':
		/*
		 * Not real ANSI I suspect, but we'll allow for the
		 * completely daft "%10%" example.
		 */
		*out++ = '%';
		break;

	    case 'd':
	    case 'i':
	    case 'u':
		/* Remember: char and short are sent as int on the stack */
		if (arg_size == -1)
		    l = (long)va_arg(ap, int);
		else if (arg_size == 1)
		    l = va_arg(ap, long); 
		else 
		    l = (long)va_arg(ap, int);

		out = append_int(out, (int)l);
		break;

	    case 'c':
		i = va_arg(ap, int);
		*out++ = i;
		break;

	    case 'f':
		d = va_arg(ap, double);
		out += sprintf(out, "%f", d);
		break;

	    case 'e':
	    case 'E':
	    case 'g':
	    case 'G':
		d = va_arg(ap, double);
		out += sprintf(out, "%g", d);
		break;

	    case 'p':
	    case 'x':
	    case 'X':
		l = (long)va_arg(ap, void *);
		puts("TODO");
		break;

	    case 'n':
		/* produces no output */
		break;

	    case 's': {
		char *s = (char *)va_arg(ap, char *);
		if (conv_len2) {
		    //memcpy(out, s, conv_len2);
		    //out += conv_len2;

		    //strncpy(out, s, conv_len2);
		    //out += MIN(conv_len2, strnlen(s));

		    while(conv_len2 && *s) {
		        *out++ = *s++;
		        conv_len2--;
		    }
		} else {
		    //strcpy(out, s);
		    //out += strlen(s);

		    while(*s) *out++ = *s++;
		}
		//		if (conv_len2) {
		//		    fwrite(s, conv_len2, 1, fp);
		//		} else
		//		    fputs(s, fp);
		break;
	    }

	    default:
		/* wchar_t types of 'C' and 'S' aren't supported */
		puts("TODO");
	    }
	    
	}

	case '\0':
	    break;

	default:
	    *out++ = *cp;
	}
    }

    *out = 0;
    fwrite(max_line, 1, out-max_line, fp);

    va_end(ap);
    
    return 0;
}
#endif

char *sam_aux_stringify(char *s, int len) {
    static char str[8192];
    char *cp = str, *s_end = s+len;
    int first = 1;

    //write(2, s, (int)(((char *)&b->ref) + b->blk_size - (uint8_t *)s));

    while (s < s_end) {
	if (first)
	    first = 0;
	else
	    *cp++ = '\t';

	switch (s[2]) {
	case 'A':
	    APPEND_FMT('A');
	    *cp++ = s[3];
	    s+=4;
	    break;

	case 'C':
	    APPEND_FMT('i');
	    cp = append_int(cp, *(uint8_t *)(s+3));
	    s+=4;
	    break;

	case 'c':
	    APPEND_FMT('i');
	    cp = append_int(cp, *(int8_t *)(s+3));
	    s+=4;
	    break;

	case 'S':
	    {
		uint16_t tmp;
		((char *)&tmp)[0] = s[3];
		((char *)&tmp)[1] = s[4];
		APPEND_FMT('i');
		cp = append_int(cp, tmp);
	    }
	    s+=5;
	    break;

	case 's':
	    {
		int16_t tmp;
		((char *)&tmp)[0] = s[3];
		((char *)&tmp)[1] = s[4];
		APPEND_FMT('i');
		cp = append_int(cp, tmp);
	    }
	    s+=5;
	    break;

	case 'I':
	    {
		uint32_t tmp;
		((char *)&tmp)[0] = s[3];
		((char *)&tmp)[1] = s[4];
		((char *)&tmp)[2] = s[5];
		((char *)&tmp)[3] = s[6];
		APPEND_FMT('i');
		cp = append_int(cp, tmp);
	    }
	    s+=7;
	    break;

	case 'i':
	    {
		int32_t tmp;
		((char *)&tmp)[0] = s[3];
		((char *)&tmp)[1] = s[4];
		((char *)&tmp)[2] = s[5];
		((char *)&tmp)[3] = s[6];
		APPEND_FMT('i');
		cp = append_int(cp, tmp);
	    }
	    s+=7;
	    break;

	case 'f':
	    {
		float f;
		memcpy(&f, s+3, 4);
		APPEND_FMT('f');
		cp += sprintf(cp, "%f", f);
	    }
	    s+=7;
	    break;

	case 'd':
	    {
		double d;
		memcpy(&d, s+3, 8);
		APPEND_FMT('f');
		cp += sprintf(cp, "%f", d);
	    }
	    s+=11;
	    break;

	case 'Z': case 'H':
	    APPEND_FMT(s[2]);
	    s+=3;
	    while (*s)
		*cp++ = *s++;
	    s++;
	    break;

	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	    return NULL;
	}
    }

    *cp = 0;

    return str;
}

char *bam_aux_stringify(bam_seq_t *b, int no_RG) {
    static char str[8192];
    char *s = (char *)bam_aux(b), *cp = str;
    int first = 1;
    int keep;

    no_RG = 1;

    //write(2, s, (int)(b->data + b->data_len - (uint8_t *)s));

    while ((uint8_t *)s < ((uint8_t *)&b->ref) + b->blk_size) {

	keep = (no_RG && s[0] == 'R' && s[1] == 'G') ? 0 : 1;
	if (keep) {
	    if (first)
		first = 0;
	    else
		*cp++ = '\t';
	}

	switch (s[2]) {
	case 'A':
	    if (keep) {
		APPEND_FMT('A');
		*cp++ = s[3];
	    }
	    s+=4;
	    break;

	case 'C':
	    if (keep) {
		APPEND_FMT('i');
		cp = append_int(cp, *(uint8_t *)(s+3));
	    }
	    s+=4;
	    break;

	case 'c':
	    if (keep) {
		APPEND_FMT('i');
		cp = append_int(cp, *(int8_t *)(s+3));
	    }
	    s+=4;
	    break;

	case 'S':
	    if (keep) {
		uint16_t tmp;
		((char *)&tmp)[0] = s[3];
		((char *)&tmp)[1] = s[4];
		APPEND_FMT('i');
		cp = append_int(cp, tmp);
	    }
	    s+=5;
	    break;

	case 's':
	    if (keep) {
		int16_t tmp;
		((char *)&tmp)[0] = s[3];
		((char *)&tmp)[1] = s[4];
		APPEND_FMT('i');
		cp = append_int(cp, tmp);
	    }
	    s+=5;
	    break;

	case 'I':
	    if (keep) {
		uint32_t tmp;
		((char *)&tmp)[0] = s[3];
		((char *)&tmp)[1] = s[4];
		((char *)&tmp)[2] = s[5];
		((char *)&tmp)[3] = s[6];
		APPEND_FMT('i');
		cp = append_int(cp, tmp);
	    }
	    s+=7;
	    break;

	case 'i':
	    if (keep) {
		int32_t tmp;
		((char *)&tmp)[0] = s[3];
		((char *)&tmp)[1] = s[4];
		((char *)&tmp)[2] = s[5];
		((char *)&tmp)[3] = s[6];
		APPEND_FMT('i');
		cp = append_int(cp, tmp);
	    }
	    s+=7;
	    break;

	case 'f':
	    if (keep) {
		float f;
		memcpy(&f, s+3, 4);
		APPEND_FMT('f');
		cp += sprintf(cp, "%f", f);
	    }
	    s+=7;
	    break;

	case 'd':
	    if (keep) {
		double d;
		memcpy(&d, s+3, 8);
		APPEND_FMT('f');
		cp += sprintf(cp, "%f", d);
	    }
	    s+=11;
	    break;

	case 'Z': case 'H':
	    if (keep) {
		APPEND_FMT(s[2]);
		s+=3;
		while (*s)
		    *cp++ = *s++;
		s -= 3;
	    }
	    s+=3;
	    if (!keep)
		while (*s++);
	    break;

	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	    return NULL;
	}
    }

    *cp = 0;

    return str;
}

/*
 * Filters out specific types from the sam aux records.
 * Returns a new aux record, also in the compact binary form.
 * 'len' holds the size of the returned string.
 *
 * Value returned is statically allocated. Do not free.
 */
char *bam_aux_filter(bam_seq_t *b, char **types, int ntypes, int *len) {
    static char str[8192];
    char *s = (char *)bam_aux(b), *cp = str;
    int keep, i;

    while ((uint8_t *)s < ((uint8_t *)&b->ref) + b->blk_size) {
	keep = 1;
	for (i = 0; i < ntypes; i++) {
	    if (s[0] == types[i][0] &&
		s[1] == types[i][1]) {
		keep = 0;
		break;
	    }
	}

	if (keep) {
	    *cp++ = s[0];
	    *cp++ = s[1];
	    *cp++ = s[2];
	}

	switch (s[2]) {
	case 'A':
	case 'C':
	case 'c':
	    if (keep)
		*cp++ = s[3];
	    s+=4;
	    break;

	case 'S':
	case 's':
	    if (keep) {
		*cp++ = s[3];
		*cp++ = s[4];
	    }
	    s+=5;
	    break;

	case 'I':
	case 'i':
	case 'f':
	    if (keep) {
		*cp++ = s[3];
		*cp++ = s[4];
		*cp++ = s[5];
		*cp++ = s[6];
	    }
	    s+=7;
	    break;

	case 'd':
	    if (keep) {
		*cp++ = s[3];
		*cp++ = s[4];
		*cp++ = s[5];
		*cp++ = s[6];
		*cp++ = s[7];
		*cp++ = s[8];
		*cp++ = s[9];
		*cp++ = s[10];
	    }
	    s+=11;
	    break;

	case 'Z': case 'H':
	    s+=3;
	    if (keep)
		while ((*cp++ = *s++));
	    else
		while (*s++);
	    break;

	default:
	    fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	    return NULL;
	}
    }

    *len = cp - str;

    return str;
}

/*
 * Creates a new contig and updates the bam_io_t struct.
 */
void bio_new_contig(bam_io_t *bio, int tid) {
    char *cname = bio->fp->ref[tid].name;

    printf("\n++Processing contig %d / %s\n", tid, cname);
	
    create_new_contig(bio->io, &(bio->c), cname, bio->a->merge_contigs);
    bio->n_inserts = 0;
    bio->npads = 0;
    bio->npads2 = 0;
    bio->skip = 0;

    if (bio->a->repad) {
	bio->tree = depad_consensus(bio->io, bio->c->rec);
	//padtree_dump(bio->tree);
	consensus_valid_range(bio->io, bio->c->rec, &bio->c_start, NULL);
    }
	
    bio->last_tid = tid;
}

/*
 * Samtools pileup won't iterate over unmapped reads. Therefore we have a
 * separate function to add these to the database - this one.
 * Although it shares much of the same code so is a candidate for merging
 * at some stage.
 */
int bio_add_unmapped(bam_io_t *bio, bam_seq_t *b) {
    char *LB;
    HacheItem *hi;
    HacheData hd;
    seq_t s;
    char tname[1024];
    library_t *lib = NULL;
    int new = 0;
    char *name, *suffix;
    int name_len;
    char *aux;
    int i, flags;
    tg_rec recno, bin_rec;
    int paired, is_pair = 0;
    char *filter[] = {"RG"};
    int stech;

    bio->count++;

    /* Check if it's a new contig, create if so */
    if (b->ref != bio->last_tid) {
	bio_new_contig(bio, b->ref);
    }

    /* Fetch read-group and pretend it's a library for now */
    if ((LB = bam_aux_find(b, "RG"))) {
	tag_list_t *rg_tag = bam_find_rg(bio->fp, LB);
	stech = rg_tag ? stech_str2int(rg_tag->value) : STECH_UNKNOWN;
    } else {
	LB = bio->fn;
	stech = STECH_UNKNOWN;
    }

    suffix = bam_aux_find(b, "FS");

    hd.p = NULL;
    hi = HacheTableAdd(bio->libs, (char *)LB, strlen(LB), hd, &new);
    if (new) {
	tg_rec lrec;
	printf("New library %s\n", LB);

	lrec = library_new(bio->io, (char *)LB);
	lib = get_lib(bio->io, lrec);
	lib = cache_rw(bio->io, lib);
	lib->machine = stech;
	hi->data.p = lib;
	cache_incr(bio->io, lib);
    }
    lib = hi->data.p;

    /* Construct a seq_t struct */
    name = bam_name(b);
    name_len = strlen(name);

    /*
     * Uncomment one of the two sections below if you wish to allow storing
     * sam auxillary records in gap5.
     * This is experimental and the format for these hasn't been finalised
     * yet.
     */
    aux = NULL;
    s.aux_len = 0;

    if (bio->a->sam_aux)
	aux = bam_aux_filter(b, filter, 1, &s.aux_len);

    //aux = bam_aux_stringify(b, 1);
    //s.aux_len = strlen(aux);
    
    //aux = bam1_aux(b);
    //s.aux_len = (int)(b->data + b->data_len - (uint8_t *)aux);

    s.pos = bio->npads +
	get_padded_coord(bio->tree, b->pos + 1 + bio->n_inserts
			 - bio->npads) + bio->c_start-1;
    //s.pos = b->pos+1;
    s.len = b->len;
    s.rec = 0;
    s.seq_tech = stech != STECH_UNKNOWN ? stech : stech_guess_by_name(name);
    s.flags = 0;
    s.left  = 1;
    s.right = s.len;
    s.parent_type = 0;
    s.parent_rec = 0;
    if (bio->a->data_type & DATA_NAME) {
	s.template_name_len = name_len;
	s.name_len = name_len + (suffix ? strlen(suffix) : 0);;
	s.name = (char *)malloc(s.name_len + 3 + 2*s.len + s.aux_len);
	strcpy(s.name, name);
	if (suffix)
	    strcat(s.name, suffix);
    } else {
	char *n = "";
	s.name_len = 0;
	s.template_name_len = 0;
	s.name = (char *)malloc(s.name_len + 3 + 2*s.len + s.aux_len);
	strcpy(s.name, n);
    }
    s.trace_name = s.name + s.name_len + 1;
    *s.trace_name = 0;
    s.trace_name_len = 0;
    s.alignment = s.trace_name + s.trace_name_len + 1;
    *s.alignment = 0;
    s.alignment_len = 0;
    s.seq = s.alignment + s.alignment_len+1;
    s.conf = s.seq+s.len;
    s.mapping_qual = bam_map_qual(b);
    s.format = SEQ_FORMAT_MAQ; /* pack bytes */
    s.anno = NULL;
    s.sam_aux = s.conf + s.len;

    for (i = 0; i < b->len; i++) {
	s.seq[i] = bio->a->data_type & DATA_SEQ
	    ? bam_nt16_rev_table[bam_seqi(bam_seq(b), i)]
	    : 'N';
	s.conf[i] = bio->a->data_type & DATA_QUAL
	    ? bam_qual(b)[i]
	    : 0;
    }

    if (bam_strand(b)) {
	complement_seq_t(&s);
	s.flags |= SEQ_COMPLEMENTED;
    }

    if (aux)
	memcpy(s.sam_aux, aux, s.aux_len);

    /* Create the range, save the sequence */
    paired = (bam_flag(b) & BAM_FPAIRED) ? 1 : 0;
    flags = paired ? GRANGE_FLAG_TYPE_PAIRED : GRANGE_FLAG_TYPE_SINGLE;

    if (bam_flag(b) & BAM_FREAD1)
	s.flags |= SEQ_END_FWD;

    if (bam_flag(b) & BAM_FREAD2)
	s.flags |= SEQ_END_REV;

    strcpy(tname, name);

    if (!suffix && name_len >= 2 && name[name_len-2] == '/') {
	tname[name_len-2] = 0;

	/* Check validity of name vs bit-fields */
	if ((name[name_len-1] == '1' &&
	     (s.flags & SEQ_END_MASK) != SEQ_END_FWD) ||
	    (name[name_len-1] == '2' &&
	     (s.flags & SEQ_END_MASK) != SEQ_END_REV)) {
	    fprintf(stderr, "Inconsistent read name vs flags: %s vs 0x%02x\n",
		    name, bam_flag(b));
	}
    }

    if (paired)
	flags |= (s.flags & SEQ_END_MASK) == SEQ_END_FWD
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;
    else
	/* Guess work here. For now all <--- are rev, all ---> are fwd */
	flags |= bam_strand(b)
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;

    if (bam_strand(b)) {
	flags |= GRANGE_FLAG_COMP1;
    }

    if (bio->pair) is_pair = 1;

    flags |= GRANGE_FLAG_ISUMSEQ;

    recno = save_range_sequence(bio->io, &s, s.mapping_qual, bio->pair,
    				is_pair, tname, bio->c, bio->a, flags, lib,
				&bin_rec);


#ifdef SAM_AUX_AS_TAG
    /* Make an annotation out of the sam auxillary data */
    aux = bam_aux_stringify(b, 1);
    if (aux && *aux) {
	anno_ele_t *e;
	bin_index_t *bin;
	range_t r;

	r.mqual = str2type("SAMX");
	r.start = s.pos;
	r.end = s.pos;
	r.pair_rec = recno;
	r.flags = GRANGE_FLAG_ISANNO | GRANGE_FLAG_TAG_SEQ;
	r.rec = anno_ele_new(bio->io, 0, GT_Seq, recno, 0, r.mqual,
			     ANNO_DIR_NUL, aux);
	e = (anno_ele_t *)cache_search(bio->io, GT_AnnoEle, r.rec);
	e = cache_rw(bio->io, e);
	
	bin = bin_add_to_range(bio->io, &bio->c, bin_rec, &r, NULL, NULL, 0);
	e->bin = bin->rec;
    }
#endif

    return 0;
}

static int hex[256] = {
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     0, 1, 2, 3, 4, 5, 6, 7, 8, 9,-1,-1,-1,-1,-1,-1,
    -1,10,11,12,13,14,15,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,10,11,12,13,14,15,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

/*
 * Parses a tag string and splits it into separate components returned in
 * start, end, dir, type and text. The tag is of the form:
 *   start;end;strand;type(;key=value)*
 * Right now we just treat the (;key=value)* bit as arbitrary text.
 */
static char *parse_bam_PT_tag(char *str, int *start, int *end, char *dir,
			      char **type, int *type_len,
			      char **text, int *text_len) {
    char *cp, *orig = str;
    int n;

    if (!*str)
	return NULL;

    if (3 != sscanf(str, "%d;%d;%c;%n", start, end, dir, &n))
	goto error;
    str += n;

    *type = cp = str;
    while (*str && *str != ';' && *str != '|') {
	if (*str == '%' && isxdigit(str[1]) && isxdigit(str[2])) {
	    *cp++ = (hex[str[1]]<<4) | hex[str[2]];
	    str += 3;
	} else {
	    *cp++ = *str++;
	}
    }
    *type_len = cp-*type;

    switch (*str) {
    case ';':
	str++;
	break;

    case '|':
	*text_len = 0;
	*text = NULL;
	return ++str;

    case '\0':
	*text_len = 0;
	*text = NULL;
	return str;
    }

    *text = cp = str;
    while (*str && *str != '|') {
	if (*str == '%' && isxdigit(str[1]) && isxdigit(str[2])) {
	    *cp++ = (hex[str[1]]<<4) | hex[str[2]];
	    str += 3;
	} else {
	    *cp++ = *str++;
	}
    }

    /* Gap5 doesn't support full key=value GFF syntax yet, so cheat. */
    if (0 == strncmp(*text, "Note=", 5))
	(*text) += 5;

    *text_len = cp-*text;

    return *str == '|' ? ++str : str;

 error:
    verror(ERR_WARN, "parse_bam_PT_tag", "invalid sam/bam annotation: %s",
	   orig);
    return NULL;
}

/*
 * Parses a tag string and splits it into separate components returned in
 * start, end, dir, type and text. The tag is of the form:
 * strand|type|text
 */
static char *parse_bam_CT_tag(char *str, char *dir,
			      char **type, int *type_len,
			      char **text, int *text_len) {
    char *cp, *orig = str;
    int n;

    if (!*str)
	return NULL;

    *dir = *str++;
    if (! (*str && *str == ';'))
	goto error;
    str++;

    *type = cp = str;
    while (*str && *str != ';') {
	if (*str == '%' && isxdigit(str[1]) && isxdigit(str[2])) {
	    *cp++ = (hex[str[1]]<<4) | hex[str[2]];
	    str += 3;
	} else {
	    *cp++ = *str++;
	}
    }
    *type_len = cp-*type;

    switch (*str) {
    case ';':
	str++;
	break;

    case '|':
	*text_len = 0;
	*text = NULL;
	return ++str;

    case '\0':
	*text_len = 0;
	*text = NULL;
	return str;
    }

    *text = cp = str;
    while (*str && *str != '|') {
	if (*str == '%' && isxdigit(str[1]) && isxdigit(str[2])) {
	    *cp++ = (hex[str[1]]<<4) | hex[str[2]];
	    str += 3;
	} else {
	    *cp++ = *str++;
	}
    }

    /* Gap5 doesn't support full key=value GFF syntax yet, so cheat. */
    if (0 == strncmp(*text, "Note=", 5))
	(*text) += 5;

    *text_len = cp-*text;

    return *str == '|' ? ++str : str;

 error:
    verror(ERR_WARN, "parse_bam_CT_tag", "invalid sam/bam annotation: %s",
	   orig);
    return NULL;
}

/*
 * Removes a sequence from the bam_io_t struct.
 * This actually performs the main work of adding a sequence to the gap5
 * database.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int bio_del_seq(bam_io_t *bio, pileup_t *p) {
    bio_seq_t *bs = (bio_seq_t *)p->cd;
    bam_seq_t *b;
    seq_t s;
    HacheItem *hi;
    tg_rec recno, bin_rec;
    int i, paired;
    int is_pair = 0;
    int flags;
    char tname[1024];
    library_t *lib = NULL;
    char type;
    char *LB;
    HacheData hd;
    int new = 0;
    char *name, *suffix;
    int name_len;
    char *aux;
    char *filter[] = {"RG"};
    char *handle, aux_key[2];
    int stech;
    bam_aux_t val;
    char *tags;
    int fake;

    bio->count++;

    b = p->b;
    fake = ((bam_flag(b) & BAM_FSECONDARY) &&
	    (bam_flag(b) & BAM_FQCFAIL));

    if (fake)
	goto anno_only; /* Yes I know! The code needs splitting up */

    /* Fetch read-group and pretend it's a library for now */
    if ((LB = bam_aux_find(b, "RG"))) {
	tag_list_t *rg_tag = bam_find_rg(bio->fp, LB);
	stech = rg_tag ? stech_str2int(rg_tag->value) : STECH_UNKNOWN;
    } else {
	LB = bio->fn;
	stech = STECH_UNKNOWN;
    }

    suffix = bam_aux_find(b, "FS");

    hd.p = NULL;
    hi = HacheTableAdd(bio->libs, (char *)LB, strlen(LB), hd, &new);
    if (new) {
	tg_rec lrec;
	printf("New library %s\n", LB);

	lrec = library_new(bio->io, (char *)LB);
	lib = get_lib(bio->io, lrec);
	lib = cache_rw(bio->io, lib);
	lib->machine = stech;
	hi->data.p = lib;
	cache_incr(bio->io, lib);
    }
    lib = hi->data.p;

    /*
    printf("\nSeq %d @ %6d: '%.*s' '%.*s' => nseq=%d\n",
	   snum, bs->pos, bs->seq_len, bs->seq, bs->seq_len, bs->conf,
	   bio->nseq-1);
    */

    /* Construct a seq_t struct */
    s.right = bs->seq_len;
    if (p->seq_offset+1 < b->len && p->ref_skip == 0) {
	unsigned char *b_seq  = (unsigned char *)bam_seq(p->b);
	unsigned char *b_qual = (unsigned char *)bam_qual(p->b);

	if (bs->seq_len+b->len - (p->seq_offset+1) >= bs->alloc_len) {
	    bs->alloc_len = bs->seq_len+b->len - (p->seq_offset+1) + 1;
	    if (NULL == (bs->seq  = (char *)realloc(bs->seq,  bs->alloc_len)))
		return -1;
	    if (NULL == (bs->conf = (char *)realloc(bs->conf, bs->alloc_len)))
		return -1;
	    if (NULL == (bs->pad = (int *)realloc(bs->pad,
						  bs->alloc_len*sizeof(int))))
		return -1;
	}
	for (i = p->seq_offset+1; i < b->len; i++) {
	    bs->seq [bs->seq_len] = bam_nt16_rev_table[bam_seqi(b_seq,i)];
	    bs->conf[bs->seq_len] = b_qual[i];
	    bs->pad[bs->seq_len] = bs->padded_pos++;
	    bs->seq_len++;
	}
    }
    
    name = bam_name(b);
    name_len = strlen(name);

    /*
     * Uncomment one of the two sections below if you wish to allow storing
     * sam auxillary records in gap5.
     * This is experimental and the format for these hasn't been finalised
     * yet.
     */
    aux = NULL;
    s.aux_len = 0;

    if (bio->a->sam_aux)
	aux = bam_aux_filter(b, filter, 1, &s.aux_len);
    
    //aux = bam_aux_stringify(b, 1);
    //s.aux_len = strlen(aux);

    //aux = bam_aux(b);
    //s.aux_len = (int)(b->data + b->data_len - (uint8_t *)aux);
    
    s.rec = bs->rec;
    s.pos = bs->pos;
    s.len = bs->seq_len;
    s.seq_tech = stech != STECH_UNKNOWN ? stech : stech_guess_by_name(name);
    s.flags = 0;
    s.left = bs->left+1;
    s.parent_type = 0;
    s.parent_rec = 0;
    if (bio->a->data_type & DATA_NAME) {
	s.template_name_len = name_len;
	s.name_len = name_len + (suffix ? strlen(suffix) : 0);
	s.name = (char *)malloc(s.name_len + 3 + 2*s.len + s.aux_len);
	strcpy(s.name, name);
	if (suffix)
	    strcat(s.name, suffix);
    } else {
	char *n = "";
	s.name_len = 0;
	s.template_name_len = 0;
	s.name = (char *)malloc(s.name_len + 3 + 2*s.len + s.aux_len);
	strcpy(s.name, n);
    }
    s.trace_name = s.name + s.name_len + 1;
    *s.trace_name = 0;
    s.trace_name_len = 0;
    s.alignment = s.trace_name + s.trace_name_len + 1;
    *s.alignment = 0;
    s.alignment_len = 0;
    s.seq = s.alignment + s.alignment_len+1;
    s.conf = s.seq+s.len;
    s.mapping_qual = bam_map_qual(b);
    s.format = SEQ_FORMAT_MAQ; /* pack bytes */
    s.anno = NULL;
    s.sam_aux = s.conf + s.len;

    if (bio->a->data_type & DATA_SEQ)
	memcpy(s.seq,  bs->seq,  s.len);
    else
	memset(s.seq, 'N', s.len);

    if (bio->a->data_type & DATA_QUAL)
	memcpy(s.conf, bs->conf, s.len);
    else
	memset(s.conf, 0, s.len);
    
    if (bam_strand(b)) {
	complement_seq_t(&s);
	s.flags |= SEQ_COMPLEMENTED;
    }

    if (aux)
	memcpy(s.sam_aux, aux, s.aux_len);

    /* Create the range, save the sequence */
    paired = (bam_flag(b) & BAM_FPAIRED) ? 1 : 0;
    flags = paired ? GRANGE_FLAG_TYPE_PAIRED : GRANGE_FLAG_TYPE_SINGLE;

    if (bam_flag(b) & BAM_FREAD1)
	s.flags |= SEQ_END_FWD;

    if (bam_flag(b) & BAM_FREAD2)
	s.flags |= SEQ_END_REV;

    strcpy(tname, name);

    if (!suffix && name_len >= 2 && name[name_len-2] == '/') {
	tname[name_len-2] = 0;

	/* Check validity of name vs bit-fields */
	if ((name[name_len-1] == '1' &&
	     (s.flags & SEQ_END_MASK) != SEQ_END_FWD) ||
	    (name[name_len-1] == '2' &&
	     (s.flags & SEQ_END_MASK) != SEQ_END_REV)) {
	    fprintf(stderr, "Inconsistent read name vs flags: %s vs 0x%02x\n",
		    name, bam_flag(b));
	}
    }

    if (paired)
	flags |= (s.flags & SEQ_END_MASK) == SEQ_END_FWD
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;
    else
	/* Guess work here. For now all <--- are rev, all ---> are fwd */
	flags |= bam_strand(b)
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;

    if (bam_strand(b)) {
	flags |= GRANGE_FLAG_COMP1;
    }

    if (bio->pair) is_pair = 1;

    recno = save_range_sequence(bio->io, &s, s.mapping_qual, bio->pair,
    				is_pair, tname, bio->c, bio->a, flags, lib,
				&bin_rec);


#ifdef SAM_AUX_AS_TAG
    /* Make an annotation out of the sam auxillary data */
    aux = bam_aux_stringify(b, 1);
    if (aux && *aux) {
	anno_ele_t *e;
	bin_index_t *bin;
	range_t r;

	r.mqual = str2type("SAMX");
	r.start = s.pos;
	r.end = s.pos;
	r.pair_rec = recno;
	r.flags = GRANGE_FLAG_ISANNO | GRANGE_FLAG_TAG_SEQ;
	r.rec = anno_ele_new(bio->io, 0, GT_Seq, recno, 0, r.mqual,
			     ANNO_DIR_NUL, aux);
	e = (anno_ele_t *)cache_search(bio->io, GT_AnnoEle, r.rec);
	e = cache_rw(bio->io, e);
	
	bin = bin_add_to_range(bio->io, &bio->c, bin_rec, &r, NULL, NULL, 0);
	e->bin = bin->rec;
    }
#endif

 anno_only:

    /* Add new style tags */
    if ((tags = bam_aux_find(b, "PT"))) {
	int start, end, type_len, text_len;
	char dir, *type, *text, tmp;
	char tag_type[5], *tag_text;
	tg_rec orec;
	int otype;
	range_t r;
	anno_ele_t *e;
	bin_index_t *bin;

	while (tags = parse_bam_PT_tag(tags, &start, &end, &dir,
				       &type, &type_len, &text, &text_len)) {
	    //printf("Tag %d..%d dir %c type=%.*s text=%.*s\n",
	    //	   start, end, dir, type_len, type, text_len, text);

	    if (start < 1)
		start = 1;
	    if (end > bs->seq_len)
		end = bs->seq_len;

	    start = bs->pad[start-1]+1;
	    end   = bs->pad[end-1]+1;

	    strncpy(tag_type, type, 4);
	    tag_type[4] = 0;

	    /* Create the tag */
	    if (fake) {
		orec  = r.pair_rec = bio->c->rec;
		otype = GT_Contig;
		r.flags = GRANGE_FLAG_ISANNO;

		r.start = start;
		r.end   = end;
	    } else {
		orec  = r.pair_rec = recno;
		otype = GT_Seq;
		r.flags = GRANGE_FLAG_ISANNO | GRANGE_FLAG_TAG_SEQ;

		/* FIXME: Re-pad start/end */

		if (start < 1 || end > ABS(s.len)) {
		    verror(ERR_WARN, "sam_import", "Anno. range (%d..%d) is "
			   "outside of sequence range (%d..%d)",
			   start, end, 1, ABS(s.len));
		    if (start < 1)
			start = 1;
		    if (end > ABS(s.len))
			end = ABS(s.len);
		}
		r.start = s.pos + start-1;
		r.end   = s.pos + end-1;
	    }

	    r.mqual = str2type(tag_type);

	    if (text) {
		tmp = text[text_len];
		text[text_len] = 0;
		r.rec   = anno_ele_new(bio->io, 0, otype, orec, 0, r.mqual,
				       dir, text);
		text[text_len] = tmp;
	    } else {
		r.rec   = anno_ele_new(bio->io, 0, otype, orec, 0, r.mqual,
				       dir, NULL);
	    }

	    /* Link it to a bin */
	    e = (anno_ele_t *)cache_search(bio->io, GT_AnnoEle, r.rec);
	    e = cache_rw(bio->io, e);

	    if (fake) {
		bin = bin_add_range(bio->io, &bio->c, &r, NULL, NULL, 0);
	    } else {
		bin = bin_add_to_range(bio->io, &bio->c, bin_rec, &r,
				       NULL, NULL, 0);
	    }
	    e->bin = bin->rec;
	}
    }

    if ((tags = bam_aux_find(p->b, "CT"))) {
	int start, end, type_len, text_len;
	char dir, *type, *text, tmp;
	char tag_type[5], *tag_text;
	tg_rec orec;
	int otype;
	range_t r;
	anno_ele_t *e;
	bin_index_t *bin;

	start = bs->pos;
	end   = bs->pos + bs->seq_len-1;

	while (tags = parse_bam_CT_tag(tags, &dir,
				       &type, &type_len, &text, &text_len)) {
	    strncpy(tag_type, type, 4);
	    tag_type[4] = 0;

	    /* Create the tag */
	    if (fake) {
		orec  = r.pair_rec = bio->c->rec;
		otype = GT_Contig;
		r.flags = GRANGE_FLAG_ISANNO;
	    } else {
		orec  = r.pair_rec = recno;
		otype = GT_Seq;
		r.flags = GRANGE_FLAG_ISANNO | GRANGE_FLAG_TAG_SEQ;
	    }

	    r.mqual = str2type(tag_type);
	    r.start = start;
	    r.end   = end;

	    if (text) {
		tmp = text[text_len];
		text[text_len] = 0;
		r.rec   = anno_ele_new(bio->io, 0, otype, orec, 0, r.mqual,
				       dir, text);
		text[text_len] = tmp;
	    } else {
		r.rec   = anno_ele_new(bio->io, 0, otype, orec, 0, r.mqual,
				       dir, NULL);
	    }

	    /* Link it to a bin */
	    e = (anno_ele_t *)cache_search(bio->io, GT_AnnoEle, r.rec);
	    e = cache_rw(bio->io, e);

	    if (fake) {
		bin = bin_add_range(bio->io, &bio->c, &r, NULL, NULL, 0);
	    } else {
		bin = bin_add_to_range(bio->io, &bio->c, bin_rec, &r,
				       NULL, NULL, 0);
	    }
	    e->bin = bin->rec;
	}
    }

    /* Add old style tags */
    handle = NULL;
    while (0 == bam_aux_iter(b, &handle, aux_key, &type, &val)) {
	range_t r;
	anno_ele_t *e;
	bin_index_t *bin;
	char *tokens[4], *cp, *tag_text, tag_type[5];
	int ntok;
	int tag_pos, tag_len;
	tg_rec orec;
	int otype;

	if (!(aux_key[0] == 'Z' && (aux_key[1] == 's' || aux_key[1] == 'c')))
	    continue;

	tokens[0] = val.s;
	for (ntok = 1, cp = val.s; *cp && ntok < 4; cp++) {
	    if (*cp == '|') {
		*cp = 0;
		tokens[ntok++] = cp+1;
	    }
	}

	/* Parse it */
	tag_type[0] = tag_type[1] = tag_type[2] = tag_type[3] = '-';
	tag_type[4] = 0;
	strncpy(tag_type, tokens[0], 4);
	tag_pos  = ntok >= 2 ? atoi(tokens[1]) : 0;
	tag_len  = ntok >= 3 ? atoi(tokens[2]) : 0;
	tag_text = ntok >= 4 ? (unescape_line(tokens[3]),tokens[3]) : NULL;

	/* Create the tag */
	r.mqual    = str2type(tag_type);
	r.start    = tag_pos-1 + s.pos;
	r.end      = tag_pos-1 + s.pos + tag_len-1;
	if (aux_key[1] == 'c') {
	    orec = r.pair_rec = bio->c->rec;
	    otype = GT_Contig;
	    r.flags = GRANGE_FLAG_ISANNO;
	} else {
	    orec = r.pair_rec = recno;
	    otype = GT_Seq;
	    r.flags = GRANGE_FLAG_ISANNO | GRANGE_FLAG_TAG_SEQ;
	    if (r.start < s.pos || r.end > s.pos + ABS(s.len)-1) {
		verror(ERR_WARN, "sam_import", "Anno. range (%d..%d) is "
		       "outside of sequence range (%d..%d)",
		       r.start, r.end, s.pos, s.pos + ABS(s.len)-1);
		if (r.start < s.pos)
		    r.start = s.pos;
		if (r.end > s.pos + ABS(s.len)-1)
		    r.end = s.pos + ABS(s.len)-1;
		if (r.start > r.end) {
		    int tmp = r.start;
		    r.start = r.end;
		    r.end = tmp;
		}
	    }
	}
	r.rec = anno_ele_new(bio->io, 0, otype, orec, 0, r.mqual, 
			     ANNO_DIR_NUL, tag_text);

	/* Link it to a bin */
	e = (anno_ele_t *)cache_search(bio->io, GT_AnnoEle, r.rec);
	e = cache_rw(bio->io, e);

	if (aux_key[1] == 's') {
	    bin = bin_add_to_range(bio->io, &bio->c, bin_rec, &r,
				   NULL, NULL, 0);
	} else {
	    bin = bin_add_range(bio->io, &bio->c, &r, NULL, NULL, 0);
	}
	e->bin = bin->rec;
    }

    //printf("Del seq %p / %p => %p,%p\n", p, bs, bs->seq, bs->conf);

    /* unlink */
    bs->next = bio->free_seq;
    bio->free_seq = bs;

    return 0;
}


/*
 * Called once per sequence as they're discovered
 *
 * Returns 0 to reject seq
 *         1 to accept seq
 *        -1 on error
 */
static int sam_check_unmapped(void *cd, bam_file_t *fp, pileup_t *p) {
    bam_io_t *bio = (bam_io_t *)cd;
    

    if ((++bio->total_count & 0xffff) == 0) {
	putchar('.');
	fflush(stdout);
	cache_flush(bio->io);
    }

    if (bam_flag(p->b) & BAM_FUNMAP) {
	if (!bio->a->store_unmapped)
	    return 0;

	bio_add_unmapped(bio, p->b);
    }

    if ((bam_flag(p->b) & BAM_FDUP) && bio->a->remove_dups)
	return 0;

    return 1;
}

static int sam_add_seq(void *cd, bam_file_t *fp, pileup_t *p,
		       int depth, int pos, int nth) {
    int tid, np = 0;
    bam_io_t *bio = (bam_io_t *)cd;

    if (!p)
	return 0;

    /* New contig? */
    tid = p->b->ref;
    if (tid != bio->last_tid)
	bio_new_contig(bio, tid);

    /* tg_index -g mode */
    if (bio->a->repad) {
	pos += bio->c_start-1;

	/* Pad these sequences based on existing pads in the padded contig */
	if ((np=padtree_pad_at(bio->tree, pos+bio->n_inserts-bio->npads))) {
	    int j;

	    //printf("Pos %d pads %d\n", pos, np);

	    /* Add pads to match existing consensus gaps */
	    for (j = nth; j < np; j++) {
		pileup_t *pcopy;

		if (bio->skip) {
		    bio->skip--;
		    continue;
		}
		    
		//printf("Import pads from existing consensus at %d\n",
		//       pos+bio->n_inserts-bio->npads);

		pcopy = p;
		for (; p; p = p->next) {
		    bio_seq_t *s;

		    if (p->start)
			continue;

		    s = (bio_seq_t *)p->cd;
		    if (s->seq_len + np + 1 >= s->alloc_len) {
			s->alloc_len = (s->alloc_len + 100 + np+1)*1.5;
			s->seq  = (char *)realloc(s->seq,  s->alloc_len);
			s->conf = (char *)realloc(s->conf, s->alloc_len);
			if (!s->seq || !s->conf)
			    return -1;
		    }
		    
		    if (bio->a->data_type & DATA_SEQ) {
			s->seq [s->seq_len] = '*';
		    } else {
			s->seq [s->seq_len] = 'N';
		    }
		    if (bio->a->data_type & DATA_QUAL) {
			s->conf[s->seq_len] = 2;
		    } else {
			s->conf[s->seq_len] = 0;
		    }
		    s->seq_len++;
		}
		p = pcopy;
	    }
	}

	/* And vice versa - pad existing contig based on new non-ref bases */
	if (nth) {
	    /* Check how many insertions were originally at this point */
	    np = padtree_pad_at(bio->tree, pos+1);
	    //printf("Pos %d nth %d, old pad count = %d\n", pos, nth, np);
	    
	    if (nth > np) {
		contig_insert_base(bio->io, &bio->c,
				   get_padded_coord(bio->tree, pos+1) +
				   bio->n_inserts,
				   '*', 0);

		bio->npads++;
	    } else {
		bio->n_inserts--;
		bio->skip++;
	    }
	}
    }


    /* Finally add the new column of base calls */
    for (; p; p = p->next) {
	bio_seq_t *s;

	/* New sequence */
	if (p->start) {
	    int ppos = bio->tree
		? bio->npads +
		  get_padded_coord(bio->tree,
				   pos + bio->n_inserts - bio->npads)
		: pos + bio->n_inserts;
	    p->cd = s = bio_new_seq(bio, p, ppos);
	}

	s = (bio_seq_t *)p->cd;

	/* Extend sequence */
	if (s->seq_len >= s->alloc_len) {
	    s->alloc_len = (s->alloc_len + 100)*1.5;
	    if (NULL == (s->seq  = (char *)realloc(s->seq,  s->alloc_len)))
		return -1;
	    if (NULL == (s->conf = (char *)realloc(s->conf, s->alloc_len)))
		return -1;
	    if (NULL == (s->pad = (int *)realloc(s->pad,
						 s->alloc_len * sizeof(int))))
		return -1;
	}
	if (p->seq_offset < 0 || p->first_del || (p->base == '*' && s->left == s->seq_len)) {
	    s->pos++;
	} else {
	    if (bio->a->data_type & DATA_SEQ) {
		s->seq [s->seq_len] = p->base;
	    } else {
		s->seq [s->seq_len] = 'N';
	    }
	    if (bio->a->data_type & DATA_QUAL) {
		s->conf[s->seq_len] = p->qual;
	    } else {
		s->conf[s->seq_len] = 0;
	    }

	    s->pad[s->padded_pos] = s->seq_len;
	    s->padded_pos += 1 - p->padding;
	    s->seq_len++;
	}

	/* Remove sequence */
	if (p->eof & 1) {
	    //printf("End seq %s\n", bam_name(p->b));

	    /* Either fake seq (consensus annotation) or real */
	    bio_del_seq(bio, p);
	}
    }
    if (nth)
	bio->n_inserts++;

    /*
     * Store mapping of reference to padded coordinates.
     */
    if (bio->a->store_refpos && nth) {
	range_t r;
	int ppos = bio->tree
	    ? bio->npads +
	    get_padded_coord(bio->tree,
			     pos + bio->n_inserts - bio->npads)
	    : pos + bio->n_inserts;

	//printf("%csam_seqadd(..., pos=%d, nth=%d) => {%d, %d}\n",
	//       " *"[nth > 0], pos, nth, pos, ppos);

	r.start    = ppos;
	r.end      = ppos;
	r.rec      = tid;   /* ref seq ID */
	r.pair_rec = 0;     /* size of deletion */
	r.mqual    = pos;   /* map from ppos -> pos */
	r.flags    = GRANGE_FLAG_ISREFPOS
	           | GRANGE_FLAG_REFPOS_INS
	           | GRANGE_FLAG_REFPOS_FWD;
	
	bin_add_range(bio->io, &bio->c, &r, NULL, NULL, 1);
    }

    return 0;
}

/*
 * Initialise bio->libs hache with current libraries. Useful for tg_index -a
 */
static void bio_init_libs(bam_io_t *bio) {
    GapIO *io = bio->io;
    int i;

    for (i = 0; i < io->db->Nlibraries; i++) {
	tg_rec rec = ARR(tg_rec, io->library, i);
	library_t *lib = cache_search(io, GT_Library, rec);
	HacheData hd;

	if (!lib)
	    continue;

	cache_incr(io, lib);
	hd.p = lib;
	HacheTableAdd(bio->libs, lib->name, strlen(lib->name), hd, NULL);
    }
}

int parse_sam_or_bam(GapIO *io, char *fn, tg_args *a, char *mode) {
    bam_io_t *bio = (bam_io_t*)calloc(1, sizeof(*bio));
    bam_file_t *fp;

    /* Setup bam_io_t object and create our pileup interface */
    bio->io = io;
    bio->seqs = NULL;
    bio->free_seq = NULL;
    bio->max_seq = 0;
    bio->a = a;
    bio->c = NULL;
    bio->count = 0;
    bio->total_count = 0;
    bio->fn = fn;
    bio->libs = HacheTableCreate(256, HASH_DYNAMIC_SIZE);
    bio->libs->name = "libs";
    bio->last_tid = -1;
    bio->tree = NULL;

    if (a->pair_reads) {
	bio->pair = create_pair(a->pair_queue);
    } else {
	bio->pair = NULL;
    }

    bio_init_libs(bio);

    fp = bam_open(fn, mode);
    if (!fp)
	return -1;
    bio->fp = fp;

//    if (!bio->header->dict) {
//	bio->header->dict = sam_header_parse2(bio->header->text);
//    }
//    bio->rg2pl_hash = sam_header2tbl(bio->header->dict, "RG", "ID", "PL");

    /* The main processing loop, calls sam_add_seq() */
    if (0 != pileup_loop(fp, sam_check_unmapped, sam_add_seq, bio)) {
	verror(ERR_WARN, "sam_import", "pileup failed processing line %d",
	       fp->line);
	cache_flush(io);
	bam_close(fp);
	return -1;
    }

    //pileup_loop(fp, NULL, sam_add_seq, bio);

    //    if (bio->rg2pl_hash)
    //	sam_tbl_destroy(bio->rg2pl_hash);

    cache_flush(io);
    vmessage("Loaded %d of %d sequences\n", bio->count, bio->total_count);

    if (bio->pair && !a->fast_mode) { 	
	finish_pairs(io, bio->pair);
    }
 
    /* Tidy up */
    if (fp)
	bam_close(fp);

    if (bio) {
	bio_seq_t *s, *n;

	if (bio->pair) delete_pair(bio->pair);

	if (bio->libs) {
	    /* call cache_decr on each lib too */
	    HacheItem *hi;
	    HacheIter *iter;

	    if (!(iter = HacheTableIterCreate()))
		return -1;
	
	    while (hi = HacheTableIterNext(bio->libs, iter)) {
		library_t *lib = hi->data.p;
		cache_decr(io, lib);
	    }
	    
	    HacheTableIterDestroy(iter);
	    HacheTableDestroy(bio->libs, 0);
	}

	if (bio->seqs)
	    free(bio->seqs);

	if (bio->tree)
	    depad_seq_tree_free(bio->tree);

	for (s = bio->free_seq; s; s = n) {
	    n = s->next;
	    if (s->seq)
		free(s->seq);
	    if (s->conf)
		free(s->conf);
	    if (s->pad)
		free(s->pad);
	    free(s);
	}

	if (bio->c)
	    cache_decr(io, bio->c);

	free(bio);
    }

    return 0;
}

int parse_bam(GapIO *io, char *fn, tg_args *a) {
    return parse_sam_or_bam(io, fn, a, "rb");
}

int parse_sam(GapIO *io, char *fn, tg_args *a) {
    return parse_sam_or_bam(io, fn, a, "r");
}
