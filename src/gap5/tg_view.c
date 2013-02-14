/*  Last edited: Jun 15 11:43 2009 (badger) */
/*
 * Author:         James Bonfield, Feb 2007
 *                 Wellcome Trust Sanger Institute
 *
 * g_view:
 * Loads and displays the contents of a g-library format database as created
 * by the g_index application.
 */

#include <staden_config.h>

#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#ifdef HAVE_LIBCURSES
#  include LIBCURSES_HEADER
#else
#  error "No curses library installed"
#endif
#include <signal.h>
#include <ctype.h>
#include <string.h>

#include "array.h"
#include "misc.h"
#include "tg_gio.h"

//#define TEST_MODE

#define get_seq(io, rec) ((seq_t *)cache_search((io), GT_Seq, (rec)))

/* ------------------------------------------------------------------------ */
/* Consensus generation functions */

/* ACGT to 0123 conversion */
static unsigned char lookup[256], lookup_done = 0;

/*
 * Compute a basic non-weighted consensus. We simply pick the basecall
 * most frequently used.
 *
 * FIXME: use a weighted sum based on confidence values instead?
 */
int calc_cons(GapIO *io, rangec_t *r, int nr, int xpos, int wid,
	      char *cons) {
    int i, j;
    int (*cvec)[6] = (int (*)[6])calloc(wid, 6 * sizeof(int));

    if (!lookup_done) {
	memset(lookup, 5, 256);
	lookup_done = 1;
	lookup['A'] = lookup['a'] = 0;
	lookup['C'] = lookup['c'] = 1;
	lookup['G'] = lookup['g'] = 2;
	lookup['T'] = lookup['t'] = 3;
	lookup['*'] = lookup[','] = 4;
    }

    /* Accumulate */
    for (i = 0; i < nr; i++) {
	int sp = r[i].start;
	seq_t *s = get_seq(io, r[i].rec);
	seq_t *sorig = s;
	int l = s->len > 0 ? s->len : -s->len;
	unsigned char *seq;
	int left, right;

	/* Complement data on-the-fly */
	if ((s->len < 0) ^ r[i].comp) {
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	seq = (unsigned char *)s->seq;
	left = s->left;
	right = s->right;

	if (sp < xpos) {
	    seq   += xpos - sp;
	    l     -= xpos - sp;
	    left  -= xpos - sp;
	    right -= xpos - sp;
	    sp = xpos;
	}
	if (l > wid - (sp-xpos))
	    l = wid - (sp-xpos);
	if (left < 1)
	    left = 1;

	for (j = left-1; j < right; j++) {
	    if (sp-xpos+j < wid)
		cvec[sp-xpos+j][lookup[seq[j]]]++;
	}

	if (s != sorig)
	    free(s);
    }

    memset(cons, ' ', wid);

    /* and speculate :-) */
    for (i = 0; i < wid; i++) {
	int max, max_base = 5;
	for (max = j = 0; j < 6; j++) {
	    if (max < cvec[i][j]) {
		max = cvec[i][j];
		max_base = j;
	    }
	}
	cons[i] = "ACGT*N"[max_base];
    }

    free(cvec);

    return 0;

}

/* ------------------------------------------------------------------------ */
/* Curses display and keyboard handling routines */
static WINDOW *gotowin = NULL;
static WINDOW *helpwin = NULL;

#define DISPLAY_COLOURS 1
#define DISPLAY_CUTOFFS 2
#define DISPLAY_QUAL    4
#define DISPLAY_DIFFS   8

void exit_curses(int sig) {
    endwin();
    exit(0);
}

void init_curses(void) {
    signal(SIGINT, exit_curses);

    initscr();
    keypad(stdscr, TRUE);
    start_color();
    init_pair(0, COLOR_RED,    COLOR_BLACK);
    init_pair(1, COLOR_GREEN,  COLOR_BLACK);
    init_pair(2, COLOR_BLUE,   COLOR_BLACK);
    init_pair(3, COLOR_YELLOW, COLOR_BLACK);
    init_pair(4, COLOR_WHITE,  COLOR_BLACK);
    init_pair(5, COLOR_WHITE,  COLOR_BLACK);
    clear();
    noecho();
    cbreak();

    gotowin = newwin(3, 36, 10, 5);
    helpwin = newwin(22,40, 1, 3);
}

int wgotonum(int xpos, WINDOW *win) {
    int orig = xpos;
    
    /* wborder(win, 0, 0, 0, 0, 0, 0, 0, 0); */
    wborder(win, '|', '|', '-', '-', '+', '+', '+', '+');

    mvwprintw(win, 1, 2, "Goto base:                       ", xpos);
    mvwprintw(win, 1, 13, "%-d", xpos);

    for (;;) {
	int c;
	wrefresh(win);

	switch (c = wgetch(win)) {
	case '0':
	case '1':
	case '2':
	case '3':
	case '4':
	case '5':
	case '6':
	case '7':
	case '8':
	case '9':
	    xpos = xpos * 10 + c - '0';
	    break;

	case KEY_BACKSPACE:
	case '\010':
	case '\177':
	    xpos /= 10;
	    break;

	case '\027': /* control w */
	    xpos = 0;
	    break;

	case '\033': /* escape */
	    return orig;

	case KEY_ENTER:
	case '\012':
	case '\015':
	    return xpos;
	}

	mvwprintw(win, 1, 13, "            ");
	mvwprintw(win, 1, 13, "%-d", xpos);
    }
}

char *wgotoseq(WINDOW *win) {
    static char name[1024];
    int cursor = 0;
    
    /* wborder(win, 0, 0, 0, 0, 0, 0, 0, 0); */
    wborder(win, '|', '|', '-', '-', '+', '+', '+', '+');

    mvwprintw(win, 1, 2, "Goto seq.:                       ");
    *name = 0;

    for (;;) {
	int c;

	mvwprintw(win, 1, 13, "%s", name);
	wrefresh(win);

	switch (c = wgetch(win)) {
	case KEY_BACKSPACE:
	case '\010':
	case '\177':
	    if (cursor > 0)
		name[--cursor] = 0;
	    break;

	case '\027': /* control w */
	    cursor = 0;
	    name[cursor] = 0;
	    break;

	case '\033': /* escape */
	    return "";

	case KEY_ENTER:
	case '\012':
	case '\015':
	    return name;

	default:
	    name[cursor++] = c;
	    name[cursor] = 0;
	}

	mvwprintw(win, 1, 13, "            ");
	mvwprintw(win, 1, 13, "%s", name);
    }
}

void whelp(WINDOW *win) {
    int r = 1;

    /* wborder(win, 0, 0, 0, 0, 0, 0, 0, 0); */
    wborder(win, '|', '|', '-', '-', '+', '+', '+', '+');

    mvwprintw(win, r++, 2, "        -=-    Help    -=- ");
    r++;
    mvwprintw(win, r++, 2, "?         This window");
    mvwprintw(win, r++, 2, "Arrows    Small scroll movement");
    mvwprintw(win, r++, 2, "h,j,k,l   Small scroll movement");
    mvwprintw(win, r++, 2, "H,J,K,L   Large scroll movement");
    mvwprintw(win, r++, 2, "ctrl-H    Scroll 1k left");
    mvwprintw(win, r++, 2, "ctrl-L    Scroll 1k right");
    mvwprintw(win, r++, 2, "c         Toggle cutoffs");
    mvwprintw(win, r++, 2, "C         Toggle colours");
    mvwprintw(win, r++, 2, "q         Toggle quality");
    mvwprintw(win, r++, 2, "d         Toggle differences");
    mvwprintw(win, r++, 2, "-         Decrement qual. threshold");
    mvwprintw(win, r++, 2, "+         Increment qual. threshold");
    mvwprintw(win, r++, 2, "g         Go to specific location");
    mvwprintw(win, r++, 2, "<         Go to start of contig");
    mvwprintw(win, r++, 2, ">         Go to end of contig");
    mvwprintw(win, r++, 2, "n         Next contig");
    mvwprintw(win, r++, 2, "p         Previous contig");
    mvwprintw(win, r++, 2, "x         Exit");
    wrefresh(win);
    wgetch(win);
}

static void complement_bin(GapIO *io, tg_rec bnum) {
    bin_index_t *bin = get_bin(io, bnum);
    bin->flags ^= BIN_COMPLEMENTED;
}

static void display_gap(GapIO *io, contig_t **c, int xpos, int ypos,
			int nlines, int wid, int mode, int qual_cutoff,
			int in_curses) {
    rangec_t *r;
    int i, nr, lno, y;
    char line[1024], *lp;
    char cons[1024];
    int attr;
    static int lookup_1conf[256];
    static int lookup_4conf[256];
    static int lookup_init = 0;

    if (!lookup_init) {
	for (i = 0; i < 256; i++)
	    lookup_1conf[i] = lookup_4conf[0] = 0;

	lookup_4conf['a'] = lookup_4conf['A'] = 0;
	lookup_4conf['c'] = lookup_4conf['C'] = 1;
	lookup_4conf['g'] = lookup_4conf['G'] = 2;
	lookup_4conf['t'] = lookup_4conf['T'] = 3;
    }

    wid -= MAX_NAME_LEN+2;

    //if (xpos < wid/2 + (*c)->start)
    //	xpos = wid/2 + (*c)->start;

    xpos -= wid/2;

    /* Query visible objects */
    r = contig_seqs_in_range(io, c, xpos, xpos+wid-1, CSIR_SORT_BY_X, &nr);

    /* Consensus */
    calc_cons(io, r, nr, xpos, wid, cons);
    if (in_curses) {
	clear();
	mvaddnstr(0, 1, contig_get_name(c), strlen(contig_get_name(c)));
	mvaddnstr(0, MAX_NAME_LEN+2, cons, wid);
    } else {
	printf(" %-*s %.*s\n", MAX_NAME_LEN, contig_get_name(c), wid, cons);
    }

    /* Position */
    for (lp = line, i = xpos; i < xpos+wid+19; i++) {
	if (i % 10 == 0) {
	    sprintf(lp, "%10d", i-10);
	    lp += 10;
	}
    }
    if (in_curses) {
	int m = (xpos-1)%10;
	if (m < 0) m += 10;
	mvaddnstr(1, MAX_NAME_LEN+2, line+10+m, wid);
    } else {
	printf("%*s%.*s\n", MAX_NAME_LEN+2, "", wid,
	       line+9+((xpos-1)%10));
    }


    /* Sequences */
    for (i = y = 0; i < nr && y < ypos; i++, y++);
    for (lno = 2; i < nr && lno < nlines; i++, lno++) {
	seq_t *s = get_seq(io, r[i].rec);
	seq_t *sorig = s;
	int sp = r[i].start;
	int l = s->len > 0 ? s->len : -s->len;
	unsigned char seq_a[MAX_SEQ_LEN], *seq = seq_a;
	int j, dir = '+';
	int left, right;
	int8_t *conf;
	int nc = s->format == SEQ_FORMAT_CNF4 ? 4 : 1;
	int *L = s->format == SEQ_FORMAT_CNF4 ? lookup_4conf : lookup_1conf;

	/* Complement data on-the-fly */
	if ((s->len < 0) ^ r[i].comp) {
	    dir = '-';
	    s = dup_seq(s);
	    complement_seq_t(s);
	}

	left = s->left;
	right = s->right;

	memcpy(seq, s->seq, l);
	conf = s->conf;

	if (sp < xpos) {
	    seq   += xpos - sp;
	    conf  += nc * (xpos - sp);
	    l     -= xpos - sp;
	    left  -= xpos - sp;
	    right -= xpos - sp;
	    sp = xpos;
	}
	if (l > wid - (sp-xpos))
	    l = wid - (sp-xpos);

	if (in_curses) {
	    /* Test of sequence_get_position */
	    /*
	      int c, p;
	      sequence_get_position(io, r[i].rec, &c, &p);
	      s->name_len = sprintf(s->name, ":%d-%d:", p, p+ABS(s->len)-1);
	    */
	    mvaddch(lno, 0, dir);
	    addnstr(s->name, MIN(MAX_NAME_LEN, s->name_len));
	    move(lno, MAX_NAME_LEN+2+sp-xpos);
	} else {
	    printf("%c%.*s%*s",
		   dir,
		   MIN(MAX_NAME_LEN, s->name_len), s->name,
		   MAX_NAME_LEN+1-MIN(MAX_NAME_LEN, s->name_len) +sp-xpos, "");
	}

	for (j = 0; j < l; j++) {
	    attr = (mode & DISPLAY_COLOURS) ? COLOR_PAIR(lookup[seq[j]]) : 0;

	    if (mode & DISPLAY_DIFFS
		&& sp-xpos+j < wid && seq[j] == cons[sp-xpos+j])
		seq[j] = '.';
	    if (j < left-1 || j > right-1)
		seq[j] = (mode & DISPLAY_CUTOFFS) ? tolower(seq[j]) : ' ';

	    if (conf[j*nc+L[seq[j]]] >= qual_cutoff && mode & DISPLAY_QUAL) {
		attr |= A_BOLD;
	    }

	    if (in_curses) {
		addch(seq[j] | attr);
	    } else {
		putchar(seq[j]);
	    }
	}

	if (!in_curses)
	    putchar('\n');

	if (s != sorig)
	    free(s);
    }

    /* Useful debugging code to show bin locations. */
#if 0
    free(r);
    r = contig_bins_in_range(io, c, xpos, xpos+wid-1, 0, 0, &nr);
    /* Bins */
    for (i=0; i < nr && lno < nlines; i++, lno++) {
	bin_index_t *bin = (bin_index_t *)cache_search(io, GT_Bin, r[i].rec);
	unsigned char *seq, *seqm;
	int j, dir = "+-"[r[i].comp];
	int sp = r[i].start;
	int l = ABS(r[i].end - r[i].start + 1);
	char name[100];

	sprintf(name, "bin-%"PRIrec, bin->rec);
	seqm = seq = malloc(l+1);
	memset(seq, '-', l);

	if (!(bin->start_used == 0 && bin->end_used == 0)) {
	    if (r[i].comp) {
		memset(&seq[bin->size - bin->end_used - 1], '=',
		       bin->end_used - bin->start_used + 1);
	    } else {
		memset(&seq[bin->start_used], '=',
		       bin->end_used - bin->start_used + 1);
	    }
	}

	/*
	fprintf(stderr, "Bin-%d: %d+%d %d..%d\n",
		bin->rec,
		bin->pos, bin->size,
		bin->start_used, bin->end_used);
	*/

	if (sp < xpos) {
	    seq   += xpos - sp;
	    l     -= xpos - sp;
	    sp = xpos;
	}
	if (l > wid - (sp-xpos))
	    l = wid - (sp-xpos);

	if (in_curses) {
	    mvaddch(lno, 0, dir);
	    addnstr(name, strlen(name));
	    move(lno, MAX_NAME_LEN+2+sp-xpos);
	} else {
	    printf("%c%.*s%*s",
		   dir,
		   (int)MIN(MAX_NAME_LEN, strlen(name)),
		   name,
		   (int)(MAX_NAME_LEN+1-MIN(MAX_NAME_LEN,
					    strlen(name)) +sp-xpos),
		   "");
	}

	for (j = 0; j < l; j++) {
	    if (in_curses) {
		addch(seq[j]);
	    } else {
		putchar(seq[j]);
	    }
	}

	if (!in_curses)
	    putchar('\n');

	free(seqm);
    }
#endif

    if (in_curses)
	refresh();

    free(r);
}


void next_contig(GapIO *io, contig_t **cp) {
    cache_decr(io, *cp);

    if (++io->contig_num == io->db->Ncontigs)
	io->contig_num = 0;
    
    gio_read_contig(io, io->contig_num, cp);
    cache_incr(io, *cp);
}

void prev_contig(GapIO *io, contig_t **cp) {
    cache_decr(io, *cp);

    if (io->contig_num == 0)
	io->contig_num = io->db->Ncontigs;
    io->contig_num--;

    gio_read_contig(io, io->contig_num, cp);
    cache_incr(io, *cp);
}

/* Write tests */
int edit_contig_name(GapIO *io, contig_t **cp) {
    char name[1024];
    int i;

    strcpy(name, contig_get_name(cp));
    for (i = 0; name[i]; i++)
	if (isalpha(name[i]))
	    name[i] ^= 0x20; /* change case */

    return contig_set_name(io, cp, name);
}

int save(GapIO *io) {
    return cache_flush(io);
}

void curses_loop(GapIO *io, contig_t **cp, int xpos, int mode) {
    int ypos = 0;
    int maxy, maxx;
    int qual_cutoff = 29;

#ifdef NCURSES_VERSION
    getmaxyx(stdscr, maxy, maxx);
#else
    maxx = 80;
    maxy = 40;
#endif

    if (xpos < (maxx-MAX_NAME_LEN+2)/2)
	xpos = (maxx-MAX_NAME_LEN+2)/2;

    display_gap(io, cp, xpos, ypos, maxy, maxx, mode, qual_cutoff, 1);

    for(;;) {
	switch (getch()) {
	case '!':
	    edit_contig_name(io, cp);
	    break;

	case '"': /* control-x */
	    save(io);
	    break;

	case 'd':
	    mode ^= DISPLAY_DIFFS;
	    break;

	case 'c':
	    mode ^= DISPLAY_CUTOFFS;
	    break;

	case 'C':
	    mode ^= DISPLAY_COLOURS;
	    break;

	case 'q':
	    mode ^= DISPLAY_QUAL;
	    break;

	case 'g':
	    xpos = wgotonum(xpos, gotowin);
	    break;

	case 'G': {
	    char *seq = wgotoseq(gotowin);
	    tg_rec n;
	    if (seq) {
		if ((n = sequence_index_query(io, seq)) < 0) {
		    putchar('\a');
		    fflush(stdout);
		} else {
		    tg_rec c;
		    int x;
		    sequence_get_position(io, n, &c, &x, NULL, NULL);
		    /* FIXME: check c and *cp are same contig */
		    xpos = x;
		}
	    }
	    break;
	}
	case '?':
	    whelp(helpwin);
	    break;

	case '-':
	    qual_cutoff--;
	    break;

	case '+':
	    qual_cutoff++;
	    break;

	case 'n':
	    next_contig(io, cp);
	    xpos = ypos = 0;
	    break;

	case 'p':
	    prev_contig(io, cp);
	    xpos = ypos = 0;
	    break;

	/* Orthoganol movements */
	case KEY_RIGHT:
	case 'l':
	    xpos++;
	    break;

	case KEY_SRIGHT:
	case 'L':
	    xpos+=20;
	    break;

	case '\014':
	    xpos+=1000;
	    break;

	case KEY_LEFT:
	case 'h':
	    xpos--;
	    break;

	case KEY_SLEFT:
	case 'H':
	    xpos-=20;
	    break;

	case KEY_BACKSPACE:
	case '\010':
	    xpos-=1000;
	    break;

	case KEY_DOWN:
	case 'j':
	    ypos++;
	    break;

	case 'J':
	    ypos+=20;
	    break;

	case KEY_UP:
	case 'k':
	    ypos--;
	    break;

	case 'K':
	    ypos-=20;
	    break;

	case 'x':
	    return;

	case '<':
	    xpos = (*cp)->start;
	    break;

	case '>':
	    xpos = (*cp)->end;
	    break;

	case 'f':
	    complement_bin(io, contig_get_bin(cp));
	    break;

	case 'i':
	    contig_insert_base(io, cp, xpos, '-', 20);
	    break;

	case 'o':
	    contig_delete_base(io, cp, xpos);
	    break;

#ifdef KEY_RESIZE
	case KEY_RESIZE:
	    getmaxyx(stdscr, maxy, maxx);
	    break;
#endif

	default:
	    continue;
	}

	if (xpos < (*cp)->start) xpos = (*cp)->start;
	if (xpos > (*cp)->end)   xpos = (*cp)->end;
	if (ypos < 0) ypos = 0;

	display_gap(io, cp, xpos, ypos, maxy, maxx, mode, qual_cutoff, 1);
    }
}

/*
 * Prints output to stdout without any curses control. This is a distilled
 * down display_gap() function.
 */
void print_output(GapIO *io, contig_t **c, int xpos, int width, int mode) {
    //complement_store(io, contig_get_bin(c));
    display_gap(io, c, xpos, 0, 1000000, width, mode, 0, 0);
    return;
}

#ifdef TEST_MODE
/* ------------------------------------------------------------------------ */
/* Debug functions that don't use curses - handy for valgrind testing */
static void test_mode(GapIO *io, contig_t **c, int xpos) {
    rangec_t *r;
    int nr, i;

    r = contig_seqs_in_range(io, c, xpos, xpos+79, CSIR_SORT_BY_X, &nr);
    for (i = 0; i < nr; i++) {
	seq_t *s = get_seq(io, r[i].rec);
	printf("%.*s: range %d..%d seq %d+%d st=%d en=%d %.*s\n", 
	       s->name_len, s->name,
	       r[i].start, r[i].end,
	       s->pos, s->len,
	       s->left, s->right,
	       ABS(s->len), s->seq);

	s = dup_seq(s);
	complement_seq_t(s);

	printf("%.*s: range %d..%d seq %d+%d st=%d en=%d %.*s\n", 
	       s->name_len, s->name,
	       r[i].start, r[i].end,
	       s->pos, s->len,
	       s->left, s->right,
	       ABS(s->len), s->seq);
    }

    gio_close(io);
    system("ps lx | grep g_iotest | grep -v grep");
    exit(0);
}

#define CONS_LEN 16384
static void test_mode2(GapIO *io, contig_t **c, int xpos) {
    rangec_t *r;
    int nr, i, bpv = 256;
    char cons[CONS_LEN+1];
    track_t *t;

    //    r = contig_seqs_in_range(io, c, xpos, xpos+CONS_LEN, &nr);
    //    qsort(r, nr, sizeof(*r), sort_range);
    //    calc_cons(io, r, nr, xpos, CONS_LEN, cons);
    //    printf("Cons=%.*s\n", CONS_LEN, cons);

    t = contig_get_track(io, c, xpos, xpos+CONS_LEN, TRACK_READ_DEPTH, bpv);
    for (i = 0; i < CONS_LEN/bpv; i++) {
	printf("%d\t%d\n", i*bpv, arr(int, t->data, i));
    }

    cache_flush(io);
    gio_close(io);
    exit(0);
}

static void benchmark(GapIO *io, contig_t **c) {
    int i;
    char cons[10000];

    srandom(0);
    fprintf(stderr, "=== Benchmarking ===\n");
    for (i = 0; i < 1000; i++) {
	int xpos = random() % 2000000;
	int size = random() % 1000;
	int nr;
	rangec_t *r;

	r = contig_seqs_in_range(io, c, xpos, xpos+size, 0, &nr);
	calc_cons(io, r, nr, xpos, size, cons);
	printf("%.*s\n", size, cons);
	fputc('.', stderr);
	fflush(stderr);
	free(r);
    }
    gio_close(io);
    exit(0);
}

static void test_mode3(GapIO *io, int cnum, int xpos) {
    rangec_t *r;
    contig_iterator *ci;

    ci = contig_iter_new(io, cnum, 0, CITER_FIRST, CITER_CSTART, CITER_CEND);
    while (r = contig_iter_next(io, ci)) {
	seq_t *s = get_seq(io, r->rec);
	char name[256];

	sprintf(name, "%.*s", s->name_len, s->name);
	printf("%c%-22s\t%8d..%-8d\t%.*s\n",
	       "+-"[s->len<0], name, r->start, r->end, ABS(s->len), s->seq);
    }
    contig_iter_del(ci);
    exit(0);
}
#endif
/* ------------------------------------------------------------------------ */
void usage(void) {
    fprintf(stderr, "Usage: tg_view [options] dbname [position]\n");
    fprintf(stderr, "\t-h        This help\n");
    fprintf(stderr, "\t-d        Start in highlight-disagreements mode\n");
    fprintf(stderr, "\t-C        Start with cutoffs hidden\n");
    fprintf(stderr, "\t-c        Start with cutoffs shown (default)\n");
    fprintf(stderr, "\t-l width  Output non-interactively to stdout\n");
    fprintf(stderr, "\t-e        Permit editing (default is read-only)\n");
}


int main(int argc, char **argv) {
    GapIO *io;
    int xpos = 0;
    int opt;
    int lp_mode = 0;
    int mode = DISPLAY_QUAL | DISPLAY_CUTOFFS;
    extern char *optarg;
    contig_t *c;
    int cnum = 0;
    int read_only = 1;

    while ((opt = getopt(argc, argv, "hl:dcCx:e")) != -1) {
	switch (opt) {
	case '?':
	case 'h':
	    usage();
	    return 0;

	case 'd':
	    mode |= DISPLAY_DIFFS;
	    break;

	case 'c':
	    mode |= DISPLAY_CUTOFFS;
	    break;

	case 'C':
	    mode &= ~DISPLAY_CUTOFFS;
	    break;

	case 'l':
	    lp_mode = atoi(optarg);
	    break;

	case 'x':
	    cnum = atoi(optarg)-1;
	    break;

	case 'e':
	    read_only = 0;
	    break;

	default:
	    if (opt == ':')
		fprintf(stderr, "Missing parameter\n");
	    else
		fprintf(stderr, "Unknown option '%c'\n", opt);
	    usage();
	    return 1;
	}
    }
    
    if (optind == argc) {
	usage();
	return 1;
    }

    if (NULL == (io = gio_open(argv[optind], read_only, 0))) {
	fprintf(stderr, "Unable to open db: %s\n", argv[1]);
	return 1;
    }
    optind++;

    if (optind != argc) {
	xpos = atoi(argv[optind]);
    }

    io->contig_num = cnum;
    gio_read_contig(io, cnum, &c); 
    cache_incr(io, c);

#ifdef TEST_MODE
    //test_mode(io, &c, xpos); 
    test_mode2(io, &c, xpos);
    //benchmark(io, &c);
    //test_mode3(io, arr(GCardinal, io->contig_order, cnum), xpos);
#endif

    if (lp_mode) {
	print_output(io, &c, xpos, lp_mode, mode);
	gio_close(io);
    } else {
	init_curses();
	curses_loop(io, &c, xpos, mode);
	endwin();
	
	if (io->cache && io->debug_level > 0) {
	    fputs("\n=== cache ===", stderr);
	    HacheTableStats(io->cache, stderr);
	}

	gio_close(io);
    }

    if (!lp_mode) {
	printf("\n\n\tg_view:\tShort Read Alignment Viewer, version 1.2.11"SVN_VERS"\n");
	printf("\n\tAuthor:\tJames Bonfield (jkb@sanger.ac.uk)\n");
	printf("\t\t2007-2011, Wellcome Trust Sanger Institute\n\n");
    }

    return 0;
}
