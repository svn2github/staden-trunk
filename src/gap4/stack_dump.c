#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <signal.h>

#include "misc.h"
#include "gap-error.h"
#include "g-error.h"
#include "array.h"
#include "bitmap.h"
#include "stack_dump.h"

char *prog_name = NULL; /* Needs assigning to argv[0] in main() */

/*
 * ---------------------------------------------------------------------------
 * DEC Alpha specific code
 * ---------------------------------------------------------------------------
 */
#if defined(__alpha____NOT_WORKING___)
#define DO_STACK_TRACE

/*
 * Dec Alpha stack tracing system.
 * It's very hacky and clumsy, often reporting more stack frames than exist,
 * but it does provide some information (ie better than nothing).
 *
 * We only load the symbol table for the main binary and not for the libraries.
 * (How do we know which libraries we've used?). Hence it works OK until we
 * get to the first library call, at which stage we revert to checking each
 * long word.
 */

#include <stdlib.h>
#include <signal.h>
#include <fcntl.h>
#include <c_asm.h>
#include <filehdr.h>
#include <syms.h>
#include <ldfcn.h>

#ifdef OLD_STACK_PRINTING_CODE
/*
 * Converts a pointer into a function name and an offset from the start of
 * this function. Returns NULL if it doesn't know.
 * Wildly inefficient - we ought to use a binary search.
 */
static char *get_name(unsigned long laddr, signed long *diff) {
    static LDFILE *ldp = NULL;
    unsigned long i;
    SYMR symbol;
    PDR procp;

    /* Read the LDFILE structure */
    if (!ldp) {
	if (NULL == (ldp = ldopen(prog_name, NULL))) {
	    return NULL;
	}
    }

    /* Maybe it's stripped */
    if (!SYMTAB(ldp))
	return NULL;

    /* Loop through the procedure data to find the correct address */
    for (i = 0; i < SYMHEADER(ldp).ipdMax; i++) {
	if (FAILURE == ldgetpd(ldp, i, &procp))
	    continue;
	if (laddr < procp.adr)
	    break;
    }

    if (--i < 0)
	return NULL;

    /* Read the symbol entry and then the name of this procedure */
    ldgetpd(ldp, i, &procp);
    *diff = laddr - procp.adr;
    if (FAILURE == ldtbread(ldp, procp.isym, &symbol))
	return NULL;

    return ldgetname(ldp, &symbol);
}

void stack_trace(void) {
    extern unsigned long __start;
    unsigned long *sp = (unsigned long *)asm("bis %sp,%sp,%v0");
    char *cp;
    int frame = 0, i;
    signed long diff;

    printf("Estimated stack trace (0x%lx) follows\n", sp);
    do {
	if (*sp >= 0x120000000 && *sp < 0x140000000) { /* HACK!!! */
	    if (cp = get_name(*sp, &diff)) {
		if (diff >= 0) {
		    verror(ERR_FATAL, "stack", "%2d: 0x%lx (%s + 0x%lx)\n",
			   frame++, *sp, cp, diff);
		} else {
		    /* printf("unknown frame at 0x%lx\n", *sp); */
		}
	    } else {
		verror(ERR_FATAL, "%2d: 0x%lx\n", frame++, *sp);
	    }
	}
	sp+=2;
    } while (*sp != __start + 0xe0 && frame < 20);
}
#endif

/*
 * Converts a pointer into a function name and an offset from the start of
 * this function. Returns NULL if it doesn't know.
 * Wildly inefficient - we ought to use a binary search.
 */
static char *get_name(unsigned long laddr, unsigned long *diff,
		      int *frameoffset) {
    unsigned long i;
    SYMR symbol;
    PDR procp;
    char *cp;
    static LDFILE *ldp = NULL;

    /* Read the LDFILE structure */
    if (!ldp) {
	if (NULL == (ldp = ldopen(prog_name, NULL))) {
	    return NULL;
	}
    }

    /* Maybe it's stripped */
    if (!SYMTAB(ldp))
	return NULL;

    /* Loop through the procedure data to find the correct address */
    for (i = 0; i < SYMHEADER(ldp).ipdMax; i++) {
	if (FAILURE == ldgetpd(ldp, i, &procp))
	    continue;
	if (laddr < procp.adr)
	    break;
    }

    if (--i < 0)
	return NULL;

    /* Read the symbol entry and then the name of this procedure */
    ldgetpd(ldp, i, &procp);
    *diff = laddr - procp.adr;
    if (FAILURE == ldtbread(ldp, procp.isym, &symbol))
	return NULL;

    *frameoffset = procp.frameoffset;

    cp = ldgetname(ldp, &symbol);
/*    printf("Function=%s frameoffset=%d, regoff=%d, fregoff=%d\n",
	   cp, procp.frameoffset, procp.regoffset, procp.fregoffset); */
    return cp;
}

static void stack_trace2() {
    unsigned long *sp = (unsigned long *)asm("bis %sp,%sp,%v0");
    extern unsigned long __start;
    char *cp;
    int frame = 0, i, j;
    unsigned long diff;
    int frameoffset, next_frameoffset;
    unsigned long *x;
    int counter = 0;
    char buf1[1024];

    printf("Estimated stack trace (0x%lx) follows\n", sp);

    /* Search for deadbeef long word */
    for (; *sp != 0xdeadbeefdeadbeef; sp++);
    sp--;

    /*
     * We know that this frame is for stack_trace(), which has a
     * frameoffset of 16
     */
    frameoffset = 16;

    for (i = 0; i < 200; i+=4) {
	char buf2[1024];

	sprintf(buf1, "%lx", &sp[i]);
	sprintf(buf2, "%016lx %016lx %016lx %016lx",
		sp[i], sp[i+1], sp[i+2], sp[i+3]);
	verror(ERR_FATAL, buf1, buf2);
    }

    do {
	sprintf(buf1, "%2d(%2d)", frame, counter);

	/* First word in frame is the return address. */
	if ((cp = get_name(*sp, &diff, &next_frameoffset)) &&
	    strcmp(cp, "_call_remove_gp_range") && diff <= 0x10000) {
	    verror(ERR_FATAL, buf1, "sp=0x%lx 0x%lx (%s + 0x%lx)\n",
		   sp, *sp, cp, diff);
	    frame++;
	} else {
	    frameoffset = 8; /* try again at next long word */
	    /* printf("unknown frame at 0x%lx\n", *sp); */
	}

	/* Then increment sp by frameoffset */
	sp += frameoffset / sizeof(long);
	frameoffset = next_frameoffset;
    } while ((strcmp(cp, "__start") != 0 || diff > 0xfff) && ++counter < 100);
}

void stack_trace(void) {
    unsigned long find_me = 0xdeadbeefdeadbeef;
    stack_trace2();
}

/*
 * ---------------------------------------------------------------------------
 * Solaris specific code
 * ---------------------------------------------------------------------------
 */
#elif (defined(__sun__) || defined(__sun)) && (defined(__svr4__) || defined(__SVR4))
#define DO_STACK_TRACE

/*
 * Solaris stack tracing system.
 * Very easy to do - just call the pstack program to do it all for us!
 */

#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include "text_output.h"

void stack_trace(void) {
    char cmd[1024], *cp;
    char buf[BUFSIZ];
    FILE *fp;

    sprintf(cmd, "/usr/proc/bin/pstack %d\n", (int)getpid());
    if ((fp = popen(cmd, "r")) != NULL)
	while (fgets(buf, BUFSIZ, fp) != NULL) {
	    if (cp = strchr(buf, '\n'))
		*cp = 0;
	    verror(ERR_FATAL, "stack", buf);
	}
    fclose(fp);
}

/*
 * ---------------------------------------------------------------------------
 * "Anything else" specific code
 * ---------------------------------------------------------------------------
 */
#else

void stack_trace(void) {
}

#endif
