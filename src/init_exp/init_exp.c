/*
 * Title: init_exp
 *
 * Based upon getMCH.c
 */

#include <stdio.h>
#include <string.h>

#include "Read.h"
#include "expFileIO.h"
#include "errno.h"
#include "misc.h"
#include "open_trace_file.h"
#include "mFILE.h"

double avg_qual(Read *r) {
    double aq = 0.0;
    int i;
    int i_start = 100, i_end = 200;

    if (i_start >= r->NBases)
        return 0.0;

    if (i_end > r->NBases)
        i_end = r->NBases;

    for (i = i_start ; i < i_end; i++) {
        switch(r->base[i]) {
        case 'a':
        case 'A':
            aq += r->prob_A[i];
            break;

        case 'c':
        case 'C':
            aq += r->prob_C[i];
            break;

        case 'g':
        case 'G':
            aq += r->prob_G[i];
            break;

        case 't':
        case 'T':
            aq += r->prob_T[i];
            break;
        }
    }

    aq /= (i_end - i_start);

    return aq;
}

int convert(char *file, int format, mFILE *ofp, char *name, int output_conf) {
    Read *r;
    Exp_info *e;
    char buf[50];
    double aq;

    if (format == TT_BIO) {
        if (NULL == (r = read_reading(file, format))) {
            fprintf(stderr, "%s: failed to read\n", file);
            return 1;
        }
    } else {
        FILE *infp;
        if (NULL == (infp = open_trace_file(file, NULL))) {
            perror(file);
            return 1;
        }
        if (NULL == (r = fread_reading(infp, file, format))) {
            fprintf(stderr, "%s: failed to read\n", file);
            return 1;
        }
        fclose(infp);
    }

    e = read2exp(r, name);
    if (NULL == e) {
        fprintf(stderr, "Failed to create experiment file.\n");
        read_deallocate(r);
        return 1;
    }

    sprintf(buf, "%f", aq = avg_qual(r));
    exp_set_entry(e, EFLT_AQ, buf);
    exp_print_mfile(ofp, e);

    if (output_conf && aq != 0) {
        char *cstr;
        int1 *conf;
        int i;

        conf = xmalloc(r->NBases * sizeof(*conf));
        cstr = xmalloc(5 * r->NBases+2);
        for (i = 0; i < r->NBases; i++) {
            switch (r->base[i]) {
            case 'a':
            case 'A':
                conf[i] = r->prob_A[i];
                break;
            case 'c':
            case 'C':
                conf[i] = r->prob_C[i];
                break;
            case 'g':
            case 'G':
                conf[i] = r->prob_G[i];
                break;
            case 't':
            case 'T':
                conf[i] = r->prob_T[i];
                break;
            default:
                conf[i] = (r->prob_A[i] +
                           r->prob_C[i] +
                           r->prob_G[i] +
                           r->prob_T[i]) / 4;
                break;
            }
        }

        conf2str(conf, r->NBases, cstr);
        exp_set_entry(e, EFLT_AV, cstr);

        xfree(cstr);
        xfree(conf);
    }

    read_deallocate(r);
    exp_destroy_info(e);

    mfflush(ofp);

    return 0;
}


void usage(void) {
    fprintf(stderr,
            "usage: init_exp [-(abi|alf|scf|pln)] [-output file] "
            "[-name entry_name] [-conf] file\n");
    exit(1);
}

int main(int argc, char **argv) {
    int a = 1;
    int format = TT_ANY;
    mFILE *ofp = mstdout();
    char *name = NULL, *oname = NULL;
    int output_conf = 0;

    read_sections(READ_BASES);

    while (a < argc) {
        if (strcasecmp(argv[a], "-abi") == 0)
            format = TT_ABI;

        else if (strcasecmp(argv[a], "-alf") == 0)
            format = TT_ALF;

        else if (strcasecmp(argv[a], "-scf") == 0)
            format = TT_SCF;

        else if (strcasecmp(argv[a], "-pln") == 0)
            format = TT_PLN;

        else if (strcasecmp(argv[a], "-name") == 0) {
            if (a == argc)
                usage();

            name = argv[++a];

        } else if (strcasecmp(argv[a], "-output") == 0) {
            if (a == argc)
                usage();

            ofp = mfopen(oname = argv[++a], "wb");
            if (NULL == ofp) {
                perror(argv[a]);
                return 1;
            }

        } else if (strcasecmp(argv[a], "-conf") == 0) {
            output_conf = 1;

        } else if (argv[a][0] == '-')
            usage();
        else
            break;

        a++;
    }

    if (a+1 != argc)
        usage();

    if (!name)
        name = oname;
    if (!name)
        name = argv[a];

    return convert(argv[a], format, ofp, name, output_conf);
}
