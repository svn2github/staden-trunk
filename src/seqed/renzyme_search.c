#include <tcl.h>
#include <tk.h>
#include <itcl.h>
#include <itk.h>

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "misc.h"
#include "tcl_utils.h"
#include "renzyme_box.h"
#include "renzyme_search.h"
#include "graphic_editor.h"
#include "renz_utils.h"  /*hash_dna & dna_search*/
#include "editor_reg.h"
#include "parse_feature.h"
#include "feature_table.h"
#include "read_sequence.h"
#include "editor.h"
#include "dna_utils.h"
#include "tkSeqed.h"

#define END_LEN 6

/* initialise a SITE structure */
SITE *init_site (void) {

    SITE *s;

    if ( NULL == ( s = (SITE* ) xmalloc (sizeof (SITE))))
        return NULL;
    s->pos1 = 0;
    s->pos2 = 0;
    s->name = NULL;
    s->rec_seq = NULL;
    s->direction = 0;
    return s;
}

/* free a SITE structure */
void free_site (SITE *s) {

  if (s) {
    if (s->name) xfree(s->name);
    if (s->rec_seq) xfree(s->rec_seq);
    xfree (s);
  }
}

/* initialise a SITES structure */
SITES *init_sites (void) {

    SITES *ss;
    SITE  **s;
    if ( NULL == (ss = (SITES *) xmalloc (sizeof (SITES))))
        goto err;;
    if ( NULL == (s = (SITE **) xmalloc (sizeof (SITE*))))
        goto err;

    ss->site = s;
    ss->used = 0;
    ss->capacity = 0;
    return  ss;
 err:
    if (ss) xfree(ss);
    if (s)  xfree (s);
    return NULL;
}

/* free a SITES structure */
void free_sites ( SITES *ss ) {

    int i;
    if (ss) {
        for ( i = 0; i <= ss->used; i++ ) {
            free_site (ss->site[i]);
        }
    }
    xfree (ss);
}

/* realloc the site array in the sites structure */
SITES *realloc_sites ( SITES *ss, int num_site) {

    SITE  **s;
    s = ss->site;
    if ((NULL == ( s = (SITE **)xrealloc (s, sizeof(SITE *)*(num_site+1)))))
        return NULL;
    ss->site = s;
    ss->used = num_site;
    return ss;
}




int save_fragment_in_sequence (SEQUENCE *s, char *seq) {

    SEQUENCE *sf;
    SITES *sites;
    char *tmp;
    int num_site, s1, s2;
    int type, seq_len;
    int site1, site2;

    if (s->sites == NULL) goto err;
    sites = s->sites;
    num_site = sites->used;
    type = s->type;
    seq_len = s->length;

    /*printf ("seq_len=%d\n", seq_len);*/
    type = 2;
    site1 = 1;
    site2 = 0;

    if ( !(sf = init_sequence()))
        goto err;
    s1 = sites->site[site1]->pos1;
    s2 = sites->site[site2]->pos1;
    if (sites->site[site1]->pos2 != 0)
        sf->left_end->len = sites->site[site1]->pos2 - sites->site[site1]->pos1;
    if (sites->site[site2]->pos2 != 0)
        sf->right_end->len = sites->site[site2]->pos1 - sites->site[site2]->pos2;
    if (site1 == num_site - 1)
        sf->length = seq_len - s1 + s2 + ABS(sf->right_end->len);
    else sf->length = s2 - s1 + ABS(sf->right_end->len);
    /*printf("site1=%d\n", s1);
    printf("site2=%d\n", s2);
    printf("sf->left_end_overhang=%d\n", sf->left_end_overhang);
    printf("sf->right_end_overhang=%d\n", sf->right_end_overhang);
    printf("sf->length=%d\n", sf->length);*/

    if (NULL == (tmp = (char *)malloc((sf->length + 1) * sizeof(char))))
        goto err;

    if (site1 == num_site - 1 && (type == 2 || type == 4)) {
        memmove (&tmp[0], &s->seq[s1], seq_len-s1);
        memmove (&tmp[seq_len-s1], &s->seq[0], s2+ABS(sf->right_end->len));
    } else
        strncpy (tmp, &seq[s1], sf->length);
    tmp[sf->length] = 0;
    sf->seq = tmp;

    s->fragment = sf;
    xfree (tmp);
    return 0;
 err:
    if (sf) xfree (sf);
    if (tmp) xfree (tmp);
    return -1;
}



RENZYMES *get_selected_renzyme (int num_sel, char **sel) {

    RENZYMES *sel_rs;
    int i, j, num_ren;
    int sel_id = 0;

    num_ren = renzymes->used;
    sel_rs = init_renzymes ();
    for (i = 0; i < num_ren; i++) {
        for ( j = 0; j < num_sel; j++) {

            if (!strcmp ( renzymes->renzyme[i]->name, sel[j]) ) {
                if (sel_id != 0) {
                    realloc_renzymes (sel_rs, sel_id);
                }
                sel_rs->renzyme[sel_id] = renzyme_copy (renzymes->renzyme[i]);
                sel_id++;
                sel_rs->used = sel_id;
            }
        }
    }
    return sel_rs;
}

int max_rec_seq_length (RENZYMES *sel_rs) {

    int i, num_ren;
    int rec_seq_len = 0;

    num_ren = sel_rs->used;
    for ( i = 0; i < num_ren; i++) {
        rec_seq_len = MAX (rec_seq_len, strlen(sel_rs->renzyme[i]->rec_seq));
    }
    return rec_seq_len;
}

int RenzSearchSelect(ClientData clientData, Tcl_Interp *interp, int argc, char **argv) {

    EDITOR_RECORD *er;
    RENZYMES *sel_rs;
    int num_sel;
    char **sel = NULL;
    int i, num_ren;
    int num_seq, k;
    int ed_id;
    char *seq, *subseq, *tmp;
    int seq_len;
    int seq_type;
    int total_matches, max_rec_length;
    char buf[1024];
    int group_id;

    if (argc != 3) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "ed_id\"", " selected items\"",(char*)NULL);
        return TCL_ERROR;
    }
    ed_id = atoi (argv[1]);

    /* create selected Renzyme name array */
    if (Tcl_SplitList(interp, argv[2], &num_sel, &sel) != TCL_OK)
        goto err;
    sel_rs = get_selected_renzyme ( num_sel, sel);
    num_ren = renzymes->used;

    er = GetEditor (ed_id);
    group_id = 0;
    num_seq = er->seq_group[group_id]->nmembers;
    max_rec_length = max_rec_seq_length (sel_rs);

    /*printf("sel_rs->used=%d\n", sel_rs->used);*/
    for (i = 0; i < sel_rs->used; i++) {
        Tcl_DString dstr, dstr1;
        Tcl_DStringInit(&dstr);
        Tcl_DStringInit(&dstr1);
        for (k = 1; k <= num_seq; k++) {
            R_MATCH *match;
            if (NULL == (match = (R_MATCH*)xcalloc(MAXMATCHES, sizeof(R_MATCH))))
                goto err;
            seq = er->seq_group[group_id]->members[k]->data->sequence->seq;
            seq_len = er->seq_group[group_id]->members[k]->data->sequence->length;
            seq_type = er->seq_group[group_id]->members[k]->data->sequence->type;
            /*seq_type = 1;*/
            seq_type = 2;

            if (seq_type == 2 || seq_type == 4) {
                char *subseq, *tmp;

                if (NULL == (subseq = (char *)malloc((seq_len + max_rec_length + 1 ) * sizeof(char))))
                    goto err;
                if (NULL == (tmp = (char *)malloc((max_rec_length + 1) * sizeof(char))))
                    goto err;
                strcpy (subseq, seq);
                subseq[seq_len] = 0;
                strncpy (tmp, &seq[0], max_rec_length - 1);
                tmp[max_rec_length - 1] = 0;
                strcat (subseq, tmp );
                seq = subseq;
                seq_len = seq_len + max_rec_length;
            }
            find_matches(sel_rs->renzyme[i], seq, seq_len, seq_type, i, &match, &total_matches);
            if (seq_type == 2 || seq_type == 4) {
              /*
                xfree (tmp);
                xfree (subseq);
              */
            }

            sprintf(buf, "%d", total_matches);
            Tcl_DStringAppendElement(&dstr, buf);
            xfree (match);
        }

        vTcl_DStringAppendElement(&dstr1, "%s", sel_rs->renzyme[i]->name);
        vTcl_DStringAppendElement(&dstr1, "%s", sel_rs->renzyme[i]->rec_seq_text);
        vTcl_DStringAppendElement(&dstr1, "%s", sel_rs->renzyme[i]->prototype);
        vTcl_DStringAppendElement(&dstr1, "%s", sel_rs->renzyme[i]->supplier_codes);
        vTcl_DStringAppendElement(&dstr1, "%.f", sel_rs->renzyme[i]->av_frag_size);
        vTcl_DStringAppendElement(&dstr1, "%s", Tcl_DStringValue(&dstr));
        Tcl_AppendElement(interp, Tcl_DStringValue(&dstr1));
        Tcl_DStringFree(&dstr);
        Tcl_DStringFree(&dstr1);
    }
    ckfree((char*)sel);
    return TCL_OK;
 err:
    if (sel) ckfree ((char*)sel);
    if (subseq) xfree (subseq);
    if (tmp) xfree (tmp);
    return TCL_OK;
}

int RenzymeSearchPreview (ClientData clientData, Tcl_Interp *interp, int argc, char **argv) {

    RENZYMES *sel_rs;
    int num_sel_ren, num_sel_seq;
    char **sel_ren = NULL;
    char **sel_seq = NULL;

    int i, k, num_ren;
    char *seq, *subseq, *tmp;
    int seq_len, seq_type;
    int total_matches, max_rec_length;
    char buf[1024];

    if (argc != 3) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], "selected sequence\"", " selected REnzyme\"",(char*)NULL);
        return TCL_ERROR;
    }

    /* create selected Renzyme name array */
    if (Tcl_SplitList(interp, argv[2], &num_sel_ren, &sel_ren) != TCL_OK)
        goto err;
    sel_rs = get_selected_renzyme (num_sel_ren, sel_ren);
    num_ren = renzymes->used;
    max_rec_length = max_rec_seq_length (sel_rs);

    /* create selected sequence identifier array */
    if (Tcl_SplitList (interp, argv[1], &num_sel_seq, &sel_seq) != TCL_OK)
        goto err;

    for (i = 0; i < sel_rs->used; i++) {
        Tcl_DString dstr, dstr1;
        Tcl_DStringInit(&dstr);
        Tcl_DStringInit(&dstr1);

        for (k = 0; k < num_sel_seq; k++) {
            R_MATCH *match;
            int seq_id; /* unique id for sequences */
            int ee_id;  /* unique id for editor and enzyme_map_editor */

            if (NULL == (match = (R_MATCH*)xcalloc(MAXMATCHES, sizeof(R_MATCH))))
                goto err;

            /* get id fromsequences */
            seq_id = GetSequenceIdByName (sel_seq[k]);
            /* make a copy in eden */
            ee_id = create_copy_for_editor (interp, seq_id);
            seq = GetEdenSeq (ee_id);
            seq_len = GetEdenLength (ee_id);
            seq_type = GetEdenType (ee_id);
            seq_type = 2;

            if (seq_type == 2 || seq_type == 4) {
                char *subseq, *tmp;

                if (NULL == (subseq = (char *)malloc((seq_len + max_rec_length + 1 ) * sizeof(char))))
                    goto err;
                if (NULL == (tmp = (char *)malloc((max_rec_length + 1) * sizeof(char))))
                    goto err;
                strcpy (subseq, seq);
                subseq[seq_len] = 0;
                strncpy (tmp, &seq[0], max_rec_length - 1);
                tmp[max_rec_length - 1] = 0;
                strcat (subseq, tmp );
                /*seq = subseq;*/
                seq_len = seq_len + max_rec_length;
                find_matches(sel_rs->renzyme[i], subseq, seq_len, seq_type, i, &match, &total_matches);
                xfree (tmp);
                xfree (subseq);
            } else
                find_matches(sel_rs->renzyme[i], seq, seq_len, seq_type, i, &match, &total_matches);
            sprintf(buf, "%d", total_matches);
            Tcl_DStringAppendElement(&dstr, buf);
            xfree (match);
        }

        vTcl_DStringAppendElement(&dstr1, "%s", sel_rs->renzyme[i]->name);
        vTcl_DStringAppendElement(&dstr1, "%s", sel_rs->renzyme[i]->rec_seq_text);
        vTcl_DStringAppendElement(&dstr1, "%s", sel_rs->renzyme[i]->prototype);
        vTcl_DStringAppendElement(&dstr1, "%s", sel_rs->renzyme[i]->supplier_codes);
        vTcl_DStringAppendElement(&dstr1, "%.f", sel_rs->renzyme[i]->av_frag_size);
        vTcl_DStringAppendElement(&dstr1, "%s", Tcl_DStringValue(&dstr));
        Tcl_AppendElement(interp, Tcl_DStringValue(&dstr1));
        Tcl_DStringFree(&dstr);
        Tcl_DStringFree(&dstr1);
    }
    ckfree((char*)sel_ren);
    ckfree((char*)sel_seq);
    return TCL_OK;
 err:
    if (sel_ren) ckfree ((char*)sel_ren);
    if (sel_seq) ckfree((char*)sel_seq);
    if (subseq) xfree (subseq);
    if (tmp) xfree (tmp);
    return TCL_OK;
}

void dis_fragment ( SEQUENCE *sequence) {

    int i, num_entry;

    printf("sequence %s start %d end %d length %d\n",
           sequence->seq,sequence->start,sequence->end, sequence->length);
    if (sequence->left_end) printf("left_end_hang=%s\n", sequence->left_end->hang);
    if (sequence->right_end) printf("right_end_hang=%s\n", sequence->right_end->hang);
    if(sequence->feature_table) {
        num_entry = sequence->feature_table->num_entry;
        for (i = 0; i < num_entry; i++) {
            print_entry(sequence->feature_table->entry[i]);
        }
    }
}
int dis_sequence ( SEQUENCE *sequence) {

    int i, num_entry;

    printf("sequence %s start %d end %d length %d\n",
           sequence->seq,sequence->start,sequence->end, sequence->length);
    printf("cut_site_1 %d overhang %s cut_site_2 %d overhang %s\n",
            sequence->cut_site_1->pos, sequence->cut_site_1->hang,
            sequence->cut_site_2->pos, sequence->cut_site_2->hang);
    if(sequence->feature_table) {
        num_entry = sequence->feature_table->num_entry;
        for (i = 0; i < num_entry; i++) {
            print_entry(sequence->feature_table->entry[i]);
        }
    }
    return 0;

}

int get_left_over_hang ( char *renz_start ) {

    int i, num_ren;
    int oh = 0;

    num_ren = renzymes->used;

    for ( i = 0; i < num_ren; i++) {
        if (!strcmp (renzymes->renzyme[i]->name, renz_start)) {
          if (renzymes->renzyme[i]->cut_pos_2 != 0)
              oh = strlen (renzymes->renzyme[i]->rec_seq)
                  - renzymes->renzyme[i]->cut_pos_2
                  - renzymes->renzyme[i]->cut_pos_1;
          else oh = strlen (renzymes->renzyme[i]->rec_seq)
                  - renzymes->renzyme[i]->cut_pos_1
                  - renzymes->renzyme[i]->cut_pos_1;
        }
    }
    return oh;
}

int get_right_over_hang ( char *renz_end ) {

    int i, num_ren;
    int oh = 0;

    num_ren = renzymes->used;
    for ( i = 0; i < num_ren; i++) {
        if (!strcmp (renzymes->renzyme[i]->name, renz_end)) {
          if (renzymes->renzyme[i]->cut_pos_2 != 0)
              oh = renzymes->renzyme[i]->cut_pos_1 + renzymes->renzyme[i]->cut_pos_2
                  - strlen (renzymes->renzyme[i]->rec_seq);
          else oh = renzymes->renzyme[i]->cut_pos_1 + renzymes->renzyme[i]->cut_pos_1
                  - strlen (renzymes->renzyme[i]->rec_seq);
        }
    }
    return oh;
}

char *get_hang (SEQUENCE *s, int start, int overhang) {

    char *hang;

    if ((NULL == (hang = (char * ) xmalloc ( overhang + 1) ))) return NULL;
    strncpy (hang, &s->seq[start], overhang);
    hang[overhang] = 0;
    return hang;
}

void set_ends (SEQUENCE *s, SEQUENCE *f, int start, int end, char *renz_start, char *renz_end) {

    SITE_HANG *sh_left;
    SITE_HANG *sh_right;
    int oh;

    /* set fragment left_end */
    oh = get_left_over_hang (renz_start);
    sh_left = init_site_hang ();
    sh_left->pos = start;
    sh_left->len = oh;

    if (sh_left->len == 0) {
      sh_left->hang = NULL;
    }
    /* acgt|ACGTACGTACGT */
    /*     ------        */
    /* tgcatgcat|GCATGCA */
    if (sh_left->len > 0) {
      sh_left->hang = get_hang (s, start, sh_left->len);
    }
    /* acgtacgtacgt|ACGTACGT */
    /*        ------         */
    /* tgcatgc|ATGCATGCA     */
    if (sh_left->len < 0) {
      int h = ABS(sh_left->len);
      sh_left->hang = get_hang (s, start-h, h);
    }
    f->left_end = sh_left;
    /* set fragment right_end */
    oh = get_right_over_hang (renz_end);

    sh_right = init_site_hang ();
    sh_right->len = oh;
    sh_right->pos = end;
    if (sh_right->len == 0) {
      sh_right->hang = NULL;
    }

    /* ACGTACGT|acgtacgt */
    /*         -----     */
    /* TGCATGCATGCA|tgca */
    if (sh_right->len < 0) {
      int h = ABS(sh_right->len);
      sh_right->hang = get_hang (s, end, h);
    }

    /* ACGTACGTACGTACGTACGT|acgtacgt */
    /*              --------         */
    /* TGCATGCATGCA|tgcagct          */
    if (sh_right->len > 0) {
      sh_right->hang = get_hang (s, end-sh_right->len, sh_right->len);
    }
    f->right_end = sh_right;
    /*printf ("create_fragment_from_sequence:\n");
      dis_fragment(f);*/
}

SEQUENCE *create_fragment_from_sequence (int seq_num, int start, int end, char *renz_start, char *renz_end) {

    SEQUENCE *s, *f;
    FEATURE_TABLE *ft;
    int seq_id;
    int frag_len;

    frag_len = end - start;
    if (NULL == (f = init_sequence())) goto err;
    if ((NULL == (f->seq = (char * ) xmalloc ( frag_len + 1) ))) goto err;

    if (seq_num != -1) {
        seq_id = GetEdenId (seq_num);
        s = GetEdenSequence (seq_num);
    }
    if (seq_num == -1) {
        seq_id = -1;
        s = get_editor_buffer ();
    }

    strncpy (f->seq, &s->seq[start], frag_len);
    f->seq[frag_len] = 0;
    f->start = 1;
    f->end = frag_len;
    f->length = frag_len;
    f->parent_id = seq_id;
    ft = get_feature_table (s, frag_len, start);
    if (ft) {
      change_feature_table_location ( ft, -(start) );
      f->feature_table = ft;
    }
    if (renz_start != NULL && renz_end != NULL) {
        set_ends (s, f, start, end, renz_start, renz_end);
    }
    return f;

    err:
    if (f) free_sequence (f);
    if (f->seq) xfree (f->seq);
    return NULL;
}

int ismatch (SITE_HANG *site1, SITE_HANG *site2, SITE_HANG *left,  SITE_HANG *right ) {

    int m;

    /* check the cut_site_1 of sequence and the left_end of fragment */
    if (site1 != NULL && site1->len != 0 && left != NULL && left->len != 0) {
        m = site1->len + left->len;
        if (!m && !strcmp (site1->hang, left->hang) ) {
            /* check the cut_site_2 of sequence and the right_end of fragment */
            if (site2 != NULL && site2->len != 0 && right != NULL && right->len != 0) {
              m = site2->len + right->len;
              if (!m && !strcmp (site2->hang, right->hang) ) return 1;
            }
            return -1;
        }
        return -1;
    }
    return -1;
}

/* return 0: both_end_match; (fragment in forward direction)
   return 1: both_end_match; (fragment in reverse direction)
   return 2: NO hangs at both end or matches happened at both direction;
   return 3: NO match at one of the both end in both direction */
int overhang_check (SEQUENCE *s, SEQUENCE *f) {

    int m;

    /* there are no hangs at both end */
    if (!s->cut_site_1 && !f->left_end && !s->cut_site_2 && !f->right_end) return 2;
    if (s->cut_site_1 != NULL && s->cut_site_1->len == 0
        && f->left_end != NULL && f->left_end->len == 0
        && s->cut_site_2 != NULL && s->cut_site_2->len == 0
        && f->right_end != NULL && f->right_end->len == 0) return 2;

    m = ismatch (s->cut_site_1, s->cut_site_2, f->left_end, f->right_end);
    if (m == 1) {
        /* complement fragment, continuing check */
        m = ismatch (s->cut_site_1, s->cut_site_2, f->right_end, f->left_end);
        if (m == 1) {
           return 2;
        }
        return 0;
    } else {
        m = ismatch (s->cut_site_1, s->cut_site_2, f->right_end, f->left_end);
        if (m == 1) return 1;
        /* return 3; */
    }
    return 3;
}

int insert_feature_table_in_sequence (SEQUENCE *s, FEATURE_TABLE *ft) {

    FEATURE_TABLE *f;
    int num_entry_f, num_entry_ft;
    int i;

    f = s->feature_table;
    num_entry_f = f->num_entry;
    num_entry_ft = ft->num_entry;

    for (i = 0; i < num_entry_ft; i++) {
        num_entry_f++;

        f = realloc_feature_table (f, num_entry_f);
        f->entry[num_entry_f - 1] = copy_ft_entry (ft->entry[i]);
        f->num_entry = num_entry_f;
    }
    return 0;
}

/* to break feature range into two parts when insertion is made in the feature range, */
/* and not to expand the range of the feature range, LINK??? YES!*/

int insert_fragment_modify_feature_table (SEQUENCE *s, SEQUENCE *f, int position) {

    int frag_len;
    int i, num_entry;
    ft_entry *e;
    ft_range *r, *rr;

    frag_len = f->length;
    num_entry = s->feature_table->num_entry;

    for (i = 0; i < num_entry; i++) {
        e = s->feature_table->entry[i];
        for (r = e->range; r; r = r->next) {
            if (position < r->left->min) {
                /*(a.b a.b) (a.b a.b) (a.b a.b)*/
                /*         ^                   */
                for (rr = r; rr; rr = rr->next) {
                    change_single_ft_range_insert (rr, frag_len);
                }
                break;
            } else if ((r->left->type == BASE || r->left->type == SIT) && position < r->left->max) {
                /*(a.b a.b)(a.b a.b) (a.b a.b)*/
                /*           ^                */
                ft_range *r_new;
                r_new = copy_ft_range (r);
                r_new->left->min = position + 1;
                for (rr = r_new; rr; rr = rr->next) {
                    change_single_ft_range_insert (rr, frag_len);
                }
                r->left->max = position;
                r->right = NULL;
                r->next = r_new;
                break;
            } else if (r->right && position < r->right->min) {
                /*(a.b a.b)(a.b a.b) (a.b a.b)*/
                /*             ^              */
                ft_range *r_new;
                r_new = copy_ft_range (r);
                r_new->left->min = position + 1;
                for (rr = r_new; rr; rr = rr->next) {
                    change_single_ft_range_insert (rr, frag_len);
                }
                r->right->min = position;
                r->next = r_new;
                break;
            } else if (r->right && (r->right->type == BASE || r->right->type == SIT)
                       && position < r->right->max) {
                /* (a.b a.b)(a.b a.b) (a.b a.b) */
                /*                ^             */
                ft_range *r_new;
                r_new = copy_ft_range (r);
                r_new->left->min = position + 1;
                r_new->left->max = r_new->right->max;
                r_new->right = NULL;
                for (rr = r_new; rr; rr = rr->next) {
                    change_single_ft_range_insert (rr, frag_len);
                }
                r->right->max = position;
                r->next = r_new;
                break;
            }
            /*break;*/
        }
    }
    return 0;
}

/*joint_pos jp 0: left  1: right */
char *get_joint_seq (int len, int ed_id, int member_id, SEQUENCE *f, int pos, int jp) {

    EDITOR_RECORD *er;
    char *s, *joint;
    int ll, group_id = 0;

    er = GetEditor (ed_id);
    s = er->seq_group[group_id]->members[member_id]->data->sequence->seq;
    len --;
    ll = 2*len;

    if (NULL == (joint = (char *)xmalloc( (ll + 1) * sizeof(char))))
        goto err;

    if (jp == 0) { /* left_joint */
        strncpy (&joint[0], &s[pos - len], len);
        strncpy (&joint[len], &f->seq[0], len);
        joint [ll] = 0;
    }
    if (jp == 1) { /* right_joint */
        int l = strlen (f->seq);
        strncpy (&joint[0], &f->seq[l - len], len);
        strncpy (&joint[len], &s[pos], len);
        joint [ll] = 0;
    }
    return joint;
 err: if (joint) xfree (joint);
    return NULL;
}

void cut_site_compare (R_MATCH *match_b, R_MATCH *match_l, int mb, int ml, int pos) {

    int i, j, k;

    k = 0;
    for (i = 0; i < mb; i++) {
        int flag = 0;
        for( j = 0; j < ml; j++) {
            if (match_b[i].enz_name == match_l[j].enz_name) {
                flag = 1;
                break;
            }
        }
        if (flag == 0) {
            k++;
            printf ("k=%d  renzymes->renzyme[%d]->name=%s   ",
                     k, match_b[i].enz_name,renzymes->renzyme[match_b[i].enz_name]->name);
                    printf ("cut_pos=%d\n", pos + match_b[i].cut_pos1);
        }
    }
}

int check_joint (int ed_id, int group_id, int member_id, SEQUENCE *f, int position) {

    EDITOR_RECORD *er;
    char *new_l, *new_r, *before;
    int rec_seq_len;
    int seq_type = 1;
    int ml, mr, mb;
    R_MATCH *match_l, *match_r, *match_b;
    int lb, pp, len;
    char *s;

    if (NULL == (match_l = (R_MATCH*)xcalloc(MAXMATCHES, sizeof(R_MATCH))))
                goto err;
    if (NULL == (match_b = (R_MATCH*)xcalloc(MAXMATCHES, sizeof(R_MATCH))))
                goto err;
    if (NULL == (match_r = (R_MATCH*)xcalloc(MAXMATCHES, sizeof(R_MATCH))))
                goto err;

    rec_seq_len = max_rec_seq_length (renzymes);
    pp = position - rec_seq_len;

    /* before insertion */
    lb = 2*(rec_seq_len - 1);
    if (NULL == (before = (char *)xmalloc( (lb + 1) * sizeof(char))))
        goto err;
    er = GetEditor (ed_id);
    s = er->seq_group[group_id]->members[member_id]->data->sequence->seq;
    strncpy (&before[0], &s[position - rec_seq_len + 1], lb);
    before[lb] = 0;
    find_matches_all (renzymes, before, lb, seq_type, &match_b, &mb);
    printf ("insertion_pos=%d\n", position);
    printf ("mb=%d\n", mb);
    printf ("before insertion: =%s\n", before);

    /* after insertion: LEFT*/
    new_l = get_joint_seq (rec_seq_len, ed_id, member_id, f, position, 0);
    len = strlen (new_l);
    find_matches_all (renzymes, new_l, len, seq_type, &match_l, &ml);
    printf ("ml=%d\n", ml);
    printf ("after insertion:  =%s\n", new_l);

    printf ("left joint LOST:\n");
    cut_site_compare (match_b, match_l, mb, ml, pp);
    printf ("left joint CREAT:\n");
    cut_site_compare (match_l, match_b, ml, mb, pp);

    /* after insertion: RIGHT*/
    new_r = get_joint_seq (rec_seq_len, ed_id, member_id, f, position, 1);
    len = strlen (new_r);
    find_matches_all (renzymes, new_r, len, seq_type, &match_r, &mr);
    printf ("mr=%d\n", mr);
    printf ("after insertion:  =%s\n", new_r);

    printf ("right joint LOST:\n");
    cut_site_compare (match_b, match_r, mb, mr, pp);
    printf ("righr joint CREAT:\n");
    pp = position + strlen(f->seq) - lb;
    cut_site_compare (match_r, match_b, mr, mb, pp);

    /*xfree (before);
    xfree (new_l);
    xfree (new_r);
    xfree (match_b);
    xfree (match_l);
    xfree (match_r);*/
    return 0;
 err:
    if (before) xfree (before);
    return -1;
}

int insert_fragment_in_sequence_1 (int seq_id, SEQUENCE *f, int position) {

    int ed_id, group_id, member_id;
    int err;

    ed_id = GetEdIdFromSeqId (seq_id);
    group_id = GetGroupIdFromSeqId (seq_id);
    member_id = GetMemberIdFromSeqId (seq_id);

    /*check_joint (ed_id, group_id, member_id, f, position);*/

    err = insert_update_sequence (ed_id, group_id, member_id, position, f);

    if (err != 0) return -1;
    return 0;
}

void complement_location (ft_location *l, int length) {

    int min = 0, max = 0;
    int min_lt;

    if (l->min_lt == 1) l->min_lt = -1;
    if (l->min_lt == -1) l->min_lt = 1;
    if (l->max_lt == 1) l->max_lt = -1;
    if (l->max_lt == -1) l->max_lt = 1;
    if (l->min) min = length - l->min + 1;
    if (l->max) max = length - l->max + 1;

    min_lt = l->min_lt;
    if (min != 0) {
        l->min = max;
        if (l->max_lt != 0) {
            l->min_lt = l->max_lt;
        }
    }
    if (max != 0) {
        l->max = min;
        if (min_lt != 0) {
            l->max_lt = min_lt;
        }
    }
}

void complement_seq_ft (SEQUENCE *s) {

    int i, num_entry;
    ft_entry *e;
    ft_range *r;

    if(s->feature_table) {
        num_entry = s->feature_table->num_entry;
        for (i = 0; i < num_entry; i++) {
            e = s->feature_table->entry[i];
             for (r = e->range; r; r = r->next) {
                 ft_location *fl = NULL;
                 if (r->complemented)
                     r->complemented = 0;
                 else
                     r->complemented = 1;
                 if (r->left) {
                     complement_location (r->left, s->length);
                     fl = r->left;
                 }
                 if (r->right) {
                     complement_location (r->right, s->length);
                     r->left = r->right;
                 }
                 if (fl) r->right = fl;

             }
        }
    }
}

int complement_fragment (SEQUENCE *c) {

    int len, s_len;
    char *hang;
    char *s;

    /* C|cgcatcgggttc|CG C */
    /* G GC|gtagcccaaggc|G */
    /*FIXME*/
    if (c->left_end != NULL && c->left_end->len != 0
        && c->right_end != NULL && c->right_end->len != 0) {
        if (c->left_end->len > 0 && c->right_end->len < 0) {
            s_len = c->length - c->left_end->len + ABS(c->right_end->len);
            len = c->left_end->len;
            /* complement sequence */
            if ((NULL == (s = (char * ) xmalloc (s_len + 1) ))) return -1;
            strcpy (s, &c->seq[len]);
            strcat (s, c->right_end->hang);
        }

        if (c->left_end->len < 0 && c->right_end->len > 0) {
            s_len = c->length + ABS(c->left_end->len) - c->right_end->len;
            len = c->length - c->right_end->len;
            if ((NULL == (s = (char * ) xmalloc (s_len + 1) ))) return -1;
            strcpy (s, c->left_end->hang);
            strncat (s, &c->seq[0], len);
        }
            complement_seq (s, s_len);
            xfree(c->seq);
            c->seq = s;
            c->length =s_len;
    } else {
        complement_seq (c->seq, c->length);
    }

    /* complement hang */
    if (c->left_end != NULL && c->left_end->len != 0) {
        len = c->left_end->len;
        hang = strdup (c->left_end->hang);
        complement_seq ( hang, ABS(len));
        if (c->right_end != NULL && c->right_end->len != 0) {
            c->left_end->len = -c->right_end->len;
            c->left_end->hang = strdup (c->right_end->hang);
                complement_seq (c->left_end->hang , ABS(c->left_end->len ));
            c->right_end->len = -len;
            c->right_end->hang = strdup (hang);
        } else {
            c->left_end = NULL;
            c->right_end->len = -len;
            c->right_end->hang = strdup (hang);
        }
    } else {
        if (c->right_end != NULL && c->right_end->len != 0) {
            len = c->right_end->len;
            hang = strdup (c->right_end->hang);
            complement_seq ( hang, ABS(len));
            c->right_end = NULL;
            c->left_end->len = -len;
            c->left_end->hang = strdup (hang);
        }
    }
    return 0;
}

SEQUENCE *made_complement_for_fragment (SEQUENCE *clipboard) {

    SEQUENCE *c;

    c = copy_sequence (clipboard);
    complement_fragment (c);
    complement_seq_ft (c);

    /*dis_fragment (clipboard);
    printf("complement=============================\n");
    dis_fragment (c);*/

    return c;
}

int editor_complement_shutdown (EDITOR_RECORD *er, EDIT *edit, int add_to_undo) {

    seq_reg_info info;
    int seq_num, seq_id;

    /*FIXME: in this case, edit->seq_id stores seq_num
      of the copied sequence for inserting complement fragment*/
    seq_num = edit->seq_id;
    seq_id = GetEdenId (seq_num);

    /*graphic_editor_shutdown & text_editor_shutdown, if has been opened*/
    info.job = GRAPHIC_EDITOR_QUIT;
    editor_notify(seq_num, (editor_reg_data *)&info);

    /*remove_sequence_from_registration*/
    remove_sequence_from_registration(seq_num);

    /*remove_sequence_from_sequences_list */
    remove_sequence_from_sequence_list (seq_id);
    return 0;
}

EDIT *create_insert_fragment_edit_complement (EDITOR_RECORD *er, int member_id, SEQUENCE *sc, int pos) {

    EDIT *edit;

    edit = create_edit();
    edit->seq_id    = member_id;
    edit->position  = pos;
    edit->operation = INSERT_FRAGMENT_COMPLEMENT;
    edit->sequence  = copy_sequence (sc);

    return edit;
}

int insert_fragment_in_sequence_c (int seq_id, SEQUENCE *f, int position) {

    EDIT *edit;
    EDITOR_RECORD *er;
    int ed_id, group_id, member_id;
    int err;

    ed_id = GetEdIdFromSeqId (seq_id);
    member_id = GetMemberIdFromSeqId (seq_id);
    group_id = GetGroupIdFromSeqId (seq_id);
    er = GetEditor (ed_id);

    /*check_joint (ed_id, group_id, member_id, f, position);*/

    edit = create_insert_fragment_edit_complement (er, member_id, f, position);
    err = insert_fragment_in_editor (er, edit, 1);

    if (err != 0) return -1;
    err = update_consensus (er->seq_group[group_id], 1);
    return 0;
}

int insert_update_sequence_complement (int ed_id, int member_id, int pos, SEQUENCE *cf) {

    int err, add_to_undo = 1;
    EDIT *edit;
    EDITOR_RECORD *er;
    int group_id;

    er = GetEditor (ed_id);
    pos -- ;
    edit = create_insert_fragment_edit (er, member_id, cf, pos, 0);
    err = insert_fragment_in_editor (er, edit, add_to_undo);
    if (err != 0) return -1;
    group_id = 0;
    err = update_consensus (er->seq_group[group_id], 1);

    return 0;
}

int insert_fragment_in_sequence_2 (Tcl_Interp *interp, int seq_id_c, SEQUENCE *cf, int position, int seq_id) {

    EDITOR_RECORD *er;
    EDIT *edit;
    SEQUENCE *s;
    int seq_num, ed_id, group_id, member_id;
    int err;

    /*create new editor for seq_id_c and should not inherit undo from it,
      it is because this is a new sequence*/
    ed_id = add_new_editor (interp, seq_id_c);
    seq_num = GetEdenNum(seq_id_c);
    group_id = GetGroupIdFromSeqId (seq_id_c);
    member_id = GetMemberIdFromSeqId (seq_id_c);
    er = GetEditor (ed_id);
    s = er->seq_group[group_id]->members[member_id]->data->sequence,
    position --;
    err = insert_fragment_in_sequence (s, cf, position, 0);
    if (err != 0) return -1;

    ed_id = GetEdIdFromSeqId (seq_id);
    er = GetEditor (ed_id);
    /* the seq_num of the copied sequence for insert complement fragment
       store in edit for UNDO */
    edit = create_insert_fragment_edit (er, seq_num, cf, position, 0);
    err = add_edit(er->edits, edit);
    if (err) return -1;

    return 0;
}

int insert_complement_fragment_at_cut_site (Tcl_Interp *interp,
                                            SEQUENCE *c,
                                            int seq_id_c,
                                            int cut_site,
                                            char *renz_start,
                                            char *renz_end,
                                            int seq_id) {

    char *seq_name, *renz_name;
    char seqid[10];

    /*insert complement fragment to the sequence in sequences list*/
    if (-1 == insert_fragment_in_sequence_2 (interp, seq_id_c, c, cut_site, seq_id)) return -1;

    seq_name = strdup (GetSequenceName (seq_id_c));
    sprintf (seqid, "%d", seq_id_c);

    if (strcmp (renz_start, renz_end)) {
        int len_start, len_end;
        len_start = strlen (renz_start);
        len_end = strlen (renz_end);
        if ( ( NULL == ( renz_name = (char * ) xmalloc ( len_start + len_end + 1) ))) return -1;
        strcpy (renz_name, renz_start);
        strcat (renz_name, renz_end);
    } else
        renz_name = renz_start;
    /*draw */
    if (TCL_OK != Tcl_VarEval(interp, "renzyme_draw_map ",
                              renz_name, " ", seqid, " ", seq_name, " ",  NULL)) {
        fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));
    }
    return 0;
}

int insert_fragment_in_two_direction (Tcl_Interp *interp,
                                      int seq_id,
                                      SEQUENCE *s,
                                      SEQUENCE *clipboard,
                                      int cut_site,
                                      char *renz_start,
                                      char *renz_end)
{

    SEQUENCE *c, *snew;
    int seq_id_c;

    c = made_complement_for_fragment (clipboard);
    snew = copy_sequence (s);
    seq_id_c = add_sequence_to_sequences (snew);

    /* The difference between this function and insert_fragment_at_cut_site
       is this function does not need to do check overhang but to open
       another renzyme map window with same paremeter */

    insert_complement_fragment_at_cut_site (interp, c, seq_id_c, cut_site, renz_start, renz_end, seq_id);
    insert_fragment_in_sequence_c (seq_id, clipboard, cut_site);
    return 0;
}

char *create_end_editor (Tcl_Interp *interp, int seq_id) {

    char seqid[10];
    char *op;

    sprintf (seqid, "%d", seq_id);

    if (TCL_OK != Tcl_VarEval(interp, "EndEditor ", seqid, " ", NULL)) {
        fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));
    }

    op = Tcl_GetStringResult(interp);
    /*printf("create_end_editor:op=%s\n", op);*/

    return op;
}

void end_editor_shutdown (int seq_id) {

    seq_reg_info info;
    int seq_num;

    info.job = END_EDITOR_QUIT;
    seq_num = GetEdenNum (seq_id);

    editor_notify(seq_num, (editor_reg_data *)&info);
}

/*0: success;
  1: clipboard is NULL;
  -1:failure.*/
int insert_fragment_at_cut_site (Tcl_Interp *interp, int seq_id, int cut_site,
                                 char *renz_start, char *renz_end) {

    SEQUENCE *s, *buffer;
    int seq_num, err;

    seq_num = GetEdenNum(seq_id);
    s = GetEdenSequence (seq_num);
    buffer = get_editor_buffer ();
    if (buffer == NULL) {
        bell();
        return 1;
    }

    /* chech overhang before inseration */
    err = overhang_check (s, buffer);

    if (err == 0) {
        /* made a paste at gaving position */
        /*printf ("matches in forward direction\n");*/
        if (-1 == insert_fragment_in_sequence_1 (seq_id, buffer, cut_site) ) return -1;
        return 0;
    }
    if (err == 1) {
        SEQUENCE *c;
        /* made a complement for clipboard */
        /*printf ("matches in complement direction\n");*/
        c = made_complement_for_fragment (buffer);
        if (-1 == insert_fragment_in_sequence_1 (seq_id, c, cut_site) ) return -1;
        return 0;
    }

    /* if matches happened at both direction, made copy of original sequence and
       insert fragment into both sequence with two direction */
    if (err == 2) {
        /*printf ("matches in two directions\n");*/
        if (!insert_fragment_in_two_direction (interp, seq_id, s, buffer, cut_site, renz_start, renz_end))
            return 0;
    }

    if (err == 3) {/*edit the ends of the sequence and fragment*/

        char *op;

        if (TCL_OK != Tcl_VarEval(interp, "EditHang", NULL)) {
            fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));
        }
        op = Tcl_GetStringResult(interp);
     
        if (!strcmp (op, "edit")) {
            char *oper;
            oper = create_end_editor (interp, seq_id);
            if (!strcmp (oper, "ok")) {

                int err, cut_pos;
            
                /* shut_down ends editor before do insert */
                end_editor_shutdown (seq_id);

                cut_pos = s->cut_site_1->pos;/* why cut_site_1? because after */
                                             /* trimmung and filling, */
                                             /* the cut_position has been changed */
                err = insert_fragment_at_cut_site (interp, seq_id, cut_pos, renz_start, renz_end);
                return err;
            }
            if (!strcmp (oper, "help")) {
                /*printf("press help button\n");*/
                return -1;
            }
            if (!strcmp (oper, "cancel")) {
                return -1;
            }
        }
        if (!strcmp (op, "cancel")) {
            return -1;
        }
        if (!strcmp (op, "help")) {
            return -1;
        }
    }
    return 0;
}

void set_cut_site (int seq_num, int left, int right, char *renz_left, char *renz_right) {

    SEQUENCE *s;
    SITE_HANG *sh1, *sh2;
    int oh;

    s = GetEdenSequence (seq_num);
    oh = -get_left_over_hang (renz_left);
    sh1 = init_site_hang ();
    sh1->pos = left;
    sh1->len = oh;
    if (oh < 0) {
        int h = ABS(oh);
        sh1->hang = get_hang (s, left, h);
    }
    if (oh > 0) {
        sh1->hang = get_hang (s, left-oh, oh);
    }
    s->cut_site_1 = sh1;

    oh = -get_right_over_hang ( renz_right );
    sh2 = init_site_hang ();
    sh2->pos = right;
    sh2->len = oh;
    if (oh < 0) {
        int h = ABS(oh);
        sh2->hang = get_hang (s, right-h, h);
    }
    if (oh > 0) {
        sh2->hang = get_hang (s, right, oh);
    }
    s->cut_site_2 = sh2;
    
    /*dis_sequence(s);*/
}

int sequences_notify_graphic (int seq_num, int pos) {

    seq_reg_changed sc;
    /*seq_reg_cursor_notify cn;*/
    SELECTION *sel;
    /*int ed_id, group_id, member_id, seq_id;*/

    sc.job = SEQ_CHANGED;
    sel = init_selection ();
    sel->sel_first = 0;
    sel->sel_last = 0;
    sc.selection = sel;
    editor_notify(seq_num, (editor_reg_data *)&sc);

#if 0
    pos ++;
    seq_id = GetEdenId (seq_num);
    ed_id = GetEdIdFromSeqId (seq_id);
    group_id = GetGroupIdFromSeqId (seq_id);
    member_id = GetMemberIdFromSeqId (seq_id);

    /* cn.cursor = get_editor_cursor (ed_id, member_id, NULL); */
    if (cn.cursor == NULL) {
          printf ("graphic_editor cursor move error\n");
    }
    cn.cursor->abspos = pos;
    cn.cursor->job = CURSOR_MOVE;
    cn.cursor->sent_by = -1;
    cn.job = SEQ_CURSOR_NOTIFY;
    editor_notify(seq_num, (editor_reg_data *)&cn);
#endif

    return 0;
}

int delete_fragment_from_sequence (int seq_id,
                                   int start,
                                   int end,
                                   char *renz_start,
                                   char *renz_end,
                                   SEQUENCE *fd) {

    int seq_num, ed_id, group_id, member_id;
    int err;

    seq_num = GetEdenNum(seq_id);
    ed_id = GetEdIdFromSeqId (seq_id);
    group_id = GetGroupIdFromSeqId (seq_id);
    member_id = GetMemberIdFromSeqId (seq_id);

    /*fd = create_fragment_from_sequence (seq_num, start, end, renz_start, renz_end);*/
    err = delete_update_sequence (ed_id, group_id, member_id, start, fd);

    /* set overhangs of the sequence for gaving sq_id cot_position and renzyme */
    set_cut_site (seq_num, start, end, renz_start, renz_end);

    if (err != 0) return -1;
    return 0;
}

int replace_fragment_from_sequence (Tcl_Interp *interp, int seq_id, int start, int end, char *renz_start, char *renz_end) {

    SEQUENCE *s, *fd = NULL;
    SEQUENCE *buffer;
    FEATURE_TABLE *ft;
    int seq_num, err;

    seq_num = GetEdenNum(seq_id);
    s = GetEdenSequence (seq_num);
    
    /* set overhangs of the sequence for gaving sq_id cot_position and renzyme */
    set_cut_site (seq_num, start, end, renz_start, renz_end);
    /* chech overhang before replacing */
    buffer = get_editor_buffer ();
    if (buffer == NULL) {
        bell();
        return 1; /* buffer NULL */
    } 
    err = overhang_check (s, buffer);

    if (err == 3) {/*edit the ends of the sequence and fragment*/

        char *op;

        if (TCL_OK != Tcl_VarEval(interp, "EditHang", NULL)) {
            fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));
        }
        op = Tcl_GetStringResult(interp);
        if (!strcmp (op, "edit")) {
            char *oper;
            oper = create_end_editor (interp, seq_id);
            if (!strcmp (oper, "ok")) {

                int cut_pos;
            
                /* shut_down ends editor before insertation */
                end_editor_shutdown (seq_id);
                cut_pos = s->cut_site_1->pos;/* why cut_site_1? because after */
                                             /* trimmung and filling, */
                                             /* the cut_position has been changed */
            }
            if (!strcmp (oper, "help")) {
                return -1;
            }
            if (!strcmp (oper, "cancel")) {
                return -1;
            }
        }
        if (!strcmp (op, "cancel")) {
            return -1;
        }
        if (!strcmp (op, "help")) {
            return -1;
        }
    }

    /* made a copy of replacing fragment, this information can be used in UNDO */
    fd = create_fragment_from_sequence (seq_num, start, end, renz_start, renz_end);
    /* delete fragment from sequence include feature_table */
    
    /* the returned ft for UNDO, include deleted entries, not to be used here.*/
    ft = delete_fragment_in_sequence (s, fd, start);
    err = insert_fragment_at_cut_site (interp, seq_id, start, renz_start, renz_end);
    return err;

}

int SequencesRedisplayGraphic (ClientData clientData,
                               Tcl_Interp *interp,
                               int argc,
                               char *argv[])
{
    char *operation;
    int seq_id, cut_pos, seq_num;
    int err;

    if (argc < 3) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " operation\"", "seq_id\"", (char*)NULL);
        return TCL_ERROR;
    }
    operation = argv[1];
    seq_id = atoi(argv[2]);

    if (!strcmp (operation, "paste")) {
        char *renz_name;
        cut_pos = atoi(argv[3]);

        if (argc != 5) {
            Tcl_AppendResult(interp, "wrong # args: should be \"",
                             argv[0], "operation\"", "seq_id\"", "pos\"",
                             "re_name\"",(char*)NULL);
            return -1;
        }
        renz_name = argv[4];
        seq_num = GetEdenNum(seq_id);
        /*printf ("seq_id =%d    seq_num=%d\n",  seq_id, seq_num);*/
        /* set overhangs of the sequence for gaving seq_id cut_position and renzymes */
        set_cut_site (seq_num, cut_pos, cut_pos, renz_name, renz_name);
        err = insert_fragment_at_cut_site (interp, seq_id, cut_pos, renz_name, renz_name);
    }

    if (!strcmp (operation, "cut")) {

        SEQUENCE *fd;
        int start, end;
        char *renz_start, *renz_end;

        if (argc != 7) {
            Tcl_AppendResult(interp, "wrong # args: should be \"",
                             argv[0], "operation\"","start_position\"", "end_position\"",
                             "start_enzyme_name\"","end_enzyme_name\"", (char*)NULL);
            return TCL_ERROR;
        }
        start = atoi(argv[3]);
        end = atoi(argv[4]);
        renz_start = argv[5];
        renz_end = argv[6];
        if (start > end) {
            start = atoi(argv[4]);
            end = atoi(argv[3]);
            renz_start = argv[6];
            renz_end = argv[5];
        }
        /* err==0, to notify */
        cut_pos = start;
        seq_num = GetEdenNum(seq_id);
        fd = create_fragment_from_sequence (seq_num, start, end, renz_start, renz_end);
        if (fd != NULL) set_editor_buffer (fd);
        err = delete_fragment_from_sequence (seq_id, start, end, renz_start, renz_end, fd);
    }

    if (!strcmp (operation, "replace")) {

        int start, end;
        char *renz_start, *renz_end;

        if (argc != 7) {
            Tcl_AppendResult(interp, "wrong # args: should be \"",
                             argv[0], "operation\"","start_position\"", "end_position\"",
                             "start_enzyme_name\"","end_enzyme_name\"", (char*)NULL);
            return -1;
        }
        start = atoi(argv[3]);
        end = atoi(argv[4]);
        renz_start = argv[5];
        renz_end = argv[6];

        if (start > end) {
            start = atoi(argv[4]);
            end = atoi(argv[3]);
            renz_start = argv[6];
            renz_end = argv[5];
        }
        seq_num = GetEdenNum(seq_id);
        err = replace_fragment_from_sequence (interp, seq_id, start, end, renz_start, renz_end);
        cut_pos = start;
    }

    if (!strcmp (operation, "copy")) {

        SEQUENCE *c,*c_new;
        int start, end;
        char *renz_start, *renz_end;
        int seq_num;

        if (argc != 7) {
            Tcl_AppendResult(interp, "wrong # args: should be \"",
                             argv[0],"operation\"", "seq_id\"", "start_position\"",
                             "end_position\"", "start_enzyme_name\"",
                             "end_enzyme_name\"", (char*)NULL);
            return TCL_ERROR;
        }
        start = atoi(argv[3]);
        end = atoi(argv[4]);
        renz_start = argv[5];
        renz_end = argv[6];
        if (start > end) {
            start = atoi(argv[4]);
            end = atoi(argv[3]);
            renz_start = argv[6];
            renz_end = argv[5];
        }
        seq_num = GetEdenNum(seq_id);
        c_new = create_fragment_from_sequence (seq_num, start, end, renz_start, renz_end);
        if (c_new != NULL) {
            c = get_editor_buffer();
            if (c != NULL) free_sequence (c);
            set_editor_buffer (c_new);
        }
        err = 1; /* made copy does not need to notify */
    }
    if (!strcmp (operation, "select")) {

        seq_reg_selected info;
        SELECTION *sel;
        int seq_id, seq_num, first, last;

        if (argc != 7) {
            Tcl_AppendResult(interp, "wrong # args: should be \"",
                             argv[0],"operation\"", "seq_id\"", "start_position\"",
                             "end_position\"", (char*)NULL);
            return TCL_ERROR;
        }
        seq_id = atoi(argv[2]);
        first = atoi(argv[3]);
        last = atoi(argv[4]);

        if (seq_id != -1) seq_num = GetEdenNum (seq_id);	
        else {
            /*printf (" sequence selection error\n");*/
            exit (0); /*FIXME*/
        }
        sel = init_selection ();
        if ( first <= last ) {
            sel->sel_first = first;
            sel->sel_last = last;
            sel->name_first = strdup(argv[5]);
            sel->name_last = strdup(argv[6]);
        } else {
            sel->sel_first = last;
            sel->sel_last = first;
            sel->name_first = strdup(argv[6]);
            sel->name_last = strdup(argv[5]);
        }
        set_editor_selection (sel);
        /* selection notify */
        info.job = SEQ_SELECTED;
        info.selection = sel;
        info.frame_name = ""; 
        editor_notify (seq_num, (editor_reg_data *)&info);
        err = 1;
    }
    if (err == 0)
        err = sequences_notify_graphic (seq_num, cut_pos);
    /*if (err == -1)
      goto err;*/
    return TCL_OK;
}

int REnzyme_search_Init(Tcl_Interp *interp) {

    Tcl_CreateCommand(interp, "renzyme_search_preview", RenzymeSearchPreview,
                      (ClientData) NULL,
                      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "renzyme_search_select", RenzSearchSelect,
                      (ClientData) NULL,
                      (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "sequences_redisplay_graphic", SequencesRedisplayGraphic,
                      (ClientData) NULL,
                      (Tcl_CmdDeleteProc *) NULL);
    return TCL_OK;
}


