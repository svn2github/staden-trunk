/*
 * renzyme.c to parse "renzyme_bairoch" and get information of each entry.
 *
 */

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
#include "getfile.h"
#include "renzyme_box.h"

#define LEN 1024
RENZYMES *renzymes = NULL;

/* initialise a RENZYME structure */
RENZYME *init_renzyme (void) {

    RENZYME  *r;
    if ( NULL == ( r = (RENZYME* ) xmalloc (sizeof (RENZYME))))
        return NULL;
    r->name= NULL;
    r->ID = 0;
    r->rec_seq = NULL;
    r->rec_seq_text = NULL;
    r->cut_pos_1 = 0;
    r->cut_pos_2 = 0;
    r->av_frag_size = 0;
    r->prototype = NULL;
    r->supplier_codes = NULL;
    r->stock_level = 0;
    return r;
}

/* free a RENZYME structure */
void free_renzyme (RENZYME *renzyme) {

  if (renzyme) {
    if (renzyme->name) xfree(renzyme->name);
    if (renzyme->prototype) xfree(renzyme->prototype);
    if (renzyme->supplier_codes) xfree(renzyme->supplier_codes);
    xfree (renzyme);
  }
}

/* initialise a RENZYMES structure */
RENZYMES *init_renzymes (void) {

    RENZYMES *rs;
    RENZYME  **r;
    if ( NULL == (rs = (RENZYMES *) xmalloc (sizeof (RENZYMES))))
        return  NULL;
    if ( NULL == (r = (RENZYME **) xmalloc (sizeof (RENZYME*))))
        return  NULL;

    rs->renzyme = r;
    rs->used = 0;
    rs->capacity = 0;
    rs->max_rec_seq = 0;
    return  rs;
}

/* free a RENZYMES structure */
void free_renzymes ( RENZYMES *rs ) {

    int i;
    if (rs) {
        for ( i = 0; i <= rs->used; i++ ) {
            free_renzyme (rs->renzyme[i]);
        }
    }
    /*xfree (rs);*/ /* FIXME */
}


/* realloc the renzyme array in the renzymes structure */
RENZYMES *realloc_renzymes ( RENZYMES *rs, int num_ren) {

    RENZYME  **r;
    r = rs->renzyme;
    if ((NULL == ( r = (RENZYME **)xrealloc (r, sizeof(RENZYME *)*(num_ren+1)))))
        return NULL;
    rs->renzyme = r;
    rs->used = num_ren;
    return rs;
}

int parse_rs(char *rs, float *efs, int *f, int *r, char **rss) {

  int i, rs_len = 0;
  float efss ;
  char *a, *b;
  int j;
  char *rr;

  rs_len = strlen(rs);

  if (NULL == (a = (char *)malloc((rs_len + 1) * sizeof(char))))
        goto err;
  if (NULL == (b = (char *)malloc((rs_len + 1) * sizeof(char))))
        goto err;
  if (NULL == (rr = (char *)malloc((rs_len + 1) * sizeof(char))))
        goto err;
  /*a[rs_len] = '\0';
    b[rs_len] = '\0';*/
  for (i = 0; i < rs_len; i++) {
    if (rs[i] == '?') {
      free(a);
      free(b);
      free(rr);
      return 0;
    }
  }
  j = 0;
  i = 0;
  while (rs[i] != ';') {
      if (isdigit (rs[i])) {
          a[j] = rs[i];
          j++;
      }
      i++;
  }
  a[j] = 0;
  j = 0;
  b[0] = '0';
  if ( (i+1) != rs_len) {
      i++;
      while (rs[i] != ';') {
          if (isdigit (rs[i])) {
              b[j] = rs[i-1];
              j++;
          }
          i++;
      }
      b[j] = rs[i-1];
  }
  b[j+1] = 0;
  *f = atoi(a);
  *r = atoi(b);
  i = 0;
  efss = 1;
  while( rs[i] != ',' ) {
      rr[i] = rs[i];
      if ( rs[i] == 'A' || rs[i] == 'C' || rs[i] == 'G' || rs[i] == 'T') efss = efss*4;
      else if (rs[i] == 'R' || rs[i] == 'Y' || rs[i] == 'W' || rs[i] == 'S'||
               rs[i] == 'S' || rs[i] == 'M' || rs[i] == 'k') efss = efss*2;
      else if (rs[i] == 'H' || rs[i] == 'B' || rs[i] == 'V' || rs[i] == 'D') efss = efss*4/3;
      else efss = efss;
      i++;
  }
  *efs = efss;
  rr[i]=0;
  xfree (a);
  xfree (b);
  if(efss > 200000000) {
    free(rr);
    *rss = NULL;
    return 0;
  }
  *rss = rr;
  return 1;
 err:
  if (a) xfree (a);
  if (b) xfree (b);
  if (rr) xfree (rr);
  printf("ERROR parsing '%s'\n", rs);
  return -1;
}

RENZYMES *get_enzyme(char *filename) {
    RENZYME  *r;
    RENZYMES *rs;
    FILE *pr;
    char line[LEN];
    char tmpp[1024];
    char tmp[1024];
    char exp_filename[FILENAME_MAX];
    int flag = 0;
    float efs;
    int current_num_entry = 0, num_entry = 0;
    char *id, *pt, *cr, *rss, *rst;
    int forward, reverse;

    expandpath(filename, exp_filename);

    if (NULL == (pr = fopen(exp_filename, "r"))){
        verror(ERR_WARN, "read r_enzyme file",
               "Unable to open r_enzyme file %s", exp_filename);
        return TCL_OK;
    }

    rs = init_renzymes();
    if(!rs) goto err;

    while (!feof(pr)) {
        fgets(line, LEN, pr);
        line[LEN-1]=0;
        strcpy(tmpp, line);
#if 0
        llen = strlen(line);
        while(llen >= 3 && line[llen - 2] == ',') {
            fgets(line, LEN, pr);
            tmpp[strlen(tmpp) - 1] = 0;
            strcat(tmpp, &line[5]);
        }
#endif
        if (!strncmp(&tmpp[0], "ID", 2)){

            strcpy(tmp, &tmpp[5]);
            tmp[strlen(tmp)- 1] = 0;
            if (num_entry > 0) {
                /*printf("name=%s   %s\n", rs->renzyme[num_entry-1]->name, tmp);*/
                if (!strcmp (rs->renzyme[num_entry-1]->name, tmp)){
                    id = NULL;
                } else id = strdup (tmp);
            } else id = strdup (tmp);
        }
        if (!strncmp(&tmpp[0], "RS", 2)){
            forward = 0;
            reverse = 0;
            strcpy(tmp, &tmpp[5]);
            tmp[strlen(tmp)- 1] = 0;
            rst = strdup (tmp);
            flag = parse_rs(tmp, &efs, &forward, &reverse, &rss);
            if( flag == 1) {
                current_num_entry++;
            }
        }
        if (!strncmp(&tmpp[0], "PT", 2)) {
            strcpy(tmp, &tmpp[5]);
            tmp[strlen(tmp)- 1] = 0;
            pt = strdup (tmp);
        }
        if (!strncmp(&tmpp[0], "CR", 2)){
            strcpy(tmp, &tmpp[5]);
            tmp[strlen(tmp)- 1] = 0;
            cr = strdup (tmp);
            /*printf("supplier_codes=%s\n",cr);*/
        }
        if((!strncmp(&tmpp[0], "//", 2)) && flag == 1 && id != NULL){

            r = init_renzyme();
            if (!r) goto err;
            r->name = id;
            r->rec_seq = rss;
            r->rec_seq_text = rst;
            r->prototype = pt;
            r->supplier_codes = cr;
            r->av_frag_size = efs;
            r->cut_pos_1 = forward;
            r->cut_pos_2 = reverse;
            efs = 0;
        if (num_entry != 0) {
            rs = realloc_renzymes (rs, num_entry);
        }
        rs->renzyme[num_entry] = r;
        num_entry++;
        rs->used = num_entry;
        }
    }
    fclose(pr);
    return rs;
 err:
    if (rs) free_renzymes (rs);
    if (r)  free_renzyme (r);
    fclose(pr);
    return NULL;
}

/*int GetRenzInfo(ClientData clientData, Tcl_Interp *interp, int argc, char **argv) {

    int num_entry;
    int i;
    char buf[1024];

    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], " filename\"", (char*)NULL);
        return TCL_ERROR;
    }

    if (!renzymes) {
        free_renzymes (renzymes);
    }

    renzymes = get_enzyme(argv[1]);
    printf("num_entry=%d\n", renzymes->used);

    if (!renzymes)
        return TCL_OK;

    num_entry = renzymes->used;
    Tcl_ResetResult(interp);
    for (i = 0; i < num_entry; i++) {
        sprintf(buf, "%s {%s} %s %s %.0f",renzymes->renzyme[i]->name,
                renzymes->renzyme[i]->rec_seq_text,
                renzymes->renzyme[i]->prototype,
                renzymes->renzyme[i]->supplier_codes,
                renzymes->renzyme[i]->av_frag_size);
        Tcl_AppendElement(interp, buf);
    }
    return TCL_OK;
}
*/
int GetRenzInfo(ClientData clientData, Tcl_Interp *interp, int argc, char **argv) {

    int num_entry;
    int i;

    if (argc != 2) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
            argv[0], " filename\"", (char*)NULL);
        return TCL_ERROR;
    }

    if (!renzymes) {
        free_renzymes (renzymes);
    }

    renzymes = get_enzyme(argv[1]);
    /* printf("num_entry=%d\n", renzymes->used); */
    if (!renzymes)
        return TCL_OK;

    num_entry = renzymes->used;
    Tcl_ResetResult(interp);
    for (i = 0; i < num_entry; i++) {
        Tcl_DString dstr;
        Tcl_DStringInit(&dstr);
        vTcl_DStringAppendElement(&dstr, "%s", renzymes->renzyme[i]->name);
        vTcl_DStringAppendElement(&dstr, "%s", renzymes->renzyme[i]->rec_seq_text);
        vTcl_DStringAppendElement(&dstr, "%s", renzymes->renzyme[i]->prototype);
        vTcl_DStringAppendElement(&dstr, "%s", renzymes->renzyme[i]->supplier_codes);
        vTcl_DStringAppendElement(&dstr, "%.f", renzymes->renzyme[i]->av_frag_size);
        Tcl_AppendElement(interp, Tcl_DStringValue(&dstr));

        Tcl_DStringFree(&dstr);
    }
    return TCL_OK;
}


int save_renzyme(FILE *pw, Tcl_Interp *interp, char **sel, int num_sel) {

    int i, j;
    char **item = NULL;
    int num_item;
    char item_name[5][4] = {"ID", "RS", "PT", "CR", "EFS"};

    for (i = 0; i < num_sel; i++) {
         if (Tcl_SplitList(interp, sel[i], &num_item, &item) != TCL_OK)
             return TCL_ERROR;
         for (j = 0; j < num_item; j++) {
             fprintf(pw, "%2s   %s\n", item_name[j], item[j]);
         }
        fprintf(pw, "//\n\n");
    }

    fclose(pw);
    ckfree((char*)item);
    return 0;
}

int SaveRenzInfo(ClientData clientData, Tcl_Interp *interp, int argc, char **argv) {

    FILE *pw;
    char **sel = NULL;
    int num_sel;
    if (argc != 3) {
        Tcl_AppendResult(interp, "wrong # args: should be \"",
                         argv[0], " filename selected_items\"", (char*)NULL);
        return TCL_ERROR;
    }
    if (NULL == (pw = fopen(argv[1], "w"))) {
        verror(ERR_WARN, "save personal r_enzyme file", "Unable to open file %s", argv[1]);
        return TCL_OK;
    }

    /* create selecing Renzyme name array */
    if (Tcl_SplitList(interp, argv[2], &num_sel, &sel) != TCL_OK)
     return TCL_ERROR;

    save_renzyme(pw, interp, sel, num_sel);
    return TCL_OK;
}

int REnzyme_box_Init(Tcl_Interp *interp) {

    if (Itcl_RegisterC(interp, "get_renz", GetRenzInfo, NULL, NULL) != TCL_OK) {
      return TCL_ERROR;
    }
    if (Itcl_RegisterC(interp, "save_renz", SaveRenzInfo, NULL, NULL) != TCL_OK) {
        return TCL_ERROR;
    }
    return TCL_OK;
}


