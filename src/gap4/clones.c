/*
 * File:
 * Version:
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description:
 *
 * Created:
 * Updated:
 *
 */


#include <stdlib.h>

#include "gap-defaults.h"
#include "IO.h"
#include "seqInfo.h"
#include "clones.h"
#include "io_utils.h"
#include "notes.h"
#include "tagUtils.h"




/*************************************************************
 * Find vector, clone, subclone
 *************************************************************/

static int find_vector(GapIO *io, char *V)
{
    char buf[128];		/* this should be more than enough */
    int i;
    GVectors t;

    for(i=0;i<io->db.Nvectors;i++) {
	/* read vector record */
	GT_Read(io,arr(GCardinal,io->vectors,i),&t,sizeof(t),GT_Vectors);
	/* read vector name */
	TextRead(io, t.name, buf, sizeof(buf));
	if (strcmp(buf,V)==0) return i+1;
    }

    return 0;
}

static int find_clone(GapIO *io, char *CN)
{
    char buf[128];		/* this should be more than enough */
    int i;
    GClones c;
    
    for(i=0;i<io->db.Nclones;i++) {
	/* read clone record */
	GT_Read(io,arr(GCardinal,io->clones,i),&c,sizeof(c),GT_Clones);
	/* read clone name */
	TextRead(io, c.name, buf, sizeof(buf));
	if (strcmp(buf,CN)==0) return i+1;
    }

    return 0;
}


/*************************************************************
 * Add vector, clone, subclone
 *************************************************************/



int add_vector(GapIO *io, char *V, int level)
{
    int err;
    int vector;
    int freerec;
    GVectors v;

    /* allocate vector name */
    v.name = allocate(io,GT_Text);
    err = TextWrite(io,v.name,V,strlen(V));
    v.level = level;

    /* find new vector number */
    vector = ++io->db.Nvectors;
    ArrayRef(io->vectors,vector-1);

    /* add vector to list */
    arr(GCardinal,io->vectors,vector-1) = freerec = allocate(io,GT_Vectors);
    err = GT_Write(io,freerec,&v,sizeof(v),GT_Vectors);

    /* write array */
    err = ArrayDelay(io, io->db.vectors, io->db.Nvectors, io->vectors);

    /* write db */
    DBDelayWrite(io);

    return vector;
}





int add_clone(GapIO *io, char *CN, char *CV)
{
    int err;
    int clone;
    int freerec;
    GClones c;

    /* find vector record */
    c.vector = find_vector(io,CV);
    if (!c.vector) c.vector = add_vector(io,CV,1);

    /* allocate clone name */
    c.name = allocate(io,GT_Text);
    err = TextWrite(io,c.name,CN,strlen(CN));

    /* find new clone number */
    clone = ++io->db.Nclones;
    ArrayRef(io->clones,clone-1);

    /* add clone to list */
    arr(GCardinal,io->clones,clone-1) = freerec = allocate(io,GT_Clones);
    err = GT_Write(io,freerec,&c,sizeof(c),GT_Clones);

    /* write array */
    err = ArrayDelay(io, io->db.clones, io->db.Nclones, io->clones);

    /* write db */
    DBDelayWrite(io);

    return clone;
}






int add_template(GapIO *io, char *TN, char *SV, char *ST, char *SI, int clone)
{
    int err;
    int template;
    int freerec;
    GTemplates t;

    /* find vector record */
    t.vector = find_vector(io,SV);
    if (!t.vector) t.vector = add_vector(io,SV,2);

    /* allocate template name */
    t.name = allocate(io,GT_Text);
    err = TextWrite(io,t.name,TN,strlen(TN));

    /* set strand and clone */
    t.strands = atoi(ST);
    if (t.strands < 1 || t.strands > 2)
	t.strands = 1;
    t.clone = clone;

    /* set insert length */
    t.insert_length_max = t.insert_length_max = 0;
    sscanf(SI,"%d..%d",&t.insert_length_min,&t.insert_length_max);
    if (t.insert_length_max < t.insert_length_min)
	t.insert_length_max = t.insert_length_min;

    /* find new template number */
    template = ++io->db.Ntemplates;
    ArrayRef(io->templates,template-1);

    /* add template to list */
    arr(GCardinal,io->templates,template-1) = freerec = allocate(io,GT_Templates);
    err = GT_Write(io,freerec,&t,sizeof(t),GT_Templates);

    /* write array */
    err = ArrayDelay(io, io->db.templates, io->db.Ntemplates, io->templates);

    /* write db */
    DBDelayWrite(io);

    cache_template_name(io, template, TN);

    return template;
}





/*************************************************************
 *
 *************************************************************/

static int clonestuff(GapIO *io, SeqInfo *si)
/*
 * This routine add template/clone information to the database.
 * Returns template number.
 */
{
    char *CN;			/* clone name */
    char *TN;			/* template name */
    char *CV;			/* cloning vector */
    char *SV;			/* sequencing vector*/
    char *ST;			/* strands */
    char *SI;			/* insert size - a range */

    int template;
    int clone;

    /* find clone name */
    if (exp_Nentries(si->e,EFLT_CN))
	CN = exp_get_entry(si->e,EFLT_CN);
    else
	CN = UNKNOWN;

    /* find template name */
    if (exp_Nentries(si->e,EFLT_TN))
	TN = exp_get_entry(si->e,EFLT_TN);
    else if (exp_Nentries(si->e,EFLT_EN))
	TN = exp_get_entry(si->e,EFLT_EN);
    else if (exp_Nentries(si->e,EFLT_ID))
	TN = exp_get_entry(si->e,EFLT_ID);
    else
	TN = UNKNOWN;

    /* find cloning vector */
    if (exp_Nentries(si->e,EFLT_CV))
	CV = exp_get_entry(si->e,EFLT_CV);
    else
	CV = UNKNOWN;

    /* find sequencing vector */
    if (exp_Nentries(si->e,EFLT_SV))
	SV = exp_get_entry(si->e,EFLT_SV);
    else
	SV = UNKNOWN;

    /* find strands*/
    if (exp_Nentries(si->e,EFLT_ST))
	ST = exp_get_entry(si->e,EFLT_ST);
    else
	ST = DEFAULT_ST;

    /* find strands*/
    if (exp_Nentries(si->e,EFLT_SI))
	SI = exp_get_entry(si->e,EFLT_SI);
    else
	SI = DEFAULT_SI;



    /* find clone record */
    clone = find_clone(io,CN); 
    if (!clone) clone = add_clone(io,CN,CV);

    /* find template record */
    template = template_name_to_number(io,TN);

    /*
     * If this template isn't in the database - create it.
     * Otherwise, if we're dealing with a double stranded stranded
     * template we wish to 'upgrade' the information incase it's now been
     * double stranded from an, originally, single stranded template.
     *
     * FIXME: Is this always true?
     */
    if (!template)
	template = add_template(io,TN,SV,ST,SI,clone);
    else if (strcmp(ST, "2") == 0) {
	GTemplates t;
	GT_Read(io, arr(GCardinal, io->templates, template-1),
		&t, sizeof(t), GT_Templates);
	if (t.strands < 2) {
	    t.strands = 2;
	    GT_Write(io, arr(GCardinal, io->templates, template-1),
		     &t, sizeof(t), GT_Templates);
	}
    }

    return template;
}


#if 0
static int commentstuff(GapIO *io, SeqInfo *si)
/*
 *
 */
{
    return 0;
}
#endif


int add_seq_details(GapIO *io, int N, SeqInfo *si)
{
    int err;
    GReadings r;
    char *DR;
    char *PR;

    if (N>Nreadings(io)) err = io_init_reading(io,N);

    /* read record */
    gel_read(io, N, r);

    /* update cloning information */
    r.template = clonestuff(io,si);

    /* direction */
    if (exp_Nentries(si->e,EFLT_DR))
	DR = exp_get_entry(si->e,EFLT_DR);
    else
	DR = DEFAULT_DR;

    r.strand = (*DR!='+');

    /* primer */
    if (exp_Nentries(si->e,EFLT_PR))
	PR = exp_get_entry(si->e,EFLT_PR);
    else
	PR = DEFAULT_PR;

    r.primer = atoi(PR);
    r.strand = STRAND(r);

    /* Fix up old style primer and strand information */
    if (!exp_Nentries(si->e, EFLT_PR))
	r.primer = PRIMER_TYPE(r);

    /* chemistry */
    if (exp_Nentries(si->e, EFLT_CH)) {
	exp_get_int(si->e, EFLT_CH, &r.chemistry);
    } else {
	r.chemistry = 0;
    }
    
    /* write record back */
    gel_write(io, N, r);

    return 0;
}
