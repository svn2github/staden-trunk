#ifndef _SPBIOLIMS_H
#define _SPBIOLIMS_H

/**********************************************************************
 *
 * spBiolims.h
 *************
 * This file contains the Staden Package interface to the Biolims IO library
 * No Biolims types or defines are required
 * This is intended to be used by Staden Package 'C' source code
 * it should not require ANY C++
 *
 **********************************************************************/

#define BIOLIMS_TAG "Biolims="
#define IS_BIOLIMS_PATH(x)  !strncmp(x,BIOLIMS_TAG,strlen(BIOLIMS_TAG))

/*
 * Strings for encoding Exp tags in Read structure info
 ******************************************************
 This is so we done lose info when generating an EXP via a Read Structure
 */
#define EXP_TAGLEN 5
#define EXP_CHEM "CHEM=" /* CH */
#define EXP_FEAT "FEAT="  
#define EXP_PRMR "PRMR=" /* PR */
#define EXP_VECT "VECT=" /* SV */
#define EXP_CLON "CLON=" /* CN */
#define EXP_CLOV "CLOV=" /* CV */

/*
 * Biolims String constants for conversion from EXP tags (defined in Exp.cpp)
 *******************************************************
 */

/* Strings for Biolims Chemistry types */
extern const char *CH_types[]; /* 0/1 in bit 0 of CH*/

/* Strings for Biolims Primer types */
extern const char *PR_types[];
extern const int  nPR_types;

/* Staden Package FeatureKey names
 * Staden Package annotations are stored in Biolims with the following tags
 */
#define STADEN_FKEY_LEN 11
extern const char *featCLON; /*cloning vector position */
extern const char *featQUAL; /*quality sequence position */
extern const char *featVECT; /*sequencing vector position */
extern const char *featCONS; /*consensus tag  */
extern const char *featGELR; /*gel reading tag*/
extern const char *featVECI; /*vector insertion length*/
extern const char *featTEMP; /*template name*/
extern const char *featSTRD; /*strands*/ 

/* Biolims Annotation tag types
 * Biolims features are stored as Staden Package annotations
 * with the tag below. The feature key is stored in the first line of the comment
 */
extern const char *tagBIOL;
extern const char *featKey;


/*
 * Biolims interface Declarations
 *********************************
 */

typedef void* spBiolims;    /* pointer to a Biolims connection */

#ifdef __cplusplus
extern "C" {
#endif

  /*
   * C interface to the Biolims Connection Class
   * an interface layer is provided, so that no knowledge of PE Biolims classes
   * is required by the calling source
   */
  spBiolims spBiolimsCreate(int bBrowsing);          /* create a Biolims session
						      A connection to the db is not
						      opened until the first time
						      it is required.
						      A connection to the db
						      uses the biolims config file
						      specified by BIOLIMS_CONF environment.
						      This is the same file used by Sample2DB.
						      The current collection is set to / */

  void spBiolimsDestroy(spBiolims ptr);            /* destroy a Biolims session */

  char** spCollections(spBiolims ptr);             /* return a list of collections
						      contained by current collection*/

  char** spSamples(spBiolims ptr);                 /* return a list of samples contained
						      by current collection */

  char** spAssemblies(spBiolims ptr);              /*  return a list of assemblies contained
						       by current collection */
  int spSetCollection(spBiolims ptr, char* name);  /* set current collection to be name
                                                      use NULL to specify parent of current
						      use /coll1/coll2/coll3 to specify
						      an absolute collection path */
#ifdef  _Read_h_
  Read* spReadBiolimsReading(char *name);          /* fill a staden Read structure
						      from a biolims Lane */
#endif
#ifdef _EXPFILEIO_H_
  Exp_info *spBiolims2exp(char *fname);           /* fill an Exp_info structure from a biolims lane */
#endif
#ifdef _IO_H_
#ifdef _TCL
  int spBiolimsExport(GapIO *iof, char *name,Tcl_Interp* interp); /* export from an opened GAP DB to Biolims */
#endif
  int spBiolimsImport(GapIO *iot, char *name); /* import into opened GAP DB from Biolims */ 
#endif

  int spIsLane(spBiolims ptr, char *name);         /* return true if lane name exists in
						      current collection */
  int spIsPathToLane(spBiolims ptr, char *path, char *name); /* return true if lane name
								exists in specified collection */
  
  int spBiolimsUpdateFromExp(char *name);          /* update a biolims lane from
						      experiment file specified by name */

  int biolims_found(char* retpath, char* path, char* name); /* return the full biolims path of trace name */


#ifdef __cplusplus
}
#endif

#endif










