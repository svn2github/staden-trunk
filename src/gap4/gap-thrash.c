/*
 * File: gap-thrash.c
 *
 * Author: Simon Dear
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description: routines to thrash a gap database (client based)
 *
 * Created: 29-Sep-1992
 * Updated:	   
 *
 */

#include <stdio.h> /* IMPORT: NULL */

#include "g-defs.h" /* IMPORT: G_LOCK_RW G_LOCK_RO */
#include "g-misc.h" /* IMPORT: G_Number */

#include "gap-if.h"
#include "gap-dbstruct.h"

int trace = 1;


#define thrash_3 thrash


static char *c[] = {
    "A Cuppa Tea and a Lie Down",
    "Hey Spinner",
    "Somebody Ate My Planet",
    "Peter Wang Pud",
    "Nelsh Bailter Space",
    "Tanker",
    "Thermos",
    "The Aim",
    "Robot World",
    "Compiletely",
    "Daddy's Highway",
    "The Law of Things",
    "Fear of God",
    "Dreams of Falling",
    "1st Album",
    "Time Flowing Backwards",
    "World of Sand",
    "Far from the Sun",
    "Kaleidoscope World",
    "Brave Words",
    "Submarine Bells",
    "Soft Bomb",
    "Packet",
    "Compilation",
    "Vehicle",
    "Light",
    "DR503",
    "DR503b",
    "Eusa Kills",
    "Trapdoor Fucking Exit",
    "Nerves",
    "The Fat Controller",
    "Tim Fin",
    "Tuatara",
    "In Love With These Times",
    "Pink Flying Saucers Over the Southern Alps",
    "Getting Older 1981-1991",
    "Freak the Sheep - Volume 1",
    "Freak the Sheep - Volume 2",
    "Songs from the Front Lawn",
    "On Fire",
    "Morse",
    "The Complete",
    "Gordons",
    "Collection",
    "Pure",
    "Stunt Clown",
    "Body Blow",
    "Hellmouth 66",
    "Love Songs",
    "The Size of Food",
    "Precious/Crush/Slip",
    "Bleeding Star",
    "Messages for the Cakekitchen",
    "The Last Great Challenge in a Dull World",
    "Cross Over",
    "Here Come the Cars",
    "Seizure",
    "Croaker",
    "Compilation",
    "Cactii",
    "Three Heads on a Plate",
    "Back of Her Hand",
    "Spaz Out",
    "Into the Hogger",
    "Sceptic Hagfish",
    "Fluid",
    "Alienation",
    "Ionospheres",
    "New Age Savage",
    "Into the Moon",
    "Million Lights",
    "Foaming Out",
    "Amalgam",
    "Sensible",
    "Skeptics III",
    "If I Will I Can",
    "Snapper",
    "Shotgun Blossom",
    "Send You",
    "Sentimental Education",
    "Hard Love Stories",
    "Sombretones",
    "The Very Best of",
    "Gnaw",
    "Hail",
    "Melt",
    "Done",
    "Blow",
    "Hello Cruel World",
    "The Short and Sick of It",
    "Fork Songs",
    "Weeville",
    "Disease - Living Off the Fat of Flying Nun",
    "Cul-de-sac",
    "Beard of Bees",
    "Fish Tales/Swarthy Songs for Swabs",
    "Hellzapoppin",
    "Gritt and Butts",
    "On and On with Lou Reed",
    "Juvenalia",
    "Hallelujah All the Way Home",
    "Some Disenchanted Evening",
    "Ready to Fly",
    "The Mekong Delta Blues",
    "Pile=Up - A Compilation of NZ Music",
    "Killing Capitalism with Kindness",
    "Making Losers Happy",
    "Dynamite Groove",
    "Out of the Yellow Eye"
};



void thrash_0(GapServer *s)
/*
 * for testing
 */
{
#define MAX_RECS 100
#define MAX_VIEWS 500
#define ITERATIONS 10000
    GapClient *client[G_MAX_CLIENTS];
    GView view[MAX_VIEWS];
    int state[MAX_VIEWS];
    int record[MAX_VIEWS];
    int album[MAX_RECS];
    int refs[MAX_RECS];
    int i;
    int Nc;
    int t;
    int TEST_FILE = GAP_DATABASE_FILE;

    trace = 0;
    t = 0;

    Nc = G_Number(c);

    {
	int seed;
	seed = (int)time(NULL);
	if (t) printf("seed=%d\n\n",seed);
	srandom(seed);
    }

    for(i=0;i<MAX_VIEWS;i++) state[i]=(-1);
    for(i=0;i<MAX_RECS;i++) refs[i]=0,album[i]=(-1);


    client[0] = g_connect_client(s,G_LOCK_RW);

    /* Phase 0 - precheck */
    printf("Phase 0 - precheck\n");
    for(i=0;i<MAX_RECS;i++) {
	GRecInfo info;
	char buf[1024];
	int j,found;
	
	/* length to read */
	(void) g_rec_info(client[0],TEST_FILE,i,&info);
	if (info.image != G_NO_IMAGE) {
	    (void) g_fast_read_N(client[0],TEST_FILE,i,buf,info.used);
	    buf[info.used]='\0';
	    found = 0;
	    for(j=0;j<Nc && !found;j++) found = (strcmp(buf,c[j])==0);
	    if (!found) {
		printf("Rec %d garbage - %s\n",i, buf);
		album[i] = -1;
	    } else 
		album[i] = j-1;
	}

    }

    /* Phase 1 - randomly assign albums with the intention to fragment */
    printf("Phase 1 - randomly assign albums with the intention to fragment\n");
    for(i=0;i<ITERATIONS;i++) {
	int v,r;
	GRecInfo info;
	v = random() % MAX_VIEWS;
	switch (state[v]) {
	case -1:
	    r = record[v] = random() % MAX_RECS;
	    (void) g_rec_info(client[0],TEST_FILE,r,&info);
	    if (info.lock >= G_LOCK_RW || random()&01 ) {
		/* read only when we must, and half of the other times */
		if (random()&01) {
		    view[v] = g_lock_N(client[0],TEST_FILE,r,G_LOCK_RO);
		    state[v] = -2;
		    if (t) printf("%6d: %3d %3d lock RO\n",i,v,r);
		} else {
		    char buf[1024];
		    (void) g_fast_read_N(client[0],
					 TEST_FILE,r,buf,256);
		    if (t) printf("%6d: %3d %3d fast reading - %s\n",i,v,r,buf);
		}
	    } else {
		if (random()&01) {
		    view[v] = g_lock_N(client[0],TEST_FILE,r,G_LOCK_RW);
		    state[v] = random() % 10;
		    if (t) printf("%6d: %3d %3d lock RW\n",i,v,r);
		} else {
		    /* write a random comment */
		    album[r] = random() % Nc;
		    (void) g_fast_write_N(client[0],
					  TEST_FILE,r,c[album[r]],strlen(c[album[r]]));
		    if (t) printf("%6d: %3d %3d fast writing - %s\n",i,v,r,c[album[r]]);
		    
		}
	    }
	    refs[r]++;
	    break;
	case 0:
	    r = record[v];
	    (void) g_unlock(client[0],view[v]);
	    if (t) printf("%6d: %3d %3d unlock (RW)\n",i,v,r);
	    state[v]--;
	    refs[r]--;
	    if (album[r] != G_NO_IMAGE) {
		/* check */
		/* GInfo info; */
		char buf[1024];
		/* length to read */
		(void) g_rec_info(client[0],TEST_FILE,r,&info);
		(void) g_fast_read_N(client[0],TEST_FILE,r,buf,info.used);
		buf[info.used]='\0';
		if (t) printf("%6d: %3d %3d fast reading - %s\n",i,v,r,buf);
		if (strcmp(buf,c[album[r]])!=0)
		    printf("Should be - %s\n",c[album[r]]);
	    }
	    break;
	case -2:
	    r = record[v];
	    (void) g_unlock(client[0],view[v]);
	    if (t) printf("%6d: %3d %3d unlock (RO)\n",i,v,r);
	    state[v] = -1;
	    refs[r]--;
	    break;
	default:
	    /* write a random comment */
	    r = record[v];
	    album[r] = random() % Nc;
	    (void) g_write(client[0],view[v],c[album[r]],strlen(c[album[r]]));
	    if (t) printf("%6d: %3d %3d writing - %s\n",i,v,r,c[album[r]]);
	    state[v]--;
	    break;
	}
    }
    
    /* Phase 2 - unlock everything left open */
    printf("Phase 2 - unlock everything left open\n");
    for(i=0;i<MAX_VIEWS;i++) {
	if (state[i] != -1) {
	    int r = record[i];
	    (void) g_unlock(client[0],view[i]);
	    if (t) printf("%6d: %3d %3d unlock\n",-1,i,r);
	}
    }


    /* Phase 3 - check */
    printf("Phase 3 - check\n");
    for(i=0;i<MAX_RECS;i++) {
	GRecInfo info;
	char buf[1024];
	if (album[i]==-1) continue;
	
	/* length to read */
	(void) g_rec_info(client[0],TEST_FILE,i,&info);
	(void) g_fast_read_N(client[0],TEST_FILE,i,buf,info.used);
	buf[info.used]='\0';
	printf("%d - %s\n",i,buf);
	if (strcmp(buf,c[album[i]])!=0)
	    printf("Should be - %s\n",c[album[i]]);

    }

    (void) g_disconnect_client(client[0]);
}





void thrash_1(GapServer *s)
/*
 * for testing
 */
{
    Array views;
    GapClient *client;
    int i;
    int err;
#define MAX_RECS 40000

    views = ArrayCreate(sizeof(GView),MAX_RECS);
    ArrayRef(views,MAX_RECS);

    client = g_connect_client(s,G_LOCK_EX);

    fprintf(stderr,"** About to lock records\n** ");
    for(i=0;i<MAX_RECS;i++) {
	if (i%1000==0) fprintf(stderr,"%d\n",i);
	arr(GView, views, i) = g_lock_N(client,GAP_DATABASE_FILE,i,G_LOCK_EX);
    }

    fprintf(stderr,"\n** About to unlock records\n** ");
    for(i=0;i<MAX_RECS;i++) {
	if (i%1000==0) fprintf(stderr,"%d\n",i);
	err = g_unlock(client,arr(GView,views,i));
    }
    fprintf(stderr,"\n** done\n");

    g_disconnect_client(client);

    ArrayDestroy(views);

}

void thrash_2(GapServer *s)
/*
 * for testing flush
 */
{
    GView view;
    GapClient *client;
    int err;
    char *str = "Test string";
    char *str2 = "What a load of rubbish!";
    char buf[1024];

    client = g_connect_client(s,G_LOCK_EX);
    view = g_lock_N(client,GAP_DATABASE_FILE,0,G_LOCK_EX);
    g_write(client,view,str,strlen(str));
    err = g_flush(client,view);
    g_write(client,view,str2,strlen(str2));
    /*err = g_unlock(client,view);*/
    g_disconnect_client(client);


    client = g_connect_client(s,G_LOCK_EX);
    (void) g_fast_read_N(client,GAP_DATABASE_FILE,0,buf,sizeof(buf));
    printf("wrote=%s\nread=%s\n",str,buf);
    g_disconnect_client(client);

}



void thrash_3(GapServer *s)
/*
 * for testing writev and readv
 */
{
    GView view;
    GapClient *client;
    int err;
    char *str1 = "Test string";
    char *str2 = "What a load of rubbish!";
    char buf[1024];
    GIOVec vin[10];
    GIOVec vout[10];

    client = g_connect_client(s,G_LOCK_EX);
    view = g_lock_N(client,GAP_DATABASE_FILE,0,G_LOCK_EX);

    vout[0].buf = str1; vout[0].len = strlen(str1);
    vout[1].buf = str2; vout[1].len = strlen(str2);
    err = g_writev(client,view,vout,2);

    vin[0].buf = buf; vin[0].len = 4;
    vin[1].buf = buf+100; vin[1].len = 50;
    vin[2].buf = buf+200; vin[2].len = 100;
    err = g_readv(client,view,vin,3);
    printf("vin[0].buf = %s\n",vin[0].buf);
    printf("vin[1].buf = %s\n",vin[1].buf);
    printf("vin[2].buf = %s\n",vin[2].buf);

    err = g_read(client,view,buf,100);
    printf("buf = %s\n",buf);



    err = g_unlock(client,view);
    g_disconnect_client(client);

}

