#ifndef _locks_h
#define _locks_h
/*
    Title: 	 locks

    File: 	 locks.h
    Purpose:	 General routines for locking resources
    Last update: 
*/



#ifdef _MSC_VER
#  ifdef BUILDING_MISC_DLL
#    define MISC_EXPORT __declspec(dllexport)
#  else
#    define MISC_EXPORT __declspec(dllimport)
#  endif
#else
#  define MISC_EXPORT
#endif

typedef struct {
	int count;
	int max;
	} SemaphoreStruct, *Semaphore, *Flag;


extern MISC_EXPORT Semaphore activeLock;

#ifdef noddy

#define semaphoreCreate(S,N) (S=(Semaphore)malloc(sizeof(Semaphore)))!=NULL)?S->count=0,S->max=N,S:S
#define semaphoreGrab(S) (S->count==S->max)?0:S->count++,1
#define semaphoreRelease(S) (S->count==0)?0:S->count--,1
#define semaphoreGrabN(S,N) (S->count+N>S->max)?0:(S-count+=N),1
#define semaphoreRealeaseN(S,N) (S->count<N)?0:(S->count-=N),1
#define semaphoreGrabExclusive(S) semaphoreGrabN(S,S->max)
#define semaphoreUsed(S) S->count
#define semaphoreFree(S) (S->count==0)
#define flagCreate (Flag) semaphoreCreate(1)
#define flagSet(F) semaphoreGrab((Semaphore) F)
#define flagUnset(F) semaphoreRelease((Semaphore) F)
#define flagUsed(F) F->count
#define flagFree(F) (F->count==0)

#else
extern Semaphore semaphoreCreate(int max);
extern int semaphoreGrab(Semaphore sem);
extern int semaphoreRelease(Semaphore sem);
extern int semaphoreGrabN(Semaphore sem, int n);
extern int semaphoreReleaseN(Semaphore sem, int n);
extern int semaphoreGrabExclusive(Semaphore sem);
extern int semaphoreUsed(Semaphore sem);
extern int semaphoreFree(Semaphore sem);
extern Flag flagCreate(void);
extern int flagSet(Flag flag);
extern int flagUnset(Flag flag);
extern int flagUsed(Flag flag);
extern int flagFree(Flag flag);
#endif

#endif /* _locks_h */
