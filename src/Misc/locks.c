/*
    Title: 	 locks

    File: 	 locks.c
    Purpose:	 General routines for locking resources
    Last update: 
*/

#include "locks.h"
#include <stdio.h>
#include <stdlib.h>

/* 6/1/99 johnt - must assign globals to force export under windows */
Semaphore activeLock=0;

Semaphore semaphoreCreate(int max)
{
    Semaphore sem;

    if ((sem = (Semaphore) malloc (sizeof(SemaphoreStruct)))!=NULL) {
	sem->count = 0;
	sem->max = max;
    }
    return sem;
}

int semaphoreGrab(Semaphore sem)
{
    return (sem->count==sem->max)?0:sem->count++,1;
}

int semaphoreRelease(Semaphore sem)
{
    return (sem->count==0)?0:sem->count--,1;
}

int semaphoreGrabN(Semaphore sem, int n)
{
    return (sem->count+n>sem->max)?0:(sem->count+=n),1;
}

int semaphoreReleaseN(Semaphore sem, int n)
{
    return (sem->count<n)?0:(sem->count-=n),1;
}

int semaphoreGrabExclusive(Semaphore sem)
{
    return semaphoreGrabN(sem, sem->max);
}

int semaphoreUsed(Semaphore sem)
{
    return sem->count;
}

int semaphoreFree(Semaphore sem)
{
    return (sem->count == 0);
}

Flag flagCreate(void)
{
    return (Flag) semaphoreCreate(1);
}

int flagSet(Flag flag)
{
    return semaphoreGrab((Semaphore) flag);
}

int flagUnset(Flag flag)
{
    return semaphoreRelease((Semaphore) flag);
}

int flagUsed(Flag flag)
{
    return (flag->count);
}

int flagFree(Flag flag)
{
    return (flag->count == 0);
}
