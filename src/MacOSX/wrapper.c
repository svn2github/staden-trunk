/*
Used to convert a shell script wrapper in STADENROOT/bin into a .app bundle
in the main STADENROOT dir. Handles file associations too.

cc -g wrapper.c -o gap5 -framework CoreFoundation -framework Carbon
*/

#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>
#include <signal.h>

#include <CoreFoundation/CoreFoundation.h>
#include <Carbon/Carbon.h>

#define MAX_ARGS 1000
static char *arg[MAX_ARGS+1];
static int nargs = 1;

/* Handles kAEOpenDocuments events, adding them to the arg list */
static OSErr event_handler(const AppleEvent *event,
			   AppleEvent *reply,
			   long refcon) {
    AEDescList docs;

    if (AEGetParamDesc(event, keyDirectObject, typeAEList, &docs) == noErr) {
	long n = 0;
	char strBuffer[256];
	int i;

	AECountItems(&docs, &n);

	for (i = 0; i < n && nargs < MAX_ARGS; i++) {
	    FSRef ref;
	    if (AEGetNthPtr(&docs, i + 1, typeFSRef, 0, 0,
			    &ref, sizeof(ref), 0) != noErr)
		continue;

	    if (FSRefMakePath(&ref, (UInt8 *)strBuffer, 256) == noErr) {
		arg[nargs++] = strdup(strBuffer);
	    }
	}
    }

    return noErr;
}

/* itimer callback and actually execute the program */
static void timer_callback(int sig) {
    arg[nargs] = NULL;
    execv(arg[0], arg);

    exit(0);
}

int main(int argc, char **argv) {
    char *cp, *prog;
    char buf[1024];
    struct itimerval delay;

    delay.it_interval.tv_sec = 1;
    delay.it_interval.tv_usec = 0;
    delay.it_value.tv_sec = 1;
    delay.it_value.tv_usec = 0;

    /* Parse argv[0] to get the executable we're launching */
    /* dir/foo.app/Contents/MacOS/foo -> dir/bin/foo */
    prog = strrchr(argv[0], '/')+1;
    for (cp = prog-1; *cp != '/'; cp--) ;
    for (cp--; *cp != '/'; cp--) ;
    for (cp--; *cp != '/'; cp--) ;
    for (cp--; *cp != '/'; cp--) ;
    *cp = 0;

    sprintf(buf, "%s/bin/%s", argv[0], prog);
    arg[0] = buf;
    
    /* Briefly monitor for doc open events to get extra command line args */
    AEInstallEventHandler(kCoreEventClass, kAEOpenDocuments,
			  event_handler, 0, false);
    
    signal(SIGALRM, timer_callback);
    setitimer(ITIMER_REAL, &delay, NULL);

    /* Runs forever, but the itimer above will terminate this app */
    RunApplicationEventLoop();

    return 0;
}
