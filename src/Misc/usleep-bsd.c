#include <stdio.h>
#include <signal.h>
#include <sys/time.h>

static int rung_alarm;
static void wakeup() {
    rung_alarm++;
}

/*
 * This sleeps for usecond microseconds. The implementation uses an itimer to
 * set an interrupt and wait for it. Complications arise where there is
 * already one set. Here we subtract our wait from the existing one (when
 * appropriate).
 */
int usleep(unsigned int useconds) {
    struct itimerval value, ovalue;
    struct sigvec    vec,   ovec;
    int old_mask;

    if (!useconds)
	return 0;
    
    /* Initialise our itimer to zero time */
    timerclear(&value.it_interval);
    timerclear(&value.it_value);
    
    /* Fetch the old itimer */
    if (-1 == setitimer(ITIMER_REAL, &value, &ovalue))
	return -1;

    /* Set up our timer structure */
    value.it_value.tv_sec  = useconds / 100000;
    value.it_value.tv_usec = useconds % 100000;
    
    /*
     * If it's set and is set for a time further in the future than ours then
     * we subtract our time from its. (Yuk this is a horrid macro). Otherwise
     * we set our timer to the old one (so we wake up in time ;)
     */
    if (timercmp( &ovalue.it_value, &value.it_value, > )) {
	ovalue.it_value.tv_usec -= value.it_value.tv_usec;
	ovalue.it_value.tv_sec  -= value.it_value.tv_sec;
	if (ovalue.it_value.tv_usec < 0) {
	    ovalue.it_value.tv_usec += 1000000;
	    ovalue.it_value.tv_sec--;
	} else {
	    value.it_value = ovalue.it_value;
	    ovalue.it_value.tv_sec = 0;
	    ovalue.it_value.tv_usec = 50000; /* HACK arbitrary minimal pause */
	}
    }

    /* Initialise our sigvec structure */
    vec.sv_handler = wakeup;
    vec.sv_mask = 0;
    vec.sv_onstack = 0;

    /* Create an alarm */
    (void)sigvec(SIGALRM, &vec, NULL);

    /* Wait for the alarm to go off */
    old_mask = sigblock(sigmask(SIGALRM));
    rung_alarm = 0;
    while (!rung_alarm)
	/* wait for a SIGALRM - but also allow other signals to be caught */
	sigpause(old_mask &~ sigmask(SIGALRM));

    /* tidy up */
    (void)sigvec(SIGALRM, &ovec, NULL);
    (void)sigsetmask(old_mask);
    (void)setitimer(ITIMER_REAL, &ovalue, NULL);

    return 0;
}
