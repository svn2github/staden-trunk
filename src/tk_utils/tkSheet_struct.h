    /* Tk basics */
    Tcl_Interp *interp;

    /*
     * A few publically configurable options
     * Several here, and the rest in the Sheet structure.
     */
    int relief;
    int font_width;
    int font_height;
    Tk_3DBorder border;
    XColor *fg;
    XColor *light;
    XColor *indel_fg;

    Tk_ConfigSpec *configSpecs;

    /* Internal bits and pieces */
    int flags;
    int initialised;

    /* The sheet itself */
    Sheet sw;

    /*
     * A function pointer for allowing extensions (eg tkEditor). See the
     * list of "jobs" defined above.
     */
    void (*extension)(ClientData clientData, int job, void *jobData);
    ClientData extensionData;

    /* Keeping track of how many display updates have been requested */
    /* int count; */

    /* The location of a divider bar, or 0 if none */
    int divider;

    /* Whether this is the internal widget to 'grid' the wm on */
    int grid;
