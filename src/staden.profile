#
# Setup file for Staden software.
# This setup for Bourne shell (sh) users or derivatives (eg bash).
#
# It is not normally necessary to source this file, but some cases may
# still require it. If so it should be sourced from your .profile or
# .bash_profile. E.g.
#
#     . /usr/local/staden-2.0/staden.profile
#
#
#echo 'Setting up the Staden software environment...'

#-- Check for valid root
if test "x$STADENROOT" = "x" -o ! -e "$STADENROOT/share/staden/staden.profile"
then
    echo "STADENROOT environment variable not set or is invalid" 1>&2
    echo "Please set and re-source this file." 1>&2
else


#-- Set all other paths relative to the root.
    STADLIB=$STADENROOT/lib/staden;               export STADLIB
    STADTABL=$STADENROOT/share/staden/etc;        export STADTABL
    STADTCL=$STADENROOT/share/staden/tcl;         export STADTCL
    GTAGDB=GTAGDB:$HOME/GTAGDB:$STADTABL/GTAGDB;  export GTAGDB
    
    # Set up PATHS
    [ x"$STADEN_PREPEND" != "x" ] \
    && PATH=$STADENROOT/bin:$PATH \
    || PATH=$PATH:$STADENROOT/bin
    
    if [ "`uname -s`" = "Darwin" ]
    then
        if [ "$DYLD_LIBRARY_PATH" != "" ]
        then
            [ x"$STADEN_PREPEND" != "x" ] \
    	&& DYLD_LIBRARY_PATH=$STADLIB:$DYLD_LIBRARY_PATH \
            || DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$STADLIB
        else
            DYLD_LIBRARY_PATH=$STADLIB
        fi
        export DYLD_LIBRARY_PATH
    else
        if [ "$LD_LIBRARY_PATH" != "" ]
        then
            [ x"$STADEN_PREPEND" != "x" ] \
            && LD_LIBRARY_PATH=$STADLIB:$LD_LIBRARY_PATH \
            || LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$STADLIB:$STADENROOT/lib
        else
            LD_LIBRARY_PATH=$STADLIB:$STADENROOT/lib
        fi
        export LD_LIBRARY_PATH
    fi
    
    #
    # files for gap4
    #
    # Not explicitly needed - defaults to $STADTABL/GTAGDB
    #
    if [ "$GTAGDB" = "" ]
    then
        GTAGDB="GTAGDB:$HOME/GTAGDB:$STADTABL/GTAGDB";
    fi
    export GTAGDB

fi
