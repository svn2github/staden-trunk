#
# Setup file for Staden software running on a Sun
# This setup for c-shell (csh) users
#
# This file should be source'd from your .login
# assuming the environmental variable STADENROOT has been set up
# to point to the root directory for the staden software
#
# e.g.
# setenv STADENROOT /home/BioSW/staden
# source $STADENROOT/staden.login
#
#
#echo 'Setting up the Staden software environment...'

setenv STADENROOT_1_6 "$STADENROOT"

#
# set MACHINE to one of alpha/sun/solaris/sgi/linux.
#
setenv MACHINE `uname -srm | sed 's/ /-/g;s/SunOS-4.*/sun/;s/IRIX.*/sgi/;s/SunOS-5.*/solaris/;s/OSF.*alpha/alpha/;s/Linux.*i.86/linux/;s/Linux.*ia64/linux-ia64/;s/Linux.*/linux/;s/FreeBSD.*/linux/;s/Darwin.*/macosx/'`

#
# The Digital Unix version is compiled on Digital Unix V4.0. This causes
# problems for the few people using Digital Unix 3.0 and 3.2, but these can
# be worked around by preloading the re-entrant copy of the C library.
#
if ( "`uname -sr | sed 's/\..*//'`" == "OSF1 V3" ) then
    echo "------------------------------------------------------------------"
    echo "This version of the Staden Package was built on Digital Unix 4.0."
    echo "You may experience problems on this older OS version, but I am"
    echo "setting the _RLD_LIST environment variable in an attempt to solve"
    echo "the known problems."
    echo "------------------------------------------------------------------"
    setenv _RLD_LIST /shlib/libc_r.so:DEFAULT
endif

setenv STADTABL	$STADENROOT/tables
setenv STADLIB  $STADENROOT/lib

# Set up PATHS
set path = ($STADENROOT/$MACHINE-bin $path)
if ( "$MACHINE" == "macosx" ) then
if ($?DYLD_LIBRARY_PATH) then
    setenv DYLD_LIBRARY_PATH "${STADLIB}/${MACHINE}-binaries:${DYLD_LIBRARY_PATH}"
else
    setenv DYLD_LIBRARY_PATH "${STADLIB}/${MACHINE}-binaries"
endif
else
if ($?LD_LIBRARY_PATH) then
    setenv LD_LIBRARY_PATH "${STADLIB}/${MACHINE}-binaries:${LD_LIBRARY_PATH}"
else
    setenv LD_LIBRARY_PATH "${STADLIB}/${MACHINE}-binaries"
endif
endif

#
# files for gap4
#
# Not explicitly needed - defaults to $STADTABL/GTAGDB
#
setenv GTAGDB   "GTAGDB:${HOME}/GTAGDB:${STADTABL}/GTAGDB"

#
# Find manual pages
# If MANPATH is not set then we have to get the default one somehow so
# we can append to it.
#
if (! $?MANPATH) then
    # Use the manpath program if available
    if ( -x /usr/bin/manpath ) then
        setenv MANPATH `/usr/bin/manpath`
    else
        setenv MANPATH "/usr/man:/usr/local/man:/usr/share/catman:/usr/share/man"
    endif
endif
# Finally append our own manual path
setenv MANPATH ${MANPATH}:${STADENROOT}/man
