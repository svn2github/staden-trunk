#
# Setup file for Staden software running on a Sun
# This setup for Bourne shell (sh) users
#
# This file should be source'd from your .profile
# assuming the environmental variable STADENROOT has been set up
# to point to the root directory for the staden software
#
# e.g.
# STADENROOT=/home/BioSW/staden; export STADENROOT
# . $STADENROOT/staden.profile
#
#
#echo 'Setting up the Staden software environment...'

STADENROOT_2002=$STADENROOT
export STADENROOT_2002

#
# set MACHINE to one of alpha/sun/solaris/sgi/linux/macosx.
#
MACHINE=`uname -sr | sed 's/ /-/g;s/SunOS-4.*/sun/;s/IRIX.*/sgi/;s/SunOS-5.*/solaris/;s/OSF.*/alpha/;s/Linux.*/linux/;s/FreeBSD.*/linux/;s/Darwin.*/macosx/'`
export MACHINE

#
# The Digital Unix version is compiled on Digital Unix V4.0. This causes
# problems for the few people using Digital Unix 3.0 and 3.2, but these can
# be worked around by preloading the re-entrant copy of the C library.
#
if [ "`uname -sr | sed 's/\..*//'`" = "OSF1 V3" ]
then
    echo "------------------------------------------------------------------"
    echo "This version of the Staden Package was built on Digital Unix 4.0."
    echo "You may experience problems on this older OS version, but I am"
    echo "setting the _RLD_LIST environment variable in an attempt to solve"
    echo "the known problems."
    echo "------------------------------------------------------------------"
    _RLD_LIST=/shlib/libc_r.so:DEFAULT; export _RLD_LIST
fi


STADTABL=$STADENROOT/tables;	export STADTABL
STADLIB=$STADENROOT/lib;	export STADLIB

# Set up PATHS
PATH=$STADENROOT/$MACHINE-bin:$PATH
export PATH
if [ "$MACHINE" = "macosx" ]
then
    if [ "$DYLD_LIBRARY_PATH" != "" ]
    then
        DYLD_LIBRARY_PATH=$STADLIB/$MACHINE-binaries:$DYLD_LIBRARY_PATH
    else
        DYLD_LIBRARY_PATH=$STADLIB/$MACHINE-binaries
    fi
    export DYLD_LIBRARY_PATH
else
    if [ "$LD_LIBRARY_PATH" != "" ]
    then
        LD_LIBRARY_PATH=$STADLIB/$MACHINE-binaries:$LD_LIBRARY_PATH
    else
        LD_LIBRARY_PATH=$STADLIB/$MACHINE-binaries
    fi
    export LD_LIBRARY_PATH
fi
#
# files for gap4
#
# Not explicitly needed - defaults to $STADTABL/GTAGDB
#
GTAGDB="GTAGDB:$HOME/GTAGDB:$STADTABL/GTAGDB";	export GTAGDB

#
# Find manual pages
# Linux uses /etc/man.config for listing search paths, but you cannot extend
# this - only replace it.
#
if [ "$MACHINE" = "linux" -a "$MANPATH" = "" -a -e /etc/man.config ]
then
defman=`sed -n 's/^MANPATH[ 	][ 	]*//p' /etc/man.config |tr '\012' :`
MANPATH=$STADENROOT/man:$defman
else
MANPATH=$STADENROOT/man:${MANPATH-/usr/man:/usr/local/man:/usr/share/catman};
fi
export MANPATH

#
# Sequence databases
#
#If you wish to use the embl indices you need to remove the # from the
#beginning of the next line
#. $STADTABL/libraries.config.sh

#If you wish to use the srs indices you need to source SRS/etc/prep_srs.sh 
#where SRS is the path to your installation of srs.
#. /pubseq/pubseq/srs5/srs/etc/prep_srs.sh
