# SYNOPSIS
#
#   AX_LIB_ITCL([MINIMUM-VERSION], [ACTION-IF-TRUE], [ACTION-IF-FALSE])
#
# DESCRIPTION
#
#   This macro will check for the existence of the itcl package.
#   The DIR specified should be either a tcl package path containing itcl
#   or the itcl install directory itself.
#
#   The following output variables are set using AC_SUBST:
#
#     ITCL_PATH (path including itcl<vers>)
#     ITCL_ROOT (parent dir of ITCL_PATH)
#
# LICENSE
#
#   Copyright (c) 2009 James Bonfield <jkb@sanger.ac.uk>
#
#   Copying and distribution of this file, with or without
#   modification, are permitted in any medium without royalty
#   provided the copyright notice and this notice are preserved.
#   This file is offered as-is, without any warranty.


AC_DEFUN([AX_LIB_ITCL],
[
  AC_ARG_WITH(itcl,
	      AC_HELP_STRING([--with-itcl=DIR],[look for itcl in DIR]),
	      [_itcl_with=$withval],[_itcl_with="no"])

  _ok=no

  AC_MSG_CHECKING([itcl directory])

  # Look in the place we requested and also in some standard best-guess
  # locations.
  for i in $_itcl_with $_itcl_with/itcl* /usr/share/tcl*/itcl* /usr/local/tcl*/itcl* /usr/lib/tcl*/itcl* /usr/lib/itcl* /usr/lib64/tcl*/itcl* /usr/lib64/itcl*
  do
    if test -e "$i/pkgIndex.tcl"
    then
      # Check version
      if test "x$1" != "x"
      then
        v1=`expr "$1" : '\([[0-9]]*\)'`
        v2=`expr "$1" : '[[0-9]]*\.\([[0-9]]*\)'`
        v3=`expr "$1" : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
        want_vers=`expr "${v1:-0}" "*" 1000000 + "${v2:-0}" "*" 1000 + "${v3:-0}"`
        
        # Get version from pkgIndex.tcl
	p=`egrep 'package ifneeded Itcl' "$i/pkgIndex.tcl"`
        v1=`expr "$p" : '.*Itcl  *\([[0-9]]*\)'`
        v2=`expr "$p" : '.*Itcl  *[[0-9]]*\.\([[0-9]]*\)'`
        v3=`expr "$p" : '.*Itcl  *[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
        have_vers=`expr "${v1:-0}" "*" 1000000 + "${v2:-0}" "*" 1000 + "${v3:-0}"`
	
	echo vers=/$have_vers/
        if test `expr "$have_vers" ">=" "$want_vers"` = "1"
	then
	  _ok="yes"
	else
	  _ok="no"
	fi
      else
	_ok="yes"
      fi

      if test "x$_ok" = "xyes"
      then
        ITCL_PATH=`echo "$i" | sed 's:/*$::'`
	ITCL_ROOT=`echo "$i" | sed 's:/[[^/]]*/*$::'`
	AC_SUBST([ITCL_PATH])
	AC_SUBST([ITCL_ROOT])
        break
      fi
    fi
  done

  # Execute the conditional expressions
  if test "$_ok" = "yes"
  then
     # This is the IF-YES path
     AC_MSG_RESULT([yes ($ITCL_PATH)])
     ifelse([$2],,:,[$2])
  else
     # This is the IF-NO path
     AC_MSG_RESULT([no])
     ifelse([$3],,:,[$3])
  fi
])
