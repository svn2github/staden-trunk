# SYNOPSIS
#
#   AX_LIB_ITK([MINIMUM-VERSION], [ACTION-IF-TRUE], [ACTION-IF-FALSE])
#
# DESCRIPTION
#
#   This macro will check for the existence of the itk package.
#   The DIR specified should be either a tcl package path containing itk
#   or the itk install directory itself.
#
#   The following output variables are set using AC_SUBST:
#
#     ITK_PATH (path including itk<vers>)
#     ITK_ROOT (parent dir of ITK_PATH)
#
# LICENSE
#
#   Copyright (c) 2009 James Bonfield <jkb@sanger.ac.uk>
#
#   Copying and distribution of this file, with or without
#   modification, are permitted in any medium without royalty
#   provided the copyright notice and this notice are preserved.
#   This file is offered as-is, without any warranty.


AC_DEFUN([AX_LIB_ITK],
[
  AC_ARG_WITH(itk,
	      AC_HELP_STRING([--with-itk=DIR],[look for itk in DIR]),
	      [_itk_with=$withval],[_itk_with="no"])

  _ok=no

  AC_MSG_CHECKING([itk directory])

  # Look in the place we requested and also in some standard best-guess
  # locations.
  for i in $_itk_with/itk* $_itk_with /usr/share/tcl*/itk* /usr/local/tcl*/itk* /usr/lib/tcl*/itk* /usr/lib/itk* /usr/lib64/tcl*/itk* /usr/lib64/itk*
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
	p=`egrep 'package ifneeded Itk' "$i/pkgIndex.tcl"`
        v1=`expr "$p" : '.*Itk  *\([[0-9]]*\)'`
        v2=`expr "$p" : '.*Itk  *[[0-9]]*\.\([[0-9]]*\)'`
        v3=`expr "$p" : '.*Itk  *[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
        have_vers=`expr "${v1:-0}" "*" 1000000 + "${v2:-0}" "*" 1000 + "${v3:-0}"`
	
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
        ITK_PATH=`echo "$i" | sed 's:/*$::'`
        ITK_ROOT=`echo "$i" | sed 's:/[[^/]]*/*$::'`
	AC_SUBST([ITK_PATH])
	AC_SUBST([ITK_ROOT])
        break
      fi
    fi
  done

  # Execute the conditional expressions
  if test "$_ok" = "yes"
  then
     # This is the IF-YES path
     AC_MSG_RESULT([yes ($ITK_PATH)])
     ifelse([$2],,:,[$2])
  else
     # This is the IF-NO path
     AC_MSG_RESULT([no])
     ifelse([$3],,:,[$3])
  fi
])
