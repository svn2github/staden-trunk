# SYNOPSIS
#
#   AX_LIB_IWIDGETS([MINIMUM-VERSION], [ACTION-IF-TRUE], [ACTION-IF-FALSE])
#
# DESCRIPTION
#
#   This macro will check for the existence of the iwidgets package.
#   Note this is NOT tk itself, but a miscellaneous collection of
#   extra Tk widgets written in tcl/tk.
#   The DIR specified should be either a tcl package path containing iwidgets
#   or the iwidgets install directory itself.
#
#   The following output variables are set using AC_SUBST:
#
#     IWIDGETS_PATH (path including iwidgets<vers>)
#     IWIDGETS_ROOT (parent of IWIDGETS_PATH)
#
# LICENSE
#
#   Copyright (c) 2009 James Bonfield <jkb@sanger.ac.uk>
#
#   Copying and distribution of this file, with or without
#   modification, are permitted in any medium without royalty
#   provided the copyright notice and this notice are preserved.
#   This file is offered as-is, without any warranty.


AC_DEFUN([AX_LIB_IWIDGETS],
[
  AC_ARG_WITH(iwidgets,
	      AC_HELP_STRING([--with-iwidgets=DIR],[look for iwidgets in DIR]),
	      [_iwidgets_with=$withval],[_iwidgets_with="no"])

  _ok=no

  AC_MSG_CHECKING([iwidgets directory])

  # Look in the place we requested and also in some standard best-guess
  # locations.
  for i in $_iwidgets_with/iwidgets* $_iwidgets_with /usr/share/tcl*/iwidgets* /usr/share/iwidgets* /usr/lib/iwidgets* /usr/lib64/iwidgets* /usr/local/tcl*/iwidgets*
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
	p=`egrep 'package ifneeded Iwidgets' "$i/pkgIndex.tcl"`
        v1=`expr "$p" : '.*Iwidgets  *\([[0-9]]*\)'`
        v2=`expr "$p" : '.*Iwidgets  *[[0-9]]*\.\([[0-9]]*\)'`
        v3=`expr "$p" : '.*Iwidgets  *[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
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
        IWIDGETS_PATH=`echo "$i" | sed 's:/*$::'`
        IWIDGETS_ROOT=`echo "$i" | sed 's:/[[^/]]*/*$::'`
	AC_SUBST([IWIDGETS_PATH])
	AC_SUBST([IWIDGETS_ROOT])
        break
      fi
    fi
  done

  # Execute the conditional expressions
  if test "$_ok" = "yes"
  then
     # This is the IF-YES path
     AC_MSG_RESULT([yes ($IWIDGETS_PATH)])
     ifelse([$2],,:,[$2])
  else
     # This is the IF-NO path
     AC_MSG_RESULT([no])
     ifelse([$3],,:,[$3])
  fi
])
