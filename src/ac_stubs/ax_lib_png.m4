# SYNOPSIS
#
#   AX_LIB_PNG([MINIMUM-VERSION], [ACTION-IF-TRUE], [ACTION-IF-FALSE])
#
# DESCRIPTION
#
#   This macro will check for the existence of png library.
#   It does this by checking for the header file png.h and the png library
#   object file. The location of these may be specified using the
#   --with-png=DIR command line option (eg --with-png=/usr/local),
#   using $DIR/include and $DIR/lib for the search path.
#
#   The following output variables are set using AC_SUBST:
#
#     PNG_CPPFLAGS
#     PNG_LDFLAGS
#
#   The C preprocessor symbol HAVE_PNG will be also defined with
#   AC_DEFINE if a functioning samtools is available.
#
# LICENSE
#
#   Copyright (c) 2009 James Bonfield <jkb@sanger.ac.uk>
#
#   Copying and distribution of this file, with or without
#   modification, are permitted in any medium without royalty
#   provided the copyright notice and this notice are preserved.
#   This file is offered as-is, without any warranty.


#
# Check for PNG
#
AC_DEFUN([AX_LIB_PNG],
[
  AC_ARG_WITH(png,
  	    AC_HELP_STRING([--with-png], [Look for png inc/lib in DIR]),
	    [PNG_ROOT="$withval"], [PNG_ROOT=""])
  
  png_ok=no
  if test "x$PNG_ROOT" != "x"
  then
    _cppflags=$CPPFLAGS
    CPPFLAGS="$CPPFLAGS -I${PNG_ROOT}/include"
    _ldflags=$LDFLAGS
    LDFLAGS="$LFDLAGS -L${PNG_ROOT}/lib"
    AC_LANG_PUSH([C])
    AC_CHECK_LIB(png, png_create_write_struct,
	[AC_CHECK_HEADER(png.h, png_ok=yes, png_ok=no)])
    AC_LANG_POP([C])

    if test "$png_ok" != "yes"
    then
        # Backout and whinge
        CPPFLAGS=$_cppflags
        LDFLAGS=$_ldflags
        AC_MSG_WARN("--with-png specified, but non-functioning")
    fi

  else
    # Maybe it works "out of the box"?
    AC_CHECK_LIB(png, png_create_write_struct,
	[AC_CHECK_HEADER(png.h, png_ok=yes, png_ok=no)])
  fi

  # Check version
  if test "x$1" != "x" && test "$png_ok" = "yes"
  then
      AC_MSG_CHECKING([if png version >= $1])

      for i in "$PNG_ROOT" "/usr/include" "/usr/share/include" "/usr/local/include"
      do
	  if test -f "$i/png.h"
	  then
	      PNG_VERSION=`sed -n 's/.*#define *PNG_LIBPNG_VER_STRING *"\([^"]*\)"/\1/p' "$i/png.h"`
	      break
	  fi
      done

      v1=`expr "$1" : '\([[0-9]]*\)'`
      v2=`expr "$1" : '[[0-9]]*\.\([[0-9]]*\)'`
      v3=`expr "$1" : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
      want_vers=`expr "${v1:-0}" "*" 1000000 + "${v2:-0}" "*" 1000 + "${v3:-0}"`

      v1=`expr "${PNG_VERSION:-}" : '\([[0-9]]*\)'`
      v2=`expr "${PNG_VERSION:-}" : '[[0-9]]*\.\([[0-9]]*\)'`
      v3=`expr "${PNG_VERSION:-}" : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
      have_vers=`expr "${v1:-0}" "*" 1000000 + "${v2:-0}" "*" 1000 + "${v3:-0}"`
      if test `expr "$have_vers" ">=" "$want_vers"` = "1"
      then
          AC_MSG_RESULT([yes])
          AC_SUBST([PNG_VERSION])
      else
          AC_MSG_RESULT([no])
	  png_ok="no"
      fi
  fi

  # perform substitutions
  if test "$png_ok" = "yes"
  then
      AC_DEFINE(HAVE_PNG, 1,
         [Define to 1 if you have a functional libpng.])
      if test "x$PNG_ROOT" != "x"
      then
          PNG_LDFLAGS="-L${PNG_ROOT}/lib -lpng"
	  PNG_CFLAGS="-I${PNG_ROOT}/include"
      else
          PNG_LDFLAGS="-lpng"
	  PNG_CFLAGS=
      fi
      AC_SUBST([PNG_LDFLAGS])
      AC_SUBST([PNG_CFLAGS])
  else
    AC_MSG_WARN("No functioning png found")
  fi

  # Not sure how many of these are needed, but it's belt-and-braces mode
  AH_TEMPLATE([HAVE_PNG], [Define if png is installed])
  AM_CONDITIONAL(HAVE_PNG, test "$png_ok" = "yes")


  # Execute the conditional expressions
  if test "$png_ok" = "yes"
  then
     # This is the IF-YES path
     ifelse([$2],,:,[$2])
  else
     # This is the IF-NO path
     ifelse([$3],,:,[$3])
  fi

  # Tidy up
  unset png_ok
  unset _cppflags
  unset _ldflags
])
