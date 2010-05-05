# SYNOPSIS
#
#   AX_LIB_CURSES([ACTION-IF-TRUE], [ACTION-IF-FALSE])
#
# DESCRIPTION
#
#   This macro will check for the existence of curses or ncurses library.
#   It does this by checking for the header file (n)curses.h and curses
#   library file. The location of these may be specified using the
#   --with-curses=DIR command line option (eg --with-curses=/usr/local),
#   using $DIR/include and $DIR/lib for the search path.
#
#   By default it searches for ncurses first and then the old curses library
#   second.
#
#   The following output variables are set using AC_SUBST:
#
#     CURSES_CFLAGS
#     CURSES_LDFLAGS
#     CURSES_VERSION (if MINIMUM-VERSION is not "")
#
#   The C preprocessor symbol HAVE_CURSES_H or HAVE_NCURSES_H will be
#   also defined with AC_DEFINE if a functioning curses library is available.
#   These should be used when determining whether to include <curses.h>
#   or <ncurses.h>. HAVE_LIBCURSES is also defined if either library name
#   has been detected.
#
# LICENSE
#
#   Copyright (c) 2010 James Bonfield <jkb@sanger.ac.uk>
#
#   Copying and distribution of this file, with or without
#   modification, are permitted in any medium without royalty
#   provided the copyright notice and this notice are preserved.
#   This file is offered as-is, without any warranty.

AC_DEFUN([AX_LIB_CURSES],
[
  AC_ARG_WITH(curses,
	      AC_HELP_STRING([--with-curses=DIR],[look for libcurses in DIR]),
	      [_curses_with=$withval],[_curses_with="no"])

  CURSES_ROOT=""
  if test "$_curses_with" != "no"
  then
     if test -f "$_curses_with/include/curses.h"
     then
         CURSES_ROOT=$_curses_with
     fi
  fi

  # Check if it's a working library
  curses_ok=no
  if test "$xCURSES_ROOT" != "x"
  then
    _cppflags=$CPPFLAGS
    CPPFLAGS="$CPPFLAGS -I${CURSES_ROOT}/include"
    _ldflags=$LDFLAGS
    LDFLAGS="$LFDLAGS -L${CURSES_ROOT}/lib"
    AC_LANG_PUSH([C])
    AC_CHECK_LIB(ncurses, mvwprintw, _nl=ncurses;
	[AC_CHECK_HEADER(ncurses.h, curses_ok=yes; _nh=1, 
	    [AC_CHECK_HEADER(curses.h, curses_ok=yes; _nh=0, curses_ok=no)])],
	_nl=curses;
	[AC_CHECK_LIB(curses, mvwprintw,
	    [AC_CHECK_HEADER(curses.h, curses_ok=yes; _nh=0, curses_ok=no)])])
    AC_LANG_POP([C])
    if test "$curses_ok" != "yes"
    then
        # Backout and whinge
        CPPFLAGS=$_cppflags
        LDFLAGS=$_ldflags
        AC_MSG_WARN("--with-curses specified, but non functioning")
    fi

  else
    # Maybe it works "out of the box"?
    AC_CHECK_LIB(ncurses, mvwprintw, _nl=ncurses;
	[AC_CHECK_HEADER(ncurses.h, curses_ok=yes; _nh=1, 
	    [AC_CHECK_HEADER(curses.h, curses_ok=yes; _nh=0, curses_ok=no)])],
	_nl=curses;
	[AC_CHECK_LIB(curses, mvwprintw,
	    [AC_CHECK_HEADER(curses.h, curses_ok=yes; _nh=0, curses_ok=no)])])
  fi

  # perform substitutions
  if test "$curses_ok" = "yes"
  then
      AC_DEFINE(HAVE_LIBCURSES, 1,
             [Define to 1 if you have a working curses/ncurses library.])
      if test "$_nh" = 1
      then
          AC_DEFINE(HAVE_NCURSES_H, 1,
             [Define to 1 if you have <ncurses.h> include file.])
      else
          AC_DEFINE(HAVE_CURSES_H, 1,
             [Define to 1 if you have <curses.h> include file.])
      fi
      if test "$CURSES_ROOT" != ""
      then
          CURSES_LDFLAGS="-L${CURSES_ROOT}/lib -l$_nl"
	  CURSES_CFLAGS="-I${CURSES_ROOT}/include"
      else
          CURSES_LDFLAGS="-l$_nl"
	  CURSES_CFLAGS=
      fi
      AC_SUBST([CURSES_LDFLAGS])
      AC_SUBST([CURSES_CFLAGS])

      if test "$_nl" = "ncurses"
      then
          have_ncurses=yes
      else
          have_curses=yes
      fi
  else
    AC_MSG_WARN("No functioning curses found")
  fi

  # Not sure how many of these are needed, but it's belt-and-braces mode
  AH_TEMPLATE([HAVE_NCURSES_H], [Define if you have <ncurses.h> include file])
  AH_TEMPLATE([HAVE_CURSES_H], [Define if you have <curses.h> include file])
  AM_CONDITIONAL(HAVE_LIBCURSES, test "$curses_ok" = "yes")


  # Execute the conditional expressions
  if test "$curses_ok" = "yes"
  then
     # This is the IF-YES path
     ifelse([$1],,:,[$1])
  else
     # This is the IF-NO path
     ifelse([$2],,:,[$2])
  fi

  # Tidy up
  unset curses_ok
  unset _cppflags
  unset _ldflags
  unset _nh
  unset _nl
])
