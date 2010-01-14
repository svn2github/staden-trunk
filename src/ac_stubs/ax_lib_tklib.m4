# SYNOPSIS
#
#   AX_LIB_TKLIB([MINIMUM-VERSION], [ACTION-IF-TRUE], [ACTION-IF-FALSE])
#
# DESCRIPTION
#
#   This macro will check for the existence of the tklib package.
#   Note this is NOT tk itself, but a miscellaneous collection of
#   extra Tk widgets written in tcl/tk.
#   The DIR specified should be the root of the packages, ie the directory
#   containing tablelist, tooltip, datefield etc.
#
#   The following output variables are set using AC_SUBST:
#
#     TKLIB_PATH
#
# LICENSE
#
#   Copyright (c) 2009 James Bonfield <jkb@sanger.ac.uk>
#
#   Copying and distribution of this file, with or without
#   modification, are permitted in any medium without royalty
#   provided the copyright notice and this notice are preserved.
#   This file is offered as-is, without any warranty.


AC_DEFUN([AX_LIB_TKLIB],
[
  AC_ARG_WITH(tklib,
	      AC_HELP_STRING([--with-tklib=DIR],[look for tklib in DIR]),
	      [_tklib_with=$withval],[_tklib_with="no"])

  _ok=no

  AC_MSG_CHECKING([tklib directory])

  # Look in the place we requested and also in some standard best-guess
  # locations.
  for i in "$_tklib_with" /usr/share/tcl*/tklib* /usr/share/tklib* /usr/local/tklib*
  do
    if test -d "$i/tablelist"
    then
      # Check version
      if test "x$1" != "x"
      then
        v1=`expr "$1" : '\([[0-9]]*\)'`
        v2=`expr "$1" : '[[0-9]]*\.\([[0-9]]*\)'`
        v3=`expr "$1" : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
        want_vers=`expr "${v1:-0}" "*" 1000000 + "${v2:-0}" "*" 1000 + "${v3:-0}"`
        
        # Guess version based on filename suffix
        v1=`expr "$i" : '.*tklib\([[0-9]]*\)'`
        v2=`expr "$i" : '.*tklib[[0-9]]*\.\([[0-9]]*\)'`
        v3=`expr "$i" : '.*tklib[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
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
        TKLIB_PATH=$i
	AC_SUBST([TKLIB_PATH])
        break
      fi
    fi
  done

  # Execute the conditional expressions
  if test "$_ok" = "yes"
  then
     # This is the IF-YES path
     AC_MSG_RESULT([yes ($TKLIB_PATH)])
     ifelse([$2],,:,[$2])
  else
     # This is the IF-NO path
     AC_MSG_RESULT([no])
     ifelse([$3],,:,[$3])
  fi
])
