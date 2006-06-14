# AC_CHECK_IOLIB autoconf macro to find io_lib.
# If found it defines HAVE_IOLIB and sets IOLIB_CPPFLAGS and IOLIB_LDFLAGS.

AC_DEBUF([AC_CHECK_IOLIB],
[
  AC_ARG_WITH(libiolib, AC_HELP_STRING([--with-io_lib=DIR],
	                                [Look for io_lib root in DIR]),
	      [_iolib_with=$withval],
	      [_iolib_with=ifelse([$1],,[yes],[$1])])

  # Defaults to enabled
  if test "$_iolib_with" != "no"
  then
    # Identify the location of iolib-config
    if [ -d "$_iolib_with" ]
    then
      AC_PATH_PROG([_iolib_config], ["$_iolib_with/bin/iolib-config"])
      # defaults incase iolib-config isn't found
      IOLIB_CPPFLAGS="-I$withval/include"
      IOLIB_LDFLAGS="-L$withval/lib -lread"
    else
      AC_PATH_PROG([_iolib_config], [iolib-config])
    fi

    # Configure IOLIB_CPPFLAGS and IOLIB_LDFLAGS
    if [ x$_iolib_config != "x" ]
    then
      IOLIB_CPPFLAGS=`$_iolib-config --cflags`
      IOLIB_LDFLAGS=`$_iolib-config --libs`
      AC_DEFINE(HAVE_IOLIB)
      AC_SUBST(IOLIB_CPPFLAGS)
      AC_SUBST(IOLIB_LDFLAGS)
    fi
  fi
])dnl
