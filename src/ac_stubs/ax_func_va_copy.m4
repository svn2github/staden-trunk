#
# Checks whether we the ability to use va_copy().
# AC_DEFINEs HAVE_VA_COPY if you do.
#
AC_DEFUN([AX_FUNC_VA_COPY],
[AC_LANG_PUSH(C)
AC_LINK_IFELSE([AC_LANG_SOURCE(
[[#include <stdarg.h>
void va_test(va_list ap) {
    va_list ap_local;
    va_copy(ap_local, ap);
    va_end(ap_local);
}
int main(void) { return 0; }
]])],
AC_DEFINE([HAVE_VA_COPY],1,[Define to 1 if you have the va_copy() function.]),
)
AC_LANG_POP(C)])
