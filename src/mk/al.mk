#
# Warning. Out of date
#
# If using gmake - requires version 3.64 or higher

CC	= fxc
CFLAGS += -uniproc
DEFINES	= -DNOSTRDUP
# may need strdup at any time
EXTRA_LIBS= $(MISCBIN)/misc.a

F77	= fortran
FFLAGS += -uniproc

LD	  = $(CC)
LDFLAGS  += -uniproc
CLDFLAGS += -uniproc
