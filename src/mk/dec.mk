#
# WARNING: Out of date
#
CC	  = c89
CFLAGS   += -common
INCLUDES += -I/usr/include/mit
DEFINES	  = -DNOSTRDUP
# may need strdup at any time
EXTRA_LIBS= $(MISCBIN)/misc.a

LD	  = $(CC)
FLIBS	 += -lfor -lutil -lUfor -li -lots -lm -lfor
