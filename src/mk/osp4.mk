#
# Makefile for osp (oligo selection program)
#

PROGS = subosp

SRCROOT=../..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

OSP4BIN=$(O)

#
# X VERSION compilation flag
#
VERSION	= SUBVERSION
DEFINES += -D$(VERSION)
INCLUDES = -I. $(MISC_INC)

XOSPOBJS=\
	$(OSP4BIN)/analysis.o

#
# Linking - dummy just to force the .o files to compile
#
subosp: $(XOSPOBJS)

install:

DEPEND_OBJ = $(XOSPOBJS)

include dependencies
