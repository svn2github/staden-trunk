#
# Makefile for osp (oligo selection program)
#

PROGS = subosp

SRCROOT=../..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

OSPBIN=$(O)

#
# X VERSION compilation flag
#
VERSION	= SUBVERSION
DEFINES += -D$(VERSION)

INCLUDES_E += $(XAW_INC)

XOSPOBJS=\
	$(OSPBIN)/our_allo.o\
	$(OSPBIN)/analysis.o\
	$(OSPBIN)/get_scores.o\
	$(OSPBIN)/paramIO.o\
	$(OSPBIN)/paramIOX.o\
	$(OSPBIN)/utils.o

#
# Linking - dummy just to force the .o files to compile
#
subosp: $(XOSPOBJS)

DEPEND_OBJ = $(XOSPOBJS)

install:

include dependencies
