#
# indexseqlibs
#

PROGS = \
	$(O)/addnl \
	$(O)/genbentryname1 \
	$(O)/entryname2 \
	$(O)/access4 \
	$(O)/access2 \
	$(O)/genbaccess1 \
	$(O)/title2 \
	$(O)/genbtitle1 \
	$(O)/emblentryname1 \
	$(O)/emblaccess1 \
	$(O)/embltitle1 \
	$(O)/pirentryname1 \
	$(O)/piraccess1 \
	$(O)/piraccess2 \
	$(O)/pirtitle1 \
	$(O)/pirtitle2 \
	$(O)/excludewords \
	$(O)/emblfreetext \
	$(O)/genbfreetext \
	$(O)/pirfreetext \
	$(O)/freetext2 \
	$(O)/freetext4 \
	$(O)/emblauthor \
	$(O)/genbauthor \
	$(O)/pirauthor \
	$(O)/hitNtrg \
	$(O)/gcgentryname1 \
	$(O)/gcgentryname2 \
	$(O)/gcgentryname3 \
	$(O)/gcgbigones \
	$(O)/division

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk


#
# Shared objects
#
OBJS = \
	$(O)/cdromheader.o \
	$(O)/mach-io.o

#
# Dependencies for programs
#
$(O)/addnl : $(O)/addnl.o
	$(CLD) -o $@ $(O)/addnl.o $(LIBSC)

$(O)/genbentryname1 : $(O)/genbentryname1.o
	$(CLD) -o $@ $(O)/genbentryname1.o $(LIBSC)

$(O)/entryname2 : $(O)/entryname2.o $(OBJS)
	$(CLD) -o $@ $(O)/entryname2.o $(OBJS) $(LIBSC)

$(O)/access4 : $(O)/access4.o $(OBJS)
	$(CLD) -o $@ $(O)/access4.o $(OBJS) $(LIBSC)

$(O)/access2 : $(O)/access2.o
	$(CLD) -o $@ $(O)/access2.o $(LIBSC)

$(O)/genbaccess1 : $(O)/genbaccess1.o
	$(CLD) -o $@ $(O)/genbaccess1.o $(LIBSC)

$(O)/title2 : $(O)/title2.o $(OBJS)
	$(CLD) -o $@ $(O)/title2.o $(OBJS) $(LIBSC)

$(O)/genbtitle1 : $(O)/genbtitle1.o
	$(CLD) -o $@ $(O)/genbtitle1.o $(LIBSC)

$(O)/emblentryname1 : $(O)/emblentryname1.o
	$(CLD) -o $@ $(O)/emblentryname1.o $(LIBSC)

$(O)/emblaccess1 : $(O)/emblaccess1.o
	$(CLD) -o $@ $(O)/emblaccess1.o $(LIBSC)

$(O)/embltitle1 : $(O)/embltitle1.o
	$(CLD) -o $@ $(O)/embltitle1.o $(LIBSC)

$(O)/pirentryname1 : $(O)/pirentryname1.o
	$(CLD) -o $@ $(O)/pirentryname1.o $(LIBSC)

$(O)/piraccess1 : $(O)/piraccess1.o
	$(CLD) -o $@ $(O)/piraccess1.o $(LIBSC)

$(O)/piraccess2 : $(O)/piraccess2.o
	$(CLD) -o $@ $(O)/piraccess2.o $(LIBSC)

$(O)/pirtitle1 : $(O)/pirtitle1.o
	$(CLD) -o $@ $(O)/pirtitle1.o $(LIBSC)

$(O)/pirtitle2 : $(O)/pirtitle2.o
	$(CLD) -o $@ $(O)/pirtitle2.o $(LIBSC)

$(O)/excludewords : $(O)/excludewords.o
	$(CLD) -o $@ $(O)/excludewords.o $(LIBSC)

$(O)/emblfreetext.o: freetext.c
	$(CC) $(CFLAGS) -DEMBL -c -o $@ freetext.c
$(O)/emblfreetext : $(O)/emblfreetext.o
	$(CLD) -o $@ $(O)/emblfreetext.o $(LIBSC)

$(O)/genbfreetext.o: freetext.c
	$(CC) $(CFLAGS) -DGENBANK -c -o $@ freetext.c
$(O)/genbfreetext : $(O)/genbfreetext.o
	$(CLD) -o $@ $(O)/genbfreetext.o $(LIBSC)

$(O)/pirfreetext.o: freetext.c
	$(CC) $(CFLAGS) -DPIR -c -o $@ freetext.c
$(O)/pirfreetext : $(O)/pirfreetext.o
	$(CLD) -o $@ $(O)/pirfreetext.o $(LIBSC)

$(O)/freetext2 : $(O)/freetext2.o
	$(CLD) -o $@ $(O)/freetext2.o $(LIBSC)

$(O)/freetext4 : $(O)/freetext4.o $(OBJS)
	$(CLD) -o $@ $(O)/freetext4.o $(OBJS) $(LIBSC)

$(O)/emblauthor.o : author.c
	$(CC) $(CFLAGS) -DEMBL -c -o $@ author.c
$(O)/emblauthor : $(O)/emblauthor.o
	$(CLD) -o $@ $(O)/emblauthor.o $(LIBSC)

$(O)/genbauthor.o : author.c
	$(CC) $(CFLAGS) -DGENBANK -c -o $@ author.c
$(O)/genbauthor : $(O)/genbauthor.o
	$(CLD) -o $@ $(O)/genbauthor.o $(LIBSC)

$(O)/pirauthor.o : author.c
	$(CC) $(CFLAGS) -DPIR -c -o $@ author.c
$(O)/pirauthor : $(O)/pirauthor.o
	$(CLD) -o $@ $(O)/pirauthor.o $(LIBSC)

$(O)/hitNtrg : $(O)/hitNtrg.o $(OBJS)
	$(CLD) -o $@ $(O)/hitNtrg.o $(OBJS) $(LIBSC)

$(O)/division : $(O)/division.o $(OBJS)
	$(CLD) -o $@ $(O)/division.o $(OBJS) $(LIBSC)

$(O)/gcgentryname1 : $(O)/gcgentryname1.o
	$(CLD) -o $@ $(O)/gcgentryname1.o $(LIBSC)

$(O)/gcgentryname2 : $(O)/gcgentryname2.o
	$(CLD) -o $@ $(O)/gcgentryname2.o $(LIBSC)

$(O)/gcgentryname3 : $(O)/gcgentryname3.o
	$(CLD) -o $@ $(O)/gcgentryname3.o $(LIBSC)

$(O)/gcgbigones : $(O)/gcgbigones.o
	$(CLD) -o $@ $(O)/gcgbigones.o $(LIBSC)

DEPEND_OBJ = \
	$(OBJS) \
	$(O)/addnl.o \
	$(O)/genbentryname1.o \
	$(O)/entryname2.o \
	$(O)/access4.o \
	$(O)/access2.o \
	$(O)/genbaccess1.o \
	$(O)/title2.o \
	$(O)/genbtitle1.o \
	$(O)/emblentryname1.o \
	$(O)/emblaccess1.o \
	$(O)/embltitle1.o \
	$(O)/pirentryname1.o \
	$(O)/piraccess1.o \
	$(O)/piraccess2.o \
	$(O)/pirtitle1.o \
	$(O)/pirtitle2.o \
	$(O)/excludewords.o \
	$(O)/genbfreetext.o \
	$(O)/pirfreetext.o \
	$(O)/freetext2.o \
	$(O)/freetext4.o \
	$(O)/emblauthor.o \
	$(O)/genbauthor.o \
	$(O)/pirauthor.o \
	$(O)/hitNtrg.o \
	$(O)/division.o \
	$(O)/gcgentryname1.o \
	$(O)/gcgentryname3.o \
	$(O)/gcgbigones.o \
	$(O)/gcgentryname2.o

include dependencies

install:
	-mkdir $(INSTALLSEQBIN)
	-mkdir $(INSTALLSEQSCRIPT)
	cp $(PROGS) $(INSTALLSEQBIN)
	cp *.script $(INSTALLSEQSCRIPT)
	chmod a+xr $(INSTALLSEQSCRIPT)/*
