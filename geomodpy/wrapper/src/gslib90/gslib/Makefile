BUILDOPT=build

PREFIX = ../../..
BINDIR = $(PREFIX)/bin
INCDIR = $(PREFIX)/include
LIBDIR = $(PREFIX)/lib
ifeq ($(OS), Windows_NT)
	FC = x86_64-w64-mingw32-gfortran
    SFLAG = -static
    ifeq ($(BUILDOPT), prepare)
	    BINDIR := $(BINDIR)/windows
	    INCDIR := $(INCDIR)/windows
	    LIBDIR := $(LIBDIR)/windows
    endif
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
    	FC = gfortran
        SFLAG = -static-libgcc -static-libgfortran
        ifeq ($(BUILDOPT), prepare)
		    BINDIR := $(BINDIR)/linux
		    INCDIR := $(INCDIR)/linux
	    	LIBDIR := $(LIBDIR)/linux
	    endif
    endif
    ifeq ($(UNAME_S), Darwin)
    	FC = gfortran
        SFLAG = -static-libgcc -static-libgfortran
        ifeq ($(BUILDOPT), prepare)
			UNAME_P := $(shell uname -p)
			ifeq ($(UNAME_P), i386)
				BINDIR := $(BINDIR)/macos/intel
				INCDIR := $(INCDIR)/macos/intel
				LIBDIR := $(LIBDIR)/macos/intel
			else
				BINDIR := $(BINDIR)/macos/silicon
				INCDIR := $(INCDIR)/macos/silicon
				LIBDIR := $(LIBDIR)/macos/silicon
			endif
	    endif
    endif
endif
INCDIR := $(INCDIR)/gslib

LIB = gslib.a
OBJS = \
	acorni.o \
	backtr.o \
	beyond.o \
	blue.o \
	chknam.o \
	chktitle.o \
	cova3.o \
	dlocate.o \
	dpowint.o \
	dsortem.o \
	gauinv.o \
	gcum.o \
	getindx.o \
	getz.o \
	green.o \
	hexa.o \
	ksol.o \
	ktsol.o \
	locate.o \
	nscore.o \
	numtext.o \
	ordrel.o \
	picksupr.o \
	powint.o \
	psfill.o \
	psline.o \
	pstext.o \
	rand.o \
	red.o \
	resc.o \
	scal.o \
	setrot.o \
	setsupr.o \
	sortem.o \
	sqdist.o \
	srchsupr.o \
	strlen.o


all: $(LIB)

$(LIB): $(OBJS)
	ar rcs $(LIBDIR)/$(LIB) $(addprefix $(BINDIR)/, $(OBJS))
	$(RM) $(addprefix $(BINDIR)/, $(OBJS))

.SUFFIXES: .o .for

%.o: %.for
	$(FC) -std=legacy -ffixed-line-length-none $(SFLAG) -I$(INCDIR) -J$(INCDIR) -c $< -o $(BINDIR)/$@

clean:
	$(RM) -r $(INCDIR)
	$(RM) $(LIBDIR)/$(LIB)

delete: clean
