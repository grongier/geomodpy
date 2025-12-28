BUILDOPT=prepare

PREFIX = ../..
BINDIR = $(PREFIX)/bin
INCDIR = $(PREFIX)/include
LIBDIR = $(PREFIX)/lib
ifeq ($(OS), Windows_NT)
	MAKE = mingw32-make
	FC = x86_64-w64-mingw32-gfortran
    SFLAG = -static
    ifeq ($(BUILDOPT), prepare)
	    BINDIR := $(BINDIR)/windows
	    INCDIR := $(INCDIR)/windows
	    LIBDIR := $(LIBDIR)/windows
    endif
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S), Linux)
    	MAKE = make
    	FC = gfortran
        SFLAG = -static
        ifeq ($(BUILDOPT), prepare)
		    BINDIR := $(BINDIR)/linux
		    INCDIR := $(INCDIR)/linux
	    	LIBDIR := $(LIBDIR)/linux
	    endif
    endif
    ifeq ($(UNAME_S), Darwin)
    	MAKE = make
    	FC = gfortran
        SFLAG = -static-libgfortran -static-libgcc
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
PROGS = \
	addcoord \
	anneal \
	backtr \
	bicalib \
	bigaus \
	bivplt \
	cokb3d \
	declus \
	draw \
	ellipsim \
	gam \
	gamv \
	gtsim \
	histplt \
	histsmth \
	ik3d \
	kb2d \
	kt3d \
	locmap \
	lusim \
	nscore \
	pfsim \
	pixelplt \
	plotem \
	postik \
	postsim \
	probplt \
	qpplt \
	rotcoord \
	sasim \
	scatplt \
	scatsmth \
	sgsim \
	sisim \
	sisim_gs \
	sisim_lm \
	trans \
	vargplt \
	varmap \
	vmodel


all: $(BINDIR) $(INCDIR) $(LIBDIR) $(LIB) $(PROGS)

$(BINDIR) $(INCDIR) $(LIBDIR):
	mkdir -p $@

$(LIB):
	cd gslib; $(MAKE) -e

%: %.for
	$(FC) -std=legacy -ffixed-line-length-none $(SFLAG) -I$(INCDIR) -J$(INCDIR) $< -o $(BINDIR)/$@ $(LIBDIR)/$(LIB)

clean:
	cd gslib; $(MAKE) clean; cd ..

delete: clean
	$(RM) $(addprefix $(BINDIR)/, $(PROGS))
