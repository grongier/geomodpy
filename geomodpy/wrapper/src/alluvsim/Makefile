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
INCDIR := $(INCDIR)/alluvsim

LIB = numerical_recipes.a
GSLIB = gslib.a
PROGS = \
	alluvsim \
	MAPSpp


all: $(BINDIR) $(INCDIR) $(LIB) $(PROGS)

$(BINDIR) $(INCDIR):
	mkdir -p $@

$(LIB):
	$(MAKE) -e -C numerical_recipes

%: %.for
	$(FC) -std=legacy -ffixed-line-length-none -fbounds-check $(SFLAG) -I$(INCDIR) -J$(INCDIR) $< -o $(BINDIR)/$@ $(LIBDIR)/$(LIB) $(LIBDIR)/$(GSLIB)

clean:
	$(MAKE) clean -C numerical_recipes

delete: clean
	$(RM) $(addprefix $(BINDIR)/, $(PROGS))
