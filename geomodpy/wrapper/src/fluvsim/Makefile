BUILDOPT=prepare

PREFIX = ../..
BINDIR = $(PREFIX)/bin
INCDIR = $(PREFIX)/include
ifeq ($(OS), Windows_NT)
	FC = x86_64-w64-mingw32-gfortran
    SFLAG = -static
    ifeq ($(BUILDOPT), prepare)
	    BINDIR := $(BINDIR)/windows
	    INCDIR := $(INCDIR)/windows
    endif
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S), Linux)
    	FC = gfortran
        SFLAG = -static
        ifeq ($(BUILDOPT), prepare)
		    BINDIR := $(BINDIR)/linux
		    INCDIR := $(INCDIR)/linux
	    endif
    endif
    ifeq ($(UNAME_S), Darwin)
    	FC = gfortran
        SFLAG = -static-libgfortran -static-libgcc
        ifeq ($(BUILDOPT), prepare)
			UNAME_P := $(shell uname -p)
			ifeq ($(UNAME_P), i386)
				BINDIR := $(BINDIR)/macos/intel
				INCDIR := $(INCDIR)/macos/intel
			else
				BINDIR := $(BINDIR)/macos/silicon
				INCDIR := $(INCDIR)/macos/silicon
			endif
	    endif
    endif
endif
INCDIR := $(INCDIR)/fluvsim

PROG = fluvsim


all: $(BINDIR) $(INCDIR) $(PROG)

$(BINDIR) $(INCDIR):
	mkdir -p $@

%: %.for
	$(FC) -std=legacy -ffixed-line-length-none $(SFLAG) -I$(INCDIR) -J$(INCDIR) $< -o $(BINDIR)/$@

clean:
	$(RM) -r $(INCDIR)

delete: clean
	$(RM) $(BINDIR)/$(PROG)
