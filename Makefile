export BUILDOPT=prepare

PREFIX = ./geomodpy/wrapper
SRCDIR = $(PREFIX)/src
BINDIR = $(PREFIX)/bin
INCDIR = $(PREFIX)/include
LIBDIR = $(PREFIX)/lib
ifeq ($(OS), Windows_NT)
    ifeq ($(BUILDOPT), prepare)
	    BINDIR := $(BINDIR)/windows
	    INCDIR := $(INCDIR)/windows
	    LIBDIR := $(LIBDIR)/windows
    endif
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S), Linux)
        ifeq ($(BUILDOPT), prepare)
		    BINDIR := $(BINDIR)/linux
		    INCDIR := $(INCDIR)/linux
	    	LIBDIR := $(LIBDIR)/linux
	    endif
    endif
    ifeq ($(UNAME_S), Darwin)
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

FOLDERS = \
	gslib90 \
	fluvsim \
	alluvsim \
	snesim \
	fracnet


all:
	for FOLDER in $(FOLDERS); do \
		pwd && cd $(SRCDIR)/$$FOLDER && $(MAKE) -e && pwd && cd ../../../.. && pwd; \
	done

clean:
	for FOLDER in $(FOLDERS); do \
		cd $(SRCDIR)/$$FOLDER && $(MAKE) -e clean && cd ../../../..; \
	done
ifeq ($(BUILDOPT),prepare)
	$(RM) -r $(INCDIR)
	$(RM) -r $(LIBDIR)
endif

delete:
	for FOLDER in $(FOLDERS); do \
		cd $(SRCDIR)/$$FOLDER && $(MAKE) -e delete && cd ../../../..; \
	done
ifeq ($(BUILDOPT),prepare)
	$(RM) -r $(BINDIR)
	$(RM) -r $(INCDIR)
	$(RM) -r $(LIBDIR)
endif
