# ---------------------------------------------------
# MAKEFILE FOR LIBIERS10++ LIBRARY
# TOP-LEVEL MAKEFILE
# ----------------------------------------------------

## Package-related info (set by autoconf)
## ===================================================
package     = libiers10++
version     = 1.0
tarname     = libiers10++
distdir     = $(tarname)-$(version)

## Path(s)-related info (set by autoconf)
## ===================================================
prefix      = /usr/local
exec_prefix = ${prefix}
datarootdir = ${prefix}/share
includedir  = ${prefix}/include
docdir      = ${datarootdir}/doc/${package}v${version}
libdir      = ${prefix}/lib
abs_builddir=$(shell pwd)

## Define and export a temporary build folder
## ===================================================
export builddir = $(abs_builddir)/obj

SUBDIRS = src

all: $(SUBDIRS)

info:
	@echo "======================================"
	@echo "Echoing top-level info"
	@echo "======================================"
	@echo "package     = $(package)"
	@echo "version     = $(version)"
	@echo "tarname     = $(tarname)"
	@echo "distdir     = $(distdir)"
	@echo "prefix      = $(prefix)"
	@echo "exec_prefix = $(exec_prefix)"
	@echo "datarootdir = $(datarootdir)"
	@echo "includedir  = $(includedir)"
	@echo "docdir      = $(docdir)"
	@echo "libdir      = $(libdir)"
	@echo "srcdir      = $(srcdir)"
	@echo "builddir    = $(builddir)"
	@echo "VPATH       = $(VPATH)"
	@echo "DESTDIR     = $(DESTDIR)"
	@echo "======================================"
	@echo "CXX         = $(CXX)"
	@echo "CXXFLAGS    = $(CXXFLAGS)"
	@echo "DEFS        = $(DEFS)"
	@echo "LDFLAGS     = $(LDFLAGS)"
	@echo "LIBS        = $(LIBS)"
	@echo "======================================"

dist: $(distdir).tar.gz

$(distdir).tar.gz: FORCE $(distdir)
	tar chof - $(distdir) | gzip -9 -c > $(distdir).tar.gz
	rm -rf $(distdir)

$(distdir):
	mkdir -p $(distdir)/src
	cp $(srcdir)/configure $(distdir)
	cp $(srcdir)/configure.h.in $(distdir)
	cp $(srcdir)/install-sh $(distdir)
	cp $(srcdir)/Makefile.in $(distdir)
	cp -rf $(srcdir)/src/* $(distdir)/src/

FORCE:
	-rm $(distdir).tar.gz &> /dev/null
	-rm -rf $(distdir) &> /dev/null

distcheck: $(distdir).tar.gz
	gzip -cd $+ | tar xvf -
	$(MAKE) -C $(distdir) all clean
	rm -rf $(distdir)
	@echo "*** Package $(distdir).tar.gz ready for distribution."

clean:
	$(MAKE) -C $(SUBDIRS) $(@)

install:
	$(MAKE) -C $(SUBDIRS) $(@)

check:
	$(MAKE) -C test check

$(SUBDIRS):
	$(MAKE) -C $(@)

.PHONY: $(SUBDIRS) FORCE dist clean all distcheck install info
