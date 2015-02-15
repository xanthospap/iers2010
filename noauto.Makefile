# --------------------------------------------------------------
# MAKEFILE FOR LIBIERS10++ LIBRARY
# TOP-LEVEL MAKEFILE
# @configure_input@
# --------------------------------------------------------------

export package = libiers10++
export version = 1.0
tarname = $(package)
distdir = $(tarname)-$(version)

export prefix      = /usr/local
export exec_prefix = $(prefix)
export datarootdir = $(prefix)/share
export includedir  = $(prefix)/include
export docdir      = $(datarootdir)/doc/$(package)
export libdir      = $(exec_prefix)/lib

SUBDIRS = src

all: $(SUBDIRS)

dist: $(distdir).tar.gz

$(distdir).tar.gz: FORCE $(distdir)
	tar chof - $(distdir) | gzip -9 -c > $(distdir).tar.gz
	rm -rf $(distdir)

$(distdir):
	mkdir -p $(distdir)/src
	cp Makefile $(distdir)
	cp -rf src/* $(distdir)/src/

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

$(SUBDIRS):
	$(MAKE) -C $(@)

.PHONY: $(SUBDIRS) FORCE dist clean all distcheck install
