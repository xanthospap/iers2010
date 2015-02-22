
if HAVE_DOXYGEN
ddox: doxyfile.stamp doxyapidir
	@DOXYGEN@ Doxyfile
else
ddox:
	touch doxyfile.stamp
endif

doxyfile.stamp:
	echo Timestamp > doxyfile.stamp

CLEANFILES = doxyfile.stamp refguide.stamp

all-local: doxyfile.stamp refguide.stamp

doxyapidir:
	$(MKDIR_P) api_html 2>/dev/null

#endif

if HAVE_XSLTPROC
dxsl-html: refguide.stamp refguidedir xslstylesheet
	cd xml && @XSLTPROC@ --output refguide_html $(XSLSS_CHUNK_HTML) libiers10-main.xml
	-mv xml/*.html refguide_html/ 2>/dev/null
	-mv xml/html/*.html refguide_html/ 2>/dev/null
	-rm -rf xml/html 2>/dev/null
if HAVE_FOP
dxsl-pdf: refguide.stamp refguidedir xslstylesheet
	cd xml && @XSLTPROC@ --stringparam  paper.type  A4 $(XSLSS_FO) libiers10-main.xml > libiers.fo
	@FOP@ -q -fo xml/libiers.fo -pdf refguide.pdf
	-rm xml/libiers.fo
else
dxsl-pdf:
	touch refguide.stamp
endif
else
dxsl-html: 
	touch refguide.stamp
endif

refguide.stamp:
	echo Timestamp > refguide.stamp

all: ddox dxsl-html dxsl-pdf

xslstylesheet:
	$(eval XSLSS_CHUNK_HTML = $(shell find /usr/share/ -name chunk.xsl 2>/dev/null | $(GREP) html | head -1) )
	$(eval XSLSS_FO = $(shell find /usr/share/ -name docbook.xsl 2>/dev/null | $(GREP) fo/ |head -1) )
	@echo "-----------------------------------------------------"
	@echo "XSL STYLESHEET FOR (CHUNK) HTML SET TO $(XSLSS_CHUNK_HTML)"
	@echo "XSL STYLESHEET FOR FO SET TO $(XSLSS_FO)"
	@echo "-----------------------------------------------------"

refguidedir:
	$(MKDIR_P) refguide_html 2>/dev/null