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
dxsl: refguide.stamp refguidedir xslstylesheet
	cd xml && @XSLTPROC@ --output refguide_html $(XSLSS) libiers10-main.xml
	-mv xml/*.html refguide_html/ 2>/dev/null
	-mv xml/html/*.html refguide_html/ 2>/dev/null
	-rm -rf xml/html 2>/dev/null
else
dxsl: 
	touch refguide.stamp
endif

refguide.stamp:
	echo Timestamp > refguide.stamp

all: ddox dxsl

xslstylesheet:
	$(eval XSLSS = $(shell find /usr/share/ -name chunk.xsl 2>/dev/null | grep xml/docbook/stylesheet | grep html | head -1) )
	@echo "-----------------------------------------------------"
	@echo "XSL STYLESHEET SET TO $(XSLSS)"
	@echo "-----------------------------------------------------"

refguidedir:
	$(MKDIR_P) refguide_html 2>/dev/null
#endif

## xsltproc --output admin_guide.fo /usr/share/xml/docbook/stylesheet/suse2013/fo/fo.xsl libiers10-main.xml
##
##