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
dxsl: refguide.stamp refguidedir
	cd xml && @XSLTPROC@ /usr/share/sgml/docbook/xsl-ns-stylesheets-1.78.1/html/chunk.xsl libiers10-main.xml
	mv xml/*.html refguide_html/
else
dxsl: 
	touch refguide.stamp
endif

refguide.stamp:
	echo Timestamp > refguide.stamp

all: ddox dxsl

#CLEANFILES = refguide.stamp

#all-local: refguide.stamp

refguidedir:
	$(MKDIR_P) refguide_html 2>/dev/null
#endif
