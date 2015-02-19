if HAVE_DOXYGEN

all: doxyfile.stamp doxyapidir
	@DOXYGEN@ Doxyfile

doxyfile.stamp:
	echo Timestamp > doxyfile.stamp

CLEANFILES = doxyfile.stamp

all-local: doxyfile.stamp

doxyapidir:
	$(MKDIR_P) api_html 2>/dev/null

endif
