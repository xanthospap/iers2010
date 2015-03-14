EXTRA_DIST = xml/libiers10-chap-01.xml \
			 xml/libiers10-chap-02.xml \
			 xml/libiers10-chap-03.xml \
			 xml/libiers10-chap-04.xml \
			 xml/libiers10-chap-05.xml \
			 xml/libiers10-main.xml

validate-dtd:
	xmllint --valid --noout xml/libiers10-main.xml

api_html:
	-cp @top_builddir@/README @top_builddir@/src/README.md
	$(DOXYGEN) Doxyfile-html
	-rm @top_builddir@/src/README.md

clean_api_html:
	-rm -rf api_html/

api_man:
	$(DOXYGEN) Doxyfile-man
	-rm  api_man/*.cpp*

clean_api_man:
	-rm -rf api_man/

ref_guide_html: xslstylesheet
	-$(MKDIR_P) refguide_html
	cd xml && $(XSLTPROC) --output refguide_html $(XSLSS_CHUNK_HTML) libiers10-main.xml 2>dummy.log
	-@mv xml/*.html refguide_html/ 2>dummy.log
	-@mv xml/html/*.html refguide_html/ 2>dummy.log
	-@rm -rf xml/html 2>dummy.log

clean_ref_guide_html:
	-rm -rf refguide_html

ref_guide_pdf: xslstylesheet
	cd xml && $(XSLTPROC) --stringparam  paper.type  A4 $(XSLSS_FO) libiers10-main.xml 1> libiers.fo 2>dummy.log
	$(FOP) -q -fo xml/libiers.fo -pdf refguide.pdf 2>dummy.log
	-rm xml/libiers.fo

clean_ref_guide_pdf:
	-rm -f refguide.pdf

xslstylesheet:
	$(eval XSLSS_CHUNK_HTML = $(shell find /usr/share/ -name chunk.xsl 2>/dev/null | $(GREP) html | head -1) )
	$(eval XSLSS_FO = $(shell find /usr/share/ -name docbook.xsl 2>/dev/null | $(GREP) fo/ |head -1) )
	@echo "-----------------------------------------------------"
	@echo "XSL STYLESHEET FOR (CHUNK) HTML SET TO $(XSLSS_CHUNK_HTML)"
	@echo "XSL STYLESHEET FOR FO SET TO $(XSLSS_FO)"
	@echo "-----------------------------------------------------"

if HAVE_DOXYGEN
API_HTML: api_html
API_MAN: api_man
else
API_HTML:
API_MAN:
endif

if HAVE_XSLTPROC
REF_GUIDE_HTML: ref_guide_html xslstylesheet
REF_GUIDE_PDF: ref_guide_pdf xslstylesheet
else
REF_GUIDE_HTML:
REF_GUIDE_PDF:
endif

rm_log:
	-rm dummy.log xml/dummy.log

all:

doc: REF_GUIDE_HTML REF_GUIDE_PDF API_HTML API_MAN rm_log

clean: clean_api_html clean_api_man clean_ref_guide_html clean_ref_guide_pdf rm_log

.PHONY: REF_GUIDE_HTML REF_GUIDE_PDF API_HTML API_MAN all clean rm_log

