# this is the file: Makefile.am
ACLOCAL_AMFLAGS = -I m4 
AUTOMAKE_OPTIONS = foreign

# include docygen am file
include $(top_srcdir)/doxygen.am

# the subdirectories of the project to go into
SUBDIRS = lib src interfaces tests

doc:
	$(MAKE) doxygen-doc

if DX_COND_man
# post-process man pages and install them
dist_man3_MANS = @DX_DOCDIR@/man/man3/RNAdesign.h.3

$(dist_man3_MANS): doxygen-doc

endif

docdata=
htmldocdir = $(docdir)/html
xmldocdir = $(docdir)/xml

if DX_COND_pdf
docdata += @DX_DOCDIR@/RNAdesign.pdf
endif
if DX_COND_html
dist_htmldoc_DATA = @DX_DOCDIR@/html/*

$(dist_htmldoc_DATA): doxygen-doc

endif
if DX_COND_xml
docdata += @DX_DOCDIR@/RNAdesign.tag 
dist_xmldoc_DATA = @DX_DOCDIR@/xml/*

$(dist_htmldoc_DATA): doxygen-doc

endif

dist_doc_DATA = $(docdata) 

$(dist_doc_DATA): doxygen-doc

MOSTLYCLEANFILES = $(DX_CLEANFILES)
