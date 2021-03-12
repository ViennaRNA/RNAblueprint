# this is the file: Makefile.am
ACLOCAL_AMFLAGS = -I m4 
AUTOMAKE_OPTIONS = foreign

AM_DISTCHECK_CONFIGURE_FLAGS = \
    --enable-doxygen-doc \
    --enable-doxygen-xml \
    --enable-doxygen-html \
     --disable-silent-rules


# include docygen am file
include $(top_srcdir)/doxygen.am

# the subdirectories of the project to go into
SUBDIRS = lib src interfaces tests

if DX_COND_man
# post-process man pages and install them
dist_man3_MANS = @DX_DOCDIR@/man/man3/RNAblueprint.h.3

$(dist_man3_MANS): doxygen-doc

endif

docdata=

if DX_COND_html
htmldocdir = $(docdir)/html
endif

if DX_COND_xml
xmldocdir = $(docdir)/xml
endif

if DX_COND_pdf
docdata += @DX_DOCDIR@/RNAblueprint.pdf
endif

if DX_COND_html
dist_htmldoc_DATA = @DX_DOCDIR@/html/*

$(dist_htmldoc_DATA): doxygen-doc

endif
if DX_COND_xml
docdata += @DX_DOCDIR@/RNAblueprint.tag 
dist_xmldoc_DATA = @DX_DOCDIR@/xml/*

$(dist_xmldoc_DATA): doxygen-doc

endif

dist_doc_DATA = $(docdata) 

$(dist_doc_DATA): doxygen-doc


EXTRA_DIST = \
    data

MOSTLYCLEANFILES = $(DX_CLEANFILES)
