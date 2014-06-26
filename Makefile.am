# this is the file: Makefile.am
ACLOCAL_AMFLAGS = -I m4 

# include docygen am file
include $(top_srcdir)/doxygen.am
dist_doc_DATA= \
  doc/html/* \
  doc/*.pdf \
  doc/*.ps


# the subdirectories of the project to go into
SUBDIRS = lib src tests
