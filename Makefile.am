#
# gettext support requires that subdir intl be first
#
# Eric Bowden (May 31, '07) Removing internationalization for the time being
#SUBDIRS = m4 po lui5 src util doc @GTK_DIRS@
# Eric Bowden (May 31, '07) Removing awkward GTK_DIRS conditional compilation for
# something more elegant.  (This also solves the distribution problem.)
#SUBDIRS = m4 lui5 src util doc @GTK_DIRS@
SUBDIRS = m4 po lui5 src util doc gtk
EXTRA_DIST = config.rpath  PORTING COPYRIGHT contrib convert scripts \
             userfuncs pixmaps vis5d+.pc.in

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = vis5d+.pc

#
# This target creates lib vis5dapi which contains only those functions 
# nessasary to access vis5d via the api 
#
API: 
	cd src; $(MAKE) API

# add a target to completely get rid of all the files produced by
# autoheader, automake, autoconf, and aclocal:
super-clean: maintainer-clean
	rm -f configure */Makefile.in Makefile.in INSTALL aclocal.m4 config.guess config.h.in config.sub install-sh ltconfig ltmain.sh missing mkinstalldirs stamp-h.in


ACLOCAL_AMFLAGS = -I m4
