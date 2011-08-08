## Process this file with automake to produce Makefile.in

SGML_SRC =  $(top_srcdir)/doc/chap1.sgml \
			$(top_srcdir)/doc/chap2.sgml \
			$(top_srcdir)/doc/chap3.sgml \
			$(top_srcdir)/doc/chap4.sgml \
			$(top_srcdir)/doc/chap5.sgml \
			$(top_srcdir)/doc/chap6.sgml \
			$(top_srcdir)/doc/chap7.sgml \
			$(top_srcdir)/doc/chap8.sgml \
			$(top_srcdir)/doc/chap9.sgml \
			$(top_srcdir)/doc/chap1.sgml \
			$(top_srcdir)/doc/chap1.sgml \
	   		$(top_srcdir)/doc/acknowledgments.sgml \
			$(top_srcdir)/doc/license.sgml \
           	$(top_srcdir)/doc/vis5d.sgml

EXTRA_DIST = vis5d.pdf html $(SGML_SRC)

# Yes that final suffix is od -- it has to do with a bug in the db2* programs
DB2DVI = $(top_srcdir)/missing --run db2dvi --output $(abs_builddir)
DB2HTML = $(top_srcdir)/missing --run db2html --output $(abs_builddir)/html
DB2PDF = $(top_srcdir)/missing --run db2pdf --output $(abs_builddir)
DB2PS = $(top_srcdir)/missing --run db2ps --output $(abs_builddir)
DB2RTF = $(top_srcdir)/missing --run db2rtf --output $(abs_builddir)

ps: vis5d.ps

dvi: vis5d.dvi

pdf: vis5d.pdf

rtf: vis5d.rtf

html: $(SGML_SRC)
	rm -rf html
	$(DB2HTML) $(top_srcdir)/doc/vis5d.sgml
#mv html
	cd html && ln -s -T index.html book1.html

vis5d.dvi: $(SGML_SRC)
	$(DB2DVI) $(top_srcdir)/doc/vis5d.sgml

vis5d.ps: $(SGML_SRC)
	$(DB2PS) $(top_srcdir)/doc/vis5d.sgml

vis5d.pdf: $(SGML_SRC)
	$(DB2PDF) $(top_srcdir)/doc/vis5d.sgml
#mv  $(top_srcdir)/doc/vis5d.pdf $(builddir)

vis5d.rtf: $(SGML_SRC)
	$(DB2RTF) $(top_srcdir)/doc/vis5d.sgml

clean-local:
	rm -f *.aux *.log

maintainer-clean-local:
	rm -f *.dvi *.ps *.pdf *.rtf *.tex
	rm -rf html
