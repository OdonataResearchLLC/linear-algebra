#
   BIBTEX = /usr/bin/bibtex
MAKEINDEX = /usr/bin/makeindex
 PDFLATEX = /usr/bin/pdflatex -shell-escape
   ASPELL = /usr/bin/aspell --mode=tex check

default:	linear-algebra.pdf

all:	linear-algebra.pdf linear-algebra-test.pdf

# linear-algebra
linear-algebra.pdf:	linear-algebra.tex bibliography.bib
	$(PDFLATEX)  linear-algebra && \
	$(BIBTEX)    linear-algebra && \
	$(MAKEINDEX) -s nomencl.ist -o linear-algebra.gls linear-algebra.glo && \
	$(PDFLATEX)  linear-algebra && \
	$(PDFLATEX)  linear-algebra && \
	$(PDFLATEX)  linear-algebra

# linear-algebra-test
linear-algebra-test.pdf:	linear-algebra-test.tex bibliography.bib
	$(PDFLATEX)  linear-algebra-test && \
	$(BIBTEX)    linear-algebra-test && \
	$(MAKEINDEX) -s nomencl.ist -o linear-algebra-test.gls linear-algebra-test.glo && \
	$(PDFLATEX)  linear-algebra-test && \
	$(PDFLATEX)  linear-algebra-test && \
	$(PDFLATEX)  linear-algebra-test

#
.PHONY:	clean clobber spell-check spell-check-linear-algebra-test

spell-check:
	$(ASPELL) linear-algebra.tex

spell-check-linear-algebra-test:
	$(ASPELL) linear-algebra-test.tex

clean:
	-rm *.aux
	-rm *.glo
	-rm *.gls
	-rm *.ilg
	-rm *.log
	-rm *.bbl
	-rm *.blg

clobber:
	-rm linear-algebra.pdf
	-rm linear-algebra-test.pdf
