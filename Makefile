# Commands
TANGLE = /usr/local/bin/notangle
WEAVE  = /usr/local/bin/noweave -delay
CPIF   = /usr/local/bin/cpif

# Noweb source documents
LINEAR-ALGEBRA-SRC      = noweb/linear-algebra.noweb
LINEAR-ALGEBRA-TEST-SRC = noweb/linear-algebra-test.noweb

# Common Lisp files
LINEAR-ALGEBRA += lisp/linear-algebra.asd
LINEAR-ALGEBRA += lisp/defpackage.lisp
LINEAR-ALGEBRA += lisp/auxiliary.lisp
LINEAR-ALGEBRA += lisp/fundamental-ops.lisp
LINEAR-ALGEBRA += lisp/vector.lisp
LINEAR-ALGEBRA += lisp/data-vector.lisp
LINEAR-ALGEBRA += lisp/matrix.lisp
LINEAR-ALGEBRA += lisp/identity-matrix.lisp
LINEAR-ALGEBRA += lisp/permutation-matrix.lisp
LINEAR-ALGEBRA += lisp/dense-matrix.lisp
LINEAR-ALGEBRA += lisp/square-matrix.lisp
LINEAR-ALGEBRA += lisp/hermitian-matrix.lisp
LINEAR-ALGEBRA += lisp/symmetric-matrix.lisp
LINEAR-ALGEBRA += lisp/triangular-matrix.lisp

# Unit test files
LINEAR-ALGEBRA-TEST += test/linear-algebra-test.asd
LINEAR-ALGEBRA-TEST += test/defpackage.lisp
LINEAR-ALGEBRA-TEST += test/vector.lisp
LINEAR-ALGEBRA-TEST += test/data-vector.lisp
LINEAR-ALGEBRA-TEST += test/matrix.lisp
LINEAR-ALGEBRA-TEST += test/identity-matrix.lisp
LINEAR-ALGEBRA-TEST += test/permutation-matrix.lisp
LINEAR-ALGEBRA-TEST += test/dense-matrix.lisp
LINEAR-ALGEBRA-TEST += test/square-matrix.lisp
LINEAR-ALGEBRA-TEST += test/hermitian-matrix.lisp
LINEAR-ALGEBRA-TEST += test/symmetric-matrix.lisp
LINEAR-ALGEBRA-TEST += test/triangular-matrix.lisp
LINEAR-ALGEBRA-TEST += test/auxiliary.lisp

# Documentation files
DOCUMENTATION-FILES += documentation/linear-algebra.tex
DOCUMENTATION-FILES += documentation/linear-algebra-test.tex

# Rules
default:	linear-algebra

all:	linear-algebra linear-algebra-test documentation

documentation:	$(DOCUMENTATION-FILES)

linear-algebra:	$(LINEAR-ALGEBRA)

linear-algebra-test:	$(LINEAR-ALGEBRA-TEST)

# Documentation rules

# Linear algebra rules
documentation/linear-algebra.tex:	$(LINEAR-ALGEBRA-SRC)
	$(WEAVE) $(LINEAR-ALGEBRA-SRC) | $(CPIF) $@

lisp/linear-algebra.asd:	$(LINEAR-ALGEBRA-SRC)
	$(TANGLE) -Rlinear-algebra.asd $(LINEAR-ALGEBRA-SRC) > $@

lisp/defpackage.lisp:	$(LINEAR-ALGEBRA-SRC)
	$(TANGLE) -Rdefpackage.lisp $(LINEAR-ALGEBRA-SRC) > $@

lisp/auxiliary.lisp:	$(LINEAR-ALGEBRA-SRC)
	$(TANGLE) -Rauxiliary.lisp $(LINEAR-ALGEBRA-SRC) > $@

lisp/fundamental-ops.lisp:	$(LINEAR-ALGEBRA-SRC)
	$(TANGLE) -Rfundamental-ops.lisp $(LINEAR-ALGEBRA-SRC) > $@

lisp/vector.lisp:	$(LINEAR-ALGEBRA-SRC)
	$(TANGLE) -Rvector.lisp $(LINEAR-ALGEBRA-SRC) > $@

lisp/data-vector.lisp:	$(LINEAR-ALGEBRA-SRC)
	$(TANGLE) -Rdata-vector.lisp $(LINEAR-ALGEBRA-SRC) > $@

lisp/matrix.lisp:	$(LINEAR-ALGEBRA-SRC)
	$(TANGLE) -Rmatrix.lisp $(LINEAR-ALGEBRA-SRC) > $@

lisp/identity-matrix.lisp:	$(LINEAR-ALGEBRA-SRC)
	$(TANGLE) -Ridentity-matrix.lisp $(LINEAR-ALGEBRA-SRC) > $@

lisp/permutation-matrix.lisp:	$(LINEAR-ALGEBRA-SRC)
	$(TANGLE) -Rpermutation-matrix.lisp $(LINEAR-ALGEBRA-SRC) > $@

lisp/dense-matrix.lisp:	$(LINEAR-ALGEBRA-SRC)
	$(TANGLE) -Rdense-matrix.lisp $(LINEAR-ALGEBRA-SRC) > $@

lisp/square-matrix.lisp:	$(LINEAR-ALGEBRA-SRC)
	$(TANGLE) -Rsquare-matrix.lisp $(LINEAR-ALGEBRA-SRC) > $@

lisp/hermitian-matrix.lisp:	$(LINEAR-ALGEBRA-SRC)
	$(TANGLE) -Rhermitian-matrix.lisp $(LINEAR-ALGEBRA-SRC) > $@

lisp/symmetric-matrix.lisp:	$(LINEAR-ALGEBRA-SRC)
	$(TANGLE) -Rsymmetric-matrix.lisp $(LINEAR-ALGEBRA-SRC) > $@

lisp/triangular-matrix.lisp:	$(LINEAR-ALGEBRA-SRC)
	$(TANGLE) -Rtriangular-matrix.lisp $(LINEAR-ALGEBRA-SRC) > $@

# Unit test rules
documentation/linear-algebra-test.tex:	$(LINEAR-ALGEBRA-TEST-SRC)
	$(WEAVE) $(LINEAR-ALGEBRA-TEST-SRC) | $(CPIF) $@

test/linear-algebra-test.asd:	$(LINEAR-ALGEBRA-TEST-SRC)
	$(TANGLE) -Rlinear-algebra-test.asd $(LINEAR-ALGEBRA-TEST-SRC) > $@

test/defpackage.lisp:	$(LINEAR-ALGEBRA-TEST-SRC)
	$(TANGLE) -Rdefpackage.lisp $(LINEAR-ALGEBRA-TEST-SRC) > $@

test/vector.lisp:	$(LINEAR-ALGEBRA-TEST-SRC)
	$(TANGLE) -Rvector.lisp $(LINEAR-ALGEBRA-TEST-SRC) > $@

test/data-vector.lisp:	$(LINEAR-ALGEBRA-TEST-SRC)
	$(TANGLE) -Rdata-vector.lisp $(LINEAR-ALGEBRA-TEST-SRC) > $@

test/matrix.lisp:	$(LINEAR-ALGEBRA-TEST-SRC)
	$(TANGLE) -Rmatrix.lisp $(LINEAR-ALGEBRA-TEST-SRC) > $@

test/identity-matrix.lisp:	$(LINEAR-ALGEBRA-TEST-SRC)
	$(TANGLE) -Ridentity-matrix.lisp $(LINEAR-ALGEBRA-TEST-SRC) > $@

test/permutation-matrix.lisp:	$(LINEAR-ALGEBRA-TEST-SRC)
	$(TANGLE) -Rpermutation-matrix.lisp $(LINEAR-ALGEBRA-TEST-SRC) > $@

test/dense-matrix.lisp:	$(LINEAR-ALGEBRA-TEST-SRC)
	$(TANGLE) -Rdense-matrix.lisp $(LINEAR-ALGEBRA-TEST-SRC) > $@

test/square-matrix.lisp:	$(LINEAR-ALGEBRA-TEST-SRC)
	$(TANGLE) -Rsquare-matrix.lisp $(LINEAR-ALGEBRA-TEST-SRC) > $@

test/hermitian-matrix.lisp:	$(LINEAR-ALGEBRA-TEST-SRC)
	$(TANGLE) -Rhermitian-matrix.lisp $(LINEAR-ALGEBRA-TEST-SRC) > $@

test/symmetric-matrix.lisp:	$(LINEAR-ALGEBRA-TEST-SRC)
	$(TANGLE) -Rsymmetric-matrix.lisp $(LINEAR-ALGEBRA-TEST-SRC) > $@

test/triangular-matrix.lisp:	$(LINEAR-ALGEBRA-TEST-SRC)
	$(TANGLE) -Rtriangular-matrix.lisp $(LINEAR-ALGEBRA-TEST-SRC) > $@

test/auxiliary.lisp:	$(LINEAR-ALGEBRA-TEST-SRC)
	$(TANGLE) -Rauxiliary.lisp $(LINEAR-ALGEBRA-TEST-SRC) > $@
