#|

  Linear Algebra in Common Lisp

  Copyright (c) 2011-2014, Odonata Research LLC

  Permission is hereby granted, free  of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction,  including without limitation the rights
  to use, copy, modify,  merge,  publish,  distribute,  sublicense, and/or sell
  copies of the  Software,  and  to  permit  persons  to  whom  the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and  this  permission  notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED  "AS IS",  WITHOUT  WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT  NOT  LIMITED  TO  THE  WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE  AND  NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT  HOLDERS  BE  LIABLE  FOR  ANY  CLAIM,  DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

|#

(in-package :cl-user)

(defpackage :linear-algebra
  (:use :common-lisp)
  (:use :floating-point :linear-algebra-kernel)
  ;; Fundamental operations
  (:export :norm
	   :min-vector
           :transpose :ntranspose
           :permute
           :scale :nscale
           :add :nadd
           :subtract :nsubtract
           :product
           :solve :nsolve
           :invert :ninvert)
  ;; Vector exports
  (:export :initialize-vector
           :make-vector
           :vector-in-bounds-p
           :vector-element-type
           :vector-length
           :vref
           :copy-vector
           :subvector
           :replace-vector
           :map-vector
           :map-into-vector
	   :reduce-vector
           :dovector
           :vec-equal
   	   :elem-divide
           :elem-multiply
           :vec-every
           :outer-product
           :elem-greater
   	   :distance
   :apply-rotation :napply-rotation
   )
  ;; Matrix interface
  (:export :matrix-object
           :initialize-matrix
           :make-matrix
           :matrixp
           :matrix-in-bounds-p
           :matrix-element-type
           :matrix-dimensions
           :matrix-row-dimension
           :matrix-column-dimension
           :mref
           :copy-matrix
           :mat-equal
           :trace
           :submatrix
           :replace-matrix
           :matrix-validated-range
	   :add-diagonal
	   )
  ;; Identity matrix
  (:export
   :identity-matrix
   :identity-matrix-p
   )
  ;; Permutation matrix
  (:export :permutation-matrix
           :permutation-matrix-p)
  ;; Data vector exports
  (:export :data-vector
           :row-vector
           :row-vector-p
           :column-vector
           :column-vector-p)
  ;; Dense matrix
  (:export :dense-matrix
           :dense-matrix-p)
  ;; Square matrix
  (:export :square-matrix
           :square-matrix-p)
  ;; Hermitian matrix
  (:export :hermitian-matrix
           :hermitian-matrix-p)
  ;; Symmetric matrix
  (:export :symmetric-matrix
           :symmetric-matrix-p))
