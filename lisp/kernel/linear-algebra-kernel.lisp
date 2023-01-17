#|

  Linear Algebra Kernel Functions

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

(defpackage :linear-algebra-kernel
  (:use :common-lisp :floating-point)
  ;; Utility functions
  (:export :copy-array
           :common-class-of
           :common-array-element-type
           :specific-array-element-type
           :complex-equal
           :number-equal)
  ;; Permutation
  (:export :right-permute
           :left-permute)
  ;; Unary operations
  (:export :sumsq2 :sumsq3
           :sumsq :sump
           :sumsq-row
           :sumsq-column
           :norm-vector
           :norm-array)
  ;; Binary operations
  (:export :compatible-dimensions-p
           :scaled-binary-op
           :add-vector :nadd-vector
           :subtract-vector :nsubtract-vector
           :add-array :nadd-array
           :subtract-array :nsubtract-array
           :element-multiply-vector
           :element-divide-vector
           :element-greater-vector
           :inner-product-vector
           :product-vector-array
           :product-array-vector
           :product-array-array)
  ;; Rotations
  (:export :givens-rotation
           :jacobi-rotation
           :householder-reflection)
  ;; Gauss algorithm
  (:export :gauss-solver
           :gauss-invert)
  ;; Cholesky
  (:export :symmetric-cholesky-decomposition
           :hermitian-cholesky-decomposition
           :root-free-symmetric-cholesky-decomposition
           :root-free-hermitian-cholesky-decomposition
           :symmetric-cholesky-solver
           :hermitian-cholesky-solver
           :symmetric-cholesky-invert
           :hermitian-cholesky-invert)
  ;; Conjugate gradient method
  (:export :conjugate-gradient-solver)
  ;; Tridiagonal
  (:export :tridiagonal-solver))
