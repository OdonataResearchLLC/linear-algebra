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

(in-package :linear-algebra)

;;; External Interface

(defgeneric sumsq (vector-or-matrix)
  (:documentation
   "Return the scaling parameter and the sum of the squares."))

(defgeneric sump (vector-or-matrix p)
  (:documentation
   "Return the scaling parameter and the sum of the P powers."))

(defgeneric norm (vector-or-matrix &optional measure)
  (:documentation
   "Return the norm according to measure."))

(defgeneric transpose (vector-or-matrix &optional conjugate)
  (:documentation
   "Transpose the vector or matrix."))

(defgeneric ntranspose (vector-or-matrix &optional conjugate)
  (:documentation
   "Destructively transpose the vector or matrix."))

(defgeneric permute (vector-or-matrix-1 vector-or-matrix-2)
  (:documentation
   "Permute the vector or matrix."))

(defgeneric scale (scalar vector-or-matrix)
  (:documentation
   "Scale each element by the scalar."))

(defgeneric nscale (scalar vector-or-matrix)
  (:documentation
   "Destructively scale each element by the scalar."))

(defgeneric add
    (vector-or-matrix-1 vector-or-matrix-2 &key scalar1 scalar2)
  (:documentation
   "Vector or matrix binary addition."))

(defgeneric nadd
    (vector-or-matrix-1 vector-or-matrix-2 &key scalar1 scalar2)
  (:documentation
   "Destructive vector or matrix addition."))

(defgeneric subtract
    (vector-or-matrix-1 vector-or-matrix-2 &key scalar1 scalar2)
  (:documentation
   "Vector or matrix binary subtraction."))

(defgeneric nsubtract
    (vector-or-matrix-1 vector-or-matrix-2 &key scalar1 scalar2)
  (:documentation
   "Destructive vector or matrix subtraction."))

(defgeneric product
    (vector-or-matrix-1 vector-or-matrix-2 &key scalar)
  (:documentation
   "Return the vector-vector, matrix-vector or matrix-matrix product."))

;;; Checks

(defmethod sump :before (v-or-m (p number))
  "The power must be positive."
  (declare (ignore scale sump))
  (unless (plusp p) (error "power(~A) must be positive." p)))
