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

(defgeneric norm (vector-or-matrix &optional measure)
  (:documentation
   "Return the norm according to measure."))

(defgeneric min-vector (vector)
  (:documentation
   "Return the minimum absolute value"))

(defgeneric distance (vector1 vector2 &optional measure)
  (:documentation
   "Return the distance between vector1 and vector2 using the measure"))

(defgeneric transpose (vector-or-matrix)
  (:documentation
   "Transpose the vector or matrix."))

(defgeneric ntranspose (vector-or-matrix)
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
    (vector-or-matrix-1 vector-or-matrix-2 &optional scalar)
  (:documentation
   "Return the vector-vector, matrix-vector or matrix-matrix product."))

(defgeneric solve (matrix vector)
  (:documentation
   "Return the solution to the system of equations."))

(defgeneric nsolve (matrix vector)
  (:documentation
   "Return the solution to the system of equations in-place."))

(defgeneric invert (matrix)
  (:documentation
   "Return the invert of the matrix."))

(defgeneric ninvert (matrix)
  (:documentation
   "Return the invert of the matrix with in-place decomposition."))

(defgeneric add-diagonal (scalar-or-vector matrix)
  (:documentation
   "Add a scalar or vector to every diagonal of the matrix"))

