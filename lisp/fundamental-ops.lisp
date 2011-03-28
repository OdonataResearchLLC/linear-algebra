#|

 Linear Algebra in Common Lisp

 Copyright (c) 2011, Thomas M. Hermann
 All rights reserved.

 Redistribution and  use  in  source  and  binary  forms, with or without
 modification, are permitted  provided  that the following conditions are
 met:

   o  Redistributions of  source  code  must  retain  the above copyright
      notice, this list of conditions and the following disclaimer.
   o  Redistributions in binary  form  must reproduce the above copyright
      notice, this list of  conditions  and  the  following disclaimer in
      the  documentation  and/or   other   materials  provided  with  the
      distribution.
   o  The names of the contributors may not be used to endorse or promote
      products derived from this software without  specific prior written
      permission.

 THIS SOFTWARE IS  PROVIDED  BY  THE  COPYRIGHT  HOLDERS AND CONTRIBUTORS
 "AS IS"  AND  ANY  EXPRESS  OR  IMPLIED  WARRANTIES, INCLUDING,  BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES  OF MERCHANTABILITY AND FITNESS FOR A
 PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT  NOT LIMITED TO,
 PROCUREMENT OF  SUBSTITUTE  GOODS  OR  SERVICES;  LOSS  OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER  CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER  IN  CONTRACT,  STRICT  LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR  OTHERWISE)  ARISING  IN  ANY  WAY  OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

|#

(in-package :linear-algebra)

(defgeneric sumsq (vector-or-matrix &key scale sumsq)
  (:documentation
   "Return the scaling parameter and the sum of the squares."))

(defgeneric sump (vector-or-matrix p &key scale sump)
  (:documentation
   "Return the scaling parameter and the sum of the P powers."))

(defgeneric norm (vector-or-matrix &key measure)
  (:documentation
   "Return the norm according to measure."))

(defgeneric transpose (vector-or-matrix &key conjugate)
  (:documentation
   "Transpose the vector or matrix."))

(defgeneric ntranspose (vector-or-matrix &key conjugate)
  (:documentation
   "Destructively transpose the vector or matrix."))

(defgeneric permute (vector-or-matrix vector-or-matrix)
  (:documentation
   "Permute the vector or matrix."))

(defgeneric npermute (vector-or-matrix vector-or-matrix)
  (:documentation
   "Destructively permute the vector or matrix."))

(defgeneric scale (scalar vector-or-matrix)
  (:documentation
   "Scale each element by the scalar."))

(defgeneric nscale (scalar vector-or-matrix)
  (:documentation
   "Destructively scale each element by the scalar."))

(defgeneric add (vector-or-matrix-1 vector-or-matrix-2
                 &key scalar1 scalar2)
  (:documentation
   "Vector or matrix binary addition."))

(defgeneric nadd (vector-or-matrix-1 vector-or-matrix-2
                  &key scalar1 scalar2)
  (:documentation
   "Destructive vector or matrix addition."))

(defgeneric subtract (vector-or-matrix-1 vector-or-matrix-2
                      &key scalar1 scalar2)
  (:documentation
   "Vector or matrix binary subtraction."))

(defgeneric nsubtract (vector-or-matrix-1 vector-or-matrix2
                       &key scalar1 scalar2)
  (:documentation
   "Destructive vector or matrix subtraction."))

(defgeneric product (vector-or-matrix-1 vector-or-matrix-2 &key scalar)
  (:documentation
   "Return the vector-vector, matrix-vector or matrix-matrix product."))

