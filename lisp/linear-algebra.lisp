#|

 Linear Algebra in Common Lisp

 Copyright (c) 2011-2012, Odonata Research LLC
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

(in-package :cl-user)

(defpackage :linear-algebra
  (:use :common-lisp :floating-point)
  (:use :linear-algebra-kernel)
  ;; Fundamental operations
  (:export :sumsq :sump
           :norm
           :transpose :ntranspose
           :permute
           :scale :nscale
           :add :nadd
           :subtract :nsubtract
           :product)
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
           :dovector
           :apply-rotation :napply-rotation)
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
           :submatrix
           :replace-matrix
           :matrix-validated-range)
  ;; Identity matrix
  (:export :identity-matrix
           :identity-matrix-p)
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
           :symmetric-matrix-p)
  ;; Triangular matrix
  (:export :upper-triangular-matrix
           :upper-triangular-matrix-p
           :lower-triangular-matrix
           :lower-triangular-matrix-p))
