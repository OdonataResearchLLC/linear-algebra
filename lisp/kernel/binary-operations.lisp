#|

 Linear Algebra Binary Operations Kernel

 Copyright (c) 2011-2012, Thomas M. Hermann
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

(in-package :linear-algebra-kernel)

;;; Interface

(defgeneric scaled-binary-op (op scalar1 scalar2)
  (:documentation
   "Compile and return a scaled binary operation."))

(defgeneric binary-operation
    (operation vector-or-matrix-1 vector-or-matrix-2 scalar1 scalar2)
  (:documentation
   "Perform the binary operation on the vectors or matrices."))

;;; Scaled binary operations

(defmethod scaled-binary-op (op
                             (scalar1 (eql nil))
                             (scalar2 (eql nil)))
  "Return the operation."
  op)

(defmethod scaled-binary-op ((op (eql #'+))
                             (scalar1 number)
                             (scalar2 (eql nil)))
  "Return the scaled operation."
  (lambda (n1 n2) (funcall op (* scalar1 n1) n2)))

(defmethod scaled-binary-op ((op (eql #'-))
                             (scalar1 number)
                             (scalar2 (eql nil)))
  "Return the scaled operation."
  (lambda (n1 n2) (funcall op (* scalar1 n1) n2)))

(defmethod scaled-binary-op ((op (eql #'+))
                             (scalar1 (eql nil))
                             (scalar2 number))
  "Return the scaled operation."
  (lambda (n1 n2) (funcall op n1 (* scalar2 n2))))

(defmethod scaled-binary-op ((op (eql #'-))
                             (scalar1 (eql nil))
                             (scalar2 number))
  "Return the scaled operation."
  (lambda (n1 n2) (funcall op n1 (* scalar2 n2))))

(defmethod scaled-binary-op ((op (eql #'+))
                             (scalar1 number)
                             (scalar2 number))
  "Return the scaled operation."
  (lambda (n1 n2) (funcall op (* scalar1 n1) (* scalar2 n2))))

(defmethod scaled-binary-op ((op (eql #'-))
                             (scalar1 number)
                             (scalar2 number))
  "Return the scaled operation."
  (lambda (n1 n2) (funcall op (* scalar1 n1) (* scalar2 n2))))

;;; Binary vector operations

(defun %vector-binary-operation (operation vector1 vector2)
  (map-into
   (make-array
    (length vector1)
    :element-type (common-array-element-type vector1 vector2))
   operation
   vector1 vector2))

(defun %destructive-vector-binary-operation
       (operation vector1 vector2)
  (map-into vector1 operation vector1 vector2))

(defmethod binary-operation ((operation (eql :add))
                             (vector1 vector)
                             (vector2 vector)
                             scalar1 scalar2)
  "Add the elements of the vectors and store the result in a new
vector."
  (%vector-binary-operation
   (scaled-binary-op #'+ scalar1 scalar2)
   vector1 vector2))

(defmethod binary-operation ((operation (eql :nadd))
                             (vector1 vector)
                             (vector2 vector)
                             scalar1 scalar2)
  (%destructive-vector-binary-operation
   (scaled-binary-op #'+ scalar1 scalar2)
   vector1 vector2))

(defmethod binary-operation ((operation (eql :subtract))
                             (vector1 vector)
                             (vector2 vector)
                             scalar1 scalar2)
  "Add the elements of the vectors and store the result in a new
vector."
  (%vector-binary-operation
   (scaled-binary-op #'- scalar1 scalar2)
   vector1 vector2))

(defmethod binary-operation ((operation (eql :nsubtract))
                             (vector1 vector)
                             (vector2 vector)
                             scalar1 scalar2)
  (%destructive-vector-binary-operation
   (scaled-binary-op #'- scalar1 scalar2)
   vector1 vector2))
