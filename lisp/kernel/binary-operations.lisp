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
    (operation
     vector-or-matrix-1
     vector-or-matrix-2
     scalar1 scalar2-or-flag)
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

(defmethod scaled-binary-op ((op (eql #'*))
                             (scalar (eql nil))
                             (conjugate (eql t)))
  (lambda (n1 n2) (* (conjugate n1) n2)))

;;; Binary vector operations

(defun %vector<-vector1-op-vector2 (operation vector1 vector2)
  (map-into
   (make-array
    (length vector1)
    :element-type (common-array-element-type vector1 vector2))
   operation
   vector1 vector2))

(defun %vector1<-vector1-op-vector2 (operation vector1 vector2)
  (map-into vector1 operation vector1 vector2))

(defun add-vector (vector1 vector2 scalar1 scalar2)
  (%vector<-vector1-op-vector2
   (scaled-binary-op #'+ scalar1 scalar2)
   vector1 vector2))

(defun nadd-vector (vector1 vector2 scalar1 scalar2)
  (%vector1<-vector1-op-vector2
   (scaled-binary-op #'+ scalar1 scalar2)
   vector1 vector2))

(defun subtract-vector (vector1 vector2 scalar1 scalar2)
  (%vector<-vector1-op-vector2
   (scaled-binary-op #'- scalar1 scalar2)
   vector1 vector2))

(defun nsubtract-vector (vector1 vector2 scalar1 scalar2)
  (%vector1<-vector1-op-vector2
   (scaled-binary-op #'- scalar1 scalar2)
   vector1 vector2))

(defun inner-product-vector (vector1 vector2 scalar conjugate)
  (loop with op = (scaled-binary-op #'* nil conjugate)
        for element1 across vector1
        and element2 across vector2
        sum (funcall op element1 element2) into result
        finally
        (return (if scalar (* scalar result) result))))

;;; Binary array/vector operations

(defmethod binary-operation ((operation (eql :product))
                             (vector vector) (array array)
                             scalar conjugate)
  (let* ((m-rows (array-dimension array 0))
         (n-columns (array-dimension array 1))
         (zero (coerce 0 (array-element-type vector)))
         (element)
         (result
          (make-array
           n-columns
           :element-type (array-element-type vector))))
    (dotimes (i1 n-columns result)
      (setf element zero)
      (dotimes (i0 m-rows)
        (incf element (* (aref vector i0) (aref array i0 i1))))
      (if scalar
          (setf (aref result i1) (* scalar element))
          (setf (aref result i1) element)))))

;;; Binary array operations

(defun %array<-array1-op-array2 (operation array1 array2)
  (let ((m-rows (array-dimension array1 0))
        (n-columns (array-dimension array1 1))
        (result
         (make-array
          (array-dimensions array1)
          :element-type
          (common-array-element-type array1 array2))))
    (dotimes (i0 m-rows result)
      (dotimes (i1 n-columns)
        (setf
         (aref result i0 i1)
         (funcall operation
                  (aref array1 i0 i1)
                  (aref array2 i0 i1)))))))

(defun %array1<-array1-op-array2 (operation array1 array2)
  (let ((m-rows (array-dimension array1 0))
        (n-columns (array-dimension array1 1)))
    (dotimes (i0 m-rows array1)
      (dotimes (i1 n-columns)
        (setf
         (aref array1 i0 i1)
         (funcall operation
                  (aref array1 i0 i1)
                  (aref array2 i0 i1)))))))

(defun add-array (array1 array2 scalar1 scalar2)
  (%array<-array1-op-array2
   (scaled-binary-op #'+ scalar1 scalar2)
   array1 array2))

(defun nadd-array (array1 array2 scalar1 scalar2)
  (%array1<-array1-op-array2
   (scaled-binary-op #'+ scalar1 scalar2)
   array1 array2))

(defun subtract-array (array1 array2 scalar1 scalar2)
  (%array<-array1-op-array2
   (scaled-binary-op #'- scalar1 scalar2)
   array1 array2))

(defun nsubtract-array (array1 array2 scalar1 scalar2)
  (%array1<-array1-op-array2
   (scaled-binary-op #'- scalar1 scalar2)
   array1 array2))
