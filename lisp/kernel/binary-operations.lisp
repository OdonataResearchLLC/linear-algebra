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

(defgeneric compatible-dimensions-p
    (operation vector-or-matrix-1 vector-or-matrix-2)
  (:documentation
   "Return true if the vector and matrix dimensions are compatible for
the operation."))

(defgeneric scaled-binary-op (op scalar1 scalar2)
  (:documentation
   "Compile and return a scaled binary operation."))

;;; Scaled binary operations

(defmethod scaled-binary-op
    (op (scalar1 (eql nil)) (scalar2 (eql nil)))
  "Return the operation."
  op)

(defmethod scaled-binary-op
    ((op (eql #'+)) (scalar1 number) (scalar2 (eql nil)))
  "Return the scaled operation."
  (lambda (n1 n2) (+ (* scalar1 n1) n2)))

(defmethod scaled-binary-op
    ((op (eql #'-)) (scalar1 number) (scalar2 (eql nil)))
  "Return the scaled operation."
  (lambda (n1 n2) (- (* scalar1 n1) n2)))

(defmethod scaled-binary-op
    ((op (eql #'+)) (scalar1 (eql nil)) (scalar2 number))
  "Return the scaled operation."
  (lambda (n1 n2) (+ n1 (* scalar2 n2))))

(defmethod scaled-binary-op
    ((op (eql #'-)) (scalar1 (eql nil)) (scalar2 number))
  "Return the scaled operation."
  (lambda (n1 n2) (- n1 (* scalar2 n2))))

(defmethod scaled-binary-op
    ((op (eql #'+)) (scalar1 number) (scalar2 number))
  "Return the scaled operation."
  (lambda (n1 n2) (+ (* scalar1 n1) (* scalar2 n2))))

(defmethod scaled-binary-op
    ((op (eql #'-)) (scalar1 number) (scalar2 number))
  "Return the scaled operation."
  (lambda (n1 n2) (- (* scalar1 n1) (* scalar2 n2))))

(defmethod scaled-binary-op
    ((op (eql #'*)) (scalar (eql nil)) (conjugate (eql t)))
  (lambda (n1 n2) (* (conjugate n1) n2)))

;;; Binary vector operations

(defun %vector<-vector1-op-vector2 (operation vector1 vector2)
  "Store the result of the binary operation in a new vector."
  (let ((result
         (make-array
          (length vector1)
          :element-type
          (common-array-element-type vector1 vector2))))
    (dotimes (index (length vector1) result)
      (setf
       (aref result index)
       (funcall operation
                (aref vector1 index)
                (aref vector2 index))))))

(defun %vector1<-vector1-op-vector2 (operation vector1 vector2)
  "Store the result of the binary operation in vector1."
  (dotimes (index (length vector1) vector1)
    (setf
     (aref vector1 index)
     (funcall operation
              (aref vector1 index)
              (aref vector2 index)))))

(defun add-vector (vector1 vector2 scalar1 scalar2)
  "Vector binary addition."
  (%vector<-vector1-op-vector2
   (scaled-binary-op #'+ scalar1 scalar2)
   vector1 vector2))

(defun nadd-vector (vector1 vector2 scalar1 scalar2)
  "Destructive vector binary addition."
  (%vector1<-vector1-op-vector2
   (scaled-binary-op #'+ scalar1 scalar2)
   vector1 vector2))

(defun subtract-vector (vector1 vector2 scalar1 scalar2)
  "Vector binary subtraction."
  (%vector<-vector1-op-vector2
   (scaled-binary-op #'- scalar1 scalar2)
   vector1 vector2))

(defun nsubtract-vector (vector1 vector2 scalar1 scalar2)
  "Destructive vector binary subtraction."
  (%vector1<-vector1-op-vector2
   (scaled-binary-op #'- scalar1 scalar2)
   vector1 vector2))

(defun inner-product-vector (vector1 vector2 scalar conjugate)
  "Return the vector inner product."
  (loop with op = (scaled-binary-op #'* nil conjugate)
        for element1 across vector1
        and element2 across vector2
        sum (funcall op element1 element2) into result
        finally
        (return (if scalar (* scalar result) result))))

;;; Binary array/vector operations

(defmethod compatible-dimensions-p
           ((operation (eql :product)) (vector vector) (array array))
  "Return true if the array dimensions are compatible for product."
  (and
   (= 2 (array-rank array))
   (= (length vector) (array-dimension array 0))))

(defmethod compatible-dimensions-p
           ((operation (eql :product)) (array array) (vector vector))
  "Return true if the array dimensions are compatible for product."
  (and
   (= 2 (array-rank array))
   (= (array-dimension array 1) (length vector))))

(defun %product-vector-array (vector array)
  "Return the result of the array premultiplied by the vector."
  (let* ((m-rows (array-dimension array 0))
         (n-columns (array-dimension array 1))
         (zero (coerce 0 (array-element-type vector)))
         (element)
         (result
          (make-array
           n-columns
           :element-type (array-element-type vector))))
    (dotimes (column n-columns result)
      (setq element zero)
      (dotimes (row m-rows)
        (setq
         element
         (+ element (* (aref vector row) (aref array row column)))))
      ;; Store the result
      (setf (aref result column) element))))

(defun %scaled-product-vector-array (vector array scalar)
  "Return the result of the array premultiplied by the vector and
scaled."
  (let* ((m-rows (array-dimension array 0))
         (n-columns (array-dimension array 1))
         (zero (coerce 0 (array-element-type vector)))
         (element)
         (result
          (make-array
           n-columns
           :element-type (array-element-type vector))))
    (dotimes (column n-columns result)
      (setq element zero)
      (dotimes (row m-rows)
        (setq
         element
         (+ element (* (aref vector row) (aref array row column)))))
      ;; Store the result
      (setf (aref result column) (* scalar element)))))

(defun product-vector-array (vector array scalar)
  "Return the result of the array premultiplied by the vector and
scaled."
  (if scalar
      (%scaled-product-vector-array vector array scalar)
      (%product-vector-array vector array)))

(defun %product-array-vector (array vector)
  "Return the result of the array postmultiplied by the vector."
  (let* ((m-rows (array-dimension array 0))
         (n-columns (array-dimension array 1))
         (zero (coerce 0 (array-element-type vector)))
         (element)
         (result
          (make-array
           m-rows
           :element-type (array-element-type vector))))
    (dotimes (row m-rows result)
      (setq element zero)
      (dotimes (column n-columns)
        (setq
         element
         (+ element
            (* (aref array row column) (aref vector column)))))
      ;; Store the result
      (setf (aref result row) element))))

(defun %scaled-product-array-vector (array vector scalar)
  "Return the result of the array postmultiplied by the vector and
scaled."
  (let* ((m-rows (array-dimension array 0))
         (n-columns (array-dimension array 1))
         (zero (coerce 0 (array-element-type vector)))
         (element)
         (result
          (make-array
           m-rows
           :element-type (array-element-type vector))))
    (dotimes (row m-rows result)
      (setq element zero)
      (dotimes (column n-columns)
        (setq
         element
         (+ element
            (* (aref array row column) (aref vector column)))))
      ;; Store the result
      (setf (aref result row) (* scalar element)))))

(defun product-array-vector (array vector scalar)
  "Return the result of the array postmultiplied by the vector and
scaled."
  (if scalar
      (%scaled-product-array-vector array vector scalar)
      (%product-array-vector array vector)))

;;; Binary array operations

(defun %array<-array1-op-array2 (operation array1 array2)
  (let ((m-rows (array-dimension array1 0))
        (n-columns (array-dimension array1 1))
        (result
         (make-array
          (array-dimensions array1)
          :element-type
          (common-array-element-type array1 array2))))
    (dotimes (row m-rows result)
      (dotimes (column n-columns)
        (setf
         (aref result row column)
         (funcall operation
                  (aref array1 row column)
                  (aref array2 row column)))))))

(defun %array1<-array1-op-array2 (operation array1 array2)
  (let ((m-rows (array-dimension array1 0))
        (n-columns (array-dimension array1 1)))
    (dotimes (row m-rows array1)
      (dotimes (column n-columns)
        (setf
         (aref array1 row column)
         (funcall operation
                  (aref array1 row column)
                  (aref array2 row column)))))))

(defmethod compatible-dimensions-p
           ((operation (eql :add)) (array1 array) (array2 array))
  "Return true if the array dimensions are compatible for an
addition."
  (and
   (= 2 (array-rank array1) (array-rank array2))
   (= (array-dimension array1 0) (array-dimension array2 0))
   (= (array-dimension array1 1) (array-dimension array2 1))))

(defun add-array (array1 array2 scalar1 scalar2)
  "Array binary addition."
  (%array<-array1-op-array2
   (scaled-binary-op #'+ scalar1 scalar2)
   array1 array2))

(defun nadd-array (array1 array2 scalar1 scalar2)
  "Destructive array binary addition."
  (%array1<-array1-op-array2
   (scaled-binary-op #'+ scalar1 scalar2)
   array1 array2))

(defun subtract-array (array1 array2 scalar1 scalar2)
  "Array binary subtraction."
  (%array<-array1-op-array2
   (scaled-binary-op #'- scalar1 scalar2)
   array1 array2))

(defun nsubtract-array (array1 array2 scalar1 scalar2)
  "Destructive array binary subtraction."
  (%array1<-array1-op-array2
   (scaled-binary-op #'- scalar1 scalar2)
   array1 array2))

(defmethod compatible-dimensions-p
           ((operation (eql :product)) (array1 array) (array2 array))
  "Return true if the array dimensions are compatible for product."
  (and
   (= 2 (array-rank array1) (array-rank array2))
   (= (array-dimension array1 1) (array-dimension array2 0))))

(defun %product-array-array (array1 array2)
  "Return the result of the product of 2 arrays."
  (let* ((l-columns (array-dimension array1 1))
         (m-rows (array-dimension array1 0))
         (n-columns (array-dimension array2 1))
         (zero (coerce 0 (array-element-type array1)))
         (element)
         (result
          (make-array
           (list m-rows n-columns)
           :element-type (array-element-type array1))))
    (dotimes (row m-rows result)
      (dotimes (column n-columns)
        (setq element zero)
        (dotimes (index l-columns)
          (setq element
                (+ element
                   (* (aref array1 row index)
                      (aref array2 index column)))))
        ;; Store result
        (setf (aref result row column) element)))))

(defun %scaled-product-array-array (array1 array2 scalar)
  "Return the scaled result of the product of 2 arrays."
  (let* ((l-columns (array-dimension array1 1))
         (m-rows (array-dimension array1 0))
         (n-columns (array-dimension array2 1))
         (zero (coerce 0 (array-element-type array1)))
         (element)
         (result
          (make-array
           (list m-rows n-columns)
           :element-type (array-element-type array1))))
    (dotimes (row m-rows result)
      (dotimes (column n-columns)
        (setq element zero)
        (dotimes (index l-columns)
          (setq element
                (+ element
                   (* (aref array1 row index)
                      (aref array2 index column)))))
        ;; Store result
        (setf (aref result row column) (* scalar element))))))

(defun product-array-array (array1 array2 scalar)
  "Return the scaled result of the product of 2 arrays."
  (if scalar
      (%scaled-product-array-array array1 array2 scalar)
      (%product-array-array array1 array2)))
