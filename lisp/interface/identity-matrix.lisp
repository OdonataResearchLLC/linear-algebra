#|

 Linear Algebra in Common Lisp

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

(in-package :linear-algebra)

(defclass identity-matrix (matrix-object)
  ((size
    :type    fixnum
    :initarg :size
    :reader  size
    :reader  matrix-row-dimension
    :reader  matrix-column-dimension)
   (contents
    :type    (array * (2))
    :initarg :contents
    :reader  contents))
  (:documentation
   "Identity matrix object."))

(defun identity-matrix-p (object)
  "Return true if object is an identity-matrix."
  (typep object 'identity-matrix))

(defmethod initialize-matrix ((matrix identity-matrix)
                              (data number)
                              (rows integer)
                              (columns integer)
                              element-type)
  "Initialize the identity matrix.."
  (cond
    ((not (zerop data))
     (error "Initial data is invalid for an identity matrix."))
    ((not (= rows columns))
     (error "Rows and columns are not equal."))
    (t
     (setf
      (slot-value matrix 'size) rows
      (slot-value matrix 'contents)
           (make-array 2 :element-type element-type
                       :initial-contents
                       (list (coerce 0 element-type)
                             (coerce 1 element-type))))
     ;; Return the matrix
     matrix)))

(defmethod matrix-in-bounds-p ((matrix identity-matrix)
                               (row integer) (column integer))
  "Return true if row and column do not exceed the dimensions of matrix."
  (and (<= 0 row)    (< row    (size matrix))
       (<= 0 column) (< column (size matrix))))

(defmethod matrix-element-type ((matrix identity-matrix))
  "Return the element type of the identity matrix."
  (array-element-type (contents matrix)))

(defmethod matrix-dimensions ((matrix identity-matrix))
  "Return the number of rows and columns in matrix."
  (list (size matrix) (size matrix)))

(defmethod mref ((matrix identity-matrix)
                 (row integer) (column integer))
  "Return the element of the matrix at row,column."
  (if (= row column)
      (aref (contents matrix) 1)
      (aref (contents matrix) 0)))

(defmethod copy-matrix ((matrix identity-matrix))
  "Return a copy of the matrix."
  (let ((element-type (matrix-element-type matrix)))
    (make-instance
     'identity-matrix
     :size (size matrix)
     :contents
     (make-array
      2 :element-type element-type
      :initial-contents
      (list (coerce 0 element-type)
            (coerce 1 element-type))))))

(defmethod submatrix ((matrix identity-matrix)
                      (start-row integer)
                      (start-column integer)
                      &key end-row end-column)
  "Return a matrix created from the submatrix of matrix."
  (multiple-value-bind (start-row start-column end-row end-column)
      (matrix-validated-range
       matrix start-row start-column end-row end-column)
    (let ((m-rows (- end-row start-row))
          (n-columns (- end-column start-column))
          (element-type (matrix-element-type matrix)))
      (cond
        ;; On the diagonal
        ((and (= start-row start-column) (= m-rows n-columns))
         (make-instance
          'identity-matrix :size m-rows
          :contents
          (make-array
           2 :element-type element-type
           :initial-contents
           (list (coerce 0 element-type)
                 (coerce 1 element-type)))))
        ;; Intersects the diagonal
        ((and (<= start-row end-column) (<= start-column end-row))
         (multiple-value-bind (r0 c0 size)
             (cond
               ((< start-row start-column)
                (values
                 (- start-column start-row)
                 0
                 (min n-columns (- end-row start-column))))
               ((< start-column start-row)
                (values
                 0
                 (- start-row start-column)
                 (min m-rows (- end-column start-row))))
               (t (values 0 0 (min m-rows n-columns))))
           (let ((one (coerce 1 element-type))
                 (contents
                  (make-array
                   (list m-rows n-columns)
                   :element-type element-type
                   :initial-element (coerce 0 element-type))))
             (make-instance
              (if (= m-rows n-columns) 'square-matrix 'dense-matrix)
              :contents
              (dotimes (index size contents)
                (setf
                 (aref contents (+ r0 index) (+ c0 index)) one))))))
        ;; A zero matrix
        (t
         (make-instance
          (if (= m-rows n-columns) 'square-matrix 'dense-matrix)
          :contents
          (make-array
           (list m-rows n-columns)
           :element-type element-type
           :initial-element (coerce 0 element-type))))))))
