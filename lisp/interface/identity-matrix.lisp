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

(defmethod initialize-instance :after
    ((self identity-matrix) &rest initargs
     &key dimensions element-type initial-element initial-contents)
  "Initialize the identity matrix."
  (cond
   ((slot-boundp self 'contents))
   ((or initial-element initial-contents)
    (error "Initial data is invalid for an identity matrix."))
   ((not (apply #'= dimensions))
    (error "Rows and columns are not equal."))
   (t
    (setf
     (slot-value self 'size) (first dimensions)
     (slot-value self 'contents)
     (make-array
      2 :element-type element-type
      :initial-contents
      (list (coerce 0 element-type) (coerce 1 element-type)))))))

(defmethod make-identity-matrix  (&key dimension element-type )
  (make-instance 'identity-matrix
		 :dimensions (list dimension dimension) :element-type element-type))


(defmethod matrix-in-bounds-p
    ((matrix identity-matrix) (row integer) (column integer))
  "Return true if row and column do not exceed the dimensions of matrix."
  (and
   (<= 0 row)    (< row    (size matrix))
   (<= 0 column) (< column (size matrix))))

(defmethod matrix-element-type ((matrix identity-matrix))
  "Return the element type of the identity matrix."
  (array-element-type (contents matrix)))

(defmethod matrix-dimensions ((matrix identity-matrix))
  "Return the number of rows and columns in matrix."
  (list (size matrix) (size matrix)))

(defmethod mref
    ((matrix identity-matrix) (row integer) (column integer))
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
      (list (coerce 0 element-type) (coerce 1 element-type))))))

