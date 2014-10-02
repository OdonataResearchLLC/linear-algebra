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

(defclass permutation-matrix (matrix-object)
  ((contents
    :type     (array fixnum (*))
    :initarg  :contents
    :accessor contents))
  (:documentation
   "Permutation matrix object."))

(defun permutation-matrix-p (object)
  "Return true if object is a permutation-matrix."
  (typep object 'permutation-matrix))

(defmethod initialize-instance :after
  ((self permutation-matrix) &rest initargs
   &key dimensions element-type initial-element initial-contents)
  "Verify that the element-type was not set and that rows equals
columns."
  (cond
   ((slot-boundp self 'contents))
   ((not (eq 'number element-type))
    (error
     "Cannot specify the element type of a permutation matrix."))
   ((not (apply #'= dimensions))
    (error "Number of rows must equal the number of columns."))
   (initial-element
    (error
     "Cannot specify the initial element of a permutation matrix."))
   (initial-contents
    (initialize-matrix-contents self initial-contents initargs))
   (t
    (loop
     with rows = (first dimensions)
     with contents =
     (setf (contents self) (make-array rows :element-type 'fixnum))
     for index below rows do
     (setf (aref contents index) index)))))

;;; FIXME : Use the LOOP.
(defun %initialize-permutation-matrix-with-seq (matrix data size)
  (if (= size (length data))
      (let ((contents
             (setf
              (contents matrix)
              (make-array size :element-type 'fixnum))))
        ;; Fill contents, there should be no duplicates.
        (dotimes (row size)
          (let ((data-row (elt data row)))
            (if (= size (length data-row))
                (let ((column (position 1 data-row :test #'=)))
                  (if column
                      (setf (aref contents row) column)
                      (error "Invalid permutation data.")))
                (error "Rows unequal in length."))))
        ;; FIXME : Find a better way to identify duplicates.
        ;; If duplicates, not a permutation matrix.
        (unless
            (= size (length (remove-duplicates contents :test #'=)))
          (error "Invalid permutation in data."))
        ;; Return the matrix
        matrix)
      (error "Invalid number of rows of data.")))

(defmethod initialize-matrix-contents
    ((matrix permutation-matrix) (initial-contents list) initargs)
  "Initialize the permutation matrix with a list."
  (%initialize-permutation-matrix-with-seq
   matrix initial-contents (first (getf initargs :dimensions))))

(defmethod initialize-matrix-contents
    ((matrix permutation-matrix) (initial-contents vector) initargs)
  "Initialize the permutation matrix with a list."
  (%initialize-permutation-matrix-with-seq
   matrix initial-contents (first (getf initargs :dimensions))))

(defmethod initialize-matrix-contents
    ((matrix permutation-matrix) (initial-contents array) initargs)
  "Initialize the permutation matrix with a 2D array."
  (let ((rows (first (getf initargs :dimensions)))
        (columns (second (getf initargs :dimensions))))
    (cond
     ((not (= rows (array-dimension initial-contents 0)))
      (error "Invalid number of rows of data."))
     ((not (= columns (array-dimension initial-contents 1)))
      (error "Invalid number of columns of data."))
     (t
      (let ((row -1))
        (map-into
         (setf (contents matrix) (make-array rows :element-type 'fixnum))
         (lambda ()
           (incf row)
           (do ((column 0 (1+ column)))
               ((cond
                 ((>= column columns)
                  (error "Invalid permutation data."))
                 ((= 1 (aref initial-contents row column))))
                column))))
        ;; FIXME : Find a better way to identify duplicates.
        (unless
            (=
             rows
             (length (remove-duplicates (contents matrix) :test #'=)))
          (error "Invalid permutation in data."))
        ;; Return the permutation matrix
        matrix)))))

(defmethod matrix-in-bounds-p
    ((matrix permutation-matrix) (row integer) (column integer))
  "Return true if row and column do not exceed the dimensions of matrix."
  (let ((size (length (contents matrix))))
    (and
     (<= 0 row)    (< row    size)
     (<= 0 column) (< column size))))

(defmethod matrix-element-type ((matrix permutation-matrix))
  "Element type of the permutation matrix."
  'fixnum)

(defmethod matrix-dimensions ((matrix permutation-matrix))
  "Return the number of rows and columns in matrix."
  (let ((size (length (contents matrix))))
    (list size size)))

(defmethod matrix-row-dimension ((matrix permutation-matrix))
  "Return the number of rows in matrix."
  (length (contents matrix)))

(defmethod matrix-column-dimension ((matrix permutation-matrix))
  "Return the number of columns in matrix."
  (length (contents matrix)))

(defmethod mref
    ((matrix permutation-matrix) (row integer) (column integer))
  "Return 1 if a permutation and 0 otherwise."
  (if (= column (aref (contents matrix) row)) 1 0))

(defmethod (setf mref)
    ((data (eql 1)) (matrix permutation-matrix)
     (row integer) (column integer))
  "Swap rows of the permutation matrix."
  (let ((contents (contents matrix)))
    (rotatef (aref contents row) (aref contents column))))

(defmethod copy-matrix ((matrix permutation-matrix))
  "Return a copy of the permutation matrix."
  (make-instance
   'permutation-matrix
   :contents (copy-array (contents matrix))))

(defmethod transpose ((matrix permutation-matrix))
  "Transpose the permutation matrix."
  (make-instance
   'permutation-matrix
   :contents
   (loop
    with contents = (contents matrix)
    with permuted =
    (make-array (length contents) :element-type 'fixnum)
    for column across contents
    as  row = 0 then (1+ row)
    do (setf (aref permuted column) row)
    finally (return permuted))))
