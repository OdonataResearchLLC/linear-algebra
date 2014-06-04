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

(defclass dense-matrix (matrix-object)
  ((contents
    :type     (array * (* *))
    :initarg  :contents
    :accessor contents))
  (:documentation
   "Dense matrix object."))

;;; Dense matrix interface operations

(defun dense-matrix-p (object)
  "Return true if object is a dense matrix."
  (typep object 'dense-matrix))

(defmethod initialize-matrix
    ((matrix dense-matrix) (data number)
     (rows integer) (columns integer) element-type)
  "Initialize the dense matrix with an initial element."
  (setf
   (contents matrix)
   (make-array
    (list rows columns)
    :element-type element-type
    :initial-element data))
  ;; Return the matrix
  matrix)

(defmethod initialize-matrix
    ((matrix dense-matrix) (data list)
     (rows integer) (columns integer) element-type)
  "Initialize the dense matrix with a nested sequence."
  (setf
   (contents matrix)
   (make-array
    (list rows columns)
    :element-type element-type
    :initial-contents data))
  ;; Return the matrix
  matrix)

(defmethod initialize-matrix
    ((matrix dense-matrix) (data vector)
     (rows integer) (columns integer) element-type)
  "Initialize the dense matrix with a nested sequence."
  (setf
   (contents matrix)
   (make-array
    (list rows columns)
    :element-type element-type
    :initial-contents data))
  ;; Return the matrix
  matrix)

(defmethod initialize-matrix :before
  ((matrix dense-matrix) (data array)
   (rows integer) (columns integer) element-type)
  "Verify that the size of the data is valid."
  (when (vectorp data) (return-from initialize-matrix))
  ;; Rank 2
  (unless (= 2 (array-rank data))
    (error "data rank(~D) must equal 2." (array-rank data)))
  ;; Consistent number of rows
  (unless (= rows (array-dimension data 0))
    (error
     "data rows(~D) does not equal matrix rows(~D)."
     (array-dimension data 0) rows))
  ;; Consistent number of columns
  (unless (= columns (array-dimension data 1))
    (error
     "data columns(~D) does not equal matrix columns(~D)."
     (array-dimension data 1) columns))
  ;; Consistent type of data
  (unless (subtypep
           (array-element-type data)
           (upgraded-array-element-type element-type))
    (error
     "Data type, ~A, is not of type ~A."
     (array-element-type data) element-type)))

(defmethod initialize-matrix
    ((matrix dense-matrix) (data array)
     (rows integer) (columns integer) element-type)
  "Initialize the dense matrix with a 2D array."
  (setf (contents matrix) (copy-array data))
  ;; Return the matrix
  matrix)

(defmethod matrix-in-bounds-p
    ((matrix dense-matrix) (row integer) (column integer))
  "Return true if row and column do not exceed the dimensions of matrix."
  (array-in-bounds-p (contents matrix) row column))

(defmethod matrix-dimensions ((matrix dense-matrix))
  "Return the number of rows and columns in matrix."
  (array-dimensions (contents matrix)))

(defmethod matrix-row-dimension ((matrix dense-matrix))
  "Return the number of rows in matrix."
  (array-dimension (contents matrix) 0))

(defmethod matrix-column-dimension ((matrix dense-matrix))
  "Return the number of columns in matrix."
  (array-dimension (contents matrix) 1))

(defmethod matrix-element-type ((matrix dense-matrix))
  "Return the element type of the matrix."
  (array-element-type (contents matrix)))

(defmethod mref ((matrix dense-matrix) (row integer) (column integer))
  "Return the element of matrix at row,column."
  (aref (contents matrix) row column))

(defmethod (setf mref)
    ((data number) (matrix dense-matrix)
     (row integer) (column integer))
  "Set the element of matrix at row,column."
  (setf (aref (contents matrix) row column) data))

(defmethod copy-matrix ((matrix dense-matrix))
  "Return a copy of the dense matrix."
  (make-instance
   (class-of matrix)
   :contents
   (copy-array (contents matrix))))

(defmethod submatrix
    ((matrix dense-matrix)
     (start-row integer) (start-column integer)
     &key end-row end-column)
  "Return a dense matrix created from the submatrix of a matrix."
  (multiple-value-bind (start-row start-column end-row end-column)
      (matrix-validated-range
       matrix start-row start-column end-row end-column)
    (let* ((m-rows (- end-row start-row))
           (n-columns (- end-column start-column))
           (original (contents matrix))
           (contents
            (make-array
             (list m-rows n-columns)
             :element-type (matrix-element-type matrix))))
      (make-instance
       'dense-matrix
       :contents
       (dotimes (row m-rows contents)
         (dotimes (column n-columns)
           (setf
            (aref contents row column)
            (aref original
                  (+ start-row row)
                  (+ start-column column)))))))))

(defmethod (setf submatrix)
    ((data dense-matrix) (matrix dense-matrix)
     (start-row integer) (start-column integer)
     &key end-row end-column)
  "Set the submatrix of matrix."
  (multiple-value-bind (start-row start-column end-row end-column)
      (matrix-validated-range
       matrix start-row start-column end-row end-column)
    (let ((m-rows
           (min (- end-row start-row)
                (matrix-row-dimension data)))
          (n-columns
           (min (- end-column start-column)
                (matrix-column-dimension data)))
          (mat (contents matrix))
          (dat (contents data)))
      (do ((row0 0   (1+ row0))
           (row1 start-row (1+ row1)))
          ((>= row0 m-rows) data)       ; Return the data
        (do ((column0 0      (1+ column0))
             (column1 start-column (1+ column1)))
            ((>= column0 n-columns))
          (setf (aref mat row1 column1) (aref dat row0 column0)))))))

(defmethod replace-matrix
    ((matrix1 dense-matrix) (matrix2 dense-matrix) &key
     (start-row1 0) end-row1
     (start-column1 0) end-column1
     (start-row2 0) end-row2
     (start-column2 0) end-column2)
  "Replace the elements of matrix1 with matrix2."
  (multiple-value-bind (start-row1 start-column1 end-row1 end-column1)
      (matrix-validated-range
       matrix1 start-row1 start-column1 end-row1 end-column1)
    (multiple-value-bind (start-row2 start-column2 end-row2 end-column2)
        (matrix-validated-range
         matrix2 start-row2 start-column2 end-row2 end-column2)
      (let ((m-rows
             (min (- end-row1 start-row1)
                  (- end-row2 start-row2)))
            (n-columns
             (min (- end-column1 start-column1)
                  (- end-column2 start-column2)))
            (contents1 (contents matrix1))
            (contents2 (contents matrix2)))
        (do ((row 0 (1+ row))
             (row1 start-row1 (1+ row1))
             (row2 start-row2 (1+ row2)))
            ((>= row m-rows) matrix1)   ; Return MATRIX1
          (do ((column 0 (1+ column))
               (column1 start-column1 (1+ column1))
               (column2 start-column2 (1+ column2)))
              ((>= column n-columns))
            (setf (aref contents1 row1 column1)
                  (aref contents2 row2 column2))))))))

;;; Dense matrix fundamental operations

(defmethod norm ((matrix dense-matrix) &optional (measure 1))
  "Return the norm of the matrix."
  (norm-array (contents matrix) measure))

(defmethod transpose ((matrix dense-matrix))
  "Return the transpose of the matrix."
  (make-instance
   (class-of matrix)
   :contents
   (let* ((m-rows (matrix-row-dimension matrix))
          (n-columns (matrix-column-dimension matrix))
          (original (contents matrix))
          (transposed
           (make-array
            (list n-columns m-rows)
            :element-type
            (matrix-element-type matrix))))
     (dotimes (row m-rows transposed)
       (dotimes (column n-columns)
         (setf
          (aref transposed column row)
          (aref original row column)))))))

(defmethod ntranspose ((matrix dense-matrix))
  "Replace the contents of the dense matrix with the transpose."
  (let ((m-rows (matrix-row-dimension matrix))
        (n-columns (matrix-column-dimension matrix))
        (contents (contents matrix)))
    (if (= m-rows n-columns)
        (dotimes (row m-rows matrix)
          (do ((column (1+ row) (1+ column)))
              ((>= column n-columns))
            (rotatef
             (aref contents row column) (aref contents column row))))
        (error "Rows and columns unequal."))))

(defmethod permute
    ((matrix dense-matrix) (permutation permutation-matrix))
  (if (every
       #'= (matrix-dimensions matrix) (matrix-dimensions permutation))
      (make-instance
       (class-of matrix)
       :contents
       (right-permute (contents matrix) (contents permutation)))
      (error
       "Dense matrix~A and permutation matrix~A sizes incompatible."
       (matrix-dimensions matrix)
       (matrix-dimensions permutation))))

(defmethod permute
    ((permutation permutation-matrix) (matrix dense-matrix))
  (if (every
       #'= (matrix-dimensions permutation) (matrix-dimensions matrix))
      (make-instance
       (class-of matrix)
       :contents
       (left-permute (contents permutation) (contents matrix)))
      (error
       "Permutation matrix~A and dense matrix~A sizes incompatible."
       (matrix-dimensions matrix)
       (matrix-dimensions permutation))))

(defmethod scale ((scalar number) (matrix dense-matrix))
  "Scale each element of the dense matrix."
  (make-instance
   (class-of matrix)
   :contents
   (scale scalar (contents matrix))))

(defmethod nscale ((scalar number) (matrix dense-matrix))
  "Scale each element of the dense matrix."
  (nscale scalar (contents matrix))
  matrix)

(defmethod add :before
  ((matrix1 dense-matrix) (matrix2 dense-matrix)
   &key scalar1 scalar2)
  "Audit the input data."
  (declare (ignore scalar1 scalar2))
  (unless (equal (matrix-dimensions matrix1)
                 (matrix-dimensions matrix2))
    (error "The matrix dimensions are not compatible.")))

(defmethod add
    ((matrix1 dense-matrix) (matrix2 dense-matrix)
     &key scalar1 scalar2)
  "Return the addition of the 2 matrices."
  (make-instance
   (common-class-of matrix1 matrix2 'dense-matrix)
   :contents
   (add-array
    (contents matrix1) (contents matrix2) scalar1 scalar2)))

(defmethod nadd :before
  ((matrix1 dense-matrix) (matrix2 dense-matrix)
   &key scalar1 scalar2)
  "Audit the input data."
  (declare (ignore scalar1 scalar2))
  (unless (equal (matrix-dimensions matrix1)
                 (matrix-dimensions matrix2))
    (error "The matrix dimensions are not compatible.")))

(defmethod nadd
    ((matrix1 dense-matrix) (matrix2 dense-matrix)
     &key scalar1 scalar2)
  "Return the addition of the 2 matrices."
  (nadd-array
   (contents matrix1) (contents matrix2) scalar1 scalar2)
  matrix1)

(defmethod subtract :before
  ((matrix1 dense-matrix) (matrix2 dense-matrix)
   &key scalar1 scalar2)
  "Audit the input data."
  (declare (ignore scalar1 scalar2))
  (unless (equal (matrix-dimensions matrix1)
                 (matrix-dimensions matrix2))
    (error "The matrix dimensions are not compatible.")))

(defmethod subtract
    ((matrix1 dense-matrix) (matrix2 dense-matrix)
     &key scalar1 scalar2)
  "Return the addition of the 2 matrices."
  (make-instance
   (common-class-of matrix1 matrix2 'dense-matrix)
   :contents
   (subtract-array
    (contents matrix1) (contents matrix2) scalar1 scalar2)))

(defmethod nsubtract :before
  ((matrix1 dense-matrix) (matrix2 dense-matrix)
   &key scalar1 scalar2)
  "Audit the input data."
  (declare (ignore scalar1 scalar2))
  (unless
      (equal
       (matrix-dimensions matrix1)
       (matrix-dimensions matrix2))
    (error "The matrix dimensions are not compatible.")))

(defmethod nsubtract
    ((matrix1 dense-matrix) (matrix2 dense-matrix)
     &key scalar1 scalar2)
  "Return the addition of the 2 matrices."
  (nsubtract-array
   (contents matrix1) (contents matrix2) scalar1 scalar2)
  matrix1)

(defmethod product :before
  ((vector row-vector) (matrix dense-matrix) &optional scalar)
  "Verify the inputs."
  (declare (ignore scalar))
  (unless (= (vector-length vector) (matrix-row-dimension matrix))
    (error "Row vector(~D) is incompatible with matrix~A."
           (vector-length vector) (matrix-dimensions matrix))))

(defmethod product
    ((vector row-vector) (matrix dense-matrix) &optional scalar)
  "Return a row vector generated by the pre-multiplication of a dense
matrix by a row vector."
  (make-instance
   'row-vector
   :contents
   (product-vector-array (contents vector) (contents matrix) scalar)))

(defmethod product :before
  ((matrix dense-matrix) (vector column-vector) &optional scalar)
  "Verify the input."
  (declare (ignore scalar))
  (unless (= (vector-length vector) (matrix-column-dimension matrix))
    (error "Column vector(~D) is incompatible with matrix~A."
           (vector-length vector) (matrix-dimensions matrix))))

(defmethod product
    ((matrix dense-matrix) (vector column-vector) &optional scalar)
  "Return a column vector generated by the multiplication of the dense
matrix with a column vector."
  (make-instance
   'column-vector
   :contents
   (product-array-vector (contents matrix) (contents vector) scalar)))

(defmethod product :before
  ((matrix1 dense-matrix) (matrix2 dense-matrix) &optional scalar)
  "Verify the input."
  (declare (ignore scalar))
  (unless (= (matrix-column-dimension matrix1)
             (matrix-row-dimension matrix2))
    (error "The matrix dimensions, ~A and ~A, are not compatible."
           (matrix-dimensions matrix1) (matrix-dimensions matrix2))))

(defmethod product
    ((matrix1 dense-matrix) (matrix2 dense-matrix) &optional scalar)
  "Return the product of the dense matrices."
  (make-instance
   (common-class-of matrix1 matrix2 'dense-matrix)
   :contents
   (product-array-array
    (contents matrix1) (contents matrix2) scalar)))

(defmethod compatible-dimensions-p
    ((operation (eql :solve))
     (matrix dense-matrix)
     (vector column-vector))
  "Return true if the array dimensions are compatible for product."
  (= (matrix-row-dimension matrix)
     (matrix-column-dimension matrix)
     (vector-length vector)))

(defmethod solve :before
  ((matrix dense-matrix) (vector column-vector))
  "Return the solution to the system of equations."
  (unless (compatible-dimensions-p :solve matrix vector)
    (error "Matrix~A is incompatible with column vector(~D)."
           (matrix-dimensions matrix) (vector-length vector))))

(defmethod solve ((matrix dense-matrix) (vector column-vector))
  "Return the solution to the system of equations."
  (make-instance
   'column-vector
   :contents
   (gauss-solver
    (copy-array (contents matrix)) (copy-array (contents vector)))))

(defmethod nsolve :before
  ((matrix dense-matrix) (vector column-vector))
  "Return the solution to the system of equations."
  (unless (compatible-dimensions-p :solve matrix vector)
    (error "Matrix~A is incompatible with column vector(~D)."
           (matrix-dimensions matrix) (vector-length vector))))

(defmethod nsolve ((matrix dense-matrix) (vector column-vector))
  "Return the solution to the system of equations."
  (setf
   (contents vector)
   (gauss-solver (contents matrix) (contents vector)))
  ;; Return the solution vector
  vector)

(defmethod invert ((matrix dense-matrix))
  "Return the invert of the dense matrix."
  (if (= (matrix-row-dimension matrix)
         (matrix-column-dimension matrix))
      (make-instance
       (class-of matrix)
       :contents
       (gauss-invert (copy-array (contents matrix))))
      (error "The number of rows does not equal columns.")))

(defmethod ninvert ((matrix dense-matrix))
  "Return the invert of the dense matrix."
  (if (= (matrix-row-dimension matrix)
         (matrix-column-dimension matrix))
      (make-instance
       (class-of matrix)
       :contents
       (gauss-invert (contents matrix)))
      (error "The number of rows does not equal columns.")))
