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

(defmethod initialize-matrix ((matrix dense-matrix)
                              (data number)
                              (rows integer)
                              (columns integer)
                              element-type)
  "Initialize the dense matrix with an initial element."
  (setf
   (contents matrix)
   (make-array (list rows columns)
               :element-type element-type
               :initial-element data))
  ;; Return the matrix
  matrix)

(defmethod initialize-matrix ((matrix dense-matrix)
                              (data list)
                              (rows integer)
                              (columns integer)
                              element-type)
  "Initialize the dense matrix with a nested sequence."
  (setf
   (contents matrix)
   (make-array (list rows columns)
               :element-type element-type
               :initial-contents data))
  ;; Return the matrix
  matrix)

(defmethod initialize-matrix ((matrix dense-matrix)
                              (data vector)
                              (rows integer)
                              (columns integer)
                              element-type)
  "Initialize the dense matrix with a nested sequence."
  (setf
   (contents matrix)
   (make-array (list rows columns)
               :element-type element-type
               :initial-contents data))
  ;; Return the matrix
  matrix)

(defmethod initialize-matrix :before ((matrix dense-matrix)
                                      (data array)
                                      (rows integer)
                                      (columns integer)
                                      element-type)
  "Verify that the size of the data is valid."
  (when (vectorp data) (return-from initialize-matrix))
  (unless (= 2 (array-rank data))
    (error "data rank(~D) must equal 2." (array-rank data)))
  (unless (= rows (array-dimension data 0))
    (error "data rows(~D) does not equal matrix rows(~D)."
           (array-dimension data 0) rows))
  (unless (= columns (array-dimension data 1))
    (error "data columns(~D) does not equal matrix columns(~D)."
           (array-dimension data 1) columns))
  (unless (subtypep (array-element-type data)
                    (upgraded-array-element-type element-type))
    (error "Data type, ~A, is not of type ~A."
           (array-element-type data) element-type)))

(defmethod initialize-matrix ((matrix dense-matrix)
                              (data array)
                              (rows integer)
                              (columns integer)
                              element-type)
  "Initialize the dense matrix with a 2D array."
  (let ((contents
         (setf (contents matrix)
               (make-array (list rows columns)
                           :element-type element-type))))
    (dotimes (row rows matrix)
      (dotimes (column columns)
        (setf (aref contents row column) (aref data row column))))))

(defmethod matrix-in-bounds-p ((matrix dense-matrix)
                               (row integer) (column integer))
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

(defmethod (setf mref) ((data number) (matrix dense-matrix)
                        (row integer) (column integer))
  "Set the element of matrix at row,column."
  (setf (aref (contents matrix) row column) data))

(defmethod copy-matrix ((matrix dense-matrix))
  "Return a copy of the dense matrix."
  (let ((m-rows (matrix-row-dimension matrix))
        (n-columns (matrix-column-dimension matrix))
        (original (contents matrix))
        (contents
         (make-array
          (matrix-dimensions matrix)
          :element-type (matrix-element-type matrix))))
    (make-instance
     (class-of matrix)
     :contents
     (dotimes (row m-rows contents)
       (dotimes (column n-columns)
         (setf
          (aref contents row column)
          (aref original row column)))))))

(defmethod submatrix ((matrix dense-matrix)
                      (start-row integer)
                      (start-column integer)
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

(defmethod (setf submatrix) ((data dense-matrix)
                             (matrix dense-matrix)
                             (start-row integer)
                             (start-column integer)
                             &key end-row end-column)
  "Set the submatrix of matrix."
  (multiple-value-bind (start-row start-column end-row end-column)
      (matrix-validated-range
       matrix start-row start-column end-row end-column)
    (let ((m-rows (min (- end-row start-row)
                       (matrix-row-dimension data)))
          (n-columns (min (- end-column start-column)
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

(defmethod replace-matrix ((matrix1 dense-matrix)
                           (matrix2 dense-matrix)
                           &key
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
      (let ((m-rows (min (- end-row1 start-row1)
                         (- end-row2 start-row2)))
            (n-columns (min (- end-column1 start-column1)
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

(defmethod sumsq ((matrix dense-matrix) &key (scale 0) (sumsq 1))
  "Return the scaling parameter and the sum of the squares of the matrix."
  (sumsq-array (contents matrix) scale sumsq))

(defmethod sump ((matrix dense-matrix) (p number) &key (scale 0) (sump 1))
  "Return the scaling parameter and the sum of the P powers of the matrix."
  (sump-array (contents matrix) p scale sump))

(defmethod norm ((matrix dense-matrix) &key (measure 1))
  "Return the norm of the matrix."
  (norm-array (contents matrix) measure))

(defmethod transpose ((matrix dense-matrix) &key conjugate)
  "Return the transpose of the matrix."
  (make-instance
   (class-of matrix)
   :contents
   (let* ((m-rows (matrix-row-dimension matrix))
          (n-columns (matrix-column-dimension matrix))
          (op (if conjugate #'conjugate #'identity))
          (contents (contents matrix))
          (tcontents
           (make-array
            (list n-columns m-rows)
            :element-type
            (matrix-element-type matrix))))
     (dotimes (row m-rows tcontents)
       (dotimes (column n-columns)
         (setf (aref tcontents column row)
               (funcall op (aref contents row column))))))))

(defmethod ntranspose ((matrix dense-matrix) &key conjugate)
  "Replace the contents of the dense matrix with the transpose."
  (let ((m-rows (matrix-row-dimension matrix))
        (n-columns (matrix-column-dimension matrix)))
    (if (= m-rows n-columns)
        (let ((op (if conjugate #'conjugate #'identity))
              (contents (contents matrix)))
          (dotimes (row m-rows matrix)
            ;; FIXME : Conjugate on the diagonal may not be correct.
            (setf
             (aref contents row row)
             (funcall op (aref contents row row)))
            (do ((column (1+ row) (1+ column)))
                ((>= column n-columns))
              (psetf
               (aref contents row column)
               (funcall op (aref contents column row))
               (aref contents column row)
               (funcall op (aref contents row column))))))
        (error "Rows and columns unequal."))))

(defmethod permute ((matrix dense-matrix)
                    (permutation permutation-matrix))
  (if (every
       #'= (matrix-dimensions matrix) (matrix-dimensions permutation))
      (make-instance
       (class-of matrix)
       :contents
       (right-permute-array (contents matrix) (contents permutation)))
      (error
       "Dense matrix~A and permutation matrix~A sizes incompatible."
       (matrix-dimensions matrix)
       (matrix-dimensions permutation))))

(defmethod permute ((permutation permutation-matrix)
                    (matrix dense-matrix))
  (if (every
       #'= (matrix-dimensions permutation) (matrix-dimensions matrix))
      (make-instance
       (class-of matrix)
       :contents
       (left-permute-array (contents permutation) (contents matrix)))
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

(defmethod add :before ((matrix1 dense-matrix)
                        (matrix2 dense-matrix)
                        &key scalar1 scalar2)
  "Audit the input data."
  (declare (ignore scalar1 scalar2))
  (unless (equal (matrix-dimensions matrix1)
                 (matrix-dimensions matrix2))
    (error "The matrix dimensions are not compatible.")))

(defmethod add ((matrix1 dense-matrix) (matrix2 dense-matrix)
                &key scalar1 scalar2)
  "Return the addition of the 2 matrices."
  (make-instance
   (common-class-of matrix1 matrix2 'dense-matrix)
   :contents
   (add-array
    (contents matrix1) (contents matrix2) scalar1 scalar2)))

(defmethod nadd :before ((matrix1 dense-matrix)
                         (matrix2 dense-matrix)
                         &key scalar1 scalar2)
  "Audit the input data."
  (declare (ignore scalar1 scalar2))
  (unless (equal (matrix-dimensions matrix1)
                 (matrix-dimensions matrix2))
    (error "The matrix dimensions are not compatible.")))

(defmethod nadd ((matrix1 dense-matrix) (matrix2 dense-matrix)
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

(defmethod subtract ((matrix1 dense-matrix) (matrix2 dense-matrix)
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
  (unless (equal
           (matrix-dimensions matrix1)
           (matrix-dimensions matrix2))
    (error "The matrix dimensions are not compatible.")))

(defmethod nsubtract ((matrix1 dense-matrix) (matrix2 dense-matrix)
                      &key scalar1 scalar2)
  "Return the addition of the 2 matrices."
  (nsubtract-array
   (contents matrix1) (contents matrix2) scalar1 scalar2)
  matrix1)

(defmethod product :before
  ((vector row-vector) (matrix dense-matrix) &key scalar)
  "Verify the inputs."
  (declare (ignore scalar))
  (unless (= (vector-length vector) (matrix-row-dimension matrix))
    (error "Row vector(~D) is incompatible with matrix~A."
           (vector-length vector) (matrix-dimensions matrix))))

(defmethod product ((vector row-vector)
                    (matrix dense-matrix)
                    &key scalar)
  "Return a row vector generated by the pre-multiplication of a dense
matrix by a row vector."
  (make-instance
   'row-vector
   :contents
   (product-vector-array (contents vector) (contents matrix) scalar)))

(defmethod product :before
  ((matrix dense-matrix) (vector column-vector) &key scalar)
  "Verify the input."
  (declare (ignore scalar))
  (unless (= (vector-length vector) (matrix-column-dimension matrix))
    (error "Column vector(~D) is incompatible with matrix~A."
           (vector-length vector) (matrix-dimensions matrix))))

(defmethod product ((matrix dense-matrix)
                    (vector column-vector)
                    &key scalar)
  "Return a column vector generated by the multiplication of the dense
matrix with a column vector."
  (make-instance
   'column-vector
   :contents
   (product-array-vector (contents matrix) (contents vector) scalar)))

(defmethod product :before
  ((matrix1 dense-matrix) (matrix2 dense-matrix) &key scalar)
  "Verify the input."
  (declare (ignore scalar))
  (unless (= (matrix-column-dimension matrix1)
             (matrix-row-dimension matrix2))
    (error "The matrix dimensions, ~A and ~A, are not compatible."
           (matrix-dimensions matrix1) (matrix-dimensions matrix2))))

(defmethod product ((matrix1 dense-matrix)
                    (matrix2 dense-matrix)
                    &key scalar)
  "Return the product of the dense matrices."
  (make-instance
   (common-class-of matrix1 matrix2 'dense-matrix)
   :contents
   (product-array-array
    (contents matrix1) (contents matrix2) scalar)))
