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

(defmethod initialize-matrix ((matrix dense-matrix) (data number)
                              (rows integer) (columns integer)
                              &optional (element-type t))
  "Initialize the dense matrix with an initial element."
  (setf (contents matrix)
        (make-array (list rows columns)
                    :element-type element-type
                    :initial-element data))
  ;; Return the matrix
  matrix)

(defmethod initialize-matrix ((matrix dense-matrix) (data list)
                              (rows integer) (columns integer)
                              &optional (element-type t))
  "Initialize the dense matrix with a nested sequence."
  (setf (contents matrix)
        (make-array (list rows columns)
                    :element-type element-type
                    :initial-contents data))
  ;; Return the matrix
  matrix)

(defmethod initialize-matrix ((matrix dense-matrix) (data vector)
                              (rows integer) (columns integer)
                              &optional (element-type t))
  "Initialize the dense matrix with a nested sequence."
  (setf (contents matrix)
        (make-array (list rows columns)
                    :element-type element-type
                    :initial-contents data))
  ;; Return the matrix
  matrix)

(defmethod initialize-matrix :before ((matrix dense-matrix) (data array)
                                      (rows integer) (columns integer)
                                      &optional (element-type t))
  "Verify that the size of the data is valid."
  (when (vectorp data) (return-from initialize-matrix))
  (unless (= 2 (array-rank data))
    (error "Rank of array data is invalid."))
  (unless (= rows (array-dimension data 0))
    (error "Invalid number of rows of data."))
  (unless (= columns (array-dimension data 1))
    (error "Invalid number of columns of data."))
  (unless (subtypep (array-element-type data)
                    (upgraded-array-element-type element-type))
    (error "Data type, ~A, is not of type ~A."
           (array-element-type data) element-type)))

(defmethod initialize-matrix ((matrix dense-matrix) (data array)
                              (rows integer) (columns integer)
                              &optional (element-type t))
  "Initialize the dense matrix with a 2D array."
  (let ((contents
         (setf (contents matrix)
               (make-array (list rows columns)
                           :element-type element-type))))
    (dotimes (i0 rows matrix)
      (dotimes (i1 columns)
        (setf (aref contents i0 i1) (aref data i0 i1))))))

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
  (let ((rows     (matrix-row-dimension matrix))
        (columns  (matrix-column-dimension matrix))
        (original (contents matrix))
        (contents (make-array (matrix-dimensions matrix)
                              :element-type
                              (matrix-element-type matrix))))
    (make-instance
     (class-of matrix)
     :contents
     (dotimes (i0 rows contents)
       (dotimes (i1 columns)
         (setf (aref contents i0 i1) (aref original i0 i1)))))))

(defmethod submatrix ((matrix dense-matrix)
                      (row integer) (column integer)
                      &key row-end column-end)
  "Return a dense matrix created from the submatrix of a matrix."
  (destructuring-bind (row column row-end column-end)
      (matrix-validated-range matrix row column row-end column-end)
    (let* ((numrows (- row-end row))
           (numcols (- column-end column))
           (original (contents matrix))
           (contents (make-array (list numrows numcols)
                                 :element-type
                                 (matrix-element-type matrix))))
      (make-instance
       'dense-matrix
       :contents
       (dotimes (i0 numrows contents)
         (dotimes (i1 numcols)
           (setf (aref contents i0 i1)
                 (aref original (+ row i0) (+ column i1)))))))))

(defmethod (setf submatrix) ((data dense-matrix) (matrix dense-matrix)
                             (row integer) (column integer)
                             &key row-end column-end)
  "Set the submatrix of matrix."
  (destructuring-bind (row column row-end column-end)
      (matrix-validated-range matrix row column row-end column-end)
    (let ((numrows (min (- row-end row)
                        (matrix-row-dimension data)))
          (numcols (min (- column-end column)
                        (matrix-column-dimension data)))
          (mat (contents matrix))
          (dat (contents data)))
      (do ((di0 0   (1+ di0))
           (mi0 row (1+ mi0)))
          ((>= di0 numrows) data)       ; Return the data
        (do ((di1 0      (1+ di1))
             (mi1 column (1+ mi1)))
            ((>= di1 numcols))
          (setf (aref mat mi0 mi1) (aref dat di0 di1)))))))

(defmethod replace-matrix ((matrix1 dense-matrix) (matrix2 dense-matrix)
                           &key (row1 0) row1-end (column1 0) column1-end
                           (row2 0) row2-end (column2 0) column2-end)
  "Replace the elements of matrix1 with matrix2."
  (destructuring-bind (row1 column1 row1-end column1-end)
      (matrix-validated-range matrix1 row1 column1 row1-end column1-end)
    (destructuring-bind (row2 column2 row2-end column2-end)
        (matrix-validated-range matrix2 row2 column2 row2-end column2-end)
      (let ((numrows (min (- row1-end row1) (- row2-end row2)))
            (numcols (min (- column1-end column1) (- column2-end column2)))
            (contents1 (contents matrix1))
            (contents2 (contents matrix2)))
        (do ((i0    0    (1+ i0))
             (m1-i0 row1 (1+ m1-i0))
             (m2-i0 row2 (1+ m2-i0)))
            ((>= i0 numrows) matrix1)   ; Return MATRIX1
          (do ((i1    0       (1+ i1))
               (m1-i1 column1 (1+ m1-i1))
               (m2-i1 column2 (1+ m2-i1)))
              ((>= i1 numcols))
            (setf (aref contents1 m1-i0 m1-i1)
                  (aref contents2 m2-i0 m2-i1))))))))

;;; Dense matrix fundamental operations

(defmethod sumsq ((matrix dense-matrix) &key (scale 0) (sumsq 1))
  "Return the scaling parameter and the sum of the squares of the matrix."
  (destructuring-bind (numrows numcols) (matrix-dimensions matrix)
    (let ((mat (contents matrix))
          (abs-val 0))
      (dotimes (i0 numrows (values scale sumsq))
        (dotimes (i1 numcols)
          (when (< 0 (setf abs-val (abs (aref mat i0 i1))))
            (if (< scale abs-val)
                (setf sumsq (1+ (* sumsq (expt (/ scale abs-val) 2)))
                      scale abs-val)
                (incf sumsq (expt (/ abs-val scale) 2)))))))))

(defmethod sump ((matrix dense-matrix) (p number) &key (scale 0) (sump 1))
  "Return the scaling parameter and the sum of the P powers of the matrix."
  (unless (plusp p) (error "The power(~A) must be positive." p))
  (destructuring-bind (numrows numcols) (matrix-dimensions matrix)
    (let ((mat (contents matrix))
          (abs-val 0))
      (dotimes (i0 numrows (values scale sump))
        (dotimes (i1 numcols)
          (when (< 0 (setf abs-val (abs (aref mat i0 i1))))
            (if (< scale abs-val)
                (setf sump (1+ (* sump (expt (/ scale abs-val) p)))
                      scale abs-val)
                (incf sump (expt (/ (aref mat i0 i1) scale) p)))))))))

(defun %dense-matrix-1-norm (matrix)
  "Return the 1 norm of the matrix."
  (destructuring-bind (numrows numcols) (matrix-dimensions matrix)
    (let ((mat (contents matrix))
          (zero (coerce 0 (matrix-element-type matrix)))
          (norm 0)
          (sum 0))
      (dotimes (i1 numcols norm)
        (setf sum zero)
        (dotimes (i0 numrows)
          (incf sum (abs (aref mat i0 i1))))
        (setf norm (max sum norm))))))

(defun %dense-matrix-max-norm (matrix)
  "Return the max norm of the matrix."
  (destructuring-bind (numrows numcols) (matrix-dimensions matrix)
    (let ((mat (contents matrix))
          (norm 0))
      (dotimes (i0 numrows norm)
        (dotimes (i1 numcols)
          (setf norm (max norm (abs (aref mat i0 i1)))))))))

(defun %dense-matrix-frobenius-norm (matrix)
  "Return the Frobenius norm of the matrix."
  (multiple-value-bind (scale sumsq) (sumsq matrix)
    (* scale (sqrt sumsq))))

(defun %dense-matrix-infinity-norm (matrix)
  "Return the infinity norm of the matrix."
  (destructuring-bind (numrows numcols) (matrix-dimensions matrix)
    (let ((mat (contents matrix))
          (zero (coerce 0 (matrix-element-type matrix)))
          (norm 0)
          (sum 0))
      (dotimes (i0 numrows norm)
        (setf sum zero)
        (dotimes (i1 numcols)
          (incf sum (abs (aref mat i0 i1))))
        (setf norm (max sum norm))))))

(defmethod norm ((matrix dense-matrix) &key (measure 1))
  "Return the norm of the matrix."
  (case measure
    (1          (%dense-matrix-1-norm matrix))
    (:max       (%dense-matrix-max-norm matrix))
    (:frobenius (%dense-matrix-frobenius-norm matrix))
    (:infinity  (%dense-matrix-infinity-norm matrix))
    (otherwise  (error "Unrecognized norm, ~A." measure))))

(defmethod transpose ((matrix dense-matrix) &key conjugate)
  "Return the transpose of the matrix."
  (make-instance
   (class-of matrix)
   :contents
   (destructuring-bind (numrows numcols)
       (matrix-dimensions matrix)
     (let ((op (if conjugate #'conjugate #'identity))
           (contents  (contents matrix))
           (tcontents (make-array
                       (list numcols numrows)
                       :element-type
                       (matrix-element-type matrix))))
       (dotimes (i0 numrows tcontents)
         (dotimes (i1 numcols)
           (setf (aref tcontents i1 i0)
                 (funcall op (aref contents i0 i1)))))))))

(defmethod ntranspose ((matrix dense-matrix) &key conjugate)
  "Replace the contents of the dense matrix with the transpose."
  (destructuring-bind (numrows numcols) (matrix-dimensions matrix)
    (if (= numrows numcols)
        (let ((op (if conjugate #'conjugate #'identity))
              (contents (contents matrix)))
          (dotimes (i0 numrows matrix)
            ;; FIXME : Conjugate on the diagonal may not be correct.
            (setf (aref contents i0 i0) (funcall op (aref contents i0 i0)))
            (do ((i1 (1+ i0) (1+ i1)))
                ((>= i1 numcols)) 
              (psetf
               (aref contents i0 i1) (funcall op (aref contents i1 i0))
               (aref contents i1 i0) (funcall op (aref contents i0 i1))))))
        (error "Rows and columns unequal."))))

(defmethod scale ((scalar number) (matrix dense-matrix))
  "Scale each element of the dense matrix."
  (make-instance
   (class-of matrix)
   :contents
   (destructuring-bind (numrows numcols)
       (matrix-dimensions matrix)
     (let ((contents (contents matrix))
           (scaled (make-array (list numrows numcols)
                               :element-type
                               (matrix-element-type matrix))))
       (dotimes (i0 numrows scaled)
         (dotimes (i1 numcols)
           (setf (aref scaled i0 i1)
                 (* scalar (aref contents i0 i1)))))))))

(defmethod nscale ((scalar number) (matrix dense-matrix))
  "Scale each element of the dense matrix."
  (destructuring-bind (numrows numcols)
      (matrix-dimensions matrix)
    (let ((contents (contents matrix)))
      (dotimes (i0 numrows matrix)
        (dotimes (i1 numcols)
          (setf (aref contents i0 i1)
                (* scalar (aref contents i0 i1))))))))

(defmethod add :before ((matrix1 dense-matrix) (matrix2 dense-matrix)
                        &key scalar1 scalar2)
  "Audit the input data."
  (declare (ignore scalar1 scalar2))
  (unless (equal (matrix-dimensions matrix1)
                 (matrix-dimensions matrix2))
    (error "The matrix dimensions are not compatible.")))

(defmethod add ((matrix1 dense-matrix) (matrix2 dense-matrix)
                &key scalar1 scalar2)
  "Return the addition of the 2 matrices."
  (let ((numrows (matrix-row-dimension matrix1))
        (numcols (matrix-column-dimension matrix1))
        (contents1 (contents matrix1))
        (contents2 (contents matrix2))
        (contents (make-array
                   (matrix-dimensions matrix1)
                   :element-type
                   (common-array-element-type matrix1 matrix2)))
        (op (scaled-binary-op #'+ scalar1 scalar2)))
    (make-instance
     (common-class-of matrix1 matrix2 'dense-matrix)
     :contents
     (dotimes (i0 numrows contents)
       (dotimes (i1 numcols)
         (setf (aref contents i0 i1)
               (funcall op
                        (aref contents1 i0 i1)
                        (aref contents2 i0 i1))))))))

(defmethod product :before ((vector row-vector) (matrix dense-matrix)
                            &key scalar)
  "Verify the inputs."
  (declare (ignore scalar))
  (unless (= (vector-length vector) (matrix-row-dimension matrix))
    (error "Row vector(~D) is incompatible with matrix~A."
           (vector-length vector) (matrix-dimensions matrix))))

(defmethod product ((vector row-vector) (matrix dense-matrix)
                    &key scalar)
  "Return a row vector generated by the pre-multiplication of a dense
matrix by a row vector."
  (destructuring-bind (numrows numcols) (matrix-dimensions matrix)
    (let ((vec (contents vector))
          (mat (contents matrix))
          (zero (coerce 0 (vector-element-type vector)))
          (val  nil)
          (newvec (make-array
                   numcols
                   :element-type
                   (vector-element-type vector))))
      (make-instance
       'row-vector
       :contents
       (dotimes (i1 numcols newvec)
         (setf val zero)
         (dotimes (i0 numrows)
           (incf val (* (aref vec i0) (aref mat i0 i1))))
         (if scalar
             (setf (aref newvec i1) (* scalar val))
             (setf (aref newvec i1) val)))))))

(defmethod product :before ((matrix dense-matrix) (vector column-vector)
                            &key scalar)
  "Verify the input."
  (declare (ignore scalar))
  (unless (= (vector-length vector) (matrix-column-dimension matrix))
    (error "Column vector(~D) is incompatible with matrix~A."
           (vector-length vector) (matrix-dimensions matrix))))

(defmethod product ((matrix dense-matrix) (vector column-vector)
                    &key scalar)
  "Return a column vector generated by the multiplication of the dense
matrix with a column vector."
  (destructuring-bind (numrows numcols) (matrix-dimensions matrix)
    (let ((vec (contents vector))
          (mat (contents matrix))
          (zero (coerce 0 (vector-element-type vector)))
          (val nil)
          (newvec (make-array
                   numrows
                   :element-type
                   (vector-element-type vector))))
      (make-instance
       'column-vector
       :contents
       (dotimes (i0 numrows newvec)
         (setf val zero)
         (dotimes (i1 numcols)
           (incf val (* (aref mat i0 i1) (aref vec i1))))
         (if scalar
             (setf (aref newvec i0) (* scalar val))
             (setf (aref newvec i0) val)))))))

(defmethod product :before ((matrix1 dense-matrix) (matrix2 dense-matrix)
                            &key scalar)
  "Verify the input."
  (declare (ignore scalar))
  (unless (= (matrix-column-dimension matrix1) (matrix-row-dimension matrix2))
    (error "The matrix dimensions, ~A and ~A, are not compatible."
           (matrix-dimensions matrix1) (matrix-dimensions matrix2))))

(defmethod product ((matrix1 dense-matrix) (matrix2 dense-matrix)
                    &key scalar)
  "Return the product of the dense matrices."
  (destructuring-bind (numrow1 numcol1) (matrix-dimensions matrix1)
    (let* ((mat1 (contents matrix1))
           (mat2 (contents matrix2))
           (zero (coerce 0 (matrix-element-type matrix1)))
           (val nil)
           (numcol2 (matrix-column-dimension matrix2))
           (newmat (make-array (list numrow1 numcol2)
                               :element-type
                               (matrix-element-type matrix1))))
      (make-instance
       (common-class-of matrix1 matrix2 'dense-matrix)
       :contents
       (dotimes (i0 numrow1 newmat)
         (dotimes (i2 numcol2)
           (setf val zero)
           (dotimes (i1 numcol1)
             (incf val (* (aref mat1 i0 i1) (aref mat2 i1 i2))))
           (if scalar
               (setf (aref newmat i0 i2) (* scalar val))
               (setf (aref newmat i0 i2) val))))))))
