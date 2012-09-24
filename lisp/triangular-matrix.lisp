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

(defclass upper-triangular-matrix (square-matrix)
  ()
  (:documentation
   "Upper triangular matrix object."))

(defclass lower-triangular-matrix (square-matrix)
  ()
  (:documentation
   "Lower triangular matrix object."))

;;; Triangular matrix interface operations

(defun upper-triangular-matrix-p (object)
  "Return true if object is an upper triangular matrix."
  (typep object 'upper-triangular-matrix))

(defun lower-triangular-matrix-p (object)
  "Return true if object is a lower triangular matrix."
  (typep object 'lower-triangular-matrix))

(defun %initialize-upper-triangular-matrix-with-seq (matrix data
                                                     rows columns
                                                     element-type)
  (let ((contents (setf (contents matrix)
                        (make-array (list rows columns)
                                    :element-type element-type
                                    :initial-contents data))))
    (dotimes (i0 rows matrix)
      (dotimes (i1 i0)
        (unless (zerop (aref contents i0 i1))
          (error "Data is not upper triangular."))))))

(defun %initialize-lower-triangular-matrix-with-seq (matrix data
                                                     rows columns
                                                     element-type)
  (let ((contents (setf (contents matrix)
                        (make-array (list rows columns)
                                    :element-type element-type
                                    :initial-contents data))))
    (dotimes (i1 columns matrix)
      (dotimes (i0 i1)
        (unless (zerop (aref contents i0 i1))
          (error "Data is not lower triangular."))))))

(defmethod initialize-matrix ((matrix upper-triangular-matrix) (data number)
                              (rows fixnum) (columns fixnum)
                              &optional (element-type t))
  "Initialize the upper triangular matrix with an initial element."
  (let ((contents (setf (contents matrix)
                        (make-array (list rows columns)
                                    :element-type element-type
                                    :initial-element (coerce 0 element-type)))))
    (dotimes (i1 columns matrix)
      (setf (aref contents i1 i1) data)
      (dotimes (i0 i1)
        (setf (aref contents i0 i1) data)))))

(defmethod initialize-matrix ((matrix lower-triangular-matrix) (data number)
                              (rows fixnum) (columns fixnum)
                              &optional (element-type t))
  "Initialize the lower triangular matrix with an initial element."
  (let ((contents (setf (contents matrix)
                        (make-array (list rows columns)
                                    :element-type element-type
                                    :initial-element (coerce 0 element-type)))))
    (dotimes (i0 rows matrix)
      (setf (aref contents i0 i0) data)
      (dotimes (i1 i0)
        (setf (aref contents i0 i1) data)))))

(defmethod initialize-matrix ((matrix upper-triangular-matrix) (data list)
                              (rows integer) (columns integer)
                              &optional (element-type t))
  "Initialize the upper triangular matrix with a nested sequence."
  (%initialize-upper-triangular-matrix-with-seq matrix data
                                                rows columns
                                                element-type))

(defmethod initialize-matrix ((matrix lower-triangular-matrix) (data list)
                              (rows integer) (columns integer)
                              &optional (element-type t))
  "Initialize the lower triangular matrix with a nested sequence."
  (%initialize-lower-triangular-matrix-with-seq matrix data
                                                rows columns
                                                element-type))

(defmethod initialize-matrix ((matrix upper-triangular-matrix) (data vector)
                              (rows integer) (columns integer)
                              &optional (element-type t))
  "Initialize the upper triangular matrix with a nested sequence."
  (%initialize-upper-triangular-matrix-with-seq matrix data
                                                rows columns
                                                element-type))

(defmethod initialize-matrix ((matrix lower-triangular-matrix) (data vector)
                              (rows integer) (columns integer)
                              &optional (element-type t))
  "Initialize the lower triangular matrix with a nested sequence."
  (%initialize-lower-triangular-matrix-with-seq matrix data
                                                rows columns
                                                element-type))

(defmethod initialize-matrix ((matrix upper-triangular-matrix) (data array)
                              (rows integer) (columns integer)
                              &optional (element-type t))
  "Initialize the upper triangular matrix with a 2D array."
  (let ((contents (setf (contents matrix)
                        (make-array (list rows columns)
                                    :element-type element-type
                                    :initial-element (coerce 0 element-type)))))
    (dotimes (i1 columns matrix)
      (setf (aref contents i1 i1) (aref data i1 i1))
      (dotimes (i0 i1)
        (if (zerop (aref data i1 i0))
            (setf (aref contents i0 i1) (aref data i0 i1))
            (error "Data is not upper triangular."))))))

(defmethod initialize-matrix ((matrix lower-triangular-matrix) (data array)
                              (rows integer) (columns integer)
                              &optional (element-type t))
  "Initialize the lower triangular matrix with a 2D array."
  (let ((contents (setf (contents matrix)
                        (make-array (list rows columns)
                                    :element-type element-type
                                    :initial-element (coerce 0 element-type)))))
    (dotimes (i0 rows matrix)
      (setf (aref contents i0 i0) (aref data i0 i0))
      (dotimes (i1 i0)
        (if (zerop (aref data i1 i0))
            (setf (aref contents i0 i1) (aref data i0 i1))
            (error "Data is not lower triangular."))))))

(defmethod (setf mref) ((data number) (matrix upper-triangular-matrix)
                        (row fixnum) (column fixnum))
  "Set the element of matrix at row,column to data."
  (if (<= row column)
      (setf (aref (contents matrix) row column) data)
      (error "Elements below the diagonal must equal zero.")))

(defmethod (setf mref) ((data number) (matrix lower-triangular-matrix)
                        (row fixnum) (column fixnum))
  "Set the element of matrix at row,column to data."
  (if (<= column row)
      (setf (aref (contents matrix) row column) data)
      (error "Elements above the diagonal must equal zero.")))

(defmethod copy-matrix ((matrix upper-triangular-matrix))
  "Return a copy of the matrix."
  (let ((columns  (matrix-column-dimension matrix))
        (original (contents matrix))
        (contents (make-array (matrix-dimensions matrix)
                              :element-type
                              (matrix-element-type matrix))))
    (make-instance
     'upper-triangular-matrix
     :contents
     (dotimes (i1 columns contents)
       (setf (aref contents i1 i1) (aref original i1 i1))
       (dotimes (i0 i1)
         (setf (aref contents i0 i1) (aref original i0 i1)))))))

(defmethod copy-matrix ((matrix lower-triangular-matrix))
  "Return a copy of the matrix."
  (let ((rows     (matrix-row-dimension matrix))
        (original (contents matrix))
        (contents (make-array (matrix-dimensions matrix)
                              :element-type
                              (matrix-element-type matrix))))
    (make-instance
     'lower-triangular-matrix
     :contents
     (dotimes (i0 rows contents)
       (setf (aref contents i0 i0) (aref original i0 i0))
       (dotimes (i1 i0)
         (setf (aref contents i0 i1) (aref original i0 i1)))))))

(defun %setf-upper-triangular-submatrix-on-diagonal (matrix data row numrows)
  (let ((mat (contents matrix))
        (dat (contents data)))
    (do ((di0 0   (1+ di0))
         (mi0 row (1+ mi0)))
        ((>= di0 numrows) data)         ; Return the data
      (do ((di1 di0         (1+ di1))
           (mi1 (+ row di0) (1+ mi1)))
          ((>= di1 numrows))
        (setf (aref mat mi0 mi1) (aref dat di0 di1))))))

(defun %setf-lower-triangular-submatrix-on-diagonal (matrix data row numrows)
  (let ((mat (contents matrix))
        (dat (contents data)))
    (do ((di1 0   (1+ di1))
         (mi1 row (1+ mi1)))
        ((>= di1 numrows) data)         ; Return the data
      (do ((di0 di1         (1+ di0))
           (mi0 (+ row di1) (1+ mi0)))
          ((>= di0 numrows))
        (setf (aref mat mi0 mi1) (aref dat di0 di1))))))

(defun %setf-upper-triangular-submatrix-above-diagonal (matrix data
                                                        row column
                                                        numrows numcols)
  (let ((mat (contents matrix))
        (dat (contents data)))
    (do ((di0 0   (1+ di0))
         (mi0 row (1+ mi0)))
        ((>= di0 numrows) data)         ; Return the data
      (do ((di1 0      (1+ di1))
           (mi1 column (1+ mi1)))
          ((>= di1 numcols))
        (setf (aref mat mi0 mi1) (aref dat di0 di1))))))

(defun %setf-lower-triangular-submatrix-below-diagonal (matrix data
                                                        row column
                                                        numrows numcols)
  (let ((mat (contents matrix))
        (dat (contents data)))
    (do ((di1 0      (1+ di1))
         (mi1 column (1+ mi1)))
        ((>= di1 numcols) data)         ; Return the data
      (do ((di0 0   (1+ di0))
           (mi0 row (1+ mi0)))
          ((>= di0 numrows))
        (setf (aref mat mi0 mi1) (aref dat di0 di1))))))

(defmethod submatrix ((matrix upper-triangular-matrix)
                      (row integer) (column integer)
                      &key row-end column-end)
  "Return a matrix created from the submatrix of matrix."
  (destructuring-bind (row column row-end column-end)
      (matrix-validated-range matrix row column row-end column-end)
    (let* ((numrows (- row-end row))
           (numcols (- column-end column))
           (original (contents matrix))
           (contents (make-array (list numrows numcols)
                                 :element-type
                                 (matrix-element-type matrix))))
      (make-instance
       (cond ((and (= row column) (= numrows numcols))
              'upper-triangular-matrix)
             ((= numrows numcols)
              'square-matrix)
             (t 'dense-matrix))
       :contents
       (dotimes (i0 numrows contents)
         (dotimes (i1 numcols)
           (setf (aref contents i0 i1)
                 (aref original (+ row i0) (+ column i1)))))))))

(defmethod submatrix ((matrix lower-triangular-matrix)
                      (row integer) (column integer)
                      &key row-end column-end)
  "Return a matrix created from the submatrix of matrix."
  (destructuring-bind (row column row-end column-end)
      (matrix-validated-range matrix row column row-end column-end)
    (let* ((numrows (- row-end row))
           (numcols (- column-end column))
           (original (contents matrix))
           (contents (make-array (list numrows numcols)
                                 :element-type
                                 (matrix-element-type matrix))))
      (make-instance
       (cond ((and (= row column) (= numrows numcols))
              'lower-triangular-matrix)
             ((= numrows numcols)
              'square-matrix)
             (t 'dense-matrix))
       :contents
       (dotimes (i0 numrows contents)
         (dotimes (i1 numcols)
           (setf (aref contents i0 i1)
                 (aref original (+ row i0) (+ column i1)))))))))

(defmethod (setf submatrix) ((data upper-triangular-matrix)
                             (matrix upper-triangular-matrix)
                             (row integer) (column integer)
                             &key row-end column-end)
  "Set a submatrix of matrix with data."
  (destructuring-bind (row column row-end column-end)
      (matrix-validated-range matrix row column row-end column-end)
    (let ((numrows (min (- row-end row)
                        (matrix-row-dimension data)))
          (numcols (min (- column-end column)
                        (matrix-column-dimension data))))
      (cond
        ((and (= row column) (= numrows numcols))
         (%setf-upper-triangular-submatrix-on-diagonal matrix data
                                                       row numrows))
        ((<= (+ row numrows -1) column)
         (%setf-upper-triangular-submatrix-above-diagonal matrix data
                                                          row column
                                                          numrows numcols))
        (t
         (error "Range(~D:~D,~D:~D) results in a non upper triangular matrix."
                row row-end column column-end))))))

(defmethod (setf submatrix) ((data lower-triangular-matrix)
                             (matrix lower-triangular-matrix)
                             (row integer) (column integer)
                             &key row-end column-end)
  "Set a submatrix of matrix with data."
  (destructuring-bind (row column row-end column-end)
      (matrix-validated-range matrix row column row-end column-end)
    (let ((numrows (min (- row-end row)
                        (matrix-row-dimension data)))
          (numcols (min (- column-end column)
                        (matrix-column-dimension data))))
      (cond
        ((and (= row column) (= numrows numcols))
         (%setf-lower-triangular-submatrix-on-diagonal matrix data
                                                       row numrows))
        ((<= (+ column numcols -1) row)
         (%setf-lower-triangular-submatrix-below-diagonal matrix data
                                                          row column
                                                          numrows numcols))
        (t
         (error "Range(~D:~D,~D:~D) results in a non lower triangular matrix."
                row row-end column column-end))))))

(defmethod (setf submatrix) ((data dense-matrix)
                             (matrix upper-triangular-matrix)
                             (row integer) (column integer)
                             &key row-end column-end)
  "Set the submatrix of matrix with data."
  (destructuring-bind (row column row-end column-end)
      (matrix-validated-range matrix row column row-end column-end)
    (let ((numrows (min (- row-end row)
                        (matrix-row-dimension data)))
          (numcols (min (- column-end column)
                        (matrix-column-dimension data))))
      (if (<= (+ row numrows -1) column)
          (%setf-upper-triangular-submatrix-above-diagonal matrix data
                                                           row column
                                                           numrows numcols)
          (error "Range(~D:~D,~D:~D) results in a non upper triangular matrix."
                 row row-end column column-end)))))

(defmethod (setf submatrix) ((data dense-matrix)
                             (matrix lower-triangular-matrix)
                             (row integer) (column integer)
                             &key row-end column-end)
  "Set the submatrix of matrix with data."
  (destructuring-bind (row column row-end column-end)
      (matrix-validated-range matrix row column row-end column-end)
    (let ((numrows (min (- row-end row)
                        (matrix-row-dimension data)))
          (numcols (min (- column-end column)
                        (matrix-column-dimension data))))
      (if (<= (+ column numcols -1) row)
          (%setf-lower-triangular-submatrix-below-diagonal matrix data
                                                           row column
                                                           numrows numcols)
          (error "Range(~D:~D,~D:~D) results in a non lower triangular matrix."
                 row row-end column column-end)))))

(defun %replace-upper-triangular-matrix-on-diagonal (matrix1 matrix2
                                                     row1 column1
                                                     row2 column2
                                                     numrows numcols)
  "Destructively replace a subset on the diagonal of matrix1 with
matrix2."
  (let ((contents1 (contents matrix1))
        (contents2 (contents matrix2)))
    (do ((   i0 0    (1+ i0))
         (m1-i0 row1 (1+ m1-i0))
         (m2-i0 row2 (1+ m2-i0)))
        ((>= i0 numrows) matrix1)       ; Return matrix1
      (do ((   i1 i0             (1+ i1))
           (m1-i1 (+ column1 i0) (1+ m1-i1))
           (m2-i1 (+ column2 i0) (1+ m2-i1)))
          ((>= i1 numcols))
        (setf (aref contents1 m1-i0 m1-i1) (aref contents2 m2-i0 m2-i1))))))

(defun %replace-lower-triangular-matrix-on-diagonal (matrix1 matrix2
                                                     row1 column1
                                                     row2 column2
                                                     numrows numcols)
  "Destructively replace a subset on the diagonal of matrix1 with
matrix2."
  (let ((contents1 (contents matrix1))
        (contents2 (contents matrix2)))
    (do ((   i1 0       (1+ i1))
         (m1-i1 column1 (1+ m1-i1))
         (m2-i1 column2 (1+ m2-i1)))
        ((>= i1 numcols) matrix1)       ; Return matrix1
      (do ((   i0 i1          (1+ i0))
           (m1-i0 (+ row1 i1) (1+ m1-i0))
           (m2-i0 (+ row2 i1) (1+ m2-i0)))
          ((>= i0 numrows))
        (setf (aref contents1 m1-i0 m1-i1) (aref contents2 m2-i0 m2-i1))))))

(defun %replace-upper-triangular-matrix-above-diagonal (matrix1 matrix2
                                                        row1 column1
                                                        row2 column2
                                                        numrows numcols)
  "Destructively replace a subset off the diagonal of matrix1 with
matrix2."
  (let ((contents1 (contents matrix1))
        (contents2 (contents matrix2)))
    (do ((   i0 0    (1+ i0))
         (m1-i0 row1 (1+ m1-i0))
         (m2-i0 row2 (1+ m2-i0)))
        ((>= i0 numrows) matrix1)       ; Return matrix1
      (do ((   i1 0       (1+ i1))
           (m1-i1 column1 (1+ m1-i1))
           (m2-i1 column2 (1+ m2-i1)))
          ((>= i1 numcols))
        (setf (aref contents1 m1-i0 m1-i1) (aref contents2 m2-i0 m2-i1))))))

(defun %replace-lower-triangular-matrix-below-diagonal (matrix1 matrix2
                                                        row1 column1
                                                        row2 column2
                                                        numrows numcols)
  "Destructively replace a subset off the diagonal of matrix1 with
matrix2."
  (let ((contents1 (contents matrix1))
        (contents2 (contents matrix2)))
    (do ((   i0 0    (1+ i0))
         (m1-i0 row1 (1+ m1-i0))
         (m2-i0 row2 (1+ m2-i0)))
        ((>= i0 numrows) matrix1)       ; Return matrix1
      (do ((   i1 0       (1+ i1))
           (m1-i1 column1 (1+ m1-i1))
           (m2-i1 column2 (1+ m2-i1)))
          ((>= i1 numcols))
        (setf (aref contents1 m1-i0 m1-i1) (aref contents2 m2-i0 m2-i1))))))

(defmethod replace-matrix ((matrix1 upper-triangular-matrix)
                           (matrix2 upper-triangular-matrix)
                           &key (row1 0) row1-end (column1 0) column1-end
                           (row2 0) row2-end (column2 0) column2-end)
  "Replace the elements of matrix1 with matrix2."
  (destructuring-bind (row1 column1 row1-end column1-end)
      (matrix-validated-range matrix1 row1 column1 row1-end column1-end)
    (destructuring-bind (row2 column2 row2-end column2-end)
        (matrix-validated-range matrix2 row2 column2 row2-end column2-end)
      (let ((numrows (min (- row1-end row1) (- row2-end row2)))
            (numcols (min (- column1-end column1) (- column2-end column2))))
        (cond
          ((and (= row1 column1) (= row2 column2) (= numrows numcols))
           (%replace-upper-triangular-matrix-on-diagonal matrix1 matrix2
                                                         row1 column1
                                                         row2 column2
                                                         numrows numcols))
          ((<= (+ row1 numrows -1) column1)
           (%replace-upper-triangular-matrix-above-diagonal matrix1 matrix2
                                                            row1 column1
                                                            row2 column2
                                                            numrows numcols))
          (t
           (error "Range(~D:~D,~D:~D) results in a non upper triangular matrix."
                  row1 (+ row1 numrows -1) column1 (+ column1 numcols -1))))))))

(defmethod replace-matrix ((matrix1 lower-triangular-matrix)
                           (matrix2 lower-triangular-matrix)
                           &key (row1 0) row1-end (column1 0) column1-end
                           (row2 0) row2-end (column2 0) column2-end)
  "Replace the elements of matrix1 with matrix2."
  (destructuring-bind (row1 column1 row1-end column1-end)
      (matrix-validated-range matrix1 row1 column1 row1-end column1-end)
    (destructuring-bind (row2 column2 row2-end column2-end)
        (matrix-validated-range matrix2 row2 column2 row2-end column2-end)
      (let ((numrows (min (- row1-end row1) (- row2-end row2)))
            (numcols (min (- column1-end column1) (- column2-end column2))))
        (cond
          ((and (= row1 column1) (= row2 column2) (= numrows numcols))
           (%replace-lower-triangular-matrix-on-diagonal matrix1 matrix2
                                                         row1 column1
                                                         row2 column2
                                                         numrows numcols))
          ((<= (+ column1 numcols -1) row1)
           (%replace-lower-triangular-matrix-below-diagonal matrix1 matrix2
                                                            row1 column1
                                                            row2 column2
                                                            numrows numcols))
          (t
           (error "Range(~D:~D,~D:~D) results in a non lower triangular matrix."
                  row1 (+ row1 numrows -1) column1 (+ column1 numcols -1))))))))

(defmethod replace-matrix ((matrix1 upper-triangular-matrix)
                           (matrix2 dense-matrix)
                           &key (row1 0) row1-end (column1 0) column1-end
                           (row2 0) row2-end (column2 0) column2-end)
  "Replace the elements of matrix1 with matrix2."
  (destructuring-bind (row1 column1 row1-end column1-end)
      (matrix-validated-range matrix1 row1 column1 row1-end column1-end)
    (destructuring-bind (row2 column2 row2-end column2-end)
        (matrix-validated-range matrix2 row2 column2 row2-end column2-end)
      (let ((numrows (min (- row1-end row1) (- row2-end row2)))
            (numcols (min (- column1-end column1) (- column2-end column2))))
        (if (<= (+ row1 numrows -1) column1)
            (%replace-upper-triangular-matrix-above-diagonal matrix1 matrix2
                                                             row1 column1
                                                             row2 column2
                                                             numrows numcols)
            (error "Range(~D:~D,~D:~D) results in a non upper triangular matrix."
                   row1 (+ row1 numrows -1) column1 (+ column1 numcols -1)))))))

(defmethod replace-matrix ((matrix1 lower-triangular-matrix)
                           (matrix2 dense-matrix)
                           &key (row1 0) row1-end (column1 0) column1-end
                           (row2 0) row2-end (column2 0) column2-end)
  "Replace the elements of matrix1 with matrix2."
  (destructuring-bind (row1 column1 row1-end column1-end)
      (matrix-validated-range matrix1 row1 column1 row1-end column1-end)
    (destructuring-bind (row2 column2 row2-end column2-end)
        (matrix-validated-range matrix2 row2 column2 row2-end column2-end)
      (let ((numrows (min (- row1-end row1) (- row2-end row2)))
            (numcols (min (- column1-end column1) (- column2-end column2))))
        (if (<= (+ column1 numcols -1) row1)
            (%replace-lower-triangular-matrix-below-diagonal matrix1 matrix2
                                                             row1 column1
                                                             row2 column2
                                                             numrows numcols)
            (error "Range(~D:~D,~D:~D) results in a non lower triangular matrix."
                   row1 (+ row1 numrows -1) column1 (+ column1 numcols -1)))))))
