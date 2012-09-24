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

(defclass symmetric-matrix (hermitian-matrix)
  ()
  (:documentation
   "Symmetric matrix object."))

;;; Symmetric matrix interface opterations

(defun symmetric-matrix-p (object)
  "Return true if object is a symmetric-matrix, NIL otherwise."
  (typep object 'symmetric-matrix))

(defun %initialize-symmetric-matrix-with-seq (matrix data
                                              rows columns
                                              element-type)
  "Initialize and validate a symmetric matrix with a sequence."
  (let ((contents (setf (contents matrix)
                        (make-array (list rows columns)
                                    :element-type element-type
                                    :initial-contents data))))
    (dotimes (i0 rows matrix)
      (dotimes (i1 i0)
        (unless (number-equal
                 (aref contents i0 i1)
                 (aref contents i1 i0))
          (error "The data is not symmetric."))))))

(defmethod initialize-matrix ((matrix symmetric-matrix) (data list)
                              (rows integer) (columns integer)
                              &optional (element-type t))
  "Initialize a symmetric matrix."
  (%initialize-symmetric-matrix-with-seq matrix data
                                         rows columns
                                         element-type))

(defmethod initialize-matrix ((matrix symmetric-matrix) (data vector)
                              (rows integer) (columns integer)
                              &optional (element-type t))
  "Initialize a symmetric matrix."
  (%initialize-symmetric-matrix-with-seq matrix data
                                         rows columns
                                         element-type))

(defmethod initialize-matrix ((matrix symmetric-matrix) (data array)
                              (rows integer) (columns integer)
                              &optional (element-type t))
  "Initialize a symmetric matrix."
  (let ((contents (setf (contents matrix)
                        (make-array (list rows columns)
                                    :element-type element-type))))
    (dotimes (i0 rows matrix)
      (setf (aref contents i0 i0) (aref data i0 i0))
      (dotimes (i1 i0)
        (unless (number-equal
                 (setf (aref contents i0 i1)
                       (aref data i0 i1))
                 (setf (aref contents i1 i0)
                       (aref data i1 i0)))
          (error "The data is not symmetric."))))))

(defmethod (setf mref) ((data number) (matrix symmetric-matrix)
                        (row integer) (column integer))
  "Set the element of matrix at row,column."
  (setf (aref (contents matrix) row column) data
        (aref (contents matrix) column row) data))

(defun %setf-symmetric-submatrix-on-diagonal (matrix data row numrows)
  (let ((mat (contents matrix))
        (dat (contents data)))
    (do ((di0 0   (1+ di0))
         (mi0 row (1+ mi0)))
        ((>= di0 numrows) data)         ; Return the data
      (setf (aref mat mi0 mi0) (aref dat di0 di0))
      (do ((di1 (1+ di0)      (1+ di1))
           (mi1 (+ 1 row di0) (1+ mi1)))
          ((>= di1 numrows))
        (setf (aref mat mi0 mi1) (aref dat di0 di1)
              (aref mat mi1 mi0) (aref dat di1 di0))))))

(defun %setf-symmetric-submatrix-off-diagonal (matrix data
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
        (setf (aref mat mi0 mi1) (aref dat di0 di1)
              (aref mat mi1 mi0) (aref dat di0 di1))))))

(defmethod submatrix ((matrix symmetric-matrix)
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
              'symmetric-matrix)
             ((= numrows numcols)
              'square-matrix)
             (t 'dense-matrix))
       :contents
       (dotimes (i0 numrows contents)
         (dotimes (i1 numcols)
           (setf (aref contents i0 i1)
                 (aref original (+ row i0) (+ column i1)))))))))

(defmethod (setf submatrix) ((data symmetric-matrix)
                             (matrix symmetric-matrix)
                             (row integer) (column integer)
                             &key row-end column-end)
  "Set a submatrix of the matrix."
  (destructuring-bind (row column row-end column-end)
      (matrix-validated-range matrix row column row-end column-end)
    (let ((numrows (min (- row-end row)
                        (matrix-row-dimension data)))
          (numcols (min (- column-end column)
                        (matrix-column-dimension data))))
      (cond
        ((and (= row column) (= numrows numcols))
         (%setf-symmetric-submatrix-on-diagonal matrix data
                                                row numrows))
        ((or (<= (+ row numrows -1) column)
             (<= (+ column numcols -1) row))
         (%setf-symmetric-submatrix-off-diagonal matrix data
                                                 row column
                                                 numrows numcols))
        (t
         (error "Range(~D:~D,~D:~D) results in an asymmetric matrix."
                row row-end column column-end))))))

(defmethod (setf submatrix) ((data dense-matrix)
                             (matrix symmetric-matrix)
                             (row integer) (column integer)
                             &key row-end column-end)
  "Set a submatrix of MATRIX."
  (destructuring-bind (row column row-end column-end)
      (matrix-validated-range matrix row column row-end column-end)
    (let ((numrows (min (- row-end row)
                        (matrix-row-dimension data)))
          (numcols (min (- column-end column)
                        (matrix-column-dimension data))))
      (if (or (<= (+ row numrows -1) column)
              (<= (+ column numcols -1) row))
          (%setf-symmetric-submatrix-off-diagonal matrix data
                                                  row column
                                                  numrows numcols)
          (error "Range(~D:~D,~D:~D) results in an asymmetric matrix."
                 row row-end column column-end)))))

(defun %replace-symmetric-matrix-on-diagonal (matrix1 matrix2
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
      (setf (aref contents1 m1-i0 m1-i0) (aref contents2 m2-i0 m2-i0))
      (do ((   i1 (1+ i0)             (1+ i1))
           (m1-i1 (+ 1 column1 i0) (1+ m1-i1))
           (m2-i1 (+ 1 column2 i0) (1+ m2-i1)))
          ((>= i1 numcols))
        (setf (aref contents1 m1-i0 m1-i1) (aref contents2 m2-i0 m2-i1)
              (aref contents1 m1-i1 m1-i0) (aref contents2 m2-i1 m2-i0))))))

(defun %replace-symmetric-matrix-off-diagonal (matrix1 matrix2
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
        (setf (aref contents1 m1-i0 m1-i1) (aref contents2 m2-i0 m2-i1)
              (aref contents1 m1-i1 m1-i0) (aref contents2 m2-i0 m2-i1))))))

(defmethod replace-matrix ((matrix1 symmetric-matrix) (matrix2 symmetric-matrix)
                           &key (row1 0) row1-end (column1 0) column1-end
                           (row2 0) row2-end (column2 0) column2-end)
  "Replace the elements of MATRIX1 with MATRIX2."
  (destructuring-bind (row1 column1 row1-end column1-end)
      (matrix-validated-range matrix1 row1 column1 row1-end column1-end)
    (destructuring-bind (row2 column2 row2-end column2-end)
        (matrix-validated-range matrix2 row2 column2 row2-end column2-end)
      (let ((numrows (min (- row1-end row1) (- row2-end row2)))
            (numcols (min (- column1-end column1) (- column2-end column2))))
        (cond
          ((and (= row1 column1) (= row2 column2) (= numrows numcols))
           (%replace-symmetric-matrix-on-diagonal matrix1 matrix2
                                                  row1 column1
                                                  row2 column2
                                                  numrows numcols))
          ((or (<= (+ row1 numrows -1) column1)
               (<= (+ column1 numcols -1) row1))
           (%replace-symmetric-matrix-off-diagonal matrix1 matrix2
                                                   row1 column1
                                                   row2 column2
                                                   numrows numcols))
          (t
           (error "Range(~D:~D,~D:~D) results in an asymmetric matrix."
                  row1 (+ row1 numrows -1) column1 (+ column1 numcols -1))))))))

(defmethod replace-matrix ((matrix1 symmetric-matrix) (matrix2 dense-matrix)
                           &key (row1 0) row1-end (column1 0) column1-end
                           (row2 0) row2-end (column2 0) column2-end)
  "Replace the elements of MATRIX1 with MATRIX2."
  (destructuring-bind (row1 column1 row1-end column1-end)
      (matrix-validated-range matrix1 row1 column1 row1-end column1-end)
    (destructuring-bind (row2 column2 row2-end column2-end)
        (matrix-validated-range matrix2 row2 column2 row2-end column2-end)
      (let ((numrows (min (- row1-end row1) (- row2-end row2)))
            (numcols (min (- column1-end column1) (- column2-end column2))))
        (if (or (<= (+ row1 numrows -1) column1)
                (<= (+ column1 numcols -1) row1))
            (%replace-symmetric-matrix-off-diagonal matrix1 matrix2
                                                    row1 column1
                                                    row2 column2
                                                    numrows numcols)
            (error "Range(~D:~D,~D:~D) results in an asymmetric matrix."
                   row1 (+ row1 numrows -1) column1 (+ column1 numcols -1)))))))
