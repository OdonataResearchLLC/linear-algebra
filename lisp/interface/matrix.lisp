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

;;; Matrix superclass

(defclass matrix-object ()
  ()
  (:documentation
   "A superclass for all matrices."))

;;; Matrix interface operations

(defgeneric initialize-matrix (matrix data rows columns
                               &optional element-type)
  (:documentation
   "Initialize the matrix with data."))

(defun make-matrix (rows columns &key
                    (matrix-type 'dense-matrix)
                    (element-type t)
                    (initial-element nil initial-element-p)
                    (initial-contents nil initial-contents-p))
  "Return a new matrix instance."
  (let ((new-matrix (make-instance matrix-type)))
    (cond
      ((and initial-element-p initial-contents-p)
       (error "Cannot specify both INITIAL-ELEMENT and INITIAL-CONTENTS."))
      (initial-contents-p
       (initialize-matrix
        new-matrix initial-contents rows columns element-type))
      (initial-element-p
       (initialize-matrix
        new-matrix initial-element rows columns element-type))
      (t
       (initialize-matrix
        new-matrix (coerce 0 element-type) rows columns element-type)))))

(defun matrixp (object)
  "Return true if object is a matrix, NIL otherwise."
  (typep object 'matrix-object))

(defgeneric matrix-in-bounds-p (matrix row column)
  (:documentation
   "Return true if ROW and COLUMN do not exceed the dimensions of MATRIX."))

(defgeneric matrix-element-type (matrix)
  (:documentation
   "Return the element type of MATRIX."))

(defgeneric matrix-dimensions (matrix)
  (:documentation
   "Return the number of rows and columns in MATRIX."))

(defgeneric matrix-row-dimension (matrix)
  (:documentation
   "Return the number of rows in MATRIX."))

(defgeneric matrix-column-dimension (matrix)
  (:documentation
   "Return the number of columns in MATRIX."))

(defgeneric mref (matrix row column)
  (:documentation
   "Return the matrix element at ROW,COLUMN."))

(defgeneric (setf mref) (data matrix row column)
  (:documentation
   "Set the element at row,column of matrix to data."))

(defgeneric copy-matrix (matrix)
  (:documentation
   "Return a copy of the matrix."))

(defgeneric submatrix (matrix row column &key row-end column-end)
  (:documentation
   "Return a submatrix of the matrix."))

(defgeneric (setf submatrix) (submatrix matrix row column
                              &key row-end column-end)
  (:documentation
   "Set the submatrix of the matrix."))

(defgeneric replace-matrix (matrix1 matrix2 &key
                            row1 row1-end column1 column1-end
                            row2 row2-end column2 column2-end)
  (:documentation
   "Destructively replace elements of matrix1 with matrix2."))

(defun matrix-validated-range (matrix row column &optional row-end column-end)
  "Returns a validated range of rows and columns for the matrix."
  (destructuring-bind (row-dimension column-dimension)
      (matrix-dimensions matrix)
    (let ((row-end    (or row-end row-dimension))
          (column-end (or column-end column-dimension)))
      (if (and (<= 0 row row-end row-dimension)
               (<= 0 column column-end column-dimension))
          (list row column row-end column-end)
          (error "The matrix range (~D:~D,~D:~D) is invalid."
                 row column row-end column-end)))))
