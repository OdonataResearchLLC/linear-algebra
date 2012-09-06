#|

 Linear Algebra in Common Lisp Unit Tests

 Copyright (c) 2011, Thomas M. Hermann
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

(in-package :linear-algebra-test)

;;; Matrix object predicate
(define-test matrixp
  (assert-true (linear-algebra:matrixp
                (make-instance 'linear-algebra:matrix-object)))
  (assert-false (linear-algebra:matrixp t)))

;;; Matrix bounds
(defmacro test-matrix-in-bounds-p (matrix-type &optional
                                   (initial-contents nil initial-contents-p))
  (let ((mat (gensym "MATRIX-")))
    `(let ((,mat (linear-algebra:make-matrix
                  10 10 :matrix-type ,matrix-type
                  ,@(when initial-contents-p
                          `(:initial-contents ,initial-contents)))))
      (assert-true
       (linear-algebra:matrix-in-bounds-p
        ,mat (random 10) (random 10)))
      (assert-false
       (linear-algebra:matrix-in-bounds-p ,mat -1 0))
      (assert-false
       (linear-algebra:matrix-in-bounds-p ,mat 0 -1))
      (assert-false
       (linear-algebra:matrix-in-bounds-p ,mat 10 0))
      (assert-false
       (linear-algebra:matrix-in-bounds-p ,mat 0 10)))))

;;; Matrix element type
(defmacro test-matrix-element-type (matrix-type &optional
                                    (test-real-p t) (test-complex-p t))
  (let ((numeric-types (gensym "NUMERIC-TYPES-"))
        (mtype (gensym "MTYPE-"))
        (ntype (gensym "NTYPE-")))
    `(let ((,mtype ,matrix-type)
           (,numeric-types
            '(integer fixnum bit
              short-float single-float double-float long-float)))
      (dolist (,ntype ,numeric-types)
        ;; Real
        ,(when
          test-real-p
          `(assert-true
            (subtypep
             (linear-algebra:matrix-element-type
              (linear-algebra:make-matrix
               2 2 :matrix-type ,mtype :element-type ,ntype))
             (array-element-type
              (make-array '(2 2) :element-type ,ntype)))))
        ;; Complex
        ,(when
          test-complex-p
          `(assert-true
            (subtypep
             (linear-algebra:matrix-element-type
              (linear-algebra:make-matrix
               2 2 :matrix-type ,mtype :element-type `(complex ,,ntype)))
             (array-element-type
              (make-array '(2 2) :element-type `(complex ,,ntype))))))))))

;;; Matrix dimensions
(defmacro test-matrix-dimensions (matrix-type rows columns)
  `(assert-equal
    (list ,rows ,columns)
    (linear-algebra:matrix-dimensions
     (linear-algebra:make-matrix
      ,rows ,columns :matrix-type ,matrix-type))))

;;; Matrix row dimension
(defmacro test-matrix-row-dimension (matrix-type rows columns)
  `(assert-eq ,rows
    (linear-algebra:matrix-row-dimension
     (linear-algebra:make-matrix
      ,rows ,columns :matrix-type ,matrix-type))))

;;; Matrix column dimension
(defmacro test-matrix-column-dimension (matrix-type rows columns)
  `(assert-eq ,columns
    (linear-algebra:matrix-column-dimension
     (linear-algebra:make-matrix
      ,rows ,columns :matrix-type ,matrix-type))))

;;; Return a matrix modified with (SETF SUBMATRIX)
(defmacro setf-submatrix (rows columns matrix-type submatrix-form data-form)
  "Return the matrix modified with (SETF SUBMATRIX)"
  (let ((matrix (second submatrix-form)))
    `(let ((,matrix (zero-matrix ,rows ,columns
                                 :matrix-type ,matrix-type)))
      (setf ,submatrix-form ,data-form)
      ,matrix)))

;;; Matrix validated range
(defmacro test-matrix-validated-range (matrix-type rows columns)
  "Test that matrix-validate-range executes correctly."
  (let ((matrix  (gensym "MATRIX-"))
        (row1 (gensym "ROW1-"))
        (row2 (gensym "ROW2-"))
        (col1 (gensym "COL1-"))
        (col2 (gensym "COL2-")))
    `(let ((,matrix (linear-algebra:make-matrix
                     ,rows ,columns :matrix-type ,matrix-type))
           (,row1 (random ,rows))
           (,row2 (random ,rows))
           (,col1 (random ,columns))
           (,col2 (random ,columns)))
      (assert-equal
       (list ,row1 ,col1 ,rows ,columns)
       (linear-algebra:matrix-validated-range ,matrix ,row1 ,col1))
      (assert-equal
       (list (min ,row1 ,row2) ,col1 (max ,row1 ,row2) ,columns)
       (linear-algebra:matrix-validated-range
        ,matrix (min ,row1 ,row2) ,col1 (max ,row1 ,row2)))
      (assert-equal
       (list ,row1 (min ,col1 ,col2) ,rows (max ,col1 ,col2))
       (linear-algebra:matrix-validated-range
        ,matrix ,row1 (min ,col1 ,col2) nil (max ,col1 ,col2)))
      (assert-equal
       (list (min ,row1 ,row2) (min ,col1 ,col2)
             (max ,row1 ,row2) (max ,col1 ,col2))
       (linear-algebra:matrix-validated-range
        ,matrix
        (min ,row1 ,row2) (min ,col1 ,col2)
        (max ,row1 ,row2) (max ,col1 ,col2)))
      (assert-error
       'error
       (linear-algebra:matrix-validated-range ,matrix (1+ ,rows) ,col1))
      (assert-error
       'error
       (linear-algebra:matrix-validated-range ,matrix ,row1 (1+ ,columns)))
      (assert-error
       'error
       (linear-algebra:matrix-validated-range
        ,matrix (1- ,rows) ,col1 1))
      (assert-error
       'error
       (linear-algebra:matrix-validated-range
        ,matrix ,row1 (1- ,columns) ,rows 1))
      (assert-error
       'error
       (linear-algebra:matrix-validated-range
        ,matrix
        (1- ,rows) (1- ,columns) 1 1)))))

