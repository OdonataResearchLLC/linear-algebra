#|

 Linear Algebra in Common Lisp Unit Tests

 Copyright (c) 2010-2012, Thomas M. Hermann
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

(define-test make-square-matrix
  (:tag :square-matrix :make-matrix)
  ;; A default square matrix
  (let ((matrix
         (linear-algebra:make-matrix
          10 10 :matrix-type
          'linear-algebra:square-matrix)))
    (assert-true (linear-algebra:matrixp matrix))
    (assert-true (typep matrix 'linear-algebra:square-matrix))
    (assert-rational-equal
     (make-array '(10 10) :initial-element 0)
     matrix))
  ;; Specify the square matrix element type
  (let ((matrix
         (linear-algebra:make-matrix
          10 10
          :matrix-type 'linear-algebra:square-matrix
          :element-type 'single-float)))
    (assert-true (linear-algebra:matrixp matrix))
    (assert-true (typep matrix 'linear-algebra:square-matrix))
    (assert-eq
     (array-element-type
      (linear-algebra::contents matrix))
     (array-element-type
      (make-array '(10 10) :element-type 'single-float)))
    (assert-float-equal
     (make-array '(10 10) :initial-element 0.0
                 :element-type 'single-float)
     matrix))
  ;; Specify the square matrix initial element
  (let ((matrix
         (linear-algebra:make-matrix
          10 10
          :matrix-type 'linear-algebra:square-matrix
          :initial-element 1.0)))
    (assert-true (linear-algebra:matrixp matrix))
    (assert-true (typep matrix 'linear-algebra:square-matrix))
    (assert-float-equal
     (make-array '(10 10) :initial-element 1.0)
     matrix))
  ;; Specify the square matrix contents - Nested list
  (let* ((data
          '((1.1 1.2 1.3 1.4)
            (2.1 2.2 2.3 2.4)
            (3.1 3.2 3.3 3.4)
            (4.1 4.2 4.3 4.4))) 
         (matrix
          (linear-algebra:make-matrix
           4 4
           :matrix-type 'linear-algebra:square-matrix
           :initial-contents data)))
    (assert-true (linear-algebra:matrixp matrix))
    (assert-true (typep matrix 'linear-algebra:square-matrix))
    (assert-float-equal
     (make-array '(4 4) :initial-contents data)
     matrix))
  ;; Specify the square matrix contents - Nested vector
  (let* ((data
          #(#(1.1 1.2 1.3 1.4)
            #(2.1 2.2 2.3 2.4)
            #(3.1 3.2 3.3 3.4)
            #(4.1 4.2 4.3 4.4)))
         (matrix
          (linear-algebra:make-matrix
           4 4
           :matrix-type 'linear-algebra:square-matrix
           :initial-contents data)))
    (assert-true (linear-algebra:matrixp matrix))
    (assert-true (typep matrix 'linear-algebra:square-matrix))
    (assert-float-equal
     (make-array '(4 4) :initial-contents data)
     matrix))
  ;; Specify the square matrix contents - 2D array
  (let* ((data
          (make-array
           '(4 4) :initial-contents
           '((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4))))
         (matrix
          (linear-algebra:make-matrix
           4 4
           :matrix-type 'linear-algebra:square-matrix
           :initial-contents data)))
    (assert-true (linear-algebra:matrixp matrix))
    (assert-true (typep matrix 'linear-algebra:square-matrix))
    (assert-float-equal data matrix))
  ;; Erroneous 2D array input data
  (assert-error
   'error
   (linear-algebra:make-matrix
    4 4
    :matrix-type 'linear-algebra:square-matrix
    :initial-contents
    #3A(((1.1 1.2) (2.1 2.2))
        ((3.1 3.2) (4.1 4.2))
        ((5.1 5.2) (6.1 6.2)))))
  (assert-error
   'error
   (linear-algebra:make-matrix
    3 4
    :matrix-type 'linear-algebra:square-matrix
    :initial-contents
    (coordinate-array 0 0 3 4)))
  (assert-error
   'error
   (linear-algebra:make-matrix
    4 3
    :matrix-type 'linear-algebra:square-matrix
    :initial-contents
    (coordinate-array 0 0 4 3)))
  (assert-error
   'error
   (linear-algebra:make-matrix
    3 3 :element-type 'single-float
    :matrix-type 'linear-algebra:square-matrix
    :initial-contents
    #2A((1 2 3) (4 5 6) (7 8 9))))
  ;; Specify initial element and initial contents
  (assert-error
   'error
   (linear-algebra:make-matrix
    4 4
    :matrix-type 'linear-algebra:square-matrix
    :initial-element 1.1
    :initial-contents
    (coordinate-array 0 0 4 4))))

;;; Test the square matrix predicate

(define-test square-matrix-predicate
  (:tag :square-matrix)
  (assert-true
   (linear-algebra:square-matrix-p
    (linear-algebra:make-matrix
     10 10 :matrix-type 'linear-algebra:square-matrix)))
  (assert-false
   (linear-algebra:square-matrix-p (make-array '(10 10)))))

;;; Test the square matrix bounds

(define-test square-matrix-in-bounds-p
  (:tag :square-matrix :matrix-in-bounds-p)
  (test-matrix-in-bounds-p 'linear-algebra:square-matrix))

;;; Test the square matrix element type

(define-test square-matrix-element-type
  (:tag :square-matrix :matrix-element-type)
  (test-matrix-element-type 'linear-algebra:square-matrix))

;;; Test the square matrix dimensions

(define-test square-matrix-dimensions
  (:tag :square-matrix :matrix-dimensions)
  (test-matrix-dimensions 'linear-algebra:square-matrix 9 9))

;;; Test the square matrix row dimension

(define-test square-matrix-row-dimension
  (:tag :square-matrix :matrix-row-dimension)
  (test-matrix-row-dimension 'linear-algebra:square-matrix 9 9))

;;; Test the square matrix column dimension

(define-test square-matrix-column-dimension
  (:tag :square-matrix :matrix-column-dimension)
  (test-matrix-column-dimension 'linear-algebra:square-matrix 9 9))

;;; Reference square matrix elements

(define-test square-matrix-mref
  (:tag :square-matrix :mref)
  (let* ((initial-contents
          '((1.1 1.2 1.3 1.4 1.5)
            (2.1 2.2 2.3 2.4 2.5)
            (3.1 3.2 3.3 3.4 3.5)
            (4.1 4.2 4.3 4.4 4.5)
            (5.1 5.2 5.3 5.4 5.5)))
         (rows 5) (columns 5)
         (rend (1- rows)) (cend (1- columns))
         (rowi (random-interior-index rows))
         (coli (random-interior-index columns))
         (data
          (make-array
           (list rows columns)
           :initial-contents
           initial-contents))
         (matrix
          (linear-algebra:make-matrix
           rows columns
           :matrix-type
           'linear-algebra:square-matrix
           :initial-contents
           initial-contents)))
    (assert-float-equal
     (aref data 0 0)
     (linear-algebra:mref matrix 0 0))
    (assert-float-equal
     (aref data 0 cend)
     (linear-algebra:mref matrix 0 cend))
    (assert-float-equal
     (aref data rend 0)
     (linear-algebra:mref matrix rend 0))
    (assert-float-equal
     (aref data rend cend)
     (linear-algebra:mref matrix rend cend))
    (assert-float-equal
     (aref data rowi coli)
     (linear-algebra:mref matrix rowi coli))))

;;; Set square matrix elements
(define-test square-matrix-setf-mref
  (:tag :square-matrix :setf-mref)
  (let* ((rows 5) (columns 5)
         (rend (1- rows)) (cend (1- columns))
         (rowi (random-interior-index rows))
         (coli (random-interior-index columns))
         (matrix
          (linear-algebra:make-matrix
           rows columns
           :matrix-type 'linear-algebra:square-matrix
           :initial-contents
           '((1.1 1.2 1.3 1.4 1.5)
             (2.1 2.2 2.3 2.4 2.5)
             (3.1 3.2 3.3 3.4 3.5)
             (4.1 4.2 4.3 4.4 4.5)
             (5.1 5.2 5.3 5.4 5.5)))))
    (destructuring-bind (val1 val2 val3 val4 val5)
        (make-random-list 5 1.0)
      (setf (linear-algebra:mref matrix 0 0)    val1)
      (setf (linear-algebra:mref matrix 0 cend) val2)
      (setf (linear-algebra:mref matrix rend 0) val3)
      (setf (linear-algebra:mref matrix rend cend) val4)
      (setf (linear-algebra:mref matrix rowi coli) val5)
      (assert-float-equal
       val1 (linear-algebra:mref matrix 0 0))
      (assert-float-equal
       val2 (linear-algebra:mref matrix 0 cend))
      (assert-float-equal
       val3 (linear-algebra:mref matrix rend 0))
      (assert-float-equal
       val4 (linear-algebra:mref matrix rend cend))
      (assert-float-equal
       val5 (linear-algebra:mref matrix rowi coli)))))

;;; Copy the square matrix
(define-test copy-square-matrix
  (:tag :square-matrix :copy-matrix)
  (let ((matrix
         (linear-algebra:make-matrix
          5 5
          :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          (coordinate-array 0 0 5 5))))
    (assert-true
     (linear-algebra:square-matrix-p
      (linear-algebra:copy-matrix matrix)))
    (assert-false
     (eq matrix (linear-algebra:copy-matrix matrix)))
    (assert-false
     (eq
      (linear-algebra::contents matrix)
      (linear-algebra::contents
       (linear-algebra:copy-matrix matrix))))
    (assert-float-equal
     matrix (linear-algebra:copy-matrix matrix))))

;;; Test the submatrix of a square matrix
(define-test square-submatrix
  (:tag :square-matrix :submatrix)
  (let ((matrix
         (linear-algebra:make-matrix
          10 10
          :matrix-type
          'linear-algebra:square-matrix
          :initial-contents
          (coordinate-array))))
    ;; The entire matrix
    (assert-float-equal
     (coordinate-array)
     (linear-algebra:submatrix matrix 0 0))
    ;; Start row and column to the end
    (assert-float-equal
     (coordinate-array 3 3)
     (linear-algebra:submatrix matrix 3 3))
    ;; End row and column
    (assert-float-equal
     (coordinate-array 3 4 5 5)
     (linear-algebra:submatrix
      matrix 3 4 :end-row 5 :end-column 5))
    ;; Start row exceeds dimensions
    (assert-error
     'error (linear-algebra:submatrix matrix 11 5))
    ;; Start column exceeds dimensions
    (assert-error
     'error (linear-algebra:submatrix matrix 5 11))
    ;; End row exceeds dimensions
    (assert-error
     'error (linear-algebra:submatrix matrix 5 5 :end-row 11))
    ;; End column exceeds dimensions
    (assert-error
     'error (linear-algebra:submatrix matrix 5 5 :end-column 11))
    ;; Start row exceeds end row
    (assert-error
     'error (linear-algebra:submatrix matrix 7 7 :end-row 6))
    ;; Start column exceeds end column
    (assert-error
     'error (linear-algebra:submatrix matrix 7 7 :end-column 6))))

;;; Set the submatrix of a square matrix

(define-test setf-square-submatrix
  (:tag :square-matrix :setf-submatrix)
  ;; Upper left submatrix
  (let ((array-ul
         (make-array
          '(5 5) :initial-contents
          '((1 1 0 0 0)
            (1 1 0 0 0)
            (0 0 0 0 0)
            (0 0 0 0 0)
            (0 0 0 0 0)))))
    (assert-rational-equal
     array-ul
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 0 0)
      (unit-matrix 2 2)))
    (assert-rational-equal
     array-ul
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 0 0 :end-row 2 :end-column 2)
      (linear-algebra:submatrix (unit-matrix 5 5) 0 0)))
    (assert-rational-equal
     array-ul
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 0 0)
      (linear-algebra:submatrix
       (unit-matrix 5 5) 0 0 :end-row 2 :end-column 2)))
    (assert-rational-equal
     array-ul
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 0 0 :end-row 2 :end-column 2)
      (linear-algebra:submatrix
       (unit-matrix 5 5) 2 2 :end-row 4 :end-column 4))))
  ;; Upper right submatrix
  (let ((array-ur
         (make-array
          '(5 5) :initial-contents
          '((0 0 0 1 1)
            (0 0 0 1 1)
            (0 0 0 0 0)
            (0 0 0 0 0)
            (0 0 0 0 0)))))
    (assert-rational-equal
     array-ur
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 0 3)
      (unit-matrix 2 2)))
    (assert-rational-equal
     array-ur
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 0 3)
      (linear-algebra:submatrix
       (unit-matrix 5 5) 0 3 :end-row 2 :end-column 5)))
    (assert-rational-equal
     array-ur
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 0 3 :end-row 2 :end-column 5)
      (unit-matrix 5 5)))
    (assert-rational-equal
     array-ur
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 0 3 :end-row 2 :end-column 5)
      (linear-algebra:submatrix
       (unit-matrix 5 5) 2 2 :end-row 4 :end-column 4))))
  ;; Lower left submatrix
  (let ((array-ll
         (make-array
          '(5 5) :initial-contents
          '((0 0 0 0 0)
            (0 0 0 0 0)
            (0 0 0 0 0)
            (1 1 0 0 0)
            (1 1 0 0 0)))))
    (assert-rational-equal
     array-ll
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 3 0)
      (unit-matrix 2 2)))
    (assert-rational-equal
     array-ll
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 3 0)
      (linear-algebra:submatrix
       (unit-matrix 5 5) 0 3 :end-row 2 :end-column 5)))
    (assert-rational-equal
     array-ll
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 3 0 :end-row 5 :end-column 2)
      (unit-matrix 5 5)))
    (assert-rational-equal
     array-ll
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 3 0 :end-row 5 :end-column 2)
      (linear-algebra:submatrix
       (unit-matrix 5 5) 2 2 :end-row 4 :end-column 4))))
  ;; Lower right submatrix
  (let ((array-lr
         (make-array
          '(5 5)
          :initial-contents
          '((0 0 0 0 0)
            (0 0 0 0 0)
            (0 0 0 0 0)
            (0 0 0 1 1)
            (0 0 0 1 1)))))
    (assert-rational-equal
     array-lr
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 3 3)
      (unit-matrix 2 2)))
    (assert-rational-equal
     array-lr
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 3 3)
      (linear-algebra:submatrix
       (unit-matrix 5 5) 0 3 :end-row 2 :end-column 5)))
    (assert-rational-equal
     array-lr
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 3 3 :end-row 5 :end-column 5)
      (unit-matrix 5 5)))
    (assert-rational-equal
     array-lr
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 3 3 :end-row 5 :end-column 5)
      (linear-algebra:submatrix
       (unit-matrix 5 5) 2 2 :end-row 4 :end-column 4))))
  ;; Middle submatrix
  (let ((array-mid
         (make-array
          '(5 5)
          :initial-contents
          '((0 0 0 0 0)
            (0 1 1 1 0)
            (0 1 1 1 0)
            (0 1 1 1 0)
            (0 0 0 0 0)))))
    (assert-rational-equal
     array-mid
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 1 1)
      (unit-matrix 3 3)))
    (assert-rational-equal
     array-mid
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 1 1)
      (linear-algebra:submatrix
       (unit-matrix 5 5) 1 1 :end-row 4 :end-column 4)))
    (assert-rational-equal
     array-mid
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 1 1 :end-row 4 :end-column 4)
      (unit-matrix 5 5)))
    (assert-rational-equal
     array-mid
     (setf-submatrix
      5 5 'linear-algebra:square-matrix
      (linear-algebra:submatrix matrix 1 1 :end-row 4 :end-column 4)
      (linear-algebra:submatrix
       (unit-matrix 5 5) 1 1 :end-row 4 :end-column 4)))))

;;; Replace all or part of a square matrix
(define-test square-matrix-replace
  (:tag :square-matrix :replace-matrix)
  ;; Replace the entire matrix
  (assert-rational-equal
   (unit-matrix 5 5)
   (linear-algebra:replace-matrix
    (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
    (unit-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)))
  ;; Replace the first 2 rows
  (let ((result
         (make-array
          '(5 5)
          :initial-contents
          '((1 1 1 1 1)
            (1 1 1 1 1)
            (0 0 0 0 0)
            (0 0 0 0 0)
            (0 0 0 0 0)))))
    (assert-rational-equal
     result
     (linear-algebra:replace-matrix
      (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      (unit-matrix 2 5)))
    (assert-rational-equal
     result
     (linear-algebra:replace-matrix
      (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      (unit-matrix 2 7)))
    (assert-rational-equal
     result
     (linear-algebra:replace-matrix
      (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      (unit-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      :start-row2 3))
    (assert-rational-equal
     result
     (linear-algebra:replace-matrix
      (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      (unit-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      :end-row1 2))
    (assert-rational-equal
     result
     (linear-algebra:replace-matrix
      (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      (unit-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      :end-row2 2)))
  ;; Replace the first 3 columns
  (let ((result
         (make-array
          '(5 5)
          :initial-contents
          '((1 1 1 0 0)
            (1 1 1 0 0)
            (1 1 1 0 0)
            (1 1 1 0 0)
            (1 1 1 0 0)))))
    (assert-rational-equal
     result
     (linear-algebra:replace-matrix
      (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      (unit-matrix 5 3)))
    (assert-rational-equal
     result
     (linear-algebra:replace-matrix
      (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      (unit-matrix 7 3)))
    (assert-rational-equal
     result
     (linear-algebra:replace-matrix
      (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      (unit-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      :start-column2 2))
    (assert-rational-equal
     result
     (linear-algebra:replace-matrix
      (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      (unit-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      :end-column1 3))
    (assert-rational-equal
     result
     (linear-algebra:replace-matrix
      (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      (unit-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      :end-column2 3)))
  ;; Replace the center
  (let ((result
         (make-array
          '(5 5)
          :initial-contents
          '((0 0 0 0 0)
            (0 1 1 1 0)
            (0 1 1 1 0)
            (0 1 1 1 0)
            (0 0 0 0 0)))))
    (assert-rational-equal
     result
     (linear-algebra:replace-matrix
      (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      (unit-matrix 3 3 :matrix-type 'linear-algebra:square-matrix)
      :start-row1 1 :start-column1 1))
    (assert-rational-equal
     result
     (linear-algebra:replace-matrix
      (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      (unit-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      :start-row1 1 :start-column1 1
      :end-row1 4 :end-column1 4))
    (assert-rational-equal
     result
     (linear-algebra:replace-matrix
      (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      (unit-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      :start-row1 1 :start-column1 1
      :end-row2 3 :end-column2 3))
    (assert-rational-equal
     result
     (linear-algebra:replace-matrix
      (zero-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      (unit-matrix 5 5 :matrix-type 'linear-algebra:square-matrix)
      :start-row1 1 :start-column1 1
      :start-row2 2 :start-column2 2))))

;;; Validate a range for a square matrix.

(define-test square-matrix-validated-range
  (:tag :square-matrix :matrix-validated-range)
  (test-matrix-validated-range
   'linear-algebra:square-matrix 10 10))

;;; Square matrix fundamental operations

(define-test sumsq-square-matrix
  (:tag :square-matrix :sumsq)
  (multiple-value-bind (scale sumsq)
      (linear-algebra:sumsq
       (linear-algebra:make-matrix
        4 4 :matrix-type 'linear-algebra:square-matrix
        :initial-contents
        #2A((1.1 1.2 1.3 1.4)
            (2.1 2.2 2.3 2.4)
            (3.1 3.2 3.3 3.4)
            (4.1 4.2 4.3 4.4))))
    (assert-float-equal 4.4 scale)
    (assert-float-equal 7.293389 sumsq)))

(define-test sump-square-matrix
  (:tag :square-matrix :sump)
  (multiple-value-bind (scale sump)
      (linear-algebra:sump
       (linear-algebra:make-matrix
        4 4 :matrix-type 'linear-algebra:square-matrix
        :initial-contents
        #2A((1.1 1.2 1.3 1.4)
            (2.1 2.2 2.3 2.4)
            (3.1 3.2 3.3 3.4)
            (4.1 4.2 4.3 4.4)))
       3.5)
    (assert-float-equal 4.4 scale)
    (assert-float-equal 5.3832097 sump)))

(define-test norm-square-matrix
  (:tag :square-matrix :norm)
  (let ((matrix
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A((1.1 1.2 1.3 1.4)
              (2.1 2.2 2.3 2.4)
              (3.1 3.2 3.3 3.4)
              (4.1 4.2 4.3 4.4)))))
    (assert-float-equal
     11.6 (linear-algebra:norm matrix))
    (assert-float-equal
     11.6 (linear-algebra:norm matrix :measure 1))
    (assert-float-equal
     4.4 (linear-algebra:norm matrix :measure :max))
    (assert-float-equal
     11.882762 (linear-algebra:norm matrix :measure :frobenius))
    (assert-float-equal
     17.0 (linear-algebra:norm matrix :measure :infinity))
    (assert-error
     'error
     (linear-algebra:norm matrix :measure :unknown))))

(define-test transpose-square-matrix
  (:tag :square-matrix :transpose)
  (let ((matrix
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A((1.1 1.2 1.3 1.4)
              (2.1 2.2 2.3 2.4)
              (3.1 3.2 3.3 3.4)
              (4.1 4.2 4.3 4.4))))
        (transpose
         #2A((1.1 2.1 3.1 4.1)
             (1.2 2.2 3.2 4.2)
             (1.3 2.3 3.3 4.3)
             (1.4 2.4 3.4 4.4))))
    (assert-true
     (typep
      (linear-algebra:transpose matrix)
      'linear-algebra:square-matrix))
    (assert-float-equal
     transpose (linear-algebra:transpose matrix))))

(define-test ntranspose-square-matrix
  (:tag :square-matrix :ntranspose)
  (let ((matrix
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A((1.1 1.2 1.3 1.4)
              (2.1 2.2 2.3 2.4)
              (3.1 3.2 3.3 3.4)
              (4.1 4.2 4.3 4.4))))
        (transpose
         #2A((1.1 2.1 3.1 4.1)
             (1.2 2.2 3.2 4.2)
             (1.3 2.3 3.3 4.3)
             (1.4 2.4 3.4 4.4))))
    (assert-eq matrix (linear-algebra:ntranspose matrix))
    (assert-float-equal transpose matrix)))

(define-test permute-square-matrix
  (:tag :square-matrix :permute)
  (let ((matrix
         (linear-algebra:make-matrix
          5 5 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A((1.0 1.1 1.2 1.3 1.4)
              (2.0 2.1 2.2 2.3 2.4)
              (3.0 3.1 3.2 3.3 3.4)
              (4.0 4.1 4.2 4.3 4.4)
              (5.0 5.1 5.2 5.3 5.4))))
        (pmat
         (linear-algebra:make-matrix
          5 5 :matrix-type 'linear-algebra:permutation-matrix
          :initial-contents
          '((0 0 1 0 0)
            (0 0 0 0 1)
            (1 0 0 0 0)
            (0 1 0 0 0)
            (0 0 0 1 0)))))
    (assert-float-equal
     #2A((1.2 1.3 1.0 1.4 1.1)
         (2.2 2.3 2.0 2.4 2.1)
         (3.2 3.3 3.0 3.4 3.1)
         (4.2 4.3 4.0 4.4 4.1)
         (5.2 5.3 5.0 5.4 5.1))
     (linear-algebra:permute matrix pmat))
    (assert-float-equal
     #2A((3.0 3.1 3.2 3.3 3.4)
         (5.0 5.1 5.2 5.3 5.4)
         (1.0 1.1 1.2 1.3 1.4)
         (2.0 2.1 2.2 2.3 2.4)
         (4.0 4.1 4.2 4.3 4.4))
     (linear-algebra:permute pmat matrix))))

(define-test npermute-square-matrix
  (:tag :square-matrix :npermute)
  (let ((matrix
         (linear-algebra:make-matrix
          5 5 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          (make-array
           '(5 5) :initial-contents
           '((1.0 1.1 1.2 1.3 1.4)
             (2.0 2.1 2.2 2.3 2.4)
             (3.0 3.1 3.2 3.3 3.4)
             (4.0 4.1 4.2 4.3 4.4)
             (5.0 5.1 5.2 5.3 5.4)))))
        (pmat
         (linear-algebra:make-matrix
          5 5 :matrix-type 'linear-algebra:permutation-matrix
          :initial-contents
          '((0 0 1 0 0)
            (0 0 0 0 1)
            (1 0 0 0 0)
            (0 1 0 0 0)
            (0 0 0 1 0)))))
    (assert-eq matrix (linear-algebra:npermute matrix pmat))
    (assert-float-equal
     #2A((1.2 1.3 1.0 1.4 1.1)
         (2.2 2.3 2.0 2.4 2.1)
         (3.2 3.3 3.0 3.4 3.1)
         (4.2 4.3 4.0 4.4 4.1)
         (5.2 5.3 5.0 5.4 5.1))
     matrix))
  (let ((matrix
         (linear-algebra:make-matrix
          5 5 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          (make-array
           '(5 5) :initial-contents
           '((1.0 1.1 1.2 1.3 1.4)
             (2.0 2.1 2.2 2.3 2.4)
             (3.0 3.1 3.2 3.3 3.4)
             (4.0 4.1 4.2 4.3 4.4)
             (5.0 5.1 5.2 5.3 5.4)))))
        (pmat
         (linear-algebra:make-matrix
          5 5 :matrix-type 'linear-algebra:permutation-matrix
          :initial-contents
          '((0 0 1 0 0)
            (0 0 0 0 1)
            (1 0 0 0 0)
            (0 1 0 0 0)
            (0 0 0 1 0)))))
    (assert-eq matrix (linear-algebra:npermute pmat matrix))
    (assert-float-equal
     #2A((3.0 3.1 3.2 3.3 3.4)
         (5.0 5.1 5.2 5.3 5.4)
         (1.0 1.1 1.2 1.3 1.4)
         (2.0 2.1 2.2 2.3 2.4)
         (4.0 4.1 4.2 4.3 4.4))
     matrix)))

(define-test scale-square-matrix
  (:tag :square-matrix :scale)
  (assert-float-equal
   #2A(( 3.3  3.6  3.9  4.2)
       ( 6.3  6.6  6.9  7.2)
       ( 9.3  9.6  9.9 10.2)
       (12.3 12.6 12.9 13.2))
   (linear-algebra:scale
    3.0 (linear-algebra:make-matrix
         4 4 :matrix-type 'linear-algebra:square-matrix
         :initial-contents
         #2A((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4))))))

(define-test nscale-square-matrix
  (:tag :square-matrix :nscale)
  (let ((matrix
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A((1.1 1.2 1.3 1.4)
              (2.1 2.2 2.3 2.4)
              (3.1 3.2 3.3 3.4)
              (4.1 4.2 4.3 4.4)))))
    (assert-eq matrix (linear-algebra:nscale 3.0 matrix))
    (assert-float-equal
     #2A(( 3.3  3.6  3.9  4.2)
         ( 6.3  6.6  6.9  7.2)
         ( 9.3  9.6  9.9 10.2)
         (12.3 12.6 12.9 13.2))
     matrix)))

(define-test add-square-matrix
  (:tag :square-matrix :add)
  (let ((matrix
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A((1.1 1.2 1.3 1.4)
              (2.1 2.2 2.3 2.4)
              (3.1 3.2 3.3 3.4)
              (4.1 4.2 4.3 4.4)))))
    ;; No scalar
    (assert-float-equal
     #2A(( 2.2  2.4  2.6  2.8)
         ( 4.2  4.4  4.6  4.8)
         ( 6.2  6.4  6.6  6.8)
         ( 8.2  8.4  8.6  8.8))
     (linear-algebra:add matrix matrix))
    ;; Scalar1
    (assert-float-equal
     #2A(( 3.3  3.6  3.9  4.2)
         ( 6.3  6.6  6.9  7.2)
         ( 9.3  9.6  9.9 10.2)
         (12.3 12.6 12.9 13.2))
     (linear-algebra:add matrix matrix :scalar1 2.0))
    ;; Scalar2
    (assert-float-equal
     #2A(( 3.3  3.6  3.9  4.2)
         ( 6.3  6.6  6.9  7.2)
         ( 9.3  9.6  9.9 10.2)
         (12.3 12.6 12.9 13.2))
     (linear-algebra:add matrix matrix :scalar2 2.0))
    ;; Scalar1 & Scalar2
    (assert-float-equal
     #2A(( 5.5  6.0  6.5  7.0)
         (10.5 11.0 11.5 12.0)
         (15.5 16.0 16.5 17.0)
         (20.5 21.0 21.5 22.0))
     (linear-algebra:add matrix matrix :scalar1 2.0 :scalar2 3.0))))

(define-test nadd-square-matrix
  (:tag :square-matrix :nadd)
  ;; No scalar
  (let ((matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A((1.1 1.2 1.3 1.4)
              (2.1 2.2 2.3 2.4)
              (3.1 3.2 3.3 3.4)
              (4.1 4.2 4.3 4.4)))))
    (assert-eq matrix1 (linear-algebra:nadd matrix1 matrix2))
    (assert-float-equal
     #2A((2.2 2.4 2.6 2.8)
         (4.2 4.4 4.6 4.8)
         (6.2 6.4 6.6 6.8)
         (8.2 8.4 8.6 8.8))
     matrix1))
  ;; Scalar1
  (let ((matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A((1.1 1.2 1.3 1.4)
              (2.1 2.2 2.3 2.4)
              (3.1 3.2 3.3 3.4)
              (4.1 4.2 4.3 4.4)))))
    (assert-eq
     matrix1 (linear-algebra:nadd matrix1 matrix2 :scalar1 2.0))
    (assert-float-equal
     #2A(( 3.3  3.6  3.9  4.2)
         ( 6.3  6.6  6.9  7.2)
         ( 9.3  9.6  9.9 10.2)
         (12.3 12.6 12.9 13.2))
     matrix1))
  ;; Scalar2
  (let ((matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A((1.1 1.2 1.3 1.4)
              (2.1 2.2 2.3 2.4)
              (3.1 3.2 3.3 3.4)
              (4.1 4.2 4.3 4.4)))))
    (assert-eq
     matrix1 (linear-algebra:nadd matrix1 matrix2 :scalar2 2.0))
    (assert-float-equal
     #2A(( 3.3  3.6  3.9  4.2)
         ( 6.3  6.6  6.9  7.2)
         ( 9.3  9.6  9.9 10.2)
         (12.3 12.6 12.9 13.2))
     matrix1))
  ;; Scalar1 & Scalar2
  (let ((matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A((1.1 1.2 1.3 1.4)
              (2.1 2.2 2.3 2.4)
              (3.1 3.2 3.3 3.4)
              (4.1 4.2 4.3 4.4)))))
    (assert-eq
     matrix1
     (linear-algebra:nadd
      matrix1 matrix2 :scalar1 2.0 :scalar2 3.0))
    (assert-float-equal
     #2A(( 5.5  6.0  6.5  7.0)
         (10.5 11.0 11.5 12.0)
         (15.5 16.0 16.5 17.0)
         (20.5 21.0 21.5 22.0))
     matrix1)))

(define-test subtract-square-matrix
  (:tag :square-matrix :subtract)
  (let ((*epsilon* (* 3F0 single-float-epsilon))
        (matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A(( 2.2  2.4  2.6  2.8)
              ( 4.2  4.4  4.6  4.8)
              ( 6.2  6.4  6.6  6.8)
              ( 8.2  8.4  8.6  8.8))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A((1.1 1.2 1.3 1.4)
              (2.1 2.2 2.3 2.4)
              (3.1 3.2 3.3 3.4)
              (4.1 4.2 4.3 4.4)))))
    ;; No scalar
    (assert-float-equal
     #2A((1.1 1.2 1.3 1.4)
         (2.1 2.2 2.3 2.4)
         (3.1 3.2 3.3 3.4)
         (4.1 4.2 4.3 4.4))
     (linear-algebra:subtract matrix1 matrix2))
    ;; Scalar1
    (assert-float-equal
     #2A(( 3.3  3.6  3.9  4.2)
         ( 6.3  6.6  6.9  7.2)
         ( 9.3  9.6  9.9 10.2)
         (12.3 12.6 12.9 13.2))
     (linear-algebra:subtract matrix1 matrix2 :scalar1 2.0))
    ;; Scalar2
    (assert-float-equal
     #2A((0.0 0.0 0.0 0.0)
         (0.0 0.0 0.0 0.0)
         (0.0 0.0 0.0 0.0)
         (0.0 0.0 0.0 0.0))
     (linear-algebra:subtract matrix1 matrix2 :scalar2 2.0))
    ;; Scalar1 & Scalar2
    (assert-float-equal
     #2A((1.1 1.2 1.3 1.4)
         (2.1 2.2 2.3 2.4)
         (3.1 3.2 3.3 3.4)
         (4.1 4.2 4.3 4.4))
     (linear-algebra:subtract
      matrix1 matrix2 :scalar1 2.0 :scalar2 3.0))))

(define-test nsubtract-square-matrix
  (:tag :square-matrix :nsubtract)
  ;; No scalar
  (let ((matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '(( 2.2  2.4  2.6  2.8)
             ( 4.2  4.4  4.6  4.8)
             ( 6.2  6.4  6.6  6.8)
             ( 8.2  8.4  8.6  8.8)))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A((1.1 1.2 1.3 1.4)
              (2.1 2.2 2.3 2.4)
              (3.1 3.2 3.3 3.4)
              (4.1 4.2 4.3 4.4)))))
    (assert-eq matrix1 (linear-algebra:nsubtract matrix1 matrix2))
    (assert-float-equal
     #2A((1.1 1.2 1.3 1.4)
         (2.1 2.2 2.3 2.4)
         (3.1 3.2 3.3 3.4)
         (4.1 4.2 4.3 4.4))
     matrix1))
  ;; Scalar1
  (let ((matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :initial-contents
          #2A((1.1 1.2 1.3 1.4)
              (2.1 2.2 2.3 2.4)
              (3.1 3.2 3.3 3.4)
              (4.1 4.2 4.3 4.4)))))
    (assert-eq
     matrix1
     (linear-algebra:nsubtract matrix1 matrix2 :scalar1 2.0))
    (assert-float-equal
     #2A((1.1 1.2 1.3 1.4)
         (2.1 2.2 2.3 2.4)
         (3.1 3.2 3.3 3.4)
         (4.1 4.2 4.3 4.4))
     matrix1))
  ;; Scalar2
  (let ((*epsilon* (* 4F0 single-float-epsilon))
        (matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '(( 3.3  3.6  3.9  4.2)
             ( 6.3  6.6  6.9  7.2)
             ( 9.3  9.6  9.9 10.2)
             (12.3 12.6 12.9 13.2)))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A((1.1 1.2 1.3 1.4)
              (2.1 2.2 2.3 2.4)
              (3.1 3.2 3.3 3.4)
              (4.1 4.2 4.3 4.4)))))
    (assert-eq
     matrix1
     (linear-algebra:nsubtract matrix1 matrix2 :scalar2 2.0))
    (assert-float-equal
     #2A((1.1 1.2 1.3 1.4)
         (2.1 2.2 2.3 2.4)
         (3.1 3.2 3.3 3.4)
         (4.1 4.2 4.3 4.4))
     matrix1))
  ;; Scalar1 & Scalar2
  (let ((*epsilon* (* 3F0 single-float-epsilon))
        (matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '(( 2.2  2.4  2.6  2.8)
             ( 4.2  4.4  4.6  4.8)
             ( 6.2  6.4  6.6  6.8)
             ( 8.2  8.4  8.6  8.8)))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:square-matrix
          :initial-contents
          #2A((1.1 1.2 1.3 1.4)
              (2.1 2.2 2.3 2.4)
              (3.1 3.2 3.3 3.4)
              (4.1 4.2 4.3 4.4)))))
    (assert-eq
     matrix1
     (linear-algebra:nsubtract
      matrix1 matrix2 :scalar1 2.0 :scalar2 3.0))
    (assert-float-equal
     #2A((1.1 1.2 1.3 1.4)
         (2.1 2.2 2.3 2.4)
         (3.1 3.2 3.3 3.4)
         (4.1 4.2 4.3 4.4))
     matrix1)))

(define-test product-square-matrix
  (:tag :square-matrix :product)
  ;; Row vector - dense matrix
  (assert-true
   (typep
    (linear-algebra:product
     (linear-algebra:row-vector 1.0 2.0 3.0)
     (unit-matrix 3 3 :matrix-type 'linear-algebra:square-matrix))
    'linear-algebra:row-vector))
  (assert-float-equal
   #(6.0 12.0 18.0)
   (linear-algebra:product
    (linear-algebra:row-vector 1.0 2.0 3.0)
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:square-matrix
     :initial-contents
     #2A((1.0 2.0 3.0)
         (1.0 2.0 3.0)
         (1.0 2.0 3.0)))))
  (assert-float-equal
   #(12.599999 25.199999 37.8)
   (linear-algebra:product
    (linear-algebra:row-vector 1.0 2.0 3.0)
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:square-matrix
     :initial-contents
     #2A((1.0 2.0 3.0)
         (1.0 2.0 3.0)
         (1.0 2.0 3.0)))
    :scalar 2.1))
  (assert-error
   'error
   (linear-algebra:product
    (linear-algebra:row-vector 1.0 2.0 3.0 4.0 5.0 6.0)
    (linear-algebra:make-matrix
     3 3 :initial-element 1.0
     :matrix-type 'linear-algebra:square-matrix)))
  ;; Dense matrix - column vector
  (assert-true
   (typep
    (linear-algebra:product
     (unit-matrix 3 3 :matrix-type 'linear-algebra:square-matrix)
     (linear-algebra:column-vector 1.0 2.0 3.0))
    'linear-algebra:column-vector))
  (assert-float-equal
   #(6.0 12.0 18.0)
   (linear-algebra:product
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:square-matrix
     :initial-contents
     #2A((1.0 1.0 1.0)
         (2.0 2.0 2.0)
         (3.0 3.0 3.0)))
    (linear-algebra:column-vector 1.0 2.0 3.0)))
  (assert-float-equal
   #(12.599999 25.199999 37.8)
   (linear-algebra:product
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:square-matrix
     :initial-contents
     #2A((1.0 1.0 1.0)
         (2.0 2.0 2.0)
         (3.0 3.0 3.0)))
    (linear-algebra:column-vector 1.0 2.0 3.0)
    :scalar 2.1))
  (assert-error
   'error
   (linear-algebra:product
    (linear-algebra:make-matrix
     3 3 :initial-element 1.0
     :matrix-type 'linear-algebra:square-matrix)
    (linear-algebra:column-vector 1.0 2.0 3.0 4.0 5.0 6.0)))
  ;; Dense matrix - matrix
  (assert-true
   (typep
    (linear-algebra:product
     (unit-matrix 3 3 :matrix-type 'linear-algebra:square-matrix)
     (unit-matrix 3 3 :matrix-type 'linear-algebra:square-matrix))
    'linear-algebra:square-matrix))
  (assert-float-equal
   #2A(( 6.0  6.0  6.0)
       (12.0 12.0 12.0)
       (18.0 18.0 18.0))
   (linear-algebra:product
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:square-matrix
     :initial-contents
     #2A((1.0 1.0 1.0)
         (2.0 2.0 2.0)
         (3.0 3.0 3.0)))
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:square-matrix
     :initial-contents
     #2A((1.0 1.0 1.0)
         (2.0 2.0 2.0)
         (3.0 3.0 3.0)))))
  (assert-float-equal
   #2A((12.599999 12.599999 12.599999)
       (25.199999 25.199999 25.199999)
       (37.8      37.8      37.8))
   (linear-algebra:product
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:square-matrix
     :initial-contents
     #2A((1.0 1.0 1.0)
         (2.0 2.0 2.0)
         (3.0 3.0 3.0)))
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:square-matrix
     :initial-contents
     #2A((1.0 1.0 1.0)
         (2.0 2.0 2.0)
         (3.0 3.0 3.0)))
    :scalar 2.1))
  (assert-error
   'error
   (linear-algebra:product
    (unit-matrix 3 3 :matrix-type 'linear-algebra:square-matrix)
    (unit-matrix 4 4 :matrix-type 'linear-algebra:square-matrix))))
