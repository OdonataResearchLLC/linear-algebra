#|

  Linear Algebra in Common Lisp Unit Tests

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

(in-package :linear-algebra-test)

(defun hermitian-matrix (&optional (start 0) (end 10))
  (linear-algebra:make-matrix
   (- end start) (- end start)
   :matrix-type 'linear-algebra:hermitian-matrix
   :initial-contents (hermitian-array start end)))

(defun unit-hermitian-matrix (size)
  "Return a ROWSxCOLUMNS unit Hermitian matrix."
  (linear-algebra:make-matrix
   size size
   :matrix-type 'linear-algebra:hermitian-matrix
   :initial-contents
   (let ((init (make-array (list size size))))
     (dotimes (i0 size init)
       (setf (aref init i0 i0) #C(1 0))
       (dotimes (i1 i0)
         (setf
          (aref init i0 i1) #C(1 -1)
          (aref init i1 i0) #C(1  1)))))))

(define-test make-hermitian-matrix
  (:tag :hermitian-matrix :make-matrix)
  ;; A default Hermitian matrix
  (let ((matrix
         (linear-algebra:make-matrix
          10 10
          :matrix-type 'linear-algebra:hermitian-matrix)))
    (assert-true (linear-algebra:matrixp matrix))
    (assert-true (typep matrix 'linear-algebra:hermitian-matrix))
    (assert-rational-equal
     (make-array '(10 10) :initial-element 0)
     matrix))
  ;; Specify the Hermitian matrix element type
  (let* ((data
          '((#C(1.0  0.0) #C(1.0  2.0) #C(1.0 3.0))
            (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0 3.0))
            (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0 0.0))))
         (matrix
          (linear-algebra:make-matrix
           3 3 :matrix-type 'linear-algebra:hermitian-matrix
           :element-type '(complex single-float)
           :initial-contents data)))
    (assert-true (linear-algebra:matrixp matrix))
    (assert-true (typep matrix 'linear-algebra:hermitian-matrix))
    (assert-eq
     (array-element-type
      (linear-algebra::contents matrix))
     (array-element-type
      (make-array
       '(3 3) :element-type '(complex single-float)
       :initial-contents data)))
    (assert-float-equal
     (make-array '(3 3) :element-type '(complex single-float)
                 :initial-contents data)
     matrix))
  ;; Specify the Hermitian matrix contents - Nested list
  (let* ((data
          '((#C(1  0) #C(1  2) #C(1  3) #C(1 4))
            (#C(1 -2) #C(2  0) #C(2  3) #C(2 4))
            (#C(1 -3) #C(2 -3) #C(3  0) #C(3 4))
            (#C(1 -4) #C(2 -4) #C(3 -4) #C(4 0))))
         (matrix
          (linear-algebra:make-matrix
           4 4
           :matrix-type 'linear-algebra:hermitian-matrix
           :initial-contents data)))
    (assert-true (linear-algebra:matrixp matrix))
    (assert-true (typep matrix 'linear-algebra:hermitian-matrix))
    (assert-rational-equal
     (make-array '(4 4) :initial-contents data)
     matrix))
  ;; Specify the Hermitian matrix contents - Nested vector
  (let* ((data
          #(#(#C(1  0) #C(1  2) #C(1  3) #C(1 4))
            #(#C(1 -2) #C(2  0) #C(2  3) #C(2 4))
            #(#C(1 -3) #C(2 -3) #C(3  0) #C(3 4))
            #(#C(1 -4) #C(2 -4) #C(3 -4) #C(4 0))))
         (matrix
          (linear-algebra:make-matrix
           4 4
           :matrix-type 'linear-algebra:hermitian-matrix
           :initial-contents data)))
    (assert-true (linear-algebra:matrixp matrix))
    (assert-true (typep matrix 'linear-algebra:hermitian-matrix))
    (assert-rational-equal
     (make-array '(4 4) :initial-contents data)
     matrix))
  ;; Specify the Hermitian matrix contents - 2D array
  (let* ((data
          (make-array
           '(4 4) :initial-contents
           '((#C(1  0) #C(1  2) #C(1  3) #C(1 4))
             (#C(1 -2) #C(2  0) #C(2  3) #C(2 4))
             (#C(1 -3) #C(2 -3) #C(3  0) #C(3 4))
             (#C(1 -4) #C(2 -4) #C(3 -4) #C(4 0)))))
         (matrix
          (linear-algebra:make-matrix
           4 4
           :matrix-type 'linear-algebra:hermitian-matrix
           :initial-contents data)))
    (assert-true (linear-algebra:matrixp matrix))
    (assert-true (typep matrix 'linear-algebra:hermitian-matrix))
    (assert-rational-equal data matrix))
  ;; Errors
  (assert-error
   'error
   (linear-algebra:make-matrix
    4 4 :matrix-type 'linear-algebra:hermitian-matrix
    :initial-contents #C(1 2)))
  (assert-error
   'error
   (linear-algebra:make-matrix
    4 4 :matrix-type 'linear-algebra:hermitian-matrix
    :initial-contents
    #3A(((#C(1 0) #C(1 2)) (#C(2 1) #C(2 2)))
        ((#C(3 1) #C(3 2)) (#C(4 1) #C(4 2)))
        ((#C(5 1) #C(5 2)) (#C(6 1) #C(6 2))))))
  (assert-error
   'error
   (linear-algebra:make-matrix
    3 4 :matrix-type 'linear-algebra:hermitian-matrix
    :initial-contents
    (hermitian-array 0 4)))
  (assert-error
   'error
   (linear-algebra:make-matrix
    4 3 :matrix-type 'linear-algebra:hermitian-matrix
    :initial-contents
    (hermitian-array 0 4)))
  (assert-error
   'error
   (linear-algebra:make-matrix
    5 5 :matrix-type 'linear-algebra:symmetric-matrix
    :initial-contents
    (coordinate-array 0 0 5 5)))
  ;; Specify initial element and initial contents
  (assert-error
   'error
   (linear-algebra:make-matrix
    4 4 :matrix-type 'linear-algebra:hermitian-matrix
    :initial-element 1.1
    :initial-contents
    (hermitian-array 0 4))))

;;; Test the hermitian matrix predicate
(define-test hermitian-matrix-predicate
  (:tag :hermitian-matrix)
  (assert-true
   (linear-algebra:hermitian-matrix-p
    (linear-algebra:make-matrix
     10 10 :matrix-type 'linear-algebra:hermitian-matrix
     :initial-contents (hermitian-array))))
  (assert-false
   (linear-algebra:hermitian-matrix-p (make-array '(10 10)))))

;;; Test the hermitian matrix bounds
(define-test hermitian-matrix-in-bounds-p
  (:tag :hermitian-matrix :matrix-in-bounds-p)
  (test-matrix-in-bounds-p
   'linear-algebra:hermitian-matrix (hermitian-array)))

;;; Test the hermitian matrix element type
(define-test hermitian-matrix-element-type
  (:tag :hermitian-matrix :matrix-element-type)
  (let ((numeric-types
         '(integer fixnum
           short-float single-float double-float long-float)))
    (dolist (ntype numeric-types)
      (assert-true
       (subtypep
        (linear-algebra:matrix-element-type
         (linear-algebra:make-matrix
          2 2 :matrix-type 'linear-algebra:hermitian-matrix
          :element-type `(complex ,ntype)
          :initial-contents
          (list
           (list
            (complex (coerce 1 ntype) (coerce  0 ntype))
            (complex (coerce 1 ntype) (coerce  1 ntype)))
           (list
            (complex (coerce 1 ntype) (coerce -1 ntype))
            (complex (coerce 1 ntype) (coerce  0 ntype))))))
        (upgraded-array-element-type `(complex ,ntype)))))))

;;; Test the hermitian matrix dimensions
(define-test hermitian-matrix-dimensions
  (:tag :hermitian-matrix :matrix-dimensions)
  (assert-equal
   (list 10 10)
   (linear-algebra:matrix-dimensions
    (linear-algebra:make-matrix
     10 10 :matrix-type 'linear-algebra:hermitian-matrix
     :initial-contents (hermitian-array)))))

;;; Test the hermitian matrix row dimension
(define-test hermitian-matrix-row-dimension
  (:tag :hermitian-matrix :matrix-row-dimension)
  (assert-eq
   10
   (linear-algebra:matrix-row-dimension
    (linear-algebra:make-matrix
     10 10 :matrix-type 'linear-algebra:hermitian-matrix
     :initial-contents (hermitian-array)))))

;;; Test the hermitian matrix column dimension
(define-test hermitian-matrix-column-dimension
  (:tag :hermitian-matrix :matrix-column-dimension)
  (assert-eq
   10
   (linear-algebra:matrix-column-dimension
    (linear-algebra:make-matrix
     10 10 :matrix-type 'linear-algebra:hermitian-matrix
     :initial-contents (hermitian-array)))))

;;; Reference hermitian matrix elements
(define-test hermitian-matrix-mref
  (:tag :hermitian-matrix :mref)
  (let* ((initial-contents
          '((#C(1  0) #C(1  2) #C(1  3) #C(1  4) #C(1 5))
            (#C(1 -2) #C(2  0) #C(2  3) #C(2  4) #C(2 5))
            (#C(1 -3) #C(2 -3) #C(3  0) #C(3  4) #C(3 5))
            (#C(1 -4) #C(2 -4) #C(3 -4) #C(4  0) #C(4 5))
            (#C(1 -5) #C(2 -5) #C(3 -5) #C(4 -5) #C(5 0))))
         (rows 5) (columns 5)
         (rend (1- rows)) (cend (1- columns))
         (rowi (random-interior-index rows))
         (coli
          (do ((i0 (random-interior-index columns)
                   (random-interior-index columns)))
              ((/= i0 rowi) i0)))
         (data
          (make-array
           (list rows columns)
           :initial-contents
           initial-contents))
         (matrix
          (linear-algebra:make-matrix
           rows columns
           :matrix-type
           'linear-algebra:hermitian-matrix
           :initial-contents
           initial-contents)))
    (assert-rational-equal
     (aref data 0 0)
     (linear-algebra:mref matrix 0 0))
    (assert-rational-equal
     (aref data 0 cend)
     (linear-algebra:mref matrix 0 cend))
    (assert-rational-equal
     (aref data rend 0)
     (linear-algebra:mref matrix rend 0))
    (assert-rational-equal
     (linear-algebra:mref matrix 0 cend)
     (conjugate
      (linear-algebra:mref matrix rend 0)))
    (assert-rational-equal
     (aref data rend cend)
     (linear-algebra:mref matrix rend cend))
    (assert-rational-equal
     (aref data rowi coli)
     (linear-algebra:mref matrix rowi coli))
    (assert-rational-equal
     (linear-algebra:mref matrix rowi coli)
     (conjugate
      (linear-algebra:mref matrix coli rowi)))))

;;; Set hermitian matrix elements
(define-test hermitian-matrix-setf-mref
  (:tag :hermitian-matrix :setf-mref)
  (let* ((rows 5) (columns 5)
         (rend (1- rows)) (cend (1- columns))
         (rowi (random-interior-index rows))
         (coli
          (do ((i0 (random-interior-index columns)
                   (random-interior-index columns)))
              ((/= i0 rowi) i0)))
         (matrix
          (linear-algebra:make-matrix
           rows columns
           :matrix-type 'linear-algebra:hermitian-matrix
           :initial-contents
           '((#C(1  0) #C(1  2) #C(1  3) #C(1  4) #C(1 5))
             (#C(1 -2) #C(2  0) #C(2  3) #C(2  4) #C(2 5))
             (#C(1 -3) #C(2 -3) #C(3  0) #C(3  4) #C(3 5))
             (#C(1 -4) #C(2 -4) #C(3 -4) #C(4  0) #C(4 5))
             (#C(1 -5) #C(2 -5) #C(3 -5) #C(4 -5) #C(5 0))))))
    (multiple-value-bind (val1 val2 val3 val4)
        (values
         #C(6 0) (complex-random #C(5 5))
         #C(7 0) (complex-random #C(5 5)))
      (setf (linear-algebra:mref matrix 0 0)       val1)
      (setf (linear-algebra:mref matrix 0 cend)    val2)
      (setf (linear-algebra:mref matrix rend cend) val3)
      (setf (linear-algebra:mref matrix rowi coli) val4)
      (assert-rational-equal val1 (linear-algebra:mref matrix 0 0))
      (assert-rational-equal val2 (linear-algebra:mref matrix 0 cend))
      (assert-rational-equal
       val2 (conjugate (linear-algebra:mref matrix rend 0)))
      (assert-rational-equal
       val3 (linear-algebra:mref matrix rend cend))
      (assert-rational-equal
       val4 (linear-algebra:mref matrix rowi coli))
      (assert-rational-equal
       val4 (conjugate (linear-algebra:mref matrix coli rowi))))))

;;; Copy the Hermitian matrix
(define-test copy-hermitian-matrix
  (:tag :hermitian-matrix :copy-matrix)
  (let ((matrix
         (linear-algebra:make-matrix
          5 5
          :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          (hermitian-array 0 5))))
    (assert-true
     (linear-algebra:hermitian-matrix-p
      (linear-algebra:copy-matrix matrix)))
    (assert-false
     (eq matrix (linear-algebra:copy-matrix matrix)))
    (assert-false
     (eq
      (linear-algebra::contents matrix)
      (linear-algebra::contents
       (linear-algebra:copy-matrix matrix))))
    (assert-rational-equal
     matrix (linear-algebra:copy-matrix matrix))))

;;; Test the submatrix of a Hermitian matrix
(define-test hermitian-submatrix
  (:tag :hermitian-matrix :submatrix)
  (let ((matrix
         (linear-algebra:make-matrix
          10 10
          :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents (hermitian-array)))
        (submat
         (linear-algebra:make-matrix
          10 10
          :matrix-type 'linear-algebra:dense-matrix
          :initial-contents (hermitian-array))))
    ;; The entire matrix
    (assert-rational-equal
     (hermitian-array)
     (linear-algebra:submatrix matrix 0 0))
    ;; Start row and column to the end
    (assert-rational-equal
     (hermitian-array 3)
     (linear-algebra:submatrix matrix 3 3))
    ;; End row and column
    (assert-rational-equal
     (hermitian-array 3 5)
     (linear-algebra:submatrix
      matrix 3 3 :end-row 5 :end-column 5))
    ;; Submatrix is a general matrix
    (assert-true
     (typep
      (linear-algebra:submatrix matrix 1 2)
      'linear-algebra:dense-matrix))
    (assert-rational-equal
     (linear-algebra:submatrix submat 1 2)
     (linear-algebra:submatrix matrix 1 2))
    (assert-true
     (typep
      (linear-algebra:submatrix matrix 1 1 :end-row 5)
      'linear-algebra:dense-matrix))
    (assert-rational-equal
     (linear-algebra:submatrix submat 1 1 :end-row 5)
     (linear-algebra:submatrix matrix 1 1 :end-row 5))
    (assert-true
     (typep
      (linear-algebra:submatrix matrix 1 1 :end-column 8)
      'linear-algebra:dense-matrix))
    (assert-rational-equal
     (linear-algebra:submatrix submat 1 1 :end-column 8)
     (linear-algebra:submatrix matrix 1 1 :end-column 8))
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

;;; Set the submatrix of a Hermitian matrix
(define-test setf-hermitian-submatrix
  (:tag :hermitian-matrix :setf-submatrix)
  (macrolet ((setf-submatrix (size submatrix-form data-form)
               (let ((matrix (second submatrix-form)))
                 `(let ((,matrix (unit-hermitian-matrix ,size)))
                   (setf ,submatrix-form ,data-form)
                   ,matrix))))
    ;; Upper left submatrix
    (let ((array-ul
           (make-array
            '(5 5) :initial-contents
            '((#C(1  0) #C(2  1) #C(3  1) #C(1  1) #C(1 1))
              (#C(2 -1) #C(2  0) #C(3  2) #C(1  1) #C(1 1))
              (#C(3 -1) #C(3 -2) #C(3  0) #C(1  1) #C(1 1))
              (#C(1 -1) #C(1 -1) #C(1 -1) #C(1  0) #C(1 1))
              (#C(1 -1) #C(1 -1) #C(1 -1) #C(1 -1) #C(1 0))))))
      (assert-rational-equal
       array-ul
       (setf-submatrix
        5
        (linear-algebra:submatrix matrix 0 0)
        (hermitian-matrix 0 3)))
      (assert-rational-equal
       array-ul
       (setf-submatrix
        5
        (linear-algebra:submatrix
         matrix 0 0 :end-row 3 :end-column 3)
        (hermitian-matrix))))
    ;; Lower right submatrix
    (assert-rational-equal
     (make-array
      '(5 5) :initial-contents
      '((#C(1  0) #C(1  1) #C(1  1) #C(1  1) #C(1 1))
        (#C(1 -1) #C(1  0) #C(1  1) #C(1  1) #C(1 1))
        (#C(1 -1) #C(1 -1) #C(1  0) #C(2  1) #C(3 1))
        (#C(1 -1) #C(1 -1) #C(2 -1) #C(2  0) #C(3 2))
        (#C(1 -1) #C(1 -1) #C(3 -1) #C(3 -2) #C(3 0))))
     (setf-submatrix
      5
      (linear-algebra:submatrix matrix 2 2)
      (hermitian-matrix)))
    ;; Middle submatrix
    (let ((array-mid
           (make-array
            '(5 5) :initial-contents
            '((#C(1  0) #C(1  1) #C(1  1) #C(1  1) #C(1 1))
              (#C(1 -1) #C(2  0) #C(3  2) #C(4  2) #C(1 1))
              (#C(1 -1) #C(3 -2) #C(3  0) #C(4  3) #C(1 1))
              (#C(1 -1) #C(4 -2) #C(4 -3) #C(4  0) #C(1 1))
              (#C(1 -1) #C(1 -1) #C(1 -1) #C(1 -1) #C(1 0))))))
      (assert-rational-equal
       array-mid
       (setf-submatrix
        5
        (linear-algebra:submatrix matrix 1 1)
        (hermitian-matrix 1 4)))
      (assert-rational-equal
       array-mid
       (setf-submatrix
        5
        (linear-algebra:submatrix
         matrix 1 1 :end-row 4 :end-column 4)
        (hermitian-matrix 1))))
    ;; Off diagonal submatrix
    (let ((array-off
           (make-array
            '(5 5) :initial-contents
            '((#C(1  0) #C(1  1) #C(1  1) #C(1  1) #C(1 2))
              (#C(1 -1) #C(1  0) #C(1  1) #C(2  1) #C(2 2))
              (#C(1 -1) #C(1 -1) #C(1  0) #C(1  1) #C(1 1))
              (#C(1 -1) #C(2 -1) #C(1 -1) #C(1  0) #C(1 1))
              (#C(1 -2) #C(2 -2) #C(1 -1) #C(1 -1) #C(1 0)))))
          (submatrix
           (linear-algebra:make-matrix
            3 3 :initial-contents
            '((#C(1 1) #C(1 2) #C(1 3))
              (#C(2 1) #C(2 2) #C(2 3))
              (#C(3 1) #C(3 2) #C(3 3))))))
      (assert-rational-equal
       array-off
       (setf-submatrix
        5 (linear-algebra:submatrix matrix 0 3)
        (linear-algebra:submatrix
         submatrix 0 0 :end-row 2 :end-column 2)))
      (assert-rational-equal
       array-off
       (setf-submatrix
        5 (linear-algebra:submatrix matrix 0 3 :end-row 2)
        submatrix)))
    (let ((array-off
           (make-array
            '(5 5) :initial-contents
            '((#C(1  0) #C(1  1) #C(1  1) #C(1  1) #C(2 1))
              (#C(1 -1) #C(1  0) #C(1  1) #C(1  2) #C(2 2))
              (#C(1 -1) #C(1 -1) #C(1  0) #C(1  1) #C(1 1))
              (#C(1 -1) #C(1 -2) #C(1 -1) #C(1  0) #C(1 1))
              (#C(2 -1) #C(2 -2) #C(1 -1) #C(1 -1) #C(1 0)))))
          (submatrix
           (linear-algebra:make-matrix
            3 3 :initial-contents
            '((#C(1 -1) #C(1 -2) #C(1 -3))
              (#C(2 -1) #C(2 -2) #C(2 -3))
              (#C(3 -1) #C(3 -2) #C(3 -3))))))
      (assert-rational-equal
       array-off
       (setf-submatrix
        5
        (linear-algebra:submatrix matrix 3 0 :end-column 2)
        submatrix))
      (assert-rational-equal
       array-off
       (setf-submatrix
        5
        (linear-algebra:submatrix matrix 3 0)
        (linear-algebra:submatrix submatrix 0 0 :end-column 2)))))
  ;; Non-Hermitian subsets
  (assert-error
   'error
   (setf
    (linear-algebra:submatrix
     (unit-hermitian-matrix 10) 0 1)
    (unit-matrix 5 3))))

;;; Replace all or part of a Hermitian matrix
(define-test hermitian-matrix-replace
  (:tag :hermitian-matrix :replace-matrix)
  ;; Replace the entire matrix
  (assert-rational-equal
   (hermitian-matrix)
   (linear-algebra:replace-matrix
    (unit-hermitian-matrix 10)
    (hermitian-matrix)))
  ;; Upper left submatrix
  (let ((array-ul
         (make-array
          '(5 5) :initial-contents
          '((#C(1  0) #C(2  1) #C(3  1) #C(1  1) #C(1 1))
            (#C(2 -1) #C(2  0) #C(3  2) #C(1  1) #C(1 1))
            (#C(3 -1) #C(3 -2) #C(3  0) #C(1  1) #C(1 1))
            (#C(1 -1) #C(1 -1) #C(1 -1) #C(1  0) #C(1 1))
            (#C(1 -1) #C(1 -1) #C(1 -1) #C(1 -1) #C(1 0))))))
    (assert-rational-equal
     array-ul
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5)
      (hermitian-matrix 0 3)))
    (assert-rational-equal
     array-ul
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5)
      (hermitian-matrix)
      :end-row1 3 :end-column1 3))
    (assert-rational-equal
     array-ul
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5)
      (hermitian-matrix)
      :end-row2 3 :end-column1 3))
    (assert-rational-equal
     array-ul
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5)
      (hermitian-matrix)
      :end-row1 3 :end-column2 3))
    (assert-rational-equal
     array-ul
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5)
      (hermitian-matrix)
      :end-row2 3 :end-column2 3)))
  ;; Lower right submatrix
  (assert-rational-equal
   (make-array
    '(5 5) :initial-contents
    '((#C(1  0) #C(1  1) #C(1  1) #C(1  1) #C(1 1))
      (#C(1 -1) #C(1  0) #C(1  1) #C(1  1) #C(1 1))
      (#C(1 -1) #C(1 -1) #C(1  0) #C(2  1) #C(3 1))
      (#C(1 -1) #C(1 -1) #C(2 -1) #C(2  0) #C(3 2))
      (#C(1 -1) #C(1 -1) #C(3 -1) #C(3 -2) #C(3 0))))
   (linear-algebra:replace-matrix
    (unit-hermitian-matrix 5)
    (hermitian-matrix)
    :start-row1 2 :start-column1 2))
  ;; Middle submatrix
  (let ((array-mid
         (make-array
          '(5 5) :initial-contents
          '((#C(1  0) #C(1  1) #C(1  1) #C(1  1) #C(1 1))
            (#C(1 -1) #C(2  0) #C(3  2) #C(4  2) #C(1 1))
            (#C(1 -1) #C(3 -2) #C(3  0) #C(4  3) #C(1 1))
            (#C(1 -1) #C(4 -2) #C(4 -3) #C(4  0) #C(1 1))
            (#C(1 -1) #C(1 -1) #C(1 -1) #C(1 -1) #C(1 0))))))
    (assert-rational-equal
     array-mid
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5)
      (hermitian-matrix 1 4)
      :start-row1 1 :start-column1 1))
    (assert-rational-equal
     array-mid
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5)
      (hermitian-matrix 1 4)
      :start-row1 1 :start-column1 1
      :end-row1 4 :end-column1 4))
    (assert-rational-equal
     array-mid
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5)
      (hermitian-matrix 1)
      :start-row1 1 :start-column1 1
      :end-row2 3 :end-column1 4))
    (assert-rational-equal
     array-mid
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5)
      (hermitian-matrix 1)
      :start-row1 1 :start-column1 1
      :end-row1 4 :end-column2 3))
    (assert-rational-equal
     array-mid
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5)
      (hermitian-matrix 1)
      :start-row1 1 :start-column1 1
      :end-row2 3 :end-column2 3)))
  ;; Off diagonal submatrix
  (let ((array-off
         (make-array
          '(5 5) :initial-contents
          '((#C(1  0) #C(1  1) #C(1  1) #C(1  1) #C(1 2))
            (#C(1 -1) #C(1  0) #C(1  1) #C(2  1) #C(2 2))
            (#C(1 -1) #C(1 -1) #C(1  0) #C(1  1) #C(1 1))
            (#C(1 -1) #C(2 -1) #C(1 -1) #C(1  0) #C(1 1))
            (#C(1 -2) #C(2 -2) #C(1 -1) #C(1 -1) #C(1 0)))))
        (submatrix1
         (linear-algebra:make-matrix
          2 2 :initial-contents
          '((#C(1 1) #C(1 2))
            (#C(2 1) #C(2 2)))))
        (submatrix2
         (linear-algebra:make-matrix
          4 4 :initial-contents
          '((#C(1 1) #C(1 2) #C(1 3) #C(1 4))
            (#C(2 1) #C(2 2) #C(2 3) #C(2 4))
            (#C(3 1) #C(3 2) #C(3 3) #C(3 4))
            (#C(4 1) #C(4 2) #C(4 3) #C(4 4))))))
    (assert-rational-equal
     array-off
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5) submatrix1
      :start-row1 0 :start-column1 3))
    (assert-rational-equal
     array-off
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5) submatrix2
      :start-row1 0 :start-column1 3 :end-row1 2))
    (assert-rational-equal
     array-off
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5) submatrix2
      :start-row1 0 :start-column1 3 :end-row2 2)))
  (let ((array-off
         (make-array
          '(5 5) :initial-contents
          '((#C(1  0) #C(1  1) #C(1  1) #C(1  1) #C(1 2))
            (#C(1 -1) #C(1  0) #C(1  1) #C(2  1) #C(2 2))
            (#C(1 -1) #C(1 -1) #C(1  0) #C(1  1) #C(1 1))
            (#C(1 -1) #C(2 -1) #C(1 -1) #C(1  0) #C(1 1))
            (#C(1 -2) #C(2 -2) #C(1 -1) #C(1 -1) #C(1 0)))))
        (submatrix
         (linear-algebra:make-matrix
          3 3 :initial-contents
          '((#C(1 1) #C(1 2) #C(1 3))
            (#C(2 1) #C(2 2) #C(2 3))
            (#C(3 1) #C(3 2) #C(3 3))))))
    (assert-rational-equal
     array-off
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5) submatrix
      :start-column1 3 :end-row1 2))
    (assert-rational-equal
     array-off
     (linear-algebra:replace-matrix
      (unit-hermitian-matrix 5) submatrix
      :start-column1 3 :end-row2 2)))
  ;; Non-Hermitian subsets
  (assert-error
   'error
   (linear-algebra:replace-matrix
    (unit-hermitian-matrix 5)
    (unit-matrix 5 3)
    :start-column1 1)))

;;; Validate a range for a hermitian matrix.
(define-test hermitian-matrix-validated-range
  (:tag :hermitian-matrix :matrix-validated-range)
  (let ((matrix (unit-hermitian-matrix 10))
        (row1 (random 10))
        (row2 (random 10))
        (col1 (random 10))
        (col2 (random 10)))
    (assert-equal
     (values row1 col1 10 10)
     (linear-algebra:matrix-validated-range matrix row1 col1))
    (assert-equal
     (values (min row1 row2) col1 (max row1 row2) 10)
     (linear-algebra:matrix-validated-range
      matrix (min row1 row2) col1 (max row1 row2)))
    (assert-equal
     (values row1 (min col1 col2) 10 (max col1 col2))
     (linear-algebra:matrix-validated-range
      matrix row1 (min col1 col2) nil (max col1 col2)))
    (assert-equal
     (values
      (min row1 row2) (min col1 col2)
      (max row1 row2) (max col1 col2))
     (linear-algebra:matrix-validated-range
      matrix
      (min row1 row2) (min col1 col2)
      (max row1 row2) (max col1 col2)))
    (assert-error
     'error
     (linear-algebra:matrix-validated-range matrix 11 col1))
    (assert-error
     'error
     (linear-algebra:matrix-validated-range matrix row1 11))
    (assert-error
     'error
     (linear-algebra:matrix-validated-range
      matrix 9 col1 1))
    (assert-error
     'error
     (linear-algebra:matrix-validated-range
      matrix row1 9 10 1))
    (assert-error
     'error
     (linear-algebra:matrix-validated-range
      matrix 9 9 1 1))))

(define-test norm-hermitian-matrix
  (:tag :hermitian-matrix :norm)
  (let ((matrix
         (linear-algebra:make-matrix
          5 5 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          '((#C(1  0) #C(1  2) #C(1  3) #C(1  4) #C(1 5))
            (#C(1 -2) #C(2  0) #C(2  3) #C(2  4) #C(2 5))
            (#C(1 -3) #C(2 -3) #C(3  0) #C(3  4) #C(3 5))
            (#C(1 -4) #C(2 -4) #C(3 -4) #C(4  0) #C(4 5))
            (#C(1 -5) #C(2 -5) #C(3 -5) #C(4 -5) #C(5 0))))))
    (assert-float-equal
     27.71826 (linear-algebra:norm matrix))
    (assert-float-equal
     27.71826 (linear-algebra:norm matrix 1))
    (assert-float-equal
     6.4031243 (linear-algebra:norm matrix :max))
    (assert-float-equal
     22.248597 (linear-algebra:norm matrix :frobenius))
    (assert-float-equal
     27.71826 (linear-algebra:norm matrix :infinity))
    (assert-error
     'error
     (linear-algebra:norm matrix :unknown))))

(define-test transpose-hermitian-matrix
  (:tag :hermitian-matrix :transpose)
  (let ((matrix
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          '((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
            (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
            (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
            (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0)))))
        (transpose
         #2A((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
             (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
             (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
             (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0)))))
    (assert-true
     (typep
      (linear-algebra:transpose matrix)
      'linear-algebra:hermitian-matrix))
    (assert-float-equal
     transpose (linear-algebra:transpose matrix))))

(define-test ntranspose-hermitian-matrix
  (:tag :hermitian-matrix :ntranspose)
  (let ((matrix
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          '((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
            (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
            (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
            (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0)))))
        (transpose
         #2A((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
             (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
             (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
             (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0)))))
    (assert-eq matrix (linear-algebra:ntranspose matrix))
    (assert-float-equal transpose matrix)))

(define-test permute-hermitian-matrix
  (:tag :hermitian-matrix :permute)
  (let ((matrix
         (linear-algebra:make-matrix
          5 5 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          '((#C(1  0) #C(1  2) #C(1  3) #C(1  4) #C(1 5))
            (#C(1 -2) #C(2  0) #C(2  3) #C(2  4) #C(2 5))
            (#C(1 -3) #C(2 -3) #C(3  0) #C(3  4) #C(3 5))
            (#C(1 -4) #C(2 -4) #C(3 -4) #C(4  0) #C(4 5))
            (#C(1 -5) #C(2 -5) #C(3 -5) #C(4 -5) #C(5 0)))))
        (pmat
         (linear-algebra:make-matrix
          5 5 :matrix-type 'linear-algebra:permutation-matrix
          :initial-contents
          '((0 0 1 0 0)
            (0 0 0 0 1)
            (1 0 0 0 0)
            (0 1 0 0 0)
            (0 0 0 1 0)))))
    (assert-true
     (typep
      (linear-algebra:permute matrix pmat)
      'linear-algebra:square-matrix))
    (assert-rational-equal
     #2A((#C(1  3) #C(1  4) #C(1  0) #C(1  5) #C(1  2))
         (#C(2  3) #C(2  4) #C(1 -2) #C(2  5) #C(2  0))
         (#C(3  0) #C(3  4) #C(1 -3) #C(3  5) #C(2 -3))
         (#C(3 -4) #C(4  0) #C(1 -4) #C(4  5) #C(2 -4))
         (#C(3 -5) #C(4 -5) #C(1 -5) #C(5  0) #C(2 -5)))
     (linear-algebra:permute matrix pmat))
    (assert-true
     (typep
      (linear-algebra:permute pmat matrix)
      'linear-algebra:square-matrix))
    (assert-rational-equal
     #2A((#C(1 -3) #C(2 -3) #C(3  0) #C(3  4) #C(3  5))
         (#C(1 -5) #C(2 -5) #C(3 -5) #C(4 -5) #C(5  0))
         (#C(1  0) #C(1  2) #C(1  3) #C(1  4) #C(1  5))
         (#C(1 -2) #C(2  0) #C(2  3) #C(2  4) #C(2  5))
         (#C(1 -4) #C(2 -4) #C(3 -4) #C(4  0) #C(4  5)))
     (linear-algebra:permute pmat matrix))))

(define-test scale-hermitian-matrix
  (:tag :hermitian-matrix :scale)
  (assert-float-equal
   #2A((#C(3.0   0.0) #C(3.0   6.0) #C(3.0   9.0) #C( 3.0 12.0))
       (#C(3.0  -6.0) #C(6.0   0.0) #C(6.0   9.0) #C( 6.0 12.0))
       (#C(3.0  -9.0) #C(6.0  -9.0) #C(9.0   0.0) #C( 9.0 12.0))
       (#C(3.0 -12.0) #C(6.0 -12.0) #C(9.0 -12.0) #C(12.0  0.0)))
   (linear-algebra:scale
    3.0
    (linear-algebra:make-matrix
     4 4 :matrix-type 'linear-algebra:hermitian-matrix
     :initial-contents
     '((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
       (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
       (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
       (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0)))))))

(define-test nscale-hermitian-matrix
  (:tag :hermitian-matrix :nscale)
  (let ((matrix
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          '((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
            (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
            (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
            (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0))))))
    (assert-eq matrix (linear-algebra:nscale 3.0 matrix))
    (assert-float-equal
     #2A((#C(3.0   0.0) #C(3.0   6.0) #C(3.0   9.0) #C( 3.0 12.0))
         (#C(3.0  -6.0) #C(6.0   0.0) #C(6.0   9.0) #C( 6.0 12.0))
         (#C(3.0  -9.0) #C(6.0  -9.0) #C(9.0   0.0) #C( 9.0 12.0))
         (#C(3.0 -12.0) #C(6.0 -12.0) #C(9.0 -12.0) #C(12.0  0.0)))
     matrix)))

(define-test add-hermitian-matrix
  (:tag :hermitian-matrix :add)
  (let ((matrix
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          '((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
            (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
            (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
            (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0))))))
    ;; No scalar
    (assert-float-equal
     #2A((#C(2.0  0.0) #C(2.0  4.0) #C(2.0  6.0) #C(2.0 8.0))
         (#C(2.0 -4.0) #C(4.0  0.0) #C(4.0  6.0) #C(4.0 8.0))
         (#C(2.0 -6.0) #C(4.0 -6.0) #C(6.0  0.0) #C(6.0 8.0))
         (#C(2.0 -8.0) #C(4.0 -8.0) #C(6.0 -8.0) #C(8.0 0.0)))
     (linear-algebra:add matrix matrix))
    ;; Scalar1
    (assert-float-equal
     #2A((#C(3.0   0.0) #C(3.0   6.0) #C(3.0   9.0) #C( 3.0 12.0))
         (#C(3.0  -6.0) #C(6.0   0.0) #C(6.0   9.0) #C( 6.0 12.0))
         (#C(3.0  -9.0) #C(6.0  -9.0) #C(9.0   0.0) #C( 9.0 12.0))
         (#C(3.0 -12.0) #C(6.0 -12.0) #C(9.0 -12.0) #C(12.0  0.0)))
     (linear-algebra:add matrix matrix :scalar1 2.0))
    ;; Scalar2
    (assert-float-equal
     #2A((#C(3.0   0.0) #C(3.0   6.0) #C(3.0   9.0) #C( 3.0 12.0))
         (#C(3.0  -6.0) #C(6.0   0.0) #C(6.0   9.0) #C( 6.0 12.0))
         (#C(3.0  -9.0) #C(6.0  -9.0) #C(9.0   0.0) #C( 9.0 12.0))
         (#C(3.0 -12.0) #C(6.0 -12.0) #C(9.0 -12.0) #C(12.0  0.0)))
     (linear-algebra:add matrix matrix :scalar2 2.0))
    ;; Scalar1 & Scalar2
    (assert-float-equal
     #2A((#C(5.0   0.0) #C( 5.0  10.0) #C( 5.0  15.0) #C( 5.0 20.0))
         (#C(5.0 -10.0) #C(10.0   0.0) #C(10.0  15.0) #C(10.0 20.0))
         (#C(5.0 -15.0) #C(10.0 -15.0) #C(15.0   0.0) #C(15.0 20.0))
         (#C(5.0 -20.0) #C(10.0 -20.0) #C(15.0 -20.0) #C(20.0  0.0)))
     (linear-algebra:add matrix matrix :scalar1 2.0 :scalar2 3.0))))

(define-test nadd-hermitian-matrix
  (:tag :hermitian-matrix :nadd)
  ;; No scalar
  (let ((matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
             (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
             (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
             (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0))))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          '((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
            (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
            (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
            (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0))))))
    (assert-eq matrix1 (linear-algebra:nadd matrix1 matrix2))
    (assert-float-equal
     #2A((#C(2.0  0.0) #C(2.0  4.0) #C(2.0  6.0) #C(2.0 8.0))
         (#C(2.0 -4.0) #C(4.0  0.0) #C(4.0  6.0) #C(4.0 8.0))
         (#C(2.0 -6.0) #C(4.0 -6.0) #C(6.0  0.0) #C(6.0 8.0))
         (#C(2.0 -8.0) #C(4.0 -8.0) #C(6.0 -8.0) #C(8.0 0.0)))
     matrix1))
  ;; Scalar1
  (let ((matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
             (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
             (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
             (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0))))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          '((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
            (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
            (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
            (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0))))))
    (assert-eq
     matrix1 (linear-algebra:nadd matrix1 matrix2 :scalar1 2.0))
    (assert-float-equal
     #2A((#C(3.0   0.0) #C(3.0   6.0) #C(3.0   9.0) #C( 3.0 12.0))
         (#C(3.0  -6.0) #C(6.0   0.0) #C(6.0   9.0) #C( 6.0 12.0))
         (#C(3.0  -9.0) #C(6.0  -9.0) #C(9.0   0.0) #C( 9.0 12.0))
         (#C(3.0 -12.0) #C(6.0 -12.0) #C(9.0 -12.0) #C(12.0  0.0)))
     matrix1))
  ;; Scalar2
  (let ((matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
             (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
             (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
             (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0))))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          '((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
            (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
            (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
            (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0))))))
    (assert-eq
     matrix1 (linear-algebra:nadd matrix1 matrix2 :scalar2 2.0))
    (assert-float-equal
     #2A((#C(3.0   0.0) #C(3.0   6.0) #C(3.0   9.0) #C( 3.0 12.0))
         (#C(3.0  -6.0) #C(6.0   0.0) #C(6.0   9.0) #C( 6.0 12.0))
         (#C(3.0  -9.0) #C(6.0  -9.0) #C(9.0   0.0) #C( 9.0 12.0))
         (#C(3.0 -12.0) #C(6.0 -12.0) #C(9.0 -12.0) #C(12.0  0.0)))
     matrix1))
  ;; Scalar1 & Scalar2
  (let ((matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
             (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
             (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
             (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0))))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          '((#C(1.0  0.0) #C(1.0  2.0) #C(1.0  3.0) #C(1.0 4.0))
            (#C(1.0 -2.0) #C(2.0  0.0) #C(2.0  3.0) #C(2.0 4.0))
            (#C(1.0 -3.0) #C(2.0 -3.0) #C(3.0  0.0) #C(3.0 4.0))
            (#C(1.0 -4.0) #C(2.0 -4.0) #C(3.0 -4.0) #C(4.0 0.0))))))
    (assert-eq
     matrix1
     (linear-algebra:nadd
      matrix1 matrix2 :scalar1 2.0 :scalar2 3.0))
    (assert-float-equal
     #2A((#C(5.0   0.0) #C( 5.0  10.0) #C( 5.0  15.0) #C( 5.0 20.0))
         (#C(5.0 -10.0) #C(10.0   0.0) #C(10.0  15.0) #C(10.0 20.0))
         (#C(5.0 -15.0) #C(10.0 -15.0) #C(15.0   0.0) #C(15.0 20.0))
         (#C(5.0 -20.0) #C(10.0 -20.0) #C(15.0 -20.0) #C(20.0  0.0)))
     matrix1)))

(define-test subtract-hermitian-matrix
  (:tag :hermitian-matrix :subtract)
  (let ((*epsilon* (* 3F0 single-float-epsilon))
        (matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          #2A((#C(2.0  0.0) #C(4.0  2.0) #C(6.0  2.0) #C(8.0 2.0))
              (#C(4.0 -2.0) #C(4.0  0.0) #C(6.0  4.0) #C(8.0 4.0))
              (#C(6.0 -2.0) #C(6.0 -4.0) #C(6.0  0.0) #C(8.0 6.0))
              (#C(8.0 -2.0) #C(8.0 -4.0) #C(8.0 -6.0) #C(8.0 0.0)))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0) #C(4.0 1.0))
              (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0) #C(4.0 2.0))
              (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0) #C(4.0 3.0))
              (#C(4.0 -1.0) #C(4.0 -2.0) #C(4.0 -3.0) #C(4.0 0.0))))))
    ;; No scalar
    (assert-float-equal
     #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0) #C(4.0 1.0))
         (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0) #C(4.0 2.0))
         (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0) #C(4.0 3.0))
         (#C(4.0 -1.0) #C(4.0 -2.0) #C(4.0 -3.0) #C(4.0 0.0)))
     (linear-algebra:subtract matrix1 matrix2))
    ;; Scalar1
    (assert-float-equal
     #2A((#C( 3.0  0.0) #C( 6.0  3.0) #C( 9.0  3.0) #C(12.0  3.0))
         (#C( 6.0 -3.0) #C( 6.0  0.0) #C( 9.0  6.0) #C(12.0  6.0))
         (#C( 9.0 -3.0) #C( 9.0 -6.0) #C( 9.0  0.0) #C(12.0  9.0))
         (#C(12.0 -3.0) #C(12.0 -6.0) #C(12.0 -9.0) #C(12.0  0.0)))
     (linear-algebra:subtract matrix1 matrix2 :scalar1 2.0))
    ;; Scalar2
    (assert-float-equal
     #2A((#C(0.0  0.0) #C(0.0  0.0) #C(0.0  0.0) #C(0.0 0.0))
         (#C(0.0 -0.0) #C(0.0  0.0) #C(0.0  0.0) #C(0.0 0.0))
         (#C(0.0 -0.0) #C(0.0 -0.0) #C(0.0  0.0) #C(0.0 0.0))
         (#C(0.0 -0.0) #C(0.0 -0.0) #C(0.0 -0.0) #C(0.0 0.0)))
     (linear-algebra:subtract matrix1 matrix2 :scalar2 2.0))
    ;; Scalar1 & Scalar2
    (assert-float-equal
     #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0) #C(4.0 1.0))
         (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0) #C(4.0 2.0))
         (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0) #C(4.0 3.0))
         (#C(4.0 -1.0) #C(4.0 -2.0) #C(4.0 -3.0) #C(4.0 0.0)))
     (linear-algebra:subtract
      matrix1 matrix2 :scalar1 2.0 :scalar2 3.0))))

(define-test nsubtract-hermitian-matrix
  (:tag :hermitian-matrix :nsubtract)
  ;; No scalar
  (let ((matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '((#C(2.0  0.0) #C(4.0  2.0) #C(6.0  2.0) #C(8.0 2.0))
             (#C(4.0 -2.0) #C(4.0  0.0) #C(6.0  4.0) #C(8.0 4.0))
             (#C(6.0 -2.0) #C(6.0 -4.0) #C(6.0  0.0) #C(8.0 6.0))
             (#C(8.0 -2.0) #C(8.0 -4.0) #C(8.0 -6.0) #C(8.0 0.0))))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0) #C(4.0 1.0))
              (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0) #C(4.0 2.0))
              (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0) #C(4.0 3.0))
              (#C(4.0 -1.0) #C(4.0 -2.0) #C(4.0 -3.0) #C(4.0 0.0))))))
    (assert-eq matrix1 (linear-algebra:nsubtract matrix1 matrix2))
    (assert-float-equal
     #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0) #C(4.0 1.0))
         (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0) #C(4.0 2.0))
         (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0) #C(4.0 3.0))
         (#C(4.0 -1.0) #C(4.0 -2.0) #C(4.0 -3.0) #C(4.0 0.0)))
     matrix1))
  ;; Scalar1
  (let ((matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0) #C(4.0 1.0))
             (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0) #C(4.0 2.0))
             (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0) #C(4.0 3.0))
             (#C(4.0 -1.0) #C(4.0 -2.0) #C(4.0 -3.0) #C(4.0 0.0))))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :initial-contents
          #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0) #C(4.0 1.0))
              (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0) #C(4.0 2.0))
              (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0) #C(4.0 3.0))
              (#C(4.0 -1.0) #C(4.0 -2.0) #C(4.0 -3.0) #C(4.0 0.0))))))
    (assert-eq
     matrix1
     (linear-algebra:nsubtract matrix1 matrix2 :scalar1 2.0))
    (assert-float-equal
     #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0) #C(4.0 1.0))
         (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0) #C(4.0 2.0))
         (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0) #C(4.0 3.0))
         (#C(4.0 -1.0) #C(4.0 -2.0) #C(4.0 -3.0) #C(4.0 0.0)))
     matrix1))
  ;; Scalar2
  (let ((matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '((#C( 3.0  0.0) #C( 6.0  3.0) #C( 9.0  3.0) #C(12.0  3.0))
             (#C( 6.0 -3.0) #C( 6.0  0.0) #C( 9.0  6.0) #C(12.0  6.0))
             (#C( 9.0 -3.0) #C( 9.0 -6.0) #C( 9.0  0.0) #C(12.0  9.0))
             (#C(12.0 -3.0) #C(12.0 -6.0) #C(12.0 -9.0) #C(12.0  0.0))))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0) #C(4.0 1.0))
              (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0) #C(4.0 2.0))
              (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0) #C(4.0 3.0))
              (#C(4.0 -1.0) #C(4.0 -2.0) #C(4.0 -3.0) #C(4.0 0.0))))))
    (assert-eq
     matrix1
     (linear-algebra:nsubtract matrix1 matrix2 :scalar2 2.0))
    (assert-float-equal
     #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0) #C(4.0 1.0))
         (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0) #C(4.0 2.0))
         (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0) #C(4.0 3.0))
         (#C(4.0 -1.0) #C(4.0 -2.0) #C(4.0 -3.0) #C(4.0 0.0)))
     matrix1))
  ;; Scalar1 & Scalar2
  (let ((*epsilon* (* 3F0 single-float-epsilon))
        (matrix1
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          (make-array
           '(4 4) :initial-contents
           '((#C(2.0  0.0) #C(4.0  2.0) #C(6.0  2.0) #C(8.0 2.0))
             (#C(4.0 -2.0) #C(4.0  0.0) #C(6.0  4.0) #C(8.0 4.0))
             (#C(6.0 -2.0) #C(6.0 -4.0) #C(6.0  0.0) #C(8.0 6.0))
             (#C(8.0 -2.0) #C(8.0 -4.0) #C(8.0 -6.0) #C(8.0 0.0))))))
        (matrix2
         (linear-algebra:make-matrix
          4 4 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0) #C(4.0 1.0))
              (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0) #C(4.0 2.0))
              (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0) #C(4.0 3.0))
              (#C(4.0 -1.0) #C(4.0 -2.0) #C(4.0 -3.0) #C(4.0 0.0))))))
    (assert-eq
     matrix1
     (linear-algebra:nsubtract
      matrix1 matrix2 :scalar1 2.0 :scalar2 3.0))
    (assert-float-equal
     #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0) #C(4.0 1.0))
         (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0) #C(4.0 2.0))
         (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0) #C(4.0 3.0))
         (#C(4.0 -1.0) #C(4.0 -2.0) #C(4.0 -3.0) #C(4.0 0.0)))
     matrix1)))

(define-test product-hermitian-matrix
  (:tag :hermitian-matrix :product)
  ;; Row vector - Hermitian matrix
  (assert-true
   (typep
    (linear-algebra:product
     (linear-algebra:row-vector 1.0 2.0 3.0)
     (unit-hermitian-matrix 3))
    'linear-algebra:row-vector))
  (assert-float-equal
   #(#C(14.0 -5.0) #C(15.0 -5.0) #C(18.0 5.0))
   (linear-algebra:product
    (linear-algebra:row-vector 1.0 2.0 3.0)
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:hermitian-matrix
     :initial-contents
     #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0))
         (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0))
         (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0))))))
  (assert-float-equal
   #(#C(29.399998 -10.5) #C(31.499999 -10.5) #C(37.8 10.5))
   (linear-algebra:product
    (linear-algebra:row-vector 1.0 2.0 3.0)
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:hermitian-matrix
     :initial-contents
     #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0))
         (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0))
         (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0))))
    2.1))
  (assert-error
   'error
   (linear-algebra:product
    (linear-algebra:row-vector 1.0 2.0 3.0 4.0 5.0 6.0)
    (unit-hermitian-matrix 3)))
  ;; Hermitian matrix - column vector
  (assert-true
   (typep
    (linear-algebra:product
     (unit-hermitian-matrix 3)
     (linear-algebra:column-vector 1.0 2.0 3.0))
    'linear-algebra:column-vector))
  (assert-float-equal
   #(#C(14.0 5.0) #C(15.0 5.0) #C(18.0 -5.0))
   (linear-algebra:product
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:hermitian-matrix
     :initial-contents
     #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0))
         (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0))
         (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0))))
    (linear-algebra:column-vector 1.0 2.0 3.0)))
  (assert-float-equal
   #(#C(29.399998 10.5) #C(31.499999 10.5) #C(37.8 -10.5))
   (linear-algebra:product
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:hermitian-matrix
     :initial-contents
     #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0))
         (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0))
         (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0))))
    (linear-algebra:column-vector 1.0 2.0 3.0)
    2.1))
  (assert-error
   'error
   (linear-algebra:product
    (unit-hermitian-matrix 3)
    (linear-algebra:column-vector 1.0 2.0 3.0 4.0 5.0 6.0)))
  ;; Hermitian matrix - matrix
  (assert-true
   (typep
    (linear-algebra:product
     (unit-hermitian-matrix 3) (unit-hermitian-matrix 3))
    'linear-algebra:hermitian-matrix))
  (assert-float-equal
   #2A((#C(16.0   0.0) #C(17.0  0.0) #C(16.0 11.0))
       (#C(17.0   0.0) #C(22.0  0.0) #C(22.0  9.0))
       (#C(16.0 -11.0) #C(22.0 -9.0) #C(32.0  0.0)))
   (linear-algebra:product
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:hermitian-matrix
     :initial-contents
     #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0))
         (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0))
         (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0))))
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:hermitian-matrix
     :initial-contents
     #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0))
         (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0))
         (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0))))))
  (assert-float-equal
   #2A((#C(33.6 0.0) #C(35.699997 0.0) #C(33.6 23.099999))
       (#C(35.699997 0.0) #C(46.199997 0.0) #C(46.199997 18.9))
       (#C(33.6 -23.099999) #C(46.199997 -18.9) #C(67.2 0.0)))
   (linear-algebra:product
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:hermitian-matrix
     :initial-contents
     #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0))
         (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0))
         (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0))))
    (linear-algebra:make-matrix
     3 3 :matrix-type 'linear-algebra:hermitian-matrix
     :initial-contents
     #2A((#C(1.0  0.0) #C(2.0  1.0) #C(3.0  1.0))
         (#C(2.0 -1.0) #C(2.0  0.0) #C(3.0  2.0))
         (#C(3.0 -1.0) #C(3.0 -2.0) #C(3.0  0.0))))
    2.1))
  (assert-error
   'error
   (linear-algebra:product
    (unit-hermitian-matrix 3) (unit-hermitian-matrix 4))))

(define-test solve-hermitian-matrix
  (:tag :hermitian-matrix :solve)
  (let ((vector2 (linear-algebra:column-vector 2.0 1.0))
        (vector3 (linear-algebra:column-vector 2.3 1.2 2.2))
        (matrix2
         (linear-algebra:make-matrix
          2 2 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          #2A((#C(2.0 0.0) #C(1.0 -2.0))
              (#C(1.0 2.0) #C(3.0  0.0)))))
        (matrix3
         (linear-algebra:make-matrix
          3 3 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          #2A((#C(3.31 0.0) #C(1.26 -2.0) #C(1.37 -3.0))
              (#C(1.26 2.0) #C(2.23  0.0) #C(2.31 -1.5))
              (#C(1.37 3.0) #C(2.31  1.5) #C(8.15  0.0))))))
    ;; 2x2
    (assert-float-equal
     #(#C(5.0 2.0) #C(0.0 -4.0))
     (linear-algebra:solve matrix2 vector2))
    (assert-float-equal
     #2A((#C(2.0 0.0) #C(1.0 -2.0))
         (#C(1.0 2.0) #C(3.0  0.0)))
     matrix2)
    (assert-float-equal
     #(2.0 1.0) vector2 (linear-algebra::contents vector2))
    ;; 3x3
    (assert-float-equal
     #(#C( 3.5175734   3.4673646)
       #C( 3.3198433  -4.3366637)
       #C(-0.78414906 -1.2595192))
     (linear-algebra:solve matrix3 vector3))
    (assert-float-equal
     #2A((#C(3.31 0.0) #C(1.26 -2.0) #C(1.37 -3.0))
         (#C(1.26 2.0) #C(2.23  0.0) #C(2.31 -1.5))
         (#C(1.37 3.0) #C(2.31  1.5) #C(8.15  0.0)))
     matrix3)
    (assert-float-equal
     #(2.3 1.2 2.2) vector3 (linear-algebra::contents vector3))))

(define-test nsolve-hermitian-matrix
  (:tag :hermitian-matrix :nsolve)
  (let ((vector2 (linear-algebra:column-vector 2.0 1.0))
        (vector3 (linear-algebra:column-vector 2.3 1.2 2.2))
        (matrix2
         (linear-algebra:make-matrix
          2 2 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          #2A((#C(2.0 0.0) #C(1.0 -2.0))
              (#C(1.0 2.0) #C(3.0  0.0)))))
        (matrix3
         (linear-algebra:make-matrix
          3 3 :matrix-type 'linear-algebra:hermitian-matrix
          :initial-contents
          #2A((#C(3.31 0.0) #C(1.26 -2.0) #C(1.37 -3.0))
              (#C(1.26 2.0) #C(2.23  0.0) #C(2.31 -1.5))
              (#C(1.37 3.0) #C(2.31  1.5) #C(8.15  0.0))))))
    ;; 2x2
    (assert-float-equal
     #(#C(5.0 2.0) #C(0.0 -4.0))
     (linear-algebra:nsolve matrix2 vector2))
    (assert-float-equal
     #2A((#C(2.0 0.0) #C(0.5 -1.0))
         (#C(0.5 1.0) #C(0.5  0.0)))
     matrix2)
    (assert-float-equal
     #(#C(5.0 2.0) #C(0.0 -4.0)) vector2)
    ;; 3x3
    (assert-float-equal
     #(#C( 3.5175734   3.4673646)
       #C( 3.3198433  -4.3366637)
       #C(-0.78414906 -1.2595192))
     (linear-algebra:nsolve matrix3 vector3))
    (assert-float-equal
     #2A((#C(3.31       0.0)       #C( 0.38066468 -0.6042296) #C( 0.4138973   -0.9063445))
         (#C(0.38066468 0.6042296) #C( 0.54190326  0.0)       #C(-0.044656467 -2.1882146))
         (#C(0.4138973  0.9063445) #C(-0.044656467 2.1882146) #C( 2.2680602   -3.7252904E-9)))
     matrix3)
    (assert-float-equal
     #(#C( 3.5175734   3.4673646)
       #C( 3.3198433  -4.3366637)
       #C(-0.78414906 -1.2595192))
     vector3)))
