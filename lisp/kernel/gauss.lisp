#|

  Linear Algebra Gauss Algorithm

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

 References
 [NumAlgoC] Gisela Engeln-Mullges and Frank Uhlig "Numerical
            Algorithms with C", Springer, 1996
            ISBN: 3-540-60530-4

|#

(in-package :linear-algebra-kernel)

;;; Algorithm 4.22, pg. 74
;;; Triangular matrix factorization with column pivot search

;;; Step 1
(defun initialize-pivot-selection-vector (size)
  "Return a new, initialized, pivot vector."
  (let ((pivot-selection-vector (make-array size)))
    (dotimes (index size pivot-selection-vector)
      (setf (svref pivot-selection-vector index) index))))

(defun column-pivot-search (array column)
  "Return the row index of the maximum value in the column."
  (loop
   with max-row = column
   and max-element = (abs (aref array column column))
   for row from (1+ column) below (array-dimension array 0)
   as element = (abs (aref array row column))
   when (< max-element element) do
   (setq
    max-row row
    max-element element)
   finally return max-row))

(defun swap-rows (array i0 jth)
  "Interchange the "
  (loop
   for column below (array-dimension array 1) do
   (rotatef (aref array i0 column) (aref array jth column))
   finally return array))

(defun column-pivot (array pivot-selection-vector column)
  "Return the LR pivot of the array."
  ;; Step 2.1
  (let ((i0 (column-pivot-search array column)))
    (unless (= i0 column)
      (rotatef
       (svref pivot-selection-vector i0)
       (svref pivot-selection-vector column))
      (swap-rows array i0 column)))
  ;; Check for a singular matrix
  (when (float-equal 0.0 (aref array column column))
    (error "~A is singular." array))
  ;; Step 2.2
  (loop
   with size = (array-dimension array 0)
   and pivot-column-element = (aref array column column)
   for row from (1+ column) below size do
   ;; Step 2.2.1
   (setf
    (aref array row column)
    (/ (aref array row column) pivot-column-element))
   ;; Step 2.2.2
   (loop
    for col from (1+ column) below size do
    (setf
     (aref array row col)
     (- (aref array row col)
        (* (aref array column col) (aref array row column))))))
  ;; Return the array & pivot selection vector
  (values array pivot-selection-vector))

;;; TODO: Replace or add an implicitly scaled pivot search
;;; FIXME: Improve the accuracy of the factorization

(defun gauss-factorization (array)
  "Return the Gauss factorization of the array."
  (loop
   with size = (array-dimension array 0)
   with pivot-selection-vector =
   (initialize-pivot-selection-vector size)
   for column below (1- size) do
   (column-pivot array pivot-selection-vector column)
   finally return (values array pivot-selection-vector)))

;;; Algorithm 4.23, pg. 75
;;; Gauss algorithm with column pivot search

(defun gauss-solver (array vector)
  "Gauss algorithm with column pivot search."
  (let ((size (array-dimension array 0))
        (solution
         (make-array
          (length vector)
          :element-type (array-element-type vector)
          :initial-element 0.0)))
  ;; Step 1
  (multiple-value-bind (factored pivot-selection-vector)
      (gauss-factorization array)
    ;; Step 2
    (loop
     initially
     (setf
      (aref solution 0)
      (aref vector (svref pivot-selection-vector 0)))
     for row from 1 below size do
     (setf
      (aref solution row)
      (- (aref vector (svref pivot-selection-vector row))
         (loop for col from 0 below row sum
               (* (aref factored row col) (aref solution col))))))
    ;; Step 3
    (loop
     with end = (1- size)
     initially
     (setf
      (aref solution end)
      (/ (aref solution end) (aref factored end end)))
     for row downfrom (1- end) to 0 do
     (setf
      (aref solution row)
      (/ (- (aref solution row)
            (loop for col from (1+ row) below size sum
                  (* (aref factored row col) (aref solution col))))
         (aref factored row row)))))
  ;; Return the solution
  solution))

;;; Algorithm 4.25, pg. 77
;;; Finding A^-1 via Gauss algorithm with partial column pivot search

(defun unit-pivot-value
       (pivot-selection-vector row column
        &optional (array-type *read-default-float-format*))
  "Return 1.0 if column equals the value at row of the pivot selection
vector, otherwise 0.0."
  (coerce
   (if (= column (svref pivot-selection-vector row)) 1.0 0.0)
   array-type))

(defun gauss-invert (array)
  "Find A^-1 via Gauss algorithm with partial column pivot search."
  (let* ((size (array-dimension array 0))
         (array-type (array-element-type array))
         (solution
          (make-array
           (list size size)
           :element-type array-type
           :initial-element (coerce 0.0 array-type))))
  ;; Step 1
  (multiple-value-bind (factored pivot-selection-vector)
      (gauss-factorization array)
    ;; Step 2
    (loop
     for column below size do
     (setf
      (aref solution 0 column)
      (unit-pivot-value pivot-selection-vector 0 column array-type))
     (loop
      for row from 1 below size do
      (setf
       (aref solution row column)
       (- (unit-pivot-value
           pivot-selection-vector row column array-type)
          (loop for index from 0 below row sum
                (* (aref factored row index)
                   (aref solution index column)))))))
    ;; Step 3
    (loop
     with end = (1- size)
     for column below size do
     (setf
      (aref solution end column)
      (/ (aref solution end column) (aref factored end end)))
     (loop
      for row downfrom (1- end) to 0 do
      (setf
       (aref solution row column)
       (/ (- (aref solution row column)
             (loop for index from (1+ row) below size sum
                   (* (aref factored row index)
                      (aref solution index column))))
          (aref factored row row))))))
  ;; Return the solution
  solution))
