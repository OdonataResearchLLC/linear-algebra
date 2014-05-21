#|

  Linear Algebra Cholesky Algorithm

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

;;; Algorithm 4.27, Step 1, pg. 80
;;; The standard Cholesky decomposition

(defun cholesky-decomposition (array)
  "Factor A = LL^T."
  (loop
   with size = (array-dimension array 0)
   for index-j below size do
   ;; Step 1.1
   (multiple-value-bind (scale sumsq)
       (sumsq-row array index-j :end index-j)
     (setf
      (aref array index-j index-j)
      (sqrt (- (aref array index-j index-j) (* scale scale sumsq)))))
   ;; Step 1.2
   (loop
    for index-k from (1+ index-j) below size do
    (setf
     (aref array index-k index-j)
     (/
      (- (aref array index-k index-j)
         (loop
          for index-i below index-j sum
          (* (aref array index-k index-i)
             (aref array index-j index-i))))
      (aref array index-j index-j))))
   ;; Return the factored array
   finally return array))

;;; Algorithm 4.29, pg. 82
;;; Simplified linear system solver via root-free Cholesky

(defun root-free-cholesky-decomposition (array)
  "Factor A = LDL^t."
  (loop
   with size = (array-dimension array 0)
   for index-j below size do
   (loop
    for index-i below index-j
    as var-h = (aref array index-j index-i)
    do
    (setf
     (aref array index-j index-i)
     (/ var-h (aref array index-i index-i)))
    (loop
     for index-k from (1+ index-i) upto index-j do
     (decf
      (aref array index-j index-k)
      (* var-h (aref array index-k index-i)))))
   finally return array))

(defun cholesky-solver (array vector)
  "Linear system solver for positive definite matrices using the
root-free Cholesky decomposition."
  (let ((size (array-dimension array 0)))
    ;; Step 1, decomposition
    (setq array (root-free-cholesky-decomposition array))
    ;; Step 2.1 & 2.2
    (loop
     for index-j below size do
     (loop
      for index-i below index-j do
      (decf
       (aref vector index-j)
       (* (aref array index-j index-i) (aref vector index-i)))))
    ;; Step 2.3 & 3.2
    (loop
     for index-j from (1- size) downto 0 do
     (setf
      (aref vector index-j)
      (/ (aref vector index-j) (aref array index-j index-j)))
     (loop
      for index-i from (1+ index-j) below size do
      (decf
       (aref vector index-j)
       (* (aref array index-i index-j) (aref vector index-i)))))
    ;; Return the solution
    vector))
