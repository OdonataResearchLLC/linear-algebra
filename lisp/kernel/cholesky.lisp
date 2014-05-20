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

;;; Algorithm 4.27, pg. 80
;;; Linear system solver for positive definite matrix using the
;;; standard Cholesky decomposition

(defun standard-cholesky-decomposition (array)
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

#| NOTE:

   Step 2 in algorithm 4.27 incorrectly describes the forward
   substitution algorithm. Division by the diagonal term cannot be
   performed in a separate loop from the substitution steps.

|#

(defun standard-cholesky-linear-solver (array vector)
  "Linear system solver for positive definite matrices using the
standard Cholesky decomposition."
  (let ((size (array-dimension array 0))
        (ll-array))
    ;; Step 1, decomposition
    (setq ll-array (standard-cholesky-decomposition array))
    ;; Step 2, forward substitution
    (loop
     for index-i below size do
     (setf
      (aref vector index-i)
      (/
       (- (aref vector index-i)
          (loop
           for index-j below index-i sum
           (* (aref ll-array index-i index-j) (aref vector index-j))))
       (aref array index-i index-i))))
    ;; Step 3, backward substitution
    (loop
     for index-i from (1- size) downto 0 do
     (setf
      (aref vector index-i)
      (/
       (- (aref vector index-i)
          (loop
           for index-j from (1+ index-i) below size sum
           (* (aref ll-array index-j index-i)
              (aref vector index-j))))
       (aref ll-array index-i index-i))))
    ;; Return the solution
    vector))
