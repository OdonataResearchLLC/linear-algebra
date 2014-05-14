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
   for index below size do
   ;; Step 1.1
   (multiple-value-bind (scale sumsq)
       (sumsq-row array index :end index)
     (setf
      (aref array index index)
      (sqrt (- (aref array index index) (* scale scale sumsq)))))
   ;; Step 1.2
   (loop
    for row from (1+ index) below size do
    (setf
     (aref array row index)
     (/
      (- (aref array row index)
         (loop for column below index sum
               (* (aref array row column)
                  (aref array index column))))
      (aref array index index))))
   ;; Return the factored array
   finally return array))
