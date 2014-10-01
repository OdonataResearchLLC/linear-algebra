#|

  Linear Algebra Conjugate Gradient Algorithm

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

;;; Algorithm 4.31, pg 84
;;; Conjugate gradient method

(defun %default-cg-epsilon (array vector)
  "Return a default epsilon for the conjugate gradient method."
  (case (common-array-element-type array vector)
    ((short-float single-float) (* 64 single-float-epsilon))
    ((double-float long-float) (* 64 double-float-epsilon))
    (otherwise (* 128 single-float-epsilon))))

(defun %initialize-cg-solution (array)
  "Return an initial solution vector for the conjugate gradient."
  (zero-vector (array-dimension array 0) (array-element-type array)))

(defun %initialize-cg-residual (array vector solution)
  "Return the initial residual vector for the conjugate gradient."
  (subtract-vector
   vector (product-array-vector array solution) nil nil))

(defun %negative-residual (residual)
  "Return the negative of the residual."
  (loop
   with size = (length residual)
   with result = (zero-vector size (array-element-type residual))
   for index below size do
   (setf (aref result index) (- (aref residual index)))
   finally (return result)))

(defun conjugate-gradient-solver
       (array vector &optional epsilon (limit 25))
  "Linear system solver using the conjugate gradient method."
  (loop
   with epsilon = (or epsilon (%default-cg-epsilon array vector))
   with solution = (%initialize-cg-solution array)
   with residual = (%initialize-cg-residual array vector solution)
   with -residual = (%negative-residual residual)
   for iteration below limit
   ;; Vector Adk
   as adk = (product-array-vector array residual)
   then (product-array-vector array residual nil adk)
   ;; Scalar denominator
   as denominator = (inner-product-vector residual adk nil)
   ;; Alpha
   as alpha =
   (/ (inner-product-vector residual -residual -1) denominator)
   ;; Update
   do
   (setf
    solution (nadd-vector solution residual nil alpha)
    -residual (nadd-vector -residual adk nil alpha))
   ;; Test convergence
   until (< (norm-vector -residual :infinity) epsilon) do
   (setf
    residual
    (nsubtract-vector
     residual -residual
     (/ (inner-product-vector -residual adk nil) denominator)
     nil))
   ;; Return the solution
   finally (return (values solution iteration))))
