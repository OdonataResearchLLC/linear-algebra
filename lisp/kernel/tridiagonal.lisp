#|

  Linear Algebra Tridiagonal Algorithm

  Copyright (c) 2011-2015, Odonata Research LLC

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

;;; Algorithm 4.32, pg. 91
;;; Linear equation solver for tridiagonal A

;;; Step 1

(defun tridiagonal-factorization (array)
  "Return the factorization of the tridiagonal array."
  (loop
   with end = (1- (array-dimension array 0))
   initially
   (setf (aref array 0 2) (/ (aref array 0 2) (aref array 0 1)))
   for row from 1 below end
   as alpha =
   (- (aref array row 1) (* (aref array row 0) (aref array (1- row) 2)))
   do
   (setf
    (aref array row 1) alpha
    (aref array row 2) (/ (aref array row 2) alpha))
   finally
   (decf
    (aref array end 1)
    (* (aref array end 0) (aref array (1- end) 2)))
   ;; Return the factored array.
   (return array)))

;;; Step 2

(defun tridiagonal-update (array vector)
  "Update the solution vector using the factored array."
  (loop
   with end = (1- (array-dimension array 0))
   initially
   (setf (aref vector 0) (/ (aref vector 0) (aref array 0 1)))
   for row from 1 upto end do
   (setf
    (aref vector row)
    (/ (- (aref vector row)
          (* (aref array row 0) (aref vector (1- row))))
       (aref array row 1)))
   finally (return vector)))

;;; Step 3

(defun tridiagonal-backsubstitution (array vector)
  "Perform backsubstitution to obtain the solution."
  (loop
   with end = (1- (array-dimension array 0))
   for row from (1- end) downto 0 do
   (decf
    (aref vector row)
    (* (aref array row 2) (aref vector (1+ row))))
   finally (return vector)))

(defun tridiagonal-solver (array vector)
  "Linear equation solver for a tridiagonal matrix."
  (progn
    (tridiagonal-factorization array)
    (tridiagonal-update array vector)
    (tridiagonal-backsubstitution array vector)))
