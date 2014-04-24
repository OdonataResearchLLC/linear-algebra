#|

  Linear Algebra in Common Lisp

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

(in-package :linear-algebra-kernel)

;;; Rotations

(defun givens-rotation (f g)
  "Return c,s,r defined from the Givens rotation."
  (cond
    ((zerop g)
     (values 1 0 f))
    ((zerop f)
     (values 0 (signum (conjugate g)) (abs g)))
    (t
     (let ((abs-f (abs f))
           (sqrtfg (sumsq2 f g)))
       (values
        (/ abs-f sqrtfg)
        (/ (* (signum f) (conjugate g)) sqrtfg)
        (* (signum f) sqrtfg))))))

(defun jacobi-rotation (x y z)
  "Return a, b, cos(theta) and sin(theta) terms from the Jacobi rotation."
  (let* ((yabs (abs y))
         (tau  (/ (- x z) 2.0 yabs))
         (tee  (/ (float-sign tau)
                  (+ (abs tau) (sqrt (+ 1.0 (expt tau 2))))))
         (cos-theta (/ (sqrt (+ 1.0 (expt tee 2))))) ; Invert sqrt
         (sin-theta (* cos-theta tee)))
    (values
     ;; a : first eigenvalue
     (+ (* cos-theta cos-theta x)
        (* 2.0 cos-theta sin-theta yabs)
        (* sin-theta sin-theta z))
     ;; b : second eigenvalue
     (+ (* sin-theta sin-theta x)
        (* -2.0 cos-theta sin-theta yabs)
        (* cos-theta cos-theta z))
     ;; Cosine theta
     cos-theta
     ;; Sine theta
     (* (conjugate (signum y)) sin-theta))))

(defun householder-reflection (alpha vector)
  "Return Beta, Tau and the Householder vector."
  (let* ((beta
          (- (float-sign
              (realpart alpha)
              (sumsq2 alpha (norm-vector vector 2)))))
         (tau  (- 1 (/ alpha beta))))
    (values
     beta tau
     (dotimes (index (length vector) vector)
       (setf
        (aref vector index)
        (/ (aref vector index) alpha))))))
