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

;;; Givens Rotation

(define-test givens-rotation
  (:tag :kernel :rotation)
  ;; g = 0
  (multiple-value-bind (c s r)
      (linear-algebra-kernel:givens-rotation 0 0)
    (assert-rational-equal 1 c)
    (assert-rational-equal 0 s)
    (assert-rational-equal 0 r))
  ;; Real f
  (multiple-value-bind (c s r)
      (linear-algebra-kernel:givens-rotation 1 0)
    (assert-rational-equal 1 c)
    (assert-rational-equal 0 s)
    (assert-rational-equal 1 r))
  ;; Imaginary f
  (multiple-value-bind (c s r)
      (linear-algebra-kernel:givens-rotation #C(1 1) 0)
    (assert-rational-equal 1 c)
    (assert-rational-equal 0 s)
    (assert-rational-equal #C(1 1) r))
  ;; f = 0 , negative g
  (multiple-value-bind (c s r)
      (linear-algebra-kernel:givens-rotation 0 -1)
    (assert-rational-equal  0 c)
    (assert-rational-equal -1 s)
    (assert-rational-equal  1 r))
  ;; f = 0 , Real g
  (multiple-value-bind (c s r)
      (linear-algebra-kernel:givens-rotation 0 1)
    (assert-rational-equal 0 c)
    (assert-rational-equal 1 s)
    (assert-rational-equal 1 r))
  ;; f = 0 , Imaginary g
  (multiple-value-bind (c s r)
      (linear-algebra-kernel:givens-rotation 0 #C(1 1))
    (assert-rational-equal 0 c)
    (assert-float-equal #C(0.70710677 -0.70710677) s)
    (assert-float-equal  1.4142135 r))
  ;; Rational f and g
  (multiple-value-bind (c s r)
      (linear-algebra-kernel:givens-rotation 1 2)
    (assert-float-equal 0.4472136 c)
    (assert-float-equal 0.8944272 s)
    (assert-float-equal 2.236068  r))
  ;; Float f and g
  (multiple-value-bind (c s r)
      (linear-algebra-kernel:givens-rotation 1.1 2.3)
    (assert-float-equal 0.4314555 c)
    (assert-float-equal 0.9021342 s)
    (assert-float-equal 2.5495098 r))
  ;; Complex rational f and g
  (multiple-value-bind (c s r)
      (linear-algebra-kernel:givens-rotation #C(1 2) #C(3 4))
    (assert-float-equal 0.40824828 c)
    (assert-float-equal #C(0.8981462 0.16329929) s)
    (assert-float-equal #C(2.4494898 4.8989797) r))
  ;; Complex float f and g
  (multiple-value-bind (c s r)
      (linear-algebra-kernel:givens-rotation #C(1.2 2.3) #C(3.4 4.5))
    (assert-float-equal 0.4178801 c)
    (assert-float-equal #C(0.8959895 0.15026298) s)
    (assert-float-equal #C(2.8716373 5.503971) r)))

;;; Jacobi Rotation

(define-test jacobi-rotation
  (:tag :kernel :rotation)
  ;; Symmetric test
  (multiple-value-bind (a b c s)
      (linear-algebra-kernel:jacobi-rotation 1.1 3.3 5.5)
    (let ((*epsilon* (* 7 single-float-epsilon)))
      (assert-float-equal -0.66610646 a))
    (assert-float-equal  7.266106   b)
    (assert-float-equal  0.8816746  c)
    (assert-float-equal -0.4718579  s))
  ;; Hermitian test
  (multiple-value-bind (a b c s)
      (linear-algebra-kernel:jacobi-rotation 1.1 #C(3.3 7.7) 5.5)
    (assert-float-equal -5.3614073  a)
    (assert-float-equal 11.961407   b)
    (assert-float-equal  0.79183334 c)
    (assert-float-equal #C(-0.24058115 0.561356) s)))

;;; Householder Reflection

(define-test householder-reflection
  (:tag :kernel :rotation)
  (multiple-value-bind (beta tau vector)
      (linear-algebra-kernel:householder-reflection
       #C(1.0 2.0) (vector 1.0 2.0 3.0 4.0 5.0))
    (assert-float-equal -7.745967 beta)
    (assert-float-equal #C(1.1290995 0.2581989) tau)
    (assert-float-equal
     #(#C(0.2 -0.4) #C(0.4 -0.8) #C(0.6 -1.2) #C(0.8 -1.6) #C(1.0 -2.0))
     vector)))
