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

(define-test sumsq2
  "sqrt |x|^2 + |y|^2"
  (:tag :unary :sumsq2)
  ;; Real values
  (dolist (args (cartesian-product '(-3.0 3.0) '(-4.0 4.0)))
    (assert-float-equal
     5.0 (apply #'linear-algebra-kernel:sumsq2 args)))
  ;; Complex values
  (let ((args1
         (mapcar
          (lambda (x) (apply #'complex x))
          (cartesian-product '(-1.1 1.1) '(-2.2 2.2))))
        (args2
         (mapcar
          (lambda (x) (apply #'complex x))
          (cartesian-product '(-3.3 3.3) '(-4.4 4.4)))))
    (dolist (args (cartesian-product args1 args2))
      (assert-float-equal
       6.024948 (apply #'linear-algebra-kernel:sumsq2 args)))))

(define-test sumsq3
  "sqrt |x|^2 + |y|^2 + |z|^2"
  (:tag :unary :sumsq2)
  ;; Real values
  (dolist (args (nary-product '(-2.0 2.0) '(-3.0 3.0) '(-4.0 4.0)))
    (assert-float-equal
     5.3851647 (apply #'linear-algebra-kernel:sumsq3 args)))
  ;; Complex values
  (let ((args1
         (mapcar
          (lambda (x) (apply #'complex x))
          (cartesian-product '(-1.1 1.1) '(-2.2 2.2))))
        (args2
         (mapcar
          (lambda (x) (apply #'complex x))
          (cartesian-product '(-3.3 3.3) '(-4.4 4.4))))
        (args3
         (mapcar
          (lambda (x) (apply #'complex x))
          (cartesian-product '(-5.5 5.5) '(-6.6 6.6)))))
    (dolist (args (nary-product args1 args2 args3))
      (assert-float-equal
       10.49333 (apply #'linear-algebra-kernel:sumsq3 args)))))

(define-test unary-sumsq-vector
  (:tag :unary :sumsq)
  ;; Real
  (multiple-value-bind (scale sumsq)
      (linear-algebra-kernel:sumsq-vector
       #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) 1 0)
    (assert-rational-equal 6 scale)
    (assert-rational-equal 73/18 sumsq))
  ;; Complex
  (multiple-value-bind (scale sumsq)
      (linear-algebra-kernel:sumsq-vector
       #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
         #C(-2 3) #C(-3 1) #C(-1 0)) 1 0)
    (assert-float-equal 4.0 scale)
    (assert-float-equal #C (2.75 -1.125) sumsq)))

(define-test unary-sumsq-array
  (:tag :unary :sumsq)
  (multiple-value-bind (scale sumsq)
      (linear-algebra-kernel:sumsq-array
       #2A((1.1 1.2 1.3 1.4 1.5)
           (2.1 2.2 2.3 2.4 2.5)
           (3.1 3.2 3.3 3.4 3.5)
           (4.1 4.2 4.3 4.4 4.5))
       1 0)
    (assert-float-equal 4.5 scale)
    (assert-float-equal 8.997532 sumsq)))

(define-test unary-sump-vector
  (:tag :unary :sump)
  ;; Real
  (multiple-value-bind (scale sump)
      (linear-algebra-kernel:sump-vector
       #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) 2 1 0)
    (assert-rational-equal 6 scale)
    (assert-rational-equal 73/18 sump))
  (multiple-value-bind (scale sump)
      (linear-algebra-kernel:sump-vector
       #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) 3 1 0)
    (assert-rational-equal 6 scale)
    (assert-rational-equal 1 sump))
  ;; Complex
  (multiple-value-bind (scale sump)
      (linear-algebra-kernel:sump-vector
       #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
         #C(-2 3) #C(-3 1) #C(-1 0))
       2 1 0)
    (assert-float-equal 4.0 scale)
    (assert-float-equal #C(2.75 -1.125) sump))
  (multiple-value-bind (scale sump)
      (linear-algebra-kernel:sump-vector
       #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
         #C(-2 3) #C(-3 1) #C(-1 0))
       3 1 0)
    (assert-float-equal 4.0 scale)
    (assert-float-equal #C(2.6639833 0.54687494) sump)))

(define-test unary-sumsq-array
  (:tag :unary :sumsq)
  (multiple-value-bind (scale sumsq)
      (linear-algebra-kernel:sumsq-array
       #2A((1.1 1.2 1.3 1.4 1.5)
           (2.1 2.2 2.3 2.4 2.5)
           (3.1 3.2 3.3 3.4 3.5)
           (4.1 4.2 4.3 4.4 4.5))
       1 0)
    (assert-float-equal 4.5 scale)
    (assert-float-equal 8.997532 sumsq)))

(define-test unary-sump-array
  (:tag :unary :sump)
  (multiple-value-bind (scale sump)
      (linear-algebra-kernel:sump-array
       #2A((1.1 1.2 1.3 1.4 1.5)
           (2.1 2.2 2.3 2.4 2.5)
           (3.1 3.2 3.3 3.4 3.5)
           (4.1 4.2 4.3 4.4 4.5))
       3.5 1 0)
    (assert-float-equal 4.5 scale)
    (assert-float-equal 6.540154 sump)))

;;; Norm & supporting functions

(define-test %abs-vector
  (:tag :unary :norm)
  (assert-rational-equal
   #(6 5 4 3 2 1 0 1 2 3 4 5)
   (linear-algebra-kernel::%abs-vector
    #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5))))

;;; Taxicab norm

(define-test unary-norm-1-vector
  (:tag :unary :norm)
  (assert-rational-equal
   36 (linear-algebra-kernel:norm-vector
       #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) 1))
  (assert-float-equal
   19.535658
   (linear-algebra-kernel:norm-vector
    #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
      #C(-2 3) #C(-3 1) #C(-1 0))
    1)))

;;; Euclidean norm

(define-test unary-norm-2-vector
  (:tag :unary :norm)
  (assert-float-equal
   12.083046
   (linear-algebra-kernel:norm-vector
    #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)
    2))
  (assert-float-equal
   8.0
   (linear-algebra-kernel:norm-vector
    #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
      #C(-2 3) #C(-3 1) #C(-1 0))
    2)))

;;; P-norm

(define-test unary-norm-p-vector
  (:tag :unary :norm)
  (let ((data #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5))
        (zdata #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
                 #C(-2 3) #C(-3 1) #C(-1 0))))
    ;; norm
    (assert-float-equal
     8.732892 (linear-algebra-kernel:norm-vector data 3))
    (assert-float-equal
     6.064035 (linear-algebra-kernel:norm-vector zdata 3))))

;;; Infinity norm

(define-test unary-norm-infinity-vector
  (:tag :unary :norm)
  (assert-rational-equal
   6 (linear-algebra-kernel:norm-vector
      #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)
      :infinity))
  (assert-float-equal
   4.0 (linear-algebra-kernel:norm-vector
        #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
          #C(-2 3) #C(-3 1) #C(-1 0))
        :infinity)))

(define-test unary-norm-array
  (:tag :unary :norm)
  (let ((array
         #2A((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)
             (5.1 5.2 5.3 5.4))))
    (assert-float-equal
     17.0 (linear-algebra-kernel:norm-array array 1))
    (assert-float-equal
     5.4 (linear-algebra-kernel:norm-array array :max))
    (assert-float-equal
     15.858751 (linear-algebra-kernel:norm-array array :frobenius))
    (assert-float-equal
     21.0 (linear-algebra-kernel:norm-array array :infinity))))
