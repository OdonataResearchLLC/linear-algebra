#|

 Linear Algebra in Common Lisp Unit Tests

 Copyright (c) 2011, Thomas M. Hermann
 All rights reserved.

 Redistribution and  use  in  source  and  binary  forms, with or without
 modification, are permitted  provided  that the following conditions are
 met:

   o  Redistributions of  source  code  must  retain  the above copyright
      notice, this list of conditions and the following disclaimer.
   o  Redistributions in binary  form  must reproduce the above copyright
      notice, this list of  conditions  and  the  following disclaimer in
      the  documentation  and/or   other   materials  provided  with  the
      distribution.
   o  The names of the contributors may not be used to endorse or promote
      products derived from this software without  specific prior written
      permission.

 THIS SOFTWARE IS  PROVIDED  BY  THE  COPYRIGHT  HOLDERS AND CONTRIBUTORS
 "AS IS"  AND  ANY  EXPRESS  OR  IMPLIED  WARRANTIES, INCLUDING,  BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES  OF MERCHANTABILITY AND FITNESS FOR A
 PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT  NOT LIMITED TO,
 PROCUREMENT OF  SUBSTITUTE  GOODS  OR  SERVICES;  LOSS  OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION)  HOWEVER  CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER  IN  CONTRACT,  STRICT  LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR  OTHERWISE)  ARISING  IN  ANY  WAY  OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

|#

(in-package :linear-algebra-test)

(define-test lapy2
  "sqrt |x|^2 + |y|^2"
  (:tag :utility)
  ;; Real values
  (dolist (args (cartesian-product '(-3.0 3.0) '(-4.0 4.0)))
    (assert-float-equal
     5.0 (apply #'linear-algebra-kernel:lapy2 args)))
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
       6.024948 (apply #'linear-algebra-kernel:lapy2 args)))))

(define-test lapy3
  "sqrt |x|^2 + |y|^2 + |z|^2"
  (:tag :utility)
  ;; Real values
  (dolist (args (nary-product '(-2.0 2.0) '(-3.0 3.0) '(-4.0 4.0)))
    (assert-float-equal
     5.3851647 (apply #'linear-algebra-kernel:lapy3 args)))
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
       10.49333 (apply #'linear-algebra-kernel:lapy3 args)))))

;;; Givens Rotation

(define-test givens-rotation
  (:tag :utility :rotation)
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
  (:tag :utility :rotation)
  ;; Symmetric test
  (multiple-value-bind (a b c s)
      (linear-algebra-kernel:jacobi-rotation 1.1 3.3 5.5)
    (assert-float-equal -0.66610646 a)
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
  (:tag :utility :rotation)
  (multiple-value-bind (beta tau vector)
      (linear-algebra-kernel:householder-reflection
       #C(1.0 2.0) (vector 1.0 2.0 3.0 4.0 5.0))
    (assert-float-equal -7.745967 beta)
    (assert-float-equal #C(1.1290995 0.2581989) tau)
    (assert-float-equal
     #(#C(0.2 -0.4) #C(0.4 -0.8) #C(0.6 -1.2) #C(0.8 -1.6) #C(1.0 -2.0))
     vector)))

(define-test common-class-of
  (:tag :utility)
  (let ((object-a (make-array '(3 3) :adjustable t))
        (object-v (make-array 3 :adjustable t))
        (object-l (make-list 3))
        (class-a (find-class 'array))
        (class-v (find-class 'vector))
        (class-t (find-class t)))
    (assert-eq
     class-a
     (linear-algebra-kernel:common-class-of object-a object-a))
    (assert-eq
     class-a
     (linear-algebra-kernel:common-class-of object-a object-v))
    (assert-eq
     class-a
     (linear-algebra-kernel:common-class-of object-v object-a))
    (assert-eq
     class-v
     (linear-algebra-kernel:common-class-of object-v object-v))
    (assert-eq
     class-t
     (linear-algebra-kernel:common-class-of object-v object-l t))
    (assert-error
     'error
     (linear-algebra-kernel:common-class-of object-v object-l))))

(define-test common-array-element-type
  (:tag :utility)
  (let ((array-s (make-array 0 :element-type 'single-float))
        (array-d (make-array 0 :element-type 'double-float)))
    (assert-eq
     'single-float
     (linear-algebra-kernel:common-array-element-type array-s array-s))
    (assert-eq
     'double-float
     (linear-algebra-kernel:common-array-element-type array-s array-d))
    (assert-eq
     'double-float
     (linear-algebra-kernel:common-array-element-type array-d array-s))
    (assert-eq
     'double-float
     (linear-algebra-kernel:common-array-element-type array-d array-d))))

(define-test complex-equal
  (:tag :utility :complex :equal)
  ;; complex float
  (assert-true
   (linear-algebra-kernel:complex-equal #C(1.0 2.0) #C(1.0 2.0)))
  (assert-true
   (linear-algebra-kernel:complex-equal 1.0 #C(1.0 0.0)))
  (assert-true
   (linear-algebra-kernel:complex-equal #C(1.0 0.0) 1.0))
  (assert-false
   (linear-algebra-kernel:complex-equal #C(1.0 2.0) #C(2.0 1.0)))
  (assert-false
   (linear-algebra-kernel:complex-equal 1.0 #C(0.0 1.0)))
  (assert-false
   (linear-algebra-kernel:complex-equal #C(0.0 1.0) 1.0))
  ;; complex integer
  (assert-true
   (linear-algebra-kernel:complex-equal #C(1 2) #C(1 2)))
  ;; Error
  (assert-error 'error (linear-algebra-kernel:complex-equal 1.0 1.0))
  (assert-error 'error (linear-algebra-kernel:complex-equal 1 1)))

(define-test number-equal
  (:tag :utility :equal)
  ;; float
  (assert-true (linear-algebra-kernel:number-equal 2.2 2.2))
  (assert-true (linear-algebra-kernel:number-equal 2 2.0))
  (assert-true (linear-algebra-kernel:number-equal 2.0 2))
  (assert-false (linear-algebra-kernel:number-equal 2 2.2))
  ;; rational
  (assert-true (linear-algebra-kernel:number-equal 1/3 1/3))
  (assert-true (linear-algebra-kernel:number-equal 3 3))
  (assert-false (linear-algebra-kernel:number-equal 1/3 3))
  ;; complex float
  (assert-true
   (linear-algebra-kernel:number-equal #C(1.1 2.2) #C(1.1 2.2)))
  (assert-true
   (linear-algebra-kernel:number-equal #C(1.0 2.0) #C(1 2)))
  (assert-true
   (linear-algebra-kernel:number-equal #C(1 2) #C(1.0 2.0)))
  (assert-false
   (linear-algebra-kernel:number-equal #C(1.1 2.2) #C(2.2 1.1)))
  ;; complex rational
  (assert-true
   (linear-algebra-kernel:number-equal #C(1 2) #C(1 2)))
  (assert-true
   (linear-algebra-kernel:number-equal #C(1/2 1/2) #C(1/2 1/2)))
  (assert-false
   (linear-algebra-kernel:number-equal #C(1 2) #C(1/2 1/2)))
  ;; error
  (assert-error 'error (linear-algebra-kernel:number-equal 1 t))
  (assert-error 'error (linear-algebra-kernel:number-equal t 1))
  (assert-error 'error (linear-algebra-kernel:number-equal t t)))
