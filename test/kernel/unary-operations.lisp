#|

 Linear Algebra in Common Lisp Unit Tests

 Copyright (c) 2011-2012, Odonata Research LLC
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

;;; Scaled binary operations

(define-test unary-sumsq-vector
  (:tag :unary :sumsq)
  ;; Real
  (multiple-value-bind (scale sumsq)
      (linear-algebra-kernel:sumsq-vector
       #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) 0 1)
    (assert-rational-equal 6 scale)
    (assert-rational-equal 73/18 sumsq))
  ;; Complex
  (multiple-value-bind (scale sumsq)
      (linear-algebra-kernel:sumsq-vector
       #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
         #C(-2 3) #C(-3 1) #C(-1 0)) 0 1)
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
       0 1)
    (assert-float-equal 4.5 scale)
    (assert-float-equal 8.997532 sumsq)))

(define-test unary-sump-vector
  (:tag :unary :sump)
  ;; Real
  (multiple-value-bind (scale sump)
      (linear-algebra-kernel:sump-vector
       #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) 2 0 1)
    (assert-rational-equal 6 scale)
    (assert-rational-equal 73/18 sump))
  (multiple-value-bind (scale sump)
      (linear-algebra-kernel:sump-vector
       #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) 3 0 1)
    (assert-rational-equal 6 scale)
    (assert-rational-equal 1 sump))
  ;; Complex
  (multiple-value-bind (scale sump)
      (linear-algebra-kernel:sump-vector
       #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
         #C(-2 3) #C(-3 1) #C(-1 0))
       2 0 1)
    (assert-float-equal 4.0 scale)
    (assert-float-equal #C(2.75 -1.125) sump))
  (multiple-value-bind (scale sump)
      (linear-algebra-kernel:sump-vector
       #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
         #C(-2 3) #C(-3 1) #C(-1 0))
       3 0 1)
    (assert-float-equal 4.0 scale)
    (assert-float-equal #C(2.6639833 0.54687494) sump)))

(define-test unary-sump-array
  (:tag :unary :sump)
  (multiple-value-bind (scale sump)
      (linear-algebra-kernel:sump-array
       #2A((1.1 1.2 1.3 1.4 1.5)
           (2.1 2.2 2.3 2.4 2.5)
           (3.1 3.2 3.3 3.4 3.5)
           (4.1 4.2 4.3 4.4 4.5))
       3.5 0 1)
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
