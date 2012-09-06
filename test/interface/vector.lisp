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

;;; Givens Rotation
(define-test givens-rotation
  ;; g = 0
  (multiple-value-bind (c s r)
      (linear-algebra:givens-rotation 0 0)
    (assert-rational-equal 1 c)
    (assert-rational-equal 0 s)
    (assert-rational-equal 0 r))
  ;; Real f
  (multiple-value-bind (c s r)
      (linear-algebra:givens-rotation 1 0)
    (assert-rational-equal 1 c)
    (assert-rational-equal 0 s)
    (assert-rational-equal 1 r))
  ;; Imaginary f
  (multiple-value-bind (c s r)
      (linear-algebra:givens-rotation #C(1 1) 0)
    (assert-rational-equal 1 c)
    (assert-rational-equal 0 s)
    (assert-rational-equal #C(1 1) r))
  ;; f = 0 , negative g
  (multiple-value-bind (c s r)
      (linear-algebra:givens-rotation 0 -1)
    (assert-rational-equal  0 c)
    (assert-rational-equal -1 s)
    (assert-rational-equal  1 r))
  ;; f = 0 , Real g
  (multiple-value-bind (c s r)
      (linear-algebra:givens-rotation 0 1)
    (assert-rational-equal 0 c)
    (assert-rational-equal 1 s)
    (assert-rational-equal 1 r))
  ;; f = 0 , Imaginary g
  (multiple-value-bind (c s r)
      (linear-algebra:givens-rotation 0 #C(1 1))
    (assert-rational-equal 0 c)
    (assert-float-equal #C(0.70710677 -0.70710677) s)
    (assert-float-equal  1.4142135 r))
  ;; Rational f and g
  (multiple-value-bind (c s r)
      (linear-algebra:givens-rotation 1 2)
    (assert-float-equal 0.4472136 c)
    (assert-float-equal 0.8944272 s)
    (assert-float-equal 2.236068  r))
  ;; Float f and g
  (multiple-value-bind (c s r)
      (linear-algebra:givens-rotation 1.1 2.3)
    (assert-float-equal 0.4314555 c)
    (assert-float-equal 0.9021342 s)
    (assert-float-equal 2.5495098 r))
  ;; Complex rational f and g
  (multiple-value-bind (c s r)
      (linear-algebra:givens-rotation #C(1 2) #C(3 4))
    (assert-float-equal 0.40824828 c)
    (assert-float-equal #C(0.8981462 0.16329929) s)
    (assert-float-equal #C(2.4494898 4.8989797) r))
  ;; Complex float f and g
  (multiple-value-bind (c s r)
      (linear-algebra:givens-rotation #C(1.2 2.3) #C(3.4 4.5))
    (assert-float-equal 0.4178801 c)
    (assert-float-equal #C(0.8959895 0.15026298) s)
    (assert-float-equal #C(2.8716373 5.503971) r)))

;;; Jacobi Rotation
(define-test jacobi-rotation
  ;; Symmetric test
  (multiple-value-bind (a b c s)
      (linear-algebra:jacobi-rotation 1.1 3.3 5.5)
    (assert-float-equal -0.66610646 a)
    (assert-float-equal  7.266106   b)
    (assert-float-equal  0.8816746  c)
    (assert-float-equal -0.4718579  s))
  ;; Hermitian test
  (multiple-value-bind (a b c s)
      (linear-algebra:jacobi-rotation 1.1 #C(3.3 7.7) 5.5)
    (assert-float-equal -5.3614073  a)
    (assert-float-equal 11.961407   b)
    (assert-float-equal  0.79183334 c)
    (assert-float-equal #C(-0.24058115 0.561356) s)))

