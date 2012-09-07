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

(define-test scaled-binary-op
  ;; No scalars
  (assert-float-equal
   2.2
   (funcall
    (linear-algebra-kernel:scaled-binary-op
     #'+ nil nil)
    1.1 1.1))
  (assert-float-equal
   1.1
   (funcall
    (linear-algebra-kernel:scaled-binary-op
     #'- nil nil)
    2.2 1.1))
  ;; Scalar 1
  (assert-float-equal
   3.3
   (funcall
    (linear-algebra-kernel:scaled-binary-op
     #'+ 2.0 nil)
    1.1 1.1))
  (assert-float-equal
   1.1
   (funcall
    (linear-algebra-kernel:scaled-binary-op
     #'- 2.0 nil)
    1.1 1.1))
  ;; Scalar 2
  (assert-float-equal
   3.3
   (funcall
    (linear-algebra-kernel:scaled-binary-op
     #'+ nil 2.0)
    1.1 1.1))
  (assert-float-equal
   1.1
   (funcall
    (linear-algebra-kernel:scaled-binary-op
     #'- nil 2.0)
    3.3 1.1))
  ;; Scalar 1 & 2
  (assert-float-equal
   5.5
   (funcall
    (linear-algebra-kernel:scaled-binary-op
     #'+ 2.0 3.0)
    1.1 1.1))
  (assert-float-equal
   1.1
   (funcall
    (linear-algebra-kernel:scaled-binary-op
     #'- 2.0 3.0)
    2.2 1.1)))

;;; Vector binary operations unit tests

(define-test binary-operation-add-vector
  ;; Real
  (let ((vector1 #(1.1 2.2 3.3 4.4))
        (vector2 #(1.1 2.2 3.3 4.4)))
    (assert-float-equal
     #(2.2 4.4 6.6 8.8)
     (linear-algebra-kernel:binary-operation
      :add vector1 vector2 nil nil))
    (assert-float-equal
     #(3.3 6.6 9.9 13.2)
     (linear-algebra-kernel:binary-operation
      :add vector1 vector2 2.0 nil))
    (assert-float-equal
     #(3.3 6.6 9.9 13.2)
     (linear-algebra-kernel:binary-operation
      :add vector1 vector2 nil 2.0))
    (assert-float-equal
     #(4.4 8.8 13.2 17.6)
     (linear-algebra-kernel:binary-operation
      :add vector1 vector2 2.0 2.0)))
  ;; Complex
  (let ((vector1 #(#C(1.1 2.2) #C(3.3 4.4)))
        (vector2 #(#C(1.1 2.2) #C(3.3 4.4))))
    (assert-float-equal
     #(#C(2.2 4.4) #C(6.6 8.8))
     (linear-algebra-kernel:binary-operation
      :add vector1 vector2 nil nil))
    (assert-float-equal
     #(#C(3.3 6.6) #C(9.9 13.2))
     (linear-algebra-kernel:binary-operation
      :add vector1 vector2 2.0 nil))
    (assert-float-equal
     #(#C(3.3 6.6) #C(9.9 13.2))
     (linear-algebra-kernel:binary-operation
      :add vector1 vector2 nil 2.0))
    (assert-float-equal
     #(#C(4.4 8.8) #C(13.2 17.6))
     (linear-algebra-kernel:binary-operation
      :add vector1 vector2 2.0 2.0))))

;;; Destructive vector addition

(define-test binary-operation-nadd-vector
  ;; Real
  (let ((vector1 (vector 1.1 2.2 3.3 4.4))
        (vector2 (vector 1.1 2.2 3.3 4.4)))
    (assert-eq
     vector1
     (linear-algebra-kernel:binary-operation
      :nadd vector1 vector2 nil nil))
    (assert-float-equal #(2.2 4.4 6.6 8.8) vector1)
    (assert-float-equal
     #(4.4 8.8 13.2 17.6)
     (linear-algebra-kernel:binary-operation
      :nadd vector1 vector2 nil 2.0))
    (assert-float-equal
     #(9.9 19.8 29.7 39.6)
     (linear-algebra-kernel:binary-operation
      :nadd vector1 vector2 2.0 nil))
    (assert-float-equal
     #(22.0 44.0 66.0 88.0)
     (linear-algebra-kernel:binary-operation
      :nadd vector1 vector2 2.0 2.0)))
  ;; Complex
  (let ((vector1 (vector #C(1.1 2.2) #C(3.3 4.4)))
        (vector2 (vector #C(1.1 2.2) #C(3.3 4.4))))
    (assert-eq
     vector1
     (linear-algebra-kernel:binary-operation
      :nadd vector1 vector2 nil nil))
    (assert-float-equal #(#C(2.2 4.4) #C(6.6 8.8)) vector1)
    (assert-float-equal
     #(#C(4.4 8.8) #C(13.2 17.6))
     (linear-algebra-kernel:binary-operation
      :nadd vector1 vector2 nil 2.0))
    (assert-float-equal
     #(#C(9.9 19.8) #C(29.7 39.6))
     (linear-algebra-kernel:binary-operation
      :nadd vector1 vector2 2.0 nil))
    (assert-float-equal
     #(#C(22.0 44.0) #C(66.0 88.0))
     (linear-algebra-kernel:binary-operation
      :nadd vector1 vector2 2.0 2.0))))

;;; Vector subtraction

(define-test binary-operation-subtract-vector
  ;; Real
  (let ((vector1 #(1.1 2.2 3.3 4.4))
        (vector2 #(1.1 2.2 3.3 4.4)))
    (assert-float-equal
     #(0.0 0.0 0.0 0.0)
     (linear-algebra-kernel:binary-operation
      :subtract vector1 vector2 nil nil))
    (assert-float-equal
     #(1.1 2.2 3.3 4.4)
     (linear-algebra-kernel:binary-operation
      :subtract vector1 vector2 2.0 nil))
    (assert-float-equal
     #(-1.1 -2.2 -3.3 -4.4)
     (linear-algebra-kernel:binary-operation
      :subtract vector1 vector2 nil 2.0))
    (assert-float-equal
     #(0.0 0.0 0.0 0.0)
     (linear-algebra-kernel:binary-operation
      :subtract vector1 vector2 2.0 2.0)))
  ;; Complex
  (let ((vector1 #(#C(1.1 2.2) #C(3.3 4.4)))
        (vector2 #(#C(1.1 2.2) #C(3.3 4.4))))
    (assert-float-equal
     #(#C(0.0 0.0) #C(0.0 0.0))
     (linear-algebra-kernel:binary-operation
      :subtract vector1 vector2 nil nil))
    (assert-float-equal
     #(#C(1.1 2.2) #C(3.3 4.4))
     (linear-algebra-kernel:binary-operation
      :subtract vector1 vector2 2.0 nil))
    (assert-float-equal
     #(#C(-1.1 -2.2) #C(-3.3 -4.4))
     (linear-algebra-kernel:binary-operation
      :subtract vector1 vector2 nil 2.0))
    (assert-float-equal
     #(#C(0.0 0.0) #C(0.0 0.0))
     (linear-algebra-kernel:binary-operation
      :subtract vector1 vector2 2.0 2.0))))

;;; Destructive vector subtraction

(define-test binary-operation-nsubtract-vector
  ;; Real
  (let ((vector1 (vector 1.1 2.2 3.3 4.4))
        (vector2 (vector 1.1 2.2 3.3 4.4)))
    (assert-eq
     vector1
     (linear-algebra-kernel:binary-operation
      :nsubtract vector1 vector2 nil nil))
    (assert-float-equal #(0.0 0.0 0.0 0.0) vector1)
    (assert-float-equal
     #(-2.2 -4.4 -6.6 -8.8)
     (linear-algebra-kernel:binary-operation
      :nsubtract vector1 vector2 nil 2.0))
    (assert-float-equal
     #(-5.5 -11.0 -16.5 -22.0)
     (linear-algebra-kernel:binary-operation
      :nsubtract vector1 vector2 2.0 nil))
    (assert-float-equal
     #(-13.2 -26.4 -39.6 -52.8)
     (linear-algebra-kernel:binary-operation
      :nsubtract vector1 vector2 2.0 2.0)))
  ;; Complex
  (let ((vector1 (vector #C(1.1 2.2) #C(3.3 4.4)))
        (vector2 (vector #C(1.1 2.2) #C(3.3 4.4))))
    (assert-eq
     vector1
     (linear-algebra-kernel:binary-operation
      :nsubtract vector1 vector2 nil nil))
    (assert-float-equal
     #(#C(0.0 0.0) #C(0.0 0.0)) vector1)
    (assert-float-equal
     #(#C(-2.2 -4.4) #C(-6.6 -8.8))
     (linear-algebra-kernel:binary-operation
      :nsubtract vector1 vector2 nil 2.0))
    (assert-float-equal
     #(#C(-5.5 -11.0) #C(-16.5 -22.0))
     (linear-algebra-kernel:binary-operation
      :nsubtract vector1 vector2 2.0 nil))
    (assert-float-equal
     #(#C(-13.2 -26.4) #C(-39.6 -52.8))
     (linear-algebra-kernel:binary-operation
      :nsubtract vector1 vector2 2.0 2.0))))
