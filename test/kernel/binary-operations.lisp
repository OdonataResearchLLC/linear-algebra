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

;;; Vector inner product

(define-test binary-operation-inner-product-vector
  ;; Real vectors
  (assert-rational-equal
   55
   (linear-algebra-kernel:binary-operation
    :inner-product #(1 2 3 4 5) #(1 2 3 4 5) nil nil))
  (assert-float-equal
   55F0
   (linear-algebra-kernel:binary-operation
    :inner-product
    #(1.0 2.0 3.0 4.0 5.0)
    #(1.0 2.0 3.0 4.0 5.0)
    nil nil))
  (assert-float-equal
   55D0
   (linear-algebra-kernel:binary-operation
    :inner-product
    #(1D0 2D0 3D0 4D0 5D0)
    #(1D0 2D0 3D0 4D0 5D0)
    nil nil))
  ;; Real vectors with conjugate keyword
  (assert-rational-equal
   55
   (linear-algebra-kernel:binary-operation
    :inner-product
    #(1 2 3 4 5) #(1 2 3 4 5) nil t))
  ;; Complex vectors
  (assert-rational-equal
   #C(8 18)
   (linear-algebra-kernel:binary-operation
    :inner-product
    #(#C(1 1) #C(2 1) #C(3 1))
    #(#C(1 2) #C(2 2) #C(3 2))
    nil nil))
  (assert-float-equal
   #C(8.0 18.0)
   (linear-algebra-kernel:binary-operation
    :inner-product
    #(#C(1.0 1.0) #C(2.0 1.0) #C(3.0 1.0))
    #(#C(1.0 2.0) #C(2.0 2.0) #C(3.0 2.0))
    nil nil))
  (assert-float-equal
   #C(8D0 18D0)
   (linear-algebra-kernel:binary-operation
    :inner-product
    #(#C(1D0 1D0) #C(2D0 1D0) #C(3D0 1D0))
    #(#C(1D0 2D0) #C(2D0 2D0) #C(3D0 2D0))
    nil nil))
  ;; Complex conjugate
  (assert-rational-equal
   #C(20 6)
   (linear-algebra-kernel:binary-operation
    :inner-product
    #(#C(1 1) #C(2 1) #C(3 1))
    #(#C(1 2) #C(2 2) #C(3 2))
    nil t))
  (assert-float-equal
   #C(20.0 6.0)
   (linear-algebra-kernel:binary-operation
    :inner-product
    #(#C(1.0 1.0) #C(2.0 1.0) #C(3.0 1.0))
    #(#C(1.0 2.0) #C(2.0 2.0) #C(3.0 2.0))
    nil t))
  (assert-float-equal
   #C(20D0 6D0)
   (linear-algebra-kernel:binary-operation
    :inner-product
    #(#C(1D0 1D0) #C(2D0 1D0) #C(3D0 1D0))
    #(#C(1D0 2D0) #C(2D0 2D0) #C(3D0 2D0))
    nil t)))

;;; Array addition

(define-test binary-operation-add-array
  (let ((array
         #2A((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)
             (5.1 5.2 5.3 5.4))))
    ;; No scalar
    (assert-float-equal
     #2A(( 2.2  2.4  2.6  2.8)
         ( 4.2  4.4  4.6  4.8)
         ( 6.2  6.4  6.6  6.8)
         ( 8.2  8.4  8.6  8.8)
         (10.2 10.4 10.6 10.8))
     (linear-algebra-kernel:binary-operation
      :add array array nil nil))
    ;; Scalar1
    (assert-float-equal
     #2A(( 3.3  3.6  3.9  4.2)
         ( 6.3  6.6  6.9  7.2)
         ( 9.3  9.6  9.9 10.2)
         (12.3 12.6 12.9 13.2)
         (15.3 15.6 15.9 16.2))
     (linear-algebra-kernel:binary-operation
      :add array array 2.0 nil))
    ;; Scalar2
    (assert-float-equal
     #2A(( 3.3  3.6  3.9  4.2)
         ( 6.3  6.6  6.9  7.2)
         ( 9.3  9.6  9.9 10.2)
         (12.3 12.6 12.9 13.2)
         (15.3 15.6 15.9 16.2))
     (linear-algebra-kernel:binary-operation
      :add array array nil 2.0))
    ;; Scalar1 & Scalar2
    (assert-float-equal
     #2A(( 5.5  6.0  6.5  7.0)
         (10.5 11.0 11.5 12.0)
         (15.5 16.0 16.5 17.0)
         (20.5 21.0 21.5 22.0)
         (25.5 26.0 26.5 27.0))
     (linear-algebra-kernel:binary-operation
      :add array array 2.0 3.0))))

(define-test binary-operation-nadd-array
  ;; No scalar
  (let ((array1
         (make-array
          '(5 4) :initial-contents
          '((1.1 1.2 1.3 1.4)
            (2.1 2.2 2.3 2.4)
            (3.1 3.2 3.3 3.4)
            (4.1 4.2 4.3 4.4)
            (5.1 5.2 5.3 5.4))))
        (array2
         #2A((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)
             (5.1 5.2 5.3 5.4))))
    (assert-eq
     array1
     (linear-algebra-kernel:binary-operation
      :nadd array1 array2 nil nil))
    (assert-float-equal
     #2A(( 2.2  2.4  2.6  2.8)
         ( 4.2  4.4  4.6  4.8)
         ( 6.2  6.4  6.6  6.8)
         ( 8.2  8.4  8.6  8.8)
         (10.2 10.4 10.6 10.8))
     array1))
  ;; Scalar1
  (let ((array1
         (make-array
          '(5 4) :initial-contents
          '((1.1 1.2 1.3 1.4)
            (2.1 2.2 2.3 2.4)
            (3.1 3.2 3.3 3.4)
            (4.1 4.2 4.3 4.4)
            (5.1 5.2 5.3 5.4))))
        (array2
         #2A((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)
             (5.1 5.2 5.3 5.4))))
    (assert-eq
     array1
     (linear-algebra-kernel:binary-operation
      :nadd array1 array2 2.0 nil))
    (assert-float-equal
     #2A(( 3.3  3.6  3.9  4.2)
         ( 6.3  6.6  6.9  7.2)
         ( 9.3  9.6  9.9 10.2)
         (12.3 12.6 12.9 13.2)
         (15.3 15.6 15.9 16.2))
     array1))
  ;; Scalar2
  (let ((array1
         (make-array
          '(5 4) :initial-contents
          '((1.1 1.2 1.3 1.4)
            (2.1 2.2 2.3 2.4)
            (3.1 3.2 3.3 3.4)
            (4.1 4.2 4.3 4.4)
            (5.1 5.2 5.3 5.4))))
        (array2
         #2A((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)
             (5.1 5.2 5.3 5.4))))
    (assert-eq
     array1
     (linear-algebra-kernel:binary-operation
      :nadd array1 array2 nil 2.0))
    (assert-float-equal
     #2A(( 3.3  3.6  3.9  4.2)
         ( 6.3  6.6  6.9  7.2)
         ( 9.3  9.6  9.9 10.2)
         (12.3 12.6 12.9 13.2)
         (15.3 15.6 15.9 16.2))
     array1))
  ;; Scalar1 & Scalar2
  (let ((array1
         (make-array
          '(5 4) :initial-contents
          '((1.1 1.2 1.3 1.4)
            (2.1 2.2 2.3 2.4)
            (3.1 3.2 3.3 3.4)
            (4.1 4.2 4.3 4.4)
            (5.1 5.2 5.3 5.4))))
        (array2
         #2A((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)
             (5.1 5.2 5.3 5.4))))
    (assert-eq
     array1
     (linear-algebra-kernel:binary-operation
      :nadd array1 array2 2.0 3.0))
    (assert-float-equal
     #2A(( 5.5  6.0  6.5  7.0)
         (10.5 11.0 11.5 12.0)
         (15.5 16.0 16.5 17.0)
         (20.5 21.0 21.5 22.0)
         (25.5 26.0 26.5 27.0))
     array1)))

(define-test binary-operation-subtract-array
  (let ((*epsilon* (* 3F0 single-float-epsilon))
        (array1
         #2A(( 2.2  2.4  2.6  2.8)
             ( 4.2  4.4  4.6  4.8)
             ( 6.2  6.4  6.6  6.8)
             ( 8.2  8.4  8.6  8.8)
             (10.2 10.4 10.6 10.8)))
        (array2
         #2A((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)
             (5.1 5.2 5.3 5.4))))
    ;; No scalar
    (assert-float-equal
     #2A((1.1 1.2 1.3 1.4)
         (2.1 2.2 2.3 2.4)
         (3.1 3.2 3.3 3.4)
         (4.1 4.2 4.3 4.4)
         (5.1 5.2 5.3 5.4))
     (linear-algebra-kernel:binary-operation
      :subtract array1 array2 nil nil))
    ;; Scalar1
    (assert-float-equal
     #2A(( 3.3  3.6  3.9  4.2)
         ( 6.3  6.6  6.9  7.2)
         ( 9.3  9.6  9.9 10.2)
         (12.3 12.6 12.9 13.2)
         (15.3 15.6 15.9 16.2))
     (linear-algebra-kernel:binary-operation
      :subtract array1 array2 2.0 nil))
    ;; Scalar2
    (assert-float-equal
     #2A((0.0 0.0 0.0 0.0)
         (0.0 0.0 0.0 0.0)
         (0.0 0.0 0.0 0.0)
         (0.0 0.0 0.0 0.0)
         (0.0 0.0 0.0 0.0))
     (linear-algebra-kernel:binary-operation
      :subtract array1 array2 nil 2.0))
    ;; Scalar1 & Scalar2
    (assert-float-equal
     #2A((1.1 1.2 1.3 1.4)
         (2.1 2.2 2.3 2.4)
         (3.1 3.2 3.3 3.4)
         (4.1 4.2 4.3 4.4)
         (5.1 5.2 5.3 5.4))
     (linear-algebra-kernel:binary-operation
      :subtract array1 array2 2.0 3.0))))

(define-test binary-operation-nsubtract-array
  ;; No scalar
  (let ((array1
         (make-array
          '(5 4) :initial-contents
          '(( 2.2  2.4  2.6  2.8)
            ( 4.2  4.4  4.6  4.8)
            ( 6.2  6.4  6.6  6.8)
            ( 8.2  8.4  8.6  8.8)
            (10.2 10.4 10.6 10.8))))
        (array2
         #2A((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)
             (5.1 5.2 5.3 5.4))))
    (assert-eq
     array1
     (linear-algebra-kernel:binary-operation
      :nsubtract array1 array2 nil nil))
    (assert-float-equal
     #2A((1.1 1.2 1.3 1.4)
         (2.1 2.2 2.3 2.4)
         (3.1 3.2 3.3 3.4)
         (4.1 4.2 4.3 4.4)
         (5.1 5.2 5.3 5.4))
     array1))
  ;; Scalar1
  (let ((array1
         (make-array
          '(5 4) :initial-contents
          '((1.1 1.2 1.3 1.4)
            (2.1 2.2 2.3 2.4)
            (3.1 3.2 3.3 3.4)
            (4.1 4.2 4.3 4.4)
            (5.1 5.2 5.3 5.4))))
        (array2
         #2A((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)
             (5.1 5.2 5.3 5.4))))
    (assert-eq
     array1
     (linear-algebra-kernel:binary-operation
      :nsubtract array1 array2 2.0 nil))
    (assert-float-equal
     #2A((1.1 1.2 1.3 1.4)
         (2.1 2.2 2.3 2.4)
         (3.1 3.2 3.3 3.4)
         (4.1 4.2 4.3 4.4)
         (5.1 5.2 5.3 5.4))
     array1))
  ;; Scalar2
  (let ((*epsilon* (* 4F0 single-float-epsilon))
        (array1
         (make-array
          '(5 4) :initial-contents
          '(( 3.3  3.6  3.9  4.2)
            ( 6.3  6.6  6.9  7.2)
            ( 9.3  9.6  9.9 10.2)
            (12.3 12.6 12.9 13.2)
            (15.3 15.6 15.9 16.2))))
        (array2
         #2A((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)
             (5.1 5.2 5.3 5.4))))
    (assert-eq
     array1
     (linear-algebra-kernel:binary-operation
      :nsubtract array1 array2 nil 2.0))
    (assert-float-equal
     #2A((1.1 1.2 1.3 1.4)
         (2.1 2.2 2.3 2.4)
         (3.1 3.2 3.3 3.4)
         (4.1 4.2 4.3 4.4)
         (5.1 5.2 5.3 5.4))
     array1))
  ;; Scalar1 & Scalar2
  (let ((*epsilon* (* 3F0 single-float-epsilon))
        (array1
         (make-array
          '(5 4) :initial-contents
          '(( 2.2  2.4  2.6  2.8)
            ( 4.2  4.4  4.6  4.8)
            ( 6.2  6.4  6.6  6.8)
            ( 8.2  8.4  8.6  8.8)
            (10.2 10.4 10.6 10.8))))
        (array2
         #2A((1.1 1.2 1.3 1.4)
             (2.1 2.2 2.3 2.4)
             (3.1 3.2 3.3 3.4)
             (4.1 4.2 4.3 4.4)
             (5.1 5.2 5.3 5.4))))
    (assert-eq
     array1
     (linear-algebra-kernel:binary-operation
      :nsubtract array1 array2 2.0 3.0))
    (assert-float-equal
     #2A((1.1 1.2 1.3 1.4)
         (2.1 2.2 2.3 2.4)
         (3.1 3.2 3.3 3.4)
         (4.1 4.2 4.3 4.4)
         (5.1 5.2 5.3 5.4))
     array1)))
