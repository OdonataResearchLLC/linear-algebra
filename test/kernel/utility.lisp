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

(define-test copy-array
  (:tag :kernel :utility)
  (assert-float-equal
   #(1.1 2.2 3.3 4.4 5.5)
   (linear-algebra-kernel:copy-array
    #(1.1 2.2 3.3 4.4 5.5)))
  (assert-float-equal
   #2A((1.1 1.2 1.3) (2.1 2.2 2.3) (3.1 3.2 3.3))
   (linear-algebra-kernel:copy-array
    #2A((1.1 1.2 1.3) (2.1 2.2 2.3) (3.1 3.2 3.3)))))

(defclass class-0 () ())
(defclass class-1 (class-0) ())
(defclass class-a (class-1) ())
(defclass class-sub-a (class-a) ())
(defclass class-b (class-1) ())
(defclass class-sub-b (class-b) ())

(define-test common-class-of
  (:tag :kernel :utility)
  (let ((object-a (make-instance 'class-a))
        (object-sub-a (make-instance 'class-sub-a))
        (object-b (make-instance 'class-b))
        (object-sub-b (make-instance 'class-sub-b)))
    (assert-eq
     (find-class 'class-1)
     (linear-algebra-kernel:common-class-of object-a object-b))
    (assert-eq
     (find-class 'class-1)
     (linear-algebra-kernel:common-class-of object-sub-a object-b))
    (assert-eq
     (find-class 'class-1)
     (linear-algebra-kernel:common-class-of object-a object-sub-b))))

(define-test common-array-element-type
  (:tag :kernel :utility)
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

(define-test specific-array-element-type
  (:tag :kernel :utility)
  (flet ((genarray (element-type)
           (make-array
            3 :element-type element-type
            :initial-element (coerce 1 element-type))))
    ;; Real float
    (assert-eq
     'single-float
     (linear-algebra-kernel:specific-array-element-type
      (genarray 'single-float)))
    (assert-eq
     'double-float
     (linear-algebra-kernel:specific-array-element-type
      (genarray 'double-float)))
    ;; Complex float
    (assert-equal
     (type-of (complex 1.0 0.0))
     (linear-algebra-kernel:specific-array-element-type
      (genarray '(complex single-float))))
    (assert-equal
     (type-of (complex 1D0 0D0))
     (linear-algebra-kernel:specific-array-element-type
      (genarray '(complex double-float))))))

(define-test complex-equal
  (:tag :kernel :utility :complex :equal)
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
  (:tag :kernel :utility :equal)
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
