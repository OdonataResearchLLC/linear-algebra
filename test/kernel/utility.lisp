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

(define-test common-class-of
  (:tag :kernel :utility)
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
