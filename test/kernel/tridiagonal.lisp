#|

  Linear Algebra in Common Lisp Unit Tests

  Copyright (c) 2011-2015, Odonata Research LLC

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

(define-test tridiagonal-factorization
  (:tag :kernel :tridiagonal)
  (assert-float-equal
   #2A((0.0  1.1        1.2727273)
       (1.4  0.4181819  5.499999)
       (2.3 -9.3499975 -0.3422461)
       (3.2  5.4951878  0.0))
   (linear-algebra-kernel::tridiagonal-factorization
    (make-array
     '(4 3) :initial-contents
     '((0.0 1.1 1.4) (1.4 2.2 2.3) (2.3 3.3 3.2) (3.2 4.4 0.0))))))

(define-test tridiagonal-update
  (:tag :kernel :tridiagonal)
  (assert-float-equal
   #(1.0909091 6.630434 1.3850269 -0.24240965)
   (linear-algebra-kernel::tridiagonal-update
    #2A((0.0  1.1        1.2727273)
        (1.4  0.4181819  5.499999)
        (2.3 -9.3499975 -0.3422461)
        (3.2  5.4951878  0.0))
    (make-array 4 :initial-contents '(1.2 4.3 2.3 3.1)))))

(define-test tridiagonal-backsubstitution
  (:tag :kernel :tridiagonal)
  (assert-float-equal
   #(1.7666159 -0.5309124 1.3020632 -0.24240965)
   (linear-algebra-kernel::tridiagonal-backsubstitution
    #2A((0.0  1.1        1.2727273)
        (1.4  0.4181819  5.499999)
        (2.3 -9.3499975 -0.3422461)
        (3.2  5.4951878  0.0))
    (make-array
     4 :initial-contents
     '(1.0909091 6.630434 1.3850269 -0.24240965)))))

(define-test tridiagonal-solver
  (:tag :kernel :tridiagonal)
  (assert-float-equal
   #(1.7666159 -0.5309124 1.3020632 -0.24240965)
   (linear-algebra-kernel:tridiagonal-solver
    (make-array
     '(4 3) :initial-contents
     '((0.0 1.1 1.4) (1.4 2.2 2.3) (2.3 3.3 3.2) (3.2 4.4 0.0)))
    (make-array 4 :initial-contents '(1.2 4.3 2.3 3.1)))))
