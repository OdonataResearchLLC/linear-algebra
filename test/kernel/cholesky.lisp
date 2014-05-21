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

(define-test cholesky-decomposition
  (:tag :kernel :cholesky)
  ;; 2x2
  (assert-float-equal
   #2A((1.0488088 1.2)
       (1.1441552 0.9438798))
   (linear-algebra-kernel:cholesky-decomposition
    (make-array '(2 2) :initial-contents '((1.1 1.2) (1.2 2.2)))))
  ;; 3x3
  (assert-float-equal
   #2A((1.0723806 1.26       1.37)
       (1.1749561 0.92167145 2.31)
       (1.2775316 0.8777058  0.95265186))
   (linear-algebra-kernel:cholesky-decomposition
    (make-array
     '(3 3)
     :initial-contents
     '((1.15 1.26 1.37)
       (1.26 2.23 2.31)
       (1.37 2.31 3.31)))))
  ;; 3x3 from Wikipedia
  (assert-float-equal
   #2A(( 2.0 12.0 -16.0)
       ( 6.0  1.0 -43.0)
       (-8.0  5.0   3.0))
   (linear-algebra-kernel:cholesky-decomposition
    (make-array
     '(3 3)
     :initial-contents
     '((4 12 -16) (12 37 -43) (-16 -43 98))))))

(define-test root-free-cholesky-decomposition
  (:tag :kernel :cholesky)
  ;; 2x2
  (assert-float-equal
   #2A((1.1 1.2) (1.0909091 0.8909091))
   (linear-algebra-kernel:root-free-cholesky-decomposition
    (make-array '(2 2) :initial-contents '((1.1 1.2) (1.2 2.2)))))
  ;; 3x3
  (assert-float-equal
   #2A((1.15      1.26       1.37)
       (1.0956522 0.84947825 2.31)
       (1.1913043 0.9522979  0.90754557))
   (linear-algebra-kernel:root-free-cholesky-decomposition
    (make-array
     '(3 3)
     :initial-contents
     '((1.15 1.26 1.37)
       (1.26 2.23 2.31)
       (1.37 2.31 3.31))))))

(define-test cholesky-solver
  (:tag :kernel :cholesky)
  ;; 2x2
  (assert-float-equal
   #(3.2653065 -1.3265308)
   (linear-algebra-kernel:cholesky-solver
    (make-array '(2 2) :initial-contents '((1.1 1.2) (1.2 2.2)))
    (make-array 2 :initial-contents '(2.0 1.0))))
  ;; 3x3
  (assert-float-equal
   #(3.5856622 -2.306286 0.79007966)
   (linear-algebra-kernel:cholesky-solver
    (make-array
     '(3 3)
     :initial-contents
     '((1.15 1.26 1.37)
       (1.26 2.23 2.31)
       (1.37 2.31 3.31)))
    (make-array 3 :initial-contents '(2.3 1.2 2.2)))))
