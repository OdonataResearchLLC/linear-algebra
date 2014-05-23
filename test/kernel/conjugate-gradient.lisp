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

(define-test conjugate-gradient-solver
  (:tag :kernel :conjugate-gradient)
  ;; 2x2 from NumAlgoC
  (assert-float-equal
   #(2.0 0.33333334)
   (linear-algebra-kernel:conjugate-gradient-solver
    (make-array '(2 2) :initial-contents '((2.0 0.0) (0.0 3.0)))
    (make-array 2 :initial-contents '(4.0 1.0))))
  ;; 2x2
  (assert-float-equal
   #(3.265307 -1.3265313)
   (linear-algebra-kernel:conjugate-gradient-solver
    (make-array '(2 2) :initial-contents '((1.1 1.2) (1.2 2.2)))
    (make-array 2 :initial-contents '(2.0 1.0))))
  ;; 3x3
  (assert-float-equal
   #(3.5856622 -2.3062859 0.7900801)
   (linear-algebra-kernel:conjugate-gradient-solver
    (make-array
     '(3 3)
     :initial-contents
     '((1.15 1.26 1.37)
       (1.26 2.23 2.31)
       (1.37 2.31 3.31)))
    (make-array 3 :initial-contents '(2.3 1.2 2.2)))))
