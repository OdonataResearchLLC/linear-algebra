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

(define-test symmetric-cholesky-decomposition
  (:tag :kernel :cholesky)
  ;; 2x2
  (assert-float-equal
   #2A((1.0488088 1.1441552)
       (1.1441552 0.9438798))
   (linear-algebra-kernel:symmetric-cholesky-decomposition
    (make-array '(2 2) :initial-contents '((1.1 1.2) (1.2 2.2)))))
  ;; 3x3
  (assert-float-equal
   #2A((1.0723806 1.1749561  1.2775316)
       (1.1749561 0.92167145 0.8777058)
       (1.2775316 0.8777058  0.95265186))
   (linear-algebra-kernel:symmetric-cholesky-decomposition
    (make-array
     '(3 3)
     :initial-contents
     '((1.15 1.26 1.37)
       (1.26 2.23 2.31)
       (1.37 2.31 3.31)))))
  ;; 3x3 from Wikipedia
  (assert-float-equal
   #2A(( 2.0  6.0 -8.0)
       ( 6.0  1.0  5.0)
       (-8.0  5.0  3.0))
   (linear-algebra-kernel:symmetric-cholesky-decomposition
    (make-array
     '(3 3)
     :initial-contents
     '((4 12 -16) (12 37 -43) (-16 -43 98)))))
  ;; 3x3 from Rosetta Code
  (assert-float-equal
   #2A(( 5.0 3.0 -1.0)
       ( 3.0 3.0  1.0)
       (-1.0 1.0  3.0))
   (linear-algebra-kernel:symmetric-cholesky-decomposition
    (make-array
     '(3 3)
     :initial-contents
     '((25.0 15.0 -5.0) (15.0 18.0 0.0) (-5.0 0.0 11.0)))))
  ;; 4x4 from Rosetta Code
  (assert-float-equal
   #2A(( 4.2426405 5.18545  12.727922  9.899495)
       ( 5.18545   6.565905  3.0460375 1.6245536)
       (12.727922  3.0460375 1.6497375 1.849715)
       ( 9.899495  1.6245536 1.849715  1.3926178))
   (linear-algebra-kernel:symmetric-cholesky-decomposition
    (make-array
     '(4 4)
     :initial-contents
     '(( 18.0  22.0  54.0  42.0)
       ( 22.0  70.0  86.0  62.0)
       ( 54.0  86.0 174.0 134.0)
       ( 42.0  62.0 134.0 106.0))))))

(define-test hermitian-cholesky-decomposition
  (:tag :kernel :cholesky)
  ;; 2x2
  (let ((*epsilon* (* 9 single-float-epsilon)))
    (assert-float-equal
     #2A((#C(1.4142135 0.0)       #C(0.7071068 -1.4142137))
         (#C(0.7071068 1.4142137) #C(0.7071066  0.0)))
     (linear-algebra-kernel:hermitian-cholesky-decomposition
      (make-array
       '(2 2) :initial-contents
       '((#C(2.0 0.0) #C(1.0 -2.0))
         (#C(1.0 2.0) #C(3.0  0.0)))))))
  ;; 3x3
  (let ((*epsilon* (* 3 single-float-epsilon)))
    (assert-float-equal
     #2A((#C(1.8193405  0.0)       #C( 0.69255865 -1.0992994) #C(0.75302017   -1.6489492))
         (#C(0.69255865 1.0992994) #C( 0.7361408   0.0)       #C(-0.032873448 -1.610834))
         (#C(0.75302017 1.6489492) #C(-0.032873448 1.610834)  #C(1.5060079     0.0)))
     (linear-algebra-kernel:hermitian-cholesky-decomposition
      (make-array
       '(3 3)
       :initial-contents
       '((#C(3.31 0.0) #C(1.26 -2.0) #C(1.37 -3.0))
         (#C(1.26 2.0) #C(2.23  0.0) #C(2.31 -1.5))
         (#C(1.37 3.0) #C(2.31  1.5) #C(8.15  0.0))))))))

(define-test root-free-symmetric-cholesky-decomposition
  (:tag :kernel :cholesky)
  ;; 2x2
  (assert-float-equal
   #2A((1.1 1.0909091) (1.0909091 0.8909091))
   (linear-algebra-kernel:root-free-symmetric-cholesky-decomposition
    (make-array '(2 2) :initial-contents '((1.1 1.2) (1.2 2.2)))))
  ;; 3x3
  (assert-float-equal
   #2A((1.15      1.0956522  1.1913043)
       (1.0956522 0.84947825 0.9522979)
       (1.1913043 0.9522979  0.90754557))
   (linear-algebra-kernel:root-free-symmetric-cholesky-decomposition
    (make-array
     '(3 3)
     :initial-contents
     '((1.15 1.26 1.37)
       (1.26 2.23 2.31)
       (1.37 2.31 3.31)))))
  ;; 3x3 from Wikipedia
  (assert-float-equal
   #2A(( 4.0  3.0 -4.0)
       ( 3.0  1.0  5.0)
       (-4.0  5.0  9.0))
   (linear-algebra-kernel:root-free-symmetric-cholesky-decomposition
    (make-array
     '(3 3)
     :initial-contents
     '((4 12 -16) (12 37 -43) (-16 -43 98))))))

(define-test root-free-hermitian-cholesky-decomposition
  (:tag :kernel :cholesky)
  ;; 2x2
  (assert-float-equal
   #2A((#C(2.0 0.0) #C(0.5 -1.0))
       (#C(0.5 1.0) #C(0.5  0.0)))
   (linear-algebra-kernel:root-free-hermitian-cholesky-decomposition
    (make-array
     '(2 2) :initial-contents
     '((#C(2.0 0.0) #C(1.0 -2.0))
       (#C(1.0 2.0) #C(3.0  0.0))))))
  ;; 3x3
  (let ((*epsilon* (* 3 single-float-epsilon)))
    (assert-float-equal
     #2A((#C(3.31       0.0)       #C( 0.38066468 -0.6042296) #C( 0.4138973   -0.9063445))
         (#C(0.38066468 0.6042296) #C( 0.54190326  0.0)       #C(-0.044656467 -2.1882146))
         (#C(0.4138973  0.9063445) #C(-0.044656467 2.1882146) #C( 2.2680602   -3.7252904E-9)))
     (linear-algebra-kernel:root-free-hermitian-cholesky-decomposition
      (make-array
       '(3 3)
       :initial-contents
       '((#C(3.31 0.0) #C(1.26 -2.0) #C(1.37 -3.0))
         (#C(1.26 2.0) #C(2.23  0.0) #C(2.31 -1.5))
         (#C(1.37 3.0) #C(2.31  1.5) #C(8.15  0.0))))))))

(define-test symmetric-cholesky-solver
  (:tag :kernel :cholesky)
  ;; 2x2
  (assert-float-equal
   #(3.2653065 -1.3265308)
   (linear-algebra-kernel:symmetric-cholesky-solver
    (make-array '(2 2) :initial-contents '((1.1 1.2) (1.2 2.2)))
    (make-array 2 :initial-contents '(2.0 1.0))))
  ;; 3x3
  (assert-float-equal
   #(3.5856622 -2.306286 0.79007966)
   (linear-algebra-kernel:symmetric-cholesky-solver
    (make-array
     '(3 3)
     :initial-contents
     '((1.15 1.26 1.37)
       (1.26 2.23 2.31)
       (1.37 2.31 3.31)))
    (make-array 3 :initial-contents '(2.3 1.2 2.2)))))

(define-test hermitian-cholesky-solver
  (:tag :kernel :cholesky)
  ;; 2x2
  (assert-float-equal
   #(#C(5.0 2.0) #C(0.0 -4.0))
   (linear-algebra-kernel:hermitian-cholesky-solver
    (make-array
     '(2 2)
     :initial-contents
     '((#C(2.0 0.0) #C(1.0 -2.0))
       (#C(1.0 2.0) #C(3.0  0.0))))
    (make-array 2 :initial-contents '(2.0 1.0))))
  ;; 3x3
  (let ((*epsilon* (* 5 single-float-epsilon)))
    (assert-float-equal
     #(#C( 3.5175734   3.4673646)
       #C( 3.3198433  -4.3366637)
       #C(-0.78414906 -1.2595192))
     (linear-algebra-kernel:hermitian-cholesky-solver
      (make-array
       '(3 3)
       :initial-contents
       '((#C(3.31 0.0) #C(1.26 -2.0) #C(1.37 -3.0))
         (#C(1.26 2.0) #C(2.23  0.0) #C(2.31 -1.5))
         (#C(1.37 3.0) #C(2.31  1.5) #C(8.15  0.0))))
      (make-array 3 :initial-contents '(2.3 1.2 2.2))))))

(define-test symmetric-cholesky-invert
  (:tag :kernel :cholesky)
  ;; 2x2
  (assert-float-equal
   #2A((2.2448979 -1.2244898) (-1.2244898 1.122449))
   (linear-algebra-kernel:symmetric-cholesky-invert
    (make-array '(2 2) :initial-contents '((1.1 1.2) (1.2 2.2)))))
  ;; 3x3
  (assert-float-equal
   #2A((2.3068395 -1.1345832 -0.16298579)
       (-1.1345832 2.1764503 -1.0493114)
       (-0.16298579 -1.0493114 1.101873))
   (linear-algebra-kernel:symmetric-cholesky-invert
    (make-array
     '(3 3)
     :initial-contents
     '((1.15 1.26 1.37)
       (1.26 2.23 2.31)
       (1.37 2.31 3.31))))))

(define-test hermitian-cholesky-invert
  (:tag :kernel :cholesky)
  ;; 2x2
  (assert-float-equal
   #2A((#C( 3.0 -0.0) #C(-1.0  2.0))
       (#C(-1.0 -2.0) #C( 2.0 -0.0)))
   (linear-algebra-kernel:hermitian-cholesky-invert
    (make-array
     '(2 2)
     :initial-contents
     '((#C(2.0 0.0) #C(1.0 -2.0))
       (#C(1.0 2.0) #C(3.0  0.0))))))
  ;; 3x3
  (let ((*epsilon* (* 3 single-float-epsilon)))
    (assert-float-equal
     #2A((#C( 2.602711    5.9604646E-8) #C(-0.64015717  2.808354)     #C(-0.7729426   0.04424542))
         (#C(-0.64015717 -2.808354)     #C(3.9574066   -3.7252904E-9) #C( 0.019689279 0.96479553))
         (#C(-0.7729426  -0.04424542)   #C(0.019689279 -0.96479553)   #C( 0.4409054  -7.241874E-10)))
     (linear-algebra-kernel:hermitian-cholesky-invert
      (make-array
       '(3 3)
       :initial-contents
       '((#C(3.31 0.0) #C(1.26 -2.0) #C(1.37 -3.0))
         (#C(1.26 2.0) #C(2.23  0.0) #C(2.31 -1.5))
         (#C(1.37 3.0) #C(2.31  1.5) #C(8.15  0.0))))))))
