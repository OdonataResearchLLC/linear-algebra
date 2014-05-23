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

;;; Taxicab norm

(define-test %norm-1-list
  (:tag :list :norm)
  (assert-rational-equal
   36 (linear-algebra::%norm
       '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) 1))
  (assert-rational-equal
   36 (linear-algebra:norm
       '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) 1))
  (assert-float-equal
   19.535658
   (linear-algebra::%norm
    '(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
      #C(-2 3) #C(-3 1) #C(-1 0))
    1))
  (assert-float-equal
   19.535658
   (linear-algebra:norm
    '(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
      #C(-2 3) #C(-3 1) #C(-1 0))
    1)))

;;; Euclidean norm

(define-test %norm-2-list
  (:tag :list :norm)
  (assert-float-equal
   12.083046
   (linear-algebra::%norm
    '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) 2))
  (assert-float-equal
   12.083046
   (linear-algebra:norm
    '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)
    2))
  (assert-float-equal
   8.0
   (linear-algebra::%norm
    '(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
      #C(-2 3) #C(-3 1) #C(-1 0)) 2))
  (assert-float-equal
   8.0
   (linear-algebra:norm
    '(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
      #C(-2 3) #C(-3 1) #C(-1 0))
    2)))

;;; P-norm

(define-test %norm-p-list
  (:tag :list :norm)
  (let ((data '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5))
        (zdata
         '(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
           #C(-2 3) #C(-3 1) #C(-1 0))))
    (assert-float-equal
     8.732892 (linear-algebra::%norm data 3))
    (assert-float-equal
     6.064035 (linear-algebra::%norm zdata 3))
    ;; norm
    (assert-float-equal
     8.732892 (linear-algebra:norm data 3))
    (assert-float-equal
     6.064035 (linear-algebra:norm zdata 3))))

;;; Infinity norm

(define-test %norm-infinity-list
  (:tag :list :norm)
  (assert-rational-equal
   6 (linear-algebra::%norm
      '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)
      :infinity))
  (assert-rational-equal
   6 (linear-algebra:norm
      '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)
      :infinity))
  (assert-float-equal
   4.0 (linear-algebra::%norm
        '(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
          #C(-2 3) #C(-3 1) #C(-1 0))
        :infinity))
  (assert-float-equal
   4.0 (linear-algebra:norm
        '(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
          #C(-2 3) #C(-3 1) #C(-1 0))
        :infinity)))

;;; List transpose

(define-test transpose-list
  (:tag :list :transpose)
  (assert-float-equal
   '(1.0 2.0 3.0 4.0 5.0)
   (linear-algebra:transpose
    '(1.0 2.0 3.0 4.0 5.0))))

(define-test ntranspose-list
  (:tag :list :ntranspose)
  (let ((data (list 1.0 2.0 3.0 4.0 5.0)))
    (assert-equal
     data (linear-algebra:ntranspose data))))

;;; List permutation

(define-test permute-list
  (:tag :list :permute)
  (let ((list (list 1.1 2.2 3.3 4.4 5.5))
        (pmat
         (linear-algebra:make-matrix
          5 5 :matrix-type
          'linear-algebra:permutation-matrix
          :initial-contents
          '((0 0 1 0 0)
            (0 0 0 0 1)
            (1 0 0 0 0)
            (0 1 0 0 0)
            (0 0 0 1 0)))))
    (assert-float-equal
     (list 3.3 4.4 1.1 5.5 2.2)
     (linear-algebra:permute list pmat))
    (assert-float-equal
     (list 3.3 5.5 1.1 2.2 4.4)
     (linear-algebra:permute pmat list))))

;;; List scale

(define-test scale-list
  (:tag :list :scale)
  (assert-float-equal
   '(2.0 4.0 6.0 8.0 10.0)
   (linear-algebra:scale 2.0 '(1.0 2.0 3.0 4.0 5.0)))
  (assert-float-equal
   '(#C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0) #C(4.0 4.0) #C(5.0 5.0))
   (linear-algebra:scale #C(1.0 1.0) '(1.0 2.0 3.0 4.0 5.0)))
  (assert-float-equal
   '(#C(2.0 2.0) #C(4.0 4.0) #C(6.0 6.0) #C(8.0 8.0) #C(10.0 10.0))
   (linear-algebra:scale 2.0
    '(#C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0) #C(4.0 4.0) #C(5.0 5.0))))
  (assert-float-equal
   #(#C(0.0 4.0) #C(0.0 8.0) #C(0.0 12.0) #C(0.0 16.0) #C(0.0 20.0))
   (linear-algebra:scale
    #C(2.0 2.0)
    '(#C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0) #C(4.0 4.0) #C(5.0 5.0)))))

(define-test nscale-list
  (:tag :list :nscale)
  (let ((list (list 1.0 2.0 3.0 4.0 5.0)))
    (assert-eq list (linear-algebra:nscale 2.0 list))
    (assert-float-equal '(2.0 4.0 6.0 8.0 10.0) list))
  (let ((list (list 1.0 2.0 3.0 4.0 5.0)))
    (assert-eq list (linear-algebra:nscale #C(1.0 1.0) list))
    (assert-float-equal
     '(#C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0) #C(4.0 4.0) #C(5.0 5.0))
     list))
  (let ((list (list #C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0)
                    #C(4.0 4.0) #C(5.0 5.0))))
    (assert-eq list (linear-algebra:nscale 2.0 list))
    (assert-float-equal
     '(#C(2.0 2.0) #C(4.0 4.0) #C(6.0 6.0) #C(8.0 8.0) #C(10.0 10.0))
     list))
  (let ((list (list #C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0)
                    #C(4.0 4.0) #C(5.0 5.0))))
    (assert-eq list (linear-algebra:nscale #C(2.0 2.0) list))
    (assert-float-equal
     '(#C(0.0 4.0) #C(0.0 8.0) #C(0.0 12.0) #C(0.0 16.0) #C(0.0 20.0))
     list)))

;;; List addition

(define-test add-list
  (:tag :list :add)
  ;; Real
  (let ((list1 '(1.1 2.2 3.3 4.4))
        (list2 '(1.1 2.2 3.3 4.4)))
    (assert-float-equal
     '(2.2 4.4 6.6 8.8) (linear-algebra:add list1 list2))
    (assert-float-equal
    '(3.3 6.6 9.9 13.2) (linear-algebra:add list1 list2 :scalar1 2.0))
    (assert-float-equal
     '(3.3 6.6 9.9 13.2) (linear-algebra:add list1 list2 :scalar2 2.0))
    (assert-float-equal
     '(4.4 8.8 13.2 17.6)
     (linear-algebra:add list1 list2 :scalar1 2.0 :scalar2 2.0)))
  ;; Complex
  (let ((list1 '(#C(1.1 2.2) #C(3.3 4.4)))
        (list2 '(#C(1.1 2.2) #C(3.3 4.4))))
    (assert-float-equal
     '(#C(2.2 4.4) #C(6.6 8.8)) (linear-algebra:add list1 list2))
    (assert-float-equal
     '(#C(3.3 6.6) #C(9.9 13.2))
     (linear-algebra:add list1 list2 :scalar1 2.0))
    (assert-float-equal
     '(#C(3.3 6.6) #C(9.9 13.2))
     (linear-algebra:add list1 list2 :scalar2 2.0))
    (assert-float-equal
     '(#C(4.4 8.8) #C(13.2 17.6))
     (linear-algebra:add list1 list2 :scalar1 2.0 :scalar2 2.0))))

;;; Destructive list addition

(define-test nadd-list
  (:tag :list :nadd)
  ;; Real
  (let ((list1 (list 1.1 2.2 3.3 4.4))
        (list2 (list 1.1 2.2 3.3 4.4)))
    (assert-eq list1 (linear-algebra:nadd list1 list2))
    (assert-float-equal '(2.2 4.4 6.6 8.8) list1)
    (assert-float-equal
     '(4.4 8.8 13.2 17.6)
     (linear-algebra:nadd list1 list2 :scalar2 2.0))
    (assert-float-equal
     '(9.9 19.8 29.7 39.6)
     (linear-algebra:nadd list1 list2 :scalar1 2.0))
    (assert-float-equal
     '(22.0 44.0 66.0 88.0)
     (linear-algebra:nadd list1 list2 :scalar1 2.0 :scalar2 2.0)))
  ;; Complex
  (let ((list1 (list #C(1.1 2.2) #C(3.3 4.4)))
        (list2 (list #C(1.1 2.2) #C(3.3 4.4))))
    (assert-eq list1 (linear-algebra:nadd list1 list2))
    (assert-float-equal '(#C(2.2 4.4) #C(6.6 8.8)) list1)
    (assert-float-equal
     '(#C(4.4 8.8) #C(13.2 17.6))
     (linear-algebra:nadd list1 list2 :scalar2 2.0))
    (assert-float-equal
     '(#C(9.9 19.8) #C(29.7 39.6))
     (linear-algebra:nadd list1 list2 :scalar1 2.0))
    (assert-float-equal
     '(#C(22.0 44.0) #C(66.0 88.0))
     (linear-algebra:nadd list1 list2 :scalar1 2.0 :scalar2 2.0))))

;;; List subtraction

(define-test subtract-list
  (:tag :list :subtract)
  ;; Real
  (let ((list1 '(1.1 2.2 3.3 4.4))
        (list2 '(1.1 2.2 3.3 4.4)))
    (assert-float-equal
     '(0.0 0.0 0.0 0.0)
     (linear-algebra:subtract list1 list2))
    (assert-float-equal
     '(1.1 2.2 3.3 4.4)
     (linear-algebra:subtract list1 list2 :scalar1 2.0))
    (assert-float-equal
     '(-1.1 -2.2 -3.3 -4.4)
     (linear-algebra:subtract list1 list2 :scalar2 2.0))
    (assert-float-equal
     '(0.0 0.0 0.0 0.0)
     (linear-algebra:subtract list1 list2 :scalar1 2.0 :scalar2 2.0)))
  ;; Complex
  (let ((list1 '(#C(1.1 2.2) #C(3.3 4.4)))
        (list2 '(#C(1.1 2.2) #C(3.3 4.4))))
    (assert-float-equal
     '(#C(0.0 0.0) #C(0.0 0.0))
     (linear-algebra:subtract list1 list2))
    (assert-float-equal
     '(#C(1.1 2.2) #C(3.3 4.4))
     (linear-algebra:subtract list1 list2 :scalar1 2.0))
    (assert-float-equal
     '(#C(-1.1 -2.2) #C(-3.3 -4.4))
     (linear-algebra:subtract list1 list2 :scalar2 2.0))
    (assert-float-equal
     '(#C(0.0 0.0) #C(0.0 0.0))
     (linear-algebra:subtract list1 list2 :scalar1 2.0 :scalar2 2.0))))

;;; Destructive list subtraction

(define-test nsubtract-list
  (:tag :list :nsubtract)
  ;; Real
  (let ((list1 (list 1.1 2.2 3.3 4.4))
        (list2 (list 1.1 2.2 3.3 4.4)))
    (assert-eq list1 (linear-algebra:nsubtract list1 list2))
    (assert-float-equal '(0.0 0.0 0.0 0.0) list1)
    (assert-float-equal
     '(-2.2 -4.4 -6.6 -8.8)
     (linear-algebra:nsubtract list1 list2 :scalar2 2.0))
    (assert-float-equal
     '(-5.5 -11.0 -16.5 -22.0)
     (linear-algebra:nsubtract list1 list2 :scalar1 2.0))
    (assert-float-equal
     '(-13.2 -26.4 -39.6 -52.8)
     (linear-algebra:nsubtract list1 list2 :scalar1 2.0 :scalar2 2.0)))
  ;; Complex
  (let ((list1 (list #C(1.1 2.2) #C(3.3 4.4)))
        (list2 (list #C(1.1 2.2) #C(3.3 4.4))))
    (assert-eq list1 (linear-algebra:nsubtract list1 list2))
    (assert-float-equal
     '(#C(0.0 0.0) #C(0.0 0.0)) list1)
    (assert-float-equal
     '(#C(-2.2 -4.4) #C(-6.6 -8.8))
     (linear-algebra:nsubtract list1 list2 :scalar2 2.0))
    (assert-float-equal
     '(#C(-5.5 -11.0) #C(-16.5 -22.0))
     (linear-algebra:nsubtract list1 list2 :scalar1 2.0))
    (assert-float-equal
     '(#C(-13.2 -26.4) #C(-39.6 -52.8))
     (linear-algebra:nsubtract list1 list2 :scalar1 2.0 :scalar2 2.0))))

;;; List dot product

(define-test product-list
  (:tag :list :product)
  ;; Real vectors
  (assert-rational-equal
   55 (linear-algebra:product '(1 2 3 4 5) '(1 2 3 4 5)))
  (assert-float-equal
   55F0
   (linear-algebra:product
    '(1.0 2.0 3.0 4.0 5.0) '(1.0 2.0 3.0 4.0 5.0)))
  (assert-float-equal
   55D0
   (linear-algebra:product
    '(1.0d0 2.0d0 3.0d0 4.0d0 5.0d0)
    '(1.0d0 2.0d0 3.0d0 4.0d0 5.0d0)))
  ;; Real vectors with conjugate keyword
  (assert-rational-equal
   55 (linear-algebra:product '(1 2 3 4 5) '(1 2 3 4 5)))
  ;; Complex vectors
  (assert-rational-equal
   #C(8 18)
   (linear-algebra:product
    '(#C(1 1) #C(2 1) #C(3 1))
    '(#C(1 2) #C(2 2) #C(3 2))))
  (assert-float-equal
   #C(8.0 18.0)
   (linear-algebra:product
    '(#C(1.0 1.0) #C(2.0 1.0) #C(3.0 1.0))
    '(#C(1.0 2.0) #C(2.0 2.0) #C(3.0 2.0))))
  (assert-float-equal
   #C(8.0d0 18.0d0)
   (linear-algebra:product
    '(#C(1.0d0 1.0d0) #C(2.0d0 1.0d0) #C(3.0d0 1.0d0))
    '(#C(1.0d0 2.0d0) #C(2.0d0 2.0d0) #C(3.0d0 2.0d0)))))
