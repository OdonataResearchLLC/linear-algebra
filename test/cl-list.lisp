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

(define-test sumsq-list
  ;; Real
  (let ((data '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)))
    (multiple-value-bind (scale sumsq)
        (linear-algebra:sumsq data)
      (assert-rational-equal 6 scale)
      (assert-rational-equal 73/18 sumsq)))
  ;; Complex
  (let ((data '(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
                #C(-2 3) #C(-3 1) #C(-1 0))))
    (multiple-value-bind (scale sumsq)
        (linear-algebra:sumsq data)
      (assert-float-equal 4.0 scale)
      (assert-float-equal #C(2.75 -1.125) sumsq))))

(define-test sump-list
  ;; Real
  (let ((data '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)))
    (multiple-value-bind (scale sump)
        (linear-algebra:sump data 2)
      (assert-rational-equal 6 scale)
      (assert-rational-equal 73/18 sump))
    (multiple-value-bind (scale sump)
        (linear-algebra:sump data 3)
      (assert-rational-equal 6 scale)
      (assert-rational-equal 1 sump)))
  ;; Complex
  (let ((data '(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
                #C(-2 3) #C(-3 1) #C(-1 0))))
    (multiple-value-bind (scale sump)
        (linear-algebra:sump data 2)
      (assert-float-equal 4.0 scale)
      (assert-float-equal #C(2.75 -1.125) sump))
    (multiple-value-bind (scale sump)
        (linear-algebra:sump data 3)
      (assert-float-equal 4.0 scale)
      (assert-float-equal #C(2.6639833 0.54687494) sump))))

;;; Taxicab norm

(define-test %norm-1-list
  (assert-rational-equal
   36 (linear-algebra::%norm
       '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) 1))
  (assert-rational-equal
   36 (linear-algebra:norm
       '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) :measure 1))
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
    :measure 1)))

;;; Euclidean norm

(define-test %norm-2-list
  (assert-float-equal
   12.083046
   (linear-algebra::%norm
    '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) 2))
  (assert-float-equal
   12.083046
   (linear-algebra:norm
    '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)
    :measure 2))
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
    :measure 2)))

;;; P-norm

(define-test %norm-p-list
  (let ((data '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5))
        (zdata '(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
                 #C(-2 3) #C(-3 1) #C(-1 0))))
    (assert-float-equal
     8.732892 (linear-algebra::%norm data 3))
    (assert-float-equal
     6.064035 (linear-algebra::%norm zdata 3))
    ;; norm
    (assert-float-equal
     8.732892 (linear-algebra:norm data :measure 3))
    (assert-float-equal
     6.064035 (linear-algebra:norm zdata :measure 3))))

;;; Infinity norm

(define-test %norm-infinity-list
  (assert-rational-equal
   6 (linear-algebra::%norm
      '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)
      :infinity))
  (assert-rational-equal
   6 (linear-algebra:norm
      '(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)
      :measure :infinity))
  (assert-float-equal
   4.0 (linear-algebra::%norm
        '(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
          #C(-2 3) #C(-3 1) #C(-1 0))
        :infinity))
  (assert-float-equal
   4.0 (linear-algebra:norm
        '(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
          #C(-2 3) #C(-3 1) #C(-1 0))
        :measure :infinity)))

;;; List transpose

(define-test transpose-list
  (assert-float-equal
   '(1.0 2.0 3.0 4.0 5.0)
   (linear-algebra:transpose
    '(1.0 2.0 3.0 4.0 5.0))))

(define-test ntranspose-list
  (let ((data (list 1.0 2.0 3.0 4.0 5.0)))
    (assert-equal
     data (linear-algebra:ntranspose data))))

;;; List permutation

(define-test permute-list
  (let ((list (list 1.1 2.2 3.3 4.4 5.5))
        (pmat (linear-algebra:make-matrix
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

(define-test npermute-list
  (let ((list (list 1.1 2.2 3.3 4.4 5.5))
        (pmat
         (linear-algebra:make-matrix
          5 5
          :matrix-type
          'linear-algebra:permutation-matrix
          :initial-contents
          '((0 0 1 0 0)
            (0 0 0 0 1)
            (1 0 0 0 0)
            (0 1 0 0 0)
            (0 0 0 1 0)))))
    (assert-eq list (linear-algebra:npermute list pmat))
    (assert-float-equal #(3.3 4.4 1.1 5.5 2.2) list))
  (let ((list (list 1.1 2.2 3.3 4.4 5.5))
        (pmat (linear-algebra:make-matrix
               5 5 :matrix-type
               'linear-algebra:permutation-matrix
               :initial-contents
               '((0 0 0 0 1)
                 (0 0 1 0 0)
                 (1 0 0 0 0)
                 (0 1 0 0 0)
                 (0 0 0 1 0)))))
    (assert-eq list (linear-algebra:npermute list pmat))
    (assert-float-equal #(3.3 4.4 2.2 5.5 1.1) list))
  (let ((list (list 1.1 2.2 3.3 4.4 5.5))
        (pmat (linear-algebra:make-matrix
               5 5 :matrix-type
               'linear-algebra:permutation-matrix
               :initial-contents
               '((0 0 1 0 0)
                 (0 0 0 0 1)
                 (1 0 0 0 0)
                 (0 1 0 0 0)
                 (0 0 0 1 0)))))
    (assert-eq list (linear-algebra:npermute pmat list))
    (assert-float-equal #(3.3 5.5 1.1 2.2 4.4) list))
  (let ((list (list 1.1 2.2 3.3 4.4 5.5))
        (pmat (linear-algebra:make-matrix
               5 5 :matrix-type
               'linear-algebra:permutation-matrix
               :initial-contents
               '((0 0 0 0 1)
                 (0 0 1 0 0)
                 (1 0 0 0 0)
                 (0 1 0 0 0)
                 (0 0 0 1 0)))))
    (assert-eq list (linear-algebra:npermute pmat list))
    (assert-float-equal #(5.5 3.3 1.1 2.2 4.4) list)))

;;; List scale

(define-test scale-list
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
  ;; Real
  (let ((list1 '(1.1 2.2 3.3 4.4))
        (list2 '(1.1 2.2 3.3 4.4)))
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
  (let ((list1 '(#C(1.1 2.2) #C(3.3 4.4)))
        (list2 '(#C(1.1 2.2) #C(3.3 4.4))))
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
  ;; Real
  (let ((list1 '(1.1 2.2 3.3 4.4))
        (list2 '(1.1 2.2 3.3 4.4)))
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
  (let ((list1 '(#C(1.1 2.2) #C(3.3 4.4)))
        (list2 '(#C(1.1 2.2) #C(3.3 4.4))))
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

(define-test dot-product-list
  ;; Real vectors
  (assert-rational-equal
   55 (linear-algebra:product '(1 2 3 4 5) '(1 2 3 4 5)))
  (assert-float-equal
   55F0 (linear-algebra:product
         '(1.0 2.0 3.0 4.0 5.0) '(1.0 2.0 3.0 4.0 5.0)))
  (assert-float-equal
   55D0 (linear-algebra:product
         '(1.0d0 2.0d0 3.0d0 4.0d0 5.0d0)
         '(1.0d0 2.0d0 3.0d0 4.0d0 5.0d0)))
  ;; Real vectors with conjugate keyword
  (assert-rational-equal
   55 (linear-algebra:product
       '(1 2 3 4 5) '(1 2 3 4 5) :conjugate t))
  ;; Complex vectors
  (assert-rational-equal
   #C(8 18) (linear-algebra:product
             '(#C(1 1) #C(2 1) #C(3 1))
             '(#C(1 2) #C(2 2) #C(3 2))))
  (assert-float-equal
   #C(8.0 18.0) (linear-algebra:product
                 '(#C(1.0 1.0) #C(2.0 1.0) #C(3.0 1.0))
                 '(#C(1.0 2.0) #C(2.0 2.0) #C(3.0 2.0))))
  (assert-float-equal
   #C(8.0d0 18.0d0)
   (linear-algebra:product
    '(#C(1.0d0 1.0d0) #C(2.0d0 1.0d0) #C(3.0d0 1.0d0))
    '(#C(1.0d0 2.0d0) #C(2.0d0 2.0d0) #C(3.0d0 2.0d0))))
  ;; Complex conjugate
  (assert-rational-equal
   #C(20 6)
   (linear-algebra:product
    '(#C(1 1) #C(2 1) #C(3 1)) '(#C(1 2) #C(2 2) #C(3 2))
    :conjugate t))
  (assert-float-equal
   #C(20.0 6.0)
   (linear-algebra:product
    '(#C(1.0 1.0) #C(2.0 1.0) #C(3.0 1.0))
    '(#C(1.0 2.0) #C(2.0 2.0) #C(3.0 2.0))
    :conjugate t))
  (assert-float-equal
   #C(20D0 6.0D0)
   (linear-algebra:product
    '(#C(1.0d0 1.0d0) #C(2.0d0 1.0d0) #C(3.0d0 1.0d0))
    '(#C(1.0d0 2.0d0) #C(2.0d0 2.0d0) #C(3.0d0 2.0d0))
    :conjugate t)))
