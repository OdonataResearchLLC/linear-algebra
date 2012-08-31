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

(define-test sumsq-vector
  ;; Real
  (let ((data #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)))
    (multiple-value-bind (scale sumsq)
        (linear-algebra:sumsq data)
      (assert-rational-equal 6 scale)
      (assert-rational-equal 73/18 sumsq)))
  ;; Complex
  (let ((data #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
                #C(-2 3) #C(-3 1) #C(-1 0))))
    (multiple-value-bind (scale sumsq)
        (linear-algebra:sumsq data)
      (assert-float-equal 4.0 scale)
      (assert-float-equal #C(2.75 -1.125) sumsq))))

(define-test sump-vector
  ;; Real
  (let ((data #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)))
    (multiple-value-bind (scale sump)
        (linear-algebra:sump data 2)
      (assert-rational-equal 6 scale)
      (assert-rational-equal 73/18 sump))
    (multiple-value-bind (scale sump)
        (linear-algebra:sump data 3)
      (assert-rational-equal 6 scale)
      (assert-rational-equal 1 sump)))
  ;; Complex
  (let ((data #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
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

(define-test %norm-1-vector
  (assert-rational-equal
   36 (linear-algebra::%norm
       #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) 1))
  (assert-rational-equal
   36 (linear-algebra:norm
       #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) :measure 1))
  (assert-float-equal
   19.535658
   (linear-algebra::%norm
    #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
      #C(-2 3) #C(-3 1) #C(-1 0))
    1))
  (assert-float-equal
   19.535658
   (linear-algebra:norm
    #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
      #C(-2 3) #C(-3 1) #C(-1 0))
    :measure 1)))

;;; Euclidean norm

(define-test %norm-2-vector
  (assert-float-equal
   12.083046
   (linear-algebra::%norm
    #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5) 2))
  (assert-float-equal
   12.083046
   (linear-algebra:norm
    #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)
    :measure 2))
  (assert-float-equal
   8.0
   (linear-algebra::%norm
    #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
      #C(-2 3) #C(-3 1) #C(-1 0)) 2))
  (assert-float-equal
   8.0
   (linear-algebra:norm
    #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
      #C(-2 3) #C(-3 1) #C(-1 0))
    :measure 2)))

;;; P-norm

(define-test %norm-p-vector
  (let ((data #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5))
        (zdata #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
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

(define-test %norm-infinity-vector
  (assert-rational-equal
   6 (linear-algebra::%norm
      #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)
      :infinity))
  (assert-rational-equal
   6 (linear-algebra:norm
      #(-6 -5 -4 -3 -2 -1 0 1 2 3 4 5)
      :measure :infinity))
  (assert-float-equal
   4.0 (linear-algebra::%norm
        #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
          #C(-2 3) #C(-3 1) #C(-1 0))
        :infinity))
  (assert-float-equal
   4.0 (linear-algebra:norm
        #(#C(1 0) #C(3 1) #C(2 3) #C(0 4)
          #C(-2 3) #C(-3 1) #C(-1 0))
        :measure :infinity)))

;;; Vector transpose

(define-test transpose-vector
  (assert-float-equal
   #(1.0 2.0 3.0 4.0 5.0)
   (linear-algebra:transpose
    #(1.0 2.0 3.0 4.0 5.0))))

(define-test ntranspose-vector
  (let ((data (vector 1.0 2.0 3.0 4.0 5.0)))
    (assert-equal
     data (linear-algebra:ntranspose data))))

;;; Vector permutation

(define-test permute-vector
  (let ((list (vector 1.1 2.2 3.3 4.4 5.5))
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
     (vector 3.3 4.4 1.1 5.5 2.2)
     (linear-algebra:permute list pmat))
    (assert-float-equal
     (vector 3.3 5.5 1.1 2.2 4.4)
     (linear-algebra:permute pmat list))))

(define-test npermute-vector
  (let ((vec1 (vector 1.1 2.2 3.3 4.4 5.5))
        (vec2 (vector 1.1 2.2 3.3 4.4 5.5))
        (pmat (linear-algebra:make-matrix
               5 5 :matrix-type
               'linear-algebra:permutation-matrix
               :initial-contents
               '((0 0 0 0 1)
                 (0 0 1 0 0)
                 (1 0 0 0 0)
                 (0 1 0 0 0)
                 (0 0 0 1 0)))))
    (assert-eq vec1 (linear-algebra:npermute vec1 pmat))
    (assert-float-equal #(3.3 4.4 2.2 5.5 1.1) vec1)
    (assert-eq vec2 (linear-algebra:npermute pmat vec2))
    (assert-float-equal #(5.5 3.3 1.1 2.2 4.4) vec2)))

;;; Vector scale

(define-test scale-vector
  (assert-float-equal
   #(2.0 4.0 6.0 8.0 10.0)
   (linear-algebra:scale 2.0 #(1.0 2.0 3.0 4.0 5.0)))
  (assert-float-equal
   #(#C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0) #C(4.0 4.0) #C(5.0 5.0))
   (linear-algebra:scale #C(1.0 1.0) #(1.0 2.0 3.0 4.0 5.0)))
  (assert-float-equal
   #(#C(2.0 2.0) #C(4.0 4.0) #C(6.0 6.0) #C(8.0 8.0) #C(10.0 10.0))
   (linear-algebra:scale 2.0
    #(#C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0) #C(4.0 4.0) #C(5.0 5.0))))
  (assert-float-equal
   #(#C(0.0 4.0) #C(0.0 8.0) #C(0.0 12.0) #C(0.0 16.0) #C(0.0 20.0))
   (linear-algebra:scale
    #C(2.0 2.0)
    #(#C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0) #C(4.0 4.0) #C(5.0 5.0)))))

(define-test nscale-vector
  (let ((list (list 1.0 2.0 3.0 4.0 5.0)))
    (assert-eq list (linear-algebra:nscale 2.0 list))
    (assert-float-equal #(2.0 4.0 6.0 8.0 10.0) list))
  (let ((list (list 1.0 2.0 3.0 4.0 5.0)))
    (assert-eq list (linear-algebra:nscale #C(1.0 1.0) list))
    (assert-float-equal
     #(#C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0) #C(4.0 4.0) #C(5.0 5.0))
     list))
  (let ((list (list #C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0)
                    #C(4.0 4.0) #C(5.0 5.0))))
    (assert-eq list (linear-algebra:nscale 2.0 list))
    (assert-float-equal
     #(#C(2.0 2.0) #C(4.0 4.0) #C(6.0 6.0) #C(8.0 8.0) #C(10.0 10.0))
     list))
  (let ((list (list #C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0)
                    #C(4.0 4.0) #C(5.0 5.0))))
    (assert-eq list (linear-algebra:nscale #C(2.0 2.0) list))
    (assert-float-equal
     #(#C(0.0 4.0) #C(0.0 8.0) #C(0.0 12.0) #C(0.0 16.0) #C(0.0 20.0))
     list)))

;;; Vector addition

(define-test add-vector
  ;; Real
  (let ((vector1 #(1.1 2.2 3.3 4.4))
        (vector2 #(1.1 2.2 3.3 4.4)))
    (assert-float-equal
     #(2.2 4.4 6.6 8.8) (linear-algebra:add vector1 vector2))
    (assert-float-equal
    #(3.3 6.6 9.9 13.2) (linear-algebra:add vector1 vector2 :scalar1 2.0))
    (assert-float-equal
     #(3.3 6.6 9.9 13.2) (linear-algebra:add vector1 vector2 :scalar2 2.0))
    (assert-float-equal
     #(4.4 8.8 13.2 17.6)
     (linear-algebra:add vector1 vector2 :scalar1 2.0 :scalar2 2.0)))
  ;; Complex
  (let ((vector1 #(#C(1.1 2.2) #C(3.3 4.4)))
        (vector2 #(#C(1.1 2.2) #C(3.3 4.4))))
    (assert-float-equal
     #(#C(2.2 4.4) #C(6.6 8.8)) (linear-algebra:add vector1 vector2))
    (assert-float-equal
     #(#C(3.3 6.6) #C(9.9 13.2))
     (linear-algebra:add vector1 vector2 :scalar1 2.0))
    (assert-float-equal
     #(#C(3.3 6.6) #C(9.9 13.2))
     (linear-algebra:add vector1 vector2 :scalar2 2.0))
    (assert-float-equal
     #(#C(4.4 8.8) #C(13.2 17.6))
     (linear-algebra:add vector1 vector2 :scalar1 2.0 :scalar2 2.0))))

;;; Destructive vector addition

(define-test nadd-vector
  ;; Real
  (let ((vector1 #(1.1 2.2 3.3 4.4))
        (vector2 #(1.1 2.2 3.3 4.4)))
    (assert-eq vector1 (linear-algebra:nadd vector1 vector2))
    (assert-float-equal #(2.2 4.4 6.6 8.8) vector1)
    (assert-float-equal
     #(4.4 8.8 13.2 17.6)
     (linear-algebra:nadd vector1 vector2 :scalar2 2.0))
    (assert-float-equal
     #(9.9 19.8 29.7 39.6)
     (linear-algebra:nadd vector1 vector2 :scalar1 2.0))
    (assert-float-equal
     #(22.0 44.0 66.0 88.0)
     (linear-algebra:nadd vector1 vector2 :scalar1 2.0 :scalar2 2.0)))
  ;; Complex
  (let ((vector1 #(#C(1.1 2.2) #C(3.3 4.4)))
        (vector2 #(#C(1.1 2.2) #C(3.3 4.4))))
    (assert-eq vector1 (linear-algebra:nadd vector1 vector2))
    (assert-float-equal #(#C(2.2 4.4) #C(6.6 8.8)) vector1)
    (assert-float-equal
     #(#C(4.4 8.8) #C(13.2 17.6))
     (linear-algebra:nadd vector1 vector2 :scalar2 2.0))
    (assert-float-equal
     #(#C(9.9 19.8) #C(29.7 39.6))
     (linear-algebra:nadd vector1 vector2 :scalar1 2.0))
    (assert-float-equal
     #(#C(22.0 44.0) #C(66.0 88.0))
     (linear-algebra:nadd vector1 vector2 :scalar1 2.0 :scalar2 2.0))))

;;; Vector subtraction

(define-test subtract-vector
  ;; Real
  (let ((vector1 #(1.1 2.2 3.3 4.4))
        (vector2 #(1.1 2.2 3.3 4.4)))
    (assert-float-equal
     #(0.0 0.0 0.0 0.0)
     (linear-algebra:subtract vector1 vector2))
    (assert-float-equal
     #(1.1 2.2 3.3 4.4)
     (linear-algebra:subtract vector1 vector2 :scalar1 2.0))
    (assert-float-equal
     #(-1.1 -2.2 -3.3 -4.4)
     (linear-algebra:subtract vector1 vector2 :scalar2 2.0))
    (assert-float-equal
     #(0.0 0.0 0.0 0.0)
     (linear-algebra:subtract
      vector1 vector2 :scalar1 2.0 :scalar2 2.0)))
  ;; Complex
  (let ((vector1 #(#C(1.1 2.2) #C(3.3 4.4)))
        (vector2 #(#C(1.1 2.2) #C(3.3 4.4))))
    (assert-float-equal
     #(#C(0.0 0.0) #C(0.0 0.0))
     (linear-algebra:subtract vector1 vector2))
    (assert-float-equal
     #(#C(1.1 2.2) #C(3.3 4.4))
     (linear-algebra:subtract vector1 vector2 :scalar1 2.0))
    (assert-float-equal
     #(#C(-1.1 -2.2) #C(-3.3 -4.4))
     (linear-algebra:subtract vector1 vector2 :scalar2 2.0))
    (assert-float-equal
     #(#C(0.0 0.0) #C(0.0 0.0))
     (linear-algebra:subtract
      vector1 vector2 :scalar1 2.0 :scalar2 2.0))))

;;; Destructive vector subtraction

(define-test nsubtract-vector
  ;; Real
  (let ((vector1 (vector 1.1 2.2 3.3 4.4))
        (vector2 (vector 1.1 2.2 3.3 4.4)))
    (assert-eq vector1 (linear-algebra:nsubtract vector1 vector2))
    (assert-float-equal #(0.0 0.0 0.0 0.0) vector1)
    (assert-float-equal
     #(-2.2 -4.4 -6.6 -8.8)
     (linear-algebra:nsubtract vector1 vector2 :scalar2 2.0))
    (assert-float-equal
     #(-5.5 -11.0 -16.5 -22.0)
     (linear-algebra:nsubtract vector1 vector2 :scalar1 2.0))
    (assert-float-equal
     #(-13.2 -26.4 -39.6 -52.8)
     (linear-algebra:nsubtract vector1 vector2 :scalar1 2.0 :scalar2 2.0)))
  ;; Complex
  (let ((vector1 (vector #C(1.1 2.2) #C(3.3 4.4)))
        (vector2 (vector #C(1.1 2.2) #C(3.3 4.4))))
    (assert-eq vector1 (linear-algebra:nsubtract vector1 vector2))
    (assert-float-equal
     #(#C(0.0 0.0) #C(0.0 0.0)) vector1)
    (assert-float-equal
     #(#C(-2.2 -4.4) #C(-6.6 -8.8))
     (linear-algebra:nsubtract vector1 vector2 :scalar2 2.0))
    (assert-float-equal
     #(#C(-5.5 -11.0) #C(-16.5 -22.0))
     (linear-algebra:nsubtract vector1 vector2 :scalar1 2.0))
    (assert-float-equal
     #(#C(-13.2 -26.4) #C(-39.6 -52.8))
     (linear-algebra:nsubtract vector1 vector2 :scalar1 2.0 :scalar2 2.0))))
