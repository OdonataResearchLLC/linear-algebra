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
