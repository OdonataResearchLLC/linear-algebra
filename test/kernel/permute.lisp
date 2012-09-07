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

;;; Vector permutation

(define-test right-permute-vector
  (assert-float-equal
   #(3.3 4.4 1.1 5.5 2.2)
   (linear-algebra-kernel:right-permute-vector
    (vector 1.1 2.2 3.3 4.4 5.5)
    (vector 2 4 0 1 3)))
  (assert-float-equal
   #(3.3 4.4 2.2 5.5 1.1)
   (linear-algebra-kernel:right-permute-vector
    (vector 1.1 2.2 3.3 4.4 5.5)
    (vector 4 2 0 1 3))))

(define-test left-permute-vector
  (assert-float-equal
   #(3.3 5.5 1.1 2.2 4.4)
   (linear-algebra-kernel:left-permute-vector
    (vector 2 4 0 1 3)
    (vector 1.1 2.2 3.3 4.4 5.5)))
  (assert-float-equal
   #(5.5 3.3 1.1 2.2 4.4)
   (linear-algebra-kernel:left-permute-vector
    (vector 4 2 0 1 3)
    (vector 1.1 2.2 3.3 4.4 5.5))))

(define-test right-npermute-vector
  (let ((vect (vector 1.1 2.2 3.3 4.4 5.5))
        (pmat (vector 2 4 0 1 3)))
    (assert-eq
     vect
     (linear-algebra-kernel:right-npermute-vector vect pmat))
    (assert-float-equal #(3.3 4.4 1.1 5.5 2.2) vect))
  (let ((vect (vector 1.1 2.2 3.3 4.4 5.5))
        (pmat (vector 4 2 0 1 3)))
    (assert-eq
     vect
     (linear-algebra-kernel:right-npermute-vector vect pmat))
    (assert-float-equal #(3.3 4.4 2.2 5.5 1.1) vect)))

(define-test left-npermute-vector
  (let ((vect (vector 1.1 2.2 3.3 4.4 5.5))
        (pmat (vector 2 3 0 4 1)))
    (assert-eq
     vect
     (linear-algebra-kernel:left-npermute-vector pmat vect))
    (assert-float-equal #(3.3 5.5 1.1 2.2 4.4) vect))
  (let ((vect (vector 1.1 2.2 3.3 4.4 5.5))
        (pmat (vector 2 3 1 4 0)))
    (assert-eq
     vect
     (linear-algebra-kernel:left-npermute-vector pmat vect))
    (assert-float-equal #(5.5 3.3 1.1 2.2 4.4) vect)))

;;; Array permutation

(define-test right-permute-array
  (assert-float-equal
   #2A((1.2 1.3 1.0 1.4 1.1)
       (2.2 2.3 2.0 2.4 2.1)
       (3.2 3.3 3.0 3.4 3.1)
       (4.2 4.3 4.0 4.4 4.1)
       (5.2 5.3 5.0 5.4 5.1))
   (linear-algebra-kernel:right-permute
    #2A((1.0 1.1 1.2 1.3 1.4)
        (2.0 2.1 2.2 2.3 2.4)
        (3.0 3.1 3.2 3.3 3.4)
        (4.0 4.1 4.2 4.3 4.4)
        (5.0 5.1 5.2 5.3 5.4))
    (vector 2 4 0 1 3))))

(define-test left-permute-array
  (assert-float-equal
   #2A((3.0 3.1 3.2 3.3 3.4)
       (5.0 5.1 5.2 5.3 5.4)
       (1.0 1.1 1.2 1.3 1.4)
       (2.0 2.1 2.2 2.3 2.4)
       (4.0 4.1 4.2 4.3 4.4))
   (linear-algebra-kernel:left-permute
    (vector 2 4 0 1 3)
    #2A((1.0 1.1 1.2 1.3 1.4)
        (2.0 2.1 2.2 2.3 2.4)
        (3.0 3.1 3.2 3.3 3.4)
        (4.0 4.1 4.2 4.3 4.4)
        (5.0 5.1 5.2 5.3 5.4)))))

(define-test right-npermute-array
  (let ((array
         (make-array
          '(5 5) :initial-contents
          '((1.0 1.1 1.2 1.3 1.4)
            (2.0 2.1 2.2 2.3 2.4)
            (3.0 3.1 3.2 3.3 3.4)
            (4.0 4.1 4.2 4.3 4.4)
            (5.0 5.1 5.2 5.3 5.4)))))
    (assert-eq
     array (linear-algebra-kernel:right-npermute
            array (vector 2 4 0 1 3)))
    (assert-float-equal
     #2A((1.2 1.3 1.0 1.4 1.1)
         (2.2 2.3 2.0 2.4 2.1)
         (3.2 3.3 3.0 3.4 3.1)
         (4.2 4.3 4.0 4.4 4.1)
         (5.2 5.3 5.0 5.4 5.1))
     array)))

(define-test left-npermute-array
  (let ((array
         (make-array
          '(5 5) :initial-contents
          '((1.0 1.1 1.2 1.3 1.4)
            (2.0 2.1 2.2 2.3 2.4)
            (3.0 3.1 3.2 3.3 3.4)
            (4.0 4.1 4.2 4.3 4.4)
            (5.0 5.1 5.2 5.3 5.4)))))
    (assert-eq
     array (linear-algebra-kernel:left-npermute
            (vector 2 3 0 4 1) array))
    (assert-float-equal
     #2A((3.0 3.1 3.2 3.3 3.4)
         (5.0 5.1 5.2 5.3 5.4)
         (1.0 1.1 1.2 1.3 1.4)
         (2.0 2.1 2.2 2.3 2.4)
         (4.0 4.1 4.2 4.3 4.4))
     array)))
