#|

 Linear Algebra in Common Lisp

 Copyright (c) 2011-2012, Thomas M. Hermann
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

(in-package :linear-algebra-kernel)

;;; Squared sums

(defun lapy2 (x y)
  "Return the square root of |x|^2 + |y|^2."
  (let* ((abs-x (abs x))
         (abs-y (abs y))
         (w (max abs-x abs-y))
         (v/w (/ (min abs-x abs-y) w)))
    (* w (sqrt (+ 1 (* v/w v/w))))))

(defun lapy3 (x y z)
  "Return the square root of |x|^2 + |y|^2 + |z|^2."
  (let* ((abs-x (abs x))
         (abs-y (abs y))
         (abs-z (abs z))
         (w (max abs-x abs-y abs-z))
         (x/w (/ abs-x w))
         (y/w (/ abs-y w))
         (z/w (/ abs-z w)))
    (* w (sqrt (+ (* x/w x/w) (* y/w y/w) (* z/w z/w))))))

;;; Rotations

(defun givens-rotation (f g)
  "Return c,s,r defined from the Givens rotation."
  (cond
    ((zerop g)
     (values 1 0 f))
    ((zerop f)
     (values 0 (signum (conjugate g)) (abs g)))
    (t
     (let* ((abs-f (abs f))
            (abs-g (abs g))
            (sqrtfg (sqrt (+ (* abs-f abs-f) (* abs-g abs-g)))))
       (values
        (/ abs-f sqrtfg)
        (/ (* (signum f) (conjugate g)) sqrtfg)
        (* (signum f) sqrtfg))))))

(defun jacobi-rotation (x y z)
  "Return a, b, cos(theta) and sin(theta) terms from the Jacobi rotation."
  (let* ((yabs (abs y))
         (tau  (/ (- x z) 2.0 yabs))
         (tee  (/ (float-sign tau)
                  (+ (abs tau) (sqrt (+ 1.0 (expt tau 2))))))
         (cos-theta (/ (sqrt (+ 1.0 (expt tee 2))))) ; Invert sqrt
         (sin-theta (* cos-theta tee)))
    (values
     ;; a : first eigenvalue
     (+ (* cos-theta cos-theta x)
        (* 2.0 cos-theta sin-theta yabs)
        (* sin-theta sin-theta z))
     ;; b : second eigenvalue
     (+ (* sin-theta sin-theta x)
        (* -2.0 cos-theta sin-theta yabs)
        (* cos-theta cos-theta z))
     ;; Cosine theta
     cos-theta
     ;; Sine theta
     (* (conjugate (signum y)) sin-theta))))

(defun householder-reflection (alpha vector)
  "Return Beta, Tau and the Householder vector."
  (let* ((beta (- (float-sign
                   (realpart alpha)
                   (lapy2 alpha (norm-vector vector 2)))))
         (tau  (- 1 (/ alpha beta))))
    (values
     beta tau
     (dotimes (index (length vector) vector)
       (setf
        (aref vector index)
        (/ (aref vector index) alpha))))))

;;; Class and type utilities

(defun common-class-of (object1 object2 &optional
                        (default-class nil default-class-p))
  "Return the common class of the 2 objects or default-class."
  (let ((class1 (class-of object1))
        (class2 (class-of object2)))
    (cond
      ((eq class1 class2) class1)
      ((subtypep class1 class2) class2)
      ((subtypep class2 class1) class1)
      (default-class-p (find-class default-class))
      (t (error "No common or default class.")))))

(defun common-array-element-type (array1 array2)
  "Return the array type common to both arrays."
  (let ((type1 (array-element-type array1))
        (type2 (array-element-type array2)))
    (cond
     ((eq type1 type2) type1)
     ((subtypep type1 type2) type1)
     ((subtypep type2 type1) type2)
     ((and (subtypep type1 'number) (subtypep type2 'number))
      (upgraded-array-element-type
       (type-of (+ (coerce 1 type1) (coerce 1 type2)))))
     (t))))

;;; Equality predicates

;;; (COMPLEX-EQUAL number1 number2) => true or false
(defun complex-equal (complex1 complex2 &optional (epsilon *epsilon*))
  "Return true if both numbers are complex and equal."
  (cond
    ((or (typep complex1 '(complex float))
         (typep complex2 '(complex float)))
     (float-equal complex1 complex2 epsilon))
    ((or (typep complex1 '(complex integer))
         (typep complex2 '(complex integer)))
     (= complex1 complex2))
    (t (error "Arguments are not complex."))))

;;; (NUMBER-EQUAL number1 number2) => true or false
(defun number-equal (number1 number2 &optional (epsilon *epsilon*))
  "Return true if the numbers are equal using the appropriate
comparison."
  (cond
    ((or (floatp number1) (floatp number2))
     (float-equal number1 number2 epsilon))
    ((and (rationalp number1) (rationalp number2))
     (= number1 number2))
    ((or (typep number1 '(complex float))
         (typep number2 '(complex float)))
     (float-equal number1 number2 epsilon))
    ((and (typep number1 '(complex rational))
          (typep number2 '(complex rational)))
     (= number1 number2))
    (t (error "Non-numeric arguments."))))
