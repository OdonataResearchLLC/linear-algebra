#|

 Linear Algebra in Common Lisp

 Copyright (c) 2011, Thomas M. Hermann
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

(in-package :linear-algebra)

;;; Vector interface operations

(defgeneric initialize-vector (vector data size &optional element-type)
  (:documentation
   "Initialize the vector with data."))

(defun make-vector (size &key
                    (element-type t)
                    (vector-type 'column-vector)
                    (initial-element nil initial-element-p)
                    (initial-contents nil initial-contents-p))
  "Create a 1D numeric array to represent a numeric vector."
  (let ((new-vector (make-instance vector-type)))
    (cond
      ((and initial-element-p initial-contents-p)
       (error "Cannot specify both :INITIAL-ELEMENT and :INITIAL-CONTENTS."))
      (initial-contents-p
       (initialize-vector new-vector initial-contents size element-type))
      (initial-element-p
       (initialize-vector new-vector initial-element size element-type))
      (t
       (initialize-vector
        new-vector (coerce 0 element-type) size element-type)))))

(defgeneric vector-in-bounds-p (vector index)
  (:documentation
   "Return true if index does not exceed the dimensions of vector."))

(defgeneric vector-element-type (vector)
  (:documentation
   "Return the element type of vector."))

(defgeneric vector-length (vector)
  (:documentation
   "Return the length of the vector."))

(defgeneric vref (vector index)
  (:documentation
   "Return the element of vector at index."))

(defgeneric (setf vref) (data vector index)
  (:documentation
   "Set the element of vector at index to data."))

(defgeneric copy-vector (vector)
  (:documentation
   "Return a copy of the vector."))

(defgeneric subvector (vector start &optional end)
  (:documentation
   "Return a new vector that is a subvector of the vector."))

(defgeneric (setf subvector) (subvector vector start &optional end)
  (:documentation
   "Set the subvector of the vector."))

(defgeneric replace-vector (vector1 vector2 &key
                            start1 end1 start2 end2)
  (:documentation
   "Destructively replace the elements of vector1 with vector2."))

;;; Vector iteration operations
(defgeneric map-vector (result-type function first-vector &rest more-vectors)
  (:documentation
   "Calls function on successive sets of vector objects."))

(defgeneric map-into-vector (result-vector function &rest vectors)
  (:documentation
   "Destructively modifies the result vector with the result of
applying the function to each element of the vectors."))

(defmacro dovector ((element vector &optional result)
                    &body body)
  "Iterate over vector returning result."
  (let ((pos (gensym "POS-"))
        (end (gensym "END-")))
    `(let ((,end (vector-length ,vector))
           (,element nil))
      (dotimes (,pos ,end ,result)
        (setf ,element (vref ,vector ,pos))
        ,@body))))

;;; Vector transformations
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
                   (lapy2 alpha (norm vector :measure 2)))))
         (tau  (- 1 (/ alpha beta))))
    (values beta tau (scale (/ alpha) vector))))

(defgeneric apply-rotation (vector1 vector2 cc ss)
  (:documentation
   "Return the plane rotations of vector1 and vector2 by cc and ss."))

(defgeneric napply-rotation (vector1 vector2 cc ss)
  (:documentation
   "Return the plane rotations of vector1 and vector2 by cc and ss."))

