#|

 Fundamental Common Lisp Vector Operations

 Copyright (c) 2009-2012, Odonata Research LLC
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

(defmethod sumsq ((data vector) &key (scale 0) (sumsq 1))
  "Return the scaling parameter and the sum of the squares of the
vector."
  (let ((size (length data))
        (abs-val))
    (dotimes (index size (values scale sumsq))
      (when (plusp (setf abs-val (abs (svref data index))))
        (if (< scale abs-val)
            (setf
             sumsq (1+ (* sumsq (expt (/ scale abs-val) 2)))
             scale abs-val)
            (incf sumsq (expt (/ (svref data index) scale) 2)))))))

(defmethod sump ((data vector) (p real) &key (scale 0) (sump 1))
  "Return the scaling parameter and the sum of the powers of p of the
vector."
  (let ((size (length data))
        (abs-val))
    (dotimes (index size (values scale sump))
      (when (plusp (setf abs-val (abs (svref data index))))
        (if (< scale abs-val)
            (setf
             sump (1+ (* sump (expt (/ scale abs-val) p)))
             scale abs-val)
            (incf sump (expt (/ (svref data index) scale) p)))))))

(defmethod %norm ((data vector) (measure (eql 1)))
  "Return the Taxicab norm of the list."
  (loop for element across data sum (abs element)))

(defmethod %norm ((data vector) (measure (eql 2)))
  "Return the Euclidean norm of the vector."
  (multiple-value-bind (scale sumsq)
      (sumsq (loop for val across data collect (abs val)))
    (* scale (sqrt sumsq))))

(defmethod %norm ((data vector) (measure integer))
  "Return the p-norm of the vector."
  (multiple-value-bind (scale sump)
      (sump (loop for val across data collect (abs val)) measure)
    (* scale (expt sump (/ measure)))))

(defmethod %norm ((data vector) (measure (eql :infinity)))
  "Return the infinity, or maximum, norm of vector."
  (loop for element across data maximize (abs element)))

(defmethod norm ((data vector) &key (measure 1))
  (%norm data measure))

(defmethod transpose ((data vector) &key conjugate)
  "Return a row vector."
  (map-into
   (make-array
    (length data)
    :element-type (array-element-type data))
   (if conjugate #'conjugate #'identity)
   data))

(defmethod ntranspose ((data vector) &key conjugate)
  "Return a row vector destructively."
  (if conjugate
      (map-into data #'conjugate data)
      data))

(defmethod permute ((data vector) (matrix permutation-matrix))
  "Return the permutation of the list."
  (if (= (length data) (matrix-row-dimension matrix))
      (right-permute data (contents matrix))
      (error
       "Vector(~D) and permutation matrix~A are incompatible."
       (length data) (matrix-dimensions matrix))))

(defmethod permute ((matrix permutation-matrix) (data vector))
  "Return the permutation of the list."
  (if (= (length data) (matrix-column-dimension matrix))
      (left-permute (contents matrix) data)
      (error
       "Permutation matrix~A and vector(~D) are incompatible."
       (matrix-dimensions matrix) (length data))))

(defmethod npermute ((data vector) (matrix permutation-matrix))
  "Destructively permute the vector."
  (if (= (length data) (matrix-row-dimension matrix))
      (right-npermute data (contents matrix))
      (error
       "Vector(~D) and permutation matrix~A are incompatible."
       (length data) (matrix-dimensions matrix))))

(defmethod npermute ((matrix permutation-matrix) (data vector))
  "Destructively permute the list."
  (if (= (length data) (matrix-column-dimension matrix))
      (left-npermute (contents (ntranspose matrix)) data)
      (error
       "Permutation matrix~A and vector(~D) are incompatible."
       (matrix-dimensions matrix) (length data))))

(defmethod scale ((scalar number) (data vector))
  "Return the vector scaled by scalar."
  (map-into
   (make-array (length data) :element-type
               (upgraded-array-element-type
                (type-of data)))
   (lambda (item) (* scalar item))
   data))

(defmethod nscale ((scalar number) (data vector))
  "Return the vector destructively scaled by scalar."
  (map-into data (lambda (x) (* scalar x)) data))

(defmethod add ((vector1 vector) (vector2 vector)
                &key scalar1 scalar2)
  "Return the addition of scalar1*vector1 with scalar2*vector2"
  (if (= (length vector1) (length vector2))
      (binary-operation :add vector1 vector2 scalar1 scalar2)
      (error "VECTOR1(~D) and VECTOR2(~D) are not of equal length."
             (length vector1) (length vector2))))

(defmethod nadd ((vector1 vector) (vector2 vector)
                 &key scalar1 scalar2)
  "Return the addition of scalar2*vector2 to scalar1*vector1."
  (if (= (length vector1) (length vector2))
      (binary-operation :nadd vector1 vector2 scalar1 scalar2)
      (error "VECTOR1(~D) and VECTOR2(~D) are not of equal length."
             (length vector1) (length vector2))))

(defmethod subtract ((vector1 vector) (vector2 vector)
                     &key scalar1 scalar2)
  "Return the subraction of scalar2*vector2 from scalar1*vector1."
  (if (= (length vector1) (length vector2))
      (binary-operation :subtract vector1 vector2 scalar1 scalar2)
      (error "VECTOR1(~D) and VECTOR2(~D) are not of equal length."
             (length vector1) (length vector2))))

(defmethod nsubtract ((vector1 vector) (vector2 vector)
                      &key scalar1 scalar2)
  "Return the subraction of scalar2*vector2 from scalar1*vector1."
  (if (= (length vector1) (length vector2))
      (binary-operation :nsubtract vector1 vector2 scalar1 scalar2)
      (error "VECTOR1(~D) and VECTOR2(~D) are not of equal length."
             (length vector1) (length vector2))))

(defmethod product ((vector1 vector) (vector2 vector)
                    &key scalar conjugate)
  "Return the dot product of vector1 and vector2."
  (if (= (length vector1) (length vector2))
      (binary-operation
       :inner-product vector1 vector2 scalar conjugate)
      (error "VECTOR1(~D) and VECTOR2(~D) are not of equal length."
             (length vector1) (length vector2))))
