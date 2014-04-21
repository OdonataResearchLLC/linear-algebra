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
  (sumsq-vector data scale sumsq))

(defmethod sump ((data vector) (p real) &key (scale 0) (sump 1))
  "Return the scaling parameter and the sum of the powers of p of the
vector."
  (sump-vector data p scale sump))

(defmethod norm ((data vector) &key (measure 1))
  (norm-vector data measure))

(defmethod transpose ((data vector) &key conjugate)
  "Return a row vector."
  (let ((result
         (make-array
          (length data)
          :element-type (array-element-type data))))
    (if conjugate
        (dotimes (index (length data) result)
          (setf (aref result index) (aref data index)))
        (dotimes (index (length data) result)
          (setf (aref result index) (conjugate (aref data index)))))))

(defmethod ntranspose ((data vector) &key conjugate)
  "Return a row vector destructively."
  (if conjugate
      (dotimes (index (length data) data)
        (setf (aref data index) (conjugate (aref data index))))
      data))

(defmethod permute ((data vector) (matrix permutation-matrix))
  "Return the permutation of the list."
  (if (= (length data) (matrix-row-dimension matrix))
      (right-permute-vector data (contents matrix))
      (error
       "Vector(~D) and permutation matrix~A are incompatible."
       (length data) (matrix-dimensions matrix))))

(defmethod permute ((matrix permutation-matrix) (data vector))
  "Return the permutation of the list."
  (if (= (length data) (matrix-column-dimension matrix))
      (left-permute-vector (contents matrix) data)
      (error
       "Permutation matrix~A and vector(~D) are incompatible."
       (matrix-dimensions matrix) (length data))))

(defmethod scale ((scalar number) (data vector))
  "Return the vector scaled by scalar."
  (let ((result
         (make-array
          (length data)
          :element-type (array-element-type data))))
    (dotimes (index (length data) result)
      (setf (aref result index) (* scalar (aref data index))))))

(defmethod nscale ((scalar number) (data vector))
  "Return the vector destructively scaled by scalar."
  (dotimes (index (length data) data)
    (setf (aref data index) (* scalar (aref data index)))))

(defmethod add ((vector1 vector) (vector2 vector)
                &key scalar1 scalar2)
  "Return the addition of scalar1*vector1 with scalar2*vector2"
  (if (= (length vector1) (length vector2))
      (add-vector vector1 vector2 scalar1 scalar2)
      (error "VECTOR1(~D) and VECTOR2(~D) are not of equal length."
             (length vector1) (length vector2))))

(defmethod nadd ((vector1 vector) (vector2 vector)
                 &key scalar1 scalar2)
  "Return the addition of scalar2*vector2 to scalar1*vector1."
  (if (= (length vector1) (length vector2))
      (nadd-vector vector1 vector2 scalar1 scalar2)
      (error "VECTOR1(~D) and VECTOR2(~D) are not of equal length."
             (length vector1) (length vector2))))

(defmethod subtract ((vector1 vector) (vector2 vector)
                     &key scalar1 scalar2)
  "Return the subraction of scalar2*vector2 from scalar1*vector1."
  (if (= (length vector1) (length vector2))
      (subtract-vector vector1 vector2 scalar1 scalar2)
      (error "VECTOR1(~D) and VECTOR2(~D) are not of equal length."
             (length vector1) (length vector2))))

(defmethod nsubtract ((vector1 vector) (vector2 vector)
                      &key scalar1 scalar2)
  "Return the subraction of scalar2*vector2 from scalar1*vector1."
  (if (= (length vector1) (length vector2))
      (nsubtract-vector vector1 vector2 scalar1 scalar2)
      (error "VECTOR1(~D) and VECTOR2(~D) are not of equal length."
             (length vector1) (length vector2))))

(defmethod product ((vector1 vector) (vector2 vector)
                    &key scalar conjugate)
  "Return the dot product of vector1 and vector2."
  (if (= (length vector1) (length vector2))
      (inner-product-vector vector1 vector2 scalar conjugate)
      (error "VECTOR1(~D) and VECTOR2(~D) are not of equal length."
             (length vector1) (length vector2))))
