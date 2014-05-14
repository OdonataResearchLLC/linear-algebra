#|

  Fundamental Common Lisp Vector Operations

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

(in-package :linear-algebra)

(defmethod sumsq ((data vector))
  "Return the scaling parameter and the sum of the squares of the
vector."
  (sumsq-vector data))

(defmethod sump ((data vector) (p real))
  "Return the scaling parameter and the sum of the powers of p of the
vector."
  (sump-vector data p))

(defmethod norm ((data vector) &optional (measure 1))
  (norm-vector data measure))

(defmethod transpose ((data vector) &optional conjugate)
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

(defmethod ntranspose ((data vector) &optional conjugate)
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

(defmethod product
    ((vector1 vector) (vector2 vector) &optional scalar)
  "Return the dot product of vector1 and vector2."
  (if (= (length vector1) (length vector2))
      (inner-product-vector vector1 vector2 scalar)
      (error "VECTOR1(~D) and VECTOR2(~D) are not of equal length."
             (length vector1) (length vector2))))
