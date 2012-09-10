#|

 Fundamental Array Operations

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

(defmethod sumsq ((data array) &key (scale 0) (sumsq 1))
  "Return the scaling parameter and the sum of the squares of the
array."
  (destructuring-bind (numrows numcols) (array-dimensions data)
    (let ((abs-val 0))
      (dotimes (i0 numrows (values scale sumsq))
        (dotimes (i1 numcols)
          (when (plusp (setf abs-val (abs (aref data i0 i1))))
            (if (< scale abs-val)
                (setf
                 sumsq (1+ (* sumsq (expt (/ scale abs-val) 2)))
                 scale abs-val)
                (incf sumsq (expt (/ abs-val scale) 2)))))))))

(defmethod sump ((data array) (p number) &key (scale 0) (sump 1))
  "Return the scaling parameter and the sum of the P powers of the matrix."
  (unless (plusp p) (error "The power(~A) must be positive." p))
  (destructuring-bind (numrows numcols) (array-dimensions data)
    (let ((abs-val 0))
      (dotimes (i0 numrows (values scale sump))
        (dotimes (i1 numcols)
          (when (plusp (setf abs-val (abs (aref data i0 i1))))
            (if (< scale abs-val)
                (setf
                 sump (1+ (* sump (expt (/ scale abs-val) p)))
                 scale abs-val)
                (incf sump (expt (/ (aref data i0 i1) scale) p)))))))))

(defmethod %norm ((data array) (measure (eql 1)))
  "Return the 1 norm of the array."
  (destructuring-bind (numrows numcols) (array-dimensions data)
    (let ((norm 0)
          (sum 0))
      (dotimes (i1 numcols norm)
        (setf sum 0)
        (dotimes (i0 numrows)
          (incf sum (abs (aref data i0 i1))))
        (setf norm (max sum norm))))))

(defmethod %norm ((data array) (measure (eql :max)))
  "Return the max norm of the array."
  (destructuring-bind (numrows numcols) (array-dimensions data)
    (let ((norm 0))
      (dotimes (i0 numrows norm)
        (dotimes (i1 numcols)
          (setf norm (max norm (abs (aref data i0 i1)))))))))

(defmethod %norm ((data array) (measure (eql :frobenius)))
  "Return the Frobenius norm of the array."
  (multiple-value-bind (scale sumsq) (sumsq data)
    (* scale (sqrt sumsq))))

(defmethod %norm ((data array) (measure (eql :infinity)))
  "Return the infinity norm of the array."
  (destructuring-bind (numrows numcols) (array-dimensions data)
    (let ((norm 0)
          (sum 0))
      (dotimes (i0 numrows norm)
        (setf sum 0)
        (dotimes (i1 numcols)
          (incf sum (abs (aref data i0 i1))))
        (setf norm (max sum norm))))))

(defmethod norm ((data array) &key (measure 1))
  "Return the norm of the array."
  (%norm data measure))

(defmethod transpose ((data array) &key conjugate)
  "Return the transpose of the array."
  (destructuring-bind (numrows numcols)
      (array-dimensions data)
    (let ((op (if conjugate #'conjugate #'identity))
          (result
           (make-array
            (list numcols numrows)
            :element-type (array-element-type data))))
      (dotimes (i0 numrows result)
        (dotimes (i1 numcols)
          (setf
           (aref result i1 i0)
           (funcall op (aref data i0 i1))))))))

(defmethod ntranspose ((data array) &key conjugate)
  "Replace the contents of the array with the transpose."
  (let ((m-rows (array-dimension data 0))
        (n-columns (array-dimension data 1)))
    (if (= m-rows n-columns)
        (let ((op (if conjugate #'conjugate #'identity)))
          (dotimes (i0 m-rows data)
            ;; FIXME : Conjugate on the diagonal may not be correct.
            (setf (aref data i0 i0) (funcall op (aref data i0 i0)))
            (do ((i1 (1+ i0) (1+ i1)))
                ((>= i1 n-columns))
              (psetf
               (aref data i0 i1) (funcall op (aref data i1 i0))
               (aref data i1 i0) (funcall op (aref data i0 i1))))))
        (error "Rows(~D) and columns(~D) unequal."
               m-rows n-columns))))

(defmethod permute ((data array) (matrix permutation-matrix))
  (if (every #'= (array-dimensions data) (matrix-dimensions matrix))
      (right-permute-array data (contents matrix))
      (error "Array~A and permutation matrix~A sizes incompatible."
             (array-dimensions data) (matrix-dimensions matrix))))

(defmethod permute ((matrix permutation-matrix) (data array))
  (if (every #'= (array-dimensions data) (matrix-dimensions matrix))
      (left-permute-array (contents matrix) data)
      (error "Permutation matrix~A and array~A sizes incompatible."
             (matrix-dimensions matrix) (array-dimensions data))))

(defmethod npermute ((data array) (matrix permutation-matrix))
  "Destructively permute the array."
  (if (every #'= (array-dimensions data) (matrix-dimensions matrix))
      (right-npermute-array data (contents matrix))
      (error "Array~A and permutation matrix~A sizes incompatible."
             (array-dimensions data) (matrix-dimensions matrix))))

(defmethod npermute ((matrix permutation-matrix) (data array))
  "Destructively permute the array."
  (if (every #'= (array-dimensions data) (matrix-dimensions matrix))
      (left-npermute-array (contents (ntranspose matrix)) data)
      (error "Permutation matrix~A and array~A sizes incompatible."
             (matrix-dimensions matrix) (array-dimensions data))))

(defmethod scale ((scalar number) (data array))
  "Scale each element of the array."
  (destructuring-bind (numrows numcols) (array-dimensions data)
    (let ((result
           (make-array
            (list numrows numcols)
            :element-type (array-element-type data))))
      (dotimes (i0 numrows result)
        (dotimes (i1 numcols)
          (setf
           (aref result i0 i1)
           (* scalar (aref data i0 i1))))))))

(defmethod nscale ((scalar number) (data array))
  "Scale each element of the array."
  (destructuring-bind (numrows numcols) (array-dimensions data)
    (dotimes (i0 numrows data)
      (dotimes (i1 numcols)
        (setf
         (aref data i0 i1)
         (* scalar (aref data i0 i1)))))))

(defmethod compatible-dimensions-p
           ((operation (eql :add)) (array1 array) (array2 array))
  "Return true if the array dimensions are compatible for an
addition."
  (and
   (= 2 (array-rank array1) (array-rank array2))
   (= (array-dimension array1 0) (array-dimension array2 0))
   (= (array-dimension array1 1) (array-dimension array2 1))))

(defmethod add ((array1 array) (array2 array) &key scalar1 scalar2)
  "Return the addition of the 2 arrays."
  (if (compatible-dimensions-p :add array1 array2)
      (add-array array1 array2 scalar1 scalar2)
      (error "The array dimensions, ~A,~A, are not compatible."
             (array-dimensions array1) (array-dimensions array2))))

(defmethod nadd ((array1 array) (array2 array) &key scalar1 scalar2)
  "Destructively add array2 to array1."
  (if (compatible-dimensions-p :add array1 array2)
      (nadd-array array1 array2 scalar1 scalar2)
      (error "The array dimensions, ~A,~A, are not compatible."
             (array-dimensions array1) (array-dimensions array2))))

(defmethod subtract ((array1 array) (array2 array) &key scalar1 scalar2)
  "Return the subtraction of the 2 arrays."
  (if (compatible-dimensions-p :add array1 array2)
      (subtract-array array1 array2 scalar1 scalar2)
      (error "The array dimensions, ~A,~A, are not compatible."
             (array-dimensions array1) (array-dimensions array2))))

(defmethod nsubtract ((array1 array) (array2 array) &key scalar1 scalar2)
  "Destructively subtract array2 from array1."
  (if (compatible-dimensions-p :add array1 array2)
      (nsubtract-array array1 array2 scalar1 scalar2)
      (error "The array dimensions, ~A and ~A, are not compatible."
             (array-dimensions array1) (array-dimensions array2))))

(defmethod compatible-dimensions-p
           ((operation (eql :product)) (vector vector) (array array))
  "Return true if the array dimensions are compatible for product."
  (and
   (= 2 (array-rank array))
   (= (length vector) (array-dimension array 0))))

(defmethod product ((vector vector) (array array) &key scalar)
  "Return a vector generated by the pre-multiplication of a array by a
vector."
  (if (compatible-dimensions-p :product vector array)
      (right-product-vector vector array scalar)
      (error "Vector(~D) is incompatible with array~A."
             (length vector) (array-dimensions array))))

(defmethod compatible-dimensions-p
           ((operation (eql :product)) (array array) (vector vector))
  "Return true if the array dimensions are compatible for product."
  (and
   (= 2 (array-rank array))
   (= (array-dimension array 1) (length vector))))

(defmethod product ((array array) (vector vector) &key scalar)
  "Return a vector generated by the multiplication of the array with a
vector."
  (if (compatible-dimensions-p :product array vector)
      (let* ((m-rows (array-dimension array 0))
             (n-columns (array-dimension array 1))
             (zero (coerce 0 (array-element-type vector)))
             (element)
             (vector-product
              (make-array
               m-rows
               :element-type (array-element-type vector))))
        (dotimes (i0 m-rows vector-product)
          (setf element zero)
          (dotimes (i1 n-columns)
            (incf element (* (aref array i0 i1) (aref vector i1))))
          (if scalar
              (setf (aref vector-product i0) (* scalar element))
              (setf (aref vector-product i0) element))))
      (error "Array~A is incompatible with vector(~D)."
             (array-dimensions array) (length vector))))

(defmethod compatible-dimensions-p
           ((operation (eql :product)) (array1 array) (array2 array))
  "Return true if the array dimensions are compatible for product."
  (and
   (= 2 (array-rank array1) (array-rank array2))
   (= (array-dimension array1 1) (array-dimension array2 0))))

(defmethod product ((array1 array) (array2 array) &key scalar)
  "Return the product of the arrays."
  (if (compatible-dimensions-p :product array1 array2)
      (let* ((l-columns (array-dimension array1 1))
             (m-rows (array-dimension array1 0))
             (n-columns (array-dimension array2 1))
             (zero (coerce 0 (array-element-type array1)))
             (element)
             (array-product
              (make-array
               (list m-rows n-columns)
               :element-type (array-element-type array1))))
        (dotimes (i0 m-rows array-product)
          (dotimes (i2 n-columns)
            (setf element zero)
            (dotimes (i1 l-columns)
              (incf element (* (aref array1 i0 i1) (aref array2 i1 i2))))
            (if scalar
                (setf (aref array-product i0 i2) (* scalar element))
                (setf (aref array-product i0 i2) element)))))
      (error "The array dimensions, ~A and ~A, are not compatible."
             (array-dimensions array1) (array-dimensions array2))))
