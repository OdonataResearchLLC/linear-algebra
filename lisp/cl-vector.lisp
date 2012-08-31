#|

 Fundamental Vector Operations

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
        (abs-val nil))
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
        (abs-val nil))
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
   (make-array (length data)
               :element-type
               (upgraded-array-element-type (type-of data)))
   (if conjugate #'conjugate #'identity)
   data))

(defmethod ntranspose ((data vector) &key conjugate)
  "Return a row vector destructively."
  (if conjugate
      (map-into data #'conjugate data)
      data))

(defmethod permute :before
  ((data vector) (matrix permutation-matrix))
  "Verify that the dimensions are compatible."
  (unless (= (length data) (matrix-row-dimension matrix))
    (error "Vector and permutation matrix sizes incompatible.")))

(defmethod permute ((data vector) (matrix permutation-matrix))
  "Return the permutation of the list."
  (loop with permuted =
        (make-array (length data) :element-type
                    (upgraded-array-element-type
                     (type-of data)))
        for column across (contents matrix)
        and row = 0 then (1+ row)
        do (setf (aref permuted column) (aref data row))
        finally (return permuted)))

(defmethod permute :before
  ((matrix permutation-matrix) (data vector))
  "Verify that the dimensions are compatible."
  (unless (= (length data) (matrix-column-dimension matrix))
    (error "Vector and permutation matrix sizes incompatible.")))

(defmethod permute ((matrix permutation-matrix) (data vector))
  "Return the permutation of the list."
  (loop with permuted =
        (make-array (length data) :element-type
                    (upgraded-array-element-type
                     (type-of data)))
        for column across (contents matrix)
        and row = 0 then (1+ row)
        do (setf (aref permuted row) (aref data column))
        finally (return permuted)))

(defmethod npermute :before
  ((data vector) (matrix permutation-matrix))
  "Verify that the dimensions are compatible."
  (unless (= (length data) (matrix-row-dimension matrix))
    (error "Vector and permutation matrix sizes incompatible.")))

(defmethod npermute ((data vector) (matrix permutation-matrix))
  "Destructively permute the vector."
  (multiple-value-bind (row0 skip)
      (%init-ntranspose (contents matrix))
    (loop with mat = (contents matrix)
          repeat (- (length mat) skip)
          for column = (aref mat row0) then (aref mat column)
          do (rotatef (aref data column) (aref data row0))
          finally (return data))))

(defmethod npermute :before
  ((matrix permutation-matrix) (data vector))
  "Verify that the dimensions are compatible."
  (unless (= (length data) (matrix-column-dimension matrix))
    (error "Vector and permutation matrix sizes incompatible.")))

(defmethod npermute ((matrix permutation-matrix) (data vector))
  "Destructively permute the list."
  (multiple-value-bind (row0 skip)
      (%init-ntranspose (contents matrix))
    (loop with mat = (contents matrix)
          repeat (- (length mat) skip 1)
          for row = row0 then column
          as column = (aref mat row)
          do (rotatef (aref data row) (aref data column))
          finally (return data))))

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

(defmethod add :before ((vector1 vector) (vector2 vector)
                        &key scalar1 scalar2)
  "Verify that the dimensions are equal."
  (declare (ignore scalar1 scalar2))
  (unless (= (length vector1) (length vector2))
    (error "VECTOR1 and VECTOR2 are not of equal length.")))

(defmethod add ((vector1 vector) (vector2 vector)
                &key scalar1 scalar2)
  "Return the addition of scalar1*vector1 with scalar2*vector2"
  (map-into
   (make-array (length vector1) :element-type
               (common-array-element-type vector1 vector2))
   (scaled-binary-op #'+ scalar1 scalar2)
   vector1 vector2))

(defmethod nadd :before ((vector1 vector) (vector2 vector)
                         &key scalar1 scalar2)
  "Verify that the dimensions are equal."
  (declare (ignore scalar1 scalar2))
  (unless (= (length vector1) (length vector2))
    (error
     "VECTOR1 and VECTOR2 are not of equal length in NADD-SCALED.")))

(defmethod nadd ((vector1 vector) (vector2 vector)
                 &key scalar1 scalar2)
  "Return the addition of scalar2*vector2 to scalar1*vector1."
  (map-into
   vector1
   (scaled-binary-op #'+ scalar1 scalar2)
   vector1 vector2))
