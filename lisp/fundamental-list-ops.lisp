#|

 Fundamental List Operations

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

(defmethod sumsq ((data list) &key (scale 0) (sumsq 1))
  "Return the scaling parameter and the sum of the squares of the
list."
  (let ((abs-val nil))
    (dolist (elm data (values scale sumsq))
      (when (< 0 (setf abs-val (abs elm)))
        (if (< scale abs-val)
            (setf
             sumsq (1+ (* sumsq (expt (/ scale abs-val) 2)))
             scale abs-val)
            (incf sumsq (expt (/ elm scale) 2)))))))

(defmethod sump ((data list) (p real) &key (scale 0) (sump 1))
  "Return the scaling parameter and the sum of the powers of p of the
data."
  (let ((abs-val nil))
    (dolist (elm data (values scale sump))
      (when (< 0 (setf abs-val (abs elm)))
        (if (< scale abs-val)
            (setf
             sump  (1+ (* sump (expt (/ scale abs-val) p)))
             scale abs-val)
            (incf sump (expt (/ elm scale) p)))))))

(defmethod %norm ((data list) (measure (eql 1)))
  "Return the Taxicab norm of the list."
  (loop for element in data sum (abs element)))

(defmethod %norm ((data list) (measure (eql 2)))
  "Return the Euclidean norm of the vector."
  (multiple-value-bind (scale sumsq)
      (sumsq (loop for val in data collect (abs val)))
    (* scale (sqrt sumsq))))

(defmethod %norm ((data list) (measure integer))
  "Return the p-norm of the vector."
  (multiple-value-bind (scale sump)
      (sump (loop for val in data collect (abs val)) measure)
    (* scale (expt sump (/ measure)))))

(defmethod %norm ((data list) (measure (eql :infinity)))
  "Return the infinity, or maximum, norm of vector."
  (loop for element in data maximize (abs element)))

(defmethod norm ((data list) &key (measure 1))
  (%norm data measure))

(defmethod transpose ((data list) &key conjugate)
  "Return a row vector."
  (loop with op = (if conjugate #'conjugate #'identity)
        for val in data collect (funcall op val)))

(defmethod ntranspose ((data list) &key conjugate)
  "Return a row vector destructively."
  (if conjugate
      (map-into data #'conjugate data)
      data))

(defmethod permute :before
  ((data list) (matrix permutation-matrix))
  "Verify that the dimensions are compatible."
  (unless (= (length data) (matrix-row-dimension matrix))
    (error "List and permutation matrix sizes incompatible.")))

(defmethod permute ((data list) (matrix permutation-matrix))
  "Return the permutation of the list."
  (loop with permuted = (make-list (length data))
        for column across (contents matrix)
        and row = 0 then (1+ row)
        do (setf (nth column permuted) (nth row data))
        finally (return permuted)))

(defmethod permute :before
  ((matrix permutation-matrix) (data list))
  "Verify that the dimensions are compatible."
  (unless (= (length data) (matrix-column-dimension matrix))
    (error "List and permutation matrix sizes incompatible.")))

(defmethod permute ((matrix permutation-matrix) (data list))
  "Return the permutation of the list."
  (loop with permuted = (make-list (length data))
        for column across (contents matrix)
        and row = 0 then (1+ row)
        do (setf (nth row permuted) (nth column data))
        finally (return permuted)))

(defmethod npermute :before
  ((data list) (matrix permutation-matrix))
  "Verify that the dimensions are compatible."
  (unless (= (length data) (matrix-row-dimension matrix))
    (error "List and permutation matrix sizes incompatible.")))

(defmethod npermute ((data list) (matrix permutation-matrix))
  "Destructively permute the list."
  (multiple-value-bind (row0 skip)
      (%init-ntranspose (contents matrix))
    (loop with mat = (contents matrix)
          repeat (- (length mat) skip)
          for column = (aref mat row0) then (aref mat column)
          do (rotatef (nth column data) (nth row0 data))
          finally (return data))))

(defmethod npermute :before
  ((matrix permutation-matrix) (data list))
  "Verify that the dimensions are compatible."
  (unless (= (length data) (matrix-column-dimension matrix))
    (error "Vector and permutation matrix sizes incompatible.")))

(defmethod npermute ((matrix permutation-matrix) (data list))
  "Destructively permute the list."
  (multiple-value-bind (row0 skip)
      (%init-ntranspose (contents matrix))
    (loop with mat = (contents matrix)
          repeat (- (length mat) skip 1)
          for row = row0 then column
          as column = (aref mat row)
          do (rotatef (nth row data) (nth column data))
          finally (return data))))

(defmethod scale ((scalar number) (data list))
  "Return the list scaled by scalar."
  (loop for item in data collect (* scalar item)))

(defmethod nscale ((scalar number) (data list))
  "Return the list destructively scaled by scalar."
  (map-into data (lambda (x) (* scalar x)) data))

(defmethod add :before ((list1 list) (list2 list)
                        &key scalar1 scalar2)
  "Verify that the dimensions are equal."
  (declare (ignore scalar1 scalar2))
  (unless (= (length list1) (length list2))
    (error "LIST1 and LIST2 are not of equal length.")))

(defmethod add ((list1 list) (list2 list) &key scalar1 scalar2)
  "Return the addition of scalar1*list1 with scalar2*list2"
  (loop with op = (scaled-binary-op #'+ scalar1 scalar2)
        for item1 in list1
        and item2 in list2
        collect (funcall op item1 item2)))

(defmethod nadd :before ((list1 list) (list2 list)
                         &key scalar1 scalar2)
  "Verify that the dimensions are equal."
  (declare (ignore scalar1 scalar2))
  (unless (= (length list1) (length list2))
    (error
     "LENGTH1 and LENGTH2 are not of equal length in NADD-SCALED.")))

(defmethod nadd ((list1 list) (list2 list) &key scalar1 scalar2)
  "Return the addition of scalar2*list2 to scalar1*list1."
  (map-into list1 (scaled-binary-op #'+ scalar1 scalar2) list1 list2))
