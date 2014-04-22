#|

 Fundamental Common Lisp List Operations

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
  (let ((abs-val))
    (dolist (elm data (values scale sumsq))
      (when (plusp (setq abs-val (abs elm)))
        (if (< scale abs-val)
            (setq
             sumsq (1+ (* sumsq (expt (/ scale abs-val) 2)))
             scale abs-val)
            (setq sumsq (+ sumsq (expt (/ elm scale) 2))))))))

(defmethod sump ((data list) (p real) &key (scale 0) (sump 1))
  "Return the scaling parameter and the sum of the powers of p of the
data."
  (let ((abs-val nil))
    (dolist (elm data (values scale sump))
      (when (plusp (setq abs-val (abs elm)))
        (if (< scale abs-val)
            (setq
             sump (1+ (* sump (expt (/ scale abs-val) p)))
             scale abs-val)
            (setq sump (+ sump (expt (/ elm scale) p))))))))

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
  (loop
   with op = (if conjugate #'conjugate #'identity)
   for val in data collect (funcall op val)))

(defmethod ntranspose ((data list) &key conjugate)
  "Return a row vector destructively."
  (if conjugate
      (map-into data #'conjugate data)
      data))

(defmethod permute ((data list) (matrix permutation-matrix))
  "Return the permutation of the list."
  (if (= (length data) (matrix-row-dimension matrix))
      (loop
       with permuted = (make-list (length data))
       for column across (contents matrix)
       and row = 0 then (1+ row)
       do (setf (nth column permuted) (nth row data))
       finally (return permuted))
      (error
       "List(~D) and permutation~A matrix sizes are incompatible."
       (length data) (matrix-dimensions matrix))))

(defmethod permute ((matrix permutation-matrix) (data list))
  "Return the permutation of the list."
  (if (= (length data) (matrix-row-dimension matrix))
      (loop
       with permuted = (make-list (length data))
       for column across (contents matrix)
       and row = 0 then (1+ row)
       do (setf (nth row permuted) (nth column data))
       finally (return permuted))
      (error
       "Permutation matrix~A and list(~D) sizes are incompatible."
       (matrix-dimensions matrix) (length data))))

(defmethod scale ((scalar number) (data list))
  "Return the list scaled by scalar."
  (loop for item in data collect (* scalar item)))

(defmethod nscale ((scalar number) (data list))
  "Return the list destructively scaled by scalar."
  (map-into data (lambda (x) (* scalar x)) data))

(defmethod add ((list1 list) (list2 list) &key scalar1 scalar2)
  "Return the addition of scalar1*list1 with scalar2*list2"
  (if (= (length list1) (length list2))
      (loop
       with op = (scaled-binary-op #'+ scalar1 scalar2)
       for item1 in list1
       and item2 in list2
       collect (funcall op item1 item2))
      (error "LIST1(~D) and LIST2(~D) are not of equal length."
             (length list1) (length list2))))

(defmethod nadd ((list1 list) (list2 list) &key scalar1 scalar2)
  "Return the addition of scalar2*list2 to scalar1*list1."
  (if (= (length list1) (length list2))
      (map-into
       list1 (scaled-binary-op #'+ scalar1 scalar2)
       list1 list2)
      (error "LIST1(~D) and LIST2(~D) are not of equal length."
             (length list1) (length list2))))

(defmethod subtract ((list1 list) (list2 list) &key scalar1 scalar2)
  "Return the subraction of scalar2*list2 from scalar1*list1."
  (if (= (length list1) (length list2))
      (loop
       with op = (scaled-binary-op #'- scalar1 scalar2)
       for item1 in list1
       and item2 in list2
       collect (funcall op item1 item2))
      (error "LIST1(~D) and LIST2(~D) are not of equal length."
             (length list1) (length list2))))

(defmethod nsubtract ((list1 list) (list2 list) &key scalar1 scalar2)
  "Return the subraction of scalar2*list2 from scalar1*list1."
  (if (= (length list1) (length list2))
      (map-into
       list1 (scaled-binary-op #'- scalar1 scalar2)
       list1 list2)
      (error "LIST1(~D) and LIST2(~D) are not of equal length."
             (length list1) (length list2))))

(defmethod product ((list1 list) (list2 list)
                    &key (scalar nil scalarp) conjugate)
  "Return the dot product of list1 and list2."
  (if (= (length list1) (length list2))
      (loop
       with op =
       (if conjugate
           (lambda (x y) (* (conjugate x) y))
           #'*)
       for element1 in list1
       and element2 in list2
       sum (funcall op element1 element2) into result
       finally
       (return (if scalarp (* scalar result) result)))
      (error "LIST1(~D) and LIST2(~D) are not of equal length."
             (length list1) (length list2))))
