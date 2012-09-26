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

(in-package :linear-algebra)

;;; Data vector classes

(defclass data-vector ()
  ((contents
    :type (array * (*))
    :initarg :contents
    :accessor contents))
  (:documentation
   "A data vector."))

(defclass row-vector (data-vector)
  ()
  (:documentation
   "A row vector."))

(defclass column-vector (data-vector)
  ()
  (:documentation
   "A column vector."))

;;; Data vector interface operations

(defmethod initialize-vector ((vector data-vector)
                              (data number)
                              (size integer)
                              element-type)
  "Initialize a data vector with a value."
  (setf (contents vector)
        (make-array size
                    :element-type element-type
                    :initial-element data))
  ;; Return the data vector
  vector)

(defmethod initialize-vector ((vector data-vector)
                              (data sequence)
                              (size integer)
                              element-type)
  "Initialize a data vector with a sequence."
  (setf (contents vector)
        (make-array size
                    :element-type element-type
                    :initial-contents data))
  ;; Return the data vector
  vector)

(defun row-vector (&rest numbers)
  "Create a row vector from the numbers."
  (make-vector (length numbers)
               :vector-type 'row-vector
               :initial-contents numbers))

(defun column-vector (&rest numbers)
  "Create a column vector from the numbers."
  (make-vector (length numbers)
               :vector-type 'column-vector
               :initial-contents numbers))

(defun row-vector-p (object)
  "Return true if object is a row-vector, NIL otherwise."
  (typep object 'row-vector))

(defun column-vector-p (object)
  "Return true if object is a column-vector, NIL otherwise."
  (typep object 'column-vector))

(defmethod vector-in-bounds-p ((vector data-vector) (index integer))
  "Return true if index does not exceed the dimensions of vector."
  (array-in-bounds-p (contents vector) index))

(defmethod vector-element-type ((vector data-vector))
  "Return the element type of vector."
  (array-element-type (contents vector)))

(defmethod vector-length ((vector data-vector))
  "Return the length of the vector."
  (length (contents vector)))

(defmethod vref ((vector data-vector) (index integer))
  "Return the element of vector at index."
  (aref (contents vector) index))

(defmethod (setf vref) ((data number) (vector data-vector) (index integer))
  "Set the element of vector at index to data."
  (setf (aref (contents vector) index) data))

(defmethod copy-vector ((vector data-vector))
  "Return a copy of the vector."
  (make-instance
   (class-of vector)
   :contents (copy-seq (contents vector))))

(defmethod subvector ((vector data-vector) start &optional end)
  "Return a new data vector that is a subset of vector."
  (make-instance
   (class-of vector)
   :contents (subseq (contents vector) start end)))

(defmethod (setf subvector) ((subvector data-vector)
                             (vector data-vector)
                             start &optional end)
  "Set the subvector of the data vector."
  (setf (subseq (contents vector) start end)
        (contents subvector))
  ;; Return the subvector
  subvector)

(defmethod replace-vector ((vector1 data-vector) (vector2 data-vector)
                           &key (start1 0) end1 (start2 0) end2)
  "Destructively replace the elements of vector1 with vector2."
  (replace (contents vector1) (contents vector2)
           :start1 start1 :end1 end1 :start2 start2 :end2 end2)
  ;; Return vector1
  vector1)

;;; Data vector iteration operations
(defun %map-data-vector (result-type function first-vector
                         &rest more-vectors)
  "Non-validating version of map-vector."
  (make-instance
   result-type
   :contents
   (apply #'map
          (class-of (contents first-vector))
          function
          (contents first-vector)
          (mapcar #'contents more-vectors))))

(defmethod map-vector :before (result-type (function function)
                               (first-vector data-vector)
                               &rest more-vectors)
  "Verify the arguments to map-vector."
  (declare (ignore function))
  (unless (subtypep result-type 'data-vector)
    (error "~A is not a subtype of DATA-VECTOR." result-type))
  (unless (every (lambda (x) (typep x 'data-vector)) more-vectors)
    (error "All vectors must be data vectors.")))

(defmethod map-vector (result-type (function function)
                       (first-vector data-vector)
                       &rest more-vectors)
  "Calls function on successive sets of data vectors."
  (apply #'%map-data-vector
         result-type
         function
         first-vector
         more-vectors))

(defun %map-into-data-vector (result-vector function &rest vectors)
  "Non-validating version of map-into-vector."
  (apply #'map-into
         (contents result-vector)
         function
         (mapcar #'contents vectors))
  ;; Return the result vector
  result-vector)

(defmethod map-into-vector :before ((result-vector data-vector)
                                    (function function) &rest vectors)
  "Verify the arguments to map-into-vector."
  (declare (ignore result-vector function))
  (unless (every (lambda (x) (typep x 'data-vector)) vectors)
    (error "All vectors must be data vectors.")))

(defmethod map-into-vector ((result-vector data-vector)
                            (function function) &rest vectors)
  "Destructively modifies the result vector with the result of
applying the function to each element of the vectors."
  (apply #'%map-into-data-vector
         result-vector
         function
         vectors))

;;; Data vector transformations

(defmethod apply-rotation :before ((vector1 data-vector) (vector2 data-vector) cc ss)
  "Verify the input to apply-rotation."
  (declare (ignore cc ss))
  (unless (= (vector-length vector1) (vector-length vector2))
    (error "VECTOR1 and VECTOR2 are not of equal length.")))

(defmethod apply-rotation ((vector1 data-vector) (vector2 data-vector) cc ss)
  "Return the plane rotations of vector1 and vector2 by cc and ss."
  (let ((rvec1 (make-vector (vector-length vector1)
                            :vector-type (class-of vector1)
                            :element-type (vector-element-type vector1)))
        (rvec2 (make-vector (vector-length vector2)
                            :vector-type (class-of vector2)
                            :element-type (vector-element-type vector2))))
    (dotimes (pos (vector-length vector1) (values rvec1 rvec2))
      (setf (vref rvec1 pos)
            (+ (* cc (vref vector1 pos))
               (* ss (vref vector2 pos))))
      (setf (vref rvec2 pos)
            (+ (* -1 (conjugate ss) (vref vector1 pos))
               (* cc (vref vector2 pos)))))))

(defmethod napply-rotation :before ((vector1 data-vector) (vector2 data-vector) cc ss)
  "Verify the input to napply-rotation."
  (declare (ignore cc ss))
  (unless (= (vector-length vector1) (vector-length vector2))
    (error "VECTOR1 and VECTOR2 are not of equal length.")))

(defmethod napply-rotation ((vector1 data-vector) (vector2 data-vector) cc ss)
  "Return the plane rotations of vector1 and vector2 by cc and ss."
  (dotimes (pos (vector-length vector1) (values vector1 vector2))
    (psetf (vref vector1 pos)
           (+ (* cc (vref vector1 pos))
              (* ss (vref vector2 pos)))
           (vref vector2 pos)
           (+ (* -1 (conjugate ss) (vref vector1 pos))
              (* cc (vref vector2 pos))))))

;;; Data vector fundamental operations

(defmethod sumsq ((vector data-vector) &key (scale 0) (sumsq 1))
  "Return the scaling parameter and the sum of the squares of vector."
  (sumsq-vector (contents vector) scale sumsq))

(defmethod sump ((vector data-vector) (p number) &key (scale 0) (sump 1))
  "Return the scaling parameter and the sum of the P powers of vector."
  (sump-vector (contents vector) p scale sump))

(defmethod norm ((vector data-vector) &key (measure 1))
  "Return the p-norm of the vector."
  (norm-vector (contents vector) measure))

(defmethod transpose ((vector column-vector) &key conjugate)
  "Return a row vector."
  (%map-data-vector 'row-vector
                    (if conjugate #'conjugate #'identity)
                    vector))

(defmethod transpose ((vector row-vector) &key conjugate)
  "Return a column vector."
  (%map-data-vector 'column-vector
                    (if conjugate #'conjugate #'identity)
                    vector))

(defmethod ntranspose ((vector column-vector) &key conjugate)
  "Return a row vector destructively."
  (if conjugate
      (%map-into-data-vector
       (change-class vector 'row-vector) #'conjugate vector)
      (change-class vector 'row-vector)))

(defmethod ntranspose ((vector row-vector) &key conjugate)
  "Return a column vector destructively."
  (if conjugate
      (%map-into-data-vector
       (change-class vector 'column-vector) #'conjugate vector)
      (change-class vector 'column-vector)))

(defmethod permute :before ((vector row-vector) (matrix permutation-matrix))
  "Verify that the dimensions are compatible."
  (unless (= (vector-length vector) (matrix-row-dimension matrix))
    (error "Vector and permutation matrix sizes incompatible.")))

(defmethod permute ((vector row-vector) (matrix permutation-matrix))
  "Return the permutation of the row vector."
  (make-instance
   'row-vector
   :contents
   (loop with permuted =
         (make-array (vector-length vector)
                     :element-type (vector-element-type vector))
         for column across (contents matrix)
         and row = 0 then (1+ row)
         do (setf (aref permuted column) (vref vector row))
         finally (return permuted))))

(defmethod permute :before ((matrix permutation-matrix) (vector column-vector))
  "Verify that the dimensions are compatible."
  (unless (= (vector-length vector) (matrix-column-dimension matrix))
    (error "Vector and permutation matrix sizes incompatible.")))

(defmethod permute ((matrix permutation-matrix) (vector column-vector))
  "Return the permutation of the column vector."
  (make-instance
   'row-vector
   :contents
   (loop with permuted =
         (make-array (vector-length vector)
                     :element-type (vector-element-type vector))
         for column across (contents matrix)
         and row = 0 then (1+ row)
         do (setf (aref permuted row) (vref vector column))
         finally (return permuted))))

(defmethod npermute :before ((vector row-vector) (matrix permutation-matrix))
  "Verify that the dimensions are compatible."
  (unless (= (vector-length vector) (matrix-row-dimension matrix))
    (error "Vector and permutation matrix sizes incompatible.")))

(defmethod npermute ((vector row-vector) (matrix permutation-matrix))
  "Destructively permute the row vector."
  (multiple-value-bind (row0 skip)
      (%init-ntranspose (contents matrix))
    (loop with mat = (contents matrix)
          with vec = (contents vector)
          repeat (- (length mat) skip)
          for column-swap = (aref mat row0) then column-next
          and value-swap  = (aref vec row0) then value-next
          as column-next = (aref mat column-swap)
          as value-next  = (aref vec column-swap)
          do (setf (aref vec column-swap) value-swap)
          finally (return vector))))

(defmethod npermute :before ((matrix permutation-matrix) (vector column-vector))
  "Verify that the dimensions are compatible."
  (unless (= (vector-length vector) (matrix-column-dimension matrix))
    (error "Vector and permutation matrix sizes incompatible.")))

(defmethod npermute ((matrix permutation-matrix) (vector column-vector))
  "Destructively permute the column vector."
  (multiple-value-bind (row0 skip)
      (%init-ntranspose (contents matrix))
    (loop with mat = (contents matrix)
          with vec = (contents vector)
          repeat (- (length mat) skip 1)
          for row = row0 then column
          as column = (aref mat row)
          do (rotatef (aref vec row) (aref vec column))
          finally (return vector))))

(defmethod scale ((scalar number) (vector data-vector))
  "Return the vector scaled by scalar."
  (%map-data-vector (class-of vector)
                    (lambda (x) (* scalar x))
                    vector))

(defmethod nscale ((scalar number) (vector data-vector))
  "Return the vector destructively scaled by scalar."
  (%map-into-data-vector vector
                         (lambda (x) (* scalar x))
                         vector))

(defmethod add :before ((vector1 data-vector) (vector2 data-vector)
                        &key scalar1 scalar2)
  "Verify that the dimensions are equal."
  (declare (ignore scalar1 scalar2))
  (unless (= (vector-length vector1) (vector-length vector2))
    (error "VECTOR1 and VECTOR2 are not of equal length.")))

(defmethod add ((vector1 column-vector) (vector2 column-vector)
                &key scalar1 scalar2)
  "Return the addition of scalar1*vector1 with scalar2*vector2."
  (%map-data-vector (common-class-of vector1 vector2 'column-vector)
                    (scaled-binary-op #'+ scalar1 scalar2)
                    vector1 vector2))

(defmethod add ((vector1 row-vector) (vector2 row-vector)
                &key scalar1 scalar2)
  "Return the addition of scalar1*vector1 with scalar2*vector2."
  (%map-data-vector (common-class-of vector1 vector2 'row-vector)
                    (scaled-binary-op #'+ scalar1 scalar2)
                    vector1 vector2))

(defmethod nadd :before ((vector1 data-vector) (vector2 data-vector)
                         &key scalar1 scalar2)
  "Verify that the dimensions are equal."
  (declare (ignore scalar1 scalar2))
  (unless (= (vector-length vector1) (vector-length vector2))
    (error "VECTOR1 and VECTOR2 are not of equal length in NADD-SCALED.")))

(defmethod nadd ((vector1 column-vector) (vector2 column-vector)
                 &key scalar1 scalar2)
  "Return the addition of scalar2*vector2 to scalar1*vector1."
  (%map-into-data-vector vector1
                         (scaled-binary-op #'+ scalar1 scalar2)
                         vector1 vector2))

(defmethod nadd ((vector1 row-vector) (vector2 row-vector)
                 &key scalar1 scalar2)
  "Return the addition of scalar2*vector2 to scalar1*vector1."
  (%map-into-data-vector vector1
                         (scaled-binary-op #'+ scalar1 scalar2)
                         vector1 vector2))

(defmethod subtract :before ((vector1 data-vector) (vector2 data-vector)
                             &key scalar1 scalar2)
  "Verify that the dimensions are equal."
  (declare (ignore scalar1 scalar2))
  (unless (= (vector-length vector1) (vector-length vector2))
    (error "VECTOR1 and VECTOR2 are not of equal length.")))

(defmethod subtract ((vector1 column-vector) (vector2 column-vector)
                     &key scalar1 scalar2)
  "Return the subraction of scalar2*vector2 from scalar1*vector1."
  (%map-data-vector (common-class-of vector1 vector2 'column-vector)
                    (scaled-binary-op #'- scalar1 scalar2)
                    vector1 vector2))

(defmethod subtract ((vector1 row-vector) (vector2 row-vector)
                     &key scalar1 scalar2)
  "Return the subraction of scalar2*vector2 from scalar1*vector1."
  (%map-data-vector (common-class-of vector1 vector2 'row-vector)
                    (scaled-binary-op #'- scalar1 scalar2)
                    vector1 vector2))

(defmethod nsubtract :before ((vector1 data-vector) (vector2 data-vector)
                              &key scalar1 scalar2)
  "Verify that the dimensions are equal."
  (declare (ignore scalar1 scalar2))
  (unless (= (vector-length vector1) (vector-length vector2))
    (error "VECTOR1 and VECTOR2 are not of equal length.")))

(defmethod nsubtract ((vector1 column-vector) (vector2 column-vector)
                      &key scalar1 scalar2)
  "Return the subraction of scalar2*vector2 from scalar1*vector1."
  (%map-into-data-vector vector1
                         (scaled-binary-op #'- scalar1 scalar2)
                         vector1 vector2))

(defmethod nsubtract ((vector1 row-vector) (vector2 row-vector)
                      &key scalar1 scalar2)
  "Return the subraction of scalar2*vector2 from scalar1*vector1."
  (%map-into-data-vector vector1
                         (scaled-binary-op #'- scalar1 scalar2)
                         vector1 vector2))

(defmethod product :before ((vector1 row-vector) (vector2 column-vector)
                            &key scalar conjugate)
  "Verify that the dimensions are equal."
  (declare (ignore scalar conjugate))
  (unless (= (vector-length vector1) (vector-length vector2))
    (error "VECTOR1 and VECTOR2 are not of equal length.")))

(defmethod product ((vector1 row-vector) (vector2 column-vector)
                    &key (scalar nil scalarp) conjugate)
  "Return the dot product of vector1 and vector2."
  (loop with op = (if conjugate
                      (lambda (x y) (* (conjugate x) y))
                      #'*)
        for element1 across (contents vector1)
        and element2 across (contents vector2)
        sum (funcall op element1 element2) into result
        finally
        (return (if scalarp (* scalar result) result))))
