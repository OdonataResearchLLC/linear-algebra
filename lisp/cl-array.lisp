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
  (destructuring-bind (numrows numcols) (array-dimensions data)
    (if (= numrows numcols)
        (let ((op (if conjugate #'conjugate #'identity)))
          (dotimes (i0 numrows data)
            ;; FIXME : Conjugate on the diagonal may not be correct.
            (setf (aref data i0 i0) (funcall op (aref data i0 i0)))
            (do ((i1 (1+ i0) (1+ i1)))
                ((>= i1 numcols))
              (psetf
               (aref data i0 i1) (funcall op (aref data i1 i0))
               (aref data i1 i0) (funcall op (aref data i0 i1))))))
        (error "Rows and columns unequal."))))

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

(defmethod add :before
  ((array1 array) (array2 array) &key scalar1 scalar2)
  "Audit the input data."
  (declare (ignore scalar1 scalar2))
  (unless (every #'=
                 (array-dimensions array1)
                 (array-dimensions array2))
    (error "The array dimensions are not compatible.")))

(defmethod add ((array1 array) (array2 array) &key scalar1 scalar2)
  "Return the addition of the 2 matrices."
  (destructuring-bind (numrows numcols) (array-dimensions array1)
    (let ((op (scaled-binary-op #'+ scalar1 scalar2))
          (result
           (make-array
            (list numrows numcols)
            :element-type
            (common-array-element-type array1 array2))))
      (dotimes (i0 numrows result)
        (dotimes (i1 numcols)
          (setf
           (aref result i0 i1)
           (funcall op (aref array1 i0 i1) (aref array2 i0 i1))))))))

(defmethod nadd :before
  ((array1 array) (array2 array) &key scalar1 scalar2)
  "Audit the input data."
  (declare (ignore scalar1 scalar2))
  (unless (every #'=
                 (array-dimensions array1)
                 (array-dimensions array2))
    (error "The array dimensions are not compatible.")))

(defmethod nadd ((array1 array) (array2 array) &key scalar1 scalar2)
  "Destructively add array2 to array1."
  (destructuring-bind (numrows numcols) (array-dimensions array1)
    (let ((op (scaled-binary-op #'+ scalar1 scalar2)))
      (dotimes (i0 numrows array1)
        (dotimes (i1 numcols)
          (setf
           (aref array1 i0 i1)
           (funcall op (aref array1 i0 i1) (aref array2 i0 i1))))))))
