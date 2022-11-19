#|

  Linear Algebra Unary Operations Kernel

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

(in-package :linear-algebra-kernel)

;;; Squared sums
;;; FIXME : Consider locating sumsq2 and sumsq3 in floating-point

(defun sumsq2 (x y)
  "Return the square root of |x|^2 + |y|^2."
  (let* ((abs-x (abs x))
         (abs-y (abs y))
         (w (max abs-x abs-y))
         (v/w (/ (min abs-x abs-y) w)))
    (* w (sqrt (+ 1 (* v/w v/w))))))

(defun sumsq3 (x y z)
  "Return the square root of |x|^2 + |y|^2 + |z|^2."
  (let* ((abs-x (abs x))
         (abs-y (abs y))
         (abs-z (abs z))
         (w (max abs-x abs-y abs-z))
         (x/w (/ abs-x w))
         (y/w (/ abs-y w))
         (z/w (/ abs-z w)))
    (* w (sqrt (+ (* x/w x/w) (* y/w y/w) (* z/w z/w))))))

(defgeneric sumsq (vector-or-array &optional scale sumsq)
  (:documentation
   "Return the scaling parameter and the sum of the squares."))

(defmethod sumsq ((data list) &optional (scale 1) (sumsq 0))
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

(defmethod sumsq ((data vector) &optional (scale 1) (sumsq 0))
  "Return the scaling parameter and the sum of the squares of the
vector."
  (let ((abs-val))
    (dotimes (index (length data) (values scale sumsq))
      (when (plusp (setq abs-val (abs (aref data index))))
        (if (< scale abs-val)
            (setq
             sumsq (1+ (* sumsq (expt (/ scale abs-val) 2)))
             scale abs-val)
            (setq
             sumsq
             (+ sumsq (expt (/ (aref data index) scale) 2))))))))

(defmethod sumsq ((data array) &optional (scale 1) (sumsq 0))
    "Return the scaling parameter and the sum of the squares of the
array."
  (let ((m-rows (array-dimension data 0))
        (n-columns (array-dimension data 1))
        (abs-val 0))
    (dotimes (row m-rows (values scale sumsq))
      (dotimes (column n-columns)
        (when (plusp (setq abs-val (abs (aref data row column))))
          (if (< scale abs-val)
              (setq
               sumsq (1+ (* sumsq (expt (/ scale abs-val) 2)))
               scale abs-val)
              (setq sumsq (+ sumsq (expt (/ abs-val scale) 2)))))))))

(defgeneric sump (vector-or-array p &optional scale sump)
  (:documentation
   "Return the scaling parameter and the sum of the P powers."))

(defmethod sump ((data list) (p real) &optional (scale 1) (sump 0))
  "Return the scaling parameter and the sum of the powers of p of the
data."
  (let ((abs-val))
    (dolist (elm data (values scale sump))
      (when (plusp (setq abs-val (abs elm)))
        (if (< scale abs-val)
            (setq
             sump (1+ (* sump (expt (/ scale abs-val) p)))
             scale abs-val)
            (setq sump (+ sump (expt (/ elm scale) p))))))))

(defmethod sump ((data vector) p &optional (scale 1) (sump 0))
  "Return the scaling parameter and the sum of the powers of p of the
vector."
  (let ((abs-val))
    (dotimes (index (length data) (values scale sump))
      (when (plusp (setq abs-val (abs (aref data index))))
        (if (< scale abs-val)
            (setq
             sump (1+ (* sump (expt (/ scale abs-val) p)))
             scale abs-val)
            (setq
             sump
             (+ sump (expt (/ (aref data index) scale) p))))))))

(defmethod sump ((data array) p &optional (scale 1) (sump 0))
  "Return the scaling parameter and the sum of the P powers of the
matrix."
  (unless (plusp p) (error "The power(~A) must be positive." p))
  (let ((m-rows (array-dimension data 0))
        (n-columns (array-dimension data 1))
        (abs-val 0))
    (dotimes (row m-rows (values scale sump))
      (dotimes (column n-columns)
        (when (plusp (setq abs-val (abs (aref data row column))))
          (if (< scale abs-val)
              (setq
               sump (1+ (* sump (expt (/ scale abs-val) p)))
               scale abs-val)
              (setq
               sump
               (+ sump
                  (expt (/ (aref data row column) scale) p)))))))))

(defun sumsq-row (array row &key (scale 1) (sumsq 0) start end)
  "Return the scaling parameter and the sum of the squares of the
array row."
  (loop
   with start = (or start 0)
   and end = (or end (array-dimension array 1))
   and abs-val = 0
   for column from start below end
   when (plusp (setq abs-val (abs (aref array row column)))) do
   (if (< scale abs-val)
       (setq
        sumsq (1+ (* sumsq (expt (/ scale abs-val) 2)))
        scale abs-val)
       (setq sumsq (+ sumsq (expt (/ abs-val scale) 2))))
   finally (return (values scale sumsq))))

(defun sumsq-column (array column &key (scale 1) (sumsq 0) start end)
  "Return the scaling parameter and the sum of the squares of the
array column."
  (loop
   with start = (or start 0)
   and end = (or end (array-dimension array 0))
   and abs-val = 0
   for row from start below end
   when (plusp (setq abs-val (abs (aref array row column)))) do
   (if (< scale abs-val)
       (setq
        sumsq (1+ (* sumsq (expt (/ scale abs-val) 2)))
        scale abs-val)
       (setq sumsq (+ sumsq (expt (/ abs-val scale) 2))))
   finally (return (values scale sumsq))))

;;; Norm

(defgeneric norm-vector (data measure)
  (:documentation
   "Return the norm of the vector according to the measure."))

(defgeneric norm-array (data measure)
  (:documentation
   "Return the norm of the array according to the measure."))

(defun %abs-vector (vector)
  "Return a vector containing absolute value of each element."
  (let ((result
         (make-array
          (length vector)
          :element-type (array-element-type vector))))
    (dotimes (index (length vector) result)
      (setf (aref result index) (abs (aref vector index))))))


(defmethod norm-vector ((data vector) (measure (eql 1)))
  "Return the Taxicab norm of the list."
  (loop for element across data sum (abs element)))

(defmethod norm-vector ((data vector) (measure (eql 2)))
  "Return the Euclidean norm of the vector."
  (multiple-value-bind (scale sumsq)
      (sumsq (%abs-vector data))
    (* scale (sqrt sumsq))))

(defmethod norm-vector ((data vector) (measure integer))
  "Return the p-norm of the vector."
  (multiple-value-bind (scale sump)
      (sump (%abs-vector data) measure)
    (* scale (expt sump (/ measure)))))

(defmethod norm-vector ((data vector) (measure (eql :infinity)))
  "Return the infinity, or maximum, norm of vector."
  (loop for element across data maximize (abs element)))


(defmethod norm-array ((data array) (measure (eql 1)))
  "Return the 1 norm of the array."
  (let ((m-rows (array-dimension data 0))
        (n-columns (array-dimension data 1))
        (norm 0)
        (sum 0))
    (dotimes (column n-columns norm)
      (setq sum 0)
      (dotimes (row m-rows)
        (setq sum (+ sum (abs (aref data row column)))))
      (setq norm (max sum norm)))))

(defmethod norm-array ((data array) (measure (eql :max)))
  "Return the max norm of the array."
  (let ((m-rows (array-dimension data 0))
        (n-columns (array-dimension data 1))
        (norm 0))
    (dotimes (row m-rows norm)
      (dotimes (column n-columns)
        (setq norm (max norm (abs (aref data row column))))))))

(defmethod norm-array ((data array) (measure (eql :frobenius)))
  "Return the Frobenius norm of the array."
  (multiple-value-bind (scale sumsq)
      (sumsq data)
    (* scale (sqrt sumsq))))

(defmethod norm-array ((data array) (measure (eql :infinity)))
  "Return the infinity norm of the array."
  (let ((m-rows (array-dimension data 0))
        (n-columns (array-dimension data 1))
        (norm 0)
        (sum 0))
    (dotimes (row m-rows norm)
      (setq sum 0)
      (dotimes (column n-columns)
        (setq sum (+ sum (abs (aref data row column)))))
      (setq norm (max sum norm)))))
