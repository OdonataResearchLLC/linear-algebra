#|

  Linear Algebra in Common Lisp

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

;;; Vector interface operations

(defun make-vector
       (size &key
        (vector-type 'column-vector)
        (element-type 'number)
        initial-element initial-contents)
  "Create the data structure to represent a vector."
  (make-instance
   vector-type :size size
   :element-type element-type
   :initial-element initial-element
   :initial-contents initial-contents))

(defgeneric vector-in-bounds-p (vector index)
  (:documentation
   "Return true if index does not exceed the dimensions of vector."))

(defgeneric vector-element-type (vector)
  (:documentation
   "Return the element type of vector."))

(defgeneric vector-length (vector)
  (:documentation
   "Return the length of the vector."))

(defgeneric vref (vector index)
  (:documentation
   "Return the element of vector at index."))

(defgeneric (setf vref) (data vector index)
  (:documentation
   "Set the element of vector at index to data."))

(defgeneric copy-vector (vector)
  (:documentation
   "Return a copy of the vector."))

(defgeneric subvector (vector start &optional end)
  (:documentation
   "Return a new vector that is a subvector of the vector."))

(defgeneric (setf subvector) (subvector vector start &optional end)
  (:documentation
   "Set the subvector of the vector."))

(defgeneric replace-vector
    (vector1 vector2 &key start1 end1 start2 end2)
  (:documentation
   "Destructively replace the elements of vector1 with vector2."))

;;; Vector iteration operations

(defgeneric map-vector
    (result-type function first-vector &rest more-vectors)
  (:documentation
   "Calls function on successive sets of vector objects."))

(defgeneric map-into-vector (result-vector function &rest vectors)
  (:documentation
   "Destructively modifies the result vector with the result of
applying the function to each element of the vectors."))

(defgeneric reduce-vector (function vector initial-value)
  (:documentation
   "Reduce a vector into a scalar using a function"))

(defmacro dovector ((element vector &optional result) &body body)
  "Iterate over vector returning result."
  (let ((pos (gensym "POS-"))
        (end (gensym "END-")))
    `(let ((,end (vector-length ,vector))
           (,element nil))
      (dotimes (,pos ,end ,result)
        (setf ,element (vref ,vector ,pos))
        ,@body))))

;;; Vector transformations

(defgeneric apply-rotation (vector1 vector2 cc ss)
  (:documentation
   "Return the plane rotations of vector1 and vector2 by cc and ss."))

(defgeneric napply-rotation (vector1 vector2 cc ss)
  (:documentation
   "Return the plane rotations of vector1 and vector2 by cc and ss."))

(defgeneric vec-equal (vector1 vector2)
  (:documentation
   "Test if the contents of vector1 and vector2 are the same"))

(defgeneric elem-divide (vector1 vector2)
  (:documentation
   "Return the element by element division of two vectors"))

(defgeneric elem-multiply (vector1 vector2)
  (:documentation
   "Return the element by element multiplation of two vectors"))

(defgeneric elem-greater (vector1 vector2)
  (:documentation
   "Return if vector1 is greater than vector 2 on every element"))

(defgeneric vec-every (vector predicate)
  (:documentation
   "Return true if every element in the vector satisfies
  predicate. Return nil on the first element that does not"))

(defgeneric outer-product (vec1 vec2)
  (:documentation
   "Outer product of vec1 and vec2 returns an nxm matrix"))
