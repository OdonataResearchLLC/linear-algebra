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

(defgeneric initialize-vector (vector data size element-type)
  (:documentation
   "Initialize the vector with data."))

(defun make-vector (size &key
                    (element-type 'number)
                    (vector-type 'column-vector)
                    (initial-element nil initial-element-p)
                    (initial-contents nil initial-contents-p))
  "Create a 1D numeric array to represent a numeric vector."
  (let ((new-vector (make-instance vector-type)))
    (cond
      ((and initial-element-p initial-contents-p)
       (error "Cannot specify both :INITIAL-ELEMENT and :INITIAL-CONTENTS."))
      (initial-contents-p
       (initialize-vector new-vector initial-contents size element-type))
      (initial-element-p
       (initialize-vector new-vector initial-element size element-type))
      (t
       (initialize-vector
        new-vector (coerce 0 element-type) size element-type)))))

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
