#|

  Linear Algebra Permutation Kernel

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

;;; Right permutation

(defgeneric right-permute (vector-or-array permutation)
  (:documentation
   "Permute the row vector or columns of the array."))

(defmethod right-permute ((data vector) (permutation vector))
  "Permute the row vector to create a column vector."
  (loop
   with result =
   (make-array (length data) :element-type (array-element-type data))
   for column across permutation
   and row = 0 then (1+ row)
   do (setf (aref result column) (aref data row))
   finally return result))

(defmethod right-permute ((data array) (permutation vector))
  "Permute the columns of the array."
  (loop
   with m-rows = (array-dimension data 0)
   with result =
   (make-array
    (array-dimensions data) :element-type (array-element-type data))
   for column across permutation
   and row = 0 then (1+ row)
   do
   (loop
    for irow below m-rows do
    (setf (aref result irow column) (aref data irow row)))
   finally return result))

;;; Left permutation

(defgeneric left-permute (permutation vector-or-array)
  (:documentation
   "Permute the column vector or rows of the array."))

(defmethod left-permute ((permutation vector) (data vector))
  "Permute the column vector to create a row vector."
  (loop
   with result =
   (make-array (length data) :element-type (array-element-type data))
   for column across permutation
   and row = 0 then (1+ row)
   do (setf (aref result row) (aref data column))
   finally return result))

(defmethod left-permute ((permutation vector) (data array))
  "Permute the rows of the array."
  (loop
   with n-columns = (array-dimension data 1)
   with result =
   (make-array
    (array-dimensions data) :element-type (array-element-type data))
   for column across permutation
   and row = 0 then (1+ row)
   do
   (loop
    for icol below n-columns do
    (setf (aref result row icol) (aref data column icol)))
   finally return result))
