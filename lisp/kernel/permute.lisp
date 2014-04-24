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

;;; Permute Vectors

(defun right-permute-vector (vector permutation)
  "Permute the row vector to create a column vector."
  (loop
   with result =
   (make-array
    (length vector)
    :element-type (array-element-type vector))
   for column across permutation
   and row = 0 then (1+ row)
   do (setf (aref result column) (aref vector row))
   finally return result))

(defun left-permute-vector (permutation vector)
  "Permute the column vector to create a row vector."
  (loop
   with result =
   (make-array
    (length vector)
    :element-type (array-element-type vector))
   for column across permutation
   and row = 0 then (1+ row)
   do (setf (aref result row) (aref vector column))
   finally return result))

;;; Permute Arrays

(defun right-permute-array (array permutation)
  "Permute the columns of the array."
  (loop
   with m-rows = (array-dimension array 0)
   with result =
   (make-array
    (array-dimensions array)
    :element-type (array-element-type array))
   for column across permutation
   and row = 0 then (1+ row)
   do (loop
       for irow below m-rows do
       (setf (aref result irow column) (aref array irow row)))
   finally return result))

(defun left-permute-array (permutation array)
  "Permute the rows of the array."
  (loop
   with n-columns = (array-dimension array 1)
   with result =
   (make-array
    (array-dimensions array)
    :element-type (array-element-type array))
   for column across permutation
   and row = 0 then (1+ row)
   do (loop
       for icol below n-columns do
       (setf (aref result row icol) (aref array column icol)))
   finally return result))
