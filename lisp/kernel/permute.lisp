#|

 Linear Algebra Permutation Kernel

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
   do
   (loop
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
   do
   (loop
    for icol below n-columns do
    (setf (aref result row icol) (aref array column icol)))
   finally return result))
