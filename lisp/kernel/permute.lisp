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

;;; Interface

(defgeneric right-permute (vector-or-array permutation)
  (:documentation
   "Permute the columns of the vector or array."))

(defgeneric left-permute (permutation vector-or-array)
  (:documentation
   "Permute the rows of the vector or array."))

(defgeneric right-npermute (vector-or-array permutation)
  (:documentation
   "Destructively permute teh columns of the vector or array."))

(defgeneric left-npermute (permutation vector-or-array)
  (:documentation
   "Destructively permute the rows fo the vector or array."))

;;; Permute Vectors

(defmethod right-permute ((data vector) (permutation vector))
  (loop with result =
        (make-array
         (length data)
         :element-type (array-element-type data))
        for row = 0 then (1+ row)
        and column across permutation
        do (setf (aref result column) (aref data row))
        finally (return result)))

(defmethod left-permute ((permutation vector) (data vector))
  (loop with result =
        (make-array
         (length data)
         :element-type (array-element-type data))
        for row = 0 then (1+ row)
        and column across permutation
        do (setf (aref result row) (aref data column))
        finally (return result)))

(defmethod right-npermute ((data vector) (permutation vector))
  (loop with end = (1- (length permutation))
        for row = 0 then (if (= row column) (1+ row) row)
        as column = (aref permutation row)
        until (= row end) unless (= row column) do
        (rotatef (aref data row) (aref data column))
        (rotatef (aref permutation row) (aref permutation column))
        finally (return data)))

(defmethod left-npermute ((permutation vector) (data vector))
  (loop with end = (1- (length permutation))
        for row = 0 then (if (= row column) (1+ row) row)
        as column = (aref permutation row)
        until (= row end) unless (= row column) do
        (rotatef (aref data row) (aref data column))
        (rotatef (aref permutation row) (aref permutation column))
        finally (return data)))
