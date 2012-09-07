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
  (loop with result =
        (make-array
         (length vector)
         :element-type (array-element-type vector))
        for row = 0 then (1+ row)
        and column across permutation
        do (setf (aref result column) (aref vector row))
        finally (return result)))

(defun left-permute-vector (permutation vector)
  (loop with result =
        (make-array
         (length vector)
         :element-type (array-element-type vector))
        for row = 0 then (1+ row)
        and column across permutation
        do (setf (aref result row) (aref vector column))
        finally (return result)))

(defun right-npermute-vector (vector permutation)
  (loop with end = (1- (length permutation))
        for row = 0 then (if (= row column) (1+ row) row)
        as column = (aref permutation row)
        until (= row end) unless (= row column) do
        (rotatef (aref vector row) (aref vector column))
        (rotatef (aref permutation row) (aref permutation column))
        finally (return vector)))

(defun left-npermute-vector (permutation vector)
  (loop with end = (1- (length permutation))
        for row = 0 then (if (= row column) (1+ row) row)
        as column = (aref permutation row)
        until (= row end) unless (= row column) do
        (rotatef (aref vector row) (aref vector column))
        (rotatef (aref permutation row) (aref permutation column))
        finally (return vector)))

;;; Permute Arrays

(defun right-permute-array (array permutation)
  (loop with m-rows = (array-dimension array 0)
        with result =
        (make-array
         (array-dimensions array)
         :element-type (array-element-type array))
        for row = 0 then (1+ row)
        and column across permutation
        do (loop for irow below m-rows do
                 (setf
                  (aref result irow column)
                  (aref array irow row)))
        finally (return result)))

(defun left-permute-array (permutation array)
  (loop with n-columns = (array-dimension array 1)
        with result =
        (make-array
         (array-dimensions array)
         :element-type (array-element-type array))
        for row = 0 then (1+ row)
        and column across permutation
        do (loop for icol below n-columns do
                 (setf
                  (aref result row icol)
                  (aref array column icol)))
        finally (return result)))

(defun right-npermute-array (array permutation)
  (loop with m-rows = (array-dimension array 0)
        with end = (1- (length permutation))
        for row = 0 then (if (= row column) (1+ row) row)
        as column = (aref permutation row)
        until (= row end) unless (= row column) do
        (loop for irow below m-rows do
              (rotatef (aref array irow row) (aref array irow column)))
        (rotatef (aref permutation row) (aref permutation column))
        finally (return array)))

(defun left-npermute-array (permutation array)
  (loop with n-columns = (array-dimension array 1)
        with end = (1- (length permutation))
        for row = 0 then (if (= row column) (1+ row) row)
        as column = (aref permutation row)
        until (= row end) unless (= row column) do
        (loop for icolumn below n-columns do
              (rotatef
               (aref array row icolumn)
               (aref array column icolumn)))
        (rotatef (aref permutation row) (aref permutation column))
        finally (return array)))
