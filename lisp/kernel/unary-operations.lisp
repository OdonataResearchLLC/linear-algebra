#|

 Linear Algebra Binary Operations Kernel

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

;;; Unary operations

(defun sumsq-vector (vector scale sumsq)
  "Return the scaling parameter and the sum of the squares of the
vector."
  (let ((abs-val))
    (dotimes (index (length vector) (values scale sumsq))
      (when (plusp (setq abs-val (abs (aref vector index))))
        (if (< scale abs-val)
            (progn
              (setq sumsq (1+ (* sumsq (expt (/ scale abs-val) 2))))
              (setq scale abs-val))
            (setq
             sumsq
             (+ sumsq (expt (/ (aref vector index) scale) 2))))))))

(defun sumsq-array (array scale sumsq)
    "Return the scaling parameter and the sum of the squares of the
array."
  (let ((m-rows (array-dimension array 0))
        (n-columns (array-dimension array 1))
        (abs-val 0))
    (dotimes (row m-rows (values scale sumsq))
      (dotimes (column n-columns)
        (when (plusp (setq abs-val (abs (aref array row column))))
          (if (< scale abs-val)
              (progn
                (setq sumsq (1+ (* sumsq (expt (/ scale abs-val) 2))))
                (setq scale abs-val))
              (setq sumsq (+ sumsq (expt (/ abs-val scale) 2)))))))))

(defun sump-vector (vector p scale sump)
  "Return the scaling parameter and the sum of the powers of p of the
vector."
  (let ((abs-val))
    (dotimes (index (length vector) (values scale sump))
      (when (plusp (setq abs-val (abs (aref vector index))))
        (if (< scale abs-val)
            (progn
              (setq sump (1+ (* sump (expt (/ scale abs-val) p))))
              (setq scale abs-val))
            (setq
             sump
             (+ sump (expt (/ (aref vector index) scale) p))))))))

(defun sump-array (array p scale sump)
  "Return the scaling parameter and the sum of the P powers of the
matrix."
  (unless (plusp p) (error "The power(~A) must be positive." p))
  (let ((m-rows (array-dimension array 0))
        (n-columns (array-dimension array 1))
        (abs-val 0))
    (dotimes (row m-rows (values scale sump))
      (dotimes (column n-columns)
        (when (plusp (setq abs-val (abs (aref array row column))))
          (if (< scale abs-val)
              (progn
                (setq sump (1+ (* sump (expt (/ scale abs-val) p))))
                (setq scale abs-val))
              (setq
               sump
               (+ sump
                  (expt (/ (aref array row column) scale) p)))))))))
