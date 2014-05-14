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

(in-package :linear-algebra-kernel)

;;; Copy each element of the array

(defgeneric copy-array (array)
  (:documentation
   "Return an element-wise copy of the original array."))

(defmethod copy-array ((original vector))
  "Return an element-wise copy of the original vector."
  (let* ((size (length original))
         (new-vector
          (make-array
           size :element-type (array-element-type original))))
    (dotimes (index size new-vector)
      (setf (aref new-vector index) (aref original index)))))

(defmethod copy-array ((original array))
  "Return an element-wise copy of the original array."
  (let ((new-array
         (make-array
          (array-dimensions original)
          :element-type (array-element-type original))))
    (dotimes (row (array-dimension original 0) new-array)
      (dotimes (column (array-dimension original 1))
        (setf
         (aref new-array row column)
         (aref original row column))))))

;;; Class and type utilities

(defun common-class-of (object1 object2 &optional
                        (default-class nil default-class-p))
  "Return the common class of the 2 objects or default-class."
  (let ((class1 (class-of object1))
        (class2 (class-of object2)))
    (cond
      ((eq class1 class2) class1)
      ((subtypep class1 class2) class2)
      ((subtypep class2 class1) class1)
      (default-class-p (find-class default-class))
      (t (error "No common or default class.")))))

(defun common-array-element-type (array1 array2)
  "Return the array type common to both arrays."
  (let ((type1 (array-element-type array1))
        (type2 (array-element-type array2)))
    (cond
     ((eq type1 type2) type1)
     ((subtypep type1 type2) type1)
     ((subtypep type2 type1) type2)
     ((and (subtypep type1 'number) (subtypep type2 'number))
      (upgraded-array-element-type
       (type-of (+ (coerce 1 type1) (coerce 1 type2)))))
     (t))))

;;; Equality predicates

;;; (COMPLEX-EQUAL number1 number2) => true or false
(defun complex-equal (complex1 complex2 &optional (epsilon *epsilon*))
  "Return true if both numbers are complex and equal."
  (cond
    ((or
      (typep complex1 '(complex float))
      (typep complex2 '(complex float)))
     (float-equal complex1 complex2 epsilon))
    ((or
      (typep complex1 '(complex integer))
      (typep complex2 '(complex integer)))
     (= complex1 complex2))
    (t (error "Arguments are not complex."))))

;;; (NUMBER-EQUAL number1 number2) => true or false
(defun number-equal (number1 number2 &optional (epsilon *epsilon*))
  "Return true if the numbers are equal using the appropriate
comparison."
  (cond
    ((or (floatp number1) (floatp number2))
     (float-equal number1 number2 epsilon))
    ((and (rationalp number1) (rationalp number2))
     (= number1 number2))
    ((or
      (typep number1 '(complex float))
      (typep number2 '(complex float)))
     (float-equal number1 number2 epsilon))
    ((and
      (typep number1 '(complex rational))
      (typep number2 '(complex rational)))
     (= number1 number2))
    (t (error "Non-numeric arguments."))))
