#|

  Linear Algebra in Common Lisp Unit Tests

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

(in-package :linear-algebra-test)

(define-test make-data-vector
  (:tag :data-vector :make-vector)
  ;; A default vector.
  (assert-true
   (typep
    (linear-algebra:make-vector 10)
    'linear-algebra:data-vector))
  (assert-true
   (typep
    (linear-algebra:make-vector 10)
    'linear-algebra:column-vector))
  (assert-rational-equal
   (vector 0 0 0 0 0 0 0 0 0 0)
   (linear-algebra:make-vector 10))
  ;; Specify the vector type
  (assert-true
   (typep
    (linear-algebra:make-vector
     10 :vector-type 'linear-algebra:column-vector)
    'linear-algebra:column-vector))
  (assert-true
   (typep
    (linear-algebra:make-vector
     10 :vector-type 'linear-algebra:row-vector)
    'linear-algebra:row-vector))
  ;; Specify the element type.
  (assert-true
   (subtypep
    (array-element-type
     (linear-algebra::contents
      (linear-algebra:make-vector
       10 :element-type 'single-float)))
    (upgraded-array-element-type 'single-float)))
  (assert-true
   (subtypep
    (array-element-type
     (linear-algebra::contents
      (linear-algebra:make-vector
       10 :element-type '(complex single-float))))
    (upgraded-array-element-type '(complex single-float))))
  ;; Specify the initial element.
  (assert-float-equal
   (vector 1.0 1.0 1.0 1.0 1.0)
   (linear-algebra:make-vector 5 :initial-element 1.0))
  (assert-float-equal
   (vector
    #C(1.0 2.0) #C(1.0 2.0) #C(1.0 2.0) #C(1.0 2.0))
   (linear-algebra:make-vector 4 :initial-element #C(1.0 2.0)))
  ;; Specify the initial contents.
  (let ((data (make-random-list 100 1.0)))
    (assert-float-equal
     (make-array 100 :initial-contents data)
     (linear-algebra:make-vector 100 :initial-contents data))))

(define-test column-vector
  (:tag :data-vector :column-vector)
  ;; Column vector
  (assert-true
   (typep
    (linear-algebra:column-vector 1.0 2.0 3.0 4.0 5.0)
    'linear-algebra:column-vector))
  ;; Contents
  (assert-float-equal
   (vector 1.0 2.0 3.0 4.0 5.0)
   (linear-algebra:column-vector 1.0 2.0 3.0 4.0 5.0)))

(define-test row-vector
  (:tag :data-vector :row-vector)
  ;; row vector
  (assert-true
   (typep
    (linear-algebra:row-vector 1.0 2.0 3.0 4.0 5.0)
    'linear-algebra:row-vector))
  ;; Contents
  (assert-float-equal
   (vector 1.0 2.0 3.0 4.0 5.0)
   (linear-algebra:row-vector 1.0 2.0 3.0 4.0 5.0)))

(define-test column-vector-p
  (:tag :data-vector :column-vector)
  (assert-true
   (linear-algebra:column-vector-p
    (linear-algebra:column-vector 1.0 2.0 3.0 4.0 5.0)))
  (assert-false
   (linear-algebra:column-vector-p
    (linear-algebra:row-vector 1 2 3 4 5))))

(define-test row-vector-p
  (:tag :data-vector :row-vector)
  (assert-true
   (linear-algebra:row-vector-p
    (linear-algebra:row-vector 1.0 2.0 3.0 4.0 5.0)))
  (assert-false
   (linear-algebra:row-vector-p
    (linear-algebra:column-vector 1 2 3 4 5))))

(define-test data-vector-in-bounds-p
  (:tag :data-vector :vector-in-bounds-p)
  ;; Column vector bounds
  (let ((cvec (linear-algebra:column-vector 1 2 3 4 5)))
    (assert-true
     (linear-algebra:vector-in-bounds-p cvec 0))
    (assert-true
     (linear-algebra:vector-in-bounds-p cvec 4))
    (assert-true
     (linear-algebra:vector-in-bounds-p cvec (random 5)))
    (assert-false
     (linear-algebra:vector-in-bounds-p cvec 10)))
  ;; Row vector bounds
  (let ((rvec (linear-algebra:column-vector 1 2 3 4 5)))
    (assert-true
     (linear-algebra:vector-in-bounds-p rvec 0))
    (assert-true
     (linear-algebra:vector-in-bounds-p rvec 4))
    (assert-true
     (linear-algebra:vector-in-bounds-p rvec (random 5)))
    (assert-false
     (linear-algebra:vector-in-bounds-p rvec 10))))

(define-test data-vector-element-type
  (:tag :data-vector :vector-element-type)
  ;; Column vector
  (assert-true
   (subtypep
    (linear-algebra:vector-element-type
     (linear-algebra:make-vector
      5 :element-type 'single-float
      :vector-type 'linear-algebra:column-vector))
    'single-float))
  ;; Row vector
  (assert-true
   (subtypep
    (linear-algebra:vector-element-type
     (linear-algebra:make-vector
      5 :element-type 'single-float
      :vector-type 'linear-algebra:row-vector))
    'single-float)))

(define-test data-vector-length
  (:tag :data-vector :vector-length)
  (assert-eq
   5 (linear-algebra:vector-length
      (linear-algebra:column-vector 1 2 3 4 5)))
  (assert-eq
   5 (linear-algebra:vector-length
      (linear-algebra:row-vector 1 2 3 4 5))))

(define-test data-vref
  (:tag :data-vector :vref)
  ;; Reference the element of a vector.
  (let* ((size 1000)
         (index (random-interior-index size))
         (end   (1- size))
         (list   (make-random-list size 100.0))
         (vector
          (linear-algebra:make-vector
           size :initial-contents list)))
    (assert-float-equal
     (nth 0 list) (linear-algebra:vref vector 0))
    (assert-float-equal
     (nth index list) (linear-algebra:vref vector index))
    (assert-float-equal
     (nth end list) (linear-algebra:vref vector end)))
  ;; Setting an element of a vector.
  (let* ((size 1000)
         (val0  (random 100.0))
         (vali  (random 100.0))
         (valn  (random 100.0))
         (index (random-interior-index size))
         (end   (1- size))
         (list   (make-random-list size 100.0))
         (vector
          (linear-algebra:make-vector
           size :initial-contents list)))
    (setf (linear-algebra:vref vector 0)     val0)
    (setf (linear-algebra:vref vector index) vali)
    (setf (linear-algebra:vref vector end)   valn)
    (assert-float-equal
     val0 (aref (linear-algebra::contents vector) 0))
    (assert-float-equal
     vali (aref (linear-algebra::contents vector) index))
    (assert-float-equal
     valn (aref (linear-algebra::contents vector) end)))
  ;; VREF only works for vectors.
  (let ((list (make-list 1000 :initial-element 10.0))
        (index (random 1000)))
    (assert-error
     'error (linear-algebra:vref list index))
    (assert-error
     'error (setf (linear-algebra:vref list index) 42.0))))

(define-test copy-data-vector
  (:tag :data-vector :copy-vector)
  ;; Column vector
  (let* ((vec (linear-algebra:column-vector 1 2 3 4 5))
         (new (linear-algebra:copy-vector vec)))
    (assert-true (subtypep (type-of new) (type-of vec)))
    (assert-false (eq vec new))
    (assert-false
     (eq (linear-algebra::contents vec)
         (linear-algebra::contents new)))
    (assert-rational-equal vec new))
  ;; Row vector
  (let* ((vec (linear-algebra:row-vector 1 2 3 4 5))
         (new (linear-algebra:copy-vector vec)))
    (assert-true (subtypep (type-of new) (type-of vec)))
    (assert-false (eq vec new))
    (assert-false
     (eq (linear-algebra::contents vec)
         (linear-algebra::contents new)))
    (assert-rational-equal vec new)))

(define-test data-subvector
  (:tag :data-vector :subvector)
  ;; Get subvectors
  (let* ((data (list 0 1 2 3 4 5 6 7 8 9))
         (vec
          (linear-algebra:column-vector
           0 1 2 3 4 5 6 7 8 9))
         (index (random-interior-index 10))
         (subvec0
          (linear-algebra:subvector vec index))
         (subvec1
          (linear-algebra:subvector vec 0 index)))
    (assert-true
     (subtypep (type-of vec) (type-of subvec0)))
    (assert-true
     (subtypep (type-of vec) (type-of subvec1)))
    (assert-rational-equal
     (subseq data index) subvec0)
    (assert-rational-equal
     (subseq data 0 index) subvec1))
  ;; Set subvectors
  (let ((vec (linear-algebra:make-vector 10))
        (data
         (linear-algebra:make-vector
          5 :initial-element 1)))
    (assert-rational-equal
     #(1 1 1 1 1)
     (setf (linear-algebra:subvector vec 7) data))
    (assert-rational-equal
     '(0 0 0 0 0 0 0 1 1 1) vec)))

(define-test replace-data-vector
  (:tag :data-vector :replace-vector)
  (flet ((zero-vector () (linear-algebra:make-vector 10)))
    (let ((data
           (linear-algebra:make-vector
            5 :initial-element 1)))
      (assert-rational-equal
       #(1 1 1 1 1 0 0 0 0 0)
        (linear-algebra:replace-vector (zero-vector) data))
      (assert-rational-equal
       #(0 0 1 1 1 1 1 0 0 0)
       (linear-algebra:replace-vector
        (zero-vector) data :start1 2))
      (assert-rational-equal
       #(0 0 0 0 0 0 0 0 1 1)
       (linear-algebra:replace-vector
        (zero-vector) data :start1 8))
      (assert-rational-equal
       #(1 1 1 0 0 0 0 0 0 0)
       (linear-algebra:replace-vector
        (zero-vector) data :end1 3))
      (assert-rational-equal
       #(0 0 1 1 0 0 0 0 0 0)
       (linear-algebra:replace-vector
        (zero-vector) data :start1 2 :end1 4))
      (assert-rational-equal
       #(1 1 1 0 0 0 0 0 0 0)
       (linear-algebra:replace-vector
        (zero-vector) data :start2 2))
      (assert-rational-equal
       #(1 1 0 0 0 0 0 0 0 0)
       (linear-algebra:replace-vector
        (zero-vector) data :end2 2))
      (assert-rational-equal
       #(1 1 1 0 0 0 0 0 0 0)
       (linear-algebra:replace-vector
        (zero-vector) data :start2 1 :end2 4))
      (assert-rational-equal
       #(0 0 1 1 1 0 0 0 0 0)
       (linear-algebra:replace-vector
        (zero-vector) data :start1 2 :start2 2))
      (assert-rational-equal
       #(0 0 1 1 1 0 0 0 0 0)
       (linear-algebra:replace-vector
        (zero-vector) data :start1 2 :end2 3)))))

(define-test %map-data-vector
  (:tag :data-vector :map-vector)
  (let ((col-data (linear-algebra:column-vector 1 1 1 1 1))
        (row-data (linear-algebra:row-vector 1 1 1 1 1)))
    (assert-true
     (typep
      (linear-algebra::%map-data-vector
       'linear-algebra:column-vector
       #'identity col-data)
      'linear-algebra:column-vector))
    (assert-true
     (typep
      (linear-algebra::%map-data-vector
       'linear-algebra:row-vector
       #'identity col-data)
      'linear-algebra:row-vector))
    (assert-true
     (typep
      (linear-algebra::%map-data-vector
       'linear-algebra:column-vector
       #'identity row-data)
      'linear-algebra:column-vector))
    (assert-true
     (typep
      (linear-algebra::%map-data-vector
       'linear-algebra:row-vector
       #'identity row-data)
      'linear-algebra:row-vector))
    (assert-rational-equal
     '(3 3 3 3 3)
     (linear-algebra::%map-data-vector
      'linear-algebra:column-vector
      #'+ col-data col-data col-data))
    (assert-rational-equal
     '(3 3 3 3 3)
     (linear-algebra::%map-data-vector
      'linear-algebra:row-vector
      #'+ row-data row-data row-data))
    (assert-rational-equal
     '(4 4 4 4 4)
     (linear-algebra::%map-data-vector
      'linear-algebra:column-vector
      #'+ col-data row-data col-data row-data))))

(define-test map-data-vector
  (:tag :data-vector :map-vector)
  (let ((col-data (linear-algebra:column-vector 1 1 1 1 1))
        (row-data (linear-algebra:row-vector 1 1 1 1 1)))
    (assert-error
     'error (linear-algebra:map-vector 'vector #'identity col-data))
    (assert-error
     'error
     (linear-algebra:map-vector
      'linear-algebra:column-vector #'identity col-data #(2 2 2 2 2)))
    (assert-rational-equal
     '(3 3 3 3 3)
     (linear-algebra:map-vector
      'linear-algebra:column-vector
      #'+ col-data row-data col-data))))

(define-test %map-into-data-vector
  (:tag :data-vector :map-into-vector)
  (let ((col-data (linear-algebra:column-vector 1 1 1 1 1))
        (row-data (linear-algebra:row-vector 1 1 1 1 1)))
    (assert-true
     (typep
      (linear-algebra::%map-into-data-vector
       col-data #'identity col-data)
      'linear-algebra:column-vector))
    (assert-true
     (typep
      (linear-algebra::%map-into-data-vector
       row-data #'identity col-data)
      'linear-algebra:row-vector))
    (assert-true
     (typep
      (linear-algebra::%map-into-data-vector
       col-data #'identity row-data)
      'linear-algebra:column-vector))
    (assert-true
     (typep
      (linear-algebra::%map-into-data-vector
       row-data #'identity row-data)
      'linear-algebra:row-vector))
    (assert-rational-equal
     '(3 3 3 3 3)
     (linear-algebra::%map-into-data-vector
      (linear-algebra:column-vector 0 0 0 0 0)
      #'+ col-data col-data col-data))
    (assert-rational-equal
     '(3 3 3 3 3)
     (linear-algebra::%map-into-data-vector
      (linear-algebra:row-vector 0 0 0 0 0)
      #'+ row-data row-data row-data))
    (assert-rational-equal
     '(4 4 4 4 4)
     (linear-algebra::%map-into-data-vector
      col-data #'+ col-data row-data col-data row-data))))

(define-test map-into-data-vector
  (:tag :data-vector :map-into-vector)
  (let ((col-data (linear-algebra:column-vector 1 1 1 1 1))
        (row-data (linear-algebra:row-vector 1 1 1 1 1)))
    (assert-error
     'error
     (linear-algebra:map-into-vector
      #(0 0 0 0 0) #'identity col-data))
    (assert-error
     'error
     (linear-algebra:map-vector
      (linear-algebra:column-vector 0 0 0 0 0) #'identity col-data))
    (assert-rational-equal
     '(3 3 3 3 3)
     (linear-algebra:map-into-vector
      col-data #'+ col-data row-data col-data))))

;;; Apply rotation

(define-test apply-rotation-data-vector
  (:tag :data-vector :apply-rotation)
  ;; Float
  (let ((vec1
         (linear-algebra:make-vector
          5 :element-type 'single-float
          :initial-contents '(3.0 5.0 7.0 11.0 13.0)))
        (vec2
         (linear-algebra:make-vector
          5 :element-type 'single-float
          :initial-contents '(4.0 9.0 16.0 25.0 36.0))))
    ;; Position 0
    (multiple-value-bind (cc ss rr)
        (linear-algebra-kernel:givens-rotation
         (linear-algebra:vref vec1 0) (linear-algebra:vref vec2 0))
      (multiple-value-bind (rvec1 rvec2)
          (linear-algebra:apply-rotation vec1 vec2 cc ss)
        (assert-float-equal rr (linear-algebra:vref rvec1 0))
        (assert-float-equal
         #(5.0 10.200001 17.0 26.6 36.600002) rvec1)
        (assert-float-equal
         #(0.0 1.4000001 4.0000005 6.200001 11.2) rvec2)))
    ;; Position 1
    (multiple-value-bind (cc ss rr)
        (linear-algebra-kernel:givens-rotation
         (linear-algebra:vref vec1 1) (linear-algebra:vref vec2 1))
      (multiple-value-bind (rvec1 rvec2)
          (linear-algebra:apply-rotation vec1 vec2 cc ss)
        (assert-float-equal rr (linear-algebra:vref rvec1 1))
        (assert-float-equal
         #(4.953558 10.2956295 17.386017 27.196003 37.78302) rvec1)
        (assert-float-equal
         #(-0.67990017 0.0 1.651186 2.525344 6.1191006) rvec2)))
    ;; Position 2
    (multiple-value-bind (cc ss rr)
        (linear-algebra-kernel:givens-rotation
         (linear-algebra:vref vec1 2) (linear-algebra:vref vec2 2))
      (multiple-value-bind (rvec1 rvec2)
          (linear-algebra:apply-rotation vec1 vec2 cc ss)
        (assert-float-equal rr (linear-algebra:vref rvec1 2))
        (assert-float-equal
         #(4.867086 10.249511 17.46425 27.312943 38.19231) rvec1)
        (assert-float-equal
         #(-1.1451968 -0.97341704 0.0 -0.05725956 2.519433) rvec2)))
    ;; Position 3
    (multiple-value-bind (cc ss rr)
        (linear-algebra-kernel:givens-rotation
         (linear-algebra:vref vec1 3) (linear-algebra:vref vec2 3))
      (multiple-value-bind (rvec1 rvec2)
          (linear-algebra:apply-rotation vec1 vec2 cc ss)
        (assert-float-equal rr (linear-algebra:vref rvec1 3))
        (assert-float-equal
         #(4.8694763 10.251528 17.46421 27.313 38.186943) rvec1)
        (assert-float-equal
         #(-1.1349905 -0.95192766 0.036612988 9.536743e-7 2.599495)
         rvec2)))
    ;; Position 4
    (multiple-value-bind (cc ss rr)
        (linear-algebra-kernel:givens-rotation
         (linear-algebra:vref vec1 4) (linear-algebra:vref vec2 4))
      (multiple-value-bind (rvec1 rvec2)
          (linear-algebra:apply-rotation vec1 vec2 cc ss)
        ;; FIXME : The error should not be this large.
        (let ((*epsilon* (* 32 single-float-epsilon)))
          (assert-float-equal rr (linear-algebra:vref rvec1 4))
          (assert-float-equal
           #(4.781149 10.163207 17.426374 27.249937 38.27532) rvec1)
          (assert-float-equal
           #(-1.4630839 -1.6459692 -1.1495662 -1.8549814 0.0) rvec2)))))
  ;; Complex
  (let ((vec1
         (linear-algebra:make-vector
          2 :element-type '(complex single-float)
          :initial-contents #(#C(3.0 5.0) #C(7.0 11.0))))
        (vec2
         (linear-algebra:make-vector
          2 :element-type '(complex single-float)
          :initial-contents #(#C(4.0 9.0) #C(16.0 25.0)))))
    ;; Position 0
    (multiple-value-bind (cc ss rr)
        (linear-algebra-kernel:givens-rotation
         (linear-algebra:vref vec1 0) (linear-algebra:vref vec2 0))
      (multiple-value-bind (rvec1 rvec2)
          (linear-algebra:apply-rotation vec1 vec2 cc ss)
        (assert-float-equal rr (linear-algebra:vref rvec1 0))
        (assert-float-equal
         #(#C(5.888673 9.814455) #C(19.85367 25.27784)) rvec1)
        (assert-float-equal
         #(#C(-4.7683716e-7 -4.7683716e-7) #C(3.326425 2.6071978))
         rvec2)))
    ;; Position 1
    (multiple-value-bind (cc ss rr)
        (linear-algebra-kernel:givens-rotation
         (linear-algebra:vref vec1 1) (linear-algebra:vref vec2 1))
      (multiple-value-bind (rvec1 rvec2)
          (linear-algebra:apply-rotation vec1 vec2 cc ss)
        ;; FIXME : The error should not be this large.
        (let ((*epsilon* (* 32 single-float-epsilon)))
          (assert-float-equal rr (linear-algebra:vref rvec1 1))
          (assert-float-equal
           #(#C(4.8474817 10.2603855) #C(17.405037 27.350775)) rvec1)
          (assert-float-equal
           #(#C(-1.1497686 -0.95104337) #C(0.0 0.0)) rvec2
           (linear-algebra::contents rvec2)
           :givens-rotation-1))))))

;;; Destructively apply rotation

(define-test napply-rotation-data-vector
  (:tag :data-vector :napply-rotation)
  ;; Float, position 0
  (let ((vec1
         (linear-algebra:make-vector
          5 :element-type 'single-float
          :initial-contents #(3.0 5.0 7.0 11.0 13.0)))
        (vec2
         (linear-algebra:make-vector
          5 :element-type 'single-float
          :initial-contents #(4.0 9.0 16.0 25.0 36.0))))
    (multiple-value-bind (cc ss rr)
        (linear-algebra-kernel:givens-rotation
         (linear-algebra:vref vec1 0) (linear-algebra:vref vec2 0))
      (multiple-value-bind (rvec1 rvec2)
          (linear-algebra:napply-rotation vec1 vec2 cc ss)
        (assert-eq vec1 rvec1)
        (assert-eq vec2 rvec2)
        (assert-float-equal rr (linear-algebra:vref vec1 0))
        (assert-float-equal
         #(5.0 10.200001 17.0 26.6 36.600002) vec1)
        (assert-float-equal
         #(0.0 1.4000001 4.0000005 6.200001 11.2) vec2))))
  ;; Float, position 1
  (let ((vec1
         (linear-algebra:make-vector
          5 :element-type 'single-float
          :initial-contents #(3.0 5.0 7.0 11.0 13.0)))
        (vec2
         (linear-algebra:make-vector
          5 :element-type 'single-float
          :initial-contents #(4.0 9.0 16.0 25.0 36.0))))
    (multiple-value-bind (cc ss rr)
        (linear-algebra-kernel:givens-rotation
         (linear-algebra:vref vec1 1) (linear-algebra:vref vec2 1))
      (multiple-value-bind (rvec1 rvec2)
          (linear-algebra:napply-rotation vec1 vec2 cc ss)
        (assert-eq vec1 rvec1)
        (assert-eq vec2 rvec2)
        (assert-float-equal rr (linear-algebra:vref vec1 1))
        (assert-float-equal
         #(4.953558 10.2956295 17.386017 27.196003 37.78302) vec1)
        (assert-float-equal
         #(-0.67990017 0.0 1.651186 2.525344 6.1191006) vec2))))
  ;; Float, position 2
  (let ((vec1
         (linear-algebra:make-vector
          5 :element-type 'single-float
          :initial-contents #(3.0 5.0 7.0 11.0 13.0)))
        (vec2
         (linear-algebra:make-vector
          5 :element-type 'single-float
          :initial-contents #(4.0 9.0 16.0 25.0 36.0))))
    (multiple-value-bind (cc ss rr)
        (linear-algebra-kernel:givens-rotation
         (linear-algebra:vref vec1 2) (linear-algebra:vref vec2 2))
      (multiple-value-bind (rvec1 rvec2)
          (linear-algebra:napply-rotation vec1 vec2 cc ss)
        (assert-eq vec1 rvec1)
        (assert-eq vec2 rvec2)
        (assert-float-equal rr (linear-algebra:vref vec1 2))
        (assert-float-equal
         #(4.867086 10.249511 17.46425 27.312943 38.19231) vec1)
        (assert-float-equal
         #(-1.1451968 -0.97341704 0.0 -0.05725956 2.519433) vec2))))
  ;; Float, position 3
  (let ((vec1
         (linear-algebra:make-vector
          5 :element-type 'single-float
          :initial-contents #(3.0 5.0 7.0 11.0 13.0)))
        (vec2
         (linear-algebra:make-vector
          5 :element-type 'single-float
          :initial-contents #(4.0 9.0 16.0 25.0 36.0))))
    (multiple-value-bind (cc ss rr)
        (linear-algebra-kernel:givens-rotation
         (linear-algebra:vref vec1 3) (linear-algebra:vref vec2 3))
      (multiple-value-bind (rvec1 rvec2)
          (linear-algebra:napply-rotation vec1 vec2 cc ss)
        (assert-eq vec1 rvec1)
        (assert-eq vec2 rvec2)
        (assert-float-equal rr (linear-algebra:vref vec1 3))
        (assert-float-equal
         #(4.8694763 10.251528 17.46421 27.313 38.186943) vec1)
        (assert-float-equal
         #(-1.1349905 -0.95192766 0.036612988 9.536743e-7 2.599495) vec2))))
  ;; Float, position 4
  (let ((vec1
         (linear-algebra:make-vector
          5 :element-type 'single-float
          :initial-contents #(3.0 5.0 7.0 11.0 13.0)))
        (vec2
         (linear-algebra:make-vector
          5 :element-type 'single-float
          :initial-contents #(4.0 9.0 16.0 25.0 36.0))))
    (multiple-value-bind (cc ss rr)
        (linear-algebra-kernel:givens-rotation
         (linear-algebra:vref vec1 4) (linear-algebra:vref vec2 4))
      (multiple-value-bind (rvec1 rvec2)
          (linear-algebra:napply-rotation vec1 vec2 cc ss)
        (let ((*epsilon* (* 32 single-float-epsilon)))
          (assert-eq vec1 rvec1)
          (assert-eq vec2 rvec2)
          (assert-float-equal rr (linear-algebra:vref vec1 4))
          (assert-float-equal
           #(4.781149 10.163207 17.426374 27.249937 38.27532) vec1)
          (assert-float-equal
           #(-1.4630839 -1.6459692 -1.1495662 -1.8549814 0.0) vec2)))))
  ;; Complex, position 0
  (let ((vec1
         (linear-algebra:make-vector
          2 :element-type '(complex single-float)
          :initial-contents #(#C(3.0 5.0) #C(7.0 11.0))))
        (vec2
         (linear-algebra:make-vector
          2 :element-type '(complex single-float)
          :initial-contents #(#C(4.0 9.0) #C(16.0 25.0)))))
    (multiple-value-bind (cc ss rr)
        (linear-algebra-kernel:givens-rotation
         (linear-algebra:vref vec1 0) (linear-algebra:vref vec2 0))
      (multiple-value-bind (rvec1 rvec2)
          (linear-algebra:napply-rotation vec1 vec2 cc ss)
        (assert-eq vec1 rvec1)
        (assert-eq vec2 rvec2)
        (assert-float-equal rr (linear-algebra:vref vec1 0))
        (assert-float-equal
         #(#C(5.888673 9.814455) #C(19.85367 25.27784)) vec1)
        (assert-float-equal
         #(#C(-4.7683716e-7 -4.7683716e-7) #C(3.326425 2.6071978)) vec2))))
  ;; Complex, position 1
  (let ((vec1
         (linear-algebra:make-vector
          2 :element-type '(complex single-float)
          :initial-contents #(#C(3.0 5.0) #C(7.0 11.0))))
        (vec2
         (linear-algebra:make-vector
          2 :element-type '(complex single-float)
          :initial-contents #(#C(4.0 9.0) #C(16.0 25.0)))))
    (multiple-value-bind (cc ss rr)
        (linear-algebra-kernel:givens-rotation
         (linear-algebra:vref vec1 1) (linear-algebra:vref vec2 1))
      (multiple-value-bind (rvec1 rvec2)
          (linear-algebra:napply-rotation vec1 vec2 cc ss)
        (assert-eq vec1 rvec1)
        (assert-eq vec2 rvec2)
        ;; FIXME : Why is this error so large?
        (let ((*epsilon* (* 32 single-float-epsilon)))
          (assert-float-equal rr (linear-algebra:vref vec1 1))
          (assert-float-equal
           #(#C(4.8474817 10.2603855) #C(17.405037 27.350775)) vec1)
          (assert-float-equal
           #(#C(-1.1497688 -0.9510431) #C(0.0 0.0)) vec2))))))

;;; Fundamental operations

;;; Data vector norm

(define-test norm-data-vector
  (:tag :data-vector :norm)
  (let ((data
         (linear-algebra:column-vector
          -6 -5 -4 -3 -2 -1 0 1 2 3 4 5))
        (zdata
         (linear-algebra:column-vector
          #C(1 0) #C(3 1) #C(2 3) #C(0 4)
          #C(-2 3) #C(-3 1) #C(-1 0))))
    ;; Taxicab norm
    (assert-rational-equal
     36 (linear-algebra:norm data))
    (assert-float-equal
     19.535658 (linear-algebra:norm zdata))
    ;; Euclidean norm
    (assert-float-equal
     12.083046 (linear-algebra:norm data 2))
    (assert-float-equal
     8.0 (linear-algebra:norm zdata 2))
    ;; P-norm
    (assert-float-equal
     8.732892 (linear-algebra:norm data 3))
    (assert-float-equal
     6.064035 (linear-algebra:norm zdata 3))
    ;; Data vector Infinity norm
    (assert-rational-equal
     6 (linear-algebra:norm data :infinity))
    (assert-float-equal
     4.0 (linear-algebra:norm zdata :infinity))))

;;; Data vector transpose

(define-test transpose-data-vector
  (:tag :data-vector :transpose)
  (let ((col-data
         (linear-algebra:column-vector
          1.0 2.0 3.0 4.0 5.0))
        (row-data
         (linear-algebra:row-vector
          1.0 2.0 3.0 4.0 5.0)))
    (assert-true
     (typep
      (linear-algebra:transpose col-data)
      'linear-algebra:row-vector))
    (assert-true
     (typep
      (linear-algebra:transpose row-data)
      'linear-algebra:column-vector))
    (assert-float-equal
     col-data (linear-algebra:transpose col-data))
    (assert-float-equal
     row-data (linear-algebra:transpose row-data))))

(define-test ntranspose-data-vector
  (:tag :data-vector :ntranspose)
  (let ((col-data
         (linear-algebra:column-vector
          1.0 2.0 3.0 4.0 5.0))
        (row-data
         (linear-algebra:row-vector
          1.0 2.0 3.0 4.0 5.0)))
    (assert-true
     (typep
      (linear-algebra:ntranspose col-data)
      'linear-algebra:row-vector))
    (assert-true
     (typep col-data 'linear-algebra:row-vector))
    (assert-float-equal
     #(1.0 2.0 3.0 4.0 5.0)
     (linear-algebra:transpose col-data))
    (assert-true
     (typep
      (linear-algebra:ntranspose row-data)
      'linear-algebra:column-vector))
    (assert-true
     (typep row-data 'linear-algebra:column-vector))
    (assert-float-equal
     #(1.0 2.0 3.0 4.0 5.0)
     (linear-algebra:transpose row-data))))

;;; Data vector permutation

(define-test permute-data-vector
  (:tag :data-vector :permute)
  (let ((rvec
         (linear-algebra:row-vector
          1.1 2.2 3.3 4.4 5.5))
        (cvec
         (linear-algebra:column-vector
          1.1 2.2 3.3 4.4 5.5))
        (rerr
         (linear-algebra:row-vector
          1.1 2.2 3.3 4.4 5.5 6.6))
        (cerr
         (linear-algebra:column-vector
          1.1 2.2 3.3 4.4 5.5 6.6))
        (pmat
         (linear-algebra:make-matrix
          5 5 :matrix-type
          'linear-algebra:permutation-matrix
          :initial-contents
          '((0 0 1 0 0)
            (0 0 0 0 1)
            (1 0 0 0 0)
            (0 1 0 0 0)
            (0 0 0 1 0)))))
    (assert-float-equal
     #(3.3 4.4 1.1 5.5 2.2)
     (linear-algebra:permute rvec pmat))
    (assert-float-equal
     #(3.3 5.5 1.1 2.2 4.4)
     (linear-algebra:permute pmat cvec))
    (assert-error
     'error (linear-algebra:permute rerr pmat))
    (assert-error
     'error (linear-algebra:permute pmat cerr))))

;;; Data vector scale

(define-test scale-data-vector
  (:tag :data-vector :scale)
  (assert-float-equal
   #(2.0 4.0 6.0 8.0 10.0)
   (linear-algebra:scale
    2.0 (linear-algebra:column-vector 1.0 2.0 3.0 4.0 5.0)))
  (assert-float-equal
   #(#C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0) #C(4.0 4.0) #C(5.0 5.0))
   (linear-algebra:scale
    #C(1.0 1.0) (linear-algebra:column-vector 1.0 2.0 3.0 4.0 5.0)))
  (assert-float-equal
   #(#C(2.0 2.0) #C(4.0 4.0) #C(6.0 6.0) #C(8.0 8.0) #C(10.0 10.0))
   (linear-algebra:scale
    2.0 (linear-algebra:column-vector
         #C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0) #C(4.0 4.0) #C(5.0 5.0))))
  (assert-float-equal
   #(#C(0.0 4.0) #C(0.0 8.0) #C(0.0 12.0) #C(0.0 16.0) #C(0.0 20.0))
   (linear-algebra:scale
    #C(2.0 2.0)
    (linear-algebra:column-vector
     #C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0) #C(4.0 4.0) #C(5.0 5.0)))))

(define-test nscale-data-vector
  (:tag :data-vector :nscale)
  (assert-float-equal
   #(2.0 4.0 6.0 8.0 10.0)
   (linear-algebra:nscale
    2.0 (linear-algebra:column-vector 1.0 2.0 3.0 4.0 5.0)))
  (assert-float-equal
   #(#C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0) #C(4.0 4.0) #C(5.0 5.0))
   (linear-algebra:nscale
    #C(1.0 1.0) (linear-algebra:column-vector 1.0 2.0 3.0 4.0 5.0)))
  (assert-float-equal
   #(#C(2.0 2.0) #C(4.0 4.0) #C(6.0 6.0) #C(8.0 8.0) #C(10.0 10.0))
   (linear-algebra:nscale
    2.0 (linear-algebra:column-vector
         #C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0) #C(4.0 4.0) #C(5.0 5.0))))
  (assert-float-equal
   #(#C(0.0 4.0) #C(0.0 8.0) #C(0.0 12.0) #C(0.0 16.0) #C(0.0 20.0))
   (linear-algebra:nscale
    #C(2.0 2.0)
    (linear-algebra:column-vector
     #C(1.0 1.0) #C(2.0 2.0) #C(3.0 3.0) #C(4.0 4.0) #C(5.0 5.0)))))

;;; Vector addition

(define-test add-data-vector
  (:tag :data-vector :add)
  ;; Real
  (let ((vec1 (linear-algebra:column-vector 1.1 2.2 3.3 4.4))
        (vec2 (linear-algebra:column-vector 1.1 2.2 3.3 4.4)))
    (assert-float-equal
     #(2.2 4.4 6.6 8.8)
     (linear-algebra:add vec1 vec2))
    (assert-float-equal
     #(3.3 6.6 9.9 13.2)
     (linear-algebra:add vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(3.3 6.6 9.9 13.2)
     (linear-algebra:add vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(4.4 8.8 13.2 17.6)
     (linear-algebra:add vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  (let ((vec1 (linear-algebra:row-vector 1.1 2.2 3.3 4.4))
        (vec2 (linear-algebra:row-vector 1.1 2.2 3.3 4.4)))
    (assert-float-equal
     #(2.2 4.4 6.6 8.8)
     (linear-algebra:add vec1 vec2))
    (assert-float-equal
     #(3.3 6.6 9.9 13.2)
     (linear-algebra:add vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(3.3 6.6 9.9 13.2)
     (linear-algebra:add vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(4.4 8.8 13.2 17.6)
     (linear-algebra:add vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  ;; Complex
  (let ((vec1 (linear-algebra:column-vector #C(1.1 2.2) #C(3.3 4.4)))
        (vec2 (linear-algebra:column-vector #C(1.1 2.2) #C(3.3 4.4))))
    (assert-float-equal
     #(#C(2.2 4.4) #C(6.6 8.8))
     (linear-algebra:add vec1 vec2))
    (assert-float-equal
     #(#C(3.3 6.6) #C(9.9 13.2))
     (linear-algebra:add vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(#C(3.3 6.6) #C(9.9 13.2))
     (linear-algebra:add vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(#C(4.4 8.8) #C(13.2 17.6))
     (linear-algebra:add vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  (let ((vec1 (linear-algebra:row-vector #C(1.1 2.2) #C(3.3 4.4)))
        (vec2 (linear-algebra:row-vector #C(1.1 2.2) #C(3.3 4.4))))
    (assert-float-equal
     #(#C(2.2 4.4) #C(6.6 8.8))
     (linear-algebra:add vec1 vec2))
    (assert-float-equal
     #(#C(3.3 6.6) #C(9.9 13.2))
     (linear-algebra:add vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(#C(3.3 6.6) #C(9.9 13.2))
     (linear-algebra:add vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(#C(4.4 8.8) #C(13.2 17.6))
     (linear-algebra:add vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  ;; Errors
  (assert-error
   'error
   (linear-algebra:add
    (linear-algebra:column-vector 1.1 2.2 3.3 4.4)
    (linear-algebra:row-vector 1.1 2.2 3.3 4.4)))
  (assert-error
   'error
   (linear-algebra:add
    (linear-algebra:row-vector 1.1 2.2 3.3 4.4)
    (linear-algebra:column-vector 1.1 2.2 3.3 4.4))))

;;; Destructive vector addition

(define-test nadd-data-vector
  (:tag :data-vector :nadd)
  ;; Real
  (let ((vec1 (linear-algebra:column-vector 1.1 2.2 3.3 4.4))
        (vec2 (linear-algebra:column-vector 1.1 2.2 3.3 4.4)))
    (assert-eq vec1 (linear-algebra:nadd vec1 vec2))
    (assert-float-equal
     #(2.2 4.4 6.6 8.8) vec1)
    (assert-float-equal
     #(4.4 8.8 13.2 17.6)
     (linear-algebra:nadd vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(9.9 19.8 29.7 39.6)
     (linear-algebra:nadd vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(22.0 44.0 66.0 88.0)
     (linear-algebra:nadd vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  (let ((vec1 (linear-algebra:row-vector 1.1 2.2 3.3 4.4))
        (vec2 (linear-algebra:row-vector 1.1 2.2 3.3 4.4)))
    (assert-eq vec1 (linear-algebra:nadd vec1 vec2))
    (assert-float-equal
     #(2.2 4.4 6.6 8.8) vec1)
    (assert-float-equal
     #(4.4 8.8 13.2 17.6)
     (linear-algebra:nadd vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(9.9 19.8 29.7 39.6)
     (linear-algebra:nadd vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(22.0 44.0 66.0 88.0)
     (linear-algebra:nadd vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  ;; Complex
  (let ((vec1 (linear-algebra:column-vector #C(1.1 2.2) #C(3.3 4.4)))
        (vec2 (linear-algebra:column-vector #C(1.1 2.2) #C(3.3 4.4))))
    (assert-eq vec1 (linear-algebra:nadd vec1 vec2))
    (assert-float-equal
     #(#C(2.2 4.4) #C(6.6 8.8)) vec1)
    (assert-float-equal
     #(#C(4.4 8.8) #C(13.2 17.6))
     (linear-algebra:nadd vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(#C(9.9 19.8) #C(29.7 39.6))
     (linear-algebra:nadd vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(#C(22.0 44.0) #C(66.0 88.0))
     (linear-algebra:nadd vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  (let ((vec1 (linear-algebra:row-vector #C(1.1 2.2) #C(3.3 4.4)))
        (vec2 (linear-algebra:row-vector #C(1.1 2.2) #C(3.3 4.4))))
    (assert-eq vec1 (linear-algebra:nadd vec1 vec2))
    (assert-float-equal
     #(#C(2.2 4.4) #C(6.6 8.8)) vec1)
    (assert-float-equal
     #(#C(4.4 8.8) #C(13.2 17.6))
     (linear-algebra:nadd vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(#C(9.9 19.8) #C(29.7 39.6))
     (linear-algebra:nadd vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(#C(22.0 44.0) #C(66.0 88.0))
     (linear-algebra:nadd vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  ;; Errors
  (assert-error
   'error
   (linear-algebra:nadd
    (linear-algebra:column-vector 1.1 2.2 3.3 4.4)
    (linear-algebra:row-vector 1.1 2.2 3.3 4.4)))
  (assert-error
   'error
   (linear-algebra:nadd
    (linear-algebra:row-vector 1.1 2.2 3.3 4.4)
    (linear-algebra:column-vector 1.1 2.2 3.3 4.4))))

;;; Vector subtraction

(define-test subtract-data-vector
  (:tag :data-vector :subtract)
  ;; Real
  (let ((vec1 (linear-algebra:column-vector 1.1 2.2 3.3 4.4))
        (vec2 (linear-algebra:column-vector 1.1 2.2 3.3 4.4)))
    (assert-float-equal
     #(0.0 0.0 0.0 0.0)
     (linear-algebra:subtract vec1 vec2))
    (assert-float-equal
     #(1.1 2.2 3.3 4.4)
     (linear-algebra:subtract vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(-1.1 -2.2 -3.3 -4.4)
     (linear-algebra:subtract vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(0.0 0.0 0.0 0.0)
     (linear-algebra:subtract vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  (let ((vec1 (linear-algebra:row-vector 1.1 2.2 3.3 4.4))
        (vec2 (linear-algebra:row-vector 1.1 2.2 3.3 4.4)))
    (assert-float-equal
     #(0.0 0.0 0.0 0.0)
     (linear-algebra:subtract vec1 vec2))
    (assert-float-equal
     #(1.1 2.2 3.3 4.4)
     (linear-algebra:subtract vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(-1.1 -2.2 -3.3 -4.4)
     (linear-algebra:subtract vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(0.0 0.0 0.0 0.0)
     (linear-algebra:subtract vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  ;; Complex
  (let ((vec1 (linear-algebra:column-vector #C(1.1 2.2) #C(3.3 4.4)))
        (vec2 (linear-algebra:column-vector #C(1.1 2.2) #C(3.3 4.4))))
    (assert-float-equal
     #(#C(0.0 0.0) #C(0.0 0.0))
     (linear-algebra:subtract vec1 vec2))
    (assert-float-equal
     #(#C(1.1 2.2) #C(3.3 4.4))
     (linear-algebra:subtract vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(#C(-1.1 -2.2) #C(-3.3 -4.4))
     (linear-algebra:subtract vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(#C(0.0 0.0) #C(0.0 0.0))
     (linear-algebra:subtract vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  (let ((vec1 (linear-algebra:row-vector #C(1.1 2.2) #C(3.3 4.4)))
        (vec2 (linear-algebra:row-vector #C(1.1 2.2) #C(3.3 4.4))))
    (assert-float-equal
     #(#C(0.0 0.0) #C(0.0 0.0))
     (linear-algebra:subtract vec1 vec2))
    (assert-float-equal
     #(#C(1.1 2.2) #C(3.3 4.4))
     (linear-algebra:subtract vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(#C(-1.1 -2.2) #C(-3.3 -4.4))
     (linear-algebra:subtract vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(#C(0.0 0.0) #C(0.0 0.0))
     (linear-algebra:subtract vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  ;; Errors
  (assert-error
   'error
   (linear-algebra:subtract
    (linear-algebra:column-vector 1.1 2.2 3.3 4.4)
    (linear-algebra:row-vector 1.1 2.2 3.3 4.4)))
  (assert-error
   'error
   (linear-algebra:subtract
    (linear-algebra:row-vector 1.1 2.2 3.3 4.4)
    (linear-algebra:column-vector 1.1 2.2 3.3 4.4))))

;;; Destructive vector subtraction

(define-test nsubtract-data-vector
  (:tag :data-vector :nsubtract)
  ;; Real
  (let ((vec1 (linear-algebra:column-vector 1.1 2.2 3.3 4.4))
        (vec2 (linear-algebra:column-vector 1.1 2.2 3.3 4.4)))
    (assert-eq vec1 (linear-algebra:nsubtract vec1 vec2))
    (assert-float-equal
     #(0.0 0.0 0.0 0.0) vec1)
    (assert-float-equal
     #(-2.2 -4.4 -6.6 -8.8)
     (linear-algebra:nsubtract vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(-5.5 -11.0 -16.5 -22.0)
     (linear-algebra:nsubtract vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(-13.2 -26.4 -39.6 -52.8)
     (linear-algebra:nsubtract vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  (let ((vec1 (linear-algebra:row-vector 1.1 2.2 3.3 4.4))
        (vec2 (linear-algebra:row-vector 1.1 2.2 3.3 4.4)))
    (assert-eq vec1 (linear-algebra:nsubtract vec1 vec2))
    (assert-float-equal
     #(0.0 0.0 0.0 0.0) vec1)
    (assert-float-equal
     #(-2.2 -4.4 -6.6 -8.8)
     (linear-algebra:nsubtract vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(-5.5 -11.0 -16.5 -22.0)
     (linear-algebra:nsubtract vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(-13.2 -26.4 -39.6 -52.8)
     (linear-algebra:nsubtract vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  ;; Complex
  (let ((vec1 (linear-algebra:column-vector #C(1.1 2.2) #C(3.3 4.4)))
        (vec2 (linear-algebra:column-vector #C(1.1 2.2) #C(3.3 4.4))))
    (assert-eq vec1 (linear-algebra:nsubtract vec1 vec2))
    (assert-float-equal
     #(#C(0.0 0.0) #C(0.0 0.0)) vec1)
    (assert-float-equal
     #(#C(-2.2 -4.4) #C(-6.6 -8.8))
     (linear-algebra:nsubtract vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(#C(-5.5 -11.0) #C(-16.5 -22.0))
     (linear-algebra:nsubtract vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(#C(-13.2 -26.4) #C(-39.6 -52.8))
     (linear-algebra:nsubtract vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  (let ((vec1 (linear-algebra:row-vector #C(1.1 2.2) #C(3.3 4.4)))
        (vec2 (linear-algebra:row-vector #C(1.1 2.2) #C(3.3 4.4))))
    (assert-eq vec1 (linear-algebra:nsubtract vec1 vec2))
    (assert-float-equal
     #(#C(0.0 0.0) #C(0.0 0.0)) vec1)
    (assert-float-equal
     #(#C(-2.2 -4.4) #C(-6.6 -8.8))
     (linear-algebra:nsubtract vec1 vec2 :scalar2 2.0))
    (assert-float-equal
     #(#C(-5.5 -11.0) #C(-16.5 -22.0))
     (linear-algebra:nsubtract vec1 vec2 :scalar1 2.0))
    (assert-float-equal
     #(#C(-13.2 -26.4) #C(-39.6 -52.8))
     (linear-algebra:nsubtract vec1 vec2 :scalar1 2.0 :scalar2 2.0)))
  ;; Errors
  (assert-error
   'error
   (linear-algebra:nsubtract
    (linear-algebra:column-vector 1.1 2.2 3.3 4.4)
    (linear-algebra:row-vector 1.1 2.2 3.3 4.4)))
  (assert-error
   'error
   (linear-algebra:nsubtract
    (linear-algebra:row-vector 1.1 2.2 3.3 4.4)
    (linear-algebra:column-vector 1.1 2.2 3.3 4.4))))

;;; Vector data dot product

(define-test product-data-vector
  (:tag :data-vector :product)
  ;; Real vectors
  (assert-rational-equal
   55
   (linear-algebra:product
    (linear-algebra:row-vector 1 2 3 4 5)
    (linear-algebra:column-vector 1 2 3 4 5)))
  (assert-float-equal
   55.0
   (linear-algebra:product
    (linear-algebra:row-vector 1.0 2.0 3.0 4.0 5.0)
    (linear-algebra:column-vector 1.0 2.0 3.0 4.0 5.0)))
  (assert-float-equal
   55.0d0
   (linear-algebra:product
    (linear-algebra:row-vector 1.0d0 2.0d0 3.0d0 4.0d0 5.0d0)
    (linear-algebra:column-vector 1.0d0 2.0d0 3.0d0 4.0d0 5.0d0)))
  ;; Real vectors with conjugate keyword
  (assert-rational-equal
   55
   (linear-algebra:product
    (linear-algebra:row-vector 1 2 3 4 5)
    (linear-algebra:column-vector 1 2 3 4 5)))
  ;; Complex vectors
  (assert-rational-equal
   #C(8 18)
   (linear-algebra:product
    (linear-algebra:row-vector #C(1 1) #C(2 1) #C(3 1))
    (linear-algebra:column-vector #C(1 2) #C(2 2) #C(3 2))))
  (assert-float-equal
   #C(8.0 18.0)
   (linear-algebra:product
    (linear-algebra:row-vector #C(1.0 1.0) #C(2.0 1.0) #C(3.0 1.0))
    (linear-algebra:column-vector #C(1.0 2.0) #C(2.0 2.0) #C(3.0 2.0))))
  (assert-float-equal
   #C(8.0d0 18.0d0)
   (linear-algebra:product
    (linear-algebra:row-vector
     #C(1.0d0 1.0d0) #C(2.0d0 1.0d0) #C(3.0d0 1.0d0))
    (linear-algebra:column-vector
     #C(1.0d0 2.0d0) #C(2.0d0 2.0d0) #C(3.0d0 2.0d0))))
  ;; Errors
  (assert-error
   'error
   (linear-algebra:product
    (linear-algebra:row-vector 1 2 3)
    (linear-algebra:column-vector 1 2 3 4)))
  (assert-equal
   14
   (linear-algebra:product
    (linear-algebra:column-vector 1 2 3)
    (linear-algebra:column-vector 1 2 3)))
  (assert-equal
   14
   (linear-algebra:product
    (linear-algebra:row-vector 1 2 3)
    (linear-algebra:row-vector 1 2 3)))
  (assert-equal
   14
   (linear-algebra:product
    (linear-algebra:column-vector 1 2 3)
    (linear-algebra:row-vector 1 2 3))))

;; test the inner product
(define-test inner-product-vector
  (:tag :data-vector :product)
  ;; Real vectors
  (assert-rational-equal
  14
   (linear-algebra-kernel:inner-product-vector
    (linear-algebra:row-vector 1 2 3)
    (linear-algebra:row-vector 1 2 3)
    nil))
)
