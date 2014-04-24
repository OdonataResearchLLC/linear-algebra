;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp -*-
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

(in-package :asdf)

(defsystem :linear-algebra-test
  :description "Unit Tests for Linear Algebra in Common Lisp"
  :version "0.1.0"
  :author "Thomas M. Hermann <thomas.m.hermann@odonata-research.com>"
  :license "MIT"
  :depends-on ("lisp-unit" "linear-algebra")
  :components
  ((:file "linear-algebra-test")
   ;; Linear algebra kernel tests
   (:module kernel
    :depends-on ("linear-algebra-test")
    :components
    ((:file "utility")
     (:file "permute")
     (:file "unary-operations")
     (:file "binary-operations")
     (:file "rotation")))
   ;; Linear algebra interface
   (:module interface
    :depends-on ("linear-algebra-test")
    :components
    ((:file "matrix")
     (:file "identity-matrix" :depends-on ("matrix"))
     (:file "permutation-matrix" :depends-on ("matrix"))))
   ;; Common lisp sequence tests
   (:module sequence
    :depends-on ("linear-algebra-test")
    :components
    ((:file "list")
     (:file "vector")
     (:file "array")))
   ;; Linear algebra tests
   (:file "data-vector" :depends-on ("linear-algebra-test"))
   (:file "dense-matrix" :depends-on ("interface"))
   (:file "square-matrix" :depends-on ("interface"))
   (:file "hermitian-matrix" :depends-on ("interface"))
   (:file "symmetric-matrix" :depends-on ("interface"))))
