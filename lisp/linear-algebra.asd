;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp -*-
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

(in-package :asdf)

(defsystem :linear-algebra
  :description "Linear Algebra in Common Lisp."
  :version "0.1.0"
  :author "Thomas M. Hermann <thomas.m.hermann@odonata-research.com>"
  :license "MIT"
  :depends-on ("floating-point")
  :components
  ((:file "linear-algebra" :depends-on ("kernel"))
   ;; Linear algebra kernel functions
   (:module kernel
    :components
    ((:file "linear-algebra-kernel")
     (:file "utility" :depends-on ("linear-algebra-kernel"))
     (:file "permute" :depends-on ("linear-algebra-kernel"))
     (:file "unary-operations"
      :depends-on ("linear-algebra-kernel"))
     (:file "binary-operations"
      :depends-on ("linear-algebra-kernel"))
     (:file "rotation" :depends-on ("unary-operations"))
     (:file "gauss" :depends-on ("linear-algebra-kernel"))
     (:file "cholesky" :depends-on ("unary-operations"))))
   ;; Interface
   (:module interface
    :depends-on ("linear-algebra" "kernel")
    :components
    ((:file "fundamental-ops")
     (:file "vector" :depends-on ("fundamental-ops"))
     (:file "matrix" :depends-on ("fundamental-ops"))
     (:file "identity-matrix" :depends-on ("matrix"))
     (:file "permutation-matrix" :depends-on ("matrix"))))
   ;; Common Lisp sequences
   (:module sequence
    :depends-on ("interface")
    :components
    ((:file "list")
     (:file "vector")
     (:file "array")))
   ;; Linear algebra classes and operations
   (:file "data-vector" :depends-on ("interface"))
   (:file "dense-matrix" :depends-on ("data-vector"))
   (:file "square-matrix" :depends-on ("dense-matrix"))
   (:file "hermitian-matrix" :depends-on ("square-matrix"))
   (:file "symmetric-matrix" :depends-on ("square-matrix"))))
