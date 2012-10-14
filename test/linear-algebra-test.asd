;;;; -*- Mode: Lisp; Syntax: ANSI-Common-Lisp -*-
#|

 Linear Algebra in Common Lisp Unit Tests

 Copyright (c) 2010-2012, Thomas M. Hermann
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

(in-package :asdf)

(defsystem :linear-algebra-test
  :description "Unit Tests for Linear Algebra in Common Lisp"
  :version "0.1.0"
  :author "Thomas M. Hermann <thomas.m.hermann@odonata-research.com>"
  :license "BSD"
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
     (:file "binary-operations")))
   ;; Linear algebra interface
   (:module interface
    :depends-on ("linear-algebra-test")
    :components
    ((:file "vector")
     (:file "matrix")
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
