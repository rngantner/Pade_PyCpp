#!/usr/bin/env python2
# coding: utf-8
import time
from numpy import *
from numpy.linalg import norm
from scipy.linalg import expm
import matplotlib.pyplot as plt
import sys

# add path to exmat module
sys.path.append('..')
import expmat as e

from numpy.testing import assert_almost_equal
class TestTypes:
    def helper(self,T,prec=10):
        Nvals = [10, 50, 100, 200]
        for n in Nvals:
            A = array(random.random((n,n)), dtype=T)
            A = A/norm(A)
            v = ones((n,1), dtype=T)
            # C++ version
            expA = zeros_like(A)
            e.expmc(A,expA)
            # assertion
            assert_almost_equal(expm(A),expA,prec)

    def test_float32(self):
        self.helper(float32,6)
    
    def test_float64(self):
        self.helper(float64)
    
    def test_complex64(self):
        self.helper(complex64,6)
    
    def test_complex128(self):
        self.helper(complex128)


