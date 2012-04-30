#!/usr/bin/env python2
# coding: utf-8
import time
from numpy import *
import os
import os.path
import sys

sys.path.append('..')
import expmat as e

def runtiming(T,Nvals,Tc,reps):
    """runs timing with given type T"""
    print "current type:", type_to_name[T]
    Tc[T] = []
    for n in Nvals:
        print "  n:",n,
        A = array(random.random((n,n)), dtype=T)
        ctime = inf
        # time C++ version
        for r in range(reps):
            expA = zeros_like(A)
            
            t = time.time()
            e.expmc(A,expA) # call C version with Eigen Pade approx.
            ctime = min(ctime, time.time()-t)
        print ctime
        Tc[T].append(ctime)
        del A, expA

if __name__ == '__main__':
    number_of_Nvals = 20
    Creps = 5
    Nmin = 10
    lNmin = log(Nmin)/log(10)
    #Nmax = 20000
    Nmax = 2000
    compiler = 'gcc'
    lNmax = log(Nmax)/log(10)
    Nvals = array(logspace(lNmin, lNmax, number_of_Nvals),dtype=int)
    
    type_to_name = {float32:'float32', float64:'float64', complex64:'complex64', complex128:'complex128'}
    types = [float32, float64, complex64, complex128]
    Tc = {}
    for T in types:
        # writes into Tc
        runtiming(T,Nvals,Tc,Creps)
    
    if not os.path.exists('results'):
        os.mkdir('results')
    outfilename = os.path.join('results',os.environ['OMP_NUM_THREADS']+'threads_'+compiler+'.npz')
    if os.path.exists(outfilename):
        print "overwriting",outfilename,"with minimum of current contents and new computed times"
        # load file, take minimum of already computed value and new value
        f = load(outfilename)
        Tc[float32] = minimum(Tc[float32],f['Tf32']) # component-wise minimum!
        Tc[float64] = minimum(Tc[float64],f['Tf64']) # component-wise minimum!
        Tc[complex64] = minimum(Tc[complex64],f['Tc64']) # component-wise minimum!
        Tc[complex128] = minimum(Tc[complex128],f['Tc128']) # component-wise minimum!
        f.close()
    savez(outfilename,Nvals=Nvals,Tf32=Tc[float32],Tf64=Tc[float64],Tc64=Tc[complex64],Tc128=Tc[complex128])

