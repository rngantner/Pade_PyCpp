#!/usr/bin/env python2
# coding: utf-8
import time
from numpy import *
from numpy.linalg import norm
from scipy.linalg import expm
import matplotlib.pyplot as plt
import os

import sys
sys.path.append('..')
import expmat as e
import expm as e2

if __name__ == '__main__':
    if not os.path.exists('results'):
        os.mkdir('results')
    number_of_Nvals = 20
    Creps = 3
    Preps = 1
    P2reps= 1
    Nmin = 10
    lNmin = log(Nmin)/log(10)
    #Nmax = 20000
    Nmax = 300
    compiler = 'gcc'
    lNmax = log(Nmax)/log(10)
    Nvals = array(logspace(lNmin, lNmax, number_of_Nvals),dtype=int)
    
    type_to_name = {float32:'float32', float64:'float64', complex64:'complex64', complex128:'complex128'}
    for T in [float32, float64, complex64, complex128]:
        thistime = time.time()
        print "current type:", type_to_name[T]
        Tc = []
        Tp = []
        Tp2= []
        for n in Nvals:
            print "n:",n,
            A = array(random.random((n,n)), dtype=T)
            ctime = inf; ptime = inf; p2time = inf
            # time C++ version
            for r in range(Creps):
                expA = zeros_like(A)
                
                t = time.time()
                e.expmc(A,expA) # call C version with Eigen Pade approx.
                ctime = min(ctime, time.time()-t)
            
            # time scipy version
            for r in range(Preps):
                t = time.time()
                expA2 = expm(A)
                ptime = min(ptime, time.time()-t)
            
            # time new python version
            for r in range(P2reps):
                t = time.time()
                expA3 = e2.expm(A)
                p2time = min(p2time, time.time()-t)
            
            Tc.append(ctime)
            Tp.append(ptime)
            Tp2.append(p2time)
            speedup = ptime/p2time
            cpeedup = ptime/ctime
            print "C (Eigen): %1.4f"%ctime,"  scipy: %1.4f"%ptime,"  new python: %1.4f"%p2time,"  Python speedup: %1.4f, %1.1f%%"%(speedup,100*(speedup-1)/speedup),"  C++ speedup: %1.4f, %1.1f%%"%(cpeedup,100*(cpeedup-1)/cpeedup)
            
            normA = norm(expA-expA2)/norm(expA) # "relative" norm
            normA2 =norm(expA3-expA2)/norm(expA3) # "relative" norm
            epsilon = 1e-5
            if normA > epsilon or normA2 > epsilon:
                print "A-A2 or A3-A2 is not small (%f or %f > %f)!"%(normA,normA2,epsilon)
            del A, expA, expA2, expA3
        
        savez('results/'+compiler+'_'+type_to_name[T],Tp=Tp,Tp2=Tp2,Tc=Tc)
        mn = mean(array(Tp)/Tp2)
        print "mean speedup: %1.4f %1.4f"%(mn,(mn-1)/mn)
        plt.loglog(Nvals, Tp2,'-s', label="Python new")
        plt.loglog(Nvals, Tp, '-o', label="Python scipy")
        plt.loglog(Nvals, Tc, '-^', label="C++")
        plt.legend(loc='best')
        plt.xlabel('N', fontsize=16)
        plt.ylabel('Time', fontsize=16)
        plt.title(u'Timing of Padé Approximation. Type: %s'%type_to_name[T])
        ax = plt.axis()
        plt.axis((Nmin,Nmax,ax[2],ax[3]))
        plt.draw()
        plt.savefig('results/expmtimings_'+type_to_name[T]+'.pdf')
        #plt.show()
        plt.close()
        print "time for %s: %1.3f sec"%(T,time.time()-thistime)

