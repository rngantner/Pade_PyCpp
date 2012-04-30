#-*- encoding: utf8 -*-

import matplotlib.pyplot as plt
from numpy import *
import os

type_to_name = {float32:'float32', float64:'float64', complex64:'complex64', complex128:'complex128'}
type_to_code = {float32:'f32', float64:'f64', complex64:'c64', complex128:'c128'}
types = [float32, float64, complex64, complex128]
fmt = ('-o','-^','-s','-<','->','-.','-x','-+')

plotfun = plt.loglog
basedir = 'results'

Ncores = (1,2,4)
#Ncores = (1,2,4,8,16)
compiler = 'gcc'

speedup = {}; onecoredata = {}
speedup_full = {}
for T in types:
    speedup[T] = 0
    speedup_full[T] = {}

# compute speedup; take minimum of old time and new time
for i in Ncores: # loop over number of cores
    fname = os.path.join(basedir,str(i)+'threads_'+compiler+'.npz')
    if not os.path.exists(fname):
        continue
    f = load(fname)
    Nvals = f['Nvals']
    if i == 1:
        for T in types: onecoredata[T] = f['T'+type_to_code[T]]
    j=1
    for T in types: # loop over types
        plt.figure(j)
        data = f['T'+type_to_code[T]]
        plotfun(Nvals,data,fmt[mod(i-1,len(fmt))],label=str(i)+' cores')
        # maximal speedup measurement
        speedup[T] = max(speedup[T],max(onecoredata[T]/data))
        speedup_full[T][i] = onecoredata[T]/data
        j = j + 1

for T in types: print 'maximal speedup (%s): %1.2f'%(type_to_name[T], speedup[T])

i = 1
for T in types:
    # time vs N plot
    plt.figure(i)
    plt.legend(loc='best')
    plt.xlabel('N', fontsize=16)
    plt.ylabel('Time', fontsize=16)
    plt.title(u'Padé Approximation Timing. %s, max. speedup: %1.1f'%(type_to_name[T],speedup[T]))
    
    ax = plt.axis()
    plt.axis((min(Nvals),max(Nvals),ax[2],ax[3]))
    plt.draw()
    plt.savefig(os.path.join(basedir,'expm_ncorestimings_'+type_to_name[T]+'.pdf'))
    #plt.show()
    plt.close()
    
    i = i + 1

for T in types:
    plotfun = plt.semilogx
    # speedup vs N plot
    plt.figure()
    for i in Ncores:
        plotfun(Nvals,speedup_full[T][i],fmt[mod(i-1,len(fmt))],label=str(i)+' cores')
    plt.legend(loc='best')
    plt.xlabel('N', fontsize=16)
    plt.ylabel('Speedup',fontsize=16)
    plt.title(u'Padé Approximation Speedup as Function of $N$. %s'%type_to_name[T])
    ax = plt.axis()
    plt.axis((min(Nvals),max(Nvals),ax[2],ax[3]))
    plt.draw()
    plt.savefig(os.path.join(basedir,'expm_ncoresspeedup_'+type_to_name[T]+'.pdf'))
    #plt.show()
    plt.close()

plotfun = plt.plot
plt.figure()
j=0
for T in types:
    xvals = array(Ncores)
    yvals = []
    for i in xvals:
        yvals.append(speedup_full[T][i][-1])
    plotfun(xvals,yvals,fmt[mod(j,len(fmt))],label=type_to_name[T],markersize=10)
    j = j + 1

plt.legend(loc='best')
plt.xlabel('Cores',fontsize=16)
plt.ylabel('Speedup',fontsize=16)
plt.title(r'Measured Speedup. $N=2000$',fontsize=16)
plt.draw()
plt.savefig(os.path.join(basedir,'expm_speedup_2000.pdf'))
plt.close()

