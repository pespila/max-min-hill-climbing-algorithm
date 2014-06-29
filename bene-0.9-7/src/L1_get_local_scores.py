#!/usr/bin/python
import numpy as np
import sys, math, struct
import coliche

x    = None
resf = None
sf   = None

def bic_score(y,A):
    k = len(A[0])
    n = len(A)
    residuals = np.linalg.lstsq(A, y)[1]
    return -residuals[0] - 0.5*k*math.log(n)

def aic_score(y,A):
    residuals = np.linalg.lstsq(A, y)[1]
    k = len(A[0])
    n = len(A)
    return -residuals[0] - 2*len(A[0])


scorefun = {'BIC': bic_score,
            'AIC': aic_score}

def scores(vs):
    for i in vs[1:]:
        indices = vs[:]
        indices.remove(i)
        A = np.take(x,indices,axis=1)
        y = np.take(x,[i],axis=1)
        residuals = np.linalg.lstsq(A, y)[1]
        resf.write(struct.pack('d', sf(y,A)))

def walk_varsets(vs, nof_calls):
    scores(vs);
    if len(vs) > 1:
        for i in xrange(1,nof_calls):
            nvs = vs[:]
            nvs.remove(i)
            walk_varsets(nvs, i);

def main(datfile, score, resfile):
    global x,sf,resf

    sf = scorefun[score]
    dataf = open(datfile)
    x = np.array([map(float,l.split())+[1.0]
                  for l in dataf])

    resf  = open(resfile,"wb")
    nof_vars = len(x[0]) - 1                   # -1 for constant
    walk_varsets(range(nof_vars+1),nof_vars+1) # +1 for constant
    resf.close()
    print nof_vars

coliche.che(main,"datfile; score; resfile")
