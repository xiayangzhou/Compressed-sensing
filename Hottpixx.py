# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 15:56:42 2016

@author: xiayangzhou
"""
import numpy as np
from cvxopt import matrix, solvers, spmatrix, sparse
import random

def permutation (X,i):# return a column Z verified the hypothesis z_2>...z_f
    #X is a vector, i-th column of the Matrix C
    #
    C_ii=X[i]
    a=np.max(X)+1
    X[i]=a
    index=np.argsort(X)[::-1]
    Z=np.sort(X)[::-1]
    Z[0]=C_ii
    return Z, index, C_ii
    

def proj_01(x): #projection on interval [0,1]
    if x>1 :return 1
    if x<0 :return 0
    else :  return x
    
    
def proj_on_phi_0 (Z,index):#algo 5 
    f=Z.shape[0] 
    x=np.zeros(f)
    xx = np.zeros(f)
    mu = Z[0]
    for k in range(1,f): 
        if Z[k] <= proj_01(mu):
            k_c=k
            break
        else:
            mu=(k-1)/k*mu+1./k*Z[k]
        k_c=k    
    print 'k_c=',k_c,'mu=',mu
    x[0]=proj_01(mu)
    for k in range(1,k_c):
        x[k]=proj_01(mu)
    for k in range(k_c,f): 
        x[k]=max(Z[k],0)
    #print 'x=',x    
    #to return X, after permutation inverse
    for i in range(f):
        xx[index[i]]=x[i]
    return xx

def proj_on_phi_final(C):
    f = C.shape[0]
    for i in range(f):
        X = np.squeeze(np.asarray(C[:,i]))
        print X.shape
        Z, index, C_ii=permutation(X,i)
        #there always has a problem of compatibility between matrix, array of i-th column of C
        C[:,i]=np.matrix(proj_on_phi_0(Z,index)).T
        #C[:,i]=proj_on_phi_0(Z,index)

    return C
        

def Norm_inf1 (X,W):
    # X, W are two matrix of shape respectively f*n et r*n
    f,n = X.size
    r = W.size[0]
    
    print f,n,r
    
    P = matrix(W).trans()
    F = matrix(1.0,(f,r))
    
    onesn = matrix(1.0,(n,1))
    Idn = spmatrix(1.0, range(n),range(n))
    Idr = spmatrix(1.0, range(r),range(r))
    Zrn = spmatrix(0,[r-1],[n-1])
    Zr = spmatrix(0,[r-1],[0])
    #Zn = spmatrix(0,[n-1],[0])
    A = sparse([ [P,-P,-Idr], [-Idn,-Idn,Zrn] ]) 
    
    for i in range(f):
        V =  X[i,:].trans()
        C = matrix([[V,-V,Zr]])
        e = matrix([ [Zr, onesn] ])
        solution = solvers.lp(e,A,C)['x']
        F[i,:] = solution[range(r)].trans()
    return F    
    

def Hottpixx(X,s_p,s_d,N):
    f,n = X.shape
    p = np.arange(f)[::-1]
    C = np.zeros((f,f))
    beta = 0
    GDindex=np.arange(n)
    mu = np.zeros(f)
    for j in range(f):#μj is the number of nonzeros in row j divided by n.
        mu[j] = (np.sum(X[j,:] != 0)+0.0)/n
    for t in range(N):
        for i in range(n*10):
            k=random.choice(GDindex)
            C = C + s_p * np.matrix(np.sign(X[:,k]-np.dot(C,X[:,k]))).T.dot(np.matrix(X[:,k])) - s_d*np.diag(mu*(beta-p))  
            
        #print 'C.shape=',C.shape
        print 'C=',C
        proj_on_phi_final(C)
        
        beta = beta +s_d*(np.trace(C)-r)#we assume r is known, r is a global variable
        print t,'-th iteration is done'
    print C
    diag=np.diag(C)
    print 'diag=',diag
    index=np.argsort(diag)[::-1]
    I=index[:r]
    #print I
    W = X[I]
    
    #convert X,W into matrix type defined in cvxopt in order to apply the argmin function
    W1=matrix(W)
    X1=matrix(X)
    F = Norm_inf1(X1,W1)
    #convert again F into numpy matrix
    F1=np.array(F)
    return F1,W

