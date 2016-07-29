# -*- coding: utf-8 -*-

from cvxopt import matrix, solvers, spmatrix, sparse

#============================================================
# distFromCon: L1 distance between a vector and a convex hull
#------------------------------------------------------------
# params : 
#    V vector 
#    P matrix which rows define the convex hull
#============================================================

def distFromConv(V, P):
    # Make sure that the parameters are supported by cvxopt
    V = matrix(V)
    P = matrix(P).trans()
    # Dimension of the vector space
    m = len(V)
    # Number of lements defining the hull
    n = P.size[1]
    # Vectors vith only ones as entry
    onesn = matrix(1.0,(n,1))
    onesm = matrix(1.0,(m,1))
    # Identity matrices
    Idm = spmatrix(1.0, range(m),range(m))
    Idn = spmatrix(1.0, range(n),range(n))
    # Zero matrix
    Znm = spmatrix(0,[n-1],[m-1])
    Zn = spmatrix(0,[n-1],[0])
    Zm = spmatrix(0,[m-1],[0])
    # Expressing the problem as an LP problem
    # Inequality constraints
    A = sparse([ [P,-P,-Idn], [-Idm,-Idm,Znm] ]) 
    C = matrix([[V,-V,Zn]])
    e = matrix([ [Zn, onesm] ])
    # Equality constraint
    G = matrix([onesn,Zm]).trans()
    # Computing the solution
    solution = solvers.lp(e,A,C,G,matrix(1.0))['x']
    # Computing the L1-closest point in the hull
    u = solution[range(n)]
    V_hull = P*u
    return sum(abs(V-V_hull))

#============================================================
# AGKM: Approximably Separable non-negative matrix 
# factorization
#------------------------------------------------------------
# params : 
#    X matrix to be factorized 
#    alpha simplicial robustness if hott topics (see paper)
#    epsilon  approximation control parameter (see paper)
#============================================================

def AGKM(X, alpha, epsilon):
    # Make sure that X is supported by cvxopt library
    X = matrix(X)
    # Number of features
    nrows = X.size[0]
    # Pairwise L1 distances between rows 
    D = matrix(0.0, (nrows,nrows))
    for i in range(nrows):
        for j in range(i+1,nrows):
            D[i,j] = sum(abs(X[i,:]-X[j,:]))
            D[j,i] = D[i,j]
    # Finding the set R of rows that are simplicial
    R = []
    # Computing minimal distance between elements of R
    minDist = 5*epsilon/alpha+2*epsilon
    for k in range(nrows):    
        # Set of indices of rows far enough
        Nk = []
        for j in range(nrows):
            if (D[k,j] >= minDist and k!=j): Nk = Nk+[j]
        # Matrix which rows are the ones indexed by Nk
        P = X[Nk,:]
        # Computing the L1 distance between Xk and the convex 
        # hull defined by the rows of Nk
        deltak = distFromConv(X[k,:].trans(),P)
        # Add the k-th row to R if deltak is large enough
        if deltak > 2*epsilon: R = R+[k]
    # Cluster the elements of R ot obtain the rows of 
    # W (due to the particular geometry of the problem,
    # the clustering is done in O(n), see paper)
    # L1 diameter of the clusters
    diam = 10*epsilon/alpha+epsilon
    clusterR = []
    for i in R:
        if filter(lambda x: x<= diam, D[i,clusterR]) == []:
            clusterR = clusterR+[i]
    # W Matrix of canonical rows
    W = X[clusterR,:]
    
    F = Norm_inf1(X,W)
    return (F,W)
        
        
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
    
            
