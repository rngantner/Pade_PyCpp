from numpy.linalg import norm, solve
from numpy import asarray, eye, dot, ceil, log2

def expm(A):
    """Compute the matrix exponential using Pade approximation.
    
    Parameters
    ----------
    A : array, shape(M,M)
        Matrix to be exponentiated

    Returns
    -------
    expA : array, shape(M,M)
        Matrix exponential of A
    
    See: "The Scaling and Squaring Method for the Matrix
    Exponential Revisited", Nicholas Higham, 2005.
    
    """
    A = asarray(A)
    A_L1 = norm(A,1)
    n_squarings = 0
    
    if A.dtype == 'float64' or A.dtype == 'complex128':
        if A_L1 < 1.495585217958292e-002:
            U,V = _pade3(A)
        elif A_L1 < 2.539398330063230e-001:
            U,V = _pade5(A)
        elif A_L1 < 9.504178996162932e-001:
            U,V = _pade7(A)
        elif A_L1 < 2.097847961257068e+000:
            U,V = _pade9(A)
        else:
            maxnorm = 5.371920351148152
            n_squarings = max(0, int(ceil(log2(A_L1 / maxnorm))))
            A = A / 2**n_squarings
            U,V = _pade13(A)
    elif A.dtype == 'float32' or A.dtype == 'complex64':
        if A_L1 < 4.258730016922831e-001:
            U,V = _pade3(A)
        elif A_L1 < 1.880152677804762e+000:
            U,V = _pade5(A)
        else:
            maxnorm = 3.925724783138660
            n_squarings = max(0, int(ceil(log2(A_L1 / maxnorm))))
            A = A / 2**n_squarings
            U,V = _pade7(A)
    else:
        print "invalid type!",A.dtype
        return
    
    P = U + V  # p_m(A) : numerator
    Q = -U + V # q_m(A) : denominator
    R = solve(Q,P)
    # squaring step to undo scaling
    for i in range(n_squarings):
        R = dot(R,R)
    
    return R

# implementation of Pade approximations of various degree using the
# algorithm presented in [Higham 2005]

def _pade3(A):
    b = (120., 60., 12., 1.)
    ident = eye(*A.shape, dtype=A.dtype)
    A2 = dot(A,A)
    U = dot(A , (b[3]*A2 + b[1]*ident))
    V = b[2]*A2 + b[0]*ident
    return U,V

def _pade5(A):
    b = (30240., 15120., 3360., 420., 30., 1.)
    ident = eye(*A.shape, dtype=A.dtype)
    A2 = dot(A,A)
    A4 = dot(A2,A2)
    U = dot(A, b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

def _pade7(A):
    b = (17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1.)
    ident = eye(*A.shape, dtype=A.dtype)
    A2 = dot(A,A)
    A4 = dot(A2,A2)
    A6 = dot(A4,A2)
    U = dot(A, b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

def _pade9(A):
    b = (17643225600., 8821612800., 2075673600., 302702400., 30270240.,
                2162160., 110880., 3960., 90., 1.)
    ident = eye(*A.shape, dtype=A.dtype)
    A2 = dot(A,A)
    A4 = dot(A2,A2)
    A6 = dot(A4,A2)
    A8 = dot(A6,A2)
    U = dot(A, b[9]*A8 + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[8]*A8 + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

def _pade13(A):
    b = (64764752532480000., 32382376266240000., 7771770303897600.,
    1187353796428800., 129060195264000., 10559470521600., 670442572800.,
    33522128640., 1323241920., 40840800., 960960., 16380., 182., 1.)
    ident = eye(*A.shape, dtype=A.dtype)
    A2 = dot(A,A)
    A4 = dot(A2,A2)
    A6 = dot(A4,A2)
    U = dot(A,dot(A6, b[13]*A6 + b[11]*A4 + b[9]*A2) + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = dot(A6, b[12]*A6 + b[10]*A4 + b[8]*A2) + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

if __name__ == '__main__':
    import time
    from numpy import random, inf
    A = random.random((200,200))
    ttot = inf
    for i in range(10):
        t = time.time()
        B = expm(A)
        ttot = min(ttot,time.time()-t)
    print "time used:",ttot

