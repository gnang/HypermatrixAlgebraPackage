#*************************************************************************#
#       Copyright (C) 2014 Edinah K. Gnang <kgnang@gmail.com>,            #
#                          Ori Parzanchevski                              #
#                          Yuval Filmus                                   #
#                                                                         #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                                                                         #
#    This code is distributed in the hope that it will be useful,         #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    #
#    General Public License for more details.                             #
#                                                                         #
#  The full text of the GPL is available at:                              #
#                                                                         #
#                  http://www.gnu.org/licenses/                           #
#*************************************************************************#

# Definition of the functions 
def MatrixGenerate(nr, nc, c):
    """
    Generates a list of lists associated with a symbolic nr x nc
    matrix using the input character c followed by indices.

    EXAMPLES:
    ::
        sage: M = MatrixGenerate(2, 2, 'm'); M
        [[m00, m01], [m10, m11]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_q_rows = nr
    n_q_cols = nc

    # Test for dimension match
    if n_q_rows > 0 and n_q_cols > 0:
        # Initialization of the hypermatrix
        q = []
        for i in range(n_q_rows):
            q.append([])
        for i in range(len(q)):
            for j in range(n_q_cols):
                # Filling up the matrix
                (q[i]).append(var(c+str(i)+str(j)))
        return q

    else :
        raise ValueError, "Input dimensions "+\
str(nr)+" and "+str(nc)+" must both be non-zero positive integers."

def SymMatrixGenerate(nr, c):
    """
    Generates a list of lists associated with a symbolic nr x nc
    symmetric matrix using the input character c followed by
    indices.

    EXAMPLES:
    ::
        sage: M = SymMatrixGenerate(2, 'm'); M
        [[m00, m01], [m10, m11]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_q_rows = nr
    n_q_cols = nr

    # Test for dimension match
    if n_q_rows > 0 and n_q_cols > 0:
        # Initialization of the hypermatrix
        q = []
        for i in range(n_q_rows):
            q.append([])
        for i in range(len(q)):
            for j in range(n_q_cols):
                # Filling up the matrix
                (q[i]).append(var(c+str(min(i,j))+str(max(i,j))))
        return q

    else :
        raise ValueError, "Input dimensions "+\
str(nr)+" must be a non-zero positive integers."

def HypermatrixGenerate(*args):
    """
    Generates a list of lists associated with a symbolic arbitrary
    order hypematrix of size specified by the input using and the
    entries are determined by the last input character

    EXAMPLES:
    ::
        sage: M = HypermatrixGenerate(2, 2, 2, 'm'); M


    AUTHORS:
    - Edinah K. Gnang, Ori Parzanchevski and Yuval Filmus
    """
    if len(args) == 1:
        return var(args[0])
    return [apply(HypermatrixGenerate, args[1:-1] + (args[-1] + str(i),)) for i in range(args[0])]

def HypermatrixGenerateAllOne(*args):
    """
    Generates a list of lists associated with a symbolic arbitrary
    order hypematrix of size specified by the input using and the
    entries are determined by the last input character

    EXAMPLES:
    ::
        sage: M = HypermatrixGenerateAllOne(2, 2, 2); M


    AUTHORS:
    - Edinah K. Gnang, Ori Parzanchevski and Yuval Filmus
    """
    if len(args) == 1:
        return [1 for i in range(args[0])]
    return [apply(HypermatrixGenerateAllOne, args[1:] ) for i in range(args[0])]

def HypermatrixGenerateAllZero(*args):
    """
    Generates a list of lists associated with a symbolic arbitrary
    order hypematrix of size specified by the input using and the
    entries are determined by the last input character

    EXAMPLES:
    ::
        sage: M = HypermatrixGenerateAllZero(2, 2, 2); M


    AUTHORS:
    - Edinah K. Gnang, Ori Parzanchevski and Yuval Filmus
    """
    if len(args) == 1:
        return [0 for i in range(args[0])]
    return [apply(HypermatrixGenerateAllZero, args[1:] ) for i in range(args[0])]

def SymHypermatrixGenerate(nr, c):
    """
    Generates a list of lists associated with a symbolic nr x nc x nd
    third order hypematrix using the input character c followed by
    indices.

    EXAMPLES:
    ::
        sage: M = SymHypermatrixGenerate(2, 'm'); M


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_q_rows = nr
    n_q_cols = nr
    n_q_dpts = nr

    # Test for dimension match
    if n_q_rows > 0 and n_q_cols > 0 and n_q_dpts >0:
        # Initialization of the hypermatrix
        q = []
        for i in range(n_q_rows):
            q.append([])
        for i in range(len(q)):
            for j in range(n_q_cols):
                (q[i]).append([])
        for i in range(len(q)):
            for j in range(len(q[i])):
                for k in range(n_q_dpts):
                    if i==j or i==k or j==k:
                        (q[i][j]).append(\
var(c+str(min(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))+str(max(i,j,k))))
                    else:
                        if i == min(i,j,k) and k == max(i,j,k):
                            (q[i][j]).append(\
var(c+str(min(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))+str(max(i,j,k))))
                        elif k == min(i,j,k) and j == max(i,j,k):
                            (q[i][j]).append(\
var(c+str(min(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))+str(max(i,j,k))))
                        elif i == max(i,j,k) and j == min(i,j,k):
                            (q[i][j]).append(\
var(c+str(min(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))+str(max(i,j,k))))
                        else:
                            (q[i][j]).append(\
var(c+str(i+j+k-min(i,j,k)-max(i,j,k))+str(min(i,j,k))+str(max(i,j,k))))
        return q

    else :
        raise ValueError, "Input dimensions "+\
str(nr)+" must be a non-zero positive integer."

def HypermatrixVectorize(A):
    """
    Outputs our canonical vectorization a list
    the input hypermatrices A.

    EXAMPLES:
    ::
        sage: M = HypermatrixVectorize(A); M


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_q_rows = len(A)
    n_q_cols = len(A[0])
    n_q_dpts = len(A[0][0])

    # Test for dimension match
    if n_q_rows>0 and n_q_cols>0 and n_q_dpts>0:
        # Initialization of the hypermatrix
        q = []
        for i in range(n_q_rows):
            for j in range(n_q_cols):
                for k in range(n_q_dpts):
                    q.append(A[i][j][k])
        return q

    else :
        raise ValueError, "The Dimensions non zero."

def HypermatrixAdd(A, B):
    """
    Outputs a list of lists associated with the addtion of
    the two input hypermatrices A and B

    EXAMPLES:
    ::
        sage: M = HypermatrixAdd(A, B); M


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_q_rows = len(B)
    n_q_cols = len(B[0])
    n_q_dpts = len(B[0][0])

    # Test for dimension match
    if n_q_rows==len(A) and n_q_cols==len(A[0]) and n_q_dpts==len(A[0][0]):
        # Initialization of the hypermatrix
        q = []
        for i in range(n_q_rows):
            q.append([])
        for i in range(len(q)):
            for j in range(n_q_cols):
                (q[i]).append([])
        for i in range(len(q)):
            for j in range(len(q[i])):
                for k in range(n_q_dpts):
                    (q[i][j]).append(A[i][j][k]+B[i][j][k])
        return q

    else :
        raise ValueError, "The Dimensions of the input hypermatrices must match."

def HypermatrixHadamardProduct(A, B):
    """
    Outputs a list of lists associated with the addtion of
    the two input hypermatrices A and B

    EXAMPLES:
    ::
        sage: M = HypermatrixHadamardProduct(A, B); M


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_q_rows = len(A)
    n_q_cols = len(A[0])
    n_q_dpts = len(A[0][0])

    # Test for dimension match
    if n_q_rows==len(A) and n_q_cols==len(A[0]) and n_q_dpts==len(A[0][0]):
        # Initialization of the hypermatrix
        q = []
        for i in range(n_q_rows):
            q.append([])
        for i in range(len(q)):
            for j in range(n_q_cols):
                (q[i]).append([])
        for i in range(len(q)):
            for j in range(len(q[i])):
                for k in range(n_q_dpts):
                    (q[i][j]).append(A[i][j][k]*B[i][j][k])
        return q

    else :
        raise ValueError, "The Dimensions of the input hypermatrices must match."

def HypermatrixScale(A, s):
    """
    Outputs a list of lists associated with product of the
    scalar s with the hypermatrix A.

    EXAMPLES:
    ::
        sage: M = HypermatrixScale(A, 3); M


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_q_rows = len(A)
    n_q_cols = len(A[0])
    n_q_dpts = len(A[0][0])

    # Initialization of the hypermatrix
    q = []
    for i in range(n_q_rows):
        q.append([])
    for i in range(len(q)):
        for j in range(n_q_cols):
            (q[i]).append([])
    for i in range(len(q)):
        for j in range(len(q[i])):
            for k in range(n_q_dpts):
                (q[i][j]).append(A[i][j][k]*s)
    return q

def HypermatrixEntryExponent(A, s):
    """
    Outputs a list of lists associated with product of the
    scalar s with the hypermatrix A.

    EXAMPLES:
    ::
        sage: M = HypermatrixEntryExponent(A, 3); M


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_q_rows = len(A)
    n_q_cols = len(A[0])
    n_q_dpts = len(A[0][0])

    # Initialization of the hypermatrix
    q = []
    for i in range(n_q_rows):
        q.append([])
    for i in range(len(q)):
        for j in range(n_q_cols):
            (q[i]).append([])
    for i in range(len(q)):
        for j in range(len(q[i])):
            for k in range(n_q_dpts):
                (q[i][j]).append(A[i][j][k]^s)
    return q

def HypermatrixEntryExponentB(s, A):
    """
    Outputs a list of lists associated with product of the
    scalar s with the hypermatrix A.

    EXAMPLES:
    ::
        sage: M = HypermatrixEntryExponentB(3,A); M


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_q_rows = len(A)
    n_q_cols = len(A[0])
    n_q_dpts = len(A[0][0])

    # Initialization of the hypermatrix
    q = []
    for i in range(n_q_rows):
        q.append([])
    for i in range(len(q)):
        for j in range(n_q_cols):
            (q[i]).append([])
    for i in range(len(q)):
        for j in range(len(q[i])):
            for k in range(n_q_dpts):
                (q[i][j]).append(s^(A[i][j][k]))
    return q

def HypermatrixProduct(A, B, C):
    """
    Outputs a list of lists associated with the ternary
    product of the input hypermatrices A, B and C.

    EXAMPLES:
    ::
        sage: M = HypermatrixProduct(A, B, C); M


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_a_rows = len(A)
    n_a_cols = len(A[0])
    n_a_dpts = len(A[0][0])

    n_b_rows = len(B)
    n_b_cols = len(B[0])
    n_b_dpts = len(B[0][0])

    n_c_rows = len(C)
    n_c_cols = len(C[0])
    n_c_dpts = len(C[0][0])

    # Test for dimension match
    if n_a_rows==n_b_rows and n_b_cols==n_c_cols and n_c_dpts==n_a_dpts and \
n_a_cols==n_b_dpts and n_b_dpts==n_c_rows:
        # Initialization of the hypermatrix
        q = []
        for i in range(n_a_rows):
            q.append([])
        for i in range(len(q)):
            for j in range(n_b_cols):
                (q[i]).append([])
        for i in range(len(q)):
            for j in range(len(q[i])):
                for k in range(n_c_dpts):
                    (q[i][j]).append(\
sum([A[i][l][k]*B[i][j][l]*C[l][j][k] for l in range(n_a_cols)]))
        return q

    else :
        raise ValueError, "Hypermatrix dimension mismatch."

def HypermatrixLogProduct(A, B, C):
    """
    Outputs a list of lists associated with the ternary
    product of the input hypermatrices A, B and C.

    EXAMPLES:
    ::
        sage: M = HypermatrixProduct(A, B, C); M


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_a_rows = len(A)
    n_a_cols = len(A[0])
    n_a_dpts = len(A[0][0])

    n_b_rows = len(B)
    n_b_cols = len(B[0])
    n_b_dpts = len(B[0][0])

    n_c_rows = len(C)
    n_c_cols = len(C[0])
    n_c_dpts = len(C[0][0])

    # Test for dimension match
    if n_a_rows==n_b_rows and n_b_cols==n_c_cols and n_c_dpts==n_a_dpts and \
n_a_cols==n_b_dpts and n_b_dpts==n_c_rows:
        # Initialization of the hypermatrix
        q = []
        for i in range(n_a_rows):
            q.append([])
        for i in range(len(q)):
            for j in range(n_b_cols):
                (q[i]).append([])
        for i in range(len(q)):
            for j in range(len(q[i])):
                for k in range(n_c_dpts):
                    (q[i][j]).append(\
sum([A[i][l][k]+B[i][j][l]+C[l][j][k] for l in range(n_a_cols)]))
        return q

    else :
        raise ValueError, "Hypermatrix dimension mismatch."

def HypermatrixKroneckerProduct(A, B, C):
    """
    Outputs a list of lists associated with the ternary
    product of the input hypermatrices A, B and C.

    EXAMPLES:
    ::
        sage: X=HM(2,1,1,'x'); Y=HM(1,2,1,'y'); Z=HM(1,1,2,'z') 
        sage: T = HM(HypermatrixKroneckerProduct(X.listHM(),Y.listHM(),Z.listHM()))
        [[[x000*y000*z000, x000*y000*z001], [x000*y010*z000, x000*y010*z001]], [[x100*y000*z000, x100*y000*z001], [x100*y010*z000, x100*y010*z001]]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_a_rows = len(A)
    n_a_cols = len(A[0])
    n_a_dpts = len(A[0][0])

    n_b_rows = len(B)
    n_b_cols = len(B[0])
    n_b_dpts = len(B[0][0])

    n_c_rows = len(C)
    n_c_cols = len(C[0])
    n_c_dpts = len(C[0][0])

    # Test for zero dimension
    if n_a_rows*n_b_rows*n_c_rows>0 and n_a_cols*n_b_cols*n_c_cols>0 and n_a_dpts*n_b_dpts*n_c_dpts>0:
        # Initialization of the hypermatrix
        q = []
        for i in range(n_a_rows*n_b_rows*n_c_rows):
            q.append([])
        for i in range(len(q)):
            for j in range(n_a_cols*n_b_cols*n_c_cols):
                (q[i]).append([])
        for i in range(len(q)):
            for j in range(len(q[i])):
                for k in range(n_a_dpts*n_b_dpts*n_c_dpts):
                    (q[i][j]).append(0)
        for i0 in range(n_a_rows):
            for i1 in range(n_b_rows):
                for i2 in range(n_c_rows):
                    for j0 in range(n_a_cols):
                        for j1 in range(n_b_cols):
                            for j2 in range(n_c_cols):
                                for k0 in range(n_a_dpts):
                                    for k1 in range(n_b_dpts):
                                        for k2 in range(n_c_dpts):
                                            q[n_a_rows*n_b_rows*i2+n_a_rows*i1+i0][n_b_cols*n_c_cols*j0+n_b_cols*j2+j1][n_c_cols*n_a_dpts*k1+n_c_dpts*k0+k2]=A[i0][i1][i2]*B[j0][j1][j2]*C[k0][k1][k2]
        return q

    else :
        raise ValueError, "Hypermatrix dimension mismatch."

def HypermatrixProductB(A, B, C, D):
    """
    Outputs a list of lists associated with the ternary
    product the input hypermatrices A, B and C with
    background hypermatrix D.

    EXAMPLES:
    ::
        sage: M = HypermatrixProductB(A, B, C, D); M


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_a_rows = len(A)
    n_a_cols = len(A[0])
    n_a_dpts = len(A[0][0])

    n_b_rows = len(B)
    n_b_cols = len(B[0])
    n_b_dpts = len(B[0][0])

    n_c_rows = len(C)
    n_c_cols = len(C[0])
    n_c_dpts = len(C[0][0])

    n_d_rows = len(D)
    n_d_cols = len(D[0])
    n_d_dpts = len(D[0][0])

    # Test for dimension match
    if \
n_a_rows==n_b_rows and n_b_cols==n_c_cols and n_c_dpts==n_a_dpts and \
n_a_cols==n_b_dpts and n_b_dpts==n_c_rows and n_a_cols==n_d_rows and \
n_a_cols==n_d_cols and n_a_cols==n_d_dpts:
        # Initialization of the hypermatrix
        q = []
        for i in range(n_a_rows):
            q.append([])
        for i in range(len(q)):
            for j in range(n_b_cols):
                (q[i]).append([])
        for i in range(len(q)):
            for j in range(len(q[i])):
                for k in range(n_c_dpts):
                    (q[i][j]).append(\
sum([A[i][l0][k]*B[i][j][l1]*C[l2][j][k]*D[l0][l1][l2] for l0 in range(n_d_rows)\
for l1 in range(n_d_cols) for l2 in range(n_d_dpts)]))
        return q

    else :
        raise ValueError, "Hypermatrix dimension mismatch."

def HypermatrixCyclicPermute(A):
    """
    Outputs a list of lists associated with the hypermatrix
    with entries index cycliclly permuted.

    EXAMPLES:
    ::
        sage: M = HypermatrixCyclicPermute(A); M


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_q_rows = len(A[0])
    n_q_cols = len(A[0][0])
    n_q_dpts = len(A)

    # Initialization of the hypermatrix
    q = []
    for i in range(n_q_rows):
        q.append([])
    for i in range(len(q)):
        for j in range(n_q_cols):
            (q[i]).append([])
    for i in range(len(q)):
        for j in range(len(q[i])):
            for k in range(n_q_dpts):
                (q[i][j]).append(A[k][i][j])
    return q

def HypermatrixKroneckerDelta(nr):
    """
    Generates a list of lists associated with the nr x nr x nr
    Kronecker Delta hypermatrix.

    EXAMPLES:
    ::
        sage: M = HypermatrixKroneckerDelta(2); M


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_q_rows = nr
    n_q_cols = nr
    n_q_dpts = nr

    # Test for dimension match
    if n_q_rows > 0 and n_q_cols > 0 and n_q_dpts >0:
        # Initialization of the hypermatrix
        q = []
        for i in range(n_q_rows):
            q.append([])
        for i in range(len(q)):
            for j in range(n_q_cols):
                (q[i]).append([])
        for i in range(len(q)):
            for j in range(len(q[i])):
                for k in range(n_q_dpts):
                    if i==j and i==k:
                        (q[i][j]).append(1)
                    else:
                        (q[i][j]).append(0)
        return q

    else :
        raise ValueError, "Input dimensions "+\
str(nr)+" must be a non-zero positive integer."

def HypermatrixPermutation(s):
    """
    Generates a list of lists associated with the permutation
    hypermatrix deduced from sigma.

    EXAMPLES:
    ::
        sage: M = HypermatrixPermutation([0,2,1]); M


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    n = len(s)
    # Setting the dimensions parameters.
    n_q_rows = n
    n_q_cols = n
    n_q_dpts = n

    # Test for dimension match
    if n_q_rows > 0 and n_q_cols > 0 and n_q_dpts >0:
        # Initialization of the hypermatrix
        q = []
        T = HypermatrixKroneckerDelta(n)
        U = HypermatrixGenerateAllOne(n,n,n)
        Id= HypermatrixProduct(U,U,T)
        Id= HypermatrixCyclicPermute(Id)
        for i in range(n):
            q.append(Id[s[i]])
        return HypermatrixCyclicPermute(HypermatrixCyclicPermute(q))

    else :
        raise ValueError, "Input dimensions "+\
str(n)+" must be a non-zero positive integer."

def DiagonalHypermatrix(Mtrx):
    """
    Outputs a diagonal third order hypermatrix
    constructed using the input square matrix
    to enforce the symmetry constraint we will
    only take entry from the lower triangular
    part of the input matrix.

     EXAMPLES:
    ::
        sage: var('a00, a11, a01')
        sage: Mtrx = Matrix(Sr,[[a00,a01],[a01,a11]])
        sage: d = DiagonalHypermatrix(Mtrx)

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initialization of the dimensions
    n = min(Mtrx.nrows(),Mtrx.ncols())
    n_d_rows = n
    n_d_cols = n
    n_d_dpts = n

    # Initialization of the identity permutations hypermatrix
    D = HypermatrixPermutation(range(n))

    # Filling up the entries of the hypermatrix.
    for i in range(n_d_rows):
        for j in range(n_d_cols):
            for k in range(n_d_dpts):
                if (D[i][j][k] != 0):
                    D[i][j][k] = Mtrx[min(i,k),max(i,k)]
    return D

def Orthogonal2x2x2Hypermatrix(t):
    """
    Outputs an orthogonal third order hypermatrix
    of size 2 by 2 by 2.

     EXAMPLES:
    ::
        sage: t=var('t')
        sage: Orthogonal2x2x2Hypermatrix(t)
        [[[cos(t)^(2/3), sin(t)^(2/3)], [sin(t)^(2/3), cos(t)^(2/3)]], [[-sin(t)^(2/3), cos(t)^(2/3)], [sin(t)^(2/3), sin(t)^(2/3)]]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    return [[[cos(t)^(2/3),sin(t)^(2/3)],[sin(t)^(2/3), cos(t)^(2/3)]],\
[[-sin(t)^(2/3),cos(t)^(2/3)],[sin(t)^(2/3),sin(t)^(2/3)]]]
def Orthogonal3x3x3Hypermatrix(t1,t2):
    """
    Outputs an orthogonal third order hypermatrix
    of size 3 by 3 by 3.

     EXAMPLES:
    ::
        sage: t1,t2=var('t1,t2')
        sage: Orthogonal3x3x3Hypermatrix(t1,t2)

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    c1=cos(t1)^(2/3)
    s1=sin(t1)^(2/3)
    c2=cos(t2)^(2/3)
    s2=sin(t2)^(2/3)
    return [[[c1,s1*c2,0],[s1*c2,s1*s2,0],[s1*s2,exp(-I*2*pi/3)*c1,0]],\
[[s1*s2,c1,exp(-I*2*pi/3)*s1*c2],[exp(I*2*pi/3)*c1,s1*c2,s1*s2],\
[s1*c2,s1*s2,c1]],[[0,s1*s2,c1],[0,c1,s1*c2],[0,exp(I*2*pi/3)*s1*c2,s1*s2]]]


#@cached_function
# I should find a way to hash my lists
def HypermatrixCayleyHamiltonList(A,n):
    """
    Outpts a list of hypermatrices of all product
    composition of order n.

     EXAMPLES:
    ::
        sage: A = HypermatrixGenerate(2,2,2,'a')
        sage: L = HypermatrixCayleyHamiltonList(A,3)

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    if n == 1:
        return [A]
    else:
        gu = []
        for i in range(1,n,2):
            for j in range(1,n-i,2):
                gu = gu + [HypermatrixProduct(g1,g2,g3) \
for g1 in HypermatrixCayleyHamiltonList(A,i) \
for g2 in HypermatrixCayleyHamiltonList(A,j) \
for g3 in HypermatrixCayleyHamiltonList(A,n-(i+j))]
        return gu

def HypermatrixCayleyHamiltonListII(A,n):
    """
    Outpts a list of hypermatrices of all product
    composition of order n.

     EXAMPLES:
    ::
        sage: A = HM(2,2,2,'a')
        sage: L = HypermatrixCayleyHamiltonListII(A,3)

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    if A.order()==3:
        if n == 1:
            return [A]
        else:
            gu = []
            for i in range(1,n,2):
                for j in range(1,n-i,2):
                    #gu = gu + [HypermatrixProduct(g1,g2,g3) for g1 in HypermatrixCayleyHamiltonListII(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,n-(i+j))]
                    gu = gu + [GeneralHypermatrixProduct(g1,g2,g3) for g1 in HypermatrixCayleyHamiltonListII(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,n-(i+j))]
            return gu

    elif A.order()==4:
        if n == 1:
            return [A]
        else:
            gu = []
            for i in range(1,n,2):
                for j in range(1,n-i,2):
                    for k in range(1,n-i-j,2):
                        gu = gu + [GeneralHypermatrixProduct(g1,g2,g3,g4) for g1 in HypermatrixCayleyHamiltonListII(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,k) for g4 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k))]
    else :
        raise ValueError, "Not supported for order > 3 and for non cube hypermpatrix of order 3 "
   
def HypermatrixCayleyHamiltonListIII(n):
    """
    Outpts a list of hypermatrices of all product
    composition of order n.

     EXAMPLES:
    ::
        sage: L = HypermatrixCayleyHamiltonListIII(A,3)

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    if n == 1:
        return ['A.listHM()']
    else:
        gu = []
        for i in range(1,n,2):
            for j in range(1,n-i,2):
                gu = gu + ['HypermatrixProduct('+str(g1)+','+str(g2)+','+str(g3)+')' \
for g1 in HypermatrixCayleyHamiltonListIII(i) \
for g2 in HypermatrixCayleyHamiltonListIII(j) \
for g3 in HypermatrixCayleyHamiltonListIII(n-(i+j))]
        # Creating the string corresponding to the file name
        filename = 'output.txt'
        # Opening the file
        f = open(filename,'w')
        f.write(str(gu).replace("'",""))
        f.close()
        return gu

def ConstraintFormator(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs matrix
    and the right hand side vector associate
    with the matrix formulation of the constraints.

    EXAMPLES:
    ::
        sage: x, y = var('x,y')
        sage: CnstrLst = [x+y==1, x-y==2]
        sage: VrbLst = [x, y]
        sage: [A,b] = ConstraintFormator(CnstrLst, VrbLst)

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initializing the Matrix
    A=Matrix(CC,len(CnstrLst),len(VrbLst),zero_matrix(len(CnstrLst),len(VrbLst)))
    b=vector(CC, [eq.rhs() for eq in CnstrLst]).column()
    for r in range(len(CnstrLst)):
        for c in range(len(VrbLst)):
            A[r,c] = diff((CnstrLst[r]).lhs(),VrbLst[c])
    return [A,b]

def ConstraintFormatorII(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs matrix
    and the right hand side vector associate
    with the matrix formulation of the constraints.

    EXAMPLES:
    ::
        sage: x,y = var('x,y')
        sage: CnstrLst = [x+y==1, x-y==2]
        sage: VrbLst = [x, y]
        sage: [A,b] = ConstraintFormatorII(CnstrLst, VrbLst)

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initializing the Matrix
    A=Matrix(SR,len(CnstrLst),len(VrbLst),zero_matrix(len(CnstrLst),len(VrbLst)))
    b=vector(SR, [eq.rhs() for eq in CnstrLst]).column()
    for r in range(len(CnstrLst)):
        for c in range(len(VrbLst)):
            A[r,c] = diff((CnstrLst[r]).lhs(),VrbLst[c])
    return [A,b]

def HypermatrixPseudoInversePairs(A,B):
    """
     Outputs the pseudo inverse pairs associated with the input pairs of matrices

    EXAMPLES:
    ::
        sage: A1=[[[0.1631135370902057,0.11600112072013125],[0.9823708115400902,0.39605960486710756]]\
,[[0.061860929755424676,0.2325542810173995],[0.39111210957450926,0.2019809359102137]]]
        sage: A2=[[[0.15508921433883183,0.17820377184410963],[0.48648171594508205,0.01568017636082064]]\
,[[0.8250247759993575,0.1938307874191597],[0.23867299119274843,0.3935578730402869]]]
        sage: [B1,B2]=HypermatrixPseudoInversePairs(A1,A2)

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    sz = len(A)

    # Initializing the list of linear constraints
    CnstrLst = []

    # Initilizing the variable list
    Vrbls  = [var('ln_al'+str(i)+str(j)+str(k)) \
for i in range(sz) for j in range(sz) for k in range(sz)]+\
[var('ln_bt'+str(i)+str(j)+str(k)) for i in range(sz) for j in range(sz) \
for k in range(sz)]

    for m in range(sz):
        for p in range(sz):
            for n in range(sz):
                V=Matrix(CC, sz, sz, [(A[m][k1][k0])*(B[k0][k1][p]) \
for k0 in range(sz) for k1 in range(sz)]).inverse()
                CnstrLst=CnstrLst+[\
var('ln_al'+str(m)+str(n)+str(k1))+var('ln_bt'+str(k1)+str(n)+str(p))==\
ln(V[k1,n])  for k1 in range(sz)]
    [A,b]=ConstraintFormator(CnstrLst,Vrbls)

    # Importing the Numerical Python package
    # for computing the matrix pseudo inverse
    import numpy

    sln = matrix(numpy.linalg.pinv(A))*b
    R1 = HypermatrixGenerateAllZero(sz,sz,sz)
    for i in range(sz):
        for j in range(sz):
            for k in range(sz):
                R1[i][j][k] = exp(sln[i*sz^2+j*sz^1+k*sz^0,0])
    R2 = HypermatrixGenerateAllZero(sz, sz, sz)
    for i in range(sz):
        for j in range(sz):
            for k in range(sz):
                R2[i][j][k] = exp(sln[sz^3+i*sz^2+j*sz^1+k*sz^0,0])
    return [R1,R2]

@cached_function
def C(n):
    """
    Counts the number of product composition involving the input
    hypermatrix n times.

    EXAMPLES:
    The input n must be greater than 0
    ::
        sage: C(3)
        1

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    if n == 1 :
        return 1
    else :
        return sum([C(i)*C(j)*C(n-i-j) for i in range(1,n,2) for j in range(1,n-i,2)])

# Deinfition of the third order hypermatrix class
class HM:
    """HM class"""
    def __init__(self,*args):
        if len(args) == 1:
            inp = args[0]
            if type(inp)==type(Matrix(SR,2,1,[var('x'),var('y')])) or \
type(inp)==type(Matrix(RR,2,1,[1,2])) or type(inp)==type(Matrix(CC,2,1,[1,1])):
                self.hm=DiagonalHypermatrix(inp)
            elif type(inp) == list:
                self.hm = inp
            elif type(inp) == type("abcd"):
                # Importing the numerical library
                import pylab, numpy
                from scipy import misc
                # Defining the input Pictures
                X  = pylab.imread(inp)
                sz = max(len(X),len(X[0]),len(X[0][0]))
                args = (len(X), len(X[0]), len(X[0][0]))
                #print args 
                T = apply(HypermatrixGenerateAllZero, args)
                # Filling up the image Hypermatrix
                for i in range(len(X)):
                    for j in range(len(X[0])):
                        for k in range(len(X[0][0])):
                            T[i][j][k] = X[i,j,k]
                self.hm = T

            else:
                raise ValueError, "Expected list as input when generating diagonal hypermatrix"
            return
    	s = args[-1]
        dims = args[:-1]
        if s == 'one':
            self.hm = apply(HypermatrixGenerateAllOne, dims)
        elif s == 'zero':
            self.hm = apply(HypermatrixGenerateAllZero, dims)
        elif s == 'ortho':
            if len(dims) == 1:
                self.hm=Orthogonal2x2x2Hypermatrix(dims[0])
            elif len(dims) == 2:
                self.hm=Orthogonal3x3x3Hypermatrix(dims[0],dims[1])
            else:
                raise ValueError, "ortho not supported for order %d tensors" % len(dims)
        elif s == 'perm':
            self.hm=HypermatrixPermutation(dims[0])
        elif s == 'kronecker':
            self.hm=HypermatrixKroneckerDelta(dims[0])
        elif s == 'sym':
            if len(dims) == 2:
                self.hm=SymHypermatrixGenerate(dims[0],dims[1])
            else:
                raise ValueError, "kronecker not supported for order %d tensors" % len(dims)
        else:
            self.hm=apply(HypermatrixGenerate, args)

    def __repr__(self):
        return `self.hm`

    def __add__(self, other):
        return GeneralHypermatrixAdd(self,other)

    def __neg__(self):
        return GeneralHypermatrixScale(self.hm,-1)

    def __sub__(self, other):
        return GeneralHypermatrixAdd(self, GeneralHypermatrixScale(other,-1))

    def __mul__(self, other):
        if other.__class__.__name__=='HM':
            return HM(GeneralHypermatrixHadamardProduct(self,other))
        elif other.__class__.__name__=='tuple':
            # This function takes a a list as intput
            l = other
            return GeneralHypermatrixProduct(self,*l)
        else: 
            return GeneralHypermatrixScale(self,other)

    def __rmul__(self, a):
        return self*a

    def __getitem__(self,i):
        if i.__class__.__name__=='tuple':
            tmp = self.hm
            for j in i:
                tmp = tmp[j]
            return tmp

    def __setitem__(self, i, v):
        if   i.__class__.__name__=='tuple':
	    tmp = self.hm
            while len(i)>1:
                tmp = tmp[i[0]]
                i = i[1:]
            tmp[i[0]] = v

    def __call__(self, *inpts):
        # This function takes a a list as intput
        return GeneralHypermatrixProduct(self, *inpts)

    def hprod(self,*inpts):
        # This function takes a a list as intput
        return GeneralHypermatrixProduct(self,*inpts)

    def hprod3b(self, b, c, t):
        return HM(HypermatrixProductB(self.hm, b.hm, c.hm, t.hm))

    def elementwise_product(self,B):
        return GeneralHypermatrixHadamardProduct(self,B)

    def elementwise_exponent(self,s):
        return GeneralHypermatrixExponent(self,s)

    def elementwise_base_exponent(self,s):
        return GeneralHypermatrixBaseExponent(self,s)

    def expand(self):
        return GeneralHypermatrixExpand(self)

    def simplify(self):
        return GeneralHypermatrixSimplify(self)

    def subs(self, Dct):
        return GeneralHypermatrixSubstitute(self, Dct)

    def subsn(self, Dct):
        return GeneralHypermatrixSubstituteN(self, Dct)

    def transpose(self, i=1):
        t = Integer(mod(i, self.order()))
        A = self 
        for i in range(t):
            A = GeneralHypermatrixCyclicPermute(A)
        return A

    def nrows(self):
        return len(self.hm)

    def ncols(self):
        return len(self.hm[0])

    def ndpts(self):
        return len(self.hm[0][0])

    def Print(self):
        if self.order() == 3:
            L = self.listHM()
            for m in L:
                print '\n'+Matrix(m).str()
        else:
            raise ValueError, "not supported for order %d tensors" % len(dims)
            
    def n(self,i):
        if i == 0:
            return self.nrows()
        elif i == 1:
            return self.ncols()
        elif i == 2:
            return self.ndpts()
        else:
            tmp = self.listHM()
            for j in range(i):
                tmp = tmp[0]
            return len(tmp)

    def list(self):
        lst = []
        l = [self.n(i) for i in range(self.order())]
        # Main loop canonicaly listing the elements
        for i in range(prod(l)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entry = [mod(i,l[0])]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            # Appending to the list
            lst.append(self[tuple(entry)])
        return lst

    def listHM(self):
        return self.hm

    def cayley_hamilton_list(self,n):
        tmp = HypermatrixCayleyHamiltonListII(self,n)
        return [h for h in tmp]

    def cayley_hamilton_mtrxII(self,itr,bnd):
        tmp = []
        for i in range(itr):
            tmp = tmp + HypermatrixCayleyHamiltonList(self.hm, 2*i+1)
        return Matrix([HM(h).list() for h in tmp[0:bnd]]).transpose()

    def cayley_hamilton_mtrxI(self):
        mn = min(self.nrows(), self.ncols(), self.ndpts())
        mx = max(self.nrows(), self.ncols(), self.ndpts())
        if self.order() == 3 and mn == mx:
            tmp = []
            i = 0
            while len(tmp) < mn^3:
                tmp = tmp + HypermatrixCayleyHamiltonList(self.hm, 2*i+1)
                i = i+1
            return Matrix([HM(h).list() for h in tmp[0:mn^3]]).transpose()
        else :
            raise ValueError, "Not supported for order > 3 and for non cube hypermpatrix of order 3 "

    def cayley_hamilton_coef(self, v):
        # Computing the matrix of powers
        M = self.cayley_hamilton_mtrxI()
        dtrm = Deter(M)
        L = []
        for i in range(M.nrows()):
            T = copy(M)
            T[:,i] = v[:,0]
            L.append(Deter(T)/dtrm)
        return L
 
    def cayley_hamilton_coefN(self, v):
        # Computing the matrix of powers
        M = Matrix(CC,self.cayley_hamilton_mtrxI())
        dtrm = M.det()
        L = []
        for i in range(M.nrows()):
            T = copy(M)
            T[:,i] = v[:,0]
            L.append(T.det()/dtrm)
        return L

    def cayley_hamilton_coefNI(self):
        # Initializing thr size of the hypermatrix
        sz = min([self.n(i) for i in range(self.order())])
        # Creating the list of hypermatrix powers
        TmpL = HypermatrixCayleyHamiltonList(self.listHM(), 1+sz^3)
        # Filling up the left hand side
        v = Matrix(CC,sz^3,1,HM(TmpL[len(TmpL)-1]).list())
        # Filling up the matrix
        M = Matrix(CC, zero_matrix(sz^3,sz^3))
        for l in range(sz^3):
            for i in range(sz):
                for j in range(sz):
                    for k in range(sz):
                        M[i*sz^2+j*sz+k,l] = (TmpL[l])[i][j][k] 
        dtrm = M.det()
        L = []
        for i in range(M.nrows()):
            T = copy(M)
            T[:,i] = v[:,0]
            L.append(T.det()/dtrm)
        return L

    def order(self):
        cnt = 0
        H = self.listHM()
        while type(H) == type([]):
            H = H[0]
            cnt = cnt+1
        return cnt

    def conj(self, k):
        Tmp = self.listHM()
        for r in range(len(Tmp)):
            for c in range(len(Tmp[0])):
                for d in range(len(Tmp[0][0])): 
                    Tmp[r][c][d]=conj3(Tmp[r][c][d], k)
        return HM(Tmp)

    def zero_padd(self):
        sz  = max(self.nrows(), self.ncols(), self.ndpts())
        Tmp = HM(sz,sz,sz,'zero') 
        for i in range(self.nrows()):
            for j in range(self.ncols()):
                for k in range(self.ndpts()):
                    Tmp[i,j,k]=self.hm[i][j][k]
        return Tmp

    def fill_with(self,T):
        if T.nrows()>=self.nrows() or T.ncols()>=self.ncols or T.ndpts()>=self.ndpts():
            for r in range(self.nrows()):
                for c in range(self.ncols()):
                    for d in range(self.ndpts()):
                        self.hm[r][c][d]=T[r,c,d]
        else:
            raise ValueError, "Expected the input 3 hypermatrix to have larger dimensions in all directions"
    
    def show(self):
        import pylab, numpy
        from scipy import misc
        # This line of code corresponds
        # to the lazy way of getting
        # the image size
        #X = pylab.imread(image_name)
        X = misc.toimage(pylab.array(self.listHM()))
        g = graphics_array([[matrix_plot(X)]])
        g.show()

    def save(self,filename):
        import pylab, numpy
        from scipy import misc
        # This line of code corresponds
        # to the lazy way of getting
        # the image size
        #X = pylab.imread(image_name)
        X = misc.toimage(pylab.array(self.listHM()))
        X.save(filename)


# Implementing the General version of the Hypermatrix product 
# I now need to test this implementation
def GeneralHypermatrixProduct(*args):
    # Initialization of the list specifying the dimensions of the output
    l = [(args[i]).n(i) for i in range(len(args))]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the assignement
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # computing the Hypermatrix product
        if len(args)<2:
            raise ValueError, "The number of operands must be >= 2"
        elif len(args) >= 2:
            Rh[tuple(entry)]=sum(\
[prod([args[s][tuple(entry[0:Integer(mod(s+1,len(args)))]+[t]+entry[Integer(mod(s+2,len(args))):])] for s in range(len(args)-2)]+\
[args[len(args)-2][tuple(entry[0:len(args)-1]+[t])]]+[args[len(args)-1][tuple([t]+entry[1:])]]) for t in range((args[0]).n(1))])
    return Rh

# Implementing the General version of the Hypermatrix product 
# I now need to test this implementation
def GeneralHypermatrixProductB(*args):
    # Initialization of the list specifying the dimensions of the output
    l = [(args[i]).n(i) for i in range(len(args)-1)]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the assignement
    # Initializing the background hypermatrix
    B = args[len(args)-1]
    args = tuple([args[id] for id in range(len(args)-1)])
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # computing the Hypermatrix product
        if len(args) < 2:
            raise ValueError, "The number of operands must be >= 2"
        elif len(args) >= 2:
            Rh[tuple(entry)] = 0
            l2 = [B.n(sz) for sz in range(B.order())]
            for j in range(prod(l2)):
                # Turning the index j into an hypermatrix array location using the decimal encoding trick
                entry2 = [mod(j,l2[0])]
                sm2 = Integer(mod(j,l2[0]))
                for z in range(len(l2)-1):
                    entry2.append(Integer(mod(Integer((j-sm2)/prod(l2[0:z+1])),l2[z+1])))
                    sm2 = sm2+prod(l2[0:z+1])*entry2[len(entry2)-1]
                Rh[tuple(entry)]= Rh[tuple(entry)]+prod([args[s][tuple(entry[0:Integer(mod(s+1,len(args)))]+[entry2[s]]+entry[Integer(mod(s+2,len(args))):])] for s in range(len(args)-2)]+[args[len(args)-2][tuple(entry[0:len(args)-1]+[entry2[len(entry2)-1]])]]+[args[len(args)-1][tuple([entry2[0]]+entry[1:])]])*B[tuple(entry2)]
    return Rh

# Implementing the General version of the Hypermatrix logarithmic product
# I now need to test this implementation
def GeneralHypermatrixLogProduct(*args):
    # Initialization of the list specifying the dimensions of the output
    l = [(args[i]).n(i) for i in range(len(args))]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the assignement
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # computing the Hypermatrix product
        if len(args)<2:
            raise ValueError, "The number of operands must be >= 2"
        elif len(args) >= 2:
            Rh[tuple(entry)]=sum(\
[sum([args[s][tuple(entry[0:Integer(mod(s+1,len(args)))]+[t]+entry[Integer(mod(s+2,len(args))):])] for s in range(len(args)-2)]+\
[args[len(args)-2][tuple(entry[0:len(args)-1]+[t])]]+[args[len(args)-1][tuple([t]+entry[1:])]]) for t in range((args[0]).n(1))])
    return Rh

def GeneralHypermatrixCyclicPermute(A):
    # Initialization of the list specifying the dimensions of the output
    l = [A.n(i) for i in range(A.order())]
    l = l[1:]+[l[0]]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # Performing the transpose
        Rh[tuple(entry)]=A[tuple([entry[len(entry)-1]]+entry[:len(entry)-1])]
    return Rh

def GeneralHypermatrixScale(A,s):
    # Initialization of the list specifying the dimensions of the output
    l = [A.n(i) for i in range(A.order())]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=s*A[tuple(entry)]
    return Rh

def GeneralHypermatrixExponent(A,s):
    # Initialization of the list specifying the dimensions of the output
    l = [A.n(i) for i in range(A.order())]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=(A[tuple(entry)])^s
    return Rh

def GeneralHypermatrixBaseExponent(A,s):
    # Initialization of the list specifying the dimensions of the output
    l = [A.n(i) for i in range(A.order())]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=s^(A[tuple(entry)])
    return Rh

def GeneralHypermatrixExpand(A):
    # Initialization of the list specifying the dimensions of the output
    l = [A.n(i) for i in range(A.order())]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=(A[tuple(entry)]).expand()
    return Rh

def GeneralHypermatrixSimplify(A):
    """
    Performs the symbolic simplification of the expressions
    associated with the hypermatrix entries. 

    EXAMPLES:
    ::
        sage: x,y = var('x,y') 
        sage:((x+y)^2*HM(2,2,2,'one')).simplify() 
        [[[x^2+2*x*y+y^2,x^2+2*x*y+y^2],[x^2+2*x*y+y^2,x^2+2*x*y+y^2]],[[x^2+2*x*y+y^2,x^2+2*x*y+y^2],[x^2+2*x*y+y^2,x^2+2*x*y+y^2]]]
    AUTHORS:
    - Edinah K. Gnang
    """

    # Initialization of the list specifying the dimensions of the output
    l = [A.n(i) for i in range(A.order())]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=(A[tuple(entry)]).simplify_full()
    return Rh

def GeneralHypermatrixSubstitute(A, Dct):
    """
    Procedure for computes the substitution in the Hypermatrix entries
    the inputs are the corresponding Hypermatric and a dictionary 
    datastructure.

    EXAMPLES:
    ::
        sage: A = HM(2,'a','sym')
        sage: R = GeneralHypermatrixAdd(A,dict([(a011,var('x')),(a001,var('y')),(a000,var('z')),(a111,var('t'))]));R
        [[[z, y], [y, x]], [[y, x], [x, t]]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initialization of the list specifying the dimensions of the output
    l = [A.n(i) for i in range(A.order())]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=(A[tuple(entry)]).subs(Dct)
    return Rh

def GeneralHypermatrixSubstituteN(A, Dct):
    """
    Procedure for computes the substitution in the Hypermatrix entries
    the inputs are the corresponding Hypermatric and a dictionary 
    datastructure.

    EXAMPLES:
    ::
        sage: A = HM(2,'a','sym')
        sage: R = GeneralHypermatrixAdd(A,dict([(a011,var('x')),(a001,var('y')),(a000,var('z')),(a111,var('t'))]));R
        [[[z, y], [y, x]], [[y, x], [x, t]]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initialization of the list specifying the dimensions of the output
    l = [A.n(i) for i in range(A.order())]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=N((A[tuple(entry)]).subs(Dct))
    return Rh

def GeneralHypermatrixCopy(A):
    """
    Procedure for computing Hypermatrix Hadamard addition.

    EXAMPLES:
    ::
        sage: A = HM(2,2,2,'a');B = HM(2,2,2,'b')
        sage: R = GeneralHypermatrixCopy(M, A+B);R
        [[[a000+b000,a001+b001],[a010+b010,a011+b011]],[[a100+b100,a101+b101],[a110+b110,a111+b111]]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initialization of the list specifying the dimensions of the output
    l = [A.n(i) for i in range(A.order())]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=A[tuple(entry)]
    return Rh

def GeneralHypermatrixAdd(A,B):
    """
    Procedure for computing Hypermatrix Hadamard addition.

    EXAMPLES:
    ::
        sage: A = HM(2,2,2,'a');B = HM(2,2,2,'b')
        sage: R = GeneralHypermatrixAdd(A,B);R
        [[[a000*b000,a001*b001],[a010*b010,a011*b011]],[[a100*b100,a101*b101],[a110*b110,a111*b111]]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initialization of the list specifying the dimensions of the output
    l = [A.n(i) for i in range(A.order())]
    s = [B.n(i) for i in range(B.order())]
    # Testing the dimensions 
    x = var('x')
    if(sum([l[i]*x^i for i in range(len(l))])==sum([s[i]*x^i for i in range(len(s))])):
        # Initializing the input for generating a symbolic hypermatrix
        inpts = l+['zero']
        # Initialization of the hypermatrix
        Rh = HM(*inpts)
        # Main loop performing the transposition of the entries
        for i in range(prod(l)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entry = [mod(i,l[0])]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            Rh[tuple(entry)]=A[tuple(entry)]+B[tuple(entry)]
        return Rh
    else:
        raise ValueError, "The Dimensions of the input hypermatrices must match."
 
def GeneralHypermatrixHadamardProduct(A,B):
    """
    Procedure for computing Hypermatrix Hadamard products.

    EXAMPLES:
    ::
        sage: A = HM(2,2,2,'a');B = HM(2,2,2,'b')
        sage: R = GeneralHypermatrixHadamardProduct(A,B);R
        [[[a000*b000,a001*b001],[a010*b010,a011*b011]],[[a100*b100,a101*b101],[a110*b110,a111*b111]]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initialization of the list specifying the dimensions of the output
    l = [A.n(i) for i in range(A.order())]
    s = [B.n(i) for i in range(B.order())]
    # Testing the dimensions 
    x = var('x')
    if(sum([l[i]*x^i for i in range(len(l))])==sum([s[i]*x^i for i in range(len(s))])):
        # Initializing the input for generating a symbolic hypermatrix
        inpts = l+['zero']
        # Initialization of the hypermatrix
        Rh = HM(*inpts)
        # Main loop performing the transposition of the entries
        for i in range(prod(l)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entry = [mod(i,l[0])]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            Rh[tuple(entry)]=A[tuple(entry)]*B[tuple(entry)]
        return Rh
    else:
        raise ValueError, "The Dimensions of the input hypermatrices must match."

def GeneralHypermatrixKroneckerDelta(od, sz):
    # Initialization of the list specifying the dimensions of the output
    l = [sz for i in range(od)] 
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        if len(Set(entry)) == 1:
            Rh[tuple(entry)] = 1
    return Rh

def GeneralOrthogonalHypermatrixU(od):
    """
    Generates an orthogonal hypermatrix of the appropriate order
    for which each one of the dimensions are equal to 2.
    The vectors are not normalized.

    EXAMPLES:
    ::
        sage: Q = GeneralOrthogonalHypermatrix(3); Q
        [[[e^(-r1+r3+r6),e^r4],[e^r6,e^r2]],[[-e^(r1+r2-r3-r4+r5),e^r3],[e^r5,e^r1]]]

    AUTHORS:
    - Edinah K. Gnang, Ori Parzanchevski and Yuval Filmus
    """
    # Initializing the hypermatrix
    Q = apply(HM,[2 for i in range(od)]+['q'])
    # Initilizing the list of variable
    VrbLst = Q.list()
    # Reinitializing of Q by exponentiation 
    Q = Q.elementwise_base_exponent(e)
    # Computing the product
    Eq = apply(GeneralHypermatrixProduct, [Q.transpose(j) for j in range(od,0,-1)])
    # Writting up the constraints
    LeQ = (Set(Eq.list())).list()
    # Removing the normalization constraints
    LeQ.remove(e^(od*var('q'+''.join(['0' for i in range(od)])))+\
e^(od*var('q01'+''.join(['0' for i in range(od-2)]))))
    LeQ.remove( e^(od*var('q10'+''.join(['1' for i in range(od-2)])))+\
e^(od*var('q'+''.join(['1' for i in range(od)]))))
    # Filling up the linear constraints
    CnstrLst= [] 
    for f in LeQ:
        CnstrLst.append(ln((f.operands())[0]).simplify_exp()-I*pi-ln((f.operands())[1]).simplify_exp()==0)
    # Directly solving the constraints
    Sl = solve(CnstrLst,VrbLst)
    # Main loop performing the substitution of the entries
    Lr = [var('r'+str(i)) for i in range(1,2^od+1)]
    l = [Q.n(i) for i in range(Q.order())]
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Q[tuple(entry)]=Q[tuple(entry)].subs(dict(map(lambda eq: (eq.lhs(),eq.rhs()), Sl[0]))).simplify_exp()
    return Q

def GeneralOrthogonalHypermatrix(od):
    """
    Generates an orthogonal hypermatrix of the appropriate order
    for which each one of the dimensions are equal to 2.
    The vectors are not normalized.

    EXAMPLES:
    ::
        sage: Q = GeneralOrthogonalHypermatrix(3); Q
        [[[e^(-r1+r3+r6),e^r4],[e^r6,e^r2]],[[-e^(r1+r2-r3-r4+r5),e^r3],[e^r5,e^r1]]]

    AUTHORS:
    - Edinah K. Gnang, Ori Parzanchevski and Yuval Filmus
    """
    if od == 2:
        nrm = sqrt(exp(2*var('r1'))+exp(2*var('r2')))
        return HM([[exp(var('r1'))/nrm, exp(var('r2'))/nrm],[-exp(var('r2'))/nrm, exp(var('r1'))/nrm]])
    else :
        # Initializing the hypermatrix
        Q = apply(HM,[2 for i in range(od)]+['q'])
    
        # Initilizing the list of variable
        VrbLst = Q.list()
    
        # Reinitializing of Q by exponentiation 
        Q = Q.elementwise_base_exponent(e)
    
        # Computing the product
        Eq = apply(GeneralHypermatrixProduct, [Q.transpose(j) for j in range(od,0,-1)])
    
        # Writting up the constraints
        LeQ = (Set(Eq.list())).list()
        
        # Removing the normalization constraints
        LeQ.remove(e^(od*var('q'+''.join(['0' for i in range(od)])))+ e^(od*var('q01'+''.join(['0' for i in range(od-2)]))))
        LeQ.remove( e^(od*var('q10'+''.join(['1' for i in range(od-2)])))+ e^(od*var('q'+''.join(['1' for i in range(od)]))))
    
        # Filling up the linear constraints
        CnstrLst= [] 
        for f in LeQ:
            CnstrLst.append(ln((f.operands())[0]).simplify_exp()-I*pi-ln((f.operands())[1]).simplify_exp()==0)
    
        # Directly solving the constraints
        Sl = solve(CnstrLst,VrbLst)
    
        # Main loop performing the substitution of the entries
        Lr = [var('r'+str(i)) for i in range(1,2^od+1)]
        l = [Q.n(i) for i in range(Q.order())]
        for i in range(prod(l)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entry = [mod(i,l[0])]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            Q[tuple(entry)]=Q[tuple(entry)].subs(dict(map(lambda eq: (eq.lhs(),eq.rhs()), Sl[0]))).simplify_exp()

        # Initialization of the output hypermatrix
        U = GeneralHypermatrixCopy(Q)
        # first row to normalize 
        entry = [0 for i in range(Q.order())]
        U[tuple(entry)]=Q[tuple(entry)]/sum([ Q[tuple([entry[0]]+[j]+entry[2:])]^Q.order() for j in range(2) ])^(1/Q.order())
        entry[1] = 1
        U[tuple(entry)]=Q[tuple(entry)]/sum([ Q[tuple([entry[0]]+[j]+entry[2:])]^Q.order() for j in range(2) ])^(1/Q.order())

        # last row to normalize 
        entry = [1 for i in range(Q.order())]
        U[tuple(entry)]=Q[tuple(entry)]/sum([ Q[tuple([entry[0]]+[j]+entry[2:])]^Q.order() for j in range(2) ])^(1/Q.order())
        entry[1] = 0
        U[tuple(entry)]=Q[tuple(entry)]/sum([ Q[tuple([entry[0]]+[j]+entry[2:])]^Q.order() for j in range(2) ])^(1/Q.order())
    return U

def DFT_image_resizer(sz, dm):
    if mod(sz,dm) == 0:
        # Initializing the identity matrix of the appropriate size
        Idm = identity_matrix(Integer(sz/dm))
        # Computing the Kronecker product with the hadamard 2x2 matrix
        Rs = Idm.tensor_product(Matrix(SR, dm, dm, [exp(I*2*pi*u*v/dm) for u in range(dm) for v in range(dm)]))
        # Permuting the colum of the matrix in order to put the resized
        # image in the top left corner of the image
        for i in range(1,Integer(sz/dm)):
            tmp = Rs[:,i]
            Rs[:,i] = Rs[:,dm*i]
            Rs[:,dm*i] = tmp
        return [HM([[(Rs[i,:]).list() for i in range(sz)] for j in range(3)]).transpose(), HM([[((Rs.transpose())[i,:]).list() for i in range(sz)] for j in range(3)]).transpose()]
    else:
        print 'Dimension mismatch !!'  

def channel_product(A,B):
    P0 = HM(A.n(0),A.n(1),'zero')
    for i in range(A.n(0)):
        for j in range(A.n(1)):
            P0[i,j]=A[i,j,0]
    P1 = HM(A.n(0),A.n(1),'zero')
    for i in range(A.n(0)):
        for j in range(A.n(1)):
            P1[i,j]=A[i,j,1]
    P2 = HM(A.n(0),A.n(1),'zero')
    for i in range(A.n(0)):
        for j in range(A.n(1)):
            P2[i,j]=A[i,j,2]
    Q0 = HM(B.n(0),B.n(1),'zero')
    for i in range(A.n(0)):
        for j in range(A.n(1)):
            Q0[i,j]=B[i,j,0]
    Q1 = HM(B.n(0),B.n(1),'zero')
    for i in range(A.n(0)):
        for j in range(A.n(1)):
            Q1[i,j]=B[i,j,1]
    Q2 = HM(B.n(0),B.n(1),'zero')
    for i in range(A.n(0)):
        for j in range(A.n(1)):
            Q2[i,j]=B[i,j,2]
    R0 = Matrix(SR,P0.listHM())*Matrix(SR,Q0.listHM())
    R1 = Matrix(SR,P1.listHM())*Matrix(SR,Q1.listHM())
    R2 = Matrix(SR,P2.listHM())*Matrix(SR,Q2.listHM())
    return HM([[(R0[i,:]).list() for i in range(A.n(1))],[(R1[i,:]).list() for i in range(A.n(1))],[(R2[i,:]).list() for i in range(A.n(1))]]).transpose()

# first order Lagrange interpolation
def lagrange1t(u0, m0):
    var('x0')
    f = 1
    for m in range(m0):
        if(m!=u0):
            f = f*(exp(I*2*pi*(x0/m0+0/2)) + exp(I*2*pi*(m/m0+1/2)))/(exp(I*2*pi*(u0/m0+0/2))+ exp(I*2*pi*(m/m0+1/2)))
    return f

# Secom2 order Lagrange interpolation
def lagrange2t(u0, m0, u1, m1):
    var('x0,x1')
    f = 1
    for m in range(m0):
        for n in range(m1):
            if(m!=u0 or n!=u1):
                 f = f*(exp(I*2*pi*((x0/m0)^2+(n/m1)+0/3)) + exp(I*2*pi*((m/m0)^2+(x1/m1)+1/3)) + exp(I*2*pi*((m/m0)^2+(n/m1)+2/3)))/(exp(I*2*pi*((u0/m0)^2+(n/m1)+0/3))+ exp(I*2*pi*((m/m0)^2+(u1/m1)+1/3)) + exp(I*2*pi*((m/m0)^2+(n/m1)+2/3)))
    return f

# Third order Lagrange interpolation
def lagrange3t(u0, m0, u1, m1, u2, m2):
    var('x0,x1,x2')
    f = 1
    for m in range(m0):
        for n in range(m1):
            for p in range(m2):
                if(m!=u0 or n!=u1 or p!=u2):
                    f = f*((exp(I*2*pi*((x0/m0)^3+(n/m1)^2+(p/m2)^1+0/4))+exp(I*2*pi*((m/m0)^3+(x1/m1)^2+(p/m2)^1+1/4))+exp(I*2*pi*((m/m0)^3+(n/m1)^2+(x2/m2)^1+2/4))+exp(I*2*pi*((m/m0)^3+(n/m1)^2+(p/m2)^1+3/4))))/((exp(I*2*pi*((u0/m0)^3+(n/m1)^2+(p/m2)^1+0/4))+exp(I*2*pi*((m/m0)^3+(u1/m1)^2+(p/m2)^1+1/4))+exp(I*2*pi*((m/m0)^3+(n/m1)^2+(u2/m2)^1+2/4))+exp(I*2*pi*((m/m0)^3+(n/m1)^2+(p/m2)^1+3/4))))
    return f

# Fourth order Lagrange interpolation
def lagrange4t(u0, m0, u1, m1, u2, m2, u3, m3):
    var('x0,x1,x2,x3')
    f = 1
    for m in range(m0):
        for n in range(m1):
            for p in range(m2):
                for q in range(m3):
                    if(m!=u0 or n!=u1 or p!=u2 or q!=u3):
                        f = f*(exp(I*2*pi*((x0/m0)^4+(n/m1)^3+(p/m2)^2+(q/m3)^1+0/5))+exp(I*2*pi*((m/m0)^4+(x1/m1)^3+(p/m2)^2+(q/m3)^1+1/5))+exp(I*2*pi*((m/m0)^4+(n/m1)^3+(x2/m2)^2+(q/m3)^1+2/5))+exp(I*2*pi*((m/m0)^4+(n/m1)^3+(p/m2)^2+(x3/m3)^1+3/5))+ exp(I*2*pi*((m/m0)^4+(n/m1)^3+(p/m2)^2+(q/m3)^1+4/5)))/(exp(I*2*pi*((u0/m0)^4+(n/m1)^3+(p/m2)^2+(q/m3)^1+0/5))+exp(I*2*pi*((m/m0)^4+(u1/m1)^3+(p/m2)^2+(q/m3)^1+1/5))+exp(I*2*pi*((m/m0)^4+(n/m1)^3+(u2/m2)^2+(q/m3)^1+2/5))+exp(I*2*pi*((m/m0)^4+(n/m1)^3+(p/m2)^2+(u3/m3)^1+3/5))+ exp(I*2*pi*((m/m0)^4+(n/m1)^3+(p/m2)^2+(q/m3)^1+4/5)))
    return f

# Fourth order Lagrange interpolation
def lagrange5t(u0, m0, u1, m1, u2, m2, u3, m3, u4, m4):
    var('x0,x1,x2,x3,x4')
    f = 1
    for m in range(m0):
        for n in range(m1):
            for p in range(m2):
                for q in range(m3):
                    for r in range(m4):
                        if(m!=u0 or n!=u1 or p!=u2 or q!=u3 or r!=u4):
                            f = f*(exp(I*2*pi*((x0/m0)^5+(n/m1)^4+(p/m2)^3+(q/m3)^2+(r/m4)^1+0/6))+exp(I*2*pi*((m/m0)^5+(x1/m1)^4+(p/m2)^3+(q/m3)^2+(r/m4)^1+1/6))+exp(I*2*pi*((m/m0)^5+(n/m1)^4+(x2/m2)^3+(q/m3)^2+(r/m4)^1+2/6))+exp(I*2*pi*((m/m0)^5+(n/m1)^4+(p/m2)^3+(x3/m3)^2+(r/m4)^1+3/6))+exp(I*2*pi*((m/m0)^5+(n/m1)^4+(p/m2)^3+(q/m3)^2+(x4/m4)^1+4/6))+exp(I*2*pi*((m/m0)^5+(n/m1)^4+(p/m2)^3+(q/m3)^2+(r/m4)^1+5/6)))/(exp(I*2*pi*((u0/m0)^5+(n/m1)^4+(p/m2)^3+(q/m3)^2+(r/m4)^1+0/6))+exp(I*2*pi*((m/m0)^5+(u1/m1)^4+(p/m2)^3+(q/m3)^2+(r/m4)^1+1/6))+exp(I*2*pi*((m/m0)^5+(n/m1)^4+(u2/m2)^3+(q/m3)^2+(r/m4)^1+2/6))+exp(I*2*pi*((m/m0)^5+(n/m1)^4+(p/m2)^3+(u3/m3)^2+(r/m4)+3/6))+exp(I*2*pi*((m/m0)^5+(n/m1)^4+(p/m2)^3+(q/m3)^2+(u4/m4)+4/6))+exp(I*2*pi*((m/m0)^5+(n/m1)^4+(p/m2)^3+(q/m3)^2+(r/m4)+5/6)))
    return f

# Matrix0 Lagrange interpolation
def matrix_lagrange_polynomial(A,order):
    f = -1
    if order == 1:
        # Initialix2ation of the polx1nomial
        f = 0
        for u0 in range(len(A)):
            f = f + A[u0]*lagrange1t(u0, Integer(len(A)))

    elif order == 2:
        # Initialix2ation of the polx1nomial
        f = 0
        for u0 in range(len(A)):
            for u1 in range(len(A[0])):
                f = f + A[u0][u1]*lagrange2t(u0, Integer(len(A)), u1, Integer(len(A[0])))

    elif order == 3 :
        # Initialix2ation of the polx1nomial
        f = 0
        for u0 in range(len(A)):
            for u1 in range(len(A[0])):
                for u2 in range(len(A[0][0])):
                    f = f + A[u0][u1][u2]*lagrange3t(u0, Integer(len(A)), u1, Integer(len(A[0])), u2, Integer(len(A[0][0])))

    elif order == 4 :
        # Initialix2ation of the polx1nomial
        f = 0
        for u0 in range(len(A)):
            for u1 in range(len(A[0])):
                for u2 in range(len(A[0][0])):
                    for u3 in range(len(A[0][0][0])):
                        f = f + A[u0][u1][u2][u3]*lagrange4t(u0, Integer(len(A)), u1, Integer(len(A[0])), u2, Integer(len(A[0][0])), u3, Integer(len(A[0][0][0])))
 
    elif order == 5 :
        # Initialix2ation of the polx1nomial
        f = 0
        for u0 in range(len(A)):
            for u1 in range(len(A[0])):
                for u2 in range(len(A[0][0])):
                    for u3 in range(len(A[0][0][0])):
                        for u4 in range(len(A[0][0][0][0])):
                            f = f + A[u0][u1][u2][u3][u4]*lagrange5t(u0, Integer(len(A)), u1, Integer(len(A[0])), u2, Integer(len(A[0][0])), u3, Integer(len(A[0][0][0])), u4, Integer(len(A[0][0][0][0])))
 
    else:
        print 'The interpolation for order', order,'-tensors is not implemented'
 
    return f



def Deter(A):
    """
    Computes symbolically the determinant of a square matrix
    using the sum over permutation formula.

    EXAMPLES:
    ::
        sage: M = Matrix(SR, MatrixGenerate(2, 2, 'm')); det(M)
        -m01*m10 + m00*m11

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the permutations
    P = Permutations(range(A.nrows()))
    return sum([Permutation([p[i]+1 for i in range(len(p))]).signature()*\
prod([A[k][p[k]] for k in range(A.nrows())]) for p in P])

def conj3(z, k):
    """
    Computes the k-th conjugate corresponding to the input complex number

    EXAMPLES:
    ::
        sage: conj3(1+3*I,1)
        0.999999999999998 - 3.00000000000000*I
    AUTHORS:
    - Edinah K. Gnang
    """
    # Definition of the variables
    a, b = var('a, b')
    if 0 < CC(pi+z).arg() and CC(pi+z).arg() <= N(2*pi/3):
        S = solve([a+b*CC(exp(I*2*pi/3)).real()==CC(z).real(), b*CC(exp(I*2*pi/3)).imag()==CC(z).imag()], a, b)
        if Integer(mod(k,3))==0:
            return 1.0*S[0][0].rhs()+1.0*CC(exp(I*2*pi/3))*S[0][1].rhs()
        elif Integer(mod(k,3))==1:
            return 1.0*S[0][0].rhs()+1.0*CC(exp(I*4*pi/3))*S[0][1].rhs()
        elif Integer(mod(k,3))==2:
            return 1.0*S[0][0].rhs()+1.0*S[0][1].rhs()
    elif N(2*pi/3) < CC(pi+z).arg() and CC(pi+z).arg() <= N(4*pi/3):
        S = solve([a*CC(exp(I*2*pi/3)).real()+b*CC(exp(I*4*pi/3)).real()==CC(z).real(), b*CC(exp(I*2*pi/3)).imag()+a*CC(exp(I*4*pi/3)).imag()==CC(z).imag()], a, b)
        if Integer(mod(k,3))==0:
            return 1.0*CC(exp(I*2*pi/3))*S[0][0].rhs()+1.0*CC(exp(I*4*pi/3))*S[0][1].rhs()
        elif Integer(mod(k,3))==1:
            return 1.0*CC(exp(I*4*pi/3))*S[0][0].rhs()+1.0*CC(exp(I*2*pi/3))*S[0][1].rhs()
        elif Integer(mod(k,3))==2:
            return 1.0*S[0][0].rhs()+1.0*S[0][1].rhs()
    else :
        S = solve([a+b*CC(exp(I*4*pi/3)).real()==CC(z).real(), b*CC(exp(I*4*pi/3)).imag()==CC(z).imag()], a, b)
        if Integer(mod(k,3))==0:
            return 1.0*S[0][0].rhs()+1.0*CC(exp(I*4*pi/3))*S[0][1].rhs()
        elif Integer(mod(k,3))==1:
            return 1.0*S[0][0].rhs()+1.0*CC(exp(I*2*pi/3))*S[0][1].rhs()
        elif Integer(mod(k,3))==2:
            return 1.0*S[0][0].rhs()+1.0*S[0][1].rhs()
  
def MeanApproximation(Im):
    import pylab, numpy
    from scipy import misc
    # Initialization of the Slices.
    A = HM(len(Im), 1         , len(Im[0][0]), 'zero')
    B = HM(len(Im), len(Im[0]), 1            , 'zero')
    Ch = HM(1      , len(Im[0]), len(Im[0][0]), 'zero')
    # Initialization of the Hypermatrix associated with the image.
    T = HM(len(Im),len(Im[0]),len(Im[0][0]),'zero')
    # Filling up the image Hypermatrix
    for i in range(T.nrows()):
        for j in range(T.ncols()):
            for k in range(T.ndpts()):
                T[i,j,k] = Im[i,j,k]
    # Computing the mean of row depth slice
    for u in range(A.nrows()):
        for v in range(A.ndpts()):
            A[u,0,v] = mean([T[u,i,v] for i in range(T.ncols())])
    # Computing the mean row column slice
    for u in range(B.nrows()):
        for v in range(B.ncols()):
            B[u,v,0] = mean([T[u,v,i] for i in range(T.ndpts())])
    # Computing the mean column depth slice
    for u in range(Ch.ncols()):
        for v in range(Ch.ndpts()):
            Ch[0,u,v] = mean([T[i,u,v] for i in range(T.nrows())])
    # Computing the outer-product of the mean slices.
    return [T, A, B, Ch]

def MedianApproximation(Im):
    import pylab, numpy
    from scipy import misc
    # Initialization of the Slices.
    A = HM(len(Im), 1         , len(Im[0][0]), 'zero')
    B = HM(len(Im), len(Im[0]), 1            , 'zero')
    Ch = HM(1      , len(Im[0]), len(Im[0][0]), 'zero')
    # Initialization of the Hypermatrix associated with the image.
    T = HM(len(Im),len(Im[0]),len(Im[0][0]),'zero')
    # Filling up the image Hypermatrix
    for i in range(T.nrows()):
        for j in range(T.ncols()):
            for k in range(T.ndpts()):
                T[i,j,k] = Im[i,j,k]
    # Computing the mean of row depth slice
    for u in range(A.nrows()):
        for v in range(A.ndpts()):
            A[u,0,v] = median([T[u,i,v] for i in range(T.ncols())])
    # Computing the mean row column slice
    for u in range(B.nrows()):
        for v in range(B.ncols()):
            B[u,v,0] = median([T[u,v,i] for i in range(T.ndpts())])
    # Computing the mean column depth slice
    for u in range(Ch.ncols()):
        for v in range(Ch.ndpts()):
            Ch[0,u,v] = median([T[i,u,v] for i in range(T.nrows())])
    # Computing the outer-product of the mean slices.
    return [T, A, B, Ch]


def GrdDcnt(T, A, B, Ch, p=10, nb_stp=100, stp_sz=0.250):    
    import pylab, numpy
    from scipy import misc
    # Initialization of the variable slices.
    X  = HM(T.n(0), 1     , T.n(2), 'x'); Rx = A
    Y  = HM(T.n(0), T.n(1), 1     , 'y'); Ry = B
    Z  = HM(1     , T.n(1), T.n(2), 'z'); Rz = Ch
    
    # Computing the p norm associated
    f = sum([((T-X*(Y,Z)).elementwise_exponent(p))[i,j,k] for i in range(T.n(0)) for j in range(T.n(1)) for k in range(T.n(2))])
    
    # Computing symbolicaly the gradient vector
    Gv = []
    for i in range(X.n(0)):
        for j in range(X.n(1)):
            for k in range(X.n(2)):
                Gv.append(f.diff(X[i,j,k]))
    for i in range(Y.n(0)):
        for j in range(Y.n(1)):
            for k in range(Y.n(2)):
                Gv.append(f.diff(Y[i,j,k]))
    for i in range(Z.n(0)):
        for j in range(Z.n(1)):
            for k in range(Z.n(2)):
                Gv.append(f.diff(Z[i,j,k]))
    # Performing the Gradient descent
    for stp in range(nb_stp):
        # Updating the Slice Rx
        indx = 0
        for i in range(X.n(0)):
            for j in range(X.n(1)):
                for k in range(X.n(2)):
                    Rx[i,j,k] = N(A[i,j,k]-stp_sz*(Gv[indx]/sum([Gv[s]^p for s in range(len(Gv))])^(1/p)).subs(dict([(X[u,v,w],A[u,v,w]) for u in range(X.n(0)) for v in range(X.n(1)) for w in range(X.n(2))]+[(Y[u,v,w],B[u,v,w]) for u in range(Y.n(0)) for v in range(Y.n(1)) for w in range(Y.n(2))]+[(Z[u,v,w],Ch[u,v,w]) for u in range(Z.n(0)) for v in range(Z.n(1)) for w in range(Z.n(2))])))
                    indx = indx+1
        # Updating the Slice Ry
        for i in range(Y.n(0)):
            for j in range(Y.n(1)):
                for k in range(Y.n(2)):
                    Ry[i,j,k] = N(B[i,j,k]-stp_sz*(Gv[indx]/sum([Gv[s]^p for s in range(len(Gv))])^(1/p)).subs(dict([(X[u,v,w],A[u,v,w]) for u in range(X.n(0)) for v in range(X.n(1)) for w in range(X.n(2))]+[(Y[u,v,w],B[u,v,w]) for u in range(Y.n(0)) for v in range(Y.n(1)) for w in range(Y.n(2))]+[(Z[u,v,w],Ch[u,v,w]) for u in range(Z.n(0)) for v in range(Z.n(1)) for w in range(Z.n(2))])))
                    indx = indx+1
        # Updating the Slice Rz
        for i in range(Z.n(0)):
            for j in range(Z.n(1)):
                for k in range(Z.n(2)):
                    Rz[i,j,k] = N(Ch[i,j,k]-stp_sz*(Gv[indx]/sum([Gv[s]^p for s in range(len(Gv))])^(1/p)).subs(dict([(X[u,v,w],A[u,v,w]) for u in range(X.n(0)) for v in range(X.n(1)) for w in range(X.n(2))]+[(Y[u,v,w],B[u,v,w]) for u in range(Y.n(0)) for v in range(Y.n(1)) for w in range(Y.n(2))]+[(Z[u,v,w],Ch[u,v,w]) for u in range(Z.n(0)) for v in range(Z.n(1)) for w in range(Z.n(2))])))
                    indx = indx+1
        # Updating the gradient points
        A = Rx; B = Ry; Ch = Rz
        # Printing the current error
        print "At iteration",stp," the error is ",N(f.subs(dict([(X[u,v,w],A[u,v,w]) for u in range(X.n(0)) for v in range(X.n(1)) for w in range(X.n(2))]+[(Y[u,v,w],B[u,v,w]) for u in range(Y.n(0)) for v in range(Y.n(1)) for w in range(Y.n(2))]+[(Z[u,v,w],Ch[u,v,w]) for u in range(Z.n(0)) for v in range(Z.n(1)) for w in range(Z.n(2))])))
    return [A,B,Ch]

def ZeroPadding(A):
    # Initializing the size parameter
    sz = max([A.n(i) for i in range(A.order())])
    # Initializing the Hypermatrix
    T = HM(sz, sz, sz, 'zero')
    # Filling up the Hypermatrix
    for r in A.n(0):
        for c in A.n(1):
            for d in A.n(2):
                T[r,c,d]=A[r,c,d]
    return T

def GenerateUnitLpNormVector(n,p = 2,indx=0):
    if n == 1:
        return [1]
    else :
        X = []
        X.append(cos(var('t'+str(indx)))^(2/p))
        for i in range(1,n-1):
            X.append((prod([sin(var('t'+str(j+indx))) for j in range(i)])*cos(var('t'+str(i+indx))))^(2/p))
        X.append((prod([sin(var('t'+str(j+indx))) for j in range(n-1)]))^(2/p))
        return X

def ProbabilityMatrix(n, xi=0):
    # Initializing the matrix to be filled
    M = Matrix(SR, zero_matrix(n,n))
    # Initialixzing the variable index
    indx=xi
    for c in range(n-1):
        # Initializing the probability vector associated with the row c
        La = GenerateUnitLpNormVector(n-c,1,indx)
        # Updating the variable index
        indx = indx+len(La)-1
        # Initializing the probability vector associated with the column c
        Lb = GenerateUnitLpNormVector(n-c-1,1,indx)
        # Updating the variable index
        indx = indx+len(Lb)-1
        # Loop which fills up the Matrix
        for i in range(c, c+len(La)):
            if c > 0:
                # Filling up the row c of the Matrix M
                M[c,i] = (1-sum([M[c,j] for j in range(c)]))*La[i-c]
                if i > c:
                    # Filling up the column c of the Matrix M
                    M[i,c] = (1-sum([M[j,c] for j in range(c+1)]))*Lb[i-c-1]
            else:
                # Filling up the row c of the Matrix M
                M[c,i] = La[i-c]
                if i > c:
                    # Filling up the column c of the Matrix M
                    M[i,c] = (1-sum([M[j,c] for j in range(c+1)]))*Lb[i-c-1]
    M[n-1,n-1]=1-sum([M[j,n-1] for j in range(n-1)])
    return M

def ProbabilitySymMatrix(n, xi=0):
    # Initializing the matrix to be filled
    M = Matrix(SR, zero_matrix(n,n))
    # Initialixzing the variable index
    indx=xi
    for c in range(n-1):
        # Initializing the probability vector associated with the row c
        La = GenerateUnitLpNormVector(n-c,1,indx)
        # Updating the variable index
        indx = indx+len(La)-1
        # Loop which fills up the Matrix
        for i in range(c, c+len(La)):
            if c > 0:
                # Filling up the row c of the Matrix M
                M[c,i] = (1-sum([M[c,j] for j in range(c)]))*La[i-c]
                if i > c:
                    # Filling up the column c of the Matrix M
                    M[i,c] = M[c,i] 
            else:
                # Filling up the row c of the Matrix M
                M[c,i] = La[i-c]
                if i > c:
                    # Filling up the column c of the Matrix M
                    M[i,c] = M[c,i]
    M[n-1,n-1]=1-sum([M[j,n-1] for j in range(n-1)])
    return M

def HypermatrixPseudoInversePairsII(A,B):
    """
     Outputs the pseudo inverse pairs associated with the input pairs of matrices

    EXAMPLES:
    ::
        sage: A1=HM([[[0.1631135370902057,0.11600112072013125],[0.9823708115400902,0.39605960486710756]]\
,[[0.061860929755424676,0.2325542810173995],[0.39111210957450926,0.2019809359102137]]])
        sage: A2=HM([[[0.15508921433883183,0.17820377184410963],[0.48648171594508205,0.01568017636082064]]\
,[[0.8250247759993575,0.1938307874191597],[0.23867299119274843,0.3935578730402869]]])
        sage: [B1,B2]=HypermatrixPseudoInversePairsII(A1,A2)

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    sz = len(A.listHM())

    # Initializing the list of linear constraints
    CnstrLst = []

    # Initilizing the variable list
    Vrbls  = [var('ln_al'+str(i)+str(j)+str(k)) \
for i in range(sz) for j in range(sz) for k in range(sz)]+\
[var('ln_bt'+str(i)+str(j)+str(k)) for i in range(sz) for j in range(sz) \
for k in range(sz)]

    for m in range(sz):
        for p in range(sz):
            for n in range(sz):
                V=Matrix(CC, sz, sz, [(A[m,k1,k0])*(B[k0,k1,p]) \
for k0 in range(sz) for k1 in range(sz)]).inverse()
                CnstrLst=CnstrLst+[\
var('ln_al'+str(m)+str(n)+str(k1))+var('ln_bt'+str(k1)+str(n)+str(p))==\
ln(V[k1,n])  for k1 in range(sz)]
    [A,b]=ConstraintFormator(CnstrLst,Vrbls)

    # Importing the Numerical Python package
    # for computing the matrix pseudo inverse
    import numpy
    sln = matrix(numpy.linalg.pinv(A))*b
    R1 = HM(sz,sz,sz,'zero')
    for i in range(sz):
        for j in range(sz):
            for k in range(sz):
                R1[i,j,k] = exp(sln[i*sz^2+j*sz^1+k*sz^0,0])
    R2 = HM(sz, sz, sz,'zero')
    for i in range(sz):
        for j in range(sz):
            for k in range(sz):
                R2[i,j,k] = exp(sln[sz^3+i*sz^2+j*sz^1+k*sz^0,0])
    return [R1,R2]

def HypermatrixPseudoInversePairsIII(A,B,p=2,nb_stp=100,stp_sz=0.25):
    # Defining the variables
    #X=HM(2,2,2,'x');Rx=(HM(HypermatrixPermutation([0,1,2])));Mx=Rx 
    #Y=HM(2,2,2,'y');Ry=(HM(HypermatrixPermutation([0,1,2]))).transpose();My=Ry
    X=HM(2,2,2,'x');Rx=A
    Y=HM(2,2,2,'y');Ry=B

    # Computing the associated p norm
    T=HM(2,2,2,'one')-((HM(2,2,2,'one'))*(A,B))*(X,Y)
    f=sum([((HM(2,2,2,'one')-((HM(2,2,2,'one'))*(A,B))*(X,Y)).elementwise_exponent(p))[i,j,k] for i in range(T.n(0)) for j in range(T.n(1)) for k in range(T.n(2))])
    
    # Computing symbolicaly the gradient vector
    Gv = []
    for i in range(X.n(0)):
        for j in range(X.n(1)):
            for k in range(X.n(2)):
                Gv.append(f.diff(X[i,j,k]))
    for i in range(Y.n(0)):
        for j in range(Y.n(1)):
            for k in range(Y.n(2)):
                Gv.append(f.diff(Y[i,j,k]))
    # Performing the Gradient descent
    for stp in range(nb_stp):
        # Updating the Slice Rx
        indx = 0
        for i in range(X.n(0)):
            for j in range(X.n(1)):
                for k in range(X.n(2)):
                    #Rx[i,j,k] = N(Mx[i,j,k]-stp_sz*(Gv[indx]/sum([Gv[s]^p for s in range(len(Gv))])^(1/p)).subs(dict([(X[u,v,w],Mx[u,v,w]) for u in range(X.n(0)) for v in range(X.n(1)) for w in range(X.n(2))]+[(Y[u,v,w],My[u,v,w]) for u in range(Y.n(0)) for v in range(Y.n(1)) for w in range(Y.n(2))])))
                    Rx[i,j,k] = N(Mx[i,j,k]-stp_sz*Gv[indx].subs(dict([(X[u,v,w],Mx[u,v,w]) for u in range(X.n(0)) for v in range(X.n(1)) for w in range(X.n(2))]+[(Y[u,v,w],My[u,v,w]) for u in range(Y.n(0)) for v in range(Y.n(1)) for w in range(Y.n(2))]))/sum([(abs(Gv[l].subs(dict([(X[u,v,w],Mx[u,v,w]) for u in range(X.n(0)) for v in range(X.n(1)) for w in range(X.n(2))]+[(Y[u,v,w],My[u,v,w]) for u in range(Y.n(0)) for v in range(Y.n(1)) for w in range(Y.n(2))]))))^p for l in range(len(Gv))])^(1/p))
                    indx = indx+1
        # Updating the Slice Ry
        for i in range(Y.n(0)):
            for j in range(Y.n(1)):
                for k in range(Y.n(2)):
                    #Ry[i,j,k] = N(My[i,j,k]-stp_sz*(Gv[indx]/sum([Gv[s]^p for s in range(len(Gv))])^(1/p)).subs(dict([(X[u,v,w],Mx[u,v,w]) for u in range(X.n(0)) for v in range(X.n(1)) for w in range(X.n(2))]+[(Y[u,v,w],My[u,v,w]) for u in range(Y.n(0)) for v in range(Y.n(1)) for w in range(Y.n(2))])))
                    Ry[i,j,k] = N(My[i,j,k]-stp_sz*(Gv[indx]).subs(dict([(X[u,v,w],Mx[u,v,w]) for u in range(X.n(0)) for v in range(X.n(1)) for w in range(X.n(2))]+[(Y[u,v,w],My[u,v,w]) for u in range(Y.n(0)) for v in range(Y.n(1)) for w in range(Y.n(2))])))
                    Ry[i,j,k] = N(My[i,j,k]-stp_sz*Gv[indx].subs(dict([(X[u,v,w],Mx[u,v,w]) for u in range(X.n(0)) for v in range(X.n(1)) for w in range(X.n(2))]+[(Y[u,v,w],My[u,v,w]) for u in range(Y.n(0)) for v in range(Y.n(1)) for w in range(Y.n(2))]))/sum([(abs(Gv[l].subs(dict([(X[u,v,w],Mx[u,v,w]) for u in range(X.n(0)) for v in range(X.n(1)) for w in range(X.n(2))]+[(Y[u,v,w],My[u,v,w]) for u in range(Y.n(0)) for v in range(Y.n(1)) for w in range(Y.n(2))]))))^p for l in range(len(Gv))])^(1/p))
                    indx = indx+1
        # Updating the gradient points
        Mx = Rx; My = Ry
        # Printing the current error
        print "At iteration",stp," the error is ",N(f.subs(dict([(X[u,v,w],Mx[u,v,w]) for u in range(X.n(0)) for v in range(X.n(1)) for w in range(X.n(2))]+[(Y[u,v,w],My[u,v,w]) for u in range(Y.n(0)) for v in range(Y.n(1)) for w in range(Y.n(2))])))
    return [Mx,My]

def HypermatrixPseudoInversePairsUnsplit(A,B):
    """
     Outputs the pseudo inverse pairs associated with the input pairs of matrices

    EXAMPLES:
    ::
        sage: A1=HM([[[0.1631135370902057,0.11600112072013125],[0.9823708115400902,0.39605960486710756]]\
,[[0.061860929755424676,0.2325542810173995],[0.39111210957450926,0.2019809359102137]]])
        sage: A2=HM([[[0.15508921433883183,0.17820377184410963],[0.48648171594508205,0.01568017636082064]]\
,[[0.8250247759993575,0.1938307874191597],[0.23867299119274843,0.3935578730402869]]])
        sage: B1B2=HypermatrixPseudoInversePairsUnsplit(A1,A2)

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initializing the size
    sz = A.nrows()
    # Initializing the symbolic matrix
    T = HM(sz,sz,sz,'one')
    # We adopt the entry convention XY_mnjp=X_mnj*Y_jnp
    XY = HM(sz,sz,sz,sz, 'xy')
    # Initialization of the variables list
    Vrbls = XY.list()
    Tmp = T*(A,B)
    # Initializing the list of linear constraints
    CnstrLst = [sum([Tmp[m,j,p]*XY[m,n,j,p] for j in range(sz)]) == T[m,n,p] for m in range(sz) for n in range(sz) for p in range(sz)]
    # Formatting the constraints
    [M,b] = ConstraintFormator(CnstrLst,Vrbls)
    # Importing the Numerical Python package for computing the matrix pseudo inverse
    import numpy; sln = matrix(numpy.linalg.pinv(M))*b
    # Filling up the hypermartix XY
    # Initialization of the list specifying the dimensions of the output
    l = [XY.n(i) for i in range(XY.order())]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        XY[tuple(entry)]=sln[i][0]
    return XY

def HypermatrixPseudoInversePairAction(T, A, B):
    # Computing the unsplit Inverse pairs
    XY = HypermatrixPseudoInversePairsUnsplit(A,B)
    Rs = HM(A.nrows(), A.ncols(), A.ndpts(),'zero')
    for m in range(Rs.nrows()):
        for n in range(Rs.ncols()):
            for p in range(Rs.ndpts()):
                Rs[m,n,p] = sum([T[m,j,p]*XY[m,n,j,p] for j in range(A.nrows())])
    return Rs 

def GenerateRandomHypermatrix(*l):
    if prod(list(l)) != 0:
        # Initializing the input for generating a symbolic hypermatrix
        inpts = list(l)+['zero']
        # Initialization of the hypermatrix
        Rh = HM(*inpts)
        # Main loop performing the transposition of the entries
        for i in range(prod(l)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entry = [mod(i,l[0])]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            Rh[tuple(entry)]=random()
        return Rh
    else :
        raise ValueError, "The Dimensions must all be non-zero."

def nearest_orthogonal2x2x2(A, p=2, stp_sz=0.25, nb_stp=25):
    # Orthogonal matrix parametrization
    r1,r2,r3,r4,r5,r6=var('r1,r2,r3,r4,r5,r6')
    # Complete parametric description of orthogonal hypermatrix
    Q = HM([[[e^(-r1 + r3 + r6)/(e^(-3*r1 + 3*r3 + 3*r6)+e^(3*r6))^(1/3), e^r4],[e^r6/(e^(-3*r1+3*r3+3*r6)+e^(3*r6))^(1/3), e^r2]],[[-e^(r1+r2-r3-r4+r5),e^r3/(e^(3*r1)+e^(3*r3))^(1/3)],[e^r5, e^r1/(e^(3*r1)+e^(3*r3))^(1/3)]]])
    # Initializing the values in the orthogonal hypermatrix parametrization
    Rv = [0, 0, 0, 0, 0, 0]
    # Initializing the function
    f = sum([((A-Q).elementwise_exponent(p))[i,j,k] for i in range(2) for j in range(2) for k in range(2)])
    # Initializing the gradient descent vector
    Gf = []
    for i in range(1,7):
        Gf.append(f.diff(var('r'+str(i))))
        #print "Gf["+str(i-1)+"] = ",Gf[i-1]
    # Updating the orthogonality parametrization
    for stp in range(nb_stp):
        Tmp = copy(Rv)
        Rv  = [N(Tmp[i]-stp_sz*(Gf[i].subs(dict([(var('r'+str(j+1)),Tmp[j]) for j in range(6)]))/sum([abs(Gf[s].subs(dict([(var('r'+str(k+1)),Tmp[k]) for k in range(6)])))^p for s in range(len(Gf))])^(1/p))) for i in range(6)]
        #print "Rv = ",Rv
    # Outputing of the final computation        
    return Q.subs(dict([(var('r'+str(j+1)), Rv[j]) for j in range(6)]))

def nearest_orthogonal2x2x2II(A, D, p=2, stp_sz=0.25, nb_stp=25):
    # Orthogonal matrix parametrization
    r1,r2,r3,r4,r5,r6=var('r1,r2,r3,r4,r5,r6')
    # Complete parametric description of orthogonal hypermatrix
    Q=HM([[[e^(-r1 + r3 + r6)/(e^(-3*r1 + 3*r3 + 3*r6)+e^(3*r6))^(1/3), e^r4],[e^r6/(e^(-3*r1+3*r3+3*r6)+e^(3*r6))^(1/3), e^r2]],[[-e^(r1+r2-r3-r4+r5),e^r3/(e^(3*r1)+e^(3*r3))^(1/3)],[e^r5, e^r1/(e^(3*r1)+e^(3*r3))^(1/3)]]])
    # Initializing the values in the orthogonal hypermatrix parametrization
    Rv = [0, 0, 0, 0, 0, 0]
    # Initializing the function
    Qs=Q*(D,D.transpose())
    f=sum([((A-Qs*(Qs.transpose(2),Qs.transpose())).elementwise_exponent(p))[i,j,k] for i in range(2) for j in range(2) for k in range(2)])
    # Initializing the gradient descent vector
    Gf=[]
    for i in range(1,7):
        Gf.append(f.diff(var('r'+str(i))))
    # Updating the orthogonality parametrization
    for stp in range(nb_stp):
        Tmp = copy(Rv)
        Rv  = [N(Tmp[i]-stp_sz*(Gf[i].subs(dict([(var('r'+str(j+1)),Tmp[j]) for j in range(6)]))/sum([abs(Gf[s].subs(dict([(var('r'+str(k+1)),Tmp[k]) for k in range(6)])))^p for s in range(len(Gf))])^(1/p))) for i in range(6)]
    # Outputing of the final computation        
    return Q.subs(dict([(var('r'+str(j+1)), Rv[j]) for j in range(6)]))

def best_diagonal_fit(A, Q, p=2, stp_sz=0.25, nb_stp=50):
    # Initializing the symbolic diagonal hypermatrix
    D = HM(Matrix(SR, SymMatrixGenerate(2,'w')))
    # Initializing the values in the orthogonal hypermatrix parametrization
    v00 = 1; v01 = 1; v11 = 1
    # Initializing the function
    #f = sum([((A-Q*(D,D.transpose())).elementwise_exponent(p))[i,j,k] for i in range(A.nrows()) for j in range(A.ncols()) for k in range(A.ndpts())])
    f = sum([((  (Q*(D,D.transpose()))*((Q*(D,D.transpose())).transpose(2),(Q*(D,D.transpose())).transpose())-A*(A.transpose(2),A.transpose())  ).elementwise_exponent(p))[i,j,k] for i in range(A.nrows()) for j in range(A.ncols()) for k in range(A.ndpts())])
    #print f
    # Initializing the gradient descent vector
    Gf = [];Gf.append(f.diff(var(w00)));Gf.append(f.diff(var(w01)));Gf.append(f.diff(var(w11)))
    # Updating the diagonality fit
    for stp in range(nb_stp):
        v00 = N(v00-stp_sz*(Gf[0].subs(w00=v00,w01=v01,w11=v11)/sum([abs(Gf[i].subs(w00=v00,w01=v01,w11=v11))^p for i in range(3)])^(1/p)))
        #print "v00 = ",v00
        v01 = N(v01-stp_sz*(Gf[1].subs(w00=v00,w01=v01,w11=v11)/sum([abs(Gf[i].subs(w00=v00,w01=v01,w11=v11))^p for i in range(3)])^(1/p)))
        #print "v01 = ",v01
        v11 = N(v11-stp_sz*(Gf[2].subs(w00=v00,w01=v01,w11=v11)/sum([abs(Gf[i].subs(w00=v00,w01=v01,w11=v11))^p for i in range(3)])^(1/p)))
        #print "v11 = ",v11
    # Outputing of the final computation 
    return D.subs(dict([(w00,v00), (w01,v01), (w11,v11)]))

def best_diagonal_fitII(A, Q, p=2, stp_sz=0.25, nb_stp=50):
    # Initializing the symbolic diagonal hypermatrix
    D = HM(Matrix(SR, SymMatrixGenerate(2,'w')))
    # Initializing the values in the orthogonal hypermatrix parametrization
    v00 = 1; v01 = 1; v11 = 1
    # Initializing the function
    f = sum([(((Q*(D,D.transpose()))*((Q*(D,D.transpose())).transpose(2),(Q*(D,D.transpose())).transpose())-A).elementwise_exponent(p))[i,j,k] for i in range(A.nrows()) for j in range(A.ncols()) for k in range(A.ndpts())])
    #print f
    # Initializing the gradient descent vector
    Gf = [];Gf.append(f.diff(var(w00)));Gf.append(f.diff(var(w01)));Gf.append(f.diff(var(w11)))
    # Updating the diagonality fit
    for stp in range(nb_stp):
        v00 = N(v00-stp_sz*(Gf[0].subs(w00=v00,w01=v01,w11=v11)/sum([abs(Gf[i].subs(w00=v00,w01=v01,w11=v11))^p for i in range(3)])^(1/p)))
        #print "v00 = ",v00
        v01 = N(v01-stp_sz*(Gf[1].subs(w00=v00,w01=v01,w11=v11)/sum([abs(Gf[i].subs(w00=v00,w01=v01,w11=v11))^p for i in range(3)])^(1/p)))
        #print "v01 = ",v01
        v11 = N(v11-stp_sz*(Gf[2].subs(w00=v00,w01=v01,w11=v11)/sum([abs(Gf[i].subs(w00=v00,w01=v01,w11=v11))^p for i in range(3)])^(1/p)))
        #print "v11 = ",v11
    # Outputing of the final computation 
    return D.subs(dict([(w00,v00), (w01,v01), (w11,v11)]))


def SpectralDecomposition(A,nb_itr):
    # Initial conditions
    Qtmp = A
    """
    First iteration for setting the initial benchmarks
    """
    # Initializing the hypermatrix
    D = HM(2,2,2,'zero');D[0,0,0]=1;D[1,0,0]=1;D[0,1,1]=1;D[1,1,1]=1
    # Updating Qtmp to the nearest orthogonal matrix
    Q = nearest_orthogonal2x2x2II(Qtmp, D, 2, 0.15, 100)
    # Initialization of the diagonal matrix
    D = best_diagonal_fitII(A, Q, 2, 0.10, 100)
    # Inverting the non-zero entries of D in Di
    Di=D
    if Di[0,0,0]!=0:
        Di[0,0,0]=1/Di[0,0,0]
    if Di[1,0,0]!=0:
        Di[1,0,0]=1/Di[1,0,0]
    if Di[0,1,1]!=0:
        Di[0,1,1]=1/Di[0,1,1]
    if Di[1,1,1]!=0:
        Di[1,1,1]=1/Di[1,1,1]
    # Updating the temporary orthogonal matrix
    Qs = Q*(D,D.transpose())
    print "At iteration "+str(0)+" the error is ",sum([abs(((A-Qs*(Qs.transpose(2),Qs.transpose())).elementwise_exponent(2))[i,j,k]) for i in range(A.nrows()) for j in range(A.ncols()) for k in range(A.ndpts())])
    err = sum([abs(((A-Qs*(Qs.transpose(2),Qs.transpose())).elementwise_exponent(2))[i,j,k]) for i in range(A.nrows()) for j in range(A.ncols()) for k in range(A.ndpts())])
    Qtmp=(HypermatrixPseudoInversePairAction(A,Qs.transpose(2),Qs.transpose()))*(Di,Di.transpose())
    # Variables for storing the final results.
    ResQ = Q
    ResD = D
    # Main Loop
    for itr in range(1,nb_itr):
        # Updating Qtmp to the nearest orthogonal matrix
        Q = nearest_orthogonal2x2x2II(Qtmp, 2, 0.15, 100)
        # Initialization of the diagonal matrix
        D = best_diagonal_fitII(A, Q, 6, 0.10, 100)
        # Inverting the non-zero entries of D in Di
        Di = D
        if Di[0,0,0]!=0:
            Di[0,0,0]=1/Di[0,0,0]
        if Di[1,0,0]!=0:
            Di[1,0,0]=1/Di[1,0,0]
        if Di[0,1,1]!=0:
            Di[0,1,1]=1/Di[0,1,1]
        if Di[1,1,1]!=0:
            Di[1,1,1]=1/Di[1,1,1]
        # Updating the temporary orthogonal matrix
        Qs = Q*(D,D.transpose())
        print "At iteration "+str(itr)+" the error is ",sum([abs(((A-Qs*(Qs.transpose(2),Qs.transpose())).elementwise_exponent(2))[i,j,k]) for i in range(A.nrows()) for j in range(A.ncols()) for k in range(A.ndpts())])
        # Checking if the current run has beaten the past established record
        if sum([abs(((A-Qs*(Qs.transpose(2),Qs.transpose())).elementwise_exponent(2))[i,j,k]) for i in range(A.nrows()) for j in range(A.ncols()) for k in range(A.ndpts())]) < err:
            err=sum([abs(((A-Qs*(Qs.transpose(2),Qs.transpose())).elementwise_exponent(2))[i,j,k]) for i in range(A.nrows()) for j in range(A.ncols()) for k in range(A.ndpts())])
            ResQ=Q
            ResD=D
        Qtmp=(HypermatrixPseudoInversePairAction(A,Qs.transpose(2),Qs.transpose()))*(Di,Di.transpose())
    return [ResQ,ResD]

def CoskewnessSpectralDecomposition(A, nb_itr):
    # Initializing the corresponding symmetric hypermatrix
    AAttAt = A*(A.transpose(2),A.transpose())
    print "\nA*(A.transpose(2),A.transpose()) = \n", AAttAt.listHM()
    # Initial conditions
    Qtmp = A
    """
    First iteration for setting the initial benchmarks
    """
    # Updating Qtmp to the nearest orthogonal matrix
    Q = nearest_orthogonal2x2x2(Qtmp, 2, 0.15, 100)
    # Initialization of the diagonal matrix
    D = best_diagonal_fit(A, Q, 2, 0.10, 100)
    # Inverting the non-zero entries of D in Di
    Di=D
    if Di[0,0,0]>10^(-20):
        Di[0,0,0]=1/Di[0,0,0]
    if Di[1,0,0]>10^(-20):
        Di[1,0,0]=1/Di[1,0,0]
    if Di[0,1,1]>10^(-20):
        Di[0,1,1]=1/Di[0,1,1]
    if Di[1,1,1]>10^(-20):
        Di[1,1,1]=1/Di[1,1,1]
    # Updating the temporary orthogonal matrix
    Qs = Q*(D, D.transpose())
    print "At iteration "+str(0)+" the error is ",sum([abs(((AAttAt-Qs*(Qs.transpose(2),Qs.transpose())).elementwise_exponent(2))[i,j,k]) for i in range(A.nrows()) for j in range(A.ncols()) for k in range(A.ndpts())])
    err = sum([abs(((AAttAt-Qs*(Qs.transpose(2),Qs.transpose())).elementwise_exponent(2))[i,j,k]) for i in range(A.nrows()) for j in range(A.ncols()) for k in range(A.ndpts())])
    Qtmp=(HypermatrixPseudoInversePairAction(AAttAt,Qs.transpose(2),Qs.transpose()))*(Di,Di.transpose())
    # Variables for storing the final results.
    ResQ = Q
    ResD = D
    # Main Loop
    for itr in range(1,nb_itr):
        # Updating Qtmp to the nearest orthogonal matrix
        #Q = nearest_orthogonal2x2x2(Qtmp, 2, 0.15, 100)
        Q = nearest_orthogonal2x2x2II(A, D, 2, 0.15, 100)
        # Initialization of the diagonal matrix
        D = best_diagonal_fit(A, Q, 6, 0.10, 100)
        # Inverting the non-zero entries of D in Di
        Di = D
        if Di[0,0,0]>10^(-20):
            Di[0,0,0]=1/Di[0,0,0]
        if Di[1,0,0]>10^(-20):
            Di[1,0,0]=1/Di[1,0,0]
        if Di[0,1,1]>10^(-20):
            Di[0,1,1]=1/Di[0,1,1]
        if Di[1,1,1]>10^(-20):
            Di[1,1,1]=1/Di[1,1,1]
        # Updating the temporary orthogonal matrix
        Qs = Q*(D,D.transpose())
        print "At iteration "+str(itr)+" the error is ",sum([abs(((AAttAt-Qs*(Qs.transpose(2),Qs.transpose())).elementwise_exponent(2))[i,j,k]) for i in range(A.nrows()) for j in range(A.ncols()) for k in range(A.ndpts())])
        # Checking if the current run has beaten the past established record
        if sum([abs(((AAttAt-Qs*(Qs.transpose(2),Qs.transpose())).elementwise_exponent(2))[i,j,k]) for i in range(A.nrows()) for j in range(A.ncols()) for k in range(A.ndpts())]) < err:
            err = sum([abs(((AAttAt-Qs*(Qs.transpose(2),Qs.transpose())).elementwise_exponent(2))[i,j,k]) for i in range(A.nrows()) for j in range(A.ncols()) for k in range(A.ndpts())])
            ResQ = Q
            ResD = D
        Qtmp=(HypermatrixPseudoInversePairAction(AAttAt,Qs.transpose(2),Qs.transpose()))*(Di,Di.transpose())
    return [ResQ, ResD]

def GeneralStochasticHypermatrix(t, od):
    """
    Generates an stochastic hypermatrix of the appropriate order
    for which each one of the dimensions are equal to 2.
    The vectors are not normalized.

    EXAMPLES:
    ::
        sage: Q = GeneralStochasticHypermatrix(var('t'), 2); Q
        [[cos(t)^2, sin(t)^2], [sin(t)^2, cos(t)^2]]

    AUTHORS:
    - Edinah K. Gnang, Ori Parzanchevski
    """
    if od < 2 and type(od)==type(1):
        raise ValueError, "The order must be greater or equal to 2."
    elif od == 2:
        return HM([[cos(t)^2, sin(t)^2], [sin(t)^2, cos(t)^2]])
    elif od > 2 and type(od) == type(1):
        dms = [2,2]
        B = GeneralStochasticHypermatrix(t,2)
        for z in range(2,od):
            A = B 
            dms.append(2) 
            B = apply(HM, dms+['zero'])
            l = [A.n(i) for i in range(A.order())]
            # Main loop performing the transposition of the entries
            for i in range(prod(l)):
                # Turning the index i into an hypermatrix array location using the decimal encoding trick
                entry = [mod(i,l[0])]
                sm = Integer(mod(i,l[0]))
                for k in range(len(l)-1):
                    entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                    sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                B[tuple(entry+[0])]=A[tuple(entry)]
                if A[tuple(entry)] == cos(t)^2:
                    B[tuple(entry+[1])]=sin(t)^2
                else :
                    B[tuple(entry+[1])]=cos(t)^2
        return B
    else :
        raise ValueError, "The order must be a positive integer"

def MatrixOrthogonalizationParametrization(sz):
    # Initializing the hadamard matrix
    H = hadamard_matrix(sz)
    # Initialization of the variables associated 
    # with the logarithm of the entries of the
    # orthogonal matrices
    LnQ = HM(sz,sz,'lnq')
    # Initialization of the parametrization variables
    C = HM(sz,sz,sz,'c')
    # Initialization of the list of variables
    VrbLst = LnQ.list()
    # Initialization of the list storing the constraints
    CnstrLst = []
    for j in range(sz):
        for i0 in range(sz):
            for i1 in range(sz):
                if i0!=i1:
                    CnstrLst.append(LnQ[i0,j]+LnQ[i1,j] == ln(sum([C[i0,i1,k]*H[j,k] for k in range(1,sz)])))
    CnstrLst = Set(CnstrLst).list()
    # Formating constraints in standard form
    [A,b] = ConstraintFormatorII(CnstrLst, VrbLst)
    # Solving for the consrtraints 
    import numpy
    sln = matrix(numpy.linalg.pinv(A))*b
    Q = HM(sz,sz,'zero')
    for j in range(sz):
        for k in range(sz):
            Q[j,k] = exp(sln[k*sz^1+j*sz^0,0])
    return Q   

def HypermatrixOrthogonalizationParametrization(sz):
    # Initializing the hadamard matrix
    H = hadamard_matrix(sz)
    # Initialization of the variables associated 
    # with the logarithm of the entries of the
    # orthogonal matrices
    LnQ = HM(sz,sz,sz,'lnq')
    # Initialization of the parametrization variables
    C = HM(sz,sz,sz,sz,'c')
    # Initialization of the list of variables
    VrbLst = LnQ.list()
    # Initialization of the list storing the constraints
    CnstrLst = []
    for j in range(sz):
        for i0 in range(sz):
            for i1 in range(sz):
                for i2 in range(sz):
                    if i0!=i1 or i1!=i2 or i0!=i2:
                        CnstrLst.append(LnQ[i0,j,i2]+LnQ[i1,j,i0]+LnQ[i2,j,i1] == ln(sum([C[i0,i1,i2,k]*H[j,k] for k in range(1,sz)])))
    CnstrLst = Set(CnstrLst).list()
    # Formating constraints in standard form
    [A,b] = ConstraintFormatorII(CnstrLst, VrbLst)
    # Solving for the consrtraints 
    import numpy
    sln = matrix(numpy.linalg.pinv(A))*b
    Q = HM(sz,sz,sz,'zero')
    for i in range(sz):
        for j in range(sz):
            for k in range(sz):
                Q[i,j,k] = exp(sln[k*sz^2+j*sz^1+i*sz^0,0])
    return Q 

def HypermatrixSliceKroneckerProduct(U, V):
    if (min([U.n(i) for i in range(U.order())])==max([U.n(i) for i in range(U.order())])) and (min([V.n(i) for i in range(V.order())])==max([V.n(i) for i in range(V.order())])):
        # Getting the sizes
        m = U.n(0); n = V.n(0)
        # Initializing the hypermatrices
        Q = HM(m*n, m*n, m*n, 'zero')
        # The buffering matrices
        M0 = matrix(SR, m, m, [0 for i in range(m^2)]); M1 = matrix(SR, n, n, [0 for i in range(n^2)])
        for u in range(m):
            for v in range(n):
                # Filling up the matrices
                for i in range(m):
                    for j in range(m):
                        M0[i,j]=U[i,u,j]
                for i in range(n):
                    for j in range(n):
                        M1[i,j]=V[i,v,j]
                # Computing the matrix outer product
                M = M0.tensor_product(M1)
                for a in range(m*n):
                    for b in range(m*n):
                        Q[a,m*v+u,b]=M[a,b]
        return Q
    else :
        raise ValueError, "Input dimensions must all be the same "

@cached_function
def Ca3(n):
    """
    Outputs the number of formula-binary trees only using fan-in three addition gates.
    This is a special  case of the Fuss-Catalan sequence.

    EXAMPLES:
    The input n must be greater than 0
    ::
        sage: Ca3(3)
        1


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n == 1 :
        return 1
    else :
        return sum([Ca3(i)*Ca3(j)*Ca3(n-i-j) for i in range(1,n,2) for j in range(1,n-i,2)])


