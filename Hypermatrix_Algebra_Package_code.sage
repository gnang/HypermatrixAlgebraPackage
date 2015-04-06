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

# Definition of the hypermatrix class HM.
class HM:
    """Hypermatrix class"""
    def __init__(self,*args):
        if len(args) == 1:
            inp = args[0]
            if type(inp)==type(Matrix(SR,2,1,[var('x'),var('y')])) or type(inp)==type(Matrix(RR,2,1,[1,2])) or type(inp)==type(Matrix(CC,2,1,[1,1])):
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
        # Obtaining the last argument
    	s = args[-1]
        # Initializing the dimension parameters
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
            self.hm=apply(GeneralHypermatrixKroneckerDelta, args[:-1]).listHM()
        elif s == 'sym':
            if len(dims) == 2:
                self.hm=SymHypermatrixGenerate(dims[0],dims[1])
            else:
                raise ValueError, "kronecker not supported for order %d tensors" % len(dims)
        elif type(s) == list:
            self.hm=(apply(List2Hypermatrix, args)).listHM()
        else:
            self.hm=apply(HypermatrixGenerate, args)
    def __repr__(self):
        return `self.hm`
    def __add__(self, other):
        return GeneralHypermatrixAdd(self,other)
    def __radd__(self, other):
        return GeneralHypermatrixAdd(self,other)
    def __neg__(self):
        return GeneralHypermatrixScale(self.hm,-1)
    def __sub__(self, other):
        return GeneralHypermatrixAdd(self, GeneralHypermatrixScale(other,-1))
    def __mul__(self, other):
        if other.__class__.__name__=='HM':
            return HM(GeneralHypermatrixHadamardProduct(self,other))
        elif other.__class__.__name__=='tuple':
            # This function takes a list as intput
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
        return GeneralHypermatrixHadamardProduct(self, B)
    def elementwise_exponent(self,s):
        return GeneralHypermatrixExponent(self, s)
    def elementwise_base_exponent(self, s):
        return GeneralHypermatrixBaseExponent(self, s)
    def elementwise_base_logarithm(self, s):
        return GeneralHypermatrixLogarithm(self, s)
    def slicekroneckerproduct(self, V):
        if  self.order()==2:
            return SecondOrderSliceKroneckerProduct(self, V)
        elif self.order()==3:
            return ThirdOrderSliceKroneckerProduct(self, V)
        elif self.order()==4:
            return FourthOrderSliceKroneckerProduct(self, V)
        elif self.order()==5:
            return FifthOrderSliceKroneckerProduct(self, V)
        else :
            raise ValueError,"not supported for order %d tensors" % len(dims)
    def block_sum(self, V):
        return GeneralHypermatrixKroneckerSum(self, V)
    def expand(self):
        return GeneralHypermatrixExpand(self)
    def simplify(self):
        return GeneralHypermatrixSimplify(self)
    def subs(self,Dct):
        return GeneralHypermatrixSubstitute(self, Dct)
    def subsn(self,Dct):
        return GeneralHypermatrixSubstituteN(self, Dct)
    def transpose(self, i=1):
        t = Integer(mod(i, self.order()))
        A = self 
        for i in range(t):
            A = GeneralHypermatrixCyclicPermute(A)
        return A
    def dagger(self, i=1):
        t = Integer(mod(i, self.order()))
        A = self 
        for i in range(t):
            A = GeneralHypermatrixCyclicPermute(A)
        if Integer(mod(i,2))==0:
            return A
        else:
            return GeneralHypermatrixConjugate(A)
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
        if i==0:
            return self.nrows()
        elif i==1:
            return self.ncols()
        elif i==2:
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
    def order(self):
        cnt = 0
        H = self.listHM()
        while type(H) == type([]):
            H = H[0]
            cnt = cnt+1
        return cnt
    def dimensions(self):
        return [self.n(i) for i in range(self.order())]
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
        # Obtaining the size of the image
        X = misc.toimage(pylab.array(self.listHM()))
        g = graphics_array([[matrix_plot(X)]])
        g.show()
    def save(self,filename):
        import pylab, numpy
        from scipy import misc
        # obtaining the image size
        X = misc.toimage(pylab.array(self.listHM()))
        X.save(filename)
    def copy(self):
        return GeneralHypermatrixCopy(self)
    def norm(self,p=2):
        return sum([abs(i)^p for i in self.list()])
    def det(self):
        if self.order()==2:
            return Deter(Matrix(SR,self.listHM()))
        elif self.order()==3:
            return ThirdOrderHyperdeterminant(self)
        elif self.order()==4:
            return FourthOrderHyperdeterminant(self)
        elif self.order()==5:
            return FifthOrderHyperdeterminant(self)
        elif self.order()==6:
            return SixthOrderHyperdeterminant(self)
        else:
            raise ValueError, "The hyperdeterminant is not implemented for order "+str(self.order())+"."
    def is_zero(self):
        if Set(self.list()).list()==[0]:
            return True
        else:
            return False
    def is_symmetric(self):
        return (self-self.transpose()).is_zero()

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
        raise ValueError, "Input dimensions "+str(nr)+" and "+str(nc)+" must both be non-zero positive integers."

def SymMatrixGenerate(nr, c):
    """
    Generates a list of lists associated with a symbolic nr x nr
    symmetric matrix using the input character c followed by
    indices.

    EXAMPLES:
    ::
        sage: M = SymMatrixGenerate(2, 'm'); M
        [[m00, m01], [m01, m11]]

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
        raise ValueError, "Input dimensions "+str(nr)+" must be a non-zero positive integers."

def HypermatrixGenerate(*args):
    """
    Generates a list of lists associated with a symbolic arbitrary
    order hypematrix of size specified by the input using and the
    entries are determined by the last input character

    EXAMPLES:
    ::
        sage: M = HypermatrixGenerate(2, 2, 2, 'm'); M
        [[[m000, m001], [m010, m011]], [[m100, m101], [m110, m111]]]


    AUTHORS:
    - Edinah K. Gnang, Ori Parzanchevski and Yuval Filmus
    """
    if len(args) == 1:
        return var(args[0])
    return [apply(HypermatrixGenerate, args[1:-1]+(args[-1]+str(i),)) for i in range(args[0])]

def HypermatrixGenerateAllOne(*args):
    """
    Generates a list of lists associated with a symbolic arbitrary
    order hypematrix of size specified by the input using and the
    entries are determined by the last input character

    EXAMPLES:
    ::
        sage: M = HypermatrixGenerateAllOne(2, 2, 2); M
        [[[1, 1], [1, 1]], [[1, 1], [1, 1]]]


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
        [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]

    AUTHORS:
    - Edinah K. Gnang, Ori Parzanchevski and Yuval Filmus
    """
    if len(args) == 1:
        return [0 for i in range(args[0])]
    return [apply(HypermatrixGenerateAllZero, args[1:] ) for i in range(args[0])]

def SymHypermatrixGenerate(nr, c):
    """
    Generates a list of lists associated with a symbolic third order hypermatrix of size
    nr x nc x nd third order hypematrix using the input character c followed by indices.

    EXAMPLES:
    ::
        sage: M = SymHypermatrixGenerate(2, 'm'); M
        [[[m000, m001], [m001, m011]], [[m001, m011], [m011, m111]]]

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
                        (q[i][j]).append(var(c+str(min(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))+str(max(i,j,k))))
                    else:
                        if i == min(i,j,k) and k == max(i,j,k):
                            (q[i][j]).append(var(c+str(min(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))+str(max(i,j,k))))
                        elif k == min(i,j,k) and j == max(i,j,k):
                            (q[i][j]).append(var(c+str(min(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))+str(max(i,j,k))))
                        elif i == max(i,j,k) and j == min(i,j,k):
                            (q[i][j]).append(var(c+str(min(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))+str(max(i,j,k))))
                        else:
                            (q[i][j]).append(var(c+str(i+j+k-min(i,j,k)-max(i,j,k))+str(min(i,j,k))+str(max(i,j,k))))
        return q

    else :
        raise ValueError, "Input dimensions "+str(nr)+" must be a non-zero positive integer."

def HypermatrixVectorize(A):
    """
    Outputs our canonical vectorization list associated with
    the input third order hypermatrices A.

    EXAMPLES:
    ::
        sage: A = SymHypermatrixGenerate(2, 'a') 
        sage: Lst = HypermatrixVectorize(A); Lst
        [a000, a001, a001, a011, a001, a011, a011, a111]


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
    the two input third order hypermatrices A and B

    EXAMPLES:
    ::
        sage: A = HypermatrixGenerate(2,2,2,'a'); B = HypermatrixGenerate(2,2,2,'b')
        sage: Rslt = HypermatrixAdd(A, B); Rslt
        [[[a000 + b000, a001 + b001], [a010 + b010, a011 + b011]], [[a100 + b100, a101 + b101], [a110 + b110, a111 + b111]]]


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
    Outputs a list of lists associated with the Hadamard product of
    the two input third order hypermatrices A and B

    EXAMPLES:
    ::
        sage: A = HypermatrixGenerate(2,2,2,'a'); B = HypermatrixGenerate(2,2,2,'b')
        sage: Rslt = HypermatrixHadamardProduct(A, B); Rslt
        [[[a000 * b000, a001 * b001], [a010 * b010, a011 * b011]], [[a100 * b100, a101 * b101], [a110 * b110, a111 * b111]]]


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
    scalar s with the third order hypermatrix A.

    EXAMPLES:
    ::
        sage: A = HypermatrixGenerate(2,2,2,'a')
        sage: Rslt = HypermatrixScale(A, 3); Rslt
        [[[3*a000, 3*a001], [3*a010, 3*a011]], [[3*a100, 3*a101], [3*a110, 3*a111]]]

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
    Outputs a list of lists associated with raising every entry of the 
    third order input hypermatrix A by the scalar s.

    EXAMPLES:
    ::
        sage: A = HypermatrixGenerate(2,2,2,'a')
        sage: Rslt = HypermatrixEntryExponent(A, 3); Rslt
        [[[a000^3, a001^3], [a010^3, a011^3]], [[a100^3, a101^3], [a110^3, a111^3]]]

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
    Outputs a list of lists associated with the exponentiation by the
    scalar s of every entry of the third order input hypermatrix A.

    EXAMPLES:
    ::
        sage: A = HypermatrixGenerate(2,2,2,'a')
        sage: Rslt = HypermatrixEntryExponentB(3,A); Rslt
        [[[3^a000, 3^a001], [3^a010, 3^a011]], [[3^a100, 3^a101], [3^a110, 3^a111]]]


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
    Bhattacharya-Mesner product of the input third order 
    hypermatrices A, B and C.
    The code writen here handles both list of list data
    structures as well as the Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: M = HypermatrixProduct(A, B, C); M


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    typ = 'list'
    if type(A) == type(HM(1,1,1,'one')) and type(B) == type(HM(1,1,1,'one')) and type(C)==type(HM(1,1,1,'one')):
        A = A.listHM(); B = B.listHM(); C = C.listHM()
        typ = 'HM'
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
    if n_a_rows==n_b_rows and n_b_cols==n_c_cols and n_c_dpts==n_a_dpts and n_a_cols==n_b_dpts and n_b_dpts==n_c_rows:
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
                    (q[i][j]).append(sum([A[i][l][k]*B[i][j][l]*C[l][j][k] for l in range(n_a_cols)]))
        if typ=='list':
            return q
        else:
            return HM(q)

    else :
        raise ValueError, "Hypermatrix dimension mismatch."

def HypermatrixLogProduct(A, B, C):
    """
    Outputs a list of lists associated with the ternary
    product of the input hypermatrices A, B and C.

    EXAMPLES:
    ::
        sage: Ha = HypermatrixGenerate(2, 1, 2, 'a')
        sage: Hb = HypermatrixGenerate(2, 2, 1, 'b')
        sage: Hc = HypermatrixGenerate(1, 2, 2, 'c')
        sage: Rslt = HypermatrixProduct(Ha, Hb, Hc); Rslt


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    typ = 'list'
    if type(A) == type(HM(1,1,1,'one')) and type(B) == type(HM(1,1,1,'one')) and type(C)==type(HM(1,1,1,'one')):
        A = A.listHM(); B = B.listHM(); C = C.listHM()
        typ = 'HM'
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
    if n_a_rows==n_b_rows and n_b_cols==n_c_cols and n_c_dpts==n_a_dpts and n_a_cols==n_b_dpts and n_b_dpts==n_c_rows:
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
                    (q[i][j]).append(sum([A[i][l][k]+B[i][j][l]+C[l][j][k] for l in range(n_a_cols)]))
        if typ == 'list':
            return q
        else:
            return HM(q)
    else :
        raise ValueError, "Hypermatrix dimension mismatch."

def HypermatrixKroneckerProduct(A, B, C):
    """
    Outputs a list of lists associated with the ternary
    analog of the Kronecker product for the inputh third
    order hypermatrices A, B and C.

    EXAMPLES:
    ::
        sage: X=HypermatrixGenerate(2,1,1,'x'); Y=HypermatrixGenerate(1,2,1,'y'); Z=HypermatrixGenerate(1,1,2,'z') 
        sage: T=HypermatrixKroneckerProduct(X, Y, Z); T
        [[[x000*y000*z000, x000*y000*z001], [x000*y010*z000, x000*y010*z001]], [[x100*y000*z000, x100*y000*z001], [x100*y010*z000, x100*y010*z001]]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    typ = 'list'
    if type(A) == type(HM(1,1,1,'one')) and type(B) == type(HM(1,1,1,'one')) and type(C)==type(HM(1,1,1,'one')):
        A = A.listHM(); B = B.listHM(); C = C.listHM()
        typ = 'HM'
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
        if typ=='list':
            return q
        else:
            return HM(q)
    else :
        raise ValueError, "Hypermatrix dimension mismatch."

def HypermatrixProductB(A, B, C, D):
    """
    Outputs a list of lists associated with the ternary
    product the input hypermatrices A, B and C with
    background hypermatrix D.

    EXAMPLES:
    ::
        sage: Ha=HypermatrixGenerate(2,2,2,'a')
        sage: Hb=HypermatrixGenerate(2,2,2,'b')
        sage: Hc=HypermatrixGenerate(2,2,2,'c')
        sage: Hd=HypermatrixGenerate(2,2,2,'d')
        sage: Rslt=HypermatrixProductB(Ha,Hb,Hc,Hd); Rslt
        [[[a000*b000*c000*d000 + a000*b000*c100*d001 + a000*b001*c000*d010 + a000*b001*c100*d011 + a010*b000*c000*d100 + a010*b000*c100*d101 + a010*b001*c000*d110 + a010*b001*c100*d111, a001*b000*c001*d000 + a001*b000*c101*d001 + a001*b001*c001*d010 + a001*b001*c101*d011 + a011*b000*c001*d100 + a011*b000*c101*d101 + a011*b001*c001*d110 + a011*b001*c101*d111], [a000*b010*c010*d000 + a000*b010*c110*d001 + a000*b011*c010*d010 + a000*b011*c110*d011 + a010*b010*c010*d100 + a010*b010*c110*d101 + a010*b011*c010*d110 + a010*b011*c110*d111, a001*b010*c011*d000 + a001*b010*c111*d001 + a001*b011*c011*d010 + a001*b011*c111*d011 + a011*b010*c011*d100 + a011*b010*c111*d101 + a011*b011*c011*d110 + a011*b011*c111*d111]], [[a100*b100*c000*d000 + a100*b100*c100*d001 + a100*b101*c000*d010 + a100*b101*c100*d011 + a110*b100*c000*d100 + a110*b100*c100*d101 + a110*b101*c000*d110 + a110*b101*c100*d111, a101*b100*c001*d000 + a101*b100*c101*d001 + a101*b101*c001*d010 + a101*b101*c101*d011 + a111*b100*c001*d100 + a111*b100*c101*d101 + a111*b101*c001*d110 + a111*b101*c101*d111], [a100*b110*c010*d000 + a100*b110*c110*d001 + a100*b111*c010*d010 + a100*b111*c110*d011 + a110*b110*c010*d100 + a110*b110*c110*d101 + a110*b111*c010*d110 + a110*b111*c110*d111, a101*b110*c011*d000 + a101*b110*c111*d001 + a101*b111*c011*d010 + a101*b111*c111*d011 + a111*b110*c011*d100 + a111*b110*c111*d101 + a111*b111*c011*d110 + a111*b111*c111*d111]]]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    typ = 'list'
    if type(A)==type(HM(1,1,1,'one')) and type(B)==type(HM(1,1,1,'one')) and type(C)==type(HM(1,1,1,'one')) and type(D)==type(HM(1,1,1,'one')):
        A = A.listHM(); B = B.listHM(); C = C.listHM(); D = D.listHM()
        typ = 'HM'
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
    if n_a_rows==n_b_rows and n_b_cols==n_c_cols and n_c_dpts==n_a_dpts and n_a_cols==n_b_dpts and n_b_dpts==n_c_rows and n_a_cols==n_d_rows and n_a_cols==n_d_cols and n_a_cols==n_d_dpts:
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
                    (q[i][j]).append(sum([A[i][l0][k]*B[i][j][l1]*C[l2][j][k]*D[l0][l1][l2] for l0 in range(n_d_rows) for l1 in range(n_d_cols) for l2 in range(n_d_dpts)]))
        if typ=='list':
            return q
        else:
            return HM(q)
    else :
        raise ValueError, "Hypermatrix dimension mismatch."

def HypermatrixDualProductB(A, B, C, D):
    """
    Outputs a list of lists associated with the ternary
    dual product the third order input hypermatrices A, B
    and C with background hypermatrix D which relates to 
    the product introduced by Richard Kerner.

    EXAMPLES:
    ::
        sage: Ha=HypermatrixGenerate(2,2,2,'a')
        sage: Hb=HypermatrixGenerate(2,2,2,'b')
        sage: Hc=HypermatrixGenerate(2,2,2,'c')
        sage: Hd=HypermatrixGenerate(2,2,2,'d')
        sage: Rslt = HypermatrixDualProductB(Ha, Hb, Hc, Hd); Rslt
        [[[a000*b000*c000*d000 + a100*b100*c000*d000 + a001*b000*c001*d000 + a101*b100*c001*d000 + a000*b010*c010*d000 + a100*b110*c010*d000 + a001*b010*c011*d000 + a101*b110*c011*d000, a000*b000*c100*d001 + a100*b100*c100*d001 + a001*b000*c101*d001 + a101*b100*c101*d001 + a000*b010*c110*d001 + a100*b110*c110*d001 + a001*b010*c111*d001 + a101*b110*c111*d001], [a000*b001*c000*d010 + a100*b101*c000*d010 + a001*b001*c001*d010 + a101*b101*c001*d010 + a000*b011*c010*d010 + a100*b111*c010*d010 + a001*b011*c011*d010 + a101*b111*c011*d010, a000*b001*c100*d011 + a100*b101*c100*d011 + a001*b001*c101*d011 + a101*b101*c101*d011 + a000*b011*c110*d011 + a100*b111*c110*d011 + a001*b011*c111*d011 + a101*b111*c111*d011]], [[a010*b000*c000*d100 + a110*b100*c000*d100 + a011*b000*c001*d100 + a111*b100*c001*d100 + a010*b010*c010*d100 + a110*b110*c010*d100 + a011*b010*c011*d100 + a111*b110*c011*d100, a010*b000*c100*d101 + a110*b100*c100*d101 + a011*b000*c101*d101 + a111*b100*c101*d101 + a010*b010*c110*d101 + a110*b110*c110*d101 + a011*b010*c111*d101 + a111*b110*c111*d101], [a010*b001*c000*d110 + a110*b101*c000*d110 + a011*b001*c001*d110 + a111*b101*c001*d110 + a010*b011*c010*d110 + a110*b111*c010*d110 + a011*b011*c011*d110 + a111*b111*c011*d110, a010*b001*c100*d111 + a110*b101*c100*d111 + a011*b001*c101*d111 + a111*b101*c101*d111 + a010*b011*c110*d111 + a110*b111*c110*d111 + a011*b011*c111*d111 + a111*b111*c111*d111]]]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    typ = 'list'
    if type(A)==type(HM(1,1,1,'one')) and type(B)==type(HM(1,1,1,'one')) and type(C)==type(HM(1,1,1,'one')) and type(D)==type(HM(1,1,1,'one')):
        A = A.listHM(); B = B.listHM(); C = C.listHM(); D = D.listHM()
        typ = 'HM'
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
    if n_a_rows==n_b_rows and n_b_cols==n_c_cols and n_c_dpts==n_a_dpts and n_a_cols==n_b_dpts and n_b_dpts==n_c_rows and n_a_cols==n_d_rows and n_a_cols==n_d_cols and n_a_cols==n_d_dpts:
        # Initialization of the hypermatrix
        q = []
        for i in range(n_a_rows):
            q.append([])
        for i in range(len(q)):
            for j in range(n_b_cols):
                (q[i]).append([])
        for l0 in range(len(q)):
            for l1 in range(len(q[i])):
                for l2 in range(n_c_dpts):
                    (q[l0][l1]).append( sum([A[i][l0][k]*B[i][j][l1]*C[l2][j][k]*D[l0][l1][l2] for i in range(n_a_rows) for j in range(n_b_cols) for k in range(n_c_dpts)]))
        if typ=='list':
            return q
        else:
            return HM(q)
    else :
        raise ValueError, "Hypermatrix dimension mismatch."

def HypermatrixCyclicPermute(A):
    """
    Outputs a list of lists associated with the third order
    hypermatrix with entries index cyclicaly permuted.

    EXAMPLES:
    ::
        sage: Ha=HypermatrixGenerate(2,2,2,'a'); Ha
        [[[a000, a001], [a010, a011]], [[a100, a101], [a110, a111]]]
        sage: HypermatrixCyclicPermute(Ha)
        [[[a000, a100], [a001, a101]], [[a010, a110], [a011, a111]]]


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
    Generates a list of lists associated with the third order
    nr x nr x nr Kronecker Delta hypermatrix.

    EXAMPLES:
    ::
        sage: Dlt = HypermatrixKroneckerDelta(2); Dlt
        [[[1, 0], [0, 0]], [[0, 0], [0, 1]]]


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
        raise ValueError, "Input dimensions "+str(nr)+" must be a non-zero positive integer."

def HypermatrixPermutation(s):
    """
    Generates a list of lists associated with the permutation
    hypermatrix deduced from sigma. Note that as a result of 
    the  non associativity, permutations must be performed as
    one transposition at a time.

    EXAMPLES:
    ::
        sage: P = HypermatrixPermutation([0,2,1]); P
        [[[1, 0, 0], [0, 0, 1], [0, 1, 0]], [[1, 0, 0], [0, 0, 1], [0, 1, 0]], [[1, 0, 0], [0, 0, 1], [0, 1, 0]]]


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
        raise ValueError, "Input dimensions "+str(n)+" must be a non-zero positive integer."

def DiagonalHypermatrix(Mtrx):
    """
    Outputs a diagonal third order hypermatrix
    constructed using the input square matrix.
    We enforce the symmetry constraint by only
    taking entries from the lower triangular
    part of the input matrix.

     EXAMPLES:
    ::
        sage: var('a00, a11, a01')
        sage: Mtrx = Matrix(SR, [[a00, a01], [a01, a11]])
        sage: Dg = DiagonalHypermatrix(Mtrx); Dg
        [[[a00, 0], [0, a01]], [[a01, 0], [0, a11]]]


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
    Outputs a symbolic parametrization of third order orthogonal hypermatrix
    of size 2x2x2.

     EXAMPLES:
    ::
        sage: t=var('t')
        sage: Orthogonal2x2x2Hypermatrix(t)
        [[[cos(t)^(2/3), sin(t)^(2/3)], [sin(t)^(2/3), cos(t)^(2/3)]], [[-sin(t)^(2/3), cos(t)^(2/3)], [sin(t)^(2/3), sin(t)^(2/3)]]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    return [[[cos(t)^(2/3),sin(t)^(2/3)],[sin(t)^(2/3), cos(t)^(2/3)]], [[-sin(t)^(2/3),cos(t)^(2/3)],[sin(t)^(2/3),sin(t)^(2/3)]]]

def Orthogonal3x3x3Hypermatrix(t1,t2):
    """
    Outputs a symbolic parametrization of third order orthogonal hypermatrix
    of size 3x3x3.

     EXAMPLES:
    ::
        sage: t1, t2=var('t1, t2')
        sage: Orthogonal3x3x3Hypermatrix(t1,t2)
        [[[cos(t1)^(2/3), cos(t2)^(2/3)*sin(t1)^(2/3), 0], [cos(t2)^(2/3)*sin(t1)^(2/3), sin(t1)^(2/3)*sin(t2)^(2/3), 0], [sin(t1)^(2/3)*sin(t2)^(2/3), cos(t1)^(2/3)*e^(-2/3*I*pi), 0]], [[sin(t1)^(2/3)*sin(t2)^(2/3), cos(t1)^(2/3), cos(t2)^(2/3)*e^(-2/3*I*pi)*sin(t1)^(2/3)], [cos(t1)^(2/3)*e^(2/3*I*pi), cos(t2)^(2/3)*sin(t1)^(2/3), sin(t1)^(2/3)*sin(t2)^(2/3)], [cos(t2)^(2/3)*sin(t1)^(2/3), sin(t1)^(2/3)*sin(t2)^(2/3), cos(t1)^(2/3)]], [[0, sin(t1)^(2/3)*sin(t2)^(2/3), cos(t1)^(2/3)], [0, cos(t1)^(2/3), cos(t2)^(2/3)*sin(t1)^(2/3)], [0, cos(t2)^(2/3)*e^(2/3*I*pi)*sin(t1)^(2/3), sin(t1)^(2/3)*sin(t2)^(2/3)]]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    c1=cos(t1)^(2/3)
    s1=sin(t1)^(2/3)
    c2=cos(t2)^(2/3)
    s2=sin(t2)^(2/3)
    return [[[c1,s1*c2,0],[s1*c2,s1*s2,0],[s1*s2,exp(-I*2*pi/3)*c1,0]], [[s1*s2,c1,exp(-I*2*pi/3)*s1*c2],[exp(I*2*pi/3)*c1,s1*c2,s1*s2], [s1*c2,s1*s2,c1]],[[0,s1*s2,c1],[0,c1,s1*c2],[0,exp(I*2*pi/3)*s1*c2,s1*s2]]]

def HypermatrixCayleyHamiltonList(A,n):
    """
    Outpts a list of hypermatrices (each of which is encapsulated as a single list) of all product composition of degree n.
    This function is most adapted for constructing matrices.

     EXAMPLES:
    ::
        sage: A = HypermatrixGenerate(2,2,2,'a')
        sage: Lst = HypermatrixCayleyHamiltonList(A,3); Lst
        [[a000^3 + a001*a010*a100, a000*a001^2 + a001*a011*a101, a000*a010^2 + a010*a011*a110, a001*a010*a011 + a011^2*a111, a000*a100^2 + a100*a101*a110, a001*a100*a101 + a101^2*a111, a010*a100*a110 + a110^2*a111, a011*a101*a110 + a111^3]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    if n == 1:
        return [A.list()]
    else:
        gu = []
        for i in range(1,n,2):
            for j in range(1,n-i,2):
                gu = gu + [HypermatrixProduct(HM(A.n(0),A.n(1),A.n(2),g1), HM(A.n(0),A.n(1),A.n(2),g2), HM(A.n(0),A.n(1),A.n(2),g3)).list() for g1 in HypermatrixCayleyHamiltonList(A,i) for g2 in HypermatrixCayleyHamiltonList(A,j) for g3 in HypermatrixCayleyHamiltonList(A,n-(i+j))]
        return gu

def HypermatrixCayleyHamiltonListII(A,n):
    """
    Outpts a list of third order hypermatrices of all product composition of degree n.

     EXAMPLES:
    ::
        sage: A = HypermatrixGenerate(2,2,2,'a')
        sage: Lst = HypermatrixCayleyHamiltonList(A,3); Lst
        [[[[a000^3 + a001*a010*a100, a000*a001^2 + a001*a011*a101], [a000*a010^2 + a010*a011*a110, a001*a010*a011 + a011^2*a111]], [[a000*a100^2 + a100*a101*a110, a001*a100*a101 + a101^2*a111], [a010*a100*a110 + a110^2*a111, a011*a101*a110 + a111^3]]]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    if n == 1:
        return [A]
    else:
        gu = []
        for i in range(1,n,2):
            for j in range(1,n-i,2):
                gu = gu + [HypermatrixProduct(g1,g2,g3) for g1 in HypermatrixCayleyHamiltonList(A,i) for g2 in HypermatrixCayleyHamiltonList(A,j) for g3 in HypermatrixCayleyHamiltonList(A,n-(i+j))]
        return gu

def HypermatrixCayleyHamiltonListIII(A,n):
    """
    Outpts a list of hypermatrices of all product
    composition of degree n for hypermatrices of 
    order up to 8 tends to be slow because it uses
    the general hypermatrix product and it takes HM
    object class as the input A

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
                    gu = gu + [GeneralHypermatrixProduct(g1,g2,g3) for g1 in HypermatrixCayleyHamiltonListII(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,n-(i+j))]
            return gu
    # Case of order 4
    elif A.order()==4:
        if n == 1:
            return [A]
        else:
            gu = []
            for i in range(1,n,2):
                for j in range(1,n-i,2):
                    for k in range(1,n-i-j,2):
                        gu = gu + [GeneralHypermatrixProduct(g1,g2,g3,g4) for g1 in HypermatrixCayleyHamiltonListII(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,k) for g4 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k))]
    # Case of order 5
    elif A.order()==5:
        if n == 1:
            return [A]
        else:
            gu = []
            for i in range(1,n,2):
                for j in range(1,n-i,2):
                    for k in range(1,n-i-j,2):
                        for l in range(1,n-i-j-k,2):
                            gu = gu + [GeneralHypermatrixProduct(g1,g2,g3,g4,g5) for g1 in HypermatrixCayleyHamiltonListII(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,k) for g4 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k)) for g5 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l))]
    # Case of order 6
    elif A.order()==6:
        if n == 1:
            return [A]
        else:
            gu = []
            for i in range(1,n,2):
                for j in range(1,n-i,2):
                    for k in range(1,n-i-j,2):
                        for l in range(1,n-i-j-k,2):
                            for m in range(1,n-i-j-k-l,2):
                                gu = gu + [GeneralHypermatrixProduct(g1,g2,g3,g4,g5,g6) for g1 in HypermatrixCayleyHamiltonListII(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,k) for g4 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k)) for g5 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l)) for g6 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m))]
    # Case of order 7
    elif A.order()==7:
        if n == 1:
            return [A]
        else:
            gu = []
            for i in range(1,n,2):
                for j in range(1,n-i,2):
                    for k in range(1,n-i-j,2):
                        for l in range(1,n-i-j-k,2):
                            for m in range(1,n-i-j-k-l,2):
                                for o in range(1,n-i-j-k-l,2):
                                    gu = gu + [GeneralHypermatrixProduct(g1,g2,g3,g4,g5,g6,g7) for g1 in HypermatrixCayleyHamiltonListII(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,k) for g4 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k)) for g5 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l)) for g6 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m)) for g7 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m+o))]
    # Case of order 8
    elif A.order()==8:
        if n == 1:
            return [A]
        else:
            gu = []
            for i in range(1,n,2):
                for j in range(1,n-i,2):
                    for k in range(1,n-i-j,2):
                        for l in range(1,n-i-j-k,2):
                            for m in range(1,n-i-j-k-l,2):
                                for o in range(1,n-i-j-k-l,2):
                                    for p in range(1,n-i-j-k-l-o,2):
                                        gu = gu + [GeneralHypermatrixProduct(g1,g2,g3,g4,g5,g6,g7) for g1 in HypermatrixCayleyHamiltonListII(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,k) for g4 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k)) for g5 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l)) for g6 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m)) for g7 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m+o))  for g8 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m+o+p))]
    else :
        raise ValueError, "Not supported for order > 4 and for non cube hypermpatrix of order 3 "
   
def HypermatrixSymCayleyHamiltonList(A,n):
    """
    Outpts a list of symmetric hypermatrices of all product
    composition of degree n.

     EXAMPLES:
    ::
        sage: A = HM(2,2,2,'a')
        sage: L = HypermatrixSymCayleyHamiltonList(A,3)

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
                    gu = gu + [GeneralHypermatrixProduct(g1,g2.transpose(2),g3.transpose()) for g1 in HypermatrixSymCayleyHamiltonList(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,n-(i+j))]
            return gu
    # Case of order 4
    elif A.order()==4:
        if n == 1:
            return [A]
        else:
            gu = []
            for i in range(1,n,2):
                for j in range(1,n-i,2):
                    for k in range(1,n-i-j,2):
                        gu = gu + [GeneralHypermatrixProduct(g1,g2.transpose(3),g3.transpose(2),g4.transpose()) for g1 in HypermatrixSymCayleyHamiltonList(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,k) for g4 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k))]
    # Case of order 5
    elif A.order()==5:
        if n == 1:
            return [A]
        else:
            gu = []
            for i in range(1,n,2):
                for j in range(1,n-i,2):
                    for k in range(1,n-i-j,2):
                        for l in range(1,n-i-j-k,2):
                            gu = gu + [GeneralHypermatrixProduct(g1,g2.transpose(4),g3.transpose(3),g4.transpose(2),g5.transpose()) for g1 in HypermatrixSymCayleyHamiltonList(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,k) for g4 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k)) for g5 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l))]
    # Case of order 6
    elif A.order()==6:
        if n == 1:
            return [A]
        else:
            gu = []
            for i in range(1,n,2):
                for j in range(1,n-i,2):
                    for k in range(1,n-i-j,2):
                        for l in range(1,n-i-j-k,2):
                            for m in range(1,n-i-j-k-l,2):
                                gu = gu + [GeneralHypermatrixProduct(g1,g2.transpose(5),g3.transpose(4),g4.transpose(3),g5.transpose(2),g6.transpose()) for g1 in HypermatrixSymCayleyHamiltonList(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,k) for g4 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k)) for g5 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l)) for g6 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m))]
    # Case of order 7
    elif A.order()==7:
        if n == 1:
            return [A]
        else:
            gu = []
            for i in range(1,n,2):
                for j in range(1,n-i,2):
                    for k in range(1,n-i-j,2):
                        for l in range(1,n-i-j-k,2):
                            for m in range(1,n-i-j-k-l,2):
                                for o in range(1,n-i-j-k-l,2):
                                    gu = gu + [GeneralHypermatrixProduct(g1,g2.transpose(6),g3.transpose(5),g4.transpose(4),g5.transpose(3),g6.transpose(2),g7.transpose()) for g1 in HypermatrixSymCayleyHamiltonList(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,k) for g4 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k)) for g5 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l)) for g6 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m)) for g7 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m+o))]
    # Case of order 8
    elif A.order()==8:
        if n == 1:
            return [A]
        else:
            gu = []
            for i in range(1,n,2):
                for j in range(1,n-i,2):
                    for k in range(1,n-i-j,2):
                        for l in range(1,n-i-j-k,2):
                            for m in range(1,n-i-j-k-l,2):
                                for o in range(1,n-i-j-k-l,2):
                                    for p in range(1,n-i-j-k-l-o,2):
                                        gu = gu + [GeneralHypermatrixProduct(g1,g2.transpose(7),g3.transpose(6),g4.transpose(5),g5.transpose(4),g6.transpose(3),g7.transpose(2),g8.transpose()) for g1 in HypermatrixSymCayleyHamiltonList(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,k) for g4 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k)) for g5 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l)) for g6 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m)) for g7 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m+o))  for g8 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m+o+p))]
    # Case of order 9
    elif A.order()==9:
        if n == 1:
            return [A]
        else:
            gu = []
            for i in range(1,n,2):
                for j in range(1,n-i,2):
                    for k in range(1,n-i-j,2):
                        for l in range(1,n-i-j-k,2):
                            for m in range(1,n-i-j-k-l,2):
                                for o in range(1,n-i-j-k-l,2):
                                    for p in range(1,n-i-j-k-l-o,2):
                                        for q in range(1,n-i-j-k-l-o-p,2):
                                            gu = gu + [GeneralHypermatrixProduct(g1,g2.transpose(8),g3.transpose(7),g4.transpose(6),g5.transpose(5),g6.transpose(4),g7.transpose(3),g8.transpose(2),g9.transpose()) for g1 in HypermatrixSymCayleyHamiltonList(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,k) for g4 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k)) for g5 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l)) for g6 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m)) for g7 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m+o))  for g8 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m+o+p)) for g9 in HypermatrixCayleyHamiltonListII(A,n-(i+j+k+l+m+o+p+q)) ]
    else :
        raise ValueError, "Not supported for order > 4 and for non cube hypermpatrix of order 3 "

def ConstraintFormator(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs matrix
    and the right hand side vector associate
    with the matrix formulation of the constraints.

    EXAMPLES:
    ::
        sage: x, y = var('x, y')
        sage: CnstrLst = [x + y == 1, x - y == 2]
        sage: VrbLst = [x, y]
        sage: [A, b] = ConstraintFormator(CnstrLst, VrbLst)
        sage: A
        [ 1.00000000000000  1.00000000000000]
        [ 1.00000000000000 -1.00000000000000]
        sage: b
        [1.00000000000000]
        [2.00000000000000]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initializing the Matrix
    A=Matrix(CC,len(CnstrLst),len(VrbLst),zero_matrix(len(CnstrLst),len(VrbLst)))
    b=vector(CC, [eq.rhs() for eq in CnstrLst]).column()
    for r in range(len(CnstrLst)):
        for c in range(len(VrbLst)):
            A[r,c]=(CnstrLst[r]).lhs().coefficient(VrbLst[c])
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
        sage: A
        [ 1  1]
        [ 1 -1]
        sage: b
        [1]
        [2]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initializing the Matrix
    A=Matrix(SR,len(CnstrLst),len(VrbLst),zero_matrix(len(CnstrLst),len(VrbLst)))
    b=vector(SR, [eq.rhs() for eq in CnstrLst]).column()
    for r in range(len(CnstrLst)):
        for c in range(len(VrbLst)):
            A[r,c]=(CnstrLst[r]).lhs().coefficient(VrbLst[c])
    return [A,b]

def ConstraintFormatorIII(CnstrLst, VrbLst):
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
        sage: A
        [ 1  1]
        [ 1 -1]
        sage: b
        [1]
        [2]
       

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initializing the Matrix
    A=Matrix(CC,len(CnstrLst),len(VrbLst),zero_matrix(len(CnstrLst),len(VrbLst)))
    b=vector(SR, [eq.rhs() for eq in CnstrLst]).column()
    for r in range(len(CnstrLst)):
        for c in range(len(VrbLst)):
            A[r,c]=(CnstrLst[r]).lhs().coefficient(VrbLst[c])
    return [A,b]

def Companion_matrix(p,vrbl):
    """
    Takes as input a polynomial and a variable
    and outputs the companion matrix associated
    with the polynomial in the specified variables.

    EXAMPLES:
    ::
        sage: x,y=var('x,y')
        sage: p=x^2+2*x*y+1
        sage: Companion_matrix(p,x)
        [0   1]
        [1 2*y]

    AUTHORS:
    - Edinah K. Gnang
    """
    if p.is_polynomial(vrbl):
        dg=p.degree(vrbl)
        if dg>1:
            # Initialization of the matrix
            A=Matrix(SR,HM(dg,dg,'zero').listHM())
            # Filling up the matrix
            A[0,dg-1]=-p.subs(dict([(vrbl,0)]))/p.coefficient(vrbl^dg)
            for i in range(1,dg):
                A[dg-1,i]=-p.coefficient(vrbl^i)/p.coefficient(vrbl^dg); A[i,i-1]=1
            return A
        else:
            raise ValueError, "Must be of degree at least 2."
    else:
        raise ValueError, "Must be a polynomial in the input variable."

def substitute_matrix(p, vrbl, A):
    """
    The functions takes as input a polynomial p,
    a variable vrbl, and a matrix A. The function
    outputs the polynomial in the variable.

    EXAMPLES:
    ::
        sage: x,y = var('x,y')
        sage: p=x^2+2*x*y+1
        sage: substitute_matrix(p,x,Matrix(SR,HM(2,2,'a').listHM()))
        [a00^2+a01*a10+2*a00*y+1   a00*a01+a01*a11+2*a01*y]
        [a00*a10+a10*a11+2*a10*y   a01*a10+a11^2+2*a11*y+1]

    AUTHORS:
    - Edinah K. Gnang
    """
    if A.nrows()==A.ncols():
        d=p.degree(vrbl)
        T=Matrix(SR, zero_matrix(d,d))
        T=T+p.subs(dict([(vrbl,0)]))*identity_matrix(A.nrows())
        for i in range(1,d+1):
            T=T+(A^i)*p.coefficient(vrbl^i)
        return T
    else:
        raise ValueError, "Must be a polynomial in the input variable."
     
def HypermatrixPseudoInversePairs(A,B):
    """
     Outputs the pseudo inverse pairs associated with the input pairs of matrices

    EXAMPLES:
    ::
        sage: A1=[[[0.1631135370902057,0.11600112072013125],[0.9823708115400902,0.39605960486710756]] ,[[0.061860929755424676,0.2325542810173995],[0.39111210957450926,0.2019809359102137]]]
        sage: A2=[[[0.15508921433883183,0.17820377184410963],[0.48648171594508205,0.01568017636082064]] ,[[0.8250247759993575,0.1938307874191597],[0.23867299119274843,0.3935578730402869]]]
        sage: [B1,B2]=HypermatrixPseudoInversePairs(A1,A2)

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    sz = len(A)
    # Initializing the list of linear constraints
    CnstrLst = []
    # Initilizing the variable list
    Vrbls  = [var('ln_al'+str(i)+str(j)+str(k)) for i in range(sz) for j in range(sz) for k in range(sz)]+[var('ln_bt'+str(i)+str(j)+str(k)) for i in range(sz) for j in range(sz) for k in range(sz)]
    for m in range(sz):
        for p in range(sz):
            for n in range(sz):
                V=Matrix(CC, sz, sz, [(A[m][k1][k0])*(B[k0][k1][p]) for k0 in range(sz) for k1 in range(sz)]).inverse()
                CnstrLst=CnstrLst+[var('ln_al'+str(m)+str(n)+str(k1))+var('ln_bt'+str(k1)+str(n)+str(p))==ln(V[k1,n])  for k1 in range(sz)]
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

def GeneralHypermatrixProduct(*args):
    """
    Outputs a list of lists associated with the general
    Bhattacharya-Mesner product of the input hypermatrices.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: Ha=HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c')
        sage: Rslt=GeneralHypermatrixProduct(Ha, Hb, Hc); Rslt
        [[[a000*b000*c000 + a010*b001*c100, a001*b000*c001 + a011*b001*c101], [a000*b010*c010 + a010*b011*c110, a001*b010*c011 + a011*b011*c111]], [[a100*b100*c000 + a110*b101*c100, a101*b100*c001 + a111*b101*c101], [a100*b110*c010 + a110*b111*c110, a101*b110*c011 + a111*b111*c111]]]


    AUTHORS:
    - Edinah K. Gnang
    """
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
            Rh[tuple(entry)]=sum([prod([args[s][tuple(entry[0:Integer(mod(s+1,len(args)))]+[t]+entry[Integer(mod(s+2,len(args))):])] for s in range(len(args)-2)]+[args[len(args)-2][tuple(entry[0:len(args)-1]+[t])]]+[args[len(args)-1][tuple([t]+entry[1:])]]) for t in range((args[0]).n(1))])
    return Rh

# Defining a shorter function call for the hypermatrix product implemented above.
def Prod(*args):
    """
    Outputs a list of lists associated with the general
    Bhattacharya-Mesner product of the input hypermatrices.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: Ha=HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c')
        sage: Rslt=Prod(Ha, Hb, Hc); Rslt
        [[[a000*b000*c000 + a010*b001*c100, a001*b000*c001 + a011*b001*c101], [a000*b010*c010 + a010*b011*c110, a001*b010*c011 + a011*b011*c111]], [[a100*b100*c000 + a110*b101*c100, a101*b100*c001 + a111*b101*c101], [a100*b110*c010 + a110*b111*c110, a101*b110*c011 + a111*b111*c111]]]


    AUTHORS:
    - Edinah K. Gnang
    """

    return GeneralHypermatrixProduct(*args)

def GeneralHypermatrixProductB(*args):
    """
    Outputs a list of lists associated with the general
    Bhattacharya-Mesner product of the input hypermatrices
    with non-trivial background. The code only handles the 
    Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: Ha=HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c'); Hd = HM(2,2,2,'d')
        sage: Rslt=GeneralHypermatrixProductB(Ha, Hb, Hc, Hd); Rslt
        [[[a000*b000*c000*d000 + a000*b000*c100*d001 + a000*b001*c000*d010 + a000*b001*c100*d011 + a010*b000*c000*d100 + a010*b000*c100*d101 + a010*b001*c000*d110 + a010*b001*c100*d111, a001*b000*c001*d000 + a001*b000*c101*d001 + a001*b001*c001*d010 + a001*b001*c101*d011 + a011*b000*c001*d100 + a011*b000*c101*d101 + a011*b001*c001*d110 + a011*b001*c101*d111], [a000*b010*c010*d000 + a000*b010*c110*d001 + a000*b011*c010*d010 + a000*b011*c110*d011 + a010*b010*c010*d100 + a010*b010*c110*d101 + a010*b011*c010*d110 + a010*b011*c110*d111, a001*b010*c011*d000 + a001*b010*c111*d001 + a001*b011*c011*d010 + a001*b011*c111*d011 + a011*b010*c011*d100 + a011*b010*c111*d101 + a011*b011*c011*d110 + a011*b011*c111*d111]], [[a100*b100*c000*d000 + a100*b100*c100*d001 + a100*b101*c000*d010 + a100*b101*c100*d011 + a110*b100*c000*d100 + a110*b100*c100*d101 + a110*b101*c000*d110 + a110*b101*c100*d111, a101*b100*c001*d000 + a101*b100*c101*d001 + a101*b101*c001*d010 + a101*b101*c101*d011 + a111*b100*c001*d100 + a111*b100*c101*d101 + a111*b101*c001*d110 + a111*b101*c101*d111], [a100*b110*c010*d000 + a100*b110*c110*d001 + a100*b111*c010*d010 + a100*b111*c110*d011 + a110*b110*c010*d100 + a110*b110*c110*d101 + a110*b111*c010*d110 + a110*b111*c110*d111, a101*b110*c011*d000 + a101*b110*c111*d001 + a101*b111*c011*d010 + a101*b111*c111*d011 + a111*b110*c011*d100 + a111*b110*c111*d101 + a111*b111*c011*d110 + a111*b111*c111*d111]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list specifying the dimensions of the output
    l = [(args[i]).n(i) for i in range(len(args)-1)]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Initializing the background hypermatrix
    B = (args[len(args)-1]).transpose(args[len(args)-1].order()-1)
    args = tuple([args[id] for id in range(len(args)-1)])
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # Computing the Hypermatrix product
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
                Rh[tuple(entry)] = Rh[tuple(entry)]+prod([args[s][tuple(entry[0:Integer(mod(s+1,len(args)))]+[entry2[ Integer(mod(s+1,len(args))) ]]+entry[Integer(mod(s+2,len(args))):])] for s in range(len(args)-2)]+[args[len(args)-2][tuple(entry[0:len(entry)-1]+[entry2[len(entry2)-1]])]]+[args[len(args)-1][tuple([entry2[0]]+entry[1:])]])*B[tuple(entry2)]
    return Rh

def ProdB(*args):
    """
    Outputs a list of lists associated with the general
    Bhattacharya-Mesner product of the input hypermatrices
    with non-trivial background. The code only handles the 
    Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: Ha=HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c'); Hd = HM(2,2,2,'d')
        sage: Rslt=ProdB(Ha, Hb, Hc, Hd); Rslt
        [[[a000*b000*c000*d000 + a000*b000*c100*d001 + a000*b001*c000*d010 + a000*b001*c100*d011 + a010*b000*c000*d100 + a010*b000*c100*d101 + a010*b001*c000*d110 + a010*b001*c100*d111, a001*b000*c001*d000 + a001*b000*c101*d001 + a001*b001*c001*d010 + a001*b001*c101*d011 + a011*b000*c001*d100 + a011*b000*c101*d101 + a011*b001*c001*d110 + a011*b001*c101*d111], [a000*b010*c010*d000 + a000*b010*c110*d001 + a000*b011*c010*d010 + a000*b011*c110*d011 + a010*b010*c010*d100 + a010*b010*c110*d101 + a010*b011*c010*d110 + a010*b011*c110*d111, a001*b010*c011*d000 + a001*b010*c111*d001 + a001*b011*c011*d010 + a001*b011*c111*d011 + a011*b010*c011*d100 + a011*b010*c111*d101 + a011*b011*c011*d110 + a011*b011*c111*d111]], [[a100*b100*c000*d000 + a100*b100*c100*d001 + a100*b101*c000*d010 + a100*b101*c100*d011 + a110*b100*c000*d100 + a110*b100*c100*d101 + a110*b101*c000*d110 + a110*b101*c100*d111, a101*b100*c001*d000 + a101*b100*c101*d001 + a101*b101*c001*d010 + a101*b101*c101*d011 + a111*b100*c001*d100 + a111*b100*c101*d101 + a111*b101*c001*d110 + a111*b101*c101*d111], [a100*b110*c010*d000 + a100*b110*c110*d001 + a100*b111*c010*d010 + a100*b111*c110*d011 + a110*b110*c010*d100 + a110*b110*c110*d101 + a110*b111*c010*d110 + a110*b111*c110*d111, a101*b110*c011*d000 + a101*b110*c111*d001 + a101*b111*c011*d010 + a101*b111*c111*d011 + a111*b110*c011*d100 + a111*b110*c111*d101 + a111*b111*c011*d110 + a111*b111*c111*d111]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    return GeneralHypermatrixProductB(*args) 

def GeneralHypermatrixLogProduct(*args):
    """
    Outputs a list of lists associated with the general
    Bhattacharya-Mesner Log-product of the input hypermatrices
    with non-trivial background. The code only handles the 
    Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: Ha=HM(2,1,2,'a'); Hb=HM(2,2,1,'b'); Hc=HM(1,2,2,'c')
        sage: Rslt=GeneralHypermatrixLogProduct(Ha, Hb, Hc); Rslt
        [[[a000 + b000 + c000, a001 + b000 + c001], [a000 + b010 + c010, a001 + b010 + c011]], [[a100 + b100 + c000, a101 + b100 + c001], [a100 + b110 + c010, a101 + b110 + c011]]]

    AUTHORS:
    - Edinah K. Gnang
    """
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
            Rh[tuple(entry)]=sum([sum([args[s][tuple(entry[0:Integer(mod(s+1,len(args)))]+[t]+entry[Integer(mod(s+2,len(args))):])] for s in range(len(args)-2)]+[args[len(args)-2][tuple(entry[0:len(args)-1]+[t])]]+[args[len(args)-1][tuple([t]+entry[1:])]]) for t in range((args[0]).n(1))])
    return Rh

def GeneralHypermatrixCyclicPermute(A):
    """
    Outputs a list of lists associated with the general
    transpose as defined by the cyclic permutation of indices.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: Ha=HM(2,2,2,'a')
        sage: GeneralHypermatrixCyclicPermute(Ha)
        [[[a000, a100], [a001, a101]], [[a010, a110], [a011, a111]]]

    AUTHORS:
    - Edinah K. Gnang
    """
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
    """
    Outputs a list of lists associated with the scaling of a general hypermatrix.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: Ha=HM(2,2,2,'a')
        sage: GeneralHypermatrixScale(Ha,3)
        [[[3*a000, 3*a001], [3*a010, 3*a011]], [[3*a100, 3*a101], [3*a110, 3*a111]]]

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
        Rh[tuple(entry)]=s*A[tuple(entry)]
    return Rh

def GeneralHypermatrixExponent(A,s):
    """
    Outputs a list of lists associated with the general
    whose entries are all raised to the power s.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: Ha=HM(2,2,2,'a')
        sage: GeneralHypermatrixExponent(Ha,3)
        [[[a000^3, a001^3], [a010^3, a011^3]], [[a100^3, a101^3], [a110^3, a111^3]]]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list specifying the dimensions of the output
    l = [A.n(i) for i in range(A.order())]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the computations of the entries
    for i in range(prod(l)):
        entry = [mod(i,l[0])]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        if A[tuple(entry)].is_zero():
            Rh[tuple(entry)] = 0
        else:
            Rh[tuple(entry)] = (A[tuple(entry)])^s
    return Rh

def GeneralHypermatrixBaseExponent(A,s):
    """
    Outputs a list of lists associated with the general
    whose entries are exponentiated using the input s as
    basis for the exponentiation.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: Ha=HM(2,2,2,'a')
        sage: GeneralHypermatrixBaseExponent(Ha,3)
        [[[3^a000, 3^a001], [3^a010, 3^a011]], [[3^a100, 3^a101], [3^a110, 3^a111]]]

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
        Rh[tuple(entry)]=s^(A[tuple(entry)])
    return Rh

def GeneralHypermatrixLogarithm(A,s=e):
    """
    Outputs a list of lists associated with the general
    whose entries are logarithms to the base s of the 
    original hypermatrix.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: Ha=HM(2,2,2,'a')
        sage: GeneralHypermatrixLogarithm(Ha,3)
        

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
        Rh[tuple(entry)]=log(A[tuple(entry)],s).canonicalize_radical()
    return Rh

def GeneralHypermatrixConjugate(A):
    """
    Outputs a list of lists associated with the general
    whose entries are complex conjugates of the original
    hypermatrix.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: Ha=HM(2,2,2,[exp(I*2*pi*u*v*w/4) for u in range(2) for v in range(2) for w in range(2)])
        sage: GeneralHypermatrixConjugate(Ha)
        [[[1, 1], [1, 1]], [[1, 1], [1, -I]]]

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
        Rh[tuple(entry)]=conjugate(A[tuple(entry)])
    return Rh

def GeneralHypermatrixExpand(A):
    """
    Outputs a list of lists associated with the general
    hypermatrix with expressions in the entries in their
    expanded form.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: Ha=HM(2,2,2,[(var('x')+var('y'))^(i+j+k) for i in range(2) for j in range(2) for k in range(2)])
        sage: GeneralHypermatrixExpand(Ha)
        [[[1, x + y], [x + y, x^2 + 2*x*y + y^2]], [[x + y, x^2 + 2*x*y + y^2], [x^2 + 2*x*y + y^2, x^3 + 3*x^2*y + 3*x*y^2 + y^3]]]
        

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
        Rh[tuple(entry)]=A[tuple(entry)]
    return Rh

def List2Hypermatrix(*args):
    """
    Procedure for Initializing a Hypermatrix from a size specifications and a list

    EXAMPLES:
    ::
        sage: A = List2Hypermatrix(2,2,2,HM(2,2,2,'a').list()); A
        [[[a000, a001], [a010, a011]],[[a100, a101], [a110, a111]]]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list specifying the dimensions of the output
    l = args[:-1]
    # Initialization of the list
    Lst =args[-1]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = [j for j in l]+['zero']
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
        Rh[tuple(entry)]=Lst[i]
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
    - Edinah K. Gnang
    """
    # The if statement bellow address the sum function
    if B == 0:
        tl = [A.n(i) for i in range(A.order())]+['zero']
        B = HM(*tl)
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
    - Edinah K. Gnang
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
    """
    Outputs a list of lists associated with the general
    Kronecter delta hypermatrix
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: Dlt = GeneralHypermatrixKroneckerDelta(2,2); Dlt
        [[1, 0], [0, 1]] 
        

    AUTHORS:
    - Edinah K. Gnang
    """
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

def GeneralUncorrelatedHypermatrixTupleU(od):
    """
    Generates a tuplet of hypermatrices of the appropriate order which are
    uncorrelated but I do not normalize them. Each one of the dimensions are equal to 2.

    EXAMPLES:
    ::
        sage: [A,B]=GeneralUncorrelatedHypermatrixTupleU(2)
        sage: A
        [[e^(I*pi + r1 - r2 + r6), e^r6], [e^(I*pi + r3 - r4 + r5), e^r5]]
        sage: B
        [[e^r4, e^r2], [e^r3, e^r1]]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization the alphabet list
    AlphaB = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    # Initializing the hypermatrix
    LQ = [apply(HM,[2 for i in range(od)]+[AlphaB[j]]).elementwise_base_exponent(e) for j in range(od)]
    # Initilizing the list of variable
    VrbLst = []
    for Q in LQ:
        VrbLst = VrbLst + (Q.elementwise_base_logarithm(e)).list()
    # Computing the product
    Eq = apply(GeneralHypermatrixProduct, [Q for Q in LQ])
    # Writting up the constraints
    LeQ = (Eq.list())[1:2^od-1]
    # Filling up the linear constraints
    CnstrLst= [] 
    for f in LeQ:
        CnstrLst.append(ln((f.operands())[0]).canonicalize_radical()-I*pi-ln((f.operands())[1]).canonicalize_radical()==0)
    # Directly solving the constraints
    Sl = solve(CnstrLst, VrbLst)
    # Setting up the list for the dictionary
    Dct = [(eq.lhs(),exp(eq.rhs())) for eq in Sl[0]]
    # Returning the uncorrelated tuplets
    return [apply(HM, [2 for i in range(od)]+[AlphaB[j]]).subs(dict(Dct)) for j in range(od)]

def GeneralUncorrelatedHypermatrixTuple(od):
    """
    Generates a tuplet of hypermatrices of the appropriate order which are
    uncorrelated and normalized. Each one of the dimensions are equal to 2.

    EXAMPLES:
    ::
        sage: [A,B]=GeneralUncorrelatedHypermatrixTuple(2)
        sage: A
        [[e^(I*pi + r1 - r2 + r6)/(e^(I*pi + r1 - r2 + r4 + r6) + e^(r3 + r6)), e^r6/(e^(I*pi + r1 - r2 + r4 + r6) + e^(r3 + r6))], [e^(I*pi + r3 - r4 + r5)/(e^(I*pi + r2 + r3 - r4 + r5) + e^(r1 + r5)), e^r5/(e^(I*pi + r2 + r3 - r4 + r5) + e^(r1 + r5))]]
        sage: B
        [[e^r4, e^r2], [e^r3, e^r1]]
        sage: (Prod(A,B)).simplify()
        [[1, 0], [0, 1]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the unormalized tuples
    L = GeneralUncorrelatedHypermatrixTupleU(od)
    Tp= apply(GeneralHypermatrixProduct, [h for h in L])
    Q = L[0].copy()
    # first row to normalize 
    entry = [0 for i in range(Q.order())]
    Q[tuple(entry)]=Q[tuple(entry)]/Tp[tuple([0 for i in range(Q.order())])]
    entry[1] = 1
    Q[tuple(entry)]=Q[tuple(entry)]/Tp[tuple([0 for i in range(Q.order())])]
    # last row to normalize 
    entry = [1 for i in range(Q.order())]
    Q[tuple(entry)]=Q[tuple(entry)]/Tp[tuple([1 for i in range(Q.order())])]
    entry[1] = 0
    Q[tuple(entry)]=Q[tuple(entry)]/Tp[tuple([1 for i in range(Q.order())])]
    return [Q]+[L[i] for i in range(1,len(L))]

def GeneralOrthogonalHypermatrixU(od):
    """
    Generates an orthogonal hypermatrix of the appropriate order
    for which each one of the dimensions are equal to 2.
    The vectors are not normalized.

    EXAMPLES:
    ::
        sage: Q = GeneralOrthogonalHypermatrixU(3); Q
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
    LeQ.remove(e^(od*var('q'+''.join(['0' for i in range(od)])))+e^(od*var('q01'+''.join(['0' for i in range(od-2)]))))
    LeQ.remove( e^(od*var('q10'+''.join(['1' for i in range(od-2)])))+e^(od*var('q'+''.join(['1' for i in range(od)]))))
    # Filling up the linear constraints
    CnstrLst= [] 
    for f in LeQ:
        CnstrLst.append(ln((f.operands())[0]).canonicalize_radical()-I*pi-ln((f.operands())[1]).canonicalize_radical()==0)
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
        Q[tuple(entry)]=Q[tuple(entry)].subs(dict(map(lambda eq: (eq.lhs(),eq.rhs()), Sl[0]))).canonicalize_radical()
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
            CnstrLst.append(ln((f.operands())[0]).canonicalize_radical()-I*pi-ln((f.operands())[1]).canonicalize_radical()==0)
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
            Q[tuple(entry)]=Q[tuple(entry)].subs(dict(map(lambda eq: (eq.lhs(),eq.rhs()), Sl[0]))).canonicalize_radical()
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

def GeneralUnitaryHypermatrixU(od):
    """
    Generates an unitary hypermatrix of the appropriate order
    for which each one of the dimensions are equal to 2.
    The vectors are not normalized. The order input od
    must be even for the function call to be meaningful.

    EXAMPLES:
    ::
        sage: [A, Ac]=GeneralUnitaryHypermatrixU(4)
        sage: (Prod(A,Ac.transpose(3),A.transpose(2),Ac.transpose())).simplify()
        [[[[(e^(4*r18+4*r19+4*r21+4*r24) + e^(4*r19+4*r20+4*r22+4*r23))*e^(-4*r18-4*r21-4*r24),0],[0,0]],[[0,0],[0,0]]],[[[0,0],[0,0]], [[0,0],[0,(e^(4*r17+4*r18+4*r21+4*r24) + e^(4*r17+4*r20+4*r22+4*r23))*e^(-4*r20-4*r22-4*r23)]]]]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the variable playing the role of sqrt(-1)
    z=var('z')
    # Initializing the real part and the imaginary part
    X=apply(HM,[2 for i in range(od)]+['x']); Y=apply(HM,[2 for i in range(od)]+['y'])
    # Initialization of the list
    Lh = [(X+z*Y).elementwise_base_exponent(e), (X-z*Y).elementwise_base_exponent(e)]
    # Computation of the Product.
    B = apply(Prod,[Lh[Integer(mod(i,2))].transpose(i) for i in range(od,0,-1)])
    B[tuple([0 for i in range(od)])]=0; B[tuple([1 for i in range(od)])]=0
    # Initializing the list
    L=Set(B.list()).list()
    # Removing the normalization constraints
    L.remove(0)
    # Initialization of the variables.
    Vrbls =X.list()+Y.list() 
    # Initialization of the non homogeneous equations
    #Eq =[(log((l.operands())[0])).simplify_log().subs(z=0) - (log((l.operands())[1])).simplify_log().subs(z=0) == 0 for l in L]+[(log((l.operands())[0])).simplify_log().coeff(z) - (log((l.operands())[1])).simplify_log().coeff(z)==pi for l in L]
    Eq =[(log((l.operands())[0])).simplify_log().subs(z=0) - (log((l.operands())[1])).simplify_log().subs(z=0) == 0 for l in L]+[(log((l.operands())[0])).simplify_log().coefficient(z) - (log((l.operands())[1])).simplify_log().coefficient(z)==pi for l in L]
    # Calling the constraint formator
    [A,b]=ConstraintFormatorII(Eq,Vrbls)
    # Setting the signs right
    V = A.kernel().basis()
    for i in range(Integer(len(V)/2),len(V)):
        # Initilization of the index locators
        c1=-1; c2=-1
        for j in range(len(V[i])):
            if V[i][j] == 1 and c1 == -1:
                c1=j
            elif V[i][j] == 1 and c1 != -1:
                c2=j
                break
        b[c2,0] = -b[c2,0]
    # Rewriting the system
    Eq = [(A*Matrix(SR,len(Vrbls),1,Vrbls))[i,0]==b[i,0] for i in range(A.nrows())]
    # Computing the Homogeneous solution
    Sln = solve(Eq,Vrbls)[0]
    X = X.subs(dict([(s.lhs(),s.rhs()) for s in Sln]))
    Y = Y.subs(dict([(s.lhs(),s.rhs()) for s in Sln]))
    # Final result
    return [(X+I*Y).elementwise_base_exponent(e), (X-I*Y).elementwise_base_exponent(e)]

def GeneralUnitaryHypermatrix(od):
    """
    Generates an unitary hypermatrix of the appropriate order
    for which each one of the dimensions are equal to 2.
    The vectors are not normalized.

    EXAMPLES:
    ::
        sage: [A, Ac]=GeneralUnitaryHypermatrixU(4)
        sage: (Prod(A,Ac.transpose(3),A.transpose(2),Ac.transpose())).simplify()
        [[[[(e^(4*r18+4*r19+4*r21+4*r24) + e^(4*r19+4*r20+4*r22+4*r23))*e^(-4*r18-4*r21-4*r24),0],[0,0]],[[0,0],[0,0]]],[[[0,0],[0,0]], [[0,0],[0,(e^(4*r17+4*r18+4*r21+4*r24) + e^(4*r17+4*r20+4*r22+4*r23))*e^(-4*r20-4*r22-4*r23)]]]]

    AUTHORS:
    - Edinah K. Gnang
    """
    Lh=GeneralUnitaryHypermatrixU(od)
    Tp=apply(Prod,[Lh[Integer(mod(i,2))].transpose(i) for i in range(od,0,-1)])
    Q =Lh[0].copy()
    Qc=Lh[1].copy()
    # first row to normalize 
    entry=[0 for i in range(Q.order())]
    Q[tuple(entry)]=Q[tuple(entry)]/(Tp[tuple([0 for i in range(Q.order())])])^(1/od)
    Qc[tuple(entry)]=Qc[tuple(entry)]/(Tp[tuple([0 for i in range(Q.order())])])^(1/od)
    entry[1]=1
    Q[tuple(entry)]=Q[tuple(entry)]/(Tp[tuple([0 for i in range(Q.order())])])^(1/od)
    Qc[tuple(entry)]=Qc[tuple(entry)]/(Tp[tuple([0 for i in range(Q.order())])])^(1/od)
    # last row to normalize 
    entry=[1 for i in range(Q.order())]
    Q[tuple(entry)]=Q[tuple(entry)]/(Tp[tuple([1 for i in range(Q.order())])])^(1/od)
    Qc[tuple(entry)]=Qc[tuple(entry)]/(Tp[tuple([1 for i in range(Q.order())])])^(1/od)
    entry[1]=0
    Q[tuple(entry)]=Q[tuple(entry)]/(Tp[tuple([1 for i in range(Q.order())])])^(1/od)
    Qc[tuple(entry)]=Qc[tuple(entry)]/(Tp[tuple([1 for i in range(Q.order())])])^(1/od)
    # Final result
    return [Q,Qc]

def DFT_image_resizer(sz, dm):
    """
    Generates a third order hypermatrix of 3 slices
    for performing third reduction on the size in each color
    channel.

    EXAMPLES:
    ::
        sage: DFT_image_resizer(4,2)
        [[[[1, 1, 1], [0, 0, 0], [1, 1, 1], [0, 0, 0]], [[1, 1, 1], [0, 0, 0], [-1, -1, -1], [0, 0, 0]], [[0, 0, 0], [1, 1, 1], [0, 0, 0], [1, 1, 1]], [[0, 0, 0], [1, 1, 1], [0, 0, 0], [-1, -1, -1]]],
 [[[1, 1, 1], [1, 1, 1], [0, 0, 0], [0, 0, 0]], [[0, 0, 0], [0, 0, 0], [1, 1, 1], [1, 1, 1]], [[1, 1, 1], [-1, -1, -1], [0, 0, 0], [0, 0, 0]], [[0, 0, 0], [0, 0, 0], [1, 1, 1], [-1, -1, -1]]]]


    AUTHORS:
    - Edinah K. Gnang
    """
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
    """
    Performs channel specific matrix multiplication

    EXAMPLES:
    ::
        sage: Ha=HM(2,2,3,'a'); Hb=HM(2,2,3,'b')
        sage:  channel_product(Ha,Hb)
        [[[a000*b000 + a010*b100, a001*b001 + a011*b101, a002*b002 + a012*b102], [a000*b010 + a010*b110, a001*b011 + a011*b111, a002*b012 + a012*b112]], [[a100*b000 + a110*b100, a101*b001 + a111*b101, a102*b002 + a112*b102], [a100*b010 + a110*b110, a101*b011 + a111*b111, a102*b012 + a112*b112]]]

    AUTHORS:
    - Edinah K. Gnang
    """
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

def Deter(A):
    """
    Computes symbolically the determinant of a square matrix
    using the sum over permutation formula.

    EXAMPLES:
    ::
        sage: M = Matrix(SR, MatrixGenerate(2, 2, 'm')); Deter(M)
        -m01*m10 + m00*m11

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the permutations
    P = Permutations(range(A.nrows()))
    return sum([Permutation([p[i]+1 for i in range(len(p))]).signature()*prod([A[k][p[k]] for k in range(A.nrows())]) for p in P])

 
def MeanApproximation(T):
    """
    Computes  the mean slice approximation. This is mostly used for images
    as a way to get a hold of the background.

    EXAMPLES:
    ::
        sage: MeanApproximation(HM(2,2,2,'a'))
        [[[[1/2*a000 + 1/2*a010, 1/2*a001 + 1/2*a011]], [[1/2*a100 + 1/2*a110, 1/2*a101 + 1/2*a111]]], [[[1/2*a000 + 1/2*a001], [1/2*a010 + 1/2*a011]], [[1/2*a100 + 1/2*a101], [1/2*a110 + 1/2*a111]]], [[[1/2*a000 + 1/2*a100, 1/2*a001 + 1/2*a101], [1/2*a010 + 1/2*a110, 1/2*a011 + 1/2*a111]]]]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the Slices.
    Ha = HM(T.n(0), 1     , T.n(2), 'zero')
    Hb = HM(T.n(0), T.n(1), 1,      'zero')
    Hc = HM(1     , T.n(1), T.n(2), 'zero')
    # Computing the mean of row depth slice
    for u in range(Ha.nrows()):
        for v in range(Ha.ndpts()):
            Ha[u,0,v] = mean([T[u,i,v] for i in range(T.ncols())])
    # Computing the mean row column slice
    for u in range(Hb.nrows()):
        for v in range(Hb.ncols()):
            Hb[u,v,0] = mean([T[u,v,i] for i in range(T.ndpts())])
    # Computing the mean column depth slice
    for u in range(Hc.ncols()):
        for v in range(Hc.ndpts()):
            Hc[0,u,v] = mean([T[i,u,v] for i in range(T.nrows())])
    # Computing the outer-product of the mean slices.
    return [Ha, Hb, Hc]

def ZeroPadding(A):
    """
    outputs the zero padding into a cube of the cuboid third order hypermatrix.

    EXAMPLES:
    ::
        sage: ZeroPadding(HM(1,1,2,'a'))
        [[[a000, a001], [0, 0]], [[0, 0], [0, 0]]]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the size parameter
    sz = max(A.dimensions())
    # Initializing the Hypermatrix
    T = HM(sz, sz, sz, 'zero')
    # Filling up the Hypermatrix
    for r in range(A.n(0)):
        for c in range(A.n(1)):
            for d in range(A.n(2)):
                T[r,c,d]=A[r,c,d]
    return T

def GenerateUnitLpNormVector(n,p = 2,indx=0):
    """
    outputs a unit lp norm vector.

    EXAMPLES:
    ::
        sage: GenerateUnitLpNormVector(2) 
        [cos(t0), sin(t0)]
        

    AUTHORS:
    - Edinah K. Gnang
    """
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
    """
    outputs the symbolic parametrization of a doubly stochastic matrix

    EXAMPLES:
    ::
        sage: ProbabilityMatrix(2)
        [     cos(t0)^2      sin(t0)^2]
        [-cos(t0)^2 + 1 -sin(t0)^2 + 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
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
    """
    outputs the symbolic parametrization of a symetric doubly stochastic matrix

    EXAMPLES:
    ::
        sage: ProbabilitySymMatrix(2)
        [     cos(t0)^2      sin(t0)^2]
        [-cos(t0)^2 + 1 -sin(t0)^2 + 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
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
        sage: A1=HM([[[0.1631135370902057,0.11600112072013125],[0.9823708115400902,0.39605960486710756]] ,[[0.061860929755424676,0.2325542810173995],[0.39111210957450926,0.2019809359102137]]])
        sage: A2=HM([[[0.15508921433883183,0.17820377184410963],[0.48648171594508205,0.01568017636082064]] ,[[0.8250247759993575,0.1938307874191597],[0.23867299119274843,0.3935578730402869]]])
        sage: [B1,B2]=HypermatrixPseudoInversePairsII(A1,A2)


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    sz = len(A.listHM())
    # Initializing the list of linear constraints
    CnstrLst = []
    # Initilizing the variable list
    Vrbls  = [var('ln_al'+str(i)+str(j)+str(k))  for i in range(sz) for j in range(sz) for k in range(sz)]+[var('ln_bt'+str(i)+str(j)+str(k)) for i in range(sz) for j in range(sz) for k in range(sz)]

    for m in range(sz):
        for p in range(sz):
            for n in range(sz):
                V=Matrix(CC, sz, sz, [(A[m,k1,k0])*(B[k0,k1,p]) for k0 in range(sz) for k1 in range(sz)]).inverse()
                CnstrLst=CnstrLst+[var('ln_al'+str(m)+str(n)+str(k1))+var('ln_bt'+str(k1)+str(n)+str(p))==ln(V[k1,n])  for k1 in range(sz)]
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
    """
     Outputs the pseudo inverse pairs associated with the input pairs of matrices

    EXAMPLES:
    ::
        sage: A1=HM([[[0.1631135370902057,0.11600112072013125],[0.9823708115400902,0.39605960486710756]] ,[[0.061860929755424676,0.2325542810173995],[0.39111210957450926,0.2019809359102137]]])
        sage: A2=HM([[[0.15508921433883183,0.17820377184410963],[0.48648171594508205,0.01568017636082064]] ,[[0.8250247759993575,0.1938307874191597],[0.23867299119274843,0.3935578730402869]]])
        sage: [B1,B2]=HypermatrixPseudoInversePairsIII(A1,A2)

    AUTHORS:
    - Edinah K. Gnang
    """
    # Defining the variables
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
                    Rx[i,j,k] = N(Mx[i,j,k]-stp_sz*Gv[indx].subs(dict([(X[u,v,w],Mx[u,v,w]) for u in range(X.n(0)) for v in range(X.n(1)) for w in range(X.n(2))]+[(Y[u,v,w],My[u,v,w]) for u in range(Y.n(0)) for v in range(Y.n(1)) for w in range(Y.n(2))]))/sum([(abs(Gv[l].subs(dict([(X[u,v,w],Mx[u,v,w]) for u in range(X.n(0)) for v in range(X.n(1)) for w in range(X.n(2))]+[(Y[u,v,w],My[u,v,w]) for u in range(Y.n(0)) for v in range(Y.n(1)) for w in range(Y.n(2))]))))^p for l in range(len(Gv))])^(1/p))
                    indx = indx+1
        # Updating the Slice Ry
        for i in range(Y.n(0)):
            for j in range(Y.n(1)):
                for k in range(Y.n(2)):
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
        sage: A1=HM([[[0.1631135370902057,0.11600112072013125],[0.9823708115400902,0.39605960486710756]], [[0.061860929755424676,0.2325542810173995],[0.39111210957450926,0.2019809359102137]]])
        sage: A2=HM([[[0.15508921433883183,0.17820377184410963],[0.48648171594508205,0.01568017636082064]], [[0.8250247759993575,0.1938307874191597],[0.23867299119274843,0.3935578730402869]]])
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
    """
     Outputs the pseudo inverse pairs associated with the input pairs of matrices

    EXAMPLES:
    ::
        sage: A1=HM([[[0.1631135370902057,0.11600112072013125],[0.9823708115400902,0.39605960486710756]], [[0.061860929755424676,0.2325542810173995],[0.39111210957450926,0.2019809359102137]]])
        sage: A2=HM([[[0.15508921433883183,0.17820377184410963],[0.48648171594508205,0.01568017636082064]], [[0.8250247759993575,0.1938307874191597],[0.23867299119274843,0.3935578730402869]]])
        sage: B=HypermatrixPseudoInversePairAction(HM(2,2,2,'one'),A1,A2)

    AUTHORS:
    - Edinah K. Gnang
    """
    # Computing the unsplit Inverse pairs
    XY = HypermatrixPseudoInversePairsUnsplit(A,B)
    Rs = HM(A.nrows(), A.ncols(), A.ndpts(),'zero')
    for m in range(Rs.nrows()):
        for n in range(Rs.ncols()):
            for p in range(Rs.ndpts()):
                Rs[m,n,p] = sum([T[m,j,p]*XY[m,n,j,p] for j in range(A.nrows())])
    return Rs 

def GenerateRandomHypermatrix(*l):
    """
     Outputs a random hypermatrix

    EXAMPLES:
    ::
        sage: A=GenerateRandomHypermatrix(2,2,2); A
        [[[0.5381896123570098, 0.17614064161822574], [0.002681324221582626, 0.3889148087368086]], [[0.35090401651983916, 0.5172460210635879], [0.5838421765928398, 0.3732279225915548]]]

    AUTHORS:
    - Edinah K. Gnang
    """
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

def GenerateRandomIntegerHypermatrix(*l):
    """
     Outputs a random hypermatrix

    EXAMPLES:
    ::
        sage: A=GenerateRandomIntegerHypermatrix(3,3,3); A
        [[[-2, 1, 4], [1, 2, 1], [3, 1, -4]], [[-5, 1, 0], [2, -2, 0], [0, -1, -6]], [[0, 2, 1], [-2, -8, -1], [-2, -2, 0]]]

    AUTHORS:
    - Edinah K. Gnang
    """
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
            Rh[tuple(entry)]=ZZ.random_element()
        return Rh
    else :
        raise ValueError, "The Dimensions must all be non-zero."


def GeneralStochasticHypermatrix(t, od):
    """
    Generates an stochastic hypermatrix of the appropriate order
    for which each one of the dimensions are equal to 2.
    The vectors are not normalized.

    EXAMPLES:
    ::
        sage: Q = GeneralStochasticHypermatrix(var('t'), 2); Q
        [[cos(t)^2, sin(t)^2],[sin(t)^2, cos(t)^2]]

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
    """
    Generates orthogonalization constraints.

    EXAMPLES:
    ::
        sage: MatrixOrthogonalizationParametrization(2)
        [[e^(-(1.96261557335e-17)*log(-c011) + 0.25*log(c011) + (1.96261557335e-17)*log(-c101) + 0.25*log(c101)), e^(0.25*log(-c011) + (3.92523114671e-17)*log(c011) + 0.25*log(-c101))], [e^(-(5.88784672006e-17)*log(-c011) + 0.25*log(c011) - (1.96261557335e-17)*log(-c101) + 0.25*log(c101)), e^(0.25*log(-c011) + (3.92523114671e-17)*log(c011) + 0.25*log(-c101))]]

    AUTHORS:
    - Edinah K. Gnang, Ori Parzanchevski
    """
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
    """
    Generates hypermatrix orthogonalization constraints.

    EXAMPLES:
    ::
        sage: f = HypermatrixOrthogonalizationParametrization(2)
        

    AUTHORS:
    - Edinah K. Gnang
    """
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

def SecondOrderSliceKroneckerProduct(Ha, Hb):
    """
    Computes the Kronecker Product for the two input second order hypermatrices.

    EXAMPLES:
    ::
        sage: A=SecondOrderSliceKroneckerProduct(HM(2,2,'a'),HM(2,2,'b'))
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the hypermatrix
    Hc = HM(Ha.n(0)*Hb.n(0), Ha.n(1)*Hb.n(1), 'zero')
    for i0 in range(Ha.n(0)):
        for i1 in range(Ha.n(1)):
            for j0 in range(Hb.n(0)):
                for j1 in range(Hb.n(1)):
                    Hc[Hb.n(0)*i0+j0,Hb.n(1)*i1+j1]=Ha[i0,i1]*Hb[j0,j1]
    return Hc

def ThirdOrderSliceKroneckerProduct(Ha, Hb):
    """
    Computes the Kronecker Product for the two input third order hypermatrices.

    EXAMPLES:
    ::
        sage: A=SecondOrderSliceKroneckerProduct(HM(2,2,'a'),HM(2,2,'b'))
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the hypermatrix
    Hc = HM(Ha.n(0)*Hb.n(0), Ha.n(1)*Hb.n(1), Ha.n(2)*Hb.n(2), 'zero')
    for i0 in range(Ha.n(0)):
        for i1 in range(Ha.n(1)):
            for i2 in range(Ha.n(2)):
                for j0 in range(Hb.n(0)):
                    for j1 in range(Hb.n(1)):
                        for j2 in range(Hb.n(2)):
                            Hc[Hb.n(0)*i0+j0, Hb.n(1)*i1+j1, Hb.n(2)*i2+j2]=Ha[i0,i1,i2]*Hb[j0,j1,j2]
    return Hc

def FourthOrderSliceKroneckerProduct(Ha, Hb):
    """
    Computes the Kronecker Product for the two input fourth order hypermatrices.

    EXAMPLES:
    ::
        sage: A=FourthOrderSliceKroneckerProduct(HM(2,2,'a'),HM(2,2,'b'))
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the hypermatrix
    Hc = HM(Ha.n(0)*Hb.n(0), Ha.n(1)*Hb.n(1), Ha.n(2)*Hb.n(2), Ha.n(3)*Hb.n(3) , 'zero')
    for i0 in range(Ha.n(0)):
        for i1 in range(Ha.n(1)):
            for i2 in range(Ha.n(2)):
                for i3 in range(Ha.n(3)):
                    for j0 in range(Hb.n(0)):
                        for j1 in range(Hb.n(1)):
                            for j2 in range(Hb.n(2)):
                                for j3 in range(Hb.n(3)):
                                    Hc[Hb.n(0)*i0+j0,Hb.n(1)*i1+j1,Hb.n(2)*i2+j2,Hb.n(3)*i3+j3]=Ha[i0,i1,i2,i3]*Hb[j0,j1,j2,j3]
    return Hc

def FifthOrderSliceKroneckerProduct(Ha, Hb):
    """
    Computes the Kronecker Product for the two input fifth order hypermatrices.

    EXAMPLES:
    ::
        sage: A=FourthOrderSliceKroneckerProduct(HM(2,2,'a'),HM(2,2,'b'))
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the hypermatrix
    Hc = HM(Ha.n(0)*Hb.n(0), Ha.n(1)*Hb.n(1), Ha.n(2)*Hb.n(2), Ha.n(3)*Hb.n(3), Ha.n(4)*Hb.n(4), 'zero')
    for i0 in range(Ha.n(0)):
        for i1 in range(Ha.n(1)):
            for i2 in range(Ha.n(2)):
                for i3 in range(Ha.n(3)):
                    for i4 in range(Ha.n(4)):
                        for j0 in range(Hb.n(0)):
                            for j1 in range(Hb.n(1)):
                                for j2 in range(Hb.n(2)):
                                    for j3 in range(Hb.n(3)):
                                        for j4 in range(Hb.n(4)):
                                            Hc[Hb.n(0)*i0+j0,Hb.n(1)*i1+j1,Hb.n(2)*i2+j2,Hb.n(3)*i3+j3,Hb.n(4)*i4+j4]=Ha[i0,i1,i2,i3,i4]*Hb[j0,j1,j2,j3,j4]
    return Hc

def ThirdOrderSliceKroneckerProductII(U, V):
    if (min(U.dimensions())==max(U.dimensions())) and (min(V.dimensions())==max(V.dimensions())):
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

def HypermatrixSliceKroneckerPower(U,n):
    """
    Computes the repeated Kronecker Product.

    EXAMPLES:
    ::
        sage: A=HypermatrixSliceKroneckerPower(HM(2,2,'a'),2)
        [[a00^2, a00*a01, a00*a01, a01^2], [a00*a10, a00*a11, a01*a10, a01*a11], [a00*a10, a01*a10, a00*a11, a01*a11], [a10^2, a10*a11, a10*a11, a11^2]]

    AUTHORS:
    - Edinah K. Gnang
    """
    T = U.copy()
    for i in range(n-1):
        T = T.slicekroneckerproduct(U)
    return T

def HypermatrixKroneckerSum(A,B):
    """
    Computes the  Kronecker sum for third order hypermatrix.

    EXAMPLES:
    ::
        sage: Rslt=HypermatrixKroneckerSum(A,B)

    AUTHORS:
    - Edinah K. Gnang
    """
    if A.order()==B.order() and A.order()==3:
        # Initializing the hypermatrix
        T = HM(A.n(0)+B.n(0), A.n(1)+B.n(1), A.n(2)+B.n(2), 'zero')
        # Loop filling up the top part of the hypermatrix
        for i in range(A.n(0)):
            for j in range(A.n(1)):
                for k in range(A.n(2)):
                    T[i,j,k]=A[i,j,k]
        # Loop filling up the bottom part of the hypermatrix
        for i in range(B.n(0)):
            for j in range(B.n(1)):
                for k in range(B.n(2)):
                    T[A.n(0)+i,A.n(1)+j,A.n(2)+k] = B[i,j,k]
        return T
    elif A.order()==B.order() and A.order()==4:
        # Initializing the hypermatrix
        T=HM(A.n(0)+B.n(0), A.n(1)+B.n(1), A.n(2)+B.n(2), A.n(3)+B.n(3), 'zero')
        # Loop filling up the top part of the hypermatrix
        for i in range(A.n(0)):
            for j in range(A.n(1)):
                for k in range(A.n(2)):
                    for l in range(A.n(3)):
                        T[i,j,k,l]=A[i,j,k,l]
        # Loop filling up the bottom part of the hypermatrix
        for i in range(B.n(0)):
            for j in range(B.n(1)):
                for k in range(B.n(2)):
                    for l in range(B.n(3)):
                        T[A.n(0)+i,A.n(1)+j,A.n(2)+k,A.n(3)+l]=B[i,j,k,l]
        return T
    elif A.order()==B.order() and A.order()==5:
        # Initializing the hypermatrix
        T=HM(A.n(0)+B.n(0), A.n(1)+B.n(1), A.n(2)+B.n(2), A.n(3)+B.n(3), A.n(4)+B.n(4), 'zero')
        # Loop filling up the top part of the hypermatrix
        for i in range(A.n(0)):
            for j in range(A.n(1)):
                for k in range(A.n(2)):
                    for l in range(A.n(3)):
                        for m in range(A.n(4)):
                            T[i,j,k,l,m]=A[i,j,k,l,m]
        # Loop filling up the bottom part of the hypermatrix
        for i in range(B.n(0)):
            for j in range(B.n(1)):
                for k in range(B.n(2)):
                    for l in range(B.n(3)):
                        for m in range(B.n(4)):
                            T[A.n(0)+i,A.n(1)+j,A.n(2)+k,A.n(3)+l,A.n(4)+m]=B[i,j,k,l,m]
        return T
    else:
        raise ValueError, "The order of the input hypermatrices must match."

def GeneralHypermatrixKroneckerSum(A,B):
    """
    Computes the  Kronecker sum for third order hypermatrix.

    EXAMPLES:
    ::
        sage: A = HM(2,2,'a'); B = HM(2,2,'b')
        sage: Rslt=GeneralHypermatrixKroneckerSum(A,B)
        [[a00, a01, 0, 0], [a10, a11, 0, 0], [0, 0, b00, b01], [0, 0, b10, b11]]

    AUTHORS:
    - Edinah K. Gnang
    """
    if A.order()==B.order():
        # Initialization of the output hypermatrix
        T=apply(HM,[A.n(i)+B.n(i) for i in range(A.order())]+['zero'])
        # Initialization of the list specifying the dimensions of A
        l = [A.n(i) for i in range(A.order())]
        # Main loop performing the transposition of the entries
        for i in range(prod(l)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entry = [mod(i,l[0])]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            T[tuple(entry)]=A[tuple(entry)]
        # Initialization of the list specifying the dimensions of B
        l = [B.n(i) for i in range(A.order())]
        # Main loop performing the transposition of the entries
        for i in range(prod(l)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entry = [mod(i,l[0])]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            T[tuple((Matrix(ZZ,entry)+Matrix(ZZ,A.dimensions())).list())]=B[tuple(entry)]
        return T
    else:
        raise ValueError, "The order of the input hypermatrices must match."
 
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

def GeneralDualHypermatrixProductB(*args):
    """
    Outputs a list of lists associated with the general
    Dual to the Bhattacharya-Mesner product of the input 
    hypermatrices with non trivial background.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:
    ::
        sage: Ha=HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c')
        sage: Rslt=GeneralDualHypermatrixProductB(Ha,Hb,Hc,HM(2,2,2,'b')); Rslt
        [[[a000*b000^2*c000 + a100*b001*b100*c000 + a001*b000*b010*c001 + a101*b011*b100*c001 + a000*b010*b100*c010 + a100*b101*b110*c010 + a001*b010*b110*c011 + a101*b110*b111*c011, a000*b000*b001*c000 + a100*b001*b101*c000 + a001*b001*b010*c001 + a101*b011*b101*c001 + a000*b011*b100*c010 + a100*b101*b111*c010 + a001*b011*b110*c011 + a101*b111^2*c011], [a010*b000^2*c000 + a110*b001*b100*c000 + a011*b000*b010*c001 + a111*b011*b100*c001 + a010*b010*b100*c010 + a110*b101*b110*c010 + a011*b010*b110*c011 + a111*b110*b111*c011, a010*b000*b001*c000 + a110*b001*b101*c000 + a011*b001*b010*c001 + a111*b011*b101*c001 + a010*b011*b100*c010 + a110*b101*b111*c010 + a011*b011*b110*c011 + a111*b111^2*c011]], [[a000*b000^2*c100 + a100*b001*b100*c100 + a001*b000*b010*c101 + a101*b011*b100*c101 + a000*b010*b100*c110 + a100*b101*b110*c110 + a001*b010*b110*c111 + a101*b110*b111*c111, a000*b000*b001*c100 + a100*b001*b101*c100 + a001*b001*b010*c101 + a101*b011*b101*c101 + a000*b011*b100*c110 + a100*b101*b111*c110 + a001*b011*b110*c111 + a101*b111^2*c111], [a010*b000^2*c100 + a110*b001*b100*c100 + a011*b000*b010*c101 + a111*b011*b100*c101 + a010*b010*b100*c110 + a110*b101*b110*c110 + a011*b010*b110*c111 + a111*b110*b111*c111, a010*b000*b001*c100 + a110*b001*b101*c100 + a011*b001*b010*c101 + a111*b011*b101*c101 + a010*b011*b100*c110 + a110*b101*b111*c110 + a011*b011*b110*c111 + a111*b111^2*c111]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list specifying the dimensions of the output
    l = [(args[i]).n(i) for i in range(len(args)-1)]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the assignement
    # Initializing the background hypermatrix
    B = (args[len(args)-1]).transpose(args[len(args)-1].order()-1)
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
                Rh[tuple(entry)] = Rh[tuple(entry)]+prod([args[s][tuple(entry2[0:Integer(mod(s+1,len(args)))]+[entry[Integer(mod(s+1,len(args)))]]+entry2[Integer(mod(s+2,len(args))):])] for s in range(len(args)-2)]+[args[len(args)-2][tuple(entry2[0:len(entry2)-1]+[entry[len(entry2)-1]])]]+[args[len(args)-1][tuple([entry[0]]+entry2[1:])]])*B[tuple(entry2)]
    return Rh

def PathAdjcencyHypermatrix(A, pthl):
    """
    Procedure for Generating a (k-1)-Path adjacency hypermatrix

    EXAMPLES:
    ::
        sage: M = HM(2,2,'a')
        sage: A = PathAdjcencyHypermatrix(M, k)
        [[[a000*b000,a001*b001],[a010*b010,a011*b011]],[[a100*b100,a101*b101],[a110*b110,a111*b111]]]

    AUTHORS:
    - Edinah K. Gnang
    """
    if A.order() == 2 and pthl == 1:
        return A
    elif A.order() == 2 and pthl > 1:
        # Initializing the number of vertices in the graph.
        sz = max(A.n(0), A.n(1))
        # Initialization of the list specifying the dimensions of the output
        l = [sz for i in range(pthl)] 
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
            Rh[tuple(entry)] = prod([A[tuple([entry[i],entry[i+1]])] for i in range(pthl-1)])
    else:
        raise ValueError, "Input hypermatrix must be order 2 and the path length must be an integer greater then 0"
    return Rh

# program for extending to hypermatrices the background
# approach to the Cayley-Hamilton Theorem.
def  GeneralHypermatrixCayleyHamiltonB(A, t):
    """
    Implements the background hypermatrix approach to the Cayley-Hamilton theorem.

    EXAMPLES:
    ::
        sage: GeneralHypermatrixCayleyHamiltonB(HM(2,2,'a'), 0) 
        [  1   0   0   1]
        [a00 a10 a01 a11]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the hypermatrix order.
    od = A.order()

    # Verifying that the input hypermatrix is cubic
    if len(Set([A.n(i) for i in range(od)]).list()) == 1:
        # initial conditions for the recurrence.
        A0 = GeneralHypermatrixKroneckerDelta(od, A.n(0))
        A1 = A
        # Initializing the list assoaciated with the first two rows of the output matrix
        L = [A0.list(), A1.list()]
        # Loop filling up the remaining lists which make up the rows of the matrix
        for j in range(t):
            # Computing the two next matrices of the recurence.
            A0 = apply(GeneralHypermatrixProductB, [A for k in range(od)]+[A0])
            A1 = apply(GeneralHypermatrixProductB, [A for k in range(od)]+[A1])
            # Append the result to the list
            L.append(A0.list()); L.append(A1.list())
        return Matrix(SR,L)
    else:
        # return the error message if the input hypermatrix is cubic
        raise ValueError, "The input hypermpatrix must be cubic"

def Rank1ApproximationI(A):
    # Initializing the list of variables
    VbL=[var('ln_u'+str(i)) for i in range(A.n(0))]+[var('ln_v'+str(j)) for j in range(A.n(1))]+[var('ln_w'+str(k)) for k in range(A.n(2))]
    # Initialization of the constraints
    CnL=[var('ln_u'+str(i))+var('ln_v'+str(j))+var('ln_w'+str(k))==ln(A[i,j,k]) for i in range(A.n(0)) for j in range(A.n(1)) for k in range(A.n(2))]
    # Initialization of the matrix and right hand side
    [M,b] = ConstraintFormatorII(CnL,VbL)
    # Importing the numpy library
    import numpy
    # computing the least square solution to the problem
    sln = matrix(numpy.linalg.pinv(M))*b
    # Filling up the hypermatrix
    U = HM(A.n(0), 1     , 1     , 'zero')
    V = HM(1     , A.n(1), 1     , 'zero')
    W = HM(1     , 1     , A.n(2), 'zero')    
    for i in range(A.n(0)):
        U[i,0,0]=sln[0 + i, 0]
    for j in range(A.n(1)):
        V[0,j,0]=sln[A.n(0) + j, 0]
    for k in range(A.n(2)):
        W[0,0,k]=sln[A.n(0) + A.n(1) + k, 0]
    return [M,VbL,b,U.elementwise_base_exponent(e),V.elementwise_base_exponent(e),W.elementwise_base_exponent(e)]

def Rank1ApproximationII(A):
    # Initializing the list of variables
    VbL=[var('ln_u'+str(i)+'0'+str(k)) for i in range(A.n(0)) for k in range(A.n(2))]+[var('ln_v'+str(i)+str(j)+'0') for i in range(A.n(0)) for j in range(A.n(1))]+[var('ln_w'+'0'+str(j)+str(k)) for j in range(A.n(1)) for k in range(A.n(2))]
    # Initialization of the constraints
    CnL=[var('ln_u'+str(i)+'0'+str(k))+var('ln_v'+str(i)+str(j)+'0')+var('ln_w'+'0'+str(j)+str(k))==ln(A[i,j,k]) for i in range(A.n(0)) for j in range(A.n(1)) for k in range(A.n(2))]
    [M,b] = ConstraintFormatorII(CnL,VbL)
    # Importing the numpy library
    import numpy
    # computing the least square solution to the problem
    sln = matrix(numpy.linalg.pinv(M))*b
    # Filling up the hypermatrix
    U = HM(A.n(0), 1     , A.n(2), 'zero')
    V = HM(A.n(0), A.n(1), 1     , 'zero')
    W = HM(1     , A.n(1), A.n(2), 'zero')
    for i in range(A.n(0)):
        for k in range(A.n(2)):
            U[i,0,k]=sln[0 + A.n(2)*i+k,0]
    for i in range(A.n(0)):
        for j in range(A.n(1)):
            V[i,j,0]=sln[A.n(0)*A.n(2) + A.n(1)*i+j,0]
    for j in range(A.n(1)):
        for k in range(A.n(2)):
            W[0,j,k]=sln[A.n(0)*A.n(2) + A.n(0)*A.n(1) + A.n(2)*j+k,0]
    return [M,VbL,b,U.elementwise_base_exponent(e),V.elementwise_base_exponent(e),W.elementwise_base_exponent(e)]

def Rank1MulMtrxApprox(A):
    # Initializing the list of variables
    VbL=[var('ln_u'+'0'+str(i)) for i in range(A.n(0))]+[var('ln_v'+str(j)+'0') for j in range(A.n(1))]
    # Initialization of the constraints
    CnL=[var('ln_u'+'0'+str(i))+var('ln_v'+str(j)+'0')==ln(A[i,j]) for i in range(A.n(0)) for j in range(A.n(1))]
    # Initializing the constraints
    [M, b] = ConstraintFormatorII(CnL, VbL)
    # Importing the numpy library
    import numpy
    # computing the least square solution to the problem
    sln = matrix(numpy.linalg.pinv(M))*b
    # Filling up the hypermatrix
    U = HM(A.n(0), 1     , 'zero')
    V = HM(1     , A.n(1), 'zero')
    for i in range(A.n(0)):
        U[i,0]=sln[0 + i, 0]
    for j in range(A.n(1)):
        V[0,j]=sln[A.n(0) + j, 0]
    return [M,VbL,b,U.elementwise_base_exponent(e), V.elementwise_base_exponent(e)]

def fast_reduce(f, monom, subst):
    """
    computes the reduction by monomial substitution
    
    EXAMPLES:
 
    ::  
        sage: fast_reduce(x1^3+x2+x3^3,[x3^3],[1])
        x1^3+x2+1

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    if len(monom) == len(subst):
        s = str(f)
        for i in range(len(monom)):
            s = s.replace(str(monom[i]), str(subst[i]))
        return expand((SR(s)).simplify_full())
    else:
        print 'Error the monomial list and the substitution list must have the same length'

def ThirdOrderHyperdeterminant(H):
    """
    computes third order hypermatrix determinant
    
    EXAMPLES:
 
    ::  
        sage: A = HM(2,2,2,'a')
        sage: ThirdOrderHyperdeterminant(A)
        -a000*a011*a101*a110 + a001*a010*a100*a111

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Testing to see that the hypermatrix is indeed a cube
    if len(Set(H.dimensions()).list())==1 and H.order()==3:
        # Initializing the matrix for the mnemonic construction
        A = Matrix(SR, H.n(0), H.n(0), [var('x'+str(i)+str(j)) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0))])
        # Computing the mnemonique polynomial
        P = Permutations(range(A.nrows()))
        L = expand(sum([Permutation([p[i]+1 for i in range(len(p))]).signature()*prod([ (A[k][p[k]]) for k in range(A.nrows())]) for p in P])*prod([sum([prod([ (A[k][p[k]])^j for k in range(A.nrows())]) for p in P]) for j in range(2,H.n(0)+1)])).operands()
        # Computing the polynomial
        f = sum([l for l in L if len((l^2).operands())==(H.n(0))^2])
        # Loop performing the umbral expression
        for k in range(H.n(0),0,-1):
            f = fast_reduce(f, [var('x'+str(i)+str(j))^k for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0))], [var('a'+str(i)+str(j)+str(k)) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0))])
        return f.subs(dict([(var('a'+str(i)+str(j)+str(k)),H[i-1,j-1,k-1]) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0))]))
    else :
        # Print an error message indicating that the matrix must be a cube.
        raise ValueError, "The hypermatrix must be a third order cube hypermatrix."

def FourthOrderHyperdeterminant(H):
    """
    computes third order hypermatrix determinant
    
    EXAMPLES:
 
    ::  
        sage: A = HM(2,2,2,2,'a')
        sage: FourthOrderHyperdeterminant(A)
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Testing to see that the hypermatrix is indeed a cube
    if len(Set(H.dimensions()).list())==1 and H.order()==4:
        # Initializing the matrix for the mnemonic construction
        A=HM(H.n(0),H.n(0),H.n(0), [var('x'+str(i)+str(j)+str(k)) for k in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for i in range(1,1+H.n(0))])
        # Computing the polynomial
        L=expand(ThirdOrderHyperdeterminant(A)*prod([sum([sqrt(g^2).simplify_radical() for g in ThirdOrderHyperdeterminant(A.elementwise_exponent(j)).operands()]) for j in range(2,1+A.n(0))])).operands()
        # Computing the polynomial
        f = sum([l for l in L if len((l^2).operands())==(H.n(0))^3])
        # Loop performing the umbral expression
        for l in range(H.n(0),0,-1):
            f = fast_reduce(f,[var('x'+str(i)+str(j)+str(k))^l for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0))],[var('a'+str(i)+str(j)+str(k)+str(l)) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0))])
        return f.subs(dict([(var('a'+str(i)+str(j)+str(k)+str(l)),H[i-1,j-1,k-1,l-1]) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for l in range(1,1+H.n(0))]))
    else :
        # Print an error message indicating that the matrix must be a cube.
        raise ValueError, "The hypermatrix must be a fourth order hypercube hypermatrix."

def FifthOrderHyperdeterminant(H):
    """
    computes third order hypermatrix determinant
    
    EXAMPLES:
 
    ::  
        sage: A = HM(2,2,2,2,2,'a')
        sage: FifthOrderHyperdeterminant(A)
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Testing to see that the hypermatrix is indeed a cube
    if len(Set(H.dimensions()).list())==1 and H.order()==5:
        # Initializing the matrix for the mnemonic construction
        A=HM(H.n(0),H.n(0),H.n(0),H.n(0),[var('x'+str(i)+str(j)+str(k)+str(l)) for l in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for i in range(1,1+H.n(0))])
        # Computing the polynomial
        L = expand(FourthOrderHyperdeterminant(A)*prod([sum([sqrt(g^2).simplify_radical() for g in FourthOrderHyperdeterminant(A.elementwise_exponent(j)).operands()]) for j in range(2,1+A.n(0))])).operands()
        f = sum([l for l in L if len((l^2).operands())==(H.n(0))^4])
        # Loop performing the umbral expression
        for m in range(H.n(0),0,-1):
            f = fast_reduce(f,[var('x'+str(i)+str(j)+str(k)+str(l))^m for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for l in range(1,1+H.n(0))],[var('a'+str(i)+str(j)+str(k)+str(l)+str(m)) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for l in range(1,1+H.n(0))])
        return f.subs(dict([(var('a'+str(i)+str(j)+str(k)+str(l)+str(m)),H[i-1,j-1,k-1,l-1,m-1]) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for l in range(1,1+H.n(0)) for m in range(1,1+H.n(0))]))
    else :
        # Print an error message indicating that the matrix must be a cube.
        raise ValueError, "The hypermatrix must be a fifth order hypercube hypermatrix."

def SixthOrderHyperdeterminant(H):
    """
    computes third order hypermatrix determinant
    
    EXAMPLES:
 
    ::  
        sage: A = HM(2,2,2,2,2,2,'a')
        sage: SixthOrderHyperdeterminant(A)
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Testing to see that the hypermatrix is indeed a cube
    if len(Set(H.dimensions()).list())==1 and H.order()==6:
        # Initializing the matrix for the mnemonic construction
        A=HM(H.n(0),H.n(0),H.n(0),H.n(0),H.n(0),[var('x'+str(i)+str(j)+str(k)+str(l)+str(m)) for m in range(1,1+H.n(0)) for l in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for i in range(1,1+H.n(0))])
        # Computing the polynomial
        L = expand(FifthOrderHyperdeterminant(A)*prod([sum([sqrt(g^2).simplify_radical() for g in FifthOrderHyperdeterminant(A.elementwise_exponent(j)).operands()]) for j in range(2,1+A.n(0))])).operands()
        f = sum([l for l in L if len((l^2).operands())==(H.n(0))^5])
        # Loop performing the umbral expression
        for p in range(H.n(0),0,-1):
            f = fast_reduce(f,[var('x'+str(i)+str(j)+str(k)+str(l)+str(m))^p for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for l in range(1,1+H.n(0)) for m in range(1,1+H.n(0))],[var('a'+str(i)+str(j)+str(k)+str(l)+str(m)+str(p)) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for l in range(1,1+H.n(0)) for m in range(1,1+H.n(0))])
        return f.subs(dict([(var('a'+str(i)+str(j)+str(k)+str(l)+str(m)+str(p)),H[i-1,j-1,k-1,l-1,m-1,p-1]) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for l in range(1,1+H.n(0)) for m in range(1,1+H.n(0)) for p in range(1,1+H.n(0))]))
    else :
        # Print an error message indicating that the matrix must be a cube.
        raise ValueError, "The hypermatrix must be a sixth order hypercube hypermatrix."

def GeneralSidelength2Hyperdeterminant(od,c):
    """
    outputs the symbolic expression with the determinant of hypermatrices of arbitrary orders.
    but every size of the hypermatrix must be equal to two. It ouputs an equality derived via
    the rank one argument.

    EXAMPLES:
 
    ::  
        sage: GeneralSidelength2Hyperdeterminant(2,'a'):
        a00*a11-a01*a10

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization the alphabet list
    AlphaB=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    # Initializing the list of hypermatrices
    L=[apply(HM,[2 for i in range(j)]+[1]+[2 for i in range(j+1,od)]+[AlphaB[j]]).elementwise_base_exponent(e) for j in range(1,od)+[0]]
    # Initilizing the list of variable
    VrbLst=[]
    for Q in L:
        VrbLst = VrbLst + (Q.elementwise_base_logarithm(e)).list()
    Eq=apply(GeneralHypermatrixProduct, [Q for Q in L])+apply(HM,[2 for i in range(od)]+[c])
    # Writting up the constraints
    Le=Eq.list()
    # Filling up the linear constraints
    CnstrLst=[] 
    for f in Le:
        CnstrLst.append(ln((f.operands())[1]).canonicalize_radical() == ln((f.operands())[0]).canonicalize_radical())
    [Mtr, b]=ConstraintFormatorII(CnstrLst, VrbLst)
    # Returning the determinantal equality
    return exp((Matrix(SR, ((Mtr.kernel()).basis()[0]).list())*b)[0,0]).canonicalize_radical().numerator()-exp((Matrix(SR, ((Mtr.kernel()).basis()[0]).list())*b)[0,0]).canonicalize_radical().denominator()

def HypermatrixDodgsonCondensation(A):
    if len(Set(A.dimensions()).list())==1:
        if A.n(0)==2:
            # Computing the determinant
            f=GeneralSidelength2Hyperdeterminant(A.order(),'xi')
            La=A.list(); Lb=apply(HM,A.dimensions()+['xi']).list()
            return(f.subs(dict([(Lb[i],La[i]) for i in range(len(La))])))
        elif A.n(0)>2:
            # Initializing the copy
            Bold=A.copy()
            Bnew=apply(HM,[i-1 for i in Bold.dimensions()]+['zero'])
            while Bold.n(0)>2:
                l=Bnew.dimensions()
                # Main loop performing the transposition of the entries
                for i in range(prod(l)):
                    # Turning the index i into an hypermatrix array location using the decimal encoding trick
                    entry = [mod(i,l[0])]
                    sm = Integer(mod(i,l[0]))
                    for k in range(len(l)-1):
                        entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                        sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                    # Initialization of the summand
                    Bl=[]
                    l2=[2 for j in range(A.order())]
                    for j in range(prod(l2)):
                        ent=[mod(j,l2[0])]
                        ms=Integer(mod(j,l2[0]))
                        for t in range(len(l2)-1):
                            ent.append(Integer(mod(Integer((j-ms)/prod(l2[0:t+1])),l2[t+1])))
                            ms = ms+prod(l2[0:t+1])*ent[len(ent)-1]
                        Bl.append((Matrix(ZZ,entry)+Matrix(ZZ,ent)).list())
                    Bnew[tuple(entry)]=HypermatrixDodgsonCondensation(apply(HM,[2 for j in range(A.order())]+[[Bold[tuple(entry2)] for entry2 in Bl ]]))
                # Performing the update
                Bold=Bnew.copy()
                Bnew=apply(HM,[i-1 for i in Bold.dimensions()]+['zero'])
            return HypermatrixDodgsonCondensation(Bold)
        else:
            raise ValueError, "The input hypermatrix side length greater then 1."
    else:
        raise ValueError, "The input hypermatrix must be hypercubic"

def GeneralHypermatrixRank1Parametrization(sz,od):
    """
    Outputs the symbolic constraints associated with rank one parametrization
    of hypermatrices.

    EXAMPLES:
 
    ::  
        sage: GeneralHypermatrixRank1Parametrization(2,2) 
        [1 0 1 0]
        [0 1 1 0]
        [1 0 0 1]
        [0 1 0 1]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization the alphabet list
    AlphaB=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    # Initializing the list of hypermatrices
    L=[apply(HM,[sz for i in range(j)]+[1]+[sz for i in range(j+1,od)]+[AlphaB[j]]) for j in range(1,od)+[0]]
    # Initilizing the list of variable
    VrbLst=[]
    for Q in L:
        VrbLst = VrbLst+Q.list()
    Eq=apply(GeneralHypermatrixLogProduct, [Q for Q in L])
    CnstrLst=[eq==0 for eq in Eq.list()]
    return ConstraintFormatorII(CnstrLst, VrbLst)[0]

def SecondOrderHadamardBlockU(l):
    """
    outputs the  direct sum construction of hadamard block matrices.
    the vectors of the matrices are not normalized.

    EXAMPLES:
 
    ::  
        sage: SecondOrderHadamardBlockU(2)
        [[1, 1], [1, -1]]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    bns = l.str(2)
    Szl = [(2^(len(bns)-1-i),bns[i]) for i in range(len(bns)) if bns[i]!='0']
    H = Matrix(QQ,hadamard_matrix(Szl[0][0]))
    for i in range(1,len(Szl)):
        if Szl[i][0]==1:
            H = H.block_sum(Matrix(QQ,1,1,[1]))
        else:
            H = H.block_sum(Matrix(QQ,hadamard_matrix(Szl[i][0])))
    return HM(l,l,H.list())

def SecondOrderHadamardBlock(l):
    """
    outputs the  direct sum construction of hadamard block matrices.
    the vectors of the matrices are not normalized.

    EXAMPLES:
 
    ::  
        sage: SecondOrderHadamardBlock(2)
        [[1/2*sqrt(2), 1/2*sqrt(2)], [1/2*sqrt(2), -1/2*sqrt(2)]]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    bns = l.str(2)
    Szl = [(2^(len(bns)-1-i),bns[i]) for i in range(len(bns)) if bns[i]!='0']
    H = (1/sqrt(Szl[0][0]))*Matrix(QQ,hadamard_matrix(Szl[0][0]))
    for i in range(1,len(Szl)):
        if Szl[i][0]==1:
            H = H.block_sum(Matrix(QQ,1,1,[1]))
        else:
            H = H.block_sum((1/sqrt(Szl[i][0]))*Matrix(QQ,hadamard_matrix(Szl[i][0])))
    return HM(l,l,H.list())

def ThirdOrderHadamardBlockU(l):
    """
    outputs the  direct sum construction of hadamard block hypermatrices.
    the vectors of the matrices are not normalized.

    EXAMPLES:
 
    ::  
        sage: ThirdOrderHadamardBlockU(2)
        [[[1, 1], [1, 1]], [[-1, 1], [1, 1]]]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing the 2x2x2 hadamard Hypermatrix.
    Hd = HM([[[1, 1], [1, 1]], [[-1, 1], [1, 1]]]) 
    # Obtaining the binary encoding of the input integer
    bns = l.str(2)
    Szl = [len(bns)-1-i for i in range(len(bns)) if bns[i] != '0']
    Lh = [HypermatrixSliceKroneckerPower(Hd, i) for i in Szl if i > 0]
    H = Lh[0]
    if Integer(mod(l,2)) == 0:
        for i in range(1,len(Lh)):
            H = H.block_sum(Lh[i])
    else :
        for i in range(1,len(Lh)):
            H = H.block_sum(Lh[i])
        H = H.block_sum(HM(1,1,1,'one')) 
    return H

def ThirdOrderHadamardBlock(l):
    """
    outputs the  direct sum construction of hadamard block hypermatrices.
    the vectors of the matrices are not normalized.

    EXAMPLES:
 
    ::  
        sage: ThirdOrderHadamardBlock(2)
        [[[(1/2)^(1/3), (1/2)^(1/3)], [(1/2)^(1/3), (1/2)^(1/3)]], [[-(1/2)^(1/3), (1/2)^(1/3)], [(1/2)^(1/3), (1/2)^(1/3)]]]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing the 2x2x2 hadamard Hypermatrix.
    Hd = (1/2)^(1/3)*HM([[[1, 1], [1, 1]], [[-1, 1], [1, 1]]]) 
    # Obtaining the binary encoding of the input integer
    bns = l.str(2)
    Szl = [len(bns)-1-i for i in range(len(bns)) if bns[i] != '0']
    Lh = [HypermatrixSliceKroneckerPower(Hd, i) for i in Szl if i > 0]
    H = Lh[0]
    if Integer(mod(l,2)) == 0:
        for i in range(1,len(Lh)):
            H = H.block_sum(Lh[i])
    else :
        for i in range(1,len(Lh)):
            H = H.block_sum(Lh[i])
        H = H.block_sum(HM(1,1,1,'one')) 
    return H

# Fifth order Hadamard block.
def FifthOrderHadamardBlockU(l):
    """
    outputs the  direct sum construction of hadamard block hypermatrices.
    the vectors of the matrices are not normalized.

    EXAMPLES:
 
    ::  
        sage: FifthOrderHadamardBlockU(2)
        [[[[[1,1],[1,1]],[[1,1],[1,1]]],[[[1,1],[1,1]],[[1,1],[1,1]]]],[[[[1,1],[-1,-1]],[[1,1],[1,1]]],[[[-1,1],[1,1]],[[1,1],[1,1]]]]]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing the 2x2x2x2x2 hadamard hypermatrix.
    Hd=HM([[[[[1,1],[1,1]],[[1,1],[1,1]]],[[[1,1],[1,1]],[[1,1],[1,1]]]],[[[[1,1],[-1,-1]],[[1,1],[1,1]]],[[[-1,1],[1,1]],[[1,1],[1,1]]]]])
    # Obtaining the binary encoding of the input integer
    bns = l.str(2)
    Szl = [len(bns)-1-i for i in range(len(bns)) if bns[i] != '0']
    Lh = [HypermatrixSliceKroneckerPower(Hd, i) for i in Szl if i > 0]
    H = Lh[0]
    if Integer(mod(l,2))==0:
        for i in range(1,len(Lh)):
            H = H.block_sum(Lh[i])
    else :
        for i in range(1,len(Lh)):
            H = H.block_sum(Lh[i])
        H = H.block_sum(HM(1,1,1,1,1,'one')) 
    return H

def FifthOrderHadamardBlock(l):
    """
    outputs the  direct sum construction of hadamard block hypermatrices.
    the vectors of the matrices are not normalized.

    EXAMPLES:
 
    ::  
        sage: FifthOrderHadamardBlock(2)

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing the 2x2x2x2x2 hadamard Hypermatrix.
    Hd = (1/2)^(1/5)*HM([[[[[1,1],[1,1]],[[1,1],[1,1]]],[[[1,1],[1,1]],[[1,1],[1,1]]]],[[[[1,1],[-1,-1]],[[1,1],[1,1]]],[[[-1,1],[1,1]],[[1,1],[1,1]]]]])
    # Obtaining the binary encoding of the input integer
    bns = l.str(2)
    Szl = [len(bns)-1-i for i in range(len(bns)) if bns[i] != '0']
    Lh = [HypermatrixSliceKroneckerPower(Hd, i) for i in Szl if i > 0]
    H = Lh[0]
    if Integer(mod(l,2)) == 0:
        for i in range(1,len(Lh)):
            H = H.block_sum(Lh[i])
    else :
        for i in range(1,len(Lh)):
            H = H.block_sum(Lh[i])
        H = H.block_sum(HM(1,1,1,1,1,'one')) 
    return H

def RandomTransposition(n):
    """
    outputs  the random transposition permutation of n elements.
    The fucntion is used for randominzing the Hadamard block construction
    

    EXAMPLES:
 
    ::  
        sage: RandomTransposition(3)
        [1, 0, 2] 

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the vector
    p = [0 .. n-1]
    idx = ZZ.random_element(n)
    jdx = idx
    # avoiding the trivial transposition
    while jdx == idx:
        jdx = ZZ.random_element(n)
    # recording the transposition
    tmp = p[idx]; p[idx] = p[jdx]; p[jdx] = tmp
    return p

def RandomColumnSlicePermutation(H):
    """
    Outputs a permutation of the slices of the input third order hypermatrix.

    EXAMPLES:
 
    ::  
        sage: RandomColumnSlicePermutation(HM(2,2,2,'a'))
        [[[a010, a011], [a000, a001]], [[a110, a111], [a100, a101]]] 
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    sz = H.n(2)
    for i in range(Integer(ceil(sqrt(sz)))+1):
        # Initialization of the permutation
        P = HypermatrixPermutation(RandomTransposition(sz))
        H = Prod(H, HM(P), HM(P).transpose())
    return H

def ThirdOrderHypermatrixResolutionPartition(U, V, W, Ha, Hb, Hc, NbPrts=2):
    """
    outputs the spliting of a third order hypermatrix into the pieces
    as suggested by the resolution of identity. The first three input
    hypermatrices are uncorrelated tuples.
    the last three imputs correspond to the factors for the spliting.

    EXAMPLES:
 
    ::  
        sage: [U,V,W]=GeneralUncorrelatedHypermatrixTuple(3)
        sage: L=ThirdOrderHypermatrixResolutionPartition(U,V,W,HM(2,2,2,'a'),HM(2,2,2,'b'),HM(2,2,2,'c'))
        sage: len(L)
        2
        sage: sum(L).simplify()
        [[[a000*b000*c000 + a010*b001*c100, a001*b000*c001 + a011*b001*c101], [a000*b010*c010 + a010*b011*c110, a001*b010*c011 + a011*b011*c111]], [[a100*b100*c000 + a110*b101*c100, a101*b100*c001 + a111*b101*c101], [a100*b110*c010 + a110*b111*c110, a101*b110*c011 + a111*b111*c111]]]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    L = []
    for i in range(U.n(1)):
        # Filling up the slice
        Tp0 = HM(U.n(0), 1, U.n(2), 'zero')
        for u in range(Tp0.n(0)):
            for w in range(Tp0.n(2)):
                Tp0[u,0,w] = U[u,i,w]
        # Filling up the slice
        Tp1 = HM(V.n(0), V.n(1), 1, 'zero')
        for u in range(Tp1.n(0)):
            for v in range(Tp1.n(1)):
                Tp1[u,v,0] = V[u,v,i]
        # Filling up the slice
        Tp2 = HM(1, W.n(1), W.n(2), 'zero')
        for v in range(Tp2.n(1)):
            for w in range(Tp2.n(2)):
                Tp2[0,v,w] = W[i,v,w]
        # Appending the components to the list
        #L.append(ProdB(Ha, Hb, Hc, Prod(Tp0, Tp1, Tp2)))
        L.append(Prod(Tp0, Tp1, Tp2))
    #return [sum(L[j*NbPrts:min((j+1)*NbPrts,len(L))]) for j in range(ceil(len(L)/NbPrts))]
    return [ProdB(Ha,Hb,Hc,sum(L[j*NbPrts:min((j+1)*NbPrts,len(L))])) for j in range(ceil(len(L)/NbPrts))]

def FourthOrderHypermatrixResolutionPartition(Q, U, V, W, Ha, Hb, Hc, Hd):
    """
    outputs the spliting of a fourth order hypermatrix into the pieces
    as suggested by the resolution of identity. The first three input
    hypermatrices are uncorrelated tuples.
    the last three imputs correspond to the factors for the spliting.

    EXAMPLES:
 
    ::  
        sage: [Q, U, V, W] = GeneralUncorrelatedHypermatrixTuple(3)
        sage: L = FourthOrderHypermatrixResolutionPartition(Q,U,V,W,HM(2,2,2,2,'a'),HM(2,2,2,2,'b'),HM(2,2,2,2,'c'), HM(2,2,2,2,'d'))
        sage: len(L)
        2
        sage: sum(L).simplify()
        [[[[a0000*b0000*c0000*d0000 + a0100*b0010*c0001*d1000, a0001*b0001*c0000*d0001 + a0101*b0011*c0001*d1001], [a0010*b0000*c0010*d0010 + a0110*b0010*c0011*d1010, a0011*b0001*c0010*d0011 + a0111*b0011*c0011*d1011]], [[a0000*b0100*c0100*d0100 + a0100*b0110*c0101*d1100, a0001*b0101*c0100*d0101 + a0101*b0111*c0101*d1101], [a0010*b0100*c0110*d0110 + a0110*b0110*c0111*d1110, a0011*b0101*c0110*d0111 + a0111*b0111*c0111*d1111]]], [[[a1000*b1000*c1000*d0000 + a1100*b1010*c1001*d1000, a1001*b1001*c1000*d0001 + a1101*b1011*c1001*d1001], [a1010*b1000*c1010*d0010 + a1110*b1010*c1011*d1010, a1011*b1001*c1010*d0011 + a1111*b1011*c1011*d1011]], [[a1000*b1100*c1100*d0100 + a1100*b1110*c1101*d1100, a1001*b1101*c1100*d0101 + a1101*b1111*c1101*d1101], [a1010*b1100*c1110*d0110 + a1110*b1110*c1111*d1110, a1011*b1101*c1110*d0111 + a1111*b1111*c1111*d1111]]]]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    L = []
    for i in range(Q.n(1)):
        # Filling up the slice
        Tp0 = HM(Q.n(0), 1, Q.n(2), Q.n(3), 'zero')
        for q in range(Tp0.n(0)):
            for v in range(Tp0.n(2)):
                for w in range(Tp0.n(3)):
                    Tp0[q,0,v,w] = Q[q,i,v,w]
        # Filling up the slice
        Tp1 = HM(U.n(0), U.n(1), 1, U.n(3), 'zero')
        for q in range(Tp1.n(0)):
            for u in range(Tp1.n(1)):
                for w in range(Tp1.n(3)):
                    Tp1[q,u,0,w] = U[q,u,i,w]
        # Filling up the slice
        Tp2 = HM(V.n(0), V.n(1), V.n(2), 1, 'zero')
        for q in range(Tp2.n(0)):
            for u in range(Tp2.n(1)):
                for v in range(Tp2.n(2)):
                    Tp2[q,u,v,0] = V[q,u,v,i]
        # Filling up the slice
        Tp3 = HM(1, W.n(1), W.n(2), W.n(3), 'zero')
        for u in range(Tp3.n(1)):
            for v in range(Tp3.n(2)):
                for w in range(Tp3.n(3)):
                    Tp3[0,u,v,w] = W[i,u,v,w]
        # Appending the components to the list
        L.append(ProdB(Ha, Hb, Hc, Hd, Prod(Tp0, Tp1, Tp2, Tp3)))
    return L

def CanonicalThirdOrderHypermatrixFactors(A):
    """
    Outputs the canonical factorization for third order hypermatrices. 
    The canonical factorization results from the fact that any sum of
    outer products is a hypermatrix product.

    EXAMPLES:
 
    ::  
        sage: [Q,U,V] = CanonicalThirdOrderHypermatrixFactors(HM(2,2,2,'a')):
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the size parameters
    sz = A.dimensions()
    # Initializing the list of column vectors
    e_i=[identity_matrix(sz[0])[:,i] for i in range(sz[0])]
    e_j=[identity_matrix(sz[1])[:,j] for j in range(sz[1])]
    e_k=[identity_matrix(sz[2])[:,k] for k in range(sz[2])]
    # Initializing the first hypermatrix
    T0 = HM(sz[0],prod(sz),sz[2],'zero')
    # Initializing the second hypermatrix
    T1 = HM(sz[0],sz[1],prod(sz),'zero')
    # Filling up of the third hypermatrix
    T2 = HM(prod(sz),sz[1],sz[2],'zero')
    cnt = 0
    for i in range(len(e_i)):
        for j in range(len(e_j)):
            for k in range(len(e_k)):
                M0 = e_i[i]*e_k[k].transpose()
                M1 = e_i[i]*e_j[j].transpose()
                M2 = e_j[j]*e_k[k].transpose()
                for r in range(sz[0]):
                    for d in range(sz[2]):
                        T0[r,cnt,d]=(A[i,j,k])^(1/3)*M0[r,d]
                for r in range(sz[0]):
                    for c in range(sz[1]):
                        T1[r,c,cnt]=(A[i,j,k])^(1/3)*M1[r,c]
                for c in range(sz[1]):
                    for d in range(sz[2]):
                        T2[cnt,c,d]=(A[i,j,k])^(1/3)*M2[c,d]
                cnt = cnt+1
    return [T0,T1,T2]

def GeneralHypermatrixCompositionMatrix(*args):
    """
    Outputs the matrix constraints associated with the product 
    this came up when I tried to use system of linear algebra
    and resoltion of identity to express a composition rule.
    This was not overall successfull.

    EXAMPLES:
 
    ::  
        sage: [Eqnts,Mtr]=GeneralHypermatrixCompositionMatrix(HM(2,1,'a'),HM(1,2,'b'))
        sage: Eqnts
        [a00+b00==0,a10+b00==0,a00+b01==0,a10+b01==0]
        sage: Mtr
        [1 0 1 0]
        [0 1 1 0]
        [1 0 0 1]
        [0 1 0 1]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the list specifying the dimensions of the output
    l = [(args[i]).n(i) for i in range(len(args))]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = []
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
            Rh = Rh + [sum([args[s][tuple(entry[0:Integer(mod(s+1,len(args)))]+[t]+entry[Integer(mod(s+2,len(args))):])] for s in range(len(args)-2)] + [args[len(args)-2][tuple(entry[0:len(args)-1]+[t])]]+[args[len(args)-1][tuple([t]+entry[1:])]]) == 0 for t in range((args[0]).n(1))]
    # Initialization of the variables.
    Vrbls = []
    for i in range(len(args)):
        Vrbls = Vrbls + (args[i]).list()
    return [Rh, ConstraintFormatorII(Rh, Vrbls)[0]]

def GeneralHypermatrixLogProductTermList(*args):
    """
    Outputs the matrix constraints associated with the Logproduct 
    this came up when I tried to use system of linear algebra
    and resoltion of identity to express a composition rule.
    This was not overall successfull.

    EXAMPLES:
 
    ::  
        sage: [Eqnts,Mtr]=GeneralHypermatrixCompositionMatrix(HM(2,1,'a'),HM(1,2,'b'))
        sage: Eqnts
        [a00+b00==0,a10+b00==0,a00+b01==0,a10+b01==0]
        sage: Mtr
        [1 0 1 0]
        [0 1 1 0]
        [1 0 0 1]
        [0 1 0 1]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the list specifying the dimensions of the output
    l = [(args[i]).n(i) for i in range(len(args))]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = []
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
            Rh = Rh + [sum([args[s][tuple(entry[0:Integer(mod(s+1,len(args)))]+[t]+entry[Integer(mod(s+2,len(args))):])] for s in range(len(args)-2)] + [args[len(args)-2][tuple(entry[0:len(args)-1]+[t])]]+[args[len(args)-1][tuple([t]+entry[1:])]]) for t in range((args[0]).n(1))]
    return Rh

def gaussian_elimination(Cf, rs):
    """
    Outputs the row echelon form of the input matrix and the right hand side.

    EXAMPLES:
 
    ::  
        sage: [RefA, c] = gaussian_elimination(Matrix(SR,HM(2,2,'a').listHM()), Matrix(SR,HM(2,1,'b').listHM()))
        sage: RefA
        [      1 a01/a00]
        [      0       1]
        sage: c
        [                                b00/a00]
        [(a10*b00/a00 - b10)/(a01*a10/a00 - a11)]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    A=copy(Cf); b=copy(rs)
    """
    # Zero padding the matrix if necessary.
    if A.nrows() < A.ncols():
        # Temporary matrix for A
        Ta = Matrix(SR, zero_matrix(A.ncols(), A.ncols()))
        Ta[:A.nrows(), :A.ncols()] = copy(A)
        # Temporary matrix for b
        Tb = Matrix(SR, zero_matrix(A.ncols(), b.ncols())) 
        Tb[:b.nrows(), :b.ncols()] = copy(b)
        # replacing the matrix with the zero apdding.
        A = Ta; b = Tb
    """
    # Initialization of the row and column index
    i = 0; j = 0;
    while i < A.nrows() and j < A.ncols():
        while (A[i:,j]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        if (A[i:,:].is_zero()) == False:
            while (A[i,j]).is_zero(): 
                Ta=A[i:,:]
                Tb=b[i:,:]
                # Initializing the cyclic shift permutation matrix
                Id=identity_matrix(Ta.nrows())
                P=sum([Id[:,k]*Id[mod(k+1,Ta.nrows()),:] for k in range(Ta.nrows())])
                Ta=P*Ta; Tb=P*Tb
                A[i:,:]=Ta
                b[i:,:]=Tb
            # Performing the row operations.
            b[i,:]=(1/A[i,j])*b[i,:]
            A[i,:]=(1/A[i,j])*A[i,:]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if A[r,:].is_zero():
                    r=r+1
                else:
                    b[r,:]=-A[r,j]*b[i,:]+b[r,:]
                    A[r,:]=-A[r,j]*A[i,:]+A[r,:]
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return [A, b]

def gauss_jordan_elimination(Cf,rs):
    """
    Outputs the reduced row echelon form of the input matrix and the right hand side.

    EXAMPLES:
 
    ::  
        sage: [RefA, c] = gauss_jordan_elimination(Matrix(SR,HM(2,2,'a').listHM()), Matrix(SR,HM(2,1,'b').listHM()))
        sage: RefA
        [1 0]
        [0 1]
        sage: c
        [-a01*(a10*b00/a00 - b10)/(a00*(a01*a10/a00 - a11)) + b00/a00]
        [                     (a10*b00/a00 - b10)/(a01*a10/a00 - a11)]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    [A, b] = gaussian_elimination(Cf,rs)
    # Initialization of the row and column index
    i=A.nrows()-1; j=0
    while i>0 or j>0:
        if (A[i,:]).is_zero():
            # decrementing the row index and initializing the column index
            i=i-1; j=0
        else :
            while (A[i,j]).is_zero():
                # Incrementing the column index
                j = j + 1
            # performing row operations
            for r in range(i-1,-1,-1):
                b[r,:] = -A[r,j]*b[i,:]+b[r,:]
                A[r,:] = -A[r,j]*A[i,:]+A[r,:]
            i=i-1; j=0
    return [A,b]

def linear_solver(A,b,c):
    """
    Outputs the reduced row echelon form of the input matrix and the right hand side.

    EXAMPLES:
 
    ::  
        sage: sz=2; Eq=[var('x'+str(i))+var('x'+str(sz+j))==var('a'+str(i)+str(j)) for i in range(sz) for j in range(sz)]
        sage: [A,b]=ConstraintFormatorII(Eq,[var('x'+str(i)) for i in range(2*sz)])
        sage: [Ap,bp]=gauss_jordan_elimination(A.transpose()*A,A.transpose()*b)
        sage: linear_solver(Ap,bp,'t')
        [ 1/4*a00 + 3/4*a01 - 1/4*a10 + 1/4*a11 - t30]
        [-1/4*a00 + 1/4*a01 + 1/4*a10 + 3/4*a11 - t30]
        [ 1/2*a00 - 1/2*a01 + 1/2*a10 - 1/2*a11 + t30]
        [                                           0]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the reduced echelon form.
    [Ap,bp]=gauss_jordan_elimination(A,b)
    # Filling up the matrix to a square.
    B=Matrix(SR,zero_matrix(max(A.nrows(),A.ncols()),max(Ap.nrows(),Ap.ncols())))
    B[0:A.nrows(),0:A.ncols()]=Ap
    # returning the result
    return bp+((B.elementwise_product(identity_matrix(B.nrows()))-B)*Matrix(HM(B.nrows(),1,[var(c+str(k)) for k in range(B.nrows())]).listHM()))[0:A.nrows()]

def multiplicative_gaussian_elimination(Cf,rs,jndx=0):
    """
    Outputs the row echelon form of the input matrix and the right hand side.

    EXAMPLES:
 
    ::  
        sage: [EfA,c,indx,Lst]=multiplicative_gaussian_elimination(Matrix(SR,HM(2,2,'a').listHM()), Matrix(SR,HM(2,1,'b').listHM()))
        sage: EfA
        [1             a01/a00]
        [0                   1]
        sage: c
        [                                                                    (b00*e^(2*I*pi*k0))^(1/a00)]
        [(((b00*e^(2*I*pi*k0))^(1/a00)*e^(2*I*pi*k1))^(-a10)*b10*e^(2*I*pi*k2))^(-1/(a01*a10/a00 - a11))]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    A = copy(Cf); b = copy(rs)
    """
    # Zero padding the matrix if necessary.
    if A.nrows()<A.ncols():
        # Temporary matrix for A
        Ta = Matrix(SR, zero_matrix(A.ncols(), A.ncols()))
        Ta[:A.nrows(), :A.ncols()] = copy(A)
        # Temporary matrix for b
        Tb = Matrix(SR, zero_matrix(A.ncols(), b.ncols())) 
        Tb[:b.nrows(), :b.ncols()] = copy(b)
        # replacing the matrix with the zero apdding.
        A=Ta; b=Tb
    """
    # Initialization of the row and column index
    i=0; j=0; indx=jndx; Lst = []
    while i<A.nrows() and j<A.ncols():
        while (A[i:,j]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        if A[i:,:].is_zero()==False:
            while A[i,j].is_zero():
                Ta=A[i:,:]
                Tb=b[i:,:]
                # Initializing the cyclic shift permutation matrix
                Id=identity_matrix(Ta.nrows())
                P =sum([Id[:,k]*Id[mod(k+1,Ta.nrows()),:] for k in range(Ta.nrows())])
                Ta=P*Ta; Tb=P*Tb
                A[i:,:]=Ta
                b[i:,:]=Tb 
            # Performing the row operations.
            b[i,0]=(b[i,0]*exp(I*2*pi*var('k'+str(indx))))^(1/A[i,j])
            indx = indx+1
            Lst.append(A[i,j])
            A[i,:]=(1/A[i,j])*A[i,:]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if A[r,:].is_zero():
                    r=r+1
                else:
                    b[r,0]=(b[i,0]*exp(I*2*pi*var('k'+str(indx))))^(-A[r,j])*b[r,0]
                    if (A[r,j]).is_zero()==False:
                        indx = indx+1
                        Lst.append(A[r,j])
                    A[r,:]=-A[r,j]*A[i,:]+A[r,:]
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return [A, b, indx, Lst]

def multiplicative_gauss_jordan_elimination(Cf,rs,jndx=0):
    """
    Outputs the reduced row echelon form of the input matrix and the right hand side.

    EXAMPLES:
 
    ::  
        sage: [RefA, c, indx] = multiplicative_gauss_jordan_elimination(Matrix(SR,HM(2,2,'a').listHM()), Matrix(SR,HM(2,1,'b').listHM()))
        sage: RefA
        [1 0]
        [0 1]
        sage: c
        [(b00*e^(2*I*pi*k0))^(1/a00)*((((b00*e^(2*I*pi*k0))^(1/a00)*e^(2*I*pi*k1))^(-a10)*b10*e^(2*I*pi*k2))^(-1/(a01*a10/a00 - a11))*e^(2*I*pi*k3))^(-a01/a00)]
        [                                                       (((b00*e^(2*I*pi*k0))^(1/a00)*e^(2*I*pi*k1))^(-a10)*b10*e^(2*I*pi*k2))^(-1/(a01*a10/a00 - a11))]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    [A, b, indx, Lst] = multiplicative_gaussian_elimination(Cf,rs,jndx)
    # Initialization of the row and column index
    i = A.nrows()-1; j = 0
    while i > 0 or j > 0:
        if (A[i,:]).is_zero():
            # decrementing the row index and initializing the column index
            i =i-1; j=0
        else :
            while (A[i,j]).is_zero():
                # Incrementing the column index
                j=j+1
            # performing row operations
            for r in range(i-1, -1, -1):
                b[r,0]=(b[i,0]*exp(I*2*pi*var('k'+str(indx))))^(-A[r,j])*b[r,0]
                if (A[r,j]).is_zero()==False:
                    indx = indx+1
                    Lst.append(A[r,j])
                A[r,:] = -A[r,j]*A[i,:]+A[r,:]
            i = i - 1; j = 0
    return [A, b, indx, Lst]

def log_gaussian_elimination(Cf, rs, jndx=0):
    """
    Outputs the row echelon form of the input matrix and the right hand side.

    EXAMPLES:
 
    ::  
        sage: [RefA, c] = log_gaussian_elimination(Matrix(SR,HM(2,2,'a').listHM()), Matrix(SR,HM(2,1,'b').listHM()))
        sage: RefA
        [      1 a01/a00]
        [      0       1]
        sage: c
        [                                                          (2*I*pi*k0 + b00)/a00]
        [((2*I*pi*k1 + (2*I*pi*k0 + b00)/a00)*a10 - 2*I*pi*k2 - b10)/(a01*a10/a00 - a11)]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    A = copy(Cf); b = copy(rs)
    """
    # Zero padding the matrix if necessary.
    if A.nrows() < A.ncols():
        # Temporary matrix for A
        Ta = Matrix(SR, zero_matrix(A.ncols(), A.ncols()))
        Ta[:A.nrows(), :A.ncols()] = copy(A)
        # Temporary matrix for b
        Tb = Matrix(SR, zero_matrix(A.ncols(), b.ncols())) 
        Tb[:b.nrows(), :b.ncols()] = copy(b)
        # replacing the matrix with the zero apdding.
        A = Ta; b = Tb
    """
    # Initialization of the row and column index
    i = 0; j = 0; indx = jndx
    while i < A.nrows() and j < A.ncols():
        while (A[i:,j]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        if (A[i:,:].is_zero()) == False:
            while (A[i,j]).is_zero():
                Ta = A[i:,:]
                Tb = b[i:,:]
                # Initializing the cyclic shift permutation matrix
                Id = identity_matrix(Ta.nrows())
                P  = sum([Id[:,k]*Id[mod(k+1,Ta.nrows()),:] for k in range(Ta.nrows())])
                Ta = P*Ta; Tb = P*Tb
                A[i:,:] = Ta
                b[i:,:] = Tb
            # Performing the row operations.
            b[i,0] = (1/A[i,j])*(b[i,0] + I*2*pi*var('k'+str(indx)))
            indx = indx + 1
            A[i,:] = (1/A[i,j])*A[i,:]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if A[r,:].is_zero():
                    r=r+1
                else:
                    b[r,0] = -A[r,j]*( b[i,0] + I*2*pi*var('k'+str(indx)) ) + b[r,0]
                    indx = indx+1
                    A[r,:] = -A[r,j]*A[i,:]+A[r,:]
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return [A, b, indx]

def log_gauss_jordan_elimination(Cf, rs, jndx=0):
    """
    Outputs the reduced row echelon form of the input matrix and the right hand side.

    EXAMPLES:
 
    ::  
        sage: [RefA, c] = log_gauss_jordan_elimination(Matrix(SR,HM(2,2,'a').listHM()), Matrix(SR,HM(2,1,'b').listHM()))
        sage: RefA
        [1 0]
        [0 1]
        sage: c
        [-(2*I*pi*k3 + ((2*I*pi*k1 + (2*I*pi*k0 + b00)/a00)*a10 - 2*I*pi*k2 - b10)/(a01*a10/a00 - a11))*a01/a00 + (2*I*pi*k0 + b00)/a00]
        [                                               ((2*I*pi*k1 + (2*I*pi*k0 + b00)/a00)*a10 - 2*I*pi*k2 - b10)/(a01*a10/a00 - a11)]
        
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    [A, b, indx] = log_gaussian_elimination(Cf,rs,jndx)
    # Initialization of the row and column index
    i = A.nrows()-1; j = 0
    while i > 0 or j > 0:
        if (A[i,:]).is_zero():
            # decrementing the row index and initializing the column index
            i = i - 1; j = 0
        else :
            while (A[i,j]).is_zero():
                # Incrementing the column index
                j = j + 1
            # performing row operations
            for r in range(i-1, -1, -1):
                b[r,0] = -A[r,j]*(b[i,0] + I*2*pi*var('k'+str(indx))) + b[r,0]
                indx = indx+1
                A[r,:] = -A[r,j]*A[i,:]+A[r,:]
            i = i - 1; j = 0
    return [A, b, indx]

def SecondOrderHadamardFactorization(A,B):
    """
    Outputs the matrix factorization induced by hadamard matrices.

    EXAMPLES:
 
    ::  
        sage: [U,V]=SecondOrderHadamardFactorization(HM(2,2,'a'),HM(2,2,'b'))
        sage: Prod(U,V).simplify()
        [[-a00*b00-a01*b10-a02*b20-a03*b30, a00*b01+a01*b11+a02*b21+a03*b31, a00*b02+a01*b12+a02*b22 +a03*b32, a00*b03+a01*b13+a02*b23+a03*b33], [a10*b00+a11*b10+a12*b20+a13*b30, -a10*b01-a11*b11-a12*b21-a13*b31, a10*b02+a11*b12+a12*b22+a13*b32, a10*b03+a11*b13+a12*b23+a13*b33], [a20*b00+a21*b10+a22*b20+a23*b30, a20*b01+a21*b11+a22*b21+a23*b31, -a20*b02-a21*b12-a22*b22-a23*b32, a20*b03+a21*b13+a22*b23+a23*b33], [a30*b00+a31*b10+a32*b20+a33*b30, a30*b01+a31*b11+a32*b21+a33*b31, a30*b02+a31*b12+a32*b22+a33*b32, -a30*b03-a31*b13-a32*b23-a33*b33]]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    if A.n(0)==B.n(1) and A.n(1)==B.n(0) and (log(Integer(A.n(0)),2)).is_integer():
        # Initializing the hadamard matrix
        H = SecondOrderHadamardBlockU(Integer(A.n(0)))
        L = [GeneralUncorrelatedHypermatrixTuple(2) for i in range(log(Integer(A.n(0)),2))]
        # Initializing the temporary hypermatrices.
        Tp0 = L[0][0]; Tp1 = L[0][1]
        for i in range(1,len(L)):
            Tp0 = Tp0.slicekroneckerproduct(L[i][0])
            Tp1 = Tp1.slicekroneckerproduct(L[i][1])
        Tp0 = sqrt(Tp0.n(0))*Tp0
        Tp1 = sqrt(Tp1.n(0))*Tp1
        # initializing the extened matrices
        M0 = HM(A.n(0), 2*A.n(0)-1, 'zero')
        for j in range(1,A.n(0)):
            for i in range(A.n(0)):
                M0[i,j-1] = H[i,j]
        for j in range(A.n(0), 2*A.n(0)):
            for i in range(A.n(0)):
                M0[i,j-1] = Tp0[i,j-A.n(0)]
        M1 = HM(2*A.n(0)-1,A.n(0), 'zero')
        for j in range(1,A.n(0)):
            for i in range(A.n(0)):
                M1[j-1,i] = -H[i,j]
        for j in range(A.n(0), 2*A.n(0)):
            for i in range(A.n(0)):
                M1[j-1,i] = Tp1[j-A.n(0),i]
        # Filling up the hypermatrices.
        U = HM(A.n(0), (2*A.n(0)-1)*A.n(1),'zero')
        for i in range(A.n(1)):
            for j in range(M0.n(1)):
                for k in range(A.n(0)):
                    U[k, A.n(1)*j+i] = A[k,i]*M0[k,j]
        V = HM((2*A.n(0)-1)*A.n(1), A.n(0),'zero')
        for i in range(B.n(0)):
            for j in range(M1.n(0)):
                for k in range(B.n(1)):
                    V[A.n(1)*j+i,k] = B[i,k]*M1[j,k]
        return [U,V]
    else:
        # return the error message if the input hypermatrix is cubic
        raise ValueError, "The input hypermpatrix are of inapropriate sizes"

def ThirdOrderHadamardFactorization(Ha, Hb, Hc):
    """
    Outputs the matrix factorization induced by Hadamard third order hypermatrices.
    The slices of third order matrices involved in outer-product computations must
    be square and must have sides whose length is a power of 2.   

    EXAMPLES:
 
    ::  
        sage: [U,V,W]=ThirdOrderHadamardFactorization(HM(4,2,4,'a'),HM(4,4,2,'b'), HM(2,4,4,'c'))
        sage: Prod(U, V, W).simplify()
[[[a000*b000*c000+a010*b001*c100,a001*b000*c001+a011*b001*c101,a002*b000*c002+a012*b001*c102,a003*b000*c003+a013*b001*c103],[a000*b010*c010+a010*b011*c110,a001*b010*c011+a011*b011*c111,a002*b010*c012+a012*b011*c112,a003*b010*c013+a013*b011*c113],[a000*b020*c020+a010*b021*c120,a001*b020*c021+a011*b021*c121,a002*b020*c022+a012*b021*c122,a003*b020*c023+a013*b021*c123],[a000*b030*c030+a010*b031*c130,a001*b030*c031+a011*b031*c131,a002*b030*c032+a012*b031*c132,a003*b030*c033+a013*b031*c133]],[[a100*b100*c000+a110*b101*c100,a101*b100*c001+a111*b101*c101,a102*b100*c002+a112*b101*c102,a103*b100*c003+a113*b101*c103],[a100*b110*c010+a110*b111*c110,a101*b110*c011+a111*b111*c111,a102*b110*c012+a112*b111*c112,a103*b110*c013+a113*b111*c113],[a100*b120*c020+a110*b121*c120,a101*b120*c021+a111*b121*c121,a102*b120*c022+a112*b121*c122,a103*b120*c023+a113*b121*c123],[a100*b130*c030+a110*b131*c130,a101*b130*c031+a111*b131*c131,a102*b130*c032+a112*b131*c132,a103*b130*c033+a113*b131*c133]],[[a200*b200*c000+a210*b201*c100,a201*b200*c001+a211*b201*c101,a202*b200*c002+a212*b201*c102,a203*b200*c003+a213*b201*c103],[a200*b210*c010+a210*b211*c110,a201*b210*c011+a211*b211*c111,a202*b210*c012+a212*b211*c112,a203*b210*c013+a213*b211*c113],[a200*b220*c020+a210*b221*c120,a201*b220*c021+a211*b221*c121,a202*b220*c022+a212*b221*c122,a203*b220*c023+a213*b221*c123],[a200*b230*c030+a210*b231*c130,a201*b230*c031+a211*b231*c131,a202*b230*c032+a212*b231*c132,a203*b230*c033+a213*b231*c133]],[[a300*b300*c000+a310*b301*c100,a301*b300*c001+a311*b301*c101,a302*b300*c002+a312*b301*c102,a303*b300*c003+a313*b301*c103],[a300*b310*c010+a310*b311*c110,a301*b310*c011+a311*b311*c111,a302*b310*c012+a312*b311*c112,a303*b310*c013+a313*b311*c113],[a300*b320*c020+a310*b321*c120,a301*b320*c021+a311*b321*c121,a302*b320*c022+a312*b321*c122,a303*b320*c023+a313*b321*c123],[a300*b330*c030+a310*b331*c130,a301*b330*c031+a311*b331*c131,a302*b330*c032+a312*b331*c132,a303*b330*c033+a313*b331*c133]]]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    if Ha.n(1)==Hb.n(2) and Hb.n(2)==Hc.n(0) and (log(Integer(Ha.n(0)),2)).is_integer():
        # Initializing the hadamard matrix
        H = ThirdOrderHadamardBlockU(Integer(Ha.n(0)))
        # Initializing the hypermatrices for the slice Kronecker product.
        L = [GeneralUncorrelatedHypermatrixTuple(3) for i in range(log(Integer(Ha.n(0)),2))]
        # Initializing the temporary hypermatrices.
        Tp0 = L[0][0]; Tp1 = L[0][1]; Tp2 = L[0][2]
        # Computing the slice kronecker product
        for i in range(1,len(L)):
            Tp0 = Tp0.slicekroneckerproduct(L[i][0])
            Tp1 = Tp1.slicekroneckerproduct(L[i][1])
            Tp2 = Tp2.slicekroneckerproduct(L[i][2])
        Tp0 = (Tp0.n(0))^(1/3)*Tp0
        Tp1 = (Tp1.n(0))^(1/3)*Tp1
        Tp2 = (Tp2.n(0))^(1/3)*Tp2
        print 'Prod(H, H.transpose(2), H.transpose())= ', Prod(H,H.transpose(2),H.transpose())
        # Initializing the extended third order hypermatrices
        M0 = HM(Ha.n(0), 2*Ha.n(0)-1, Ha.n(0),'zero')
        for j in range(Ha.n(0)-1):
            for k in range(Ha.n(0)):
                for i in range(Ha.n(0)):
                    M0[i,j,k] = -H[i,j,k]
        for j in range(Ha.n(0), 2*Ha.n(0)):
            for k in range(Ha.n(0)):
                for i in range(Ha.n(0)):
                    M0[i,j-1,k] = Tp0[i,j-Ha.n(0),k]
        M1 = HM(Ha.n(0), Ha.n(0), 2*Ha.n(0)-1, 'zero')
        for j in range(Ha.n(0)-1):
            for k in range(Ha.n(0)):
                for i in range(Ha.n(0)):
                    M1[i,k,j] = (H.transpose(2))[i,k,j]
        for j in range(Ha.n(0), 2*Ha.n(0)):
            for k in range(Ha.n(0)):
                for i in range(Ha.n(0)):
                    M1[i,k,j-1] = Tp1[i,k,j-Ha.n(0)]
        M2 = HM(2*Ha.n(0)-1, Ha.n(0), Ha.n(0),'zero')
        for j in range(Ha.n(0)-1):
            for k in range(Ha.n(0)):
                for i in range(Ha.n(0)):
                    M2[j,i,k] = (H.transpose())[j,i,k]
        for j in range(Ha.n(0), 2*Ha.n(0)):
            for k in range(Ha.n(0)):
                for i in range(Ha.n(0)):
                    M2[j-1,i,k] = Tp2[j-Ha.n(0),i,k]
        # Filling up the three thirdhypermatrices to be outputed.
        U = HM(Ha.n(0), (2*Ha.n(0)-1)*Ha.n(1), Ha.n(0), 'zero')
        for i in range(Ha.n(1)):
            for j in range(M0.n(1)):
                for k in range(Ha.n(0)):
                    for l in range(Ha.n(2)):
                        U[k,Ha.n(1)*j+i,l] = Ha[k,i,l]*M0[k,j,l]
        V = HM(Ha.n(0), Ha.n(0), (2*Ha.n(0)-1)*Ha.n(1), 'zero')
        for i in range(Hb.n(2)):
            for j in range(M1.n(2)):
                for k in range(Hb.n(0)):
                    for l in range(Hb.n(1)):
                        V[k,l,Ha.n(1)*j+i] = Hb[k,l,i]*M1[k,l,j]
        W = HM((2*Ha.n(0)-1)*Ha.n(1), Ha.n(0), Ha.n(0),'zero')
        for i in range(Hc.n(0)):
            for j in range(M2.n(0)):
                for k in range(Hc.n(1)):
                    for l in range(Hc.n(2)):
                        W[Ha.n(1)*j+i,k,l] = Hc[i,k,l]*M2[j,k,l]
        return [U, V, W]
    else:
        # return the error message if the input hypermatrix is cubic
        raise ValueError, "The input hypermpatrix are of inapropriate sizes"

