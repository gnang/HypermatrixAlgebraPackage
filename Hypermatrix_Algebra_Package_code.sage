#*************************************************************************#
#    Copyright (C) 2017, 18, 19, 20 Edinah K.Gnang <kgnang@gmail.com>,    #
#                          Ori Parzanchevski,                             #
#                          Yuval Filmus,                                  #
#                          Doron Zeilberger,                              #
#                                                                         #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                                                                         #
#    This code is distributed in the hope that it will be useful,         #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU     #
#    General Public License for more details.                             #
#                                                                         #
#  The full text of the GPL is available at:                              #
#                                                                         #
#                  http://www.gnu.org/licenses/                           #
#*************************************************************************#

# Definition of the hypermatrix class HM.
class HM:
    """
    Describes the Hypermatrix constructors
    
    INPUT:

    -  ``nrows`` - int, the number of rows
    -  ``ncols`` - int, the number of columns
    -  ``ndpts`` - int, the number of depths
    -  ``a`` - char, the symbolic entry prefix

    EXAMPLES::
    
        sage: HM(2,2,2,'a') # Initialization of a third order hypermatrix
        [[[a000, a001], [a010, a011]], [[a100, a101], [a110, a111]]]
        sage: od=2; sz=2; HM(od,sz,'kronecker') # Creating a Kronecker delta or identity matrix
        [[1, 0], [0, 1]]
        sage: od=2; HM(od,[0,1,2],'perm')
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        sage: od=2; sz=2; HM(od,sz,'a','sym') # Creating a symmetric matrix
        [[a00, a01], [a01, a11]]
        sage: od=2; sz=2; HM(od,sz,'a','skewsym') # Creating a symmetric matrix
        [[0, a01], [-a01, 0]]
        sage: sz=2; HM(sz,sz,sz,'a','shift') # Creating a hypermatrix where the index are shifted by one
        [[[a111, a112], [a121, a122]], [[a211, a212], [a221, a222]]]
        sage: HM(2,2,'one') # creating a 2 x 2 matrix whose entries are all equal to one
        [[1, 1], [1, 1]]
        sage: HM(2,2,'zero') # creating a 2 x 2 matrix whose entries are all equal to one
        [[0, 0], [0, 0]]
        sage: HM(HM(2,2,'a').matrix()) # Creating a scaling hypermatrix
        [[[a00, 0], [0, a01]], [[a10, 0], [0, a11]]]
        sage: HM(3,'x').list() # creatling a list of variables
        [x0, x1, x2]
        sage: sz=3; var_list('x',sz) # alternative approach to creatling a list of variables
        [x0, x1, x2]
        sage: od=2; sz=2; HM(od,var_list('a',sz),'diag') # diagonal matrix
        [[a0, 0], [0, a1]]
        sage: La=HM(2, 2, 2, 'a').list(); Lb=HM(2, 2, 2, 'b').list(); Lc=HM(2, 2, 2, 'c').list() # Initialization of the variables
        sage: F = FreeAlgebra(QQ, 24, La + Lb + Lc)
        sage: F.<a000,a100,a010,a110,a001,a101,a011,a111,b000,b100,b010,b110,b001,b101,b011,b111,c000,c100,c010,c110,c001,c101,c011,c111> = FreeAlgebra(QQ,24)
        sage: Ha=HM(2, 2, 2, [a000,a100,a010,a110,a001,a101,a011,a111])
        sage: Hb=HM(2, 2, 2, [b000,b100,b010,b110,b001,b101,b011,b111])
        sage: Hc=HM(2, 2, 2, [c000,c100,c010,c110,c001,c101,c011,c111])
        sage: Prod(Ha,Hb,Hc)[0,0,0]
        a000*b000*c000 + a010*b001*c100
        sage: Prod(Hb,Ha,Hc)[0,0,0]
        b000*a000*c000 + b010*a001*c100
        sage: od=2; sz=3; DltL=GeneralHypermatrixKroneckerDeltaL(od,sz)
        sage: DltL[0].printHM()
        [:, :]=
        [1 0 0]
        [0 0 0]
        [0 0 0]
        sage: DltL[1].printHM()
        [:, :]=
        [0 0 0]
        [0 1 0]
        [0 0 0]
        sage: DltL[2].printHM()
        [:, :]=
        [0 0 0]
        [0 0 0]
        [0 0 1]
        sage: Ha=HM(2,1,2,'a');Hb=HM(2,1,2,'b');Hc=HM(2,2,1,'c');Hd=HM(2,2,1,'d');Hf=HM(1,2,2,'f');Hg=HM(1,2,2,'g')
        sage: A=HM(1,2,1,[Ha,Hb]); B=HM(1,1,2,[Hc,Hd]); C=HM(2,1,1,[Hf,Hg])
        sage: BlockProd(A, B, C)[0,0,0].printHM() # Computing the block product
        [:, :, 0]=
        [a000*c000*f000 + b000*d000*g000 a000*c010*f010 + b000*d010*g010]
        [a100*c100*f000 + b100*d100*g000 a100*c110*f010 + b100*d110*g010]
        <BLANKLINE>
        [:, :, 1]=
        [a001*c000*f001 + b001*d000*g001 a001*c010*f011 + b001*d010*g011]
        [a101*c100*f001 + b101*d100*g001 a101*c110*f011 + b101*d110*g011]
        sage: HM(3,3,3,'a').slice([2], 'dpt').printHM() # Illustrating the slicing
        [:, :, 0]=
        [a002 a012 a022]
        [a102 a112 a122]
        [a202 a212 a222] 
        <BLANKLINE>
    """
    def __init__(self,*args):
        if len(args) == 1:
            inp = args[0]
            if type(inp)==type(Matrix(SR,2,1,[var('xxi'),var('yyj')])) or type(inp)==type(Matrix(RR,2,1,[1,2])) or type(inp)==type(Matrix(CC,2,1,[1,1])):
                self.hm=DiagonalHypermatrix(inp)
            elif type(inp) == list:
                self.hm = inp
            # Case associated with the initialization of a picture
            # the inlustration of how images are turned into a hypermatrix
            # A=HM("I.png")
            elif type(inp) == type("abcd"):
                # Importing the numerical library
                import pylab, numpy
                from scipy import misc
                # Defining the input Pictures
                X  = pylab.imread(inp)
                sz = max(len(X),len(X[0]),len(X[0][0]))
                args = (len(X), len(X[0]), len(X[0][0]))
                #T = apply(HypermatrixGenerateAllZero, args)
                AtmpL=args
                T = HypermatrixGenerateAllZero(*AtmpL)
                # Filling up the image Hypermatrix
                for i in range(len(X)):
                    for j in range(len(X[0])):
                        for k in range(len(X[0][0])):
                            T[i][j][k] = X[i,j,k]
                self.hm = T
            else:
                raise ValueError("Expected list as input when generating diagonal hypermatrix")
            return
        # Obtaining the last argument
        s = args[-1]
        # Initializing the dimension parameters
        dims = args[:-1]
        if s == 'one':
            AtmpL=dims
            self.hm = HypermatrixGenerateAllOne(*AtmpL)
        elif s == 'zero':
            AtmpL=dims
            self.hm = HypermatrixGenerateAllZero(*AtmpL)
        elif s == 'shift':
            AtmpL=dims
            self.hm = HypermatrixGenerateII(*AtmpL)
        elif s == 'perm':
            if dims[0]==2:
                # Initialization of the size parameter
                sz=len(dims[1])
                Id=HM(2,sz,'kronecker')
                # The permutation acts on columns
                self.hm=sum(HM(sz,1,[Id[i,dims[1][k]] for i in range(sz)])*HM(1,sz,[Id[k,j] for j in range(sz)]) for k in range(sz)).transpose().listHM()
            else:
                raise ValueError("not supported for order %d hypermatrices" % len(dims))
        elif s == 'kronecker':
            AtmpL=args[:-1]
            self.hm=GeneralHypermatrixKroneckerDelta(*AtmpL).listHM()
        elif s == 'diag':
            AtmpL=args[:-1]
            self.hm=GeneralHypermatrixMainDiag(*AtmpL).listHM()
        elif s == 'sym':
            if (len(dims) == 3) and (dims[0]==2):
                self.hm=SymMatrixGenerate(dims[1],dims[2])
            elif (len(dims) == 3) and (dims[0]==3):
                self.hm=SymHypermatrixGenerate(dims[1],dims[2])
            else:
                raise ValueError("SymHypermatrixGenerate not supported for order %d hypermatrices" % dims[0])
        elif s == 'skewsym':
            if (len(dims) == 3) and (dims[0]==2):
                self.hm=SkewSymMatrixGenerate(dims[1],dims[2])
            elif (len(dims) == 3) and (dims[0]==3):
                self.hm=SkewSymHypermatrixGenerate(dims[1],dims[2])
            else:
                raise ValueError("SkewSymHypermatrixGenerate not supported for order %d hypermatrices" % dims[0])
        elif type(s) == list:
            self.hm=(List2Hypermatrix(*args)).listHM()
        else:
            self.hm=HypermatrixGenerate(*args)
    def __repr__(self):
        return self.hm.__repr__()
    def __pow__(self, other):
        """
        This method returns the exponentiation
        operation. I uses the convention I have
        introduced for exponentiation which work
        for coformable block partition and relates
        multiplicative constraints. The functions
        should not be confused with exp(B*ln(A))


        EXAMPLES:

        ::

            sage: A=HM(2,2,'a'); B=HM(2,2,'b')
            sage: (A^B).printHM()
            [:, :]=
            [a00^b00*a10^b01 a01^b00*a11^b01]
            [a00^b10*a10^b11 a01^b10*a11^b11]


        AUTHORS:

        - Edinah K. Gnang
        """
        if self.order()==2 and other.order()==2 and self.n(0)==other.n(1):
            #return mprod(other, self)
            return GProdIII([other,self],prod,BaseExp)
        elif self.order()==2 and self.is_hypercolumn() and prod(self.dimensions())>1:
            return vec_exp(self, other)
        elif self.order()==2 and other==-1 and self.is_cubical():
            return self.inverse()
        elif self.order()==2 and other==0 and self.is_cubical():
            return HM(self.order(), self.n(0), 'kronecker')
        elif self.order()==2 and type(other)==type(5) and self.is_cubical():
            if other > 0:
                return self*(self^(other-1))
            elif other < 0:
                return self.inverse()^(-other)
        elif self.dimensions()==other.dimensions():
            return GeneralHypermatrixHadamardExponent(self,other)
        else:
            raise ValueError("Operation is not supported !")
    def __add__(self, other):
        return GeneralHypermatrixAdd(self,other)
    def __truediv__(other,self):
        """
        This method returns the product self
        with the inverse of the input other.
        It only works for second order hypermatrices.


        EXAMPLES:

        ::

            sage: A=HM(2,2,'a'); B=HM(2,var_list('b',2),'diag')
            sage: (A/B).printHM()
            [:, :]=
            [a00/b0 a01/b1]
            [a10/b0 a11/b1]


        AUTHORS:

        - Edinah K. Gnang
        """
        if self.order()==2:
            return other*self.inverse()
    def __radd__(self, other):
        return GeneralHypermatrixAdd(self,other)
    def __neg__(self):
        return GeneralHypermatrixScale(self,-1)
    def __sub__(self, other):
        return GeneralHypermatrixAdd(self, GeneralHypermatrixScale(other,-1))
    def __mul__(self, other):
        """
        This method returns the BM product, it
        is a short cuc of sort which avoid having
        to write Prod all the time.


        EXAMPLES:

        ::

            sage: A=HM(2,2,'a'); B=HM(2,2,'b') # The matrix setting
            sage: (A*B).printHM()
            [:, :]=
            [a00*b00 + a01*b10 a00*b01 + a01*b11]
            [a10*b00 + a11*b10 a10*b01 + a11*b11]
            sage: A=HM(2,2,2,'a'); B=HM(2,2,2,'a'); C=HM(2,2,2,'c')
            sage: D=A*(B,C); D.printHM()
            [:, :, 0]=
            [   a000^2*c000 + a001*a010*c100 a000*a010*c010 + a010*a011*c110]
            [   a100^2*c000 + a101*a110*c100 a100*a110*c010 + a110*a111*c110]
            <BLANKLINE>
            [:, :, 1]=
            [a000*a001*c001 + a001*a011*c101    a001*a010*c011 + a011^2*c111]
            [a100*a101*c001 + a101*a111*c101    a101*a110*c011 + a111^2*c111]
            <BLANKLINE>


        AUTHORS:

        - Edinah K. Gnang
        """
        if other.__class__.__name__=='HM':
            return Prod(self,other)
        elif other.__class__.__name__=='tuple':
            l=other
            return Prod(self,*l)
        else: 
            return GeneralHypermatrixScaleRight(self,other)
    def __rmul__(self, a):
        return GeneralHypermatrixScale(self,a)
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
        """
        This method uses the call functionality
        to return the general BM product inspired
        the notation the self is the background
        hypermatrix while inpts refers to the 
        conformable hypermatrices to be multiplied.


        EXAMPLES:

        ::

            sage: A=HM(2,2,2,'a'); B=HM(2,2,2,'b'); C=HM(2,2,2,'c'); D=HM(2,2,2,'d')
            sage: D(A,B,C)-ProdB(A,B,C,D)
            [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]


        AUTHORS:

        - Edinah K. Gnang
        """
        AtmpL=[inp for inp in inpts]+[self]
        return ProdB(*AtmpL)
    def __eq__(self, other):
        return (isinstance(other, self.__class__) and self.list() == other.list())
    def __ne__(self, other):
        return not self.__eq__(other)

    def elementwise_product(self, B):
        """
        Returns the elementwise product of two 
        hypermatrices.
        This routine assumes that ``self`` and ``B``
        are both hypermatrices, with identical
        sizes. It is "unsafe" in the sense that these
        conditions are not checked and no sensible 
        errors are raised.

        EXAMPLES:

        ::

            sage: A=HM(2,2,2,'a'); B=HM(2,2,2,'b')
            sage: A.elementwise_product(B)
            [[[a000*b000, a001*b001], [a010*b010, a011*b011]], [[a100*b100, a101*b101], [a110*b110, a111*b111]]]


        AUTHORS:

        - Edinah K. Gnang
        """
        return GeneralHypermatrixHadamardProduct(self, B)

    def hadamard_product(self, B):
        """
        Returns the elementwise product of two 
        hypermatrices.
        This routine assumes that ``self`` and ``B``
        are both hypermatrices, with identical
        sizes. It is "unsafe" in the sense that these
        conditions are not checked and no sensible 
        errors are raised.

        EXAMPLES:

        ::

            sage: A=HM(2,2,2,'a'); B=HM(2,2,2,'b')
            sage: A.hadamard_product(B)
            [[[a000*b000, a001*b001], [a010*b010, a011*b011]], [[a100*b100, a101*b101], [a110*b110, a111*b111]]]


        AUTHORS:

        - Edinah K. Gnang
        """
        return GeneralHypermatrixHadamardProduct(self, B)

    def elementwise_exponent(self,s):
        """
        Returns the elementwise exponent of the entries
        ``self`` and ``B`` are both hypermatrices, with
        identical sizes. It is "unsafe" in the sense that
        these conditions are not checked and no sensible 
        errors are raised.

        EXAMPLES:

        ::

            sage: A=HM(2,2,2,'a')
            sage: A.elementwise_exponent(3)
            [[[a000^3, a001^3], [a010^3, a011^3]], [[a100^3, a101^3], [a110^3, a111^3]]]


        AUTHORS:

        - Edinah K. Gnang
        """
        return GeneralHypermatrixExponent(self, s)

    def hadamard_exponent(self, B):
        """
        Returns the elementwise of two hypermatrices.
        This routine assumes that ``self`` and ``B``
        are both hypermatrices, with identical
        sizes. It is "unsafe" in the sense that these
        conditions are not checked and no sensible 
        errors are raised.

        EXAMPLES:

        ::

            sage: A=HM(2,2,2,'a'); B=HM(2,2,2,'b')
            sage: A.hadamard_exponent(B)
            [[[a000^b000, a001^b001], [a010^b010, a011^b011]], [[a100^b100, a101^b101], [a110^b110, a111^b111]]]


        AUTHORS:

        - Edinah K. Gnang
        """
        return GeneralHypermatrixHadamardExponent(self, B)

    def elementwise_base_exponent(self, s):
        """
        Returns a hypermatrix whose entries are all
        raised by s.
        This routine assumes that we are not dividing
        by zero. It is "unsafe" in the sense that these
        conditions are not checked and no sensible 
        errors are raised.

        EXAMPLES:

        ::

            sage: A=HM(2,2,2,'a')
            sage: A.elementwise_base_exponent(3)
            [[[3^a000, 3^a001], [3^a010, 3^a011]], [[3^a100, 3^a101], [3^a110, 3^a111]]]


        AUTHORS:

        - Edinah K. Gnang
        """
        return GeneralHypermatrixBaseExponent(self, s)

    def elementwise_base_exponentN(self, s, dgts=50):
        """
        Returns a hypermatrix whose entries are all
        raised by s.
        This routine assumes that we are not dividing
        by zero. It is "unsafe" in the sense that these
        conditions are not checked and no sensible 
        errors are raised.

        EXAMPLES:

        ::

            sage: A=HM(2,2,[1,2,3,4]); A.elementwise_base_exponentN(3)
            [[3.0000000000000, 27.000000000000], [9.0000000000000, 81.000000000000]]
            


        AUTHORS:

        - Edinah K. Gnang
        """
        return GeneralHypermatrixBaseExponentN(self, s, dgts)

    def elementwise_base_logarithm(self, s):
        """
        Outputs a list of lists associated with the general
        whose entries are logarithms to the base s of the 
        original hypermatrix.
        The code only handles the Hypermatrix HM class objects.

        EXAMPLES:

        ::

            sage: Ha=HM(2,2,2,'a')
            sage: Ha.elementwise_base_logarithm(3)
            [[[log(a000)/log(3), log(a001)/log(3)], [log(a010)/log(3), log(a011)/log(3)]], [[log(a100)/log(3), log(a101)/log(3)], [log(a110)/log(3), log(a111)/log(3)]]]


        AUTHORS:
        - Edinah K. Gnang
        """
        return GeneralHypermatrixLogarithm(self, s)

    def elementwise_base_logarithmN(self, s, dgts=50):
        """
        Outputs a list of lists associated with the general
        whose entries are numerical logarithms to the base s of the 
        original hypermatrix.
        The code only handles the Hypermatrix HM class objects.

        EXAMPLES:

        ::

            sage: Ha=HM(2,2,[1,2,3,4]); Ha.elementwise_base_logarithmN(3)
            [[0.00000000000000, 1.0000000000000], [0.63092975357146, 1.2618595071429]]


        AUTHORS:
        - Edinah K. Gnang
        """
        return GeneralHypermatrixLogarithmN(self, s, dgts)

    def apply_map(self, phi):
        """
        Apply the given map phi (an arbitrary Python function or callable
        object) to this hypermatrix.
        
         
        EXAMPLES::

        ::

            sage: A = HM(2,2,'a')
            sage: phi = lambda x: sin(x)
            sage: A.apply_map(phi).printHM()
            [:, :]=
            [sin(a00) sin(a01)]
            [sin(a10) sin(a11)]

        """
        return GeneralHypermatrixApplyMap(self, phi)

    def map(self, phi):
        """
        Apply the given map phi (an arbitrary Python function or callable
        object) to this hypermatrix.
        
         
        EXAMPLES::

        ::

            sage: A = HM(2,2,'a')
            sage: phi = lambda x: sin(x)
            sage: A.map(phi).printHM()
            [:, :]=
            [sin(a00) sin(a01)]
            [sin(a10) sin(a11)]

        """
        return GeneralHypermatrixApplyMap(self, phi)

    def tensor_product(self, V):
        """
        Computes the  Kronecker product of arbitrary hypermatrices A, B of the same order.

        EXAMPLES:

        ::

            sage: A=HM(2,2,'a'); B=HM(2,2,'b')
            sage: A.tensor_product(B)
            [[a00*b00, a00*b01, a01*b00, a01*b01], [a00*b10, a00*b11, a01*b10, a01*b11], [a10*b00, a10*b01, a11*b00, a11*b01], [a10*b10, a10*b11, a11*b10, a11*b11]]
        

        AUTHORS:
        - Edinah K. Gnang
        """
        if  self.order()==2:
            return SecondOrderSliceKroneckerProduct(self, V)
        elif self.order()==3:
            return ThirdOrderSliceKroneckerProduct(self, V)
        elif self.order()==4:
            return FourthOrderSliceKroneckerProduct(self, V)
        elif self.order()==5:
            return FifthOrderSliceKroneckerProduct(self, V)
        elif self.order()==6:
            return SixthOrderSliceKroneckerProduct(self, V)
        elif self.order()==7:
            return SeventhOrderSliceKroneckerProduct(self, V)
        elif self.order()>7:
            return GeneralHypermatrixKroneckerProduct(self, V)
        else :
            raise ValueError("not supported for order %d hypermatrices" % self.order())

    def tensor_power(self, m):
        """
        Computes the Kronecker product with itself of arbitrary hypermatrices
        m times.

        EXAMPLES:

        ::

            sage: A=HM(2,2,'a')
            sage: A.tensor_power(2)
            [[a00^2, a00*a01, a00*a01, a01^2], [a00*a10, a00*a11, a01*a10, a01*a11], [a00*a10, a01*a10, a00*a11, a01*a11], [a10^2, a10*a11, a10*a11, a11^2]]
        

        AUTHORS:
        - Edinah K. Gnang
        """
        if m == 1:
            return self.copy()
        elif m in ZZ and m > 1 :
            Tmp=self.copy()
            for i in range(m-1):
                Tmp=Tmp.tensor_product(self)  
            return Tmp
        else:
            raise ValueError("Expected a positive integer as argument")

    def block_sum(self, V):
        """
        Returns the block sum of two hypermatrices.
        This routine assumes that ``self`` and ``B``
        are both hypermatrices, with identical
        sizes. It is "unsafe" in the sense that these
        conditions are not checked and no sensible 
        errors are raised.

        EXAMPLES:

        ::

            sage: HM(2,2,'a').block_sum(HM(2,2,'b'))
            [[a00, a01, 0, 0], [a10, a11, 0, 0], [0, 0, b00, b01], [0, 0, b10, b11]]

        AUTHORS:
        - Edinah K. Gnang
        """
        return GeneralHypermatrixKroneckerSum(self, V)

    def expand(self):
        """
        Outputs a list of lists associated with the general
        hypermatrix with expressions in the entries in their
        expanded form.
        The code only handles the Hypermatrix HM class objects.

        EXAMPLES:

        ::

            sage: Ha=HM(2,2,2,[(var('x')+var('y'))^(i+j+k) for i in range(2) for j in range(2) for k in range(2)])
            sage: Ha.expand()
            [[[1, x + y], [x + y, x^2 + 2*x*y + y^2]], [[x + y, x^2 + 2*x*y + y^2], [x^2 + 2*x*y + y^2, x^3 + 3*x^2*y + 3*x*y^2 + y^3]]]
        

        AUTHORS:
        - Edinah K. Gnang
        """
        return GeneralHypermatrixExpand(self)

    def factor(self):
        """
        Outputs a list of lists associated with the general
        hypermatrix with expressions in the entries in their
        factored form.
        The code only handles the Hypermatrix HM class objects.

        EXAMPLES:

        ::

            sage: sz=2; vx=HM(sz,'x').list(); vy=HM(sz,'y').list(); vz=HM(sz,'z').list()
            sage: X=HM(sz, 1, sz, [vx[u] for v in range(sz) for t in range(1) for u in range(sz)])
            sage: Y=HM(sz, 1, sz, [vy[u] for v in range(sz) for t in range(1) for u in range(sz)])
            sage: Z=HM(sz, 1, sz, [vz[u] for v in range(sz) for t in range(1) for u in range(sz)])
            sage: A=Prod(X, Y.transpose(2), Z.transpose())
            sage: Prod(A,A,A).factor()
            [[[(x0*y0*z0 + x1*y1*z1)*x0^2*y0^2*z0^2, (x0*y0*z0 + x1*y1*z1)*x0^2*y0^2*z1^2], [(x0*y0*z0 + x1*y1*z1)*x0^2*y1^2*z0^2, (x0*y0*z0 + x1*y1*z1)*x0^2*y1^2*z1^2]], [[(x0*y0*z0 + x1*y1*z1)*x1^2*y0^2*z0^2, (x0*y0*z0 + x1*y1*z1)*x1^2*y0^2*z1^2], [(x0*y0*z0 + x1*y1*z1)*x1^2*y1^2*z0^2, (x0*y0*z0 + x1*y1*z1)*x1^2*y1^2*z1^2]]]

        AUTHORS:
        - Edinah K. Gnang
        """
        return GeneralHypermatrixFactor(self)

    def simplify(self):
        """
        Performs the symbolic simplification of the expressions
        associated with the hypermatrix entries. 

        EXAMPLES:

        ::

            sage: x,y=var('x,y'); ((x+y)^2*HM(2,2,2,'one')).simplify()
            [[[(x + y)^2, (x + y)^2], [(x + y)^2, (x + y)^2]], [[(x + y)^2, (x + y)^2], [(x + y)^2, (x + y)^2]]]
 

        AUTHORS:

        - Edinah K. Gnang
        """
        return GeneralHypermatrixSimplify(self)

    def simplify_full(self):
        """
        Performs the symbolic simplification of the expressions
        associated with the hypermatrix entries. 

        EXAMPLES:

        ::

            sage: x,y=var('x,y'); ((x+y)^2*HM(2,2,2,'one')).simplify_full()
            [[[x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2], [x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2]], [[x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2], [x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2]]]
 

        AUTHORS:

        - Edinah K. Gnang
        """
        return GeneralHypermatrixSimplifyFull(self)

    def canonicalize_radical(self):
        """
        Performs the symbolic simplification of the expressions
        associated with the hypermatrix entries. 

        EXAMPLES:

        ::

            sage: x,y = var('x,y') 
            sage: ((x^2+2*x*y+y^2)*HM(2,2,2,'one')).canonicalize_radical() 
            [[[x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2], [x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2]], [[x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2], [x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2]]]

        AUTHORS:

        - Edinah K. Gnang
        """
        return GeneralHypermatrixCanonicalizeRadical(self)

    def numerical(self, dgts=15):
        """
        Performs the symbolic simplification of the expressions
        associated with the hypermatrix entries. 

        EXAMPLES:

        ::

            sage: Ha=HM(3,3,[exp(I*2*pi*u*v/3) for v in range(3) for u in range(3)]) 
            sage: Ha.numerical()
            [[1.00000000000000, 1.00000000000000, 1.00000000000000], [1.00000000000000, -0.500000000000000 + 0.866025403784439*I, -0.500000000000000 - 0.866025403784439*I], [1.00000000000000, -0.500000000000000 - 0.866025403784439*I, -0.500000000000000 + 0.866025403784439*I]] 
         

        AUTHORS:

        - Edinah K. Gnang
        """
        return GeneralHypermatrixNumerical(self, dgts)

    def subs(self, *args, **kwds):
        """
        Substitute values to the variables in the hypermatrix.
        All the arguments are transmitted unchanged to the method ``subs`` of
        the coefficients.

        EXAMPLES:

        ::

            sage: var('a,b,d,e')
            (a, b, d, e)
            sage: m=HM([[a,b], [d,e]])
            sage: m.subs(a=1).printHM()
            [:, :]=
            [1 b]
            [d e]
            sage: m.subs(a=b, b=d).printHM()
            [:, :]=
            [b d]
            [d e]
            sage: m.subs({a: 3, b:2, d:1, e:-1}).printHM()
            [:, :]=
            [ 3  2]
            [ 1 -1]
            sage: m.subs([a==3, b==2, d==1, e==-1]).printHM()
            [:, :]=
            [ 3  2]
            [ 1 -1]
        """
        return GeneralHypermatrixSubstituteII(self, *args, **kwds)

    def subsn(self,Dct):
        """
        Procedure for computes the substitution in the Hypermatrix entries
        the inputs are the corresponding Hypermatric and a dictionary 
        datastructure.

        EXAMPLES:

        ::

            sage: HM(3,2,'a','sym').subsn(dict([(var('a011'),1),(var('a001'),2),(var('a000'),3),(var('a111'),4)]))
            [[[3.00000000000000, 2.00000000000000], [2.00000000000000, 1.00000000000000]], [[2.00000000000000, 1.00000000000000], [1.00000000000000, 4.00000000000000]]]


        AUTHORS:

        - Edinah K. Gnang
        """
        return GeneralHypermatrixSubstituteN(self, Dct)

    def substituteHMinto(self, poly, vrbl):
        """
        Procedure for computing substitution of the Hypermatrix into
        the input polynomial

        EXAMPLES:

        ::

            sage: x = var('x'); A=HM(2,2,'a')
            sage: p=x^2 - A.trace()*x + A.det()
            sage: A.substituteHMinto(p, x).expand()
            [[0, 0], [0, 0]]

        AUTHORS:

        - Edinah K. Gnang
        """
        return substituteHM(poly, vrbl, self)

    def substituteInHM(self, vrbl, M):
        """
        Procedure for computing substitution into the Hypermatrix 
        polynomial entries of a given input matrix for the input
        variable.


        EXAMPLES:

        ::

            sage: x = var('x'); Ha=HM(2,1,[x+1,x^2+1])
            sage: Y=var_list('y',2); rM=Ha.substituteInHM(x,HM(2,2,[Y[0],0,0,Y[1]])); rM
            [[[[y0 + 1, 0], [0, y1 + 1]]], [[[y0^2 + 1, 0], [0, y1^2 + 1]]]]


        AUTHORS:

        - Edinah K. Gnang
        """
        return GeneralHypermatrixSubstituteInMatrix(self,vrbl,M)

    def transpose(self, i=1):
        """
        Outputs a list of lists associated with the general
        transpose as defined by the cyclic permutation of indices.
        The code only handles the Hypermatrix HM class objects.

        EXAMPLES:

        ::

            sage: Ha=HM(2,2,2,'a')
            sage: Ha.transpose()
            [[[a000, a100], [a001, a101]], [[a010, a110], [a011, a111]]]

        AUTHORS:
        - Edinah K. Gnang
        """
        t = Integer(mod(i, self.order()))
        A = self 
        for i in range(t):
            A = GeneralHypermatrixCyclicPermute(A)
        return A

    def conjugate(self):
        """
        Outputs the conjugate of the hypermatrix. 

        EXAMPLES:

        ::

            sage: Ha=HM(2,1,[1+3*I, I])
            sage: Ha.conjugate()
            [[-3*I + 1], [-I]]


        AUTHORS:

        - Edinah K. Gnang
        """ 
        A = self 
        return GeneralHypermatrixConjugate(A)

    def conjugate_transpose(self, i=1):
        """
        Outputs the conjugate transpose of the hypermatrix. 

        EXAMPLES:

        ::

            sage: Ha=HM(2,1,[1+3*I, I])
            sage: Ha.conjugate_transpose()
            [[-3*I + 1, -I]]

        AUTHORS:

        - Edinah K. Gnang
        """ 
        t = Integer(mod(i, self.order()))
        A = self 
        for i in range(1,t+1):
            A = GeneralHypermatrixCyclicPermute(A)
        if Integer(mod(i,2))==0:
            return A
        else:
            return GeneralHypermatrixConjugate(A)

    def tumble(self, i=1):
        """
        Outputs a list of lists associated with the tumble
        transpose. 
        The code only handles second order HM class objects.

        EXAMPLES:

        ::

            sage: sz=5; (HM(sz,sz,'a') - HM(sz,sz,'a').tumble()).printHM()
            [:, :]=
            [ a00 - a40  a01 - a30  a02 - a20  a03 - a10 -a00 + a04]
            [ a10 - a41  a11 - a31  a12 - a21 -a11 + a13 -a01 + a14]
            [ a20 - a42  a21 - a32          0 -a12 + a23 -a02 + a24]
            [ a30 - a43  a31 - a33 -a23 + a32 -a13 + a33 -a03 + a34]
            [ a40 - a44 -a34 + a41 -a24 + a42 -a14 + a43 -a04 + a44]


        AUTHORS:
        - Edinah K. Gnang
        """
        if self.order()==2 and self.is_cubical():
            t = Integer(mod(i,4))
            tA = self.copy() 
            for itr in range(t):
                B=HM(tA.n(0),tA.n(1),'zero')
                for i in rg(tA.n(0)):
                    for j in rg(tA.n(1)):
                        B[i,j]=tA[tA.n(0)-1-j,i]
                tA=B.copy()
            return tA
        else:
            raise ValueError("Expected a cubical second order hypermatrix")

    def index_rotation(self, T):
        """
        Outputs the index rotation by angle T of self.
        The method is only implemented to second and third
        order hypermatrices. The rotation is performed
        clockwise by multiples of 2*pi/4 for second order hypermatrices as follows
        [i,j] -> [(i-floor(sz/2))*cos(T)+(-j+floor(sz/2))*sin(T)+floor(sz/2), (i-floor(sz/2))*sin(T)-(-j+floor(sz/2))*cos(T)+floor(sz/2)] if sz is odd
        [i,j] -> [(i-(sz-1)/2)*cos(T)+(-j+(sz-1)/2)*sin(T)+(sz-1)/2, (i-(sz-1)/2)*sin(T)-(-j+(sz-1)/2)*cos(T)+(sz-1)/2] if sz is even


        EXAMPLES:

        ::

            sage: sz=5; (HM(sz,sz,'a') - HM(sz,sz,'a').index_rotation(pi/2)).printHM()
            [:, :]=
            [ a00 - a40  a01 - a30  a02 - a20  a03 - a10 -a00 + a04]
            [ a10 - a41  a11 - a31  a12 - a21 -a11 + a13 -a01 + a14]
            [ a20 - a42  a21 - a32          0 -a12 + a23 -a02 + a24]
            [ a30 - a43  a31 - a33 -a23 + a32 -a13 + a33 -a03 + a34]
            [ a40 - a44 -a34 + a41 -a24 + a42 -a14 + a43 -a04 + a44]
            sage: sz=2; A=HM(sz, sz, sz, 'a')
            sage: A.index_rotation([2*pi/4, 0, 0]).printHM()
            [:, :, 0]=
            [a100 a110]
            [a101 a111]
            <BLANKLINE>
            [:, :, 1]=
            [a000 a010]
            [a001 a011]
            <BLANKLINE>


        AUTHORS:
        - Edinah K. Gnang and Fan Tian
        """
        if self.order()==2:
            return SecondOrderIndexRotation(self, T)
        elif self.order()==3:
            return ThirdOrderIndexRotation(self, T)
        else :
            raise ValueError("not supported for order %d hypermatrices" % self.order())

    def select_index_rotation(self, T, EntryList):
        """
        Outputs a list of lists associated with the tumble
        transpose. 
        The code only handles second order HM class objects.

        EXAMPLES:

        ::

            sage: sz=5; (HM(sz,sz,'a') - HM(sz,sz,'a').select_index_rotation(pi/2, rg(sz))).printHM()
            [:, :]=
            [ a00 - a40  a01 - a30  a02 - a20  a03 - a10 -a00 + a04]
            [ a10 - a41  a11 - a31  a12 - a21 -a11 + a13 -a01 + a14]
            [ a20 - a42  a21 - a32          0 -a12 + a23 -a02 + a24]
            [ a30 - a43  a31 - a33 -a23 + a32 -a13 + a33 -a03 + a34]
            [ a40 - a44 -a34 + a41 -a24 + a42 -a14 + a43 -a04 + a44]
            sage: sz=5; (HM(sz,sz,'a') - HM(sz,sz,'a').select_index_rotation(pi/2, [1,3])).printHM()
            [:, :]=
            [         0          0          0          0          0]
            [         0  a11 - a31          0 -a11 + a13          0]
            [         0          0          0          0          0]
            [         0  a31 - a33          0 -a13 + a33          0]
            [         0          0          0          0          0]
            sage: sz=2; A=HM(sz, sz, sz, 'a')
            sage: A.select_index_rotation([2*pi/4, 0, 0], rg(sz)).printHM()
            [:, :, 0]=
            [a100 a110]
            [a101 a111]
            <BLANKLINE>
            [:, :, 1]=
            [a000 a010]
            [a001 a011]
            <BLANKLINE>


        AUTHORS:

        - Edinah K. Gnang and Fan Tian 
        - To Do: Implement the arbitrary order version
        """
        if self.order()==2:
            return SelectSecondOrderIndexRotation(self, T, EntryList)
        elif self.order()==3:
            return SelectThirdOrderIndexRotation(self, T, EntryList)
        else :
            raise ValueError("not supported for order %d hypermatrices" % self.order())

    def slice_index_rotation(self, T, Lslice):
        """
        Outputs a hypermatrix which perform index
        rotations only on the designated slice.
        This function was inspired by looking at TCO
        play with his 3 x 3 x 3 rubik's cube.
        In this context we represent a 3 x 3 x 3 rubik's cube
        as a 3 x 3 x 3 hypermatrix whose entries are 2 x 2 x 2
        rank one symbolic hypermatrices


        EXAMPLES:

        ::

            sage: sz=2; HM(sz,sz,sz,'a').slice_index_rotation([0,0,2*pi/4], [1]).printHM()
            [:, :, 0]=
            [a000 a010]
            [a100 a110]
            <BLANKLINE>
            [:, :, 1]=
            [a101 a001]
            [a111 a011]
            <BLANKLINE>
            sage: AlphaB = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
            sage: Ha = HM(2,2,2,[HM(2,2,2,AlphaB[i]) for i in rg(2^3)])
            sage: Hb = Ha.slice_index_rotation([0,0,2*pi/4], [0])
            sage: Hb[0,1,0].printHM()
            [:, :, 0]=
            [a100 a000]
            [a110 a010]
            <BLANKLINE>
            [:, :, 1]=
            [a101 a001]
            [a111 a011]
            <BLANKLINE>


        AUTHORS:

        - Edinah K. Gnang
        - To Do: Implement the arbitrary order version
        """
        if self.order()==3:
            if type(self[0,0,0]) == type(HM(2,2,2,'a')):
                # Performing the outside transformation
                Tmp=ThirdOrderSliceIndexRotation(self, T, Lslice); sz=Tmp.n(0)
                # Performing the inside transformation
                return HM(sz, sz, sz, [Tmp[i,j,k].index_rotation(T) for k in rg(sz) for j in rg(sz) for i in rg(sz)])
            else: 
                return ThirdOrderSliceIndexRotation(self, T, Lslice)
        else :
            raise ValueError("not supported for order %d hypermatrices" % self.order())

    def nrows(self):
        return len(self.hm)

    def ncols(self):
        return len(self.hm[0])

    def ndpts(self):
        return len(self.hm[0][0])

    def inverse(self):
        """
        Returns the matrix inverse. 

        EXAMPLES:

        ::

            sage: Ha=HM(2,2,'a')
            sage: Ha.inverse().printHM()
            [:, :]=
            [1/a00 - a01*a10/(a00^2*(a01*a10/a00 - a11))               a01/(a00*(a01*a10/a00 - a11))]
            [              a10/(a00*(a01*a10/a00 - a11))                      -1/(a01*a10/a00 - a11)]

        AUTHORS:

        - Edinah K. Gnang
        """
        if self.order()==2 and self.is_cubical():
            return HM(self.n(0), self.n(1), Matrix(SR, self.listHM()).inverse().transpose().list()) 
        else:
            raise ValueError("not supported for order %d hypermatrices" %self.order())

    def printHM(self):
        """
        Conveniently displays matrices and third order hypermatrices
        one depth slice at a time.

        EXAMPLES:

        ::

            sage: Ha=HM(2,2,2,'a')
            sage: Ha.printHM()
            [:, :, 0]=
            [a000 a010]
            [a100 a110]
            <BLANKLINE>
            [:, :, 1]=
            [a001 a011]
            [a101 a111]
            <BLANKLINE>


        AUTHORS:

        - Edinah K. Gnang
        """
        if self.order()==2:
            L=self.listHM()
            print('[:, :]=\n'+Matrix(SR,L).str())
        elif self.order()==3:
            L=self.listHM()
            for dpth in range(self.n(2)):
                print('[:, :, '+str(dpth)+']=\n'+Matrix(SR,self.n(0),self.n(1),[L[i][j][dpth] for i in range(self.n(0)) for j in range(self.n(1))]).str()+'\n')
        else:
            raise ValueError("not supported for order %d hypermatrices" %self.order())

    def latexHM(self):
        """
        Returns the latex string. 

        EXAMPLES:

        ::

            sage: Ha=HM(2,2,'a')
            sage: len(Ha.latexHM())
            86


        AUTHORS:

        - Edinah K. Gnang
        """
        # Initialization of the string to
        strng=""
        if self.order()==2:
            L=self.listHM()
            strng=strng+'[:, :]=\n'+latex(Matrix(SR,L))
            return strng
        elif self.order()==3:
            L=self.listHM()
            for dpth in range(self.n(2)):
                strng=strng+'[:, :, '+str(dpth)+']=\n'+latex(Matrix(SR,self.n(0),self.n(1),[L[i][j][dpth] for i in range(self.n(0)) for j in range(self.n(1))]))+'\n'
            return strng
        else:
            raise ValueError("not supported for order %d hypermatrices" %self.order())

    def n(self,i):
        """
        Returns the i-th side length

        EXAMPLES:

        ::

            sage: Ha=HM(2,3,2,'a')
            sage: Ha.n(1)
            3


        AUTHORS:

        - Edinah K. Gnang
        """
        if i==0:
            return self.nrows()
        elif i==1:
            return self.ncols()
        elif i==2:
            return self.ndpts()
        else:
            tmp=self.listHM()
            for j in range(i):
                tmp=tmp[0]
            return len(tmp)

    def list(self):
        """
        Returns the hypermatrix as a list going down by columns.

        EXAMPLES:

        ::

            sage: Ha=HM(2,3,'a'); Ha.list()
            [a00, a10, a01, a11, a02, a12]


        AUTHORS:

        - Edinah K. Gnang
        """
        lst = []
        l = [self.n(i) for i in range(self.order())]
        # Main loop canonicaly listing the elements
        for i in range(prod(l)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entry = [Integer(mod(i,l[0]))]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            # Appending to the list
            lst.append(self[tuple(entry)])
        return lst

    def listHM(self):
        """
        Returns the hypermatrix as a list of list

        EXAMPLES:

        ::

            sage: Ha=HM(2,3,'a'); Ha.listHM()
            [[a00, a01, a02], [a10, a11, a12]]


        AUTHORS:

        - Edinah K. Gnang
        """
        return self.hm

    def matrix(self):
        """
        Returns the sage matrix class conversion of a
        the second order hypermatrix associated with self

        EXAMPLES:

        ::

            sage: Ha=HM(2,3,'a'); Ha.matrix()
            [a00 a01 a02]
            [a10 a11 a12]


        AUTHORS:

        - Edinah K. Gnang
        """
        if self.order()<=2:
            return Matrix(SR,self.listHM())
        else:
            raise ValueError("not supported for order %d hypermatrices" %self.order())

    def order(self):
        """
        Returns the hypermatrix order.

        EXAMPLES:

        ::

            sage: Ha=HM(2,3,'a'); Ha.order()
            2


        AUTHORS:

        - Edinah K. Gnang
        """
        cnt = 0
        H = self.listHM()
        while type(H) == type([]):
            H = H[0]
            cnt = cnt+1
        return cnt

    def dimensions(self):
        """
        Returns the list of hypermatrix side lengths.

        EXAMPLES:

        ::

            sage: Ha=HM(2,3,'a'); Ha.dimensions()
            [2, 3]


        AUTHORS:

        - Edinah K. Gnang
        """
        return [self.n(i) for i in range(self.order())]

    def zero_pad(self, dimLst):
        """
        return the zero padding of the input hypermatrix the inputs dimLst
        is desired dimensions of the padding        

        EXAMPLES::
        
            sage: A = HM(2,3,'a'); B=A.zero_pad([4,4])
            sage: B.printHM()
            [:, :]=
            [a00 a01 a02   0]
            [a10 a11 a12   0]
            [  0   0   0   0]
            [  0   0   0   0]
            sage: U = HM(2,1,'u'); V=U.zero_pad([2,2,2])
            sage: V.printHM()
            [:, :, 0]=
            [u00   0]
            [u10   0]
            <BLANKLINE>
            [:, :, 1]=
            [0 0]
            [0 0] 
        """
        if len(dimLst) == self.order():
            l = self.dimensions()
            AtmpL=[max(dimLst[z],self.n(z)) for z in rg(self.order())]+['zero']
            Rh = HM(*AtmpL)
            # Main loop performing the transposition of the entries
            for i in range(prod(l)):
                # Turning the index i into an hypermatrix array location using the decimal encoding trick
                entry = [Integer(mod(i,l[0]))]
                sm = Integer(mod(i,l[0]))
                for k in range(len(l)-1):
                    entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                    sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                booldim = True
                for j in rg(self.order()):
                    if entry[j] >= self.n(j):
                        booldim = False
                        break
                if  booldim:
                    Rh[tuple(entry)]=self[tuple(entry)]
                else:
                    Rh[tuple(entry)]=0
            return Rh
        elif len(dimLst) > self.order():
            AtmpL=[max(dimLst[z],self.n(z)) for z in self.dimensions()]+[dimLst[z] for z in rg(self.order(),len(dimLst))]+['zero']
            Rh = HM(*AtmpL)
            l = Rh.dimensions()
            # Main loop performing the transposition of the entries
            for i in range(prod(l)):
                # Turning the index i into an hypermatrix array location using the decimal encoding trick
                entry = [Integer(mod(i,l[0]))]
                sm = Integer(mod(i,l[0]))
                for k in range(len(l)-1):
                    entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                    sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                booldim = True
                for j in rg(self.order()):
                    if entry[j] >= self.n(j):
                        booldim = False
                        break
                if entry[self.order():] == [0 for z in rg(self.order(),len(dimLst))] and booldim:
                    Rh[tuple(entry)]=self[tuple(entry[:self.order()])]
                else:
                    Rh[tuple(entry)]=0
            return Rh
        else:
            raise ValueError("The order must not be smaller the starting hypermatrix")

    def fill_with(self, T, st=0):
        """
        returns a hypermatrix whose top left corner is replaced with the entries a same order
        but smaller size input hypermatrix T. This method generalizes slightly the zero padding
        method for hypermatrices of the same order. The input st which is set by default to 0
        is the shift parameter. It adds st to all the index and has the effec of translating the
        block along the main diagonal. 


        EXAMPLES::
        
            sage: A=HM(3,3,'a'); B=HM(2,2,'b') 
            sage: A.fill_with(B).printHM()
            [:, :]=
            [b00 b01 a02]
            [b10 b11 a12]
            [a20 a21 a22]
            sage: A=HM(3,3,'a'); B=HM(2,2,'b','shift') 
            sage: A.fill_with(B,1).printHM()
            [:, :]=
            [a00 a01 a02]
            [a10 b11 b12]
            [a20 b21 b22]
 
        """
        boolsize = True
        for d in rg(self.order()):
            if T.n(d)+st > self.n(d):
                boolsize = False
                raise ValueError("Expected the input hypermatrix to have larger dimensions in all directions")
        l = self.dimensions()
        Rh = self.copy()
        # Main loop performing the transposition of the entries
        for i in range(prod(l)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entry = [Integer(mod(i,l[0]))]
            sm = Integer(mod(i, l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            booldim = True
            for j in rg(self.order()):
                if entry[j] >= T.n(j):
                    booldim = False
                    break
            if booldim:
                Rh[tuple([z+st for z in entry])]=T[tuple(entry)]
        return Rh

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

    def support(self):
        return GeneralHypermatrixSupport(self)

    def append_index(self,indx):
        return GeneralHypermatrixAppendIndex(self,indx)

    def norm(self,p=2):
        if p==Infinity:
            return max([abs(i) for i in self.list() if not i.is_zero()])
        elif p==0:
            return sum([1 for i in self.list() if not i.is_zero()])
        elif p==1:
            return sum([abs(i) for i in self.list() if not i.is_zero()])
        elif p>1:
            return sum([abs(i)^p for i in self.list()])^(1/p)
        else:
            raise ValueError("Expect a real number greater or equal to 0 or Infinity")

    def trace(self):
        if self.order()==2:
            return sum(self[i,i] for i in range(min(self.n(0),self.n(1))))
        else:
            raise ValueError("Expects a second order hypermatrix")

    def reduce(self, VrbL, Rlts):
        return GeneralHypermatrixReduce(self, VrbL, Rlts)

    def fast_reduce(self, monom, subst):
        AtmpL=self.dimensions()+[[fast_reduce(self.list()[i], monom, subst) for i in rg(prod(self.dimensions()))]]
        return HM(*AtmpL)

    def operands(self):
        # Turning the hypermatrix into a list
        L=self.list()
        # Boolean variables checking that the operator is the same
        # for all entries.
        CommonOperator=True
        # max_len keeps track of the maximum number of operands.
        max_len=len(L[0].operands())
        for f in L[1:]:
            if f.operator() != L[0].operator():
                CommonOperator=False
                break
            else:
                if max_len < len(f.operands()):
                    max_len = len(f.operands())
        if CommonOperator==False:
            raise ValueError("Expects the same operator for all entries.")
        else:
            # Initilization of the list of operands
            Lst=[f.operands()+[0 for i in rg(max_len-len(f.operands()))] for f in L]
            # Initialization of the list which will store the hypermatrix operands
            LstOprds=[]
            # Looping over operands. The operands order is the default
            # sage lexicographic order
            for i in rg(max_len):
                AtmpL=self.dimensions()+[[l[i] for l in Lst]]
                LstOprds.append(HM(*AtmpL))  
            return LstOprds

    def mod(self, m):
        AtmpL=self.dimensions()+[[Integer(mod(self.list()[i], m)) for i in rg(prod(self.dimensions()))]]
        return HM(*AtmpL )

    def adjugate(self):
        if self.order()==2:
            return Matrix2HM((self.matrix()).adjugate())
        else:
            raise ValueError("Expects a second order hypermatrix")

    def per(self):
        if self.is_cubical():
            if self.n(0)==1:
                return self.list()[0]
            elif self.order()==2:
                return Per(self)
        else:
            raise ValueError("Expects a cubical second order hypermatrix")

    def det(self):
        if self.is_cubical():
            if self.n(0)==1:
                return self.list()[0]
            elif self.n(0)==2:
                return general_side_length_2_det(self)
            elif self.order()==2:
                return Deter(self)
            elif self.order()==3:
                return ThirdOrderDeter(self)
            elif self.order()==4:
                return FourthOrderHyperdeterminant(self)
            elif self.order()==5:
                return FifthOrderHyperdeterminant(self)
            elif self.order()==6:
                return SixthOrderHyperdeterminant(self)
            else:
                return GeneralHyperdeterminant(self) 
        else:
            raise ValueError("Expects a cubical hypermatrix")

    def is_empty(self):
        return self.hm==[]

    def is_zero(self):
        if Set([f.is_zero() for f in self.list()])==Set([True]):
            return True
        else:
            return False

    def is_unit(self):
        if self.n(0)==self.n(1) and self.n(0)==1 and self.list()==[1]:
            return True
        elif self.order()==2 and self.n(0)==self.n(1) and self==HM(self.order(),self.n(1),'kronecker'):
            return True
        else:
            return False

    def is_symmetric(self):
        return (self-self.transpose()).is_zero()

    def is_cubical(self):
        return len(Set(self.dimensions()).list())==1

    def is_hypercolumn(self):
        return Set(self.dimensions()[1:])==Set([1])

    def ref(self):
        if self.order()==2:
            return gaussian_eliminationHMIV(self)
        else:
            raise ValueError("Expected a second order hypermatrix")

    def refII(self):
        if self.order()==2:
            return gaussian_eliminationHMIII(self, HM(self.n(0),1,'zero'))[0]
        else:
            raise ValueError("Expected a second order hypermatrix")

    def rref(self):
        if self.order()==2:
            return gauss_jordan_eliminationHMII(self, HM(self.n(0),1,'zero'))[0]
        elif self.order()==3 and self.is_cubical():
            # Initializing the size parameter
            sz=self.n(0)
            # initializing the list second order slices
            TmpL=[HM(sz, sz, [self[i,j,k] for j in range(sz) for i in range(sz)]).rref() for k in range(sz)]
            # Initilizing the final second order hypermatrix
            M=HM(sz, sz, [TmpL[k][i,i] for i in range(sz) for k in range(sz)]).rref()
            # Initializing and filling up the resulting hypermatrix
            Rslt=HM(sz,sz,sz,'zero')
            for i in range(sz):
                for k in range(sz):
                    Rslt[i,i,k]=M[k,i]
            return Rslt
        else:
            raise ValueError("Expected a second order hypermatrix or a cubic third order hypermatrix")

    def rrefII(self):
        if self.order()==2:
            # Initiallization of the size parameter
            sz=max(self.dimensions())
            # Initialization of the identity matrix
            Id=HM(2,sz,'kronecker')
            # Initialization of the permutation matrix
            Hp=sum(HM(sz,1,[Id[i,sz-1-k] for i in range(sz)])*HM(1,sz,[Id[k,j] for j in range(sz)]) for k in range(sz))
            return (Hp.transpose()*((Hp.transpose()*ZeroPadding(self.refII())*Hp).refII())*Hp).refII()
        else:
            raise ValueError("Expected a second order hypermatrix")

    def rank(self):
        if self.order()==2:
            if self.is_zero():
                return 0
            else:
                cnt=0
                TmpHm=self.ref().copy()
                for i in range(min(TmpHm.dimensions())):
                    if not HM(1,self.n(1),[TmpHm[i,j] for j in range(self.n(1))]).is_zero():
                        cnt=cnt+1
                return cnt
        else:
            raise ValueError("Expected a second order hypermatrix")

    def matrix_from_rows_and_columns(self, rows, columns):
        """
        Return the matrix constructed from self from the given rows and
        columns.
    
        EXAMPLES::
    
            sage: A = HM(3,3,rg(3^2)); A.printHM()
            [:, :]=
            [0 3 6]
            [1 4 7]
            [2 5 8]
            sage: A.matrix_from_rows_and_columns([1], [0,2]).printHM()
            [:, :]=
            [1 7]
            sage: A.matrix_from_rows_and_columns([1,2], [1,2]).printHM()
            [:, :]=
            [4 7]
            [5 8]
    
        Note that row and column indices can be reordered or repeated::
    
            sage: A.matrix_from_rows_and_columns([2,1], [2,1]).printHM()
            [:, :]=
            [8 5]
            [7 4]
    
        For example here we take from row 1 columns 2 then 0 twice, and do
        this 3 times.
    
        ::
    
            sage: A.matrix_from_rows_and_columns([1,1,1],[2,0,0]).printHM()
            [:, :]=
            [7 1 1]
            [7 1 1]
            [7 1 1]

 
        AUTHORS:
    
        - Edinah K. Gnang
        """
        tMp=self.matrix().matrix_from_rows_and_columns(rows, columns)
        return HM(tMp.nrows(), tMp.ncols(), tMp.transpose().list())
 
    def third_order_hypermatrix_from_rows_columns_and_depths(self, rows, columns, depths):
        """
        Return the third order hypermatrix constructed from self from the given rows
        columns and depths.

        EXAMPLES::

            sage: A = HM(3,3,3,rg(3^3)); A.printHM()
            [:, :, 0]=
            [0 3 6]
            [1 4 7]
            [2 5 8]
            <BLANKLINE>
            [:, :, 1]=
            [ 9 12 15]
            [10 13 16]
            [11 14 17]
            <BLANKLINE>
            [:, :, 2]=
            [18 21 24]
            [19 22 25]
            [20 23 26]
            <BLANKLINE>
            sage: A.third_order_hypermatrix_from_rows_columns_and_depths([1], [0,2], [0,2]).printHM()
            [:, :, 0]=
            [1 7]
            <BLANKLINE>
            [:, :, 1]=
            [19 25]
            <BLANKLINE>
            sage: A.third_order_hypermatrix_from_rows_columns_and_depths([1,2], [1,2], [1,2]).printHM()
            [:, :, 0]=
            [13 16]
            [14 17]
            <BLANKLINE>
            [:, :, 1]=
            [22 25]
            [23 26]
            <BLANKLINE>

        Note that row and column indices can be reordered or repeated::

            sage: A.third_order_hypermatrix_from_rows_columns_and_depths([2,1], [2,1], [2,1]).printHM()
            [:, :, 0]=
            [26 23]
            [25 22]
            <BLANKLINE>
            [:, :, 1]=
            [17 14]
            [16 13]
            <BLANKLINE>

        For example here we take from row 1 columns 2 then 0 twice, and do
        this 3 times.

        ::

            sage: A.third_order_hypermatrix_from_rows_columns_and_depths([1,1,1],[2,0,0],[2,0,0]).printHM()
            [:, :, 0]=
            [25 19 19]
            [25 19 19]
            [25 19 19]
            <BLANKLINE>
            [:, :, 1]=
            [7 1 1]
            [7 1 1]
            [7 1 1]
            <BLANKLINE>
            [:, :, 2]=
            [7 1 1]
            [7 1 1]
            [7 1 1]
            <BLANKLINE>

        AUTHORS:

        - Edinah K. Gnang
        """
        if self.order()==3:
            if not type(rows)==list:
                raise TypeError("rows must be a list of integers")
            if not type(columns)==list:
                raise TypeError("columns must be a list of integers")
            if not type(depths)==list:
                raise TypeError("depths must be a list of integers")
            A = HM(len(rows), len(columns), len(depths), 'zero')
            r = 0 
            for i in rows:
                k = 0
                for j in columns:
                    t = 0
                    for l in depths:
                        A[r,k,t]=self.hm[i][j][l]
                        t += 1
                    k += 1
                r += 1
            return A
        else:
            raise ValueError("Expected a third order hypermatrix")

    def side_length_subhypermatrices(self, k):
        r"""
        Return the list of all side length hypermatrices.
        Only supported for second and third order hypermatrices. 
    
        In the case of matrix The returned list is sorted in lexicographical
        row major ordering, e.g., if A is a `3 \times 3` matrix then the minors
        returned are with these rows/columns: [ [0, 1]x[0, 1], [0, 1]x[0, 2], [0, 1]x[1, 2],
        [0, 2]x[0, 1], [0, 2]x[0, 2], [0, 2]x[1, 2], [1, 2]x[0, 1], [1, 2]x[0, 2], [1, 2]x[1, 2] ].
    
        INPUT:
    
        - ``k`` -- integer
    
        EXAMPLES::
    
            sage: A = HM(2,3,[1,2,3,4,5,6]); A.printHM()
            [:, :]=
            [1 3 5]
            [2 4 6]
            sage: [m.det() for m in A.side_length_subhypermatrices(2)]
            [-2, -4, -2]
            sage: [m.det() for m in A.side_length_subhypermatrices(1)]
            [1, 3, 5, 2, 4, 6]


        AUTHORS:

        - Edinah K. Gnang (2016-12-11)

        """
        if self.order() == 2:
            from sage.combinat.combination import Combinations
            all_rows = range(self.nrows())
            all_cols = range(self.ncols())
            m = []
            for rows in Combinations(all_rows,k):
                for cols in Combinations(all_cols,k):
                    m.append(self.matrix_from_rows_and_columns(rows,cols))
            return m
        elif self.order() == 3:
            from sage.combinat.combination import Combinations
            all_rows = range(self.nrows())
            all_cols = range(self.ncols())
            all_dpts = range(self.ndpts())
            m = []
            for rows in Combinations(all_rows,k):
                for cols in Combinations(all_cols,k):
                    for dpts in Combinations(all_dpts,k):
                        m.append(self.third_order_hypermatrix_from_rows_columns_and_depths(rows, cols, dpts))
            return m
        else :
            raise ValueError("Not supported for hypermatrices of this order")

    def matrix_minor(self,u,v):
        """
        Outputs a second order minor hypermatrix
        where the row u and column v is removed

        EXAMPLES:
 
        ::

            sage: HM(3,3,'a').matrix_minor(0,0).printHM()
            [:, :]=
            [a11 a12]
            [a21 a22]

        AUTHORS:
        - Edinah K. Gnang
        - To Do: 
        """
        if self.order()==2 and self.n(0)>1 and self.n(1)>1:
            return HM(self.n(0)-1,self.n(1)-1, [self[i,j] for j in rg(self.n(1)) for i in rg(self.n(0)) if j!=v and i!=u])
        else :
            raise ValueError("Not supported for hypermatrices of order > 2 and evry dimension must exceed 1")

    def slice(self,L,strg):
        """
        Outputs the result of the slicifing 


        EXAMPLES:
 
        ::

            sage: HM(3,3,'a').slice([0], 'row').printHM()
            [:, :]=
            [a00 a01 a02]
            <BLANKLINE>
            sage: HM(3,3,'a').slice([1], 'col').printHM()
            [:, :]=
            [a01]
            [a11]
            [a21]
            <BLANKLINE>
            sage: HM(3,3,3,'a').slice([2], 'dpt').printHM()
            [:, :, 0]=
            [a002 a012 a022]
            [a102 a112 a122]
            [a202 a212 a222] 
            <BLANKLINE>
  

        AUTHORS:
        - Edinah K. Gnang
        - To Do: 
        """
        return GeneralHypermatrixSlicer(self, L, strg)

    def flatten(self, Rg, ord):
        """
        Outputs a lower order flattened hypermatrix of the higher order input.             
        The first input is the hypermatrix to be flattened. The second input 
        corresponds to the indices that will remain unchanged after flattening. 
        The third input is the order of the desired output hypermatrix. When 
        flatten an higher order input hypermatrix to order 1, the result is a list.
        
        In the example of A.flatten([0, 1], 2) for flattening 
        the order 3 hypermatrix A into a matrix, the input [0, 1] tells us that the
        column and the row indices will be left unchanged, and we are stacking depth
        slices along the row direction. 
        
        The flattening of hypermatrices with ord>2 is obtained recursively. 


        EXAMPLES:
 
        ::


            sage: sz=3; A=HM(sz, sz, sz, 'a')
            sage: A.flatten([0, 1], 2).dimensions()
            [3, 9]
            sage: A.flatten([0, 1], 2).printHM()
            [:, :]=
            [a000 a010 a020 a001 a011 a021 a002 a012 a022]
            [a100 a110 a120 a101 a111 a121 a102 a112 a122]
            [a200 a210 a220 a201 a211 a221 a202 a212 a222]
            sage: sz=2; B=HM(sz,sz,sz,'b')
            sage: B.flatten([2], 1)
            [b000, b001, b100, b101, b010, b011, b110, b111]
  

        AUTHORS:
        - Edinah K. Gnang and Fan Tian
        - To Do: 
        """
        return GeneralHypermatrixFlatten(self, Rg, ord)

def MatrixGenerate(nr, nc, c):
    """
    Generates a list of lists associated with a symbolic nr x nc
    matrix by indexing the input character c by the indices.

    EXAMPLES:

    ::

        sage: M = MatrixGenerate(2, 2, 'm'); M
        [[m00, m01], [m10, m11]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_q_rows = nr; n_q_cols = nc
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
        raise ValueError("Input dimensions "+str(nr)+" and "+str(nc)+" must both be non-zero positive integers.")

def SymMatrixGenerate(nr, c):
    """
    Generates a list of lists associated with a symbolic nr x nr
    symmetric matrix by indexing the input character c by indices.

    EXAMPLES:

    ::

        sage: M = SymMatrixGenerate(2, 'm'); M
        [[m00, m01], [m01, m11]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_q_rows = nr; n_q_cols = nr
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
        raise ValueError("Input dimensions "+str(nr)+" must be a non-zero positive integers.")

def SkewSymMatrixGenerate(nr, c):
    """
    Generates a list of lists associated with a symbolic nr x nr
    symmetric matrix by indexing the input character c by indices.

    EXAMPLES:

    ::

        sage: M = SkewSymMatrixGenerate(2, 'm'); M
        [[0, m01], [-m01, 0]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Setting the dimensions parameters.
    n_q_rows = nr; n_q_cols = nr
    # Test for dimension match
    if n_q_rows > 0 and n_q_cols > 0:
        # Initialization of the hypermatrix
        q = []
        for i in range(n_q_rows):
            q.append([])
        for i in range(len(q)):
            for j in range(n_q_cols):
                # Filling up the matrix
                if i==j:
                    (q[i]).append(0)
                elif i>j:
                    (q[i]).append(-var(c+str(min(i,j))+str(max(i,j))))
                else:
                    (q[i]).append(var(c+str(min(i,j))+str(max(i,j))))
        return q
    else :
        raise ValueError("Input dimensions "+str(nr)+" must be a non-zero positive integers.")

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
    return [HypermatrixGenerate(*(args[1:-1]+(args[-1]+str(i),))) for i in range(args[0])]

def HypermatrixGenerateII(*args):
    """
    Generates a list of lists associated with a symbolic arbitrary
    order hypematrix of size specified by the input using and the
    entries are determined by the last input character. The difference
    between this function and the previous one is the fact that the index
    start from 1 as apposed to starting from zero. This was motivated by
    general determinant implementation.

    EXAMPLES:

    ::

        sage: M = HypermatrixGenerateII(2, 2, 2, 'm'); M
        [[[m111, m112], [m121, m122]], [[m211, m212], [m221, m222]]]


    AUTHORS:
    - Edinah K. Gnang, Ori Parzanchevski and Yuval Filmus
    """
    if len(args) == 1:
        return var(args[0])
    return [HypermatrixGenerateII( *(args[1:-1]+(args[-1]+str(i),)) ) for i in range(1,1+args[0])]

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
        return [SR(1) for i in rg(args[0])]
    AtmpL=args[1:]
    return [HypermatrixGenerateAllOne(*AtmpL) for i in rg(args[0])]

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
        return [SR(0) for i in range(args[0])]
    AtmpL=args[1:]
    return [HypermatrixGenerateAllZero(*AtmpL) for i in range(args[0])]

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
        raise ValueError("Input dimensions "+str(nr)+" must be a non-zero positive integer.")

def SkewSymHypermatrixGenerate(nr, c):
    """
    Generates a list of lists associated with a symbolic third order hypermatrix of size
    nr x nc x nd third order hypematrix using the input character c followed by indices.

    EXAMPLES:

    ::

        sage: M = SkewSymHypermatrixGenerate(3, 'm'); M
        [[[0, 0, 0], [0, 0, m012], [0, m021, 0]],
         [[0, 0, m021*w^2], [0, 0, 0], [m012*w, 0, 0]],
         [[0, m012*w^2, 0], [m021*w, 0, 0], [0, 0, 0]]]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initialization of the variable to be subsituted
    w=var('w')
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
                        (q[i][j]).append(0*var(c+str(min(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))+str(max(i,j,k))))
                    else:
                        if i == min(i,j,k) and k == max(i,j,k):
                            (q[i][j]).append(w^0*var(c+str(min(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))+str(max(i,j,k))))
                        elif k == min(i,j,k) and j == max(i,j,k):
                            (q[i][j]).append(w^1*var(c+str(min(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))+str(max(i,j,k))))
                        elif i == max(i,j,k) and j == min(i,j,k):
                            (q[i][j]).append(w^2*var(c+str(min(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))+str(max(i,j,k))))

                        elif i == min(i,j,k) and j == max(i,j,k):
                            (q[i][j]).append(w^0*var(c+str(min(i,j,k))+str(max(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))))
                        elif k == min(i,j,k) and i == max(i,j,k):
                            (q[i][j]).append(w^1*var(c+str(min(i,j,k))+str(max(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))))
                        elif k == max(i,j,k) and j == min(i,j,k):
                            (q[i][j]).append(w^2*var(c+str(min(i,j,k))+str(max(i,j,k))+str(i+j+k-min(i,j,k)-max(i,j,k))))
        return q
    else :
        raise ValueError("Input dimensions "+str(nr)+" must be a non-zero positive integer.")

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
        raise ValueError("The Dimensions non zero.")

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
        raise ValueError("The Dimensions of the input hypermatrices must match.")

def HypermatrixHadamardProduct(A, B):
    """
    Outputs a list of lists associated with the Hadamard product of
    the two input third order hypermatrices A and B

    EXAMPLES:

    ::

        sage: A = HypermatrixGenerate(2,2,2,'a'); B = HypermatrixGenerate(2,2,2,'b')
        sage: Rslt = HypermatrixHadamardProduct(A, B); Rslt
        [[[a000*b000, a001*b001], [a010*b010, a011*b011]],
         [[a100*b100, a101*b101], [a110*b110, a111*b111]]]


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
        raise ValueError("The Dimensions of the input hypermatrices must match.")

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

        sage: HypermatrixProduct(HM(2,2,2,'a'), HM(2,2,2,'b'), HM(2,2,2,'c'))
        [[[a000*b000*c000 + a010*b001*c100, a001*b000*c001 + a011*b001*c101], [a000*b010*c010 + a010*b011*c110, a001*b010*c011 + a011*b011*c111]], [[a100*b100*c000 + a110*b101*c100, a101*b100*c001 + a111*b101*c101], [a100*b110*c010 + a110*b111*c110, a101*b110*c011 + a111*b111*c111]]]


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
        raise ValueError("Hypermatrix dimension mismatch.")

def HypermatrixProductII(A, B, C, support):
    """
    Outputs a list of lists associated with the ternary
    Bhattacharya-Mesner product of the input third order 
    hypermatrices A, B and C. But we restrict the sum in products
    to the support 
    The code writen here handles both list of list data
    structures as well as the Hypermatrix HM class objects.

    EXAMPLES:

    ::

        sage: M = HypermatrixProductII(HM(2,2,2,'a'), HM(2,2,2,'b'), HM(2,2,2,'c'), [0]); M
        [[[a000*b000*c000, a001*b000*c001], [a000*b010*c010, a001*b010*c011]], [[a100*b100*c000, a101*b100*c001], [a100*b110*c010, a101*b110*c011]]]

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
                    (q[i][j]).append(sum([A[i][l][k]*B[i][j][l]*C[l][j][k] for l in support]))
        if typ=='list':
            return q
        else:
            return HM(q)

    else :
        raise ValueError("Hypermatrix dimension mismatch.")

def HypermatrixLogProduct(A,B,C):
    """
    Outputs a list of lists associated with the ternary
    product of the input hypermatrices A, B and C.

    EXAMPLES:

    ::

        sage: Ha=HM(2,1,2,'a'); Hb=HM(2,2,1,'b'); Hc=HM(1,2,2,'c'); HypermatrixLogProduct(Ha,Hb,Hc)
        [[[a000 + b000 + c000, a001 + b000 + c001], [a000 + b010 + c010, a001 + b010 + c011]], [[a100 + b100 + c000, a101 + b100 + c001], [a100 + b110 + c010, a101 + b110 + c011]]]

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
        raise ValueError("Hypermatrix dimension mismatch.")

def HypermatrixKroneckerProduct(A, B, C):
    """
    Outputs a list of lists associated with the ternary
    analog of the Kronecker product for the inputh third
    order hypermatrices A, B and C.

    EXAMPLES:

    ::

        sage: X=HM(2, 1, 1, HM(2, 'x').list()); Y=HM(1, 2, 1, HM(2, 'y').list()); Z=HM(1, 1, 2, HM(2, 'z').list())
        sage: T=HypermatrixKroneckerProduct(X, Y, Z); T
        [[[x0*y0*z0, x0*y0*z1], [x0*y1*z0, x0*y1*z1]], [[x1*y0*z0, x1*y0*z1], [x1*y1*z0, x1*y1*z1]]]

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
        raise ValueError("Hypermatrix dimension mismatch.")

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
        raise ValueError("Hypermatrix dimension mismatch.")

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
    - Edinah K. Gnang
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
        raise ValueError("Hypermatrix dimension mismatch.")

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
        raise ValueError("Input dimensions "+str(nr)+" must be a non-zero positive integer.")

def Vandermonde(l):
    """
    Constructs a Vandermonde matrix from the input list
    assumed to be either numbers or symbolic variables
    nothing breaks however if one presents as input a list of
    hypermatrices.

    EXAMPLES:

    ::

        sage: Vandermonde(var_list('x',2))
        [[1, 1], [x0, x1]]


    AUTHORS:
    - Edinah K. Gnang
    """
    return HM(len(l),len(l),[l[j]^i for j in range(len(l)) for i in range(len(l))])

def var_list(c,sz):
    """
    Returns a variable list of size sz indexing the input character c.

    EXAMPLES:

    ::

        sage: var_list('x', 3)
        [x0, x1, x2]


    AUTHORS:
    - Edinah K. Gnang
    """
    return HM(sz,c).list()

def Lagrange(X):
    """
    Constructs a Lagrange matrix from the input list
    assumed to be either numbers or symbolic variables
    nothing breaks however if one presents as input a 
    list of hypermatrices.

    EXAMPLES:

    ::

        sage: Lagrange(HM(3,'x').list())
        [[x0 - x2, x0 - x1], [x1 - x2, 1]]


    AUTHORS:
    - Edinah K. Gnang
    """
    sz=len(X)-1; L=[]
    for i in range(sz):
        tmpl=range(sz,i,-1)
        L.append([X[i]-X[j] for j in tmpl]+[1 for k in range(sz-len(tmpl))]) 
    return HM(L) 

def HypermatrixPermutation(s):
    """
    Generates a list of lists associated with a transposition 
    hypermatrix deduced from sigma. Note that as a result of 
    the  non associativity, permutations must be performed as
    one transposition at a time. This is one way of implementing
    permutation hypermatrices in this setting only non overlaping
    cycles can be combined into a single permutation hypermatrix


    EXAMPLES:


    ::


        sage: P = HypermatrixPermutation([0,2,1]); P
        [[[1, 0, 0], [0, 0, 1], [0, 1, 0]], [[1, 0, 0], [0, 0, 1], [0, 1, 0]], [[1, 0, 0], [0, 0, 1], [0, 1, 0]]]
        sage: P = HM(HypermatrixPermutation([1,0])); A=HM(2,2,2,'a')
        sage: Prod(P,A,P.transpose(2)).printHM()
        [:, :, 0]=
        [a001 a011]
        [a101 a111]
        <BLANKLINE>
        [:, :, 1]=
        [a000 a010]
        [a100 a110]


    AUTHORS:

    - Edinah K. Gnang and Ori Parzanchevski
    """
    sz = len(s)
    # Setting the dimensions parameters.
    n_q_rows = sz; n_q_cols = sz; n_q_dpts = sz
    # Test for dimension match
    if n_q_rows > 0 and n_q_cols > 0 and n_q_dpts >0:
        # Initialization of the hypermatrix
        q = []
        T = HypermatrixKroneckerDelta(sz)
        U = HypermatrixGenerateAllOne(sz,sz,sz)
        Id= HypermatrixProduct(U,U,T)
        Id= HypermatrixCyclicPermute(Id)
        for i in range(sz):
            q.append(Id[s[i]])
        return HypermatrixCyclicPermute(HypermatrixCyclicPermute(q))
    else :
        raise ValueError("Input dimensions "+str(sz)+" must be a non-zero positive integer.")

def DiagonalHypermatrix(Mtrx):
    """
    Outputs a diagonal third order hypermatrix
    constructed using the input square matrix.
    We enforce the symmetry constraint by only
    taking entries from the lower triangular
    part of the input matrix.

     EXAMPLES:

    ::

        sage: a00, a11, a01=var('a00, a11, a01'); Mtrx = Matrix(SR, [[a00, a01], [a01, a11]])
        sage: Dg = DiagonalHypermatrix(Mtrx); Dg
        [[[a00, 0], [0, a01]], [[a01, 0], [0, a11]]]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initialization of the dimensions
    sz = min(Mtrx.nrows(), Mtrx.ncols())
    n_d_rows = sz; n_d_cols = sz; n_d_dpts = sz
    # Initialization of the identity permutations hypermatrix
    D = HypermatrixPermutation(range(sz))
    # Filling up the entries of the hypermatrix.
    for i in range(n_d_rows):
        for j in range(n_d_cols):
            for k in range(n_d_dpts):
                if D[i][j][k] != 0:
                    D[i][j][k] = Mtrx[i,k]
    return D

def Orthogonal2x2x2Hypermatrix(t,x,y):
    """
    Outputs a symbolic parametrization of third order orthogonal hypermatrix
    of size 2x2x2.

     EXAMPLES:

    ::

        sage: t,x,y=var('t,x,y')
        sage: Orthogonal2x2x2Hypermatrix(t,x,y)
        [[[sin(t)^(2/3), -x*cos(t)^(2/3)], [cos(t)^(2/3), y*sin(t)^(2/3)]], [[1/x, sin(t)^(2/3)], [1/y, cos(t)^(2/3)]]]
        sage: U=Orthogonal2x2x2Hypermatrix(t,x,y); Prod(U,U.transpose(2),U.transpose())
        [[[cos(t)^2 + sin(t)^2, 0], [0, 0]], [[0, 0], [0, cos(t)^2 + sin(t)^2]]]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    return HM([[[sin(t)^(2/3), -x*cos(t)^(2/3)], [cos(t)^(2/3), y*sin(t)^(2/3)]], [[1/x, sin(t)^(2/3)], [1/y, cos(t)^(2/3)]]])

def Orthogonal2x2x2HypermatrixII(t,x,y):
    """
    Outputs a symbolic parametrization of third order orthogonal hypermatrix
    of size 2x2x2.

     EXAMPLES:

    ::

        sage: t,x,y=var('t,x,y')
        sage: Orthogonal2x2x2HypermatrixII(t,x,y)
        [[[sin(t)^(2/3), y*cos(t)^(2/3)], [cos(t)^(2/3), -x*sin(t)^(2/3)]], [[1/y, sin(t)^(2/3)], [1/x, cos(t)^(2/3)]]]
        sage: U=Orthogonal2x2x2HypermatrixII(t,x,y); Prod(U,U.transpose(2),U.transpose())
        [[[cos(t)^2 + sin(t)^2, 0], [0, 0]], [[0, 0], [0, cos(t)^2 + sin(t)^2]]]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initialization of the Permutation hypermatrix
    P,Q=HypermatrixSn([1,0])
    return Prod(HM([[[cos(t)^(2/3), -x*sin(t)^(2/3)], [sin(t)^(2/3), y*cos(t)^(2/3)]], [[1/x, cos(t)^(2/3)], [1/y, sin(t)^(2/3)]]]), Q.transpose(), P.transpose())

def Orthogonal2x2x2HypermatrixIII(t,x,y):
    """
    Outputs a symbolic parametrization of third order orthogonal hypermatrix
    of size 2x2x2.

     EXAMPLES:

    ::

        sage: t,x,y=var('t,x,y')
        sage: Orthogonal2x2x2HypermatrixIII(t,x,y)
        [[[sin(t)^(2/3), x*cos(t)^(2/3)], [cos(t)^(2/3), -y*sin(t)^(2/3)]], [[1/x, sin(t)^(2/3)], [1/y, cos(t)^(2/3)]]]
        sage: U=Orthogonal2x2x2HypermatrixIII(t,x,y); Prod(U,U.transpose(2),U.transpose())
        [[[cos(t)^2 + sin(t)^2, 0], [0, 0]], [[0, 0], [0, cos(t)^2 + sin(t)^2]]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    return HM([[[sin(t)^(2/3), x*cos(t)^(2/3)], [cos(t)^(2/3), -y*sin(t)^(2/3)]], [[1/x, sin(t)^(2/3)], [1/y, cos(t)^(2/3)]]])


def Orthogonal3x3x3Hypermatrix(t1,t2):
    """
    Outputs a symbolic parametrization of third order orthogonal hypermatrix
    of size 3x3x3.

     EXAMPLES:

    ::

        sage: t1, t2=var('t1, t2')
        sage: Orthogonal3x3x3Hypermatrix(t1,t2)
        [[[cos(t1)^(2/3), cos(t2)^(2/3)*sin(t1)^(2/3), 0],
          [cos(t2)^(2/3)*sin(t1)^(2/3), sin(t1)^(2/3)*sin(t2)^(2/3), 0],
          [sin(t1)^(2/3)*sin(t2)^(2/3), -1/2*(I*sqrt(3) + 1)*cos(t1)^(2/3), 0]],
         [[sin(t1)^(2/3)*sin(t2)^(2/3),
           cos(t1)^(2/3),
           -1/2*(I*sqrt(3) + 1)*cos(t2)^(2/3)*sin(t1)^(2/3)],
          [-1/2*(-I*sqrt(3) + 1)*cos(t1)^(2/3),
           cos(t2)^(2/3)*sin(t1)^(2/3),
           sin(t1)^(2/3)*sin(t2)^(2/3)],
          [cos(t2)^(2/3)*sin(t1)^(2/3), sin(t1)^(2/3)*sin(t2)^(2/3), cos(t1)^(2/3)]],
         [[0, sin(t1)^(2/3)*sin(t2)^(2/3), cos(t1)^(2/3)],
          [0, cos(t1)^(2/3), cos(t2)^(2/3)*sin(t1)^(2/3)],
          [0,
           -1/2*(-I*sqrt(3) + 1)*cos(t2)^(2/3)*sin(t1)^(2/3),
           sin(t1)^(2/3)*sin(t2)^(2/3)]]]

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    c1=cos(t1)^(2/3)
    s1=sin(t1)^(2/3)
    c2=cos(t2)^(2/3)
    s2=sin(t2)^(2/3)
    return [[[c1,s1*c2,0],[s1*c2,s1*s2,0],[s1*s2,exp(-I*2*pi/3)*c1,0]], [[s1*s2,c1,exp(-I*2*pi/3)*s1*c2],[exp(I*2*pi/3)*c1,s1*c2,s1*s2], [s1*c2,s1*s2,c1]],[[0,s1*s2,c1],[0,c1,s1*c2],[0,exp(I*2*pi/3)*s1*c2,s1*s2]]]

def HypermatrixCayleyHamiltonStringList(n,c):
    """
    Outpts a list of strings describint hypermatrix powers of degree n. 

     EXAMPLES:

    ::

        sage: HypermatrixCayleyHamiltonStringList(5,'A')
        ['Prod(A, A, Prod(A, A, A))',
         'Prod(A, Prod(A, A, A), A)',
         'Prod(Prod(A, A, A), A, A)']

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    if n == 1:
        return [c]
    else:
        gu = []
        for i in range(1,n,2):
            for j in range(1,n-i,2):
                gu = gu + ["Prod("+g1+", "+g2+", "+g3+")" for g1 in HypermatrixCayleyHamiltonStringList(i,c) for g2 in HypermatrixCayleyHamiltonStringList(j,c) for g3 in HypermatrixCayleyHamiltonStringList(n-(i+j),c)]
        return gu

def HypermatrixCayleyHamiltonList(A,n):
    """
    Outpts a list of hypermatrices (each of which is encapsulated as a single list) of all product composition of degree n.
    This function is most adapted for constructing matrices.

     EXAMPLES:

    ::

        sage: A=HM(2,2,2,'a'); HypermatrixCayleyHamiltonList(A,3)
        [[a000^3 + a001*a010*a100,
          a000*a100^2 + a100*a101*a110,
          a000*a010^2 + a010*a011*a110,
          a010*a100*a110 + a110^2*a111,
          a000*a001^2 + a001*a011*a101,
          a001*a100*a101 + a101^2*a111,
          a001*a010*a011 + a011^2*a111,
          a011*a101*a110 + a111^3]]

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
        sage: Lst = HypermatrixCayleyHamiltonListII(A,3); Lst
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
                gu = gu + [HypermatrixProduct(g1,g2,g3) for g1 in HypermatrixCayleyHamiltonListII(A,i) for g2 in HypermatrixCayleyHamiltonListII(A,j) for g3 in HypermatrixCayleyHamiltonListII(A,n-(i+j))]
        return gu

def HypermatrixCayleyHamiltonListIII(A,n):
    """
    Outputs a list of hypermatrices of all product
    composition of degree n for hypermatrices of 
    order up to 8 tends to be slow because it uses
    the general hypermatrix product and it takes HM
    object class as the input A

     EXAMPLES:

    ::

        sage: A = HM(2,2,2,'a')
        sage: L = HypermatrixCayleyHamiltonListIII(A,3)

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
                    gu = gu + [GeneralHypermatrixProduct(g1,g2,g3) for g1 in HypermatrixCayleyHamiltonListIII(A,i) for g2 in HypermatrixCayleyHamiltonListIII(A,j) for g3 in HypermatrixCayleyHamiltonListIII(A,n-(i+j))]
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
                        gu = gu + [GeneralHypermatrixProduct(g1,g2,g3,g4) for g1 in HypermatrixCayleyHamiltonListIII(A,i) for g2 in HypermatrixCayleyHamiltonListIII(A,j) for g3 in HypermatrixCayleyHamiltonListIII(A,k) for g4 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k))]
            return gu
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
                            gu = gu + [GeneralHypermatrixProduct(g1,g2,g3,g4,g5) for g1 in HypermatrixCayleyHamiltonListIII(A,i) for g2 in HypermatrixCayleyHamiltonListIII(A,j) for g3 in HypermatrixCayleyHamiltonListIII(A,k) for g4 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k)) for g5 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k+l))]
            return gu
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
                                gu = gu + [GeneralHypermatrixProduct(g1,g2,g3,g4,g5,g6) for g1 in HypermatrixCayleyHamiltonListIII(A,i) for g2 in HypermatrixCayleyHamiltonListIII(A,j) for g3 in HypermatrixCayleyHamiltonListIII(A,k) for g4 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k)) for g5 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k+l)) for g6 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k+l+m))]
            return gu
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
                                    gu = gu + [GeneralHypermatrixProduct(g1,g2,g3,g4,g5,g6,g7) for g1 in HypermatrixCayleyHamiltonListIII(A,i) for g2 in HypermatrixCayleyHamiltonListIII(A,j) for g3 in HypermatrixCayleyHamiltonListIII(A,k) for g4 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k)) for g5 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k+l)) for g6 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k+l+m)) for g7 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k+l+m+o))]
            return gu
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
                                        gu = gu + [GeneralHypermatrixProduct(g1,g2,g3,g4,g5,g6,g7) for g1 in HypermatrixCayleyHamiltonListIII(A,i) for g2 in HypermatrixCayleyHamiltonListIII(A,j) for g3 in HypermatrixCayleyHamiltonListIII(A,k) for g4 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k)) for g5 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k+l)) for g6 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k+l+m)) for g7 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k+l+m+o))  for g8 in HypermatrixCayleyHamiltonListIII(A,n-(i+j+k+l+m+o+p))]
            return gu
    else :
        raise ValueError("Not supported for order > 4 and for non cube hypermpatrix of order 3 ")
   
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
        raise ValueError("Not supported for order > 4 and for non cube hypermpatrix of order 3 ")

def ConstraintFormator(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs matrix
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    working over CC for both A and b.

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
    working over SR for both A and b.

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

def ConstraintFormatorHM(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs an HM 
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    working over SR for both A and b.

    EXAMPLES:

    ::

        sage: x,y = var('x,y')
        sage: CnstrLst = [x+y==1, x-y==2]
        sage: VrbLst = [x, y]
        sage: [A,b] = ConstraintFormatorHM(CnstrLst, VrbLst)
        sage: A.printHM()
        [:, :]=
        [ 1  1]
        [ 1 -1]
        sage: b.printHM()
        [:, :]=
        [1]
        [2]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initializing the Matrix
    A=HM(len(CnstrLst),len(VrbLst),'zero')
    b=HM(len(CnstrLst), 1, [eq.rhs() for eq in CnstrLst])
    for r in range(len(CnstrLst)):
        for c in range(len(VrbLst)):
            A[r,c]=SR((CnstrLst[r]).lhs().coefficient(VrbLst[c]))
    return [A,b]

def multiplicativeConstraintFormator(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs matrix
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    working over SR for both A and b.

    EXAMPLES:

    ::

        sage: x, y = var('x, y')
        sage: CnstrLst = [x*y^2==1, x/y==2]
        sage: VrbLst = [x, y]
        sage: [A,b] = multiplicativeConstraintFormator(CnstrLst, VrbLst)
        sage: A
        [ 1  2]
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
            A[r,c]=(CnstrLst[r]).lhs().degree(VrbLst[c])
    return [A,b]

def multiplicativeConstraintFormatorII(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs matrix
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    working over SR for both A and b.

    EXAMPLES:

    ::

        sage: x, y = var('x, y')
        sage: CnstrLst = [(1/7)*x*y^2, (1/2)*x/y]
        sage: VrbLst = [x, y]
        sage: [A,b] = multiplicativeConstraintFormatorII(CnstrLst, VrbLst)
        sage: A
        [ 1  2]
        [ 1 -1]
        sage: b
        [7]
        [2]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initialization of the equations
    Eq=[f/f.subs([v==1 for v in VrbLst])==f.subs([v==1 for v in VrbLst])^(-1) for f in CnstrLst]
    return multiplicativeConstraintFormator(Eq, VrbLst)

def multiplicativeConstraintFormatorIII(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs matrix
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    working over SR for both A and b.
    The difference with the implementation above is
    that it handles symbolic exponents.

    EXAMPLES:

    ::

        sage: x, y = var('x, y'); A=HM(2,2,'a')
        sage: CnstrLst = [x^A[0,0]*y^A[0,1]==7, x^A[1,0]*y^A[1,1]==2]
        sage: VrbLst = [x, y]
        sage: [A,b] = multiplicativeConstraintFormatorIII(CnstrLst, VrbLst)
        sage: A
        [a00 a01]
        [a10 a11]
        sage: b
        [7]
        [2]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initializing the Matrix
    A=Matrix(SR,len(CnstrLst),len(VrbLst),zero_matrix(len(CnstrLst),len(VrbLst)))
    b=vector(SR, [eq.rhs() for eq in CnstrLst]).column()
    for r in range(len(CnstrLst)):
        for c in range(len(VrbLst)):
            #A[r,c]=(CnstrLst[r]).lhs().degree(VrbLst[c])
            A[r,c]=((CnstrLst[r]).lhs().diff(VrbLst[c])/(CnstrLst[r]).lhs()).subs(VrbLst[c]==1)
    return [A,b]

def multiplicativeConstraintFormatorIV(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs matrix
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    working over SR for both A and b.
    The difference with the implementation above is
    that it handles symbolic exponents.


    EXAMPLES:

    ::

        sage: x, y = var('x, y'); A=HM(2,2,'a')
        sage: CnstrLst = [(1/7)*x^A[0,0]*y^A[0,1], (1/2)*x^A[1,0]*y^A[1,1]]
        sage: VrbLst = [x, y]
        sage: [A,b] = multiplicativeConstraintFormatorIV(CnstrLst, VrbLst)
        sage: A
        [a00 a01]
        [a10 a11]
        sage: b
        [7]
        [2]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initialization of the equations
    Eq=[f/f.subs([v==1 for v in VrbLst])==f.subs([v==1 for v in VrbLst])^(-1) for f in CnstrLst]
    return multiplicativeConstraintFormatorIII(Eq, VrbLst)

def multiplicativeConstraintFormatorHM(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs HM
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    working over SR for both A and b.


    EXAMPLES:

    ::

        sage: x, y = var('x, y')
        sage: CnstrLst = [x*y^2==1, x/y==2]
        sage: VrbLst = [x, y]
        sage: [Ha,hb] = multiplicativeConstraintFormatorHM(CnstrLst, VrbLst)
        sage: Ha.printHM()
        [:, :]=
        [ 1  2]
        [ 1 -1]
        sage: hb.printHM()
        [:, :]=
        [1]
        [2]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initializing the Matrix
    A=HM(len(CnstrLst),len(VrbLst),'zero')
    b=HM(len(CnstrLst), 1, [eq.rhs() for eq in CnstrLst])
    for r in range(len(CnstrLst)):
        for c in range(len(VrbLst)):
            A[r,c]=SR((CnstrLst[r]).lhs().degree(VrbLst[c]))
    return [A,b]

def multiplicativeConstraintFormatorIIHM(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs HM
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    working over SR for both A and b.

    EXAMPLES:

    ::

        sage: x, y = var('x, y')
        sage: CnstrLst = [(1/7)*x*y^2, (1/2)*x/y]
        sage: VrbLst = [x, y]
        sage: [Ha,hb] = multiplicativeConstraintFormatorIIHM(CnstrLst, VrbLst)
        sage: Ha.printHM()
        [:, :]=
        [ 1  2]
        [ 1 -1]
        sage: hb.printHM()
        [:, :]=
        [7]
        [2]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    Eq=[f/f.subs([v==1 for v in VrbLst])==f.subs([v==1 for v in VrbLst])^(-1) for f in CnstrLst]
    return multiplicativeConstraintFormatorHM(Eq, VrbLst)

def multiplicativeConstraintFormatorIIIHM(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs HM
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    working over SR for both A and b.
    The difference with the implementation above is
    that it handles symbolic exponents.


    EXAMPLES:

    ::

        sage: x, y = var('x, y'); A=HM(2,2,'a')
        sage: CnstrLst = [x^A[0,0]*y^A[0,1]==7, x^A[1,0]*y^A[1,1]==2]
        sage: VrbLst = [x, y]
        sage: [Ha,hb] = multiplicativeConstraintFormatorIIIHM(CnstrLst, VrbLst)
        sage: Ha.printHM()
        [:, :]=
        [a00 a01]
        [a10 a11]
        sage: hb.printHM()
        [:, :]=
        [7]
        [2]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initializing the Matrix
    A=HM(len(CnstrLst),len(VrbLst),'zero')
    b=HM(len(CnstrLst), 1, [eq.rhs() for eq in CnstrLst])
    for r in range(len(CnstrLst)):
        for c in range(len(VrbLst)):
            #A[r,c]=SR((CnstrLst[r]).lhs().degree(VrbLst[c]))
            A[r,c]=((CnstrLst[r]).lhs().diff(VrbLst[c])/(CnstrLst[r]).lhs()).subs(VrbLst[c]==1)
    return [A,b]

def multiplicativeConstraintFormatorIVHM(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs matrix
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    working over SR for both A and b.
    The difference with the implementation above is
    that it handles symbolic exponents.


    EXAMPLES:

    ::

        sage: x, y = var('x, y'); A=HM(2,2,'a')
        sage: CnstrLst = [(1/7)*x^A[0,0]*y^A[0,1], (1/2)*x^A[1,0]*y^A[1,1]]
        sage: VrbLst = [x, y]
        sage: [A,b] = multiplicativeConstraintFormatorIVHM(CnstrLst, VrbLst)
        sage: A.printHM()
        [:, :]=
        [a00 a01]
        [a10 a11]
        sage: b.printHM()
        [:, :]=
        [7]
        [2]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initialization of the equations
    Eq=[f/f.subs([v==1 for v in VrbLst])==f.subs([v==1 for v in VrbLst])^(-1) for f in CnstrLst]
    return multiplicativeConstraintFormatorIIIHM(Eq, VrbLst)

def ConstraintFormatorIII(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs matrix
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    working over CC for A and SR for b

    EXAMPLES:

    ::

        sage: x,y=var('x,y'); CnstrLst=[x+y==1, x-y==2]; VrbLst=[x,y]
        sage: [A,b] = ConstraintFormatorIII(CnstrLst, VrbLst)
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

def ConstraintFormatorIV(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs matrix
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    this implementation allows for the lefthand side
    not to be specified but input variables must not
    be monomials.

    EXAMPLES:

    ::

        sage: x,y = var('x,y')
        sage: CnstrLst = [x+y-1, x-y-2]
        sage: VrbLst = [x, y]
        sage: [A,b] = ConstraintFormatorIV(CnstrLst, VrbLst)
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
    for r in range(len(CnstrLst)):
        for c in range(len(VrbLst)):
            Tmp=Set(VrbLst).difference(Set([VrbLst[c]])).list()
            A[r,c]=CnstrLst[r].subs([f==0 for f in Tmp]).coefficient(VrbLst[c])
    b=-Matrix(len(CnstrLst),1,CnstrLst).subs([f==0 for f in VrbLst])
    return [A,b]

def ConstraintFormatorIVHM(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs matrix
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    this implementation allows for the lefthand side
    not to be specified but input variables must not
    be monomials.

    EXAMPLES:

    ::

        sage: x,y = var('x,y')
        sage: CnstrLst = [x+y-1, x-y-2]
        sage: VrbLst = [x, y]
        sage: [A,b] = ConstraintFormatorIVHM(CnstrLst, VrbLst)
        sage: A.printHM()
        [:, :]=
        [ 1  1]
        [ 1 -1]
        sage: b.printHM()
        [:, :]=
        [1]
        [2]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initializing the Matrix
    A=HM(len(CnstrLst),len(VrbLst),'zero')
    for r in range(len(CnstrLst)):
        for c in range(len(VrbLst)):
            Tmp=Set(VrbLst).difference(Set([VrbLst[c]])).list()
            A[r,c]=CnstrLst[r].subs([f==0 for f in Tmp]).coefficient(VrbLst[c])
    b=-HM(len(CnstrLst),1,CnstrLst).subs([f==0 for f in VrbLst])
    return [A,b]

def MonomialConstraintFormator(L, X, MnL, Y):
    """
    Takes as input a List of polynomials, a list of
    variables used in the polynomials,the monomial list 
    a list of alternative variables to replace the monomials.
    No right hand side is given. We are implicitly working over SR.


    EXAMPLES:

    ::


        sage: x1, x2 = var('x1, x2')
        sage: X = [x1, x2]
        sage: MnL=[x1*x2, x1, x2]
        sage: L = [-2*x1*x2 + 3*x2 - 2, -38*x1*x2 - 18*x1 + 11*x2 - 8, -506*x1*x2 + 112*x1 + 121*x2 - 8, -5852*x1*x2 + 202*x1 + 1099*x2 - 8]
        sage: Y = HM(3,'y').list()
        sage: [A,b] = MonomialConstraintFormator(L, X, MnL, Y)
        sage: A
        [   -2     0     3]
        [  -38   -18    11]
        [ -506   112   121]
        [-5852   202  1099]
        sage: b
        [2]
        [8]
        [8]
        [8]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Obtaining the right hand side vector b
    tb=Matrix(SR,len(L),1,[-f.subs([v==0 for v in X]) for f in L])
    L2=copy(L)
    # Updating the list to remove the constant terms
    for i in range(len(L)):
        L2[i]=L2[i]+tb[i,0]
    # Performing the monomial substitution
    Eq=[]; Hy=HM([MnL,Y])
    cnt=0
    for g in L2:
        TmpL=[]
        for o in g.operands():
            for j in range(Hy.n(1)):
                if (o/Hy[0,j]).is_constant():
                    TmpL.append((o/Hy[0,j])*Hy[1,j]);break
        Eq.append(sum(TmpL)==tb[cnt,0]); cnt=cnt+1 
    # Ready to use the generic Constraint formator
    return ConstraintFormatorII(Eq, Y)

def MonomialConstraintFormatorII(L, X, MnL, Y):
    """
    Takes as input a List of polynomials, a list of
    variables used in the polynomials,the monomial list 
    a list of alternative variables to replace the monomials.
    No right hand side is given. We are implicitly working over SR.
    The difference with the implementation above is that it not assume
    that the coefficients are constants they can be themselves polynomials.


    EXAMPLES:

    ::


        sage: x1, x2 = var('x1, x2')
        sage: X = [x1, x2]
        sage: MnL=[x1*x2, x1, x2]
        sage: L = [-2*x1*x2 + 3*x2 - 2, -38*x1*x2 - 18*x1 + 11*x2 - 8, -506*x1*x2 + 112*x1 + 121*x2 - 8, -5852*x1*x2 + 202*x1 + 1099*x2 - 8]
        sage: Y = HM(3,'y').list()
        sage: [A,b] = MonomialConstraintFormatorII(L, X, MnL, Y)
        sage: A
        [   -2     0     3]
        [  -38   -18    11]
        [ -506   112   121]
        [-5852   202  1099]
        sage: b
        [2]
        [8]
        [8]
        [8]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Obtaining the right hand side vector b
    tb=Matrix(SR,len(L),1,[-f.subs([v==0 for v in X]) for f in L])
    L2=copy(L)
    # Updating the list to remove the constant terms
    for i in range(len(L)):
        L2[i]=L2[i]+tb[i,0]
    # Performing the monomial substitution
    Eq=[]; Hy=HM([MnL,Y])
    cnt=0
    for g in L2:
        TmpL=[]
        for o in g.operands():
            for j in range(Hy.n(1)):
                #if (o/Hy[0,j]).is_constant():
                if Set([(o/Hy[0,j]).degree(v) for v in X]).list()==[0]:
                    TmpL.append((o/Hy[0,j])*Hy[1,j]);break
        Eq.append(sum(TmpL)==tb[cnt,0]); cnt=cnt+1 
    # Ready to use the generic Constraint formator
    return ConstraintFormatorII(Eq, Y)

def MonomialConstraintFormatorIII(L, X, MnL, Y):
    """
    Takes as input a List of polynomials, a list of
    variables used in the polynomials,the monomial list 
    a list of alternative variables to replace the monomials.
    No right hand side is given. We are implicitly working over SR.
    The difference with the implementation above is
    that it handles symbolic exponents.


    EXAMPLES:

    ::


        sage: x1, x2, y1, y2, y3, y4 = var('x1, x2, y1, y2, y3, y4')
        sage: X = [x1, x2]
        sage: MnL=[x1^y1*x2^y2, x1^y3, x2^y4]
        sage: L = [-2*x1^y1*x2^y2 + 3*x2^y4 - 2, -38*x1^y1*x2^y2 - 18*x1^y3 + 11*x2^y4 - 8, -506*x1^y1*x2^y2 + 112*x1^y3 + 121*x2^y4 - 8, -5852*x1^y1*x2^y2 + 202*x1^y3 + 1099*x2^y4 - 8]
        sage: Y = HM(3,'y').list()
        sage: [A,b] = MonomialConstraintFormatorIII(L, X, MnL, Y)
        sage: A
        [   -2     0     3]
        [  -38   -18    11]
        [ -506   112   121]
        [-5852   202  1099]
        sage: b
        [                 2*0^y1*0^y2 - 3*0^y4 + 2]
        [     38*0^y1*0^y2 + 18*0^y3 - 11*0^y4 + 8]
        [  506*0^y1*0^y2 - 112*0^y3 - 121*0^y4 + 8]
        [5852*0^y1*0^y2 - 202*0^y3 - 1099*0^y4 + 8]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Obtaining the right hand side vector b
    tb=Matrix(SR,len(L),1,[-f.subs([v==0 for v in X]) for f in L])
    L2=copy(L)
    # Updating the list to remove the constant terms
    for i in range(len(L)):
        L2[i]=L2[i]+tb[i,0]
    # Performing the monomial substitution
    Eq=[]; Hy=HM([MnL,Y])
    cnt=0
    for g in L2:
        TmpL=[]
        for o in g.operands():
            for j in range(Hy.n(1)):
                #if (o/Hy[0,j]).is_constant():
                #if Set([(o/Hy[0,j]).degree(v) for v in X]).list()==[0]:
                if Set([(o/Hy[0,j]).diff(v)/((o/Hy[0,j]).subs(v==1)) for v in X]).list()==[0]:
                    TmpL.append((o/Hy[0,j])*Hy[1,j]);break
        Eq.append(sum(TmpL)==tb[cnt,0]); cnt=cnt+1 
    # Ready to use the generic Constraint formator
    return ConstraintFormatorII(Eq, Y)

def MonomialConstraintFormatorHM(L, X, MnL, Y):
    """
    Takes as input a List of polynomials, a list of
    variables used in the polynomials,the monomial list 
    a list of alternative variables to replace the monomials.
    No right hand side is given. We are implicitly working over SR.


    EXAMPLES:

    ::


        sage: x1, x2 = var('x1, x2')
        sage: X = [x1, x2]
        sage: MnL=[x1*x2, x1, x2]
        sage: L = [-2*x1*x2 + 3*x2 - 2, -38*x1*x2 - 18*x1 + 11*x2 - 8, -506*x1*x2 + 112*x1 + 121*x2 - 8, -5852*x1*x2 + 202*x1 + 1099*x2 - 8]
        sage: Y = HM(3,'y').list()
        sage: [A,b] = MonomialConstraintFormatorHM(L, X, MnL, Y)
        sage: A.printHM()
        [:, :]=
        [   -2     0     3]
        [  -38   -18    11]
        [ -506   112   121]
        [-5852   202  1099]
        sage: b.printHM()
        [:, :]=
        [2]
        [8]
        [8]
        [8]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Obtaining the right hand side vector b
    tb=Matrix(SR,len(L),1,[-f.subs([v==0 for v in X]) for f in L])
    L2=copy(L)
    # Updating the list to remove the constant terms
    for i in range(len(L)):
        L2[i]=L2[i]+tb[i,0]
    # Performing the monomial substitution
    Eq=[]; Hy=HM([MnL,Y])
    cnt=0
    for g in L2:
        TmpL=[]
        for o in g.operands():
            for j in range(Hy.n(1)):
                if (o/Hy[0,j]).is_constant():
                    TmpL.append((o/Hy[0,j])*Hy[1,j]);break
        Eq.append(sum(TmpL)==tb[cnt,0]); cnt=cnt+1 
    # Ready to use the generic Constraint formator
    return ConstraintFormatorHM(Eq, Y)

def MonomialConstraintFormatorHMII(L, X, MnL, Y):
    """
    Takes as input a List of polynomials, a list of
    variables used in the polynomials,the monomial list 
    a list of alternative variables to replace the monomials.
    No right hand side is given. We are implicitly working over SR.
    The difference with the implementation above is that this
    function deals well with symbolic coefficient matrices

    EXAMPLES:

    ::


        sage: x1, x2 = var('x1, x2')
        sage: X = [x1, x2]
        sage: MnL=[x1*x2, x1, x2]
        sage: L = [-2*x1*x2 + 3*x2 - 2, -38*x1*x2 - 18*x1 + 11*x2 - 8, -506*x1*x2 + 112*x1 + 121*x2 - 8, -5852*x1*x2 + 202*x1 + 1099*x2 - 8]
        sage: Y = HM(3,'y').list()
        sage: [A,b] = MonomialConstraintFormatorHMII(L, X, MnL, Y)
        sage: A.printHM()
        [:, :]=
        [   -2     0     3]
        [  -38   -18    11]
        [ -506   112   121]
        [-5852   202  1099]
        sage: b.printHM()
        [:, :]=
        [2]
        [8]
        [8]
        [8]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Obtaining the right hand side vector b
    tb=Matrix(SR,len(L),1,[-f.subs([v==0 for v in X]) for f in L])
    L2=copy(L)
    # Updating the list to remove the constant terms
    for i in range(len(L)):
        L2[i]=L2[i]+tb[i,0]
    # Performing the monomial substitution
    Eq=[]; Hy=HM([MnL,Y])
    cnt=0
    for g in L2:
        TmpL=[]
        for o in g.operands():
            for j in range(Hy.n(1)):
                #if (o/Hy[0,j]).is_constant():
                if Set([(o/Hy[0,j]).degree(v) for v in X]).list()==[0]:
                    TmpL.append((o/Hy[0,j])*Hy[1,j]);break
        Eq.append(sum(TmpL)==tb[cnt,0]); cnt=cnt+1 
    # Ready to use the generic Constraint formator
    return ConstraintFormatorHM(Eq, Y)

def MonomialConstraintFormatorHMIII(L, X, MnL, Y):
    """
    Takes as input a List of polynomials, a list of
    variables used in the polynomials,the monomial list 
    a list of alternative variables to replace the monomials.
    No right hand side is given. We are implicitly working over SR.
    The difference with the implementation above is that this
    function deals well with symbolic coefficient matrices

    EXAMPLES:

    ::


        sage: x1, x2, y1, y2, y3, y4 = var('x1, x2, y1, y2, y3, y4')
        sage: X = [x1, x2]
        sage: MnL=[x1^y1*x2^y2, x1^y3, x2^y4]
        sage: L = [-2*x1^y1*x2^y2 + 3*x2^y4 - 2, -38*x1^y1*x2^y2 - 18*x1^y3 + 11*x2^y4 - 8, -506*x1^y1*x2^y2 + 112*x1^y3 + 121*x2^y4 - 8, -5852*x1^y1*x2^y2 + 202*x1^y3 + 1099*x2^y4 - 8]
        sage: Y = HM(3,'y').list()
        sage: [A,b] = MonomialConstraintFormatorHMIII(L, X, MnL, Y)
        sage: A.printHM()
        [:, :]=
        [   -2     0     3]
        [  -38   -18    11]
        [ -506   112   121]
        [-5852   202  1099]
        sage: b.printHM()
        [:, :]=
        [                 2*0^y1*0^y2 - 3*0^y4 + 2]
        [     38*0^y1*0^y2 + 18*0^y3 - 11*0^y4 + 8]
        [  506*0^y1*0^y2 - 112*0^y3 - 121*0^y4 + 8]
        [5852*0^y1*0^y2 - 202*0^y3 - 1099*0^y4 + 8]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Obtaining the right hand side vector b
    tb=Matrix(SR,len(L),1,[-f.subs([v==0 for v in X]) for f in L])
    L2=copy(L)
    # Updating the list to remove the constant terms
    for i in range(len(L)):
        L2[i]=L2[i]+tb[i,0]
    # Performing the monomial substitution
    Eq=[]; Hy=HM([MnL,Y])
    cnt=0
    for g in L2:
        TmpL=[]
        for o in g.operands():
            for j in range(Hy.n(1)):
                #if (o/Hy[0,j]).is_constant():
                #if Set([(o/Hy[0,j]).degree(v) for v in X]).list()==[0]:
                if Set([(o/Hy[0,j]).diff(v)/((o/Hy[0,j]).subs(v==1)) for v in X]).list()==[0]:
                    TmpL.append((o/Hy[0,j])*Hy[1,j]);break
        Eq.append(sum(TmpL)==tb[cnt,0]); cnt=cnt+1 
    # Ready to use the generic Constraint formator
    return ConstraintFormatorHM(Eq, Y)

def exponentialConstraintFormatorHM(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs HM
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    working over SR for both A and b.
    The difference with the implementation above is
    that it handles symbolic exponents.


    EXAMPLES:

    ::

        sage: x, y = var('x, y'); A=HM(2,2,'a')
        sage: CnstrLst = [(A[0,0]^x)*(A[0,1]^y) == 7, (A[1,0]^x)*(A[1,1])^y == 2]
        sage: VrbLst = [x, y]
        sage: [Ha,hb] = exponentialConstraintFormatorHM(CnstrLst, VrbLst)
        sage: Ha.printHM()
        [:, :]=
        [a00 a01]
        [a10 a11]
        sage: hb.printHM()
        [:, :]=
        [7]
        [2]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initializing the Matrix
    A=HM(len(CnstrLst),len(VrbLst),'zero')
    b=HM(len(CnstrLst), 1, [eq.rhs() for eq in CnstrLst])
    for r in range(len(CnstrLst)):
        for c in range(len(VrbLst)):
            A[r,c]=((CnstrLst[r]).lhs().diff(VrbLst[c])/(CnstrLst[r]).lhs()).subs(VrbLst[c]==1)
    return [A.elementwise_base_exponent(e),b]

def exponentialConstraintFormatorHMII(CnstrLst, VrbLst):
    """
    Takes as input a List of linear constraints
    and a list of variables and outputs matrix
    and the right hand side vector associate
    with the matrix formulation of the constraints.
    working over SR for both A and b.
    The difference with the implementation above is
    that it handles symbolic exponents.


    EXAMPLES:

    ::

        sage: x, y = var('x, y'); A=HM(2,2,'a')
        sage: CnstrLst = [(1/7)*(A[0,0]^x)*(A[0,1]^y), (1/2)*(A[1,0]^x)*(A[1,1])^y]
        sage: VrbLst = [x, y]
        sage: [A,b] = exponentialConstraintFormatorHMII(CnstrLst, VrbLst)
        sage: A.printHM()
        [:, :]=
        [a00 a01]
        [a10 a11]
        sage: b.printHM()
        [:, :]=
        [7]
        [2]


    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    # Initialization of the equations
    Eq=[f/f.subs([v==0 for v in VrbLst])==f.subs([v==0 for v in VrbLst])^(-1) for f in CnstrLst]
    return exponentialConstraintFormatorHM(Eq, VrbLst)

def Companion_matrix(p,vrbl):
    """
    Takes as input a polynomial and a variable
    and outputs the companion matrix associated
    with the polynomial in the specified variables.

    EXAMPLES:

    ::

        sage: x=var('x'); p=sum(HM(5,'a').list()[k]*x^(k) for k in range(5))
        sage: A=Companion_matrix(p,x); A.characteristic_polynomial()
        x^4 + a3/a4*x^3 + a2/a4*x^2 + a1/a4*x + a0/a4

    AUTHORS:
    - Edinah K. Gnang
    """
    if p.is_polynomial(vrbl):
        dg=Integer(p.degree(vrbl))
        if dg>1:
            # Initialization of the matrix
            A=HM(dg,dg,'zero').matrix()
            # Filling up the matrix
            A[0,dg-1]=-p.subs(dict([(vrbl,0)]))/(p.diff(vrbl, dg)/factorial(dg))
            for i in range(1,dg):
                #A[i,dg-1]=-p.coefficient(vrbl^(i))/p.coefficient(vrbl^dg)
                A[i,dg-1]=-(p.diff(vrbl,i).subs(vrbl==0)/factorial(i))/(p.diff(vrbl,dg).subs(vrbl==0)/factorial(dg))
                A[i,i-1]=1
            return A
        elif dg==1:
            #return Matrix(SR,1,1,[p.subs(dict([(vrbl,0)]))/p.coefficient(vrbl)])
            return Matrix(SR,1,1,[p.subs(dict([(vrbl,0)]))/(diff(p,vrbl).subs(vrbl==0))])
        else:
            raise ValueError("Must be of degree at least 1.")
    else:
        raise ValueError("Must be a polynomial in the input variable.")

def CompanionHM(p,vrbl):
    """
    Takes as input a polynomial and a variable
    and outputs the companion second order hypermatrix
    associated with the polynomial in the specified variables.

    EXAMPLES:

    ::

        sage: sz=5; x=var('x'); Id=HM(2,sz-1,'kronecker')
        sage: A=CompanionHM(sum(HM(sz,'a').list()[k]*x^(k) for k in range(sz)),x); expand((x*Id-A).det())
        x^4 + a3*x^3/a4 + a2*x^2/a4 + a1*x/a4 + a0/a4

    AUTHORS:
    - Edinah K. Gnang
    """
    if p.is_polynomial(vrbl):
        dg=Integer(p.degree(vrbl))
        if dg>1:
            # Initialization of the matrix
            A=HM(dg,dg,'zero')
            # Filling up the matrix
            #A[0,dg-1]=-p.subs(dict([(vrbl,0)]))/p.coefficient(vrbl^dg)
            A[0,dg-1]=-p.subs(dict([(vrbl,0)]))/(p.diff(vrbl, dg)/factorial(dg))
            for i in range(1,dg):
                #A[i,dg-1]=-p.coefficient(vrbl^(i))/p.coefficient(vrbl^dg)
                A[i,dg-1]=-(p.diff(vrbl,i).subs(vrbl==0)/factorial(i))/(p.diff(vrbl,dg).subs(vrbl==0)/factorial(dg))
                A[i,i-1]=1
            return A
        elif dg==1:
            #return HM(1,1,[p.subs(dict([(vrbl,0)]))/p.coefficient(vrbl)])
            return HM(1,1,[p.subs(dict([(vrbl,0)]))/(diff(p,vrbl).subs(vrbl==0))])
        else:
            raise ValueError("Must be of degree at least 1.")
    else:
        raise ValueError("Must be a polynomial in the input variable.")

def Sylvester_matrix(p,q,vrbl):
    """
    Takes as input two polynomials and a variable
    and outputs the Sylvester matrix associated
    with the polynomials in the specified variables.

    EXAMPLES:

    ::

        sage: x, a0, a1, b0, b1=var('x, a0, a1, b0, b1')
        sage: p=(x-a0)*(x-a1); q=(x-b0)*(x-b1)
        sage: Sylvester_matrix(p, q, x).det().factor()
        (a0 - b0)*(a0 - b1)*(a1 - b0)*(a1 - b1)

    AUTHORS:
    - Edinah K. Gnang
    """
    if p.is_polynomial(vrbl) and q.is_polynomial(vrbl):
        dp=Integer(p.degree(vrbl)); dq=Integer(q.degree(vrbl))
        # Initialization of the matrix
        A=HM(dp+dq,dp+dq,'zero').matrix()
        # Filling up the matrix
        cp=0
        for i in range(dq):
            for j in range(dp):
                #A[i,cp+j]=p.coefficient(vrbl^(dp-j))
                A[i,cp+j]=p.diff(vrbl,Integer(dp-j)).subs(vrbl==0)/factorial(dp-j)
            A[i,cp+dp]=p.subs(dict([(vrbl,0)]))
            cp=cp+1
        cq=0
        for i in range(dp):
            for j in range(dq):
                #A[dq+i,cq+j]=q.coefficient(vrbl^(dq-j))
                A[dq+i,cq+j]=q.diff(vrbl,Integer(dq-j)).subs(vrbl==0)/factorial(dq-j)
            A[dq+i,cq+dq]=q.subs(dict([(vrbl,0)]))
            cq=cq+1
        return A
    else:
        raise ValueError("The inputs must both be polynomials in the input variable.")

def SylvesterHM(p,q,vrbl):
    """
    Takes as input two polynomials and a variable
    and outputs the Sylvester second order 
    hypermatrix associated with the polynomials
    in the specified variables.

    EXAMPLES:

    ::

        sage: x, a0, a1, b0, b1=var('x, a0, a1, b0, b1')
        sage: p=(x-a0)*(x-a1); q=(x-b0)*(x-b1)
        sage: SylvesterHM(p, q, x).ref()[3,3].factor()
        (a0 - b0)*(a0 - b1)*(a1 - b0)*(a1 - b1)

    AUTHORS:
    - Edinah K. Gnang
    """
    if p.is_polynomial(vrbl) and q.is_polynomial(vrbl):
        dp=Integer(p.degree(vrbl)); dq=Integer(q.degree(vrbl))
        # Initialization of the second order hypermatrix
        A=HM(dp+dq,dp+dq,'zero')
        # Filling up the matrix
        cp=0
        for i in range(dq):
            for j in range(dp):
                #A[i,cp+j]=p.coefficient(vrbl^(dp-j))
                A[i,cp+j]=p.diff(vrbl,Integer(dp-j)).subs(vrbl==0)/factorial(dp-j)
            A[i,cp+dp]=p.subs(dict([(vrbl,0)]))
            cp=cp+1
        cq=0
        for i in range(dp):
            for j in range(dq):
                #A[dq+i,cq+j]=q.coefficient(vrbl^(dq-j))
                A[dq+i,cq+j]=q.diff(vrbl,Integer(dq-j)).subs(vrbl==0)/factorial(dq-j)
            A[dq+i,cq+dq]=q.subs(dict([(vrbl,0)]))
            cq=cq+1
        return A
    else:
        raise ValueError("The inputs must both be polynomials in the input variable.")

def Resultant(p, q, vrbl):
    """
    Takes as input two polynomials and a variable
    and outputs the corresponding resultant.

    EXAMPLES:

    ::

        sage: x, a0, a1, b0, b1=var('x, a0, a1, b0, b1')
        sage: p=(x-a0)*(x-a1); q=(x-b0)*(x-b1)
        sage: Resultant(p, q, x).factor()
        (a0 - b0)*(a0 - b1)*(a1 - b0)*(a1 - b1)

    AUTHORS:
    - Edinah K. Gnang
    """
    return SylvesterHM(p, q, vrbl).det()

def Gmatrix(p,q,vrbl):
    """
    Takes as input two polynomials and a variable
    and outputs the G matrix associated with the 
    polynomial in the specified variables.

    EXAMPLES:

    ::

        sage: x, a0, a1, b0, b1=var('x, a0, a1, b0, b1')
        sage: p=(x-a0)*(x-a1); q=(x-b0)*(x-b1)
        sage: Gmatrix(p,q,x).det().factor()
        (a0 - b0)*(a0 - b1)*(a1 - b0)*(a1 - b1)

    AUTHORS:
    - Edinah K. Gnang
    """
    if p.is_polynomial(vrbl) and q.is_polynomial(vrbl):
        dp=Integer(p.degree(vrbl)); dq=Integer(q.degree(vrbl))
        if dp >= 1 and dq >= 1:
            return identity_matrix(dq).tensor_product(Companion_matrix(p,vrbl))-(Companion_matrix(q,vrbl)).tensor_product(identity_matrix(dp))
        else:
            raise ValueError("Both inputs must be of degree at least 2.")
    else:
        raise ValueError("Both inputs must be polynomials in the input variable.")

def GmatrixHM(p,q,vrbl):
    """
    Takes as input two polynomials and a variable
    and outputs the G matrix associated with the 
    polynomial in the specified variables.

    EXAMPLES:

    ::

        sage: x, a0, a1, b0, b1=var('x, a0, a1, b0, b1')
        sage: p=(x-a0)*(x-a1); q=(x-b0)*(x-b1)
        sage: GmatrixHM(p,q,x).ref()[3,3].factor()
        -(a0 - b0)*(a0 - b1)*(a1 - b0)*(a1 - b1)

    AUTHORS:
    - Edinah K. Gnang
    """
    if p.is_polynomial(vrbl) and q.is_polynomial(vrbl):
        dp=Integer(p.degree(vrbl)); dq=Integer(q.degree(vrbl))
        if dp >= 1 and dq >= 1:
            return HM(2,dq,'kronecker').tensor_product(CompanionHM(p,vrbl))-(CompanionHM(q,vrbl)).tensor_product(HM(2,dp,'kronecker'))
        else:
            raise ValueError("Both inputs must be of degree at least 2.")
    else:
        raise ValueError("Both inputs must be polynomials in the input variable.")

def substitute_matrix(p, vrbl, A):
    """
    The functions takes as input a polynomial p,
    a variable vrbl, and a matrix A. The function
    outputs the polynomial in the variable.

    EXAMPLES:

    ::

        sage: x,y = var('x,y')
        sage: p=x^2+2*x*y+1
        sage: substitute_matrix(p,x,HM(2,2,'a').matrix())
        [a00^2 + a01*a10 + 2*a00*y + 1   a00*a01 + a01*a11 + 2*a01*y]
        [  a00*a10 + a10*a11 + 2*a10*y a01*a10 + a11^2 + 2*a11*y + 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    if A.nrows() == A.ncols():
        T=p.subs(vrbl == 0)*identity_matrix(A.nrows())
        d=Integer(p.degree(vrbl))
        for i in rg(1,d+1):
            #T=T+(A^i)*p.coefficient(vrbl^i)
            #T=T+(A^Integer(i))*(p.diff(vrbl,i).subs(vrbl==0)/factorial(i))
            T=T+(A^i)*(p.diff(vrbl,i).subs(vrbl==0)/factorial(i))
        return T
    else:
        raise ValueError("Must be a polynomial in the input variable.")

def substituteHM(p, vrbl, A):
    """
    The functions takes as input a polynomial p,
    a variable vrbl, and a hypermatrix A of order 2.
    The function outputs the polynomial in the variable.


    EXAMPLES:

    ::

        sage: x = var('x'); A=HM(2,2,'a')
        sage: p=x^2 - A.trace()*x + A.det()
        sage: substituteHM(p,x,A).expand()
        [[0, 0], [0, 0]]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    if A.nrows()==A.ncols():
        d=Integer(p.degree(vrbl))
        T = p.subs(vrbl==0)*HM(2,A.nrows(),'kronecker')
        for i in rg(1,d+1):
            #T=T+(A^Integer(i))*(p.diff(vrbl,i).subs(vrbl==0)/factorial(i))
            T=T+(A^i)*(p.diff(vrbl,i).subs(vrbl==0)/factorial(i))
        return T
    else:
        raise ValueError("Must be a polynomial in the input variable.")

def OuterHypermatrixInversePair(U, V):
    """
    Outputs the pseudo inverse pairs associated with the input pairs of hypermatrices
    The implementation does not assume that the input third order hypermatrices are
    cubic. In fact U must be m x p x p and V is p x n x p.


    EXAMPLES:

    ::

        sage: Hu=HM(3,2,2,'u'); Hv=HM(2,4,2,'v')
        sage: [Sln, Tx, Ty]=OuterHypermatrixInversePair(Hu, Hv)[0]
        sage: Hx=Tx.subs(dict([(s.lhs(),s.rhs()) for s in Sln if s.lhs()!=1])) 
        sage: Hy=Ty.subs(dict([(s.lhs(),s.rhs()) for s in Sln if s.lhs()!=1]))
        sage: Prod(Hx, Prod(Hu, HM(3,4,2,'a'), Hv), Hy).simplify_full().list()
        [a000,
         a100,
         a200,
         ((a010*u001*u010*u201*u210*v001*v011*v100 - a010*u001*u010*u200*u211*v000*v011*v101)*v110 - (a010*u000*u011*u201*u210*v001*v010*v100 - a010*u000*u011*u200*u211*v000*v010*v101)*v111)/((u001*u010*u201*u210*v001*v011*v100 - u000*u011*u201*u210*v000*v011*v101)*v110 - (u001*u010*u200*u211*v001*v010*v100 - u000*u011*u200*u211*v000*v010*v101)*v111),
         ((a110*u101*u110*u201*u210*v001*v011*v100 - a110*u101*u110*u200*u211*v000*v011*v101)*v110 - (a110*u100*u111*u201*u210*v001*v010*v100 - a110*u100*u111*u200*u211*v000*v010*v101)*v111)/((u101*u110*u201*u210*v001*v011*v100 - u100*u111*u201*u210*v000*v011*v101)*v110 - (u101*u110*u200*u211*v001*v010*v100 - u100*u111*u200*u211*v000*v010*v101)*v111),
         a210,
         ((a020*u001*u010*u201*u210*v001*v021*v100 - a020*u001*u010*u200*u211*v000*v021*v101)*v120 - (a020*u000*u011*u201*u210*v001*v020*v100 - a020*u000*u011*u200*u211*v000*v020*v101)*v121)/((u001*u010*u201*u210*v001*v021*v100 - u000*u011*u201*u210*v000*v021*v101)*v120 - (u001*u010*u200*u211*v001*v020*v100 - u000*u011*u200*u211*v000*v020*v101)*v121),
         ((a120*u101*u110*u201*u210*v001*v021*v100 - a120*u101*u110*u200*u211*v000*v021*v101)*v120 - (a120*u100*u111*u201*u210*v001*v020*v100 - a120*u100*u111*u200*u211*v000*v020*v101)*v121)/((u101*u110*u201*u210*v001*v021*v100 - u100*u111*u201*u210*v000*v021*v101)*v120 - (u101*u110*u200*u211*v001*v020*v100 - u100*u111*u200*u211*v000*v020*v101)*v121),
         a220,
         ((a030*u001*u010*u201*u210*v001*v031*v100 - a030*u001*u010*u200*u211*v000*v031*v101)*v130 - (a030*u000*u011*u201*u210*v001*v030*v100 - a030*u000*u011*u200*u211*v000*v030*v101)*v131)/((u001*u010*u201*u210*v001*v031*v100 - u000*u011*u201*u210*v000*v031*v101)*v130 - (u001*u010*u200*u211*v001*v030*v100 - u000*u011*u200*u211*v000*v030*v101)*v131),
         ((a130*u101*u110*u201*u210*v001*v031*v100 - a130*u101*u110*u200*u211*v000*v031*v101)*v130 - (a130*u100*u111*u201*u210*v001*v030*v100 - a130*u100*u111*u200*u211*v000*v030*v101)*v131)/((u101*u110*u201*u210*v001*v031*v100 - u100*u111*u201*u210*v000*v031*v101)*v130 - (u101*u110*u200*u211*v001*v030*v100 - u100*u111*u200*u211*v000*v030*v101)*v131),
         a230,
         a001,
         a101,
         a201,
         a011,
         ((a111*u001*u010*u101*u110*v001*v011*v100 - a111*u000*u011*u101*u110*v000*v011*v101)*v110 - (a111*u001*u010*u100*u111*v001*v010*v100 - a111*u000*u011*u100*u111*v000*v010*v101)*v111)/((u001*u010*u101*u110*v001*v011*v100 - u001*u010*u100*u111*v000*v011*v101)*v110 - (u000*u011*u101*u110*v001*v010*v100 - u000*u011*u100*u111*v000*v010*v101)*v111),
         ((a211*u001*u010*u201*u210*v001*v011*v100 - a211*u000*u011*u201*u210*v000*v011*v101)*v110 - (a211*u001*u010*u200*u211*v001*v010*v100 - a211*u000*u011*u200*u211*v000*v010*v101)*v111)/((u001*u010*u201*u210*v001*v011*v100 - u001*u010*u200*u211*v000*v011*v101)*v110 - (u000*u011*u201*u210*v001*v010*v100 - u000*u011*u200*u211*v000*v010*v101)*v111),
         a021,
         ((a121*u001*u010*u101*u110*v001*v021*v100 - a121*u000*u011*u101*u110*v000*v021*v101)*v120 - (a121*u001*u010*u100*u111*v001*v020*v100 - a121*u000*u011*u100*u111*v000*v020*v101)*v121)/((u001*u010*u101*u110*v001*v021*v100 - u001*u010*u100*u111*v000*v021*v101)*v120 - (u000*u011*u101*u110*v001*v020*v100 - u000*u011*u100*u111*v000*v020*v101)*v121),
         ((a221*u001*u010*u201*u210*v001*v021*v100 - a221*u000*u011*u201*u210*v000*v021*v101)*v120 - (a221*u001*u010*u200*u211*v001*v020*v100 - a221*u000*u011*u200*u211*v000*v020*v101)*v121)/((u001*u010*u201*u210*v001*v021*v100 - u001*u010*u200*u211*v000*v021*v101)*v120 - (u000*u011*u201*u210*v001*v020*v100 - u000*u011*u200*u211*v000*v020*v101)*v121),
         a031,
         ((a131*u001*u010*u101*u110*v001*v031*v100 - a131*u000*u011*u101*u110*v000*v031*v101)*v130 - (a131*u001*u010*u100*u111*v001*v030*v100 - a131*u000*u011*u100*u111*v000*v030*v101)*v131)/((u001*u010*u101*u110*v001*v031*v100 - u001*u010*u100*u111*v000*v031*v101)*v130 - (u000*u011*u101*u110*v001*v030*v100 - u000*u011*u100*u111*v000*v030*v101)*v131),
         ((a231*u001*u010*u201*u210*v001*v031*v100 - a231*u000*u011*u201*u210*v000*v031*v101)*v130 - (a231*u001*u010*u200*u211*v001*v030*v100 - a231*u000*u011*u200*u211*v000*v030*v101)*v131)/((u001*u010*u201*u210*v001*v031*v100 - u001*u010*u200*u211*v000*v031*v101)*v130 - (u000*u011*u201*u210*v001*v030*v100 - u000*u011*u200*u211*v000*v030*v101)*v131)]

 
    AUTHORS:
    - Edinah K. Gnang
    """
    if U.n(1)==U.n(2) and V.n(0)==V.n(2) and U.n(1)==V.n(0) and U.order()==3 and V.order()==3:
        # Initialization of the size parameter
        sz0=U.n(0); sz1=V.n(1); sz2=U.n(2)
        # Initialization of the container matrix
        M = Matrix(SR,HM(sz0*sz1*sz2, sz0*sz1*sz2, 'zero').listHM())
        for i in range(sz0):
            for j in range(sz1):
                for k in range(sz2):
                    for t in range(sz2):
                        M[i*sz1*sz2+j*sz2+k, i*sz1*sz2+j*sz2+t]=U[i,t,k]*V[t,j,k]
        # Initialization of the coefficient HM
        Ha=HM(sz0*sz1, sz0*sz1, [HM(sz2,sz2,'zero') for ij in range((sz0*sz1)^2)])
        for ij in range(sz0*sz1):
            Ha[ij, ij]=HM(sz2, sz2, M[ij*sz2:ij*sz2+sz2, ij*sz2:ij*sz2+sz2].transpose().list())
        # Initialization of the RHS
        Hb=HM(sz0*sz1, 1, [HM(2,sz2,'kronecker') for ij in range(sz0*sz1)])
        # Computing the matrix pseudo inverse
        [A,b]=gauss_jordan_eliminationHM(Ha, Hb)
        # Filling up the solution with the result
        B=HM(M.nrows(), M.ncols(), 'zero').matrix()
        for ij in range(sz0*sz1):
            B[sz2*ij:sz2*ij+sz2, sz2*ij:sz2*ij+sz2]=b[ij,0].simplify_full().matrix()
        # Initializing the multiplicative constraints.
        X=HM(U.n(0), U.n(1), U.n(2), 'x'); Y=HM(V.n(0), V.n(1), V.n(2), 'y')
        Eq=[X[i,s,t]*Y[s,j,t]==B[i*sz1*sz2+j*sz2+t,i*sz1*sz2+j*sz2+s] for i in range(sz0) for j in range(sz1) for s in range(sz2) for t in range(sz2)]
        # Formating the constraints
        [A,b]=multiplicativeConstraintFormator(Eq, X.list()+Y.list())
        Mx=Matrix(SR, A.ncols(), 1, X.list()+Y.list())
        return [[multiplicative_linear_solver(A,b,Mx,Mx), X, Y], multiplicative_gauss_jordan_eliminationII(A,b)]
    else:
        raise ValueError("The input hypermatrices must be cubical third order hypermatrices of the same sizes.")

def OuterHypermatrixInversePairHM(U, V):
    """
    Outputs the pseudo inverse pairs associated with the input pairs of hypermatrices
    The implementation does not assume that the input third order hypermatrices are
    cubic. In fact U must be m x p x p and V is p x n x p.


    EXAMPLES:

    ::

        sage: Hu=HM(2,2,2,'u'); Hv=HM(2,2,2,'v')
        sage: [Sln, Tx, Ty]=OuterHypermatrixInversePairHM(Hu, Hv)[0]
        sage: Hx=Tx.subs(dict([(s.lhs(),s.rhs()) for s in Sln if s.lhs()!=1])) 
        sage: Hy=Ty.subs(dict([(s.lhs(),s.rhs()) for s in Sln if s.lhs()!=1]))
        sage: Prod(Hx, Prod(Hu, HM(3,4,2,'a'), Hv), Hy).simplify_full().list()
        [a000,
         a100,
         ((a010*u001*u010*u101*u110*v001*v011*v100 - a010*u001*u010*u100*u111*v000*v011*v101)*v110 - (a010*u000*u011*u101*u110*v001*v010*v100 - a010*u000*u011*u100*u111*v000*v010*v101)*v111)/((u001*u010*u101*u110*v001*v011*v100 - u000*u011*u101*u110*v000*v011*v101)*v110 - (u001*u010*u100*u111*v001*v010*v100 - u000*u011*u100*u111*v000*v010*v101)*v111),
         a110,
         a001,
         a101,
         a011,
         ((a111*u001*u010*u101*u110*v001*v011*v100 - a111*u000*u011*u101*u110*v000*v011*v101)*v110 - (a111*u001*u010*u100*u111*v001*v010*v100 - a111*u000*u011*u100*u111*v000*v010*v101)*v111)/((u001*u010*u101*u110*v001*v011*v100 - u001*u010*u100*u111*v000*v011*v101)*v110 - (u000*u011*u101*u110*v001*v010*v100 - u000*u011*u100*u111*v000*v010*v101)*v111)]        

 
    AUTHORS:
    - Edinah K. Gnang
    """
    if U.n(1)==U.n(2) and V.n(0)==V.n(2) and U.n(1)==V.n(0) and U.order()==3 and V.order()==3:
        # Initialization of the size parameter
        sz0=U.n(0); sz1=V.n(1); sz2=U.n(2)
        # Initialization of the container matrix
        #M = Matrix(SR,HM(sz0*sz1*sz2, sz0*sz1*sz2, 'zero').listHM())
        M = HM(sz0*sz1*sz2, sz0*sz1*sz2, 'zero')
        for i in range(sz0):
            for j in range(sz1):
                for k in range(sz2):
                    for t in range(sz2):
                        M[i*sz1*sz2+j*sz2+k, i*sz1*sz2+j*sz2+t]=U[i,t,k]*V[t,j,k]
        # Initialization of the coefficient HM
        Ha=HM(sz0*sz1, sz0*sz1, [HM(sz2,sz2,'zero') for ij in range((sz0*sz1)^2)])
        for ij in range(sz0*sz1):
            #Ha[ij, ij]=HM(sz2, sz2, M[ij*sz2:ij*sz2+sz2, ij*sz2:ij*sz2+sz2].transpose().list())
            Ha[ij, ij]=M.slice(rg(ij*sz2,ij*sz2+sz2),0).slice(rg(ij*sz2,ij*sz2+sz2),1)
        # Initialization of the RHS
        Hb=HM(sz0*sz1, 1, [HM(2,sz2,'kronecker') for ij in range(sz0*sz1)])
        # Computing the matrix pseudo inverse
        [A,b]=gauss_jordan_eliminationHM(Ha, Hb)
        # Filling up the solution with the result
        #B=HM(M.nrows(), M.ncols(), 'zero').matrix()
        #for ij in range(sz0*sz1):
        #    B[sz2*ij:sz2*ij+sz2, sz2*ij:sz2*ij+sz2]=b[ij,0].simplify_full().matrix()
        B=HM(M.nrows(), M.ncols(), 'zero')
        for ij in range(sz0*sz1):
            for ii in rg(sz2*ij,sz2*ij+sz2):
                for jj in rg(sz2*ij,sz2*ij+sz2):
                    B[ii,jj]=(b[ij,0])[ii-sz2*ij,jj-sz2*ij].canonicalize_radical()
        # Initializing the multiplicative constraints.
        X=HM(U.n(0), U.n(1), U.n(2), 'x'); Y=HM(V.n(0), V.n(1), V.n(2), 'y')
        #Eq=[X[i,s,t]*Y[s,j,t]==B[i*sz1*sz2+j*sz2+t,i*sz1*sz2+j*sz2+s] for i in range(sz0) for j in range(sz1) for s in range(sz2) for t in range(sz2)]
        Eq=[X[i,s,t]*Y[s,j,t]==B[i*sz1*sz2+j*sz2+t,i*sz1*sz2+j*sz2+s] for i in rg(sz0) for j in rg(sz1) for s in rg(sz2) for t in rg(sz2)]
        # Formating the constraints
        #[A,b]=multiplicativeConstraintFormator(Eq, X.list()+Y.list())
        [A,b]=multiplicativeConstraintFormatorHM(Eq, X.list()+Y.list())
        #Mx=Matrix(SR, A.ncols(), 1, X.list()+Y.list())
        Mx=HM(A.ncols(), 1, X.list()+Y.list())
        #return [[multiplicative_linear_solver(A,b,Mx,Mx), X, Y], multiplicative_gauss_jordan_eliminationII(A,b)]
        return [[multiplicative_linear_solverHM(A,b,Mx,Mx), X, Y], multiplicative_gauss_jordan_eliminationHMII(A,b)]
    else:
        raise ValueError("The input hypermatrices must be cubical third order hypermatrices of the same sizes.")

def InnerHypermatrixInversePair(X, Y):
    """
    Outputs the pseudo inverse pairs associated with the input pairs of hypermatrices
    The implementation does not assume that the input third order hypermatrices are 
    cubic. In fact X is m x p x p and Y is p x n x p.


    EXAMPLES:

    ::

        sage: Hx=HM(2,2,2,'x'); Hy=HM(2,2,2,'y')
        sage: [Sln, Tu, Tv]=InnerHypermatrixInversePair(Hx, Hy)[0]
        sage: Hu=Tu.subs(dict([(s.lhs(),s.rhs()) for s in Sln if s.lhs() !=1])) 
        sage: Hv=Tv.subs(dict([(s.lhs(),s.rhs()) for s in Sln if s.lhs() !=1]))
        sage: Prod(Hx, Prod(Hu, HM(2,2,2,'a'), Hv), Hy).simplify_full().list()
        [a000,
         a100,
         ((a010*x001*x010*x101*x110*y001*y011*y100 - a010*x001*x010*x100*x111*y000*y011*y101)*y110 - (a010*x000*x011*x101*x110*y001*y010*y100 - a010*x000*x011*x100*x111*y000*y010*y101)*y111)/((x001*x010*x101*x110*y001*y011*y100 - x000*x011*x101*x110*y000*y011*y101)*y110 - (x001*x010*x100*x111*y001*y010*y100 - x000*x011*x100*x111*y000*y010*y101)*y111),
         a110,
         a001,
         a101,
         ((a011*x001*x010*x101*x110*y001*y011*y100 - a011*x001*x010*x100*x111*y000*y011*y101)*y110 - (a011*x000*x011*x101*x110*y001*y010*y100 - a011*x000*x011*x100*x111*y000*y010*y101)*y111)/((x001*x010*x101*x110*y001*y011*y100 - x000*x011*x101*x110*y000*y011*y101)*y110 - (x001*x010*x100*x111*y001*y010*y100 - x000*x011*x100*x111*y000*y010*y101)*y111),
         a111]


    AUTHORS:
    - Edinah K. Gnang
    """
    if X.n(1)==X.n(2) and Y.n(0)==Y.n(2) and X.n(1)==Y.n(0) and X.order()==3 and Y.order()==3:
        # Initialization of the size parameter
        sz0=X.n(0); sz1=Y.n(1); sz2=Y.n(2)
        # Initialization of the container matrix
        M = Matrix(SR,HM(sz0*sz1*sz2, sz0*sz1*sz2, 'zero').listHM())
        for i in range(sz0):
            for j in range(sz1):
                for s in range(sz2):
                    for t in range(sz2):
                        M[i*sz1*sz2+j*sz2+t,sz1*sz2*i+sz2*j+s]=X[i,s,t]*Y[s,j,t]
        # Initialization of the coefficient HM
        Ha=HM(sz0*sz1, sz0*sz1, [HM(sz2,sz2,'zero') for ij in range((sz0*sz1)^2)])
        for ij in range(sz0*sz1):
            Ha[ij, ij]=HM(sz2, sz2, M[ij*sz2:ij*sz2+sz2, ij*sz2:ij*sz2+sz2].transpose().list())
        # Initialization of the RHS
        Hb=HM(sz0*sz1, 1, [HM(2,sz2,'kronecker') for ij in range(sz0*sz1)])
        # Computing the matrix pseudo inverse
        [A,b]=gauss_jordan_eliminationHM(Ha, Hb)
        # Filling up the solution with the result
        B=HM(M.nrows(), M.ncols(), 'zero').matrix()
        for ij in range(sz0*sz1):
            B[sz2*ij:sz2*ij+sz2, sz2*ij:sz2*ij+sz2]=b[ij,0].simplify_full().matrix()
        # Initializing the multiplicative constraints.
        U=HM(X.n(0),X.n(1),X.n(2),'u'); V=HM(Y.n(0),Y.n(1),Y.n(2),'v')
        Eq=[U[i,t,k]*V[t,j,k]==B[i*sz1*sz2+j*sz2+k,i*sz1*sz2+j*sz2+t] for i in range(sz0) for j in range(sz1) for k in range(sz2) for t in range(sz2)]
        # Formating the constraints
        [A,b]=multiplicativeConstraintFormator(Eq, U.list()+V.list())
        Mx=Matrix(SR, A.ncols(), 1, U.list()+V.list())
        return [[multiplicative_linear_solver(A,b,Mx,Mx), U, V], multiplicative_gauss_jordan_eliminationII(A,b)]
    else:
        raise ValueError("The input hypermatrices must be cubical third order hypermatrices of the same sizes.")

def InnerHypermatrixInversePairHM(X, Y):
    """
    Outputs the pseudo inverse pairs associated with the input pairs of hypermatrices
    The implementation does not assume that the input third order hypermatrices are 
    cubic. In fact X is m x p x p and Y is p x n x p.


    EXAMPLES:

    ::

        sage: Hx=HM(2,2,2,'x'); Hy=HM(2,2,2,'y')
        sage: [Sln, Tu, Tv]=InnerHypermatrixInversePairHM(Hx, Hy)[0]
        sage: Hu=Tu.subs(dict([(s.lhs(),s.rhs()) for s in Sln if s.lhs() !=1])) 
        sage: Hv=Tv.subs(dict([(s.lhs(),s.rhs()) for s in Sln if s.lhs() !=1]))
        sage: Prod(Hx, Prod(Hu, HM(2,2,2,'a'), Hv), Hy).simplify_full().list()
        [a000,
         a100,
         ((a010*x001*x010*x101*x110*y001*y011*y100 - a010*x001*x010*x100*x111*y000*y011*y101)*y110 - (a010*x000*x011*x101*x110*y001*y010*y100 - a010*x000*x011*x100*x111*y000*y010*y101)*y111)/((x001*x010*x101*x110*y001*y011*y100 - x000*x011*x101*x110*y000*y011*y101)*y110 - (x001*x010*x100*x111*y001*y010*y100 - x000*x011*x100*x111*y000*y010*y101)*y111),
         a110,
         a001,
         a101,
         ((a011*x001*x010*x101*x110*y001*y011*y100 - a011*x001*x010*x100*x111*y000*y011*y101)*y110 - (a011*x000*x011*x101*x110*y001*y010*y100 - a011*x000*x011*x100*x111*y000*y010*y101)*y111)/((x001*x010*x101*x110*y001*y011*y100 - x000*x011*x101*x110*y000*y011*y101)*y110 - (x001*x010*x100*x111*y001*y010*y100 - x000*x011*x100*x111*y000*y010*y101)*y111),
         a111] 


    AUTHORS:
    - Edinah K. Gnang
    """
    if X.n(1)==X.n(2) and Y.n(0)==Y.n(2) and X.n(1)==Y.n(0) and X.order()==3 and Y.order()==3:
        # Initialization of the size parameter
        sz0=X.n(0); sz1=Y.n(1); sz2=Y.n(2)
        # Initialization of the container matrix
        #M = Matrix(SR,HM(sz0*sz1*sz2, sz0*sz1*sz2, 'zero').listHM())
        M = HM(sz0*sz1*sz2, sz0*sz1*sz2, 'zero')
        for i in range(sz0):
            for j in range(sz1):
                for s in range(sz2):
                    for t in range(sz2):
                        M[i*sz1*sz2+j*sz2+t,sz1*sz2*i+sz2*j+s]=X[i,s,t]*Y[s,j,t]
        # Initialization of the coefficient HM
        Ha=HM(sz0*sz1, sz0*sz1, [HM(sz2,sz2,'zero') for ij in range((sz0*sz1)^2)])
        for ij in range(sz0*sz1):
            #Ha[ij, ij]=HM(sz2, sz2, M[ij*sz2:ij*sz2+sz2, ij*sz2:ij*sz2+sz2].transpose().list())
            Ha[ij, ij]=M.slice(rg(ij*sz2,ij*sz2+sz2),0).slice(rg(ij*sz2,ij*sz2+sz2),1)
        # Initialization of the RHS
        Hb=HM(sz0*sz1, 1, [HM(2,sz2,'kronecker') for ij in range(sz0*sz1)])
        # Computing the matrix pseudo inverse
        [A,b]=gauss_jordan_eliminationHM(Ha, Hb)
        # Filling up the solution with the result
        #B=HM(M.nrows(), M.ncols(), 'zero').matrix()
        #for ij in range(sz0*sz1):
        #    B[sz2*ij:sz2*ij+sz2, sz2*ij:sz2*ij+sz2]=b[ij,0].simplify_full().matrix()
        B=HM(M.nrows(), M.ncols(), 'zero')
        for ij in range(sz0*sz1):
            for ii in rg(sz2*ij,sz2*ij+sz2):
                for jj in rg(sz2*ij,sz2*ij+sz2):
                    B[ii,jj]=(b[ij,0])[ii-sz2*ij,jj-sz2*ij].canonicalize_radical()
        # Initializing the multiplicative constraints.
        U=HM(X.n(0),X.n(1),X.n(2),'u'); V=HM(Y.n(0),Y.n(1),Y.n(2),'v')
        #Eq=[U[i,t,k]*V[t,j,k]==B[i*sz1*sz2+j*sz2+k,i*sz1*sz2+j*sz2+t] for i in range(sz0) for j in range(sz1) for k in range(sz2) for t in range(sz2)]
        Eq=[U[i,t,k]*V[t,j,k]==B[i*sz1*sz2+j*sz2+k,i*sz1*sz2+j*sz2+t] for i in rg(sz0) for j in rg(sz1) for k in rg(sz2) for t in rg(sz2)]
        # Formating the constraints
        #[A,b]=multiplicativeConstraintFormator(Eq, U.list()+V.list())
        [A,b]=multiplicativeConstraintFormatorHM(Eq, U.list()+V.list())
        #Mx=Matrix(SR, A.ncols(), 1, U.list()+V.list())
        Mx=HM(A.ncols(), 1, U.list()+V.list())
        #return [[multiplicative_linear_solver(A,b,Mx,Mx), U, V], multiplicative_gauss_jordan_eliminationII(A,b)]
        return [[multiplicative_linear_solverHM(A,b,Mx,Mx), U, V], multiplicative_gauss_jordan_eliminationHMII(A,b)]
    else:
        raise ValueError("The input hypermatrices must be cubical third order hypermatrices of the same sizes.")

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

def CountCompositions(n):
    """
    Counts the number of product composition involving the input
    hypermatrix n times.

    EXAMPLES:
    The input n must be greater than 0
    ::

        sage: CountCompositions(3)
        1

    AUTHORS:
    - Edinah K. Gnang and Ori Parzanchevski
    """
    if n == 1 :
        return 1
    else :
        return sum([CountCompositions(i)*CountCompositions(j)*CountCompositions(n-i-j) for i in range(1,n,2) for j in range(1,n-i,2)])

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
        sage: Ha=HM(2,2,'a'); Hb=HM(2,2,'b'); GeneralHypermatrixProduct(Ha, Hb)
        [[a00*b00 + a01*b10, a00*b01 + a01*b11], [a10*b00 + a11*b10, a10*b01 + a11*b11]]


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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # computing the Hypermatrix product
        if len(args)<2:
            raise ValueError("The number of operands must be >= 2")
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
        sage: Prod(HM(2,2,'a'),HM(2,2,'b'))
        [[a00*b00 + a01*b10, a00*b01 + a01*b11], [a10*b00 + a11*b10, a10*b01 + a11*b11]]


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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # Computing the Hypermatrix product
        if len(args) < 2:
            raise ValueError("The number of operands must be >= 2")
        elif len(args) >= 2:
            Rh[tuple(entry)] = 0
            l2 = [B.n(sz) for sz in range(B.order())]
            for j in range(prod(l2)):
                # Turning the index j into an hypermatrix array location using the decimal encoding trick
                entry2 = [Integer(mod(j,l2[0]))]
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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # computing the Hypermatrix product
        if len(args)<2:
            raise ValueError("The number of operands must be >= 2")
        elif len(args) >= 2:
            Rh[tuple(entry)]=sum([sum([args[s][tuple(entry[0:Integer(mod(s+1,len(args)))]+[t]+entry[Integer(mod(s+2,len(args))):])] for s in range(len(args)-2)]+[args[len(args)-2][tuple(entry[0:len(args)-1]+[t])]]+[args[len(args)-1][tuple([t]+entry[1:])]]) for t in range((args[0]).n(1))])
    return Rh

def GeneralHypermatrixBlockProduct(*args):
    """  
    Outputs a list of lists associated with the general
    Bhattacharya-Mesner block product of the input hypermatrices.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:

    ::   

        sage: Ha=HM(2,2,2,'a');Hb=HM(2,2,2,'b');Hc=HM(2,2,2,'c');Hd=HM(2,2,2,'d');Hf=HM(2,2,2,'f');Hg=HM(2,2,2,'g')
        sage: A=HM(1,2,1,[Ha,Hb]); B=HM(1,1,2,[Hc,Hd]); C=HM(2,1,1,[Hf,Hg])
        sage: Rslt=GeneralHypermatrixBlockProduct(A, B, C); Rslt[0,0,0].printHM()
        [:, :, 0]=
        [a000*c000*f000 + a010*c001*f100 + b000*d000*g000 + b010*d001*g100 a000*c010*f010 + a010*c011*f110 + b000*d010*g010 + b010*d011*g110]
        [a100*c100*f000 + a110*c101*f100 + b100*d100*g000 + b110*d101*g100 a100*c110*f010 + a110*c111*f110 + b100*d110*g010 + b110*d111*g110]
        <BLANKLINE>
        [:, :, 1]=
        [a001*c000*f001 + a011*c001*f101 + b001*d000*g001 + b011*d001*g101 a001*c010*f011 + a011*c011*f111 + b001*d010*g011 + b011*d011*g111]
        [a101*c100*f001 + a111*c101*f101 + b101*d100*g001 + b111*d101*g101 a101*c110*f011 + a111*c111*f111 + b101*d110*g011 + b111*d111*g111]


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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # computing the Hypermatrix product
        if len(args)<2:
            raise ValueError("The number of operands must be >= 2")
        elif len(args) >= 2:
            Rh[tuple(entry)]=sum([Prod( *([args[s][tuple(entry[0:Integer(mod(s+1,len(args)))]+[t]+entry[Integer(mod(s+2,len(args))):])] for s in range(len(args)-2)]+[args[len(args)-2][tuple(entry[0:len(args)-1]+[t])]]+[args[len(args)-1][tuple([t]+entry[1:])]]) ) for t in range((args[0]).n(1))])
    return Rh

def BlockProd(*args):
    """  
    Outputs a list of lists associated with the general
    Bhattacharya-Mesner block product of the input hypermatrices.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:

    ::   

        sage: Ha=HM(2,2,2,'a');Hb=HM(2,2,2,'b');Hc=HM(2,2,2,'c');Hd=HM(2,2,2,'d');Hf=HM(2,2,2,'f');Hg=HM(2,2,2,'g')
        sage: A=HM(1,2,1,[Ha,Hb]); B=HM(1,1,2,[Hc,Hd]); C=HM(2,1,1,[Hf,Hg])
        sage: Rslt=BlockProd(A, B, C); Rslt[0,0,0].printHM()
        [:, :, 0]=
        [a000*c000*f000 + a010*c001*f100 + b000*d000*g000 + b010*d001*g100 a000*c010*f010 + a010*c011*f110 + b000*d010*g010 + b010*d011*g110]
        [a100*c100*f000 + a110*c101*f100 + b100*d100*g000 + b110*d101*g100 a100*c110*f010 + a110*c111*f110 + b100*d110*g010 + b110*d111*g110]
        <BLANKLINE>
        [:, :, 1]=
        [a001*c000*f001 + a011*c001*f101 + b001*d000*g001 + b011*d001*g101 a001*c010*f011 + a011*c011*f111 + b001*d010*g011 + b011*d011*g111]
        [a101*c100*f001 + a111*c101*f101 + b101*d100*g001 + b111*d101*g101 a101*c110*f011 + a111*c111*f111 + b101*d110*g011 + b111*d111*g111]
        sage: Ha=HM(2,1,2,'a');Hb=HM(2,1,2,'b');Hc=HM(2,2,1,'c');Hd=HM(2,2,1,'d');Hf=HM(1,2,2,'f');Hg=HM(1,2,2,'g')
        sage: A=HM(1,2,1,[Ha,Hb]); B=HM(1,1,2,[Hc,Hd]); C=HM(2,1,1,[Hf,Hg])
        sage: BlockProd(A, B, C)[0,0,0].printHM() # Computing the block product
        [:, :, 0]=
        [a000*c000*f000 + b000*d000*g000 a000*c010*f010 + b000*d010*g010]
        [a100*c100*f000 + b100*d100*g000 a100*c110*f010 + b100*d110*g010]
        <BLANKLINE>
        [:, :, 1]=
        [a001*c000*f001 + b001*d000*g001 a001*c010*f011 + b001*d010*g011]
        [a101*c100*f001 + b101*d100*g001 a101*c110*f011 + b101*d110*g011]


    AUTHORS:
    - Edinah K. Gnang
    """
    return GeneralHypermatrixBlockProduct(*args) 

def GeneralHypermatrixProductII(Lh, Op, F):
    """
    Outputs an HM which is a list of lists associated with the
    construct approach to the Bhattacharya-Mesner product of the input
    hypermatrices. Here Op is the combinator and F is the composer.
    This implementation in theory captures the full scope of
    the construct products subsequently implemented here but in 
    practice is hard for to specify arbitrary composers for F.
    The code only handles the Hypermatrix HM class objects.


    EXAMPLES:

    ::

        sage: Ha=HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c')
        sage: Rslt=GeneralHypermatrixProductII([Ha, Hb, Hc], prod, sum); Rslt.printHM()
        [:, :, 0]=
        [(a000 + b000 + c000)*(a010 + b001 + c100) (a000 + b010 + c010)*(a010 + b011 + c110)]
        [(a100 + b100 + c000)*(a110 + b101 + c100) (a100 + b110 + c010)*(a110 + b111 + c110)]
        <BLANKLINE>
        [:, :, 1]=
        [(a001 + b000 + c001)*(a011 + b001 + c101) (a001 + b010 + c011)*(a011 + b011 + c111)]
        [(a101 + b100 + c001)*(a111 + b101 + c101) (a101 + b110 + c011)*(a111 + b111 + c111)]        
        
        sage: Rslt=GeneralHypermatrixProductII([Ha, Hb, Hc], sum, prod); Rslt.printHM()
        [:, :, 0]=
        [a000*b000*c000 + a010*b001*c100 a000*b010*c010 + a010*b011*c110]
        [a100*b100*c000 + a110*b101*c100 a100*b110*c010 + a110*b111*c110]
        <BLANKLINE>
        [:, :, 1]=
        [a001*b000*c001 + a011*b001*c101 a001*b010*c011 + a011*b011*c111]
        [a101*b100*c001 + a111*b101*c101 a101*b110*c011 + a111*b111*c111]

        sage: Ha=HM(2,2,'a'); Hb=HM(2,2,'b'); GeneralHypermatrixProductII([Ha, Hb], sum, prod)
        [[a00*b00 + a01*b10, a00*b01 + a01*b11], [a10*b00 + a11*b10, a10*b01 + a11*b11]]
        sage: A=HM([[59, -3, 2], [-6, 1, 1], [1, -1, 1]]); B=HM([[-1, 1, 0], [-1, 1, 0], [0, -43, 1]])
        sage: GeneralHypermatrixProductII([A, B], sum, prod)-A*B # One way of recovering matrix multiplication
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        sage: GeneralHypermatrixProductII([A, B], min, sum) # Recovering the Min-plus matrix multiplication
        [[-4, -41, -3], [-7, -42, -6], [-2, -42, -1]]
        sage: MSta=HM([[Set([1,2]), Set([1,3,2])], [Set([1]), Set([2,3])]])
        sage: MStb=HM([[Set([1,2,3]), Set([2])], [Set([1,3]), Set([1,3])]])
        sage: GeneralHypermatrixProductII([MSta, MStb], SetIntersection, SetUnion) 
        [[{1, 2, 3}, {1, 2}], [{1, 2, 3}, {1, 2}]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list specifying the dimensions of the output
    l = [(Lh[i]).n(i) for i in range(len(Lh))]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the assignement
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # computing the Hypermatrix product
        if len(Lh)<2:
            raise ValueError("The number of operands must be >= 2")
        elif len(Lh) >= 2:
            #Rh[tuple(entry)]=apply(Op, [[apply(F, [[Lh[s][tuple(entry[0:Integer(mod(s+1,len(Lh)))]+[t]+entry[Integer(mod(s+2,len(Lh))):])] for s in range(len(Lh)-2)]+[Lh[len(Lh)-2][tuple(entry[0:len(Lh)-1]+[t])]]+[Lh[len(Lh)-1][tuple([t]+entry[1:])]]]) for t in range((Lh[0]).n(1))]])
            Rh[tuple(entry)]=Op(*[[F( *[[Lh[s][tuple(entry[0:Integer(mod(s+1,len(Lh)))]+[t]+entry[Integer(mod(s+2,len(Lh))):])] for s in range(len(Lh)-2)]+[Lh[len(Lh)-2][tuple(entry[0:len(Lh)-1]+[t])]]+[Lh[len(Lh)-1][tuple([t]+entry[1:])]]] ) for t in range((Lh[0]).n(1))]])
    return Rh

def GeneralHypermatrixProductIII(Lh, Op, Lv):
    """
    Outputs a list of lists associated with the composition
    based Bhattacharya-Mesner product of the input hypermatrices.
    The entries of the hypermatrices are taken to be functions
    so that while performing the product we compose with the entries
    of the first of the list of inputs. This implementation is aesthetically
    more pleasing then the previous one because it explicitly articulate the
    preference for composition as our defacto product operation. Hoever it is
    theoretically less general then the previous one. Both these implementations
    are inspired by initial exposure to ideas from Category theory, the implementation
    also make painfully obvious some of the programming constraints imposed by Python.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:

    ::

        sage: x,y=var('x,y'); Ha=x*y*HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c')
        sage: Rslt=GeneralHypermatrixProductIII([Ha,Hb,Hc], sum, [x,y]); Rslt.printHM()
        [:, :, 0]=
        [a000*b000*c000 + a010*b001*c100 a000*b010*c010 + a010*b011*c110]
        [a100*b100*c000 + a110*b101*c100 a100*b110*c010 + a110*b111*c110]
        <BLANKLINE>
        [:, :, 1]=
        [a001*b000*c001 + a011*b001*c101 a001*b010*c011 + a011*b011*c111]
        [a101*b100*c001 + a111*b101*c101 a101*b110*c011 + a111*b111*c111]

        sage: x=var('x'); Ha=x*HM(2,2,'a'); Hb=HM(2,2,'b'); GeneralHypermatrixProductIII([Ha, Hb], sum, [x])
        [[a00*b00 + a01*b10, a00*b01 + a01*b11], [a10*b00 + a11*b10, a10*b01 + a11*b11]]
        sage: Ha=HM(2,2,'a').elementwise_exponent(x); Hb=HM(2,2,'b'); GeneralHypermatrixProductIII([Ha, Hb], prod, [x])
        [[a00^b00*a01^b10, a00^b01*a01^b11], [a10^b00*a11^b10, a10^b01*a11^b11]]



    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list specifying the dimensions of the output
    l = [(Lh[i]).n(i) for i in range(len(Lh))]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the assignement
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # computing the Hypermatrix product
        if len(Lh)<2:
            raise ValueError("The number of operands must be >= 2")
        elif len(Lh) >= 2:
            #Rh[tuple(entry)]=apply(Op, [ [([Lh[s][tuple(entry[0:Integer(mod(s+1,len(Lh)))]+[t]+entry[Integer(mod(s+2,len(Lh))):])] for s in range(len(Lh)-2)]+[Lh[len(Lh)-2][tuple(entry[0:len(Lh)-1]+[t])]]+[Lh[len(Lh)-1][tuple([t]+entry[1:])]])[0].subs([Lv[z-1]==([Lh[s][tuple(entry[0:Integer(mod(s+1,len(Lh)))]+[t]+entry[Integer(mod(s+2,len(Lh))):])] for s in range(len(Lh)-2)]+[Lh[len(Lh)-2][tuple(entry[0:len(Lh)-1]+[t])]]+[Lh[len(Lh)-1][tuple([t]+entry[1:])]])[z] for z in rg(1,len(Lh))]) for t in range((Lh[0]).n(1))] ])
            Rh[tuple(entry)]=Op( *[ [([Lh[s][tuple(entry[0:Integer(mod(s+1,len(Lh)))]+[t]+entry[Integer(mod(s+2,len(Lh))):])] for s in range(len(Lh)-2)]+[Lh[len(Lh)-2][tuple(entry[0:len(Lh)-1]+[t])]]+[Lh[len(Lh)-1][tuple([t]+entry[1:])]])[0].subs([Lv[z-1]==([Lh[s][tuple(entry[0:Integer(mod(s+1,len(Lh)))]+[t]+entry[Integer(mod(s+2,len(Lh))):])] for s in range(len(Lh)-2)]+[Lh[len(Lh)-2][tuple(entry[0:len(Lh)-1]+[t])]]+[Lh[len(Lh)-1][tuple([t]+entry[1:])]])[z] for z in rg(1,len(Lh))]) for t in range((Lh[0]).n(1))] ])
    return Rh

def GeneralHypermatrixProductIV(Lh, Op, Lv):
    """
    Outputs a list of lists associated with the composition
    based Bhattacharya-Mesner product of the input hypermatrices.
    The entries of the hypermatrices are taken to be functions
    so that while performing the product we compose with the entries
    of the first of the list of inputs. This implementation is aesthetically
    more pleasing then the previous one because it explicitly articulate the
    preference for composition as our defacto product operation. Hoever it is
    theoretically less general then the previous one. Both these implementations
    are inspired by initial exposure to ideas from Category theory, the implementation
    also make painfully obvious some of the programming constraints imposed by Python.
    The code only handles the Hypermatrix HM class objects. This new variant is to
    accomodate substitution with free variables, it is mainly a work around a sage
    bug the substitute function behaves very differently for free field.
    This implmentation suggest that free variables which support inversions would
    be a great feature to have. I will keep this in mind as I get acquainted with
    sage developers.


    EXAMPLES:

    ::

        sage: sz=3; l=2; Lx=var_list('x',l); La=HM(sz,l,'a').list(); Lb=HM(sz,l,'b').list(); z=var('z')
        sage: F=FreeAlgebra(QQ,len(La+Lx+Lb+[z]),La+Lx+Lb+[z])
        sage: F.<a00, a10, a20, a01, a11, a21, x0, x1, b00, b10, b20, b01, b11, b21, z>=FreeAlgebra(QQ,len(La+Lx+Lb+[z]))
        sage: Ha=HM(sz,l,[a00, a10, a20, a01, a11, a21]); Hx=HM(l,1,[x0, x1])
        sage: Hb=HM(sz,l,[b00, b10, b20, b01, b11, b21])
        sage: Hr=Ha.elementwise_product(z*Hb); GeneralHypermatrixProductIV([Hr, Hx], sum, [z]).printHM()
        [:, :]=
        [a00*x0*b00 + a01*x1*b01]
        [a10*x0*b10 + a11*x1*b11]
        [a20*x0*b20 + a21*x1*b21]
        sage: x,y=var('x,y'); Ha=x*y*HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c')
        sage: Rslt=GeneralHypermatrixProductIV([Ha,Hb,Hc], sum, [x,y]); Rslt.printHM()
        [:, :, 0]=
        [a000*b000*c000 + a010*b001*c100 a000*b010*c010 + a010*b011*c110]
        [a100*b100*c000 + a110*b101*c100 a100*b110*c010 + a110*b111*c110]
        <BLANKLINE>
        [:, :, 1]=
        [a001*b000*c001 + a011*b001*c101 a001*b010*c011 + a011*b011*c111]
        [a101*b100*c001 + a111*b101*c101 a101*b110*c011 + a111*b111*c111]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Testing conformability by initiailizing the common dimension
    el=Lh[len(Lh)-1].n(0)
    # Loop running throug the dimentsion of the input
    for indx in rg(len(Lh)-1):
        if Lh[indx].n(indx+1) != el:
            print('WARNING !!! The input hypermatrices are not conformable. Truncating to allow the product.')
            break
    # Initialization of the list specifying the dimensions of the output
    l = [(Lh[i]).n(i) for i in range(len(Lh))]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the assignement
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # computing the Hypermatrix product
        if len(Lh)<2:
            raise ValueError("The number of operands must be >= 2")
        elif len(Lh) >= 2:
            Ltmp=[]
            for t in range((Lh[0]).n(1)):
                Tmp=([Lh[s][tuple(entry[0:Integer(mod(s+1,len(Lh)))]+[t]+entry[Integer(mod(s+2,len(Lh))):])] for s in range(len(Lh)-2)]+[Lh[len(Lh)-2][tuple(entry[0:len(Lh)-1]+[t])]]+[Lh[len(Lh)-1][tuple([t]+entry[1:])]])[0]
                #print 'Tmp before', Tmp
                #print 'Ltmp before', Ltmp
                for vz in rg(1,len(Lh)):
                    #print {Lv[vz-1]:([Lh[s][tuple(entry[0:Integer(mod(s+1,len(Lh)))]+[t]+entry[Integer(mod(s+2,len(Lh))):])] for s in range(len(Lh)-2)]+[Lh[len(Lh)-2][tuple(entry[0:len(Lh)-1]+[t])]]+[Lh[len(Lh)-1][tuple([t]+entry[1:])]])[vz]}
                    Tmp=Tmp.subs({Lv[vz-1]:([Lh[s][tuple(entry[0:Integer(mod(s+1,len(Lh)))]+[t]+entry[Integer(mod(s+2,len(Lh))):])] for s in range(len(Lh)-2)]+[Lh[len(Lh)-2][tuple(entry[0:len(Lh)-1]+[t])]]+[Lh[len(Lh)-1][tuple([t]+entry[1:])]])[vz]})
                    #print 'Tmp after', Tmp
                Ltmp.append(Tmp)
                #print 'Ltmp after', Ltmp
            #print apply(Op, [Ltmp])
            #Rh[tuple(entry)]=apply(Op, [Ltmp])
            Rh[tuple(entry)]=Op(*[Ltmp])
    return Rh

def GeneralHypermatrixProductV(Lh, Op, Lv, indx):
    """
    Outputs a list of lists associated with the composition
    based Bhattacharya-Mesner product of the input hypermatrices.
    The entries of the hypermatrices are taken to be functions
    so that while performing the product we compose with the entries
    of the first of the list of inputs. This implementation is aesthetically
    more pleasing then the previous one because it explicitly articulate the
    preference for composition as our defacto product operation. Hoever it is
    theoretically less general then the previous one. Both these implementations
    are inspired by initial exposure to ideas from Category theory, the implementation
    also make painfully obvious some of the programming constraints imposed by Python.
    The code only handles the Hypermatrix HM class objects. This new variant is to
    accomodate substitution with free variables, it is mainly a work around a sage
    bug the substitute function behaves very differently for free field.
    This implmentation suggest that free variables which support inversions would
    be a great feature to have. I will keep this in mind as I get acquainted with
    sage developers. The difference with the previous function is that this support
    an indexing parameter for the hypermatrix input in which the substitions will be
    performed. This gives a tiny bit more flexibility with setting the composer.


    EXAMPLES:

    ::

        sage: sz=3; l=2; Lx=var_list('x',l); La=HM(sz,l,'a').list(); Lb=HM(sz,l,'b').list(); z=var('z')
        sage: F=FreeAlgebra(QQ,len(La+Lx+Lb+[z]),La+Lx+Lb+[z])
        sage: F.<a00, a10, a20, a01, a11, a21, x0, x1, b00, b10, b20, b01, b11, b21, z>=FreeAlgebra(QQ,len(La+Lx+Lb+[z]))
        sage: Ha=HM(sz,l,[a00, a10, a20, a01, a11, a21]); Hx=HM(l,1,[x0, x1])
        sage: Hb=HM(sz,l,[b00, b10, b20, b01, b11, b21])
        sage: Hr=Ha.elementwise_product(z*Hb); GeneralHypermatrixProductV([Hr, Hx], sum, [z], 0).printHM()
        [:, :]=
        [a00*x0*b00 + a01*x1*b01]
        [a10*x0*b10 + a11*x1*b11]
        [a20*x0*b20 + a21*x1*b21]
        sage: x,y=var('x,y'); Ha=HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c')
        sage: GeneralHypermatrixProductV([x*y*Ha,Hb,Hc], sum, [x,y], 0).printHM()
        [:, :, 0]=
        [a000*b000*c000 + a010*b001*c100 a000*b010*c010 + a010*b011*c110]
        [a100*b100*c000 + a110*b101*c100 a100*b110*c010 + a110*b111*c110]
        <BLANKLINE>
        [:, :, 1]=
        [a001*b000*c001 + a011*b001*c101 a001*b010*c011 + a011*b011*c111]
        [a101*b100*c001 + a111*b101*c101 a101*b110*c011 + a111*b111*c111]
        sage: x,y=var('x,y'); Ha=HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c')
        sage: GeneralHypermatrixProductV([Ha,x*y*Hb,Hc], sum, [x,y], 1).printHM()
        [:, :, 0]=
        [a000*b000*c000 + a010*b001*c100 a000*b010*c010 + a010*b011*c110]
        [a100*b100*c000 + a110*b101*c100 a100*b110*c010 + a110*b111*c110]
        <BLANKLINE>
        [:, :, 1]=
        [a001*b000*c001 + a011*b001*c101 a001*b010*c011 + a011*b011*c111]
        [a101*b100*c001 + a111*b101*c101 a101*b110*c011 + a111*b111*c111]
        sage: GeneralHypermatrixProductV([x*y*Ha,Hb,Hc], sum, [x,y], 0).printHM()
        [:, :, 0]=
        [a000*b000*c000 + a010*b001*c100 a000*b010*c010 + a010*b011*c110]
        [a100*b100*c000 + a110*b101*c100 a100*b110*c010 + a110*b111*c110]
        <BLANKLINE>
        [:, :, 1]=
        [a001*b000*c001 + a011*b001*c101 a001*b010*c011 + a011*b011*c111]
        [a101*b100*c001 + a111*b101*c101 a101*b110*c011 + a111*b111*c111]


    AUTHORS:
    - Edinah K. Gnang
    """
    if indx == 0:
        return GeneralHypermatrixProductIV(Lh, Op, Lv)
    else:
        # Initialization of the list specifying the dimensions of the output
        l = [(Lh[i]).n(i) for i in range(len(Lh))]
        # Initializing the input for generating a symbolic hypermatrix
        inpts = l+['zero']
        # Initialization of the hypermatrix
        Rh = HM(*inpts)
        # Main loop performing the assignement
        for i in range(prod(l)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entry = [Integer(mod(i,l[0]))]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            # computing the Hypermatrix product
            if len(Lh)<2:
                raise ValueError("The number of operands must be >= 2")
            elif len(Lh) >= 2:
                Ltmp=[]
                for t in range((Lh[0]).n(1)):
                    Tmp=([Lh[s][tuple(entry[0:Integer(mod(s+1,len(Lh)))]+[t]+entry[Integer(mod(s+2,len(Lh))):])] for s in range(len(Lh)-2)]+[Lh[len(Lh)-2][tuple(entry[0:len(Lh)-1]+[t])]]+[Lh[len(Lh)-1][tuple([t]+entry[1:])]])[indx]
                    #print 'Tmp before', Tmp
                    #print 'Ltmp before', Ltmp
                    LindX1 = rg(indx)
                    #print 'indx=', indx
                    #print 'LindX1=', LindX1
                    for vz in LindX1:
                        #print {Lv[vz]:([Lh[s][tuple(entry[0:Integer(mod(s+1,len(Lh)))]+[t]+entry[Integer(mod(s+2,len(Lh))):])] for s in range(len(Lh)-2)]+[Lh[len(Lh)-2][tuple(entry[0:len(Lh)-1]+[t])]]+[Lh[len(Lh)-1][tuple([t]+entry[1:])]])[vz]}
                        Tmp=Tmp.subs({Lv[vz]:([Lh[s][tuple(entry[0:Integer(mod(s+1,len(Lh)))]+[t]+entry[Integer(mod(s+2,len(Lh))):])] for s in range(len(Lh)-2)]+[Lh[len(Lh)-2][tuple(entry[0:len(Lh)-1]+[t])]]+[Lh[len(Lh)-1][tuple([t]+entry[1:])]])[vz]})
                        #print 'Tmp after', Tmp
                    #print 'LindX2=', LindX2
                    LindX2 = rg(1+indx,1+len(Lh)-indx)
                    for vz in LindX2:
                        #print {Lv[vz-1]:([Lh[s][tuple(entry[0:Integer(mod(s+1,len(Lh)))]+[t]+entry[Integer(mod(s+2,len(Lh))):])] for s in range(len(Lh)-2)]+[Lh[len(Lh)-2][tuple(entry[0:len(Lh)-1]+[t])]]+[Lh[len(Lh)-1][tuple([t]+entry[1:])]])[vz]}
                        Tmp=Tmp.subs({Lv[vz-1]:([Lh[s][tuple(entry[0:Integer(mod(s+1,len(Lh)))]+[t]+entry[Integer(mod(s+2,len(Lh))):])] for s in range(len(Lh)-2)]+[Lh[len(Lh)-2][tuple(entry[0:len(Lh)-1]+[t])]]+[Lh[len(Lh)-1][tuple([t]+entry[1:])]])[vz]})
                        #print 'Tmp after', Tmp
                    Ltmp.append(Tmp)
                    #print 'Ltmp after', Ltmp
                #print apply(Op, [Ltmp])
                #Rh[tuple(entry)]=apply(Op, [Ltmp])
                Rh[tuple(entry)]=Op(*[Ltmp])
        return Rh

def GProd(Lh, Op, Lv):
    """
    Outputs a list of lists associated with the composition
    based Bhattacharya-Mesner product of the input hypermatrices.
    The entries of the hypermatrices are taken to be functions
    so that while performing the product we compose with the entries
    of the first of the list of inputs. This implementation is aesthetically
    more pleasing then the previous one because it explicitly articulate the
    preference for composition as our defacto product operation. However it is
    theoretically less general then the previous one. Both these implementations
    are inspired by initial exposure to ideas from Category theory, the implementation
    also make painfully obvious some of the programming constraints imposed by Python.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:

    ::

        sage: x,y=var('x,y'); Ha=x*y*HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c')
        sage: Rslt=GProd([Ha,Hb,Hc], sum, [x,y]); Rslt.printHM()
        [:, :, 0]=
        [a000*b000*c000 + a010*b001*c100 a000*b010*c010 + a010*b011*c110]
        [a100*b100*c000 + a110*b101*c100 a100*b110*c010 + a110*b111*c110]
        <BLANKLINE>
        [:, :, 1]=
        [a001*b000*c001 + a011*b001*c101 a001*b010*c011 + a011*b011*c111]
        [a101*b100*c001 + a111*b101*c101 a101*b110*c011 + a111*b111*c111]
        <BLANKLINE>
        sage: x=var('x'); Ha=x*HM(2,2,'a'); Hb=HM(2,2,'b'); GProd([Ha, Hb], sum, [x])
        [[a00*b00 + a01*b10, a00*b01 + a01*b11], [a10*b00 + a11*b10, a10*b01 + a11*b11]]
        sage: Ha=HM(2,2,'a').elementwise_exponent(x); Hb=HM(2,2,'b'); GProd([Ha, Hb], prod, [x])
        [[a00^b00*a01^b10, a00^b01*a01^b11], [a10^b00*a11^b10, a10^b01*a11^b11]]
        sage: sz=3; l=2; Lx=var_list('x',l); La=HM(sz,l,'a').list(); Lb=HM(sz,l,'b').list(); Lc=var_list('c',sz); z=var('z')
        sage: F=FreeAlgebra(QQ,len(La+Lx+Lb+Lc+[z]),La+Lx+Lb+Lc+[z])
        sage: F.<a00, a10, a20, a01, a11, a21, x0, x1, b00, b10, b20, b01, b11, b21, c0, c1, c2, z>=FreeAlgebra(QQ,len(La+Lx+Lb+Lc+[z]))
        sage: Ha=HM(sz,l,[a00, a10, a20, a01, a11, a21]); Hx=HM(l,1,[x0, x1])
        sage: Hb=HM(sz,l,[b00, b10, b20, b01, b11, b21])
        sage: Hr=Ha.elementwise_product(z*Hb)-HM(sz,1,[c0, c1, c2])*HM(1,l,[QQ(1/2) for i in rg(l)])
        sage: GProd([Hr, Hx], sum, [z]).printHM()
        [:, :]=
        [-c0 + a00*x0*b00 + a01*x1*b01]
        [-c1 + a10*x0*b10 + a11*x1*b11]
        [-c2 + a20*x0*b20 + a21*x1*b21]


    AUTHORS:
    - Edinah K. Gnang
    """
    return GeneralHypermatrixProductIV(Lh, Op, Lv)

def GProdII(Lh, Op, Lv, indx):
    """
    Outputs an HM which is list of lists associated with the composition
    based Bhattacharya-Mesner product of the input hypermatrices.
    The entries of the hypermatrices are taken to be functions
    so that while performing the product we compose with the entries
    of the indx-th of the list of inputs. This implementation is aesthetically
    more pleasing then the previous one because it explicitly articulate the
    preference for a particular composition operation. However both this
    and the previous implementation are special instances of ways to specify the composer.
    Both these implementations are inspired by initial exposure to ideas from Category theory,
    the implementation also make painfully obvious some of the programming constraints imposed by Python.
    This code only handles the Hypermatrix HM class objects.

    EXAMPLES:

    ::

        sage: x,y=var('x,y'); Ha=x*y*HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c')
        sage: GProdII([Ha,Hb,Hc], sum, [x,y], 0).printHM()
        [:, :, 0]=
        [a000*b000*c000 + a010*b001*c100 a000*b010*c010 + a010*b011*c110]
        [a100*b100*c000 + a110*b101*c100 a100*b110*c010 + a110*b111*c110]
        <BLANKLINE>
        [:, :, 1]=
        [a001*b000*c001 + a011*b001*c101 a001*b010*c011 + a011*b011*c111]
        [a101*b100*c001 + a111*b101*c101 a101*b110*c011 + a111*b111*c111]
        <BLANKLINE>
        sage: x=var('x'); Ha=x*HM(2,2,'a'); Hb=HM(2,2,'b'); GProdII([Ha, Hb], sum, [x], 0)
        [[a00*b00 + a01*b10, a00*b01 + a01*b11], [a10*b00 + a11*b10, a10*b01 + a11*b11]]
        sage: Ha=HM(2,2,'a').elementwise_exponent(x); Hb=HM(2,2,'b'); GProdII([Ha, Hb], prod, [x], 0)
        [[a00^b00*a01^b10, a00^b01*a01^b11], [a10^b00*a11^b10, a10^b01*a11^b11]]
        sage: sz=3; l=2; Lx=var_list('x',l); La=HM(sz,l,'a').list(); Lb=HM(sz,l,'b').list(); Lc=var_list('c',sz); z=var('z')
        sage: F=FreeAlgebra(QQ,len(La+Lx+Lb+Lc+[z]),La+Lx+Lb+Lc+[z])
        sage: F.<a00, a10, a20, a01, a11, a21, x0, x1, b00, b10, b20, b01, b11, b21, c0, c1, c2, z>=FreeAlgebra(QQ,len(La+Lx+Lb+Lc+[z]))
        sage: Ha=HM(sz,l,[a00, a10, a20, a01, a11, a21]); Hx=HM(l,1,[x0, x1])
        sage: Hb=HM(sz,l,[b00, b10, b20, b01, b11, b21])
        sage: Hr=Ha.elementwise_product(z*Hb)-HM(sz,1,[c0, c1, c2])*HM(1,l,[QQ(1/2) for i in rg(l)])
        sage: GProdII([Hr, Hx], sum, [z], 0).printHM()
        [:, :]=
        [-c0 + a00*x0*b00 + a01*x1*b01]
        [-c1 + a10*x0*b10 + a11*x1*b11]
        [-c2 + a20*x0*b20 + a21*x1*b21]


    AUTHORS:
    - Edinah K. Gnang
    """
    return GeneralHypermatrixProductV(Lh, Op, Lv, indx)

def GProdIII(Lh, Op, F):
    """
    Outputs an HM which is a list of lists associated with the
    construct approach to the Bhattacharya-Mesner product of the input
    hypermatrices. Here Op is the combinator and F is the composer.
    This implementation in theory captures the full scope of
    the construct products subsequently implemented here but in 
    practice is hard for to specify arbitrary composers for F.
    The code only handles the Hypermatrix HM class objects.
    This implementation comes in handy for combinatorial constructs
    such as shortest path problems and set valued constructs.


    EXAMPLES:

    ::

        sage: Ha=HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c')
        sage: Rslt=GProdIII([Ha, Hb, Hc], prod, sum); Rslt.printHM()
        [:, :, 0]=
        [(a000 + b000 + c000)*(a010 + b001 + c100) (a000 + b010 + c010)*(a010 + b011 + c110)]
        [(a100 + b100 + c000)*(a110 + b101 + c100) (a100 + b110 + c010)*(a110 + b111 + c110)]
        <BLANKLINE>
        [:, :, 1]=
        [(a001 + b000 + c001)*(a011 + b001 + c101) (a001 + b010 + c011)*(a011 + b011 + c111)]
        [(a101 + b100 + c001)*(a111 + b101 + c101) (a101 + b110 + c011)*(a111 + b111 + c111)]        
        <BLANKLINE>
        sage: Rslt=GProdIII([Ha, Hb, Hc], sum, prod); Rslt.printHM()
        [:, :, 0]=
        [a000*b000*c000 + a010*b001*c100 a000*b010*c010 + a010*b011*c110]
        [a100*b100*c000 + a110*b101*c100 a100*b110*c010 + a110*b111*c110]
        <BLANKLINE>
        [:, :, 1]=
        [a001*b000*c001 + a011*b001*c101 a001*b010*c011 + a011*b011*c111]
        [a101*b100*c001 + a111*b101*c101 a101*b110*c011 + a111*b111*c111]
        <BLANKLINE>
        sage: Ha=HM(2,2,'a'); Hb=HM(2,2,'b'); GProdIII([Ha, Hb], sum, prod)
        [[a00*b00 + a01*b10, a00*b01 + a01*b11], [a10*b00 + a11*b10, a10*b01 + a11*b11]]
        sage: Ha=HM(2,2,'a'); Hb=HM(2,2,'b'); GProdIII([Ha, Hb], prod, Exp)
        [[a00^b00*a01^b10, a00^b01*a01^b11], [a10^b00*a11^b10, a10^b01*a11^b11]]
        sage: Ha=HM(2,2,'a'); Hb=HM(2,2,'b'); GProdIII([Ha, Hb], prod, BaseExp)
        [[b00^a00*b10^a01, b01^a00*b11^a01], [b00^a10*b10^a11, b01^a10*b11^a11]]
        sage: A=HM([[59, -3, 2], [-6, 1, 1], [1, -1, 1]]); B=HM([[-1, 1, 0], [-1, 1, 0], [0, -43, 1]])
        sage: GProdIII([A, B], sum, prod)-A*B # One way of recovering matrix multiplication
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        sage: GProdIII([A, B], min, sum) # Recovering the Min-plus matrix multiplication
        [[-4, -41, -3], [-7, -42, -6], [-2, -42, -1]]
        sage: MSta=HM([[Set([1,2]), Set([1,3,2])], [Set([1]), Set([2,3])]])
        sage: MStb=HM([[Set([1,2,3]), Set([2])], [Set([1,3]), Set([1,3])]])
        sage: GProdIII([MSta, MStb], SetIntersection, SetUnion) 
        [[{1, 2, 3}, {1, 2}], [{1, 2, 3}, {1, 2}]]
        sage: GProdIII([MSta, MStb], SetUnion, SetIntersection)
        [[{1, 2, 3}, {1, 2, 3}], [{1, 3}, {3}]]
        sage: sz=3; l=2; Lx=var_list('x',l); La=HM(sz,l,'a').list(); Lb=HM(sz,l,'b').list(); Lc=var_list('c',sz); z=var('z')
        sage: F=FreeAlgebra(QQ,len(La+Lx+Lb+Lc+[z]),La+Lx+Lb+Lc+[z])
        sage: F.<a00, a10, a20, a01, a11, a21, x0, x1, b00, b10, b20, b01, b11, b21, c0, c1, c2, z>=FreeAlgebra(QQ,len(La+Lx+Lb+Lc+[z]))
        sage: Ha=HM(1,l,sz,[a00, a10, a20, a01, a11, a21])
        sage: Hx=HM(1,1,l,[x0, x1]) 
        sage: Hb=HM(l,1,sz,[b00, b10, b20, b01, b11, b21])
        sage: GProdIII([Ha, Hx, Hb], sum, prod).printHM()
        [:, :, 0]=
        [a00*x0*b00 + a10*x1*b10]
        <BLANKLINE>
        [:, :, 1]=
        [a20*x0*b20 + a01*x1*b01]
        <BLANKLINE>
        [:, :, 2]=
        [a11*x0*b11 + a21*x1*b21]
        <BLANKLINE>
        sage: Ca=HM(2,2,[HM(1,1,[var('a00')]), HM(1,1,[var('a10')]), HM(1,1,[var('a01')]), HM(1,1,[var('a11')])])
        sage: Cb=HM(2,2,[HM(1,1,[var('b00')]), HM(1,1,[var('b10')]), HM(1,1,[var('b01')]), HM(1,1,[var('b11')])])
        sage: GProdIII([Ca, Cb], DirectSum, TensorProduct)[0,0].printHM()
        [:, :]=
        [a00*b00       0]
        [      0 a01*b10]
        sage: GProdIII([Ca, Cb], DirectSum, TensorProduct)[0,1].printHM()
        [:, :]=
        [a00*b01       0]
        [      0 a01*b11]
        sage: GProdIII([Ca, Cb], DirectSum, TensorProduct)[1,0].printHM()
        [:, :]=
        [a10*b00       0]
        [      0 a11*b10]
        sage: GProdIII([Ca, Cb], DirectSum, TensorProduct)[1,1].printHM()
        [:, :]=
        [a10*b01       0]
        [      0 a11*b11]


    AUTHORS:
    - Edinah K. Gnang
    """
    return GeneralHypermatrixProductII(Lh, Op, F)

def CProd(Lh, Op, F):
    """
    Outputs an HM which is a list of lists associated with the
    construct approach to the Bhattacharya-Mesner product of the input
    hypermatrices. Here Op is the combinator and F is the composer.
    This implementation in theory captures the full scope of
    the construct products subsequently implemented here but in 
    practice is hard for to specify arbitrary composers for F.
    The code only handles the Hypermatrix HM class objects.
    This implementation comes in handy for combinatorial constructs
    such as shortest path problems and set valued constructs.


    EXAMPLES:

    ::

        sage: Ha=HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c')
        sage: Rslt=CProd([Ha, Hb, Hc], prod, sum); Rslt.printHM()
        [:, :, 0]=
        [(a000 + b000 + c000)*(a010 + b001 + c100) (a000 + b010 + c010)*(a010 + b011 + c110)]
        [(a100 + b100 + c000)*(a110 + b101 + c100) (a100 + b110 + c010)*(a110 + b111 + c110)]
        <BLANKLINE>
        [:, :, 1]=
        [(a001 + b000 + c001)*(a011 + b001 + c101) (a001 + b010 + c011)*(a011 + b011 + c111)]
        [(a101 + b100 + c001)*(a111 + b101 + c101) (a101 + b110 + c011)*(a111 + b111 + c111)]        
        <BLANKLINE>
        sage: Rslt=CProd([Ha, Hb, Hc], sum, prod); Rslt.printHM()
        [:, :, 0]=
        [a000*b000*c000 + a010*b001*c100 a000*b010*c010 + a010*b011*c110]
        [a100*b100*c000 + a110*b101*c100 a100*b110*c010 + a110*b111*c110]
        <BLANKLINE>
        [:, :, 1]=
        [a001*b000*c001 + a011*b001*c101 a001*b010*c011 + a011*b011*c111]
        [a101*b100*c001 + a111*b101*c101 a101*b110*c011 + a111*b111*c111]
        <BLANKLINE>
        sage: Ha=HM(2,2,'a'); Hb=HM(2,2,'b'); CProd([Ha, Hb], sum, prod)
        [[a00*b00 + a01*b10, a00*b01 + a01*b11], [a10*b00 + a11*b10, a10*b01 + a11*b11]]
        sage: Ha=HM(2,2,'a'); Hb=HM(2,2,'b'); CProd([Ha, Hb], prod, Exp)
        [[a00^b00*a01^b10, a00^b01*a01^b11], [a10^b00*a11^b10, a10^b01*a11^b11]]
        sage: Ha=HM(2,2,'a'); Hb=HM(2,2,'b'); CProd([Ha, Hb], prod, BaseExp)
        [[b00^a00*b10^a01, b01^a00*b11^a01], [b00^a10*b10^a11, b01^a10*b11^a11]]
        sage: A=HM([[59, -3, 2], [-6, 1, 1], [1, -1, 1]]); B=HM([[-1, 1, 0], [-1, 1, 0], [0, -43, 1]])
        sage: CProd([A, B], sum, prod)-A*B # One way of recovering matrix multiplication
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        sage: CProd([A, B], min, sum) # Recovering the Min-plus matrix multiplication
        [[-4, -41, -3], [-7, -42, -6], [-2, -42, -1]]
        sage: MSta=HM([[Set([1,2]), Set([1,3,2])], [Set([1]), Set([2,3])]])
        sage: MStb=HM([[Set([1,2,3]), Set([2])], [Set([1,3]), Set([1,3])]])
        sage: CProd([MSta, MStb], SetIntersection, SetUnion) 
        [[{1, 2, 3}, {1, 2}], [{1, 2, 3}, {1, 2}]]
        sage: CProd([MSta, MStb], SetUnion, SetIntersection)
        [[{1, 2, 3}, {1, 2, 3}], [{1, 3}, {3}]]
        sage: sz=3; l=2; Lx=var_list('x',l); La=HM(sz,l,'a').list(); Lb=HM(sz,l,'b').list(); Lc=var_list('c',sz); z=var('z')
        sage: F=FreeAlgebra(QQ,len(La+Lx+Lb+Lc+[z]),La+Lx+Lb+Lc+[z])
        sage: F.<a00, a10, a20, a01, a11, a21, x0, x1, b00, b10, b20, b01, b11, b21, c0, c1, c2, z>=FreeAlgebra(QQ,len(La+Lx+Lb+Lc+[z]))
        sage: Ha=HM(1,l,sz,[a00, a10, a20, a01, a11, a21])
        sage: Hx=HM(1,1,l,[x0, x1]) 
        sage: Hb=HM(l,1,sz,[b00, b10, b20, b01, b11, b21])
        sage: CProd([Ha, Hx, Hb], sum, prod).printHM()
        [:, :, 0]=
        [a00*x0*b00 + a10*x1*b10]
        <BLANKLINE>
        [:, :, 1]=
        [a20*x0*b20 + a01*x1*b01]
        <BLANKLINE>
        [:, :, 2]=
        [a11*x0*b11 + a21*x1*b21]
        <BLANKLINE>
        sage: Ca=HM(2,2,[HM(1,1,[var('a00')]), HM(1,1,[var('a10')]), HM(1,1,[var('a01')]), HM(1,1,[var('a11')])])
        sage: Cb=HM(2,2,[HM(1,1,[var('b00')]), HM(1,1,[var('b10')]), HM(1,1,[var('b01')]), HM(1,1,[var('b11')])])
        sage: CProd([Ca, Cb], DirectSum, TensorProduct)[0,0].printHM()
        [:, :]=
        [a00*b00       0]
        [      0 a01*b10]
        sage: CProd([Ca, Cb], DirectSum, TensorProduct)[0,1].printHM()
        [:, :]=
        [a00*b01       0]
        [      0 a01*b11]
        sage: CProd([Ca, Cb], DirectSum, TensorProduct)[1,0].printHM()
        [:, :]=
        [a10*b00       0]
        [      0 a11*b10]
        sage: CProd([Ca, Cb], DirectSum, TensorProduct)[1,1].printHM()
        [:, :]=
        [a10*b01       0]
        [      0 a11*b11]


    AUTHORS:
    - Edinah K. Gnang
    """
    return GeneralHypermatrixProductII(Lh, Op, F)

def CProdII(Lh, Op, Lv):
    """
    Outputs a list of lists associated with the composition
    based Bhattacharya-Mesner product of the input hypermatrices.
    The entries of the hypermatrices are taken to be functions
    so that while performing the product we compose with the entries
    of the first of the list of inputs. This implementation is aesthetically
    more pleasing then the previous one because it explicitly articulate the
    preference for composition as our defacto product operation. However it is
    theoretically less general then the previous one. Both these implementations
    are inspired by initial exposure to ideas from Category theory, the implementation
    also make painfully obvious some of the programming constraints imposed by Python.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:

    ::

        sage: x,y=var('x,y'); Ha=x*y*HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c')
        sage: Rslt=CProdII([Ha,Hb,Hc], sum, [x,y]); Rslt.printHM()
        [:, :, 0]=
        [a000*b000*c000 + a010*b001*c100 a000*b010*c010 + a010*b011*c110]
        [a100*b100*c000 + a110*b101*c100 a100*b110*c010 + a110*b111*c110]
        <BLANKLINE>
        [:, :, 1]=
        [a001*b000*c001 + a011*b001*c101 a001*b010*c011 + a011*b011*c111]
        [a101*b100*c001 + a111*b101*c101 a101*b110*c011 + a111*b111*c111]
        <BLANKLINE>
        sage: x=var('x'); Ha=x*HM(2,2,'a'); Hb=HM(2,2,'b'); CProdII([Ha, Hb], sum, [x])
        [[a00*b00 + a01*b10, a00*b01 + a01*b11], [a10*b00 + a11*b10, a10*b01 + a11*b11]]
        sage: Ha=HM(2,2,'a').elementwise_exponent(x); Hb=HM(2,2,'b'); CProdII([Ha, Hb], prod, [x])
        [[a00^b00*a01^b10, a00^b01*a01^b11], [a10^b00*a11^b10, a10^b01*a11^b11]]
        sage: sz=3; l=2; Lx=var_list('x',l); La=HM(sz,l,'a').list(); Lb=HM(sz,l,'b').list(); Lc=var_list('c',sz); z=var('z')
        sage: F=FreeAlgebra(QQ,len(La+Lx+Lb+Lc+[z]),La+Lx+Lb+Lc+[z])
        sage: F.<a00, a10, a20, a01, a11, a21, x0, x1, b00, b10, b20, b01, b11, b21, c0, c1, c2, z>=FreeAlgebra(QQ,len(La+Lx+Lb+Lc+[z]))
        sage: Ha=HM(sz,l,[a00, a10, a20, a01, a11, a21]); Hx=HM(l,1,[x0, x1])
        sage: Hb=HM(sz,l,[b00, b10, b20, b01, b11, b21])
        sage: Hr=Ha.elementwise_product(z*Hb)-HM(sz,1,[c0, c1, c2])*HM(1,l,[QQ(1/2) for i in rg(l)])
        sage: CProdII([Hr, Hx], sum, [z]).printHM()
        [:, :]=
        [-c0 + a00*x0*b00 + a01*x1*b01]
        [-c1 + a10*x0*b10 + a11*x1*b11]
        [-c2 + a20*x0*b20 + a21*x1*b21]


    AUTHORS:
    - Edinah K. Gnang
    """
    return GeneralHypermatrixProductIV(Lh, Op, Lv)

def GProdII(Lh, Op, Lv, indx):
    """
    Outputs an HM which is list of lists associated with the composition
    based Bhattacharya-Mesner product of the input hypermatrices.
    The entries of the hypermatrices are taken to be functions
    so that while performing the product we compose with the entries
    of the indx-th of the list of inputs. This implementation is aesthetically
    more pleasing then the previous one because it explicitly articulate the
    preference for a particular composition operation. However both this
    and the previous implementation are special instances of ways to specify the composer.
    Both these implementations are inspired by initial exposure to ideas from Category theory,
    the implementation also make painfully obvious some of the programming constraints imposed by Python.
    This code only handles the Hypermatrix HM class objects.

    EXAMPLES:

    ::

        sage: x,y=var('x,y'); Ha=x*y*HM(2,2,2,'a'); Hb=HM(2,2,2,'b'); Hc=HM(2,2,2,'c')
        sage: GProdII([Ha,Hb,Hc], sum, [x,y], 0).printHM()
        [:, :, 0]=
        [a000*b000*c000 + a010*b001*c100 a000*b010*c010 + a010*b011*c110]
        [a100*b100*c000 + a110*b101*c100 a100*b110*c010 + a110*b111*c110]
        <BLANKLINE>
        [:, :, 1]=
        [a001*b000*c001 + a011*b001*c101 a001*b010*c011 + a011*b011*c111]
        [a101*b100*c001 + a111*b101*c101 a101*b110*c011 + a111*b111*c111]
        <BLANKLINE>
        sage: x=var('x'); Ha=x*HM(2,2,'a'); Hb=HM(2,2,'b'); GProdII([Ha, Hb], sum, [x], 0)
        [[a00*b00 + a01*b10, a00*b01 + a01*b11], [a10*b00 + a11*b10, a10*b01 + a11*b11]]
        sage: Ha=HM(2,2,'a').elementwise_exponent(x); Hb=HM(2,2,'b'); GProdII([Ha, Hb], prod, [x], 0)
        [[a00^b00*a01^b10, a00^b01*a01^b11], [a10^b00*a11^b10, a10^b01*a11^b11]]
        sage: sz=3; l=2; Lx=var_list('x',l); La=HM(sz,l,'a').list(); Lb=HM(sz,l,'b').list(); Lc=var_list('c',sz); z=var('z')
        sage: F=FreeAlgebra(QQ,len(La+Lx+Lb+Lc+[z]),La+Lx+Lb+Lc+[z])
        sage: F.<a00, a10, a20, a01, a11, a21, x0, x1, b00, b10, b20, b01, b11, b21, c0, c1, c2, z>=FreeAlgebra(QQ,len(La+Lx+Lb+Lc+[z]))
        sage: Ha=HM(sz,l,[a00, a10, a20, a01, a11, a21]); Hx=HM(l,1,[x0, x1])
        sage: Hb=HM(sz,l,[b00, b10, b20, b01, b11, b21])
        sage: Hr=Ha.elementwise_product(z*Hb)-HM(sz,1,[c0, c1, c2])*HM(1,l,[QQ(1/2) for i in rg(l)])
        sage: GProdII([Hr, Hx], sum, [z], 0).printHM()
        [:, :]=
        [-c0 + a00*x0*b00 + a01*x1*b01]
        [-c1 + a10*x0*b10 + a11*x1*b11]
        [-c2 + a20*x0*b20 + a21*x1*b21]


    AUTHORS:
    - Edinah K. Gnang
    """
    return GeneralHypermatrixProductV(Lh, Op, Lv, indx)


def GProdB(Lh, Op, F):
    """
    Outputs an HM whose entries are themselves hypermatrices
    providing a construct approach to the general Bhattacharya-Mesner
    product with non-trivial background hypermatrix. This code emphasizes
    the outer product picture. The function is currently implemented
    having the sum as the combinator input noted  Op  and the composer
    input noted F is the general BM product.
    Note that one must be careful here with the combinator because of
    the conformablity issue. The sum combinator is currently the safe
    bet.
    This implementation in theory captures the full scope of
    the construct products subsequently implemented here but in 
    practice is hard for to specify arbitrary composers for F.
    The code only handles the Hypermatrix HM class objects.


    EXAMPLES:

    ::

        sage: AlphaB=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
        sage: sz=2; A=HM(sz,1,[HM(1,sz,var_list(AlphaB[i],sz)) for i in rg(sz)]); B=HM(1,sz,[HM(sz,1,var_list(AlphaB[sz+j],sz)) for j in rg(sz)])
        sage: C=HM(sz,sz,AlphaB[2*sz]) # Initialization of the Background hypermatrix
        sage: GProdB([A, B, C], sum, ProdB)
        [[[[a0*c0*e00 + a0*c1*e01 + a1*c0*e10 + a1*c1*e11]], [[a0*d0*e00 + a0*d1*e01 + a1*d0*e10 + a1*d1*e11]]], [[[b0*c0*e00 + b0*c1*e01 + b1*c0*e10 + b1*c1*e11]], [[b0*d0*e00 + b0*d1*e01 + b1*d0*e10 + b1*d1*e11]]]]
        sage: Hr=GProdB([A, B, C], sum, ProdB); D=HM(sz,sz,[Hr[i,j][0,0] for j in rg(sz) for i in rg(sz)]); D.printHM()
        [:, :]=
        [a0*c0*e00 + a0*c1*e01 + a1*c0*e10 + a1*c1*e11 a0*d0*e00 + a0*d1*e01 + a1*d0*e10 + a1*d1*e11]
        [b0*c0*e00 + b0*c1*e01 + b1*c0*e10 + b1*c1*e11 b0*d0*e00 + b0*d1*e01 + b1*d0*e10 + b1*d1*e11]
        sage: sz=2; A=HM(sz,1,[HM(1,sz,var_list(AlphaB[i],sz)) for i in rg(sz)]); B=HM(1,sz,[HM(sz,1,var_list(AlphaB[sz+j],sz)) for j in rg(sz)])
        sage: Hr=GProdB([A, B, HM(2,sz,'kronecker')], sum, ProdB); D=HM(sz,sz,[Hr[i,j][0,0] for j in rg(sz) for i in rg(sz)]); D.printHM()
        [:, :]=
        [a0*c0 + a1*c1 a0*d0 + a1*d1]
        [b0*c0 + b1*c1 b0*d0 + b1*d1]
        sage: sz=2 # Initialization of the size parameter
        sage: A=HM(sz,1,sz,[HM(1,sz,1,var_list(AlphaB[0*sz^2+sz*i+j],sz)) for j in rg(sz) for i in rg(sz)])
        sage: B=HM(sz,sz,1,[HM(1,1,sz,var_list(AlphaB[1*sz^2+sz*i+j],sz)) for j in rg(sz) for i in rg(sz)])
        sage: C=HM(1,sz,sz,[HM(sz,1,1,var_list(AlphaB[2*sz^3+sz*i+j],sz)) for j in rg(sz) for i in rg(sz)])
        sage: D=HM(sz,sz,sz,'d') # Initialization of the background matrix
        sage: Ht=GProdB([A, B, C, D], sum, ProdB); E=HM(sz,sz,sz,[Ht[i,j,k][0,0,0] for k in rg(sz) for j in rg(sz) for i in rg(sz)]); E.printHM()
        [:, :, 0]=
        [a0*d000*e0*q0 + a1*d100*e0*q0 + a0*d010*e1*q0 + a1*d110*e1*q0 + a0*d001*e0*q1 + a1*d101*e0*q1 + a0*d011*e1*q1 + a1*d111*e1*q1 a0*d000*f0*s0 + a1*d100*f0*s0 + a0*d010*f1*s0 + a1*d110*f1*s0 + a0*d001*f0*s1 + a1*d101*f0*s1 + a0*d011*f1*s1 + a1*d111*f1*s1]
        [c0*d000*g0*q0 + c1*d100*g0*q0 + c0*d010*g1*q0 + c1*d110*g1*q0 + c0*d001*g0*q1 + c1*d101*g0*q1 + c0*d011*g1*q1 + c1*d111*g1*q1 c0*d000*h0*s0 + c1*d100*h0*s0 + c0*d010*h1*s0 + c1*d110*h1*s0 + c0*d001*h0*s1 + c1*d101*h0*s1 + c0*d011*h1*s1 + c1*d111*h1*s1]
        <BLANKLINE>
        [:, :, 1]=
        [b0*d000*e0*r0 + b1*d100*e0*r0 + b0*d010*e1*r0 + b1*d110*e1*r0 + b0*d001*e0*r1 + b1*d101*e0*r1 + b0*d011*e1*r1 + b1*d111*e1*r1 b0*d000*f0*t0 + b1*d100*f0*t0 + b0*d010*f1*t0 + b1*d110*f1*t0 + b0*d001*f0*t1 + b1*d101*f0*t1 + b0*d011*f1*t1 + b1*d111*f1*t1]
        [d0*d000*g0*r0 + d1*d100*g0*r0 + d0*d010*g1*r0 + d1*d110*g1*r0 + d0*d001*g0*r1 + d1*d101*g0*r1 + d0*d011*g1*r1 + d1*d111*g1*r1 d0*d000*h0*t0 + d1*d100*h0*t0 + d0*d010*h1*t0 + d1*d110*h1*t0 + d0*d001*h0*t1 + d1*d101*h0*t1 + d0*d011*h1*t1 + d1*d111*h1*t1]
        <BLANKLINE>
        sage: AlphaB=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
        sage: sz=2 # Initialization of the size parameter
        sage: A=HM(sz,1,sz,[HM(1,sz,1,var_list(AlphaB[0*sz^2+sz*i+j],sz)) for j in rg(sz) for i in rg(sz)])
        sage: B=HM(sz,sz,1,[HM(1,1,sz,var_list(AlphaB[1*sz^2+sz*i+j],sz)) for j in rg(sz) for i in rg(sz)])
        sage: C=HM(1,sz,sz,[HM(sz,1,1,var_list(AlphaB[2*sz^3+sz*i+j],sz)) for j in rg(sz) for i in rg(sz)])
        sage: Ht=GProdB([A, B, C, HM(3, sz, 'kronecker')], sum, ProdB); E=HM(sz,sz,sz,[Ht[i,j,k][0,0,0] for k in rg(sz) for j in rg(sz) for i in rg(sz)]); E.printHM()
        [:, :, 0]=
        [a0*e0*q0 + a1*e1*q1 a0*f0*s0 + a1*f1*s1]
        [c0*g0*q0 + c1*g1*q1 c0*h0*s0 + c1*h1*s1]
        <BLANKLINE>
        [:, :, 1]=
        [b0*e0*r0 + b1*e1*r1 b0*f0*t0 + b1*f1*t1]
        [d0*g0*r0 + d1*g1*r1 d0*h0*t0 + d1*h1*t1]
        <BLANKLINE>


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list specifying the dimensions of the output
    l = [(Lh[i]).n(i) for i in range(len(Lh)-1)]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the assignement
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # computing the Hypermatrix product
        if len(Lh)<2:
            raise ValueError("The number of operands must be >= 2")
        elif len(Lh) >= 2:
            t=0
            #Rh[tuple(entry)]=Op(*[[F(*[ [Lh[s][tuple(entry[0:Integer(mod(s+1,len(l)))]+[t]+entry[Integer(mod(s+2,len(l))):])] for s in range(len(l)-2)]+[Lh[len(l)-2][tuple(entry[0:len(l)-1]+[t])]]+[Lh[len(l)-1][tuple([t]+entry[1:])]]+[Lh[len(Lh)-1]] ] ) for t in range((Lh[0]).n(1))]])
            LsT=[[F( *([Lh[s][tuple(entry[0:Integer(mod(s+1,len(l)))]+[t]+entry[Integer(mod(s+2,len(l))):])] for s in range(len(l)-2)]+[Lh[len(l)-2][tuple(entry[0:len(l)-1]+[t])]]+[Lh[len(l)-1][tuple([t]+entry[1:])]]+[Lh[len(Lh)-1]]) ) for t in range((Lh[0]).n(1))]]
            Rh[tuple(entry)]=Op(*LsT)
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
        entry = [Integer(mod(i,l[0]))]
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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=s*A[tuple(entry)]
    return Rh

def GeneralHypermatrixScaleRight(A,s):
    """
    Outputs a list of lists associated with the scaling of a general hypermatrix.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:

    ::

        sage: Ha=HM(2,2,2,'a')
        sage: GeneralHypermatrixScaleRight(Ha,3)
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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=A[tuple(entry)]*s
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
        entry = [Integer(mod(i,l[0]))]
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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=s^(A[tuple(entry)])
    return Rh

def GeneralHypermatrixBaseExponentN(A,s, dgts=50):
    """
    Outputs a list of lists associated with the general
    whose entries are exponentiated using the input s as
    basis for the numerical exponentiation.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:

    ::

        sage: Ha=HM(2,2,[1,2,3,4]); GeneralHypermatrixBaseExponentN(Ha,3)
        [[3.0000000000000, 27.000000000000], [9.0000000000000, 81.000000000000]]
        

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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=ComplexField(dgts)(s)^ComplexField(dgts)(A[tuple(entry)])
    return Rh

def GeneralHypermatrixLogarithm(A,s=e):
    """
    Outputs a list of lists associated with the general
    whose entries are logarithms to the base s of the 
    original hypermatrix.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:

    ::

        sage: Ha=HM(2,2,2,'a'); GeneralHypermatrixLogarithm(Ha,3)
        [[[log(a000)/log(3), log(a001)/log(3)], [log(a010)/log(3), log(a011)/log(3)]], [[log(a100)/log(3), log(a101)/log(3)], [log(a110)/log(3), log(a111)/log(3)]]]


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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=log(A[tuple(entry)],s).canonicalize_radical()
    return Rh

def GeneralHypermatrixLogarithmN(A,s=e, dgts=50):
    """
    Outputs a list of lists associated with the general
    whose entries are numerical logarithms to the base s of the 
    original hypermatrix.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:

    ::

        sage: Ha=HM(2,2,[1,2,3,4]); GeneralHypermatrixLogarithmN(Ha,3)
        [[0.00000000000000, 1.0000000000000], [0.63092975357146, 1.2618595071429]]


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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=ComplexField(dgts)(log(A[tuple(entry)],s))
    return Rh

def GeneralHypermatrixApplyMap(A, phi):
    """
    Apply the given map phi (an arbitrary Python function or callable
    object) to this hypermatrix.
    

    INPUT:

    -  ``phi`` - arbitrary Python function or callable object


    OUTPUT: a symbolic hypermatrix

    EXAMPLES::

        sage: A = HM(2,2,'a')
        sage: phi = lambda x: sin(x)
        sage: GeneralHypermatrixApplyMap(A, phi).matrix()
        [sin(a00) sin(a01)]
        [sin(a10) sin(a11)]

    """
    # Initialization of the list specifying the dimensions of the output
    l = [A.n(i) for i in range(A.order())]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the computations of the entries
    for i in range(prod(l)):
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        if A[tuple(entry)].is_zero():
            Rh[tuple(entry)] = 0
        else:
            Rh[tuple(entry)] = phi(A[tuple(entry)])
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
        entry = [Integer(mod(i,l[0]))]
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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=(A[tuple(entry)]).expand()
    return Rh

def GeneralHypermatrixFactor(A):
    """
    Outputs a list of lists associated with the general
    hypermatrix with expressions in the entries in their
    factored form.
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:

    ::

        sage: sz=2; vx=HM(sz,'x').list(); vy=HM(sz,'y').list(); vz=HM(sz,'z').list()
        sage: X=HM(sz, 1, sz, [vx[u] for v in range(sz) for t in range(1) for u in range(sz)])
        sage: Y=HM(sz, 1, sz, [vy[u] for v in range(sz) for t in range(1) for u in range(sz)])
        sage: Z=HM(sz, 1, sz, [vz[u] for v in range(sz) for t in range(1) for u in range(sz)])
        sage: A=Prod(X, Y.transpose(2), Z.transpose())
        sage: GeneralHypermatrixFactor(Prod(A,A,A))
        [[[(x0*y0*z0 + x1*y1*z1)*x0^2*y0^2*z0^2, (x0*y0*z0 + x1*y1*z1)*x0^2*y0^2*z1^2], [(x0*y0*z0 + x1*y1*z1)*x0^2*y1^2*z0^2, (x0*y0*z0 + x1*y1*z1)*x0^2*y1^2*z1^2]], [[(x0*y0*z0 + x1*y1*z1)*x1^2*y0^2*z0^2, (x0*y0*z0 + x1*y1*z1)*x1^2*y0^2*z1^2], [(x0*y0*z0 + x1*y1*z1)*x1^2*y1^2*z0^2, (x0*y0*z0 + x1*y1*z1)*x1^2*y1^2*z1^2]]]

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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=(A[tuple(entry)]).factor()
    return Rh

def GeneralHypermatrixSimplifyFull(A):
    """
    Performs the symbolic simplification of the expressions
    associated with the hypermatrix entries. 

    EXAMPLES:

    ::

        sage: x,y=var('x,y'); ((x+y)^2*HM(2,2,2,'one')).simplify_full()
        [[[x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2], [x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2]], [[x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2], [x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2]]]
 

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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        if type(Integer(0)) != type(A[tuple(entry)]):
            Rh[tuple(entry)]=(A[tuple(entry)]).simplify_full()
        else:
            Rh[tuple(entry)]=A[tuple(entry)]
    return Rh

def GeneralHypermatrixSimplify(A):
    """
    Performs the symbolic simplification of the expressions
    associated with the hypermatrix entries. 

    EXAMPLES:

    ::

        sage: x,y=var('x,y'); GeneralHypermatrixSimplify((x+y)^2*HM(2,2,2,'one'))
        [[[(x + y)^2, (x + y)^2], [(x + y)^2, (x + y)^2]], [[(x + y)^2, (x + y)^2], [(x + y)^2, (x + y)^2]]]
 

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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        if type(Integer(0)) != type(A[tuple(entry)]):
            Rh[tuple(entry)]=(A[tuple(entry)]).simplify()
        else:
            Rh[tuple(entry)]=A[tuple(entry)]
    return Rh

def GeneralHypermatrixCanonicalizeRadical(A):
    """
    Performs the symbolic simplification of the expressions
    associated with the hypermatrix entries. 

    EXAMPLES:

    ::

        sage: x,y = var('x,y') 
        sage: GeneralHypermatrixCanonicalizeRadical((x+y)^2*HM(2,2,2,'one'))
        [[[x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2], [x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2]], [[x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2], [x^2 + 2*x*y + y^2, x^2 + 2*x*y + y^2]]]

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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        if type(Integer(0)) != type(A[tuple(entry)]):
            Rh[tuple(entry)]=(A[tuple(entry)]).canonicalize_radical()
        else:
            Rh[tuple(entry)]=A[tuple(entry)]
    return Rh

def GeneralHypermatrixNumerical(A, dgts=15):
    """
    Performs the symbolic simplification of the expressions
    associated with the hypermatrix entries. 


    EXAMPLES:

    ::

        sage: Ha=HM(3,3,[exp(I*2*pi*u*v/3) for v in range(3) for u in range(3)]) 
        sage: GeneralHypermatrixNumerical(Ha)
        [[1.00000000000000, 1.00000000000000, 1.00000000000000], [1.00000000000000, -0.500000000000000 + 0.866025403784439*I, -0.500000000000000 - 0.866025403784439*I], [1.00000000000000, -0.500000000000000 - 0.866025403784439*I, -0.500000000000000 + 0.866025403784439*I]] 
         

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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=N(A[tuple(entry)], digits=dgts)
    return Rh

def GeneralHypermatrixSubstitute(A, Dct):
    """
    Procedure for computes the substitution in the Hypermatrix entries
    the inputs are the corresponding Hypermatric and a dictionary 
    datastructure.

    EXAMPLES:

    ::

        sage: GeneralHypermatrixSubstitute(HM(3,2,'a','sym'), dict([(var('a011'),var('x')),(var('a001'),var('y')),(var('a000'),var('z')),(var('a111'),var('t'))]))
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
        entry = [Integer(mod(i,l[0]))]
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

        sage: GeneralHypermatrixSubstituteN(HM(3,2,'a','sym'),dict([(var('a011'),1),(var('a001'),2),(var('a000'),3),(var('a111'),4)]))
        [[[3.00000000000000, 2.00000000000000], [2.00000000000000, 1.00000000000000]], [[2.00000000000000, 1.00000000000000], [1.00000000000000, 4.00000000000000]]]
        


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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=N((A[tuple(entry)]).subs(Dct))
    return Rh

def GeneralHypermatrixSubstituteII(A, *args, **kwds):
    """
    Procedure for computes the substitution in the Hypermatrix entries
    the inputs are the corresponding Hypermatric and a dictionary 
    datastructure.

    EXAMPLES:

    ::

        sage: a,b,d,e=var('a,b,d,e'); M=HM([[a,b],[d,e]])
        sage: GeneralHypermatrixSubstituteII(M,a=1)
        [[1, b], [d, e]]
        sage: GeneralHypermatrixSubstituteII(M, a=b, b=d)
        [[b, d], [d, e]]
        sage: GeneralHypermatrixSubstituteII(M, {a: 3, b:2, d:1, e:-1})
        [[3, 2], [1, -1]]
        

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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=(A[tuple(entry)]).subs(*args, **kwds)
    return Rh

def GeneralHypermatrixCopy(A):
    """
    Procedure for performing a Hypermatrix copy.

    EXAMPLES:

    ::

        sage: A=HM(2,2,2,'a'); B=HM(2,2,2,'b'); GeneralHypermatrixCopy(A+B)
        [[[a000 + b000, a001 + b001], [a010 + b010, a011 + b011]], [[a100 + b100, a101 + b101], [a110 + b110, a111 + b111]]]

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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=A[tuple(entry)]
    return Rh

def GeneralHypermatrixAppendIndex(A,indx):
    """
    Procedure for computing Hypermatrix Hadamard addition.

    EXAMPLES:

    ::

        sage: GeneralHypermatrixAppendIndex(HM(2,2,'a','shift'),1)
        [[a111, a121], [a211, a221]]

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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=var(str(A[tuple(entry)])+str(indx))
    return Rh

def List2Hypermatrix(*args):
    """
    Procedure for Initializing a Hypermatrix from a size specifications and a list

    EXAMPLES:

    ::

        sage: A = List2Hypermatrix(2,2,2,HM(2,2,2,'a').list()); A
        [[[a000, a001], [a010, a011]], [[a100, a101], [a110, a111]]]

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
        entry = [Integer(mod(i,l[0]))]
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

        sage: A=HM(2,2,2,'a'); B=HM(2,2,2,'b'); GeneralHypermatrixAdd(A,B)
        [[[a000 + b000, a001 + b001], [a010 + b010, a011 + b011]], [[a100 + b100, a101 + b101], [a110 + b110, a111 + b111]]]


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
            entry = [Integer(mod(i,l[0]))]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            Rh[tuple(entry)]=A[tuple(entry)]+B[tuple(entry)]
        return Rh
    else:
        raise ValueError("The Dimensions of the input hypermatrices must match.")
 
def GeneralHypermatrixHadamardProduct(A,B):
    """
    Procedure for computing Hypermatrix Hadamard products.

    EXAMPLES:

    ::

        sage: A=HM(2,2,2,'a'); B=HM(2,2,2,'b')
        sage: GeneralHypermatrixHadamardProduct(A,B)
        [[[a000*b000, a001*b001], [a010*b010, a011*b011]], [[a100*b100, a101*b101], [a110*b110, a111*b111]]]


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
            entry = [Integer(mod(i,l[0]))]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            Rh[tuple(entry)]=A[tuple(entry)]*B[tuple(entry)]
        return Rh
    else:
        raise ValueError("The Dimensions of the input hypermatrices must match.")

def GeneralHypermatrixHadamardExponent(A,B):
    """
    Procedure for computing Hypermatrix elementwise exponent of two input hypermatrices.

    EXAMPLES:

    ::

        sage: A = HM(2,2,2,'a');B = HM(2,2,2,'b')
        sage: GeneralHypermatrixHadamardExponent(A,B) 
        [[[a000^b000, a001^b001], [a010^b010, a011^b011]], [[a100^b100, a101^b101], [a110^b110, a111^b111]]]


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
            entry = [Integer(mod(i,l[0]))]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            Rh[tuple(entry)]=A[tuple(entry)]^B[tuple(entry)]
        return Rh
    else:
        raise ValueError("The Dimensions of the input hypermatrices must match.")

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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        if len(Set(entry)) == 1:
            Rh[tuple(entry)] = 1
    return Rh

def GeneralHypermatrixMainDiag(od, Lv):
    """
    Outputs a list of lists associated with the general
    Kronecter delta type hypermatrix using the entries of
    list Lv as elements to be placed on the main diagonal
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:

    ::

        sage: Dlt = GeneralHypermatrixMainDiag(2, HM(2,'x').list()); Dlt
        [[x0, 0], [0, x1]] 
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the size parameter
    sz=len(Lv)
    # Initialization of the list specifying the dimensions of the output
    l = [sz for i in range(od)] 
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Initializing the index
    Indx=0
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        if len(Set(entry)) == 1:
            Rh[tuple(entry)] = Lv[Indx]
            Indx=Indx+1
    return Rh

def GeneralHypermatrixKroneckerDeltaL(od, sz):
    """
    Outputs a list of lists associated with the general
    Kronecter delta hypermatrix
    The code only handles the Hypermatrix HM class objects.

    EXAMPLES:

    ::

        sage: L = GeneralHypermatrixKroneckerDeltaL(2,2); L
        [[[1, 0], [0, 0]], [[0, 0], [0, 1]]]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    L=[]
    for t in range(sz):
        # Initialization of the list specifying the dimensions of the output
        l = [sz for i in range(od)] 
        # Initializing the input for generating a symbolic hypermatrix
        inpts = l+['zero']
        # Initialization of the hypermatrix
        Rh = HM(*inpts)
        # Main loop performing the transposition of the entries
        for i in range(prod(l)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entry = [Integer(mod(i,l[0]))]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            if Set(entry).list() == [t]:
                Rh[tuple(entry)] = 1
        L.append(Rh)
    return L

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
    LQ = [HM(*([2 for i in range(od)]+[AlphaB[j]])).elementwise_base_exponent(e) for j in range(od)]
    # Initilizing the list of variable
    VrbLst = []
    for Q in LQ:
        VrbLst = VrbLst + (Q.elementwise_base_logarithm(e)).list()
    # Computing the product
    Eq = GeneralHypermatrixProduct(*[Q for Q in LQ])
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
    return [HM(*([2 for i in range(od)]+[AlphaB[j]])).subs(dict(Dct)) for j in range(od)]

def GeneralUncorrelatedHypermatrixTuple(od):
    """
    Generates a tuplet of hypermatrices of the appropriate order which are
    uncorrelated and normalized. Each one of the dimensions are equal to 2.

    EXAMPLES:

    ::

        sage: [A,B]=GeneralUncorrelatedHypermatrixTuple(2); Prod(A,B).simplify_full()
        [[1, 0], [0, 1]]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initializing the unormalized tuples
    L = GeneralUncorrelatedHypermatrixTupleU(od)
    #Tp= apply(GeneralHypermatrixProduct, [h for h in L])
    Tp= GeneralHypermatrixProduct(*[h for h in L])
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

        sage: Q=GeneralOrthogonalHypermatrixU(3); Rs=Prod(Q,Q.transpose(2),Q.transpose())
        sage: [Rs[i,j,k] for k in range(2) for j in range(2) for i in range(2) if i!=j or j!=k or i!=k]
        [0, 0, 0, 0, 0, 0]


    AUTHORS:

    - Edinah K. Gnang, Ori Parzanchevski and Yuval Filmus
    """
    # Initializing the hypermatrix
    #Q=apply(HM,[2 for i in range(od)]+['q'])
    Q=HM(*[2 for i in range(od)]+['q'])
    # Initilizing the list of variable
    VrbLst=Q.list()
    # Reinitializing of Q by exponentiation 
    Q=Q.elementwise_base_exponent(e)
    # Computing the product
    Eq=GeneralHypermatrixProduct(*[Q.transpose(j) for j in range(od,0,-1)])
    # Writting up the constraints
    LeQ=(Set(Eq.list())).list()
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
        entry = [Integer(mod(i,l[0]))]
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

        sage: Q=GeneralOrthogonalHypermatrix(3); Prod(Q,Q.transpose(2),Q.transpose()).simplify_full()
        [[[1, 0], [0, 0]], [[0, 0], [0, 1]]]


    AUTHORS:

    - Edinah K. Gnang, Ori Parzanchevski and Yuval Filmus
    """
    if od == 2:
        nrm = sqrt(exp(2*var('r1'))+exp(2*var('r2')))
        return HM([[exp(var('r1'))/nrm, exp(var('r2'))/nrm],[-exp(var('r2'))/nrm, exp(var('r1'))/nrm]])
    else :
        # Initializing the hypermatrix
        #Q = apply(HM,[2 for i in range(od)]+['q'])
        Q = HM(*[2 for i in range(od)]+['q'])
        # Initilizing the list of variable
        VrbLst = Q.list()
        # Reinitializing of Q by exponentiation 
        Q = Q.elementwise_base_exponent(e)
        # Computing the product
        Eq = GeneralHypermatrixProduct(*[Q.transpose(j) for j in range(od,0,-1)])
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
            entry = [Integer(mod(i,l[0]))]
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

        sage: [A, Ac]=GeneralUnitaryHypermatrixU(4); Rs=(Prod(A,Ac.transpose(3),A.transpose(2),Ac.transpose())).simplify()
        sage: [Rs[i,j,k,l] for l in range(2) for k in range(2) for j in range(2) for i in range(2) if i!=j or j!=k or k!=l or i!=k or i!=l or j!=k]
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the variable playing the role of sqrt(-1)
    z=var('z')
    # Initializing the real part and the imaginary part
    #X=apply(HM,[2 for i in range(od)]+['x']); Y=apply(HM,[2 for i in range(od)]+['y'])
    X=HM(*([2 for i in range(od)]+['x'])); Y=HM(*([2 for i in range(od)]+['y']))
    # Initialization of the list
    Lh = [(X+z*Y).elementwise_base_exponent(e), (X-z*Y).elementwise_base_exponent(e)]
    # Computation of the Product.
    #B = apply(Prod,[Lh[Integer(mod(i,2))].transpose(i) for i in range(od,0,-1)])
    B = Prod(*[Lh[Integer(mod(i,2))].transpose(i) for i in range(od,0,-1)])
    B[tuple([0 for i in range(od)])]=0; B[tuple([1 for i in range(od)])]=0
    # Initializing the list
    L=Set(B.list()).list()
    # Removing the normalization constraints
    L.remove(0)
    # Initialization of the variables.
    Vrbls =X.list()+Y.list() 
    # Initialization of the non homogeneous equations
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

        sage: [A, Ac]=GeneralUnitaryHypermatrix(4)
        sage: (Prod(A, Ac.transpose(3), A.transpose(2), Ac.transpose())).simplify_full()
        [[[[1, 0], [0, 0]], [[0, 0], [0, 0]]], [[[0, 0], [0, 0]], [[0, 0], [0, 1]]]]

    AUTHORS:
    - Edinah K. Gnang
    """
    Lh=GeneralUnitaryHypermatrixU(od)
    Tp=Prod(*[Lh[Integer(mod(i,2))].transpose(i) for i in range(od,0,-1)])
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

def GeneralHypermatrixReduce(A, VrbL, Rlts):
    """
    Outputs a list of lists associated with the general
    hypermatrix with expressions in the entries reduced 
    modulo the input relations in the list Rlts on the
    variables VrbL. The relation are assume to be monic.


    EXAMPLES:

    ::
        
        sage: VrbL=[var('x'), var('y')]
        sage: Ha=HM(2,2,2,[(VrbL[0]+VrbL[1])^(i+j+k) for i in range(2) for j in range(2) for k in range(2)])
        sage: GeneralHypermatrixReduce(Ha, VrbL, [VrbL[0]^2-5, VrbL[1]^3-7])
        [[[1, x + y], [x + y, 2*x*y + y^2 + 5]], [[x + y, 2*x*y + y^2 + 5], [2*x*y + y^2 + 5, 3*x*y^2 + 5*x + 15*y + 7]]]
        

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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # Initialization of the function
        f=(A[tuple(entry)]).expand()
        # performing the reduction
        for v in range(len(VrbL)):
            for d in range(f.degree(VrbL[v])-Rlts[v].degree(VrbL[v]),-1,-1):
                f=expand(fast_reduce(f,[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))],[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))-expand(Rlts[v]*VrbL[v]^d)])) 
        Rh[tuple(entry)]=f
    return Rh

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
    if Integer(mod(sz,dm)) == 0:
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
        print('Dimension mismatch !!')

def channel_product(A,B):
    """
    Performs channel specific matrix multiplication

    EXAMPLES:

    ::

        sage: Ha=HM(2,2,3,'a'); Hb=HM(2,2,3,'b')
        sage: channel_product(Ha,Hb)
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

def Matrix2HM(A):
    """
    Converts a matrix to a hypermatrix.


    EXAMPLES:

    ::

        sage: Lv=var_list('a',5)
        sage: Matrix2HM(Matrix(SR,[Lv[1:3],Lv[3:]]))
        [[a1, a2], [a3, a4]]


    AUTHORS:
    - Edinah K. Gnang
    """
    return HM(A.nrows(), A.ncols(), [SR(v) for v in A.transpose().list()])

def Deter(A):
    """
    Computes symbolically the determinant of a square matrix
    using the sum over permutation formula.

    EXAMPLES:

    ::

        sage: Deter(HM(2,2,'m'))
        -m01*m10 + m00*m11

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the permutations
    P = Permutations(range(A.nrows()))
    return sum([Permutation([p[i]+1 for i in range(len(p))]).signature()*prod([A[k,p[k]] for k in range(A.nrows())]) for p in P])

def DeterHadamard(A):
    """
    Computes symbolically the determinant of a square matrix
    using the sum over permutation formula.

    EXAMPLES:

    ::

        sage: DeterHadamard(HM(2,2,[HM(2,2,'a'),HM(2,2,'c'),HM(2,2,'b'),HM(2,2,'d')])).printHM()
        [:, :]=
        [-b00*c00 + a00*d00 -b01*c01 + a01*d01]
        [-b10*c10 + a10*d10 -b11*c11 + a11*d11]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the permutations
    P = Permutations(range(A.nrows()))
    return sum([Permutation([p[i]+1 for i in range(len(p))]).signature()*list_elementwise_product([A[k,p[k]] for k in range(A.nrows())]) for p in P])

def ThirdOrderDeter(H):
    """
    computes third order hypermatrix determinant using the recursive determinant construction.
    Where the base case is the matrix case. The third order hyperdeterminant are only defined
    for cubic third order hypermatrices.
    
    EXAMPLES:
 
    ::

        sage: ThirdOrderDeter(HM(2,2,2,'a'))
        a000*a011*a101*a110 - a001*a010*a100*a111
        

    AUTHORS:

    - Edinah K. Gnang
    - To Do: Implement a faster and more generic version.
    """
    # Testing to see that the hypermatrix is indeed a cube
    if len(Set(H.dimensions()).list())==1 and H.order()==3:
        # Initializing the matrix for the mnemonic construction
        A=HM(H.n(0),H.n(0),[var('x'+str(i)+str(j)) for j in range(1,1+H.n(0)) for i in range(1,1+H.n(0))])
        # Computing the mnemonique polynomial
        P=Permutations(range(A.nrows()))
        L=expand(Deter(A)*prod([sum([sqrt(g^2).canonicalize_radical() for g in Deter(A.elementwise_exponent(j)).operands()]) for j in range(2,1+A.n(0))])).operands()
        # Computing the polynomial
        f=sum([l for l in L if len((l^2).operands())==(H.n(0))^2])
        # Loop performing the umbral expression
        for k in range(H.n(0),0,-1):
            f=fast_reduce(f,[var('x'+str(i)+str(j))^k for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0))],[var('a'+str(i)+str(j)+str(k)) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0))])
        return f.subs(dict([(var('a'+str(i)+str(j)+str(k)),H[i-1,j-1,k-1]) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0))]))
    else :
        # Print an error message indicating that the matrix must be a cube.
        raise ValueError("The hypermatrix must be a third order cube hypermatrix.")
 
def Per(A):
    """
    Computes symbolically the permanent of a square matrix
    using the sum over permutation formula.

    EXAMPLES:

    ::

        sage: M = HM(2, 2, 'm'); Per(M)
        m01*m10 + m00*m11

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the permutations
    P = Permutations(rg(A.nrows()))
    return sum([prod([A[k,p[k]] for k in rg(A.nrows())]) for p in P])

def PerII(A):
    """
    Computes symbolically the permanent of a square matrix
    using the sum over permutation formula. This function 
    follows a maple implementation suggested


    EXAMPLES:

    ::

        sage: M = HM(2, 2, 'm'); PerII(M)
        m01*m10 + m00*m11


    AUTHORS:
    - Edinah K. Gnang, Harry Crane
    """
    # Initializing the symbolic matrix
    sz = min(A.nrows(), A.ncols())
    P = GeneratePartition(sz)
    Pml = []
    for part in P:
        Pml.append((-1)^(max(part))*factorial(max(part))*((Partition2HM(part)).elementwise_product(A)).det())
    return expand(sum(Pml))*(-1)^A.nrows()
 
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
    # Initialization of the list specifying the dimensions of the output
    l = [A.n(i) for i in range(A.order())]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = [max(A.dimensions()) for i in range(A.order())]+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=A[tuple(entry)]
    return Rh

def Permutation_to_PseudoTuple(M):
    """
    Returns list of edge tuple desctiption associated with the 
    input permutation. The encoding is based on the diagonals.
    The row of the input nxn matrix M determines where the loop
    edge is placed. The second row of the input matrix M determines
    where the edge whose edge weight equals one. This continues
    up utill we reach the edge of weight (n-1) for there is only
    one possible choice 


    EXAMPLES:
    ::
        sage: Permutation_to_PseudoTuple(HM([[1,0,0],[0,1,0],[0,0,1]]))
        [[(0, 0)], [(0, 1), (1, 0)], [(0, 2), (2, 0)]]
        sage: Permutation_to_PseudoTuple(HM([[0,0,1],[0,1,0],[1,0,0]]))
        [[(2, 2)], [(1, 2), (2, 1)], [(0, 2), (2, 0)]]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of running edge weight parameter and the Tuple list
    edgw=0; tp=[]
    while M.n(0) > 1:
        for i in rg(M.n(0)):
            if M[0,i] == 1:
                M=M.slice([k for k in rg(M.n(1)) if k != i], 'col').slice([j  for j in rg(1,M.n(0))], 'row')
                if edgw == 0:
                    tp.append([(i, i+edgw)])
                else:
                    tp.append([(i, i+edgw), (i+edgw, i)])
                edgw=edgw+1
                break
    tp.append([(0, edgw), (edgw, 0)])
    return tp

def Tuple_to_Adjacency(T):
    """
    The method returns the adjacency matrix of input edge tuple
    description of the input directed graph.


    EXAMPLES:

    ::

        sage: Tuple_to_Adjacency([(0, 1), (1, 2), (2, 0), (3, 3)]).printHM()
        [:, :]= 
        [0 1 0 0]
        [0 0 1 0]
        [1 0 0 0]
        [0 0 0 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=HM(2,len(T),'kronecker')
    return (sum([Id.slice([t[0]],'col')*Id.slice([t[1]],'row') for t in T]))

def Tuple_to_PseudoTuple(tp):
    """
    Returns list of unidrected edge tuple desctiption associated with the 
    input tuple.


    EXAMPLES:
    ::
        sage: Tuple_to_PseudoTuple([(0, 0), (1, 0), (2, 0)])
        [[(0, 0)], [(0, 1), (1, 0)], [(0, 2), (2, 0)]]
        sage: Tuple_to_PseudoTuple([(0, 2), (1, 2), (2, 2)])
        [[(2, 2)], [(1, 2), (2, 1)], [(0, 2), (2, 0)]]



    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the pseudo tuple list
    psdT = []
    for i in rg(len(tp)):
        if tp[i][0]==tp[i][1]:
            psdT.append([(tp[i][0],tp[i][1])])
        else:
            psdT.append([(min(tp[i][0], tp[i][1]), max(tp[i][0], tp[i][1])),   (max(tp[i][1], tp[i][0]), min(tp[i][1], tp[i][0]))])
    T=[]
    for i in rg(len(tp)):
        for j in rg(len(tp)):
            if abs(psdT[j][0][0]-psdT[j][0][1]) == i:
                T.append(psdT[j])
    return T

def PseudoTuple_to_Permutation(psdT):
    """
    Returns a permutation matrix associated with list of edge specified as a pseudo tuple edge list.


    EXAMPLES:
    ::
        sage: PseudoTuple_to_Permutation([[(0, 0)], [(0, 1), (1, 0)], [(2, 0), (0, 2)]])
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        sage: PseudoTuple_to_Permutation([[(2, 2)], [(1, 2), (2, 1)], [(2, 0), (0, 2)]])
        [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
        sage: PseudoTuple_to_Permutation([[(1, 1)], [(1, 2), (2, 1)], [(2, 0), (0, 2)]])
        [[0, 1, 0], [0, 0, 1], [1, 0, 0]]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the input matrix
    M=HM(len(psdT),len(psdT),'zero'); M[0,psdT[0][0][0]]=1
    # Initialization of the forbiden index
    allowed_row_index=rg(1,len(psdT))
    allowed_col_index=[i for i in rg(len(psdT)) if i != psdT[0][0][0]]
    # Initialization fo the counter
    for indx in rg(1,len(psdT)):
        M[indx, allowed_col_index[min(psdT[indx][0])]]=1
        # Updating the allowable indices
        allowed_row_index.remove(indx)
        allowed_col_index.remove(allowed_col_index[min(psdT[indx][0])])
    return M

def Permutation_to_InsertionPattern(T):
    """
    Returns list of the list of insertion patterns associated
    with the bijection from permutation to gracefully labeled
    undirected graph having a loop edge (of weight 0). The 
    input is tuple descrition of a permutation. The function 
    checks to see that it is indeed a permutation.
    The output is a list of insertion patterns reflecting the
    edge choice made in deacreasin order of subtractive edge
    weight magnitudes.


    EXAMPLES:
    ::
        sage: Permutation_to_InsertionPattern([(0, 0), (1, 1), (2, 2)])
        [0, 0, 0]
        sage: Permutation_to_InsertionPattern([(0, 2), (1, 1), (2, 0)])
        [0, 1, 2]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the size and order parameters
    sz=len(T); od=2
    # Initialization of the identity matrix
    Id=HM(od, sz, 'kronecker')
    P=sum(Id.slice([T[i][0]],'col')*Id.slice([T[i][1]],'row') for i in rg(sz))
    if (P*P.transpose()-Id).is_zero():
        # Initialization of the list of options for each subtractive edge weights
        Lopt=[[(i,i) for i in rg(sz)]]+[[] for i in rg(sz-1)]
        for i in rg(1,sz):
            for j in rg(i):
                Lopt[abs(i-j)].append((j,i))
        # Initialization of the list of dictionaries for edge weight options
        LDct=[{Lopt[i][j] : j for j in rg(len(Lopt[i]))} for i in rg(sz)]
        # Converting the permutation matrix into a Pseudotuple list
        # to exploit the correspondence between permutations and 
        # gracefully labeled undirected graph having a loop edge.
        pT=Permutation_to_PseudoTuple(P)
        # Returing the reveserd difference this ensure that the insertion
        # entry never exceed the length of the partial list. For instance
        # the first entry is always zero to account for the fact that we
        # start from the constant 1.
        return [LDct[pT[sz-1-i][0][1]-pT[sz-1-i][0][0]][(pT[sz-1-i][0][0],pT[sz-1-i][0][1])] for i in rg(sz)]
    else:
        raise ValueError("The input must be a spanning union of cycles as input.")

def GenerateUnitLpNormVector(n, p = 2, indx=0):
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

def GenerateUnitLpNormVectorII(sz,p=2,indx=0):
    """
    outputs a unit lp norm vector.

    EXAMPLES:

    ::

        sage: GenerateUnitLpNormVectorII(2) 
        [1/2*e^(I*t0) + 1/2*e^(-I*t0), -1/2*I*e^(I*t0) + 1/2*I*e^(-I*t0)]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    Cos(x)=exp(I*x)/(2*1) + exp(-I*x)/(2*1); Sin(x)=exp(I*x)/(2*I) - exp(-I*x)/(2*I)
    if n == 1:
        return [1]
    else :
        X = []
        X.append( Cos(var('t'+str(indx)))^(2/p) )
        for i in range(1,sz-1):
            X.append( prod(Sin( var('t'+str(j+indx)) ) for j in range(i))*Cos(var('t'+str(i+indx)))^(2/p) )
        X.append( prod(Sin(var('t'+str(j+indx)))^(2/p) for j in range(sz-1)) )
        return X

def GenerateUnitLpNormVectorIII(sz, p, indx, T):
    """
    outputs a unit lp norm vector. The implmentation
    here uses exponential polynomial expression of 
    of trigonometric polynomials. The input sz
    specifies the length of the vector. The input 
    p specifies the norm which we choose. This implementation
    differs from the ones above in the fact that accounts
    for the insertion pattern specified by the input tuple
    description of the a permutation T. The insertion pattern
    is devised from T by using the correspondence between 
    permutations and gracefully labeled undirected graph having
    a loop edge.
    Crucially T must a permutation of sz-1 elements.


    EXAMPLES:

    ::

        sage: GenerateUnitLpNormVectorIII(2, 2, 0, [(0,0)]) 
        [-1/2*I*e^(I*t0) + 1/2*I*e^(-I*t0), 1/2*e^(I*t0) + 1/2*e^(-I*t0)]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    if n == 1:
        return [1]
    else :
        # Initialization of the exponential encoding of trigonometric functions
        Cos(x) = exp(I*x)/(2*1) + exp(-I*x)/(2*1); Sin(x) = exp(I*x)/(2*I) - exp(-I*x)/(2*I)
        # Initialization of the insertion of sequence
        L = Permutation_to_InsertionPattern(T)
        #X = [Cos(var('t'+str(indx)))^(2/p), Sin(var('t'+str(indx)))^(2/p)]
        X = [1]
        #for i in range(1,sz-1):
        for i in range(sz-1):
            tmp = X[L[i]]
            X[L[i]] = tmp*Cos(var('t'+str(i+indx)))^(2/p)
            X.insert(L[i],tmp*Sin(var('t'+str(i+indx)))^(2/p))
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

def ProbabilityHM(sz, xi=0):
    """
    outputs the symbolic parametrization of a doubly stochastic matrix

    EXAMPLES:

    ::

        sage: ProbabilityHM(2).printHM()
        [:, :]=
        [     cos(t0)^2      sin(t0)^2]
        [-cos(t0)^2 + 1 -sin(t0)^2 + 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the matrix to be filled
    M = HM(sz,sz, 'zero')
    # Initialixzing the variable index
    indx=xi
    for c in rg(sz-1):
        # Initializing the probability vector associated with the row c
        La = GenerateUnitLpNormVector(sz-c,1,indx)
        # Updating the variable index
        indx = indx+len(La)-1
        # Initializing the probability vector associated with the column c
        Lb = GenerateUnitLpNormVector(sz-c-1,1,indx)
        # Updating the variable index
        indx = indx+len(Lb)-1
        # Loop which fills up the Matrix
        for i in range(c, c+len(La)):
            if c > 0:
                # Filling up the row c of the Matrix M
                M[c,i] = (1-sum([M[c,j] for j in rg(c)]))*La[i-c]
                if i > c:
                    # Filling up the column c of the Matrix M
                    M[i,c] = (1-sum([M[j,c] for j in rg(c+1)]))*Lb[i-c-1]
            else:
                # Filling up the row c of the Matrix M
                M[c,i] = La[i-c]
                if i > c:
                    # Filling up the column c of the Matrix M
                    M[i,c] = (1-sum([M[j,c] for j in rg(c+1)]))*Lb[i-c-1]
    M[sz-1,sz-1]=1-sum([M[j,sz-1] for j in rg(sz-1)])
    return M

def ProbabilityHMII(sz, xi=0):
    """
    outputs the symbolic parametrization of a doubly stochastic matrix

    EXAMPLES:

    ::

        sage: ProbabilityHMII(2).printHM()
        [:, :]=
        [                                 1/4*(e^(I*t0) + e^(-I*t0))^2                         (-1/2*I*e^(I*t0) + 1/2*I*e^(-I*t0))^2]
        [-1/16*((e^(I*t0) + e^(-I*t0))^2 - 4)*(e^(I*t1) + e^(-I*t1))^2                    -(-1/2*I*e^(I*t0) + 1/2*I*e^(-I*t0))^2 + 1]
 

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the matrix to be filled
    M = HM(sz,sz,'zero')
    # Initialixzing the variable index
    indx=xi
    for c in rg(sz-1):
        # Initializing the probability vector associated with the row c
        La = GenerateUnitLpNormVectorII(sz-c,1,indx)
        # Updating the variable index
        indx = indx+len(La)-1
        # Initializing the probability vector associated with the column c
        Lb = GenerateUnitLpNormVectorII(sz-c-1,1,indx)
        # Updating the variable index
        indx = indx+len(Lb)-1
        # Loop which fills up the Matrix
        for i in range(c, c+len(La)):
            if c > 0:
                # Filling up the row c of the Matrix M
                M[c,i] = (1-sum([M[c,j] for j in rg(c)]))*La[i-c]
                if i > c:
                    # Filling up the column c of the Matrix M
                    M[i,c] = (1-sum([M[j,c] for j in rg(c+1)]))*Lb[i-c-1]
            else:
                # Filling up the row c of the Matrix M
                M[c,i] = La[i-c]
                if i > c:
                    # Filling up the column c of the Matrix M
                    M[i,c] = (1-sum([M[j,c] for j in rg(c+1)]))*Lb[i-c-1]
    M[sz-1,sz-1]=1-sum([M[j,sz-1] for j in rg(sz-1)])
    return M

def ProbabilitySymHMIII(sz, xi, Lt):
    """
    Outputs a symbolic parametrization of a symetric doubly stochastic matrix.
    where the points on the unit sphere a constructed according to the permutation
    patterns specified in the list Lt. Lt must be carefully constructed to account
    for every call of the function which generates symbolic vectors of l1 norm 1


    EXAMPLES:

    ::

        sage: ProbabilitySymHMIII(2,0,[[(0,0)]]).printHM()
        [:, :]=
        [(-1/2*I*e^(I*t0) + 1/2*I*e^(-I*t0))^2          1/4*(e^(I*t0) + e^(-I*t0))^2]
        [         1/4*(e^(I*t0) + e^(-I*t0))^2     -1/4*(e^(I*t0) + e^(-I*t0))^2 + 1]

    AUTHORS:

    - Edinah K. Gnang
    """
    # Initializing the matrix to be filled
    M = HM(sz, sz, 'zero')
    # Initialixzing the variable index
    indx=xi
    # Initialization of the index which keeps track of the permutation
    # specifiying the isertion patterns
    indxp=0
    for c in rg(sz-1):
        # Initializing the probability vector associated with the row c
        La = GenerateUnitLpNormVectorIII(sz-c, 1, indx, Lt[indxp])
        # Updating the insertion pattern index
        indxp = indxp+1
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
    M[sz-1,sz-1]=1-sum([M[j,sz-1] for j in range(sz-1)])
    return M

def ProbabilityHMIII(sz, xi, Lt):
    """
    Outputs a symbolic parametrization of a symetric doubly stochastic matrix.
    where the points on the unit sphere a constructed according to the permutation
    patterns specified in the list Lt


    EXAMPLES:

    ::

        sage: ProbabilityHMIII(2,0,[[(0,0)],[(0,0)]]).printHM()
        [:, :]=
        [                                    (-1/2*I*e^(I*t0) + 1/2*I*e^(-I*t0))^2                                              1/4*(e^(I*t0) + e^(-I*t0))^2]
        [-1/4*((-1/2*I*e^(I*t0) + 1/2*I*e^(-I*t0))^2 - 1)*(e^(I*t1) + e^(-I*t1))^2                                         -1/4*(e^(I*t0) + e^(-I*t0))^2 + 1] 


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the matrix to be filled
    M = HM(sz,sz,'zero')
    # Initialixzing the variable index
    indx=xi
    # Initialization of the index which keeps track of the permutation
    # specifiying the isertion patterns
    indxp=0
    for c in rg(sz-1):
        # Initializing the probability vector associated with the row c
        La = GenerateUnitLpNormVectorIII(sz-c, 1, indx, Lt[indxp])
        # Updating the insertion pattern index
        indxp = indxp+1
        # Updating the variable index
        indx = indx+len(La)-1
        # Initializing the probability vector associated with the column c
        Lb = GenerateUnitLpNormVectorII(sz-c-1,1,indx)
        # Updating the variable index
        indx = indx+len(Lb)-1
        # Loop which fills up the Matrix
        for i in range(c, c+len(La)):
            if c > 0:
                # Filling up the row c of the Matrix M
                M[c,i] = (1-sum([M[c,j] for j in rg(c)]))*La[i-c]
                if i > c:
                    # Filling up the column c of the Matrix M
                    M[i,c] = (1-sum([M[j,c] for j in rg(c+1)]))*Lb[i-c-1]
            else:
                # Filling up the row c of the Matrix M
                M[c,i] = La[i-c]
                if i > c:
                    # Filling up the column c of the Matrix M
                    M[i,c] = (1-sum([M[j,c] for j in rg(c+1)]))*Lb[i-c-1]
    M[sz-1,sz-1]=1-sum([M[j,sz-1] for j in rg(sz-1)])
    return M

def ProbabilitySymMatrix(n, xi=0):
    """
    outputs the symbolic parametrization of a symetric doubly stochastic matrix

    EXAMPLES:

    ::

        sage: ProbabilitySymMatrix(2)
        [     cos(t0)^2      sin(t0)^2]
        [     sin(t0)^2 -sin(t0)^2 + 1]
        

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

def ProbabilitySymHM(sz, xi=0):
    """
    outputs the symbolic parametrization of a symetric doubly stochastic matrix

    EXAMPLES:

    ::

        sage: ProbabilitySymHM(2).printHM()
        [:, :]=
        [     cos(t0)^2      sin(t0)^2]
        [     sin(t0)^2 -sin(t0)^2 + 1]
        

    AUTHORS:

    - Edinah K. Gnang
    """
    # Initializing the matrix to be filled
    M = HM(sz, sz, 'zero')
    # Initialixzing the variable index
    indx=xi
    for c in rg(sz-1):
        # Initializing the probability vector associated with the row c
        La = GenerateUnitLpNormVector(sz-c,1,indx)
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
    M[sz-1,sz-1]=1-sum([M[j,sz-1] for j in range(sz-1)])
    return M

def ProbabilitySymHMII(sz, xi=0):
    """
    outputs the symbolic parametrization of a symetric doubly stochastic matrix

    EXAMPLES:

    ::

        sage: ProbabilitySymHMII(2).printHM()
        [:, :]=
        [              1/4*(e^(I*t0) + e^(-I*t0))^2      (-1/2*I*e^(I*t0) + 1/2*I*e^(-I*t0))^2]
        [     (-1/2*I*e^(I*t0) + 1/2*I*e^(-I*t0))^2 -(-1/2*I*e^(I*t0) + 1/2*I*e^(-I*t0))^2 + 1]        

    AUTHORS:

    - Edinah K. Gnang
    """
    # Initializing the matrix to be filled
    M = HM(sz, sz, 'zero')
    # Initialixzing the variable index
    indx=xi
    for c in rg(sz-1):
        # Initializing the probability vector associated with the row c
        La = GenerateUnitLpNormVectorII(sz-c,1,indx)
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
    M[sz-1,sz-1]=1-sum([M[j,sz-1] for j in range(sz-1)])
    return M

def ProbabilityHMIV(sz,c):
    """
    Outputs a symbolic parametrization of a doubly stochastic matrix.
    Using the analog of the adjoint matrix based upon the permanent.
    The inputs correspond respectively to the size and the caracter
    to use for the symbolic parametrization.


    EXAMPLES:

    ::

        sage: ProbabilityHMIV(2,'a').printHM()
        [:, :]=
        [a00*a11/(a01*a10 + a00*a11) a01*a10/(a01*a10 + a00*a11)]
        [a01*a10/(a01*a10 + a00*a11) a00*a11/(a01*a10 + a00*a11)]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the symbolic matrix
    A = HM(sz, sz, c)
    # Computing the permanent of the matrix
    F = PerII(A)
    # Initialization of the desired matrix
    return (1/F)*HM(sz, sz, [A[i,j]*diff(F,A[i,j]) for j in rg(sz) for i in rg(sz)])

def GenerateRandomOneZeroHM(sz, k):
    """
    Outputs uniformly randomly chosen 0,1 vector having k ones and sz-k zeros
    with the index fi being forbidden

    EXAMPLES:

    ::

        sage: GenerateRandomOneZeroHM(2,2)
        [1, 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    if type(k)==type(2) and k > 0 :
        # Initialization of the set of powers of 2
        S=Set([2^i for i in rg(sz)])
        # Initialization of the list of integers
        L=[sum(s) for s in S.subsets(k)]
        # Chosing uniformly at random a memeber of L
        if len(L)==0 :
            return HM(sz,'zero').list()
        else :
            rn=L[randint(0,len(L)-1)]
            # obtaining the bits of the chose number
            bL=rn.bits()
            return bL+[0 for i in rg(len(bL),sz)]
    else :
        return HM(sz,'zero').list()

def Random_d_regular_adjacencyHM(sz, d):
    """
    Outputs a complementary pair of directed graphs where 
    the first graphs has vertex in-degree and out-degree 
    both equal to d for every vertex. While the second has
    vertex in-degree and out-degree both equal to sz-d-1
    for all of it's vertices.


    EXAMPLES:

    ::

        sage: sz=5; [A,B]=Random_d_regular_adjacencyHM(sz,3)
        sage: (A*HM(sz,1,'one')).printHM()
        [:, :]=
        [3]
        [3]
        [3]
        [3]
        [3]
        sage: (HM(1,sz,'one')*A).printHM()
        [:, :]=
        [3 3 3 3 3]
        

    AUTHORS:
    - Edinah K. Gnang, Yan Jiang
    """
    # Initializing the matrix to be filled
    M = HM(sz,sz, 'zero'); M[sz-1,sz-1]=1
    # List of regular graph on 2 vertices
    while M[sz-1,sz-1]!=0 or M*HM(sz,1,'one') != HM(sz,1,'one')*d:
        for c in rg(sz-1):
            # Initializing the vector associated with the row c
            La = [0]+GenerateRandomOneZeroHM(sz-c-1, d-sum(M[c,j] for j in rg(c)))
            # Initializing the vector associated with the column c
            Lb = [ ]+GenerateRandomOneZeroHM(sz-c-1, d-sum(M[j,c] for j in rg(c+1)))
            # Loop which fills up the Matrix
            for i in range(c, c+len(La)):
                    M[c,i] = La[i-c]
                    if i > c:
                        # Filling up the column c of the Matrix M
                        M[i,c] = Lb[i-c-1]
        M[sz-1,sz-1]=d-sum([M[j,sz-1] for j in rg(sz-1)])
    return [M, HM(sz,sz,'one')-M-HM(2,sz,'kronecker')]

def Random_d_regular_undirected_adjacencyHM(sz, d):
    """
    Outputs a complementary pair of undirected graphs where
    the first graphs has vertex in-degree and out-degree
    both equal to d for every vertex. While the second has
    vertex in-degree and out-degree both equal to sz-d-1
    for all of it's vertices.
    Note that if sz is odd, then d must be even because
    2*|E| = sz * d in case both d and sz are odd we have
    a contradictions
    
 
    EXAMPLES:
    
    ::
    
    sage: sz=5; [A,B]=Random_d_regular_undirected_adjacencyHM(sz,2)
    sage: A.printHM()
    [:, :]=
    [0 1 0 0 1]
    [1 0 0 1 0]
    [0 0 0 1 1]
    [0 1 1 0 0]
    [1 0 1 0 0]
    sage: (A*HM(sz,1,'one')).printHM()
    [:, :]=
    [2]
    [2]
    [2]
    [2]
    [2]
    sage: (HM(1,sz,'one')*A).printHM()
    [:, :]=
    [2 2 2 2 2]
    sage: (A-A.transpose()).printHM()
    [:, :]=
    [0 0 0 0 0]
    [0 0 0 0 0]
    [0 0 0 0 0]
    [0 0 0 0 0]
    [0 0 0 0 0]
    
    
    AUTHORS:
    - Yan Jiang, Edinah K. Gnang
    """
    # Initializing the matrix to be filled
    M = HM(sz,sz, 'zero'); M[sz-1,sz-1]=1
    # List of regular graph on 2 vertices
    while M[sz-1,sz-1]!=0 or M*HM(sz,1,'one') != HM(sz,1,'one')*d:
        for c in rg(sz-1):
            # Initializing the vector associated with the row c
            La = [0]+GenerateRandomOneZeroHM(sz-c-1,d-sum(M[c,j] for j in rg(c)))
            # Loop which fills up the Matrix
            for i in range(c, c+len(La)):
                # Filling up the row c of the Matrix M
                M[c,i] = La[i-c]
                if i > c:
                    # Filling up the column c of the Matrix M
                    M[i,c] = La[i-c]
        M[sz-1,sz-1]=d-sum([M[j,sz-1] for j in rg(sz-1)])
    return [M, HM(sz,sz,'one')-M-HM(2,sz,'kronecker')]

def Random_3_regular_directed_adjacencyHM(sz):
    """
    Outputs a complementary pair of undirected graphs where
    the first graphs has vertex in-degree and out-degree
    both equal to 3 for every vertex. While the second has
    vertex in-degree and out-degree both equal to sz-3-1
    for all of it's vertices.
    Note that if sz is odd, then d must be even because
    2*|E| = sz * 3 in case both d and sz are odd we have
    a contradictions
    
 
    EXAMPLES:
    
    ::
    
    sage: sz=12; [A,B]=Random_3_regular_directed_adjacencyHM(sz)
    sage: (A*HM(sz,1,'one')).printHM()
    [:, :]=
    [3]
    [3]
    [3]
    [3]
    [3]
    [3]
    [3]
    [3]
    [3]
    [3]
    [3]
    [3]
    
    
    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the random permutations
    Sigma=RandomDerangementTuple(sz); Gamma=RandomDerangementTuple(sz)
    while prod(Sigma[i][1]-Gamma[i][1] for i in rg(sz)).is_zero() or not (Tuple_to_Adjacency(Sigma) + Tuple_to_Adjacency(Gamma) + HM(2,sz,'kronecker')).index_rotation(pi/2).trace().is_zero():
        Sigma=RandomDerangementTuple(sz); Gamma=RandomDerangementTuple(sz)
    return [(Tuple_to_Adjacency(Sigma) + Tuple_to_Adjacency(Gamma) + HM(2,sz,'kronecker')).index_rotation(pi/2),HM(sz,sz,'one')-(Tuple_to_Adjacency(Sigma) + Tuple_to_Adjacency(Gamma) + HM(2,sz,'kronecker')).index_rotation(pi/2)]

def BipInflateHM(sz):
    """
    outputs bipartite where vertex in-degree and out-degree 
    both equal to (sz-1)/2 for every vertex. The partition
    are made up of complementary adjacency matrices.
    The number of vertices sz must be odd. 

    EXAMPLES:

    ::

        sage: sz=5; A= BipInflateHM(sz)
        sage: (A*HM(2*sz,1,'one')).printHM()
        [:, :]=
        [2]
        [2]
        [2]
        [2]
        [2]
        [2]
        [2]
        [2]
        [2]
        [2]
        sage: (HM(1,2*sz,'one')*A).printHM()
        [:, :]=
        [2 2 2 2 2 2 2 2 2 2]

    AUTHORS:
    - Edinah K. Gnang, James Murphy
    """
    if Integer(mod(sz,2))==1:
        d=Integer((sz-1)/2)
        [A,B]=Random_d_regular_adjacencyHM(sz, d)
        return HM([[0,1],[0,0]]).tensor_product(A)+HM([[0,0],[1,0]]).tensor_product(B)
    else:
        raise ValueError("The number of vertices must be odd")

def HypermatrixPseudoInversePairsII(A,B):
    """
     Outputs the pseudo inverse pairs associated with the input pairs of hypermatrices

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
        entry = [Integer(mod(i,l[0]))]
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

        sage: A=GenerateRandomHypermatrix(2,2,2); A.dimensions()
        [2, 2, 2]

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
            entry = [Integer(mod(i,l[0]))]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            Rh[tuple(entry)]=random()
        return Rh
    else :
        raise ValueError("The Dimensions must all be non-zero.")

def GenerateRandomRationalHypermatrix(*l):
    """
     Outputs a random hypermatrix

    EXAMPLES:

    ::

        sage: A=GenerateRandomHypermatrix(2,2,2); A.dimensions()
        [2, 2, 2]

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
            entry = [Integer(mod(i,l[0]))]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            Rh[tuple(entry)]=QQ(random())
        return Rh
    else :
        raise ValueError("The Dimensions must all be non-zero.")

def GenerateRandomIntegerHypermatrix(*l):
    """
     Outputs a random hypermatrix

    EXAMPLES:

    ::

        sage: A=GenerateRandomIntegerHypermatrix(3,3,3); A.dimensions()
        [3, 3, 3]

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
            entry = [Integer(mod(i,l[0]))]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            Rh[tuple(entry)]=ZZ.random_element()
        return Rh
    else :
        raise ValueError("The Dimensions must all be non-zero.")

def GenerateRandomBinaryHypermatrix(*l):
    """
     Outputs a random hypermatrix

    EXAMPLES:

    ::

        sage: A=GenerateRandomBinaryHypermatrix(3,3,3); A.dimensions()
        [3, 3, 3]

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
            entry = [Integer(mod(i,l[0]))]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            Rh[tuple(entry)]=Integer(mod(ZZ.random_element(),2))
        return Rh
    else :
        raise ValueError("The Dimensions must all be non-zero.")

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
        raise ValueError("The order must be greater or equal to 2.")
    elif od == 2:
        return HM([[cos(t)^2, sin(t)^2], [sin(t)^2, cos(t)^2]])
    elif od > 2 and type(od) == type(1):
        dms = [2,2]
        B = GeneralStochasticHypermatrix(t,2)
        for z in range(2,od):
            A = B 
            dms.append(2) 
            #B = apply(HM, dms+['zero'])
            B = HM(*(dms+['zero']))
            l = [A.n(i) for i in range(A.order())]
            # Main loop performing the transposition of the entries
            for i in range(prod(l)):
                # Turning the index i into an hypermatrix array location using the decimal encoding trick
                entry = [Integer(mod(i,l[0]))]
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
        raise ValueError("The order must be a positive integer")

def GeneralStochasticHypermatrixII(t, od):
    """
    Generates an stochastic hypermatrix of the appropriate order
    for which each one of the dimensions are equal to 2.
    The vectors are not normalized.

    EXAMPLES:

    ::

        sage: GeneralStochasticHypermatrixII(var('t'), 2)
        [[1/4*(e^(I*t) + e^(-I*t))^2, (-1/2*I*e^(I*t) + 1/2*I*e^(-I*t))^2], [(-1/2*I*e^(I*t) + 1/2*I*e^(-I*t))^2, 1/4*(e^(I*t) + e^(-I*t))^2]]


    AUTHORS:
    - Edinah K. Gnang, Ori Parzanchevski
    """
    if od < 2 and type(od)==type(1):
        raise ValueError("The order must be greater or equal to 2.")
    elif od == 2:
        return HM([[(exp(I*t)/2+exp(-I*t)/2)^2, (exp(I*t)/(2*I)-exp(-I*t)/(2*I))^2], [(exp(I*t)/(2*I)-exp(-I*t)/(2*I))^2, (exp(I*t)/2+exp(-I*t)/2)^2]])
    elif od > 2 and type(od) == type(1):
        dms = [2,2]
        B = GeneralStochasticHypermatrixII(t,2)
        for z in range(2,od):
            A = B 
            dms.append(2) 
            #B = apply(HM, dms+['zero'])
            B = HM(*(dms+['zero']))
            l = [A.n(i) for i in range(A.order())]
            # Main loop performing the transposition of the entries
            for i in range(prod(l)):
                # Turning the index i into an hypermatrix array location using the decimal encoding trick
                entry = [Integer(mod(i,l[0]))]
                sm = Integer(mod(i,l[0]))
                for k in range(len(l)-1):
                    entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                    sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                B[tuple(entry+[0])]=A[tuple(entry)]
                if A[tuple(entry)] == (exp(I*t)/2+exp(-I*t)/2)^2:
                    B[tuple(entry+[1])]=(exp(I*t)/(2*I)-exp(-I*t)/(2*I))^2
                else :
                    B[tuple(entry+[1])]=(exp(I*t)/2+exp(-I*t)/2)^2
        return B
    else :
        raise ValueError("The order must be a positive integer")

def SecondOrderSliceKroneckerProduct(Ha, Hb):
    """
    Computes the Kronecker Product for the two input second order hypermatrices.

    EXAMPLES:

    ::

        sage: SecondOrderSliceKroneckerProduct(HM(2,2,'a'),HM(2,2,'b'))
        [[a00*b00, a00*b01, a01*b00, a01*b01], [a00*b10, a00*b11, a01*b10, a01*b11], [a10*b00, a10*b01, a11*b00, a11*b01], [a10*b10, a10*b11, a11*b10, a11*b11]]

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

        sage: ThirdOrderSliceKroneckerProduct(HM(2,2,2,'a'),HM(3,3,3,'b'))
        [[[a000*b000, a000*b001, a000*b002, a001*b000, a001*b001, a001*b002], [a000*b010, a000*b011, a000*b012, a001*b010, a001*b011, a001*b012], [a000*b020, a000*b021, a000*b022, a001*b020, a001*b021, a001*b022], [a010*b000, a010*b001, a010*b002, a011*b000, a011*b001, a011*b002], [a010*b010, a010*b011, a010*b012, a011*b010, a011*b011, a011*b012], [a010*b020, a010*b021, a010*b022, a011*b020, a011*b021, a011*b022]], [[a000*b100, a000*b101, a000*b102, a001*b100, a001*b101, a001*b102], [a000*b110, a000*b111, a000*b112, a001*b110, a001*b111, a001*b112], [a000*b120, a000*b121, a000*b122, a001*b120, a001*b121, a001*b122], [a010*b100, a010*b101, a010*b102, a011*b100, a011*b101, a011*b102], [a010*b110, a010*b111, a010*b112, a011*b110, a011*b111, a011*b112], [a010*b120, a010*b121, a010*b122, a011*b120, a011*b121, a011*b122]], [[a000*b200, a000*b201, a000*b202, a001*b200, a001*b201, a001*b202], [a000*b210, a000*b211, a000*b212, a001*b210, a001*b211, a001*b212], [a000*b220, a000*b221, a000*b222, a001*b220, a001*b221, a001*b222], [a010*b200, a010*b201, a010*b202, a011*b200, a011*b201, a011*b202], [a010*b210, a010*b211, a010*b212, a011*b210, a011*b211, a011*b212], [a010*b220, a010*b221, a010*b222, a011*b220, a011*b221, a011*b222]], [[a100*b000, a100*b001, a100*b002, a101*b000, a101*b001, a101*b002], [a100*b010, a100*b011, a100*b012, a101*b010, a101*b011, a101*b012], [a100*b020, a100*b021, a100*b022, a101*b020, a101*b021, a101*b022], [a110*b000, a110*b001, a110*b002, a111*b000, a111*b001, a111*b002], [a110*b010, a110*b011, a110*b012, a111*b010, a111*b011, a111*b012], [a110*b020, a110*b021, a110*b022, a111*b020, a111*b021, a111*b022]], [[a100*b100, a100*b101, a100*b102, a101*b100, a101*b101, a101*b102], [a100*b110, a100*b111, a100*b112, a101*b110, a101*b111, a101*b112], [a100*b120, a100*b121, a100*b122, a101*b120, a101*b121, a101*b122], [a110*b100, a110*b101, a110*b102, a111*b100, a111*b101, a111*b102], [a110*b110, a110*b111, a110*b112, a111*b110, a111*b111, a111*b112], [a110*b120, a110*b121, a110*b122, a111*b120, a111*b121, a111*b122]], [[a100*b200, a100*b201, a100*b202, a101*b200, a101*b201, a101*b202], [a100*b210, a100*b211, a100*b212, a101*b210, a101*b211, a101*b212], [a100*b220, a100*b221, a100*b222, a101*b220, a101*b221, a101*b222], [a110*b200, a110*b201, a110*b202, a111*b200, a111*b201, a111*b202], [a110*b210, a110*b211, a110*b212, a111*b210, a111*b211, a111*b212], [a110*b220, a110*b221, a110*b222, a111*b220, a111*b221, a111*b222]]]


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

        sage: Hc=FourthOrderSliceKroneckerProduct(HM(2,2,2,2,'a'),HM(2,2,2,2,'b'))
        sage: Hc.dimensions()
        [4, 4, 4, 4]

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

        sage: Hc=FifthOrderSliceKroneckerProduct(HM(2,2,2,2,2,'a'),HM(2,2,2,2,2,'b'))
        sage: Hc.dimensions()
        [4, 4, 4, 4, 4]


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

def SixthOrderSliceKroneckerProduct(Ha, Hb):
    """
    Computes the Kronecker Product for the two input fifth order hypermatrices.

    EXAMPLES:

    ::

        sage: Hc=SixthOrderSliceKroneckerProduct(HM(2,2,2,2,2,2,'a'),HM(2,2,2,2,2,2,'b'))
        sage: Hc.dimensions()
        [4, 4, 4, 4, 4, 4]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the hypermatrix
    Hc = HM(Ha.n(0)*Hb.n(0), Ha.n(1)*Hb.n(1), Ha.n(2)*Hb.n(2), Ha.n(3)*Hb.n(3), Ha.n(4)*Hb.n(4), Ha.n(5)*Hb.n(5), 'zero')
    for i0 in range(Ha.n(0)):
        for i1 in range(Ha.n(1)):
            for i2 in range(Ha.n(2)):
                for i3 in range(Ha.n(3)):
                    for i4 in range(Ha.n(4)):
                        for i5 in range(Ha.n(5)):
                            for j0 in range(Hb.n(0)):
                                for j1 in range(Hb.n(1)):
                                    for j2 in range(Hb.n(2)):
                                        for j3 in range(Hb.n(3)):
                                            for j4 in range(Hb.n(4)):
                                                for j5 in range(Hb.n(5)):
                                                    Hc[Hb.n(0)*i0+j0,Hb.n(1)*i1+j1,Hb.n(2)*i2+j2,Hb.n(3)*i3+j3,Hb.n(4)*i4+j4,Hb.n(5)*i5+j5]=Ha[i0,i1,i2,i3,i4,i5]*Hb[j0,j1,j2,j3,j4,j5]
    return Hc

def SeventhOrderSliceKroneckerProduct(Ha, Hb):
    """
    Computes the Kronecker Product for the two input fifth order hypermatrices.

    EXAMPLES:

    ::

        sage: A=SeventhOrderSliceKroneckerProduct(HM(2,2,2,2,2,2,2,'a'),HM(2,2,2,2,2,2,2,'b'))
        sage: A.dimensions()
        [4, 4, 4, 4, 4, 4, 4]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the hypermatrix
    Hc = HM(Ha.n(0)*Hb.n(0), Ha.n(1)*Hb.n(1), Ha.n(2)*Hb.n(2), Ha.n(3)*Hb.n(3), Ha.n(4)*Hb.n(4), Ha.n(5)*Hb.n(5), Ha.n(6)*Hb.n(6), 'zero')
    for i0 in range(Ha.n(0)):
        for i1 in range(Ha.n(1)):
            for i2 in range(Ha.n(2)):
                for i3 in range(Ha.n(3)):
                    for i4 in range(Ha.n(4)):
                        for i5 in range(Ha.n(5)):
                            for i6 in range(Ha.n(6)):
                                for j0 in range(Hb.n(0)):
                                    for j1 in range(Hb.n(1)):
                                        for j2 in range(Hb.n(2)):
                                            for j3 in range(Hb.n(3)):
                                                for j4 in range(Hb.n(4)):
                                                    for j5 in range(Hb.n(5)):
                                                        for j6 in range(Hb.n(6)):
                                                            Hc[Hb.n(0)*i0+j0,Hb.n(1)*i1+j1,Hb.n(2)*i2+j2,Hb.n(3)*i3+j3,Hb.n(4)*i4+j4,Hb.n(5)*i5+j5,Hb.n(6)*i6+j6]=Ha[i0,i1,i2,i3,i4,i5,i6]*Hb[j0,j1,j2,j3,j4,j5,j6]
    return Hc

def HypermatrixSliceKroneckerPower(U,n):
    """
    Computes the repeated Kronecker Product.

    EXAMPLES:

    ::

        sage: HypermatrixSliceKroneckerPower(HM(2,2,'a'),2)
        [[a00^2, a00*a01, a00*a01, a01^2], [a00*a10, a00*a11, a01*a10, a01*a11], [a00*a10, a01*a10, a00*a11, a01*a11], [a10^2, a10*a11, a10*a11, a11^2]]


    AUTHORS:
    - Edinah K. Gnang
    """
    T = U.copy()
    for i in range(n-1):
        T = T.tensor_product(U)
    return T

def ThirdOrderHypermatrixKroneckerSum(A,B):
    """
    Computes the  Kronecker sum for third order hypermatrix.

    EXAMPLES:

    ::

        sage: ThirdOrderHypermatrixKroneckerSum(HM(2,2,2,'a'),HM(2,2,2,'b'))
        [[[a000, a001, 0, 0], [a010, a011, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], [[a100, a101, 0, 0], [a110, a111, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, b000, b001], [0, 0, b010, b011]], [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, b100, b101], [0, 0, b110, b111]]]
        

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
        raise ValueError("The order of the input hypermatrices must match and be equal to 3.")

def GeneralHypermatrixKroneckerSum(A,B):
    """
    Computes the  Kronecker sum for arbitrary order hypermatrix.

    EXAMPLES:

    ::

        sage: GeneralHypermatrixKroneckerSum(HM(2,2,'a'), HM(2,2,'b'))
        [[a00, a01, 0, 0], [a10, a11, 0, 0], [0, 0, b00, b01], [0, 0, b10, b11]]

    AUTHORS:
    - Edinah K. Gnang
    """
    if A.order()==B.order():
        # Initialization of the output hypermatrix
        #T=apply(HM,[A.n(i)+B.n(i) for i in range(A.order())]+['zero'])
        T=HM(*([A.n(i)+B.n(i) for i in range(A.order())]+['zero']))
        # Initialization of the list specifying the dimensions of A
        la = [A.n(i) for i in range(A.order())]
        # Main loop performing the transposition of the entries
        for i in range(prod(la)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entry = [Integer(mod(i,la[0]))]
            sm = Integer(mod(i,la[0]))
            for k in range(len(la)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(la[0:k+1])),la[k+1])))
                sm = sm+prod(la[0:k+1])*entry[len(entry)-1]
            T[tuple(entry)]=A[tuple(entry)]
        # Initialization of the list specifying the dimensions of B
        lb = [B.n(j) for j in range(A.order())]
        # Main loop performing the transposition of the entries
        for j in range(prod(lb)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entry = [Integer(mod(j,lb[0]))]
            sm = Integer(mod(j,lb[0]))
            for k in range(len(lb)-1):
                entry.append(Integer(mod(Integer((j-sm)/prod(lb[0:k+1])),lb[k+1])))
                sm = sm+prod(lb[0:k+1])*entry[len(entry)-1]
            T[tuple((Matrix(ZZ,entry)+Matrix(ZZ,A.dimensions())).list())]=B[tuple(entry)]
        return T
    else:
        raise ValueError("The order of the input hypermatrices must match.")

def GeneralHypermatrixKroneckerProduct(A,B):
    """
    Computes the  Kronecker product of arbitrary hypermatrices A, B of the same order.

    EXAMPLES:

    ::

        sage: GeneralHypermatrixKroneckerProduct(HM(2,2,'a'), HM(2,2,'b'))
        [[a00*b00, a00*b01, a01*b00, a01*b01], [a00*b10, a00*b11, a01*b10, a01*b11], [a10*b00, a10*b01, a11*b00, a11*b01], [a10*b10, a10*b11, a11*b10, a11*b11]]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    if A.order()==B.order():
        # Initialization of the output hypermatrix
        T=HM(*([A.n(i)*B.n(i) for i in range(A.order())]+['zero']))
        # Initialization of the list specifying the dimensions of A and B
        la = [A.n(i) for i in range(A.order())]; lb = [B.n(j) for j in range(B.order())]
        # Main loop performing the transposition of the entries
        for i in range(prod(la)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entrya = [Integer(mod(i,la[0]))]
            sma = Integer(mod(i,la[0]))
            for ka in range(len(la)-1):
                entrya.append(Integer(mod(Integer((i-sma)/prod(la[0:ka+1])),la[ka+1])))
                sma = sma+prod(la[0:ka+1])*entrya[len(entrya)-1]
            # Main loop performing the transposition of the entries
            for j in range(prod(lb)):
                # Turning the index i into an hypermatrix array location using the decimal encoding trick
                entryb = [Integer(mod(j,lb[0]))]
                smb = Integer(mod(j,lb[0]))
                for kb in range(len(lb)-1):
                    entryb.append(Integer(mod(Integer((j-smb)/prod(lb[0:kb+1])),lb[kb+1])))
                    smb = smb+prod(lb[0:kb+1])*entryb[len(entryb)-1]
                T[tuple([lb[z]*Integer(entrya[z])+Integer(entryb[z]) for z in range(A.order())])]=A[tuple(entrya)]*B[tuple(entryb)]
        return T
    else:
        raise ValueError("The order of the input hypermatrices must match.")

def T2Pre(expr):
    """
    Converts formula written in the bracket tree encoding to the Prefix string encoding notation
    the symbol m will stand for the input -1


    EXAMPLES:

    The function implemented here tacitly assume that the input is valid

    ::

        sage: T2Pre(['+',1,1])
        '+11'


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    s = str(expr)
    return ((((s.replace("[","")).replace("]","")).replace(",","")).replace("'","")).replace(" ","").replace('-1','m')


def T2P(expr):
    """
    Converts the formula tree to Postfix notation


    EXAMPLES:
    The tacitly assume that the input is valid

    ::

        sage: T2P(['+',1,1])
        '11+'


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    s = str(expr)
    return ((((s.replace("[","")).replace("]","")).replace(",","")).replace("'","")).replace(" ","")[::-1].replace('-1','m')

def RollLD(L):
    """
    Given an Loaded die, L, the procedures rolls it
    the function output numbers from 1 to len(L)
    each weighted by the 


    EXAMPLES:
    The tacitly assume that the input list is made up of positve integers

    ::

        sage: RollLD([1])
        1


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    N = sum(L)
    r = randint(1,N)
    for i in rg(len(L)):
        if sum(L[:i+1]) >= r:
            return 1+i

@cached_function
def FaT(n):
    """
    The list of formula-binary trees only using addition gates
    which evaluates to the input integer n.

    EXAMPLES:
    The input n must be greater than 0

    ::

        sage: FaT(3)
        [['+', 1, ['+', 1, 1]], ['+', ['+', 1, 1], 1]]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n == 1:
        return [1]
    else :
        gu = []
        for i in range(1,n):
            gu = gu + [['+', g1, g2] for g1 in FaT(i) for g2 in FaT(n-i)]
        return gu

@cached_function
def FaPre(n):
    """
    The list of formula only using addition gates
    which evaluates to the input integer n in prefix notation.

    EXAMPLES:
    The input n must be greater than 0

    ::

        sage: FaPre(3)
        ['+1+11', '++111']


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return [T2Pre(g) for g in FaT(n)]

@cached_function
def FaP(n):
    """
    The list of formula only using addition gates
    which evaluates to the input integer n in postfix notation.

    EXAMPLES:
    The input n must be greater than 0

    ::

        sage: FaP(3)
        ['11+1+', '111++']


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return [T2P(g) for g in FaT(n)]

@cached_function
def Ca(n):
    """
    Outputs the number of formula-binary trees only using addition gates.

    EXAMPLES:
    The input n must be greater than 0

    ::

        sage: Ca(3)
        2


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n == 1:
        return 1
    else :
        return sum([Ca(i)*Ca(n-i) for i in range(1,n)])

@cached_function
def LopFaT(n):
    """
    Outputs all the formula-binary trees only using addition
    but the first term of the addition is >= the second term

    EXAMPLES:
    The input n must be greater than 0

    ::

        sage: LopFaT(3)
        [['+', ['+', 1, 1], 1]]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n == 1:
        return [1]
    else :
        gu = []
        for i in range(1, 1+floor(n/2)):
            gu = gu + [['+', g1, g2] for g1 in LopFaT(n-i) for g2 in LopFaT(i)]
        return gu

@cached_function
def LopCa(n):
    """
    Outputs the number of formula-binary trees only using addition gates
    such that the first term of the addition is >= the second term.
    EXAMPLES:
    The input n must be greater than 0

    ::

        sage: LopCa(3)
        1


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n == 1:
        return 1
    else :
        return sum([LopCa(i)*LopCa(n-i) for i in range(1,1+floor(n/2))])

def RaFaT(n):
    """
    Outputs a uniformly randomly chosen formula-binary tree
    which evaluate to the input integer n.
    EXAMPLES:
    The input n must be greater than 0

    ::

        sage: RaFaT(2)
        ['+', 1, 1]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n == 1:
        return 1
    else :
        # Rolling the Loaded Die.
        j = RollLD([Ca(i)*Ca(n-i) for i in range(1,n+1)])
        return ['+', RaFaT(j), RaFaT(n-j)]

def RaFaPre(n):
    """
    Outputs a uniformly randomly chosen formula-binary tree
    which evaluate to the input integer n in Prefix notation.

    EXAMPLES:
    The input n must be greater than 0

    ::

        sage: RaFaPre(2)
        '+11'


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return(T2Pre(RaFaT(n)))

def RaFaP(n):
    """
    Outputs a uniformly randomly chosen formula-binary tree
    which evaluate to the input integer n in Prefix notation.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaFaP(2)
        '11+'


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return(T2P(RaFaT(n)))

def RaLopFaT(n):
    """
    Outputs a uniformly randomly chosen formula-binary tree
    which evaluate to the input integer n such that the first
    term of the addition is >= the second term.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaLopFaT(2)
        ['+', 1, 1]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n == 1:
        return 1
    else:
        # Rolling the Loaded Die.
        j = RollLD([LopCa(i)*LopCa(n-i) for i in range(1,2+floor(n/2))])
        return ['+', RaLopFaT(n-j), RaLopFaT(j)]

def RaLopFaPre(n):
    """
    Outputs a uniformly randomly chosen formula-binary tree
    which evaluate to the input integer n such that the first
    term of the addition is >= the second term in prefix notation.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaLopFaPre(2)
        '+11'


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return T2Pre(RaLopFaT(n))

def RaLopFaP(n):
    """
    Outputs a uniformly randomly chosen formula-binary tree
    which evaluate to the input integer n such that the first
    term of the addition is >= the second term in poistfix notation.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaLopFaP(2)
        '11+'


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return T2P(RaLopFaT(n))

def LopFaPre(n):
    """
    Outputs all the formula-binary tree
    which evaluate to the input integer n such that the first
    term of the addition is >= the second term in prefix notation.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: LopFaPre(2)
        ['+11']


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return [T2Pre(f) for f in LopFaT(n)]

def LopFaP(n):
    """
    Outputs all the formula-binary tree
    which evaluate to the input integer n such that the first
    term of the addition is >= the second term in prefix notation.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: LopFaP(2)
        ['11+']


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return [T2P(f) for f in LopFaT(n)]

def FamTa(n):
    """
    The list of formula-binary trees only using addition  and
    multiplication gates with the top gate being an addition
    gate which evaluates to the input integer n.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: FamTa(2)
        [['+', 1, 1]]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n == 1:
        return [1]
    else :
        gu = []
        for i in range(1,n):
            gu = gu + [['+', g1, g2] for g1 in FamT(i) for g2 in FamT(n-i)]
        return gu

@cached_function
def FamTm(n):
    """
    The list of formula-binary trees only using addition  and
    multiplication gates with the top gate being a multiplication
    gate which evaluates to the input integer n.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: FamTm(4)
        [['*', ['+', 1, 1], ['+', 1, 1]]]
        

    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n == 1:
        return []
    else :
        gu = []
        for i in range(2,1+floor(n/2)):
            if mod(n,i) == 0:
                gu = gu + [['*', g1, g2] for g1 in FamT(i) for g2 in FamT(n/i)]
        return gu

def FamT(n):
    """
    The list of formula-binary trees only using addition and
    multiplication gates.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: FamT(3)
        [['+', 1, ['+', 1, 1]], ['+', ['+', 1, 1], 1]]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return FamTa(n) + FamTm(n)

@cached_function
def Cama(n):
    """
    Output the size of the set of formulas produced by the procedure FamTa(n).

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: Cama(6)
        48


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n==1:
        return 1
    else:
        return sum([Cam(i)*Cam(n-i) for i in range(1,n)])

@cached_function
def Camm(n):
    """
    Output the size of the set of formulas produced by the procedure FamTa(n).

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: Camm(6)
        4


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return sum([Cam(i)*Cam(n/i) for i in range(2,1+floor(n/2)) if mod(n,i)==0])

@cached_function
def Cam(n):
    """
    Output the size of the set of formulas produced by the procedure FamTa(n).

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: Cam(2)
        1


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return Cama(n)+Camm(n)

def RaFamTa(n):
    """
    Outputs a formula-binary tree sampled uniformly at random
    which evaluates to the input integer n using only addition
    and multiplication gates.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaFamT(2)
        ['+', 1, 1]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n==1:
        return 1
    else:
        j = RollLD([Cam(i)*Cam(n-i) for i in range(1,n)])
        return ['+', RaFamT(j), RaFamT(n-j)]

def RaFamTm(n):
    """
    Outputs a formula-binary tree sampled uniformly at random
    which evaluates to the input integer n using only addition
    and multiplication gates.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaFamTm(4)
        ['*', ['+', 1, 1], ['+', 1, 1]]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if not is_prime(n):
        lu = []
        L  = []
        for i in range(2,1+floor(n/2)):
            if mod(n,i)==0:
                lu.append(i)
                L.append(Cam(i)*Cam(n/i))
        j = RollLD(L)
        return ['*', RaFamT(lu[j-1]), RaFamT(n/lu[j-1])]

def RaFamT(n):
    """
    Outputs a formula-binary tree sampled uniformly at random
    which evaluates to the input integer n using only addition
    and multiplication gates.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaFamT(2)
        ['+', 1, 1]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n==1:
        return 1
    else:
        i = RollLD([Cama(n),Camm(n)])
        if i==1:
            return RaFamTa(n)
        else :
            return RaFamTm(n)

@cached_function
def FamP(n):
    """
    Outputs the set of formula-binary tree written in postfix notation
    which evaluates to the input integer n using only addition
    and multiplication gates.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: FamP(3)
        ['11+1+', '111++']


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return [T2P(f) for f in FamT(n)]

@cached_function
def FamPre(n):
    """
    Outputs the set of formula-binary tree written in prefix notation
    which evaluates to the input integer n using only addition
    and multiplication gates.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: FamPre(3)
        ['+1+11', '++111']


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return [T2Pre(f) for f in FamT(n)]

def RaFamP(n):
    """
    Outputs a uniformly randomly sample formula-binary tree written
    in postfix notation which evaluates to the input integer n using
    only addition and multiplication gates.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaFamP(2)
        '11+'


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return T2P(RaFamT(n))

def RaFamPre(n):
    """
    Outputs a uniformly randomly sample formula-binary tree written
    in prefix notation which evaluates to the input integer n using
    only addition and multiplication gates.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaFamPre(2)
        '+11'


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return T2Pre(RaFamT(n))

@cached_function
def FameTa(n):
    """
    The list of formula-binary trees only using addition,
    multiplication, and exponentiation gates. The top gate
    being an addition gate and and the formula evaluates to
    the input integer n.

    EXAMPLES:

    The input n must be greater than 0

    ::
        sage: FameTa(3)
        [['+', 1, ['+', 1, 1]], ['+', ['+', 1, 1], 1]]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n == 1:
        return [1]
    else:
        gu = []
        for i in range(1,n):
            gu = gu + [['+', g1, g2] for g1 in FameT(i) for g2 in FameT(n-i)]
        return gu

@cached_function
def FameTm(n):
    """
    The list of formula-binary trees only using addition.
    multiplication and exponentiation gates with the top
    gate being a multiplication gate which evaluates to the
    input integer n.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: FameTm(6)
        [['*', ['+', 1, 1], ['+', 1, ['+', 1, 1]]],
         ['*', ['+', 1, 1], ['+', ['+', 1, 1], 1]],
         ['*', ['+', 1, ['+', 1, 1]], ['+', 1, 1]],
         ['*', ['+', ['+', 1, 1], 1], ['+', 1, 1]]]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n == 1:
        return []
    else :
        gu = []
        for i in range(2,1+floor(n/2)):
            if mod(n,i) == 0:
                gu = gu + [['*', g1, g2] for g1 in FameT(i) for g2 in FameT(n/i)]
        return gu

@cached_function
def FameTe(n):
    """
    The list of formula-binary trees only using addition.
    multiplication and exponentiation gates with the top
    gate being an exponetiation gate which evaluates to the
    input integer n.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: FameTe(4)
        [['^', ['+', 1, 1], ['+', 1, 1]]]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n == 1:
        return []
    else :
        gu = []
        for i in range(2,2+floor(log(n)/log(2))):
            if floor(n^(1/i)) == ceil(n^(1/i)):
                gu = gu + [['^', g1, g2] for g1 in FameT(i) for g2 in FameT(n^(1/i))]
        return gu

@cached_function
def FameT(n):
    """
    The list of formula-binary trees only using addition.
    multiplication and exponentiation gates which evaluates to the
    input integer n.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: FameT(3)
        [['+', 1, ['+', 1, 1]], ['+', ['+', 1, 1], 1]]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return FameTa(n) + FameTm(n) + FameTe(n)

@cached_function
def Camea(n):
    """
    Output the size of the set of formulas produced by the procedure FameTa(n).

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: Camea(6)
        54


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n==1:
        return 1
    else:
        return sum([Came(i)*Came(n-i) for i in range(1,n)])

@cached_function
def Camem(n):
    """
    Output the size of the set of formulas produced by the procedure FameTm(n).

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: Camem(6)
        4


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return sum([Came(i)*Came(n/i) for i in range(2,1+floor(n/2)) if mod(n,i)==0])

@cached_function
def Camee(n):
    """
    Output the size of the set of formulas produced by the procedure FameTe(n).

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: Camee(9)
        2


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return sum([Came(i)*Came(n^(1/i)) for i in range(2,2+floor(log(n)/log(2))) if floor(n^(1/i)) == ceil(n^(1/i))])

@cached_function
def Came(n):
    """
    Output the size of the set of formulas produced by the procedure FameT(n).
    Which counts all monotone formula encodings evaluating using a combination
    of fanin two addition, multiplication and exponentiation gates which evaluates
    to the input integer n

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: Came(6)
        58


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return Camea(n)+Camem(n)+Camee(n)

def RaFameTa(n):
    """
    Output a Random Formula Tres chosen uniformly at random amoung all
    Tree representation of the positive integer n which use a combination
    of addition, multiplication and exponentiation gates with the top gate
    being the addition gate

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaFameTa(2)
        ['+', 1, 1]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n == 1:
        return 1
    else:
        j = RollLD([Came(i)*Came(n-i) for i in range(1,n)])
        return ['+', RaFameT(j), RaFameT(n-j)]
        
def RaFameTm(n):
    """
    Outputs a formula-binary tree sampled uniformly at random
    which evaluates to the input integer n using only addition
    and multiplication and exponentiaiton gates, with the top
    gate being a multiplication gate.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaFameT(2)
        ['+', 1, 1]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if not is_prime(n):
        lu = []
        L  = []
        for i in range(2,1+floor(n/2)):
            if mod(n,i)==0:
                lu.append(i)
                L.append(Came(i)*Came(n/i))
        j = RollLD(L)
        return ['*', RaFameT(lu[j-1]), RaFameT(n/lu[j-1])]

def RaFameTe(n):
    """
    Outputs a formula-binary tree sampled uniformly at random
    which evaluates to the input integer n using only addition
    and multiplication and exponentiaiton gates, with the top
    gate being an exponentiation gate.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaFameTe(4)
        ['^', ['+', 1, 1], ['+', 1, 1]]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if not is_prime(n) and n>1:
        lu = []
        L  = []
        for i in range(2,1+floor(n/2)):
            if floor(n^(1/i)) == ceil(n^(1/i)):
                lu.append(i)
                L.append(Came(i)*Came(n^(1/i)))
        j = RollLD(L)
        return ['^', RaFameT( n^(1/lu[j-1]) ), RaFameT(lu[j-1])]

def RaFameT(n):
    """
    Outputs a formula-binary tree sampled uniformly at random
    which evaluates to the input integer n using only addition
    multiplication and exponentiation gates.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaFameT(2)
        ['+', 1, 1]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n==1:
        return 1
    else:
        i = RollLD([Camea(n),Camem(n),Camee(n)])
        if i==1:
            return RaFameTa(n)
        elif i==2:
            return RaFameTm(n)
        else :
            return RaFameTe(n)

def RaFameP(n):
    """
    Outputs a uniformly randomly sample formula-binary tree written
    in postfix notation which evaluates to the input integer n using
    only addition, multiplication and exponentiation gates.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaFameP(2)
        '11+'


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return T2P(RaFameT(n))

def RaFamePre(n):
    """
    Outputs a uniformly randomly sample formula-binary tree written
    in prefix notation which evaluates to the input integer n using
    only addition, multiplication  and exponentiation gates.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: RaFamePre(2)
        '+11'


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return T2Pre(RaFameT(n))

def Tsize(T):
    """
    Outputs the size of the Tree associated with the formula.
 
    EXAMPLES:
    
    ::

        sage: Tsize(['+', 1, 1])
        3


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if T==1:
        return 1
    elif T==-1:
        return 1
    else:
        return 1+Tsize(T[1])+Tsize(T[2])

@cached_function
def ShortestTame(n):
    """
    Outputs the length and an example of the smallest binary-tree
    formula using addition, multiplication and exponentiation

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: ShortestTame(2)
        [3, ['+', 1, 1]]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n==1:
        return [1,1]
    else:
        aluf = []
        si = 2*n
        for i in range(1,n):
            T1 = ShortestTame(i)
            T2 = ShortestTame(n-i)
            if (T1[0]+T2[0]+1) < si:
                si = T1[0]+T2[0]+1
                if EvalT(T1[1]) <= EvalT(T2[1]):
                    aluf = ['+', T1[1], T2[1]]
                else:
                    aluf = ['+', T2[1], T1[1]]
        for i in range(2,floor(n/2)):
            if mod(n,i)==0:
                T1 = ShortestTame(i)
                T2 = ShortestTame(n/i)
                if (T1[0]+T2[0]+1) < si:
                    si = T1[0]+T2[0]+1
                    if EvalT(T1[1]) <= EvalT(T2[1]):
                        aluf = ['*', T1[1], T2[1]]
                    else:
                        aluf = ['*', T2[1], T1[1]]
        for i in range(2,2+floor(log(n)/log(2))):
            if floor(n^(1/i)) == ceil(n^(1/i)):
                T1 = ShortestTame(n^(1/i))
                T2 = ShortestTame(i)
                if (T1[0]+T2[0]+1) < si:
                    si = T1[0]+T2[0]+1
                    aluf = ['^', T1[1], T2[1]]
        return [si, aluf]

@cached_function
def ShortestTameList(n):
    """
    Outputs the list of the smallest binary-tree
    formula using addition, multiplication and exponentiation

    EXAMPLES:

    ::

        sage: ShortestTameList(4)
        [['+', 1, ['+', 1, ['+', 1, 1]]], ['+', 1, ['+', ['+', 1, 1], 1]], ['+', ['+', 1, 1], ['+', 1, 1]], ['+', ['+', 1, ['+', 1, 1]], 1], ['+', ['+', ['+', 1, 1], 1], 1], ['*', ['+', 1, 1], ['+', 1, 1]], ['^', ['+', 1, 1], ['+', 1, 1]]]


    AUTHORS:
    - Edinah K. Gnang

    To Do :
    - Try to implement faster version of this procedure

    """
    # Obtaining the minimal size
    si=ShortestTame(n)[0]
    return [f for f in FameT(n) if Tsize(f)==si]

def get_permutation(la,lb):
    """
    Obtains a permutation list from two lists of the same size.
    No check is performed here the user must be very carefull
    to input lists of the same size
    
    EXAMPLES:

    ::

        sage: get_permutation([1, 2, 3, 4, 6, 8],[1, 2, 4, 8, 3, 6])
        [0, 1, 4, 2, 5, 3]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing the output
    L = list()

    # Loop performing the evaluation.
    for i1 in range(len(la)):
        for i2 in range(len(lb)):
            if la[i1] == lb[i2]:
                L.append(i2)
                break
    return L

def permute(l,p):
    """
    Permutes the entries of the list l according to the permutation p
    No check is performed here the user must be very carefull
    to input lists of the same size
 
    EXAMPLES:

    ::

        sage: permute([1, x, x^x, x^(x + 1), x + 1, (x + 1)*x],[0, 1, 4, 2, 5, 3])
        [1, x, x + 1, x^x, (x + 1)*x, x^(x + 1)]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing the output
    L = list()
    # Loop performing the evaluation.
    for i in range(len(l)):
        L.append(l[p[i]])
    return L

@cached_function
def NaiveZetaT(nbit):
    """
    Produces Tree associated with the Second Canonical Forms.

    EXAMPLES:

    ::

        sage: NaiveZetaT(1)
        [[['+', 1, 1]], [1, ['+', 1, 1]]]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 

    """

    # Initial conditions
    Pk = [['+',1,1]]
    Nk = [1] + Pk
    for it in range(1,nbit):
        L = []
        for p in Pk:
            if L == []:
                L = [1] + [p] + [['^',p,Nk[i]] for i in range(1,len(Nk))]

            else:
                Lp = [p] + [['^',p,Nk[i]] for i in range(1,len(Nk))]
                L  = L + [['*',L[i],n] for i in range(1,len(L)) for n in Lp] + Lp

        # Sorting the list
        Va = [EvalT(l) for l in L]
        Vb = copy(Va)
        Vb.sort()
        perm = get_permutation(Vb, Va)
        # Reinitialization of the list Nk
        Nk = permute(L, perm)
        # Set completion
        l = len(Nk)
        i = 0
        while i < l-1:
            if EvalT(Nk[i+1]) - EvalT(Nk[i]) == 2 :
                Pk.append(['+',1,Nk[i]])
                Nk.insert(i+1, ['+',1,Nk[i]])
                l = l+1
            else:
                i = i+1
    return [Pk, Nk]


@cached_function
def Goodstein(number_of_iterations=1):
    """
    Produces the set of symbolic expressions associated with the
    the first canonical form. In all the expressions the symbolic
    variable x stands for a short hand notation for the formula (1+1).

    EXAMPLES:

    ::

        sage: Goodstein(1)[:3]
        [1, x, x + 1]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initial condition of Initial set
    Ng0 = [1, x]
    # Main loop performing the iteration
    for iteration in range(number_of_iterations):
        # Implementation of the set recurrence
        Ng0 = [1] + [x^n for n in Ng0]
        # Initialization of a buffer list Ng1
        # which will store updates to Ng0
        Ng1 = []
        for n in Set(Ng0).subsets():
            if n.cardinality() > 0:
                Ng1.append(sum(n))
        Ng0 = list(Ng1)
    Nf = []
    for i in range(len(Ng0)):
        Nf.append([])
    for i in range(len(Nf)):
        Nf[(Ng0[i]).subs(x=2)-1].append(Ng0[i])
    Ng0=[]
    for i in range(len(Nf)):
        Ng0.append(Nf[i][0])
    return Ng0

@cached_function
def GoodsteinT(number_of_iterations=1):
    """
    Produces Tree associated with Goodstein Trees.

    ::

        sage: GoodsteinT(1)
        [1,
         ['+', 1, 1],
         ['+', 1, ['+', 1, 1]],
         ['^', ['+', 1, 1], ['+', 1, 1]],
         ['+', 1, ['^', ['+', 1, 1], ['+', 1, 1]]],
         ['+', ['+', 1, 1], ['^', ['+', 1, 1], ['+', 1, 1]]],
         ['+', ['+', 1, ['+', 1, 1]], ['^', ['+', 1, 1], ['+', 1, 1]]]]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initial condition of Initial set
    Ng0 = [1, ['+',1,1]]
    # Main loop performing the iteration
    for iteration in range(number_of_iterations):
        # Implementation of the set recurrence
        Tmp = copy(Ng0)
        Tmp.pop(0)
        Ng0 = [1] + [['+',1,1]] + [['^',['+',1,1], m] for m in Tmp]
        # Initialization of a buffer list Ng1
        # which will store updates to Ng0
        Ng1 = []
        for n in Set(range(len(Ng0))).subsets():
            if n.cardinality() == 1:
                Ng1.append(Ng0[n[0]])
            elif n.cardinality() > 1:
                T = Ng0[n[0]]
                for j in range(1,n.cardinality()):
                    T = ['+', T, Ng0[n[j]]]
                Ng1.append(T)
        Ng0 = copy(Ng1)
    # Sorting the obtained list
    Nf = []
    for i in range(len(Ng0)):
        Nf.append([])
    for i in range(len(Nf)):
        Nf[EvalT(Ng0[i])-1].append(Ng0[i])
    Ng0=[]
    for i in range(len(Nf)):
        Ng0.append(Nf[i][0])
    return Ng0

def list_eval(L):
    """
    Perform the evaluation of the list to integers
 
    EXAMPLES:

    ::

        sage: x=var('x'); list_eval([1, x, x+1, x^x])
        [1, 2, 3, 4]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    return [i.substitute(x=2) for i in L]

def get_permutation(la,lb):
    """
    Obtains a permutation list from two lists of the same size.
    No check is performed here the user must be very carefull
    to input lists of the same size
    
 
    EXAMPLES:

    ::

        sage: get_permutation([1, 2, 3, 4, 6, 8], [1, 2, 4, 8, 3, 6])
        [0, 1, 4, 2, 5, 3]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing the output
    L = list()
    # Loop performing the evaluation.
    for i1 in range(len(la)):
        for i2 in range(len(lb)):
            if la[i1] == lb[i2]:
                L.append(i2)
                break
    return L

def permute(l,p):
    """
    Permutes the entries of the list l according to the permutation p
    No check is performed here the user must be very carefull
    to input lists of the same size
 
    EXAMPLES:

    ::

        sage: permute([1, x, x^x, x^(x + 1), x + 1, (x + 1)*x],[0, 1, 4, 2, 5, 3])
        [1, x, x + 1, x^x, (x + 1)*x, x^(x + 1)]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing the output
    L = list()
    # Loop performing the evaluation.
    for i in range(len(l)):
        L.append(l[p[i]])
    return L

@cached_function
def base2expansion(n):
    """
    Returns the polynomial encoding the binary expansion
     
    EXAMPLES:

    The function is a crucial first step to the recursive encoding
    
    ::

        sage: p = base2expansion(10)
        sage: p
        x^3 + x

    AUTHORS:
    - Edinah K. Gnang
    -  
    """
    x = var('x')
    # polynomial
    p = 0
    k = 2
    if n == 1:
        return 1
    elif n > 1:
        while k < n:
            k = k^2
        if k == n:
            return x^(log(k,2))
        elif k > n:
            k = sqrt(k)
            while k < n:
                k = 2*k
                if k == n:
                    return x^(log(k,2))
                elif k > n:
                    p = x^(floor(log(k/2,2))) + base2expansion(n-k/2)
    return p

@cached_function
def base2expansionT(n):
    """
    Returns the polynomial encoding the binary expansion
     
    EXAMPLES: 
    The function is a crucial first step to the recursive encoding
    
    ::

        sage: base2expansionT(2)
        ['^', ['+', 1, 1], 1]


    AUTHORS:
    - Edinah K. Gnang
    -  
    """
    x = var('x')
    # polynomial
    p = 0; k = 2
    if n == 1:
        return 1
    elif n > 1:
        while k < n:
            k = k^2
        if k == n:
            return ['^',['+',1,1],log(k,2)]
        elif k > n:
            k = sqrt(k)
            while k < n:
                k = 2*k
                if k == n:
                    return ['^',['+',1,1],log(k,2)]
                elif k > n:
                    p = ['+',['^',['+',1,1],(floor(log(k/2,2)))],base2expansionT(n-k/2)]
    return p

@cached_function
def recurse_base2expansion(n):
    """
    Returns the Goodstein encoding of the input integer n
     
    EXAMPLES: 
    The function builds the crucial base2expansion function
    however it should be noted that it's a horrible idea
    to use this function to compute the recursive encoding
    for a list of consecutive integer the number of recursive
    call is unmanageable
    
    ::

        sage: recurse_base2expansion(4)
        x^x
        

    AUTHORS:
    - Edinah K. Gnang
    -  
    """
    p = 0; k = 2
    if n == 1:
        return 1
    elif n > 1:
        while k < n:
            k = k^2
        if k == n:
            return x^recurse_base2expansion(log(k,2))
        elif k > n:
            k = sqrt(k)
            while k < n:
                k = 2*k
                if k == n:
                    return x^recurse_base2expansion(log(k,2))
                elif k > n:
                    p = x^recurse_base2expansion(floor(log(k/2,2))) + recurse_base2expansion(n-k/2)
    return p

@cached_function
def recurse_base2expansionT(n):
    """
    Returns the Goodstein Tree encoding of the input integer n
     
    EXAMPLES: 
    The function builds the crucial base2expansion function
    however it should be noted that it's a horrible idea
    to use this function to compute the recursive encoding
    for a list of consecutive integer the number of recursive
    call is unmanageable
    
    ::

        sage: recurse_base2expansionT(2)
        ['^', ['+', 1, 1], 1]


    AUTHORS:
    - Edinah K. Gnang
    -  
    """
    p = 0; k = 2
    if n == 1:
        return 1
    elif n > 1:
        while k < n:
            k = k^2
        if k == n:
            return ['^',['+',1,1],recurse_base2expansionT(log(k,2))]
        elif k > n:
            k = sqrt(k)
            while k < n:
                k = 2*k
                if k == n:
                    return ['^',['+',1,1],recurse_base2expansionT(log(k,2))]
                elif k > n:
                    p = ['+',['^',['+',1,1],recurse_base2expansionT(floor(log(k/2,2)))],recurse_base2expansionT(n-k/2)]
    return p

@cached_function
def Fa3T(n):
    """
    The list of formula-binary trees only using fan-in three addition gates
    which evaluates to the input integer n.

    EXAMPLES:

    The input n must be greater than 0

    ::

        sage: Fa3T(3)
        [['+', 1, 1, 1]]

    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n == 1:
        return [1]
    else :
        gu = []
        for i in range(1,n,2):
            for j in range(1,n-i,2):
                gu = gu + [['+', g1, g2, g3] for g1 in Fa3T(i) for g2 in Fa3T(j) for g3 in Fa3T(n-i-j)]
        return gu

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

def Zeta(nbitr):
    """
    Produces Tree associated with the Second Canonical Forms.
    Implements an improved version of the zeta recurrence and
    the combinatorial tower sieve.

    EXAMPLES:

    ::

        sage: Zeta(1)
        The current iteration will uncover 1.00000000000000 new primes in the range [2.00000000000000, 4.00000000000000]
        [[x, x + 1], [1, x, x + 1, x^x]]

    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger
    
    To Do :
    - Try to implement faster version of this procedure

    """
    x = var('x')
    # Pr corresponds to the initial list of primes
    Pr = [x]
    # Nu corresponds to the initial list of integer
    NuC  = [1,x]; TNuC = [1,x]
    # Initializing the upper and lower bound
    upr_bnd = 2^2; lwr_bnd = 2
    # Computing the set recurrence
    for itr in range(nbitr):
        for jtr in range(log(upr_bnd,2)-log(lwr_bnd,2)):
            TpNu = [1]
            for p in Pr:
                TpNu=TpNu+[m*pn for m in TpNu for pn in [p^n for n in NuC if (p^n).subs(x=2)<=2^(log(lwr_bnd,2)+jtr+1)] if (m*pn).subs(x=2)<=2^(log(lwr_bnd,2)+jtr+1)]
            # Keeping only the elements within the range of the upper and lower bound
            Nu = [f for f in TpNu if (2^(log(lwr_bnd,2)+jtr)<f.subs(x=2) and f.subs(x=2)<=2^(log(lwr_bnd,2)+jtr+1))]
            print('The current iteration will uncover '+str(2^(N(log(lwr_bnd,2))+jtr+1)-2^(N(log(lwr_bnd,2))+jtr)-len(Nu))+' new primes in the range ['+str(2^(N(log(lwr_bnd,2))+jtr))+', '+str(2^(N(log(lwr_bnd,2))+jtr+1))+']')
            # Obtaining the corresponding sorted integer list
            la = [f.subs(x=2) for f in Nu]; lb = copy(la); lb.sort()
            # Obtaining the sorting permutation
            perm = []
            for i1 in range(len(la)):
                for i2 in range(len(lb)):
                    if lb[i1]==la[i2]:
                        perm.append(i2)
                        break
            # Sorting the list using the obtained permutation
            Nu = [Nu[perm[j]] for j in range(len(Nu))]
            # Computing the set completion
            TNuC = TNuC + Nu
            l = len(TNuC)
            i = 2^(log(lwr_bnd,2)+jtr-1)
            while i<l-1:
                if(TNuC[i+1].subs(x=2)-TNuC[i].subs(x=2)==2):
                    Pr.append(TNuC[i]+1)
                    TNuC.insert(i+1,TNuC[i]+1)
                    l=l+1
                else:
                    i=i+1
        # Updating the list of integers
        NuC = TNuC
        # Updating the upper and lower bound
        lwr_bnd = upr_bnd
        upr_bnd = 2^upr_bnd
    return [Pr,NuC]        

@cached_function
def ZetaT(nbitr):
    """
    Produces Tree associated with the Second Canonical Forms.
    Implements an improved version of the zeta recurrence and
    the combinatorial tower sieve.

    EXAMPLES:

    ::

        sage: ZetaT(1)
        The current iteration will uncover 1.00000000000000 new primes in the range [2.00000000000000, 4.00000000000000]
        [[['+', 1, 1], ['+', 1, ['+', 1, 1]]],
         [1, ['+', 1, 1], ['+', 1, ['+', 1, 1]], ['^', ['+', 1, 1], ['+', 1, 1]]]]

    AUTHORS:
    - Edinah K. Gnang

    To Do :
    - Try to implement faster version of this procedure

    """
    # Pr corresponds to the initial list of primes
    Pr = [ ['+',1,1] ]
    # Nu corresponds to the initial list of integer
    NuC  = [1] + Pr
    TNuC = [1] + Pr
    # Initializing the upper and lower bound
    upr_bnd = 2^2
    lwr_bnd = 2
    # Computing the set recurrence
    for itr in range(nbitr):
        for jtr in range(log(upr_bnd,2)-log(lwr_bnd,2)):
            TpNu = [1]
            for p in Pr:
                TpNu = TpNu+[pn for pn in [['^',p,n] for n in NuC if EvalT(['^',p,n])<=2^(log(lwr_bnd,2)+jtr+1)]]+[['*',m,pn] for m in TpNu[1:] for pn in [['^',p,n] for n in NuC if EvalT(['^',p,n])<=2^(log(lwr_bnd,2)+jtr+1)]  if EvalT(['*',m,pn])<=2^(log(lwr_bnd,2)+jtr+1)]
            # Keeping only the elements within the range of the upper and lower bound
            Nu = [f for f in TpNu if (2^(log(lwr_bnd,2)+jtr)<EvalT(f) and EvalT(f)<=2^(log(lwr_bnd,2)+jtr+1))] 
            print('The current iteration will uncover '+str(2^(N(log(lwr_bnd,2))+jtr+1)-2^(N(log(lwr_bnd,2))+jtr)-len(Nu))+' new primes in the range ['+str(2^(N(log(lwr_bnd,2))+jtr))+', '+str(2^(N(log(lwr_bnd,2))+jtr+1))+']')
            # Obtaining the corresponding sorted integer list
            la = [EvalT(f) for f in Nu]; lb = copy(la); lb.sort()
            # Obtaining the sorting permutation
            perm = []
            for i1 in range(len(la)):
                for i2 in range(len(lb)):
                    if lb[i1]==la[i2]:
                        perm.append(i2)
                        break
            # Sorting the list using the obtained permutation
            Nu = [Nu[perm[j]] for j in range(len(Nu))]
            # Perfoming the set completion
            TNuC = TNuC + Nu
            l = len(TNuC)
            i = 2^(log(lwr_bnd,2)+jtr-1)
            while i<l-1:
                if(EvalT(TNuC[i+1])-EvalT(TNuC[i])==2):
                    Pr.append(['+',1,TNuC[i]])
                    TNuC.insert(i+1,['+',1,TNuC[i]])
                    l=l+1
                else:
                    i=i+1
        # Updating the list of integers
        NuC = TNuC
        # Updating the upper and lower bound
        lwr_bnd = upr_bnd
        upr_bnd = 2^upr_bnd
    return [Pr,NuC]        

def Horner(nbitr):
    """
    Produces list of symbolic expressions associated with recursive horner encoding.

    EXAMPLES:

    ::

        sage: Horner(2)
        [1, x, x + 1, x^x, (x + 1)*x, (x + 1)*x^x, x^(x^x), x^(x + 1), x^x + 1, (x^x + 1)*x, (x^x + 1)*x^x, (x^x + 1)*x^(x^x), x^(x + 1)*(x^x + 1), x^((x + 1)*x), x^((x + 1)*x^x), x^(x^(x^x)), x^(x^(x + 1)), x^(x^x + 1), (x + 1)*x + 1, (x + 1)*x^x + 1, x^(x^x) + 1, x^(x + 1) + 1]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 

    """
    x = var('x')
    Nk  = [1, x, 1+x, x^x]
    # Initialization of the lists
    LEk = [x^x]; LOk = [1+x]; LPk = [x, x^x]
    # Main loop computing the encoding
    for i in range(nbitr):
        # Updating the list
        LEkp1 = [lp*lo for lp in LPk for lo in LOk] + [x^m for m in LEk+LOk]
        LOkp1 = [n+1 for n in LEk]
        LPkp1 = LPk + [x^m for m in LEk+LOk]
        # The New replaces the old
        Nk = Nk + LEkp1+LOkp1
        LEk = LEkp1; LOk = LOkp1; LPk = LPkp1
    return Nk

def HornerT(nbitr):
    """
    Produces Tree associated with the recursive Horner encoding

    EXAMPLES:

    ::

        sage: HornerT(1)
        [1,
         ['+', 1, 1],
         ['+', 1, ['+', 1, 1]],
         ['^', ['+', 1, 1], ['+', 1, 1]],
         ['*', ['+', 1, 1], ['+', 1, ['+', 1, 1]]],
         ['*', ['^', ['+', 1, 1], ['+', 1, 1]], ['+', 1, ['+', 1, 1]]],
         ['^', ['+', 1, 1], ['^', ['+', 1, 1], ['+', 1, 1]]],
         ['^', ['+', 1, 1], ['+', 1, ['+', 1, 1]]],
         ['+', 1, ['^', ['+', 1, 1], ['+', 1, 1]]]]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 

    """
    # Initial set
    Nk  = [ 1, ['+',1,1], ['+',1,['+',1,1]], ['^',['+',1,1],['+',1,1]] ]
    # Initialization of the lists
    LEk = [ ['^',['+',1,1],['+',1,1]] ]
    LOk = [ ['+',1,['+',1,1]] ]
    LPk = [ ['+',1,1], ['^',['+',1,1],['+',1,1]] ]
    # Main loop computing the recursive horner encoding
    for i in range(nbitr):
        # Updating the list
        LEkp1 = [['*',lp,lo] for lp in LPk for lo in LOk] + [['^',['+',1,1],m] for m in LEk+LOk]
        LOkp1 = [['+',1,n] for n in LEk]
        LPkp1 = LPk + [['^',['+',1,1],m] for m in LEk+LOk]
        # The New replaces the old
        Nk = Nk+LEkp1+LOkp1
        LEk = LEkp1; LOk = LOkp1; LPk = LPkp1
    return Nk

def EvalT(T):
    """
    Outputs the evaluation value of a tree.

    EXAMPLES:

    ::

        sage: EvalT(['+', ['+', 1, ['+', 1, 1]], ['+', ['+', 1, 1], 1]])
        6

    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger
    """
    if T == 1:
        return 1
    elif T == -1:
        return -1
    elif T[0] == '+':
        return EvalT(T[1]) + EvalT(T[2])
    elif T[0] == '*':
        return EvalT(T[1]) * EvalT(T[2])
    elif T[0] == '^':
        return EvalT(T[1]) ^ EvalT(T[2])
    else:
        print('IMPROPER INPUT !!!')

@cached_function
def MonotoneFormula(n):
    """
    Outputs Monotone formula encodings of length at most n.

    EXAMPLES:

    ::

        sage: MonotoneFormula(3)
        [[], [1], [], []]

    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n<=3:
        return [[], [1], [], []]
    elif n>3:
        # Initialization of the list of formula.
        A=[[], [1]] + [[] for t in range(n-1)]
        # Main loop.
        for sz in range(3,n+1):
            # Initialization of the fifth entry
            for o in ['+', '*', '^']:
                for i in range(1,sz-1):
                        A[sz]=A[sz]+[[o,s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0)]
        return A 

@cached_function
def ReducedMonotoneFormula(n):
    """
    Outputs monotone formula encodings of length at most n.

    EXAMPLES:

    ::

        sage: ReducedMonotoneFormula(3)
        [[], [1], [], []]

    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n<=3:
        return [[], [1], [], []]
    elif n>3:
        # Initialization of the list of formula.
        A=[[], [1]] + [[] for t in range(n-1)]
        # Main loop.
        for sz in range(3,n+1):
            # Initialization of the fifth entry
            for i in range(1,sz-1):
                A[sz]=A[sz]+[['+',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and not EvalT(s).is_zero() and not EvalT(t).is_zero() and not EvalT(['+',s,t]).is_zero()]
            for i in range(1,sz-1):
                A[sz]=A[sz]+[['*',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and EvalT(s)!=1 and EvalT(t)!=1 and EvalT(['*',s,t])!=1]
            for i in range(1,sz-1):
                A[sz]=A[sz]+[['^',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and EvalT(s)!=1 and EvalT(t)!=1 and EvalT(['^',s,t])!=1]
        return A

def ReducedMonotoneFormulaSets(sz):
    """
    Outputs set of numbers associated with monotone encodings
    of complexity less then the size input parameter sz.

    EXAMPLES:

    ::

        sage: ReducedMonotoneFormulaSets(5)
        [{}, {1}, {}, {2}, {}, {3}]

    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    Lt=ReducedMonotoneFormula(sz)
    # Filling up the result list
    Rslt=[]
    for n in range(sz+1):
        L=[]; i=0
        while i<len(Lt[n]):
            L.append(Lt[n][i])
            i=i+1
        Rslt.append(Set([EvalT(L[i]) for i in range(len(L))]))
    # Cleaning up the list
    for i in range(1,len(Rslt)):
        for j in range(i):
            Rslt[i]=Rslt[i].difference(Rslt[j])
    return Rslt

@cached_function
def NonMonotoneFormula(n):
    """
    Outputs non-monotone formula encodings of length at most n.

    EXAMPLES:

    ::

        sage: NonMonotoneFormula(3)
        [[], [1, -1], [], []]

    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n<=3:
        return [[], [1,-1], [], []]
    elif n>3:
        # Initialization of the list of formula.
        A=[[], [1,-1]] + [[] for t in range(n-1)]
        # Main loop.
        for sz in range(3,n+1):
            # Initialization of the fifth entry
            for o in ['+', '*', '^']:
                for i in range(1,sz-1):
                        A[sz]=A[sz]+[[o,s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0)]
        return A

@cached_function
def ReducedNonMonotoneFormula(n):
    """
    Outputs non-monotone formula encodings of length at most n.

    EXAMPLES:

    ::

        sage: ReducedNonMonotoneFormula(3)
        [[], [1, -1], [], []]

    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n<=3:
        return [[], [1,-1], [], []]
    elif n>3:
        # Initialization of the list of formula.
        A=[[], [1,-1]] + [[] for t in range(n-1)]
        # Main loop.
        for sz in range(3,n+1):
            # Initialization of the fifth entry
            for i in range(1,sz-1):
                A[sz]=A[sz]+[['+',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and not EvalT(s).is_zero() and not EvalT(t).is_zero() and not EvalT(['+',s,t]).is_zero()]
            for i in range(1,sz-1):
                A[sz]=A[sz]+[['*',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and EvalT(s)!=1 and EvalT(t)!=1 and EvalT(['*',s,t])!=1]
            for i in range(1,sz-1):
                A[sz]=A[sz]+[['^',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and EvalT(s)!=1 and EvalT(t)!=1 and EvalT(['^',s,t])!=1]
        return A

def ReducedNonMonotoneFormulaSets(sz):
    """
    Outputs set of numbers associated with non monotone encodings
    of complexity less then the size input parameter sz.

    EXAMPLES:

    ::

        sage: ReducedNonMonotoneFormulaSets(5)
        [{}, {1, -1}, {}, {2, -2}, {}, {3, -1/2, -3, 1/2}]

    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    Lt=ReducedNonMonotoneFormula(sz)
    # Filling up the result list
    Rslt=[]
    for n in range(sz+1):
        L=[]; i=0
        while i<len(Lt[n]):
            L.append(Lt[n][i])
            i=i+1
        Rslt.append(Set([EvalT(L[i]) for i in range(len(L))]))
    # Cleaning up the list
    for i in range(1,len(Rslt)):
        for j in range(i):
            Rslt[i]=Rslt[i].difference(Rslt[j])
    return Rslt

@cached_function
def NonMonotoneShortestTameList(m):
    """
    Outputs the list of the smallest binary-tree
    formula using addition, multiplication and exponentiation

    EXAMPLES:

    ::

        sage: NonMonotoneShortestTameList(10)
        [['+', 1, ['^', ['+', 1, ['+', 1, 1]], ['+', 1, 1]]],
         ['+', 1, ['^', ['+', -1, ['+', -1, -1]], ['+', 1, 1]]],
         ['+', 1, ['^', ['+', ['+', 1, 1], 1], ['+', 1, 1]]],
         ['+', 1, ['^', ['+', ['+', -1, -1], -1], ['+', 1, 1]]],
         ['+', ['^', ['+', 1, ['+', 1, 1]], ['+', 1, 1]], 1],
         ['+', ['^', ['+', -1, ['+', -1, -1]], ['+', 1, 1]], 1],
         ['+', ['^', ['+', ['+', 1, 1], 1], ['+', 1, 1]], 1],
         ['+', ['^', ['+', ['+', -1, -1], -1], ['+', 1, 1]], 1]]
        

    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    # Obtaining the length of the shortes monotone formula
    n=ShortestTame(m)[0]
    # Initialization of the stop flag
    StopFlag=False
    if n<=3:
        return m
    elif n>3:
        # Initialization of the list of formula.
        A=[[], [1,-1]] + [[] for t in range(n-1)]
        # Main loop.
        for sz in range(3,n+1):
            # Initialization of the fifth entry
            for i in range(1,sz-1):
                A[sz]=A[sz]+[['+',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and not EvalT(s).is_zero() and not EvalT(t).is_zero() and not EvalT(['+',s,t]).is_zero()]
            for i in range(1,sz-1):
                A[sz]=A[sz]+[['*',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and EvalT(s)!=1 and EvalT(t)!=1 and EvalT(['*',s,t])!=1]
            for i in range(1,sz-1):
                A[sz]=A[sz]+[['^',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and EvalT(s)!=1 and EvalT(t)!=1 and EvalT(['^',s,t])!=1]
            if m in [EvalT(T) for T in A[sz]]:
                StopFlag=True
                break
        return [T for T in A[sz] if EvalT(T)==m]

def EvalTn(T, dgts):
    """
    Outputs the evaluation value of a tree using numerically truncated complex numbers
    where the number of digits is specified by the second input
    it is a real mystery how the branches are being chosen however.

    EXAMPLES:

    ::

        sage: EvalTn(['+', ['+', 1, ['+', 1, 1]], ['+', ['+', 1, 1], 1]], 50)
        6.0000000000000

    AUTHORS:
    - Edinah K. Gnang
    """
    if T == ComplexField(dgts)(1,0):
        return ComplexField(dgts)(1,0)
    elif T == ComplexField(dgts)(-1,0):
        return ComplexField(dgts)(-1,0)
    elif T[0] == '+':
        return ComplexField(dgts)(EvalTn(T[1], dgts) + EvalTn(T[2], dgts))
    elif T[0] == '*':
        return ComplexField(dgts)(EvalTn(T[1], dgts) * EvalTn(T[2], dgts))
    elif T[0] == '^':
        return ComplexField(dgts)(EvalTn(T[1], dgts) ^ EvalTn(T[2], dgts))
    else:
        print('IMPROPER INPUT !!!')

@cached_function
def nReducedNonMonotoneFormula(n, dgts):
    """
    Outputs non-monotone formula encodings of length at most n.

    EXAMPLES:

    ::

        sage: nReducedNonMonotoneFormula(3, 4)
        [[], [1, -1], [], []]

    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if n<=3:
        return [[], [1,-1], [], []]
    elif n>3:
        # Initialization of the list of formula.
        A=[[], [1,-1]] + [[] for t in range(n-1)]
        # Main loop.
        for sz in range(3,n+1):
            # Initialization of the fifth entry
            for i in range(1,sz-1):
                #A[sz]=A[sz]+[['+',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and not EvalTn(s).is_zero() and not EvalT(t).is_zero() and not EvalT(['+',s,t]).is_zero()]
                A[sz]=A[sz]+[['+',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and abs(EvalTn(s,dgts)) > ComplexField(dgts)(10^(-dgts),0) and abs(EvalTn(t,dgts)) > ComplexField(dgts)(10^(-dgts),0) and abs(EvalTn(['+',s,t],dgts)) > ComplexField(dgts)(10^(-dgts),0)]
            for i in range(1,sz-1):
                #A[sz]=A[sz]+[['*',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and EvalT(s)!=1 and EvalT(t)!=1 and EvalT(['*',s,t])!=1]
                A[sz]=A[sz]+[['*',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and abs(EvalTn(s,dgts) - ComplexField(dgts)(1,0)) > ComplexField(dgts)(10^(-dgts),0) and abs(EvalTn(t,dgts)-ComplexField(dgts)(1,0)) > ComplexField(dgts)(10^(-dgts),0) and abs(EvalTn(['*',s,t],dgts)-ComplexField(dgts)(1,0)) > ComplexField(dgts)(10^(-dgts),0)]
            for i in range(1,sz-1):
                #A[sz]=A[sz]+[['^',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and EvalT(s)!=1 and EvalT(t)!=1 and EvalT(['^',s,t])!=1]
                A[sz]=A[sz]+[['^',s,t] for s in A[i] for t in A[sz-i-1] if (len(A[i])>0) and (len(A[sz-i-1])>0) and abs(EvalTn(s,dgts)-ComplexField(dgts)(1,0)) > ComplexField(dgts)(10^(-dgts),0) and (EvalTn(t,dgts)-ComplexField(dgts)(1,0)) > ComplexField(dgts)(10^(-dgts),0) and abs(EvalTn(['^',s,t],dgts)-ComplexField(dgts)(1,0)) > ComplexField(dgts)(10^(-dgts),0)]
        return A

def nReducedNonMonotoneFormulaSets(sz, dgts):
    """
    Outputs set of numbers associated with non monotone encodings
    of complexity less then the size input parameter sz.

    EXAMPLES:

    ::

        sage: nReducedNonMonotoneFormulaSets(5,50)
        [{},
         {1.0000000000000, -1.0000000000000},
         {},
         {2.0000000000000, -2.0000000000000},
         {},
         {3.0000000000000, 1.0000000000000 - 2.4492935982947e-16*I, -3.0000000000000}]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    Lt=nReducedNonMonotoneFormula(sz,dgts)
    # Filling up the result list
    Rslt=[]
    for n in rg(sz+1):
        L=[]; i=0
        while i < len(Lt[n]):
            L.append(Lt[n][i])
            i=i+1
        Rslt.append(Set([EvalTn(L[i], dgts) for i in rg(len(L))]))
    # Cleaning up the list
    for i in rg(1,len(Rslt),2):
        for j in range(i):
            Rslt[i]=Rslt[i].difference(Rslt[j])
    return Rslt

def T2PreBool(expr):
    """
    Converts the Boolean formula written in the bracket tree encoding to the 
    Prefix string encoding notation


    EXAMPLES:

    The function implemented here tacitly assume that the input is valid

    ::

        sage: T2PreBool(['OR',var('x0'),var('x1')])
        '+x0x1'


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    s = str(expr)
    return ((((s.replace("[","")).replace("]","")).replace(",","")).replace("'","")).replace(" ","").replace('NOT','-').replace('OR','+').replace('AND','*')

def IncrementVariablesBool(T, incr):
    """
    Returns the same boolean formula with all variables are
    have their index incremented.


    EXAMPLES:

    The function implemented here tacitly assume that the input is valid

    ::

        sage: IncrementVariablesBool(['OR',var('x0'),var('x1')], 1)
        ['OR', x1, x2]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if type(T) == type(x):
        return var('x'+str(Integer(str(T)[1])+incr))
    elif T == 1:
        return 1
    elif T == 0:
        return 0
    elif T[0] == 'OR':
        return ['OR', IncrementVariablesBool(T[1], incr), IncrementVariablesBool(T[2], incr)]
    elif T[0] == 'AND':
        return ['AND', IncrementVariablesBool(T[1], incr), IncrementVariablesBool(T[2], incr)]
    elif T[0] == 'NOT':
        return ['NOT',IncrementVariablesBool(T[1], incr)]
    else:
        print('IMPROPER INPUT !!!')

def VariablesBool(T):
    """
    Returns the set of variables in the boolean formula


    EXAMPLES:

    The function implemented here tacitly assume that the input is valid

    ::

        sage: VariablesBool(['OR',var('x0'),var('x1')])
        {x1, x0}


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if type(T) == type(x):
        return Set([T])
    elif T == 0:
        return Set([])
    elif T == 1:
        return Set([])
    elif T[0] == 'OR':
        return VariablesBool(T[1]).union(VariablesBool(T[2]))
    elif T[0] == 'AND':
        return VariablesBool(T[1]).union(VariablesBool(T[2]))
    elif T[0] == 'NOT':
        return VariablesBool(T[1])
    else:
        print('IMPROPER INPUT !!!')

def CountVariablesBool(T):
    """
    Returns the number of variables in the boolean formula
    


    EXAMPLES:

    The function implemented here tacitly assume that the input is valid

    ::

        sage: CountVariablesBool(['OR',var('x0'),var('x1')])
        2


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    return VariablesBool(T).cardinality()

def EvalTBool(T):
    """
    Outputs the evaluation value of the output with binary assignement.

    EXAMPLES:

    ::

        sage: EvalTBool(['AND', ['OR', 1, ['OR', 1, 1]], ['AND', ['AND', 1, 1], 1]])
        1

    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger
    """
    if T == 1:
        return 1
    elif T == 0:
        return 0
    elif T[0] == 'OR':
        return EvalTBool(T[1]) or EvalTBool(T[2])
    elif T[0] == 'AND':
        return EvalTBool(T[1]) and EvalTBool(T[2])
    elif T[0] == 'NOT':
        return Integer(not EvalTBool(T[1]) )
    else:
        print('IMPROPER INPUT !!!')

def SubsTBool(T, Lx, Lv):
    """
    Outputs the evaluation value of the output with binary assignement.

    EXAMPLES:

    ::

        sage: SubsTBool(['AND', ['OR', 1, ['OR', 1, 1]], ['AND', ['AND', 1, 1], var('x0')]], [var('x0')], [1])
        1

    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger
    """
    if type(T) == type(x):
        return T.subs([Lx[i] == Lv[i] for i in rg(len(Lx))])
    elif T == 1:
        return 1
    elif T == 0:
        return 0
    elif T[0] == 'OR':
        if type(T[1]) == type(x) and type(T[2]) == type(x):
            return T[1].subs([Lx[u]==Lv[u] for u in rg(len(Lx))]) or T[2].subs([Lx[u]==Lv[u] for u in rg(len(Lx))])
        elif type(T[1]) == type(x) and type(T[2]) != type(x):
            return T[1].subs([Lx[u]==Lv[u] for u in rg(len(Lx))]) or SubsTBool(T[2], Lx, Lv)
        elif type(T[1]) != type(x) and type(T[2]) == type(x):
            return SubsTBool(T[1], Lx, Lv) or T[2].subs([Lx[u]==Lv[u] for u in rg(len(Lx))])
        else:
            return SubsTBool(T[1], Lx, Lv) or SubsTBool(T[2], Lx, Lv)
    elif T[0] == 'AND':
        if type(T[1]) == type(x) and type(T[2]) == type(x):
            return T[1].subs([Lx[u]==Lv[u] for u in rg(len(Lx))]) and T[2].subs([Lx[u]==Lv[u] for u in rg(len(Lx))])
        elif type(T[1]) == type(x) and type(T[2]) != type(x):
            return T[1].subs([Lx[u]==Lv[u] for u in rg(len(Lx))]) and SubsTBool(T[2], Lx, Lv)
        elif type(T[1]) != type(x) and type(T[2]) == type(x):
            return SubsTBool(T[1], Lx, Lv) and T[2].subs([Lx[u]==Lv[u] for u in rg(len(Lx))])
        else:
            return SubsTBool(T[1], Lx, Lv) and SubsTBool(T[2], Lx, Lv)
    elif T[0] == 'NOT':
        if type(T[1]) == type(x):
            return Integer(not T[1].subs([Lx[u]==Lv[u] for u in rg(len(Lx))]))
        else:
            return Integer(not SubsTBool(T[1], Lx, Lv))
    else:
        print('IMPROPER INPUT !!!')

def Bool2HM(T):
    """
    Outputs the hypermatrix encoding of the boolean formula
    the index correspond to the values assigned to the variables
    and the entry itself corresponds the 

    EXAMPLES:

    ::

        sage: Bool2HM(['AND', ['OR', var('x0'), ['OR', var('x0'), 1]], ['AND', ['AND', 1, var('x1')], 1]]).printHM()
        [:, :]=
        [0 1]
        [0 1]


    AUTHORS:

    - Edinah K. Gnang and Doron Zeilberger
    """
    # Initialization of the list specifying the dimensions of the output
    l = [2 for i in rg(CountVariablesBool(T))]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # Sorting the list of variables
        if VariablesBool(T).cardinality() == 1:
            LstVar=VariablesBool(T).list()
        else:
            LstVar=sum(VariablesBool(T)).operands()
        Rh[tuple(entry)]=SubsTBool(T, LstVar, entry)
    return Rh

def Bool2Integer(T):
    """
    Outputs the binary encoding of the boolean formula
    the integer is obtained by suming the output of 
    evaluation scaled by distinct powers of two

    EXAMPLES:

    ::

        sage: Bool2Integer(['AND', ['OR', var('x0'), ['OR', var('x0'), 1]], ['AND', ['AND', 1, var('x1')], 1]])
        16
        sage: X=var_list('x',2); Bool2Integer(x0)
        2
        sage: Bool2Integer(['NOT', x0])
        1
        sage: Bool2Integer(['OR', x0, ['NOT', x0]])
        3
        sage: Bool2Integer(['AND', x0, ['NOT', x0]])
        0


    AUTHORS:

    - Edinah K. Gnang and Doron Zeilberger
    """
    # Initialization of the list of evaliations
    L=Bool2HM(T).list()
    if CountVariablesBool(T) == 1:
        return sum(L[i]*2^i for i in rg(len(L)))
    else: 
        # Adding up the geomertric like sum 
        return sum(2^(2^k) for k in rg(1,CountVariablesBool(T)))+sum(L[i]*2^i for i in rg(len(L)))

def Bool2Poly(T):
    """
    Outputs the evaluation value of the output with binary assignement.

    EXAMPLES:

    ::

        sage: Bool2Poly(['AND', ['OR', var('x0'), ['OR', var('x1'), var('x2')]], ['AND', ['AND', var('x3'), var('x0')], ['NOT', var('x1')]]])
        -((x1*x2 - x1 - x2)*x0 - x1*x2 + x0 + x1 + x2)*x0*(x1 - 1)*x3

    AUTHORS:

    - Edinah K. Gnang and Doron Zeilberger
    """
    if T == 1:
        return SR(1)
    elif T == 0:
        return SR(0)
    elif type(T) == type(x):
        return T
    elif T[0] == 'OR':
        return Bool2Poly(T[1]) + Bool2Poly(T[2]) - Bool2Poly(T[1]) * Bool2Poly(T[2])
    elif T[0] == 'AND':
        return Bool2Poly(T[1]) * Bool2Poly(T[2])
    elif T[0] == 'NOT':
        return  1 - Bool2Poly(T[1])
    else:
        print('IMPROPER INPUT !!!')
 
def BoolTsize(T):
    """
    Outputs the size of the boolean formual associated with the formula.
 
    EXAMPLES:
    
    ::

        sage: BoolTsize(['AND', 0, 1])
        3


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    if T == 0 or T == 1 or type(T)==type(x):
        return 1
    elif T[0]=='NOT':
        return 1+BoolTsize(T[1])
    else:
        return 1+BoolTsize(T[1])+BoolTsize(T[2])

@cached_function
def ReducedNonMonotoneBooleanFormula(SZ):
    """
    Outputs the list of non-monotone boolean formula stratified by their
    size. The second output is the list of integer recording the tables
    which have occured thus far.

    EXAMPLES:

    ::

        sage: ReducedNonMonotoneBooleanFormula(3)[0]
        [[], [x0], [['NOT', x0]], [['AND', x0, x1], ['OR', x0, x1]]]
        sage: sz=5; X=var_list('x',sz); A=ReducedNonMonotoneBooleanFormula(sz)[0]
        sage: Lr = [[] for i in rg(len(A))]
        sage: for i in rg(len(A)):
        ....:     if len(A[i])>0:
        ....:         for T in A[i]:
        ....:             Lr[i].append((T, Bool2Integer(T)))
        ....:
        sage: Lr
        [[],
         [(x0, 2)],
         [(['NOT', x0], 1)],
         [(['AND', x0, x1], 12), (['OR', x0, x1], 18)],
         [(['AND', x0, ['NOT', x0]], 0),
          (['AND', x0, ['NOT', x1]], 6),
          (['AND', ['NOT', x0], x1], 8),
          (['OR', x0, ['NOT', x0]], 3),
          (['OR', x0, ['NOT', x1]], 15),
          (['OR', ['NOT', x0], x1], 17),
          (['NOT', ['AND', x0, x1]], 11),
          (['NOT', ['OR', x0, x1]], 5)],
         [(['AND', x0, ['AND', x1, x2]], 148),
          (['AND', x0, ['OR', x0, x1]], 14),
          (['AND', x0, ['OR', x1, x2]], 188),
          (['AND', ['OR', x0, x1], x1], 16),
          (['AND', ['OR', x0, x1], x2], 244),
          (['OR', x0, ['AND', x1, x2]], 254),
          (['OR', x0, ['OR', x1, x2]], 274),
          (['OR', ['AND', x0, x1], x2], 268)]]


    AUTHORS:

    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    # Initialization of the list of variables
    X=var_list('x',SZ)
    # Initialization of the list which store the integer encodings
    L = [Bool2Integer(X[0]), Bool2Integer(['NOT',var('x0')]), Bool2Integer(['AND', var('x0'), var('x1')]), Bool2Integer(['OR', var('x0'), var('x1')])]
    if SZ <= 3:
        return [[[], [X[0]], [['NOT',var('x0')]], [['AND', var('x0'), var('x1')], ['OR', var('x0'), var('x1')]]], L]
    elif SZ > 3:
        # Initialization of the list of formula.
        A=[[], [X[0]], [['NOT',var('x0')]], [['AND', var('x0'), var('x1')], ['OR', var('x0'), var('x1')]]] + [[] for t in range(SZ-3)]
        # Main loop.
        for sz in range(3,SZ+1):
            # Initialization of the fifth entry
            for i in rg(1,sz-1):
                for s in A[i]:
                    for t in A[sz-i-1]:
                        for j in rg(CountVariablesBool(s)+1):
                            if (len(A[i])>0) and (len(A[sz-i-1])>0) and not Bool2Integer(['AND', s, IncrementVariablesBool(t,j)]) in L:
                                A[sz].append(['AND', s, IncrementVariablesBool(t,j)])
                                L.append(Bool2Integer(['AND', s, IncrementVariablesBool(t,j)]))
            for i in range(1,sz-1):
                for s in A[i]:
                    for t in A[sz-i-1]:
                        for j in rg(CountVariablesBool(s)+1):
                            if (len(A[i])>0) and (len(A[sz-i-1])>0) and not Bool2Integer(['OR', s, IncrementVariablesBool(t,j)]) in L:
                                A[sz].append(['OR', s, IncrementVariablesBool(t,j)])
                                L.append(Bool2Integer(['OR', s, IncrementVariablesBool(t,j)]))
            for s in A[sz-1]:
                if (len(A[sz-1])>0) and not Bool2Integer(['NOT',s]) in L:
                    A[sz].append(['NOT', s])
                    L.append(Bool2Integer(['NOT', s]))
        return [A,L]

@cached_function
def ReducedNonMonotoneBooleanFormulaPoly(SZ):
    """
    Outputs the list of non-monotone boolean formula encoded as polynomial
    stratified by the size of the corresponding boolean formula. 
    The second output is the list of integer recording the tables
    which have occured thus far.

    EXAMPLES:

    ::

        sage: ReducedNonMonotoneBooleanFormulaPoly(3)[0]
        [[], [x0], [-x0 + 1], [x0*x1, -x0*x1 + x0 + x1]]
        sage: sz=5; X=var_list('x',sz); A=ReducedNonMonotoneBooleanFormula(sz)[0]
        sage: Lr = [[] for i in rg(len(A))]
        sage: for i in rg(len(A)):
        ....:     if len(A[i])>0:
        ....:         for T in A[i]:
        ....:             Lr[i].append((Bool2Poly(T), Bool2Integer(T)))
        ....:
        sage: Lr
        [[],
         [(x0, 2)],
         [(-x0 + 1, 1)],
         [(x0*x1, 12), (-x0*x1 + x0 + x1, 18)],
         [(-(x0 - 1)*x0, 0),
          (-x0*(x1 - 1), 6),
          (-(x0 - 1)*x1, 8),
          ((x0 - 1)*x0 + 1, 3),
          (x0*(x1 - 1) + x0 - x1 + 1, 15),
          ((x0 - 1)*x1 - x0 + x1 + 1, 17),
          (-x0*x1 + 1, 11),
          (x0*x1 - x0 - x1 + 1, 5)],
         [(x0*x1*x2, 148),
          (-(x0*x1 - x0 - x1)*x0, 14),
          (-(x1*x2 - x1 - x2)*x0, 188),
          (-(x0*x1 - x0 - x1)*x1, 16),
          (-(x0*x1 - x0 - x1)*x2, 244),
          (-x0*x1*x2 + x1*x2 + x0, 254),
          ((x1*x2 - x1 - x2)*x0 - x1*x2 + x0 + x1 + x2, 274),
          (-x0*x1*x2 + x0*x1 + x2, 268)]]


    AUTHORS:
    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    # Initialization of the list of variables
    X=var_list('x',SZ)
    # Initialization of the list which store the integer encodings
    L = [Bool2Integer(X[0]), Bool2Integer(['NOT',var('x0')]), Bool2Integer(['AND', var('x0'), var('x1')]), Bool2Integer(['OR', var('x0'), var('x1')])]
    if SZ <= 3:
        return [[[], [X[0]], [1-var('x0')], [x0*x1, x0 + x1- x0*x1]], L]
    elif SZ > 3:
        # Initialization of the list of formula.
        A=[[], [X[0]], [['NOT',var('x0')]], [['AND', var('x0'), var('x1')], ['OR', var('x0'), var('x1')]]] + [[] for t in range(SZ-3)]
        # Main loop.
        for sz in range(3,SZ+1):
            # Initialization of the fifth entry
            for i in rg(1,sz-1):
                for s in A[i]:
                    for t in A[sz-i-1]:
                        for j in rg(CountVariablesBool(s)+1):
                            if (len(A[i])>0) and (len(A[sz-i-1])>0) and not Bool2Integer(['AND', s, IncrementVariablesBool(t,j)]) in L:
                                A[sz].append(['AND', s, IncrementVariablesBool(t,j)])
                                L.append(Bool2Integer(['AND', s, IncrementVariablesBool(t,j)]))
            for i in range(1,sz-1):
                for s in A[i]:
                    for t in A[sz-i-1]:
                        for j in rg(CountVariablesBool(s)+1):
                            if (len(A[i])>0) and (len(A[sz-i-1])>0) and not Bool2Integer(['OR', s, IncrementVariablesBool(t,j)]) in L:
                                A[sz].append(['OR', s, IncrementVariablesBool(t,j)])
                                L.append(Bool2Integer(['OR', s, IncrementVariablesBool(t,j)]))
            for s in A[sz-1]:
                if (len(A[sz-1])>0) and not Bool2Integer(['NOT',s]) in L:
                    A[sz].append(['NOT', s])
                    L.append(Bool2Integer(['NOT', s]))
        # Returning Boolean formulas encoded as polynomials
        return [[[Bool2Poly(F) for F in A[t]] for t in rg(SZ)], L]

def ShortestNonMonotoneBooleanFormula(T):
    """
    Outputs the list of non-monotone boolean formula of minimal size
    having the same truth table as T

    EXAMPLES:

    ::

        sage: x0,x1=var('x0,x1'); ShortestNonMonotoneBooleanFormula(['OR',['OR', x0, x1],['AND',x0,x1]])
        [['OR', x0, x1]]



    AUTHORS:

    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    # Initialization of the stoping flag
    Ivl=Bool2Integer(T)
    # Initialization of the list of variables
    X=var_list('x',CountVariablesBool(T))
    # Initialization of the list which store the integer encodings
    L = [Bool2Integer(X[0]), Bool2Integer(['NOT',var('x0')]), Bool2Integer(['AND', var('x0'), var('x1')]), Bool2Integer(['OR', var('x0'), var('x1')])]
    if BoolTsize(T) <= 3:
        return T
    elif BoolTsize(T) > 3:
        # Initialization of the list of formula.
        A=[[], [X[0]], [['NOT',var('x0')]], [['AND', var('x0'), var('x1')], ['OR', var('x0'), var('x1')]]] + [[] for t in range(BoolTsize(T)-3)]
        # Main loop.
        for sz in range(3,BoolTsize(T)+1):
            # Initialization of the fifth entry
            for i in rg(1,sz-1):
                for s in A[i]:
                    for t in A[sz-i-1]:
                        for j in rg(CountVariablesBool(s)+1):
                            if (len(A[i])>0) and (len(A[sz-i-1])>0) and not Bool2Integer(['AND', s, IncrementVariablesBool(t,j)]) in L:
                                A[sz].append(['AND', s, IncrementVariablesBool(t,j)])
                                L.append(Bool2Integer(['AND', s, IncrementVariablesBool(t,j)]))
            for i in range(1,sz-1):
                for s in A[i]:
                    for t in A[sz-i-1]:
                        for j in rg(CountVariablesBool(s)+1):
                            if (len(A[i])>0) and (len(A[sz-i-1])>0) and not Bool2Integer(['OR', s, IncrementVariablesBool(t,j)]) in L:
                                A[sz].append(['OR', s, IncrementVariablesBool(t,j)])
                                L.append(Bool2Integer(['OR', s, IncrementVariablesBool(t,j)]))
            for s in A[sz-1]:
                if (len(A[sz-1])>0) and not Bool2Integer(['NOT',s]) in L:
                    A[sz].append(['NOT', s])
                    L.append(Bool2Integer(['NOT', s]))
            if Ivl in L:
                break
        return [F for F in A[sz] if Ivl==Bool2Integer(F)]

def Poly2HM(P,X):
    """
    Outputs the hypermatrix encoding of the boolean formula
    the index correspond to the values assigned to the variables
    and the entry itself corresponds the output of the boolean
    function at those particular inputs. 

    EXAMPLES:

    ::

        sage: x0,x1=var('x0,x1'); Poly2HM((x0 - 1)*(x1 - 1) + x0*x1, [x0,x1]).printHM()
        [:, :]=
        [1 0]
        [0 1]


    AUTHORS:

    - Edinah K. Gnang and Doron Zeilberger
    """
    # Initialization of the list specifying the dimensions of the output
    l = [2 for i in rg(len(X))]
    # Initializing the input for generating a symbolic hypermatrix
    inpts = l+['zero']
    # Initialization of the hypermatrix
    Rh = HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=P.subs([X[v]==entry[v] for v in rg(len(l))])
    return Rh

def Poly2Integer(P,X):
    """
    Outputs the binary encoding of the boolean formula
    the integer is obtained by suming the output of 
    evaluation scaled by distinct powers of two.
    The inputs are polynomial P and a list of variable X
    accounting for all the variables in P


    EXAMPLES:

    ::

        sage: x0,x1=var('x0,x1'); Poly2Integer((x0 - 1)*(x1 - 1) + x0*x1, [x0,x1])
        13


    AUTHORS:

    - Edinah K. Gnang and Doron Zeilberger
    """
    # Initialization of the list of evaliations
    H=Poly2HM(P,X); L=H.list()
    if H.order() == 1:
        return sum(L[i]*2^i for i in rg(len(L)))
    else: 
        # Adding up the geomertric like sum 
        return sum(2^(2^k) for k in rg(1,H.order()))+sum(L[i]*2^i for i in rg(len(L)))

def ShortestNonMonotoneBooleanFormulaII(P,X):
    """
    Outputs the list of non-monotone boolean formula of minimal size.
    having the same truth table as the input polynomial P.
    The difference with the implementation above is that it takes as
    input a polynomial


    EXAMPLES:

    ::

        sage: x0,x1=var('x0,x1'); ShortestNonMonotoneBooleanFormulaII((x0 - 1)*(x1 - 1) + x0*x1, [x0,x1])
        [['OR', ['AND', x0, x1], ['NOT', ['OR', x0, x1]]]]


    AUTHORS:

    - Edinah K. Gnang and Doron Zeilberger

    To Do :
    - Try to implement faster version of this procedure

    """
    # Initialization of the stoping flag
    Ivl=Poly2Integer(P,X)
    # Initialization of the size of the formula from the size of the infix polynomial interpolation
    SZ=len(str(P))
    # Initialization of the list which store the integer encodings
    L = [Bool2Integer(X[0]), Bool2Integer(['NOT',var('x0')]), Bool2Integer(['AND', var('x0'), var('x1')]), Bool2Integer(['OR', var('x0'), var('x1')])]
    if SZ <= 3:
        return P
    elif SZ > 3:
        # Initialization of the list of formula.
        A=[[], [X[0]], [['NOT',var('x0')]], [['AND', var('x0'), var('x1')], ['OR', var('x0'), var('x1')]]] + [[] for t in range(SZ-3)]
        # Main loop.
        for sz in range(3,SZ+1):
            # Initialization of the fifth entry
            for i in rg(1,sz-1):
                for s in A[i]:
                    for t in A[sz-i-1]:
                        for j in rg(CountVariablesBool(s)+1):
                            if (len(A[i])>0) and (len(A[sz-i-1])>0) and not Bool2Integer(['AND', s, IncrementVariablesBool(t,j)]) in L:
                                A[sz].append(['AND', s, IncrementVariablesBool(t,j)])
                                L.append(Bool2Integer(['AND', s, IncrementVariablesBool(t,j)]))
            for i in range(1,sz-1):
                for s in A[i]:
                    for t in A[sz-i-1]:
                        for j in rg(CountVariablesBool(s)+1):
                            if (len(A[i])>0) and (len(A[sz-i-1])>0) and not Bool2Integer(['OR', s, IncrementVariablesBool(t,j)]) in L:
                                A[sz].append(['OR', s, IncrementVariablesBool(t,j)])
                                L.append(Bool2Integer(['OR', s, IncrementVariablesBool(t,j)]))
            for s in A[sz-1]:
                if (len(A[sz-1])>0) and not Bool2Integer(['NOT',s]) in L:
                    A[sz].append(['NOT', s])
                    L.append(Bool2Integer(['NOT', s]))
            if Ivl in L:
                break
        return [F for F in A[sz] if Ivl==Bool2Integer(F)]

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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # computing the Hypermatrix product
        if len(args) < 2:
            raise ValueError("The number of operands must be >= 2")
        elif len(args) >= 2:
            Rh[tuple(entry)] = 0
            l2 = [B.n(sz) for sz in range(B.order())]
            for j in range(prod(l2)):
                # Turning the index j into an hypermatrix array location using the decimal encoding trick
                entry2 = [Integer(mod(j,l2[0]))]
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

        sage: PathAdjcencyHypermatrix(HM(2,2,'a'), 3)
        [[[a00^2, a00*a01], [a01*a10, a01*a11]], [[a00*a10, a01*a10], [a10*a11, a11^2]]]

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
            entry = [Integer(mod(i,l[0]))]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            Rh[tuple(entry)] = prod([A[tuple([entry[i],entry[i+1]])] for i in range(pthl-1)])
    else:
        raise ValueError("Input hypermatrix must be order 2 and the path length must be an integer greater then 0")
    return Rh

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
            #A0 = apply(GeneralHypermatrixProductB, [A for k in range(od)]+[A0])
            A0 = GeneralHypermatrixProductB( *([A for k in range(od)]+[A0]) )
            #A1 = apply(GeneralHypermatrixProductB, [A for k in range(od)]+[A1])
            A1 = GeneralHypermatrixProductB( *([A for k in range(od)]+[A1]) )
            # Append the result to the list
            L.append(A0.list()); L.append(A1.list())
        return Matrix(SR,L)
    else:
        # return the error message if the input hypermatrix is cubic
        raise ValueError("The input hypermpatrix must be cubic")

def fast_reduce(f, monom, subst):
    """
    computes the reduction by monomial substitution
    by converting the symbolic expression into a string
    of characters and performing the substition on the
    string and converts back in the end the obtained
    string back into a symbolic expression SR
    
    EXAMPLES:
 
    ::

        sage: x1,x2,x3=var('x1, x2, x3'); fast_reduce(x1^3+x2+x3^3,[x3^3],[1])
        x1^3 + x2 + 1

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    if len(monom) == len(subst):
        s = str(f)
        for i in range(len(monom)):
            s = s.replace(str(monom[i]), '('+str(subst[i])+')')
        return expand((SR(s)).simplify_full())
    else:
        print('Error the monomial list and the substitution list must have the same length')

def fast_reduce_no_expand(f, monom, subst):
    """
    computes the reduction by the symbolic expression substitution
    does not expand the resulting polynomial
    
    EXAMPLES:
 
    ::

        sage: x1,x2,x3=var('x1, x2, x3'); fast_reduce_no_expand(x1^3+x2+x3^3,[x3^3+x2],[(x2+x3)^2])
        x1^3 + (x2 + x3)^2

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    if len(monom) == len(subst):
        s = str(f)
        for i in range(len(monom)):
            s = s.replace(str(monom[i]), '('+str(subst[i])+')')
        return SR(s)
    else:
        print('Error the monomial list and the substitution list must have the same length')

def fast_reduceII(f, monom, subst, Dct):
    """
    computes the reduction by monomial substitution
    by converting the symbolic expression into a string
    of characters and performing the substition on the
    string and converts back in the end the obtained
    string back into an expression.
    The difference with the implementation above is the
    fact that this implementation handles matrix input
    of type HM

    
    EXAMPLES:
 
    ::

        sage: x1,x2,x3,A,B,C=var('x1, x2, x3, A, B, C')
        sage: fast_reduceII(x1^3+x2+x3^3,[x1,x2,x3],[A,B,C],{'A':HM(2,2,'a'),'B':HM(2,2,'b'),'C':HM(2,2,'c')}).printHM()
        [:, :]=
        [(a00^2 + a01*a10)*a00 + (a00*a10 + a10*a11)*a01 + (c00^2 + c01*c10)*c00 + (c00*c10 + c10*c11)*c01 + b00 (a00*a01 + a01*a11)*a00 + (a01*a10 + a11^2)*a01 + (c00*c01 + c01*c11)*c00 + (c01*c10 + c11^2)*c01 + b01]
        [(a00^2 + a01*a10)*a10 + (a00*a10 + a10*a11)*a11 + (c00^2 + c01*c10)*c10 + (c00*c10 + c10*c11)*c11 + b10 (a00*a01 + a01*a11)*a10 + (a01*a10 + a11^2)*a11 + (c00*c01 + c01*c11)*c10 + (c01*c10 + c11^2)*c11 + b11]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    if len(monom) == len(subst):
        s = str(f)
        for i in range(len(monom)):
            s = s.replace(str(monom[i]), str(subst[i]))
        return sage_eval(s, locals=Dct)
    else:
        print('Error the monomial list and the substitution list must have the same length')

def SecondOrderHyperdeterminant(H):
    """
    computes second order i.e. matrix determinant using the Cauchy umbral mnemonic device
    The algorithm used here is worst then sum over signed permutation monomials because it
    constructs 2^binomial(n,2) monomials in the expanded form. There is however a massive
    amount of cancelation which will leaves out only factorial(n) terms.
    
    EXAMPLES:
 
    ::

        sage: SecondOrderHyperdeterminant(HM(2,2,'a'))
        -a01*a10 + a00*a11
        

    AUTHORS:

    - Edinah K. Gnang
    - To Do: 
    """
    # Testing to see that the hypermatrix is indeed a cube
    if len(Set(H.dimensions()).list())==1 and H.order()==2:
        # Initializing the matrix for the mnemonic construction
        A=HM(H.n(0),1,[var('x'+str(i)) for i in range(1,1+H.n(0))])
        L=expand(prod([A[u,0]-A[v,0] for u in range(1,A.nrows()) for v in range(u)])*prod([A[i,0] for i in range(A.nrows())])).operands()
        # Computing the polynomial
        f=sum(L)
        # Loop performing the umbral expression
        for k in range(H.n(0),0,-1):
            f=fast_reduce(f,[var('x'+str(i))^k for i in range(1,1+H.n(0))],[var('a'+str(i)+str(k)) for i in range(1,1+H.n(0))])
        return f.subs(dict([(var('a'+str(i)+str(k)),H[i-1,k-1]) for i in range(1,1+H.n(0)) for k in range(1,1+H.n(0))]))
    else :
        # Print an error message indicating that the matrix must be a cube.
        raise ValueError("The matrix must be square.")

def ThirdOrderHyperdeterminant(H):
    """
    computes third order hypermatrix determinant using the recursive determinant construction.
    Where the base case is the matrix case. The third order hyperdeterminant are only defined
    for cubic third order hypermatrices.
    
    EXAMPLES:
 
    ::

        sage: ThirdOrderHyperdeterminant(HM(2,2,2,'a'))
        a000*a011*a101*a110 - a001*a010*a100*a111
        

    AUTHORS:

    - Edinah K. Gnang
    - To Do: Implement a faster and more generic version.
    """
    # Testing to see that the hypermatrix is indeed a cube
    if len(Set(H.dimensions()).list())==1 and H.order()==3:
        # Initializing the matrix for the mnemonic construction
        A=HM(H.n(0),H.n(0),[var('x'+str(i)+str(j)) for j in range(1,1+H.n(0)) for i in range(1,1+H.n(0))])
        # Computing the mnemonique polynomial
        L=expand(SecondOrderHyperdeterminant(A)*prod([sum([sqrt(g^2).canonicalize_radical() for g in SecondOrderHyperdeterminant(A.elementwise_exponent(j)).operands()]) for j in range(2,1+A.n(0))])).operands()
        # Computing the polynomial
        f=sum([l for l in L if len((l^2).operands())==(H.n(0))^2])
        # Loop performing the umbral expression
        for k in range(H.n(0),0,-1):
            f=fast_reduce(f,[var('x'+str(i)+str(j))^k for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0))],[var('a'+str(i)+str(j)+str(k)) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0))])
        return f.subs(dict([(var('a'+str(i)+str(j)+str(k)),H[i-1,j-1,k-1]) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0))]))
    else :
        # Print an error message indicating that the matrix must be a cube.
        raise ValueError("The hypermatrix must be a third order cube hypermatrix.")

def FourthOrderHyperdeterminant(H):
    """
    computes fourth order hyperdeterminant using the recusive construction.
    The base case of the recusrion is the matrix determinant.
    
    EXAMPLES:
 
    ::

        sage: FourthOrderHyperdeterminant(HM(2,2,2,2,'a'))
        -a0001*a0010*a0100*a0111*a1000*a1011*a1101*a1110 + a0000*a0011*a0101*a0110*a1001*a1010*a1100*a1111

    AUTHORS:

    - Edinah K. Gnang
    - To Do: 
    """
    # Testing to see that the hypermatrix is indeed a cube
    if len(Set(H.dimensions()).list())==1 and H.order()==4:
        # Initializing the matrix for the mnemonic construction
        A=HM(H.n(0),H.n(0),H.n(0), [var('x'+str(i)+str(j)+str(k)) for k in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for i in range(1,1+H.n(0))])
        # Computing the polynomial
        L=expand(ThirdOrderHyperdeterminant(A)*prod([sum([sqrt(g^2).canonicalize_radical() for g in ThirdOrderHyperdeterminant(A.elementwise_exponent(j)).operands()]) for j in range(2,1+A.n(0))])).operands()
        # Computing the polynomial
        f = sum([l for l in L if len((l^2).operands())==(H.n(0))^3])
        # Loop performing the umbral expression
        for l in range(H.n(0),0,-1):
            f = fast_reduce(f,[var('x'+str(i)+str(j)+str(k))^l for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0))],[var('a'+str(i)+str(j)+str(k)+str(l)) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0))])
        return f.subs(dict([(var('a'+str(i)+str(j)+str(k)+str(l)),H[i-1,j-1,k-1,l-1]) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for l in range(1,1+H.n(0))]))
    else :
        # Print an error message indicating that the matrix must be a cube.
        raise ValueError("The hypermatrix must be a fourth order hypercube hypermatrix.")

def FifthOrderHyperdeterminant(H):
    """
    computes the fifth order hyperderterminant using the recusrive construction using the matrix determinant as
    the base case.
    
    EXAMPLES:
 
    ::

        sage: FifthOrderHyperdeterminant(HM(2,2,2,2,2,'a'))
        a00000*a00011*a00101*a00110*a01001*a01010*a01100*a01111*a10001*a10010*a10100*a10111*a11000*a11011*a11101*a11110 - a00001*a00010*a00100*a00111*a01000*a01011*a01101*a01110*a10000*a10011*a10101*a10110*a11001*a11010*a11100*a11111

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Testing to see that the hypermatrix is indeed a cube
    if len(Set(H.dimensions()).list())==1 and H.order()==5:
        # Initializing the matrix for the mnemonic construction
        A=HM(H.n(0),H.n(0),H.n(0),H.n(0),[var('x'+str(i)+str(j)+str(k)+str(l)) for l in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for i in range(1,1+H.n(0))])
        # Computing the polynomial
        L = expand(FourthOrderHyperdeterminant(A)*prod([sum([sqrt(g^2).canonicalize_radical() for g in FourthOrderHyperdeterminant(A.elementwise_exponent(j)).operands()]) for j in range(2,1+A.n(0))])).operands()
        f = sum([l for l in L if len((l^2).operands())==(H.n(0))^4])
        # Loop performing the umbral expression
        for m in range(H.n(0),0,-1):
            f = fast_reduce(f,[var('x'+str(i)+str(j)+str(k)+str(l))^m for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for l in range(1,1+H.n(0))],[var('a'+str(i)+str(j)+str(k)+str(l)+str(m)) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for l in range(1,1+H.n(0))])
        return f.subs(dict([(var('a'+str(i)+str(j)+str(k)+str(l)+str(m)),H[i-1,j-1,k-1,l-1,m-1]) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for l in range(1,1+H.n(0)) for m in range(1,1+H.n(0))]))
    else :
        # Print an error message indicating that the matrix must be a cube.
        raise ValueError("The hypermatrix must be a fifth order hypercube hypermatrix.")

def SixthOrderHyperdeterminant(H):
    """
    computes the sixth order hypermatrix determinant using the recursive construction.
    The base case of the recursion being the matrix determinant formula
    
    EXAMPLES:
 
    ::

        sage: SixthOrderHyperdeterminant(HM(2,2,2,2,2,2,'a'))
        -a000001*a000010*a000100*a000111*a001000*a001011*a001101*a001110*a010000*a010011*a010101*a010110*a011001*a011010*a011100*a011111*a100000*a100011*a100101*a100110*a101001*a101010*a101100*a101111*a110001*a110010*a110100*a110111*a111000*a111011*a111101*a111110 + a000000*a000011*a000101*a000110*a001001*a001010*a001100*a001111*a010001*a010010*a010100*a010111*a011000*a011011*a011101*a011110*a100001*a100010*a100100*a100111*a101000*a101011*a101101*a101110*a110000*a110011*a110101*a110110*a111001*a111010*a111100*a111111
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Testing to see that the hypermatrix is indeed a cube
    if len(Set(H.dimensions()).list())==1 and H.order()==6:
        # Initializing the matrix for the mnemonic construction
        A=HM(H.n(0),H.n(0),H.n(0),H.n(0),H.n(0),[var('x'+str(i)+str(j)+str(k)+str(l)+str(m)) for m in range(1,1+H.n(0)) for l in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for i in range(1,1+H.n(0))])
        # Computing the polynomial
        L = expand(FifthOrderHyperdeterminant(A)*prod([sum([sqrt(g^2).canonicalize_radical() for g in FifthOrderHyperdeterminant(A.elementwise_exponent(j)).operands()]) for j in range(2,1+A.n(0))])).operands()
        f = sum([l for l in L if len((l^2).operands())==(H.n(0))^5])
        # Loop performing the umbral expression
        for p in range(H.n(0),0,-1):
            f = fast_reduce(f,[var('x'+str(i)+str(j)+str(k)+str(l)+str(m))^p for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for l in range(1,1+H.n(0)) for m in range(1,1+H.n(0))],[var('a'+str(i)+str(j)+str(k)+str(l)+str(m)+str(p)) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for l in range(1,1+H.n(0)) for m in range(1,1+H.n(0))])
        return f.subs(dict([(var('a'+str(i)+str(j)+str(k)+str(l)+str(m)+str(p)),H[i-1,j-1,k-1,l-1,m-1,p-1]) for i in range(1,1+H.n(0)) for j in range(1,1+H.n(0)) for k in range(1,1+H.n(0)) for l in range(1,1+H.n(0)) for m in range(1,1+H.n(0)) for p in range(1,1+H.n(0))]))
    else :
        # Print an error message indicating that the matrix must be a cube.
        raise ValueError("The hypermatrix must be a sixth order hypercube hypermatrix.")

def general_side_length_2_det(A):
    """
    outputs the symbolic expression with the determinant of hypermatrices of arbitrary orders.
    but every size of the hypermatrix must be equal to two. It ouputs an equality derived via
    the rank one argument. The difference with the function above is that this function
    takes a hypermatrix as input.

    EXAMPLES:
 
    ::

        sage: general_side_length_2_det(HM(2,2,'a'))
        -a01*a10 + a00*a11

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    if A.is_cubical() and A.n(0)==2:
        # Initialization of the list specifying the dimensions of the output
        l = [A.n(i) for i in range(A.order())]
        # Initializing the input for generating a symbolic hypermatrix
        inpts = l+['zero']
        # Initialization of the list of odd and even index list
        Lodd=[]; Leven=[]
        Rh = HM(*inpts)
        # Main loop performing the transposition of the entries
        for i in range(prod(l)):
            # Turning the index i into an hypermatrix array location using the decimal encoding trick
            entry = [Integer(mod(i,l[0]))]
            sm = Integer(mod(i,l[0]))
            for k in range(len(l)-1):
                entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
            if Integer(mod(sum(entry),2)) == 0:
                Leven.append(A[tuple(entry)])
            else:
                Lodd.append(A[tuple(entry)])
        return prod(Leven)-prod(Lodd)
    else:
        raise ValueError("The input hypermatrix must cubical with side lentgh 2.")

def DodgsonCondensation(A):
    """
    outputs the symbolic expression deduced from the Dodgson condensation applied to matrices
    and hypermatrices. The side length of the input hypermatrix must be greater then 2.

    EXAMPLES:
 
    ::

        sage: DodgsonCondensation(HM(3,3,'a')) 
        -((a02*a11 - a01*a12)*(a11*a20 - a10*a21) - (a01*a10 - a00*a11)*(a12*a21 - a11*a22))/a11

    AUTHORS:

    - Edinah K. Gnang
    - To Do: 
    """
    if A.is_cubical() and A.n(0)>2:
        # Initializing the copy
        Bold=A.copy()
        #Bnew=apply(HM,[i-1 for i in Bold.dimensions()]+['zero'])
        Bnew=HM( *([i-1 for i in Bold.dimensions()]+['zero']) )
        Temp=HM( *([i-2 for i in Bold.dimensions()]+['zero']) )
        while Bnew.n(0)>1:
            # Filling up Bnew
            l=Bnew.dimensions()
            # Main loop performing the transposition of the entries
            for i in range(prod(l)):
                # Turning the index i into an hypermatrix array location using the decimal encoding trick
                entry=[Integer(mod(i,l[0]))]
                sm=Integer(mod(i,l[0]))
                for k in range(len(l)-1):
                    entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                    sm=sm+prod(l[0:k+1])*entry[len(entry)-1]
                # Initialization of the summand
                Bl=[]
                l2=[2 for j in range(A.order())]
                for j in range(prod(l2)):
                    ent=[Integer(mod(j,l2[0]))]
                    ms=Integer(mod(j,l2[0]))
                    for t in range(len(l2)-1):
                        ent.append(Integer(mod(Integer((j-ms)/prod(l2[0:t+1])),l2[t+1])))
                        ms=ms+prod(l2[0:t+1])*ent[len(ent)-1]
                    Bl.append((Matrix(ZZ,entry)+Matrix(ZZ,ent)).list())
                Bnew[tuple(entry)]=general_side_length_2_det(HM(*([2 for j in range(A.order())]+[[Bold[tuple(entry2)] for entry2 in Bl ]])))
            # Filling up Temp
            Temp=HM(*([i-2 for i in Bold.dimensions()]+['zero']))
            l=Temp.dimensions()
            # Main loop performing the transposition of the entries
            for i in range(prod(l)):
                # Turning the index i into an hypermatrix array location using the decimal encoding trick
                entry = [Integer(mod(i,l[0]))]
                sm = Integer(mod(i,l[0]))
                for k in range(len(l)-1):
                    entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                    sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                # Initialization of the summand
                Bl=[]
                l2=[2 for j in range(A.order())]
                for j in range(prod(l2)):
                    ent=[Integer(mod(j,l2[0]))]
                    ms=Integer(mod(j,l2[0]))
                    for t in range(len(l2)-1):
                        ent.append(Integer(mod(Integer((j-ms)/prod(l2[0:t+1])),l2[t+1])))
                        ms = ms+prod(l2[0:t+1])*ent[len(ent)-1]
                    Bl.append((Matrix(ZZ,entry)+Matrix(ZZ,ent)).list())
                #Temp[tuple(entry)]=general_side_length_2_det(apply(HM,[2 for j in range(A.order())]+[[Bnew[tuple(entry2)] for entry2 in Bl ]]))/Bold[tuple((Matrix(ZZ,entry)+ones_matrix(1,len(entry))).list())]
                Temp[tuple(entry)]=general_side_length_2_det(HM(*([2 for j in range(A.order())]+[[Bnew[tuple(entry2)] for entry2 in Bl ]])))/Bold[tuple((Matrix(ZZ,entry)+ones_matrix(1,len(entry))).list())]
            # Performing the update
            if Temp.n(0)>0:
                Bold=Bnew.copy()
                Bnew=Temp.copy()
        return (Temp.list())[0]
    else:
        raise ValueError("The input hypermatrix must be hypercubic of size 2")

def GeneralHyperdeterminantExpression(od,sz):
    """
    computes the general hyperdeterminant expression using the recursive determinant
    construction. Where the base case is the matrix determinant.
    
    EXAMPLES:
 
    ::

        sage: GeneralHyperdeterminantExpression(3,2)
        x111*x122*x212*x221 - x112*x121*x211*x222
        

    AUTHORS:

    - Edinah K. Gnang
    - To Do: Implement a faster and more generic version.
    """
    # Testing to see that the hypermatrix is indeed a cube
    if od==2:
        return HM(sz,sz,'x','shift').det()
    elif od>2:
        # Initializing the base case of the recursion
        A=HM(sz,sz,'x','shift').copy()
        # Computing the matrix determinant and permanent pair
        Ldtm=[Deter(A)]+[Per(A.elementwise_exponent(j)) for j in range(2,1+sz)]
        # Main loop performing the recursion computation
        for o in range(2,od):
            # ReInitialization of the Hypermatrix A
            #A=apply(HM,[sz for i in range(o)]+['x']+['shift']).copy()
            A=HM(*([sz for i in range(o)]+['x']+['shift'])).copy()
            # Computing the mnemonique polynomial
            L=expand(prod(Ldtm)).operands()
            # Computing the polynomial
            f=sum([l for l in L if len((l^2).operands())==prod(A.dimensions())])
            # Loop performing the umbral transformation
            for k in range(sz,0,-1):
                f=fast_reduce(f,A.elementwise_exponent(k).list(), HM(*([sz for i in range(o)]+['a']+['shift'])).append_index(k).list())
            B=HM(*([sz for i in range(o+1)]+['x']+['shift'])).copy()
            L1=HM(*([sz for i in range(o+1)]+['x']+['shift'])).list();L2=HM(*([sz for i in range(o+1)]+['a']+['shift'])).list()
            f=f.subs(dict([(L2[i],L1[i]) for i in range(len(L1))]))
            Ldtm=[f]+[sum([sqrt(g^2).canonicalize_radical() for g in f.operands()]).subs(dict([(B.list()[i],B.elementwise_exponent(j).list()[i]) for i in range(prod(B.dimensions()))])) for j in range(2,1+sz)]
        return f

def GeneralHyperdeterminant(H):
    """
    computes the general hyperdeterminant using the recursive determinant construction.
    Where the base case is the matrix determinant.
    
    EXAMPLES:
 
    ::

        sage: GeneralHyperdeterminant(HM(2,2,2,'a','shift'))
        a111*a122*a212*a221 - a112*a121*a211*a222
        

    AUTHORS:

    - Edinah K. Gnang
    - To Do: Implement a faster and more generic version.
    """
    if H.is_cubical(): 
        f=GeneralHyperdeterminantExpression(H.order(),H.n(0))
        Lh=H.list()
        #Lx=apply(HM,[H.n(0) for i in range(H.order())]+['x']+['shift']).list()
        Lx=HM(*([H.n(0) for i in range(H.order())]+['x']+['shift'])).list()
        return f.subs(dict([(Lx[i],Lh[i]) for i in range(len(Lh))]))
    else:
        raise ValueError("The hypermatrix must be cubical.")

def LatinHypermatrixList(od,sz):
    """
    computes list of latin hypercubes of order od and side length sz.
    
    EXAMPLES:
 
    ::

        sage: LatinHypermatrixList(2,2)
        [[[0, x01], [x10, 0]], [[x00, 0], [0, x11]]]
        

    AUTHORS:

    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the hyperdeterminant expression.
    #Dt=GeneralHyperdeterminant(apply(HM,[sz for i in range(od)]+['x']))
    Dt=GeneralHyperdeterminant(HM(*([sz for i in range(od)]+['x'])))
    Lt = [sqrt(tm^2).canonicalize_radical() for tm in Dt.operands()]
    Sa = Set(HM(*([sz for i in range(od)]+['x'])).list())
    L = []
    for f in Lt:
        Tmp=HM(*([sz for i in range(od)]+['x']))
        Lsf=Sa.difference(Set(f.operands()))
        L.append(Tmp.subs(dict([(g,0) for g in Lsf])))
    return L

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
    L=[]
    for j in rg(1,od)+[0]:
        TmpL=[sz for i in rg(j)]+[1]+[sz for i in rg(j+1,od)]+[AlphaB[j]]
        L.append(HM(*TmpL))
    # Initilizing the list of variable
    VrbLst=[]
    for Q in L:
        VrbLst = VrbLst+Q.list()
    Eq=Prod(*[Q for Q in L])
    CnstrLst=[eq==1 for eq in Eq.list()]
    [A,b]=multiplicativeConstraintFormator(CnstrLst, VrbLst)
    return A

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
    bns = Integer(l).str(2)
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
    bns = Integer(l).str(2)
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
    bns = Integer(l).str(2)
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
    bns = Integer(l).str(2)
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

def FifthOrderHadamardBlockU(l):
    """
    outputs the  direct sum construction of hadamard block hypermatrices.
    the vectors of the matrices are not normalized.

    EXAMPLES:
 
    ::

        sage: FifthOrderHadamardBlockU(2)
        [[[[[1, 1], [1, 1]], [[1, 1], [1, 1]]], [[[1, 1], [1, 1]], [[1, 1], [1, 1]]]], [[[[1, 1], [-1, -1]], [[1, 1], [1, 1]]], [[[-1, 1], [1, 1]], [[1, 1], [1, 1]]]]]

    AUTHORS:

    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing the 2x2x2x2x2 hadamard hypermatrix.
    Hd=HM([[[[[1,1],[1,1]],[[1,1],[1,1]]],[[[1,1],[1,1]],[[1,1],[1,1]]]],[[[[1,1],[-1,-1]],[[1,1],[1,1]]],[[[-1,1],[1,1]],[[1,1],[1,1]]]]])
    # Obtaining the binary encoding of the input integer
    bns = Integer(l).str(2)
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
        [[[[[(1/2)^(1/5), (1/2)^(1/5)], [(1/2)^(1/5), (1/2)^(1/5)]], [[(1/2)^(1/5), (1/2)^(1/5)], [(1/2)^(1/5), (1/2)^(1/5)]]], [[[(1/2)^(1/5), (1/2)^(1/5)], [(1/2)^(1/5), (1/2)^(1/5)]], [[(1/2)^(1/5), (1/2)^(1/5)], [(1/2)^(1/5), (1/2)^(1/5)]]]], [[[[(1/2)^(1/5), (1/2)^(1/5)], [-(1/2)^(1/5), -(1/2)^(1/5)]], [[(1/2)^(1/5), (1/2)^(1/5)], [(1/2)^(1/5), (1/2)^(1/5)]]], [[[-(1/2)^(1/5), (1/2)^(1/5)], [(1/2)^(1/5), (1/2)^(1/5)]], [[(1/2)^(1/5), (1/2)^(1/5)], [(1/2)^(1/5), (1/2)^(1/5)]]]]]

    AUTHORS:

    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing the 2x2x2x2x2 hadamard Hypermatrix.
    Hd = (1/2)^(1/5)*HM([[[[[1,1],[1,1]],[[1,1],[1,1]]],[[[1,1],[1,1]],[[1,1],[1,1]]]],[[[[1,1],[-1,-1]],[[1,1],[1,1]]],[[[-1,1],[1,1]],[[1,1],[1,1]]]]])
    # Obtaining the binary encoding of the input integer
    bns = Integer(l).str(2)
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
    Outputs  the random transposition permutation of n elements.
    The fucntion is used for randominzing the Hadamard block construction
    

    EXAMPLES:
 
    ::

        sage: min(RandomTransposition(3))
        0 

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

def SecondOrderHypermatrixResolutionPartition(U, V, Ha, Hb, NbPrts=2):
    """
    outputs the spliting of a second order hypermatrix into the pieces
    as suggested by the resolution of identity. The first three input
    hypermatrices are uncorrelated tuples.
    the last three inputs correspond to the factors for the spliting.

    EXAMPLES:
 
    ::

        sage: [U,V]=GeneralUncorrelatedHypermatrixTuple(2)
        sage: L=SecondOrderHypermatrixResolutionPartition(U, V, HM(2,2,'a'), HM(2,2,'b'))
        sage: len(L)
        2
        sage: sum(L).simplify_full()
        [[a00*b00 + a01*b10, a00*b01 + a01*b11], [a10*b00 + a11*b10, a10*b01 + a11*b11]]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the list storing the projectors
    L = []
    for i in range(U.n(1)):
        # Filling up the slice
        Tp0 = U.slice([i],'col')
        # Filling up the slice
        Tp1 = V.slice([i],'row')
        # Appending the components to the list
        L.append(Prod(Tp0, Tp1))
    if len(L)>NbPrts:
        return [ProdB(Ha,Hb,sum(L[j*NbPrts:min((j+1)*NbPrts,len(L))])) for j in range(ceil(len(L)/NbPrts))]
    else:
        return [ProdB(Ha,Hb,L[j]) for j in range(len(L))]

def SecondOrderHypermatrixResolutionPartitionII(U, V, Ha, Hb, NbPrts=2):
    """
    outputs the spliting of a second order hypermatrix into the pieces
    as suggested by the resolution of identity. The first three input
    hypermatrices are uncorrelated tuples.
    the last three inputs correspond to the factors for the spliting.

    EXAMPLES:
 
    ::

        sage: [U,V]=GeneralUncorrelatedHypermatrixTuple(2)
        sage: L=SecondOrderHypermatrixResolutionPartitionII(U.matrix(), V.matrix(), HM(2,2,'a').matrix(), HM(2,2,'b').matrix())
        sage: len(L)
        2
        sage: sum(L).simplify_full()
        [a00*b00 + a01*b10 a00*b01 + a01*b11]
        [a10*b00 + a11*b10 a10*b01 + a11*b11]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the list storing the projectors
    L = []
    for i in rg(U.nrows()):
        # Filling up the slice
        Tp0 = U[:,i]
        # Filling up the slice
        Tp1 = V[i,:]
        # Appending the components to the list
        L.append(Tp0*Tp1)
    if len(L)>NbPrts:
        return [Ha*sum(L[j*NbPrts:min((j+1)*NbPrts,len(L))])*Hb for j in range(ceil(len(L)/NbPrts))]
    else:
        return [Ha*L[j]*Hb for j in range(len(L))]

def ThirdOrderHypermatrixResolutionPartition(U, V, W, Ha, Hb, Hc, NbPrts=2):
    """
    outputs the spliting of a third order hypermatrix into the pieces
    as suggested by the resolution of identity. The first three input
    hypermatrices are uncorrelated tuples.
    the last three inputs correspond to the factors for the spliting.

    EXAMPLES:
 
    ::

        sage: [U,V,W]=GeneralUncorrelatedHypermatrixTuple(3)
        sage: L=ThirdOrderHypermatrixResolutionPartition(U,V,W,HM(2,2,2,'a'),HM(2,2,2,'b'),HM(2,2,2,'c'))
        sage: len(L)
        2
        sage: sum(L).simplify_full()
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
        L.append(Prod(Tp0, Tp1, Tp2))
    if len(L)>NbPrts:
        return [ProdB(Ha,Hb,Hc,sum(L[j*NbPrts:min((j+1)*NbPrts,len(L))])) for j in range(ceil(len(L)/NbPrts))]
    else:
        return [ProdB(Ha,Hb,Hc,L[j]) for j in range(len(L))]

def FourthOrderHypermatrixResolutionPartition(Q, U, V, W, Ha, Hb, Hc, Hd):
    """
    outputs the spliting of a fourth order hypermatrix into the pieces
    as suggested by the resolution of identity. The first three input
    hypermatrices are uncorrelated tuples.
    the last three imputs correspond to the factors for the spliting.

    EXAMPLES:
 
    ::

        sage: [Q, U, V, W] = GeneralUncorrelatedHypermatrixTuple(4)
        sage: L = FourthOrderHypermatrixResolutionPartition(Q,U,V,W,HM(2,2,2,2,'a'),HM(2,2,2,2,'b'),HM(2,2,2,2,'c'), HM(2,2,2,2,'d'))
        sage: len(L)
        2
        sage: sum(L).simplify_full()
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

        sage: [Q,U,V] = CanonicalThirdOrderHypermatrixFactors(HM(2,2,2,'a'))
        sage: Prod(Q,U,V)-HM(2,2,2,'a')
        [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]
        

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


def GeneralHypermatrixLogProductTermList(*args):
    """
    Outputs the matrix constraints associated with the Logproduct 
    this came up when I tried to use system of linear algebra
    and resoltion of identity to express a composition rule.
    This was not overall successfull.

    EXAMPLES:
 
    ::

        sage: GeneralHypermatrixLogProductTermList(HM(2,1,'a'),HM(1,2,'b'))
        [a00 + b00, a10 + b00, a00 + b01, a10 + b01]
        

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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # computing the Hypermatrix product
        if len(args)<2:
            raise ValueError("The number of operands must be >= 2")
        elif len(args) >= 2:
            Rh = Rh + [sum([args[s][tuple(entry[0:Integer(mod(s+1,len(args)))]+[t]+entry[Integer(mod(s+2,len(args))):])] for s in range(len(args)-2)] + [args[len(args)-2][tuple(entry[0:len(args)-1]+[t])]]+[args[len(args)-1][tuple([t]+entry[1:])]]) for t in range((args[0]).n(1))]
    return Rh

def BlockSweep(A,k):
    """
    Outputs the result of the sweep operator a  block partion second order hypermatrix

    EXAMPLES:

    ::

        sage: cL=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']; sz=3; sz0=2; sz1=1 
        sage: A=HM(sz,sz,[HM(sz0,sz0,cL[0]),HM(sz1,sz0,cL[1]),HM(sz1,sz0,cL[2]), HM(sz0,sz1,cL[3]),HM(sz1,sz1,cL[4]),HM(sz1,sz1,cL[5]), HM(sz0,sz1,cL[6]),HM(sz1,sz1,cL[7]),HM(sz1,sz1,cL[8])])
        sage: B=A.copy()
        sage: for k in range(B.n(0)):
        ....:     B=BlockSweep(B,k)
        ....:
        sage: (A*B).simplify_full()
        [[[[-1, 0], [0, -1]], [[0], [0]], [[0], [0]]], [[[0, 0]], [[-1]], [[0]]], [[[0, 0]], [[0]], [[-1]]]]


    AUTHORS:
    - Edinah K. Gnang and Jeanine Gnang
    """
    # Initialization of the hypermatrix B
    B=HM(A.n(0),A.n(1), [HM(*(H.dimensions()+['zero'])) for H in A.list()])
    B[k,k] = -A[k,k].inverse()
    for i in range(A.n(0)):
        for j in range(A.n(1)):
            if i != k:
                B[i,k]= A[i,k]*A[k,k].inverse()
            if j != k:
                B[k,j]= A[k,k].inverse()*A[k,j]
            if i != k and j!= k:
                B[i,j] = A[i,j]-(A[i,k]*A[k,k].inverse()*A[k,j])
    return B.copy()

def SweepHM(A,k):
    """
    Outputs the result of the sweep operator on a block partition 
    second order hypermatrix. This version implements the BM
    product take on the sweep operator.

    EXAMPLES:

    ::

        sage: Ha=HM(2,2,'a'); Hb=HM(2,2,'b'); Hc=HM(2,2,'c'); Hd=HM(2,2,'d')
        sage: A=HM([[Ha,Hb],[Hc,Hd]]); B=A.copy()
        sage: for k in range(B.n(0)):
        ....:     B=SweepHM(B,k)[3]
        ....:
        sage: (A*B).simplify_full()
        [[[[-1, 0], [0, -1]], [[0, 0], [0, 0]]], [[[0, 0], [0, 0]], [[-1, 0], [0, -1]]]]
        sage: Ha=HM(2, 2, 'a'); A=HM(2, 2, [HM(1,1,[f]) for f in Ha.list()])
        sage: B=A.copy()
        sage: for k in range(B.n(0)):
        ....:     B=SweepHM(B,k)[3]
        ....:
        sage: (A*B).simplify_full()
        [[[[-1]], [[0]]], [[[0]], [[-1]]]]
        sage: sz=3; indx=1 # Initialization of the size  sweep index parameter
        sage: Ha=HM(sz, sz, 'a'); A=HM(sz, sz, [HM(1, 1, [f]) for f in Ha.list()]); B=A.copy()
        sage: [U, V, W, Rt]=SweepHM(A, indx) # Performing the sweep
        sage: Hu=HM(sz, 2, 1, [U[i,j,0][0,0] for j in range(2) for i in range(sz)]); Hu.printHM()
        [:, :, 0]=
        [  1 a01]
        [  0   1]
        [  1 a21]
        sage: Hv=HM(sz, sz, 2, [V[i,j,k][0,0] for k in range(2) for j in range(sz) for i in range(sz)]); Hv.printHM()
        [:, :, 0]=
        [a00 a01 a02]
        [a10 a11 a12]
        [a20 a21 a22]
        <BLANKLINE>
        [:, :, 1]=
        [-1/a11  1/a11 -1/a11]
        [ 1/a11 -1/a11  1/a11]
        [-1/a11  1/a11 -1/a11]
        <BLANKLINE>
        sage: Hw=HM(2, sz, 1, [W[i,j,0][0,0] for j in range(sz) for i in range(2)]);Hw.printHM()
        [:, :, 0]=
        [  1   0   1]
        [a10   1 a12]
        sage: Prod(Hu,Hv,Hw).printHM()
        [:, :, 0]=
        [ a00 - a01*a10/a11            a01/a11  a02 - a01*a12/a11]
        [           a10/a11             -1/a11            a12/a11]
        [ a20 - a10*a21/a11            a21/a11 -a12*a21/a11 + a22]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Checking that the input hypermatrix is of order 2
    if A.is_cubical() and A[k,k].is_cubical() and A.order()==2:
        # Initializing the hypermatrix U
        U=HM(A.n(0),2,1,[HM(A[i,j].n(0),A[i,j].n(1),'zero') for z in range(1) for j in range(2) for i in range(A.n(0))])
        for i in range(A.n(0)):
            if i!=k:
                U[i,0,0]=HM(2,A[i,k].n(0),'kronecker')
                U[i,1,0]=A[i,k]
            if i==k:
                U[i,0,0]=HM(A[k,k].n(1),A[k,k].n(0),'zero')
                U[i,1,0]=HM(2,A[k,k].n(0),'kronecker')
        # Initializing the hypermatrix W        
        W=HM(2,A.n(1),1,[HM(A[i,j].n(1),A[i,j].n(0),'zero') for z in range(1) for j in range(A.n(1)) for i in range(2)])
        for j in range(A.n(1)):
            if j!=k:
                W[0,j,0]=HM(2,A[k,j].n(1),'kronecker')
                W[1,j,0]=A[k,j]
            if j==k:
                W[0,j,0]=HM(A[k,k].n(1),A[k,k].n(0),'zero')
                W[1,j,0]=HM(2,A[k,k].n(0),'kronecker')
        # Computing the inverse of the entry A[k,k]
        if A[k,k].dimensions()==[1,1]:
            AkkI=A[k,k].inverse()
        else:
            # Part of the code responsible for computing the inverse by sweeping
            AkkI=A[k,k].copy()
            Tmp=HM(AkkI.n(0),AkkI.n(1),'zero')
            for t in range(AkkI.n(0)):
                Tmp[t,t]=-(AkkI[t,t])^(-1)
                for u in range(Tmp.n(0)):
                    for v in range(Tmp.n(1)):
                        if u!=t:
                            Tmp[u,t]=AkkI[u,t]*(AkkI[t,t]^(-1))
                        if v!=t:
                            Tmp[t,v]=(AkkI[t,t]^(-1))*AkkI[t,v]
                        if u!=t and t!=v:
                            Tmp[u,v]=AkkI[u,v]-AkkI[u,t]*(AkkI[t,t]^(-1))*AkkI[t,v]
                AkkI=Tmp.copy()
            AkkI=-AkkI.copy()
        # Initialization of the hypermatrix V
        V=HM(A.n(0),A.n(1),2,[HM(A[i,j].n(0),A[i,j].n(1),'zero') for k in range(2) for j in range(A.n(1)) for i in range(A.n(0))])
        for i in range(V.n(0)):
            for j in range(V.n(1)):
                V[i,j,0]=A[i,j]
                if i==k and j==k:
                    V[i,j,1]=-AkkI
                if i==k and j!=k:
                    V[i,j,1]=AkkI
                if i!=k and j==k:
                    V[i,j,1]=AkkI
                if i!=k and j!=k:
                    V[i,j,1]=-AkkI
        # Returning the BM product.
        return [U,V,W,HM(A.n(0),A.n(1),Prod(U,V,W).list())]
    else:
        raise ValueError("Expected square second order hypermatrix with square diagonal blocks")

def gaussian_elimination(Cf, rs):
    """
    Outputs the row echelon form of the input matrix and the right hand side.

    EXAMPLES:
 
    ::

        sage: [RefA, c] = gaussian_elimination(HM(2,2,'a').matrix(), HM(2,1,'b').matrix())
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
                P=sum([Id[:,k]*Id[Integer(mod(k+1,Ta.nrows())),:] for k in range(Ta.nrows())])
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
    return [A,b]

def gaussian_eliminationHM(Cf, rs):
    """
    Outputs the row echelon form of the input second order hypermatrix and the right hand side.
    This implementation works perfectly well if the entries are matrices or hypermatrices
    associated with appropriate block partitions.

    EXAMPLES:
 
    ::

        sage: [A,b] = gaussian_eliminationHM(HM(2,2,'a'), HM(2,1,'b'))
        sage: A.printHM()
        [:, :]=
        [      1 a01/a00]
        [      0       1]
        sage: b.printHM()
        [:, :]=
        [                                b00/a00]
        [(a10*b00/a00 - b10)/(a01*a10/a00 - a11)]
        sage: Ta=HM(2,2,'a'); Tb=HM(2,1,'b') # Initialization of the factors.
        sage: Ha=HM(2,2,[Ta[0,0]*HM(2,2,'kronecker'), Ta[1,0]*HM(2,2,'kronecker'), Ta[0,1]*HM(2,2,'kronecker'), Ta[1,1]*HM(2,2,'kronecker')])
        sage: Hb=HM(2,1,[Tb[0,0]*HM(2,2,'kronecker'), Tb[1,0]*HM(2,2,'kronecker')])
        sage: [A,b]=gaussian_eliminationHM(Ha,Hb) # performing the gaussian elimination where entries are hypermatrices.
        sage: A
        [[[[1, 0], [0, 1]], [[a01/a00, 0], [0, a01/a00]]], [[[0, 0], [0, 0]], [[1, 0], [0, 1]]]]
        sage: b
        [[[[b00/a00, 0], [0, b00/a00]]], [[[(a10*b00/a00 - b10)/(a01*a10/a00 - a11), 0], [0, (a10*b00/a00 - b10)/(a01*a10/a00 - a11)]]]]
        sage: Ta=HM(2,2,'a'); Tb=HM(2,1,HM(2,'b').list())
        sage: Ha=HM(2,2,[Ta[0,0]*HM(2,3,'kronecker'), Ta[1,0]*HM(2,3,'c'), Ta[0,1]*HM(3,2,'d'), Ta[1,1]*HM(2,2,'kronecker')])
        sage: Hb=HM(2,1,[Tb[0,0]*HM(2,3,'kronecker'), Tb[1,0]*HM(2,3,'zero')])
        sage: [A,b]=gaussian_eliminationHM(Ha,Hb)
        sage: A[0,0].printHM()
        [:, :]=
        [1 0 0]
        [0 1 0]
        [0 0 1]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing a copy of the input second order hypermatrices.
    A=Cf.copy(); b=rs.copy()
    # Initialization of the row and column index
    i=0; j=0
    while i < A.n(0) and j < A.n(1):
        while A.slice(rg(i,A.n(0)),'row').slice([j],'col').is_zero() and j < A.n(1)-1:
            # Incrementing the column index
            j=j+1
        if A.slice(rg(i,A.n(0)),'row').is_zero()==False:
            while A[i,j].is_zero(): 
                Ta=A.slice(rg(i,A.n(0)),'row')
                Tb=b.slice(rg(i,b.n(0)),'row')
                # Initializing the cyclic shift permutation matrix
                Id=HM(2, Ta.n(0), 'kronecker')
                P=Matrix2HM(sum([Id.matrix()[:,k]*Id.matrix()[Integer(mod(k+1,Ta.nrows())),:] for k in rg(Ta.n(0))]))
                Ta=P*Ta; Tb=P*Tb
                for i0 in range(Ta.n(0)):
                    for j0 in range(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
                for i0 in range(Tb.n(0)):
                    for j0 in range(Tb.n(1)):
                        b[i+i0,j0]=Tb[i0,j0]
            # Performing the row operations.
            cf1=A[i,j]
            for j0 in range(b.n(1)):
                b[i,j0]=(cf1^(-1))*b[i,j0]
            for j0 in range(A.n(1)):
                A[i,j0]=(cf1^(-1))*A[i,j0]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if HM(1,A.n(1),[A[r,j0] for j0 in range(A.n(1))]).is_zero():
                    r=r+1
                else:
                    # Initialization of the coefficient
                    cf2=A[r,j]
                    for j0 in range(b.n(1)):
                        b[r,j0]=-cf2*b[i,j0]+b[r,j0]
                    for j0 in range(A.n(1)):
                        A[r,j0]=-cf2*A[i,j0]+A[r,j0]
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return [A,b]

def gaussian_eliminationHMII(Cf, rs):
    """
    Outputs the row echelon form of the input second order hypermatrix and the right hand side.
    does not normalize the rows to ensure that the first non zero entry of non zero rows = 1
    This implementation tacitly assumes that the entries commute.
 

    EXAMPLES:
 
    ::

        sage: [A,b]=gaussian_eliminationHMII(HM(2,2,'a'), HM(2,1,'b'))
        sage: A.printHM()
        [:, :]=
        [              a00               a01]
        [                0 a01*a10 - a00*a11]
        sage: b.printHM()
        [:, :]=
        [              b00]
        [a10*b00 - a00*b10]
        sage: Ta=HM(2,2,'a'); Tb=HM(2,1,'b')
        sage: Ha=HM(2,2,[Ta[0,0]*HM(2,2,'kronecker'), Ta[1,0]*HM(2,2,'kronecker'), Ta[0,1]*HM(2,2,'kronecker'), Ta[1,1]*HM(2,2,'kronecker')])
        sage: Hb=HM(2,1,[Tb[0,0]*HM(2,2,'kronecker'), Tb[1,0]*HM(2,2,'kronecker')])
        sage: [A,b]=gaussian_eliminationHMII(Ha,Hb)
        sage: A
        [[[[a00, 0], [0, a00]], [[a01, 0], [0, a01]]], [[[0, 0], [0, 0]], [[a01*a10 - a00*a11, 0], [0, a01*a10 - a00*a11]]]]
        sage: b
        [[[[b00, 0], [0, b00]]], [[[a10*b00 - a00*b10, 0], [0, a10*b00 - a00*b10]]]]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing a copy of the input second order hypermatrices.
    A=Cf.copy(); b=rs.copy()
    # Initialization of the row and column index
    i=0; j=0
    while i < A.n(0) and j < A.n(1):
        while HM(A.n(0)-i, 1, [A[i0,j] for i0 in range(i,A.n(0))]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        if HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))]).is_zero()==False:
            while A[i,j].is_zero(): 
                Ta=HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))])
                Tb=HM(b.n(0)-i, b.n(1), [b[i0,j0] for j0 in range(b.n(1)) for i0 in range(i,b.n(0))])
                # Initializing the cyclic shift permutation matrix
                Id=HM(2, Ta.n(0), 'kronecker')
                P=sum([HM(Ta.n(0),1,[Id[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Id[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                Ta=P*Ta; Tb=P*Tb
                for i0 in range(Ta.n(0)):
                    for j0 in range(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
                for i0 in range(Tb.n(0)):
                    for j0 in range(Tb.n(1)):
                        b[i+i0,j0]=Tb[i0,j0]
            # Performing the row operations.
            cf1=A[i,j]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if HM(1,A.n(1),[A[r,j0] for j0 in range(A.n(1))]).is_zero():
                    r=r+1
                else:
                    # Initialization of the coefficient
                    cf2=A[r,j]
                    for j0 in range(b.n(1)):
                        b[r,j0]=cf2*b[i,j0]-cf1*b[r,j0]
                    for j0 in range(A.n(1)):
                        A[r,j0]=cf2*A[i,j0]-cf1*A[r,j0]
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return [A,b]

def gaussian_eliminationHMIII(Cf, rs):
    """
    Outputs the row echelon form of the input second order hypermatrix and the right hand side.
    does not normalize the rows to ensure that the first non zero entry of non zero rows = 1.
    The difference with the previous implementation is the fact that the row linear combination
    operations are performed in such a way as to not change the absolute value of the determinant.
    This implementation is skew fields or division ring friendly by inputing hypermatrices
    whose entries are themselve hypermatrices of the approprioate size. 


    EXAMPLES:
 
    ::

        sage: [A,b]=gaussian_eliminationHMIII(HM(2,2,'a'), HM(2,1,'b'))
        sage: A.printHM()
        [:, :]=
        [               a00                a01]
        [                 0 -a01*a10/a00 + a11]
        sage: b.printHM()
        [:, :]=
        [               b00]
        [-a10*b00/a00 + b10]
        sage: Ta=HM(2,2,'a'); Tb=HM(2,1,'b')
        sage: Ha=HM(2,2,[Ta[0,0]*HM(2,2,'kronecker'), Ta[1,0]*HM(2,2,'kronecker'), Ta[0,1]*HM(2,2,'kronecker'), Ta[1,1]*HM(2,2,'kronecker')])
        sage: Hb=HM(2,1,[Tb[0,0]*HM(2,2,'kronecker'), Tb[1,0]*HM(2,2,'kronecker')])
        sage: [A,b]=gaussian_eliminationHMIII(Ha,Hb)
        sage: A
        [[[[a00, 0], [0, a00]], [[a01, 0], [0, a01]]], [[[0, 0], [0, 0]], [[-a01*a10/a00 + a11, 0], [0, -a01*a10/a00 + a11]]]]
        sage: b
        [[[[b00, 0], [0, b00]]], [[[-a10*b00/a00 + b10, 0], [0, -a10*b00/a00 + b10]]]]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing a copy of the input second order hypermatrices.
    A=Cf.copy(); b=rs.copy()
    # Initialization of the row and column index
    i=0; j=0
    while i < A.n(0) and j < A.n(1):
        while HM(A.n(0)-i, 1, [A[i0,j] for i0 in range(i,A.n(0))]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        if HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))]).is_zero()==False:
            while A[i,j].is_zero(): 
                Ta=HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))])
                Tb=HM(b.n(0)-i, b.n(1), [b[i0,j0] for j0 in range(b.n(1)) for i0 in range(i,b.n(0))])
                # Initializing the cyclic shift permutation matrix
                Id=HM(2, Ta.n(0), 'kronecker')
                P=sum([HM(Ta.n(0),1,[Id[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Id[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                Ta=P*Ta; Tb=P*Tb
                for i0 in range(Ta.n(0)):
                    for j0 in range(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
                for i0 in range(Tb.n(0)):
                    for j0 in range(Tb.n(1)):
                        b[i+i0,j0]=Tb[i0,j0]
            if A.n(0)-i-1 > 0 and not (HM(A.n(0)-i-1, 1, [A[i0,j] for i0 in range(i+1,A.n(0))]).is_zero() and j <= A.ncols()-1):
                # Performing the row operations.
                cf1=A[i,j]
                for r in range(i+1,A.nrows()):
                    # Taking care of the zero row
                    if HM(1,A.n(1),[A[r,j0] for j0 in range(A.n(1))]).is_zero():
                        r=r+1
                    else:
                        # Initialization of the coefficient
                        cf2=A[r,j]
                        for j0 in range(b.n(1)):
                            b[r,j0]=-(cf2*cf1^(-1))*b[i,j0]+b[r,j0]
                        for j0 in range(A.n(1)):
                            A[r,j0]=-(cf2*cf1^(-1))*A[i,j0]+A[r,j0]
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return [A,b]

def gaussian_eliminationHMIV(Cf):
    """
    Outputs the row echelon form of the input second order hypermatrix and the right hand side.
    does not normalize the rows to ensure that the first non zero entry of non zero rows = 1
    This implementation tacitly assumes that the entries commute.
    The computation is performed in such a way that the last entry on the main diagonal is 
    the determinant.
     

    EXAMPLES:
 
    ::

        sage: A=gaussian_eliminationHMIV(HM(2,2,'a'))
        sage: A.printHM()
        [:, :]=
        [               a00                a01]
        [                 0 -a01*a10 + a00*a11]
        sage: Ta=HM(2,2,'a'); Tb=HM(2,1,'b')
        sage: Ha=HM(2,2,[Ta[0,0]*HM(2,2,'kronecker'), Ta[1,0]*HM(2,2,'kronecker'), Ta[0,1]*HM(2,2,'kronecker'), Ta[1,1]*HM(2,2,'kronecker')])
        sage: Hb=HM(2,1,[Tb[0,0]*HM(2,2,'kronecker'), Tb[1,0]*HM(2,2,'kronecker')])
        sage: A=gaussian_eliminationHMIV(Ha)
        sage: A
        [[[[a00, 0], [0, a00]], [[a01, 0], [0, a01]]], [[[0, 0], [0, 0]], [[-a01*a10 + a00*a11, 0], [0, -a01*a10 + a00*a11]]]]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing a copy of the input second order hypermatrices.
    A=Cf.copy()
    # Initialization of the row and column index
    i=0; j=0
    while i < A.n(0) and j < A.n(1):
        while HM(A.n(0)-i, 1, [A[i0,j] for i0 in range(i,A.n(0))]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        if HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))]).is_zero()==False:
            while A[i,j].is_zero(): 
                Ta=HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))])
                # Initializing the cyclic shift permutation matrix
                Id=HM(2, Ta.n(0), 'kronecker')
                P=sum([HM(Ta.n(0),1,[Id[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Id[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                Ta=P*Ta
                for i0 in range(Ta.n(0)):
                    for j0 in range(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
            # Performing the row operations.
            cf1=A[i,j]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if HM(1,A.n(1),[A[r,j0] for j0 in range(A.n(1))]).is_zero():
                    r=r+1
                else:
                    # Initialization of the coefficient
                    cf2=A[r,j]
                    if r==j+1:
                        for j0 in range(j,A.n(1)):
                            A[r,j0]=-cf2*A[i,j0]+cf1*A[r,j0]
                    else:
                        for j0 in range(j,A.n(1)):
                            A[r,j0]=(-cf2*A[i,j0]+cf1*A[r,j0])/cf1
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return A

def gaussian_elimination_ReductionHM(Cf, rs, VrbL, Rlts):
    """
    Outputs the row echelon form of the input second order hypermatrix and the right hand side.
    does not normalize the rows to ensure that the first non zero entry of non zero rows = 1
    The algorithm perform the reduction assuming that the the leading term in each relations
    is a monic powers in a distinct variable as illustrated in the example bellow.

    EXAMPLES:
 
    ::

        sage: x1, x2=var('x1, x2')
        sage: Cf=HM([[-2, -2*x1 + 3], [-12*x1 + 10, -2*x1 + 3]])
        sage: rs=HM(2,1,'zero')
        sage: VrbL=[x1, x2]
        sage: Rlts=[x1^2 - 3*x1 + 2, x2^2 - 3*x2 + 2]
        sage: [A,b]=gaussian_elimination_ReductionHM(Cf, rs, VrbL, Rlts)
        sage: A.printHM()
        [:, :]=
        [        -2  -2*x1 + 3]
        [         0 12*x1 - 12]
        sage: b.printHM()
        [:, :]=
        [0]
        [0]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing a copy of the input second order hypermatrices.
    A=Cf.copy(); b=rs.copy()
    # Initialization of the row and column index
    i=0; j=0
    while i < A.n(0) and j < A.n(1):
        while HM(A.n(0)-i, 1, [A[i0,j] for i0 in range(i,A.n(0))]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        if HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))]).is_zero()==False:
            while A[i,j].is_zero(): 
                Ta=HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))])
                Tb=HM(b.n(0)-i, b.n(1), [b[i0,j0] for j0 in range(b.n(1)) for i0 in range(i,b.n(0))])
                # Initializing the cyclic shift permutation matrix
                Id=HM(2, Ta.n(0), 'kronecker')
                P=sum([HM(Ta.n(0),1,[Id[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Id[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                Ta=P*Ta; Tb=P*Tb
                for i0 in range(Ta.n(0)):
                    for j0 in range(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
                for i0 in range(Tb.n(0)):
                    for j0 in range(Tb.n(1)):
                        b[i+i0,j0]=Tb[i0,j0]
            # Performing the row operations.
            cf1=A[i,j]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if HM(1,A.n(1),[A[r,j0] for j0 in range(A.n(1))]).is_zero():
                    r=r+1
                else:
                    # Initialization of the coefficient
                    cf2=A[r,j]
                    for j0 in range(b.n(1)):
                        # Performing the reduction
                        #b[r,j0]=cf2*b[i,j0]-cf1*b[r,j0]
                        f=expand(cf2*b[i,j0]-cf1*b[r,j0])
                        for v in range(len(VrbL)):
                            for d in range(f.degree(VrbL[v])-Rlts[v].degree(VrbL[v]),-1,-1):
                                f=expand(fast_reduce(f,[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))],[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))-expand(Rlts[v]*VrbL[v]^d)])) 
                        b[r,j0]=f
                    for j0 in range(A.n(1)):
                        # Performing the reduction
                        #A[r,j0]=cf2*A[i,j0]-cf1*A[r,j0]
                        g=expand(cf2*A[i,j0]-cf1*A[r,j0])
                        for v in range(len(VrbL)):
                            for d in range(g.degree(VrbL[v])-Rlts[v].degree(VrbL[v]),-1,-1):
                                g=expand(fast_reduce(g,[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))],[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))-expand(Rlts[v]*VrbL[v]^d)]))
                        A[r,j0]=g
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return [A,b]

def gaussian_elimination_ReductionHMII(Cf, VrbL, Rlts):
    """
    Outputs the row echelon form of the input second order hypermatrix and the right hand side.
    does not normalize the rows to ensure that the first non zero entry of non zero rows = 1
    The algorithm perform the reduction assuming that the the leading term in each relations
    is a monic powers in a distinct variable as illustrated in the example bellow.
    The computation is perform in such a way that the last diagonal entry holds the determinant
    of the whole matrix it is also true that each diagonal entry corresponds to the determinant
    of the corresponding top diagonal block of the matrix. Note that the relation in Rlts are
    assumed to be univariate leading term. This implementaion also gets the sign right for the
    determinant. This implementation is considerably more efficient then the previous one above
    because the extra multiplicative factor is kept minimal.


    EXAMPLES:
 
    ::

        sage: x1, x2=var('x1, x2')
        sage: Cf=HM([[-2, -2*x1 + 3], [-12*x1 + 10, -2*x1 + 3]])
        sage: VrbL=[x1, x2]
        sage: Rlts=[x1^2 - 3*x1 + 2, x2^2 - 3*x2 + 2]
        sage: A=gaussian_elimination_ReductionHMII(Cf, VrbL, Rlts)
        sage: A.printHM()
        [:, :]=
        [         -2   -2*x1 + 3]
        [          0 -12*x1 + 12]
        sage: od=2; sz=3 # Initialization of the order and size parameter
        sage: A=HM(od,sz,'a','sym'); X=HM(sz,sz,[x^(sz^abs(j-i)) for j in rg(sz) for i in rg(sz)])
        sage: Hb=HM(sz,binomial(sz,2),'zero'); clidx=0 # Initialization of the incidence matrix
        sage: for i in rg(sz):
        ....:     for j in rg(sz):
        ....:         if i < j:
        ....:             Hb[i,clidx]=-sqrt(A[i,j]*X[i,j])
        ....:             Hb[j,clidx]=+sqrt(A[i,j]*X[i,j])
        ....:             clidx=clidx+1
        ....:
        sage: t=0; B=HM(sz-1,Hb.n(1),[Hb[i,j] for j in rg(Hb.n(1)) for i in rg(Hb.n(0)) if i!=t]) # Grounding at the vertex t
        sage: M=(HM(od,[x*A[t,t]]+[1 for i in rg(sz-2)],'diag')*(B*B.transpose())).expand() # Initialization of the fundamental matrix
        sage: d=1+sum(sz^k for k in rg(1+floor((sz-1)/2),sz))+sum(sz^(sz-1-k) for k in rg(1,1+floor((sz-1)/2))) # Initializing the max degree
        sage: Rh=gaussian_elimination_ReductionHMII(M,[x],[x^(1+d)])
        sage: Rh[1,1]
        a00*a01*a02*x^13 + a00*a02*a12*x^13 + a00*a01*a12*x^7


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing a copy of the input second order hypermatrices.
    A=Cf.copy()
    # Initialization of the row and column index
    i=0; j=0
    while i < A.n(0) and j < A.n(1):
        while HM(A.n(0)-i, 1, [A[i0,j] for i0 in range(i,A.n(0))]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        if HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))]).is_zero()==False:
            while A[i,j].is_zero(): 
                Ta=HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))])
                # Initializing the cyclic shift permutation matrix
                Id=HM(2, Ta.n(0), 'kronecker')
                P=sum([HM(Ta.n(0),1,[Id[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Id[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                Ta=P*Ta
                for i0 in range(Ta.n(0)):
                    for j0 in range(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
            # Performing the row operations.
            cf1=A[i,j]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if HM(1,A.n(1),[A[r,j0] for j0 in range(A.n(1))]).is_zero():
                    r=r+1
                else:
                    # Initialization of the coefficient
                    cf2=A[r,j]
                    if r==j+1:
                        for j0 in range(j,A.n(1)):
                            # Performing the reduction
                            f=(-cf2*A[i,j0]+cf1*A[r,j0]).numerator()
                            for v in range(len(VrbL)):
                                #f=f.maxima_methods().divide(Rlts[v])[1]
                                for d in range(f.degree(VrbL[v])-Rlts[v].degree(VrbL[v]),-1,-1):
                                    f=expand(fast_reduce(f,[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))],[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))-expand(Rlts[v]*VrbL[v]^d)])) 
                            A[r,j0]=f
                    else:
                        for j0 in range(j,A.n(1)):
                            # Performing the reduction
                            g=expand((-cf2*A[i,j0]+cf1*A[r,j0])/cf1)
                            A[r,j0]=g
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return A

def gaussian_elimination_ReductionHMIII(Cf, VrbL, Rlts, RltsII):
    """
    Outputs the row echelon form of the input second order hypermatrix and the right hand side.
    does not normalize the rows to ensure that the first non zero entry of non zero rows = 1
    The algorithm perform the reduction assuming that the the leading term in each relations
    is a monic powers in a distinct variable as illustrated in the example bellow.
    The computation is perform in such a way that the last diagonal entry holds the determinant
    of the whole matrix it is also true that each diagonal entry corresponds to the determinant
    of the corresponding top diagonal block of the matrix. Note that the relation in Rlts are
    assumed to be univariate leading term. This implementaion also gets the sign right for the
    determinant. This implementation is considerably more efficient then the previous one above
    because the extra multiplicative factor is kept minimal.


    EXAMPLES:
 
    ::

        sage: x1, x2=var('x1, x2')
        sage: Cf=HM([[-2, -2*x1 + 3], [-12*x1 + 10, -2*x1 + 3]])
        sage: VrbL=[x1, x2]
        sage: Rlts=[x1^2 - 3*x1 + 2, x2^2 - 3*x2 + 2]
        sage: A=gaussian_elimination_ReductionHMIII(Cf, VrbL, Rlts, Rlts)
        sage: A.printHM()
        [:, :]=
        [         -2   -2*x1 + 3]
        [          0 -12*x1 + 12]
        sage: od=2; sz=3 # Initialization of the order and size parameter
        sage: A=HM(od,sz,'a','sym'); X=HM(sz,sz,[x^(sz^abs(j-i)) for j in rg(sz) for i in rg(sz)])
        sage: Hb=HM(sz,binomial(sz,2),'zero'); clidx=0 # Initialization of the incidence matrix
        sage: for i in rg(sz):
        ....:     for j in rg(sz):
        ....:         if i < j:
        ....:             Hb[i,clidx]=-sqrt(A[i,j]*X[i,j])
        ....:             Hb[j,clidx]=+sqrt(A[i,j]*X[i,j])
        ....:             clidx=clidx+1
        ....:
        sage: t=0; B=HM(sz-1,Hb.n(1),[Hb[i,j] for j in rg(Hb.n(1)) for i in rg(Hb.n(0)) if i!=t]) # Grounding at the vertex t
        sage: M=(HM(od,[x*A[t,t]]+[1 for i in rg(sz-2)],'diag')*(B*B.transpose())).expand() # Initialization of the fundamental matrix
        sage: d=1+sum(sz^k for k in rg(1+floor((sz-1)/2),sz))+sum(sz^(sz-1-k) for k in rg(1,1+floor((sz-1)/2))) # Initializing the max degree
        sage: Rh=gaussian_elimination_ReductionHMIII(M,[x],[x^(1+d)],[x^(1+d)])
        sage: Rh[1,1]
        a00*a01*a02*x^13 + a00*a02*a12*x^13 + a00*a01*a12*x^7
        sage: od=2; sz=4 # Initialization of the order and size parameter
        sage: Lx=var_list('x',sz) # Initialization of the list of vraiables
        sage: A=HM(od,sz,'a','sym') # Initialization of the adjacency matrix
        sage: X=HM(sz,sz,[Lx[abs(i-j)] for j in rg(sz) for i in rg(sz)]) # initialization of the edge weight matrix
        sage: Hb=HM(sz,binomial(sz,2),'zero'); clidx=0 # Initialization of the incidence matrix
        sage: for i in rg(sz):
        ....:     for j in rg(sz):
        ....:         if i < j:
        ....:             Hb[i,clidx]=-sqrt(A[i,j]*X[i,j])
        ....:             Hb[j,clidx]=+sqrt(A[i,j]*X[i,j])
        ....:             clidx=clidx+1
        ....:
        sage: t=0; B=HM(sz-1,Hb.n(1),[Hb[i,j] for j in rg(Hb.n(1)) for i in rg(Hb.n(0)) if i!=t]) # Grounding at the vertex t
        sage: M=(HM(od,[Lx[0]*A[t,t]]+[1 for i in rg(sz-2)],'diag')*(B*B.transpose())).expand() # Initialization of the fundamental matrix
        sage: d0=3; d1=2; Rh=gaussian_elimination_ReductionHMIII(M,Lx,[v^d0 for v in Lx],[v^d1 for v in Lx])
        sage: Rh[sz-2,sz-2]
        a00*a01*a02*a03*x0*x1*x2*x3 + a00*a02*a03*a12*x0*x1*x2*x3 + a00*a03*a12*a13*x0*x1*x2*x3 + a00*a03*a13*a23*x0*x1*x2*x3
        sage: od=2; sz=4 # Initialization of the order and size parameter
        sage: Lx=var_list('x',sz) # Initialization of the list of variables
        sage: A=HM(sz,sz,'a') # Initialization of part of the adjacency matrix
        sage: X=HM(sz,sz,[Lx[abs(j-i)] for j in rg(sz) for i in rg(sz)]) # Initialization of the symbolic edge weights
        sage: Hb=Matrix2HM((HM(od,(A.elementwise_product(X)*HM(sz,1,'one')).list(),'diag')-A.elementwise_product(X)).matrix()[1:,1:]) # Directed Laplacian
        sage: d0=3; d1=2; Rh2=A[0,0]*Lx[0]*gaussian_elimination_ReductionHMIII(Hb,Lx,[v^3 for v in Lx],[v^2 for v in Lx])
        sage: Rh2[sz-2,sz-2]
        (a10*a20*a30*x1*x2*x3 + a12*a20*a30*x1*x2*x3 + a13*a21*a30*x1*x2*x3 + a13*a23*a30*x1*x2*x3)*a00*x0


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing a copy of the input second order hypermatrices.
    A=Cf.copy()
    # Initialization of the row and column index
    i=0; j=0
    while i < A.n(0) and j < A.n(1):
        while HM(A.n(0)-i, 1, [A[i0,j] for i0 in range(i,A.n(0))]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        if HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))]).is_zero()==False:
            while A[i,j].is_zero(): 
                Ta=HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))])
                # Initializing the cyclic shift permutation matrix
                Id=HM(2, Ta.n(0), 'kronecker')
                P=sum([HM(Ta.n(0),1,[Id[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Id[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                Ta=P*Ta
                for i0 in range(Ta.n(0)):
                    for j0 in range(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
            # Performing the row operations.
            cf1=A[i,j]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if HM(1,A.n(1),[A[r,j0] for j0 in range(A.n(1))]).is_zero():
                    r=r+1
                else:
                    # Initialization of the coefficient
                    cf2=A[r,j]
                    if r==j+1 and r!=A.n(0)-1:
                        for j0 in range(j,A.n(1)):
                            # Performing the reduction
                            f=(-cf2*A[i,j0]+cf1*A[r,j0]).numerator()
                            for v in range(len(VrbL)):
                                #f=f.maxima_methods().divide(Rlts[v])[1]
                                for d in range(f.degree(VrbL[v])-Rlts[v].degree(VrbL[v]),-1,-1):
                                    f=expand(fast_reduce(f,[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))],[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))-expand(Rlts[v]*VrbL[v]^d)])) 
                            A[r,j0]=f
                    elif r==j+1 and r==A.n(0)-1:
                        for j0 in range(j,A.n(1)):
                            # Performing the reduction
                            f=(-cf2*A[i,j0]+cf1*A[r,j0]).numerator()
                            for v in range(len(VrbL)):
                                #f=f.maxima_methods().divide(Rlts[v])[1]
                                for d in range(f.degree(VrbL[v])-RltsII[v].degree(VrbL[v]),-1,-1):
                                    f=expand(fast_reduce(f,[VrbL[v]^(d+RltsII[v].degree(VrbL[v]))],[VrbL[v]^(d+RltsII[v].degree(VrbL[v]))-expand(RltsII[v]*VrbL[v]^d)])) 
                            A[r,j0]=f
                    else:
                        for j0 in range(j,A.n(1)):
                            # Performing the reduction
                            g=expand((-cf2*A[i,j0]+cf1*A[r,j0])/cf1)
                            A[r,j0]=g
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return A

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

def gauss_jordan_eliminationHM(Cf,rs):
    """
    Outputs the reduced row echelon form of the input matrix and the right hand side.
    This implementation is skew field friendly as illustrated in some of the examples
    below. We do not assume that the input entries commute. 


    EXAMPLES:
 
    ::

        sage: [RefA, c] = gauss_jordan_eliminationHM(HM(2,2,'a'), HM(2,1,'b'))
        sage: RefA.printHM()
        [:, :]=
        [1 0]
        [0 1]
        sage: c.printHM()
        [:, :]=
        [-a01*(a10*b00/a00 - b10)/(a00*(a01*a10/a00 - a11)) + b00/a00]
        [                     (a10*b00/a00 - b10)/(a01*a10/a00 - a11)]
        sage: Ta=HM(2,2,'a'); Tb=HM(2,1,HM(2,'b').list()) # Initialization of the factors.
        sage: Ha=HM(2,2,[Ta[0,0]*HM(2,2,'kronecker'), Ta[1,0]*HM(2,2,'kronecker'), Ta[0,1]*HM(2,2,'kronecker'), Ta[1,1]*HM(2,2,'kronecker')])
        sage: Hb=HM(2,1,[Tb[0,0]*HM(2,2,'kronecker'), Tb[1,0]*HM(2,2,'kronecker')])
        sage: [A,b]=gauss_jordan_eliminationHM(Ha,Hb) # performing the gaussian elimination where entries are hypermatrices.
        sage: A
        [[[[1, 0], [0, 1]], [[0, 0], [0, 0]]], [[[0, 0], [0, 0]], [[1, 0], [0, 1]]]]
        sage: b
        [[[[-a01*(a10*b0/a00 - b1)/(a00*(a01*a10/a00 - a11)) + b0/a00, 0], [0, -a01*(a10*b0/a00 - b1)/(a00*(a01*a10/a00 - a11)) + b0/a00]]], [[[(a10*b0/a00 - b1)/(a01*a10/a00 - a11), 0], [0, (a10*b0/a00 - b1)/(a01*a10/a00 - a11)]]]]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    [A, b] = gaussian_eliminationHM(Cf, rs)
    # Initialization of the row and column index
    i=A.nrows()-1; j=0
    while i>0 or j>0:
        #if (A[i,:]).is_zero():
        if HM(1,A.n(1),[A[i,j0] for j0 in range(A.n(1))]).is_zero():
            # decrementing the row index and initializing the column index
            i=i-1; j=0
        else :
            while (A[i,j]).is_zero():
                # Incrementing the column index
                j = j + 1
            # performing row operations
            for r in range(i-1,-1,-1):
                #b[r,:] = -A[r,j]*b[i,:]+b[r,:]
                #Tb0=HM(1, b.n(1), [b[i,j0] for j0 in range(b.n(1))])
                #Tb1=HM(1, b.n(1), [b[r,j0] for j0 in range(b.n(1))])
                Trb=-HM(1, b.n(1), [A[r,j]*b[i,j0] for j0 in range(b.n(1))]) + HM(1, b.n(1), [b[r,j0] for j0 in range(b.n(1))])
                for j0 in range(b.n(1)):
                    b[r,j0]=Trb[0,j0]
                #A[r,:] = -A[r,j]*A[i,:]+A[r,:]
                #Ta0=HM(1, A.n(1), [A[i,j0] for j0 in range(A.n(1))])
                #Ta1=HM(1, A.n(1), [A[r,j0] for j0 in range(A.n(1))])
                Tra=-HM(1, A.n(1), [A[r,j]*A[i,j0] for j0 in range(A.n(1))]) + HM(1, A.n(1), [A[r,j0] for j0 in range(A.n(1))])
                for j0 in range(A.n(1)):
                    A[r,j0]=Tra[0,j0]
            i=i-1; j=0
    return [A,b]

def gauss_jordan_eliminationHMII(Cf,rs):
    """
    Outputs the reduced row echelon form of the input matrix and the right hand side.
    This implementation assumes that the input entries commute and is therefore NOT
    skew field friendly


    EXAMPLES:
 
    ::

        sage: [RefA, c] = gauss_jordan_eliminationHMII(HM(2,2,'a'), HM(2,1,'b'))
        sage: RefA.printHM()
        [:, :]=
        [-(a01*a10 - a00*a11)*a00                        0]
        [                       0        a01*a10 - a00*a11]
        sage: c.printHM()
        [:, :]=
        [(a10*b00 - a00*b10)*a01 - (a01*a10 - a00*a11)*b00]
        [                                a10*b00 - a00*b10] 
        sage: Ta=HM(2,2,'a'); Tb=HM(2,1,'b')
        sage: Ha=HM(2,2,[Ta[0,0]*HM(2,2,'kronecker'), Ta[1,0]*HM(2,2,'kronecker'), Ta[0,1]*HM(2,2,'kronecker'), Ta[1,1]*HM(2,2,'kronecker')])
        sage: Hb=HM(2,1,[Tb[0,0]*HM(2,2,'kronecker'), Tb[1,0]*HM(2,2,'kronecker')])
        sage: [A,b]=gauss_jordan_eliminationHMII(Ha,Hb)
        sage: A
        [[[[-(a01*a10 - a00*a11)*a00, 0], [0, -(a01*a10 - a00*a11)*a00]], [[0, 0], [0, 0]]], [[[0, 0], [0, 0]], [[a01*a10 - a00*a11, 0], [0, a01*a10 - a00*a11]]]]
        sage: b
        [[[[(a10*b00 - a00*b10)*a01 - (a01*a10 - a00*a11)*b00, 0], [0, (a10*b00 - a00*b10)*a01 - (a01*a10 - a00*a11)*b00]]], [[[a10*b00 - a00*b10, 0], [0, a10*b00 - a00*b10]]]]



    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    [A, b] = gaussian_eliminationHMII(Cf,rs)
    # Initialization of the row and column index
    i=A.nrows()-1; j=0
    while i>0 or j>0:
        #if (A[i,:]).is_zero():
        if HM(1,A.n(1),[A[i,j0] for j0 in range(A.n(1))]).is_zero():
            # decrementing the row index and initializing the column index
            i=i-1; j=0
        else :
            while (A[i,j]).is_zero():
                # Incrementing the column index
                j = j + 1
            # performing row operations
            cf1=A[i,j]
            for r in range(i-1,-1,-1):
                #b[r,:] = -A[r,j]*b[i,:]+b[r,:]
                cf2=A[r,j]
                for j0 in range(b.n(1)):
                    b[r,j0]=cf2*b[i,j0]-cf1*b[r,j0]
                #A[r,:] = -A[r,j]*A[i,:]+A[r,:]
                for j0 in range(A.n(1)):
                    A[r,j0]=cf2*A[i,j0]-cf1*A[r,j0]
            i=i-1; j=0
    return [A,b]

def gauss_jordan_elimination_ReductionHM(Cf, rs, VrbL, Rlts):
    """
    Outputs the reduced row echelon form of the input matrix and the right hand side.

    EXAMPLES:
 
    ::

        sage: x1, x2=var('x1, x2')
        sage: Cf=HM([[-2, -2*x1 + 3], [-12*x1 + 10, -2*x1 + 3]])
        sage: rs=HM(2,1,'zero')
        sage: VrbL=[x1, x2]
        sage: Rlts=[x1^2 - 3*x1 + 2, x2^2 - 3*x2 + 2]
        sage: [A,b]=gauss_jordan_elimination_ReductionHM(Cf, rs, VrbL, Rlts)
        sage: A.printHM()
        [:, :]=
        [24*x1 - 24          0]
        [         0 12*x1 - 12]        
        sage: b.printHM()
        [:, :]=
        [0]
        [0]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    [A, b] = gaussian_elimination_ReductionHM(Cf, rs, VrbL, Rlts)
    # Initialization of the row and column index
    i=A.nrows()-1; j=0
    while i>0 or j>0:
        #if (A[i,:]).is_zero():
        if HM(1,A.n(1),[A[i,j0] for j0 in range(A.n(1))]).is_zero():
            # decrementing the row index and initializing the column index
            i=i-1; j=0
        else :
            while (A[i,j]).is_zero():
                # Incrementing the column index
                j = j + 1
            # performing row operations
            cf1=A[i,j]
            for r in range(i-1,-1,-1):
                #b[r,:] = -A[r,j]*b[i,:]+b[r,:]
                cf2=A[r,j]
                for j0 in range(b.n(1)):
                    #b[r,j0]=cf2*b[i,j0]-cf1*b[r,j0]
                    f=expand(cf2*b[i,j0]-cf1*b[r,j0])
                    for v in range(len(VrbL)):
                        for d in range(f.degree(VrbL[v])-Rlts[v].degree(VrbL[v]),-1,-1):
                            f=expand(fast_reduce(f,[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))],[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))-expand(Rlts[v]*VrbL[v]^d)])) 
                    b[r,j0]=f
                #A[r,:] = -A[r,j]*A[i,:]+A[r,:]
                for j0 in range(A.n(1)):
                    #A[r,j0]=cf2*A[i,j0]-cf1*A[r,j0]
                    g=expand(cf2*A[i,j0]-cf1*A[r,j0])
                    for v in range(len(VrbL)):
                        for d in range(g.degree(VrbL[v])-Rlts[v].degree(VrbL[v]),-1,-1):
                            g=expand(fast_reduce(g,[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))],[VrbL[v]^(d+Rlts[v].degree(VrbL[v]))-expand(Rlts[v]*VrbL[v]^d)]))
                    A[r,j0]=g
            i=i-1; j=0
    return [A,b]

def gauss_jordan_eliminationHMIII(Cf,rs):
    """
    Outputs the reduced row echelon form of the input second order hypermatrix and the right hand side.
    does not normalize the rows to ensure that the first non zero entry of non zero rows = 1.
    The difference with gauss_jordan_eliminationHMII is the fact that the row linear combination
    operations are performed in such a way as to not change the absolute value of the determinant.
    

    EXAMPLES:
 
    ::

        sage: [RefA, c] = gauss_jordan_eliminationHMIII(HM(2,2,'a'), HM(2,1,'b'))
        sage: RefA.printHM()
        [:, :]=
        [               a00                  0]
        [                 0 -a01*a10/a00 + a11]
        sage: c.printHM()
        [:, :]=
        [-a01*(a10*b00/a00 - b10)/(a01*a10/a00 - a11) + b00]
        [                                -a10*b00/a00 + b10]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    [A, b] = gaussian_eliminationHMIII(Cf,rs)
    # Initialization of the row and column index
    i=A.nrows()-1; j=0
    while i>0 or j>0:
        #if (A[i,:]).is_zero():
        if HM(1,A.n(1),[A[i,j0] for j0 in range(A.n(1))]).is_zero():
            # decrementing the row index and initializing the column index
            i=i-1; j=0
        else :
            while (A[i,j]).is_zero():
                # Incrementing the column index
                j = j + 1
            # performing row operations
            cf1=A[i,j]
            for r in range(i-1,-1,-1):
                #b[r,:] = -A[r,j]*b[i,:]+b[r,:]
                cf2=A[r,j]
                for j0 in range(b.n(1)):
                    b[r,j0]=-(cf2*cf1^(-1))*b[i,j0]+b[r,j0]
                #A[r,:] = -A[r,j]*A[i,:]+A[r,:]
                for j0 in range(A.n(1)):
                    A[r,j0]=-(cf2*cf1^(-1))*A[i,j0]+A[r,j0]
            i=i-1; j=0
    return [A,b]

def multiplicative_gaussian_elimination(Cf,rs,jndx=0):
    """
    Outputs the row echelon form of the input matrix and the right hand side.

    EXAMPLES:
 
    ::

        sage: [EfA,c,indx,Lst]=multiplicative_gaussian_elimination(Matrix(SR,HM(2,2,'a').listHM()), Matrix(SR,HM(2,1,'b').listHM()))
        sage: EfA
        [      1 a01/a00]
        [      0       1]
        sage: c
        [                                                                    (b00*e^(2*I*pi*k0))^(1/a00)]
        [1/((b10*e^(2*I*pi*k2)/((b00*e^(2*I*pi*k0))^(1/a00)*e^(2*I*pi*k1))^a10)^(1/(a01*a10/a00 - a11)))]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    A = copy(Cf); b = copy(rs)
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
                P =sum([Id[:,k]*Id[Integer(mod(k+1,Ta.nrows())),:] for k in range(Ta.nrows())])
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

        sage: [RefA, c, indx, Lst] = multiplicative_gauss_jordan_elimination(Matrix(SR,HM(2,2,'a').listHM()), Matrix(SR,HM(2,1,'b').listHM()))
        sage: RefA
        [1 0]
        [0 1]
        sage: c
        [(b00*e^(2*I*pi*k0))^(1/a00)/(e^(2*I*pi*k3)/(b10*e^(2*I*pi*k2)/((b00*e^(2*I*pi*k0))^(1/a00)*e^(2*I*pi*k1))^a10)^(1/(a01*a10/a00 - a11)))^(a01/a00)]
        [                                                  1/((b10*e^(2*I*pi*k2)/((b00*e^(2*I*pi*k0))^(1/a00)*e^(2*I*pi*k1))^a10)^(1/(a01*a10/a00 - a11)))]
        

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

def multiplicative_gaussian_eliminationII(Cf,rs,jndx=0):
    """
    Outputs the row echelon form of the input matrix and the right hand side.
    The solver here differs from the one above in the fact that it assumes
    that the entries of the Cf matrix are not symbolic and checks during
    the elimination steps whether or we are indeed adding new branches.

    EXAMPLES:
 
    ::

        sage: [EfA,c,indx,Lst]=multiplicative_gaussian_eliminationII(HM(2,2,'a').matrix(), HM(2,1,'b').matrix())
        sage: EfA
        [      1 a01/a00]
        [      0       1]
        sage: c
        [                                                                    (b00*e^(2*I*pi*k0))^(1/a00)]
        [1/((b10*e^(2*I*pi*k2)/((b00*e^(2*I*pi*k0))^(1/a00)*e^(2*I*pi*k1))^a10)^(1/(a01*a10/a00 - a11)))]


    AUTHORS:

    - Edinah K. Gnang
    - To Do: 
    """
    A = copy(Cf); b = copy(rs)
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
                P =sum([Id[:,k]*Id[Integer(mod(k+1,Ta.nrows())),:] for k in range(Ta.nrows())])
                Ta=P*Ta; Tb=P*Tb
                A[i:,:]=Ta
                b[i:,:]=Tb 
            # Performing the row operations.
            #if A[i,j]==-1 or A[i,j]==1:
            if (A[i,j].numerator()==1 or A[i,j].numerator()==-1) and A[i,j].denominator().is_integer():
                # In the case where no branching is introduced by log scaling the pivot
                for j0 in rg(b.ncols()):
                    b[i,j0]=b[i,j0]^(1/A[i,j])
            else:
                # In the case where some branching is introduced by log scaling the pivot
                for j0 in rg(b.ncols()):
                    b[i,j0]=(b[i,j0]*exp(I*2*pi*var('k'+str(indx))))^(1/A[i,j])
                indx = indx+1
                Lst.append(A[i,j])
            A[i,:]=(1/A[i,j])*A[i,:]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if A[r,:].is_zero():
                    r=r+1
                else:
                    #if A[r,j]==-1 or A[r,j]==1:
                    if A[r,j].is_integer():
                        # In the case where no branching is introduced by log scaling the pivot
                        for j0 in rg(b.ncols()):
                            b[r,j0]=b[i,j0]^(-A[r,j])*b[r,j0]
                    else:
                        # In the case where some branching is introduced by log scaling the pivot
                        for j0 in rg(b.ncols()):
                            b[r,j0]=(b[i,j0]*exp(I*2*pi*var('k'+str(indx))))^(-A[r,j])*b[r,j0]
                        if (A[r,j]).is_zero()==False:
                            indx = indx+1
                            Lst.append(1/A[r,j])
                    A[r,:]=-A[r,j]*A[i,:]+A[r,:]
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return [A, b, indx, Lst]

def multiplicative_gauss_jordan_eliminationII(Cf,rs,jndx=0):
    """
    Outputs the reduced row echelon form of the input matrix and the right hand side.
    The solver here differs from the one above in the fact that it assumes
    that the entries of the Cf matrix are not symbolic and checks during
    the elimination steps whether or we are indeed adding new branches.

    EXAMPLES:
 
    ::

        sage: [RefA, c, indx, L] = multiplicative_gauss_jordan_eliminationII(HM(2,2,'a').matrix(), HM(2,1,'b').matrix())
        sage: RefA
        [1 0]
        [0 1]
        sage: c
        [(b00*e^(2*I*pi*k0))^(1/a00)/(e^(2*I*pi*k3)/(b10*e^(2*I*pi*k2)/((b00*e^(2*I*pi*k0))^(1/a00)*e^(2*I*pi*k1))^a10)^(1/(a01*a10/a00 - a11)))^(a01/a00)]
        [                                                  1/((b10*e^(2*I*pi*k2)/((b00*e^(2*I*pi*k0))^(1/a00)*e^(2*I*pi*k1))^a10)^(1/(a01*a10/a00 - a11)))]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    [A, b, indx, Lst] = multiplicative_gaussian_eliminationII(Cf,rs,jndx)
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
                #if A[r,j]==-1 or A[r,j]==1:
                if A[r,j].is_integer():
                    # In the case where no branching is introduced by log scaling the pivot
                    for j0 in rg(b.ncols()):
                        b[r,j0]=b[i,j0]^(-A[r,j])*b[r,j0]
                else:
                    # In the case where some branching is introduced by log scaling the pivot
                    for j0 in rg(b.ncols()):
                        b[r,j0]=(b[i,j0]*exp(I*2*pi*var('k'+str(indx))))^(-A[r,j])*b[r,j0]
                    if (A[r,j]).is_zero()==False:
                        indx = indx+1
                        Lst.append(1/A[r,j])
                A[r,:] = -A[r,j]*A[i,:]+A[r,:]
            i = i - 1; j = 0
    return [A, b, indx, Lst]

def multiplicative_gaussian_eliminationHM(Cf,rs,jndx=0):
    """
    Outputs the row echelon form of the input matrix and the right hand side.
    The solver here differs from the one above in the fact that it assumes
    that the entries of the Cf HM are not symbolic and checks during
    the elimination steps whether or we are indeed adding new branches.

    EXAMPLES:
 
    ::

        sage: [EfA,c,indx,Lst]=multiplicative_gaussian_eliminationHM(HM(2,2,'a'), HM(2,1,'b'))
        sage: EfA.printHM()
        [:, :]=
        [      1 a01/a00]
        [      0       1]
        sage: c.printHM()
        [:, :]=
        [                                                                    (b00*e^(2*I*pi*k0))^(1/a00)]
        [1/((b10*e^(2*I*pi*k2)/((b00*e^(2*I*pi*k0))^(1/a00)*e^(2*I*pi*k1))^a10)^(1/(a01*a10/a00 - a11)))]


    AUTHORS:

    - Edinah K. Gnang
    - To Do: 
    """
    A = Cf.copy(); b = rs.copy()
    # Initialization of the row and column index
    i=0; j=0; indx=jndx; Lst = []
    #while i<A.nrows() and j<A.ncols():
    while i < A.n(0) and j < A.n(1):
        #while (A[i:,j]).is_zero() and j < A.ncols()-1:
        while A.slice(rg(i,A.n(0)),'row').slice([j],'col').is_zero() and j < A.n(1)-1:
            # Incrementing the column index
            j=j+1
        #if A[i:,:].is_zero()==False:
        if A.slice(rg(i,A.n(0)),'row').is_zero()==False:
            while A[i,j].is_zero():
                #Ta=A[i:,:]
                Ta=A.slice(rg(i,A.n(0)),'row')
                #Tb=b[i:,:]
                Tb=b.slice(rg(i,b.n(0)),'row')
                # Initializing the cyclic shift permutation matrix
                #Id=identity_matrix(Ta.nrows())
                Id=HM(2,Ta.n(0),'kronecker')
                #P=Matrix2HM(sum([Id.matrix()[:,k]*Id.matrix()[Integer(mod(k+1,Ta.nrows())),:] for k in rg(Ta.n(0))]))
                P=HM(2, [Integer(mod(fi-1, Ta.n(0))) for fi in rg(Ta.n(0))], 'perm')
                Ta=P*Ta; Tb=P*Tb
                #A[i:,:]=Ta
                for u in rg(i,A.n(0)):
                    for v in rg(A.n(1)):
                        A[u,v]=Ta[u-i,v]
                #b[i:,:]=Tb 
                for u in rg(i,b.n(0)):
                    for v in rg(b.n(1)):
                        b[u,v]=Tb[u-i,v]
            # Performing the row operations.
            if (A[i,j].numerator()==1 or A[i,j].numerator()==-1) and A[i,j].denominator().is_integer():
                for j0 in rg(b.ncols()):
                    b[i,j0]=b[i,j0]^(1/A[i,j])
            else:
                for j0 in rg(b.ncols()):
                    b[i,j0]=(b[i,j0]*exp(I*2*pi*var('k'+str(indx))))^(1/A[i,j])
                indx = indx+1
                Lst.append(A[i,j])
            #A[i,:]=(1/A[i,j])*A[i,:]
            tpv=A[i,j]
            for v in rg(A.n(1)):
                A[i,v]=(1/tpv)*A[i,v]
            for r in rg(i+1,A.nrows()):
                # Taking care of the zero row
                #if A[r,:].is_zero():
                if A.slice([r],'row').is_zero():
                    r=r+1
                else:
                    if A[r,j].is_integer():
                        for j0 in rg(b.ncols()):
                            b[r,j0]=b[i,j0]^(-A[r,j])*b[r,j0]
                    else:
                        for j0 in rg(b.ncols()):
                            b[r,j0]=(b[i,j0]*exp(I*2*pi*var('k'+str(indx))))^(-A[r,j])*b[r,j0]
                        if (A[r,j]).is_zero()==False:
                            indx = indx+1
                            Lst.append(1/A[r,j])
                    #A[r,:]=-A[r,j]*A[i,:]+A[r,:]
                    tpv=-A[r,j]
                    for v in rg(A.n(1)):
                        A[r,v]=tpv*A[i,v]+A[r,v]
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return [A, b, indx, Lst]

def multiplicative_gauss_jordan_eliminationHM(Cf,rs):
    """
    Outputs the reduced row echelon form of the input second order hypermatrix and the right hand side.
    does not normalize the rows to ensure that the first non zero entry of non zero rows = 1.
    The difference with previous implementations is the fact that the row linear combination
    operations are performed in such a way as to not change the absolute value of the determinant.
    

    EXAMPLES:
 
    ::

        sage: [A, b, indx, Lst] = multiplicative_gauss_jordan_eliminationHM(HM(2,2,'a'), HM(2,1,'b'))
        sage: A.printHM()
        [:, :]=
        [1 0]
        [0 1]
        sage: b.printHM()
        [:, :]=
        [(b00*e^(2*I*pi*k0))^(1/a00)/(e^(2*I*pi*k3)/(b10*e^(2*I*pi*k2)/((b00*e^(2*I*pi*k0))^(1/a00)*e^(2*I*pi*k1))^a10)^(1/(a01*a10/a00 - a11)))^(a01/a00)]
        [                                                  1/((b10*e^(2*I*pi*k2)/((b00*e^(2*I*pi*k0))^(1/a00)*e^(2*I*pi*k1))^a10)^(1/(a01*a10/a00 - a11)))]
        sage: indx
        4
        sage: Lst
        [a00, 1/a10, -a01*a10/a00 + a11, -a01/a00]



    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    [A, b, indx, Lst] = multiplicative_gaussian_eliminationHM(Cf,rs)
    # Initialization of the row and column index
    i=A.nrows()-1; j=0
    while i>0 or j>0:
        #if (A[i,:]).is_zero():
        if HM(1,A.n(1),[A[i,j0] for j0 in range(A.n(1))]).is_zero():
            # decrementing the row index and initializing the column index
            i=i-1; j=0
        else :
            while (A[i,j]).is_zero():
                # Incrementing the column index
                j = j + 1
            # performing row operations
            cf1=A[i,j]
            for r in range(i-1,-1,-1):
                #b[r,:] = -A[r,j]*b[i,:]+b[r,:]
                cf2=A[r,j]
                # Performing the row operations.
                if (cf2/cf1).is_integer():
                    for j0 in range(b.n(1)):
                        b[r,j0] = b[i,j0]^(-cf2/cf1)*b[r,j0]
                    for j0 in range(A.n(1)):
                        A[r,j0]=A[i,j0]*(-cf2/cf1)+A[r,j0]
                else:
                    for j0 in range(b.n(1)):
                        b[r,j0]=(b[i,j0]*exp(I*2*pi*var('k'+str(indx))))^(-cf2/cf1)*b[r,j0]
                    for j0 in range(A.n(1)):
                        A[r,j0]=A[i,j0]*(-cf2/cf1)+A[r,j0]
                    indx = indx+1
                    Lst.append(-cf2/cf1)
            i=i-1; j=0
    return [A, b, indx, Lst]

def multiplicative_gaussian_eliminationHMII(Cf, rs, jndx=0):
    """
    Outputs the row echelon form of the input second order hypermatrix and the right hand side.
    does not normalize the rows to ensure that the first non zero entry of non zero rows = 1.
    The difference with the previous implementation is the fact that the row linear combination
    operations are performed in such a way as to not change the absolute value of the determinant.
    This implementation is skew fields or division ring friendly by inputing hypermatrices
    whose entries are themselve hypermatrices of the approprioate size. 


    EXAMPLES:
 
    ::

        sage: [A, b, indx, Lst] = multiplicative_gaussian_eliminationHMII(HM(2,2,'a'), HM(2,1,'b'))
        sage: A.printHM()
        [:, :]=
        [               a00                a01]
        [                 0 -a01*a10/a00 + a11]
        sage: b.printHM()
        [:, :]=
        [                              b00]
        [b10/(b00*e^(2*I*pi*k0))^(a10/a00)]
        sage: indx
        1
        sage: Lst
        [-a10/a00]



    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing a copy of the input second order hypermatrices.
    A=Cf.copy(); b=rs.copy()
    # Initialization of the row and column index
    i=0; j=0; indx=jndx; Lst = []
    while i < A.n(0) and j < A.n(1):
        while HM(A.n(0)-i, 1, [A[i0,j] for i0 in range(i,A.n(0))]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        if HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))]).is_zero()==False:
            while A[i,j].is_zero(): 
                Ta=HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))])
                Tb=HM(b.n(0)-i, b.n(1), [b[i0,j0] for j0 in range(b.n(1)) for i0 in range(i,b.n(0))])
                # Initializing the cyclic shift permutation matrix
                Id=HM(2, Ta.n(0), 'kronecker')
                #P=sum([HM(Ta.n(0),1,[Id[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Id[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                P=HM(2, [Integer(mod(fi-1, Ta.n(0))) for fi in rg(Ta.n(0))], 'perm')
                Ta=P*Ta; Tb=P*Tb
                for i0 in range(Ta.n(0)):
                    for j0 in range(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
                for i0 in range(Tb.n(0)):
                    for j0 in range(Tb.n(1)):
                        b[i+i0,j0]=Tb[i0,j0]
            if A.n(0)-i-1 > 0 and not (HM(A.n(0)-i-1, 1, [A[i0,j] for i0 in range(i+1,A.n(0))]).is_zero() and j <= A.ncols()-1):
                # Performing the row operations.
                cf1=A[i,j]
                for r in range(i+1,A.nrows()):
                    # Taking care of the zero row
                    if HM(1,A.n(1),[A[r,j0] for j0 in range(A.n(1))]).is_zero():
                        r=r+1
                    else:
                        # Initialization of the coefficient
                        cf2=A[r,j]
                        # Performing the row operations.
                        if (cf2/cf1).is_integer():
                            for j0 in range(b.n(1)):
                                b[r,j0] = b[i,j0]^(-cf2/cf1)*b[r,j0]
                            for j0 in range(A.n(1)):
                                #A[r,j0] = A[i,j0]^(-cf2/cf1)*A[r,j0]
                                A[r,j0]=A[i,j0]*(-cf2/cf1)+A[r,j0]
                        else:
                            for j0 in range(b.n(1)):
                                b[r,j0]=(b[i,j0]*exp(I*2*pi*var('k'+str(indx))))^(-cf2/cf1)*b[r,j0]
                            for j0 in range(A.n(1)):
                                #A[r,j0]=A[i,j0]^(-cf2/cf1)*A[r,j0]
                                A[r,j0]=A[i,j0]*(-cf2/cf1)+A[r,j0]
                            indx = indx+1
                            Lst.append(-cf2/cf1)
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return [A, b, indx, Lst]

def multiplicative_gauss_jordan_eliminationHMII(Cf,rs):
    """
    Outputs the reduced row echelon form of the input second order hypermatrix and the right hand side.
    does not normalize the rows to ensure that the first non zero entry of non zero rows = 1.
    The difference with previous implementations is the fact that the row linear combination
    operations are performed in such a way as to not change the absolute value of the determinant.
    

    EXAMPLES:
 
    ::

        sage: [A, b, indx, Lst] = multiplicative_gauss_jordan_eliminationHMII(HM(2,2,'a'), HM(2,1,'b'))
        sage: A.printHM()
        [:, :]=
        [               a00                  0]
        [                 0 -a01*a10/a00 + a11]
        sage: b.printHM()
        [:, :]=
        [b00*(b10*e^(2*I*pi*k1)/(b00*e^(2*I*pi*k0))^(a10/a00))^(a01/(a01*a10/a00 - a11))]
        [                                              b10/(b00*e^(2*I*pi*k0))^(a10/a00)]
        sage: indx
        2
        sage: Lst
        [-a10/a00, a01/(a01*a10/a00 - a11)]



    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    [A, b, indx, Lst] = multiplicative_gaussian_eliminationHMII(Cf,rs)
    # Initialization of the row and column index
    i=A.nrows()-1; j=0
    while i>0 or j>0:
        #if (A[i,:]).is_zero():
        if HM(1,A.n(1),[A[i,j0] for j0 in range(A.n(1))]).is_zero():
            # decrementing the row index and initializing the column index
            i=i-1; j=0
        else :
            while (A[i,j]).is_zero():
                # Incrementing the column index
                j = j + 1
            # performing row operations
            cf1=A[i,j]
            for r in range(i-1,-1,-1):
                #b[r,:] = -A[r,j]*b[i,:]+b[r,:]
                cf2=A[r,j]
                # Performing the row operations.
                if (cf2/cf1).is_integer():
                    for j0 in range(b.n(1)):
                        b[r,j0] = b[i,j0]^(-cf2/cf1)*b[r,j0]
                    for j0 in range(A.n(1)):
                        A[r,j0]=A[i,j0]*(-cf2/cf1)+A[r,j0]
                else:
                    for j0 in range(b.n(1)):
                        b[r,j0]=(b[i,j0]*exp(I*2*pi*var('k'+str(indx))))^(-cf2/cf1)*b[r,j0]
                    for j0 in range(A.n(1)):
                        A[r,j0]=A[i,j0]*(-cf2/cf1)+A[r,j0]
                    indx = indx+1
                    Lst.append(-cf2/cf1)
            i=i-1; j=0
    return [A, b, indx, Lst]

def multiplicative_matrix_product(A,B):
    """
    Outputs the result of the multiplicative product of the
    two input matrices.

    EXAMPLES:
 
    ::

        sage: multiplicative_matrix_product(Matrix(SR,HM(2,2,'a').listHM()), Matrix(SR,HM(2,2,'b').listHM()))
        [b00^a00*b10^a01 b01^a00*b11^a01]
        [b00^a10*b10^a11 b01^a10*b11^a11]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    Rslt=Matrix(SR,zero_matrix(A.nrows(), B.ncols()))
    for i in range(A.nrows()):
        for k in range(B.ncols()):
            Rslt[i,k]=prod([B[j,k]^A[i,j] for j in range(A.ncols())])
    return Rslt

def multiplicative_matrix_productHM(A,B):
    """
    Outputs the result of the multiplicative product of the
    two input matrices.

    EXAMPLES:
 
    ::

        sage: multiplicative_matrix_productHM(HM(2,2,'a'), HM(2,2,'b')).printHM()
        [:, :]=
        [b00^a00*b10^a01 b01^a00*b11^a01]
        [b00^a10*b10^a11 b01^a10*b11^a11]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    Rslt=HM(A.nrows(), B.ncols(), 'zero')
    for i in range(A.nrows()):
        for k in range(B.ncols()):
            Rslt[i,k]=prod([B[j,k]^A[i,j] for j in range(A.ncols())])
    return Rslt

def vec_exp(vx,A):
    """
    Outputs the result of the multiplicative product of a
    matrix by a vector. Does not check that the left input
    actually is a vector.

    EXAMPLES:
 
    ::

        sage: vec_exp(HM(2,1,HM(2,'x').list()), HM(2,2,'a')).printHM()
        [:, :]=
        [x0^a00*x1^a01]
        [x0^a10*x1^a11]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    #return mprod(A,x)
    return GProdIII([A,vx],prod,BaseExp)

def linear_solver(A,b,x,v):
    """
    Outputs the Reduced Row Echelon Form of the input matrix and the right hand side.
    where A denotes the input matrix, b denotes the right-hand side vector, x denotes
    the variable vector coming from the original system of equations, and v denotes 
    the free variable vector.

    EXAMPLES:
 
    ::

        sage: sz=2; Eq=[var('x'+str(i))+var('x'+str(sz+j))==var('a'+str(i)+str(j)) for i in range(sz) for j in range(sz)]
        sage: [A,b]=ConstraintFormatorII(Eq,[var('x'+str(i)) for i in range(2*sz)])
        sage: Mx=Matrix(SR,A.ncols(),1,[var('x'+str(i)) for i in range(A.ncols())])
        sage: Mv=Matrix(SR,A.ncols(),1,[var('t'+str(i)) for i in range(A.ncols())])
        sage: linear_solver(A,b,Mx,Mv)
        [x0 == a00 - a10 + a11 - t3,
         x1 == a11 - t3,
         x2 == a10 - a11 + t3,
         0  == -a00 + a01 + a10 - a11]

    AUTHORS:
    - Initial implementation by Edinah K. Gnang updates to the doc string by Jeanine S. Gnang
    - To Do: 
    """
    # Initialization of the reduced echelon form.
    [Ap,bp]=gauss_jordan_elimination(A,b)
    Id1 = identity_matrix(Ap.nrows())
    Id2 = identity_matrix(Ap.ncols())
    # Obtainin the list of pivot variables.
    Pm=Matrix(SR,zero_matrix(Ap.nrows(),Ap.ncols()))
    for i in range(Ap.nrows()):
        if not Ap[i,:].is_zero():
            for j in range(Ap.ncols()):
                if Ap[i,j]==1:
                    break
            Pm=Pm+Id1[:,i]*Id2[j,:]
    # Expressing the solutions
    tp1=Pm*x; tp2=bp-(Ap-Pm)*v
    return [tp1[i,0]==tp2[i,0] for i in range(tp1.nrows())]

def linear_solverHM(A,b,x,v):
    """
    Outputs the Reduced Row Echelon Form of the input matrix and the right hand side.
    where A denotes the input matrix, b denotes the right-hand side vector, x denotes
    the variable vector coming from the original system of equations, and v denotes 
    the free variable vector.

    EXAMPLES:
 
    ::

        sage: sz=2; Eq=[var('x'+str(i))+var('x'+str(sz+j))==var('a'+str(i)+str(j)) for i in range(sz) for j in range(sz)]
        sage: [A,b]=ConstraintFormatorHM(Eq,[var('x'+str(i)) for i in range(2*sz)])
        sage: linear_solverHM(A,b,HM(A.n(1),1,var_list('x',A.n(1))),HM(A.n(1),1,var_list('t',A.n(1))))
        [x0 == a00 - a10 + a11 - t3,
         x1 == a11 - t3,
         x2 == a10 - a11 + t3,
         0  == -a00 + a01 + a10 - a11]

    AUTHORS:
    - Initial implementation by Edinah K. Gnang updates to the doc string by Jeanine S. Gnang
    - To Do: 
    """
    # Initialization of the reduced echelon form.
    [Ap,bp]=gauss_jordan_eliminationHM(A,b)
    Id1 = HM(2, Ap.n(0), 'kronecker')
    Id2 = HM(2, Ap.n(1), 'kronecker')
    # Obtainin the list of pivot variables.
    Pm=HM(Ap.n(0), Ap.n(1), 'zero')
    for i in range(Ap.n(0)):
        if not HM(1, Ap.n(1), [SR(Ap[i,u]) for u in range(Ap.n(1))]).is_zero():
            for j in range(Ap.n(1)):
                if Ap[i,j]==1:
                    break
            Pm=Pm+HM(Id1.n(0), 1, [SR(Id1[s,i]) for s in range(Id1.n(0))])*HM(1,Id2.n(1),[SR(Id2[j,t]) for t in range(Id2.n(1))])
    # Expressing the solutions
    tp1=Pm*x; tp2=bp-(Ap-Pm)*v
    return [tp1[i,0]==tp2[i,0] for i in range(tp1.n(0))]

def default_linear_solver(EqL, Lv, Lf):
    """
    Formats the constraints and outputs the solution of the system of linear constraints.
    The input EqL is the list of constraints. The input Lv corresponds to the list of
    variable constraints appearing in EqL. The input Lf corresponds to the free variable
    name each of which is put in correspondent with the entries of Lv.


    EXAMPLES:
 
    ::

        sage: sz=2; EqL=[var('x'+str(i))+var('x'+str(sz+j))==var('a'+str(i)+str(j)) for i in range(sz) for j in range(sz)]
        sage: default_linear_solver(EqL, var_list('x', 2*sz), var_list('t', 2*sz))
        [x0 == a00 - a10 + a11 - t3,
         x1 == a11 - t3,
         x2 == a10 - a11 + t3,
         0  == -a00 + a01 + a10 - a11]

    AUTHORS:
    - Initial implementation by Edinah K. Gnang updates to the doc string by Jeanine S. Gnang
    - To Do: 
    """
    # Formating the constraints
    [A,b]=ConstraintFormatorHM(EqL,Lv)
    # Initialization of the variable vectors
    x=HM(A.ncols(), 1, Lv); v=HM(A.ncols(), 1, Lf)
    # Initialization of the reduced echelon form.
    [Ap,bp]=gauss_jordan_eliminationHM(A,b)
    Id1 = HM(2, Ap.n(0), 'kronecker')
    Id2 = HM(2, Ap.n(1), 'kronecker')
    # Obtainin the list of pivot variables.
    Pm=HM(Ap.n(0), Ap.n(1), 'zero')
    for i in range(Ap.n(0)):
        if not HM(1, Ap.n(1), [SR(Ap[i,u]) for u in range(Ap.n(1))]).is_zero():
            for j in range(Ap.n(1)):
                if Ap[i,j]==1:
                    break
            Pm=Pm+HM(Id1.n(0), 1, [SR(Id1[s,i]) for s in range(Id1.n(0))])*HM(1,Id2.n(1),[SR(Id2[j,t]) for t in range(Id2.n(1))])
    # Expressing the solutions
    tp1=Pm*x; tp2=bp-(Ap-Pm)*v
    return [tp1[i,0]==tp2[i,0] for i in range(tp1.n(0))]

def multiplicative_linear_solver(A,b,x,v):
    """
    Outputs the solution to a multiplicative linear system of equations.

    EXAMPLES:
 
    ::

        sage: sz=2; Eq=[var('x'+str(i))*var('x'+str(sz+j))==var('a'+str(i)+str(j)) for i in range(sz) for j in range(sz)]
        sage: [A,b]=multiplicativeConstraintFormator(Eq,[var('x'+str(i)) for i in range(2*sz)])
        sage: Mx=Matrix(SR,A.ncols(),1,[var('x'+str(i)) for i in range(A.ncols())])
        sage: Mv=Matrix(SR,A.ncols(),1,[var('t'+str(i)) for i in range(A.ncols())])
        sage: multiplicative_linear_solver(A,b,Mx,Mv)
        [x0 == a00*a11/(a10*t3),
         x1 == a11/t3,
         x2 == a10*t3/a11,
         1 == a01*a10/(a00*a11)]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the reduced echelon form.
    [Ap,bp]=multiplicative_gauss_jordan_eliminationII(A,b)[:2]
    Id1=identity_matrix(Ap.nrows())
    Id2=identity_matrix(Ap.ncols())
    # Obtainin the list of pivot variables.
    Pm=Matrix(SR,zero_matrix(Ap.nrows(),Ap.ncols()))
    for i in range(Ap.nrows()):
        if not Ap[i,:].is_zero():
            for j in range(Ap.ncols()):
                if Ap[i,j]==1:
                    break
            Pm=Pm+Id1[:,i]*Id2[j,:]
    # Expressing the solutions
    tp1=multiplicative_matrix_product(Pm,x)
    tp2=multiplicative_matrix_product((Ap-Pm),v)
    return [tp1[i,0]==bp[i,0]/tp2[i,0] for i in range(tp1.nrows())]

def multiplicative_linear_solverHM(A,b,x,v):
    """
    Outputs the solution to a multiplicative linear system of equations.

    EXAMPLES:
 
    ::

        sage: sz=2; Eq=[var('x'+str(i))*var('x'+str(sz+j))==var('a'+str(i)+str(j)) for i in range(sz) for j in range(sz)]
        sage: [A,b]=multiplicativeConstraintFormatorHM(Eq,var_list('x',2*sz))
        sage: Mx=HM(A.ncols(),1,var_list('x',A.ncols()))
        sage: Mv=HM(A.ncols(),1,var_list('t',A.ncols()))
        sage: multiplicative_linear_solverHM(A,b,Mx,Mv)
        [x0 == a00*a11/(a10*t3),
         x1 == a11/t3,
         x2 == a10*t3/a11,
         1 == a01*a10/(a00*a11)]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the reduced echelon form.
    [Ap,bp]=multiplicative_gauss_jordan_eliminationII(A.matrix(),b.matrix())[:2]
    Id1=identity_matrix(Ap.nrows())
    Id2=identity_matrix(Ap.ncols())
    # Obtainin the list of pivot variables.
    Pm=Matrix(SR,zero_matrix(Ap.nrows(),Ap.ncols()))
    for i in range(Ap.nrows()):
        if not Ap[i,:].is_zero():
            for j in range(Ap.ncols()):
                if Ap[i,j]==1:
                    break
            Pm=Pm+Id1[:,i]*Id2[j,:]
    # Expressing the solutions
    tp1=multiplicative_matrix_product(Pm,x)
    tp2=multiplicative_matrix_product((Ap-Pm),v)
    return [tp1[i,0]==bp[i,0]/tp2[i,0] for i in range(tp1.nrows())]

def default_multiplicative_linear_solver(EqL, Lv, Lf):
    """
    Formats the constraints performs and solves the multiplicatively linear constraints.
    Outputs the solutions. The input EqL corresponds to a list of constraints. The input Lv
    corresponds to the list of variables appearing in the constraints. The input Lf corresponds
    to the list of free varaibles each taken in correspondence with the entries of Lv. This 
    implementation tacitly assumes that the  the input constraints are indeed multiplicatively linear.


    EXAMPLES:
 
    ::

        sage: sz=2; EqL=[var('x'+str(i))*var('x'+str(sz+j))==var('a'+str(i)+str(j)) for i in range(sz) for j in range(sz)]
        sage: default_multiplicative_linear_solver(EqL, var_list('x', 2*sz), var_list('t', 2*sz))
        [x0 == a00*a11/(a10*t3),
         x1 == a11/t3,
         x2 == a10*t3/a11,
         1 == a01*a10/(a00*a11)]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Formatting the constraints
    [A,b]=multiplicativeConstraintFormatorHM(EqL, Lv)
    # Initialization of the variables
    x=HM(A.ncols(), 1, Lv); v=HM(A.ncols(), 1, Lf)
    # Initialization of the reduced echelon form.
    [Ap,bp]=multiplicative_gauss_jordan_eliminationII(A.matrix(),b.matrix())[:2]
    #Id1=identity_matrix(Ap.nrows())
    Id1=HM(2, Ap.nrows(), 'kronecker')
    #Id2=identity_matrix(Ap.ncols())
    Id2=HM(2, Ap.ncols(), 'kronecker')
    # Obtainin the list of pivot variables.
    #Pm=Matrix(SR,zero_matrix(Ap.nrows(),Ap.ncols()))
    Pm=HM(Ap.nrows(), Ap.ncols(), 'zero')
    for i in range(Ap.nrows()):
        if not Ap[i,:].is_zero():
            for j in range(Ap.ncols()):
                if Ap[i,j]==1:
                    break
            Pm=Pm+HM(Ap.nrows(),1,[Id1[f,i] for f in rg(Ap.nrows())])*HM(1,Ap.ncols(),[Id2[j,g] for g in rg(Ap.ncols())])
    # Expressing the solutions
    tp1=x^Pm; tp2=v^(Matrix2HM(Ap)-Pm)
    return [tp1[i,0]==bp[i,0]/tp2[i,0] for i in range(tp1.nrows())]

def multiplicative_least_square_linear_solver(A,b,x,v):
    """
    Outputs the solution to the multiplicative least square problem

    EXAMPLES:
 
    ::

        sage: sz=2; Eq=[var('x'+str(i))+var('x'+str(sz+j))==var('a'+str(i)+str(j)) for i in range(sz) for j in range(sz)]
        sage: [A,b]=ConstraintFormatorII(Eq,[var('x'+str(i)) for i in range(2*sz)])
        sage: Mx=Matrix(SR,A.ncols(),1,[var('x'+str(i)) for i in range(A.ncols())])
        sage: Mv=Matrix(SR,A.ncols(),1,[var('t'+str(i)) for i in range(A.ncols())])
        sage: multiplicative_least_square_linear_solver(A,b,Mx,Mv)
        [x0 == sqrt(a00*a01)*e^(3/2*I*pi*k0 + 1/2*I*pi*k1 - I*pi*k3)/(sqrt(a00*a10/(sqrt(a00*a01)*sqrt(a10*a11)))*t3),
         x1 == sqrt(a10*a11)*e^(1/2*I*pi*k0 + 3/2*I*pi*k1 - I*pi*k2)/(sqrt(a00*a10/(sqrt(a00*a01)*sqrt(a10*a11)))*t3),
         x2 == a00*a10*t3*e^(-I*pi*k0 - I*pi*k1)/(sqrt(a00*a01)*sqrt(a10*a11)),
         1 == e^(-2*I*pi*k0 - 2*I*pi*k1)]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the reduced echelon form.
    [Ap,bp]=multiplicative_gauss_jordan_eliminationII(A.transpose()*A, multiplicative_matrix_product(A.transpose(),b))[:2]
    Id1=identity_matrix(Ap.nrows())
    Id2=identity_matrix(Ap.ncols())
    # Obtainin the list of pivot variables.
    Pm=Matrix(SR,zero_matrix(Ap.nrows(),Ap.ncols()))
    for i in range(Ap.nrows()):
        if not Ap[i,:].is_zero():
            for j in range(Ap.ncols()):
                if Ap[i,j]==1:
                    break
            Pm=Pm+Id1[:,i]*Id2[j,:]
    # Expressing the solutions
    tp1=multiplicative_matrix_product(Pm,x)
    tp2=multiplicative_matrix_product((Ap-Pm),v)
    return [tp1[i,0]==bp[i,0]/tp2[i,0] for i in range(tp1.nrows())]

def multiplicative_least_square_linear_solverHM(A,b,x,v):
    """
    Outputs the solution to the multiplicative least square problem

    EXAMPLES:
 
    ::

        sage: sz=2; Eq=[var('x'+str(i))*var('x'+str(sz+j))==var('a'+str(i)+str(j)) for i in range(sz) for j in range(sz)]
        sage: [A,b]=multiplicativeConstraintFormatorHM(Eq,[var('x'+str(i)) for i in range(2*sz)])
        sage: Mx=HM(A.ncols(),1,[var('x'+str(i)) for i in range(A.ncols())])
        sage: Mv=HM(A.ncols(),1,[var('t'+str(i)) for i in range(A.ncols())])
        sage: multiplicative_least_square_linear_solverHM(A,b,Mx,Mv)
        [x0 == sqrt(a00*a01)*e^(3/2*I*pi*k0 + 1/2*I*pi*k1 - I*pi*k3)/(sqrt(a00*a10/(sqrt(a00*a01)*sqrt(a10*a11)))*t3),
         x1 == sqrt(a10*a11)*e^(1/2*I*pi*k0 + 3/2*I*pi*k1 - I*pi*k2)/(sqrt(a00*a10/(sqrt(a00*a01)*sqrt(a10*a11)))*t3),
         x2 == a00*a10*t3*e^(-I*pi*k0 - I*pi*k1)/(sqrt(a00*a01)*sqrt(a10*a11)),
         1 == e^(-2*I*pi*k0 - 2*I*pi*k1)]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the reduced echelon form.
    [Ap,bp]=multiplicative_gauss_jordan_eliminationII((A.transpose()*A).matrix(), multiplicative_matrix_productHM(A.transpose(),b).matrix())[:2]
    Id1=identity_matrix(Ap.nrows())
    Id2=identity_matrix(Ap.ncols())
    # Obtainin the list of pivot variables.
    Pm=Matrix(SR,zero_matrix(Ap.nrows(),Ap.ncols()))
    for i in range(Ap.nrows()):
        if not Ap[i,:].is_zero():
            for j in range(Ap.ncols()):
                if Ap[i,j]==1:
                    break
            Pm=Pm+Id1[:,i]*Id2[j,:]
    # Expressing the solutions
    tp1=multiplicative_matrix_product(Pm,x)
    tp2=multiplicative_matrix_product((Ap-Pm),v)
    return [tp1[i,0]==bp[i,0]/tp2[i,0] for i in range(tp1.nrows())]

def default_multiplicative_linear_solverHM(EqL, Lv, Lf):
    """
    Formats the constraints performs and solves the multiplicatively linear constraints.
    Outputs the solutions. The input EqL corresponds to a list of constraints. The input Lv
    corresponds to the list of variables appearing in the constraints. The input Lf corresponds
    to the list of free varaibles each taken in correspondence with the entries of Lv. This 
    implementation tacitly assumes that the  the input constraints are indeed multiplicatively linear.


    EXAMPLES:
 
    ::

        sage: sz=2; EqL=[GProdIII([HM(sz,sz,'a'), HM(sz,1,var_list('x',sz))], prod, BaseExp)[i,0]==var_list('b',sz)[i] for i in rg(sz)]
        sage: default_multiplicative_linear_solverHM(EqL, var_list('x', sz), var_list('t', sz))
        [[x0 == (b0*(b1*e^(2*I*pi*k1)/(b0*e^(2*I*pi*k0))^(a10/a00))^(a01/(a01*a10/a00 - a11))*e^(2*I*pi*k2))^(1/a00),
          x1 == (b1*e^(2*I*pi*k3)/(b0*e^(2*I*pi*k0))^(a10/a00))^(1/a11)],
         4,
         [-a10/a00, a01/(a01*a10/a00 - a11), 1/a00, -1/(a01*a10/a00 - a11)]]
        sage: sz=2; Eq=[var('x'+str(i))*var('x'+str(sz+j))==var('a'+str(i)+str(j)) for i in range(sz) for j in range(sz)]
        sage: sz=2; default_multiplicative_linear_solverHM(Eq, var_list('x',2*sz), var_list('t',2*sz))
        [[x0 == a01, x1 == a11, x2 == a00/a01, 1 == a01*a10/(a00*a11)], 0, []]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Formatting the constraints
    [A, b]=multiplicativeConstraintFormatorIIIHM(EqL, Lv)
    # Initialization of the variables
    Hx = HM(A.ncols(), 1, Lv); Hv = HM(A.ncols(), 1, Lf)
    # Initialization of the reduced echelon form.
    [Ap, bp, indx, Lst] = multiplicative_gauss_jordan_eliminationHMII(A, b)
    Id1=HM(2, Ap.nrows(), 'kronecker'); Id2=HM(2, Ap.ncols(), 'kronecker')
    # Obtainin the list of pivot variables.
    Pm=HM(Ap.nrows(), Ap.ncols(), 'zero'); Qm=HM(Ap.nrows(), Ap.ncols(), 'zero')
    for i in range(Ap.nrows()):
        if not Ap.slice([i],'row').is_zero():
            for j in range(Ap.ncols()):
                if Ap[i,j] != 0:
                    break
            Pm=Pm+HM(Ap.nrows(),1,[Id1[f,i] for f in rg(Ap.nrows())])*HM(1,Ap.ncols(),[Ap[j,g] for g in rg(Ap.ncols())])
            Qm=Qm+HM(Ap.nrows(),1,[Id1[f,i] for f in rg(Ap.nrows())])*HM(1,Ap.ncols(),[Id2[j,g] for g in rg(Ap.ncols())])
    # Expressing the solutions
    tp1=GProdIII([Pm, Hx], prod, BaseExp); tp2=GProdIII([Ap-Pm, Hv], prod, BaseExp); tq=GProdIII([Qm, Hx], prod, BaseExp)
    # Initialization of the variables
    [Af, bf]=multiplicativeConstraintFormatorIIIHM([tp1[gi,0] == bp[gi,0]/tp2[gi,0] for gi in range(tp1.nrows())], [vb for vb in tq.list() if vb in Lv])
    for r in rg(min(Af.dimensions())):
        if (1/Af[r,r]).is_integer():
            for j0 in rg(bf.n(1)):
                bf[r,j0]=(bf[r,j0])^(1/Af[r,r])
        else:
            for j0 in rg(bf.n(1)):
                bf[r,j0]=(bf[r,j0]*exp(I*2*pi*var('k'+str(indx))))^(1/A[r,r])
            indx = indx+1
            Lst.append(1/Af[r,r])
    return [[tq[i,0] == bf[i,0] for i in range(tq.nrows())], indx, Lst]

def exponential_linear_solverHM(Ha,b,x,v):
    """
    Outputs the solution to a multiplicative linear system of equations.

    EXAMPLES:
 
    ::

        sage: sz=2; X=var_list('x',sz); A=HM(sz,sz,'a'); Eq=[(A[0,0]^X[0])*(A[0,1]^X[1])==7, (A[1,0]^X[0])*(A[1,1]^X[1])==2]
        sage: [A,b]=exponentialConstraintFormatorHM(Eq,X); Mx=HM(sz,1,X); Mv=HM(sz,1,var_list('t',sz))
        sage: exponential_linear_solverHM(A,b,Mx,Mv)
        [e^x0 == (7*e^(2*I*pi*k0))^(1/log(a00))/(e^(2*I*pi*k3)/(2*e^(2*I*pi*k2)/((7*e^(2*I*pi*k0))^(1/log(a00))*e^(2*I*pi*k1))^log(a10))^(1/(log(a01)*log(a10)/log(a00) - log(a11))))^(log(a01)/log(a00)),
         e^x1 == (1/((2*e^(2*I*pi*k2)/((7*e^(2*I*pi*k0))^(1/log(a00))*e^(2*I*pi*k1))^log(a10))^(1/(log(a01)*log(a10)/log(a00) - log(a11)))))]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    phi = lambda x: ln(x)
    #A = apply(HM, Ha.dimensions()+[ Ha.apply_map(phi).list() ])
    A = HM(*(Ha.dimensions()+[ Ha.apply_map(phi).list() ]))
    # Initialization of the reduced echelon form.
    [Ap,bp] = multiplicative_gauss_jordan_eliminationII(A.matrix(),b.matrix())[:2]
    Id1 = identity_matrix(Ap.nrows()); Id2 = identity_matrix(Ap.ncols())
    # Obtainin the list of pivot variables.
    Pm = Matrix(SR,zero_matrix(Ap.nrows(),Ap.ncols()))
    for i in range(Ap.nrows()):
        if not Ap[i,:].is_zero():
            for j in range(Ap.ncols()):
                if Ap[i,j] == 1:
                    break
            Pm = Pm+Id1[:,i]*Id2[j,:]
    # Expressing the solutions
    tp1 = multiplicative_matrix_product(Pm,x)
    tp2 = multiplicative_matrix_product((Ap-Pm),v)
    return [exp(tp1[i,0]) == bp[i,0]/tp2[i,0] for i in range(tp1.nrows())]

def SecondOrderHadamardFactorization(A,B):
    """
    Outputs the matrix factorization induced by hadamard matrices.

    EXAMPLES:
 
    ::

        sage: [U,V]=SecondOrderHadamardFactorization(HM(2,2,'a'),HM(2,2,'b'))
        sage: Prod(U,V).simplify_full()
        [[a00*b00 + a01*b10, a00*b01 + a01*b11], [a10*b00 + a11*b10, a10*b01 + a11*b11]]


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
        raise ValueError("The input hypermpatrix are of inapropriate sizes")

def RealRow_Gram_Schmidt(M):
    """
    Implements the naive Gram-Schmidt algorithm for rows.

    EXAMPLES:

    ::

        sage: Q=RealRow_Gram_Schmidt(Matrix(SR,HM(2,2,'a').listHM())); Matrix(SR,[[(Q*Q.transpose())[i,j].simplify_full() for i in range(2)] for j in range(2)])
        [                                                  a00^2 + a01^2                                                               0]
        [                                                              0 (a01^2*a10^2 - 2*a00*a01*a10*a11 + a00^2*a11^2)/(a00^2 + a01^2)]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the resulting matrix
    Rs = Matrix(SR, zero_matrix(M.nrows(),M.ncols()))
    # Initializing the first vector
    # the implementation assumes that the first vector is non zero
    Rs[0,:]=M[0,:]
    for i in range(1,M.nrows()):
        v = M[i,:]
        v = v-sum([(v*Rs[j,:].transpose())[0,0]/sum(Rs[j,s]^2 for s in range(Rs.ncols()))*Rs[j,:] for j in range(i) if not sum(Rs[j,s]^2 for s in range(Rs.ncols())).is_zero()])
        Rs[i,:] = v
    return Rs

def RealRow_Gram_SchmidtHM(Hm):
    """
    Implements the naive Gram-Schmidt algorithm for rows.
    The implementation here is a naive varian of the 
    implementation above.

    EXAMPLES:

    ::

        sage: Q=RealRow_Gram_SchmidtHM(HM(2,2,'a')); (Q*Q.transpose()).simplify_full().printHM()
        [:, :]=
        [                                                  a00^2 + a01^2                                                               0]
        [                                                              0 (a01^2*a10^2 - 2*a00*a01*a10*a11 + a00^2*a11^2)/(a00^2 + a01^2)]
        sage: A = HM([[-2*I + 3, -2], [5*I - 1, -I + 2]]); od=2; sz = 2
        sage: B=HM(2*sz, 2*sz,'zero') # Conversion from complex to real
        sage: for i in range(sz):
        ....:     for j in range(sz):
        ....:         Tmp=HM(sz, sz,'zero'); Tmp[i,j]=1
        ....:         B=B+Tmp.tensor_product(HM([[A[i,j].real(),-A[i,j].imag()],[A[i,j].imag(),A[i,j].real()]]))
        ....:
        sage: Q=RealRow_Gram_SchmidtHM(B) # Performing the orthogonalization
        sage: U=HM(sz, sz,'zero') # Conversion from real back to complex
        sage: for i in range(0, 2*sz, 2):
        ....:     for j in range(0, 2*sz, 2):
        ....:         U[i/2,j/2]=Q[i,j]+I*Q[i+1,j]
        ....:
        sage: U.printHM()
        [:, :]=
        [     -2*I + 3            -2]
        [6/17*I + 4/17       13/17*I]
        sage: U*U.conjugate_transpose()
        [[17, 0], [0, 13/17]]
        sage: A = HM([[-2*I + 3, -2], [5*I - 1, -I + 2]]); U = RR2CC_deflate(RealRow_Gram_SchmidtHM(CC2RR_inflate(A))); U
        [[-2*I + 3, -2], [6/17*I + 4/17, 13/17*I]]
        sage: U*U.conjugate_transpose()
        [[17, 0], [0, 13/17]]


    AUTHORS:
    - Edinah K. Gnang
    """
    M=Hm.matrix()
    # Initialization of the resulting matrix
    Rs = Matrix(SR, zero_matrix(M.nrows(),M.ncols()))
    # Initializing the first vector
    # the implementation assumes that the first vector is non zero
    Rs[0,:]=M[0,:]
    for i in range(1,M.nrows()):
        v = M[i,:]
        v = v-sum([(v*Rs[j,:].transpose())[0,0]/sum(Rs[j,s]^2 for s in range(Rs.ncols()))*Rs[j,:] for j in range(i) if not sum(Rs[j,s]^2 for s in range(Rs.ncols())).is_zero()])
        Rs[i,:] = v
    return HM(M.nrows(),M.ncols(),Rs.transpose().list())

def RealColumn_Gram_Schmidt(M):
    """
    Implements the naive Gram-Schmidt algorithm for the columns.

    EXAMPLES:

    ::

        sage: Q=RealColumn_Gram_Schmidt(Matrix(SR,HM(2,2,'a').listHM())); Matrix(SR,[[(Q.transpose()*Q)[i,j].simplify_full() for i in range(2)] for j in range(2)])
        [                                                  a00^2 + a10^2                                                               0]
        [                                                              0 (a01^2*a10^2 - 2*a00*a01*a10*a11 + a00^2*a11^2)/(a00^2 + a10^2)]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the symbolic matrix to be used as output
    Rs = Matrix(SR, zero_matrix(M.nrows(), M.ncols()))
    # This implementation assumes that the first column is a non zero column.
    # Initializing the first vector to a unit vector.
    Rs[:,0]=M[:,0]
    # Loop removing the bad components for all the remaining vectors
    for i in range(1, M.ncols()):
        v = M[:,i]
        v = v-sum([(v.transpose()*Rs[:,j])[0,0]/sum(Rs[s,j]^2 for s in range(Rs.nrows()))*Rs[:,j] for j in range(i) if not sum(Rs[s,j]^2 for s in range(Rs.nrows())).is_zero()])
        Rs[:,i] = v
    return Rs

def RealColumn_Gram_SchmidtHM(Hm):
    """
    Implements the naive Gram-Schmidt algorithm for the columns.
    This implementation is a naive variant of the implementation
    above.

    EXAMPLES:

    ::

        sage: Q=RealColumn_Gram_SchmidtHM(HM(2,2,'a')); (Q.transpose()*Q).simplify_full().printHM()
        [:, :]=
        [                                                  a00^2 + a10^2                                                               0]
        [                                                              0 (a01^2*a10^2 - 2*a00*a01*a10*a11 + a00^2*a11^2)/(a00^2 + a10^2)]
        sage: A = HM([[-2*I + 3, -2], [5*I - 1, -I + 2]]); od=2; sz = 2
        sage: B=HM(2*sz, 2*sz,'zero') # Conversion from complex to real
        sage: for i in range(sz):
        ....:     for j in range(sz):
        ....:         Tmp=HM(sz, sz,'zero'); Tmp[i,j]=1
        ....:         B=B+Tmp.tensor_product(HM([[A[i,j].real(),-A[i,j].imag()],[A[i,j].imag(),A[i,j].real()]]))
        ....:
        sage: Q=RealColumn_Gram_SchmidtHM(B) # Performing the orthogonalization
        sage: U=HM(sz, sz,'zero') # Conversion from real back to complex
        sage: for i in range(0, 2*sz, 2):
        ....:     for j in range(0, 2*sz, 2):
        ....:         U[i/2,j/2]=Q[i,j]+I*Q[i+1,j]
        ....:
        sage: U.printHM()
        [:, :]=
        [   -2*I + 3 1/3*I - 1/3]
        [    5*I - 1       1/3*I]
        sage: U.conjugate_transpose()*U
        [[39, 0], [0, 1/3]]


    AUTHORS:
    - Edinah K. Gnang
    """
    M=Hm.matrix()
    # Initialization of the symbolic matrix to be used as output
    Rs = Matrix(SR, zero_matrix(M.nrows(), M.ncols()))
    # This implementation assumes that the first column is a non zero column.
    # Initializing the first vector to a unit vector.
    Rs[:,0]=M[:,0]
    # Loop removing the bad components for all the remaining vectors
    for i in range(1, M.ncols()):
        v = M[:,i]
        v = v-sum([(v.transpose()*Rs[:,j])[0,0]/sum(Rs[s,j]^2 for s in range(Rs.nrows()))*Rs[:,j] for j in range(i) if not sum(Rs[s,j]^2 for s in range(Rs.nrows())).is_zero()])
        Rs[:,i] = v
    return HM(M.nrows(), M.ncols(), Rs.transpose().list())

def All_RealRow_Gram_Schmidt(M):
    """
    Takes a square matrix as input and returns an orthognal martix with minimal row distance.

    EXAMPLES:

    ::

        sage: All_RealRow_Gram_Schmidt(Matrix(SR,HM(2,2,'a').listHM()))
        [
        [                                           a00                                            a01]  [a00 - (a00*a10 + a01*a11)*a10/(a10^2 + a11^2) a01 - (a00*a10 + a01*a11)*a11/(a10^2 + a11^2)]
        [-(a00*a10 + a01*a11)*a00/(a00^2 + a01^2) + a10 -(a00*a10 + a01*a11)*a01/(a00^2 + a01^2) + a11], [                                          a10                                           a11]
        ]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initilization of the list of Permutations
    P = Permutations(range(M.nrows()))
    # Initialization of the list
    L=[]
    # Initialization of the final matrix
    Q = Matrix(SR, zero_matrix(M.nrows(), M.ncols()))
    for p in P:
        # Initialization of the temporary matrix.
        Rs = Matrix(SR, zero_matrix(M.nrows(), M.ncols()))
        # Initialization of the first vector.
        # the implementation assumes that the first vector is non zero
        Rs[p[0],:] = M[p[0],:]
        for i in range(1,M.nrows()):
            v = M[p[i],:]
            v = v-sum([(v*Rs[p[j],:].transpose())[0,0]/sum(Rs[p[j],s]^2 for s in range(Rs.ncols()))*Rs[p[j],:] for j in range(i) if not sum(Rs[p[j],s]^2 for s in range(Rs.ncols())).is_zero()])
            Rs[p[i],:] = v
        L.append(copy(Rs))
    return L

def All_RealColumn_Gram_Schmidt(M):
    """
    Takes a square matrix as input and returns an orthognal martix with minimal row distance.

    EXAMPLES:

    ::

        sage: All_RealColumn_Gram_Schmidt(Matrix(SR,HM(2,2,'a').listHM()))
        [
        [                                           a00 -(a00*a01 + a10*a11)*a00/(a00^2 + a10^2) + a01]  [a00 - (a00*a01 + a10*a11)*a01/(a01^2 + a11^2)                                           a01]
        [                                           a10 -(a00*a01 + a10*a11)*a10/(a00^2 + a10^2) + a11], [a10 - (a00*a01 + a10*a11)*a11/(a01^2 + a11^2)                                           a11]
        ]

    AUTHORS:
    - Edinah K. Gnang
    """
    M = M.transpose()
    # Initilization of the list of Permutations
    P = Permutations(range(M.nrows()))
    # Initialization of the list
    L=[]
    # Initialization of the final matrix
    Q = Matrix(SR, zero_matrix(M.nrows(), M.ncols()))
    for p in P:
        # Initialization of the temporary matrix.
        Rs = Matrix(SR, zero_matrix(M.nrows(), M.ncols()))
        # Initialization of the first vector.
        # the implementation assumes that the first vector is non zero
        Rs[p[0],:] = M[p[0],:]
        for i in range(1,M.nrows()):
            v = M[p[i],:]
            v = v-sum([(v*Rs[p[j],:].transpose())[0,0]/sum(Rs[p[j],s]^2 for s in range(Rs.ncols()))*Rs[p[j],:] for j in range(i) if not sum(Rs[p[j],s]^2 for s in range(Rs.ncols())).is_zero()])
            Rs[p[i],:] = v
        L.append(copy(Rs).transpose())
    return L

def Nearest_RealRow_Gram_Schmidt(M):
    """
    Takes a square matrix as input and returns an orthognal martix with minimal row distance.

    EXAMPLES:

    ::

        sage: Nearest_RealRow_Gram_Schmidt(Matrix(QQ,2,2,[1, -1/2, -1, 0]))
        Distance from the zero matrix  9/4
        The current error for p= [0, 1] is 4/5
        The current error for p= [1, 0] is 1
        [
        [   1 -1/2]  [   1 -1/2]     
        [  -1    0], [-1/5 -2/5], 4/5
        ]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Normalizing the Rows of the matrix
    for i in range(M.nrows()):
        M[i,:]=M[i,:]
    # Initilization of the list of Permutations
    P = Permutations(range(M.nrows()))
    # Initialization of the current error.
    cur_err=sum(nt^2 for nt in M.list())
    print("Distance from the zero matrix "+str(cur_err))
    # Initialization of the final matrix
    Q = Matrix(QQ, zero_matrix(M.nrows(), M.ncols()))
    for p in P:
        # Initialization of the temporary matrix.
        Rs = Matrix(QQ, zero_matrix(M.nrows(), M.ncols()))
        # Initialization of the first vector.
        # the implementation assumes that the first vector is non zero
        Rs[0,:]=M[p[0],:]
        for i in range(1,M.nrows()):
            v=M[p[i],:]
            v = v-sum([(v*Rs[p[j],:].transpose())[0,0]/sum(Rs[p[j],s]^2 for s in range(Rs.ncols()))*Rs[p[j],:] for j in range(i) if not sum(Rs[p[j],s]^2 for s in range(Rs.ncols())).is_zero()])
            Rs[p[i],:] = v
        tmp_err=sum(nt^2 for nt in (Rs-M).list())
        print('The current error for p= '+str(p)+' is '+str(tmp_err))
        if(tmp_err < cur_err):
            cur_err=tmp_err
            Q[:,:] = Rs[:,:]
    return [M, Q, cur_err]

def Nearest_RealColumn_Gram_Schmidt(M):
    """
    Takes a square matrix as input and returns an orthognal martix with minimal column distance.
    This function simply uses the Nearest Real Row Gram-Schmidt function.

    EXAMPLES:

    ::

        sage: Nearest_RealColumn_Gram_Schmidt(Matrix(QQ,2,2,[1, -1/2, -1, 0]))
        Distance from the zero matrix  9/4
        The current error for p= [0, 1] is 1/8
        The current error for p= [1, 0] is 1/4
        [
        [   1 -1/2]  [   1 -1/4]     
        [  -1    0], [  -1 -1/4], 1/8
        ]


    AUTHORS:
    - Edinah K. Gnang
    """
    [M, Q, cur_err]=Nearest_RealRow_Gram_Schmidt(M.transpose())
    return [M.transpose(), Q.transpose(), cur_err]

def GeneralHypermatrixConstrainedOrthogonalization(H, X):
    """
    Implements the general hypermatrix constrained orthogonalization algorithm.
    The function simply returns the solution to the corresponding constraints.


    EXAMPLES:

    ::

        sage: od=2; sz=2; Sln=GeneralHypermatrixConstrainedOrthogonalization(HM(*([sz for i in range(od)]+['h'])), HM(*([sz for i in range(od)]+['x']))); Sln
        [x00 == 1/2*(h00*h10 - h01*h11)/x10, x01 == -1/2*(h00*h10 - h01*h11)/x11]
        sage: H=HM(*([sz for i in range(od)]+['x'])).subs(dict([(s.lhs(),s.rhs()) for s in Sln]))
        sage: Prod(*[H.transpose(od-i) for i in range(od)])
        [[1/4*(h00*h10 - h01*h11)^2/x10^2 + 1/4*(h00*h10 - h01*h11)^2/x11^2, 0], [0, x10^2 + x11^2]]


    AUTHORS:
    - Edinah K. Gnang
    """
    if H.order() > 1:
        # Initializing the order and the size
        od = H.order(); szL = H.dimensions()
        # Constrained Orthogonalization procedure.
        DltL=GeneralHypermatrixKroneckerDeltaL(od,szL[0])
        Dlt=sum(DltL)
        # Loop initializing the hypermartrix enrtry list 
        Lx=[]; Lh=[]
        for t in range(H.n(1)):
            Lx=Lx+((HM(*(szL+['one']))-Dlt).elementwise_product(ProdB(*([X.transpose(od-j) for j in range(od)]+[DltL[t]])))).list()
            Lh=Lh+((HM(*(szL+['one']))-Dlt).elementwise_product(ProdB(*([H.transpose(od-j) for j in range(od)]+[DltL[t]]))-(1/H.n(1))*Prod(*[H.transpose(od-j) for j in range(od)]))).list()
        # Initialization of the equation
        EqL=Set([Lx[i]==Lh[i] for i in range(len(Lx)) if not Lx[i].is_zero()]).list()
        # Formating the constraints
        if len(X.list())==len(Set(X.list()).list()):
            VrbL=X.list()
        else:
            VrbL=Set(X.list()).list()
        [A, b]=multiplicativeConstraintFormator(EqL, VrbL)
        # Initialization of the vector of variables
        v=Matrix(SR, A.ncols(), 1, VrbL)
        # returning the solutions to the system obtained via Gauss-Jordan elimination
        return multiplicative_linear_solver(A, b, v, v)
    else :
        raise ValueError("The input hypermatrix must be order greater then 1")

def GeneralHypermatrixConstrainedOrthogonalizationII(H, X, sz):
    """
    Implements the general hypermatrix constrained orthogonalization algorithm.
    The difference with the implementation above is that we can change here the
    size of the support of the inner product. This enables orthogonalization
    of a smaller hypermatrix block. The function returns the solutions


    EXAMPLES:

    ::

        sage: od=2; sz=2; Sln=GeneralHypermatrixConstrainedOrthogonalizationII(HM(*([sz for i in range(od)]+['h'])), HM(*([sz for i in range(od)]+['x'])), 2); Sln
        [x00 == 1/2*(h00*h10 - h01*h11)/x10, x01 == -1/2*(h00*h10 - h01*h11)/x11]
        sage: H=HM(*([sz for i in range(od)]+['x'])).subs(dict([(s.lhs(),s.rhs()) for s in Sln]))
        sage: Prod(*[H.transpose(od-i) for i in range(od)])
        [[1/4*(h00*h10 - h01*h11)^2/x10^2 + 1/4*(h00*h10 - h01*h11)^2/x11^2, 0], [0, x10^2 + x11^2]]
        sage: Ha=HM(3,2,'a'); A=ZeroPadding(Ha)
        sage: Hx=HM(3,2,'x'); X=ZeroPadding(Hx)
        sage: k0,k1=var('k0,k1')
        sage: Sln=GeneralHypermatrixConstrainedOrthogonalizationII(A, X, 2)
        sage: H=Hx.subs(dict([(s.lhs(),s.rhs()) for s in Sln])).subs(dict([(k0,0), (k1,0)]))
        sage: Prod(H,H.transpose())
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]]


    AUTHORS:
    - Edinah K. Gnang
    """
    if H.order() > 1:
        # Initializing the order and the size
        od = H.order(); szL = H.dimensions()
        # Constrained Orthogonalization procedure.
        DltL=GeneralHypermatrixKroneckerDeltaL(od,szL[0])
        Dlt=sum(DltL)
        # Loop initializing the hypermartrix enrtry list 
        Lx=[]; Lh=[]
        for t in range(H.n(1)):
            Lx=Lx+((HM(*(szL+['one']))-Dlt).elementwise_product(ProdB(*([X.transpose(od-j) for j in range(od)]+[DltL[t]])))).list()
            Lh=Lh+((HM(*(szL+['one']))-Dlt).elementwise_product(ProdB(*([H.transpose(od-j) for j in range(od)]+[DltL[t]]))-(1/sz)*Prod(*[H.transpose(od-j) for j in range(od)]))).list()
        # Initialization of the equation
        EqL=Set([Lx[i]==Lh[i] for i in range(len(Lx)) if not (Lx[i].is_zero() and Lh[i].is_zero())]).list()
        # Formating the constraints
        if len(X.list())==len(Set(X.list()).list()):
            VrbL=X.list()
        else:
            VrbL=Set(X.list()).list()
        [A, b]=multiplicativeConstraintFormator(EqL, VrbL)
        # Initialization of the vector of variables
        v=Matrix(SR, A.ncols(), 1, VrbL)
        # returning the solutions to the system obtained via Gauss-Jordan elimination
        return multiplicative_linear_solver(A, b, v, v)
    else :
        raise ValueError("The input hypermatrix must be order greater then 1")

def Filmus_GnangConstrainedOrthogonalizationHM(Hm):
    """
    Implements Filmus-Gnang orthogonalization procedure.


    EXAMPLES:

    ::

        sage: Q=Filmus_GnangConstrainedOrthogonalizationHM(HM(2,2,'a')); od=Q.order() 
        sage: Prod(*[Q.transpose(od-i) for i in range(od)]).printHM()
        [:, :]=
        [1/4*(a00*a10 - a01*a11)^2/x10^2 + 1/4*(a00*a10 - a01*a11)^2/x11^2                                                                 0]
        [                                                                0                                                     x10^2 + x11^2]
        sage: Q=Filmus_GnangConstrainedOrthogonalizationHM(HM(3,3,'a')).subs(k0=0,k1=0,k2=0); od=Q.order()
        sage: Prod(*[Q.transpose(od-i) for i in range(od)]).elementwise_product(HM(3,3,'one')-HM(2,3,'kronecker'))
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]]


    AUTHORS:
    - Edinah K. Gnang
    """
    if Hm.is_cubical():
        # Initialization of the parameters
        od = Hm.order()
        sz = Hm.n(0)
        Sln=GeneralHypermatrixConstrainedOrthogonalizationII(Hm, HM(*([sz for i in range(od)]+['x'])), sz)
        return HM(*([sz for i in range(od)]+['x'])).subs(dict([(s.lhs(),s.rhs()) for s in Sln if s.lhs() !=1 ]))
    else:
        raise ValueError("Expected a cubical  hypermatrix")

def GeneralHypermatrixConstrainedUncorrelatedTuples(Hl, Xl):
    """
    Implements the general hypermatrix constrained uncorrelated tuple  algorithm. 


    EXAMPLES:

    ::

        sage: Hl=[HM(2,2,'u'),HM(2,2,'v')]; Xl=[HM(2,2,'x'),HM(2,2,'y')] 
        sage: Sln=GeneralHypermatrixConstrainedUncorrelatedTuples(Hl, Xl); Sln
        [x00 == 1/2*(u00*v01 - u01*v11)/y01,
         x10 == 1/2*(u10*v00 - u11*v10)/y00,
         x01 == -1/2*(u00*v01 - u01*v11)/y11,
         x11 == -1/2*(u10*v00 - u11*v10)/y10]
        sage: Prod(Xl[0].subs(dict([(s.lhs(),s.rhs()) for s in Sln])),Xl[1].subs(dict([(s.lhs(),s.rhs()) for s in Sln])))
        [[1/2*(u00*v01 - u01*v11)*y00/y01 - 1/2*(u00*v01 - u01*v11)*y10/y11, 0], [0, 1/2*(u10*v00 - u11*v10)*y01/y00 - 1/2*(u10*v00 - u11*v10)*y11/y10]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the order and the size
    od = Hl[0].order(); dimL=[Hl[i].dimensions() for i in range(len(Hl))]
    # Constrained Orthogonalization procedure.
    DltL=GeneralHypermatrixKroneckerDeltaL(od, Hl[0].n(1))
    Dlt=sum(DltL)
    # Loop initializing the hypermartrix enrtry list 
    Lx=[]; Lh=[]
    for t in range(Hl[0].n(1)):
        #Lx=Lx+((apply(HM,dimL[t]+['one'])-Dlt).elementwise_product(apply(ProdB,[X for X in Xl]+[DltL[t]]))).list()
        Lx=Lx+((HM(*(dimL[t]+['one']))-Dlt).elementwise_product(ProdB(*([X for X in Xl]+[DltL[t]])))).list()
        #Lh=Lh+((apply(HM,dimL[t]+['one'])-Dlt).elementwise_product(apply(ProdB,[H for H in Hl]+[DltL[t]])-(1/Hl[0].n(1))*apply(Prod,[H for H in Hl]))).list()
        Lh=Lh+((HM(*(dimL[t]+['one']))-Dlt).elementwise_product(ProdB(*([H for H in Hl]+[DltL[t]]))-(1/Hl[0].n(1))*Prod(*[H for H in Hl]))).list()
    # Initialization of the equation
    EqL=Set([Lx[i]==Lh[i] for i in range(len(Lx)) if not Lx[i].is_zero()]).list()
    # Formating the constraints
    LstX=[]
    for x in Xl:
        LstX=LstX+x.list()
    if len(Set(LstX).list())==len(LstX):
        VrbL=LstX
    else:
        VrbL=Set(LstX).list()
    [A,b]=multiplicativeConstraintFormator(EqL, VrbL)
    # Initialization of the vector of variables
    v=Matrix(SR, A.ncols(), 1, VrbL)
    # returning the solutions to the system obtained via Gauss-Jordan elimination
    return multiplicative_linear_solver(A, b, v, v)

def GeneralHypermatrixConstrainedUncorrelatedTuplesII(Hl, Xl, sz):
    """
    Implements the general hypermatrix constrained uncorrelated tuple  algorithm. 


    EXAMPLES:

    ::

        sage: Hl=[HM(2,2,'u'),HM(2,2,'v')]; Xl=[HM(2,2,'x'),HM(2,2,'y')] 
        sage: Sln=GeneralHypermatrixConstrainedUncorrelatedTuplesII(Hl, Xl, 2); Sln
        [x00 == 1/2*(u00*v01 - u01*v11)/y01,
         x10 == 1/2*(u10*v00 - u11*v10)/y00,
         x01 == -1/2*(u00*v01 - u01*v11)/y11,
         x11 == -1/2*(u10*v00 - u11*v10)/y10]
        sage: Prod(Xl[0].subs(dict([(s.lhs(),s.rhs()) for s in Sln])),Xl[1].subs(dict([(s.lhs(),s.rhs()) for s in Sln])))
        [[1/2*(u00*v01 - u01*v11)*y00/y01 - 1/2*(u00*v01 - u01*v11)*y10/y11, 0], [0, 1/2*(u10*v00 - u11*v10)*y01/y00 - 1/2*(u10*v00 - u11*v10)*y11/y10]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the order and the size
    od = Hl[0].order(); dimL=[Hl[i].dimensions() for i in range(len(Hl))]
    # Constrained Orthogonalization procedure.
    DltL=GeneralHypermatrixKroneckerDeltaL(od, Hl[0].n(1))
    Dlt=sum(DltL)
    # Loop initializing the hypermartrix enrtry list 
    Lx=[]; Lh=[]
    for t in range(Hl[0].n(1)):
        #Lx=Lx+((apply(HM,dimL[t]+['one'])-Dlt).elementwise_product(apply(ProdB,[X for X in Xl]+[DltL[t]]))).list()
        Lx=Lx+((HM(*(dimL[t]+['one']))-Dlt).elementwise_product(ProdB(*([X for X in Xl]+[DltL[t]])))).list()
        #Lh=Lh+((apply(HM,dimL[t]+['one'])-Dlt).elementwise_product(apply(ProdB,[H for H in Hl]+[DltL[t]])-(1/sz)*apply(Prod,[H for H in Hl]))).list()
        Lh=Lh+((HM(*(dimL[t]+['one']))-Dlt).elementwise_product(ProdB(*([H for H in Hl]+[DltL[t]]))-(1/sz)*Prod(*[H for H in Hl]))).list()
    # Initialization of the equation
    EqL=Set([Lx[i]==Lh[i] for i in range(len(Lx)) if not (Lx[i].is_zero() and Lh[i].is_zero())]).list()
    # Formating the constraints
    LstX=[]
    for x in Xl:
        LstX=LstX+x.list()
    if len(Set(LstX).list())==len(LstX):
        VrbL=LstX
    else:
        VrbL=Set(LstX).list()
    [A,b]=multiplicativeConstraintFormator(EqL, VrbL)
    # Initialization of the vector of variables
    v=Matrix(SR, A.ncols(), 1, VrbL)
    # returning the solutions to the system obtained via Gauss-Jordan elimination
    return multiplicative_linear_solver(A, b, v, v)

def GeneralHypermatrixPerturbedUncorrelatedTuples(Hl, Xl, A):
    """
    Implements the general hypermatrix perturbation algorithm.


    EXAMPLES:

    ::

        sage: Hl=[HM(2,2,'u'), HM(2,2,'u').inverse()]; Xl=[HM(2,2,'x'), HM(2,2,'y')] 
        sage: Sln=GeneralHypermatrixPerturbedUncorrelatedTuples(Hl, Xl, HM(2,2,'a')); Sln
        [x00 == 1/2*(a01 + 2*u01/(u01*u10/u00 - u11))/y01,
         x10 == 1/2*(2*u10*(1/u00 - u01*u10/(u00^2*(u01*u10/u00 - u11))) + a10)/y00,
         x01 == 1/2*(a01 - 2*u01/(u01*u10/u00 - u11))/y11,
         x11 == 1/2*(a10 + 2*u10*u11/(u00*(u01*u10/u00 - u11)))/y10] 
        sage: Prod(Xl[0].subs(dict([(s.lhs(),s.rhs()) for s in Sln])),Xl[1].subs(dict([(s.lhs(),s.rhs()) for s in Sln])))[0,1].simplify_rational()
        a01 
        sage: Prod(Xl[0].subs(dict([(s.lhs(),s.rhs()) for s in Sln])),Xl[1].subs(dict([(s.lhs(),s.rhs()) for s in Sln])))[1,0].simplify_rational()
        a10 
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the order and the size
    od = Hl[0].order(); dimL=[Hl[i].dimensions() for i in range(len(Hl))]
    # Constrained Orthogonalization procedure.
    DltL=GeneralHypermatrixKroneckerDeltaL(od,Hl[0].n(1))
    Dlt=sum(DltL)
    # Loop initializing the hypermartrix enrtry list 
    Lx=[]; Lh=[]
    for t in range(Hl[0].n(1)):
        #Lx=Lx+((apply(HM, dimL[t]+['one'])-Dlt).elementwise_product(apply(ProdB,[X for X in Xl]+[DltL[t]]))).list()
        Lx=Lx+((HM(*(dimL[t]+['one']))-Dlt).elementwise_product(ProdB(*([X for X in Xl]+[DltL[t]])))).list()
        #Lh=Lh+((apply(HM, dimL[t]+['one'])-Dlt).elementwise_product(apply(ProdB,[H for H in Hl]+[DltL[t]])+(1/Hl[0].n(1))*A)).list()
        Lh=Lh+((HM(*(dimL[t]+['one']))-Dlt).elementwise_product(ProdB(*([H for H in Hl]+[DltL[t]]))+(1/Hl[0].n(1))*A)).list()
    # Initialization of the equation
    EqL=Set([Lx[i]==Lh[i] for i in range(len(Lx)) if not Lx[i].is_zero()]).list()
    # Formating the constraints
    LstX=[]
    for x in Xl:
        LstX=LstX+x.list()
    if len(Set(LstX).list())==len(LstX):
        VrbL=LstX
    else:
        VrbL=Set(LstX).list()
    [A,b]=multiplicativeConstraintFormator(EqL, VrbL)
    # Initialization of the vector of variables
    v=Matrix(SR, A.ncols(), 1, VrbL)
    # returning the solutions to the system obtained via Gauss-Jordan elimination
    return multiplicative_linear_solver(A, b, v, v)

def GeneralHypermatrixPerturbedUncorrelatedTuplesII(Hl, Xl, A, sz):
    """
    Implements the general hypermatrix perturbation algorithm.


    EXAMPLES:

    ::

        sage: Hl=[HM(2,2,'u'), HM(2,2,'u').inverse()]; Xl=[HM(2,2,'x'), HM(2,2,'y')] 
        sage: Sln=GeneralHypermatrixPerturbedUncorrelatedTuplesII(Hl, Xl, HM(2,2,'a'), 2); Sln
        [x00 == 1/2*(a01 + 2*u01/(u01*u10/u00 - u11))/y01,
         x10 == 1/2*(2*u10*(1/u00 - u01*u10/(u00^2*(u01*u10/u00 - u11))) + a10)/y00,
         x01 == 1/2*(a01 - 2*u01/(u01*u10/u00 - u11))/y11,
         x11 == 1/2*(a10 + 2*u10*u11/(u00*(u01*u10/u00 - u11)))/y10] 
        sage: Prod(Xl[0].subs(dict([(s.lhs(),s.rhs()) for s in Sln])),Xl[1].subs(dict([(s.lhs(),s.rhs()) for s in Sln])))[0,1].simplify_rational()
        a01 
        sage: Prod(Xl[0].subs(dict([(s.lhs(),s.rhs()) for s in Sln])),Xl[1].subs(dict([(s.lhs(),s.rhs()) for s in Sln])))[1,0].simplify_rational()
        a10 
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the order and the size
    od = Hl[0].order(); dimL=[Hl[i].dimensions() for i in range(len(Hl))]
    # Constrained Orthogonalization procedure.
    DltL=GeneralHypermatrixKroneckerDeltaL(od,Hl[0].n(1))
    Dlt=sum(DltL)
    # Loop initializing the hypermartrix enrtry list 
    Lx=[]; Lh=[]
    for t in range(Hl[0].n(1)):
        #Lx=Lx+((apply(HM, dimL[t]+['one'])-Dlt).elementwise_product(apply(ProdB,[X for X in Xl]+[DltL[t]]))).list()
        Lx=Lx+((HM(*(dimL[t]+['one']))-Dlt).elementwise_product(ProdB(*([X for X in Xl]+[DltL[t]])))).list()
        #Lh=Lh+((apply(HM, dimL[t]+['one'])-Dlt).elementwise_product(apply(ProdB,[H for H in Hl]+[DltL[t]])+(1/sz)*A)).list()
        Lh=Lh+((HM(*(dimL[t]+['one']))-Dlt).elementwise_product(ProdB(*([H for H in Hl]+[DltL[t]]))+(1/sz)*A)).list()
    # Initialization of the equation
    EqL=Set([Lx[i]==Lh[i] for i in range(len(Lx)) if not (Lx[i].is_zero() and Lh[i].is_zero())]).list()
    # Formating the constraints
    LstX=[]
    for x in Xl:
        LstX=LstX+x.list()
    if len(Set(LstX).list())==len(LstX):
        VrbL=LstX
    else:
        VrbL=Set(LstX).list()
    [A,b]=multiplicativeConstraintFormator(EqL, VrbL)
    # Initialization of the vector of variables
    v=Matrix(SR, A.ncols(), 1, VrbL)
    # returning the solutions to the system obtained via Gauss-Jordan elimination
    return multiplicative_linear_solver(A, b, v, v) 

def TriangulationCompositionStringList(n, c):
    """
    Outputs a composition list of strings associated with triangulations

     EXAMPLES:

    ::

        sage: sz=4; A=HM(sz,sz,'a').elementwise_product(HM(sz,sz,'one')-HM(2,sz,'kronecker'))
        sage: for i0 in range(1,sz):
        ....:   for i1 in range(i0):
        ....:       A[i0,i1]=0
        sage: L = TriangulationCompositionStringList(sz-1,'A'); L
        ['A.elementwise_product(Prod(A, A.elementwise_product(Prod(A, A))))',
         'A.elementwise_product(Prod(A.elementwise_product(Prod(A, A)), A))'] 
        sage: eval(L[0])
        [[0, 0, 0, a01*a03*a12*a13*a23], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
        sage: eval(L[1])
        [[0, 0, 0, a01*a02*a03*a12*a23], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

    AUTHORS:

    - Edinah K. Gnang
    """
    if n == 1:
        return [c]
    else:
        gu = []
        for i in range(1,n):
            gu = gu+[c+".elementwise_product(Prod("+g1+", "+g2+"))" for g1 in TriangulationCompositionStringList(i,c) for g2 in TriangulationCompositionStringList(n-i,c)]
        return gu

def Triangulations(A,Ha,n,sz):
    """
    Outputs a list of second order hypermatrices each of which have a single nonzero symbolic entry which
    describes a triangulation of a regular polygon on n vertices. The input matrix is meant to be 
    upper-triangular matrices.

     EXAMPLES:

    ::

        sage: sz=4
        sage: A=HM(sz,sz,'a').elementwise_product(HM(sz,sz,'one')-HM(2,sz,'kronecker'))
        sage: for i0 in range(1,sz):
        ....:   for i1 in range(i0):
        ....:       A[i0,i1]=0
        sage: L=Triangulations(A,A,sz-1,sz)
        sage: L[0]
        [[0, 0, 0, a01*a03*a12*a13*a23], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

    AUTHORS:

    - Edinah K. Gnang
    """
    if n == 1:
        return [A]
    else:
        gu = []
        for i in range(1,n):
            gu = gu+[Ha.elementwise_product(Prod(g1,g2)).expand() for g1 in Triangulations(A,Ha,i,sz) for g2 in Triangulations(A,Ha,n-i,sz)]
        return gu

def TriangulationsII(A,Ha,n,sz):
    """
    Outputs a list of second order hypermatrices each of which have a single nonzero symbolic entry which
    describes a triangulation of a regular polygon on n vertices. The input matrix is meant to be 
    upper-triangular matrices. The difference with the implementaiton above is that we do not expand the
    polynomials. This implementation is well suited for non commuting variables defined over free fields.

     EXAMPLES:

    ::

        sage: sz=4
        sage: A=HM(sz,sz,'a').elementwise_product(HM(sz,sz,'one')-HM(2,sz,'kronecker'))
        sage: for i0 in range(1,sz):
        ....:   for i1 in range(i0):
        ....:       A[i0,i1]=0
        sage: L=TriangulationsII(A,A,sz-1,sz)
        sage: L[0]
        [[0, 0, 0, a01*a03*a12*a13*a23], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

    AUTHORS:

    - Edinah K. Gnang
    """
    if n == 1:
        return [A]
    else:
        gu = []
        for i in range(1,n):
            gu = gu+[Prod(g1,g2).elementwise_product(Ha) for g1 in TriangulationsII(A,Ha,i,sz) for g2 in TriangulationsII(A,Ha,n-i,sz)]
        return gu

def generate_triangulation_script(sz):
    """
    Creates a sage file which corresponds to a script
    which computes triangulation using non-commutative
    variables. The script starts with a stricly upper-
    triangular symbolic adjacency matrix whose non-
    zero entries defined as free variables.


    EXAMPLES:

    ::

        sage: generate_triangulation_script(4)
        sage: load('triangulation_4.sage')
        sage: L[0].printHM()
        [:, :]=
        [                  0                   0                   0 a01*a12*a23*a13*a03]
        [                  0                   0                   0                   0]
        [                  0                   0                   0                   0]
        [                  0                   0                   0                   0]
        sage: L[1].printHM()
        [:, :]=
        [                  0                   0                   0 a01*a12*a02*a23*a03]
        [                  0                   0                   0                   0]
        [                  0                   0                   0                   0]
        [                  0                   0                   0                   0]
        sage: from subprocess import call
        sage: call("rm triangulation_4.sage", shell=True)
        0

        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the hypermatrices.
    A=HM(sz,sz,'a')
    for i in rg(sz):
        for j in rg(sz):
            if i >= j:
                A[i,j]=0
    Al=[HM(sz,sz,'a')[i,j] for i in rg(sz) for j in rg(sz) if i<j]
    # Creating the file name string.
    filename='triangulation_'+str(sz)+'.sage'
    # Opening the file
    f=open(filename,'w')
    #f.write('# Loading the Hypermatrix Package\n')
    #f.write("load('./Hypermatrix_Algebra_tst.sage')\n\n")
    f.write('# Initializing the number of constraints and the number of variableas\n')
    f.write('sz='+str(sz)+'\n\n')
    f.write('# Initialization of the variables\n')
    f.write("La=[HM(sz,sz,'a')[i,j] for i in rg(sz) for j in rg(sz) if i<j]\n")
    f.write('# Initializing the free variables\n')
    f.write('F=FreeAlgebra(QQ,len(La),La)\n')
    f.write('F.<'+str(Al)[1:len(str(Al))-1]+'>=FreeAlgebra(QQ,len(La))\n\n')
    f.write('# Initialization of the hypermatrices with symbolic variable entries which do not commute\n')
    f.write('Ha=HM(sz,sz,'+str(A.list())+')\n')
    f.write('L=TriangulationsII(Ha,Ha,sz-1,sz)')
    # Closing the file
    f.close()

def TriangulationGraphsEdgeList(sz):
    """
    Takes as input the size paramater which corresponds to the number of vertices of the graph
    and outputs list of triangulation of the convex regular polygon. Each graph in the list
    is describe by as a list of edges.
    
     EXAMPLES:

    ::

        sage: TriangulationGraphsEdgeList(4)
        [[a01, a03, a12, a13, a23], [a01, a02, a03, a12, a23]]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the hypermatrix 
    A=HM(sz,sz,'a').elementwise_product(HM(sz,sz,'one')-HM(2,sz,'kronecker'))
    for i0 in range(1,sz):
        for i1 in range(i0):
            A[i0,i1]=0
    # Computing the list of triangulations
    L=Triangulations(A,A,sz-1,sz)
    # Initializing the list which will store triangulations
    # as a list of edges  
    list_of_graphs=[]
    for h in L:
        list_of_graphs.append((Set(h.list()).list())[1].operands())
    return list_of_graphs

def TriangulationGraphsAdjacencyMatrix(sz):
    """
    Takes as input the size paramater which corresponds to the number of vertices of the graph
    and outputs list of triangulation of the convex regular polygon. Each graph in the list
    is describe by its adjacency matrix.
    
     EXAMPLES:

    ::

        sage: TriangulationGraphsAdjacencyMatrix(4)
        [
        [0 1 0 1]  [0 1 1 1]
        [1 0 1 1]  [1 0 1 0]
        [0 1 0 1]  [1 1 0 1]
        [1 1 1 0], [1 0 1 0]
        ]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Obtaining the list of triangles from edge list
    L=TriangulationGraphsEdgeList(sz)
    # Initializing the list of graphs
    list_of_graphs=[]
    for l in L:
        Ha=HM(sz,sz,'a'); Ht=HM(sz,sz,'zero')
        for i in range(1,sz):
            for j in range(i):
                if Ha[j,i] in l:
                    Ht[j,i]=1; Ht[i,j]=1
        # Initialization of the result
        list_of_graphs.append(Matrix(SR,Ht.listHM()))
    return list_of_graphs

def TriangulationGraphsTriangleList(sz):
    """
    Takes as input the size paramater which corresponds to the number of vertices of the graph
    and outputs list of triangulation of the convex regular polygon. Each graph in the list
    is describe by as a list of triangle specified by their edges.
    
     EXAMPLES:

    ::

        sage: len(TriangulationGraphsTriangleList(6))
        14


    AUTHORS:

    - Edinah K. Gnang
    """
    # Obtaining the list of triangles from edge list
    L=TriangulationGraphsEdgeList(sz)
    # Initializing the list of graphs
    list_of_graphs=[]
    for l in L:
        # Initializing the color variables
        Ha=HM(sz,sz,'a')
        # Initialzing the colored hypermatrix
        Ht=Ha
        for i in range(sz):
            for j in range(sz):
                if Ha[i,j] not in l:
                    Ht[i,j]=0
        # Initialization of the result
        Hr=GeneralHypermatrixHadamardProduct(Prod(Ht,Ht),Ht)
        list_of_graphs.append([f.operands() for f in Set(Hr.list()).difference(Set([0])).list()])
    return list_of_graphs

def TriangulationGraphsDualAdjacencyMatrixList(sz):
    """
    Takes as input the size paramater which corresponds to the number of vertices of the graph
    and outputs list of descriptions of triangulation of the convex regular polygon. Each list
    is made up of a dual adjancency matrix followed by list of triangle specified by their edges.
    
     EXAMPLES:

    ::

        sage: len(TriangulationGraphsDualAdjacencyMatrixList(5))
        5
        
    AUTHORS:
    - Edinah K. Gnang
    """
    # Obtaining the list of graphs as list of triangles
    L = TriangulationGraphsTriangleList(sz)
    L2= TriangulationGraphsEdgeList(sz)
    list_of_graphs=[]
    for l in range(len(L)):
        # Initialization of the temporary adjacency matrix
        TmpA = Matrix(SR,HM(sz-2,sz-2,'zero').listHM())
        for i in range(1,len(L[l])):
            for j in range(i):
                if not Set(L[l][i]).intersection(Set(L[l][j])).is_empty():
                    TmpA[i,j]=1; TmpA[j,i]=1
        list_of_graphs.append([copy(TmpA),L2[l],L[l]]) 
    return list_of_graphs

def TriangulationListing(A,B,m,sz):
    """
    Outputs a list of third order hypermatrices of size sz x sz x 1 each describing a
    triangulation of a regular polygon on n vertices. The input hypermatrix A of size
    sz x sz x 1 is meant to be a symbolic upper-triangular with free variable entries.
    defined over free fields. This implements is construct centric.


     EXAMPLES:

    ::

        sage: sz=4;a01,a02,a03,a12,a13,a23,b01,b02,b03,b12,b13,b23,c01,c02,c03,c12,c13,c23,z0,z1=\
        ....: var('a01,a02,a03,a12,a13,a23,b01,b02,b03,b12,b13,b23,c01,c02,c03,c12,c13,c23,z0,z1')
        sage: F=FreeAlgebra(QQ,20,[a01,a02,a03,a12,a13,a23,b01,b02,b03,b12,b13,b23,c01,c02,c03,c12,c13,c23,z0,z1])
        sage: F.<a01,a02,a03,a12,a13,a23,b01,b02,b03,b12,b13,b23,c01,c02,c03,c12,c13,c23,z0,z1>=FreeAlgebra(QQ,20)
        sage: A=HM(sz,sz,1,[QQ(0),QQ(0),QQ(0),QQ(0),a01,QQ(0),QQ(0),QQ(0),a02,a12,QQ(0),QQ(0),a03,a13,a23,0])
        sage: Tb=HM(sz,sz,[QQ(0),QQ(0),QQ(0),QQ(0),b01,QQ(0),QQ(0),QQ(0),b02,b12,QQ(0),QQ(0),b03,b13,b23,0])
        sage: B=HM(sz,sz,sz,[Tb[i,j] for k in rg(sz) for j in rg(sz) for i in rg(sz)])*z0*z1
        sage: len(TriangulationListing(A,B,sz-1,sz))
        2


    AUTHORS:

    - Edinah K. Gnang
    """
    if m == 1:
        return [A]
    else:
        gu = []
        for i in range(1,m):
            gu = gu+[GProdII([g1,B,g2],sum,[z0,z1],1) for g1 in TriangulationListing(A,B,i,sz) for g2 in TriangulationListing(A,B,m-i,sz)]
        return gu

def generate_triangulation_scriptII(sz):
    """
    Creates a sage file which corresponds to a script
    which computes triangulation using non-commutative
    variables. The script starts with a stricly upper-
    triangular symbolic adjacency matrix whose non-
    zero entries are defined as free variables.
    this implementation is construct centric


    EXAMPLES:

    ::

        sage: generate_triangulation_scriptII(4)
        sage: load('triangulation_4.sage')
        sage: L[0].printHM()
        [:, :, 0]=
        [                  0                   0                   0 b03*a01*b13*a12*a23]
        [                  0                   0                   0                   0]
        [                  0                   0                   0                   0]
        [                  0                   0                   0                   0]
        sage: L[1].printHM()
        [:, :, 0]=
        [                  0                   0                   0 b03*b02*a01*a12*a23]
        [                  0                   0                   0                   0]
        [                  0                   0                   0                   0]
        [                  0                   0                   0                   0]
        sage: from subprocess import call
        sage: call("rm triangulation_4.sage", shell=True)
        0
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the hypermatrices and the construct variables
    A=HM(sz,sz,'a'); TpB=HM(sz,sz,'b'); z0,z1=var('z0,z1')
    for i in rg(sz):
        for j in rg(sz):
            if i >= j:
                A[i,j]=0; TpB[i,j]=0
    B=HM(sz,sz,sz,[TpB[i,j] for k in rg(sz) for j in rg(sz) for i in rg(sz)])
    # Initialization of the list which will store the variable to be set free
    Al=[HM(sz,sz,'a')[i,j] for i in rg(sz) for j in rg(sz) if i<j]
    Bl=[HM(sz,sz,'b')[i,j] for i in rg(sz) for j in rg(sz) if i<j]
    Zl=[z0,z1]
    # Creating the file name string.
    filename='triangulation_'+str(sz)+'.sage'
    # Opening the file
    f=open(filename,'w')
    f.write('# Loading the Hypermatrix Package\n')
    f.write("load('./Hypermatrix_Algebra_0.0.1.sage')\n\n")
    f.write('# Initializing the size and order parameters\n')
    f.write('sz='+str(sz)+'; od=2\n\n')
    f.write('# Initialization of the variables\n')
    f.write("z0,z1=var('z0,z1')\n")
    f.write("La=[HM(od,sz,'a','sym')[i,j] for i in rg(sz) for j in rg(sz) if i<j]\n")
    f.write("Lb=[HM(od,sz,'b','sym')[i,j] for i in rg(sz) for j in rg(sz) if i<j]\n")
    f.write("Lz=[z0,z1]\n\n")
    f.write('# Initializing the free variables\n')
    f.write('F=FreeAlgebra(QQ,len(La+Lb+Lz),La+Lb+Lz)\n')
    f.write('F.<'+str(Al+Bl+Zl)[1:len(str(Al+Bl+Zl))-1]+'>=FreeAlgebra(QQ,len(La+Lb+Lz))\n\n')
    f.write('# Initialization of the hypermatrices with symbolic variable entries which do not commute\n')
    f.write('Ha=HM(sz,sz,1,'+str(A.list())+')\n')
    f.write('Hb=HM(sz,sz,sz,'+str(B.list())+')*z0*z1\n\n')
    f.write('# Obtaining the triangulations\n')
    f.write('L=TriangulationListing(Ha,Hb,sz-1,sz)')
    # Closing the file
    f.close()

def generate_triangulation_scriptIII(sz):
    """
    Creates a sage file which corresponds to a script
    which computes triangulation using non-commutative
    variables. The script starts with a stricly upper-
    triangular symbolic adjacency matrix whose non-
    zero entries are defined as free variables.
    this implementation is construct centric


    EXAMPLES:

    ::

        sage: generate_triangulation_scriptIII(4)
        sage: load('triangulation_4.sage')
        sage: L[0].printHM()
        [:, :, 0]=
        [                  0                   0                   0 a01*a12*a23*b13*b03]
        [                  0                   0                   0                   0]
        [                  0                   0                   0                   0]
        [                  0                   0                   0                   0]
        sage: L[1].printHM()
        [:, :, 0]=
        [                  0                   0                   0 a01*a12*b02*a23*b03]
        [                  0                   0                   0                   0]
        [                  0                   0                   0                   0]
        [                  0                   0                   0                   0]
        sage: from subprocess import call
        sage: call("rm triangulation_4.sage", shell=True)
        0
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the hypermatrices and the construct variables
    A=HM(sz,sz,'a'); TpB=HM(sz,sz,'b'); z0,z1=var('z0,z1')
    for i in rg(sz):
        for j in rg(sz):
            if i >= j:
                A[i,j]=0; TpB[i,j]=0
    B=HM(sz,sz,sz,[TpB[i,j] for k in rg(sz) for j in rg(sz) for i in rg(sz)])
    # Initialization of the list which will store the variable to be set free
    Al=[HM(sz,sz,'a')[i,j] for i in rg(sz) for j in rg(sz) if i<j]
    Bl=[HM(sz,sz,'b')[i,j] for i in rg(sz) for j in rg(sz) if i<j]
    Zl=[z0,z1]
    # Creating the file name string.
    filename='triangulation_'+str(sz)+'.sage'
    # Opening the file
    f=open(filename,'w')
    f.write('# Loading the Hypermatrix Package\n')
    f.write("load('./Hypermatrix_Algebra_0.0.1.sage')\n\n")
    f.write('# Initializing the size and order parameters\n')
    f.write('sz='+str(sz)+'; od=2\n\n')
    f.write('# Initialization of the variables\n')
    f.write("z0,z1=var('z0,z1')\n")
    f.write("La=[HM(od,sz,'a','sym')[i,j] for i in rg(sz) for j in rg(sz) if i<j]\n")
    f.write("Lb=[HM(od,sz,'b','sym')[i,j] for i in rg(sz) for j in rg(sz) if i<j]\n")
    f.write("Lz=[z0,z1]\n\n")
    f.write('# Initializing the free variables\n')
    f.write('F=FreeAlgebra(QQ,len(La+Lb+Lz),La+Lb+Lz)\n')
    f.write('F.<'+str(Al+Bl+Zl)[1:len(str(Al+Bl+Zl))-1]+'>=FreeAlgebra(QQ,len(La+Lb+Lz))\n\n')
    f.write('# Initialization of the hypermatrices with symbolic variable entries which do not commute\n')
    f.write('Ha=HM(sz,sz,1,'+str(A.list())+')\n')
    f.write('Hb=z0*z1*HM(sz,sz,sz,'+str(B.list())+')\n\n')
    f.write('# Obtaining the triangulations\n')
    f.write('L=TriangulationListing(Ha,Hb,sz-1,sz)')
    # Closing the file
    f.close()

def Tetrahedralizations(A,B,n,sz):
    """
    Outpts a list of hypermatrices whoes nonzero symbolic entries describes
    a tetrahedral partition of a regular convex polytope on n vertices.
    In order to avoid degenerate hyperedges in the triangulation induced
    by the BM algebra, the input hypermatrix should be of the type illustrated
    in the example bellow.

     EXAMPLES:

    ::

        sage: sz=4; S=HM(sz,sz,sz,'zero')
        sage: for i in range(sz):
        ....:   for j in range(sz):
        ....:       for k in range(sz):
        ....:           if i<j and j<k:
        ....:               S[i,j,k]=1; S[i,k,j]=1
        sage: A=HM(sz,sz,sz,'a').elementwise_product(S)
        sage: L=Tetrahedralizations(A,A,sz-1,sz)
        sage: len(L)
        1

    AUTHORS:
    - Edinah K. Gnang
    """
    if n == 1:
        return [A]
    else:
        gu = []
        for i in range(1,n,2):
            for j in range(1,n-i,2):
                gu=gu+[B.elementwise_product(Prod(g1,g2,g3)).expand() for g1 in Tetrahedralizations(A,B,i,sz) for g2 in Tetrahedralizations(A,B,j,sz) for g3 in Tetrahedralizations(A,B,n-(i+j),sz)]
        return gu

def RandomTriangulation(A,Ha,n,sz):
    """
    Outpts a list of second order hypermatrices each of which have a single nonzero symbolic entry which
    describes a triangulation of a regular polygon on n vertices. The input matrix is meant to be 
    upper-triangular matrices.

     EXAMPLES:

    ::

        sage: sz=3
        sage: A=HM(sz,sz,'a').elementwise_product(HM(sz,sz,'one')-HM(2,sz,'kronecker'))
        sage: for i0 in range(1,sz):
        ....:   for i1 in range(i0):
        ....:       A[i0,i1]=0
        sage: L=RandomTriangulation(A,A,sz-1,sz)
        

    AUTHORS:
    - Edinah K. Gnang
    """
    if n == 1:
        return [A]
    else:
        gu = []
        j = RollLD([Ca(i)*Ca(n-i) for i in range(1,n+1)])
        return [Ha.elementwise_product(Prod(g1,g2)).expand() for g1 in RandomTriangulation(A,Ha,j,sz) for g2 in RandomTriangulation(A,Ha,n-j,sz)]

def SecondOrderDFT(m,n):
    """
    outputs second order uncorrelated pair of second order
    hypermatrices classically associated with the Discrete
    Fourier Transform.
    

    EXAMPLES:
 
    ::

        sage: [Ha, Hb]=SecondOrderDFT(2,2)
        sage: Ha
        [[1, 1], [1, -1]]
        sage: Hb
        [[1, 1], [1, -1]]
        sage: Prod(Ha,Hb).simplify()
        [[2, 0], [0, 2]]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the hypermatrices
    Ha=HM(n,m,[exp( I*2*pi*u*t/m) for t in range(m) for u in range(n)])
    Hb=HM(m,n,[exp(-I*2*pi*t*v/m) for v in range(n) for t in range(m)])
    return [Ha, Hb]

def SecondOrderDFT2(n):
    """
    outputs second order uncorrelated pair of second order
    hypermatrices classically associated with the Discrete
    Fourier Transform.
    

    EXAMPLES:
 
    ::

        sage: [Ha, Hb]=SecondOrderDFT2(2)
        sage: Prod(Ha,Hb).printHM()
        [:, :]=
        [[2] [0]]
        [[0] [2]]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    if n==2:
        # Initialization of the hypermatrices
        Ha=HM(n,n,[Matrix(SR,1,1,exp(I*2*pi*u*t/n)) for t in range(n) for u in range(n)])
        Hb=HM(n,n,[Matrix(SR,1,1,exp(I*2*pi*u*t/n)) for t in range(n) for u in range(n)])
    elif n>2:
        # Initialization of the hypermatrices
        Ha=HM(n,n,[Matrix(SR,[[cos(2*pi*u*t/n),-sin(2*pi*u*t/n)],[ sin(2*pi*u*t/n),cos(2*pi*u*t/n)]]) for t in range(n) for u in range(n)])
        Hb=HM(n,n,[Matrix(SR,[[cos(2*pi*u*t/n), sin(2*pi*u*t/n)],[-sin(2*pi*u*t/n),cos(2*pi*u*t/n)]]) for t in range(n) for u in range(n)])
    return [Ha, Hb]

def ThirdOrderDFT(m,n):
    """
    outputs third order uncorrelated tuple of third order
    hypermatrices generalizing the classical Discrete
    Fourier Transform matrix. The input m is associated
    with the order of the root of unity and the input n
    correspond to the size of the matrix.
    

    EXAMPLES:
 
    ::

        sage: [Ha, Hb, Hc]=ThirdOrderDFT(3,2)
        sage: Ha
        [[[1, 1], [1, 1/2*I*sqrt(3) - 1/2], [1, -1/2*I*sqrt(3) - 1/2]], [[1, 1], [1/2*I*sqrt(3) - 1/2, 1], [-1/2*I*sqrt(3) - 1/2, 1]]]
        sage: Hb
        [[[1, 1, 1], [1, 1/2*I*sqrt(3) - 1/2, -1/2*I*sqrt(3) - 1/2]], [[1, 1/2*I*sqrt(3) - 1/2, -1/2*I*sqrt(3) - 1/2], [1, 1, 1]]]
        sage: Hc
        [[[1, 1], [1, 1]], [[1, 1/2*I*sqrt(3) - 1/2], [1/2*I*sqrt(3) - 1/2, 1]], [[1, -1/2*I*sqrt(3) - 1/2], [-1/2*I*sqrt(3) - 1/2, 1]]]
        sage: Prod(Ha,Hb,Hc).simplify_full()
        [[[3, 0], [0, 0]], [[0, 0], [0, 3]]]
        sage: [Ha, Hb, Hc]=ThirdOrderDFT(5,5)
        sage: Prod(Ha,Hb,Hc).simplify_full()
        [[[5, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]], [[0, 0, 0, 0, 0], [0, 5, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]], [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 5, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]], [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 5, 0], [0, 0, 0, 0, 0]], [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 5]]]
        
    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the hypermatrices
    Ha=HM(n, m, n, [exp(I*2*pi*t*(u-w)^2/m) for w in range(n) for t in range(m) for u in range(n)])
    Hb=HM(n, n, m, [exp(I*2*pi*t*(u-v)^2/m) for t in range(m) for v in range(n) for u in range(n)])
    Hc=HM(m, n, n, [exp(I*2*pi*t*(v-w)^2/m) for w in range(n) for v in range(n) for t in range(m)])
    return [Ha, Hb, Hc]

def ThirdOrderDFT2(n):
    """
    outputs third order uncorrelated tuple of third order
    hypermatrices generalizing the classical Discrete
    Fourier Transform matrix. The input m is associated
    with the order of the root of unity and the input n
    correspond to the size of the matrix.
    

    EXAMPLES:
 
    ::

        sage: [Ha, Hb, Hc]=ThirdOrderDFT2(5)
        sage: Prod(Ha,Hb,Hc).expand().list()
        [
        [5 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]
        [0 5], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0],
        <BLANKLINE>
        [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]
        [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0],
        <BLANKLINE>
        [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]
        [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0],
        <BLANKLINE>
        [0 0]  [5 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]
        [0 0], [0 5], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0],
        <BLANKLINE>
        [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]
        [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0],
        <BLANKLINE>
        [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]
        [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0],
        <BLANKLINE>
        [0 0]  [0 0]  [5 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]
        [0 0], [0 0], [0 5], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0],
        <BLANKLINE>
        [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]
        [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0],
        <BLANKLINE>
        [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]
        [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0],
        <BLANKLINE>
        [0 0]  [0 0]  [0 0]  [5 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]
        [0 0], [0 0], [0 0], [0 5], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0],
        <BLANKLINE>
        [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]
        [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0],
        <BLANKLINE>
        [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]
        [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0], [0 0],
        <BLANKLINE>
        [0 0]  [0 0]  [0 0]  [0 0]  [5 0]
        [0 0], [0 0], [0 0], [0 0], [0 5]
        ]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the hypermatrices
    Ha=HM(n,n,n,[Matrix(SR,[[cos(2*pi*t*(u-w)^2/n), -sin(2*pi*t*(u-w)^2/n)], [sin(2*pi*t*(u-w)^2/n), cos(2*pi*t*(u-w)^2/n)]]) for w in range(n) for t in range(n) for u in range(n)])
    Hb=HM(n,n,n,[Matrix(SR,[[cos(2*pi*t*(u-v)^2/n), -sin(2*pi*t*(u-v)^2/n)], [sin(2*pi*t*(u-v)^2/n), cos(2*pi*t*(u-v)^2/n)]]) for t in range(n) for v in range(n) for u in range(n)])
    Hc=HM(n,n,n,[Matrix(SR,[[cos(2*pi*t*(v-w)^2/n), -sin(2*pi*t*(v-w)^2/n)], [sin(2*pi*t*(v-w)^2/n), cos(2*pi*t*(v-w)^2/n)]]) for w in range(n) for v in range(n) for t in range(n)])
    return [Ha, Hb, Hc]

def ThirdOrdercharpoly(A, c):
    """
    Outpts the solution associated with the linear dependence of powers.

     EXAMPLES:

    ::

        sage: x0,x1,x2,x3,x4,x5,x6,x7,x8=var('x0,x1,x2,x3,x4,x5,x6,x7,x8')
        sage: A=HM([[[1, 1], [-1, 2]], [[-1, 4], [1, 17]]])
        sage: [Sln, s]=ThirdOrdercharpoly(A,'A')
        sage: Tp=eval(s)
        sage: Sln
        [x0 == 186641252224723742786128720548/302131816499269505043899*x8,
         x1 == -16334511509638140306907278508/302131816499269505043899*x8,
         x2 == -39916709921235256189474933621/302131816499269505043899*x8,
         x3 == 444869993198324697111572893/302131816499269505043899*x8,
         x4 == -241741381030800897895387320/302131816499269505043899*x8,
         x5 == -122250819499746431538540833/302131816499269505043899*x8,
         x6 == 120487356897906588291809840/302131816499269505043899*x8,
         x7 == 139067504458363202088911118/302131816499269505043899*x8]
        sage: Tp.subs(dict([(s.lhs(),s.rhs()) for s in Sln]))
        [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of powers
    L=[]; Ls=[]; i=0
    while len(L) < 1+prod(A.dimensions()):
        L=L+HypermatrixCayleyHamiltonList(A, 2*i+1)
        Ls=Ls+HypermatrixCayleyHamiltonStringList(2*i+1, c)
        i=i+1
    # Initialization of the boolean variables which assert the status of the search.
    Fnd=False
    n=1+prod(A.dimensions())
    while (not Fnd) and n>1:
        # Initializing the index variables
        Indx = Set(range(len(L))).subsets(n)
        # Looping through all the subset of the appropriate size.
        for index in Indx:
            M=Matrix(SR, [L[i] for i in index]).transpose()
            #print 'Indexing set =', index
            if M.rank()==n-1:
                Fnd=True
                break
        # Initialization the result
        if M.rank()==n-1:
            b=Matrix(SR, M.nrows(), 1, HM(M.nrows(), 1, 'zero').list())
            x=Matrix(SR, M.ncols(), 1, [var('x'+str(i)) for i in range(M.ncols())])
            return [linear_solver(M, b, x, x) ," + ".join([str(x[i,0])+'*'+Ls[index[i]] for i in range(M.ncols())])]
        n=n-1
    if Fnd==False:
        return []
 
def llsa(Ha):
    """
    Outputs the rank 1 Logarithmic Least Square Approximation of 
    the input Hypermatrix Ha. The function uses the multiplicative
    least square solver to ensure that there always exist a solution.

    EXAMPLES:
 
    ::

        sage: [[Sln, U, V], MulSln]=llsa(Prod(HM(2,1,'a'), HM(1,2,'b')))
        sage: Hu = U.subs(dict([(s.lhs(),s.rhs()) for s in Sln[:3]]))
        sage: Hv = V.subs(dict([(s.lhs(),s.rhs()) for s in Sln[:3]]))
        sage: Hr = Prod(Hu,Hv)
        sage: Matrix(2,2,[Hr[i,j].subs(k0=0,k1=0,k2=0,k3=0,k4=0,k5=0,k6=0,k7=0).canonicalize_radical() for i in range(2) for j in range(2)])
        [a00*b00 a00*b01]
        [a10*b00 a10*b01]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the dimension list
    dm = Ha.dimensions()
    # Initialization of the list of Hypermatrices
    L=[]
    for i in range(Ha.order()):
        dmt=copy(dm); dmt[Integer(mod(i+1, Ha.order()))]=1
        #L.append(HM(apply(HypermatrixGenerate, dmt+['x'+str(i)+'y'])))
        L.append(HM(HypermatrixGenerate(*(dmt+['x'+str(i)+'y']))))
    # Initializing the constraints
    #Lh=apply(GeneralHypermatrixLogProduct, L).list(); Rh=Ha.list()
    Lh=GeneralHypermatrixLogProduct(*L).list(); Rh=Ha.list()
    CnstrLst=[Lh[j]==Rh[j] for j in range(prod(Ha.dimensions()))]
    # Initializing the list of variables
    VrbLst=[]
    for i in range(len(L)):
        VrbLst=VrbLst+L[i].list()
    # Formatting the constraints
    [Atmp, btmp]=ConstraintFormatorII(CnstrLst, VrbLst)
    # Initializing the least square constraints
    A=Atmp.transpose()*Atmp; b=multiplicative_matrix_product(Atmp.transpose(),btmp)
    # Solving the least square equations
    Mx=Matrix(SR,A.ncols(),1,VrbLst)
    Sln=multiplicative_least_square_linear_solver(A.transpose()*A, multiplicative_matrix_product(A.transpose(),b), Mx, Mx)
    return [[Sln]+L, multiplicative_gauss_jordan_eliminationII(A.transpose()*A, multiplicative_matrix_product(A.transpose(),b))] 

def LogarithmicRankOneApproximation(Ha):
    """
    Outputs the logarithmic BM-rank one approximationassociated with
    the input Hypermatrix Ha

    EXAMPLES:
 
    ::

        sage: [L0,L1]=LogarithmicRankOneApproximation(Prod(HM(2,1,'a'), HM(1,2,'b')))
        sage: Hu=(L0[1]).copy().subs(dict([(s.lhs(),s.rhs()) for s in L0[0]]))
        sage: Hv=(L0[2]).copy().subs(dict([(s.lhs(),s.rhs()) for s in L0[0]]))
        sage: Prod(Hu,Hv).subs(dict([(var('k'+str(s)),0) for s in range(8)])).canonicalize_radical()
        [[a00*b00, a00*b01], [a10*b00, a10*b01]]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the dimension list
    dm = Ha.dimensions()
    # Initialization of the list of Hypermatrices
    L=[]
    for i in range(Ha.order()):
        dmt=copy(dm); dmt[Integer(mod(i+1, Ha.order()))]=1
        #L.append(HM(apply(HypermatrixGenerate, dmt+['x'+str(i)+'y'])))
        L.append(HM(HypermatrixGenerate(*(dmt+['x'+str(i)+'y']))))
    # Initializing the constraints
    #Lh=apply(GeneralHypermatrixLogProduct, L).list(); Rh=Ha.list()
    Lh=GeneralHypermatrixLogProduct(*L).list(); Rh=Ha.list()
    CnstrLst=[Lh[j]==Rh[j] for j in range(prod(Ha.dimensions()))]
    # Initializing the list of variables
    VrbLst=[]
    for i in range(len(L)):
        VrbLst=VrbLst+L[i].list()
    # Formatting the constraints
    [A, b]=ConstraintFormatorII(CnstrLst, VrbLst)
    # Solving the least square equations
    Mx=Matrix(SR,A.ncols(),1,VrbLst)
    #Sln=multiplicative_linear_solver(A, b, Mx, Mx)
    Sln=multiplicative_least_square_linear_solver(A.transpose()*A, multiplicative_matrix_product(A.transpose(),b), Mx, Mx)
    #return [[Sln]+L, multiplicative_gauss_jordan_eliminationII(A,b)] 
    return [[Sln]+L, multiplicative_gauss_jordan_eliminationII(A.transpose()*A, multiplicative_matrix_product(A.transpose(),b))] 

def UncorrelatedSideLength2Triple(c1,c2,c3):
    """
    Generates a symbolic parametrization of an uncorrelated hypermatrix
    triplet Each having side length equal to 2. The inputs are three
    characters to used in the symbolic parametrization of the hyperma
    trix entries.

    EXAMPLES:

    ::

        sage: [Sln1, Sln2, Sln3]=UncorrelatedSideLength2Triple('a','b','c')
        sage: Ha=HM(2,2,2,'a').subs(dict([(s.lhs(),s.rhs()) for s in Sln1]))
        sage: Hb=HM(2,2,2,'b').subs(dict([(s.lhs(),s.rhs()) for s in Sln1]))
        sage: Hc=HM(2,2,2,'c').subs(dict([(s.lhs(),s.rhs()) for s in Sln1]))
        sage: Hd=Prod(Ha,Hb,Hc).simplify_full(); Hd
        [[[-(a010*b001*c000*c011*c101*c110 - a010*b001*c001*c010*c100*c111)/(c001*c010*c111), 0], [0, 0]], [[0, 0], [0, -(a111*b111*c000*c011*c101*c110 - a111*b111*c001*c010*c100*c111)/(c001*c010*c100)]]]
        sage: Ah=Ha.copy(); Bh=Hb.copy(); Ch=Hc.copy()
        sage: Ch[0,0,0]=Ch[0,0,0]/Hd[0,0,0]; Ch[1,0,0]=Ch[1,0,0]/Hd[0,0,0]
        sage: Ch[0,1,1]=Ch[0,1,1]/Hd[1,1,1]; Ch[1,1,1]=Ch[1,1,1]/Hd[1,1,1]
        sage: Prod(Ah,Bh,Ch).simplify_full()
        [[[1, 0], [0, 0]], [[0, 0], [0, 1]]]
        sage: Ha=HM(2,2,2,'a').subs(dict([(s.lhs(),s.rhs()) for s in Sln2]))
        sage: Hb=HM(2,2,2,'b').subs(dict([(s.lhs(),s.rhs()) for s in Sln2]))
        sage: Hc=HM(2,2,2,'c').subs(dict([(s.lhs(),s.rhs()) for s in Sln2]))
        sage: Hd=Prod(Ha,Hb,Hc).simplify_full(); Hd
        [[[-(a010*b001*c000*c011*c101*c110 - a010*b001*c001*c010*c100*c111)/(c001*c010*c111), 0], [0, 0]], [[0, 0], [0, -(a111*b111*c000*c011*c101*c110 - a111*b111*c001*c010*c100*c111)/(c001*c010*c100)]]]
        sage: Ah=Ha.copy(); Bh=Hb.copy(); Ch=Hc.copy()
        sage: Ch[0,0,0]=Ch[0,0,0]/Hd[0,0,0]; Ch[1,0,0]=Ch[1,0,0]/Hd[0,0,0]
        sage: Ch[0,1,1]=Ch[0,1,1]/Hd[1,1,1]; Ch[1,1,1]=Ch[1,1,1]/Hd[1,1,1]
        sage: Prod(Ah,Bh,Ch).simplify_full()
        [[[1, 0], [0, 0]], [[0, 0], [0, 1]]]
        sage: Ha=HM(2,2,2,'a').subs(dict([(s.lhs(),s.rhs()) for s in Sln3]))
        sage: Hb=HM(2,2,2,'b').subs(dict([(s.lhs(),s.rhs()) for s in Sln3]))
        sage: Hc=HM(2,2,2,'c').subs(dict([(s.lhs(),s.rhs()) for s in Sln3]))
        sage: Hd=Prod(Ha,Hb,Hc).simplify_full(); Hd
        [[[1, 0], [0, 0]], [[0, 0], [0, 1]]]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the hypermatrices.
    Ha=HM(2,2,2,c1); Hb=HM(2,2,2,c2); Hc=HM(2,2,2,c3)
    # Initialization fo the list of variables.
    VrbLst=Ha.list()+Hb.list()+Hc.list()
    # Initialization of the first list of constraints.
    CnstrLst1=[Ha[i,0,k]+Hb[i,j,0]+Hc[0,j,k]-Ha[i,1,k]-Hb[i,j,1]-Hc[1,j,k]==-1 for i in range(2) for j in range(2) for k in range(2) if (i!=j or j!=k or i!=k)]
    # Initialization of the second list of constraints.
    CnstrLst2=[Ha[i,1,k]+Hb[i,j,1]+Hc[1,j,k]-Ha[i,0,k]-Hb[i,j,0]-Hc[0,j,k]==-1 for i in range(2) for j in range(2) for k in range(2) if (i!=j or j!=k or i!=k)]
    # Formating the first set of constraints.
    [A1,b1]=ConstraintFormatorII(CnstrLst1, VrbLst)
    Mx1=Matrix(SR,A1.ncols(),1,VrbLst); Sln1=multiplicative_linear_solver(A1, b1, Mx1, Mx1)
    # Formating the second set of constraints.
    [A2,b2]=ConstraintFormatorII(CnstrLst2, VrbLst)
    Mx2=Matrix(SR,A2.ncols(),1,VrbLst); Sln2=multiplicative_linear_solver(A2, b2, Mx2, Mx2)
    # Initializing the product
    Htmp=Prod(Ha,Hb,Hc)
    # Initializing the Kroenker delta
    Kdlt=HM(3,2,'kronecker')
    CnstrLst3=[Htmp[i,j,k]==Kdlt[i,j,k] for i in range(2) for j in range(2) for k in range(2)]
    # Formating the last list of constraints.
    [A3, b3]=ConstraintFormatorII(CnstrLst3, Hc.list())
    Mx3=Matrix(SR,A3.ncols(),1,Hc.list()); Sln3=linear_solver(A3, b3, Mx3, Mx3)
    return [Sln1, Sln2, Sln3]

def GeneralHypermatrixTransform(Hl, X):
    """
    Implements the general hypermatrix Transform

    EXAMPLES:

    ::

        sage: sz=2; t=var('t'); Q=HM(sz,sz,[cos(t),sin(t),-sin(t), cos(t)])
        sage: Y=GeneralHypermatrixTransform([Q.transpose(),Q], HM(2,1,var_list('x',sz))); Y.canonicalize_radical()
        [[-x0*cos(t) + x1*sin(t)], [x1*cos(t) + x0*sin(t)]]
        sage: GeneralHypermatrixTransform([Q.transpose(),Q], Y).canonicalize_radical()
        [[(cos(t)^2 + sin(t)^2)*x0], [(cos(t)^2 + sin(t)^2)*x1]]
        sage: sz=2; A=HM(2,2,'a')
        sage: GeneralHypermatrixTransform([A.transpose(),A], HM(2,1,var_list('x',sz))).canonicalize_radical()
        [[a00*x0 + a01*x1], [a10*x0 + a11*x1]]
        sage: sz=2; P0=HM(sz,1,var_list('x',sz)); A=HM(sz,sz,'a'); B=i2x2(A)
        sage: P1=GeneralHypermatrixTransform([A,B],P0).canonicalize_radical()
        sage: sum(f^2 for f in P1.list()).canonicalize_radical()
        x0^2 + x1^2


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the order and the size
    od = Hl[0].order(); sz = Hl[0].n(0)
    # Constrained Orthogonalization procedure.
    DltL=GeneralHypermatrixKroneckerDeltaL(od,sz)
    Dlt=sum(DltL)
    # Loop initializing the hypermartrix enrtries
    Ly=[]
    for t in range(sz):
        #Ly.append((apply(ProdB,[X.transpose(i) for i in range(od-1,-1,-1)]+[apply(ProdB,[H for H in Hl]+[DltL[t]])])).list()[0]^(1/od))
        Ly.append((ProdB(*([X.transpose(i) for i in range(od-1,-1,-1)]+[ ProdB(*([H for H in Hl]+[DltL[t]])) ]))).list()[0]^(1/od))
    return HM(*([sz]+[1 for i in range(od-1)]+[Ly])) 

def weighted_hypermatrix_minor(A):
    """
    Outputs a list of second order minor hypermatrices
    which add up to the original hypermatrices.

    EXAMPLES:
 
    ::

        sage: L=weighted_hypermatrix_minor(HM(3,3,'a')); L
        [[[0, 0, 0], [0, 1/2*a11, a12], [0, a21, 1/2*a22]],
         [[1/2*a00, 0, a02], [0, 0, 0], [a20, 0, 1/2*a22]],
         [[1/2*a00, a01, 0], [a10, 1/2*a11, 0], [0, 0, 0]]]
        sage: sum(L)
        [[a00, a01, a02], [a10, a11, a12], [a20, a21, a22]]

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Checking that the matrix is square
    if A.is_cubical() and A.order()==2 and A.n(0)>2:
        # Initialization of the size parameter
        sz=A.n(0)
        # Initialization of the Kronecker delta hypermatrices
        Id=HM(2,sz,'kronecker')
        # Initialization of the list
        L=[HM(sz,1,[Id[i,t] for i in range(sz)]) for t in range(sz)]
        # Initialization of the list of matrices
        MtrL=[]
        for t in range(sz):
            MtrL.append((sz-2)^(-1)*(Prod(sum([L[i] for i in range(sz) if i!=t]),sum([L[i] for i in range(sz) if i!=t]).transpose())-sum([Prod(L[i],L[i].transpose()) for i in range(sz) if i!=t]))+(sz-1)^(-1)*sum([Prod(L[i],L[i].transpose()) for i in range(sz) if i!=t]))
        return [A.elementwise_product(MtrL[i]) for i in range(sz)]
    elif A.is_cubical() and A.order()==3 and A.n(0)>3:
        # Initialization of the size parameter
        sz=A.n(0)
        # Initialization of the list
        L=[]
        for t in range(sz):
            # Initialization of a temporary HM
            Tp=HM(sz,sz,sz,'zero')
            for i in range(sz):
                for j in range(sz):
                    for k in range(sz):
                        if i!=t and j!=t and k!=t:
                            if i==j and j==k:
                                Tp[i,j,k]=A[i,j,k]/(sz-1)
                            if (i==j and j!=k) or (i==k and i!=j) or (j==k and i!=k):
                                Tp[i,j,k]=A[i,j,k]/(sz-2)
                            elif i!=j and i!=k and j!=k:
                                Tp[i,j,k]=A[i,j,k]/(sz-3)
            L.append(Tp.copy())
        return L
    else:
        raise ValueError("The input hypermpatrix must be cubic and order 3  and lower.")

def SecondOrderCharpolyI(A, U, mu, nu):
    """
    Outputs second order hypermatrix characterisitic polynomial.
    Or generator of the first spectral elimination ideal.
    The inputs are matrix A for which we want to compute
    the spectral decompostion. The matrices U and V which are
    associated with the e-vectors and the diagonal matrices
    mu, nu.

    EXAMPLES:
 
    ::

        sage: Ha=HM(2,2,'a'); Hu=HM(2,2,'u'); Hd0=HM(2,2,diagonal_matrix(HM(2,'x').list()).list()); Hd1=HM(2,2,diagonal_matrix(HM(2,'y').list()).list())
        sage: SecondOrderCharpolyI(Ha, Hu, Hd0, Hd1)[0].printHM()
        [:, :]=
        [u00*x0*y0*(1/u00 - u01*u10/(u00^2*(u01*u10/u00 - u11))) + u01*u10*x1*y1/(u00*(u01*u10/u00 - u11))                                     u01*x0*y0/(u01*u10/u00 - u11) - u01*x1*y1/(u01*u10/u00 - u11)]
        [u10*x0*y0*(1/u00 - u01*u10/(u00^2*(u01*u10/u00 - u11))) + u10*u11*x1*y1/(u00*(u01*u10/u00 - u11))                           u01*u10*x0*y0/(u00*(u01*u10/u00 - u11)) - u11*x1*y1/(u01*u10/u00 - u11)] 
        sage: x,y=var('x,y')
        sage: Ha=HM(2,2,'a'); Hu=HM(2,2,'u'); Hd0=HM(2,2,diagonal_matrix([x,x]).list()); Hd1=HM(2,2,diagonal_matrix([y,y]).list())
        sage: (SecondOrderCharpolyI(Ha, Hu, Hd0, Hd1)[1]).factor()
        x^2*y^2 - a00*x*y - a11*x*y - a01*a10 + a00*a11

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Checking that the matrix is square and that mu and nu are diagonal matrices
    if A.is_cubical() and A.order()==2 and (mu.elementwise_exponent(2)-Prod(mu.transpose(),mu)).is_zero() and (nu.elementwise_exponent(2)-Prod(nu.transpose(),nu)).is_zero():
        # Initialization of the hypermatrix V
        V=U.inverse().transpose()
        return [Prod(Prod(U,mu),Prod(V,nu).transpose()), (A-Prod(Prod(U,mu),Prod(V,nu).transpose())).det()]
    else:
        raise ValueError("Not supported for the input hypermatrices.")

def SecondOrderCharpolyII(A, U, Dg):
    """
    Outputs second order hypermatrix characterisitic polynomial.
    Or generator of the first spectral elimination ideal.
    The inputs are matrix A for which we want to compute
    the spectral decompostion. The matrices U and V which are
    associated with the e-vectors and the diagonal matrices
    mu, nu.

    EXAMPLES:
 
    ::

        sage: A=HM(2,2,'a'); U=HM(2,2,'u'); Dg=HM(2,2,diagonal_matrix(HM(2,'x').list()).list())
        sage: SecondOrderCharpolyII(A, U, Dg)[2].simplify_rational().numerator().factor()
        (a10*u00^2 - a00*u00*u10 + a11*u00*u10 - a01*u10^2)*(a10*u01^2 - a00*u01*u11 + a11*u01*u11 - a01*u11^2)*(a01*a10 - a00*a11)
        sage: SecondOrderCharpolyII(A, U, Dg)[2].simplify_rational().denominator().factor()
        (a10*u00*u01 - a00*u01*u10 + a11*u00*u11 - a01*u10*u11)*(a10*u00*u01 + a11*u01*u10 - a00*u00*u11 - a01*u10*u11)

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Checking that the matrix is square and that mu and nu are diagonal matrices
    if A.is_cubical() and A.order() == 2 and (Dg.elementwise_exponent(2)-Prod(Dg.transpose(),Dg)).is_zero():
        # Initialization of the hypermatrix V
        V=U.inverse().transpose()
        # Initialization of the Kronecker delta list
        L=GeneralHypermatrixKroneckerDeltaL(2,A.n(0))
        # Initializing the equations
        Eq=[(A-ProdB(U,Prod(Dg,V.transpose()),L[i])).det()==0 for i in range(A.n(0))]
        # Initializing the list of equations
        Sln=solve(Eq, [Dg[i,i] for i in range(Dg.n(0))])[0]
        return [(A-ProdB(U,Prod(Dg,V.transpose()),L[i])).det() for i in range(A.n(0))]+[(A-Prod(U,Prod(Dg,V.transpose()))).det().subs(Sln)]
    else:
        raise ValueError("Not supported for the input hypermatrices.")

def GeneralUncorrelatedComposition(Lu, La, Lf): 
    """
    Implements the generic composition procedure to be used for the spectral
    decomposition procedure. The hypermatrix inputs are three lists.


    EXAMPLES:

    ::

        sage: Lu=[HM(2,2,'u'), HM(2,2,'v')]
        sage: Sln=GeneralUncorrelatedComposition(Lu, [HM(2,2,'a'),HM(2,2,'a').inverse()], [HM(2,2,'f'),HM(2,2,'f').inverse()])
        sage: U=HM(2,2,'u').subs(dict([(s.lhs(),s.rhs()) for s in Sln if s.lhs()!=1]))
        sage: V=HM(2,2,'v').subs(dict([(s.lhs(),s.rhs()) for s in Sln if s.lhs()!=1]))
        sage: Prod(U,V)[1,0].is_zero()
        True
        sage: Prod(U,V)[0,1].is_zero()
        True
        sage: Hl1=[HM(2,2,'a'),HM(2,2,'b')]; Xl1=[HM(2,2,'s'),HM(2,2,'t')]
        sage: Hl2=[HM(2,2,'f'),HM(2,2,'g')]; Xl2=[HM(2,2,'x'),HM(2,2,'y')]
        sage: Sln1=GeneralHypermatrixConstrainedUncorrelatedTuples(Hl1,Xl1)
        sage: Sln2=GeneralHypermatrixConstrainedUncorrelatedTuples(Hl2,Xl2)
        sage: La=[Hx.subs(dict([(s.lhs(),s.rhs()) for s in Sln1 if s.lhs()!=1])) for Hx in Xl1]
        sage: Tmp1=HM(2,2,'one')
        sage: Tmp1[0,0]=1/Prod(La[0],La[1])[0,0];Tmp1[0,1]=1/Prod(La[0],La[1])[0,0]
        sage: Tmp1[1,0]=1/Prod(La[0],La[1])[1,1];Tmp1[1,1]=1/Prod(La[0],La[1])[1,1]
        sage: Lf=[Hx.subs(dict([(s.lhs(),s.rhs()) for s in Sln2 if s.lhs()!=1])) for Hx in Xl2]
        sage: Tmp2=HM(2,2,'one')
        sage: Tmp2[0,0]=1/Prod(Lf[0],Lf[1])[0,0];Tmp2[0,1]=1/Prod(Lf[0],Lf[1])[0,0]
        sage: Tmp2[1,0]=1/Prod(Lf[0],Lf[1])[1,1];Tmp2[1,1]=1/Prod(Lf[0],Lf[1])[1,1]
        sage: Lu=[HM(2,2,'u'), HM(2,2,'v')]
        sage: Sln=GeneralUncorrelatedComposition(Lu, [La[0].elementwise_product(Tmp1),La[1]], [Lf[0].elementwise_product(Tmp2),Lf[1]])
        sage: Prod(U,V)[1,0].is_zero()
        True
        sage: Prod(U,V)[0,1].is_zero()
        True
        sage: Lu=[HM(2,2,2,'u'), HM(2,2,2,'v'), HM(2,2,2,'w')]
        sage: # Initialization of the first uncorrelated tuple
        sage: [Sln1, Sln2, Sln3]=UncorrelatedSideLength2Triple('a','b','c')
        sage: Ha=HM(2,2,2,'a').subs(dict([(s.lhs(),s.rhs()) for s in Sln1]))
        sage: Hb=HM(2,2,2,'b').subs(dict([(s.lhs(),s.rhs()) for s in Sln1]))
        sage: Hc=HM(2,2,2,'c').subs(dict([(s.lhs(),s.rhs()) for s in Sln1]))
        sage: Hd=Prod(Ha,Hb,Hc).simplify()
        sage: Ah=Ha.copy(); Bh=Hb.copy(); Ch=Hc.copy()
        sage: Ch[0,0,0]=Ch[0,0,0]/Hd[0,0,0]; Ch[1,0,0]=Ch[1,0,0]/Hd[0,0,0]
        sage: Ch[0,1,1]=Ch[0,1,1]/Hd[1,1,1]; Ch[1,1,1]=Ch[1,1,1]/Hd[1,1,1]
        sage: La=[Ah, Bh, Ch]
        sage: # Initialization of the second uncorrelated tuple
        sage: [Tln1, Tln2, Tln3]=UncorrelatedSideLength2Triple('f','g','h')
        sage: Hf=HM(2,2,2,'f').subs(dict([(s.lhs(),s.rhs()) for s in Tln1]))
        sage: Hg=HM(2,2,2,'g').subs(dict([(s.lhs(),s.rhs()) for s in Tln1]))
        sage: Hh=HM(2,2,2,'h').subs(dict([(s.lhs(),s.rhs()) for s in Tln1]))
        sage: He=Prod(Hf,Hg,Hh).simplify()
        sage: Fh=Hf.copy(); Gh=Hg.copy(); Hh=Hh.copy()
        sage: Hh[0,0,0]=Hh[0,0,0]/He[0,0,0]; Hh[1,0,0]=Hh[1,0,0]/He[0,0,0]
        sage: Hh[0,1,1]=Hh[0,1,1]/He[1,1,1]; Hh[1,1,1]=Hh[1,1,1]/He[1,1,1]
        sage: Lf=[Fh, Gh, Hh]
        sage: # Performing the composition
        sage: Sln=GeneralUncorrelatedComposition(Lu, La, Lf)
        sage: Hu=HM(2,2,2,'u').subs(dict([(s.lhs(),s.rhs()) for s in Sln if s.lhs()!=1]))
        sage: Hv=HM(2,2,2,'v').subs(dict([(s.lhs(),s.rhs()) for s in Sln if s.lhs()!=1]))
        sage: Hw=HM(2,2,2,'w').subs(dict([(s.lhs(),s.rhs()) for s in Sln if s.lhs()!=1]))
        sage: Hx=Prod(Hu,Hv,Hw).simplify()
        sage: Uh=Hu.copy(); Vh=Hv.copy(); Wh=Hw.copy()
        sage: Wh[0,0,0]=Wh[0,0,0]/Hx[0,0,0]; Wh[1,0,0]=Wh[1,0,0]/Hx[0,0,0]
        sage: Wh[0,1,1]=Wh[0,1,1]/Hx[1,1,1]; Wh[1,1,1]=Wh[1,1,1]/Hx[1,1,1]
        sage: Prod(Uh,Vh,Wh).simplify_full()
        [[[1, 0], [0, 0]], [[0, 0], [0, 1]]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the order paramter
    od=Lu[0].order()
    # Initializing the size parameter
    sz=Lu[0].n(0)
    # Initializing the list of hypermatrices
    DltL=GeneralHypermatrixKroneckerDeltaL(od,sz)
    # We seek to solve for Hu,Hv,Hw for input hypermatrices Ha,Hb,Hc, He,Hf,Hg
    L1=[]; L2=[]
    for i in range(sz): 
        #L1=L1+apply(ProdB,[Hu for Hu in Lu]+[DltL[i]]).list()
        L1=L1+ProdB(*([Hu for Hu in Lu]+[DltL[i]])).list()
        #L2=L2+apply(ProdB,[Ha for Ha in La]+[apply(ProdB,[Hf for Hf in Lf]+[DltL[i]])]).list()
        L2=L2+ProdB(*([Ha for Ha in La]+[ProdB(*([Hf for Hf in Lf]+[DltL[i]]))])).list()
    # Initialization of the equation
    EqL=[L1[i]==L2[i] for i in range(len(L1))]
    LstX=[]
    for x in Lu:
        LstX=LstX+x.list()
    if len(Set(LstX).list())==len(LstX):
        VrbL=LstX
    else:
        VrbL=Set(LstX).list()
    # Formating the constraints
    [A,b]=multiplicativeConstraintFormator(EqL, VrbL)
    # Initialization of the vector of variables
    v=Matrix(SR, A.ncols(), 1, VrbL)
    # returning the solutions to the system obtained via Gauss-Jordan elimination
    return multiplicative_linear_solver(A, b, v, v)

def ThirdOrdeCharpolyI(A, Mu, Mv, Mw):
    """
    Outputs the third order hypermatrix characterisitic polynomial.
    Can alternatively generator of the first spectral elimination ideal.
    Assumes that the hypermatrix U, V, W arises from a parametrization
    of uncorrelated tuples but checks that the the inputs hypermatrices
    Du, Dv, Dw are diagonal hypremartices.


    EXAMPLES:
 
    ::

        sage: Mu=HM(2,2,'f','sym'); Mv=HM(2,2,'g','sym'); Mw=HM(2,2,'h','sym')
        sage: ThirdOrdeCharpolyI(HM(2,2,2,'a'), Mu, Mv, Mw)
        (f01^2*g01^2*h01^2 - a111)*a001*a010*a100 - (f00^2*g00^2*h00^2 - a000)*a011*a101*a110 
 

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Checking that the matrix is square and that mu and nu are diagonal matrices
    if A.is_cubical() and A.order() == 3 and A.n(0)==2:
        # Initializing the hypermatrix B
        B=HM(A.n(0),A.n(1),A.n(2),'zero')
        for i in range(A.n(0)):
            for j in range(A.n(1)):
                for k in range(A.n(2)):
                    if i==j and j==k:
                        B[i,j,k]=A[i,j,k]-(Mu[0,i]*Mv[0,i]*Mw[0,i])^2
                    else:
                        B[i,j,k]=A[i,j,k]
        return B.det()
    else:
        raise ValueError("Not supported for the input hypermatrices.")

def ThirdOrdeCharpolyII(A, U, V, W, Du, Dv, Dw):
    """
    Outputs the third order hypermatrix characterisitic polynomial.
    Can alternatively generator of the first spectral elimination ideal.
    Assumes that the hypermatrix U, V, W arises from a parametrization
    of uncorrelated tuples but checks that the the inputs hypermatrices
    Du, Dv, Dw are diagonal hypremartices.


    EXAMPLES:
 
    ::

        sage: Du=HM(Matrix(2,2,HM(2,2,'f','sym').listHM())); Dv=HM(Matrix(2,2,HM(2,2,'g','sym').listHM())); Dw=HM(Matrix(2,2,HM(2,2,'h','sym').listHM()))
        sage: [Sln1, Sln2, Sln3]=UncorrelatedSideLength2Triple('u','v','w')
        sage: Hu=HM(2,2,2,'u').subs(dict([(s.lhs(),s.rhs()) for s in Sln1]))
        sage: Hv=HM(2,2,2,'v').subs(dict([(s.lhs(),s.rhs()) for s in Sln1]))
        sage: Hw=HM(2,2,2,'w').subs(dict([(s.lhs(),s.rhs()) for s in Sln1]))
        sage: Hd=Prod(Hu,Hv,Hw).simplify_full()
        sage: Uh=Hu.copy(); Vh=Hv.copy(); Wh=Hw.copy()
        sage: Wh[0,0,0]=Wh[0,0,0]/Hd[0,0,0]; Wh[1,0,0]=Wh[1,0,0]/Hd[0,0,0]
        sage: Wh[0,1,1]=Wh[0,1,1]/Hd[1,1,1]; Wh[1,1,1]=Wh[1,1,1]/Hd[1,1,1]
        sage: ThirdOrdeCharpolyII(HM(2,2,2,'a'), Uh, Vh, Wh, Du, Dv, Dw)[0]
        -(f01^2*g01^2*h01^2*u111*v111*w000*w011*w101*w110/(u111*v111*w000*w011*w101*w110 - u111*v111*w001*w010*w100*w111) - a111)*(f00*f01*g00*g01*h00^2*u110*v101*w001*w010*w100*w111/(u010*v001*w000*w011*w101*w110 - u010*v001*w001*w010*w100*w111) - a100)*(f00*f01*g00^2*h00*h01*u011*v001*w101 + a001)*(f00^2*g00*g01*h00*h01*u010*v011*w110 + a010) + (f00^2*g00^2*h00^2*u010*v001*w000*w011*w101*w110/(u010*v001*w000*w011*w101*w110 - u010*v001*w001*w010*w100*w111) - a000)*(f00*f01*g00*g01*h01^2*u011*v011*w001*w010*w100*w111/(u111*v111*w000*w011*w101*w110 - u111*v111*w001*w010*w100*w111) - a011)*(f01^2*g00*g01*h00*h01*u111*v101*w101 + a101)*(f00*f01*g01^2*h00*h01*u110*v111*w110 + a110)


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Checking that the matrix is square and that mu and nu are diagonal matrices
    if A.order()==3 and A.n(0)==2 and (Prod(Du.transpose(),Du.transpose(2),Du)-Du.elementwise_exponent(3)).is_zero() and (Prod(Dv.transpose(),Dv.transpose(2),Dv)-Dv.elementwise_exponent(3)).is_zero() and (Prod(Dw.transpose(),Dw.transpose(2),Dw)-Dw.elementwise_exponent(3)).is_zero():
        DltL=GeneralHypermatrixKroneckerDeltaL(A.order(),A.n(0))
        return [(A-ProdB(Prod(U,Du,Du.transpose()), Prod(Dv,V,Dv.transpose(2)), Prod(Dw.transpose(),Dw.transpose(2),W), DltL[t])).det() for t in range(2)]
    else:
        raise ValueError("Not supported for the input hypermatrices.") 

def Additive_Determinant_Matrix(f, g, t):
    """
    Outputs a symbolic second order hypermatrix whose determinant is the sum of
    the two polynomials.

    EXAMPLES:

    ::

        sage: x=var('x'); f=x^2+sum(HM(2,'a').list()[i]*x^i for i in range(2)); g=x^3+sum(HM(3,'b').list()[i]*x^i for i in range(3))
        sage: Ha=Additive_Determinant_Matrix(f, g, x); Hb=Ha.copy()
        sage: for k in range(Hb.n(0)):
        ....:     Hb=BlockSweep(Hb,k)
        ....:
        sage: Hb[0,0][0,0].canonicalize_radical()
        -(a1*x + x^2 + a0)/((b2 + 1)*x^2 + x^3 + (a1 + b1)*x + a0 + b0)
        sage: (Ha*Hb).simplify_full()
        [[[[-1]], [[0, 0]], [[0, 0, 0]]], [[[0], [0]], [[-1, 0], [0, -1]], [[0, 0, 0], [0, 0, 0]]], [[[0], [0], [0]], [[0, 0], [0, 0], [0, 0]], [[-1, 0, 0], [0, -1, 0], [0, 0, -1]]]]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the size parameter
    sz0=f.degree(t); sz1=g.degree(t)
    # Initialization of the companion matrices
    A=CompanionHM(f,t)-t*HM(2,sz0,'kronecker')
    B=CompanionHM(g,t)-t*HM(2,sz1,'kronecker')
    Tp=HM(sz1,sz0,'zero');Tp[0,sz0-1]=1
    # Initialization of the first row
    M00=HM(1,1,'one'); M01=HM(1,sz0,'zero'); M02=HM(1,sz1,[B[0,j] for j in range(B.n(1))])
    # Performing the row substitution
    for j in range(B.n(1)):
        B[0,j]=HM(2,sz1,'kronecker')[sz1-1,j]
    # Initialization of the  second row
    M10=HM(sz0,1,[HM(2,sz0,'kronecker')[i,0] for i in range(sz0)]); M11=A; M12=HM(sz0,sz1,'zero')
    # Initialization of the last row
    M20=HM(sz1,1,'zero'); M21=Tp.copy(); M22=B
    # Initialization of the hypermatrix
    return HM([[M00,M01,M02],[M10,M11,M12],[M20,M21,M22]])

def GeneralHypermatrixRankOnePartition(B, Hl, Xl):
    """
    Returns the constraints associated with linearizations of the rank one decomposition
    constraints.The polynomial returned come from zero rows. The function takes as inputs
    a hypermatrix and two lists. The input Xl is the list of hypermatrices associated with
    the variables. The input Hl is the list of hypermatrix associated with the decompositions.
    The first list element of Hl is the Hypermatrix to be deocomposed the other elements 
    correspond to the parameters of the decomposition. 

    EXAMPLES:

    ::

        sage: B=HM(2,2,'b'); Hl=[HM(2,2,'c')]; Xl=[HM(2,2,'x'), HM(2,2,'y')] 
        sage: Sln=GeneralHypermatrixRankOnePartition(B, Hl, Xl); Sln
        [-(b01 - c01)*(b10 - c10) + (b00 - c00)*(b11 - c11),
         -(b01 + c01)*(b10 + c10) + (b00 + c00)*(b11 + c11)]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the primitive root of unity
    w=exp(I*2*pi/Xl[0].n(1))
    # Initializing the order
    od=Xl[0].order()
    # Initialization the Kronecker slice selectors
    DltL=GeneralHypermatrixKroneckerDeltaL(od, Xl[0].n(1))
    # Loop initializing the hypermartrix enrtry lists associaed with constraints 
    Lx=[]; Lh=[]
    for t in range(Xl[0].n(1)):
        #Lx=Lx+apply(ProdB,[X for X in Xl]+[DltL[t]]).list()
        Lx=Lx+ProdB(*([X for X in Xl]+[DltL[t]])).list()
        Lh=Lh+((1/Xl[0].n(1))*sum([B]+[Hl[j]*w^(t*(j+1)) for j in rg(len(Hl))])).list()
    # Initialization of the equation
    EqL=[Lx[i]==Lh[i] for i in rg(len(Lx))]
    # Formating the constraints
    LstX=[]
    for x in Xl:
        LstX=LstX+x.list()
    if len(Set(LstX).list())==len(LstX):
        VrbL=LstX
    else:
        VrbL=Set(LstX).list()
    [A,b]=multiplicativeConstraintFormator(EqL, VrbL)
    # Initialization of the vector of variables
    v=Matrix(SR, A.ncols(), 1, VrbL)
    # computing the solutions to the system obtained via Gauss-Jordan elimination
    Sln=multiplicative_linear_solver(A, b, v, v)
    # returning the polynomial conditions eliminating the variables
    return [f.rhs().numerator()-f.rhs().denominator() for f in Sln if f.lhs()==1]

def Form2TotallySymmetricHypermatrix(f, od, Vrbls):
    """
    Procedure for extracting a Hypermatrix from a multivariate
    homogeneous from the the order corresponds to the degree
    of the homogeneous form and the side length is determined
     by the number of variables

    EXAMPLES:

    ::

        sage: sz=2; od=2; X=HM(sz,sz,var_list('x',sz^2)); f=X.det()
        sage: H=Form2TotallySymmetricHypermatrix(f, 2, X.list()); H.printHM()
        [:, :]=
        [   0    0    0  1/2]
        [   0    0 -1/2    0]
        [   0 -1/2    0    0]
        [ 1/2    0    0    0]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the list 
    l=[len(Vrbls) for i in range(od)]
    # Initialization of the hypermatrix 
    inpts=l+['zero']; Rh=HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry=[Integer(mod(i,l[0]))]
        sm=Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm=sm+prod(l[0:k+1])*entry[len(entry)-1]
        Rh[tuple(entry)]=f.diff([Vrbls[v] for v in entry])/factorial(od)
    return Rh

def Form2Hypermatrix(f, Vrbls):
    """
    Procedure for extracting a side length sz Hypermatrix from
    a non homogeneous multivariate from. The order input is od
    the side length is sz. This is inpired by the theory of 
    system of linear equations in terms of side length two 
    hypermatrices. The list Vrbls stores all the variables.
    In the underlying form each variables is used to construct
    a vector of powers of the variable.


    EXAMPLES:

    ::

        sage: sz=2; od=2; A=HM(sz,sz,'a'); B=HM(sz,sz,'b'); X=HM(od,'x').list()
        sage: Hv0=HM(sz, 1, [X[0]^0, X[0]^1]); Hv1=HM(sz, 1, [X[1]^0,X[1]^1])
        sage: f=ProdB(Hv0.transpose(),Hv1,A)[0,0]; f
        a11*x0*x1 + a10*x0 + a01*x1 + a00
        sage: Form2Hypermatrix(f, X).printHM()
        [:, :]=
        [a00 a01]
        [a10 a11]
        sage: g=ProdB(Hv0.transpose(),Hv1,B)[0,0]; g
        b11*x0*x1 + b10*x0 + b01*x1 + b00
        sage: Form2Hypermatrix(g, X).printHM()
        [:, :]=
        [b00 b01]
        [b10 b11]
        sage: sz=2; od=3; A=HM(sz, sz, sz,'a'); X=HM(od, 'x').list()
        sage: Hv0=HM(sz, 1, 1, [1,X[0]]); Hv1=HM(sz, 1, 1, [1,X[1]]); Hv2=HM(sz, 1, 1, [1,X[2]])
        sage: f=ProdB(Hv0.transpose(2), Hv1.transpose(), Hv2, A)[0,0,0]; f
        a111*x0*x1*x2 + a110*x0*x1 + a101*x0*x2 + a011*x1*x2 + a100*x0 + a010*x1 + a001*x2 + a000
        sage: Form2Hypermatrix(f, X).printHM()
        [:, :, 0]=
        [a000 a010]
        [a100 a110]
        <BLANKLINE> 
        [:, :, 1]=
        [a001 a011]
        [a101 a111]
        sage: sz=3; A=HM(sz,sz,'a'); X=var_list('x',2)
        sage: Hv0=HM(sz,1,[X[0]^i for i in range(sz)]); Hv1=HM(sz,1,[X[1]^j for j in range(sz)])
        sage: f=ProdB(Hv0.transpose(),Hv1,A)[0,0]; f
        a22*x0^2*x1^2 + a21*x0^2*x1 + a12*x0*x1^2 + a20*x0^2 + a11*x0*x1 + a02*x1^2 + a10*x0 + a01*x1 + a00
        sage: Form2Hypermatrix(f, X).printHM()
        [:, :]=
        [a00 a01 a02]
        [a10 a11 a12]
        [a20 a21 a22]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the size and order parameters
    od=len(Vrbls); sz=1+max([f.degree(v) for v in Vrbls])
    # Initialization of the list 
    l=[sz for i in range(od)]
    # Initialization of the hypermatrix 
    inpts=l+['zero']; Rh=HM(*inpts)
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry=[Integer(mod(i,l[0]))]
        sm=Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm=sm+prod(l[0:k+1])*entry[len(entry)-1]
        Tmplst=[]
        for z in range(len(entry)):
            Tmplst=Tmplst+[Vrbls[z],entry[z]]
        #print [factorial(Tmplst[b]) for b in range(1,len(Tmplst),2)]
        Rh[tuple(entry)]=f.diff(*Tmplst).subs([v==0 for v in Vrbls])/prod(factorial(Tmplst[b]) for b in range(1,len(Tmplst),2))
    return Rh

def diagonal_coef_gaussian_eliminationHM(Cf1, Vx, Cf2, rs):
    """
    Outputs the row echelon form of the input coefficient third order hypermatrices 
    and the corresponding right hand side. This implementation assumes that the
    inputs are third order hypermatrices whose entries are themselves diagonal
    second order hypermatrices. All entries of Cf1 are matrices of the size 
    m x m and all entries of Cf2 are matrices of the size n x n (this is not 
    checked). Consequently the entries of Vx and rs are all m x n matrices.
    This implementation avoids division and would work on other inputs so
    long as the entries of the coefficien matrices commutes among themeselve.


    EXAMPLES:
    ::

       
        sage: A00=HM(2,HM(2,'a').list(),'diag'); A01=HM(2,HM(2,'b').list(),'diag'); A10=HM(2,HM(2,'c').list(),'diag'); A11=HM(2,HM(2,'d').list(),'diag')
        sage: B00=HM(2,HM(2,'e').list(),'diag'); B01=HM(2,HM(2,'f').list(),'diag'); B10=HM(2,HM(2,'g').list(),'diag'); B11=HM(2,HM(2,'h').list(),'diag')
        sage: Cf1=HM([[[A00, A10], [A01, A11]]]); Cf1.dimensions() # Initialization of the left coefficient matrix
        [1, 2, 2]
        sage: Cf2=HM([[[B00, B01]], [[B10, B11]]]); Cf2.dimensions() # Initialization of the right coefficient matrix
        [2, 1, 2]
        sage: Vx=HM([[[HM(2,2,'x'), HM(2,2,'y')]]])
        sage: rs=HM([[[HM(2,2,'m'), HM(2,2,'n')]]])
        sage: [A, X, B, C]=diagonal_coef_gaussian_eliminationHM(Cf1, Vx, Cf2, rs)
        sage: A[0,0,0].printHM()
        [:, :]=
        [a0  0  0  0]
        [ 0 a1  0  0]
        [ 0  0 a0  0]
        [ 0  0  0 a1]
        sage: A[0,1,0].printHM()
        [:, :]=
        [b0  0  0  0]
        [ 0 b1  0  0]
        [ 0  0 b0  0]
        [ 0  0  0 b1]
        sage: A[0,0,1].printHM()
        [:, :]=
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        sage: A[0,1,1].printHM()
        [:, :]=
        [-b0*c0      0      0      0]
        [     0 -b1*c1      0      0]
        [     0      0  a0*d0      0]
        [     0      0      0  a1*d1]
        sage: B[0,0,0].printHM()
        [:, :]=
        [e0  0  0  0]
        [ 0 e1  0  0]
        [ 0  0 e0  0]
        [ 0  0  0 e1]
        sage: B[0,0,1].printHM()
        [:, :]=
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        sage: B[1,0,0].printHM()
        [:, :]=
        [g0  0  0  0]
        [ 0 g1  0  0]
        [ 0  0 g0  0]
        [ 0  0  0 g1]
        sage: B[1,0,1].printHM()
        [:, :]=
        [f0*g0     0     0     0]
        [    0 f1*g1     0     0]
        [    0     0 e0*h0     0]
        [    0     0     0 e1*h1]
        sage: A00=HM(2,HM(2,'a').list(),'diag'); A01=HM(2,HM(2,'b').list(),'diag'); A10=HM(2,HM(2,'c').list(),'diag')
        sage: A11=HM(2,HM(2,'d').list(),'diag'); A20=HM(2,HM(2,'e').list(),'diag'); A21=HM(2,HM(2,'f').list(),'diag')
        sage: Cf1=HM([[[A00, A10, A20], [A01, A11, A21]]])
        sage: B00=HM(2,HM(2,'g').list(),'diag'); B01=HM(2,HM(2,'h').list(),'diag'); B02=HM(2,HM(2,'i').list(),'diag')
        sage: B10=HM(2,HM(2,'j').list(),'diag'); B11=HM(2,HM(2,'k').list(),'diag'); B12=HM(2,HM(2,'l').list(),'diag')
        sage: Cf2=HM([[[B00, B01, B02]], [[B10, B11, B12]]]) 
        sage: Vx=HM([[[HM(2,2,'x'), HM(2,2,'y')]]])
        sage: rs=HM([[[HM(2,2,'m'), HM(2,2,'n'), HM(2,2,'p')]]])
        sage: [A, X, B, C]=diagonal_coef_gaussian_eliminationHM(Cf1, Vx, Cf2, rs)
        sage: A[0,0,0].printHM()
        [:, :]=
        [a0  0  0  0  0  0  0  0]
        [ 0 a1  0  0  0  0  0  0]
        [ 0  0 a0  0  0  0  0  0]
        [ 0  0  0 a1  0  0  0  0]
        [ 0  0  0  0 a0  0  0  0]
        [ 0  0  0  0  0 a1  0  0]
        [ 0  0  0  0  0  0 a0  0]
        [ 0  0  0  0  0  0  0 a1]
        sage: A00=HM(2,HM(2,'a').list(),'diag'); A01=HM(2,HM(2,'b').list(),'diag'); A10=HM(2,HM(2,'c').list(),'diag')
        sage: A11=HM(2,HM(2,'d').list(),'diag'); A20=HM(2,HM(2,'e').list(),'diag'); A21=HM(2,HM(2,'f').list(),'diag')
        sage: Cf1=HM([[[A00, HM(2,2,'zero'), HM(2,2,'zero')], [A01, A11, HM(2,2,'zero')]]])
        sage: Cf1
        [[[[[a0, 0], [0, a1]], [[0, 0], [0, 0]], [[0, 0], [0, 0]]], [[[b0, 0], [0, b1]], [[d0, 0], [0, d1]], [[0, 0], [0, 0]]]]]
        sage: B00=HM(2,HM(2,'g').list(),'diag'); B01=HM(2,HM(2,'h').list(),'diag'); B02=HM(2,HM(2,'i').list(),'diag')
        sage: B10=HM(2,HM(2,'j').list(),'diag'); B11=HM(2,HM(2,'k').list(),'diag'); B12=HM(2,HM(2,'l').list(),'diag')
        sage: Cf2=HM([[[B00, HM(2,2,'zero'), HM(2,2,'zero')]], [[B10, B11, HM(2,2,'zero')]]])
        sage: Cf2
        [[[[[g0, 0], [0, g1]], [[0, 0], [0, 0]], [[0, 0], [0, 0]]]], [[[[j0, 0], [0, j1]], [[k0, 0], [0, k1]], [[0, 0], [0, 0]]]]]
        sage: Vx=HM([[[HM(2,2,'x'), HM(2,2,'y')]]])
        sage: rs=HM([[[HM(2,2,'m'), HM(2,2,'n'), HM(2,2,'p')]]])
        sage: [A, X, B, C]=diagonal_coef_gaussian_eliminationHM(Cf1, Vx, Cf2, rs)
        sage: A
        [[[[[a0, 0], [0, a1]], [[0, 0], [0, 0]], [[0, 0], [0, 0]]], [[[b0, 0], [0, b1]], [[d0, 0], [0, d1]], [[0, 0], [0, 0]]]]]
        sage: B
        [[[[[g0, 0], [0, g1]], [[0, 0], [0, 0]], [[0, 0], [0, 0]]]], [[[[j0, 0], [0, j1]], [[k0, 0], [0, k1]], [[0, 0], [0, 0]]]]]

 
    AUTHORS:
    - Edinah K. Gnang
    - To Do:
    """
    if Cf1.n(1)==Vx.n(2) and Vx.n(2)==Cf2.n(0) and Cf1.n(2)==Cf2.n(2) and 1==Vx.n(0) and Vx.n(0)==Vx.n(1):
        # Initialization of the variable index
        vindx=0
        # Initializing copies of the input hypermatrices.
        A=Cf1.copy(); X=Vx.copy(); B=Cf2.copy(); C=rs.copy()
        # Initialization of the row and column index
        i=0; j=0
        while i < A.n(2) and j < A.n(1):
            while (HM(1,1,A.n(2)-i,[A[0,j,i0] for i0 in range(i,A.n(2))]).is_zero() and j < A.n(1)-1) or (HM(1,1,B.n(0)-i,[B[j,0,i0] for i0 in range(i,B.n(0))]).is_zero() and j < B.n(0)-1):
                # Incrementing the column index
                j=j+1
            if (HM(1,A.n(1),A.n(2)-i,[A[0,j0,i0] for i0 in range(i,A.n(2)) for j0 in range(A.n(1))]).is_zero()==False) and (HM(B.n(0),1,B.n(2)-i,[B[j0,0,i0] for i0 in range(i,B.n(2)) for j0 in range(B.n(0))]).is_zero()==False) and j < A.n(1):
                while A[0,j,i].is_zero() or B[j,0,i].is_zero():
                    # Initialization of the matrices
                    Ta=HM(A.n(2)-i,A.n(1),[A[0,j0,i0] for j0 in range(A.n(1)) for i0 in range(i,A.n(2))])
                    Tb=HM(B.n(2)-i,B.n(0),[B[j0,0,i0] for j0 in range(B.n(0)) for i0 in range(i,B.n(2))])
                    Tc=HM(C.n(2)-i,1,[C[0,0,i0] for i0 in range(i,C.n(2))])
                    # Inflating the entries of the identity matrix
                    idta=HM(2, A[0,0,0].n(0),'kronecker')
                    Ida=HM(2, Ta.n(0), 'kronecker')
                    idtb=HM(2, B[0,0,0].n(0),'kronecker')
                    Idb=HM(2, Tb.n(0), 'kronecker')
                    for u in range(Ta.n(0)):
                        for v in range(Ta.n(1)):
                            Ida[u,v]=idta*Ida[u,v]
                    for u in range(Tb.n(0)):
                        for v in range(Tb.n(1)):
                            Idb[u,v]=idtb*Idb[u,v]
                    # Initialization of the cyclic shift permutation matrix
                    Pa=sum([HM(Ta.n(0),1,[Ida[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Ida[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                    Pb=sum([HM(Tb.n(0),1,[Idb[i0,k] for i0 in range(Tb.n(0))])*HM(1,Tb.n(0),[Idb[Integer(mod(k+1,Tb.n(0))),j0] for j0 in range(Tb.n(0))]) for k in range(Tb.n(0))])
                    # Performing the shift
                    Ta=Pa*Ta; Tb=Pb*Tb; Tc=Pa*Tc
                    for i0 in range(Ta.n(0)):
                        for j0 in range(Ta.n(1)):
                            A[0,j0,i+i0]=Ta[i0,j0]
                    for i0 in range(Ta.n(0)):
                        for j0 in range(Ta.n(1)):
                            B[j0,0,i+i0]=Tb[i0,j0]
                    for i0 in range(b.n(0)):
                        C[0,0,i+i0]=Tc[i0,0]
                # Main part
                if (A.n(2)-i-1>0) and not ((HM(1,1,A.n(2)-i-1,[A[0,j,i0] for i0 in range(i+1,A.n(2))]).is_zero() and j <= A.n(1)-1) and (HM(1,1,B.n(2)-i-1,[B[j,0,i0] for i0 in range(i+1,B.n(2))]).is_zero() and j <= B.n(0)-1)):
                    # Performing the row operations.
                    cf1a = A[0,j,i]; cf1b = B[j,0,i]
                    # Backing up the variable prior to the inflation
                    Xold=X.copy()
                    # Updating the variables
                    for t in range(X.n(2)):
                        X[0,0,t]=HM(2,2,'kronecker').tensor_product(X[0,0,t])
                    # Updating the right hand side
                    for r in range(i+1,A.n(2)):
                        # Taking care of the zero row
                        if (HM(1,A.n(1),1,[A[0,j0,r] for j0 in range(A.n(1))]).is_zero()) or (HM(B.n(0),1,1,[B[j0,0,r] for j0 in range(B.n(0))]).is_zero()):
                            r=r+1
                        else:
                            # Initialization of the coefficient
                            cf2a=A[0,j,r]; cf2b=B[j,0,r]
                            # Updating the right hand side.
                            n0=C[0,0,0].n(0); n1=C[0,0,0].n(1)
                            U=HM(n0,n1,[var('z'+str(vindx+t)) for t in range(n0*n1)])
                            # Incrementing the free variable index
                            vindx=(n0*n1)+vindx
                            C[0,0,r]=(U-cf2a*C[0,0,i]*cf2b).block_sum(cf1a*C[0,0,r]*cf1b-U)
                            # Updating the constraints
                            for j0 in range(A.n(1)):
                                if (-cf2a*A[0,j0,i]*Xold[0,0,j0]*B[j0,0,i]*cf2b + cf1a*A[0,j0,r]*Xold[0,0,j0]*B[j0,0,r]*cf1b).is_zero():
                                    A[0,j0,r]=HM(2,2,'zero').tensor_product(A[0,j0,r])
                                    B[j0,0,r]=HM(2,2,'zero').tensor_product(B[j0,0,r])
                                else:
                                    A[0,j0,r]=(-cf2a*A[0,j0,i]).block_sum(cf1a*A[0,j0,r])
                                    B[j0,0,r]=( B[j0,0,i]*cf2b).block_sum(B[j0,0,r]*cf1b)
                    for r in range(i+1):
                        # Updating the other entries.
                        for j0 in range(A.n(1)):
                            A[0,j0,r]=HM(2,2,'kronecker').tensor_product(A[0,j0,r])
                        for j0 in range(B.n(0)):
                            B[j0,0,r]=HM(2,2,'kronecker').tensor_product(B[j0,0,r])
                        C[0,0,r]=HM(2,2,'kronecker').tensor_product(C[0,0,r])
            # Incrementing the row and column index
            i=i+1; j=j+1
        return [A, X, B, C]
    else:
        raise ValueError("Incorrect inputs")

def diagonal_coef_gauss_jordan_eliminationHM(Cf1, Vx, Cf2, rs):
    """
    Outputs the reduced row echelon form of the input coefficient third order 
    hypermatrices and the corresponding right hand side. This implementation 
    assumes that the inputs are third order hypermatrices whose entries are 
    themselves diagonal second order hypermatrices (this is not checked).
    All entries of Cf1 are matrices of the size m x m and all entries of Cf2 
    are matrices of the size n x n. Consequently the entries of Vx and rs are 
    all m x n matrices. The implementation avoids division and consequently,
    does normalize the pivots to the units. The procedure would on other inputs so
    long as the entries of the coefficien matrices commutes among themeselve.


    EXAMPLES:
    ::

       
        sage: A00=HM(2,HM(2,'a').list(),'diag'); A01=HM(2,HM(2,'b').list(),'diag'); A10=HM(2,HM(2,'c').list(),'diag'); A11=HM(2,HM(2,'d').list(),'diag')
        sage: B00=HM(2,HM(2,'e').list(),'diag'); B01=HM(2,HM(2,'f').list(),'diag'); B10=HM(2,HM(2,'g').list(),'diag'); B11=HM(2,HM(2,'h').list(),'diag')
        sage: Cf1=HM([[[A00, A10], [A01, A11]]]); Cf1.dimensions() # Initialization of the left coefficient matrix
        [1, 2, 2] 
        sage: Cf2=HM([[[B00, B01]], [[B10, B11]]]); Cf2.dimensions() # Initialization of the right coefficient matrix
        [2, 1, 2]
        sage: Vx=HM([[[HM(2,2,'x'), HM(2,2,'y')]]])
        sage: rs=HM([[[HM(2,2,'m'), HM(2,2,'n')]]])
        sage: [A, X, B, C]=diagonal_coef_gauss_jordan_eliminationHM(Cf1, Vx, Cf2, rs)
        sage: A[0,0,0].printHM()
        [:, :]=
        [-a0*b0*c0         0         0         0]
        [        0 -a1*b1*c1         0         0]
        [        0         0   a0^2*d0         0]
        [        0         0         0   a1^2*d1]
        sage: A[0,1,0].printHM()
        [:, :]=
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        sage: A[0,0,1].printHM()
        [:, :]=
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        sage: A[0,1,1].printHM()
        [:, :]=
        [-b0*c0      0      0      0]
        [     0 -b1*c1      0      0]
        [     0      0  a0*d0      0]
        [     0      0      0  a1*d1]
        sage: B[0,0,0].printHM()
        [:, :]=
        [e0*f0*g0        0        0        0]
        [       0 e1*f1*g1        0        0]
        [       0        0  e0^2*h0        0]
        [       0        0        0  e1^2*h1]
        sage: B[0,0,1].printHM()
        [:, :]=
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        sage: B[1,0,0].printHM()
        [:, :]=
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        sage: B[1,0,1].printHM()
        [:, :]=
        [f0*g0     0     0     0]
        [    0 f1*g1     0     0]
        [    0     0 e0*h0     0]
        [    0     0     0 e1*h1]
        sage: A00=HM(2,HM(2,'a').list(),'diag'); A01=HM(2,HM(2,'b').list(),'diag'); A10=HM(2,HM(2,'c').list(),'diag')
        sage: A11=HM(2,HM(2,'d').list(),'diag'); A20=HM(2,HM(2,'e').list(),'diag'); A21=HM(2,HM(2,'f').list(),'diag')
        sage: Cf1=HM([[[A00, HM(2,2,'zero'), HM(2,2,'zero')], [A01, A11, HM(2,2,'zero')]]])
        sage: B00=HM(2,HM(2,'g').list(),'diag'); B01=HM(2,HM(2,'h').list(),'diag'); B02=HM(2,HM(2,'i').list(),'diag')
        sage: B10=HM(2,HM(2,'j').list(),'diag'); B11=HM(2,HM(2,'k').list(),'diag'); B12=HM(2,HM(2,'l').list(),'diag')
        sage: Cf2=HM([[[B00, HM(2,2,'zero'), HM(2,2,'zero')]], [[B10, B11, HM(2,2,'zero')]]]) 
        sage: Vx=HM([[[HM(2,2,'x'), HM(2,2,'y')]]])
        sage: rs=HM([[[HM(2,2,'m'), HM(2,2,'n'), HM(2,2,'p')]]])
        sage: [A, X, B, C]=diagonal_coef_gauss_jordan_eliminationHM(Cf1, Vx, Cf2, rs)
        sage: A
        [[[[[a0*d0, 0], [0, a1*d1]], [[0, 0], [0, 0]], [[0, 0], [0, 0]]], [[[0, 0], [0, 0]], [[d0, 0], [0, d1]], [[0, 0], [0, 0]]]]]
        sage: C
        [[[[[d0*k0*m00 - b0*j0*n00, d0*k1*m01 - b0*j1*n01], [d1*k0*m10 - b1*j0*n10, d1*k1*m11 - b1*j1*n11]], [[n00, n01], [n10, n11]], [[p00, p01], [p10, p11]]]]]

 
    AUTHORS:
    - Edinah K. Gnang
    - To Do:
    """
    [A, X, B, C]=diagonal_coef_gaussian_eliminationHM(Cf1, Vx, Cf2, rs)
    # Initialization of the row and column index
    i=A.n(2)-1; j=0
    while i > 0 or j > 0:
        if HM(1,A.n(1),1,[A[0,j0,i] for j0 in range(A.n(1))]).is_zero() or HM(B.n(0),1,1,[B[j0,0,i] for j0 in range(B.n(0))]).is_zero():
            # decrementing the row index and initializing the column index
            i=i-1; j=0
        else :
            while A[0,j,i].is_zero() or B[j,0,i].is_zero():
                # Incrementing the column index
                j = j + 1
            # Performing row operations
            cf1a=A[0,j,i]; cf1b=B[j,0,i]
            for r in range(i-1,-1,-1):
                cf2a=A[0,j,r]; cf2b=B[j,0,r]
                # Updating the right hand side
                C[0,0,r]=-cf2a*C[0,0,i]*cf2b + cf1a*C[0,0,r]*cf1b
                # Updating the coefficients
                for j0 in range(A.n(1)):
                    A[0,j0,r]=-cf2a*A[0,j0,i] + cf1a*A[0,j0,r]
                    B[j0,0,r]=-B[j0,0,i]*cf2b + B[j0,0,r]*cf1b
            i=i-1; j=0
    return [A, X, B, C]

def general_gaussian_eliminationHM(Cf1, Vx, Cf2, rs):
    """
    Outputs the row echelon form of the input coefficient third order hypermatrices 
    and the corresponding right hand side. This implementation assumes that the
    inputs are third order hypermatrices whose entries are themselves second order
    hypermatrices (square matrices to be precise). All entries of Cf1 are matrices
    of the size m x m and all entries of Cf2 are matrices of the size n x n.
    Consequently the entries of Vx and rs are all m x n matrices.


    EXAMPLES:
    ::

       
        sage: A00=HM(2,HM(2,'a').list(),'diag'); A01=HM(2,HM(2,'b').list(),'diag'); A10=HM(2,HM(2,'c').list(),'diag'); A11=HM(2,HM(2,'d').list(),'diag')
        sage: B00=HM(2,HM(2,'e').list(),'diag'); B01=HM(2,HM(2,'f').list(),'diag'); B10=HM(2,HM(2,'g').list(),'diag'); B11=HM(2,HM(2,'h').list(),'diag')
        sage: Cf1=HM([[[A00, A10], [A01, A11]]]); Cf1.dimensions() # Initialization of the left coefficient matrix
        [1, 2, 2] 
        sage: Cf2=HM([[[B00, B01]], [[B10, B11]]]); Cf2.dimensions() # Initialization of the right coefficient matrix
        [2, 1, 2]
        sage: Vx=HM([[[HM(2,2,'x'), HM(2,2,'y')]]])
        sage: rs=HM([[[HM(2,2,'m'), HM(2,2,'n')]]])
        sage: [A, X, B, C]=general_gaussian_eliminationHM(Cf1, Vx, Cf2, rs)
        sage: A[0,0,0].printHM()
        [:, :]=
        [a0  0  0  0]
        [ 0 a1  0  0]
        [ 0  0 a0  0]
        [ 0  0  0 a1]
        sage: A[0,1,0].printHM()
        [:, :]=
        [b0  0  0  0]
        [ 0 b1  0  0]
        [ 0  0 b0  0]
        [ 0  0  0 b1]
        sage: A[0,0,1].printHM()
        [:, :]=
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        sage: A[0,1,1].printHM()
        [:, :]=
        [-b0*c0/a0         0         0         0]
        [        0 -b1*c1/a1         0         0]
        [        0         0        d0         0]
        [        0         0         0        d1]
        sage: B[0,0,0].printHM()
        [:, :]=
        [e0  0  0  0]
        [ 0 e1  0  0]
        [ 0  0 e0  0]
        [ 0  0  0 e1]
        sage: B[0,0,1].printHM()
        [:, :]=
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        sage: B[1,0,0].printHM()
        [:, :]=
        [g0  0  0  0]
        [ 0 g1  0  0]
        [ 0  0 g0  0]
        [ 0  0  0 g1]
        sage: B[1,0,1].printHM()
        [:, :]=
        [f0*g0/e0        0        0        0]
        [       0 f1*g1/e1        0        0]
        [       0        0       h0        0]
        [       0        0        0       h1] 
        sage: A00=HM(2,HM(2,'a').list(),'diag'); A01=HM(2,HM(2,'b').list(),'diag'); A10=HM(2,HM(2,'c').list(),'diag')
        sage: A11=HM(2,HM(2,'d').list(),'diag'); A20=HM(2,HM(2,'e').list(),'diag'); A21=HM(2,HM(2,'f').list(),'diag')
        sage: Cf1=HM([[[A00, A10, A20], [A01, A11, A21]]])
        sage: B00=HM(2,HM(2,'g').list(),'diag'); B01=HM(2,HM(2,'h').list(),'diag'); B02=HM(2,HM(2,'i').list(),'diag')
        sage: B10=HM(2,HM(2,'j').list(),'diag'); B11=HM(2,HM(2,'k').list(),'diag'); B12=HM(2,HM(2,'l').list(),'diag')
        sage: Cf2=HM([[[B00, B01, B02]], [[B10, B11, B12]]]) 
        sage: Vx=HM([[[HM(2,2,'x'), HM(2,2,'y')]]])
        sage: rs=HM([[[HM(2,2,'m'), HM(2,2,'n'), HM(2,2,'p')]]])
        sage: [A, X, B, C]=general_gaussian_eliminationHM(Cf1, Vx, Cf2, rs)
        sage: A[0,0,0].printHM()
        [:, :]=
        [a0  0  0  0  0  0  0  0]
        [ 0 a1  0  0  0  0  0  0]
        [ 0  0 a0  0  0  0  0  0]
        [ 0  0  0 a1  0  0  0  0]
        [ 0  0  0  0 a0  0  0  0]
        [ 0  0  0  0  0 a1  0  0]
        [ 0  0  0  0  0  0 a0  0]
        [ 0  0  0  0  0  0  0 a1]
        sage: A00=HM(2,2,'a'); A01=HM(2,2,'b'); A10=HM(2,2,'c')
        sage: A11=HM(2,2,'d'); A20=HM(2,2,'e'); A21=HM(2,2,'f')
        sage: Cf1=HM([[[A00, A10, A20], [A01, A11, A21]]])
        sage: B00=HM(2,2,'g'); B01=HM(2,2,'h'); B02=HM(2,2,'i')
        sage: B10=HM(2,2,'j'); B11=HM(2,2,'k'); B12=HM(2,2,'l')
        sage: Cf2=HM([[[B00, B01, B02]], [[B10, B11, B12]]]) 
        sage: Vx=HM([[[HM(2,2,'x'), HM(2,2,'y')]]])
        sage: rs=HM([[[HM(2,2,'m'), HM(2,2,'n'), HM(2,2,'p')]]])
        sage: [A, X, B, C]=general_gaussian_eliminationHM(Cf1, Vx, Cf2, rs)
        sage: A[0,0,0].printHM()
        [:, :]=
        [a00 a01   0   0   0   0   0   0]
        [a10 a11   0   0   0   0   0   0]
        [  0   0 a00 a01   0   0   0   0]
        [  0   0 a10 a11   0   0   0   0]
        [  0   0   0   0 a00 a01   0   0]
        [  0   0   0   0 a10 a11   0   0]
        [  0   0   0   0   0   0 a00 a01]
        [  0   0   0   0   0   0 a10 a11]
        sage: A00=HM(2,HM(2,'a').list(),'diag'); A01=HM(2,HM(2,'b').list(),'diag'); A10=HM(2,HM(2,'c').list(),'diag')
        sage: A11=HM(2,HM(2,'d').list(),'diag'); A20=HM(2,HM(2,'e').list(),'diag'); A21=HM(2,HM(2,'f').list(),'diag')
        sage: Cf1=HM([[[A00, HM(2,2,'zero'), HM(2,2,'zero')], [A01, A11, HM(2,2,'zero')]]])
        sage: B00=HM(2,HM(2,'g').list(),'diag'); B01=HM(2,HM(2,'h').list(),'diag'); B02=HM(2,HM(2,'i').list(),'diag')
        sage: B10=HM(2,HM(2,'j').list(),'diag'); B11=HM(2,HM(2,'k').list(),'diag'); B12=HM(2,HM(2,'l').list(),'diag')
        sage: Cf2=HM([[[B00, HM(2,2,'zero'), HM(2,2,'zero')]], [[B10, B11, HM(2,2,'zero')]]]) 
        sage: Vx=HM([[[HM(2,2,'x'), HM(2,2,'y')]]])
        sage: rs=HM([[[HM(2,2,'m'), HM(2,2,'n'), HM(2,2,'p')]]])
        sage: [A, X, B, C]=general_gaussian_eliminationHM(Cf1, Vx, Cf2, rs)
        sage: A[0,0,0].printHM()
        [:, :]=
        [a0  0]
        [ 0 a1]

 
    AUTHORS:
    - Edinah K. Gnang
    - To Do:
    """
    if Cf1.n(1)==Vx.n(2) and Vx.n(2)==Cf2.n(0) and Cf1.n(2)==Cf2.n(2) and 1==Vx.n(0) and Vx.n(0)==Vx.n(1):
        # Initialization of the variable index
        vindx=0
        # Initializing copies of the input hypermatrices.
        A=Cf1.copy(); X=Vx.copy(); B=Cf2.copy(); C=rs.copy()
        # Initialization of the row and column index
        i=0; j=0
        while i < A.n(2) and j < A.n(1):
            while (HM(1,1,A.n(2)-i,[A[0,j,i0] for i0 in range(i,A.n(2))]).is_zero() and j < A.n(1)-1) or (HM(1,1,B.n(2)-i,[B[j,0,i0] for i0 in range(i,B.n(2))]).is_zero() and j < B.n(2)-1):
                # Incrementing the column index
                j=j+1
            if (HM(1,A.n(1),A.n(2)-i,[A[0,j0,i0] for i0 in range(i,A.n(2)) for j0 in range(A.n(1))]).is_zero()==False) and (HM(B.n(0),1,B.n(2)-i,[B[j0,0,i0] for i0 in range(i,B.n(2)) for j0 in range(B.n(0))]).is_zero()==False) and j < A.n(1):
                while A[0,j,i].is_zero() or B[j,0,i].is_zero():
                    # Initialization of the matrices
                    Ta=HM(A.n(2)-i,A.n(1),[A[0,j0,i0] for j0 in range(A.n(1)) for i0 in range(i,A.n(2))])
                    Tb=HM(B.n(2)-i,B.n(0),[B[j0,0,i0] for j0 in range(B.n(0)) for i0 in range(i,B.n(2))])
                    Tc=HM(C.n(2)-i,1,[C[0,0,i0] for i0 in range(i,C.n(2))])
                    # Inflating the entries of the identity matrix
                    idta=HM(2, A[0,0,0].n(0),'kronecker')
                    Ida=HM(2, Ta.n(0), 'kronecker')
                    idtb=HM(2, B[0,0,0].n(0),'kronecker')
                    Idb=HM(2, Tb.n(0), 'kronecker')
                    for u in range(Ta.n(0)):
                        for v in range(Ta.n(1)):
                            Ida[u,v]=idta*Ida[u,v]
                    for u in range(Tb.n(0)):
                        for v in range(Tb.n(1)):
                            Idb[u,v]=idtb*Idb[u,v]
                    # Initialization of the cyclic shift permutation matrix
                    Pa=sum([HM(Ta.n(0),1,[Ida[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Ida[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                    Pb=sum([HM(Tb.n(0),1,[Idb[i0,k] for i0 in range(Tb.n(0))])*HM(1,Tb.n(0),[Idb[Integer(mod(k+1,Tb.n(0))),j0] for j0 in range(Tb.n(0))]) for k in range(Tb.n(0))])
                    # Performing the shift
                    Ta=Pa*Ta; Tb=Pb*Tb; Tc=Pa*Tc
                    for i0 in range(Ta.n(0)):
                        for j0 in range(Ta.n(1)):
                            A[0,j0,i+i0]=Ta[i0,j0]
                    for i0 in range(Ta.n(0)):
                        for j0 in range(Ta.n(1)):
                            B[j0,0,i+i0]=Tb[i0,j0]
                    for i0 in range(b.n(0)):
                        C[0,0,i+i0]=Tc[i0,0]
                # Here we mean business
                if (A.n(2)-i-1> 0) and not ((HM(1,1,A.n(2)-i-1,[A[0,j,i0] for i0 in range(i+1,A.n(2))]).is_zero() and j <= A.n(1)-1) and (HM(1,1,B.n(2)-i-1,[B[j,0,i0] for i0 in range(i+1,B.n(2))]).is_zero() and j <= B.n(0)-1)):
                    # Performing the row operations.
                    cf1a = A[0,j,i]; cf1b = B[j,0,i]
                    # Backing up the variable prior to the inflation
                    Xold=X.copy()
                    # Updating the variables
                    for t in range(X.n(2)):
                        X[0,0,t]=HM(2,2,'kronecker').tensor_product(X[0,0,t])
                    # Updating the right hand side
                    for r in range(i+1,A.n(2)):
                        # Taking care of the zero row
                        if (HM(1,A.n(1),1,[A[0,j0,r] for j0 in range(A.n(1))]).is_zero()) or (HM(B.n(0),1,1,[B[j0,0,r] for j0 in range(B.n(0))]).is_zero()):
                            r=r+1
                        else:
                            # Initialization of the coefficient
                            cf2a=A[0,j,r]; cf2b=B[j,0,r]
                            # Updating the right hand side.
                            n0=C[0,0,0].n(0); n1=C[0,0,0].n(1)
                            U=HM(n0,n1,[var('z'+str(vindx+t)) for t in range(n0*n1)])
                            # Incrementing the free variable index
                            vindx=(n0*n1)+vindx
                            C[0,0,r]=(U-(cf2a*cf1a^(-1))*C[0,0,i]*(cf1b^(-1)*cf2b)).block_sum(C[0,0,r]-U)
                            # Updating the constraints
                            for j0 in range(A.n(1)):
                                if (-cf2a*A[0,j0,i]*Xold[0,0,j0]*B[j0,0,i]*cf2b + cf1a*A[0,j0,r]*Xold[0,0,j0]*B[j0,0,r]*cf1b).is_zero():
                                    A[0,j0,r]=HM(2,2,'zero').tensor_product(A[0,j0,r])
                                    B[j0,0,r]=HM(2,2,'zero').tensor_product(B[j0,0,r])
                                else:
                                    A[0,j0,r]=(-(cf2a*cf1a^(-1))*A[0,j0,i]).block_sum(A[0,j0,r])
                                    B[j0,0,r]=(B[j0,0,i]*(cf1b^(-1)*cf2b)).block_sum(B[j0,0,r])
                    for r in range(i+1):
                        # Updating the other entries.
                        for j0 in range(A.n(1)):
                            A[0,j0,r]=HM(2,2,'kronecker').tensor_product(A[0,j0,r])
                            B[j0,0,r]=HM(2,2,'kronecker').tensor_product(B[j0,0,r])
                        C[0,0,r]=HM(2,2,'kronecker').tensor_product(C[0,0,r])
            # Incrementing the row and column index
            i=i+1; j=j+1
        return [A, X, B, C]
    else:
        raise ValueError("Incorrect inputs")

def general_gauss_jordan_eliminationHM(Cf1, Vx, Cf2, rs):
    """
    Outputs the reduced row echelon form of the input coefficient third order hypermatrices 
    and the corresponding right hand side. This implementation assumes that the
    inputs are third order hypermatrices whose entries are themselves second order
    hypermatrices (square matrices to be precise). All entries of Cf1 are matrices
    of the size m x m and all entries of Cf2 are matrices of the size n x n.
    Consequently the entries of Vx and rs are all m x n matrices.


    EXAMPLES:
    ::


        sage: A00=HM(2,HM(2,'a').list(),'diag'); A01=HM(2,HM(2,'b').list(),'diag'); A10=HM(2,HM(2,'c').list(),'diag'); A11=HM(2,HM(2,'d').list(),'diag')
        sage: B00=HM(2,HM(2,'e').list(),'diag'); B01=HM(2,HM(2,'f').list(),'diag'); B10=HM(2,HM(2,'g').list(),'diag'); B11=HM(2,HM(2,'h').list(),'diag')
        sage: Cf1=HM([[[A00, A10], [A01, A11]]]); Cf1.dimensions() # Initialization of the left coefficient matrix
        [1, 2, 2] 
        sage: Cf2=HM([[[B00, B01]], [[B10, B11]]]); Cf2.dimensions() # Initialization of the right coefficient matrix
        [2, 1, 2]
        sage: Vx=HM([[[HM(2,2,'x'), HM(2,2,'y')]]])
        sage: rs=HM([[[HM(2,2,'m'), HM(2,2,'n')]]])
        sage: [A, X, B, C]=general_gauss_jordan_eliminationHM(Cf1, Vx, Cf2, rs)
        sage: A[0,0,0].printHM()
        [:, :]=
        [a0  0  0  0]
        [ 0 a1  0  0]
        [ 0  0 a0  0]
        [ 0  0  0 a1]
        sage: A[0,1,0].printHM()
        [:, :]=
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        sage: A[0,0,1].printHM()
        [:, :]=
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        sage: A[0,1,1].printHM()
        [:, :]=
        [-b0*c0/a0         0         0         0]
        [        0 -b1*c1/a1         0         0]
        [        0         0        d0         0]
        [        0         0         0        d1]
        sage: B[0,0,0].printHM()
        [:, :]=
        [e0  0  0  0]
        [ 0 e1  0  0]
        [ 0  0 e0  0]
        [ 0  0  0 e1]
        sage: B[0,0,1].printHM()
        [:, :]=
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        sage: B[1,0,0].printHM()
        [:, :]=
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        sage: B[1,0,1].printHM()
        [:, :]=
        [f0*g0/e0        0        0        0]
        [       0 f1*g1/e1        0        0]
        [       0        0       h0        0]
        [       0        0        0       h1]
        sage: A00=HM(2,HM(2,'a').list(),'diag'); A01=HM(2,HM(2,'b').list(),'diag'); A10=HM(2,HM(2,'c').list(),'diag')
        sage: A11=HM(2,HM(2,'d').list(),'diag'); A20=HM(2,HM(2,'e').list(),'diag'); A21=HM(2,HM(2,'f').list(),'diag')
        sage: Cf1=HM([[[A00, HM(2,2,'zero'), HM(2,2,'zero')], [A01, A11, HM(2,2,'zero')]]])
        sage: B00=HM(2,HM(2,'g').list(),'diag'); B01=HM(2,HM(2,'h').list(),'diag'); B02=HM(2,HM(2,'i').list(),'diag')
        sage: B10=HM(2,HM(2,'j').list(),'diag'); B11=HM(2,HM(2,'k').list(),'diag'); B12=HM(2,HM(2,'l').list(),'diag')
        sage: Cf2=HM([[[B00, HM(2,2,'zero'), HM(2,2,'zero')]], [[B10, B11, HM(2,2,'zero')]]]) 
        sage: Vx=HM([[[HM(2,2,'x'), HM(2,2,'y')]]])
        sage: rs=HM([[[HM(2,2,'m'), HM(2,2,'n'), HM(2,2,'p')]]])
        sage: [A, X, B, C]=general_gauss_jordan_eliminationHM(Cf1, Vx, Cf2, rs)
        sage: A
        [[[[[a0, 0], [0, a1]], [[0, 0], [0, 0]], [[0, 0], [0, 0]]], [[[0, 0], [0, 0]], [[d0, 0], [0, d1]], [[0, 0], [0, 0]]]]]        
        sage: C
        [[[[[m00 - b0*j0*n00/(d0*k0), m01 - b0*j1*n01/(d0*k1)], [m10 - b1*j0*n10/(d1*k0), m11 - b1*j1*n11/(d1*k1)]], [[n00, n01], [n10, n11]], [[p00, p01], [p10, p11]]]]]

 
    AUTHORS:
    - Edinah K. Gnang
    - To Do:
    """
    [A, X, B, C]=general_gaussian_eliminationHM(Cf1, Vx, Cf2, rs)
    # Initialization of the row and column index
    i=A.n(2)-1; j=0
    while i > 0 or j > 0:
        if HM(1,A.n(1),1,[A[0,j0,i] for j0 in range(A.n(1))]).is_zero() or HM(B.n(0),1,1,[B[j0,0,i] for j0 in range(B.n(0))]).is_zero():
            # decrementing the row index and initializing the column index
            i=i-1; j=0
        else :
            while A[0,j,i].is_zero() or B[j,0,i].is_zero():
                # Incrementing the column index
                j = j + 1
            # Performing row operations
            cf1a=A[0,j,i]; cf1b=B[j,0,i]
            for r in range(i-1,-1,-1):
                cf2a=A[0,j,r]; cf2b=B[j,0,r]
                # Updating the right hand side
                C[0,0,r]=-(cf2a*cf1a^(-1))*C[0,0,i]*(cf1b^(-1)*cf2b) + C[0,0,r]
                # Updating the coefficients
                for j0 in range(A.n(1)):
                    A[0,j0,r]=-(cf2a*cf1a^(-1))*A[0,j0,i] + A[0,j0,r]
                for j0 in range(B.n(0)):
                    B[j0,0,r]=-B[j0,0,i]*(cf1b^(-1)*cf2b) + B[j0,0,r]
            i=i-1; j=0
    return [A, X, B, C]

def ThirdOrderDepthCyclicShift(A, s=1):
    """ 
    This function performs a cyclic shift to the order of
    the depth slices of an input third order hypermatrix A.

    EXAMPLES:
    ::

  
        sage: A = HM(3,3,3,'a'); ThirdOrderDepthCyclicShift(A).printHM()
        [:, :, 0]= 
        [a002 a012 a022]
        [a102 a112 a122]
        [a202 a212 a222]
        <BLANKLINE>
        [:, :, 1]=
        [a000 a010 a020]
        [a100 a110 a120]
        [a200 a210 a220]
        <BLANKLINE>
        [:, :, 2]=
        [a001 a011 a021]
        [a101 a111 a121]
        [a201 a211 a221]
        <BLANKLINE>


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the output hypermatrix to be filled
    #B = apply(HM, A.dimensions()+['zero'])
    B = HM(*(A.dimensions()+['zero']))
    for i in range(A.n(0)):
        for j in range(A.n(1)):
            for k in range(A.n(2)):
                B[i,j,k]=A[i,j,Integer(mod(k-s,A.n(2)))]
    return B

def hadamard_gaussian_eliminationHM(Cf, rs):
    """
    Outputs the row echelon form of the input second order hypermatrix and the right hand side.
    does not normalize the rows to ensure that the first non zero entry of non zero rows = 1
    This implementation tacitly assumes that the entries commute. As a result this implementation
    is NOT skew field friendly. Furthermore the implementation assumes that the inputs are
    second order hypermatrices whose entries are themselves second order hypermatrices.
 

    EXAMPLES:
 
    ::

        sage: Ha=HM(2,2,[HM(2,2,'a'), HM(2,2,'c'), HM(2,2,'b'), HM(2,2,'d')])
        sage: Hb=HM(2,1,[HM(2,2,'f'), HM(2,2,'g')])
        sage: [A,b]=hadamard_gaussian_eliminationHM(Ha,Hb)
        sage: A
        [[[[a00, a01], [a10, a11]], [[b00, b01], [b10, b11]]], [[[0, 0], [0, 0]], [[b00*c00 - a00*d00, b01*c01 - a01*d01], [b10*c10 - a10*d10, b11*c11 - a11*d11]]]]
        sage: b
        [[[[f00, f01], [f10, f11]]], [[[c00*f00 - a00*g00, c01*f01 - a01*g01], [c10*f10 - a10*g10, c11*f11 - a11*g11]]]]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing a copy of the input second order hypermatrices.
    A=Cf.copy(); b=rs.copy()
    # Initialization of the row and column index
    i=0; j=0
    while i < A.n(0) and j < A.n(1):
        while HM(A.n(0)-i, 1, [A[i0,j] for i0 in range(i,A.n(0))]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        if HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))]).is_zero()==False:
            while A[i,j].is_zero(): 
                Ta=HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))])
                Tb=HM(b.n(0)-i, b.n(1), [b[i0,j0] for j0 in range(b.n(1)) for i0 in range(i,b.n(0))])
                # Initializing the cyclic shift permutation matrix
                Id=HM(2, Ta.n(0), 'kronecker')
                P=sum([HM(Ta.n(0),1,[Id[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Id[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                Ta=P*Ta; Tb=P*Tb
                for i0 in range(Ta.n(0)):
                    for j0 in range(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
                for i0 in range(Tb.n(0)):
                    for j0 in range(Tb.n(1)):
                        b[i+i0,j0]=Tb[i0,j0]
            # Performing the row operations.
            cf1=A[i,j]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if HM(1,A.n(1),[A[r,j0] for j0 in range(A.n(1))]).is_zero():
                    r=r+1
                else:
                    # Initialization of the coefficient
                    cf2=A[r,j]
                    for j0 in range(b.n(1)):
                        b[r,j0]=cf2.elementwise_product(b[i,j0])-cf1.elementwise_product(b[r,j0])
                    for j0 in range(A.n(1)):
                        A[r,j0]=cf2.elementwise_product(A[i,j0])-cf1.elementwise_product(A[r,j0])
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return [A,b]

def hadamard_gauss_jordan_eliminationHM(Cf,rs):
    """
    Outputs the reduced row echelon form of the input matrix and the right hand side.
    This implementation assumes that the input entries commute and is therefore NOT
    skew field friendly. Furthermore the implementation assumes that the inputs are
    second order hypermatrices whose entries are themselves second order hypermatrices.


    EXAMPLES:
 
    ::

        sage: Ha=HM(2,2,[HM(2,2,'a'), HM(2,2,'c'), HM(2,2,'b'), HM(2,2,'d')])
        sage: Hb=HM(2,1,[HM(2,2,'f'), HM(2,2,'g')])
        sage: [A,b]=hadamard_gauss_jordan_eliminationHM(Ha,Hb)
        sage: A
        [[[[-(b00*c00 - a00*d00)*a00, -(b01*c01 - a01*d01)*a01], [-(b10*c10 - a10*d10)*a10, -(b11*c11 - a11*d11)*a11]], [[0, 0], [0, 0]]], [[[0, 0], [0, 0]], [[b00*c00 - a00*d00, b01*c01 - a01*d01], [b10*c10 - a10*d10, b11*c11 - a11*d11]]]]
        sage: b
        [[[[(c00*f00 - a00*g00)*b00 - (b00*c00 - a00*d00)*f00, (c01*f01 - a01*g01)*b01 - (b01*c01 - a01*d01)*f01], [(c10*f10 - a10*g10)*b10 - (b10*c10 - a10*d10)*f10, (c11*f11 - a11*g11)*b11 - (b11*c11 - a11*d11)*f11]]], [[[c00*f00 - a00*g00, c01*f01 - a01*g01], [c10*f10 - a10*g10, c11*f11 - a11*g11]]]]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    [A, b] = hadamard_gaussian_eliminationHM(Cf,rs)
    # Initialization of the row and column index
    i=A.nrows()-1; j=0
    while i>0 or j>0:
        #if (A[i,:]).is_zero():
        if HM(1,A.n(1),[A[i,j0] for j0 in range(A.n(1))]).is_zero():
            # decrementing the row index and initializing the column index
            i=i-1; j=0
        else :
            while (A[i,j]).is_zero():
                # Incrementing the column index
                j = j + 1
            # performing row operations
            cf1=A[i,j]
            for r in range(i-1,-1,-1):
                #b[r,:] = -A[r,j]*b[i,:]+b[r,:]
                cf2=A[r,j]
                for j0 in range(b.n(1)):
                    b[r,j0]=cf2.elementwise_product(b[i,j0])-cf1.elementwise_product(b[r,j0])
                #A[r,:] = -A[r,j]*A[i,:]+A[r,:]
                for j0 in range(A.n(1)):
                    A[r,j0]=cf2.elementwise_product(A[i,j0])-cf1.elementwise_product(A[r,j0])
            i=i-1; j=0
    return [A,b]

def list_elementwise_product(L):
    """
    Procedure for computing Hypermatrix Hadamard products of list elementwise.

    EXAMPLES:

    ::

        sage: A=HM(2,2,'a'); B=HM(2,2,'b'); C=HM(2,2,'c')
        sage: list_elementwise_product([A,B,C]).printHM()
        [:, :]=
        [a00*b00*c00 a01*b01*c01]
        [a10*b10*c10 a11*b11*c11]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the hypermatrix
    Tmp=L[0]
    for h in L[1:]:
        Tmp=Tmp.elementwise_product(h)  
    return Tmp

def GeneralHypermatrixProduct_with_elementwise_product(*args):
    """
    Outputs a list of lists associated with the general
    Bhattacharya-Mesner product of the input hypermatrices
    whose entries are themselves hypermatrices of the same
    size. The BM product is performed as usual with the 
    exception that we perform elementwise product performing
    elementwise product of the entries


    EXAMPLES:

    ::

        sage: A=HM(2,2,[HM(2,2,'a'),HM(2,2,'c'),HM(2,2,'b'),HM(2,2,'d')])
        sage: B=HM(2,2,[HM(2,2,'e'),HM(2,2,'f'),HM(2,2,'g'),HM(2,2,'h')])
        sage: GeneralHypermatrixProduct_with_elementwise_product(A,B)[0,0].printHM()
        [:, :]=
        [a00*e00 + b00*f00 a01*e01 + b01*f01]
        [a10*e10 + b10*f10 a11*e11 + b11*f11]
        sage: GeneralHypermatrixProduct_with_elementwise_product(A,B)[0,1].printHM()
        [:, :]=
        [a00*g00 + b00*h00 a01*g01 + b01*h01]
        [a10*g10 + b10*h10 a11*g11 + b11*h11]
        sage: GeneralHypermatrixProduct_with_elementwise_product(A,B)[1,0].printHM()
        [:, :]=
        [c00*e00 + d00*f00 c01*e01 + d01*f01]
        [c10*e10 + d10*f10 c11*e11 + d11*f11]
        sage: GeneralHypermatrixProduct_with_elementwise_product(A,B)[1,1].printHM()
        [:, :]=
        [c00*g00 + d00*h00 c01*g01 + d01*h01]
        [c10*g10 + d10*h10 c11*g11 + d11*h11]


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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        # computing the Hypermatrix product
        if len(args)<2:
            raise ValueError("The number of operands must be >= 2")
        elif len(args) >= 2:
            Rh[tuple(entry)]=sum([list_elementwise_product([args[s][tuple(entry[0:Integer(mod(s+1,len(args)))]+[t]+entry[Integer(mod(s+2,len(args))):])] for s in range(len(args)-2)]+[args[len(args)-2][tuple(entry[0:len(args)-1]+[t])]]+[args[len(args)-1][tuple([t]+entry[1:])]]) for t in range((args[0]).n(1))])
    return Rh

def hadamard_linear_solverHM(A,b,x,v):
    """
    Outputs the Reduced Row Echelon Form of the input matrix and the right hand side.
    where A denotes the input matrix, b denotes the right-hand side vector, x denotes
    the variable vector coming from the original system of equations, and v denotes 
    the free variable vector.

    EXAMPLES:
 
    ::

        sage: sz=2; Eq=[var('x'+str(i))+var('x'+str(sz+j))==var('a'+str(i)+str(j)) for i in range(sz) for j in range(sz)]
        sage: [A,b]=ConstraintFormatorHM(Eq,[var('x'+str(i)) for i in range(2*sz)])
        sage: Mx=HM(A.ncols(),1,[HM(1,1,[var('x'+str(i))]) for i in range(A.ncols())])
        sage: Mv=HM(A.ncols(),1,[HM(1,1,[var('t'+str(i))]) for i in range(A.ncols())])
        sage: Ha=HM(A.n(0),A.n(1),'zero'); Hb=HM(b.n(0),b.n(1),'zero')
        sage: for i in range(A.n(0)):
        ....:     Hb[i,0]=HM(1,1,[b[i,0]])
        ....:     for j in range(A.n(1)):
        ....:         Ha[i,j]=HM(1,1,[A[i,j]])
        ....:
        sage: Sln=hadamard_linear_solverHM(Ha,Hb,Mx,Mv); Sln
        [[[[x0]], [[a00 - a10 + a11 - t3]]],
         [[[x1]], [[a11 - t3]]],
         [[[x2]], [[a10 - a11 + t3]]],
         [[[0]], [[-a00 + a01 + a10 - a11]]]] 
        sage: A=HM(2,2,'a'); b=HM(2,1,HM(2,'b').list()); X=HM(2,1,HM(2,'x').list())
        sage: for i in range(2):
        ....:     b[i,0]=HM(1,1,[b[i,0]]); X[i,0]=HM(1,1,[X[i,0]])
        ....:     for j in range(2):
        ....:         A[i,j]=HM(1,1,[A[i,j]])
        ....:
        sage: Sln=hadamard_linear_solverHM(A,b,X,X); Sln
        [[[[-(a01*a10 - a00*a11)*a00*x0]],
          [[(a10*b0 - a00*b1)*a01 - (a01*a10 - a00*a11)*b0]]],
         [[[(a01*a10 - a00*a11)*x1]], [[a10*b0 - a00*b1]]]]

        
    AUTHORS:
    - Initial implementation by Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the reduced echelon form.
    [Ap,bp]=hadamard_gauss_jordan_eliminationHM(A,b)
    Id1=HM(2,Ap.n(0),'kronecker')
    for i in range(Id1.n(0)):
        for j in range(Id1.n(1)):
            Id1[i,j]=HM(2,Ap[0,0].n(0),'kronecker')*Id1[i,j]
    Id2=HM(2,Ap.n(1),'kronecker')
    for i in range(Id2.n(0)):
        for j in range(Id2.n(1)):
            Id2[i,j]=HM(2,Ap[0,0].n(1),'kronecker')*Id2[i,j]
    # Obtainin the list of pivot variables.
    Pm=HM(Ap.n(0),Ap.n(1),[HM(Ap[0,0].n(0),Ap[0,0].n(1),'zero') for cnt in range(Ap.n(0)*Ap.n(1))])
    for i in range(Ap.n(0)):
        if not HM(1,Ap.n(1),[Ap[i,u] for u in range(Ap.n(1))]).is_zero():
            for j in range(Ap.n(1)):
                #if Ap[i,j].is_unit():
                if not Ap[i,j].is_zero():
                    break
            #Pm=Pm+HM(Id1.n(0),1,[Id1[s,i] for s in range(Id1.n(0))])*HM(1,Id2.n(1),[Id2[j,t]*Ap[i,j] for t in range(Id2.n(1))])
            Pm[i,j]=Ap[i,j].copy()
    # Expressing the solutions
    tp1=GeneralHypermatrixProduct_with_elementwise_product(Pm,x)
    tp2=bp-GeneralHypermatrixProduct_with_elementwise_product((Ap-Pm),v)
    return [[tp1[i,0],tp2[i,0]] for i in range(tp1.n(0))]

def eulerian_eliminationHM(PolyLst, VrbLst):
    """
    Outputs list of contraints whose degree matrix is in row echelon form.
    The general problem of determining the existence of solutions to a
    system of polynomial equations having at most finitely many solutions
    is NP hard. This implementation should therefore be used with caution.


    EXAMPLES:
 
    ::

        sage: sz=3; VrbLst=HM(sz,'x').list(); Ha=HM(sz,sz,'a'); Hb=HM(sz,1,HM(sz,'b').list())
        sage: CnstrLst=(Ha*HM(sz,1,VrbLst)-Hb).list()
        sage: Lf=eulerian_eliminationHM(CnstrLst, VrbLst)
        sage: Lf
        [a00*x0 + a01*x1 + a02*x2 - b0,
         -(a11*x1 + a12*x2 - b1)*a00 + (a01*x1 + a02*x2 - b0)*a10,
         ((a22*x2 - b2)*a00 - (a02*x2 - b0)*a20)*(a01*a10 - a00*a11) - ((a12*x2 - b1)*a00 - (a02*x2 - b0)*a10)*(a01*a20 - a00*a21)]
        sage: degree_matrix(Lf, var_list('x',sz)).printHM()
        [:, :]=
        [1 1 1]
        [0 1 1]
        [0 0 1]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    CnstrLst=copy(PolyLst)
    # Initializing the degree matrix.
    A=HM([[SR(CnstrLst[indx].degree(VrbLst[jndx])) for jndx in range(len(VrbLst))] for indx in range(len(CnstrLst))])
    #A.printHM()
    # Initialization of the row and column index
    i=0; j=0
    while i < A.n(0) and j < A.n(1):
        #while (A[i:,j]).is_zero() and j < A.ncols()-1:
        while HM(A.n(0)-i, 1, [A[i0,j] for i0 in range(i,A.n(0))]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        #if (A[i:,:].is_zero())==False:
        if HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))]).is_zero()==False:
            while A[i,j].is_zero(): 
                #Ta=A[i:,:]
                Ta=HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))])
                # Initializing the cyclic shift permutation matrix
                #Id=identity_matrix(Ta.nrows())
                Id=HM(2, Ta.n(0), 'kronecker')
                #P=sum([Id[:,k]*Id[mod(k+1,Ta.nrows()),:] for k in range(Ta.nrows())])
                P=sum([HM(Ta.n(0),1,[Id[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Id[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                Ta=P*Ta; CnstrLst=(P*HM(len(CnstrLst), 1, CnstrLst)).list()
                #A[i:,:]=Ta
                for i0 in range(Ta.n(0)):
                    for j0 in range(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
            # Performing the row operations.
            cf1=A[i,j]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if HM(1, A.n(1), [A[r,j0] for j0 in range(A.n(1))]).is_zero():
                    r=r+1
                else:
                    if (CnstrLst[r].degree(VrbLst[j]))*(CnstrLst[i].degree(VrbLst[j]))>0 and not SylvesterHM(CnstrLst[r], CnstrLst[i], VrbLst[j]).is_empty():
                        if not SylvesterHM(CnstrLst[r], CnstrLst[i], VrbLst[j]).det().is_zero():
                            CnstrLst[r]=SylvesterHM(CnstrLst[r], CnstrLst[i], VrbLst[j]).det()
                            #print 'i=', i,'j=', j,' r=', r
                            #print 'CnstrLst=', CnstrLst
                            #degree_matrix(CnstrLst,VrbLst).printHM()
                            A=HM([[SR(CnstrLst[indx].degree(VrbLst[jndx])) for jndx in range(len(VrbLst))] for indx in range(len(CnstrLst))])
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return CnstrLst

def euler_sylvester_eliminationHM(PolyLst, VrbLst):
    """
    Outputs list of contraints whose degree matrix is in reduced row echelon form.
    The general problem of determining the existence of solutions to a
    system of polynomial equations having at most finitely many solutions
    is NP hard. This implementation should therefore be used with caution.


    EXAMPLES:
 
    ::

        sage: sz=3; VrbLst=HM(sz,'x').list(); Ha=HM(sz,sz,'a'); Hb=HM(sz,1,HM(sz,'b').list())
        sage: CnstrLst=(Ha*HM(sz,1,VrbLst)-Hb).list()
        sage: Lf=euler_sylvester_eliminationHM(CnstrLst, VrbLst)
        sage: Lf
        [-((a02*a10 - a00*a12)*(a01*a20 - a00*a21) - (a01*a10 - a00*a11)*(a02*a20 - a00*a22))*(((a02*a10 - a00*a12)*(a01*a20 - a00*a21) - (a01*a10 - a00*a11)*(a02*a20 - a00*a22))*(a00*x0 - b0) + ((a01*a20 - a00*a21)*(a10*b0 - a00*b1) - (a01*a10 - a00*a11)*(a20*b0 - a00*b2))*a02)*(a01*a10 - a00*a11) + (((a01*a20 - a00*a21)*(a10*b0 - a00*b1) - (a01*a10 - a00*a11)*(a20*b0 - a00*b2))*(a02*a10 - a00*a12) - ((a02*a10 - a00*a12)*(a01*a20 - a00*a21) - (a01*a10 - a00*a11)*(a02*a20 - a00*a22))*(a10*b0 - a00*b1))*((a02*a10 - a00*a12)*(a01*a20 - a00*a21) - (a01*a10 - a00*a11)*(a02*a20 - a00*a22))*a01,
         ((a02*a10 - a00*a12)*(a01*a20 - a00*a21) - (a01*a10 - a00*a11)*(a02*a20 - a00*a22))*((a11*x1 - b1)*a00 - (a01*x1 - b0)*a10) - ((a01*a20 - a00*a21)*(a10*b0 - a00*b1) - (a01*a10 - a00*a11)*(a20*b0 - a00*b2))*(a02*a10 - a00*a12),
         ((a22*x2 - b2)*a00 - (a02*x2 - b0)*a20)*(a01*a10 - a00*a11) - ((a12*x2 - b1)*a00 - (a02*x2 - b0)*a10)*(a01*a20 - a00*a21)]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    CnstrLst=copy(eulerian_eliminationHM(PolyLst, VrbLst))
    # Initializing the degree matrix.
    A=HM(len(CnstrLst), len(VrbLst), [SR(CnstrLst[i].degree(VrbLst[j])) for j in range(len(VrbLst)) for i in range(len(CnstrLst))])
    # Initialization of the row and column index
    i=A.nrows()-1; j=0
    while i>0 or j>0:
        #if (A[i,:]).is_zero():
        if HM(1,A.n(1),[A[i,j0] for j0 in range(A.n(1))]).is_zero():
            # decrementing the row index and initializing the column index
            i=i-1; j=0
        else :
            while (A[i,j]).is_zero():
                # Incrementing the column index
                j = j + 1
            # performing row operations
            for r in range(i-1,-1,-1):
                if (CnstrLst[r].degree(VrbLst[j]))*(CnstrLst[i].degree(VrbLst[j]))>0 and not SylvesterHM(CnstrLst[r], CnstrLst[i], VrbLst[j]).is_empty():
                    if not SylvesterHM(CnstrLst[r], CnstrLst[i], VrbLst[j]).det().is_zero():
                        CnstrLst[r]=SylvesterHM(CnstrLst[r], CnstrLst[i], VrbLst[j]).det()
                        A=HM([[SR(CnstrLst[indx].degree(VrbLst[jndx])) for jndx in range(len(VrbLst))] for indx in range(len(CnstrLst))])
            i=i-1; j=0
    return CnstrLst

def eulerian_elimination_reductionHM(PolyLst, VrbLst, Rlts):
    """
    Outputs list of contraints whose degree matrix is in row echelon form.
    The general problem of determining the existence of solutions to a
    system of polynomial equations having at most finitely many solutions
    is NP hard. This implementation should therefore be used with caution.
    The polynomial expressions obtained are reduced modulo the single variable
    variables relations inputed in Rlts. The input list of polynomials must
    be given in their expanded form


    EXAMPLES:
 
    ::

        sage: sz=3; VrbLst=var_list('x',sz); f=sum(var_list('a',sz)[k]*x^(k-1) for k in rg(1,sz))
        sage: Rlts=[prod(VrbLst[i]-j for j in rg(sz-1)) for i in rg(sz)]
        sage: CnstrLst=Rlts[:sz-1]+[sum(f.subs(x==VrbLst[k]) for k in rg(sz))]
        sage: Lf=eulerian_elimination_reductionHM(CnstrLst, VrbLst, Rlts)
        sage: Lf
        [(x0 - 1)*x0,
         (x1 - 1)*x1,
         48*(135*a1^4*a2^14*x2 + 270*a1^3*a2^15*x2 + 201*a1^2*a2^16*x2 + 66*a1*a2^17*x2 + 8*a2^18*x2 + 81*a1^5*a2^13 + 135*a1^4*a2^14 + 81*a1^3*a2^15 + 21*a1^2*a2^16 + 2*a1*a2^17)*(6*a1*a2^11*x2 + 3*a2^12*x2 + 9*a1^2*a2^10 + 6*a1*a2^11 + a2^12)*(a2^3*x2 + 3*a1*a2^2 + a2^3)*a2^3]
        sage: degree_matrix(Lf, VrbLst).printHM()
        [:, :]=
        [2 0 0]
        [0 2 0]
        [0 0 3]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    CnstrLst=copy(PolyLst)
    # Initializing the degree matrix.
    A=HM([[SR(CnstrLst[indx].degree(VrbLst[jndx])) for jndx in range(len(VrbLst))] for indx in range(len(CnstrLst))])
    #A.printHM()
    # Initialization of the row and column index
    i=0; j=0
    while i < A.n(0) and j < A.n(1):
        #while (A[i:,j]).is_zero() and j < A.ncols()-1:
        while HM(A.n(0)-i, 1, [A[i0,j] for i0 in range(i,A.n(0))]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        #if (A[i:,:].is_zero())==False:
        if HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))]).is_zero()==False:
            while A[i,j].is_zero(): 
                #Ta=A[i:,:]
                Ta=HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))])
                # Initializing the cyclic shift permutation matrix
                #Id=identity_matrix(Ta.nrows())
                Id=HM(2, Ta.n(0), 'kronecker')
                #P=sum([Id[:,k]*Id[mod(k+1,Ta.nrows()),:] for k in range(Ta.nrows())])
                P=sum([HM(Ta.n(0),1,[Id[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Id[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                Ta=P*Ta; CnstrLst=(P*HM(len(CnstrLst), 1, CnstrLst)).list()
                #A[i:,:]=Ta
                for i0 in range(Ta.n(0)):
                    for j0 in range(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
            # Performing the row operations.
            cf1=A[i,j]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if HM(1, A.n(1), [A[r,j0] for j0 in range(A.n(1))]).is_zero():
                    r=r+1
                else:
                    if (CnstrLst[r].degree(VrbLst[j]))*(CnstrLst[i].degree(VrbLst[j])) > 0 and not SylvesterHM(CnstrLst[r], CnstrLst[i], VrbLst[j]).is_empty():
                        #if not SylvesterHM(CnstrLst[r], CnstrLst[i], VrbLst[j]).det().is_zero():
                        Cf=SylvesterHM(CnstrLst[r], CnstrLst[i], VrbLst[j])
                        rs=HM(Cf.n(0),1,'zero')
                        #print prod(gaussian_elimination_ReductionHM(Cf, rs, VrbLst, Rlts)[0][z,z] for z in range(Cf.n(0)))
                        if not prod(gaussian_elimination_ReductionHM(Cf, rs, VrbLst, Rlts)[0][z,z] for z in range(Cf.n(0))).is_zero():
                            #CnstrLst[r]=SylvesterHM(CnstrLst[r], CnstrLst[i], VrbLst[j]).det()
                            CnstrLst[r]=prod(gaussian_elimination_ReductionHM(Cf, rs, VrbLst, Rlts)[0][z,z] for z in range(Cf.n(0)))
                            #print 'i=', i,'j=', j,' r=', r
                            #print 'CnstrLst=', CnstrLst
                            A=HM([[SR(CnstrLst[indx].degree(VrbLst[jndx])) for jndx in range(len(VrbLst))] for indx in range(len(CnstrLst))])
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return CnstrLst

def degree_matrix(EqL, VrbL):
    """
    Outputs the degree matrix associated with the in input system 
    relative to the input list of variables.


    EXAMPLES:
 
    ::

        sage: sz=3; VrbL=HM(sz,'x').list(); Ha=HM(sz,sz,'a'); Hb=HM(sz,1,HM(sz,'b').list())
        sage: EqL=(Ha*HM(sz,1,VrbL)-Hb).list()
        sage: Ha=degree_matrix(EqL, VrbL)
        sage: Ha.printHM()
        [:, :]=
        [1 1 1]
        [1 1 1]
        [1 1 1]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    return HM([[EqL[i].degree((VrbL)[j]) for j in range(len(VrbL))] for i in range(len(EqL))])

def i2x2(A):
    """
    Outputs the symbolic inverse of a 2x2 matrix.


    EXAMPLES:
 
    ::

        sage: i2x2(HM(2,2,'a')).printHM()
        [:, :]=
        [-a11/(a01*a10 - a00*a11)  a01/(a01*a10 - a00*a11)]
        [ a10/(a01*a10 - a00*a11) -a00/(a01*a10 - a00*a11)]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    return HM([\
[ A[1,1]/(A[0,0]*A[1,1]-A[0,1]*A[1,0]), -A[0,1]/(A[0,0]*A[1,1]-A[0,1]*A[1,0])],\
[-A[1,0]/(A[0,0]*A[1,1]-A[0,1]*A[1,0]),  A[0,0]/(A[0,0]*A[1,1]-A[0,1]*A[1,0])]])

def inxn(A):
    """
    Outputs the symbolic inverse of a 2x2 matrix.


    EXAMPLES:
 
    ::

        sage: inxn(HM(2,2,'a')).printHM()
        [:, :]=
        [-a11/(a01*a10 - a00*a11)  a01/(a01*a10 - a00*a11)]
        [ a10/(a01*a10 - a00*a11) -a00/(a01*a10 - a00*a11)]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the computation of the determinant
    G=Deter(A); sz=A.n(0)
    return HM(sz, sz, [diff(ln(G),A[j,i]) for j in rg(sz) for i in rg(sz)])

def LeftRightDiagonalDependence3x3x3(A):
    """
    Outputs the a pair of solutions to the left right diagonal dependence problem.


    EXAMPLES:
 
    ::

        sage: A=HM([[[0,1,1],[6,-5,-1],[1,0,19]],[[1,2,3],[3,-2,-2],[1,1,0]],[[-6,0,2],[1,-3,-1],[-3,2,1]]])
        sage: [[Xfa,Yfa], [Xfb,Yfb]]=LeftRightDiagonalDependence3x3x3(A)
        sage: sz=3; sum(HM(sz,sz,[Xfa[i,t]*A[i,j,t]*Yfa[t,j] for j in range(sz) for i in range(sz)]) for t in range(sz)).canonicalize_radical()
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        sage: sz=3; sum(HM(sz,sz,[Xfb[i,t]*A[i,j,t]*Yfb[t,j] for j in range(sz) for i in range(sz)]) for t in range(sz)).canonicalize_radical()
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        sage: (HM(2,Xfa.matrix()[:,2].list(),'diag')*HM(2,Yfa.matrix()[2,:].list(),'diag')).det().is_zero() # Testing the first slice 
        False
        sage: (HM(2,Xfb.matrix()[:,2].list(),'diag')*HM(2,Yfb.matrix()[2,:].list(),'diag')).det().is_zero() # Testing the first slice 
        False


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    sz=A.nrows()
    # Initialization of the variable associated with the diagonal dependence.
    X=HM(sz,sz,'x'); Y=HM(sz,sz,'y')
    # Initializing the permutations
    P=Permutations(range(sz))
    # Initialization of constraints
    Eq=[sum(Permutation([p[indx]+1 for indx in range(len(p))]).signature()*\
    prod([X[i,p[i]]*A[i,k,p[i]] for i in range(sz)]) for p in P)==0 for k in range(sz)]
    # Formating the constraints using the first column of X as variables
    [F,g]=ConstraintFormatorHM(Eq, X.list()[-sz:])
    # Obtatining the determinant of the matrix F
    dtF=F.det()
    # Obtaining the parametrization by solving in the variable X[sz-1, sz-2].
    v=X[sz-1, sz-2]
    a0=dtF.subs(v==0)
    a1=diff(dtF, v, 1).subs(v==0)/factorial(1)
    a2=diff(dtF, v, 2).subs(v==0)/factorial(2)
    # Initialization of the solution of the quadratic equation
    Sln1=[\
    v == -1/2*( a1 + sqrt(a1^2 - 4*a0*a2) )/a2,\
    v == -1/2*( a1 - sqrt(a1^2 - 4*a0*a2) )/a2]
    # Performing the substitution
    Fa=F.subs(Sln1[0]).matrix()
    Sln2a=[X[i,sz-1]==-(X[sz-1,sz-1]*i2x2(Fa[:sz-1,:sz-1]).matrix()*Fa[:sz-1,sz-1])[i,0] for i in range(sz-1)]+Sln1[:1]
    Fb=F.subs(Sln1[1]).matrix()
    Sln2b=[X[i,sz-1]==-(X[sz-1,sz-1]*i2x2(Fb[:sz-1,:sz-1]).matrix()*Fb[:sz-1,sz-1])[i,0] for i in range(sz-1)]+Sln1[1:]
    # Initialization of the list of Mks 
    Xfa=X.subs(Sln2a)
    Lma=[HM(sz,sz,[Xfa[i,j]*A[i,k,j] for j in range(sz) for i in range(sz)]) for k in range(sz)]
    Xfb=X.subs(Sln2b)
    Lmb=[HM(sz,sz,[Xfb[i,j]*A[i,k,j] for j in range(sz) for i in range(sz)]) for k in range(sz)]
    # Initialization of solutions in Y
    Sln3a=[Y[j,k]==\
    -(Y[sz-1,k]*i2x2((Lma[k].matrix())[:sz-1,:sz-1]).matrix()*(Lma[k].matrix())[:sz-1,sz-1])[j,0] for j in range(sz-1) for k in range(sz)]
    Sln3b=[Y[j,k]==\
    -(Y[sz-1,k]*i2x2((Lmb[k].matrix())[:sz-1,:sz-1]).matrix()*(Lmb[k].matrix())[:sz-1,sz-1])[j,0] for j in range(sz-1) for k in range(sz)]
    Yfa=Y.subs(Sln3a)
    # Final verification
    #Rsa=sum(HM(sz,sz,[Xfa[i,t]*A[i,j,t]*Yfa[t,j] for j in range(sz) for i in range(sz)]) for t in range(sz)).canonicalize_radical()
    Yfb=Y.subs(Sln3b)
    # Final verification
    #Rsb=sum(HM(sz,sz,[Xfb[i,t]*A[i,j,t]*Yfb[t,j] for j in range(sz) for i in range(sz)]) for t in range(sz)).canonicalize_radical()
    return [[Xfa,Yfa], [Xfb,Yfb]]

def Reduced3x3x3CanonicalFactorization(A, Xfa, Yfa, indx):
    """
    Outputs the a pair of solutions to the left right diagonal dependence problem.


    EXAMPLES:
 
    ::

        sage: A=HM([[[0,1,1],[6,-5,-1],[1,0,19]],[[1,2,3],[3,-2,-2],[1,1,0]],[[-6,0,2],[1,-3,-1],[-3,2,1]]])
        sage: [[Xfa, Yfa], [Xfb,Yfb]]=LeftRightDiagonalDependence3x3x3(A)
        sage: (HM(2,Xfa.matrix()[:,2].list(),'diag')*HM(2,Yfa.matrix()[2,:].list(),'diag')).det().is_zero() # Testing the first slice 
        False
        sage: [Uf, Vf, Wf]=Reduced3x3x3CanonicalFactorization(A, Xfa, Yfa, 2)
        sage: Prod(Uf,Vf,Wf).canonicalize_radical() # Checking the factorization
        [[[0, 1, 1], [6, -5, -1], [1, 0, 19]], [[1, 2, 3], [3, -2, -2], [1, 1, 0]], [[-6, 0, 2], [1, -3, -1], [-3, 2, 1]]] 
        sage: Prod(Uf,Vf,Wf).canonicalize_radical()==A
        True


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the size and order parameters
    sz=3; od=3 
    # Checking that the last column of X and the last row of Y have no zero entries
    #print "(HM(2, Xfa.matrix()[:,2].list(),'diag')*HM(2, Yfa.matrix()[2,:].list(),'diag')).det().is_zero() is ",\
    #(HM(2, Xfa.matrix()[:,2].list(),'diag')*HM(2, Yfa.matrix()[2,:].list(),'diag')).det().is_zero()
    if (HM(2, Xfa.matrix()[:,2].list(),'diag')*HM(2, Yfa.matrix()[2,:].list(),'diag')).det().is_zero():
        raise ValueError("The input index has zero entries")
    else:
        # A convenien way to Initialize of the identity pairs
        J1=Prod(HM(sz,sz,sz,'one'), HM(sz,sz,sz,'one'), HM(od,sz,'kronecker')); U=J1.copy()
        J2=Prod(HM(od,sz,'kronecker'), HM(sz,sz,sz,'one'), HM(sz,sz,sz,'one')); W=J2.copy()
        # Initialization of the hypermatrix
        Ha=A.copy()
        # Updating the last slice Hypermatrix
        for i in range(sz):
            for j in range(sz):
                Ha[i,j,indx]=0
        # Initialization of the left weight coefficient
        XXfa=Xfa.copy()
        for i in range(sz):
            for j in range(sz):
                if j!=indx:
                    XXfa[i,j]=-Xfa[i,j]/Xfa[i,indx]
        for i in range(sz):
            XXfa[i,indx]=1
        # Initialization of the right weight coefficient
        YYfa=Yfa.copy()
        for i in range(sz):
            for j in range(sz):
                if i!=indx:
                    YYfa[i,j]=Yfa[i,j]/Yfa[indx,j]
        for j in range(sz):
            YYfa[indx,j]=1
        # Updating the slices of U
        for i in range(sz):
            for t in range(sz):
                for k in range(sz):
                    if t!=indx:
                        U[i,t,k]=XXfa[i,t]*J1[i,indx,k]+J1[i,t,k]
        # Updating the slices of W
        for j in range(sz):
            for t in range(sz-1):
                for k in range(sz):
                    if t!=indx:
                        W[t,j,k]=J2[t,j,k]+J2[indx,j,k]*YYfa[t,j]
        # Printing the difference
        #print '\n'
        #(A-Prod(U,Ha,W)).canonicalize_radical().printHM()
        # Obtaining the finall output
        Uf=HM(sz,sz-1,sz,[U[i,j,k] for k in range(sz) for j in range(sz-1) for i in range(sz)])
        Vf=HM(sz,sz,sz-1,[A[i,j,k] for k in range(sz-1) for j in range(sz) for i in range(sz)])
        Wf=HM(sz-1,sz,sz,[W[i,j,k] for k in range(sz) for j in range(sz) for i in range(sz-1)])
        #print '\n'
        #(A-Prod(Uf,Vf,Wf)).canonicalize_radical().printHM()
        return [Uf, Vf, Wf]

def Remnant(p, q, vrbl):
    """
    Takes as input two polynomials and a variable
    and outputs the corresponding remnant with the
    modular arithmetic method.
 

    EXAMPLES:

    ::

        sage: x, a0, a1, a2, b0, b1, b2=var('x, a0, a1, a2, b0, b1, b2')
        sage: p=(x-a0)*(x-a1); q=(x-b0)*(x-b1)
        sage: factor(Remnant(p, q, x).det())
        (a0 - b0)*(a0 - b1)*(a1 - b0)*(a1 - b1)
 

    AUTHORS:
    - Edinah K. Gnang
    """
    # Updating the firt input polynomial to make it monic in vrbl
    p=p/(p.diff(vrbl,Integer(p.degree(vrbl))).subs(vrbl==0)/factorial(Integer(p.degree(vrbl))))
    #print 'p=',p
    # Initialization of the list
    L=[]; f=expand(q)
    #print 'Initial f=',f
    for d in range(Integer(f.degree(vrbl)-p.degree(vrbl)),-1,-1):
        f=expand(fast_reduce(f,[vrbl^(d+p.degree(vrbl))],[vrbl^(d+p.degree(vrbl))-expand(p*vrbl^d)]))
        #print '    f=',f
    L.append(f)
    while len(L) < p.degree(vrbl):
        # Initialization of the update of q
        #f=expand(q*L[len(L)-1])
        f=expand(vrbl*L[len(L)-1])
        for d in range(f.degree(vrbl)-p.degree(vrbl),-1,-1):
            f=expand(fast_reduce(f,[vrbl^(d+p.degree(vrbl))],[vrbl^(d+p.degree(vrbl))-expand(p*vrbl^d)]))
        L.append(f)
    # Initialisation of the matrix
    return HM(p.degree(vrbl), p.degree(vrbl),[diff(L[i],vrbl,j).subs(vrbl==0)/factorial(j) for j in range(p.degree(vrbl)) for i in range(len(L))])

def RemnantII(p, q, vrbl):
    """
    Takes as input two polynomials and a variable
    and outputs the corresponding remnant with the
    modular arithmetic method.
 

    EXAMPLES:

    ::

        sage: x, a0, a1, a2, b0, b1, b2=var('x, a0, a1, a2, b0, b1, b2')
        sage: p=(x-a0)*(x-a1); q=(x-b0)*(x-b1)
        sage: factor(RemnantII(p, q, x).det())
        (a0 + a1 - b0 - b1)*(a0 - b0)*(a0 - b1)*(a1 - b0)*(a1 - b1)
 

    AUTHORS:
    - Edinah K. Gnang
    """
    # Updating the firt input polynomial to make it monic in vrbl
    p=p/(p.diff(vrbl,Integer(p.degree(vrbl))).subs(vrbl==0)/factorial(Integer(p.degree(vrbl))))
    #print 'p=',p
    # Initialization of the list
    L=[]; f=expand(q)
    #print 'Initial f=',f
    for d in range(Integer(f.degree(vrbl)-p.degree(vrbl)),-1,-1):
        f=expand(fast_reduce(f,[vrbl^(d+p.degree(vrbl))],[vrbl^(d+p.degree(vrbl))-expand(p*vrbl^d)]))
        #print '    f=',f
    L.append(f)
    while len(L) < p.degree(vrbl):
        # Initialization of the update of q
        f=expand(q*L[len(L)-1])
        for d in range(f.degree(vrbl)-p.degree(vrbl),-1,-1):
            f=expand(fast_reduce(f,[vrbl^(d+p.degree(vrbl))],[vrbl^(d+p.degree(vrbl))-expand(p*vrbl^d)]))
        L.append(f)
    # Initialisation of the matrix
    return HM(p.degree(vrbl), p.degree(vrbl),[diff(L[i],vrbl,j).subs(vrbl==0)/factorial(j) for j in range(p.degree(vrbl)) for i in range(len(L))]) 

def modular_eliminationHM(PolyLst, VrbLst):
    """
    Outputs list of contraints whose degree matrix is in row echelon form.
    The general problem of determining the existence of solutions to a
    system of polynomial equations having at most finitely many solutions
    is NP hard. This implementation should therefore be used with caution.


    EXAMPLES:
 
    ::

        sage: sz=3; VrbLst=var_list('x',sz); Ha=Vandermonde(range(1,sz+1)); Hb=HM(sz,1,var_list('b',sz))
        sage: CnstrLst=(Ha*HM(sz,1,VrbLst)-Hb).list()
        sage: Lf=modular_eliminationHM(CnstrLst, VrbLst); Lf
        [-b0 + x0 + x1 + x2, -b0 + b1 - x1 - 2*x2, -2/3*b0 + b1 - 1/3*b2 + 2/3*x2]
        sage: degree_matrix(Lf, var_list('x',sz)).printHM()
        [:, :]=
        [1 1 1]
        [0 1 1]
        [0 0 1]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    CnstrLst=copy(PolyLst)
    # Initializing the degree matrix.
    A=HM([[SR(CnstrLst[indx].degree(VrbLst[jndx])) for jndx in range(len(VrbLst))] for indx in range(len(CnstrLst))])
    #A.printHM()
    # Initialization of the row and column index
    i=0; j=0
    while i < A.n(0) and j < A.n(1):
        #while (A[i:,j]).is_zero() and j < A.ncols()-1:
        while HM(A.n(0)-i, 1, [A[i0,j] for i0 in range(i,A.n(0))]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        #if (A[i:,:].is_zero())==False:
        if HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))]).is_zero()==False:
            while A[i,j].is_zero(): 
                #Ta=A[i:,:]
                Ta=HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))])
                # Initializing the cyclic shift permutation matrix
                #Id=identity_matrix(Ta.nrows())
                Id=HM(2, Ta.n(0), 'kronecker')
                #P=sum([Id[:,k]*Id[mod(k+1,Ta.nrows()),:] for k in range(Ta.nrows())])
                P=sum([HM(Ta.n(0),1,[Id[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Id[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                Ta=P*Ta; CnstrLst=(P*HM(len(CnstrLst), 1, CnstrLst)).list()
                #A[i:,:]=Ta
                for i0 in range(Ta.n(0)):
                    for j0 in range(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
            # Performing the row operations.
            cf1=A[i,j]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if HM(1, A.n(1), [A[r,j0] for j0 in range(A.n(1))]).is_zero():
                    r=r+1
                else:
                    if (CnstrLst[r].degree(VrbLst[j]))*(CnstrLst[i].degree(VrbLst[j]))>0 and not Remnant(CnstrLst[r], CnstrLst[i], VrbLst[j]).is_empty():
                        if not Remnant(CnstrLst[r], CnstrLst[i], VrbLst[j]).det().is_zero():
                            CnstrLst[r]=Remnant(CnstrLst[r], CnstrLst[i], VrbLst[j]).det()
                            #print 'i=', i,'j=', j,' r=', r
                            #print 'CnstrLst=', CnstrLst
                            A=HM([[SR(CnstrLst[indx].degree(VrbLst[jndx])) for jndx in range(len(VrbLst))] for indx in range(len(CnstrLst))])
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return CnstrLst

def complete_modular_eliminationHM(PolyLst, VrbLst):
    """
    Outputs list of contraints whose degree matrix is in reduced row echelon form.
    The general problem of determining the existence of solutions to a
    system of polynomial equations having at most finitely many solutions
    is NP hard. This implementation should therefore be used with caution.


    EXAMPLES:
 
    ::

        sage: x, a0, a1, a2, b0, b1, b2=var('x, a0, a1, a2, b0, b1, b2')
        sage: p=(x-a0)*(x-a1); q=(x-b0)*(x-b1)
        sage: [f.factor() for f in complete_modular_eliminationHM([p,q],[x])]
        [(a0 - x)*(a1 - x), (a0 - b0)*(a0 - b1)*(a1 - b0)*(a1 - b1)]
        sage: sz=3; VrbLst=var_list('x',sz); Ha=Vandermonde(range(1,sz+1)); Hb=HM(sz,1,var_list('b',sz))
        sage: CnstrLst=(Ha*HM(sz,1,VrbLst)-Hb).list()
        sage: Lf=complete_modular_eliminationHM(CnstrLst, VrbLst)
        sage: Lf
        [-b0 + 5/6*b1 - 1/6*b2 + 1/3*x0,
         -b0 + 4/3*b1 - 1/3*b2 - 1/3*x1,
         -2/3*b0 + b1 - 1/3*b2 + 2/3*x2]        
        sage: degree_matrix(Lf, var_list('x',sz)).printHM()
        [:, :]=
        [1 0 0]
        [0 1 0]
        [0 0 1]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    CnstrLst=copy(modular_eliminationHM(PolyLst, VrbLst))
    # Initializing the degree matrix.
    A=HM(len(CnstrLst), len(VrbLst), [SR(CnstrLst[i].degree(VrbLst[j])) for j in range(len(VrbLst)) for i in range(len(CnstrLst))])
    # Initialization of the row and column index
    i=A.nrows()-1; j=0
    while i>0 or j>0:
        #if (A[i,:]).is_zero():
        if HM(1,A.n(1),[A[i,j0] for j0 in range(A.n(1))]).is_zero():
            # decrementing the row index and initializing the column index
            i=i-1; j=0
        else :
            while (A[i,j]).is_zero():
                # Incrementing the column index
                j = j + 1
            # performing row operations
            for r in range(i-1,-1,-1):
                if (CnstrLst[r].degree(VrbLst[j]))*(CnstrLst[i].degree(VrbLst[j]))>0 and not Remnant(CnstrLst[r], CnstrLst[i], VrbLst[j]).is_empty():
                    if not Remnant(CnstrLst[r], CnstrLst[i], VrbLst[j]).det().is_zero():
                        CnstrLst[r]=Remnant(CnstrLst[r], CnstrLst[i], VrbLst[j]).det()
                        A=HM([[SR(CnstrLst[indx].degree(VrbLst[jndx])) for jndx in range(len(VrbLst))] for indx in range(len(CnstrLst))])
            i=i-1; j=0
    return CnstrLst

def modular_eliminationHMII(PolyLst, VrbLst):
    """
    Outputs list of contraints whose degree matrix is in row echelon form.
    The general problem of determining the existence of solutions to a
    system of polynomial equations having at most finitely many solutions
    is NP hard. This implementation should therefore be used with caution.


    EXAMPLES:
 
    ::

        sage: sz=3; VrbLst=var_list('x',sz); Ha=Vandermonde(range(1,sz+1)); Hb=HM(sz,1,var_list('b',sz))
        sage: CnstrLst=(Ha*HM(sz,1,VrbLst)-Hb).list()
        sage: Lf=modular_eliminationHMII(CnstrLst, VrbLst); Lf
        [-b0 + x0 + x1 + x2, -b0 + b1 - x1 - 2*x2, -2/3*b0 + b1 - 1/3*b2 + 2/3*x2]
        sage: degree_matrix(Lf, var_list('x',sz)).printHM()
        [:, :]=
        [1 1 1]
        [0 1 1]
        [0 0 1]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    CnstrLst=copy(PolyLst)
    # Initializing the degree matrix.
    A=HM([[SR(CnstrLst[indx].degree(VrbLst[jndx])) for jndx in range(len(VrbLst))] for indx in range(len(CnstrLst))])
    #A.printHM()
    # Initialization of the row and column index
    i=0; j=0
    while i < A.n(0) and j < A.n(1):
        #while (A[i:,j]).is_zero() and j < A.ncols()-1:
        while HM(A.n(0)-i, 1, [A[i0,j] for i0 in range(i,A.n(0))]).is_zero() and j < A.ncols()-1:
            # Incrementing the column index
            j=j+1
        #if (A[i:,:].is_zero())==False:
        if HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))]).is_zero()==False:
            while A[i,j].is_zero(): 
                #Ta=A[i:,:]
                Ta=HM(A.n(0)-i, A.n(1), [A[i0,j0] for j0 in range(A.n(1)) for i0 in range(i,A.n(0))])
                # Initializing the cyclic shift permutation matrix
                #Id=identity_matrix(Ta.nrows())
                Id=HM(2, Ta.n(0), 'kronecker')
                #P=sum([Id[:,k]*Id[mod(k+1,Ta.nrows()),:] for k in range(Ta.nrows())])
                P=sum([HM(Ta.n(0),1,[Id[i0,k] for i0 in range(Ta.n(0))])*HM(1,Ta.n(0),[Id[Integer(mod(k+1,Ta.n(0))),j0] for j0 in range(Ta.n(0))]) for k in range(Ta.n(0))])
                Ta=P*Ta; CnstrLst=(P*HM(len(CnstrLst), 1, CnstrLst)).list()
                #A[i:,:]=Ta
                for i0 in range(Ta.n(0)):
                    for j0 in range(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
            # Performing the row operations.
            cf1=A[i,j]
            for r in range(i+1,A.nrows()):
                # Taking care of the zero row
                if HM(1, A.n(1), [A[r,j0] for j0 in range(A.n(1))]).is_zero():
                    r=r+1
                else:
                    if (CnstrLst[r].degree(VrbLst[j]))*(CnstrLst[i].degree(VrbLst[j]))>0 and not RemnantII(CnstrLst[r], CnstrLst[i], VrbLst[j]).is_empty():
                        if not RemnantII(CnstrLst[r], CnstrLst[i], VrbLst[j]).det().is_zero():
                            CnstrLst[r]=RemnantII(CnstrLst[r], CnstrLst[i], VrbLst[j]).det()
                            #print 'i=', i,'j=', j,' r=', r
                            #print 'CnstrLst=', CnstrLst
                            A=HM([[SR(CnstrLst[indx].degree(VrbLst[jndx])) for jndx in range(len(VrbLst))] for indx in range(len(CnstrLst))])
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return CnstrLst

def complete_modular_eliminationHMII(PolyLst, VrbLst):
    """
    Outputs list of contraints whose degree matrix is in reduced row echelon form.
    The general problem of determining the existence of solutions to a
    system of polynomial equations having at most finitely many solutions
    is NP hard. This implementation should therefore be used with caution.


    EXAMPLES:
 
    ::

        sage: x, a0, a1, a2, b0, b1, b2=var('x, a0, a1, a2, b0, b1, b2')
        sage: p=(x-a0)*(x-a1); q=(x-b0)*(x-b1)
        sage: [f.factor() for f in complete_modular_eliminationHMII([p,q],[x])]
        [(a0 - x)*(a1 - x),
         -(a0 + a1 - b0 - b1)*(a0 - b0)*(a0 - b1)*(a1 - b0)*(a1 - b1)]
        sage: sz=3; VrbLst=var_list('x',sz); Ha=Vandermonde(range(1,sz+1)); Hb=HM(sz,1,var_list('b',sz))
        sage: CnstrLst=(Ha*HM(sz,1,VrbLst)-Hb).list()
        sage: Lf=complete_modular_eliminationHMII(CnstrLst, VrbLst)
        sage: Lf
        [-b0 + 5/6*b1 - 1/6*b2 + 1/3*x0,
         -b0 + 4/3*b1 - 1/3*b2 - 1/3*x1,
         -2/3*b0 + b1 - 1/3*b2 + 2/3*x2]        
        sage: degree_matrix(Lf, var_list('x',sz)).printHM()
        [:, :]=
        [1 0 0]
        [0 1 0]
        [0 0 1]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    CnstrLst=copy(modular_eliminationHMII(PolyLst, VrbLst))
    # Initializing the degree matrix.
    A=HM(len(CnstrLst), len(VrbLst), [SR(CnstrLst[i].degree(VrbLst[j])) for j in range(len(VrbLst)) for i in range(len(CnstrLst))])
    # Initialization of the row and column index
    i=A.nrows()-1; j=0
    while i>0 or j>0:
        #if (A[i,:]).is_zero():
        if HM(1,A.n(1),[A[i,j0] for j0 in range(A.n(1))]).is_zero():
            # decrementing the row index and initializing the column index
            i=i-1; j=0
        else :
            while (A[i,j]).is_zero():
                # Incrementing the column index
                j = j + 1
            # performing row operations
            for r in range(i-1,-1,-1):
                if (CnstrLst[r].degree(VrbLst[j]))*(CnstrLst[i].degree(VrbLst[j]))>0 and not RemnantII(CnstrLst[r], CnstrLst[i], VrbLst[j]).is_empty():
                    if not RemnantII(CnstrLst[r], CnstrLst[i], VrbLst[j]).det().is_zero():
                        CnstrLst[r]=RemnantII(CnstrLst[r], CnstrLst[i], VrbLst[j]).det()
                        A=HM([[SR(CnstrLst[indx].degree(VrbLst[jndx])) for jndx in range(len(VrbLst))] for indx in range(len(CnstrLst))])
            i=i-1; j=0
    return CnstrLst

def outerdeterminant(A, B):
    """
    Computes symbolically the outer-product expansions of the determinant
    using the sum over permutation formula. The inputs should be second 
    order hypermatrices.

    EXAMPLES:

    ::

        sage: sz=2; outerdeterminant(HM(sz,sz,'a'), HM(sz,sz,'b')).printHM()
        [:, :]= 
        [-a01*a10*b01*b10 + a00*a11*b01*b10 -a01*a10*b01*b11 + a00*a11*b01*b11]
        [ a01*a10*b00*b10 - a00*a11*b00*b10  a01*a10*b00*b11 - a00*a11*b00*b11]
        sage: outerdeterminant(HM(2,2,'a'), HM(2,2,'b')).trace().factor()
        -(a01*a10 - a00*a11)*(b01*b10 - b00*b11)


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing the permutations
    P = Permutations(range(A.nrows()))
    return sum([Permutation([p[i]+1 for i in range(len(p))]).signature()*prod([\
HM(A.nrows(),1,[A[i,p[k]] for i in range(A.nrows())])*\
HM(1,B.ncols(),[B[k,j] for j in range(B.ncols())]) for k in range(A.nrows())]) for p in P])

def CC2RR_inflate(A):
    """ 
    Outputs the inflated matrix obtained by replacing
    complex entries of A by their canonical 2x2 
    real matrix representations.


    EXAMPLES:

    ::  

        sage: A = HM([[-2*I + 3, -2], [5*I - 1, -I + 2]]); CC2RR_inflate(A).printHM()
        [:, :]=
        [ 3  2 -2  0]
        [-2  3  0 -2]
        [-1 -5  2  1]
        [ 5 -1 -1  2]


    AUTHORS:
    - Edinah K. Gnang
    """
    if A.order()==2:
        # Initialization of the matrix
        B=HM(2*A.n(0), 2*A.n(1), 'zero') # Conversion from complex to real
        for i in range(A.n(0)):
            for j in range(A.n(1)):
                Tmp=HM(A.n(0), A.n(1), 'zero'); Tmp[i,j]=1
                B=B+Tmp.tensor_product(HM([[A[i,j].real(),-A[i,j].imag()],[A[i,j].imag(),A[i,j].real()]]))
        return B
    else:
        raise ValueError("Expected a second order hypermatrix")
        
def RR2CC_deflate(Q):
    """ 
    Outputs the deflated matrix obtained by replacing
    the canonical 2x2 real matrix representations by
    complex numbers.


    EXAMPLES:

    ::  

        sage: A = HM([[3, 2, -2, 0], [-2, 3, 0, -2], [-1, -5, 2, 1], [5, -1, -1, 2]]); RR2CC_deflate(A).printHM()
        [:, :]=
        [-2*I + 3       -2]
        [ 5*I - 1   -I + 2]


    AUTHORS:
    - Edinah K. Gnang
    """
    if Q.order()==2:
        # Initialization of the matrix
        U=HM(Q.n(0)/2, Q.n(1)/2, 'zero') # Conversion from complex to real
        for i in range(0, Q.n(0), 2):
            for j in range(0, Q.n(1), 2):
                Tmp=HM(Q.n(0), Q.n(1), 'zero'); Tmp[i,j]=1
                U[i/2, j/2] = Q[i, j]+I*Q[i+1, j]
        return U
    else:
        raise ValueError("Expected a second order hypermatrix")
 
def quaternion_2x2_CC_rep(VrbL): 
    """ 
    Outputs the canonical 2x2 complex matrix 
    representation of quaternions using the input 
    variables. The input VrbL is a list of pairs 
    associated with real part and imaginary part
    of the input complex numbers 


    EXAMPLES:

    ::  

        sage: VrbLx=var_list('x',2); VrbLy=var_list('y',2)
        sage: quaternion_2x2_CC_rep([(VrbLx[i], VrbLy[i]) for i in range(2)]).printHM()
        [:, :]=
        [ x1 + I*y1  x0 + I*y0]
        [-x0 + I*y0  x1 - I*y1]


    AUTHORS:
    - Edinah K. Gnang
    """ 
    return HM(2,2,[VrbL[1][0]+I*VrbL[1][1], -VrbL[0][0]+I*VrbL[0][1],  VrbL[0][0]+I*VrbL[0][1], VrbL[1][0]-I*VrbL[1][1]])

def quaternion_2x2_unit_abs_CC(VrbL): 
    """ 
    Out puts the canonical 2x2 matrix representation
    using the input variables. Assumes that the 
    2 variables in VrbL are associated with complex
    numbers lying on the unit circle


    EXAMPLES:

    ::  

        sage: VrbL=var_list('x',2); quaternion_2x2_unit_abs_CC(VrbL).printHM()
        [:, :]=
        [   x1    x0]
        [-1/x0  1/x1]

    AUTHORS:
    - Edinah K. Gnang
    """ 
    return HM(2,2,[VrbL[1],-1/VrbL[0],  VrbL[0],1/VrbL[1]])

def rg(*args): 
    """ 
    Adapts the range function to our purposes.
    This function makes sure to return a list of
    integers and not object of type int.


    EXAMPLES:

    ::  

        sage: rg(0,10,2)
        [0, 2, 4, 6, 8]


    AUTHORS:
    - Edinah K. Gnang
    """ 
    return [Integer(i) for i in range(*args)]

def multivariate_leading_term(f, Xv, Pp):
    """
    Takes as input a polynomial in the variables
    specified in the list Xv and outputs the 
    leading term of the input polynomial.
    The last input is the list of primes
    assigned to each variable to determine
    the monomial ordering
 

    EXAMPLES:

    ::

        sage: sz=3; Xv=var_list('x',sz); P=Primes(); Pp=[P.unrank(i) for i in rg(sz)] 
        sage: f = 5*Xv[0]*Xv[1]^3*Xv[2]^2 + 4*Xv[1] + 2*Xv[1]*Xv[2]^7 + 5
        sage: multivariate_leading_term(f, Xv, Pp)
        2*x1*x2^7
 

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter
    sz=len(Xv)
    # Expression used for specifying the type of the operation.
    cst = 2
    add = var('x0') + var('x1')
    mul = var('x0') * var('x1')
    xpo = var('x0') ^ var('x1')
    if f.operator() == add.operator():
        # Collecting the terms
        L=f.operands()
        # Collecting the terms striped from their coefficients
        L_strpd = list()
        for i in rg(len(L)):
            if L[i].arguments() != ():
                if (L[i].operator()==mul.operator() or (L[i]).operator()==xpo.operator()):
                    cst_fctr = 1
                    lst_i = L[i].operands()
                    for j in rg(len(lst_i)):
                        if lst_i[j].arguments()==():
                            cst_fctr = lst_i[j]
                    L_strpd.append((expand(L[i]/cst_fctr), i))
                elif L[i] in Xv:
                    L_strpd.append((L[i], i))
        # Storing the integer and the index associated with the term
        tmp_value = L_strpd[0][0].subs([Xv[i]==Pp[i] for i in rg(sz)])
        idx = L_strpd[0][1]
        for k in rg(len(L_strpd)):
            if L_strpd[k][0].subs([Xv[i]==Pp[i] for i in rg(sz)]) > tmp_value:
                tmp_value = L_strpd[k][0].subs([Xv[i]==Pp[i] for i in rg(sz)])
                idx = L_strpd[k][1]
        return L[idx]
    elif f.operator() == mul.operator():
        return f
    elif f.operator() == xpo.operator():
        return f 
    else :
        return f

def multivariate_division(f, List, Xv, Pp):
    """
    Takes as input a polynomial f the list of 
    polynomials and performs the multivariable
    division algorithm and the list of variables
    used. The last input is a list of primes which
    determines the monomial ordering.
 

    EXAMPLES:

    ::

        sage: sz=3; Xv=var_list('x',sz); P=Primes(); Pp=[P.unrank(i) for i in rg(sz)] 
        sage: f = 5*Xv[0]*Xv[1]^3*Xv[2]^2 + 4*Xv[1] + 2*Xv[1]*Xv[2]^7 + 5
        sage: List = [Xv[0]^2+5*Xv[1]^3+2, Xv[1]^3*Xv[0]^2+5*Xv[0]+1]
        sage: multivariate_division(f, List, Xv, Pp)
        [x0*x2^2, 0, 2*x1*x2^7 - x0^3*x2^2 - 2*x0*x2^2 + 4*x1 + 5]
 

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter
    sz=len(Xv)
    # Initializing the output List
    L = []
    for j in rg(len(List)+1):
        L.append(0)
    # Initializing the Polynomial
    p = f
    while not p.is_zero():
        i = 0
        division_occured = False
        while (i in rg(len(List))) and (not division_occured):
            if List[i] == 0:
                i = i+1
            else:
                # Getting the Leading term of fi
                Lt_fi = multivariate_leading_term(List[i], Xv, Pp)
                # Getting the leading Monomial of fi
                Lm_fi = Lt_fi/Lt_fi.subs([Xv[j]==1 for j in rg(sz)])
                # Getting the Leading term of p
                Lt_p = multivariate_leading_term(p, Xv, Pp)
                # Getting the leading Monomial of p
                Lm_p = Lt_p/Lt_p.subs([Xv[j]==1 for j in rg(sz)])
                m_p  = Lm_p.subs([Xv[j]==Pp[j] for j in rg(sz)])
                m_fi = Lm_fi.subs([Xv[j]==Pp[j] for j in rg(sz)])
                if gcd(m_p, m_fi) == m_fi or gcd(m_p, m_fi) == -m_fi:
                    L[i] = expand(L[i] + Lt_p/Lt_fi)
                    p = expand(p - List[i] * Lt_p/Lt_fi)
                    division_occured = True
                else :
                    i = i+1
        if division_occured == False:
            L[len(List)] = L[len(List)] + multivariate_leading_term(p, Xv, Pp)
            p = p - multivariate_leading_term(p, Xv, Pp)
    return L

def multivariate_monomial_lcm(t1, t2, Xv, Pp):
    """
    Takes as input two terms t1 and t2 and
    returns a the least common multiple monomial.
    the function does not check the inputs.
    The last imput is a list of prime which
    determines the monomial ordering. 


    EXAMPLES:

    ::

        sage: sz=3; Xv=var_list('x',sz); P=Primes(); Pp=[P.unrank(i) for i in rg(sz)] 
        sage: multivariate_monomial_lcm(Xv[0]*Xv[1]^2*Xv[2], Xv[0]*Xv[2]^2, Xv, Pp)
        x0*x1^2*x2^2


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter
    sz=len(Xv)
    # Initialization of the Prime variable dictionary datastructure
    PpXv = dict([(Pp[i],Xv[i]) for i in rg(sz)])
    # These 2 lines of code get rid of the coefficient of the leading terms
    m1 = t1/t1.subs([Xv[j]==1 for j in rg(sz)])
    m2 = t2/t2.subs([Xv[j]==1 for j in rg(sz)])
    # The following computes the lcm value associated with the monomial we seek
    monomial_lcm_value = lcm(m1.subs([Xv[i]==Pp[i] for i in rg(sz)]), m2.subs([Xv[i]==Pp[i] for i in rg(sz)]))
    # The next section of line of codes recovers the 
    # monomial in question from the computed lcm integer.
    prime_factors = factor(Integer(monomial_lcm_value))
    factor_list = list(prime_factors)
    # Initialization of the monomial
    m = 1
    for i in rg(len(factor_list)):
        tmp_list = list(factor_list[i])
        m = m * PpXv[tmp_list[0]]^Integer(tmp_list[1])
    return m

def multivariate_S_polynomials(List, Xv, Pp):
    """
    Takes as input a list of polynomials
    returns the list of substracted polynomials.
    The last input is a list of primes which
    determines the monomial ordering.


    EXAMPLES:

    ::

        sage: sz=3; Xv=var_list('x',sz); P=Primes(); Pp=[P.unrank(i) for i in rg(sz)] 
        sage: multivariate_S_polynomials([Xv[0]*Xv[1]^2*Xv[2]+Xv[1]*Xv[2]+1, Xv[0]*Xv[2]^2+1], Xv, Pp)
        [x1*x2^2 - x1^2 + x2]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter
    sz=len(Xv)
    L = []
    for i in rg(len(List)-1):
        for j in rg(i+1,len(List)):
            Lt_fi = multivariate_leading_term(List[i], Xv, Pp)
            Lt_fj = multivariate_leading_term(List[j], Xv, Pp)
            monomial_lcm = multivariate_monomial_lcm(Lt_fi,Lt_fj, Xv, Pp)
            Sij = expand(List[i]*(monomial_lcm/Lt_fi) - List[j]*(monomial_lcm/Lt_fj))
            L.append(Sij)
    return L

def multivariate_reduce_polynomial(f, Xv, Pp):
    """
    Takes as input a polynomial f and removes any redundant monomial
    common factors between the terms of f. The reduction referes to 
    the reduced grobner bases. The last input is a list of primes
    which determines the variable ordering.


    EXAMPLES:

    ::

        sage: sz=3; Xv=var_list('x',sz); P=Primes(); Pp=[P.unrank(i) for i in rg(sz)] 
        sage: multivariate_reduce_polynomial(Xv[1]*Xv[2]^2 - Xv[1]^2 + Xv[1]*Xv[2], Xv, Pp)
        x2^2 - x1 + x2


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter
    sz=len(Xv)
    # Initialization of the Prime variable dictionary datastructure
    PpXv = dict([(Pp[i],Xv[i]) for i in rg(sz)])
    # Checks to see if there is a constant term 
    # in which case no reduction is needed
    if not f.subs([Xv[j]==0 for j in rg(sz)]).is_zero():
        return f
    elif f == 0:
        return f
    else:
        L = list(f.iterator())
        # The next piece of code determines i
        # if the expression is a monomial
        prd = 1
        for j in rg(len(L)):
            prd = prd*L[j]
        if prd == f:
            return 1
        elif (len(L) == 2) and (L[0]^L[1] == f):
            return f
        for i in rg(len(L)):
            # The next line of code gets rid of 
            # the coefficients in the list
            L[i] = L[i]/(L[i].subs([Xv[j]==1 for j in rg(sz)]))
            L[i] = Integer(L[i].subs([Xv[j]==Pp[j] for j in rg(sz)]))
        # Computing the greatest common divisior
        cmn_fctr = gcd(L)
        if cmn_fctr == 1 :
            return f
        else :
            # The next section of line of codes recover the monomial
            # in question from the computed gcd integer.
            prime_factors = factor(cmn_fctr)
            factor_list = list(prime_factors)
            m = 1
            for i in rg(len(factor_list)):
                tmp_list = list(factor_list[i])
                m = m * PpXv[tmp_list[0]]^Integer(tmp_list[1])
            g = expand(f/m)
            return g

def prime_induced_grobner_basis(Idl, Xv, Pp):
    """
    Takes as input a polynomial f and removes any redundant monomial
    common factors between the terms of f. The reduction referes to 
    the reduced grobner bases. The last input is a list of primes 
    which determines the variable ordering.


    EXAMPLES:

    ::

        sage: sz=3; Xv=var_list('x',sz); P=Primes(); Pp=[P.unrank(i) for i in rg(sz)] 
        sage: Idl=[expand((Xv[0]+2*Xv[1])*(2*Xv[2]))-6, Xv[2]^2-Xv[2], Xv[1]^2-Xv[1], Xv[0]^2-Xv[0]] 
        sage: prime_induced_grobner_basis(Idl, Xv, Pp)[0]
        2*x0*x2 + 4*x1*x2 - 6


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter
    sz=len(Xv)
    # Initialization step 
    I_curr = list()
    for i in rg(len(Idl)):
        I_curr.append(Idl[i])
    l_old = 0; l_new = len(I_curr)
    # Boolean variable tracking contradictions
    finished = False
    while l_old != l_new and finished == False:
        # Computes the single pass of the substraction polynomials
        S  = multivariate_S_polynomials(I_curr, Xv, Pp)
        #print '\n\n The subtraction polynomials yield'
        for i in rg(len(S)):
            S[i] = multivariate_reduce_polynomial(S[i], Xv, Pp)
            #print 'S[',i,']= ',S[i]
        # The instruction bellow is the lazy way of getting rid of the duplicates.
        St = Set(S)
        S = list(St)
        # Recording the size of the Ideal generator set before the division
        l_old = len(I_curr)
        for i in rg(len(S)):
            tmp_list = multivariate_division(S[i], I_curr, Xv, Pp)
            if tmp_list[len(tmp_list)-1]!=0:
                I_curr.append(tmp_list[len(tmp_list)-1])
        # Printing the result of the first pass of the Buchberger algorithm.
        #print '\n\n The Current generator for the Ideal is given by'
        for i in rg(len(I_curr)):
            #print I_curr[i]
            if I_curr[i] == I_curr[i].subs([Xv[j]==0 for j in rg(sz)]) and not I_curr[i].subs([Xv[j]==0 for j in rg(sz)]).is_zero():
                finished = True
        # Recording the size of the generator set after the division
        l_new = len(I_curr)
    return I_curr

def generate_general_linear_constraints(sz,l):
    """
    Creates a sage file which intializes a general linear 
    system of sz constraints in l variables.

    EXAMPLES:

    ::

        sage: generate_general_linear_constraints(3,2)
        sage: load('general_linear_system_3_2.sage')
        [:, :, 0]=
        [a00*x0*b00 + a01*x1*b10]
        <BLANKLINE>
        [:, :, 1]=
        [a10*x0*b01 + a11*x1*b11]
        <BLANKLINE>
        [:, :, 2]=
        [a20*x0*b02 + a21*x1*b12]
        sage: from subprocess import call
        sage: call("rm general_linear_system_3_2.sage", shell=True)
        0
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the hypermatrices.
    Al=HM(sz,l,'a').list(); Bl=HM(l,sz,'b').list(); Xl=var_list('x',l)
    # Creating the file name string.
    filename='general_linear_system_'+str(sz)+'_'+str(l)+'.sage'
    # Opening the file
    f=open(filename,'w')
    #f.write('# Loading the Hypermatrix Package\n')
    #f.write("load('./Hypermatrix_Algebra_tst.sage')\n\n")
    f.write('# Initializing the number of constraints and the number of variableas\n')
    f.write('sz='+str(sz)+'; l='+str(l)+'\n\n')
    f.write('# Initialization of the variables\n')
    f.write("Lx=var_list('x',l)\n")
    f.write("La=HM(sz,l,'a').list()\n")
    f.write("Lb=HM(l,sz,'b').list()\n\n")
    f.write('# Initializing the free variables\n')
    f.write('F=FreeAlgebra(QQ,len(La+Lx+Lb),La+Lx+Lb)\n')
    f.write('F.<'+str(Al+Xl+Bl)[1:len(str(Al+Xl+Bl))-1]+'>=FreeAlgebra(QQ,len(La+Lx+Lb))\n\n')
    f.write('# Initialization of the hypermatrices with symbolic variable entries which do not commute\n')
    f.write('# associated with the map\n')
    f.write('Ha=HM(1,l,sz,HM(sz,l,'+str(Al)+').transpose().list())\n')
    f.write('Hx=HM(1,1,l,'+str(Xl)+')\n')
    f.write('Hb=HM(l,1,sz,'+str(Bl)+')\n\n')
    f.write('# Initialization of the product\n')
    f.write('Hr=Prod(Ha,Hx,Hb)\n')
    f.write('Hr.printHM()\n')
    # Closing the file
    f.close()

def GeneralHypermatrixSlicer(A, Rg, indx):
    """
    Outputs slices specified by index list L.
    the last string input is either row or col
    and determines the slices to be collected
    into a hypermatrix in the specified order
    by the input list Rg.


    EXAMPLES:
 
    ::

        sage: sz=3; A=HM(sz,sz,'a')
        sage: GeneralHypermatrixSlicer(A, [0], 'row').printHM()
        [:, :]=
        [a00 a01 a02]
        sage: GeneralHypermatrixSlicer(A, [1], 'col').printHM()
        [:, :]=
        [a01]
        [a11]
        [a21]
        sage: GeneralHypermatrixSlicer(A, [1,0], 'col').printHM()
        [:, :]=
        [a01 a00]
        [a11 a10]
        [a21 a20]
        sage: GeneralHypermatrixSlicer(A, [1,0], 'row').printHM()
        [:, :]=
        [a10 a11 a12]
        [a00 a01 a02]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the desired permutation
    Pmt=get_permutation(sorted(Rg),Rg)
    # Initialization of a number for the test
    nb=2
    if type(indx) == type(nb):
        if indx < A.order() :
            if len(Rg) <= A.n(indx):
                # Initialization of the hypermatrix which stores the result
                dms=A.dimensions(); dms[indx]=len(Rg); 
                # Initializing the list of entries
                Lst=[]
                # Initialization of the list specifying the dimensions of the output
                l = [A.n(i) for i in range(A.order())]
                # Main loop performing the transposition of the entries
                for i in range(prod(l)):
                    # Turning the index i into an hypermatrix array location using the decimal encoding trick
                    entry = [Integer(mod(i,l[0]))]
                    sm = Integer(mod(i,l[0]))
                    for k in range(len(l)-1):
                        entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                        sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                    if entry[indx] in Rg:
                        entry[indx] = Rg[Pmt[Rg.index(entry[indx])]]
                        Lst.append(A[tuple(entry)])
                #return apply(HM,dms+[Lst])
                return HM(*(dms+[Lst]))
            else:
                raise ValueError("The range must be smaller then corresponding index range")
        else:
            raise ValueError("The index must be smaller then the order of the Hypermatrix")
    elif type(indx) == type('tst'):
        if indx == 'row':
            indx=0
            if len(Rg) <= A.n(indx):
                # Initialization of the hypermatrix which stores the result
                dms=A.dimensions(); dms[indx]=len(Rg); 
                # Initializing the list of entries
                Lst=[]
                # Initialization of the list specifying the dimensions of the output
                l = [A.n(i) for i in range(A.order())]
                # Main loop performing the transposition of the entries
                for i in range(prod(l)):
                    # Turning the index i into an hypermatrix array location using the decimal encoding trick
                    entry = [Integer(mod(i,l[0]))]
                    sm = Integer(mod(i,l[0]))
                    for k in range(len(l)-1):
                        entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                        sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                    if entry[indx] in Rg:
                        entry[indx] = Rg[Pmt[Rg.index(entry[indx])]]
                        Lst.append(A[tuple(entry)])
                #return apply(HM,dms+[Lst])
                return HM(*(dms+[Lst]))
            else:
                raise ValueError("The range must be smaller then corresponding index range")
        elif indx == 'col':
            indx=1
            if len(Rg) <= A.n(indx):
                # Initialization of the hypermatrix which stores the result
                dms=A.dimensions(); dms[indx]=len(Rg); 
                # Initializing the list of entries
                Lst=[]
                # Initialization of the list specifying the dimensions of the output
                l = [A.n(i) for i in range(A.order())]
                # Main loop performing the transposition of the entries
                for i in range(prod(l)):
                    # Turning the index i into an hypermatrix array location using the decimal encoding trick
                    entry = [Integer(mod(i,l[0]))]
                    sm = Integer(mod(i,l[0]))
                    for k in range(len(l)-1):
                        entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                        sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                    if entry[indx] in Rg:
                        entry[indx] = Rg[Pmt[Rg.index(entry[indx])]]
                        Lst.append(A[tuple(entry)])
                #return apply(HM,dms+[Lst])
                return HM(*(dms+[Lst]))
            else:
                raise ValueError("The range must be smaller then corresponding index range")
        elif indx == 'dpt':
            indx=2
            if len(Rg) <= A.n(indx):
                # Initialization of the hypermatrix which stores the result
                dms=A.dimensions(); dms[indx]=len(Rg); 
                # Initializing the list of entries
                Lst=[]
                # Initialization of the list specifying the dimensions of the output
                l = [A.n(i) for i in range(A.order())]
                # Main loop performing the transposition of the entries
                for i in range(prod(l)):
                    # Turning the index i into an hypermatrix array location using the decimal encoding trick
                    entry = [Integer(mod(i,l[0]))]
                    sm = Integer(mod(i,l[0]))
                    for k in range(len(l)-1):
                        entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                        sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                    if entry[indx] in Rg:
                        entry[indx] = Rg[Pmt[Rg.index(entry[indx])]]
                        Lst.append(A[tuple(entry)])
                #return apply(HM,dms+[Lst])
                return HM(*(dms+[Lst]))
            else:
                raise ValueError("The range must be smaller then corresponding index range")
        elif indx == 'tme':
            indx=3
            if len(Rg) <= A.n(indx):
                # Initialization of the hypermatrix which stores the result
                dms=A.dimensions(); dms[indx]=len(Rg); 
                # Initializing the list of entries
                Lst=[]
                # Initialization of the list specifying the dimensions of the output
                l = [A.n(i) for i in range(A.order())]
                # Main loop performing the transposition of the entries
                for i in range(prod(l)):
                    # Turning the index i into an hypermatrix array location using the decimal encoding trick
                    entry = [Integer(mod(i,l[0]))]
                    sm = Integer(mod(i,l[0]))
                    for k in range(len(l)-1):
                        entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                        sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                    if entry[indx] in Rg:
                        entry[indx] = Rg[Pmt[Rg.index(entry[indx])]]
                        Lst.append(A[tuple(entry)])
                #return apply(HM,dms+[Lst])
                return HM(*(dms+[Lst]))
            else:
                raise ValueError("The range must be smaller then corresponding index range")
        else:
            raise ValueError("The string must be one of the following 4 choices row, col, dpt or tme")
    else:
        raise ValueError("The index must be a string or an integer")

def GeneralHypermatrixSlicerII(A, Rg, indx):
    """
    Outputs slices specified by index list L.
    the last string input is either row or col
    and determines the slices to be collected
    into a hypermatrix. The difference with
    the implementation above is that slice
    are placed in increasing order independently
    of the specification in the input list Rg


    EXAMPLES:
 
    ::

        sage: sz=3; A=HM(sz,sz,'a')
        sage: GeneralHypermatrixSlicerII(A, [0], 'row').printHM()
        [:, :]=
        [a00 a01 a02]
        sage: GeneralHypermatrixSlicerII(A, [1], 'col').printHM()
        [:, :]=
        [a01]
        [a11]
        [a21]
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of a number for the test
    nb=2
    if type(indx) == type(nb):
        if indx < A.order() :
            if len(Rg) <= A.n(indx):
                # Initialization of the hypermatrix which stores the result
                dms=A.dimensions(); dms[indx]=len(Rg); 
                # Initializing the list of entries
                Lst=[]
                # Initialization of the list specifying the dimensions of the output
                l = [A.n(i) for i in range(A.order())]
                # Main loop performing the transposition of the entries
                for i in range(prod(l)):
                    # Turning the index i into an hypermatrix array location using the decimal encoding trick
                    entry = [Integer(mod(i,l[0]))]
                    sm = Integer(mod(i,l[0]))
                    for k in range(len(l)-1):
                        entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                        sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                    if entry[indx] in Rg:
                        Lst.append(A[tuple(entry)])
                #return apply(HM,dms+[Lst])
                return HM(*(dms+[Lst]))
            else:
                raise ValueError("The range must be smaller then corresponding index range")
        else:
            raise ValueError("The index must be smaller then the order of the Hypermatrix")
    elif type(indx) == type('tst'):
        if indx == 'row':
            indx=0
            if len(Rg) <= A.n(indx):
                # Initialization of the hypermatrix which stores the result
                dms=A.dimensions(); dms[indx]=len(Rg); 
                # Initializing the list of entries
                Lst=[]
                # Initialization of the list specifying the dimensions of the output
                l = [A.n(i) for i in range(A.order())]
                # Main loop performing the transposition of the entries
                for i in range(prod(l)):
                    # Turning the index i into an hypermatrix array location using the decimal encoding trick
                    entry = [Integer(mod(i,l[0]))]
                    sm = Integer(mod(i,l[0]))
                    for k in range(len(l)-1):
                        entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                        sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                    if entry[indx] in Rg:
                        Lst.append(A[tuple(entry)])
                #return apply(HM,dms+[Lst])
                return HM(*(dms+[Lst]))
            else:
                raise ValueError("The range must be smaller then corresponding index range")
        elif indx == 'col':
            indx=1
            if len(Rg) <= A.n(indx):
                # Initialization of the hypermatrix which stores the result
                dms=A.dimensions(); dms[indx]=len(Rg); 
                # Initializing the list of entries
                Lst=[]
                # Initialization of the list specifying the dimensions of the output
                l = [A.n(i) for i in range(A.order())]
                # Main loop performing the transposition of the entries
                for i in range(prod(l)):
                    # Turning the index i into an hypermatrix array location using the decimal encoding trick
                    entry = [Integer(mod(i,l[0]))]
                    sm = Integer(mod(i,l[0]))
                    for k in range(len(l)-1):
                        entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                        sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                    if entry[indx] in Rg:
                        Lst.append(A[tuple(entry)])
                #return apply(HM,dms+[Lst])
                return HM(*(dms+[Lst]))
            else:
                raise ValueError("The range must be smaller then corresponding index range")
        elif indx == 'dpt':
            indx=2
            if len(Rg) <= A.n(indx):
                # Initialization of the hypermatrix which stores the result
                dms=A.dimensions(); dms[indx]=len(Rg); 
                # Initializing the list of entries
                Lst=[]
                # Initialization of the list specifying the dimensions of the output
                l = [A.n(i) for i in range(A.order())]
                # Main loop performing the transposition of the entries
                for i in range(prod(l)):
                    # Turning the index i into an hypermatrix array location using the decimal encoding trick
                    entry = [Integer(mod(i,l[0]))]
                    sm = Integer(mod(i,l[0]))
                    for k in range(len(l)-1):
                        entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                        sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                    if entry[indx] in Rg:
                        Lst.append(A[tuple(entry)])
                #return apply(HM,dms+[Lst])
                return HM(*(dms+[Lst]))
            else:
                raise ValueError("The range must be smaller then corresponding index range")
        elif indx == 'tme':
            indx=3
            if len(Rg) <= A.n(indx):
                # Initialization of the hypermatrix which stores the result
                dms=A.dimensions(); dms[indx]=len(Rg); 
                # Initializing the list of entries
                Lst=[]
                # Initialization of the list specifying the dimensions of the output
                l = [A.n(i) for i in range(A.order())]
                # Main loop performing the transposition of the entries
                for i in range(prod(l)):
                    # Turning the index i into an hypermatrix array location using the decimal encoding trick
                    entry = [Integer(mod(i,l[0]))]
                    sm = Integer(mod(i,l[0]))
                    for k in range(len(l)-1):
                        entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
                        sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
                    if entry[indx] in Rg:
                        Lst.append(A[tuple(entry)])
                #return apply(HM,dms+[Lst])
                return HM(*(dms+[Lst]))
            else:
                raise ValueError("The range must be smaller then corresponding index range")
        else:
            raise ValueError("The string must be one of the following 4 choices row, col, dpt or tme")
    else:
        raise ValueError("The index must be a string or an integer")

def KroneckerResultant(L, vrbl, VrbLp, VrbLq):
    """
    Takes as input a list of polynomials, a variable,
    two lists of dummy variables to be used in the linear
    combinations of constraints and outputs the corresponding
    Kronecker resultant list of polynomials.


    EXAMPLES:

    ::

        sage: sz=2;La=var_list('a',sz); Lb=var_list('b',sz); Lc=var_list('c',sz)
        sage: p=expand(prod(x-La[i] for i in rg(sz)))
        sage: q=expand(prod(x-La[i] for i in rg(floor(sz/2)))*prod(x-Lb[1] for i in rg(floor(sz/2))))
        sage: h=expand(prod((x-Lc[i]) for i in rg(sz)))
        sage: L=[p, q, h]
        sage: len(Set([factor(v) for v in KroneckerResultant(L, x, var_list('u',len(L)), var_list('v',len(L)))]).list())
        10        
        sage: X=var_list('x',4); L=[X[1] + X[2] + X[3] - 5, 7*X[1]*X[2] + 4*X[3] - 2*X[2] - 8, 10*X[2]+5*X[1]*X[2]-2*X[2]*X[3]+1]
        sage: len(Set([factor(v) for v in KroneckerResultant(L, X[1], var_list('u',len(L)), var_list('v',len(L)))]).list())
        6 


    AUTHORS:
    - Edinah K. Gnang, Jonathan Earl
    """
    # Initialization of the symbolic linear combination
    p=expand(sum(L[i]*VrbLp[i] for i in rg(len(L)))); q=expand(sum(L[i]*VrbLq[i] for i in rg(len(L))))
    # Initialization of the determinant
    f=expand(SylvesterHM(p, q, vrbl).det())
    # Initialization of the list of monomial
    Lm=[mn/mn.subs([v==1 for v in VrbLp+VrbLq]) for mn in f.operands()]
    return [f.coefficient(mn) for mn in Lm]

def KroneckerResultantII(L, vrbl, VrbLp, VrbLq):
    """
    Takes as input a list of polynomials, a variable,
    two lists of dummy variables to be used in the linear
    combinations of constraints and outputs the corresponding
    Kronecker resultant list of polynomials.


    EXAMPLES:

    ::

        sage: sz=2;La=var_list('a',sz); Lb=var_list('b',sz); Lc=var_list('c',sz)
        sage: p=expand(prod(x-La[i] for i in rg(sz)))
        sage: q=expand(prod(x-La[i] for i in rg(floor(sz/2)))*prod(x-Lb[1] for i in rg(floor(sz/2))))
        sage: h=expand(prod((x-Lc[i]) for i in rg(sz)))
        sage: L=[p, q, h]
        sage: len(Set([factor(v) for v in KroneckerResultant(L, x, var_list('u',len(L)), var_list('v',len(L)))]).list())
        10


    AUTHORS:
    - Jonathan Earl
    """
    lenL = len(L)
    var1, var2 = str(VrbLp[0])[0], str(VrbLq[0])[0]
    M = SylvesterHM(sum(VrbLp[i]*L[i] for i in range(lenL)), sum(VrbLq[i]*L[i] for i in range(lenL)), vrbl).matrix()
    f = expand(M.det())
    appender = {}
    for i in f.operands():
        parts = str(i).replace('-'+var1, '-1*'+var1).replace('-'+var2, '-1*'+var2).split('*')
        myKey = '*'.join([j for j in parts if j.startswith(var1) or j.startswith(var2)])
        myVal = '*'.join([j for j in parts if j not in myKey])
        myVal = 1 if not bool(myVal) else SR(myVal)
        appender[myKey] = appender.get(myKey, 0) + myVal
    return appender.values()

def kroneckerian_elimination(L, VrbL):
    """
    Takes as input a list of polynomials and a list of variable
    and outputs the corresponding resultant. Performs the Kroneckerian
    elimination algorithm

    EXAMPLES:

    ::

        sage: sz=3; VrbL=var_list('x',sz); L=(Vandermonde(rg(1,1+sz))*HM(sz,1,VrbL)-HM(sz,1,rg(sz))).list()
        sage: degree_matrix(L,VrbL)
        [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
        sage: len(kroneckerian_elimination(L, VrbL))
        3
        sage: sz=2; VrbL=var_list('x',sz); A=HM(sz,sz,'a'); b=HM(sz,1,var_list('b',sz))
        sage: L=(A*HM(sz,1,VrbL)-b).list()
        sage: kroneckerian_elimination(L, VrbL)
        [[a00*x0 + a01*x1 - b0, a10*x0 + a11*x1 - b1],
         [-a01*a10*x1 + a00*a11*x1 + a10*b0 - a00*b1,
          a01*a10*x1 - a00*a11*x1 - a10*b0 + a00*b1]]        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the resulting list
    RsL=[L]
    # Initialization of the loop eliminating all the variables
    for i in rg(len(VrbL)-1):
        RsL.append(Set(KroneckerResultant(RsL[len(RsL)-1], VrbL[i], var_list('u',len(RsL[len(RsL)-1])), var_list('v',len(RsL[len(RsL)-1])))).list()) 
    return RsL

def GeneralHypermatrixSubstituteInMatrix(A,vrbl,M):
    """
    Outputs a hypermatrix whose polynomial entries
    have been substited in the input matrix for the
    input variable.


    EXAMPLES:

    ::

        sage: x=var('x'); Ha=HM(2,1,[x+1,x^2+1])
        sage: Y=var_list('y',2); rM=GeneralHypermatrixSubstituteInMatrix(Ha,x,HM(2,2,[Y[0],0,0,Y[1]]))
        sage: rM[0,0].printHM()
        [:, :]=
        [y0 + 1      0]
        [     0 y1 + 1]
        sage: rM[1,0].printHM()
        [:, :]=
        [y0^2 + 1        0]
        [       0 y1^2 + 1]


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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        if A[tuple(entry)].is_zero():
            Rh[tuple(entry)] = 0
        else:
            Rh[tuple(entry)] = substituteHM(A[tuple(entry)],vrbl,M).expand()
    return Rh

def naught_eliminationHM(Cf):
    """
    Outputs the row echelon form of a multiplicative linear constraints where the RHS is zero.
    This implementation assumes that there is not division. This assumption incurs no loss of 
    generality at all since we can collect the denominators to make up a new multiplicative 
    system. Solve it independently and check whether they have non overlapping solutions.
    The corresponding problem is striking by its combinatorial flavor.
    

    EXAMPLES:
 
    ::

        sage: A = naught_eliminationHM(HM(2,2,'a'))
        sage: A.printHM()
        [:, :]=
        [a00 a01]
        [  0   0]
        sage: Ta=HM(2,2,'a') # Initialization of the coefficient matrix.
        sage: Ha=HM(2,2,[Ta[0,0]*HM(2,2,'kronecker'), Ta[1,0]*HM(2,2,'kronecker'), Ta[0,1]*HM(2,2,'kronecker'), Ta[1,1]*HM(2,2,'kronecker')])
        sage: A=naught_eliminationHM(Ha) # performing the gaussian elimination where entries are hypermatrices.
        sage: A
        [[[[a00, 0], [0, a00]], [[a01, 0], [0, a01]]], [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]]
        sage: sz=2; A=HM(sz,sz,'a'); B=HM(sz,sz,'b')
        sage: A00=HM([[A[0,0],-B[0,0]],[B[0,0],A[0,0]]]); A01=HM([[A[0,1],-B[0,1]],[B[0,1],A[0,1]]])
        sage: A10=HM([[A[1,0],-B[1,0]],[B[1,0],A[1,0]]]); A11=HM([[A[1,1],-B[1,1]],[B[1,1],A[1,1]]])
        sage: M=naught_eliminationHM(HM([[A00,A01],[A10,A11]]))
        sage: M
        [[[[a00, -b00], [b00, a00]], [[a01, -b01], [b01, a01]]], [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initializing a copy of the input second order hypermatrices.
    A=Cf.copy()
    # Initialization of the row and column index
    i=0; j=0
    while i < A.n(0) and j < A.n(1):
        while A.slice(rg(i,A.n(0)),'row').slice([j],'col').is_zero() and j < A.n(1)-1:
            # Incrementing the column index
            j=j+1
        if A.slice(rg(i,A.n(0)),'row').is_zero()==False:
            while A[i,j].is_zero(): 
                Ta=A.slice(rg(i,A.n(0)),'row')
                # Initializing the cyclic shift permutation matrix
                Id=HM(2,Ta.n(0), 'kronecker')
                P=Matrix2HM(sum([Id.matrix()[:,k]*Id.matrix()[Integer(mod(k+1,Ta.nrows())),:] for k in rg(Ta.n(0))]))
                Ta=P*Ta
                for i0 in rg(Ta.n(0)):
                    for j0 in rg(Ta.n(1)):
                        A[i+i0,j0]=Ta[i0,j0]
            # Performing the row operations.
            for r in rg(i+1,A.nrows()):
                # Taking care of the zero row
                if HM(1,A.n(1),[A[r,j0] for j0 in range(A.n(1))]).is_zero():
                    r=r+1
                else:
                    # Initialization of the coefficient
                    if A[r,j].is_zero() == False:
                        for j0 in rg(A.n(1)):
                            A[r,j0]=0*A[r,j0]
        # Incrementing the row and column index.
        i=i+1; j=j+1
    return A

def naught_reduced_eliminationHM(Cf):
    """
    Outputs the reduced row echelon form associated with the naught elimination.
    This implementation assumes that the input entries commute. 


    EXAMPLES:
 
    ::

        sage: RefA = naught_reduced_eliminationHM(HM(2,2,'a'))
        sage: RefA.printHM()
        [:, :]=
        [a00 a01]
        [  0   0]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    A=naught_eliminationHM(Cf)
    # Initialization of the row and column index
    i=A.nrows()-1; j=0
    while i>0 or j>0:
        #if (A[i,:]).is_zero():
        #if HM(1,A.n(1),[A[i,j0] for j0 in rg(A.n(1))]).is_zero():
        if A.slice([i],'row').is_zero():
            # decrementing the row index and initializing the column index
            i=i-1; j=0
        else :
            while (A[i,j]).is_zero():
                # Incrementing the column index
                j = j + 1
            # performing row operations
            for r in rg(i-1,-1,-1):
                #A[r,:] = -A[r,j]*A[i,:]+A[r,:]
                Tra=HM(1, A.n(1), 'zero')
                for j0 in rg(A.n(1)):
                    if j0 == j:
                        Tra[0,j0]=-A[r,j0]*A[i,j0]+A[r,j0]
                    else:
                        Tra[0,j0]=A[r,j0]
                for j0 in rg(A.n(1)):
                    A[r,j0]=Tra[0,j0]
            i=i-1; j=0
    return A

def default_naught_solver(Eq, La, Lf):
    """
    Formats the constraints performs and solves the multiplicatively linear constraints
    where the right hand side equals zero. This function outputs the solutions. The input
    EqL corresponds to a list of constraints. The input Lv corresponds to the list of 
    variables appearing in the constraints. The input Lf corresponds to the list of free
    varaibles each taken in correspondence with the entries of Lv. This implementation 
    tacitly assumes that the  the input constraints are indeed multiplicatively linear.
    This implementation performs cyclic permutations to the equations and all permutations
    to the variables in the first constraints.


    EXAMPLES:
 
    ::

        sage: sz=2; len(default_naught_solver([var('x'+str(i))*var('x'+str(sz+j)) for i in range(sz) for j in range(sz)], var_list('x', 2*sz), var_list('t', 2*sz)))
        4


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the list which stores the solutions
    Sln=[]
    for itr in rg(len(Eq)):
        # Performing a cyclic permutation of the constraints
        Eq.insert(0,Eq.pop())

        # Obtaining the variables in the first equations
        [TmpHa,hb]=multiplicativeConstraintFormatorIIHM(Eq[:1],La)
        La1=(TmpHa*HM(len(La),1,La))[0,0].operands()
        La2=Set(La).difference(Set(La1)).list()

        # Initializing the permutation of the variables in the first constraints
        P=Permutations(len(La1))
        for p in P:
            q=[p[i]-1 for i in rg(len(La1))]
            # Updating the ordering of the variables
            tLa=[La1[q[i]] for i in rg(len(La1))]+La2
            # Formatting the constraints to obtain the coefficient matrix
            [Ha, hb]=multiplicativeConstraintFormatorIIHM(Eq, tLa)

            # Performing the Gaussian elimination procedure
            tA=naught_reduced_eliminationHM(Ha)

            # Identifying the non zero rows of the matrix
            r=1
            while HM(tA.n(0),tA.n(1),'zero').fill_with(tA.slice(rg(r),'row')) != tA:
                r=r+1
        
            # Taking only the nonzero rows
            Ca=tA.slice(rg(r),'row')

            # Initialization of the vector
            #vA=HM(len(tLa),1,tLa)
        
            # Obtaining the resulting constraints
            #qE=(vA^Ca).list(); print 'Eq =', Eq,'qE =', qE, 'tLa =', tLa
        
            # Obtaining a solution to the system
            Mx=HM(Ca.n(1),1,tLa); Mv=HM(Ca.n(1),1,Lf) # Initialization of the pivot and free variables
            tmpSln=multiplicative_linear_solverHM(Ca,HM(Ca.n(0),1,'zero'),Mx,Mv)
            if (Set(tmpSln).list() in Sln) == False:
                Sln.append(Set(tmpSln).list())
    return Sln 

def naught_solver(EqL, La, Lf):
    """
    Formats the constraints performs and solves the multiplicatively linear constraints
    where the right hand side equals zero. This function outputs the solutions. The input
    EqL corresponds to a list of constraints. The input Lv corresponds to the list of 
    variables appearing in the constraints. The input Lf corresponds to the list of free
    varaibles each taken in correspondence with the entries of Lv. This implementation 
    tacitly assumes that the  the input constraints are indeed multiplicatively linear.
    This implementation performs all permutations of the variables and all permutations
    of the equations.


    EXAMPLES:
 
    ::

        sage: sz=2; len(naught_solver([var('x'+str(i))*var('x'+str(sz+j)) for i in range(sz) for j in range(sz)], var_list('x', 2*sz), var_list('t', 2*sz)))
        9


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Initialization of the permutations
    P=Permutations(len(La)); Q=Permutations(len(EqL))
    # Initialization of the list which stores the solutions
    Sln=[]
    for jtr in Q:
        itr=[jtr[i]-1 for i in rg(len(EqL))]
        # Performing a cyclic permutation of the constraints
        Eq=[EqL[itr[i]] for i in rg(len(EqL))]

        # Obtaining the variables in the first equations
        [TmpHa,hb]=multiplicativeConstraintFormatorIIHM(Eq[:1],La)
        La1=(TmpHa*HM(len(La),1,La))[0,0].operands()
        La2=Set(La).difference(Set(La1)).list()

        for p in P:
            q=[p[i]-1 for i in rg(len(La))]
            # Updating the ordering of the variables
            tLa=[La[q[i]] for i in rg(len(La))]
            # Formatting the constraints to obtain the coefficient matrix
            [Ha, hb]=multiplicativeConstraintFormatorIIHM(Eq, tLa)

            # Performing the Gaussian elimination procedure
            tA=naught_reduced_eliminationHM(Ha)

            # Identifying the non zero rows of the matrix
            r=1
            while HM(tA.n(0),tA.n(1),'zero').fill_with(tA.slice(rg(r),'row')) != tA:
                r=r+1
        
            # Taking only the nonzero rows
            Ca=tA.slice(rg(r),'row')

            # Initialization of the vector
            #vA=HM(len(tLa),1,tLa)
        
            # Obtaining the resulting constraints
            #qE=(vA^Ca).list(); print 'Eq =', Eq,'qE =', qE, 'tLa =', tLa
        
            # Obtaining a solution to the system
            Mx=HM(Ca.n(1),1,tLa); Mv=HM(Ca.n(1),1,Lf) # Initialization of the pivot and free variables
            tmpSln=multiplicative_linear_solverHM(Ca,HM(Ca.n(0),1,'zero'),Mx,Mv)
            if (Set(tmpSln).list() in Sln) == False:
                Sln.append(Set(tmpSln).list())
    return Sln
 
def SecondOrderIndexRotation(Ha, T):
    """
    The function perform the rotation of angle T for the indices.
    Ha is input second order hypermatrices. The rotation is performed
    clockwise by multiples of 2*pi/4.
    [i,j] -> [(i-floor(sz/2))*cos(T)+(-j+floor(sz/2))*sin(T)+floor(sz/2), (i-floor(sz/2))*sin(T)-(-j+floor(sz/2))*cos(T)+floor(sz/2)] if sz is odd
    [i,j] -> [(i-(sz-1)/2)*cos(T)+(-j+(sz-1)/2)*sin(T)+(sz-1)/2, (i-(sz-1)/2)*sin(T)-(-j+(sz-1)/2)*cos(T)+(sz-1)/2] if sz is even


    EXAMPLES:

    ::

        sage: sz=5; Ha=HM(sz,sz,'a') # Initialization of the input Hypermatrix
        sage: (Ha.tumble()-SecondOrderIndexRotation(Ha, 2*pi/4)).is_zero()
        True
        sage: sz=6; Ha=HM(sz,sz,'a') # Initialization of the input Hypermatrix
        sage: (Ha.tumble()-SecondOrderIndexRotation(Ha, 2*pi/4)).is_zero()
        True


    AUTHORS:
    - Edinah K. Gnang
    """
    if Ha.is_cubical() and Integer(mod(Ha.n(0),2)) == 1:
        # Initialization of the matrix
        sz=Ha.n(0); B=HM(sz,sz,'zero')
        for i in rg(sz):
            for j in rg(sz):
                B[i,j]=Ha[(i-floor(sz/2))*cos(T)+(-j+floor(sz/2))*sin(T)+floor(sz/2), (i-floor(sz/2))*sin(T)-(-j+floor(sz/2))*cos(T)+floor(sz/2)]
        return B
    elif Ha.is_cubical() and Integer(mod(Ha.n(0),2)) == 0:
        # Initialization of the matrix
        sz=Ha.n(0); B=HM(sz,sz,'zero')
        for i in rg(sz):
            for j in rg(sz):
                B[i,j]=Ha[(i-(sz-1)/2)*cos(T)+(-j+(sz-1)/2)*sin(T)+(sz-1)/2, (i-(sz-1)/2)*sin(T)-(-j+(sz-1)/2)*cos(T)+(sz-1)/2]
        return B
    else:
        raise ValueError("The input matrices must be square.")

def ThirdOrderIndexRotation(A, Langle):
    """
    The function takes a 3rd order hypermatrix and performs
    an index rotation around the axis specified in the Row,
    Column and Depth order specified by the input list angles.
    This implement only handles cubic hypermatrices. In case
    the hypermatrix is not cubic zeropadd to a cubic hypermatrix.

    EXAMPLES:
 
    ::

        sage: sz=2; A=HM(sz, sz, sz, 'a')
        sage: Langle=[2*pi/4, 0, 0]
        sage: ThirdOrderIndexRotation(A, Langle).printHM()
        [:, :, 0]=
        [a100 a110]
        [a101 a111]
        <BLANKLINE>
        [:, :, 1]=
        [a000 a010]
        [a001 a011]
        <BLANKLINE>


    AUTHORS:
    - Fan Tian and Edinah K. Gnang
    - To Do: Implement the arbitrary order version
    """
    if A.is_cubical():
        sz=A.n(0)
        # Initializing the output
        B=A.copy()
        # First performs a rotation around the Row ( or x ) axis performs a rotation of the column slices
        # Second performs a rotation around the Column ( or y ) axis performs a rotation of the  row slices
        # Third performs a rotation around the Depth ( or z ) axis performs a rotation of the depth slices
        axes = [1, 0, 2]
        for i in rg(len(axes)):
            axis = axes[i]
            for j in rg(sz):
                M=HM(sz, sz, B.slice([j], axis).list()).index_rotation(Langle[i])
                for u in rg(sz):
                    for v in rg(sz):
                        # Performing the index rotation relative to the row axis
                        if axis == 1:
                            B[u,j,v]=M[u,v]
                        # Performing the index rotation relative to the col axis
                        elif axis == 0:
                            B[j,u,v]=M[u,v]
                        # Performing the index rotation relative to the dpt axis
                        elif axis == 2:
                            B[u,v,j]=M[u,v]
        return B
    else:
        raise ValueError("Expected a cubic hypermatrix")

def SelectThirdOrderIndexRotation(A, Langle, EntryList):
    """
    The function takes a 3rd order hypermatrix and performs
    an index rotation  of the select indices specified by the 
    EntryList input around the axis specified in the Row,
    Column and Depth order specified by the input list angles.
    This implement only handles cubic hypermatrices. In case
    the hypermatrix is not cubic zeropadd to a cubic hypermatrix.


    EXAMPLES:
 
    ::

        sage: sz=2; A=HM(sz, sz, sz, 'a')
        sage: Langle=[2*pi/4, 0, 0]
        sage: SelectThirdOrderIndexRotation(A, Langle, rg(sz)).printHM()
        [:, :, 0]=
        [a100 a110]
        [a101 a111]
        <BLANKLINE>
        [:, :, 1]=
        [a000 a010]
        [a001 a011]
        <BLANKLINE>
        sage: sz=4; A=HM(sz, sz, sz, 'a')
        sage: Langle=[2*pi/4, 0, 0]
        sage: SelectThirdOrderIndexRotation(A, Langle, [0,1]).printHM()
        [:, :, 0]=
        [a100 a110 a020 a030]
        [a101 a111 a120 a130]
        [a200 a210 a220 a230]
        [a300 a310 a320 a330]
        <BLANKLINE>
        [:, :, 1]=
        [a000 a010 a021 a031]
        [a001 a011 a121 a131]
        [a201 a211 a221 a231]
        [a301 a311 a321 a331]
        <BLANKLINE>
        [:, :, 2]=
        [a002 a012 a022 a032]
        [a102 a112 a122 a132]
        [a202 a212 a222 a232]
        [a302 a312 a322 a332]
        <BLANKLINE>
        [:, :, 3]=
        [a003 a013 a023 a033]
        [a103 a113 a123 a133]
        [a203 a213 a223 a233]
        [a303 a313 a323 a333]
        <BLANKLINE>


    AUTHORS:
    - Edinah K. Gnang and Fan Tian 
    - To Do: Implement the arbitrary order version
    """
    if A.is_cubical():
        # Sorting the EntryList
        EntryList.sort()
        sz=A.n(0); TmpB=HM(len(EntryList),len(EntryList),len(EntryList),[A[i,j,k] for k in EntryList for j in EntryList for i in EntryList]).index_rotation(Langle)
        B=A.copy()
        for i in rg(len(EntryList)):
            for j in rg(len(EntryList)):
                for k in rg(len(EntryList)):
                    B[EntryList[i], EntryList[j], EntryList[k]]=TmpB[i,j,k]
        return B
    else:
        raise ValueError("Expected a cubic hypermatrix and the entry list must be smaller then the side length.")

def ThirdOrderSliceIndexRotation(A, Langle, Lslice):
    """
    The function takes a 3rd order hypermatrix and performs
    an index rotation around the axis specified in the Row,
    Column and Depth order specified by the input list angles.
    The operation is performed only to the slices which specified
    in the input Lslice.
    This implement only handles cubic hypermatrices. In case
    the hypermatrix is not cubic zeropadd to a cubic hypermatrix.

    EXAMPLES:
 
    ::

        sage: sz=2; A=HM(sz, sz, sz, 'a')
        sage: Langle=[0, 0, 2*pi/4]; Lslice=[0]
        sage: ThirdOrderSliceIndexRotation(A, Langle, Lslice).printHM()
        [:, :, 0]=
        [a100 a110]
        [a101 a111]
        <BLANKLINE>
        [:, :, 1]=
        [a000 a010]
        [a001 a011]
        <BLANKLINE>


    AUTHORS:
    - Fan Tian and Edinah K. Gnang
    - To Do: Implement the arbitrary order version
    """
    if A.is_cubical():
        sz=A.n(0)
        # Initializing the output
        B=A.copy()
        # First performs a rotation around the Row ( or x ) axis performs a rotation of the column slices
        # Second performs a rotation around the Column ( or y ) axis performs a rotation of the  row slices
        # Third performs a rotation around the Depth ( or z ) axis performs a rotation of the depth slices
        axes = [1, 0, 2]
        for i in rg(len(axes)):
            axis = axes[i]
            for j in Lslice:
                M=HM(sz, sz, B.slice([j], axis).list()).index_rotation(Langle[i])
                for u in rg(sz):
                    for v in rg(sz):
                        # Performing the index rotation relative to the row axis
                        if axis == 1:
                            B[u,j,v]=M[u,v]
                        # Performing the index rotation relative to the col axis
                        elif axis == 0:
                            B[j,u,v]=M[u,v]
                        # Performing the index rotation relative to the dpt axis
                        elif axis == 2:
                            B[u,v,j]=M[u,v]
        return B
    else:
        raise ValueError("Expected a cubic hypermatrix")

def SelectSecondOrderIndexRotation(Ha, T, EntryList):
    """
    The function perform the rotation of angle T for the 
    selected indices specified by the EntryList input.
    Ha is input second order hypermatrices.


    EXAMPLES:

    ::

        sage: sz=5; Ha=HM(sz,sz,'a') # Initialization of the input Hypermatrix
        sage: (Ha.tumble()-SelectSecondOrderIndexRotation(Ha, 2*pi/4, rg(sz))).is_zero()
        True
        sage: sz=6; Ha=HM(sz,sz,'a') # Initialization of the input Hypermatrix
        sage: (Ha.tumble()-SelectSecondOrderIndexRotation(Ha, 2*pi/4, rg(sz))).is_zero()
        True
        sage: sz=5; Ha=HM(sz,sz,'a')
        sage: (Ha-SelectSecondOrderIndexRotation(Ha, 2*pi/4, [0,1])).printHM()
        [:, :]=
        [ a00 - a10 -a00 + a01          0          0          0]
        [ a10 - a11 -a01 + a11          0          0          0]
        [         0          0          0          0          0]
        [         0          0          0          0          0]
        [         0          0          0          0          0]


    AUTHORS:
    - Edinah K. Gnang
    """
    if Ha.is_cubical() and len(EntryList) <= Ha.n(0):
        # Sorting the EntryList
        EntryList.sort()
        # Initialization of the matrix
        sz=Ha.n(0); TmpB=HM(len(EntryList),len(EntryList),[Ha[i,j] for j in EntryList for i in EntryList]).index_rotation(T)
        B=Ha.copy()
        for i in rg(len(EntryList)):
            for j in rg(len(EntryList)):
                B[EntryList[i], EntryList[j]]=TmpB[i,j]
        return B
    else:
        raise ValueError("The input matrices must be square and the entry list must be smaller then the side length.")

def SecondOrderIndexMap(A):
    """
    The function perform a very special index map to the the the indices.
    Ha is input second order hypermatrices and is assumed to be ?-diagonal


    EXAMPLES:

    ::

        sage: sz=5; Ha=HM(sz, sz, 'a') # Initialization of the input Hypermatrix
        sage: SecondOrderIndexMap(Ha).printHM()
        [:, :]=
        [      a00       a11       a22       a33       a44]
        [        0 a01 + a10 a12 + a21 a23 + a32 a34 + a43]
        [        0         0 a02 + a20 a13 + a31 a24 + a42]
        [        0         0         0 a03 + a30 a14 + a41]
        [        0         0         0         0 a04 + a40]


    AUTHORS:
    - Edinah K. Gnang
    """
    if A.is_cubical():
        # Initialization of the matrix
        sz=A.n(0); B=HM(sz,sz,'zero')
        for i in rg(sz):
            for j in rg(i+1):
                if i == A.n(0)-1: 
                    B[i,j]=A[i-j,A.n(1)-1-j]
                else:
                    B[i,j]=A[i-j,A.n(1)-1-j]+A[A.n(1)-1-j,i-j]
        return B.index_rotation(2*(2*pi/4))
    else:
        raise ValueError("The input matrices must be square.")

def geometric_mean(L):
    """
    The function computes the geometric mean
    of the input list L.


    EXAMPLES:

    ::

        sage: sz0=2; sz1=3; Ha=HM(sz0, sz1, 'a')
        sage: geometric_mean(Ha.list())
        (a00*a01*a02*a10*a11*a12)^(1/6)


    AUTHORS:
    - Edinah K. Gnang, Jeanine S. Gnang
    """
    # Initialization of the construct variable
    z=var('z')
    # Initialization of the size determined by the length of the list
    sz=len(L)
    # Return the geometric mean viewed as an inner-product of sorts
    return GProd([HM(1,sz,L), HM(sz,1,'one')], prod, [z])[0,0]^(1/sz)

def geometric_meanII(L):
    """
    The function computes the geometric mean
    of the input list L. To be used for numerical
    computation


    EXAMPLES:

    ::

        sage: X=[0.03658233053840871, 0.022693452948493835, 0.09667167574741581, 0.1354454280082533, 0.02906998361722122, 0.09068172967913375, 0.18808247869583106, 0.06596997458401005, 0.043486280252139436, 0.6087259605232003, 0.14746917844106938, 0.0024856908759630526, 0.11116106918595019, 0.02867537872396155, 0.15283808208604527]
        sage: geometric_meanII(X)
        0.0653181757705423


    AUTHORS:
    - Edinah K. Gnang, Jeanine S. Gnang
    """
    return prod(L)^(1.0/len(L))

def arithmetic_mean(L):
    """
    The function computes the arithmetic mean
    of the input list L.


    EXAMPLES:

    ::

        sage: sz0=2; sz1=3; Ha=HM(sz0, sz1, 'a')
        sage: arithmetic_mean(Ha.list())
        1/6*a00 + 1/6*a01 + 1/6*a02 + 1/6*a10 + 1/6*a11 + 1/6*a12


    AUTHORS:
    - Edinah K. Gnang, Jeanine S. Gnang
    """
    # Initialization of the construct variable
    z=var('z')
    # Initialization of the size determined by the length of the list
    sz=len(L)
    # Return the arithmetic mean viewed as an inner-product of sorts
    return GProd([HM(1,sz,L), HM(sz,1,'one')], sum, [z])[0,0]*(1/sz)

def arithmetic_meanII(L):
    """
    The function computes the arithmetic mean
    of the input list L. To be used for numerical
    computation


    EXAMPLES:

    ::

        sage: X=[0.03658233053840871, 0.022693452948493835, 0.09667167574741581, 0.1354454280082533, 0.02906998361722122, 0.09068172967913375, 0.18808247869583106, 0.06596997458401005, 0.043486280252139436, 0.6087259605232003, 0.14746917844106938, 0.0024856908759630526, 0.11116106918595019, 0.02867537872396155, 0.15283808208604527]
        sage: arithmetic_mean(X)
        0.117335912927140


    AUTHORS:
    - Edinah K. Gnang, Jeanine S. Gnang
    """
    return sum(L)*(1.0/len(L))

def GeneratePartition(sz):
    """
    Creates Partition to be used for computing the permanent.
    This function follows a maple implementation suggested
    Harry Crane


    EXAMPLES:
    ::


        sage: GeneratePartition(3)
        [[1, 1, 1], [1, 1, 2], [1, 2, 1], [1, 2, 2], [1, 2, 3]]
        


    AUTHORS:
    - Edinah K. Gnang, Harry Crane
    """
    N = [[1]]
    if sz > 1:
        for i in rg(1,sz):
            M = N
            r = len(M)
            N = []
            for j in rg(r):
                mx = max(M[j])
                for k in rg(1,mx+2):
                    N = N + [M[j]+[k]]
    return N

def Partition2HM(part):
    """
    Converts a partition into matrices. This function follows a maple implementation suggested
    Harry Crane


    EXAMPLES:
    ::


        sage: Partition2HM([1, 2, 1]).printHM()
        [:, :]=
        [1 0 1]
        [0 1 0]
        [1 0 1]


    AUTHORS:
    - Edinah K. Gnang, Harry Crane
    """
    M = max(part); N = len(part)
    B = HM(N,N,'zero')
    for i in rg(1,M+1):
        d = HM(N,1,'zero')
        for j in rg(N):
            if part[j] == i:
                d[j,0] = 1
        B = B + d*d.transpose()
    return B

def SetIntersection(L):
    """
    Outputs the intersection of the input list of Sets.
    This implementation does not check the validity of the inputs


    EXAMPLES:
 
    ::


        sage: SetIntersection([Set([1,2,3]),Set([1,2])])
        {1, 2}

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    St=L[0]
    for i in rg(1,len(L)):
        St=St.intersection(L[i])
    return St

def SetUnion(L):
    """
    Outputs the union of the input list of Sets.
    This implementation does not check the validity of the inputs


    EXAMPLES:
 
    ::


        sage: SetUnion([Set([1,2,3]),Set([1,2])])
        {1, 2, 3}

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    St=L[0]
    for i in rg(1,len(L)):
        St=St.union(L[i])
    return St

def DirectSum(L):
    """
    Outputs the direct sum of the input list of Hypermatrices.
    This implementation does not check the validity of the inputs


    EXAMPLES:
 
    ::


        sage: DirectSum([HM(2,2,'a'), HM(2,2,'b')]).printHM()
        [:, :]=
        [a00 a01   0   0]
        [a10 a11   0   0]
        [  0   0 b00 b01]
        [  0   0 b10 b11]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    TmpH=L[0]
    for i in rg(1,len(L)):
        TmpH=TmpH.block_sum(L[i])
    return TmpH

def TensorProduct(L):
    """
    Outputs the tensor product of the input list of Hypermatrices.
    This implementation does not check the validity of the inputs


    EXAMPLES:
 
    ::


        sage: TensorProduct([HM(2,2,'a'), HM(2,2,'b')]).printHM()
        [:, :]=
        [a00*b00 a00*b01 a01*b00 a01*b01]
        [a00*b10 a00*b11 a01*b10 a01*b11]
        [a10*b00 a10*b01 a11*b00 a11*b01]
        [a10*b10 a10*b11 a11*b10 a11*b11]


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    TmpH=L[0]
    for i in rg(1,len(L)):
        TmpH=TmpH.tensor_product(L[i])
    return TmpH

def Exp(L):
    """
    Outputs the exponentiation of the input list of two elements.
    The function checks that the list has only two elements
    When used with GProdIII it can only handle second order constructs


    EXAMPLES:
 
    ::


        sage: Exp(var_list('a',2))
        a0^a1


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    if len(L) == 2:
        return L[0]^L[1]
    else:
        raise ValueError("Expected list of two elements")

def BaseExp(L):
    """
    Outputs the base exponentiation of the input list
    This implementation check the validity of the inputs
    When used with GProdIII it can only handle second order constructs


    EXAMPLES:
 
    ::


        sage: Exp(var_list('a',2))
        a0^a1


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    if len(L) == 2:
        return L[1]^L[0]
    else:
        raise ValueError("Expected list of two elements")

def ExpN(L, dgts=50):
    """
    Outputs the exponentiation of the input list of two elements
    The function checks that the list has only two elements


    EXAMPLES:
 
    ::


        sage: ExpN([2.000000, 3.000000])
        8.0000000000000


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    if len(L) == 2:
        return ComplexField(dgts)(exp(L[1]*ln(L[0])))
    else:
        raise ValueError("Expected list of two elements")

def BaseExpN(L,dgts=50):
    """
    Outputs the base exponentiation of the input list
    This implementation check the validity of the inputs


    EXAMPLES:
 
    ::


        sage: BaseExpN([2.000000, 3.000000])
        9.0000000000000


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    if len(L) == 2:
        return ComplexField(dgts)(exp(L[0]*ln(L[1])))
    else:
        raise ValueError("Expected list of two elements")

def Sigmoid(L):
    """
    Outputs the sgmoid of the input list
    This implementation check the validity of the inputs
    When used with GProdIII it can only handle second order constructs


    EXAMPLES:
 
    ::


        sage: Exp(var_list('a',2))
        a0^a1


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    if len(L) == 2:
        return L[1]^L[0]/(L[1]^L[0]+1)
    else:
        raise ValueError("Expected list of two elements")

def And(L):
    """
    Outputs the conjuction of the input list of boolean values.
    This implementation does not check the validity of the inputs


    EXAMPLES:
 
    ::


        sage: And([True, False])
        False

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    St=L[0]
    for i in rg(1,len(L)):
        St = St and L[i]
    return St

def Or(L):
    """
    Outputs the disjunction of the input list of boolean values.
    This implementation does not check the validity of the inputs


    EXAMPLES:
 
    ::


        sage: Or([False, True])
        True

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    St=L[0]
    for i in rg(1,len(L)):
        St = St or L[i]
    return St

def MinSum(L):
    """
    Outputs the minimum between the sum of entries of the input
    list L and the unit integer 1.
    This implementation does not check the validity of the inputs


    EXAMPLES:
 
    ::


        sage: MinSum([2,3])
        1
        sage: MinSum([2, -3, 1])
        0

    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    return min(1, sum(L))

def is_Tree(A):
    """
    Returns an boolean value determining if the input unweighted
    adjacency matrix is associated with a tree. The implementation
    is based on a implementation of the matrix tree theorem.


    EXAMPLES:
    ::
        sage: is_Tree(Matrix([[0, 1, 0], [1, 0, 1], [0, 1, 0]]))
        True

    AUTHORS:
    - Edinah K. Gnang
    """
    if 1==((diagonal_matrix((A*ones_matrix(A.nrows(),1)).list())-A)[0:A.nrows()-1,0:A.ncols()-1]).det():
        return True
    else:
        return False

def Tuple2DiGraph(T,sz):
    """
    The method returns a directed graph object associated with 
    with the tuple list description of the directed graph

    EXAMPLES:

    ::

        sage: Tuple2DiGraph([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],5).degree_sequence()
        [4, 2, 2, 1, 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=HM(2,sz,'kronecker')
    return DiGraph(sum([Id.slice([t[0]],'col')*Id.slice([t[1]],'row') for t in T]).matrix())

def TupleFunctionList(sz):
    """
    Returns a list of edge tuple desctiption for all 
    functional directed graphs on sz vertices.


    EXAMPLES:
    ::
        sage: TupleFunctionList(3)
        [[(0, 0), (1, 0), (2, 0)],
         [(0, 1), (1, 0), (2, 0)],
         [(0, 2), (1, 0), (2, 0)],
         [(0, 0), (1, 1), (2, 0)],
         [(0, 1), (1, 1), (2, 0)],
         [(0, 2), (1, 1), (2, 0)],
         [(0, 0), (1, 2), (2, 0)],
         [(0, 1), (1, 2), (2, 0)],
         [(0, 2), (1, 2), (2, 0)],
         [(0, 0), (1, 0), (2, 1)],
         [(0, 1), (1, 0), (2, 1)],
         [(0, 2), (1, 0), (2, 1)],
         [(0, 0), (1, 1), (2, 1)],
         [(0, 1), (1, 1), (2, 1)],
         [(0, 2), (1, 1), (2, 1)],
         [(0, 0), (1, 2), (2, 1)],
         [(0, 1), (1, 2), (2, 1)],
         [(0, 2), (1, 2), (2, 1)],
         [(0, 0), (1, 0), (2, 2)],
         [(0, 1), (1, 0), (2, 2)],
         [(0, 2), (1, 0), (2, 2)],
         [(0, 0), (1, 1), (2, 2)],
         [(0, 1), (1, 1), (2, 2)],
         [(0, 2), (1, 1), (2, 2)],
         [(0, 0), (1, 2), (2, 2)],
         [(0, 1), (1, 2), (2, 2)],
         [(0, 2), (1, 2), (2, 2)]]



    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Appending the function to the list
        Lf.append([(i,f[i]) for i in rg(sz)])
    return Lf

def RepresentativeTupleFunctionList(sz):
    """
    Returns a list of edge tuple desctiption for all 
    functional directed graphs. Outputing one per
    isomorphism class the graph isomorphism routine
    is doing the heavy lifting here.


    EXAMPLES:
    ::
        sage: RepresentativeTupleFunctionList(3)
        [[(0, 0), (1, 0), (2, 0)],
         [(0, 1), (1, 0), (2, 0)],
         [(0, 0), (1, 1), (2, 0)],
         [(0, 1), (1, 1), (2, 0)],
         [(0, 2), (1, 1), (2, 0)],
         [(0, 1), (1, 2), (2, 0)],
         [(0, 0), (1, 1), (2, 2)]]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the lists
    L=TupleFunctionList(sz)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for tp in L:
        nwT=True
        for i in range(len(cL)):
            if Tuple2DiGraph(tp,sz).is_isomorphic(Tuple2DiGraph(cL[i],sz)):
                nwT=False
                break
        if nwT==True:
            cL.append(tp)
    return cL

def PermutationFunctionList(sz):
    """
    Returns a list of edge tuple descriptions associated 
    with permutations.


    EXAMPLES:
    ::
        sage: PermutationFunctionList(3)
        [[(0, 0), (1, 1), (2, 2)],
         [(0, 0), (1, 2), (2, 1)],
         [(0, 1), (1, 0), (2, 2)],
         [(0, 1), (1, 2), (2, 0)],
         [(0, 2), (1, 0), (2, 1)],
         [(0, 2), (1, 1), (2, 0)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list of permutations of elements from 1 to (n-1).
    P=Permutations(sz)
    return [[(i,p[i]-1) for i in rg(sz)] for p in P]

def RepresentativePermutationFunctionList(sz):
    """
    Returns a list of edge tuple descriptions associated 
    with permutations.


    EXAMPLES:
    ::
        sage: RepresentativePermutationFunctionList(3)
        [[(0, 0), (1, 1), (2, 2)], [(0, 0), (1, 2), (2, 1)], [(0, 1), (1, 2), (2, 0)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    L=PermutationFunctionList(sz)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for tp in L:
        nwT=True
        for i in range(len(cL)):
            if Tuple2DiGraph(tp,sz).is_isomorphic(Tuple2DiGraph(cL[i],sz)):
                nwT=False
                break
        if nwT==True:
            cL.append(tp)
    return cL

def is_permutation(T):
    """
    Dertermines whether or not the input tuple describes a permutation


    EXAMPLES:

    ::


        sage: is_permutation([(0, 0), (1, 2), (2, 3), (3, 4), (4, 5), (5, 1), (6, 7), (7, 6)])
        True


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the adjacency matrix
    sz=len(T); A=HM(sz,sz,'zero')
    for t in T:
        A[t[0],t[1]]=1
    if (A*HM(sz,1,'one')-HM(sz,1,'one')).is_zero() and (A.transpose()*HM(sz,1,'one')-HM(sz,1,'one')).is_zero():
        return True
    else:
        return False

def RandomPermutationTuple(sz):
    """
    Outputs a permutation described in tuple notation
    chosen uniformly at random


    EXAMPLES:

    ::

        sage: RandomPermutationTuple(1)
        [(0, 0)]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the Tuple and the list which
    # will store values that were already used
    L=rg(sz); T=[]
    for i in rg(sz):
        T.append((i, L[randint(0,len(L)-1)]))
        L.remove(T[len(T)-1][1])
    return T

def RandomDerangementTuple(sz):
    """
    Outputs a derangement described in tuple notation
    chosen uniformly at random


    EXAMPLES:

    ::

        sage: RandomDerangementTuple(2)
        [(0, 1), (1, 0)]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the randomly chosen permutation
    T = RandomPermutationTuple(sz)
    while prod(T[i][1]-T[i][0] for i in rg(sz)).is_zero():
        T = RandomPermutationTuple(sz) 
    return T

def is_functional_tree(T):
    """
    Dertermines whether or not the input tuple describes a permutation


    EXAMPLES:

    ::


        sage: is_functional_tree([(0, 0), (1, 2), (2, 3), (3, 4), (4, 5), (5, 1), (6, 7), (7, 6)])
        False
        sage: is_functional_tree([(0, 3), (1, 3), (2, 3), (3, 3)])
        True


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the adjacency matrix
    sz=len(T); A=HM(sz,sz,'zero')
    for t in T:
        A[t[0],t[1]]=1
    if 1==sum(A[i,i]*((HM(2,(A*HM(sz,1,'one')).list(),'diag')-A).slice([j for j in rg(sz) if j!=i],'row').slice([k for k in rg(sz) if k!=i],'col')).det() for i in rg(sz)):
        return True
    else:
        return False

def RootedTupleTreeFunctionList(sz):
    """
    Goes through all the functions and determines which ones
    are associated with trees.

    EXAMPLES:
    ::
        sage: RootedTupleTreeFunctionList(3)
        [[(0, 0), (1, 0), (2, 0)],
         [(0, 1), (1, 1), (2, 0)],
         [(0, 0), (1, 2), (2, 0)],
         [(0, 0), (1, 0), (2, 1)],
         [(0, 1), (1, 1), (2, 1)],
         [(0, 2), (1, 1), (2, 1)],
         [(0, 2), (1, 0), (2, 2)],
         [(0, 1), (1, 2), (2, 2)],
         [(0, 2), (1, 2), (2, 2)]]
        sage: sz=3; X=var_list('x',sz+1); Y=var_list('y',sz+1); Z=var_list('z',sz+1); Ta=HM(sz+1,sz+1,'a'); A=HM(sz+1,sz+1,'zero'); DgA=[Ta[i,i] for i in rg(sz+1)]
        sage: for i in rg(sz+1):
        ....:     for j in rg(sz+1):
        ....:         A[i,j]=Ta[i,j]*X[abs(j-i)]*prod(Y[u] for u in rg(min(i,j), max(i,j)))*prod(Z[u] for u in rg(1+min(i,j), 1+max(i,j)))
        ....:
        sage: sum(prod(A[t[0],t[1]] for t in tp) for tp in RootedTupleTreeFunctionList(sz))
        a00*a10*a20*x0*x1*x2*y0^2*y1*z1^2*z2 + a01*a11*a20*x0*x1*x2*y0^2*y1*z1^2*z2 + a02*a10*a22*x0*x1*x2*y0^2*y1*z1^2*z2 + a00*a12*a20*x0*x1*x2*y0*y1^2*z1*z2^2 + a02*a11*a21*x0*x1*x2*y0*y1^2*z1*z2^2 + a02*a12*a22*x0*x1*x2*y0*y1^2*z1*z2^2 + a00*a10*a21*x0*x1^2*y0*y1*z1*z2 + a01*a11*a21*x0*x1^2*y0*y1*z1*z2 + a01*a12*a22*x0*x1^2*y0*y1*z1*z2
        sage: sz=5; X=HM(sz,1,var_list('x',sz)); Y=HM(1,sz,var_list('y',sz)); A=X*Y # Quick demonstration of Cayley's Theorem
        sage: factor(sum(prod(A[t[0],t[1]] for t in tp) for tp in RootedTupleTreeFunctionList(sz)))
        (y0^2 + y1^2 + y2^2 + y3^2 + y4^2)*x0*x1*x2*x3*x4*(y0 + y1 + y2 + y3 + y4)^3
        sage: sz=5; A=HM(sz,sz,'a'); F=sum(prod(A[t[0],t[1]] for t in tp) for tp in RootedTupleTreeFunctionList(sz))


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Initialization of the adjacency martrix of the functional graph
        A = sum([Id[:,j]*Id[f[j],:] for j in range(sz)])
        # Initialization of the directed Laplacian
        lA= diagonal_matrix((A*ones_matrix(A.nrows(),1)).list())-A
        # Initialization of the list of sumbratrices
        Lmtr=[Matrix(ZZ,sz-1,sz-1,[lA[u,v] for u in range(sz) for v in range(sz) if u!=t if v!=t]) for t in range(sz)]
        # Testing treeness
        if sum(A[t,t]*Lmtr[t].det() for t in range(sz)) == 1:
            # Appending the function to the list
            Lf.append([(i, f[i]) for i in range(sz)])
    return Lf

def RepresentativeRootedTupleTreeFunctionList(sz):
    """
    Goes through all the functions and determines which ones
    are associated with trees. The method returns only one
    per isomorphism equivalence class


    EXAMPLES:
    ::
        sage: RepresentativeRootedTupleTreeFunctionList(3)
        [[(0, 0), (1, 0), (2, 0)], [(0, 1), (1, 1), (2, 0)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the list
    L = RootedTupleTreeFunctionList(sz)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for tp in L:
        nwT=True
        for i in range(len(cL)):
            if Tuple2DiGraph(tp,sz).is_isomorphic(Tuple2DiGraph(cL[i],sz)):
                nwT=False
                break
        if nwT==True:
            cL.append(tp)
    return cL

def RootedTupleInducedTreeFunctionList(sz, induced_edge_label_sequence):
    """
    Goes through all the functions and determines which ones
    are associated with a given edge labeled sequence sorted
    in non decreasing order

    EXAMPLES:
    ::
        sage: RootedTupleInducedTreeFunctionList(3, [0, 1, 2])
        [[(0, 0), (1, 0), (2, 0)],
         [(0, 1), (1, 1), (2, 0)],
         [(0, 0), (1, 2), (2, 0)],
         [(0, 2), (1, 1), (2, 1)],
         [(0, 2), (1, 0), (2, 2)],
         [(0, 2), (1, 2), (2, 2)]]
        sage: sz=3; X=var_list('x',sz+1); Y=var_list('y',sz+1); Z=var_list('z',sz+1); Ta=HM(sz+1,sz+1,'a'); A=HM(sz+1,sz+1,'zero'); DgA=[Ta[i,i] for i in rg(sz+1)]
        sage: for i in rg(sz+1):
        ....:     for j in rg(sz+1):
        ....:         A[i,j]=Ta[i,j]*X[abs(j-i)]*prod(Y[u] for u in rg(min(i,j), max(i,j)))*prod(Z[u] for u in rg(1+min(i,j), 1+max(i,j)))
        ....:
        sage: sum(prod(A[t[0],t[1]] for t in tp) for tp in RootedTupleInducedTreeFunctionList(sz,[0,1,2]))
        a00*a10*a20*x0*x1*x2*y0^2*y1*z1^2*z2 + a01*a11*a20*x0*x1*x2*y0^2*y1*z1^2*z2 + a02*a10*a22*x0*x1*x2*y0^2*y1*z1^2*z2 + a00*a12*a20*x0*x1*x2*y0*y1^2*z1*z2^2 + a02*a11*a21*x0*x1*x2*y0*y1^2*z1*z2^2 + a02*a12*a22*x0*x1*x2*y0*y1^2*z1*z2^2


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Initialization of the adjacency martrix of the functional graph
        A = sum([Id[:,j]*Id[f[j],:] for j in range(sz)])
        # Initialization of the directed Laplacian
        lA= diagonal_matrix((A*ones_matrix(A.nrows(),1)).list())-A
        # Initialization of the list of sumbratrices
        Lmtr=[Matrix(ZZ,sz-1,sz-1,[lA[u,v] for u in range(sz) for v in range(sz) if u!=t if v!=t]) for t in range(sz)]
        # Testing treeness
        EdgLblSeq=[abs(f[i]-i) for i in range(sz)]; EdgLblSeq.sort()
        if sum(A[t,t]*Lmtr[t].det() for t in range(sz)) == 1 and  EdgLblSeq == induced_edge_label_sequence:
            # Appending the function to the list
            Lf.append([(i,f[i]) for i in range(sz)])
    return Lf

def NonDecreasingFunctionList(sz):
    """
    Goes through all the functions and determines which ones
    are associated with pointwise non decreasing functions.
    They capture the isomorphism class of all spanning union
    of functional trees.

    EXAMPLES:
    ::
        sage: NonDecreasingFunctionList(3)[0]
        [[(0, 0), (1, 1), (2, 2)],
         [(0, 1), (1, 1), (2, 2)],
         [(0, 2), (1, 1), (2, 2)],
         [(0, 0), (1, 2), (2, 2)],
         [(0, 1), (1, 2), (2, 2)],
         [(0, 2), (1, 2), (2, 2)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]; Lg=[]
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Testing whether the candidate fucntion is non decreasing 
        # Setting the boolean function
        Decreasing = True
        for i in rg(sz):
            if f[i] < i:
                Decreasing = False
                break
        # Appending the function to the list
        if Decreasing == True:
            Lf.append([(i, f[i]) for i in range(sz)])
        else:
            Lg.append([(i, f[i]) for i in range(sz)])
    return [Lf, Lg]

def DecreasingFunctionList(sz):
    """
    Goes through all the functions and determines which ones
    are associated with pointwise decreasing functions.

    EXAMPLES:
    ::
        sage: DecreasingFunctionList(3)
        [[(0, 0), (1, 0), (2, 0)], [(0, 0), (1, 0), (2, 1)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    A=HM(sz,sz,'a')
    if sz == 1:
        return [(0, 0)]
    elif sz == 2:
        return [(0, 0), (1, 0)]
    elif sz > 2:
        return [Monomial2Tuple(mnm, A.list(), sz) for mnm in expand(A[0,0]*prod(sum(A[i,j] for j in rg(i)) for i in rg(1,sz))).operands()]
    else:
        raise ValueError("Expected the size to be an integer >=1")

def NonIncreasingFunctionList(sz):
    """
    Goes through all the functions and determines which ones
    are associated with non increasing functions.
    They capture the isomorphism class of all spanning
    unions of functional trees.

    EXAMPLES:
    ::
        sage: NonIncreasingFunctionList(3)[0]
        [[(0, 0), (1, 0), (2, 0)],
         [(0, 0), (1, 1), (2, 0)],
         [(0, 0), (1, 0), (2, 1)],
         [(0, 0), (1, 1), (2, 1)],
         [(0, 0), (1, 0), (2, 2)],
         [(0, 0), (1, 1), (2, 2)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]; Lg=[]
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Testing whether the candidate fucntion is non decreasing 
        # Setting the boolean function
        Increasing = True
        for i in rg(sz):
            if f[i] > i:
                Increasing = False
                break
        # Appending the function to the list
        if Increasing == True:
            Lf.append([(i, f[i]) for i in range(sz)])
        else:
            Lg.append([(i, f[i]) for i in range(sz)])
    return [Lf, Lg]

def TupleComplementaryLabel(T):
    """
    Returns the tuple encoding of the involuted labeling
    This implementation does not assume that the tuple is rooted at 0.

    EXAMPLES:

    ::

        sage: T=[(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)]
        sage: TupleComplementaryLabel(T)
        [(0, 4), (1, 4), (2, 0), (3, 2), (4, 4)]


    AUTHORS:
    - Edinah K. Gnang
    """
    sz=len(T)
    tp=[(sz-1-t[0], sz-1-t[1]) for t in T]
    tp.sort()
    return tp

def IncreasingFunctionList(sz):
    """
    Goes through all the functions and determines which ones
    are associated with pointwise increasing functions.

    EXAMPLES:
    ::
        sage: IncreasingFunctionList(3)
        [[(0, 1), (1, 2), (2, 2)], [(0, 2), (1, 2), (2, 2)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    A=HM(sz,sz,'a')
    if sz == 1:
        return [(0, 0)]
    elif sz == 2:
        return [(0, 1), (1, 1)]
    elif sz > 2:
        return [Monomial2Tuple(mnm, A.list(), sz) for mnm in expand(A[sz-1,sz-1]*prod(sum(A[sz-1-i,sz-1-j] for j in rg(i)) for i in rg(1,sz))).operands()]
    else:
        raise ValueError("Expected the size to be an integer >=1")

def RootedTupleTreeNonIncreasingFunctionList(sz):
    """
    Goes through all the functions and determines which ones
    are associated with trees and separates then according
    to whether or they are non increasing

    EXAMPLES:
    ::
        sage: RootedTupleTreeNonIncreasingFunctionList(3)[0]
        [[(0, 0), (1, 0), (2, 0)], [(0, 0), (1, 0), (2, 1)]] 



    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]; Lg=[]
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Initialization of the adjacency martrix of the functional graph
        A = sum([Id[:,j]*Id[f[j],:] for j in range(sz)])
        # Initialization of the dierected Laplacian
        lA= diagonal_matrix((A*ones_matrix(A.nrows(),1)).list())-A
        # Initialization of the list of sumbratrices
        Lmtr=[Matrix(ZZ,sz-1,sz-1,[lA[u,v] for u in range(sz) for v in range(sz) if u!=t if v!=t]) for t in range(sz)]
        # Testing treeness
        if sum(A[t,t]*Lmtr[t].det() for t in range(sz)) == 1:
            # Setting the boolean function
            Increasing = True
            for i in rg(sz):
                if f[i] > i:
                    Increasing = False
                    break
            # Appending the function to the list
            if Increasing == True:
                Lf.append([(i, f[i]) for i in range(sz)])
            else:
                Lg.append([(i, f[i]) for i in range(sz)])
    return [Lf, Lg]

def RootedTupleTreeNonDecreasingFunctionList(sz):
    """
    Goes through all the functions and determines which ones
    are associated with non decreasing functional trees.

    EXAMPLES:
    ::
        sage: RootedTupleTreeNonDecreasingFunctionList(3)[0]
        [[(0, 1), (1, 2), (2, 2)], [(0, 2), (1, 2), (2, 2)]]



    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[sz for i in range(sz)]; Lf=[]; Lg=[]
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Initialization of the adjacency martrix of the functional graph
        A = sum([Id[:,j]*Id[f[j],:] for j in range(sz)])
        # Initialization of the dierected Laplacian
        lA= diagonal_matrix((A*ones_matrix(A.nrows(),1)).list())-A
        # Initialization of the list of sumbratrices
        Lmtr=[Matrix(ZZ,sz-1,sz-1,[lA[u,v] for u in range(sz) for v in range(sz) if u!=t if v!=t]) for t in range(sz)]
        # Testing treeness
        if sum(A[t,t]*Lmtr[t].det() for t in range(sz)) == 1:
            # Setting the boolean function
            Decreasing = True
            for i in rg(sz):
                if f[i] < i:
                    Decreasing = False
                    break
            # Appending the function to the list
            if Decreasing == True:
                Lf.append([(i, f[i]) for i in range(sz)])
            else:
                Lg.append([(i, f[i]) for i in range(sz)])
    return [Lf, Lg]

def RepresentativeTupleSpanningFunctionalTreeList(sz):
    """
    Returns a list of edge tuple desctiption for all 
    spanning unions of functional trees. Outputing one
    per isomorphism class the graph isomorphism routine
    is doing the heavy lifting here. The number of
    vertices size must be at least 2.


    EXAMPLES:
    ::
        sage: RepresentativeTupleSpanningFunctionalTreeList(4)
        [[(0, 0), (1, 0), (2, 0), (3, 0)],
         [(0, 0), (1, 1), (2, 0), (3, 0)],
         [(0, 0), (1, 0), (2, 1), (3, 0)],
         [(0, 0), (1, 1), (2, 1), (3, 0)],
         [(0, 0), (1, 1), (2, 2), (3, 0)],
         [(0, 0), (1, 0), (2, 1), (3, 1)],
         [(0, 0), (1, 0), (2, 2), (3, 1)],
         [(0, 0), (1, 0), (2, 1), (3, 2)],
         [(0, 0), (1, 1), (2, 2), (3, 3)]]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the lists
    A=HM(sz,sz,'a')
    L=[Monomial2Tuple(mnm, A.list(), sz) for mnm in expand(A[0,0]*prod(sum(A[i,j] for j in rg(i+1)) for i in rg(1,sz))).operands()]
    # Initialization of the list storing the equivalence class of trees.
    cL=[ L[0] ]
    # Loop perfomring the isomorphism binning.
    for tp in L:
        nwT=True
        for i in range(len(cL)):
            if Tuple2DiGraph(tp,sz).is_isomorphic(Tuple2DiGraph(cL[i],sz)):
                nwT=False
                break
        if nwT==True:
            cL.append(tp)
    return cL

def TupleSpanningFunctionalTreeList(sz):
    """
    Returns a list of edge tuple desctiption for all 
    spanning unions of functional trees. There are
    (sz+1)^(sz-1) of them. 



    EXAMPLES:
    ::
        sage: TupleSpanningFunctionalTreeList(3)
        [[(0, 0), (1, 0), (2, 0)],
         [(0, 0), (1, 1), (2, 0)],
         [(0, 1), (1, 1), (2, 0)],
         [(0, 0), (1, 2), (2, 0)],
         [(0, 0), (1, 0), (2, 1)],
         [(0, 0), (1, 1), (2, 1)],
         [(0, 1), (1, 1), (2, 1)],
         [(0, 2), (1, 1), (2, 1)],
         [(0, 0), (1, 0), (2, 2)],
         [(0, 2), (1, 0), (2, 2)],
         [(0, 0), (1, 1), (2, 2)],
         [(0, 1), (1, 1), (2, 2)],
         [(0, 2), (1, 1), (2, 2)],
         [(0, 0), (1, 2), (2, 2)],
         [(0, 1), (1, 2), (2, 2)],
         [(0, 2), (1, 2), (2, 2)]] 


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the lists
    A=HM(sz, sz, 'a')
    L =TupleFunctionList(sz)
    # Initialization of the list storing the equivalence class of trees.
    cL=[ ]
    cL0=RepresentativeTupleSpanningFunctionalTreeList(sz)
    # Loop perfomring the isomorphism binning.
    for tp in L:
        nwT=False
        for tq in cL0:
            if Tuple2DiGraph(tp,sz).is_isomorphic(Tuple2DiGraph(tq,sz)):
                nwT=True
                break
        if nwT==True:
            cL.append(tp)
    return cL

def TupleSpanningFunctionalTreeListII(sz):
    """
    Returns a list of edge tuple desctiption for all 
    spanning unions of functional trees. There are
    (sz+1)^(sz-1) of them. This implementation is based
    on the matrix tree theorem proof argument using
    functional trees on sz+1 vertices rooted at the
    vertex labeled sz. The listing of spanning unions
    of functional trees on sz vertices is obtained by
    contracting every subgraph described by the edge
    monomial A[i,sz]*A[sz,sz] into the self loop edge
    A[i,i] where i is in rg(sz).


    EXAMPLES:
    ::
        sage: TupleSpanningFunctionalTreeListII(3)
        [[(0, 0), (1, 0), (2, 0)],
         [(0, 0), (1, 1), (2, 0)],
         [(0, 1), (1, 1), (2, 0)],
         [(0, 0), (1, 2), (2, 0)],
         [(0, 0), (1, 0), (2, 1)],
         [(0, 0), (1, 1), (2, 1)],
         [(0, 1), (1, 1), (2, 1)],
         [(0, 2), (1, 1), (2, 1)],
         [(0, 0), (1, 0), (2, 2)],
         [(0, 2), (1, 0), (2, 2)],
         [(0, 0), (1, 1), (2, 2)],
         [(0, 1), (1, 1), (2, 2)],
         [(0, 2), (1, 1), (2, 2)],
         [(0, 0), (1, 2), (2, 2)],
         [(0, 1), (1, 2), (2, 2)],
         [(0, 2), (1, 2), (2, 2)]] 


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the lists
    A=HM(sz+1, sz+1, 'a'); As=HM(sz, sz, 'a')
    # Initialization of the directed Laplacian
    LpmA=HM(2,(A*HM(sz+1,1,'one')).list(),'diag')-A
    # Initialization of the symbolic listing
    tmpF=expand(Deter(LpmA.slice(rg(sz),'row').slice(rg(sz),'col')))
    F=fast_reduce(tmpF, [A[i,sz] for i in rg(sz)], [A[i,i] for i in rg(sz)])
    return [Monomial2Tuple(mnm, As.list(), sz) for mnm in F.operands()]

def TreeFunctionList(n):
    """
    Goes through all the functions and determines which ones
    are associated with trees. One of think of thes trees as
    rooted at 0.

    EXAMPLES:
    ::
        sage: TreeFunctionList(4)
        [[0, 0, 0],
        [2, 0, 0],
        [3, 0, 0],
        [0, 1, 0],
        [3, 1, 0],
        [0, 3, 0],
        [2, 3, 0],
        [3, 3, 0],
        [0, 0, 1],
        [2, 0, 1],
        [0, 1, 1],
        [0, 3, 1],
        [0, 0, 2],
        [2, 0, 2],
        [3, 0, 2],
        [0, 1, 2]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the lists
    l=[n for i in range(n-1)]; Lf=[]
    # Initialization of the identity matrix
    Id=identity_matrix(n)
    # Main loop performing the collecting the functions.
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        f=entry
        # Appending the function to the list
        if is_Tree(sum([Id[:,j]*Id[f[j-1],:]+Id[:,f[j-1]]*Id[j,:] for j in range(1,n)])):
            Lf.append(f)
    return Lf

def Tuple2DiGraphII(T,sz):
    """
    The method returns a directed graph object associated with 
    with the tuple list description of the directed graph

    EXAMPLES:

    ::

        sage: Tuple2DiGraph([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)],5).degree_sequence()
        [4, 2, 2, 1, 1]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the identity matrix
    Id=identity_matrix(sz)
    return DiGraph(sum([Id[:,t[0]]*Id[t[1],:] for t in T]))

def mPer(A,tq):
    """
    Computes symbolically the partial modified permanent by summing only over
    representatives of the coeset.


    EXAMPLES:

    ::

        sage: sz=4; tp=[(0, 0), (1, 2), (2, 0), (3, 0)]
        sage: A=HM(sz,sz,'a').elementwise_product(HM(sz,sz,[x^(sz^abs(j-i)) for j in rg(sz) for i in rg(sz)]))
        sage: mPer(A,tp).coefficient(x^((sz^sz-1)/(sz-1)))
        a00*a12*a20*a30 + a02*a12*a22*a30 + a03*a11*a21*a31 + a03*a13*a21*a33
        sage: sz=4; mPer(HM(sz,sz,'a'), [(i,0) for i in rg(sz)])
        a00*a10*a20*a30 + a01*a11*a21*a31 + a02*a12*a22*a32 + a03*a13*a23*a33
       

    AUTHORS:
    - Edinah K. Gnang
    """
    sz=A.n(0)
    tp=[(1+tq[i][0], 1+tq[i][1]) for i in rg(len(tq))]
    # Initializing the permutations
    P = Permutations(sz); S=SymmetricGroup(sz)
    # Initializing the graph
    grph=Tuple2DiGraph(tp, sz+1)
    # Initializing the automorphism group
    AutGrp=grph.automorphism_group()
    # Initializing representatives of Left coset as strings
    Lcstr=[CstL[0].cycle_string() for CstL in S.cosets(AutGrp)]
    # Loop enumerating the number of graceful labelings
    fctr=0
    # Initialization of the function
    f=0
    for p in P:
        if p.cycle_string() in Lcstr:
            # Initializing the inverse
            pinv=p.inverse()
            # fixing the permutations index
            q =[p[i]-1 for i in rg(sz)]
            qi=[pinv[i]-1 for i in rg(sz)]
            f=f+prod([A[j, q[tp[qi[j]][1]-1]] for j in range(sz)])
    return f

def Monomial2Tuple(mnm, VrbL, sz):
    """
    Outputs the tuple edge list description of the input monomial mnm.
    The variables are obtain by listing the entries of a hypermatrix.
    directed graph case.


    EXAMPLES:
 
    ::

        sage: sz=4; A=HM(sz,sz,'a'); Tp=Monomial2Tuple(prod(A[0,i] for i in rg(sz)), A.list(), sz)
        sage: Tp
        [(0, 0), (0, 1), (0, 2), (0, 3)]
        


    AUTHORS:
    - Edinah K. Gnang
    - To Do: 
    """
    # Getting rid of the coefficient
    f=mnm/(mnm.subs([v==1 for v in VrbL]))
    # Initialization of the dictionary
    EdgeDct=dict([(VrbL[j*sz+i],(i,j)) for j in rg(sz) for i in rg(sz)])
    if f in VrbL:
        return [EdgeDct[f]] 
    else:
        Tp=[EdgeDct[g] for g in f.operands()]; Tp.sort()
    return Tp

def RootedTupleClassTreeFunctionList(tp):
    """
    Computes determines the set all functional
    trees with the same skeleton unalbeled tree.


    EXAMPLES:

    ::

        sage: sz=3; RootedTupleClassTreeFunctionList([(0, 0)]+[(i, i-1) for i in rg(1,sz)])
        [[(0, 0), (1, 0), (2, 0)],
         [(0, 1), (1, 1), (2, 0)],
         [(0, 0), (1, 2), (2, 0)],
         [(0, 0), (1, 0), (2, 1)],
         [(0, 1), (1, 1), (2, 1)],
         [(0, 2), (1, 1), (2, 1)],
         [(0, 2), (1, 0), (2, 2)],
         [(0, 1), (1, 2), (2, 2)],
         [(0, 2), (1, 2), (2, 2)]] 
       
 

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter
    sz=len(tp)
    # Initialization of a symbolic matrix
    A=HM(sz,sz,'a')
    # Computing the sum over modified permanents
    F=sum(mPer(A,T) for T in [switch_sink(tp,i) for i in rg(sz)])
    # converting term to tuple description of functional directed trees
    return [Monomial2Tuple(mnm, A.list(), sz) for mnm in F.operands()]

def BasicLagrangeInterpolation(L, x):
    """
    Implements the basic lagrange interpolation.
    The functions take as input a list of tuples
    and outputs a polynomial in the variable x


    EXAMPLES:
    ::


        sage: x=var('x'); BasicLagrangeInterpolation([(0,0), (1,1), (2,2), (3,3)], x)
        1/2*(x - 1)*(x - 2)*x - (x - 1)*(x - 3)*x + 1/2*(x - 2)*(x - 3)*x


    AUTHORS:
    - Edinah K. Gnang
    """
    # L is a list of tuples
    # Initialized the lenght of the list
    n = len(L)
    # Initialization of the function
    f = 0
    # Code for building the parts 
    for idx in range(len(L)):
        fk = 1
        for j in [i for i in range(len(L)) if i != idx]:
            fk = fk*((x-L[j][0])/(L[idx][0]-L[j][0]))
        f = f + L[idx][1]*fk
    return f

def composition(f, x, k, sz):
    """
    This function performs the composition of the function f in the 
    variable x with itself k times. The implementation assumes that
    f is the functions defined by a functional directed graph on sz
    vertices in order to ensure that the function has degree at most
    sz-1 


    EXAMPLES:
    ::


        sage: f=1/2*(x - 1)*(x - 2)*x - (x - 1)*(x - 3)*x + 1/2*(x - 2)*(x - 3)*x
        sage: composition(f, x, 1, 4)
        1/2*(x - 1)*(x - 2)*x - (x - 1)*(x - 3)*x + 1/2*(x - 2)*(x - 3)*x
        sage: sz=4; tq=[(0, 1), (1, 2), (2, 3), (3, 3)]
        sage: f=BasicLagrangeInterpolation(tq,var('z'))
        sage: tq == [(i, f.subs(z==i)) for i in rg(sz)]
        True
        sage: L=[composition(f, var('z'), i, sz) for i in rg(sz)]; L
        [z,
         -1/6*(z - 1)*(z - 2)*(z - 3) + 1/2*(z - 1)*(z - 2)*z - 3/2*(z - 1)*(z - 3)*z + (z - 2)*(z - 3)*z,
         -1/3*(z - 1)*(z - 2)*(z - 3) + 1/2*(z - 1)*(z - 2)*z - 3/2*(z - 1)*(z - 3)*z + 3/2*(z - 2)*(z - 3)*z,
         -1/2*(z - 1)*(z - 2)*(z - 3) + 1/2*(z - 1)*(z - 2)*z - 3/2*(z - 1)*(z - 3)*z + 3/2*(z - 2)*(z - 3)*z]
        sage: [[(i, g.subs(var('z')==i)) for i in rg(sz)] for g in [composition(f, var('z'), i, sz) for i in rg(sz)]]
        [[(0, 0), (1, 1), (2, 2), (3, 3)],
         [(0, 1), (1, 2), (2, 3), (3, 3)],
         [(0, 2), (1, 3), (2, 3), (3, 3)],
         [(0, 3), (1, 3), (2, 3), (3, 3)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    if k==0:
        return x
    elif k==1:
        return f
    elif k > 1:
        # Initialization of the function
        g=f
        for t in rg(k-1):
            # Initialization of the function
            TmpF=0
            # Code for building the parts 
            for idx in range(sz):
                fk=1
                for j in [i for i in rg(sz) if i != idx]:
                    fk=fk*(x-j)/(idx-j)
                TmpF=TmpF + g.subs(x==f).subs(x==idx)*fk
            g=TmpF
        return g

def compose_with(f, g, sz, vrbl):
    """
    This function performs the composition of the function f in the 
    variable x with itself k times. The implementation assumes that
    f is the functions defined by a functional directed graph on sz
    vertices in order to ensure that the function has degree at most
    sz-1 


    EXAMPLES:
    ::


        sage: f0=BasicLagrangeInterpolation([(0,1), (1,2), (2,0)],x); f1=BasicLagrangeInterpolation([(0,2), (1,0), (2,1)],x); f2=x
        sage: f3=0; f4=BasicLagrangeInterpolation([(0,0), (1,2), (2,1)],x)
        sage: expand(compose_with(f0, f1, 3, x))
        x
        sage: compose_with(f4, f4, 3, x) - composition(f4, x, 2, 3)
        0


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the function
    TmpF=0
    # Code for building the parts 
    for idx in range(sz):
        fk=1
        for j in [i for i in rg(sz) if i != idx]:
            fk=fk*(vrbl-j)/(idx-j)
        TmpF=TmpF + f.subs(vrbl==g).subs(vrbl==idx)*fk
    g=TmpF
    return g

def compose_tuple(T0, T1):
    """
    This function performs the composition of the input tuple functions.


    EXAMPLES:
    ::


        sage: T0=[(0,1), (1,2), (2,0)]; T1=[(0,2), (1,0), (2,1)]
        sage: compose_tuple(T0, T1)
        [(0, 0), (1, 1), (2, 2)]


    AUTHORS:
    - Edinah K. Gnang
    """
    if len(T0) != len(T1):
        raise ValueError("Expected tuples of the same size.")
    else:
        # Initialization of the polynomial encoding
        f0=BasicLagrangeInterpolation(T0, x); f1=BasicLagrangeInterpolation(T1, x)
        # Computing the polynomial enconding of the composition
        g=compose_with(f0, f1, len(T0), x)
        return [(i,g.subs(x==i)) for i in rg(len(T0))]

def invert_tuple(T):
    """
    This function computes the composition inverse of 
    the input tuple.


    EXAMPLES:
    ::


        sage: sz=3; T0=[(0, 1), (1, 2), (2, 0)]; T1=invert_tuple(T0)
        sage: compose_tuple(T0, T1)
        [(0, 0), (1, 1), (2, 2)]
        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of a tuple description for f
    TpL=[(T[i][1], i) for i in rg(len(T))]; TpL.sort()
    if [t[0] for t in TpL] == rg(len(TpL)):
        return TpL
    else:
        raise ValueError("The input function must be invertible")

def compose_tuple_pow(T, k):
    """
    This function performs the composition of the input tuple.
    The implementation assumes that T is a functional directed graph.


    EXAMPLES:
    ::


        sage: T=[(0, 1), (1, 2), (2, 3), (3, 3)]; L=[compose_tuple_pow(T, i) for i in rg(len(T))]; L
        [[(0, 0), (1, 1), (2, 2), (3, 3)],
         [(0, 1), (1, 2), (2, 3), (3, 3)],
         [(0, 2), (1, 3), (2, 3), (3, 3)],
         [(0, 3), (1, 3), (2, 3), (3, 3)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    if k==0:
        return [(i, i) for i in rg(len(T))]
    elif k==1:
        return T
    elif k > 1:
        # Initialization of the function
        f=BasicLagrangeInterpolation(T, x)
        g=f
        for t in rg(k-1):
            # Initialization of the function
            TmpF=0
            # Code for building the parts 
            for idx in range(len(T)):
                fk=1
                for j in [i for i in rg(len(T)) if i != idx]:
                    fk=fk*(x-j)/(idx-j)
                TmpF=TmpF + g.subs(x==f).subs(x==idx)*fk
            g=TmpF
        return [(i, g.subs(x==i)) for i in rg(len(T))]

def composition_inverse(f, sz, vrbl):
    """
    This function computes the composition inverse of f.
    The implementation test whether f is invertible.


    EXAMPLES:
    ::


        sage: sz=3; f0=BasicLagrangeInterpolation([(0,1), (1,2), (2,0)], x)
        sage: f1=composition_inverse(f0, sz, x)
        sage: expand(compose_with(f0, f1, sz, x))
        x


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of a tuple description for f
    TpL=[(f.subs(vrbl==i), i) for i in rg(sz)]; TpL.sort()
    if [t[0] for t in TpL] == rg(sz):
        g = BasicLagrangeInterpolation(TpL, vrbl) 
    else:
        raise ValueError("The input function must be invertible")
    return g

def find_sink(tp):
    """
    Finds the sink (also called root of a functional tree).


    EXAMPLES:

    ::


        sage: find_sink([(0, 0), (1, 3), (2, 1), (3, 0), (4, 0)])
        0
        

    AUTHORS:
    - Edinah K. Gnang
    """
    T=[(tp[i][0],tp[i][1]) for i in tpl_image_set(tp)]
    for i in rg(len(tp)):
        T=[(tp[i][0],tp[i][1]) for i in tpl_image_set(T)]
    return T[0][0]

def FindTreeTupleComponentsII(T):
    """
    Returns a tuple list each of which corresponds to a pair
    made up of a vertex and the associated sink vertex.
    This implementation assume the input is a tree.


    EXAMPLES:

    ::


        sage: FindTreeTupleComponentsII([(0, 0), (1, 2), (2, 4), (3, 7), (4, 4), (5, 0), (6, 0), (7, 0)])
        [(0, 0), (1, 4), (2, 4), (3, 0), (4, 4), (5, 0), (6, 0), (7, 0)]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Classifying the vertices
    C=[(find_sink(T),find_sink(T))]
    Rg=rg(len(T)); Rg.remove(find_sink(T))
    for j in Rg:
        i=j; cnt=0
        while T[i][1] != T[i][0] and cnt <= len(T)+1:
            i=T[i][1]; cnt=cnt+1
        if T[i][1] == T[i][0]:
            C.append((j,T[i][0]))
        else:
            raise ValueError("Expected a tree")
    C.sort()
    return C

def switch_sink(Tn, alpha):
    """
    Switches the sink from of the input functional tree to alpha.


    EXAMPLES:

    ::


        sage: switch_sink([(0, 0), (1, 3), (2, 1), (3, 0), (4, 0)], 1)
        [(0, 3), (1, 1), (2, 1), (3, 1), (4, 0)]



    AUTHORS:
    - Edinah K. Gnang
    """
    a = find_sink(Tn)
    if a == alpha:
        return Tn
    else:
        # The new sink will be alpha instead of a
        Tp=[]
        # Obtaining the spine made up of the path from alpha to a
        Snk=[]; i=alpha; Tp.append((alpha, alpha))
        while Tn[i][0] != a and Tn[i][1] != a:
            Tp.append((Tn[i][1], Tn[i][0])); Snk.append(Tn[i][0])
            i=Tn[i][1]
            #print 'Tp=',Tp
        Tp.append((Tn[i][1], Tn[i][0])); Snk.append(Tn[i][0])
        Snk.append(Tn[i][1])
        #print 'Snk=',Snk
        # Partioning the vertices of Tn by components
        C=FindTreeTupleComponentsII(Tn)
        #print 'C=',C
        # Initialization of the list of vertices in the partition a
        Sa=Set([C[i][0] for i in range(len(Tn)) if C[i][1]==a])
        #print 'Sa=',Sa
        # Correcting the orrientation of the remaining edges
        while Set(Snk) != Sa:
            for i in Sa:
                if (Tn[i][0] in Snk) and (Tn[i][1] not in Snk):
                    Tp.append((Tn[i][1], Tn[i][0])); Snk.append(Tn[i][1])
                elif (Tn[i][0] not in Snk) and (Tn[i][1] in Snk):
                    Tp.append((Tn[i][0], Tn[i][1])); Snk.append(Tn[i][0])
        # Sorting the list
        Tp.sort()
        return Tp

def tpl_pre_image_set(tp, i):
    """
    returns the list of vertex pre-images of the input vertex.
    The input graphs is assumes to have no isolated vertices.

    EXAMPLES:

    ::

        sage: tpl_pre_image_set([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)], 0)
        [0, 3, 4]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    return Set([t[0] for t in tp if t[1]==i]).list()

def tpl_image_set(tp):
    """
    returns the list of vertex image of the domain of the input
    function associated with the list of tuples.
    The input graphs is assumes to have no isolated vertices.


    EXAMPLES:

    ::

        sage: tpl_image_set([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)])
        [0, 2, 4]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing and sorting the list
    L=Set([t[1] for t in tp]).list(); L.sort()
    return L

def tpl_pre_image_set_function(tp):
    """
    returns the list of vertex pre-images of every vertex
    the vertices being sorted in increasing order.
    The input graphs is assumes to have no isolated vertices.


    EXAMPLES:

    ::

        sage: tpl_pre_image_set_function([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)])
        [[0, 3, 4], [], [1], [], [2]]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    return [tpl_pre_image_set(tp, i) for i in rg(len(tp))]

def tpl_leaf_set(tp):
    """
    returns the list of leaf vertex set of the functional directed graph
    function associated with the list of tuples.


    EXAMPLES:

    ::

        sage: tpl_leaf_set([(0, 0), (1, 2), (2, 4), (3, 0), (4, 0)])
        [1, 3]
        

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initializing and sorting the list
    L0=Set([t[1] for t in tp]).list(); L0.sort()
    L1=[tp[i][0] for i in rg(len(tp)) if not tp[i][0] in L0]; L1.sort()
    return L1

def gcomp(tp, tq):
    """
    Performs the absolute difference composition of functions
    the inputs to the function are two tuple lists of the same
    size. The current version does not check that the list are
    in fact of the same size.


    EXAMPLES:

    ::

        sage: tp=[(0,1), (1,2), (2,3), (3,0)]; tq=[(0,0), (1,0), (2,0), (3,0)]
        sage: gcomp(tp, tq)
        [(0, 1), (1, 2), (2, 3), (3, 0)]
       

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter
    sz=len(tq); tr=[(i,i) for i in rg(sz)]
    # Updating the composition
    for i in rg(sz):
        tr[abs(tq[i][1]-tq[i][0])] = (abs(tq[i][0]-tq[i][1]), tp[abs(tq[i][0]-tq[i][1])][1])
    return tr

def functionaldigraphincidenceHM(T):
    """
    Returns the transpose of the conventional incidence matrix associated with the input
    tree specified by a list of tuples associated with the function. 


    EXAMPLES:

    ::

        sage: functionaldigraphincidenceHM([(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)])
        [[0, 0, 0, 0, 0], [1, -1, 0, 0, 0], [1, 0, -1, 0, 0], [1, 0, 0, -1, 0], [1, 0, 0, 0, -1]]


    AUTHORS:
    - Edinah K. Gnang
    """
    sz = len(T)
    # Initialization of the list of variables
    X = var_list('x',sz); Y = var_list('y',sz)
    # Initialization of the list of constraints
    Eq=[X[T[i][1]]-X[T[i][0]] == Y[i] for i in rg(sz)]
    return ConstraintFormatorHM(Eq, X)[0]

def functionaldigraphunsignedincidenceHM(T):
    """
    Returns the transpose of the conventional incidence matrix associated with the input
    tree specified by a list of tuples associated with the function. 


    EXAMPLES:

    ::

        sage: functionaldigraphunsignedincidenceHM([(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)])
        [[2, 0, 0, 0, 0], [1, 1, 0, 0, 0], [1, 0, 1, 0, 0], [1, 0, 0, 1, 0], [1, 0, 0, 0, 1]]


    AUTHORS:
    - Edinah K. Gnang
    """
    sz = len(T)
    # Initialization of the list of variables
    X = var_list('x',sz); Y = var_list('y',sz)
    # Initialization of the list of constraints
    Eq=[X[T[i][1]]+X[T[i][0]] == Y[i] for i in rg(sz)]
    return ConstraintFormatorHM(Eq, X)[0]

def EuclidsGCD(a, b):
    """

    This function implements Euclid's GCD algorithm.
    The function checks that the inputs are non-zero
    integers and returns as output the matrix which
    describes all the iterations of the algorithm.


    EXAMPLES:
    ::


        sage: G=EuclidsGCD(89, 55); G.printHM()
        [:, :]=
        [89 55  1 34]
        [55 34  1 21]
        [34 21  1 13]
        [21 13  1  8]
        [13  8  1  5]
        [ 8  5  1  3]
        [ 5  3  1  2]
        [ 3  2  1  1]
        [ 2  1  2  0]


    AUTHORS:
    - Edinah K. Gnang
    """
    if type(a)==type(Integer(1)) and type(b)==type(Integer(1)) and a*b !=0:
        # Initialization of the matrix Data Structure.
        G = HM(1, 4, 'zero')
        # Initialization of the initial conditions
        G[0, 0] = a; G[0, 1] = b
        G[0, 3] = Integer(mod(G[0, 0], G[0, 1]))
        G[0, 2] = (G[0, 0] - G[0, 3])/G[0, 1]
        # Initialization of the index
        indx = 0
        while G[indx, 3] != 0:
            # Updating the size of G
            G=G.zero_pad([G.n(0)+1, G.n(1)])
            # Incrementing the index
            indx=indx+1
            G[indx, 0] = G[indx-1, 1]; G[indx, 1] = G[indx-1, 3]
            G[indx, 3] = Integer(mod(G[indx, 0], G[indx, 1]))
            G[indx, 2] = (G[indx, 0] - G[indx, 3])/G[indx, 1]
        return G
    else:
        raise ValueError("Expected non zero integer inputs.")

def Modulo(f, VrbL, Rlts):
    """
    Outputs the quotient and the remainder of the simultaneous
    Euclidean division. The algorithm takes as input
    a multivariate polynomial f in the variables in the
    input list VrbL and a list of monic univariate polynomial
    which correspond to the relations we are moding by.


    EXAMPLES:

    ::
        
        sage: VrbL = var_list('x', 2); f = expand((VrbL[0]+VrbL[1])^10); Rlts = [VrbL[0]^2-5, VrbL[1]^3-7]; Modulo(f, VrbL, Rlts)[1]
        44100*x0*x1^2 + 35650*x0*x1 + 39150*x1^2 + 108430*x0 + 184093*x1 + 260375


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the degree matrix
    Dm = degree_matrix(Rlts, VrbL)
    for i in rg(min(Dm.n(0), Dm.n(1))):
        for j in rg(min(Dm.n(0), Dm.n(1))):
            Dm[i,j]=0
    if Dm.is_zero():
        # Initialization of the quotient
        # and the initial remainder.
        q = 0; r = f
        for v in rg(len(VrbL)):
            for d in rg(f.degree(VrbL[v])-Rlts[v].degree(VrbL[v]), -1, -1):
                q = q + (VrbL[v]^d)*r.coefficient(VrbL[v]^(d+Rlts[v].degree(VrbL[v])))*Rlts[v]
                r = expand(fast_reduce(r, [VrbL[v]^(d+Rlts[v].degree(VrbL[v]))], [VrbL[v]^(d+Rlts[v].degree(VrbL[v]))-expand(Rlts[v]*VrbL[v]^d)])) 
        return [q, r]
    else:
        raise ValueError("Expected univariate algebraic relations.")

def ModuloII(f, VrbL, Rlts):
    """
    Outputs the remainder of the simultaneous
    Euclidean division. The algorithm takes as input
    a multivariate polynomial f in the variables in the
    input list VrbL and a list of monic univariate polynomial
    which correspond to the relations we are moding by.


    EXAMPLES:

    ::
        
        sage: VrbL = var_list('x', 2); f = expand((VrbL[0]+VrbL[1])^10); Rlts = [VrbL[0]^2-5, VrbL[1]^3-7]
        sage: ModuloII(f, VrbL, Rlts)
        44100*x0*x1^2 + 35650*x0*x1 + 39150*x1^2 + 108430*x0 + 184093*x1 + 260375
        

    AUTHORS:
    - Edinah K. Gnang
    - To Do: Implment faster version
    """
    # Initialization of the degree matrix
    Dm = degree_matrix(Rlts, VrbL)
    # Updating only the diagonal entries to zero
    for i in rg(min(Dm.n(0), Dm.n(1))):
        for j in rg(min(Dm.n(0), Dm.n(1))):
            Dm[i,j]=0
    if Dm.is_zero():
        # Initialization of the initial remainder.
        r = f
        for v in rg(len(VrbL)):
            for d in rg(f.degree(VrbL[v])-Rlts[v].degree(VrbL[v]), -1, -1):
                r = expand(fast_reduce(r, [VrbL[v]^(d+Rlts[v].degree(VrbL[v]))], [VrbL[v]^(d+Rlts[v].degree(VrbL[v]))-expand(Rlts[v]*VrbL[v]^d)])) 
        return r
    else:
        raise ValueError("Expected univariate algebraic relations.")

def remainder_via_lagrange_interpolation(f, Lr, X):
    """
    Returns the canonical representative of the residue class
    f modulo Ld. The input f is an arbitrary multivariate polynomial
    in the variables stored in the list input X. Lr is is the list of distinct
    roots of the monic univariate polynomial associated with each variable.
    This implementation does not require the input polunomial to be in 
    factored form


    EXAMPLES:
    ::

        sage: sz=3; X=var_list('x',sz); f=expand(var('x0')^2*var('x1')+var('x2')^4); Lr=rg(sz)
        sage: expand(remainder_via_lagrange_interpolation(f, Lr, X))
        x0^2*x1 + 7*x2^2 - 6*x2


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization the number of variables
    sz=len(X)
    return sum(f.subs([X[i] == Lr[t[i][1]] for i in rg(sz)])*prod(prod((X[k]-jk)/(t[k][1]-jk) for jk in Lr if jk !=t[k][1]) for k in rg(sz) ) for t in TupleFunctionList(sz))

def GeneralHypermatrixFlatten(A, Rg, ord):
    """
    Outputs a lower order flattened hypermatrix of the higher order input. 
    The first input is the hypermatrix to be flattened. The second input 
    corresponds to the indices that will remain unchanged after flattening. 
    The third input is the order of the desired output hypermatrix. When 
    flatten an higher order input hypermatrix to order 1, the result is a list.
    
    In the example of GeneralHypermatrixFlatten(A, [0, 1], 2) for flattening 
    the order 3 hypermatrix A into a matrix, the input [0, 1] tells us that the
    column and the row indices will be left unchanged, and we are stacking depth 
    slices along the row direction. 
    
    The flattening of hypermatrices with ord>2 is obtained recursively.

    EXAMPLES:

    ::
        sage: sz=3; A=HM(sz,sz,sz,'a')
        sage: GeneralHypermatrixFlatten(A, [0, 1], 2).dimensions()
        [3, 9]
        sage: GeneralHypermatrixFlatten(A, [0, 1], 2).printHM()
        [:, :]=
        [a000 a010 a020 a001 a011 a021 a002 a012 a022]
        [a100 a110 a120 a101 a111 a121 a102 a112 a122]
        [a200 a210 a220 a201 a211 a221 a202 a212 a222]
        sage: sz=2; B=HM(sz,sz,sz,'b')
        sage: GeneralHypermatrixFlatten(B, [2], 1)
        [b000, b001, b100, b101, b010, b011, b110, b111]
 
 

    AUTHORS:
    - Edinah K. Gnang and Fan Tian
    - To Do: 
    """
    if A.order() <= ord:
        return A
    elif ord == 1:
        return HM(prod(A.dimensions()),A.transpose(Rg[0]).list())
    elif ord == 2:
        return HM([GeneralHypermatrixFlatten(A.slice([t],Rg[0]),[Rg[1]],1).list() for t in rg(A.n(Rg[0]))])
    elif ord == 3:
        return HM([GeneralHypermatrixFlatten(A.slice([t],Rg[0]),[Rg[1],Rg[2]],2).listHM() for t in rg(A.n(Rg[0]))])
    elif ord == 4:
        return HM([GeneralHypermatrixFlatten(A.slice([t],Rg[0]),[Rg[1],Rg[2],Rg[3]],3).listHM() for t in rg(A.n(Rg[0]))])
    elif ord == 5:
        return HM([GeneralHypermatrixFlatten(A.slice([t],Rg[0]),[Rg[1],Rg[2],Rg[3],Rg[4]],4).listHM() for t in rg(A.n(Rg[0]))])

def KroneckerVectorOuterProduct(*args):
    """
    Outputs the hypermatrix resulting from the outer product
    of properly embeded and orriented vectors. The input
    are column vector of any order. This operation is motivated
    by the classical vector outer-product.


    EXAMPLES:
 
    ::

        sage: U=HM(2,1,var_list('u',3)); V=HM(3,1,var_list('v',3)); W=HM(4,1,var_list('w',4))
        sage: KroneckerVectorOuterProduct(U, V, W).printHM()
        [:, :, 0]=
        [u0*v0*w0 u0*v1*w0 u0*v2*w0]
        [u1*v0*w0 u1*v1*w0 u1*v2*w0]
        <BLANKLINE>
        [:, :, 1]=
        [u0*v0*w1 u0*v1*w1 u0*v2*w1]
        [u1*v0*w1 u1*v1*w1 u1*v2*w1]
        <BLANKLINE>
        [:, :, 2]=
        [u0*v0*w2 u0*v1*w2 u0*v2*w2]
        [u1*v0*w2 u1*v1*w2 u1*v2*w2]
        <BLANKLINE>
        [:, :, 3]=
        [u0*v0*w3 u0*v1*w3 u0*v2*w3]
        [u1*v0*w3 u1*v1*w3 u1*v2*w3]
        sage: U=HM(2,1,var_list('u',3)); V=HM(3,1,var_list('v',3))
        sage: KroneckerVectorOuterProduct(U, V).printHM()
        [:, :]=
        [u0*v0 u0*v1 u0*v2]
        [u1*v0 u1*v1 u1*v2]


    AUTHORS: Edinah Gnang and Fan Tian
    - To Do: 
    """
    if len(args) == 2:
        return Prod(HM((args[0]).n(0),1,(args[0]).list()), HM((args[1]).n(0),1,(args[1]).list()).transpose())
    else:
        # Initialization of the canonical embeding into vectors of the appropriate order.
        #Lv = [apply(HM,[(args[i]).n(0)]+[1 for j in rg(len(args)-1)]+[(args[i]).list()]) for i in rg(len(args))]
        Lv = [HM(*([(args[i]).n(0)]+[1 for j in rg(len(args)-1)]+[(args[i]).list()])) for i in rg(len(args))]
        # Initialization of the list hypermatrix slices of the appropriate orders.
        #Lh=[Lv[j].tensor_product(apply(HM,[1,1]+[(args[i]).n(0) for i in rg(2,len(args))]+['one'])) for j in rg(len(args))]
        Lh=[Lv[j].tensor_product(HM(*([1,1]+[(args[i]).n(0) for i in rg(2,len(args))]+['one']))) for j in rg(len(args))]
        # returning the BM product of the appropriately orriented slices
        return Prod(*[Lh[i].transpose(len(args)-i) for i in rg(len(args))])

def HypermatrixList(bnd, DmsL):
    """
    Outputs a list of all hypermatrix with entries
    in {0,1,2,...,bnd-1} of size specified by the list DmsL
    This fucntion is motivated by the study of hypermatrices
    with elements from a finite field.


    EXAMPLES:
 
    ::

        sage: L=HypermatrixList(2,[2,2,2])
        sage: L[0].printHM()
        [:, :, 0]=
        [0 0]
        [0 0]
        <BLANKLINE>
        [:, :, 1]=
        [0 0]
        [0 0]


    AUTHORS: Edinah Gnang 
    - To Do: 
    """
    # Initialization of the list specifying the dimensions of the output
    l=[bnd for i in rg(prod(DmsL))]
    # Initialization of the list of Hypermatrices
    Lh=[]
    # Main loop performing the transposition of the entries
    for i in range(prod(l)):
        # Turning the index i into an hypermatrix array location using the decimal encoding trick
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        #Lh.append(apply(HM,DmsL+[entry]))
        Lh.append(HM(*(DmsL+[entry])))
    return Lh

def HypermatrixSn(L):
    """
    Generates the permutation pair of third order hypermatrices.
    hypermatrix deduced from sigma. Note that as a result of 
    the  non associativity, permutations must be performed as
    one transposition at a time.

    EXAMPLES:

    ::

        sage: P0,P1 = HypermatrixSn([1,2,0])
        sage: Prod(P0,HM(3,3,3,'a'),P1).printHM()
        [:, :, 0]=
        [a001 a011 a021]
        [a101 a111 a121]
        [a201 a211 a221]
        <BLANKLINE>
        [:, :, 1]=
        [a002 a012 a022]
        [a102 a112 a122]
        [a202 a212 a222]
        <BLANKLINE>
        [:, :, 2]=
        [a000 a010 a020]
        [a100 a110 a120]
        [a200 a210 a220]


    AUTHORS:
    - Edinah K. Gnang, Ori Parzanchevski, Fan Tian 
    """
    sz=len(L)
    # Test for dimension match
    if sz > 0:
        # Initialization of the first input (P0)
        P0=HM([HM(2,L,'perm').transpose().listHM() for i in rg(sz)])
        # Initialization of the second input (P1)
        P1=HM([HM(2,L,'perm').listHM() for i in rg(sz)]).transpose(2)
        return [P0, P1]
    else :
        raise ValueError("Input list must me non empty ")

def layer_rotation(M, v):
    """
    Performs a layered rotation which amounts to viewing the entries
    as forming a spanning union of cycles whose cycles are made up
    of layers of entries. The input matrix M must be square, the second
    input v must be a character will will be used as the base symbol 
    variable to be used to fill up the template matrix. It is 
    important that indexings of v do not appear in M. In this current
    implementation it will cause bugs.


    EXAMPLES:

    ::


        sage: [U, V]=layer_rotation(HM(4,4,'c'), 'a')
        sage: U.printHM()
        [:, :]=
        [c10 c00 c01 c02]
        [c20 c21 c11 c03]
        [c30 c22 c12 c13]
        [c31 c32 c33 c23]
        sage: V.printHM()
        [:, :]=
        [c01 c02 c03 c13]
        [c00 c12 c22 c23]
        [c10 c11 c21 c33]
        [c20 c30 c31 c32]        


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the size parameter
    sz=min(M.dimensions())
    # Iinitialization of the symbolic template matrix
    # This codes collects all the layers into a list of lists
    A=HM(sz, sz, v); L=[[] for i in rg(ceil(sz/2))]
    for l in rg(ceil(sz/2)): 
        TmpA = A.slice(rg(l,sz-l),0).slice(rg(l,sz-l),1)
        for i in rg(3):
            hA = TmpA.index_rotation(-2*i*pi/4)
            tmpL=hA.listHM()[0]
            for j in rg(hA.n(1)):
                if not tmpL[j] in L[l]:
                    L[l].append(tmpL[j])
        tmpL=TmpA.transpose().listHM()[0]
        tmpL.reverse()
        for j in rg(hA.n(1)):
            if not tmpL[j] in L[l]:
                L[l].append(tmpL[j])
    # Performing the cyclic permutation
    Ln=[l[1:]+[l[0]] for l in L]
    # Concatenating the lists
    L0=L[0];L1=Ln[0]
    for i in rg(1,len(L)):
        L0=L0+L[i]; L1=L1+Ln[i]
    # Performing the substition for the forward mapping
    B=A.subs([L1[i]==L0[i] for i in rg(len(L0))])
    # Performing the substition for the backward mapping
    C=A.subs([L0[i]==L1[i] for i in rg(len(L0))])
    return [B.subs([A[i,j]==M[i,j] for j in rg(sz) for i in rg(sz)]), C.subs([A[i,j]==M[i,j] for j in rg(sz) for i in rg(sz)])]
 
def ThirdOrderIdentityPair(DimList):
    """
    Generates the left right identity pair for a midel third
    order hypermatrix whose dimensions corresponds to the input
    

    EXAMPLES:

    ::

        sage: A=HM(2,3,4,'a')
        sage: I0,I1 = ThirdOrderIdentityPair(A.dimensions())
        sage: Prod(I0, A, I1).printHM()
        [:, :, 0]=
        [a000 a010 a020]
        [a100 a110 a120]
        <BLANKLINE>
        [:, :, 1]=
        [a001 a011 a021]
        [a101 a111 a121]
        <BLANKLINE>
        [:, :, 2]=
        [a002 a012 a022]
        [a102 a112 a122]
        <BLANKLINE>
        [:, :, 3]=
        [a003 a013 a023]
        [a103 a113 a123]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the first of the two identity hypermatrix pairs
    I0=HM(DimList[0], DimList[2], DimList[2], 'zero')
    for i in rg(I0.n(0)):
        for t in rg(I0.n(1)):
            I0[i,t,t]=1
    # Initialization of the second of the two identity hypermatrix pairs
    I1=HM(DimList[2], DimList[1], DimList[2], 'zero')
    for j in rg(I1.n(1)):
        for t in rg(I1.n(0)):
            I1[t,j,t]=1
    return [I0, I1]

def GeneralHypermatrixSupport(A):
    """
    Procedure for determining the support
    of non-zero entries

    EXAMPLES:

    ::

        sage: GeneralHypermatrixSupport(HM(2,2,2,'a')).printHM()
        [:, :, 0]=
        [1 1]
        [1 1]
        <BLANKLINE>
        [:, :, 1]=
        [1 1]
        [1 1]


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
        entry = [Integer(mod(i,l[0]))]
        sm = Integer(mod(i,l[0]))
        for k in range(len(l)-1):
            entry.append(Integer(mod(Integer((i-sm)/prod(l[0:k+1])),l[k+1])))
            sm = sm+prod(l[0:k+1])*entry[len(entry)-1]
        if A[tuple(entry)].is_zero():
            Rh[tuple(entry)]=Integer(0)
        else:
            Rh[tuple(entry)]=Integer(1)
    return Rh

def EliminationHMI(Ha,i):
    """
    Procedure for performing a row linear combination which cancels the first entry
    of non-zero entries. This implementation can serve as the basis for a parallel
    solver


    EXAMPLES:

    ::

        sage: sz=3; Ha=HM(sz,sz,'a'); Hb=EliminationHMI(Ha,1).canonicalize_radical()
        sage: Hb.printHM()
        [:, :]=
        [                                                                              a00                                                                               a01                                                                               a02]
        [                                                                                0 -1/2*a00*a10*a21*(I*sqrt(3) + 1) - 1/2*(a00*a11*(-I*sqrt(3) + 1) - 2*a01*a10)*a20 -1/2*a00*a10*a22*(I*sqrt(3) + 1) - 1/2*(a00*a12*(-I*sqrt(3) + 1) - 2*a02*a10)*a20]
        [                                                                              a20                                                                               a21                                                                               a22]
        sage: sz=3; Ha=HM(sz,sz,'a'); Hb=EliminationHMI(Ha,2).canonicalize_radical()
        sage: Hb.printHM()
        [:, :]=
        [                                                                             a00                                                                              a01                                                                              a02]
        [                                                                             a10                                                                              a11                                                                              a12]
        [                                                                               0 1/2*a00*a10*a21*(I*sqrt(3) - 1) + 1/2*(a00*a11*(-I*sqrt(3) - 1) + 2*a01*a10)*a20 1/2*a00*a10*a22*(I*sqrt(3) - 1) + 1/2*(a00*a12*(-I*sqrt(3) - 1) + 2*a02*a10)*a20] 
        sage: Hb=EliminationHMI(Ha,1); HM(sz,sz,[Hb[i,j].is_zero() for j in rg(sz) for i in rg(sz)]).printHM()
        [:, :]=
        [0 0 0]
        [1 0 0]
        [0 0 0]
        sage: Hb=EliminationHMI(Ha,2); HM(sz,sz,[Hb[i,j].is_zero() for j in rg(sz) for i in rg(sz)]).printHM()
        [:, :]=
        [0 0 0]
        [0 0 0]
        [1 0 0]


    AUTHORS:
    - Edinah K. Gnang
    """
    A=Ha.copy()
    #Sm=sum([exp(I*2*pi*i*k/A.n(0))*prod([A[j,0] for j in rg(A.n(0)) if j!=k])*A.slice([k],'row') for k in rg(A.n(0))])
    Sm=sum([exp(I*2*pi*i*k/A.n(0))*GeneralHypermatrixScaleRight(GeneralHypermatrixScale(A.slice([k],'row'),prod(A[j,0] for j in rg(A.n(0)) if j<k)), prod([A[j,0] for j in rg(A.n(0)) if j>k])) for k in rg(A.n(0))])
    for jndx in rg(A.n(1)):
        A[i,jndx]=Sm[0,jndx]
    return A

def EliminationHMII(Ha):
    """
    Procedure for clearing bellow the first pivot


    EXAMPLES:

    ::

        sage: sz=3; Ha=HM(sz,sz,'a'); Hc=EliminationHMII(Ha).canonicalize_radical()
        sage: Hc.printHM()
        [:, :]=
        [                                                                              a00                                                                               a01                                                                               a02]
        [                                                                                0 -1/2*a00*a10*a21*(I*sqrt(3) + 1) - 1/2*(a00*a11*(-I*sqrt(3) + 1) - 2*a01*a10)*a20 -1/2*a00*a10*a22*(I*sqrt(3) + 1) - 1/2*(a00*a12*(-I*sqrt(3) + 1) - 2*a02*a10)*a20]
        [                                                                                0  1/2*a00*a10*a21*(I*sqrt(3) - 1) + 1/2*(a00*a11*(-I*sqrt(3) - 1) + 2*a01*a10)*a20  1/2*a00*a10*a22*(I*sqrt(3) - 1) + 1/2*(a00*a12*(-I*sqrt(3) - 1) + 2*a02*a10)*a20]
        sage: Deter(Hc).is_zero()
        False
        sage: Ha = HM([[HM(2,2,'a'),HM(2,2,'b'),HM(2,2,'c')],[HM(2,2,'d'),HM(2,2,'e'),HM(2,2,'f')],[HM(2,2,'g'),HM(2,2,'h'),HM(2,2,'i')]])
        sage: Hb=EliminationHMII(Ha); HM(sz,sz,[Hb[i,j].is_zero() for j in rg(sz) for i in rg(sz)]).printHM()
        [:, :]=
        [0 0 0]
        [1 0 0]
        [1 0 0]


    AUTHORS:
    - Edinah K. Gnang
    """
    A=Ha.copy(); B=Ha.copy()
    for i in rg(1,A.n(0)):
        # Zeroing out the first entry of row i
        Tmp=EliminationHMI(A,i)
        for jndx in rg(A.n(1)):
            B[i,jndx]=Tmp[i,jndx]
    return B

def GaussEliminationHM(Ha):
    """
    Procedure for obtaining the Row Echelon Form (REF)


    EXAMPLES:

    ::

        sage: sz=3; Ha=HM(sz,sz,'a'); Hc=GaussEliminationHM(Ha).canonicalize_radical()
        sage: Hc.printHM()
        [:, :]=
        [                                                                                                                                                                                        a00                                                                                                                                                                                         a01                                                                                                                                                                                         a02]
        [                                                                                                                                                                                          0                                                                                                           -1/2*a00*a10*a21*(I*sqrt(3) + 1) - 1/2*(a00*a11*(-I*sqrt(3) + 1) - 2*a01*a10)*a20                                                                                                           -1/2*a00*a10*a22*(I*sqrt(3) + 1) - 1/2*(a00*a12*(-I*sqrt(3) + 1) - 2*a02*a10)*a20]
        [                                                                                                                                                                                          0                                                                                                                                                                                           0 (-I*sqrt(3)*a00*a02*a10*a11 + I*sqrt(3)*a00*a01*a10*a12)*a20^2 + (I*sqrt(3)*a00*a02*a10^2 - I*sqrt(3)*a00^2*a10*a12)*a20*a21 + (-I*sqrt(3)*a00*a01*a10^2 + I*sqrt(3)*a00^2*a10*a11)*a20*a22]
        sage: Ha = HM([[HM(2,2,'a'),HM(2,2,'b'),HM(2,2,'c')],[HM(2,2,'d'),HM(2,2,'e'),HM(2,2,'f')],[HM(2,2,'g'),HM(2,2,'h'),HM(2,2,'i')]])
        sage: Hb=GaussEliminationHM(Ha); HM(sz,sz,[Hb[i,j].is_zero() for j in rg(sz) for i in rg(sz)]).printHM()
        [:, :]=
        [0 0 0]
        [1 0 0]
        [1 1 0]


    AUTHORS:
    - Edinah K. Gnang
    """
    A=Ha.copy()
    for i in rg(A.n(0)-1):
        A = A.fill_with(EliminationHMII(A.slice(rg(i,A.n(0)),'row').slice(rg(i,A.n(1)),'col')),i)
    return A

def GaussJordanEliminationHM(Ha):
    """
    Procedure for obtaining the Reduced Row Echelon Form (RREF)


    EXAMPLES:

    ::

        sage: sz=2; Ha=HM(sz,sz,'a'); Hd=GaussJordanEliminationHM(Ha).canonicalize_radical()
        sage: Hd.printHM()
        [:, :]=
        [-a00*a01*a10 + a00^2*a11                        0]
        [                       0        a01*a10 - a00*a11]
        sage: Hb=GaussJordanEliminationHM(Ha); Hb[0,1].is_zero()
        True


    AUTHORS:
    - Edinah K. Gnang
    """
    A=Ha.copy()
    B=GaussEliminationHM(A).index_rotation(2*2*pi/4)
    for i in rg(B.n(0)-1):
        B=B.fill_with(GaussEliminationHM(B.slice(rg(i,A.n(0)),'row').slice(rg(i,A.n(1)),'col')),i)
    return B.index_rotation(2*2*pi/4)

def general_derivative(f, v, od):
    """
    Returns the vector of generalized derivatives for the input function
    f in the variable v of order od.


    EXAMPLES:
    ::

        sage: f=x^2
        sage: general_derivative(f, x, 2).printHM()
        [:, :]=
        [2*h^2 + 4*h*x + 2*x^2]
        [                4*h*x]
        sage: general_derivative(f, x, 3).printHM()
        [:, :]=
        [3*h^2 + 6*h*x + 3*x^2]
        [                3*h^2]
        [                3*h^2]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the step variable
    h=var('h')
    # Initialization of the DFT matrix
    Ha=SecondOrderDFT(od,od)[0]
    return HM(od, 1, [sum(Ha[i,j]*f.subs(v==v+Ha[i,j]*h) for j in rg(od)) for i in rg(od)]).expand()

def Grph_EliminationHMII(Ha):
    """
    Procedure for clearing bellow the first pivot.
    This implementation is motivated by functional directed graph listing


    EXAMPLES:

    ::

        sage: sz=2; Ha=HM(sz,sz,'a'); Hc=Grph_EliminationHMII(Ha)
        sage: Hc.printHM()
        [:, :]=
        [        2*a00*a10 a01*a10 + a00*a11]
        [                0 a01*a10 - a00*a11]
        sage: sz=3; Ha=HM(sz,sz,'a'); Hc=Grph_EliminationHMII(Ha)
        sage: Hc.printHM()
        [:, :]=
        [                                                                    3*a00*a10*a20                                           a01*a10*a20 + a00*a11*a20 + a00*a10*a21                                           a02*a10*a20 + a00*a12*a20 + a00*a10*a22]
        [-1/2*a00*a10*a20*(I*sqrt(3) + 1) - 1/2*a00*a10*a20*(-I*sqrt(3) + 1) + a00*a10*a20 -1/2*a00*a10*a21*(I*sqrt(3) + 1) - 1/2*a00*a11*a20*(-I*sqrt(3) + 1) + a01*a10*a20 -1/2*a00*a10*a22*(I*sqrt(3) + 1) - 1/2*a00*a12*a20*(-I*sqrt(3) + 1) + a02*a10*a20]
        [-1/2*a00*a10*a20*(I*sqrt(3) + 1) - 1/2*a00*a10*a20*(-I*sqrt(3) + 1) + a00*a10*a20 -1/2*a00*a11*a20*(I*sqrt(3) + 1) - 1/2*a00*a10*a21*(-I*sqrt(3) + 1) + a01*a10*a20 -1/2*a00*a12*a20*(I*sqrt(3) + 1) - 1/2*a00*a10*a22*(-I*sqrt(3) + 1) + a02*a10*a20]


    AUTHORS:
    - Edinah K. Gnang
    """
    A=Ha.copy(); B=Ha.copy()
    #for i in rg(1,A.n(0)):
    for i in rg(A.n(0)):
        # Zeroing out the first entry of row i
        Tmp=EliminationHMI(A,i)
        for jndx in rg(A.n(1)):
            B[i,jndx]=Tmp[i,jndx]
    return B

def Grph_GaussEliminationHM(Ha):
    """
    Procedure for obtaining the Row Echelon Form (REF)
    This implementation is motivated by functional directed graph listing

    EXAMPLES:

    ::

        sage: sz=3; Ha=HM(sz,sz,'a'); Hc=Grph_GaussEliminationHM(Ha)
        sage: Hc.printHM()
        [:, :]=
        [                                                                                                                                                                                                                                                                                                                 3*a00*a10*a20                                                                                                                                                                                                                                                                                        a01*a10*a20 + a00*a11*a20 + a00*a10*a21                                                                                                                                                                                                                                                                                        a02*a10*a20 + a00*a12*a20 + a00*a10*a22]
        [                                                                                                                                                                                                                                             -1/2*a00*a10*a20*(I*sqrt(3) + 1) - 1/2*a00*a10*a20*(-I*sqrt(3) + 1) + a00*a10*a20                                                                                                                                                                  1/2*(a00*a11*a20*(I*sqrt(3) + 1) + a00*a10*a21*(-I*sqrt(3) + 1) - 2*a01*a10*a20)*(a00*a10*a21*(I*sqrt(3) + 1) + a00*a11*a20*(-I*sqrt(3) + 1) - 2*a01*a10*a20)  1/4*(a00*a12*a20*(I*sqrt(3) + 1) + a00*a10*a22*(-I*sqrt(3) + 1) - 2*a02*a10*a20)*(a00*a10*a21*(I*sqrt(3) + 1) + a00*a11*a20*(-I*sqrt(3) + 1) - 2*a01*a10*a20) + 1/4*(a00*a11*a20*(I*sqrt(3) + 1) + a00*a10*a21*(-I*sqrt(3) + 1) - 2*a01*a10*a20)*(a00*a10*a22*(I*sqrt(3) + 1) + a00*a12*a20*(-I*sqrt(3) + 1) - 2*a02*a10*a20)]
        [                                                                                                                                                                                                                                             -1/2*a00*a10*a20*(I*sqrt(3) + 1) - 1/2*a00*a10*a20*(-I*sqrt(3) + 1) + a00*a10*a20                                                                                                                                                                                                                                                                                                                              0 -1/4*(a00*a12*a20*(I*sqrt(3) + 1) + a00*a10*a22*(-I*sqrt(3) + 1) - 2*a02*a10*a20)*(a00*a10*a21*(I*sqrt(3) + 1) + a00*a11*a20*(-I*sqrt(3) + 1) - 2*a01*a10*a20) + 1/4*(a00*a11*a20*(I*sqrt(3) + 1) + a00*a10*a21*(-I*sqrt(3) + 1) - 2*a01*a10*a20)*(a00*a10*a22*(I*sqrt(3) + 1) + a00*a12*a20*(-I*sqrt(3) + 1) - 2*a02*a10*a20)]



    AUTHORS:
    - Edinah K. Gnang
    """
    A=Ha.copy()
    for i in rg(A.n(0)-1):
        A = A.fill_with(Grph_EliminationHMII(A.slice(rg(i,A.n(0)),'row').slice(rg(i,A.n(1)),'col')),i)
    return A

def FindCycleTupleComponents(T):
    """
    Returns the cycle decomposition of the input
    permutation.


    EXAMPLES:

    ::


        sage: FindCycleTupleComponents([(0, 0), (1, 2), (2, 3), (3, 4), (4, 5), (5, 1), (6, 7), (7, 6)])
        [[(0, 0)], [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)], [(6, 7), (7, 6)]]


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of the temporary list
    Tp=[T[i] for i in rg(len(T))]
    # Initialization of the list of components
    cL=[]
    while len(Tp) > 0:
        # Testing for the cycle length equal to one
        if Tp[0][0] == Tp[0][1]:
            # Updating the temporary list
            cL.append([Tp[0]]); Tp.remove(Tp[0])
        # case where the cycle length equal > one
        else:
            tmp_vL=[Tp[0][0]]; tmp_cL=[Tp[0]]
            nv=tmp_cL[len(tmp_cL)-1][1]
            Tp.remove( (Tp[0][0], Tp[0][1]) )
            while nv not in tmp_vL:
                tmp_vL.append(tmp_cL[len(tmp_cL)-1][1]); tmp_cL.append(T[nv])
                Tp.remove(T[nv])
                nv=tmp_cL[len(tmp_cL)-1][1]
            cL.append(tmp_cL)
    return cL

def signf(T):
    """
    Returns the sign of a permutation inputed as
    a list of tuples of the form [(i, f(i))  for i in rg(sz)].


    EXAMPLES:
    ::
        sage: signf([(0, 0), (1, 1), (2, 2)])
        1
        sage: signf([(0, 0), (1, 2), (2, 3), (3, 4), (4, 5), (5, 1), (6, 7), (7, 6)])
        -1 


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the cylce decomposition of the permutation
    cL=FindCycleTupleComponents(T)
    return prod((-1)^Integer(mod(1+len(lst),2)) for lst in cL)

def mod_coefficients(F, d, Lx):
    """
    Returns the multivariate polynomial where all coefficients of the
    input multivariate polynomial F in the variables in the input list Lx
    are moded by in the integer input d.The implementation of this function
    is rather naive.
    This implementation expects the input polynomial to be in expanded form.


    EXAMPLES:
    ::
        sage: sz=3; A=HM(sz,sz,'a'); F=expand(prod(A[j,1]+A[i,1] for i in rg(sz) for j in rg(sz) if i<j))
        sage: Lv =[HM(*([sz]+[1 for k in rg(sz-1)]+[[A[sz-1-i,1]^j for j in rg(sz)]])) for i in rg(sz)]
        sage: prm = expand(prod(A[u,0] for u in rg(sz))*\
        ....: fast_reduce(F, [A[i,1]^(sz-j-1) for j in rg(sz-1) for i in rg(sz)], \
        ....: [A[i,sz-j-1]/A[i,0] for j in rg(sz-1) for i in rg(sz)]))
        sage: mod_coefficients(prm, 2, A.list())
        a02*a11*a20 + a01*a12*a20 + a02*a10*a21 + a00*a12*a21 + a01*a10*a22 + a00*a11*a22


    AUTHORS:

    - Edinah K. Gnang
    """
    add = var('x0') + var('x1')
    mul = var('x0') * var('x1')
    if F.operator() == add.operator():
        return sum(Integer(mod(mnm.subs([v==1 for v in Lx]),2))*(mnm/mnm.subs([v==1 for v in Lx])) for mnm in F.operands())
    elif F.operator() == mul.operator():
        return Integer(mod(F.subs([v==1 for v in Lx]),2))*(F/F.subs([v==1 for v in Lx]))

def HM2Poly(A, Lr):
    """
    Returns multivariate polynomial whose evaluation over 
    the grid specified by the cartesian product of Lr with
    itself yields entries of the hypermatrix.


    EXAMPLES:
    ::

        sage: A=HM(2,2,[1,0,0,1])
        sage: HM2Poly(A,rg(2))
        (x0 - 1)*(x1 - 1) + x0*x1

    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the list of variables and the size parameter
    X=var_list('x', A.order()); sz=len(X)
    return sum(A[tuple([Lr[t[i][1]] for i in rg(sz)])]*prod(prod((X[k]-jk)/(t[k][1]-jk) for jk in Lr if jk !=t[k][1]) for k in rg(sz) ) for t in TupleFunctionList(sz))

def svd_numeric2x2(A, dgts=15):
    """
    Returns the matrix of left singular vectors and the matrix of
    right singular vectors both with the diagonal matrix
    whose entries are [1, sigma1/sigma0].
    The format is [[U,Dg0],[V,Dg1]]
    This code is only needed when A is not rank deficient    


    EXAMPLES:
    ::

        sage: A=HM(2,2,[1,2,1,4]); svd_numeric2x2(A)
        [[[[1.35353252429997, -4.47041424613582], [4.47041424613582, 1.35353252429997]],
          [[1, 0], [0, 0.00840395484417606]]],
         [[[-0.377523729112191, 0.202044568264666], [-0.202044568264666, -0.377523729112191]],
          [[1, 0], [0, 118.991596045157]]]]


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the size parameter
    sz=2; z=var('z')
    # Initialization of the first product of transposes
    AAt=A*A.transpose()
    # Initialization of the variables
    x00,x10=var('x00,x10')
    # Obtaining the constraints for the e-vectors
    EqL0=[\
    (AAt-HM(sz,1,[x00,x10])*HM(1,sz,[x00,x10])-z*HM(sz,1,[-x10,x00])*HM(1,sz,[-x10,x00]))[0,0],\
    (AAt-HM(sz,1,[x00,x10])*HM(1,sz,[x00,x10])-z*HM(sz,1,[-x10,x00])*HM(1,sz,[-x10,x00]))[0,1],\
    (AAt-HM(sz,1,[x00,x10])*HM(1,sz,[x00,x10])-z*HM(sz,1,[-x10,x00])*HM(1,sz,[-x10,x00]))[1,1]]
    # Putting the system in row echelon form
    L0=[expand(p) for p in eulerian_eliminationHM(EqL0,[z,x00,x10])]
    # Initialization of the degree matrix
    M0=degree_matrix(L0, [z,x00,x10])
    # Solution to the equation in x10
    Sln_x10=[N(eq.rhs(), digits=dgts) for eq in solve(L0[2]==0, x10) if eq.rhs()!=0]
    Sln_x00=[]
    for u10 in Sln_x10:
        Sln_x00=Sln_x00+[N(eq.rhs(), digits=dgts) for eq in solve(L0[1].subs(x10==u10)==0, x00) if eq.rhs()!=0]
    Sln_z0=[]
    for u10 in Sln_x10:
        for u00 in Sln_x00:
            Sln_z0=Sln_z0+[N(eq.rhs(), digits=dgts) for eq in solve(L0[0].subs([x10==u10, x00==u00])==0, z) if eq.rhs()!=0]
    min_norm0=AAt.norm()
    # Performing the checking
    for u00 in Sln_x00:
        for u10 in Sln_x10:
            for z0 in Sln_z0:
                if (AAt-HM(sz,1,[u00,u10])*HM(1,sz,[u00,u10])-z0*HM(sz,1,[-u10,u00])*HM(1,sz,[-u10,u00])).norm() < min_norm0:
                    min_norm0=(AAt-HM(sz,1,[u00,u10])*HM(1,sz,[u00,u10])-z0*HM(sz,1,[-u10,u00])*HM(1,sz,[-u10,u00])).norm()
                    U=HM([[u00, -u10], [u10, u00]]); Dg0=HM(2,[1,z0],'diag')
    # Initialization of the second product of transposes
    AtA=A.transpose()*A
    # Initialization of the variables
    y00,y01=var('y00,y01')
    # Obtaining the constraints for the e-vectors
    EqL1=[\
    (AtA-HM(sz,1,[y00,y01])*HM(1,sz,[y00,y01])-z*HM(sz,1,[-y01,y00])*HM(1,sz,[-y01,y00]))[0,0],\
    (AtA-HM(sz,1,[y00,y01])*HM(1,sz,[y00,y01])-z*HM(sz,1,[-y01,y00])*HM(1,sz,[-y01,y00]))[0,1],\
    (AtA-HM(sz,1,[y00,y01])*HM(1,sz,[y00,y01])-z*HM(sz,1,[-y01,y00])*HM(1,sz,[-y01,y00]))[1,1]]
    # Putting the system in row echelon form
    L1=[expand(p) for p in eulerian_eliminationHM(EqL1,[z,y00,y01])]
    # Initialization of the degree matrix
    M1=degree_matrix(L1, [z,y00,y01])
    # Solution to the equation in y01
    Sln_y01=[N(eq.rhs(), digits=dgts) for eq in solve(L1[2]==0, y01) if eq.rhs()!=0]
    Sln_y00=[]
    for v01 in Sln_y01:
        Sln_y00=Sln_y00+[N(eq.rhs(), digits=dgts) for eq in solve(L1[1].subs(y01==v01)==0, y00) if eq.rhs()!=0]
    Sln_z1=[]
    for v01 in Sln_y01:
        for v00 in Sln_y00:
            Sln_z1=Sln_z1+[N(eq.rhs(), digits=dgts) for eq in solve(L1[0].subs([y01==v01, y00==v00])==0, z) if eq.rhs()!=0]
    min_norm1=AtA.norm()
    # Performing the checking
    for v00 in Sln_y00:
        for v01 in Sln_y01:
            for z1 in Sln_z1:
                if (AtA-HM(sz,1,[v00,v01])*HM(1,sz,[v00,v01])-z1*HM(sz,1,[-v01,v00])*HM(1,sz,[-v01,v00])).norm() < min_norm1:
                    min_norm1=(AtA-HM(sz,1,[v00,v01])*HM(1,sz,[v00,v01])-z1*HM(sz,1,[-v01,v00])*HM(1,sz,[-v01,v00])).norm()
                    V=HM([[v00, v01], [-v01, v00]]); Dg1=HM(2,[1,z1],'diag')
    return [[U,Dg0],[V,Dg1]]

def svd_symbolic2x2(A):
    """
    Returns the solutions to the algebraic constraints associated with
    the SVD decomposition of the 2x2 matrix.
    


    EXAMPLES:
    ::

        sage: A=HM(2,2,[1,2,1,4]); svd_symbolic2x2(A)[0][0]
        x00 == 1/39*(sqrt(13)*(5*sqrt(13) + 18) - 13)*(-I*sqrt(3) + 1)/(-3*13^(1/4)*sqrt(2/13)*sqrt(5*sqrt(13) + 18) + 1/39*sqrt(2/3)*sqrt(4*sqrt(13)*(5*sqrt(13) + 18)^3 - 156*(5*sqrt(13) + 18)^2 + 156*sqrt(13)*(5*sqrt(13) + 18) + 56862*sqrt(13) + 204659))^(1/3) - 1/2*(-3*13^(1/4)*sqrt(2/13)*sqrt(5*sqrt(13) + 18) + 1/39*sqrt(2/3)*sqrt(4*sqrt(13)*(5*sqrt(13) + 18)^3 - 156*(5*sqrt(13) + 18)^2 + 156*sqrt(13)*(5*sqrt(13) + 18) + 56862*sqrt(13) + 204659))^(1/3)*(I*sqrt(3) + 1)


    AUTHORS:

    - Edinah K. Gnang
    """
    # Initialization of the size parameter
    sz=2; z=var('z')
    # Initialization of the first product of transposes
    AAt=A*A.transpose()
    # Initialization of the variables
    x00,x10=var('x00,x10')
    # Obtaining the constraints for the e-vectors
    EqL0=[\
    (AAt-HM(sz,1,[x00,x10])*HM(1,sz,[x00,x10])-z*HM(sz,1,[-x10,x00])*HM(1,sz,[-x10,x00]))[0,0],\
    (AAt-HM(sz,1,[x00,x10])*HM(1,sz,[x00,x10])-z*HM(sz,1,[-x10,x00])*HM(1,sz,[-x10,x00]))[0,1],\
    (AAt-HM(sz,1,[x00,x10])*HM(1,sz,[x00,x10])-z*HM(sz,1,[-x10,x00])*HM(1,sz,[-x10,x00]))[1,1]]
    # Putting the system in row echelon form
    L0=[expand(p) for p in eulerian_eliminationHM(EqL0,[z,x00,x10])]
    # Initialization of the degree matrix
    M0=degree_matrix(L0, [z,x00,x10])
    # Solution to the equation in x10
    Sln_x10=[eq for eq in solve(L0[2]==0, x10) if eq.rhs()!=0]
    Sln_x00=[]
    for u10 in Sln_x10:
        Sln_x00=Sln_x00+[eq for eq in solve(L0[1].subs(x10==u10.rhs())==0, x00) if eq.rhs()!=0]
    Sln_z0=[]
    for u10 in Sln_x10:
        for u00 in Sln_x00:
            Sln_z0=Sln_z0+[eq for eq in solve(L0[0].subs([x10==u10.rhs(), x00==u00.rhs()])==0, z) if eq.rhs()!=0]
    # Initialization of the second product of transposes
    AtA=A.transpose()*A
    # Initialization of the variables
    y00,y01=var('y00,y01')
    # Obtaining the constraints for the e-vectors
    EqL1=[\
    (AtA-HM(sz,1,[y00,y01])*HM(1,sz,[y00,y01])-z*HM(sz,1,[-y01,y00])*HM(1,sz,[-y01,y00]))[0,0],\
    (AtA-HM(sz,1,[y00,y01])*HM(1,sz,[y00,y01])-z*HM(sz,1,[-y01,y00])*HM(1,sz,[-y01,y00]))[0,1],\
    (AtA-HM(sz,1,[y00,y01])*HM(1,sz,[y00,y01])-z*HM(sz,1,[-y01,y00])*HM(1,sz,[-y01,y00]))[1,1]]
    # Putting the system in row echelon form
    L1=[expand(p) for p in eulerian_eliminationHM(EqL1,[z,y00,y01])]
    # Initialization of the degree matrix
    M1=degree_matrix(L1, [z,y00,y01])
    # Solution to the equation in y01
    Sln_y01=[eq for eq in solve(L1[2]==0, y01) if eq.rhs()!=0]
    Sln_y00=[]
    for v01 in Sln_y01:
        Sln_y00=Sln_y00+[eq for eq in solve(L1[1].subs(y01==v01.rhs())==0, y00) if eq.rhs()!=0]
    Sln_z1=[]
    for v01 in Sln_y01:
        for v00 in Sln_y00:
            Sln_z1=Sln_z1+[eq for eq in solve(L1[0].subs([y01==v01.rhs(), y00==v00.rhs()])==0, z) if eq.rhs()!=0]
    return [Sln_x00,Sln_x10,Sln_z0,Sln_y00,Sln_y01,Sln_z1]


def EuclidsPolynomialGCD(a, b, v):
    """

    This function implements Euclid's GCD algorithm for
    polynomials.
    The function checks that the inputs are not degree 0
    polynomials and returns as output the matrix which
    describes all the iterations of the algorithm.


    EXAMPLES:
    ::


        sage: p=x^5+1; d=x^2+1; G=EuclidsPolynomialGCD(p, d, x); G.printHM()
        [:, :]=
        [x^5 + 1 x^2 + 1 x^3 - x   x + 1]
        [x^2 + 1   x + 1   x - 1       2]
        sage: p=x^5-1; d=x^2-1; G=EuclidsPolynomialGCD(p, d, x); G.printHM()
        [:, :]=
        [x^5 - 1 x^2 - 1 x^3 + x   x - 1]
        [x^2 - 1   x - 1   x + 1       0]


    AUTHORS:
    - Edinah K. Gnang
    """
    if a.degree(v)>0 and b.degree(v)>0:
        # Initialization of the matrix Data Structure.
        G = HM(1, 4, 'zero')
        # Initialization of the initial conditions
        G[0, 0] = a; G[0, 1] = b; [q, r] = Division(G[0, 0], G[0, 1], v)
        G[0, 3] = r; G[0, 2] = q
        # Initialization of the index
        indx = 0
        while G[indx, 3].degree(v) > 0:
            # Updating the size of G
            G=G.zero_pad([G.n(0)+1, G.n(1)])
            # Incrementing the index
            indx=indx+1
            G[indx, 0] = G[indx-1, 1]; G[indx, 1] = G[indx-1, 3]; [q, r]=Division(G[indx, 0], G[indx, 1], v)
            G[indx, 3] = r; G[indx, 2] = q
        return G
    else:
        raise ValueError("Expected inputs of degree >= 1.")

def CompositionalDivision(p, d,  v):
    """
    Outputs the quotient list and the remainder of the composition version
    the Euclidean division. The algorithm takes as input
    two univariate polynomials in the input variable v,  p and d respectively
    associated with the dividend and the divisor. 


    EXAMPLES:

    ::
        
        sage: p = 4*x^3+3*x^2+2*x+1; d = 5*x^2+3*x+7; [Lq, r] = CompositionalDivision(p, d, x); [Lq, r] 
        [[2*sqrt(1/5)*x^(3/2), sqrt(3/5)*x],
         -1/5*(3*sqrt(5)*sqrt(3) - 10)*x - 6/5*sqrt(5)*x^(3/2) - 13]
        sage: (sum(d.subs(x==q) for q in Lq)+r).canonicalize_radical() # Checking the computation
        4*x^3 + 3*x^2 + 2*x + 1
        sage: p=13*x^7; d=5*x^2+3*x+7; [Lq,r] = CompositionalDivision(p, d, x); [Lq, r]
        [[sqrt(13/5)*x^(7/2), 1/5*sqrt(3)*sqrt(-sqrt(13)*sqrt(5))*x^(7/4)],
         -3/5*5^(1/4)*sqrt(3)*x^(7/4)*sqrt(-sqrt(13)) - 14]
        sage: (sum(d.subs(x==q) for q in Lq)+r).canonicalize_radical() # Checking the computation
        13*x^7
        sage: p=x^7; d=5*x^2+3*x+7; [Lq,r] = CompositionalDivision(p, d, x); [Lq,r]
        [[sqrt(1/5)*x^(7/2), 1/5*sqrt(3)*x^(7/4)*sqrt(-sqrt(5))],
         -3/5*sqrt(3)*x^(7/4)*sqrt(-sqrt(5)) - 14]
        sage: (sum(d.subs(x==q) for q in Lq)+r).canonicalize_radical() # Checking the computation
        x^7
        sage: p=x^7; d=5*x^2; [Lq,r] = CompositionalDivision(p, d, x); [Lq,r]
        [[sqrt(1/5)*x^(7/2)], 0]
        sage: sum(d.subs(x==q) for q in Lq)+r # Checking the computation
        x^7
        sage: p=x*(var('t')+1)*(var('t')+2); d=x*(var('t')+1); [Lq,r] = CompositionalDivision(p, d, x); [Lq,r]
        [[(t + 2)*x], 0]

    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of remainder
    r = p
    # Initialization of the Quotient list
    Lq = []
    # Checking that neither p or d are monomials.
    if p.degree(v)>=1 and d.degree(v)>=1:
        # Obtaining the leading term of d
        ltd = d.coefficients(v)[len(d.coefficients(v))-1][0]* v^d.coefficients(v)[len(d.coefficients(v))-1][1]
        # Main loop
        while r.degree(v) >= d.degree(v):
            # Obtaining the leading term of r
            ltr = r.coefficients(v)[len(r.coefficients(v))-1][0]* v^r.coefficients(v)[len(r.coefficients(v))-1][1]
            # Updating the quotient list
            Lq.append((ltr.subs(v==1)/ltd.subs(v==1))^(1/ltd.degree(v))*v^(ltr.degree(v)/ltd.degree(v)))
            #print 'Lq = ',Lq
            # Updating the remainder
            #r = r-d.subs(v==Lq[len(Lq)-1])
            #r = (r-d.subs(v==Lq[len(Lq)-1])).canonicalize_radical()
            r = sum( l[0]*v^l[1] for l in (r-d.subs(v==Lq[len(Lq)-1])).canonicalize_radical().coefficients() )
            #print 'r = ',r
        return [Lq, r]
    else:
        raise ValueError("Expected inputs of degree >= 1 in the input variable.")

def Division(p, d,  v):
    """
    Outputs the quotient and the remainder of the Euclidean division.
    The algorithm takes as input two univariate polynomials in the
    input variable v,  p and d respectively associated with the 
    dividend and the divisor. 


    EXAMPLES:

    ::
        
        sage: p = 4*x^3+3*x^2+2*x+1; d = 5*x^2+3*x+7; [q, r] = Division(p, d, x); [q, r]
        [4/5*x + 3/25, -99/25*x + 4/25]
        sage: expand(d*q+r) # Checking the computation
        4*x^3 + 3*x^2 + 2*x + 1
        sage: p=13*x^7; d=5*x^2+3*x+7; [q,r] = Division(p, d, x); [q, r]
        [13/5*x^5 - 39/25*x^4 - 338/125*x^3 + 2379/625*x^2 + 4693/3125*x - 97344/15625,
         127777/15625*x + 681408/15625]
        sage: expand(d*q+r) # Checking the computation
        13*x^7
        sage: p=x^7; d=5*x^2+3*x+7; [q,r] = Division(p, d, x); [q,r]
        [1/5*x^5 - 3/25*x^4 - 26/125*x^3 + 183/625*x^2 + 361/3125*x - 7488/15625,
         9829/15625*x + 52416/15625]
        sage: expand(d*q+r) # Checking the computation
        x^7
        sage: p=x^7; d=5*x^2; [q,r] = Division(p, d, x); [q,r]
        [1/5*x^5, 0]
        sage: expand(d*q+r) # Checking the computation
        x^7


    AUTHORS:
    - Edinah K. Gnang
    """
    # Initialization of remainder
    r = p
    # Initialization of the Quotient list
    Lq = []
    # Checking that neither p or d are monomials.
    if p.degree(v)>=1 and d.degree(v)>=1:
        # Obtaining the leading term of d
        ltd = d.coefficients(v)[len(d.coefficients(v))-1][0]* v^d.coefficients(v)[len(d.coefficients(v))-1][1]
        # Main loop
        while r.degree(v) >= d.degree(v):
            # Obtaining the leading term of r
            ltr = r.coefficients(v)[len(r.coefficients(v))-1][0]* v^r.coefficients(v)[len(r.coefficients(v))-1][1]
            # Updating the quotient list
            Lq.append((ltr.subs(v==1)/ltd.subs(v==1))*v^(ltr.degree(v)-ltd.degree(v)))
            #print 'Lq = ',Lq
            # Updating the remainder
            r = expand(r-d*Lq[len(Lq)-1])
            #print 'r = ',r
        return [sum(Lq), r]
    else:
        raise ValueError("Expected inputs of degree >= 1 in the input variable.")

def EuclidsCompositionalGCD(a, b, v):
    """

    This function implements the composition version
    of the Euclid's GCD algorithm.
    The function checks that the inputs are polynomials
    of degree >1 and returns as output the matrix which
    describes all the iterations of the algorithm.


    EXAMPLES:
    ::


        sage: p=x^5+1; d=x^2+1; G=EuclidsCompositionalGCD(p, d, x); G.slice([0,1,3],'col').printHM()
        [:, :]=
        [x^5 + 1 x^2 + 1       0]
        sage: p=x^5-1; d=x^2-1; G=EuclidsCompositionalGCD(p, d, x); G.slice([0,1,3],'col').printHM()
        [:, :]=
        [x^5 - 1 x^2 - 1       0]
        sage: p=x^5+1; d=x^3+2*x^2+1; G=EuclidsCompositionalGCD(p, d, x); G.slice([0,1,3],'col').printHM() 
        [:, :]=
        [                           x^5 + 1                    x^3 + 2*x^2 + 1 -2*2^(2/3)*(-1)^(2/3)*x^(20/9) - 1]
        [                   x^3 + 2*x^2 + 1 -2*2^(2/3)*(-1)^(2/3)*x^(20/9) - 1                          2*x^2 + 2]
        [-2*2^(2/3)*(-1)^(2/3)*x^(20/9) - 1                          2*x^2 + 2                                 -3]


    AUTHORS:
    - Edinah K. Gnang
    """
    if a.degree(v)>1 and b.degree(v)>1:
        # Initialization of the matrix Data Structure.
        G = HM(1, 4, 'zero')
        # Initialization of the initial conditions
        G[0, 0] = a; G[0, 1] = b
        # Performing the compositional division
        [Lq, r] = CompositionalDivision(G[0, 0], G[0, 1], v)
        G[0, 3] = r.canonicalize_radical(); G[0, 2] = Lq
        #print G.slice([0],'row')
        # Initialization of the index
        indx = 0
        while G[indx, 3].degree(v) > 1:
            #print G.slice([G.n(0)-1],'row')
            # Updating the size of G
            G=G.zero_pad([G.n(0)+1, G.n(1)])
            # Incrementing the index
            indx=indx+1
            G[indx, 0] = G[indx-1, 1]; G[indx, 1] = G[indx-1, 3]
            #print G.slice([G.n(0)-1],'row')
            # Performing the compositional division
            [Lq, r] = CompositionalDivision(G[indx, 0], G[indx, 1], v)
            G[indx, 3] = r.canonicalize_radical(); G[indx, 2] = Lq
            #print G.slice([G.n(0)-1],'row')
        return G
    else:
        raise ValueError("Expected inputs of degree > 1.")

