=======
# Hypermatrix Algebra Package

We provide here a sagemath implementation of the Bhattacharya-Mesner algebra as well as the algebra of constructs.

The `Hypermatrix Algebra Package` is a symbolic hypermatrix package designed to experimentally investigate symbolic
structural and combinatorial properties of the Bhattacharya-Mesner algebra and the algebra of constructs more generally.

# Installation 

A working install of [sage](http://sagemath.org/) is the only prerequisite to using the hypermatrix package. To use the hypermatrix
algebra package, simply download the [hypermatrix algebra package sage file](https://github.com/gnang/HypermatrixAlgebraPackage/blob/master/Hypermatrix_Algebra_Package_code.sage) into your working directory and load the file into your SageMath [interactive shell](http://doc.sagemath.org/html/en/tutorial/interactive_shell.html) session via the command:

```python
sage: load('Hypermatrix_Algebra_Package_code.sage')
```

The hypermatrix package passes all tests on SageMath version 9.3 To run the hypermatrix tests on your local sage, 
use the following command in a [bash](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) shell session

```bash
n: g$ sage -t Hypermatrix_Algebra_Package_code.sage
```

To get started with SageMath, the authors of the hypermatrix package highly recommend reading 
[Mathematical Computation with SageMath] (http://sagebook.gforge.inria.fr/)

The python file can be created from the sage file in the bash terminal as follows

```bash
n: g$ sage -preparse Hypermatrix_Algebra_Package_code.sage; mv Hypermatrix_Algebra_Package_code.sage.py Hypermatrix_Algebra_Package_code.py
```

The functions implemented in the hypermatrix package can also be imported as follows

```python
sage: from Hypermatrix_Algebra_Package_code import *
```


# Usage

To create a symbolic hypermatrix of size say 2 by 3 by 4, use

```python
sage: Ha=HM(2,3,4,'a')
sage: Ha.printHM()
[:, :, 0]=
[a000 a010 a020]
[a100 a110 a120]

[:, :, 1]=
[a001 a011 a021]
[a101 a111 a121]

[:, :, 2]=
[a002 a012 a022]
[a102 a112 a122]

[:, :, 3]=
[a003 a013 a023]
[a103 a113 a123]
```

Alternatively, hypermatrices can be initialized from an arbitrary Python list as follows

```python
sage: rng=range(1,3); Lst=[var('x{0}{1}{2}'.format(i,j,k)) for k in rng for j in rng for i in rng]
sage: Hx=HM(2,2,2,Lst)
sage: Hx.printHM()
[:, :, 0]=
[x111 x121]
[x211 x221]

[:, :, 1]=
[x112 x122]
[x212 x222]
```

Note that hypermatrix entries are not restricted to numbers and symbolic expressions. The hypermatrix 
entries may be of other sage object types or mixture of different sage object types including list, Matrix,
Boolean or even hypermatrices as illustrated below

```python
sage: rng = range(2)
sage: Hh=HM(2,2,2,[(i*2^2+j*2+k)*hadamard_matrix(2) for k in rng for j in rng for i in rng])
```

Hypermatrix entries are accessed the same way matrix entries are usually accessed in sage.
For example, the (0,1,0) entry of the hypermatrix constructed in the previous instruction
is accessed as follows

```python
sage: Hh[0,1,0]
[ 2  2]
[ 2 -2]
```

# Basic hypermatrix operation

The [BM product](http://arxiv.org/abs/1411.6270) is implemented for conformable hypermatrices
of all orders. The BM product is computed as follows

```python
sage: H1=Prod(HM(3,2,'a'), HM(2,2,'c')) # Matrix product
sage: H1.printHM()
[:, :]=
[a00*c00 + a01*c10 a00*c01 + a01*c11]
[a10*c00 + a11*c10 a10*c01 + a11*c11]
[a20*c00 + a21*c10 a20*c01 + a21*c11]

sage: H2=Prod(HM(2,3,4,'a'), HM(2,2,3,'b'), HM(3,2,4,'c')) # Third order product 
sage: H2.printHM()
[:, :, 0]=
[a000*b000*c000 + a010*b001*c100 + a020*b002*c200 a000*b010*c010 + a010*b011*c110 + a020*b012*c210]
[a100*b100*c000 + a110*b101*c100 + a120*b102*c200 a100*b110*c010 + a110*b111*c110 + a120*b112*c210]

[:, :, 1]=
[a001*b000*c001 + a011*b001*c101 + a021*b002*c201 a001*b010*c011 + a011*b011*c111 + a021*b012*c211]
[a101*b100*c001 + a111*b101*c101 + a121*b102*c201 a101*b110*c011 + a111*b111*c111 + a121*b112*c211]

[:, :, 2]=
[a002*b000*c002 + a012*b001*c102 + a022*b002*c202 a002*b010*c012 + a012*b011*c112 + a022*b012*c212]
[a102*b100*c002 + a112*b101*c102 + a122*b102*c202 a102*b110*c012 + a112*b111*c112 + a122*b112*c212]

[:, :, 3]=
[a003*b000*c003 + a013*b001*c103 + a023*b002*c203 a003*b010*c013 + a013*b011*c113 + a023*b012*c213]
[a103*b100*c003 + a113*b101*c103 + a123*b102*c203 a103*b110*c013 + a113*b111*c113 + a123*b112*c213]
```

As illustrated, the product of conformable matrices recovers the matrix product.
Other basic operations including addition and multiplication by scalars
are analogous to their matrix counterpart in sage. The hypermatrix transpose 
generalizes the matrix transpose and performs a cyclic permutation of the entry
indices and is performed as follows 

```python
sage: HM(2,2,2,'a').printHM()
[:, :, 0]=
[a000 a010]
[a100 a110]

[:, :, 1]=
[a001 a011]
[a101 a111]

sage: HM(2,2,2,'a').transpose().printHM()
[:, :, 0]=
[a000 a001]
[a010 a011]

[:, :, 1]=
[a100 a101]
[a110 a111]

sage: HM(2,2,2,'a').transpose().transpose().printHM()
[:, :, 0]=
[a000 a100]
[a001 a101]

[:, :, 1]=
[a010 a110]
[a011 a111]
```

To perform two or more consecutive transposes, simply provide as argument
the number of times we wish the transpose to be performed as follows

```python
sage: At=HM(2,2,2,2,'a').transpose()
sage: Att=HM(2,2,2,2,'a').transpose(2)
sage: Attt=HM(2,2,2,2,'a').transpose(3)
sage: (HM(2,2,2,2,'a').transpose(4)-HM(2,2,2,2,'a')).is_zero()
True
```

Hypermatrix sums can be taken over a lists as follows.

```python
sage: A=sum(HM(2,2,2,'a').transpose(i) for i in range(3))
sage: A.is_symmetric()
True
```

The hypermatrix package provides implementations of the Kronecker product for hypermatrices
of all orders. The Kronecker product of two hypermatrices are obtained as follows

```python
sage: Ha=HM(2,2,2,'a'); Hb=HM(3,3,3,'b')
sage: Hc=Ha.tensor_product(Hb)
sage: Hc.dimensions()
[6, 6, 6]
sage: Ha=HM(2,2,2,'a'); Hb=HM(1,1,1,'b')
sage: Ha.tensor_product(Hb).printHM()
[:, :, 0]=
[a000*b000 a010*b000]
[a100*b000 a110*b000]

[:, :, 1]=
[a001*b000 a011*b000]
[a101*b000 a111*b000]
```

Similarly the block diagonal hypermatrix sum is performed as follows

```python
sage: Ha=HM(2,2,'a'); Hb=HM(3,3,'b')
sage: Ha.block_sum(Hb).printHM()
[:, :]=
[a00 a01   0   0   0]
[a10 a11   0   0   0]
[  0   0 b00 b01 b02]
[  0   0 b10 b11 b12]
[  0   0 b20 b21 b22]

sage: Ha=HM(2,2,2,'a'); Hb=HM(2,2,2,'b')
sage: Ha.block_sum(Hb).printHM()
[:, :, 0]=
[a000 a010    0    0]
[a100 a110    0    0]
[   0    0    0    0]
[   0    0    0    0]

[:, :, 1]=
[a001 a011    0    0]
[a101 a111    0    0]
[   0    0    0    0]
[   0    0    0    0]

[:, :, 2]=
[   0    0    0    0]
[   0    0    0    0]
[   0    0 b000 b010]
[   0    0 b100 b110]

[:, :, 3]=
[   0    0    0    0]
[   0    0    0    0]
[   0    0 b001 b011]
[   0    0 b101 b111]
```

Let X, Y denote 3 by 1 hypermatrices, the product of the transpose of X with Y
is obtained as follows

```python
sage: sz=3; X=HM(sz, 1, HM(sz,'x').list()); Y=HM(sz, 1, HM(sz,'y').list()); L=[X, Y]
sage: Prod(X.transpose(),Y)[0,0]
x0*y0 + x1*y1 + x2*y2
sage: Prod(*[L[i].transpose(i) for i in range(X.order()-1,-1,-1)])[0,0]
x0*y0 + x1*y1 + x2*y2
```

Here is an illustration of a similar product in the hypermatrix case for columns X, Y, Z
each of size 3 by 1 by 1

```python
sage: sz=3; X=HM(sz,1,1,HM(sz,'x').list()); Y=HM(sz,1,1,HM(sz,'y').list()); Z=HM(sz,1,1,HM(sz,'z').list()); L=[X, Y, Z]
sage: Prod(X.transpose(2), Y.transpose(), Z)[0,0,0]
x0*y0*z0 + x1*y1*z1 + x2*y2*z2
sage: Prod(*[L[i].transpose(i) for i in range(X.order()-1,-1,-1)])[0,0,0]
x0*y0*z0 + x1*y1*z1 + x2*y2*z2
```

Here is a similar command in a more general set up

```python
sage: AlphaB=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
sage: od=4; sz=3 # Initialization of the order and size paramaters
sage: L=[HM(*([sz]+[1 for i in range(od-1)]+[HM(sz,AlphaB[j-od]).list()])) for j in range(od)]
sage: Prod(*[L[i].transpose(i) for i in range(od-1,-1,-1)])
[[[[w0*x0*y0*z0 + w1*x1*y1*z1 + w2*x2*y2*z2]]]]
```


Multilinear forms associated with an input hypermatrix is conveniently expressed
using the general BM product as follows

```python
sage: X=HM(2,1,HM(2,'x').list()); Y=HM(2,1,HM(2,'y').list()); A=HM(2,2,'a')
sage: ProdB(X.transpose(), Y, A)[0,0]
a00*x0*y0 + a10*x1*y0 + a01*x0*y1 + a11*x1*y1
```

The corresponding commands for third order hypermatrices is

```python
sage: X=HM(2,1,1,HM(2,'x').list()); Y=HM(2,1,1,HM(2,'y').list()); Z=HM(2,1,1,HM(2,'z').list())
sage: A=HM(2,2,2,'a')
sage: ProdB(X.transpose(2), Y.transpose(), Z, A)[0,0,0]
a000*x0*y0*z0 + a100*x1*y0*z0 + a010*x0*y1*z0 + a110*x1*y1*z0 + a001*x0*y0*z1 + a101*x1*y0*z1 + a011*x0*y1*z1 + a111*x1*y1*z1
```

Here is the same command in a more general setup

```python
sage: AlphaB=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
sage: od=3; sz=2 # Initialization of the order and size paramaters
sage: L=[HM(*([sz]+[1 for i in range(od-1)]+[HM(sz,AlphaB[j-od]).list()])) for j in range(od)]
sage: A=HM(2,2,2,'a') # Initialization of the background hypermatrix
sage: ProdB(*([L[i].transpose(i) for i in range(od-1,-1,-1)]+[A]))[0,0,0]
a000*x0*y0*z0 + a100*x1*y0*z0 + a010*x0*y1*z0 + a110*x1*y1*z0 + a001*x0*y0*z1 + a101*x1*y0*z1 + a011*x0*y1*z1 + a111*x1*y1*z1
```

Hypermatrices can be extracted from multilinear forms as follows

```python
sage: sz=2; od=2; X=HM(sz,sz,HM(sz^2,'x').list()); f=X.det()
sage: H=Form2TotallySymmetricHypermatrix(f, 2, X.list()); H.printHM()
[:, :]=
[   0    0    0  1/2]
[   0    0 -1/2    0]
[   0 -1/2    0    0]
[ 1/2    0    0    0]
```

# Symbolic parametrization of some special hypermatrices

We further illustrate here the connection between matrices and hypermatrices.

## Kronecker delta hypermatrices
The [Kronecker delta](http://en.wikipedia.org/wiki/Kronecker_delta#Properties_of_generalized_Kronecker_delta)
hypermatrices are generalizations of identity matrices in the sense that all entries are zero 
except for the entries located on the main diagonal as illustated below 
```python
sage: sz=2; od=3; rng=range(2) # Initialization of the size and order parameters
sage: Dlt=HM(od, sz, 'kronecker')
sage: A=HM(sz,sz,sz,'d')
sage: He=HM(sz,sz,sz,[A[i,j,k]==Dlt[i,j,k] for k in rng for j in rng for i in rng])
sage: He.list() # Printing constraints associated with the Kronecker delta
[d000 == 1,
 d100 == 0,
 d010 == 0,
 d110 == 0,
 d001 == 0,
 d101 == 0,
 d011 == 0,
 d111 == 1]
```

## Orthogonal hypermatrices

Orthogonal hypermatrices are analogous to orthogonal matrices in the fact that the product of the transposes 
yields the Kronecker delta hypermatrix. A symbolic parametetrization of 2 by 2 by 2 orthogonal hypermatrices is 
obtained as follows

```python
sage: Q=GeneralOrthogonalHypermatrix(3)
sage: (Prod(*[Q.transpose(i) for i in range(Q.order(),0,-1)]).simplify()-HM(Q.order(),2,'kronecker')).is_zero()
True
```

Similarly a symbolic parametrization of orthogonal 2 by 2 by 2 by 2 hypermatrices is obtained
as follows

```python
sage: Q=GeneralOrthogonalHypermatrix(4)
sage: (Prod(*[Q.transpose(i) for i in range(Q.order(),0,-1)]).simplify()-HM(Q.order(),2,'kronecker')).is_zero()
True
```

The parameters which appear in the parametrization can be assigned arbitrary complex values. The hypermatrix 
package also provides constructions of third order hypermatrix analog of Hadamard matrices, 
```python
sage: H=ThirdOrderHadamardBlockU(4); H
[[[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]], [[-1, 1, -1, 1], [1, 1, 1, 1], [-1, 1, -1, 1], [1, 1, 1, 1]], [[-1, -1, 1, 1], [-1, -1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]], [[1, -1, -1, 1], [-1, -1, 1, 1], [-1, 1, -1, 1], [1, 1, 1, 1]]]
sage: (Prod(*[H.transpose(i) for i in range(H.order(),0,-1)])-4*HM(H.order(),4,'kronecker')).is_zero()
True
```


## Unitary hypermatrices

The hypermatrix algebra package also implements a symbolic parametrization of even order unitary hypermatrices having 
of sides length equal to 2. As suggested by the [Gowers norm](http://en.wikipedia.org/wiki/Gowers_norm) operation, an 
even order hypermatrix is said to be unitary if the product of the conjugate transposes equals the Kronecker delta.
Symbolic parametrization of 2 by 2 unitary matrices is obtained as follows

```python
sage: [U,Uc]=GeneralUnitaryHypermatrix(2)
sage: U
[[e^(-I*pi-I*r1+I*r2+I*r3-r4+r5+r6)/sqrt(e^(-2*r4+2*r5+2*r6)+e^(2*r6)), e^(I*r2+r6)/sqrt(e^(-2*r4+2*r5+2*r6)+e^(2*r6))], [e^(I*r3+r4)/sqrt(e^(2*r4)+e^(2*r5)), e^(I*r1+r5)/sqrt(e^(2*r4)+e^(2*r5))]]
sage: Uc # Symbolic complex conjugate of U
[[e^(I*pi + I*r1 - I*r2 - I*r3 - r4 + r5 + r6)/sqrt(e^(-2*r4 + 2*r5 + 2*r6) + e^(2*r6)), e^(-I*r2 + r6)/sqrt(e^(-2*r4 + 2*r5 + 2*r6) + e^(2*r6))], [e^(-I*r3 + r4)/sqrt(e^(2*r4) + e^(2*r5)), e^(-I*r1 + r5)/sqrt(e^(2*r4) + e^(2*r5))]]
sage: Prod(U,Uc.transpose()).simplify()
[[1, 0], [0, 1]]
```

Similarly a symbolic paramterization of a fourth order unitary hypermatrix of side length 2 is obtained as follows 

```python
sage: [U,Uc]=GeneralUnitaryHypermatrix(4) 
sage: Prod(U, Uc.transpose(3), U.transpose(2), Uc.transpose()).simplify()
[[[[1, 0], [0, 0]], [[0, 0], [0, 0]]], [[[0, 0], [0, 0]], [[0, 0], [0, 1]]]]
```

Note that the function returns a symbolic parametrization of unitary hypermatrix and it's complex conjugate.
It is therefore important to emphasize that all the parameters in the parametrization are to be assigned only
real numbers.

## Uncorrelated tuples: Generalization of matrix inverse pair
The [general linear group](http://en.wikipedia.org/wiki/General_linear_group) plays a crucial role in many areas of
mathematics. The BM algebra suggests a natural generalization for the general linear group, for which very little is known.
Recall that a pair of square matrices are said to form an inverse pair if their product yields the Kronecker delta of the
same size and order. Similarly a tuple of cubic hypermatrices are said to form an uncorrelated tuple if their product
yield the Kronecker delta of the same order and size. The hypermatrix package implements a symbolic parametrization for
uncorrelated tuples of arbitrary order but having side length equal to 2 obtained as follows

```python
sage: [Ha,Hb]=GeneralUncorrelatedHypermatrixTuple(2); Prod(Ha,Hb).simplify() 
[[1, 0], [0, 1]]
sage: [Ha,Hb,Hc]=GeneralUncorrelatedHypermatrixTuple(3); Prod(Ha,Hb,Hc).simplify()
[[[1, 0], [0, 0]], [[0, 0], [0, 1]]]
sage: [Ha,Hb,Hc,Hd]=GeneralUncorrelatedHypermatrixTuple(4); Prod(Ha,Hb,Hc,Hd).simplify()
[[[[1, 0], [0, 0]], [[0, 0], [0, 0]]], [[[0, 0], [0, 0]], [[0, 0], [0, 1]]]]
```

## Multistochastic hypermatrices
The hypermatrix package also provide and symbolic parametrization of multistochastic hypermartices of arbitrary orders
and side length equal to 2 obtained as follows

```python
sage: GeneralStochasticHypermatrix(var('x'), 2).printHM()
[:, :]=
[cos(x)^2 sin(x)^2]
[sin(x)^2 cos(x)^2] 
sage: GeneralStochasticHypermatrix(var('x'), 3).printHM()
[:, :, 0]=
[cos(x)^2 sin(x)^2]
[sin(x)^2 cos(x)^2]

[:, :, 1]=
[sin(x)^2 cos(x)^2]
[cos(x)^2 sin(x)^2]

sage: GeneralStochasticHypermatrix(var('x'), 4)
[[[[cos(x)^2, sin(x)^2], [sin(x)^2, cos(x)^2]], [[sin(x)^2, cos(x)^2], [cos(x)^2, sin(x)^2]]], [[[sin(x)^2, cos(x)^2], [cos(x)^2, sin(x)^2]], [[cos(x)^2, sin(x)^2], [sin(x)^2, cos(x)^2]]]]
```


## Diagonal third order hypermatrices.
The package provides a routine for consrtructing third order generalization of diagonal matrices. Recall that using
sage diagonal matrices are constructed in sage using vectors as follows

```python
sage: Dg=diagonal_matrix(SR,HM(2,'d').list()); Dg
[d0  0]
[ 0 d1]

```

A defining property of diagonal matrices is the identity

```python
sage: (Dg.transpose()*Dg-Dg.elementwise_product(Dg)).is_zero()
True
```

Diagonal third order hypermatrices are constructed from symmetric matrices as follows

```python
sage: Mtr=Matrix(SR,SymMatrixGenerate(2,'d')); Mtr
[d00 d01]
[d01 d11]
sage: Dg=HM(Mtr); Dg.printHM()
[:, :, 0]=
[d00   0]
[d01   0]

[:, :, 1]=
[  0 d01]
[  0 d11]
```

The defining property of diagonal third order hypermatrices is the identity

```python
sage: (Prod(Dg.transpose(),Dg.transpose(2),Dg)-Dg.elementwise_product(Dg).elementwise_product(Dg)).is_zero()
True
```


## Constructing general rank one hypermatrices

Recall from the algebra of matrices that rank one matrices of size 2 by 3 are obtained as follows

```python
sage: Prod(HM(2,1,'a'),HM(1,3,'b')).printHM()
[:, :]=
[a00*b00 a00*b01 a00*b02]
[a10*b00 a10*b01 a10*b02]
```

Quite similarly a general rank one third order hypermatrix of size 2 by 3 by 4 is obtained as follows

```python
sage: Prod(HM(2,1,4,'a'),HM(2,3,1,'b'),HM(1,3,4,'c')).printHM()
[:, :, 0]=
[a000*b000*c000 a000*b010*c010 a000*b020*c020]
[a100*b100*c000 a100*b110*c010 a100*b120*c020]

[:, :, 1]=
[a001*b000*c001 a001*b010*c011 a001*b020*c021]
[a101*b100*c001 a101*b110*c011 a101*b120*c021]

[:, :, 2]=
[a002*b000*c002 a002*b010*c012 a002*b020*c022]
[a102*b100*c002 a102*b110*c012 a102*b120*c022]

[:, :, 3]=
[a003*b000*c003 a003*b010*c013 a003*b020*c023]
[a103*b100*c003 a103*b110*c013 a103*b120*c023]
```

It is possible to extract outer product summands from the any BM product by using some Kronecker delta type matrices as follows

```python
sage: sz=3;A=HM(sz,sz,'a'); B=HM(sz,sz,'b')
sage: od=2; Dlt0=HM(od,[1,0,0],'diag')
sage: ProdB(A,B,Dlt0).printHM()
[:, :]=
[a00*b00 a00*b01 a00*b02]
[a10*b00 a10*b01 a10*b02]
[a20*b00 a20*b01 a20*b02]
```

For third order hypermatrices we proceed as follows

```python
sage: sz=2;A=HM(sz,sz,sz,'a'); B=HM(sz,sz,sz,'b'); C=HM(sz,sz,sz,'c')
sage: od=3; Dlt0=HM(od,[1,0],'diag')
sage: ProdB(A,B,C,Dlt0).printHM()
[:, :, 0]=
[a000*b000*c000 a000*b010*c010]
[a100*b100*c000 a100*b110*c010]

[:, :, 1]=
[a001*b000*c001 a001*b010*c011]
[a101*b100*c001 a101*b110*c011]
```


## Revisiting Hyperdeterminants

The BM algebra of hypermatrices suggests a natural generalization to the classical matrix determinant.
determinant of side length two hypermatrices are computed as follows 

```python
sage: HM(2,2,'m').det()
-m01*m10 + m00*m11
sage: HM(2,2,2,'m').det()
m000*m011*m101*m110 - m001*m010*m100*m111
sage: HM(2,2,2,2,'m').det()
-m0001*m0010*m0100*m0111*m1000*m1011*m1101*m1110 + m0000*m0011*m0101*m0110*m1001*m1010*m1100*m1111
sage: HM(2,2,2,2,2,'m').det()
m00000*m00011*m00101*m00110*m01001*m01010*m01100*m01111*m10001*m10010*m10100*m10111*m11000*m11011*m11101*m11110 - m00001*m00010*m00100*m00111*m01000*m01011*m01101*m01110*m10000*m10011*m10101*m10110*m11001*m11010*m11100*m11111
```

# Transposition hypermatrices.
Transposition hypermatrices are analog of permutation matrices. We illustrate here their use
for transposing the first and second row slices of a 3x3x3 hypermatrix

```python
sage: A=HM(3,3,3,'a'); A.list()
[a000,a100,a200,a010,a110,a210,a020,a120,a220,a001,a101,a201,a011,a111,a211,a021,a121,a221,a002,a102,a202,a012,a112,a212,a022,a122,a222]
sage: P=HM(3,[1,0,2],'perm'); P
[[[0, 1, 0], [1, 0, 0], [0, 0, 1]], [[0, 1, 0], [1, 0, 0], [0, 0, 1]], [[0, 1, 0], [1, 0, 0], [0, 0, 1]]]
sage: Prod(P.transpose(),P.transpose(2),A).list()
[a100,a000,a200,a110,a010,a210,a120,a020,a220,a101,a001,a201,a111,a011,a211,a121,a021,a221,a102,a002,a202,a112,a012,a212,a122,a022,a222]
```

For transposing the first two column slices we use the following instructions

```python
sage: A=HM(3,3,3,'a'); A
sage: P=HM(3,[1,0,2],'perm'); P
sage: Prod(A,P,P.transpose())
```

For transposing the first two depth slices we use the following instructions

```python
sage: A=HM(3,3,3,'a'); A
sage: P=HM(3,[1,0,2],'perm'); P
sage: Prod(P,A,P.transpose(2))
```

Transposition can composed to perform an arbitrary permutation as illustrated below

```python
sage: A=HM(3,3,3,'a'); A
sage: P=HM(3,[1,0,2],'perm'); P
sage: Q=HM(3,[2,1,0],'perm'); Q
sage: Prod(Q, Prod(P,A,P.transpose(2)), Q.transpose(2))
```

# Slice linear combination

The package also allows us to extend familiar row/column operations to hypermatrices.
We describe how to perform row/column/depth slice operations on hypermartrices.
The general principle is best illustrated for hypermatrices of size 2 by 2 by 2.

```python
sage: A=HM(2,2,2,'a'); A
[[[a000, a001], [a010, a011]], [[a100, a101], [a110, a111]]]
```
The basic slice operation is the transposition of the two slices via transposition
hypermatrices.

```python
sage: Prod(HM(3,[1,0],'perm'),A,HM(3,[1,0],'perm').transpose(2))
[[[a001, a000], [a011, a010]], [[a101, a100], [a111, a110]]]
```

The last operation is the slice linear combination operation computed in two steps.
The first step amount to slice selection and scaling.

```python
sage: Idj=HM(3,range(2),'perm'); Idj[0,0,0]=0; Idj[1,0,0]=0 
sage: A+var('c')*Prod(HM([1,0],'perm'), Prod(Idj,A,Idj.transpose(2)), HM([1,0],'perm').transpose(2))
[[[a001*c + a000, a001], [a011*c + a010, a011]], [[a101*c + a100, a101], [a111*c + a110, a111]]]
```

# Action of hypermatrices.
Recall from linear algebra that the action of 2 by 3 matrix on a row vector of size
1 by 2 is defined by the following product

```python
sage: Prod(HM(1,2,'x'),HM(2,3,'a')).printHM()
[:,:]=
[a00*x00 + a10*x01 a01*x00 + a11*x01 a02*x00 + a12*x01]
```

Also recall that the action of a 2 by 3 matrix on a column vector of size 3 by 1
is prescribed by the product

```python
sage: Prod(HM(2,3,'a'),HM(3,1,'x')).printHM()
[a00*x00 + a01*x10 + a02*x20]
[a10*x00 + a11*x10 + a12*x20]
```

Quite similarly the a action of a hypermatrices of size 2,3,4 on a pair 
depth slice matrices of respective sizes 2 by 4 and 4 by 3 is prescribed
by the product resulting in a matrix slice of size 2 by 3

```python
sage: Prod(HM(2,4,1,'x'),HM(2,3,4,'a'),HM(4,3,1,'y'))[0,0,0]
a000*x000*y000 + a001*x010*y100 + a002*x020*y200 + a003*x030*y300
sage: Prod(HM(2,4,1,'x'),HM(2,3,4,'a'),HM(4,3,1,'y'))[1,0,0]
a100*x100*y000 + a101*x110*y100 + a102*x120*y200 + a103*x130*y300
```
By transposing the result we the product we obtain derivation of the action of a third order hypermatrix on other
pair of matrix slices.

# Triangulations and Tetrahedralizations

Recall that matrix algebra of can be used derive edge an edge list description of triangulations of convex polygon.
To list all the triangulation of a convex polygon on 5 vertices for instance we proceed as follows

```python
sage: TriangulationGraphs(4)
[[a01, a04, a12, a14, a23, a24, a34],
 [a01, a04, a12, a13, a14, a23, a34],
 [a01, a02, a04, a12, a23, a24, a34],
 [a01, a03, a04, a12, a13, a23, a34],
 [a01, a02, a03, a04, a12, a23, a34]]
```

If instead we wanted to sample uniformly from such triangulations we proceed by initializing the adjacency matrix
of the graph which embeds all the triangulations

```python
sage: sz=4
sage: A=HM(sz,sz,'a').elementwise_product(HM(sz,sz,'one')-HM(2,sz,'kronecker'))
sage: for i0 in range(1,sz):
    :   for i1 in range(i0):
    :       A[i0,i1]=0
```

Once the adjancency matrix A has been initialized, the sampling is performed as follows

```python
sage: Set(RandomTriangulation(A,A,sz-1,sz)[0].list()).difference(Set([0])).list()[0]
```

Similarly the BM algebra can be used to list tetrahedralization embeded in some hypergraph
as follows. First we initialize the following Hypergraph on 5 vertices
```python
sage: sz=5; S=HM(sz,sz,sz,'zero') 
sage: for i in range(sz): 
....:   for j in range(sz):
....:       for k in range(sz):       
....:           if i<j and j<k:
....:               S[i,j,k]=1; S[i,k,j]=1
sage: A=HM(sz,sz,sz,'a').elementwise_product(S)
```

The hypermatrix A denotes the symbolic adjacency hypermatrix of the hypergraph. The
tetrahedralization is obtained as follows

 ```python
sage: L=Tetrahedralizations(A,A,sz,sz)
sage: [Set(f.list()).difference(Set([0])).list() for f in L]
[[a013*a041*a043*a123*a142*a143*a243, a014*a031*a034*a124*a132*a134*a234],
 [a012*a023*a041*a042*a043*a142*a243, a012*a024*a031*a032*a034*a132*a234],
 [a014*a021*a024*a032*a034*a124*a234, a013*a021*a023*a042*a043*a123*a243]]
```

# Misc

We collect here description of functions in the package that the authors found handy.
In many situations one is presented with a set of linear constraints presented as
list of symbolic expressions and one wants to extract a matrix and right hand side.
This can be done as follows
 
 ```python
sage: sz=2; od=2; X=HM(sz,sz,var_list('x',sz^2)); f=X.det()
sage: H=Form2TotallySymmetricHypermatrix(f, 2, X.list()); H.printHM()
[:, :]=
[   0    0    0  1/2]
[   0    0 -1/2    0]
[   0 -1/2    0    0]
[ 1/2    0    0    0]
sage: sz=2; od=2; A=HM(sz,sz,'a'); X=HM(od,'x').list()
sage: Hv0=HM(sz,1,[1,X[0]]); Hv1=HM(sz,1,[1,X[1]])
sage: f=ProdB(Hv0.transpose(),Hv1,A)[0,0]; f
a11*x0*x1 + a10*x0 + a01*x1 + a00
sage: Form2Hypermatrix(f, X).printHM()
[:, :]=
[a00 a01]
[a10 a11]
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
sage: sz=2; EqL=[var('x'+str(i))+var('x'+str(sz+j))==var('a'+str(i)+str(j)) for i in range(sz) for j in range(sz)]
sage: default_linear_solver(EqL, var_list('x', 2*sz), var_list('t', 2*sz))
[x0 == a00 - a10 + a11 - t3,
 x1 == a11 - t3,
 x2 == a10 - a11 + t3,
 0  == -a00 + a01 + a10 - a11]
```

To compute resultant via Sylvester matrix construction we proceed as follows

 ```python
sage: x, a0, a1, b0, b1=var('x, a0, a1, b0, b1')
sage: p=expand((x-a0)*(x-a1))
sage: q=expand((x-b0)*(x-b1))
sage: Sylvester_matrix(p, q, x).det().factor()
(a0 - b0)*(a0 - b1)*(a1 - b0)*(a1 - b1)
```

To substitute a matrix in polynomial expression we proceed as follows

 ```python
sage: x,y = var('x,y')
sage: p=x^2+2*x*y+1
sage: substitute_matrix(p,x,Matrix(SR,HM(2,2,'a').listHM()))
[a00^2 + a01*a10 + 2*a00*y + 1   a00*a01 + a01*a11 + 2*a01*y]
[  a00*a10 + a10*a11 + 2*a10*y a01*a10 + a11^2 + 2*a11*y + 1]
```

To loading a PNG image into a hypermatrix ( warning: very slow !)

```python
sage: A=HM("pic.png")
```

The slicer is given by 

```python
sage: sz=3; A=HM(sz,sz,sz,'a')
sage: ThirdOrderSlicer(A, [0], 'row').printHM()
[:, :, 0]=
[a000 a010 a020]
<BLANKLINE>
[:, :, 1]=
[a001 a011 a021]
<BLANKLINE>
[:, :, 2]=
[a002 a012 a022]
<BLANKLINE>
sage: ThirdOrderSlicer(A, [1], 'col').printHM()
[:, :, 0]=
[a010]
[a110]
[a210]
<BLANKLINE>
[:, :, 1]=
[a011]
[a111]
[a211]
<BLANKLINE>
[:, :, 2]=
[a012]
[a112]
[a212]
<BLANKLINE>
sage: ThirdOrderSlicer(A, [2], 'dpt').printHM()
[:, :, 0]=
[a002 a012 a022]
[a102 a112 a122]
[a202 a212 a222] 
<BLANKLINE> 
```

Code for Discrete Fourier matrices

```python
sage: [Ha, Hb]=SecondOrderDFT(2,2)
sage: Ha
[[1, 1], [1, -1]]
sage: Hb
[[1, 1], [1, -1]]
sage: Prod(Ha,Hb).simplify()
[[2, 0], [0, 2]]
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
```



# Bug report

Please reporting bugs is greatly appreciated.

