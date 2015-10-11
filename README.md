=======
# Hypermatrix Algebra Package

We provide here a sagemath implementation of the Bhattacharya-Mesner(BM) algebra as well as the general BM algebra and its
dual.

The `Hypermatrix Algebra Package` is a symbolic hypermatrix package designed to experimentally investigate structural and
combinatorial properties of the algebra of hypermatrices.

# Installation 

A properly working install of [sage](http://sagemath.org/) is the only prerequisite to using the 
hypermatrix package. The hypermatrix algebra package has been tested on SageMath Version 6.8.
To get started with SageMath, the authors of the hypermatrix algebra package highly recommend reading 
[Calcul mathématique avec Sage] (http://sagebook.gforge.inria.fr/)

To use the hypermatrix algebra package, simply download the [hypermatrix algebra package sage file](https://github.com/gnang/HypermatrixAlgebraPackage/blob/master/Hypermatrix_Algebra_Package_code.sage) into your working directory and load the file into your SageMath [interactive shell](http://doc.sagemath.org/html/en/tutorial/interactive_shell.html) session using the command:

```python
sage: %runfile("Hypermatrix_Algebra_Package_code.sage")
```

# Usage

To create a symbolic hypermatrix instance of size say 2x3x4, use the instructions

```python
sage: Ha=HM(2,3,4,'a')
sage: Ha 
[[[a000,a001,a002,a003],[a010,a011,a012,a013],[a020,a021,a022,a023]],[[a100,a101,a102,a103],[a110,a111,a112,a113],[a120,a121,a122,a123]]]
```

Alternatively, hypermatrices can be initialized from an arbitrary Python list as follows

```python
sage: rng=range(1,3)
sage: Lst=[var('x{0}{1}{2}'.format(i,j,k)) for k in rng for j in rng for i in rng]
sage: Hx=HM(2,2,2,Lst)
sage: Hx
[[[x111, x112], [x121, x122]], [[x211, x212], [x221, x222]]]
```

Note that hypermatrix entries are not restricted to numbers and symbolic expressions. The hypermatrix 
entries may be of other sage object types or mixture of different sage object types including list, Matrix,
Boolean or even hypermatrices as illustrated below

```python
sage: rng = range(2)
sage: Hh=HM(2,2,2,[(i*2^2+j*2+k)*hadamard_matrix(2) for k in rng for j in rng for i in rng])
```

Hypermatrix entries are accessed quite similarly to the way in which matrix entries are usually
accessed in sage. For example, the (0,1,0) entry of the hypermatrix constructed in the previous 
instruction is accessed as follows

```python
sage: Hh[0,1,0]
[ 2  2]
[ 2 -2]
```

# Basic hypermatrix operation

The [BM product](http://arxiv.org/abs/1411.6270) is implemented for hypermatrices of all orders 
having compatible sizes. The BM hypermatrix product is computed as follows

```python
sage: H1=Prod(HM(2,3,'a'), HM(3,4,'c'))
sage: H2=Prod(HM(2,3,4,'a'), HM(2,2,3,'b'), HM(3,2,4,'c'))
sage: H3=Prod(HM(2,1,2,2,'a'), HM(2,2,1,2,'b'), HM(2,2,2,1,'c'), HM(1,2,2,2,'d'))
```

As illustrated by the previous instructions, the product of second order hypermatrices with compatible
sizes corresponds to the usual matrix product. Other basic hypermatrix operations including addition and 
multiplication by scalars are analogous to their matrix counterpart in sage. The hypermatrix transpose 
generalizes the matrix transpose and amounts to a cyclic permutation of the entry indices and is performed
as follows 

```python
sage: HM(2,2,2,'a')
[[[a000, a001], [a010, a011]], [[a100, a101], [a110, a111]]]
sage: HM(2,2,2,'a').transpose()
[[[a000, a100], [a001, a101]], [[a010, a110], [a011, a111]]]
```

To perform two or more consecutive transposes, we pass on to the transpose function an integer argument as 
follows

```python
sage: At=HM(2,2,2,2,'a').transpose( )
sage: Att=HM(2,2,2,2,'a').transpose(2)
sage: Attt=HM(2,2,2,2,'a').transpose(3)
sage: (HM(2,2,2,2,'a').transpose(4)-HM(2,2,2,2,'a')).is_zero()
True
```

The sum of hypermatrices can taken over a list of hypermatrices.
We illustrate this by constructing a symbolic symmetric hypermatrix as follows.

```python
sage: A=sum(HM(2,2,2,'a').transpose(i) for i in range(3))
sage: A.is_symmetric()
True
```

Many of the properties of the special hypermatrices that we describe subsequently are preserved by
the Kroencker product. The hypermatrix package provides implementations of the Kronecker product 
for hypermatrices of orders less or equal to 5 (later versions of the package will address this limitation).
The Kronecker product of two hypermatrices are obtained
as follows

```python
sage: Ha=HM(2,2,2,'a'); Hb=HM(3,3,3,'b')
sage: Hc=Ha.tensor_product(Hb)
sage: Hc.dimensions()
[6, 6, 6]
```

An additional basic hypermatrix operations implemented in the package are the inner-product and 
multilinear forms induced by hypermatrices of various order. We recall from matrix algebra that
the innerproduct of a pairs column vectors X and Y each of size 3 by 1, is defined by the product

```python
sage: X=HM(3, 1, HM(3,'x').list()); Y=HM(3, 1, HM(3,'y').list())
sage: Prod(X.transpose(),Y)[0,0]
x0*y0 + x1*y1 + x2*y2
sage: apply(Prod, [X.transpose(i) for i in range(X.order()-1,-1,-1)])[0,0]
```

Similarly for third order hypermatrices the innerproduct of a triplet of column vectors X, Y, Z
each of size 3 by 1 by 1, is defined by the product

```python
sage: X=HM(3,1,1,HM(3,'x').list()); Y=HM(3,1,1,HM(3,'y').list()); Z=HM(3,1,1,HM(3,'z').list())
sage: Prod(X.transpose(2), Y.transpose(), Z)[0,0,0]
x0*y0*z0 + x1*y1*z1 + x2*y2*z2
sage: apply(Prod, [X.transpose(i) for i in range(X.order()-1,-1,-1)])[0,0,0]
```

and for fourth order hypermatrices the innerproduct of a quadruplet of column vectors W, X, Y, Z
each of size 3 by 1 by 1 by 1, is defined by the product

```python
sage: W=HM(3,1,1,1,HM(3,'w').list())
sage: X=HM(3,1,1,1,HM(3,'x').list())
sage: Y=HM(3,1,1,1,HM(3,'y').list())
sage: Z=HM(3,1,1,1,HM(3,'z').list())
sage: Prod(W.transpose(3), X.transpose(2), Y.transpose(), Z)[0,0,0,0]
w0*x0*y0*z0 + w1*x1*y1*z1 + w2*x2*y2*z2
sage: apply(Prod, [X.transpose(i) for i in range(X.order()-1,-1,-1)])[0,0,0,0]
```

Furthermore the bilinear form associated with a 2 by 2 matrix A is expressed as follows

```python
sage: X=HM(2,1,HM(2,'x').list()); Y=HM(2,1,HM(2,'y').list()); A=HM(2,2,'a')
sage: ProdB(X.transpose(), Y, A)[0,0]
a00*x0*y0 + a10*x1*y0 + a01*x0*y1 + a11*x1*y1
```

Finally a trilinear form associated with a 2 by 2 by 2 hypermatrix A is expressed as follows

```python
sage: X=HM(2,1,1,HM(2,'x').list()); Y=HM(2,1,1,HM(2,'y').list()); Z=HM(2,1,1,HM(2,'z').list())
sage: A=HM(2,2,2,'a')
sage: ProdB(X.transpose(2), Y.transpose(), Z, A)[0,0,0]
```


# Symbolic parametrization of some special hypermatrices

To emphasize that hypermatrices naturally extend the algebra of matrices, we describe instructions
for obtaining symbolic parametrization of special families of hypermatrices which generalizes well
known symbolic parametrization of special families of matrices .

## Kronecker delta hypermatrices
The [Kronecker delta](http://en.wikipedia.org/wiki/Kronecker_delta#Properties_of_generalized_Kronecker_delta)
hypermatrices are generalization of the identity matrices in the sense that all of it's entries are zero 
except for the entries on the main diagonal as illustated below 
```python
sage: rng=range(2)
sage: Dlt=HM(3,2,'kronecker')
sage: A=HM(2,2,2,'a')
sage: HM(2,2,2,[A[i,j,k]==Dlt[i,j,k] for k in rng for j in rng for i in rng])
[[[a000==1, a001==0], [a010==0, a011==0]], [[a100==0, a101==0], [a110==0, a111==1]]]
```

## Orthogonal hypermatrices

Orthogonal hypermatrices are analogous to orthogonal matrices in the fact that the product of the transposes 
yields the Kronecker delta hypermatrix. A symbolic parametetrization of 2x2x2 orthogonal hypermatrices is 
obtained via the following instructions

```python
sage: Q=GeneralOrthogonalHypermatrix(3)
sage: (apply(Prod,[Q.transpose(i) for i in range(Q.order(),0,-1)]).simplify()-HM(Q.order(),2,'kronecker')).is_zero()
True
```

Similarly a symbolic parametrization of orthogonal hypermatrices of size 2x2x2x2 hypermatrices is obtained
obtained by the command

```python
sage: Q=GeneralOrthogonalHypermatrix(4)
sage: (apply(Prod,[Q.transpose(i) for i in range(Q.order(),0,-1)]).simplify()-HM(Q.order(),2,'kronecker')).is_zero()
True
```

The parameters which appear in the parametrization can be assigned arbitrary complex values. The hypermatrix 
package also provides construction of cubic third order Hadamard hypermatrices, which are analogous to 
Hadamard matrices in the fact that the product of the transposes is equal to the Kronecker delta scaled by
the side length.
It is easy to see that the classical Hadamard matrix [construction](http://en.wikipedia.org/wiki/Hadamard_matrix#Sylvester.27s_construction) first proposed by James Joseph Sylvester, and going back to 1867, extends to all
hypermatrices of prime order (recall that the order refers to the number of indices associated with the 
hypermatrix entries). Third order Hadamard hypermatrices whose side length are powers of 2 and in the 
particular case of size 4x4x4 is obtained via the following instructions

```python
sage: H=ThirdOrderHadamardBlockU(4); H
[[[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]], [[-1, 1, -1, 1], [1, 1, 1, 1], [-1, 1, -1, 1], [1, 1, 1, 1]], [[-1, -1, 1, 1], [-1, -1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]], [[1, -1, -1, 1], [-1, -1, 1, 1], [-1, 1, -1, 1], [1, 1, 1, 1]]]
sage: (apply(Prod,[H.transpose(i) for i in range(H.order(),0,-1)])-4*HM(H.order(),4,'kronecker')).is_zero()
True
```

One can not help but wonder if the reach of Hadamard matrix conjecture extendeds to hypermatrices of prime order.

## Unitary hypermatrices

The hypermatrix algebra package also implements a symbolic parametrization of even order unitary hypermatrices having 
of sides length equal to 2. As suggested by the [Gowers norm](http://en.wikipedia.org/wiki/Gowers_norm), an even order
hypermatrix is said to be unitary if the product of the conjugate transposes yield equals the Kronecker delta.
Consequently, second order unitary hypermatrices correspond to unitary matrices and a symbolic parametrization
of 2x2 unitary matrices is obtained from the hypermatrix package via the commands

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
real number values in order for the hypermatrix unitary parametrization to be valid.

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
sage: GeneralStochasticHypermatrix(var('x'), 2)
[[cos(x)^2, sin(x)^2], [sin(x)^2, cos(x)^2]] 
sage: GeneralStochasticHypermatrix(var('x'), 3)
[[cos(x)^2, sin(x)^2], [sin(x)^2, cos(x)^2]], [[sin(x)^2, cos(x)^2], [cos(x)^2, sin(x)^2]]]
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
and diagonal matrices are characterized by the identity

```python
sage: (Dg.transpose()*Dg-Dg.elementwise_product(Dg)).is_zero()
True
```
Quite similarly, diagonal third order hypermatrices are constructed from symmetric matrices.
```python
sage: Mtr=Matrix(SR,SymMatrixGenerate(2,'d')); Mtr
[d00 d01]
[d01 d11]
sage: Dg=HM(Mtr); Dg
[[[d00, 0], [0, d01]], [[d01, 0], [0, d11]]]
```

Futhermore, the defining property of diagonal third order hypermatrices is the identity

```python
sage: (Prod(Dg.transpose(),Dg.transpose(2),Dg)-Dg.elementwise_product(Dg).elementwise_product(Dg)).is_zero()
True
```


## Constructing general rank one hypermatrices
Recall from the algebra of matrices that rank one matrices of size 2 by 3 are obtained as follows
```python
sage: Matrix(SR,Prod(HM(2,1,'a'),HM(1,3,'b')).listHM())
[a00*b00 a00*b01 a00*b02]
[a10*b00 a10*b01 a10*b02]
```
Quite similarly a general rank one third order hypermatrix of size 2 by 3 by 4 is obtained as follows
```python
sage: Prod(HM(2,1,4,'a'),HM(2,3,1,'b'),HM(1,3,4,'c'))
[[[a000*b000*c000, a001*b000*c001, a002*b000*c002, a003*b000*c003], [a000*b010*c010, a001*b010*c011, a002*b010*c012, a003*b010*c013], [a000*b020*c020, a001*b020*c021, a002*b020*c022, a003*b020*c023]], [[a100*b100*c000, a101*b100*c001, a102*b100*c002, a103*b100*c003], [a100*b110*c010, a101*b110*c011, a102*b110*c012, a103*b110*c013], [a100*b120*c020, a101*b120*c021, a102*b120*c022, a103*b120*c023]]]
```


## Revisiting Hyperdeterminants

The BM algebra of hypermatrices suggests a natural generalization to the classical matrix determinant which is also 
implemented in the hypermatrix package. The user must be warn that while it is relatively straight forward to compute
determinant of hypermatrices of arbitrary order of with side length equal to 2 via the following expressions 

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
however computing the determinant of hypermatrices of arbitrary order with side length greater then 3 is considerably 
more difficult.


# Transposition hypermatrices.
Transposition hypermatrices implement the analog of permutation matrices. The package implements the hypermatrix 
analog of permutation matrices. We illustrate here their use for transposing the first and second row slices of a
3x3x3 hypermatrix

```python
sage: A=HM(3,3,3,'a'); A.list()
[a000,a100,a200,a010,a110,a210,a020,a120,a220,a001,a101,a201,a011,a111,a211,a021,a121,a221,a002,a102,a202,a012,a112,a212,a022,a122,a222]
sage: P=HM([1,0,2],'perm'); P
[[[0, 1, 0], [1, 0, 0], [0, 0, 1]], [[0, 1, 0], [1, 0, 0], [0, 0, 1]], [[0, 1, 0], [1, 0, 0], [0, 0, 1]]]
sage: Prod(P.transpose(),P.transpose(2),A).list()
[a100,a000,a200,a110,a010,a210,a120,a020,a220,a101,a001,a201,a111,a011,a211,a121,a021,a221,a102,a002,a202,a112,a012,a212,a122,a022,a222]
```

For transposing the first two column slices we use the following instructions

```python
sage: A=HM(3,3,3,'a'); A
sage: P=HM([1,0,2],'perm'); P
sage: Prod(A,P,P.transpose())
```

For transposing the first two depth slices we use the following instructions

```python
sage: A=HM(3,3,3,'a'); A
sage: P=HM([1,0,2],'perm'); P
sage: Prod(P,A,P.transpose(2))
```

Transposition can composed to perform an arbitrary permutation as illustrated below

```python
sage: A=HM(3,3,3,'a'); A
sage: P=HM([1,0,2],'perm'); P
sage: Q=HM([2,1,0],'perm'); Q
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
sage: Prod(HM([1,0],'perm'),A,HM([1,0],'perm').transpose(2))
[[[a001, a000], [a011, a010]], [[a101, a100], [a111, a110]]]
```
The last operation is the slice linear combination operation computed in two steps.
The first step amount to slice selection and scaling.
```python
sage: Idj=HM(range(2),'perm'); Idj[0,0,0]=0; Idj[1,0,0]=0 
sage: A+var('c')*Prod(HM([1,0],'perm'), Prod(Idj,A,Idj.transpose(2)), HM([1,0],'perm').transpose(2))
[[[a001*c + a000, a001], [a011*c + a010, a011]], [[a101*c + a100, a101], [a111*c + a110, a111]]]
```

# Action of hypermatrices.
We recall from linear algebra that the action of 2 by 3 matrix on a row vector of size
1 by 2 is defined by the following product
```python
sage: Matrix(SR,Prod(HM(1,2,'x'),HM(2,3,'a')).listHM())
[a00*x00 + a10*x01 a01*x00 + a11*x01 a02*x00 + a12*x01]
```
We also recall that the action of a 2 by 3 matrix on a column vector of size 3 by 1
is prescribed by the product
```python
sage: Matrix(SR,Prod(HM(2,3,'a'),HM(3,1,'x')).listHM())
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
The algebra of matrices can be used to produce edge lists for triangulations of convex polygon
To list all the triangulation of a convex polygon on 5 vertices for instance we proceed as follows
```python
sage: TriangulationGraphs(4)
[[a01, a04, a12, a14, a23, a24, a34],
 [a01, a04, a12, a13, a14, a23, a34],
 [a01, a02, a04, a12, a23, a24, a34],
 [a01, a03, a04, a12, a13, a23, a34],
 [a01, a02, a03, a04, a12, a23, a34]]
```
For sampling such triangulations uniformly at random we proceed by initializing the adjacency matrix
of the graph which embeds all the triangulations

```python
sage: sz=4
sage: A=HM(sz,sz,'a').elementwise_product(HM(sz,sz,'one')-HM(2,sz,'kronecker'))
sage: for i0 in range(1,sz):
    :   for i1 in range(i0):
    :       A[i0,i1]=0
```
Once the adjancency matrix A has been initialized, we sample as follows
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

# Bug report

Please report any bugs, it is greatly appreciated.
