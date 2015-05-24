=======
# Hypermatrix Algebra Package

We provide here a sagemath implementation of the Bhattacharya-Mesner(BM) algebra and it's
generalizations.

**UPDATE 2015-02-15** Major changes to the codebase, as reflected by the current post.

The `Hypermatrix Algebra Package` is a symbolic hypermatrix package designed to
investigate structural and combinatorial properties of hypermatrices.

# Installation 

A properly working install of [sage](http://sagemath.org/) is a prerequisite to using the hypermatrix 
package. Download the Hypermatrix sage file into your working directory. Load the hypermatrix package 
into a sage terminal session using the following command:

```python
sage: %runfile("Hypermatrix_Algebra_Package_code.sage")
```

# Usage

To create a symbolic hypermatrix instance of size 2x3x4 for example, we use the instructions

```python
sage: Ha=HM(2,3,4,'a')
sage: Ha 
[[[a000,a001,a002,a003],[a010,a011,a012,a013],[a020,a021,a022,a023]],[[a100,a101,a102,a103],[a110,a111,a112,a113],[a120,a121,a122,a123]]]
```

Alternatively hypermatrices can be initialized from a list as follows

```python
sage: rng=range(1,3)
sage: Lst=[var('x{0}{1}{2}'.format(i,j,k)) for k in rng for j in rng for i in rng]
sage: Hx=HM(2,2,2,Lst)
sage: Hx
[[[x111, x112], [x121, x122]], [[x211, x212], [x221, x222]]]
```

The hypermatrix entries are not restricted to numbers and symbolic expression, they may be matrices or more generally
can even be other hypermatrices as illustrated below

```python
sage: rng = range(2)
sage: Hh=HM(2,2,2,[(i*2^2+j*2+k)*hadamard_matrix(2) for k in rng for j in rng for i in rng])
```

Hypermatrix entries are accessed quite similarly to the way in which matrix entries are usually accessed in sage.
For example, the (0,1,0) entry of the hypermatrix constructed in the previous instruction is accessed as follows

```python
sage: Hh[0,1,0]
[ 2  2]
[ 2 -2]
```

# Basic hypermatrix operation

The [BM product](http://arxiv.org/abs/1411.6270) is implemented for hypermatrices of all orders and all 
compatible sizes. The hypermatrix product is performed as follows

```python
sage: Prod(HM(2,3,'a'), HM(3,4,'c'))
sage: Prod(HM(2,3,4,'a'), HM(2,2,3,'b'), HM(3,2,4,'c'))
```

As illustrated by the previous instructions, the product of compatible second order hypermatrices recovers the 
usual matrix product. Other basic hypermatrix operations including addition and multiplication by scalars are
quite similar to their matrix counter part in sage. The hypermatrix transpose amounts to a cyclic permutation of the
entry indices and is performed as follows 

```python
sage: HM(2,2,2,'a')
[[[a000, a001], [a010, a011]], [[a100, a101], [a110, a111]]]
sage: HM(2,2,2,'a').transpose()
[[[a000, a100], [a001, a101]], [[a010, a110], [a011, a111]]]
```

In order to perform two or more consecutive transposes we use the following instructions

```python
sage: HM(2,2,2,2,'a').transpose( )
sage: HM(2,2,2,2,'a').transpose(2)
sage: HM(2,2,2,2,'a').transpose(3)
sage: HM(2,2,2,2,'a').transpose(4)-HM(2,2,2,2,'a')
```

The sum of hypermatrices can taken over a list of hypermatrices.
We illustrate this by constructing a symbolic symmetric hypermatrix as follows.

```python
sage: sum([HM(2,2,2,'a').transpose(i) for i in range(3)])
```

Many of the properties of special hypermatrices that we describe subsequently are preserved by the Kroencker product.
The hypermatrix package provides an implementation of the Kronecker product for hypermatrices of orders going up to 5
(later versions of the package will address this limitation). The Kronecker product of two hypermatrices are obtained
as follows

```python
sage: A=HM(2,2,2,'a'); B=HM(3,3,3,'b')
sage: A.slicekroneckerproduct(B)
```

An additional basic hypermatrix operations implemented in the package are the inner-product and 
multilinear forms induced by hypermatrices of various order. We recall from matrix algebra that
the innerproduct of a pairs column vectors X and Y each of size 3 by 1, is defined by the product

```python
sage: X=HM(3,1,'x'); Y=HM(3,1,'y')
sage: Prod(X.transpose(),Y)[0,0]
```

Similarly for third order hypermatrices the innerproduct of a triplet of column vectors X, Y, Z
each of side 3 by 1 by 1, is defined by the product

```python
sage: X=HM(3,1,1,'x'); Y=HM(3,1,1,'y'); Z=HM(3,1,1,'z')
sage: Prod(X.transpose(2), Y.transpose(), Z)[0,0,0]
```

Furthermore the bilinear form associated with a 2 by 2 matrix A is expressed as follows

```python
sage: X=HM(2,1,'x'); Y=HM(2,1,'y'); A=HM(2,2,'a')
sage: ProdB(X.transpose(), Y, A)[0,0]
a00*x00*y00+a10*x10*y00+a01*x00*y10+a11*x10*y10
```

Finally a trilinear form associated with a 2 by 2 by 2 hypermatrix A is expressed as follows

```python
sage: X=HM(2,1,1,'x'); Y=HM(2,1,1,'y'); Z=HM(2,1,1,'z'); A=HM(2,2,2,'a')
sage: ProdB(X.transpose(2), Y.transpose(), Z, A)[0,0,0]
a000*x000*y000*z000+a100*x100*y000*z000+a010*x000*y100*z000+a110*x100*y100*z000+a001*x000*y000*z100+a101*x100*y000*z100+a011*x000*y100*z100+a111*x100*y100*z100
```


# Symbolic parametrization of some special hypermatrices
To emphasize that hypermatrices naturally extend the algebra of matrices, we describe here instructions
for obtaining some symbolic parametrization of special families of hypermatrices.

## Kronecker delta hypermatrices
The [Kronecker delta](http://en.wikipedia.org/wiki/Kronecker_delta#Properties_of_generalized_Kronecker_delta) hypermatrices
are generalization of the identity matrices in the sense that all of it's entries are zero except for the entries on the
main diagonal as illustated below 
```python
sage: Dlt=GeneralHypermatrixKroneckerDelta(3, 2)
sage: A=HM(2,2,2,'a')
sage: rng=range(2)
sage: HM(2,2,2,[A[i,j,k]==Dlt[i,j,k] for k in rng for j in rng for i in rng])
[[[a000==1, a001==0], [a010==0, a011==0]], [[a100==0, a101==0], [a110==0, a111==1]]]
```

## Orthogonal hypermatrices

Orthogonal hypermatrices are analogous to Orthogonal matrices in the fact that the product of the transposes equals
the Kronecker delta hypermatrix of the same size and order. A symbolic parametetrization of 2x2x2 orthogonal hypermatrices is 
obtained via the following instructions

```python
sage: Q=GeneralOrthogonalHypermatrix(3); Q
[[e^(-r1+r3+r6)/(e^(-3*r1+3*r3+3*r6)+e^(3*r6))^(1/3), e^r4], [e^r6/(e^(-3*r1+3*r3+3*r6) + e^(3*r6))^(1/3), e^r2]], [[-e^(r1+r2-r3-r4+r5), e^r3/(e^(3*r1)+e^(3*r3))^(1/3)], [e^r5, e^r1/(e^(3*r1)+e^(3*r3))^(1/3)]]]
sage: Prod(Q,Q.transpose(2),Q.transpose()).simplify()
[[[1, 0], [0, 0]], [[0, 0], [0, 1]]]
```

Similarly a symbolic parametrization of orthogonal parametrization of size 2x2x2x2 hypermatrices is obtained via
the following instructions

```python
sage: Q=GeneralOrthogonalHypermatrix(4); Q
sage: Prod(Q,Q.transpose(3),Q.transpose(2),Q.transpose()).simplify()
```

The parameters which appear in the parametrization can be assigned arbitrary complex values The hypermatrix package also
provides construction of cubic third order Hadamard hypermatrices, which are quite analogous
to Hadamard matrices in the fact that the product of the transposes equals to the Kronecker delta (of the same size and 
order) scaled by the side length of the hypermatrix. It is easy to see that the classical Hadamard matrix 
[construction](http://en.wikipedia.org/wiki/Hadamard_matrix#Sylvester.27s_construction) proposed by 
James Joseph Sylvester, and going back to 1867, extends to all hypermatrices of prime order (recall that the order refers
to the number of indices associated with the hypermatrix entries, for instance matrices are second order hypermatrices 
because every matrix entry has a row and column index). Third order Hadamard hypermatrices whose side length are powers 
of 2 and in the particular case of size 4x4x4 is obtained via the following instructions

```python
sage: H=ThirdOrderHadamardBlockU(4); H
[[[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]], [[-1, 1, -1, 1], [1, 1, 1, 1], [-1, 1, -1, 1], [1, 1, 1, 1]], [[-1, -1, 1, 1], [-1, -1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]], [[1, -1, -1, 1], [-1, -1, 1, 1], [-1, 1, -1, 1], [1, 1, 1, 1]]]
sage: Prod(H, H.transpose(2), H.transpose())
[[[4, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 4, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 4, 0], [0, 0, 0, 0]], [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 4]]]
```

One can not help but wonder if the reach of the famous Hadamard matrix conjecture extendeds to hypermatrices of prime order.

## Unitary hypermatrices
The hypermatrix package also provides a symbolic parametrization of even order unitary hypermatrices with sides length 
equal to 2. As suggested by the [Gowers norm](http://en.wikipedia.org/wiki/Gowers_norm), an even order hypermatrix is
said to be unitary if the product of the conjugate transposes yield the Kronecker delta of the same order and size.
Consequently, second order unitary hypermatrices correspond to the usual unitary matrices and a symbolic parametrization
of 2x2 unitary matrices is obtained from the hypermatrix package via the following commands

```python
sage: [U,Uc]=GeneralUnitaryHypermatrix(2); U
[[e^(-I*pi-I*r1+I*r2+I*r3-r4+r5+r6)/sqrt(e^(-2*r4+2*r5+2*r6)+e^(2*r6)), e^(I*r2+r6)/sqrt(e^(-2*r4+2*r5+2*r6)+e^(2*r6))], [e^(I*r3+r4)/sqrt(e^(2*r4)+e^(2*r5)), e^(I*r1+r5)/sqrt(e^(2*r4)+e^(2*r5))]] 
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
It is therefore important to emphasize that all the parameters in the parametrization are to be assigned real
values for the hypermatrix to be unitary.

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
sage diagonal matrices are constructed using vectors as follows
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
Transposition hypermatrices implement the analog of permutation matrices. The package implements the hypermatrix  analog of
permutation matrices. We illustrate here their use for transposing the first and second row slices
of a 3x3x3 hypermatrix

```python
sage: A=HM(3,3,3,'a'); A
[[[a000, a001, a002], [a010, a011, a012], [a020, a021, a022]], [[a100, a101, a102], [a110, a111, a112], [a120, a121, a122]], [[a200, a201, a202], [a210, a211, a212], [a220, a221, a222]]]
sage: P=HM([1,0,2],'perm'); P
[[[0, 1, 0], [1, 0, 0], [0, 0, 1]], [[0, 1, 0], [1, 0, 0], [0, 0, 1]], [[0, 1, 0], [1, 0, 0], [0, 0, 1]]]
sage: Prod(P.transpose(),P.transpose(2),A)
[[[a100, a101, a102], [a110, a111, a112], [a120, a121, a122]], [[a000, a001, a002], [a010, a011, a012], [a020, a021, a022]], [[a200, a201, a202], [a210, a211, a212], [a220, a221, a222]]]
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
sage: Prod(HM(2,4,1,'x'), HM(2,3,4,'a'), HM(4,3,1,'y'))
[[[a000*x000*y000 + a001*x010*y100 + a002*x020*y200 + a003*x030*y300], [a010*x000*y010 + a011*x010*y110 + a012*x020*y210 + a013*x030*y310], [a020*x000*y020 + a021*x010*y120 + a022*x020*y220 + a023*x030*y320]], [[a100*x100*y000 + a101*x110*y100 + a102*x120*y200 + a103*x130*y300], [a110*x100*y010 + a111*x110*y110 + a112*x120*y210 + a113*x130*y310], [a120*x100*y020 + a121*x110*y120 + a122*x120*y220 + a123*x130*y320]]]
```
By transposing the result we the product we obtain derivation of the action of a third order hypermatrix on other
pair of matrix slices.

# Bug report

Please report any bugs, it is greatly appreciated.
