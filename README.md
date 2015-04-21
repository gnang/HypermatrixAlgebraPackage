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
sage: Ha = HM(2,3,4,'a')
sage: Ha 
[[[a000, a001, a002, a003], [a010, a011, a012, a013], [a020, a021, a022, a023]], [[a100, a101, a102, a103], [a110, a111, a112, a113], [a120, a121, a122, a123]]]
```

Alternatively hypermatrices can be initialized from a list as follows

```python
sage: rng = range(1,3)
sage: Lst = [var('x{0}{1}{2}'.format(i,j,k)) for k in rng for j in rng for i in rng]
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
sage: Prod(HM(2,3,'a'),HM(3,4,'c'))
sage: Prod(HM(2,3,4,'a'),HM(2,2,3,'b'),HM(3,2,4,'c'))
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
sage: HM(2,2,2,2,'a').transpose()
sage: HM(2,2,2,2,'a').transpose(2)
sage: HM(2,2,2,2,'a').transpose(3)
sage: HM(2,2,2,2,'a').transpose(4) - HM(2,2,2,2,'a')
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


# Symbolic parametrization of some special hypermatrices

To emphasize that hypermatrices naturally extend the algebra of matrices, we describe here instructions for obtaining
some symbolic parametrization of special families of hypermatrices.

## Kronecker delta hypermatrices
The [Kronecker delta](http://en.wikipedia.org/wiki/Kronecker_delta#Properties_of_generalized_Kronecker_delta) hypermatrices
are generalization of the identity matrices in the sense that all of it's entries are zero except for the entries on the
main diagonal as illustated below 
```python
sage: Dlt=GeneralHypermatrixKroneckerDelta(3, 2)
sage: A=HM(2,2,2,'a')
sage: rng=range(2)
sage: HM(2,2,2,[A[i,j,k]==Dlt[i,j,k] for k in rng for j in rng for i in rng])
[[[a000 == 1, a001 == 0], [a010 == 0, a011 == 0]], [[a100 == 0, a101 == 0], [a110 == 0, a111 == 1]]]
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


## Revisiting Hyperdeterminants

The BM algebra of hypermatrices suggests a natural generalization to the classical matrix determinant which is also 
implemented in the hypermatrix package. The user must be warn that while it is relatively straight forward to compute
determinant of hypermatrices of arbitrary order of with side length equal to 2 via the following expressions 

```python
sage: GeneralHypermatrixDeterminant(2)
-m01*m10 + m00*m11
sage: GeneralHypermatrixDeterminant(3)
m000*m011*m101*m110 - m001*m010*m100*m111
sage: GeneralHypermatrixDeterminant(4)
-m0001*m0010*m0100*m0111*m1000*m1011*m1101*m1110 + m0000*m0011*m0101*m0110*m1001*m1010*m1100*m1111
sage: GeneralHypermatrixDeterminant(5)
m00000*m00011*m00101*m00110*m01001*m01010*m01100*m01111*m10001*m10010*m10100*m10111*m11000*m11011*m11101*m11110 - m00001*m00010*m00100*m00111*m01000*m01011*m01101*m01110*m10000*m10011*m10101*m10110*m11001*m11010*m11100*m11111
```

however computing the determinant of hypermatrices of arbitrary order with side length greater then 3 is considerably 
more difficult as the expression of the proposed generalization of the determinant involves summing over all 
[latin hypercubes](http://en.wikipedia.org/wiki/Latin_hypercube_sampling) of the corresponding order and side length.
We perform a brute forces search for the latin Hypercubes (suggestions for improvements are more then welcome),
consequently the package provides an implementation of the hyperdeterminant of hypermatrices with order less than 7.
The computation of the determinant of a hypermatrix is obtained as follows

```python
sage: Ha=HM(3,3,3,'a')
sage: Ha.det()
-a001*a010*a022*a100*a112*a121*a202*a211*a220 + a000*a012*a021*a101*a110*a122*a202*a211*a220 + a000*a011*a022*a102*a110*a121*a201*a212*a220 - a002*a010*a021*a100*a111*a122*a201*a212*a220 - a000*a011*a022*a101*a112*a120*a202*a210*a221 + a001*a012*a020*a100*a111*a122*a202*a210*a221 + a001*a010*a022*a102*a111*a120*a200*a212*a221 - a002*a011*a020*a101*a110*a122*a200*a212*a221 - a000*a012*a021*a102*a111*a120*a201*a210*a222 + a002*a011*a020*a100*a112*a121*a201*a210*a222 + a002*a010*a021*a101*a112*a120*a200*a211*a222 - a001*a012*a020*a102*a110*a121*a200*a211*a222
```

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


# Bug report

Please report any bugs, it is greatly appreciated.
