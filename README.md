---
layout: page-mathjax
---
=======
# Hypermatrix Algebra Package

We provide here a sagemath implementation of the Bhattacharya-Mesner(BM) algebra and it's
generalizations.

**UPDATE 2015-02-15** Major changes to the codebase, as reflected by the current post.

The `Hypermatrix Algebra Package` is a symbolic hypermatrix package designed to
investigate structural and combinatorial properties of hypermatrices.

# Installation 

A properly working install of [Sage](http://sagemath.org/) is a prerequisite to using the hypermatrix 
package. In a sage terminal session, the hypermatrix package is to be loaded via the following 
command:

```python
sage: %runfile("Hypermatrix_Algebra_Package_code.sage")
```

# Usage

To create a symbolic hypermatrix instance of size 2x3x4 for example, we use the instruction

```python
sage: Ha = HM(2,3,4,'a'); Ha
```

Alternatively hypermatrices can be initialized from a list as follows

```python
sage: Hx=HM(2,2,2,[var('x'+str(i)+str(j)+str(k)) for k in range(1,3) for j in range(1,3) for i in range(1,3)]); Hx
```

The hypermatrix entries are not restricted to numbers and symbolic expression, they may be matrices or more generally
can be other hypermatrices as illustrated below

```python
sage: Hh=HM(2,2,2,[(i*2^2+j*2+k)*hadamard_matrix(2) for k in range(2) for j in range(2) for i in range(2)]); Hh
```

Hypermatrix entries are accessed quite similarly to the way in which matrix entries are usually accessed in sage.
For example, the (0,1,0) entry of the hypermatrix constructed in the previous instruction is accessed as follows

```python
sage: Hh[0,1,0]
```

# Basic hypermatrix operation

The [BM product](http://arxiv.org/abs/1411.6270) is implemented for hypermatrices of all orders and all 
compatible sizes. The hypermatrix product is performed as follows

```python
sage: Prod(HM(2,3,'a'),HM(3,4,'c'))
sage: Prod(HM(2,3,4,'a'),HM(2,2,3,'b'),HM(3,2,4,'c'))
```

As illustrated by the previous instructions, the product of compatible second order hypermatrices recover the 
usual matrix product. Other basic hypermatrix operations including addition and multiplication by scalars are
quite similar to their matrix counter part. The hypermatrix transpose amounts to a clyclic permutation of the
entry indices and is performed as follows 

```python
sage: HM(2,2,2,'a').transpose()
```

In order to perform two consecutive transposes we use the following instruction

```python
sage: HM(2,2,2,2,'a').transpose()
sage: HM(2,2,2,2,'a').transpose(2)
sage: HM(2,2,2,2,'a').transpose(3)
sage: HM(2,2,2,2,'a').transpose(4) - HM(2,2,2,2,'a')
```

The sum of hypermatrices can taken over a list of hypermatrices.
We ilustrate this by constructing a symbolic symmetric hypermatrix.

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

To emphasize that hypermatricess naturally extend the algebra of matrices, we describe here instructions for obtaining
some symbolic parametrization of special families of hypermatrices.

## Orthogonal hypermatrices

Orthogonal hypermatrices are analogous to Orthogonal matrices in the fact that the product of their transposes yields
the [Kronecker delta](http://en.wikipedia.org/wiki/Kronecker_delta#Properties_of_generalized_Kronecker_delta) of the 
same size and order (Kronecker delta is the hypermatrix whose non zero entries equal 1 and correspond the entries whose
indices are all equal as illustrated by the identity matrices). A symbolic parametetrization of 2x2x2 hypermatrices is 
given by

```python
sage: Q=GeneralOrthogonalHypermatrix(3); Q
sage: Prod(Q,Q.transpose(2),Q.transpose()).simplify()
```

Similarly a symbolic parametrization of orthogonal parametrization of 2x2x2x2 hypermatrices is given
by

```python
sage: Q=GeneralOrthogonalHypermatrix(4); Q
sage: Prod(Q,Q.transpose(3),Q.transpose(2),Q.transpose()).simplify()
```

The parameters which appear in the parametrization can be assigned arbitrary complex values (one must be mindful of 
possible division by zeros when assigning values to the parameters).
The hypermatrix package also provides construction of cubic third order Hadamard hypermatrices, which are quite analogous
to Hadamard matrices in the fact that the product of the transposes equals to the Kronecker delta ( of the same size and 
order ) scaled by the side lenght of the hypermatrix. It is infact easy to see that the classical Hadamard matrix 
[construction](http://en.wikipedia.org/wiki/Hadamard_matrix#Sylvester.27s_construction) proposed by 
James Joseph Sylvester and going back to 1867 extends to all hypermatrices of prime order ( recall that the order refers
to the number of indices associated with the hypermatrix entries, for instance matrices are second order hypermatrices 
because every matrix entry has a row and column index). Third order Hadamard hypermatrices whose side length are powers 
of 2 and in the particular case of size 4x4x4 is obtained as follows

```python
sage: H=ThirdOrderHadamardBlockU(4); H
sage: Prod(H, H.transpose(2), H.transpose())
```

One can not help but wonder if the reach of the famous Hadamard matrix conjecture should be extended to all hypermatrices 
of prime order.


## Unitary hypermatrices
The hypermatrix package also provides a symbolic parametrization of even order unitary hypermatrix with sides length 
equal to 2. As suggested by the [Gowers norm](http://en.wikipedia.org/wiki/Gowers_norm), an even order hypermatrix is
said to be unitary if the product of the conjugate transposes yield the Kronecker delta of the same order and size.
Consequently, second order unitary hypermatrices correspond to the usual unitary matrices and a symbolic parametrization
of 2x2 unitary matrices is obtained from the hypermatrix package via the commands

```python
sage: [U,Uc]=GeneralUnitaryHypermatrix(2); U 
sage: Prod(U,Uc.transpose()).simplify()
```

Similarly a symblic paramterization of a fourth order unitary hypermatrix of side length 2 is obtained as follows 

```python
sage: [U,Uc]=GeneralUnitaryHypermatrix(4); U 
sage: Prod(U, Uc.transpose(3), U.transpose(2), Uc.transpose()).simplify()
```

Note that the symbolic parametrization returns the unitary hypermatrix and it's complex conjugate. It is therefore 
important to emphasize that all the parameters in the parametrization are to be assigned real values for the hypermatrix
to be unitary.

## Uncorrelated tuples: Generalization of matrix inverse pair
The [General linear group](http://en.wikipedia.org/wiki/General_linear_group) plays a crucial role in many areas of
mathematics. The BM algebra suggest natural generalization for the General linear group, for which very little is known.
Recall that a pair of square matrices are said to form an inverse pair if their product yields the Kronecker delta of the
same size and order. Similarly a tuple of cubic hypermatrices are said to form an uncorrelated tuple if their product
yield the Kronecker delta of the same order and size. The hypermatrix package implements a symbolic parametrization for
uncorrelated tuples of arbitrary order but of side lenght equal to 2 and are obtained as follows

```python
sage: [Ha,Hb]=GeneralUncorrelatedHypermatrixTuple(2); Prod(Ha,Hb).simplify() 
sage: [Ha,Hb,Hc]=GeneralUncorrelatedHypermatrixTuple(3); Prod(Ha,Hb,Hc).simplify()
sage: [Ha,Hb,Hc,Hd]=GeneralUncorrelatedHypermatrixTuple(4); Prod(Ha,Hb,Hc,Hd).simplify() 
```

## Multistochastic hypermatrices
The hypermatrix package also provide and symbolic parametrization of multistochastic hypermartices of arbitrary orders
and side lenght equal to 2 obtained as follows

```python
sage: GeneralStochasticHypermatrix(var('x'), 2) 
sage: GeneralStochasticHypermatrix(var('x'), 3) 
sage: GeneralStochasticHypermatrix(var('x'), 4) 
```


## Revisiting Hyperdeterminants
The BM algebra of hypermatrices suggests very natural generalization to the matrix determinant which has also been 
implemented in the hypermatrix package. The user must be warn that while it is relatively straight forward to compute
determinant of hypermatrices of arbitrary order of with side length equal to 2 via the following expressions 

```python
sage: GeneralHypermatrixDeterminant(2)
sage: GeneralHypermatrixDeterminant(3)
sage: GeneralHypermatrixDeterminant(4) 
sage: GeneralHypermatrixDeterminant(5) 
```

however computing determinant of hypermatrices of arbitrary order with side length greater then 3 is considerably 
more difficult as the expression of the proposed generalization of the determinant involves all 
[latin hypercubes](http://en.wikipedia.org/wiki/Latin_hypercube_sampling) of the corresponding order and side lenght.
Unefortunately as our implementation suggestis we perform a brute forces search for the latin Hypercubes (suggestions for
improvemements are more then welcome). Consequently the package provides an implementation of the hyperdeterminant of 
hypermatrices of arbitrary size but order less than 7. The computation of the determinant of a hypermatrix is obtained 
as follows

```python
sage: Hc=HM(3,3,3,'c')
sage: Hc.det()
```

# Bug report
Please report any bugs, it will be greatly appreciated.
