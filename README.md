---
layout: page-mathjax
---
=======
# Hypermatrix Algebra Package

We provide here a sagemath package implementation for the Mesner-Bhattacharya Hypermatrix Algebra.

**UPDATE 2015-02-15** Major changes to the codebase, as reflected by the current post.

`Hypermatrix Algebra Package` is a symbolic hypermatrix package for investigating
structural and combinatorial properties of the Bhattacharya Mesner Algebra of 
hypermatrices.

# Installation 

A properly working install of [Sage](http://sagemath.org/) is a pre-requisite to using the Hypermatrix Algebra package.
In a sage terminal session, the hypermatrix package is to be loaded into sage via the following command:

```python
sage: %runfile("Hypermatrix_Algebra_Package_code.sage")
```

# Usage

To creates a symbolic Hypermatrix instance of size \(2\times 2 \times 2\), we use the instruction

```python
sage: Ha = HM(2,3,4,'a')
```

alternatively Hypermatrices can be initialized from a list as follows

```python
sage: Hx = HM(2,2,2,[var('x'+str(i)+str(j)+str(k)) for k in range(1,3) for j in range(1,3) for i in range(1,3)])
```

Hypermatrix entries are not restricted to numbers and symbolic expression, they can be matrices and even hypermatrices
as illustrated below

```python
sage: Hh = HM(2,2,2, [ (i*2^2+j*2+k)*hadamard_matrix(2) for k in range(2) for j in range(2) for k in range(2)])
```

Hypermatrix entries are accessed quite similarly to the way Matrix entries are accessed in sage. For example,
the 0,1,0 entry of the hypermatrix above is accessed as follows

```python
sage: Hh[0,1,0]
```

# Basic hypermatrix operation

The [BM product](http://arxiv.org/abs/1411.6270) is implemented for hypermatrices of all orders
and all sizes. The hypermatrix product is performed as follows

```python
sage: Prod(HM(2,3,4,'a'), HM(2,2,3,'b'), HM(3,2,4,'c'))
```

The product of second order hypermatrices recovers the usual matrix product. Other basic hypermatrix operations
including addition and multiplication by scalars is quite similar to their matrix counter part. The hypermatrix
transpose amounts to a clyclic permutation of the entry indices and is performed as follows 

```python
sage: HM(2,2,2,'a').transpose()
```

In order to perform two consecutive transposes we use the following instruction

```python
sage: HM(2,2,2,'a').transpose(2)
```

The sum function call can be used to sum over a list of hypermatrices. We use this fact
to illustrate the construction of a symmetric hypermatrix.

```python
sage: sum([ HM(2,2,2,'a').transpose(i) for i in range(3)])
```

Many of the properties of the special hypermatrices that we describe subsequently are preserved by the Kroencker product.
The package provides an implementation of the Kronecker product for hypermatrices of orders going up to 5 (possibly later
verision of the package will address this limitation). The Kroneckr product of two hypermatrices are obtained as follows

```python
sage: A = HM(2,2,2,'a'); B = HM(3,3,3,'b')
sage: A.slicekroneckerproduct(B)
```


# Symbolic parametrization of some special hypermatrices

To emphasize that hypermatricess naturally extend the algebra of matrices, we describe here instructions for obtaining
some symbolic parametrization of special families of hypermatrices.

## Orthogonal matrix

Orthogonal hypermatrices are analogous to Orthogonal matrices in the fact that the product of their transposes yields
the [Kronecker delta](http://en.wikipedia.org/wiki/Kronecker_delta#Properties_of_generalized_Kronecker_delta) of the 
same size and order. A symbolic parametetrization of $2\times 2\times 2$ hypermatrices is given by

```python
sage: Q = GeneralOrthogonalHypermatrix(3)
sage: Prod(Q,Q.transpose(2),Q.transpose()).simplify()
```

Similarly a symbolic parametrization of orthogonal parametrization of $2\times 2\times 2\times 2$ hypermatrices is given
by

```python
sage: Q = GeneralOrthogonalHypermatrix(4)
sage: Prod(Q,Q.transpose(3),Q.transpose(2),Q.transpose()).simplify()
```

The parameters which appear in the parametrization can be assigned arbitrary values (be mindful of division by zeros).
The Hypermatrix package provides construction of third order Hadamard Hypermatrices, which are analogous to Hadamard 
matrices in the fact that the product of the transposes equals to the Kronecker delta of the same order multiplied
by the side lenght of the cubic hypermatrix. It is infact easy to see that the classical Hadamard matrix 
[construction](http://en.wikipedia.org/wiki/Hadamard_matrix#Sylvester.27s_construction) due to James Joseph Sylvester and
going back to 1867 extends to all Hypermatrices of prime order ( recall that the order refers to the number of indices 
associated with the hypermatrix entries for instance matrices are second order hypermatrices because every entry has a 
row and column index). Third order Hadamard hypermatrices whose side length are powers of 2 and in particular in our 
example of size $4\times 4\times 4$

```python
sage: H = ThirdOrderHadamardBlockU(4)
sage: Prod(H,H.transpose(2),H.transpose())
```

One can not help but wonder if the reach of the famous Hadamard conjecture should be extended to hypermatrices of prime 
order.


## Unitary hypermatrices
The Hypermatrix package also provides a symbolic parametrization of even order unitary Hypermatrix with sides length 
equal to 2. As suggested by [Gowers norm](http://en.wikipedia.org/wiki/Gowers_norm) an even order hypermatrix is said
to be unitary if the product of the conjugate transposes yield the Kronecker delta of the same order. Consequently,
second order Unitary matrices correspond to the usual unitary matrices and a symbolic parametrization of $2\times 2$
unitary matrices is obtained via the commands

```python
sage: [U,Uc] = GeneralUnitaryHypermatrix(2) 
sage: Prod(U,Uc.transpose()).simplify()
```

Similarly a symblic paramterization of a fourth order unitary hypermatrix of side length 2 is given by 

```python
sage: [U,Uc] = GeneralUnitaryHypermatrix(4) 
sage: Prod(U, Uc.transpose(3), U.transpose(2), Uc.transpose()).simplify()
```

Note that the symbolic parametrization return the Unitary Hypermatrix and it's complex conjugate. It is therefore important
to point out that all paramaeter are to be assign real values in the previous parametrization.

## Uncorrelated tuples: Generalization of matrix inverse pairs
The [General linear group](http://en.wikipedia.org/wiki/General_linear_group) plays a crucial role in many areas of
mathematics. The BM algebra suggest natural generalization for the General linear group, for which very little is known.
Recall that a pair of square matrices are said to form inverse pairs if their product yields the Kronecker delta of the
same size and order. Similarly a tuple of cubic Hypermatrices are said to form an uncorrelated tuple if their product
yield the Kronecker delta of the same order and size. The Hypermatrix package implements a symbolic parametrization for
uncorrelated tuples of arbitrary order but of side lenght equal to 2 as follows

```python
sage: [Ha, Hb]=GeneralUncorrelatedHypermatrixTuple(2); Prod(Ha,Hb).simplify() 
sage: [Ha, Hb, Hc]=GeneralUncorrelatedHypermatrixTuple(3); Prod(Ha,Hb,Hc).simplify()
sage: [Ha, Hb, Hc, Hd]=GeneralUncorrelatedHypermatrixTuple(4); Prod(Ha,Hb,Hc,Hd).simplify() 
```

## Revisiting Hyperdeterminant
The BM algebra of hypermatrices suggests very natural generalization to the matrix determinant which have been 
incorporated in the Hypermatrix package. The user must be warn that it is relatively straight forward to compute
determinant of hypermatrices of arbitrary order of with side length equal to two and the expression of such determinant
in terms of the hypermatrix entries is provided by the follwoing instruction

```python
sage: GeneralHypermatrixDeterminant(2)
sage: GeneralHypermatrixDeterminant(3)
sage: GeneralHypermatrixDeterminant(4) 
sage: GeneralHypermatrixDeterminant(5) 
```

However computing determinant of hypermatrices whose side length are greater then 3 is considerably more difficult
as the expression of the proposed definition of the determinant involves all 
[Latin Hypercube](http://en.wikipedia.org/wiki/Latin_hypercube_sampling) of the corresponding order and side lenght.
Unefortunately as our implementation suggest we perform a brute forces search for the latin Hypercubes (suggestions for
improvemements are more then welcome). The package provides an implementation of the hyperdeterminant of 
hypermatrices of arbitrary size but order less than 7. The computation of the determinant of a Hypermatrix is obtained as
follows

```python
sage: Hc = HM(3,3,3,'c')
sage: Hc.det()
```



