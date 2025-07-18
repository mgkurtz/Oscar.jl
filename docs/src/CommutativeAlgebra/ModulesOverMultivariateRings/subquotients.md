```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Subquotients

A subquotient  is a submodule of a quotient of a free module. In this section, the expression
*subquotient* refers to a subquotient over a ring of type `MPolyRing`, `MPolyQuoRing`,
`MPolyLocRing`, `MPolyQuoLocRing`, `ZZRing`, or  `Field`. That is, given a ring $R$ of one of these
types, a subquotient $M$ over $R$ is a module of type

$M = (\text{im } a + \text{im } b)/\text{im } b,$

where

$a:R^s ⟶R^p \;\text{ and }\; b:R^t ⟶R^p$

are two homomorphisms of free $R$-modules with the same codomain. We then refer to
- the module $M$ as the *subquotient defined by $a$ and $b$*,
- the codomain $R^p$ as the *ambient free module* of $M$,
- the images of the canonical basis vectors of $R^s$ in $R^p$ as the *ambient representatives of the generators* of $M$, and
- the images of the canonical basis vectors of $R^t$ in $R^p$ as the *relations* of $M$.

Alternatively, we speak of the *subquotient of* $\;\text{im } a\;$ *by* $\;\text{im } b\;$ or the
*subquotient defined by $A$ and $B$*, where $A$ and $B$ are the matrices representing
$a$ and $b$, respectively.

Finally, we refer to
- the quotient of $R^p$ by the submodule generated by the relations of $M$ as the *ambient module of $M$*,

and regard $M$ as a submodule of that ambient module, embedded in the natural way.

!!! note
    Recall from the section on [free modules](@ref free_modules) that by a free $R$-module we mean a free
    module of type $R^p$ , where we think of $R^p$ as a free module with a given
    basis, namely the basis of standard unit vectors. Accordingly, elements of free modules
    are represented by coordinate vectors, and homomorphisms between free modules by
    matrices. Here, by convention, vectors are row vectors, and matrices operate by
    multiplication on the right.

!!! note
    Over a graded ring $R$, we work with graded free modules $R^s$, $R^p$, $R^t$ and graded
    homomorphisms $a$, $b$. As a consequence, every module involved in the construction of
    the subquotient defined by $a$ and $b$ carries an induced grading. In particular, the
    subquotient itself carries an induced grading.

## Types

All OSCAR types for the finitely presented modules considered here belong to the
abstract type `ModuleFP{T}`, where `T` is the element type of the underlying ring.
Graded or not, the subquotients belong to the abstract subtype `AbstractSubQuo{T} <: ModuleFP{T}`,
they are modeled as objects of the concrete type `SubquoModule{T} <: AbstractSubQuo{T}`.

!!! note
    Canonical maps such us the canonical projection onto a quotient module arise in many 
    constructions in commutative algebra. The `SubquoModule` type is designed so that it allows
    for the caching of such maps when executing functions. The `tensor_product`
    function discussed in this section provides an example.

## Constructors

```@docs
subquotient(a::FreeModuleHom, b::FreeModuleHom)
```

## Data Associated to Subqotients

If `M` is a subquotient with ambient free `R`-module `F`, then

- `base_ring(M)` refers to `R`,
- `ambient_free_module(M)` to `F`,
- `gens(M)` to the generators of `M`, 
- `number_of_generators(M)` / `ngens(M)` to the number of these generators, 
- `M[i]`, `gen(M, i)` to the `i`th such generator,
- `ambient_representatives_generators(M)` to the ambient representatives of the generators of `M` in `F`,
- `relations(M)` to the relations of `M`, and
- `ambient_module(M)` to the ambient module of `M`.

##### Examples

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over R

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = subquotient(F, A, B)
Subquotient of submodule with 2 generators
  1: x*e[1]
  2: y*e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> base_ring(M)
Multivariate polynomial ring in 3 variables x, y, z
  over rational field

julia> F === ambient_free_module(M)
true

julia> gens(M)
2-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 x*e[1]
 y*e[1]

julia> number_of_generators(M)
2

julia> gen(M, 2)
y*e[1]

julia> ambient_representatives_generators(M)
2-element Vector{FreeModElem{QQMPolyRingElem}}:
 x*e[1]
 y*e[1]

julia> relations(M)
3-element Vector{FreeModElem{QQMPolyRingElem}}:
 x^2*e[1]
 y^3*e[1]
 z^4*e[1]

julia> ambient_module(M)
Subquotient of submodule with 1 generator
  1: e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]
```

In the graded case, we also have:

```@docs
 grading_group(M::SubquoModule)
```

```@docs
degrees_of_generators(M::SubquoModule)
```

## Elements of Subqotients

All OSCAR types for elements of finitely presented modules considered here belong to the
abstract type `ModuleElemFP{T}`, where `T` is the element type of the polynomial ring.
For elements of subquotients, there  are the abstract subtype `AbstractSubQuoElem{T} <: ModuleFPElem{T}`
and its concrete descendant `SubquoModuleElem{T}` which implements an element $m$ of a subquotient
$M$ over a ring $R$ as a sparse row, that is, as an object of type `SRow{T}`.
This object specifies the coefficients of an $R$-linear combination of the generators of $M$
giving $m$. To create an element, enter the coefficients as a sparse row or a vector: 

```@julia
(M::SubquoModule{T})(c::SRow{T}) where T
```

```@julia
(M::SubquoModule{T})(c::Vector{T}) where T
```

Alternatively, directly write the element as an $R$-linear combination of generators of $M$.

##### Examples

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over R

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = subquotient(F, A, B)
Subquotient of submodule with 2 generators
  1: x*e[1]
  2: y*e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> m = M(sparse_row(R, [(1,z),(2,one(R))]))
(x*z + y)*e[1]

julia> n = M([z, one(R)])
(x*z + y)*e[1]

julia> o = z*M[1] + M[2]
(x*z + y)*e[1]

julia> m == n == o
true

```

Given an element `m` of a subquotient `M` over a ring $R$ with element type `T`,
- `parent(m)` refers to `M`, 
- `coordinates(m)` to an object of type `SRow{T}` specifying the coefficients of an $R$-linear combination of the generators of $M$ which gives $m$, and
- `ambient_representative(m)` to an element of the ambient free module of `M` which represents `m`.

Given an element `f` of the ambient free module of a subquotient `M` such that `f` represents an element of `M`,
the function below creates the represented element:

```@julia
(M::SubquoModule{T})(f::FreeModElem{T}; check::Bool = true) where T
```

By default (`check = true`), it is tested whether `f` indeed represents an element of `M`.
If this is already clear, it may be convenient to omit the test (`check = false`).

##### Examples

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over R

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = subquotient(F, A, B)
Subquotient of submodule with 2 generators
  1: x*e[1]
  2: y*e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> m = z*M[1] + M[2]
(x*z + y)*e[1]

julia> parent(m)
Subquotient of submodule with 2 generators
  1: x*e[1]
  2: y*e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> coordinates(m)
Sparse row with positions [1, 2] and values QQMPolyRingElem[z, 1]

julia> fm = ambient_representative(m)
(x*z + y)*e[1]

julia> typeof(m)
SubquoModuleElem{QQMPolyRingElem}

julia> typeof(fm)
FreeModElem{QQMPolyRingElem}

julia> parent(fm) === ambient_free_module(M)
true

julia> F = ambient_free_module(M)
Free module of rank 1 over R

julia> f = x*F[1]
x*e[1]

julia> M(f)
x*e[1]

julia> typeof(f)
FreeModElem{QQMPolyRingElem}

julia> typeof(M(f))
SubquoModuleElem{QQMPolyRingElem}

```

The zero element of a subquotient is obtained as follows:

```@docs
zero(M::SubquoModule)
```

Whether a given element of a subquotient is zero can be tested as follows:

```@docs
is_zero(m::SubquoModuleElem)
```

In the graded case, we additionally have:

```@docs
is_homogeneous(m::SubquoModuleElem)
```

```@docs
degree(m::SubquoModuleElem)
```

## Tests on Subqotients

The functions [`is_graded`](@ref), [`is_standard_graded`](@ref), [`is_z_graded`](@ref),
and [`is_zm_graded`](@ref) carry over analogously to subquotients. They return `true` if the
respective property is satisfied, and `false` otherwise. In addition, we have:

```@docs
==(M::SubquoModule{T}, N::SubquoModule{T}) where T
```

```@docs
is_subset(M::SubquoModule{T}, N::SubquoModule{T}) where T
```

```@docs
is_zero(M::SubquoModule)
```

## Basic Operations on Subquotients


```@docs
:+(M::SubquoModule{T},N::SubquoModule{T}) where T
```

```@docs
sum(M::SubquoModule{T},N::SubquoModule{T}) where T
```

```@docs
intersect(M::SubquoModule{T}, N::SubquoModule{T}) where T
```

```@docs
quotient(M::SubquoModule{T}, N::SubquoModule{T}) where T
```

```@docs
quotient(M::SubquoModule, J::Ideal)
```

```@docs
saturation(M:: SubquoModule, J::Ideal)
```

```@docs
saturation_with_index(M:: SubquoModule, J::Ideal)
```


## Submodules and Quotients

```@docs
sub(M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}) where T
```

```@docs
quo(M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}; cache_morphism::Bool=false) where T
```

## Homomorphisms From Subqotients

All OSCAR types for homomorphisms of finitely presented modules considered here belong
to the abstract type `ModuleFPHom{T1, T2}`, where `T1` and `T2` are the types of domain and codomain respectively.
For homomorphisms from subquotients, OSCAR provides the concrete type `SubQuoHom{T1, T2} <: ModuleFPHom{T1, T2}`
as well as the following constructors:

```@docs
hom(M::SubquoModule{T}, N::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}) where T
```

Given a homomorphism of type `SubQuoHom`, a matrix `A` representing it
is recovered by the following function:

```@docs
matrix(a::SubQuoHom)
```

The domain and codomain of a homomorphism `a`  of type `SubQuoHom` can be
recovered by entering `domain(a)` and `codomain(a)`, respectively.

The functions below test whether a homomorphism of type
`SubQuoHom` is graded and homogeneous, respectively.

```@docs
is_graded(a:: SubQuoHom)
```

```@docs
is_homogeneous(a:: SubQuoHom)
```

In the graded case, we additionally have:

```@docs
degree(a:: SubQuoHom)
```

```@docs
grading_group(a:: SubQuoHom)
```



