```@meta
CurrentModule = Oscar
```

# Overview over Rings, Fields and Algebras in Oscar

Some rings and example elements of those rings.
Note, that if there is some canonical embedding of one ring `R` into another ring `S` and `r` is an element of `r`, usually `S(r)` embeds `r` into `S`. Especially this always works with `R` being the ring of integers `ZZ`. Likewise you can embed Julia’s native integers. So both `S(ZZ(42))` and `S(42)` work and yield the same element of `S`.

## $\mathbb Z$, $\mathbb Q$, $\bar{\mathbb Q}$, $\mathbb R$ and $\mathbb C$
Let `x` be defined by `_, x = PolynomialRing(ZZ)` or the like.

| Symbol            | Name in Oscar                           | Example Elements         |
|:------------------|:----------------------------------------|:-------------------------|
| $\mathbb Z$       | [`ZZ`](@ref ZZ)                         | $42=$`ZZ(42)`            |
| $\mathbb Q$       | [`QQ`](@ref QQ)                         | $\frac{22}7=$`QQ(22, 7)` |
| $\bar{\mathbb Q}$ | [`QQBar`](@ref QQBar)                   | $\sqrt2=$`sqrt(QQBar(2))`, $\varphi=$`roots(x^2+x+1, QQBar)[1]` |
| $\mathbb C$ [^1]  | `CC = `[`CalciumField`](@ref)`()`       | $\sqrt2=$`sqrt(CC(2))`, $\varphi=$`roots(x^2+x+1, CC)[1]`, $1+2\mathrm i=$`CC(1+2im)` |
| $\mathbb R$ [^2]  | `RR = `[`Nemo.RealField`](@ref)`()`     | $\sqrt2\approx$`sqrt(RR(2))`, $\varphi\approx$`roots(x^2+x+1, RR)[2]`, $0.1\approx$`RR("0.1")`, $\frac23\pm0.001\approx$`RR(QQ(2,3)) + RR("0 +/- 0.001")` |
| $\mathbb C$ [^2]  | `CC_ = `[`Nemo.ComplexField`](@ref)`()` | $\sqrt2\approx$`sqrt(CC_(2))`, $\varphi\approx$`roots(x^2+x+1, CC_)[2]`, $0.1+0.3\mathrm i\approx$`CC_("0.1")+CC_("0.3")*onei(CC_)`, $0.1+0.3\mathrm i\pm0.001\pm0.001\mathrm i\approx$`CC_("0.1 +/- 0.001") + CC_("0.3 +/- 0.001")*onei(CC_)` |

## $\mathbb Z/n\mathbb Z$, $\mathbb F_p$, $\mathbb F_{p^n}$, $\mathbb Z_p$, $\mathbb Q_p$, $\mathbb Z_{p^n}$, $\mathbb Q_{p^n}$
Let $p$ be a prime and $n$ a natural number. Additionaly choose some `prec::Int` as precision.

| Symbol                 | Name in Oscar                                  | Example Elements                    |
|:-----------------------|:-----------------------------------------------|:------------------------------------|
| $\mathbb Z/n\mathbb Z$ | `ZZnZ = `[`ResidueRing`](@ref Residue-rings)`(ZZ, n)` | `ZZnZ(42)`                   |
| $\mathbb F_p$          | `FFp = `[`GF`](@ref)`(p)`                      | `FFp(42)`, `FFp(ZZnZ(42))` [^3]     |
| $\mathbb F_{p^n}$      | `FFpn, o = `[`FiniteField`](@ref)`(p, n)`      | `1+o`, `o^-42`                      |
| $\mathbb Q_p$          | `QQp = `[`PadicField`](@ref)`(p, prec)`        | `QQp(QQ(2, 3))`                     |
| $\mathbb Z_p$          | `ZZp = MaximalOrder(QQp)`                      | `ZZp(QQ(p-1, p+1))`                 |
| $\mathbb Q_{p^n}$      | `QQpn, a = `[`QadicField`](@ref)`(p, n, prec)` | `QQpn(QQ(2, 3))`, `a^-42`           |
| $\mathbb Z_{p^n}$      | `ZZpn = MaximalOrder(QQpn)`                    | `ZZpn(QQ(p-1, p+1))`, `ZZpn(a)^-42` |
<!-- Shouldn’t all `ResidueRing`s be subtypes of `ResRing`?
I’m also confused, that ResidueRing itself is no type, but merely a function, that happens to be uppercase, while in other such contexts snakecase is used. Shouldn’t we settle for one case? -->

## Rings from other Rings/Algebras
In the following, let $R$ be some ring and $r1$, $r2$, $r3$ elements of $R$.
Let also $G$ be a group and $n$ a natural number.

| Symbol                 | Name in Oscar                                                  | Aliases                                                                                     | Example Elements                                                                          |
|:-----------------------|:---------------------------------------------------------------|:--------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------|
| $R[x]$                 | `Rx, x = `[`PolynomialRing`](@ref)`(R)`                        | `PolynomialRing(R, "x")`, `R["x"]`, `R[:x]`                                                 | `x^3 + 2x^2 + 3x + 4` $=$ `Rx([4, 3, 2, 1])`                                              |
| $R[x1,x2]$             | `Rx1x2, (x1, x2) = `[`PolynomialRing`](@ref)`(R, 2)`           | `PolynomialRing(R, ["x1", "x2"])`, `R["x1", "x2"]`, `R[:x1, :x2]` [^0],                     | `x1^3*x2 - 2x1*x2^3`, see documentation                                                   |
| $\operatorname{Frac}R$ | `FracR = `[`FractionField`](@ref)`(R)`                         |                                                                                             | `FracR(42)`, `FracR(r1)` $=$ `r1//1`, `FracR(r1, r2)` $=$ `r1//r2` [^4]                   |
| $R/r1$                 | `Rmodr1 = `[`ResidueRing`](@ref)`(R, r1)`                      | `quo(R, r1)[1]` [^0], `quo(R, ideal(R, r1))[1]` [^0], `quo(R, ideal(R, [r1]))[1]` [^0] [^5] | `Rmodr1(r2)`                                                                              |
| $R/(r1, r2)$           | [`Rmodr1r2 = ?`](@ref affine_algebras)                         | `quo(R, [r1, r2])[1]`, `quo(R, ideal(R, [r1, r2]))[1]` [^0] [^5]                            | `Rmodr1r2(r3)` [^6]                                                                       |
| $R[[x]]$               | `Rxx, xx = `[`PowerSeriesRing`](@ref)`(R, prec, "x")` [^7]     |                                                                                             | `42 + O(xx^2)`, `2xx - 1` $=$ `Rxx(2x - 1)` [^0] $=$ `rel_series(R, [-1, 2], 2, prec, 0)` |
| $R[x, 1/x]$            | `R_x, x_ = `[`LaurentPolynomialRing`](@ref)`(R, "x")`          |                                                                                             | `x + 2 + 3x^-1 + 4x^-2` $=$ `R_x(Rx([4, 3, 2, 1]), -2)`                                   |
| $R[[x, 1/x]]$          | `R_xx, xx_ = `[`LaurentSeriesRing`](@ref)`(R, prec, "x")` [^7] |                                                                                             | `42 + O(xx_^2)`, `2xx_^-1 - 1` $=$ `Rxx(2x_^-1 - 1)` [^0]                                 |
| $RG$                   | `RG = `[`group_algebra`](@ref)`(R, G)`                         | `AlgGrp(R, G)`, `R[G]`                                                                      | ...                                                                                       |
| $R^{n\times n}$        | `Rnxn = `[`MatrixAlgebra`](@ref)`(R, n)` [^8]                  |                                                                                             | ...                                                                                       |
| $R^{n\times n}$        | `Rnxn = `[`matrix_algebra`](@ref)`(R, n)` [^8]                 |                                                                                             | ...                                                                                       |


## Number Fields

| Symbol                     | Name in Oscar                            | Example Elements                               |
|:---------------------------|:-----------------------------------------|:-----------------------------------------------|
| $K[x]/f$                   | `L, b = `[`NumberField`](@ref)`(f, "b")` | `b^3 + 2b^2 + 3b + 4` = `L([4, 3, 2, 1])` [^9] |
| $\mathbb Z[x]/f$           | `O = `[`EquationOrder`](@ref)`(L)`       | `O(b)`                                         |
| $\operatorname{IntCls}(L)$ | `O = `[`MaximalOrder`](@ref)`(L)`        | `O(b)`                                         |
| $\operatorname{Ord}(B)$    | `O = `[`Order`](@ref)`(B)`               | `O(B[1])`                                      |

[^0]: Currently not implemented/only implemented in some cases.
[^1]: We cannot provide `CC` as a constant, since (a) `CalciumField`s are not threadsafe and (b) you can give various options to `CalciumField()` to change behaviour.
[^2]: Implements interval/rectangle arithmetic. Not provided as a constant. The precision can be given as the argument to `RealField()`/`ComplexField()`.
[^3]: In the case, this works mathematically, i. e. $p$ divides $n$.
[^4]: Caveat: The `FracR` constructor does not reduce automatically.
[^5]: The second return value is the residue map. Though it is not needed to create elements.
[^6]: Caveat: Elements will not be normalized automatically.
[^7]: Also documented [here](https://nemocas.github.io/AbstractAlgebra.jl/stable/series).
[^8]: The camel case version is from `AbstractAlgebra.Generic` while the snake case version is more specific.
[^9]: The vector length must be the degree of the number field (respectively `f`).
