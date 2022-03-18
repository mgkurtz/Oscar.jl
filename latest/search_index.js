var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Welcome to Oscar",
    "title": "Welcome to Oscar",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Welcome-to-Oscar-1",
    "page": "Welcome to Oscar",
    "title": "Welcome to Oscar",
    "category": "section",
    "text": "Oscar is a new computer algebra system under development."
},

{
    "location": "Rings/integer.html#",
    "page": "Integers",
    "title": "Integers",
    "category": "page",
    "text": "CurrentModule = Oscarusing OscarPages = [\"integer.md\"]"
},

{
    "location": "Rings/integer.html#Integers-1",
    "page": "Integers",
    "title": "Integers",
    "category": "section",
    "text": "An important design decision in Oscar.jl is to use Julia as the user language by default. This means that integers typed at the REPL are Julia integers. However, for performance reasons, Oscar has its own integer format.Julia has a number of different integer types, but the two that are most relevant here are Base.Int and Base.BigInt. All the Julia integer types belong to Base.Integer.The Base.Int type is for machine integers which are highly efficient, but can only represent integers up to a certain hardware defined size before wrapping around.The Base.BigInt type is backed by GMP multiprecision integers and can represent integers whose size is usually only limited by available memory.Oscar currently only has one integer type, Oscar.fmpz which for performance reasons scales internally from machine integers to GMP multiprecision integers. The Oscar integer type belongs to Oscar.Integer.This situation is illustrated in the following diagram.(Image: alt text)In the documentation below, we always use Base.Integer for a Julia integer and Oscar.Integer for an Oscar integer. Some functions accept only machine integers for certain arguments; in such cases, we refer to Base.Int."
},

{
    "location": "Rings/integer.html#The-ring-of-integers-1",
    "page": "Integers",
    "title": "The ring of integers",
    "category": "section",
    "text": "Every object in Oscar representing a mathematical element has a parent. This is an object encoding information about where that element belongs.The parent of an Oscar integer is the ring of integers ZZ.ZZ"
},

{
    "location": "Rings/integer.html#Integer-constructors-1",
    "page": "Integers",
    "title": "Integer constructors",
    "category": "section",
    "text": "Oscar integers are created using ZZ:ZZ(2)^100\nZZ(618970019642690137449562111)One can also construct the integer 0 with the empty constructor:ZZ()The following special constructors are also provided:zero(ZZ)\none(ZZ)zero(ZZ)\none(ZZ)Note that ZZ is not a Julia type, but the above methods of constructing Oscar integers are similar to the way that Julia integer types can be used to construct Julia integers.Int(123)\nBigInt(123456343567843598776327698374259876295438725)\nzero(BigInt)\none(Int)"
},

{
    "location": "Rings/integer.html#Limitations-1",
    "page": "Integers",
    "title": "Limitations",
    "category": "section",
    "text": "Oscar integers have the same limitations as GMP multiprecision integers, namely that they are limited by the available memory on the machine and in any case to signed integers whose absolute value does not exceed 2^37 bits.note: Note\nThe Julia Int type is either a 32 or 64 bit integer, depending on the machine architecture (usually 64 bits on most modern machines). The range of values is machine dependent, but can be found by typing typemin(Int) and typemax(Int) in Julia."
},

{
    "location": "Rings/integer.html#Julia-integers-in-Oscar-functions-1",
    "page": "Integers",
    "title": "Julia integers in Oscar functions",
    "category": "section",
    "text": "For convenience, all basic arithmetic and exact division functions in Oscar also accept Julia integers. For example:divexact(ZZ(234), 2)In this example, 2 is a Julia integer but is still valid in the call to the Oscar function divexact."
},

{
    "location": "Rings/integer.html#Predicates-1",
    "page": "Integers",
    "title": "Predicates",
    "category": "section",
    "text": "iszero(n::Oscar.Integer) -> Bool\nisone(n::Oscar.Integer) -> Bool\nisunit(n::Oscar.Integer) -> Bool\nisodd(n::Oscar.Integer) -> Bool\niseven(n::Oscar.Integer) -> Bool\nissquare(n::Oscar.Integer) -> Bool\nisprime(n::Oscar.Integer) -> Bool\nisprobable_prime(n::Oscar.Integer) -> BoolThe isprime predicate will prove primality, whereas isprobable_prime may declare a composite number to be prime with very low probability.Negative numbers, 0 and 1 are not considered prime by isprime and isprobable_prime.isone(ZZ(1))\nisunit(ZZ(-1))\nissquare(ZZ(16))\nisprobable_prime(ZZ(23))"
},

{
    "location": "Rings/integer.html#Properties-1",
    "page": "Integers",
    "title": "Properties",
    "category": "section",
    "text": "sign(n::Oscar.Integer) -> Oscar.IntegerReturns the sign of n, i.e. nn if n neq 0, or 0 otherwise.sign(ZZ(23))\nsign(ZZ(0))\nsign(ZZ(-1))abs(n::Oscar.Integer) -> Oscar.IntegerReturn the absolute value of n, i.e. n if n geq 0 and -n otherwiseabs(ZZ(-3))"
},

{
    "location": "Rings/integer.html#Basic-arithmetic-1",
    "page": "Integers",
    "title": "Basic arithmetic",
    "category": "section",
    "text": "Oscar provides the basic arithmetic operations +, - and * and comparison operators ==, !=, <, <=, >, >=, including mixed operations between Julia and Oscar integers. It also provides division and powering as described below."
},

{
    "location": "Rings/integer.html#Division-in-Oscar-1",
    "page": "Integers",
    "title": "Division in Oscar",
    "category": "section",
    "text": "Oscar distinguishes a number of different kinds of division:Exact division (divexact)\nEuclidean division (div, rem, divrem and mod)\nConstruction of fractions (a//b)\nFloating point division (a/b)\nDivisibility testing (divides)These choices have been made for maximum parsimony with the Julia language.note: Note\nIt is a common error to enter 1/2 for the fraction \'one half\' in Julia. This expression is reserved for floating point division. Instead, the double slash operator // should be used for fractions."
},

{
    "location": "Rings/integer.html#integer_exact_division-1",
    "page": "Integers",
    "title": "Exact Division",
    "category": "section",
    "text": "divexact(a::Oscar.Integer, b::Oscar.Integer) -> Oscar.IntegerOne or both arguments may be Julia integers, however if they both are, the result will be a Julia integer not an Oscar integer.Returns the quotient of a by b. The result of the exact division of two integers will always be another integer. Exact division raises an exception if the division is not exact, or if division by zero is attempted.divexact(ZZ(6), ZZ(3))\ndivexact(ZZ(6), ZZ(0))\ndivexact(ZZ(6), ZZ(5))\ndivexact(ZZ(6), 2)"
},

{
    "location": "Rings/integer.html#Powering-1",
    "page": "Integers",
    "title": "Powering",
    "category": "section",
    "text": "^(a::Oscar.Integer, b::Base.Int) -> Oscar.IntegerPowering can be accomplished naturally using the special caret infix operator:ZZ(37)^37\nZZ(1)^(-2)note: Note\nAn exception will be raised if an integer other than -1 or 1 is raised to a negative exponent.note: Note\nIn Julia 2^-2 is called a literal power. The value returned is a floating point value. To get behaviour that agrees with Oscar, one can write 2^Int(-2).The following is allowed for convenience.ZZ(0)^0note: Note\nIn Julia, 2^64 will return zero, as the Julia integer 2 is a machine integer. In Oscar, the expression ZZ(2)^64 will return the expected result, just as the Julia equivalent BigInt(2)^64 does."
},

{
    "location": "Rings/integer.html#integer_euclidean_division-1",
    "page": "Integers",
    "title": "Euclidean division",
    "category": "section",
    "text": "The ring of integers is a Euclidean domain and Oscar provides Euclidean division through the functions divrem, div and rem.Integer Euclidean division of a by b computes a quotient and remainder such thata = qb + rwith r  b."
},

{
    "location": "Rings/integer.html#Division-with-remainder-1",
    "page": "Integers",
    "title": "Division with remainder",
    "category": "section",
    "text": "divrem(a::Oscar.Integer, b::Oscar.Integer) -> (Oscar.Integer, Oscar.Integer) : division with remainder\ndiv(a::Oscar.Integer, b::Oscar.Integer) -> Oscar.Integer : quotient only\nrem(a::Oscar.Integer, b::Oscar.Integer) -> Oscar.Integer : remainder onlyOne or both arguments may be Julia integers, however if they both are, the result will be a Julia integer not an Oscar integer.Both rem and divrem compute the remainder r such that when r neq 0 the sign of r is the same as the sign of a.All three functions raise an exception if the modulus b is zero.divrem(ZZ(5), ZZ(3))\ndiv(ZZ(7), ZZ(2))\nrem(ZZ(4), ZZ(3))\ndiv(ZZ(2), ZZ(0))note: Note\nThe rem function does not provide a minimal set of representatives, e.g. rem(-2, 3) = -2 but rem(1, 3) = 1."
},

{
    "location": "Rings/integer.html#Modular-arithmetic-1",
    "page": "Integers",
    "title": "Modular arithmetic",
    "category": "section",
    "text": ""
},

{
    "location": "Rings/integer.html#Modular-reduction-1",
    "page": "Integers",
    "title": "Modular reduction",
    "category": "section",
    "text": "mod(a::Oscar.Integer, b::Oscar.Integer) -> Oscar.Integer : remainder onlyOne or both arguments may be Julia integers, however if they both are, the result will be a Julia integer not an Oscar integer.The mod function computes a remainder r such that when r neq 0 the sign of r is the same as the sign of b. Thus, if b  0 then mod(a, b) will be in the range 0 b). An exception is raised if the modulus b is zero. This is summarised in the following table.remainder division sign rounding\nrem div/divrem same as dividend towards zero\nmod  same as divisor towards -inftyThere is no function implemented to compute the quotient corresponding to the remainder given by mod.mod(ZZ(4), ZZ(3))\nmod(ZZ(2), ZZ(0)) "
},

{
    "location": "Rings/integer.html#integer_divisibility_testing-1",
    "page": "Integers",
    "title": "Divisibility testing",
    "category": "section",
    "text": "divides(a::Oscar.Integer, b::Oscar.Integer) -> (Bool, Oscar.Integer)In Oscar, we say that b divides a if there exists c in the same ring such that a = bc.The call divides(a, b) returns a tuple (flag, q) where flag is either true if b divides a in which case q will be a quotient, or flag is false if b does not divide a in which case q will be an integer whose value is not defined.divides(ZZ(6), ZZ(3))\ndivides(ZZ(5), ZZ(2))Note that for convenience we define:divides(ZZ(0), ZZ(0))"
},

{
    "location": "Rings/integer.html#Greatest-common-divisor-1",
    "page": "Integers",
    "title": "Greatest common divisor",
    "category": "section",
    "text": ""
},

{
    "location": "Rings/integer.html#Greatest-common-divisor-2",
    "page": "Integers",
    "title": "Greatest common divisor",
    "category": "section",
    "text": "gcd(a::Oscar.Integer, b::Oscar.Integer) -> Oscar.IntegerOne or both arguments may be Julia integers, however if they both are, the result will be a Julia integer not an Oscar integer.Returns the greatest common divisor of its inputs, which is by definition the largest integer dividing the two inputs, unless both inputs are zero in which case it returns zero. The result will always be non-negative and will only be zero if both inputs are zero.gcd(ZZ(34), ZZ(17))\ngcd(ZZ(3), ZZ(0))"
},

{
    "location": "Rings/integer.html#Extended-GCD-1",
    "page": "Integers",
    "title": "Extended GCD",
    "category": "section",
    "text": "gcdx(a::Oscar.Integer, b::Oscar.Integer) -> (Oscar.Integer, Oscar.Integer, Oscar.Integer)One or both arguments may be Julia integers, however if they both are, the result will be a tuple of Julia integers not Oscar integers.Returns a tuple (g s t) such that g is the greatest common divisor of a and b and g = as + bt. Normally s and t are chosen so that s  b(2g) and t  a(2g), where this uniquely defines s and t. But the following cases are handled specially:if a = b then t = bb\nif b = 0 or b = 2g then s = aa\nif a = 0 or a = 2g then t = bb"
},

{
    "location": "Rings/integer.html#Least-common-multiple-1",
    "page": "Integers",
    "title": "Least common multiple",
    "category": "section",
    "text": "lcm(a::Oscar.Integer, b::Oscar.Integer) -> Oscar.IntegerOne or both arguments may be Julia integers, however if they both are, the result will be a Julia integer not an Oscar integer.Returns the least common multiple of a and b. This is the least positive multiple of a and b, unless a = 0 or b = 0 which case we define the least common multiple to be zero.lcm(ZZ(6), ZZ(21))\nlcm(ZZ(0), ZZ(0))note: Note\nThe identity gcd(m n)mathrmlcm(m n) = mn does not hold for the definition that Oscar uses, unless both m and n are the same sign or one of them is zero."
},

{
    "location": "Rings/integer.html#Roots-1",
    "page": "Integers",
    "title": "Roots",
    "category": "section",
    "text": ""
},

{
    "location": "Rings/integer.html#Square-roots-1",
    "page": "Integers",
    "title": "Square roots",
    "category": "section",
    "text": "Julia and Oscar distinguish two kinds of square root:Integer square root (isqrt)\nFloating point square root (sqrt)We describe only the first of these here.isqrt(n::Oscar.Integer) -> Oscar.IntegerReturns the floor of the square root of its argument, i.e. the largest integer whose square does not exceed its input. An exception is raised if a negative input is passed.isqrt(ZZ(16))\nisqrt(ZZ(0))\nisqrt(ZZ(5))\nisqrt(ZZ(-3))isqrtrem(n::Oscar.Integer) -> (Oscar.Integer, Oscar.Integer)Returns the tuple (s, r) such that s is equal to isqrt(n) and n = s^2 + r.isqrtrem(ZZ(16))\nisqrtrem(ZZ(5))"
},

{
    "location": "Rings/integer.html#General-roots-1",
    "page": "Integers",
    "title": "General roots",
    "category": "section",
    "text": "root(a::Oscar.Integer, n::Base.Int) -> Oscar.IntegerReturns the value r of largest absolute value such that r^n leq a. When a is a perfect n-th power, the return value will be an n-th root of a.When n is even, the non-negative root is always returned. An exception is raised if n leq 0 or if n is even and a  0.root(ZZ(16), 4)\nroot(ZZ(5), 2)\nroot(ZZ(-5), 3)\nroot(ZZ(0), 4)\nroot(ZZ(-5), 2)\nroot(ZZ(12), -2)"
},

{
    "location": "Rings/integer.html#Conversions-1",
    "page": "Integers",
    "title": "Conversions",
    "category": "section",
    "text": "Int(n::Oscar.Integer) -> Base.Int\nBigInt(n::Oscar.Integer) -> Base.BigIntConvert the Oscar integer to the respective Julia integer.n = ZZ(123)\nInt(n)\nBigInt(n)In the case of Int, if the Oscar integer is too large to fit, an exception is raised.Int(ZZ(12348732648732648763274868732687324))fits(::Type{Int}, n::Oscar.Integer) -> BoolReturns true if the Oscar integer will fit in a Base.Int.fits(Int, ZZ(123))\nfits(Int, ZZ(12348732648732648763274868732687324))"
},

{
    "location": "Rings/integer.html#Factorisation-1",
    "page": "Integers",
    "title": "Factorisation",
    "category": "section",
    "text": "factor(n::Oscar.Integer) -> FacReturns a factorisation of the given integer. The return value is a special factorisation struct which can be manipulated using the functions below.factor(ZZ(-6000361807272228723606))\nfactor(ZZ(0))unit(F::Fac) -> Oscar.IntegerF = factor(ZZ(-12))\nunit(F)"
},

{
    "location": "Rings/integer.html#Factorisation-are-iterable-1",
    "page": "Integers",
    "title": "Factorisation are iterable",
    "category": "section",
    "text": "Once created, a factorisation is iterable:F = factor(ZZ(-60))\nfor (p, e) in F; println(\"$p^$e\"); endThe pairs (p, e) in a factorisation represent the prime power factors p^e of the non-unit part of the factorisation. They can be placed in an array using collect:F = factor(ZZ(-60))\ncollect(F)"
},

{
    "location": "Rings/integer.html#Accessing-exponents-in-a-factorisation-1",
    "page": "Integers",
    "title": "Accessing exponents in a factorisation",
    "category": "section",
    "text": "One can also determine whether a given prime is in the non-unit part of a factorisation and if so return its exponent. If the exponent of a prime that is not in a factorisation is requested, an exception is raised.For convenience, a Base.Int can be used instead of an Oscar integer for this functionality.F = factor(ZZ(-60))\n5 in F\nZZ(3) in F\n7 in F\nF[3]\nF[ZZ(7)]"
},

{
    "location": "Rings/integer.html#Combinatorial-functions-1",
    "page": "Integers",
    "title": "Combinatorial functions",
    "category": "section",
    "text": "note: Note\nThe functions in this section that take Base.Int arguments will return a Base.Int, which may overflow. Use the Oscar.Integer versions if this is not the desired behaviour."
},

{
    "location": "Rings/integer.html#Factorial-1",
    "page": "Integers",
    "title": "Factorial",
    "category": "section",
    "text": "factorial(n::Oscar.Integer) -> Oscar.IntegerReturns the factorial of n, i.e. n. An exception is raised if n  0. We define 0 = 1.note: Note\nThe function factorial is already defined in Julia for Base.Int.rising_factorial(x::Base.Int, n::Base.Int) -> Base.Int\nrising_factorial(x::Oscar.Integer, n::Base.Int) -> Oscar.Integer\nrising_factorial(x::Oscar.Integer, n::Oscar.Integer) -> Oscar.IntegerReturns x(x + 1)(x + 2)ldots(x + n - 1). An exception is raised if n  0. We define rising_factorial(x, 0) to be 1.factorial(ZZ(30))\nrising_factorial(ZZ(-30), 3)"
},

{
    "location": "Rings/integer.html#Primorial-1",
    "page": "Integers",
    "title": "Primorial",
    "category": "section",
    "text": "primorial(n::Base.Int) -> Base.Int\nprimorial(n::Oscar.Integer) -> Oscar.IntegerReturns the promorial P(n), i.e. the product of all primes less than or equal to n. An exception is raised if n  0. We define P(0) = P(1) = 1.primorial(ZZ(100))"
},

{
    "location": "Rings/integer.html#Bell-numbers-1",
    "page": "Integers",
    "title": "Bell numbers",
    "category": "section",
    "text": "bell(n::Base.Int) -> Base.Int\nbell(n::Oscar.Integer) -> Oscar.IntegerReturns the n-th Bell number B(n), i.e. the number of ways of partitioning a set of n elements. An exception is raised if n  0.bell(ZZ(20))"
},

{
    "location": "Rings/integer.html#Binomial-coefficients-1",
    "page": "Integers",
    "title": "Binomial coefficients",
    "category": "section",
    "text": "binomial(n::Oscar.Integer, k::Oscar.Integer) -> Oscar.IntegerReturns the binomial coefficient fracnk(n - k). If n k  0 or k  n we return zero.note: Note\nJulia already defines the binomial function for Base.Int.binomial(ZZ(72), ZZ(15))"
},

{
    "location": "Rings/integer.html#Integer-partitions-1",
    "page": "Integers",
    "title": "Integer partitions",
    "category": "section",
    "text": "number_of_partitions(n::Base.Int) -> Base.Int\nnumber_of_partitions(n::Oscar.Integer) -> Oscar.IntegerReturns the number of integer partitions p(n) of n, i.e. the number of distinct ways to write n as a sum of positive integers. Note that p(0) = 1, as the empty sum is counted. For n  0 we return zero.number_of_partitions(ZZ(10^6))"
},

{
    "location": "Rings/integer.html#Fibonacci-sequence-1",
    "page": "Integers",
    "title": "Fibonacci sequence",
    "category": "section",
    "text": "fibonacci(n::Base.Int) -> Base.Int\nfibonacci(n::Oscar.Integer) -> Oscar.IntegerReturns the n-th Fibonacci number F(n), defined by the recurrence relation F(1) = 1, F(2) = 1 and F(n) = F(n - 1) + F(n - 2) for n geq 3. We define F(0) = 0 and for n  0 we have F(-n) = (-1)^n+1F(n).fibonacci(ZZ(100))\nfibonacci(-2)"
},

{
    "location": "Rings/integer.html#Number-theoretic-functionality-1",
    "page": "Integers",
    "title": "Number theoretic functionality",
    "category": "section",
    "text": "note: Note\nThe functions in this section that take Base.Int arguments will return a Base.Int, which may overflow. Use the Oscar.Integer versions if this is not the desired behaviour."
},

{
    "location": "Rings/integer.html#Moebius-mu-function-1",
    "page": "Integers",
    "title": "Moebius mu function",
    "category": "section",
    "text": "moebius_mu(n::Base.Int) -> Base.Int\nmoebius_mu(n::Oscar.Integer) -> Base.Int Return the Moebius function mu(n), which is defined to be 0 if n is not squarefree and otherwise is defined to be +1 or -1 if n has an even or odd number of prime factors, respectively. Alternatively, mu(n) can be defined to be the sum of the primitive n-th roots of unity. An exception is raised if n leq 0.moebius_mu(30)"
},

{
    "location": "Rings/integer.html#Jacobi-symbols-1",
    "page": "Integers",
    "title": "Jacobi symbols",
    "category": "section",
    "text": "jacobi_symbol(m::Base.Int, n::Base.Int) -> Base.Int\njacobi_symbol(m::Oscar.Integer, n::Oscar.Integer) -> Base.IntReturn the Jacobi symbol left(fracmnright), which is defined for integers m and odd, positive integers n. If the factorisation of n is n = p_1^i_1p_2^i_2ldots p_r^i_r then we defineleft(fracmnright) = left(fracmp_1right)^i_1left(fracmp_2right)^i_2ldots left(fracmp_rright)^i_rwhere left(fracmpright) on the right hand side is the Legendre symbol, which is defined for an odd prime number p to be 0 if p divides m and otherwise +1 or -1 depending on whether m is a square modulo p or not. An exception is raised if n is even or if n leq 0.jacobi_symbol(3, 37)"
},

{
    "location": "Rings/integer.html#Sigma-function-1",
    "page": "Integers",
    "title": "Sigma function",
    "category": "section",
    "text": "divisor_sigma(m::Base.Int, n::Base.Int) -> Base.Int\ndivisor_sigma(m::Oscar.Integer, n::Base.Int) -> Oscar.Integer\ndivisor_sigma(m::Oscar.Integer, n::Oscar.Integer) -> Oscar.IntegerReturn the sum of the n-th powers of the divisors of msigma(m n) = sum_dm d^nIf m leq 0 or n  0 we raise an exception.divisor_sigma(60, 5)"
},

{
    "location": "Rings/integer.html#Euler-totient-function-1",
    "page": "Integers",
    "title": "Euler totient function",
    "category": "section",
    "text": "euler_phi(n::Base.Int) -> Base.Int\neuler_phi(n::Oscar.Integer) -> Oscar.IntegerReturn the Euler totient function varphi(n), i.e. the number of positive integers 1 leq x leq n which are coprime to n. Note that varphi(1) = 1. We raise an exception if n leq 0.euler_phi(200)"
},

{
    "location": "Rings/rational.html#",
    "page": "Rationals",
    "title": "Rationals",
    "category": "page",
    "text": "CurrentModule = Oscarusing OscarPages = [\"rational.md\"]"
},

{
    "location": "Rings/rational.html#Rationals-1",
    "page": "Rationals",
    "title": "Rationals",
    "category": "section",
    "text": "Fractions are created in Julia with the double slash operator //. If a fraction is created from Julia integers, a Julia fraction results and if either the numerator or denominator is an Oscar integer, an Oscar fraction results.Oscar has just a single rational fraction type Oscar.fmpq. The Oscar fraction type belongs to Oscar.Rational.Julia has a parameterised type Base.Rational{T} for its own fractions, where T is the integer type the fractions are constructed from, e.g. Base.Rational{Int} and Base.Rational{BigInt}.Note that conveniently, both of these Julia rational types belong to Base.Rational. Therefore in the description below we refer to Oscar.Rational for Oscar rational numbers and Base.Rational for any of the Julia rational number types.The situation is described by the following diagram.(Image: alt text)"
},

{
    "location": "Rings/rational.html#The-field-of-rationals-1",
    "page": "Rationals",
    "title": "The field of rationals",
    "category": "section",
    "text": "The parent of an Oscar rational number is the field of rationals. It can be constructed from the ring of integers ZZ using the FractionField constructor.For convenience, QQ is already defined to be the field of rational numbers.S = FractionField(ZZ)\nQQ"
},

{
    "location": "Rings/rational.html#Integer-constructors-1",
    "page": "Rationals",
    "title": "Integer constructors",
    "category": "section",
    "text": "Oscar rationals can be created using QQ. Two arguments can be passed to specify numerator and denominator. If a single argument is passed, the denominator is set to 1.For convenience, QQ also accepts Julia integers and rationals, but will always construct an Oscar rational.Naturally, Julia\'s double slash operator can also be used to construct fractions. However, unlike QQ, the double slash operator only constructs an Oscar rational if either the numerator or denominator is an Oscar integer.An exception is raised if a fraction is constructed with denominator zero.QQ(1, 2)\nQQ(5)\nZZ(3)//5\n1//ZZ(7)\nQQ(2//3)\nZZ(3)//0One can also construct the rational number 0 with the empty constructor:QQ()The following special constructors are also provided:zero(QQ)\none(QQ)zero(QQ)\none(QQ)"
},

{
    "location": "Rings/rational.html#Predicates-1",
    "page": "Rationals",
    "title": "Predicates",
    "category": "section",
    "text": "iszero(n::Oscar.Rational) -> Bool\nisone(n::Oscar.Rational) -> Bool\nisunit(n::Oscar.Rational) -> BoolThe isunit function will return true iff n neq 0.iszero(QQ())\nisone(one(QQ))\nisunit(QQ(-2, 3))"
},

{
    "location": "Rings/rational.html#Properties-1",
    "page": "Rationals",
    "title": "Properties",
    "category": "section",
    "text": "numerator(n::Oscar.Rational) -> Oscar.Integer\ndenominator(n::Oscar.Rational) -> Oscar.IntegerReturn the numerator and denominator respectively, of n.sign(n::Oscar.Rational) -> Oscar.RationalReturns the sign of n, i.e. nn if n neq 0, or 0 otherwise.sign(QQ(2, 3))\nsign(QQ())\nsign(QQ(-1))abs(n::Oscar.Rational) -> Oscar.RationalReturn the absolute value of n, i.e. n if n geq 0 and -n otherwise.abs(QQ(-3, 2))height(n::Oscar.Rational) -> Oscar.IntegerReturn the maximum of the absolute values of the numerator and denominator of n.height(QQ(324987329, -8372492324))floor(n::Oscar.Rational) -> Oscar.RationalReturn the greatest integer m such that m leq n.ceil(n::Oscar.Rational) -> Oscar.RationalReturn the least integer m such that m geq n.floor(QQ(-2, 3))\nceil(QQ(7, 2))\nceil(QQ(5))"
},

{
    "location": "Rings/rational.html#Basic-arithmetic-1",
    "page": "Rationals",
    "title": "Basic arithmetic",
    "category": "section",
    "text": "Oscar provides the basic arithmetic operations +, - and * and comparison operators ==, !=, <, <=, >, >=, including mixed operations between Julia and Oscar rationals and integers."
},

{
    "location": "Rings/rational.html#[Exact-Division]-1",
    "page": "Rationals",
    "title": "[Exact Division]",
    "category": "section",
    "text": "divexact(a::Oscar.Rational, b::Oscar.Rational) -> Oscar.Rational\ndivexact(a::Oscar.Rational, b::Oscar.Integer) -> Oscar.Rational\ndivexact(a::Oscar.Integer, b::Oscar.Rational) -> Oscar.RationalIn the first signature, one of the arguments may be a Julia rational and in the other two signatures the integers may be Julia integers.Returns the quotient of a by b. Exact division raises an exception if division by zero is attempted.divexact(QQ(2, 3), QQ(3, 5))\ndivexact(QQ(1, 3), ZZ(0))\ndivexact(QQ(3, 4), ZZ(5))\ndivexact(ZZ(6), QQ(2, 3))\ndivexact(QQ(1, 3), 5)"
},

{
    "location": "Rings/rational.html#Powering-1",
    "page": "Rationals",
    "title": "Powering",
    "category": "section",
    "text": "^(a::Oscar.Rational, b::Base.Int) -> Oscar.RationalPowering can be accomplished naturally using the special caret infix operator:QQ(5, 7)^32\nQQ(1, 2)^(-2)The following is allowed for convenience.QQ(0)^0note: Note\nIn Julia, the rational number 01 when raised to a negative power returns 10 to indicate that the value is undefined. Oscar raises an exception.QQ(0)^-2"
},

]}
