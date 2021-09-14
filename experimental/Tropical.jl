mutable struct TropicalRing{T} <: Ring
end

field(R::TropicalRing) = FlintQQ

tropical_ring(::typeof(max)) = TropicalRing{typeof(max)}()

tropical_ring(::typeof(min)) = TropicalRing{typeof(min)}()

mutable struct TropicalRingElem{T} <: RingElem
  isinf::Bool
  data::fmpq
  parent::TropicalRing{T}

  function TropicalRingElem(R::TropicalRing{T}, isinf::Bool) where {T}
    @assert isinf
    z = new{T}()
    z.isinf = true
    z.parent = R
    return z
  end

  function TropicalRingElem(R::TropicalRing{T}, x::RingElem) where {T}
    z = new{T}()
    z.isinf = false
    z.data = x
    z.parent = R
    return z
  end
end

function (R::TropicalRing{T})(u::TropicalRingElem{T}) where {T}
  @assert parent(u) === R
  return u
end

function (R::TropicalRing{T})(u::RingElem) where {T}
  v = field(R)(u)
  @assert parent(v) === field(R)
  return TropicalRingElem(R, v)
end

function (R::TropicalRing{T})(u::Integer) where {T}
  return TropicalRingElem(R, field(R)(u))
end

data(x::TropicalRingElem) = x.data

isinf(x::TropicalRingElem) = x.isinf

inf(R::TropicalRing) = TropicalRingElem(R, true)

(R::TropicalRing)() = zero(R)

parent(x::TropicalRingElem) = x.parent

iszero(x::TropicalRingElem) = isinf(x)

fun(x::TropicalRing{typeof(min)}) = min

fun(x::TropicalRing{typeof(max)}) = max

function Base.:(==)(x::TropicalRingElem, y::TropicalRingElem)
  (isinf(x) && isinf(y)) && return true
  ((isinf(x) && !isinf(y)) || (!isinf(x) && isinf(y))) && return false
  return data(x) == data(y)
end

zero(T::TropicalRing) = inf(T)

function Base.:(+)(x::TropicalRingElem{T}, y::TropicalRingElem{T}) where {T}
  if isinf(x)
    return deepcopy(y)
  else
    if isinf(y)
      return deepcopy(x)
    else
      return parent(x)(fun(parent(x))(data(x), data(y)))
    end
  end
end

AbstractAlgebra.expressify(x::TropicalRingElem; context = nothing) = isinf(x) ? "∞" : Expr(:call, "", expressify(data(x), context = context))

@enable_all_show_via_expressify TropicalRingElem

@enable_all_show_via_expressify TropicalRing

zero(x::TropicalRingElem) = zero(parent(x))

one(R::TropicalRing) = R(zero(field(R)))

function Base.:(*)(x::TropicalRingElem{T}, y::TropicalRingElem{T}) where {T}
  if isinf(x)
    return x
  else
    if isinf(y)
      return y
    else
      return parent(x)(data(x) + data(y))
    end
  end
end

#function Base.:(-)(x::TropicalRingElem)
#  error("Computer at says no!")
#end

function Base.deepcopy_internal(x::TropicalRingElem, dict::IdDict)
  if !isinf(x)
    return TropicalRingElem(x.parent, Base.deepcopy_internal(data(x), dict))
  else
    return inf(parent(x))
  end
end

Oscar.elem_type(::Type{TropicalRing{T}}) where {T} = TropicalRingElem{T}

Oscar.parent_type(::Type{TropicalRingElem{T}}) where {T} = TropicalRing{T}

function Base.:(-)(x::TropicalRingElem, y::TropicalRingElem...)
  error("Computer says no!")
end

Oscar.mul!(x::TropicalRingElem, y::TropicalRingElem, z::TropicalRingElem) = y * z

Oscar.addeq!(y::TropicalRingElem, z::TropicalRingElem) = y + z

function Base.:(^)(a::TropicalRingElem, n::Int)
  return Base.power_by_squaring(a, n)
end

Base.copy(a::TropicalRingElem) = a

function permanent(x)
  R = base_ring(x)
  S = AbstractAlgebra.SymmetricGroup(nrows(x))
  res = zero(R)
  for s in S
    o = one(R)
    for i in 1:nrows(x)
      o = o * x[i, s[i]]
    end
    res = res + o
  end
  return res
end

#function show_via_expressify(io::IO, mi::MIME, @nospecialize(obj); context = nothing)
#end

#function Base.show(io::IO, mi::MIME"text/plain", p::AbstractAlgebra.Generic.Poly{TropicalRingElem{T}}) where {T}
#   AbstractAlgebra.show_obj(io, mi, AbstractAlgebra.expressify(p))
#end
#
#function Base.show(io::IO, p::AbstractAlgebra.Generic.Poly{TropicalRingElem{T}}) where {T}
#  AbstractAlgebra.show_obj(io, MIME("text/plain"), AbstractAlgebra.expressify(p))
#end

AbstractAlgebra.expressify(R::TropicalRing{typeof(min)}; context = nothing) = "Tropical ring (min)"

AbstractAlgebra.expressify(R::TropicalRing{typeof(max)}; context = nothing) = "Tropical ring (max)"

function AbstractAlgebra.expressify(@nospecialize(a::PolyElem{<:TropicalRingElem}),
   x = var(parent(a)); context = nothing)
   if iszero(a)
     return expressify(zero(base_ring(a)), context = context)
   end
   sum = Expr(:call, :+)
   for k in degree(a):-1:0
      c = coeff(a, k)
      if !iszero(c)
        xk = k < 1 ? expressify(one(base_ring(a)), context = context) : k == 1 ? x : Expr(:call, :^, x, k)
         if isone(c)
            push!(sum.args, Expr(:call, :*, xk))
         else
            push!(sum.args, Expr(:call, :*, expressify(c, context = context), xk))
         end
      end
   end
   return sum
end

function AbstractAlgebra.expressify(a::MPolyElem{<:TropicalRingElem}, x = symbols(parent(a)); context = nothing)
   if iszero(a)
     return expressify(zero(base_ring(a)), context = context)
   end
   sum = Expr(:call, :+)
   n = nvars(parent(a))
   for (c, v) in zip(coefficients(a), exponent_vectors(a))
      prod = Expr(:call, :*)
      if !isone(c)
         push!(prod.args, expressify(c, context = context))
      end
      for i in 1:n
         if v[i] > 1
            push!(prod.args, Expr(:call, :^, x[i], v[i]))
         elseif v[i] == 1
            push!(prod.args, x[i])
         end
      end
      # Capture empty products
      if length(prod.args) == 1
        prod = expressify(one(base_ring(a)), context = context)
      end
      push!(sum.args, prod)
   end
   return sum
end

one(R::AbstractAlgebra.Generic.MPolyRing{<:TropicalRingElem}) = R(one(base_ring(R)))

zero(R::AbstractAlgebra.Generic.MPolyRing{<:TropicalRingElem}) = R(zero(base_ring(R)))

one(R::AbstractAlgebra.Generic.PolyRing{<:TropicalRingElem}) = R(one(base_ring(R)))

zero(R::AbstractAlgebra.Generic.PolyRing{TropicalRingElem{S}}) where {S} = R(zero(base_ring(R)))
