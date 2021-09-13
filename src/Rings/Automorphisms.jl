import Hecke: NfMorSet, HeckeMap, NumFieldMor, morphism_type

mutable struct PermGroupToNfMorSet{S, T} <: Map{PermGroup, NfMorSet{T}, HeckeMap, PermGroupToNfMorSet{S, T}}
  G::PermGroup
  images::Dict{PermGroupElem, S}
  preimages::Dict{S, PermGroupElem}
  header::MapHeader{PermGroup, NfMorSet{T}}

  function PermGroupToNfMorSet(G::PermGroup, images::Dict{PermGroupElem, S}, preimages::Dict{S, PermGroupElem}, s::NfMorSet{T}) where {S, T}
    z = new{S, T}()
    z.header = MapHeader(G, s)
    z.images = images
    z.preimages = preimages
    z.G = G
    return z
  end
end

function PermGroupToNfMorSet(G, K, p::Vector{PermGroupElem}, aut::Vector{S}) where S <: NumFieldMor
  return PermGroupToNfMorSet(G, Dict(q => a for (q, a) in zip(p, aut)),
                                Dict(a => q for (q, a) in zip(p, aut)),
                                NfMorSet(K))
end

#function image(f::PermGroupToNfMorSet, g::PermGroupElem)
#  @assert parent(g) == f.G
#  K = codomain(f).field
#  return f.aut[g[]]
#end
#
#function (f::PermGroupToNfMorSet)(g::PermGroupElem)
#  return image(f, g)
#end
#
#function preimage(f::PermGroupToNfMorSet{S, T}, a::S) where {S, T}
#  K = codomain(f).field
#  aut = automorphisms(K, copy = false)
#  for i in 1:length(aut)
#    if a == aut[i]
#      return domain(f)[i]
#    end
#  end
#  error("something wrong")
#end
#
#
#function evaluate(f::fmpq_poly, a::nf_elem)
#  #Base.show_backtrace(stdout, Base.stacktrace())
#  R = parent(a)
#  if iszero(f)
#    return zero(R)
#  end
#  if a == gen(R) && parent(f) == parent(parent(a).pol)
#    return R(f)
#  end
#  l = length(f) - 1
#  s = R(coeff(f, l))
#  for i in l-1:-1:0
#    #s = s*a + R(coeff(f, i))
#    mul!(s, s, a)
#    add!(s, s, coeff(f, i))
#  end
#  return s
#end
#
#Base.copy(f::NfToNfMor) = f
function automorphism_group(k::NumField)

  G, mG = Hecke.automorphism_group(k)
  H = symmetric_group(degree(k))
  gens(G) #to make sure gens are actually there...
  H = sub(H, [H(G.mult_table[i, :]) for i=G.gens])[1]
  D = Dict{PermGroupElem, morphism_type(k)}()
  DD = Dict{morphism_type(k), PermGroupElem}()

  for g in H
    m = Int[i^p for i in 1:degree(k)]
    i = Base.findfirst(x -> G.mult_table[x, :] == m, 1:degree(k))
    gg = GrpGenElem(G, i)
    a = mG(gg)
    D[gg] = a
    DD[a] = gg
  end

  return H, PermGroupNfToMorSet(H, D, DD, NfMorSet(k))
end


