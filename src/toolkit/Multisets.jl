import Base: union, keys, haskey

struct Multiset{T}
	info::Dict{T,Integer}
end

Multiset(A::Array{T}) where {T} = Multiset( Dict{T,Integer}( k=>count(==(k),A) for k in unique(A) ) )

(M::Multiset{T})(k::T) where {T} = haskey(M.info, k) ? M.info[k]::Integer : 0::Integer

keys(M::Multiset) = keys(M.info)
haskey(M::Multiset) = haskey(M.info)

union(A::Array{Multiset{T}}) where {T} = Mulitset( Dict( k=>maximum(map(a->a(k), A)) for k in unique(keys.(Set(A))...) ))
union(A::Multiset{T}...) where {T} = Mulitset( Dict( k=>maximum(map(a->a(k), A)) for k in unique(keys.(A)...) ))

intersection(A::Array{Multiset{T}}) where {T} = Mulitset( Dict( k=>minimum(map(a->a(k), A)) for k in unique(keys.(Set(A))...) ))
intersection(A::Multiset{T}...) where {T} = Mulitset( Dict( k=>minimum(map(a->a(k), A)) for k in unique(keys.(A)...) ))
