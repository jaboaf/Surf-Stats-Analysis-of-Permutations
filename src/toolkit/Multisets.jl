import Base: union, keys, haskey, ==

struct Multiset{T}
	info::Dict{T,<:Integer}
end

Multiset(A::Array{T}) where {T} = Multiset( Dict{T,Integer}( k=>count(==(k),A) for k in unique(A) ) )

(M::Multiset{T})(k::T) where {T} = haskey(M.info, k) ? M.info[k]::Integer : 0::Integer

keys(M::Multiset) = keys(M.info)
haskey(M::Multiset) = haskey(M.info)

union(A::Multiset{T}...) where {T} = Multiset( Dict( k=>maximum(map(a->a(k), A)) for k in collect(union(keys.(A)...)) ))
intersect(A::Multiset{T}...) where {T} = Multiset( Dict( k=>minimum(map(a->a(k), A)) for k in collect(union(keys.(A)...)) ))

==(A::Multiset{T},B::Multiset{T}) where {T} = (A.info == B.info)
