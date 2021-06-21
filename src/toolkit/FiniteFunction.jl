import Base.∘

struct FiniteFn{DomainT,RangeT} 
	D::Array{DomainT}
	R::Array{RangeT}
	FiniteFn(D::Array{DomainT},R::Array{RangeT}) where {DomainT,RangeT} = allunique(D) ? new{DomainT,RangeT}(D,R) : new{DomainT,RangeT}(sort(unique(D)),[  R[findall(==(x),D)] for x in sort(unique(D))])
end

FiniteFn(P::Array{Pair{U,V}}) where {U,V} = FiniteFn(first.(P),last.(P))

(f::FiniteFn{DomT,RanT})(x::DomT) where {DomT,RanT} = f.R[findfirst(==(x),f.D)]
(f::FiniteFn{DomT,RanT})(X::DomT...) where {DomT,RanT} = f.(X)

∘(f::FiniteFn{V,W},g::FiniteFn{U,V}) where {U,V,W} = FiniteFn(g.D,f.(g.R))

