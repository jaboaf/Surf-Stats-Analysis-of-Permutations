include("PanelAnalysis.jl")

# N is the number of variables
struct Monomial{D,R}
	coeff::R
	pwrs::NTuple{D,<:Integer}
end
(m::Monomial{D,R})(x::Vararg{R,D}) where {D,R<:Number} = m.coeff*prod(x.^ m.pwrs) 

struct Polynomial{D,R}
	parts::Vector{Monomial{D,R}}
end
(p::Polynomial{D,R})(x::Vararg{R,D}) where {D,R<:Number} = mapreduce(m->m(x),+,p.parts)

Monomial(I::CartesianIndex) = Monomial(1,ntuple(i->count(==(i),Tuple(I)),length(I)))
Monomial(c,I::CartesianIndex) = Monomial(c,ntuple(i->count(==(i),Tuple(I)),length(I)))
function Polynomial(T::Array{R,D}) where {R<:Number,D}
	d = ndims(T)
	supp = findall(!=(0),T)
	seen = NTuple[]
	basis = Dict()
	for I in supp
		term = ntuple(i->count(==(i),Tuple(I)),d)
		if term in seen
			basis[term] += T[I]
		else
			basis[term] = T[I]
			push!(seen,term)
		end
	end
	v = vec(map(t->Monomial(basis[t],t),seen))
	return v
end

function Polynomial(T::Array{R,D}) where {R<:Number,D}
	d = ndims(T)
	supp = findall(!=(0),T)
	seen = NTuple[]
	basis = Dict()
	for I in supp
		term = ntuple(i->count(==(i),Tuple(I)),d)
		if term in seen
			basis[term] += T[I]
		else
			basis[term] = T[I]
			push!(seen,term)
		end
	end
	v = vec(map(t->Monomial(basis[t],t),seen))
	return v
end

# Support of P
supp = findall(!=(0),P)
monomial(i::CartesianIndex) = [count(==(Int(c)), Tuple(i)) for c in JUD_ORIGS]
monomial(t::NTuple) = [ count(==(Int(c)),t) for c in JUD_ORIGS ]
monomials = unique(monomial.(supp))
basis = Dict{Array,Float64}(monomials .=> 0 )
for t in supp
	basis[monomial(t)] += P[t]
end
p(x::Vararg{<:Number,7}) = mapreduce(m->basis[m]*prod(x.^m),+,monomials)
# this is the probability generating function

for i in 1:7
	plot(LinRange(0,1,101),x->p([k==i ? x : 1 for k in 1:7]...) , label="$(ORIG(i))",hold=true)
end
plot(LinRange(0,1,101),x-> 7^-1 * sum([ p([i==k ? x : 1 for k in 1:7]...) for i in 1:7]),"x",markersize=.4, label="Avg. Marginal" ,hold=true)

