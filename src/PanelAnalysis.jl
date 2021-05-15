include("SimplifyComp.jl")

P = zeros(Float64, (7,7,7,7,7) );
Panel_λs = map(x->x[2].λ_c,WAVES)
for λ in Panel_λs
	c = prod( factorial.(length.(λ)) )
	for a in Base.product( Sym.(λ)...)
		P[ Int.(vcat(a...))... ] += 1/c
	end
end
#=
s = 0
D = zeros(Float64, (5,5,5,5,5,5,5) );
#U_(k) = 1/5^k * ones(Float64, ntuple(i->5,k));
for t in multiset_permutations([c for i in 1:4 for c in JUD_ORIGS],5)
	I = [ c in t ? findall(==(c),t) : 1:5 for c in JUD_ORIGS]
	c = prod(length.(I))
	D[I...] .+= P[Int.(t)...]/c
	s += P[Int.(t)...]
end
=#

N = 910
S = SymOp(P)
total = sum(P)
Y = Array{<:Number,5}[]
Hts = map(x->map(y->y.λ_origs,x),last.(partitionBy(:heatId)) );
for ht_λs in Hts
	Yᵢ = zeros(Float32, (7,7,7,7,7))
	for λ in ht_λs
		c = prod( factorial.(length.(λ)) )
		for a in Base.product( Sym.(λ)...)
			Yᵢ[ Int.(vcat(a...))... ] += 1/c
		end
	end
	push!(Y,Yᵢ)
end
# turn heats in to densities
Y = 