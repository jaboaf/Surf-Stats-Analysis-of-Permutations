include("SimplifyComp.jl")
using StatsBase

Panel_λs = map(x->x[2].λ_c,WAVES)
# NOTE: this is w.r.t <
P = zeros(Float64, (7,7,7,7,7) );
for λ in Panel_λs
	c = prod( factorial.(length.(λ)) )
	for a in Base.product( Sym.(λ)...)
		P[ Int.(vcat(a...))... ] += 1/c
	end
end

ℙ = P/sum(P);

# NOTE: this is w.r.t >
R = permutedims(P, [5,4,3,2,1]);
ℝ = R / sum(R);


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

# 1) Chi-sq Test of for diff from unordered partitions
N = length(Panel_λs)
Y = map(λ->[count(==(i),length.(λ)) for i in 1:5], Panel_λs)
λ_Counts = countmap(Y)
λ_Obs = Dict([x=>λ_Counts[x]/N for x in keys(λ_Counts)])
λ_Thry = Dict(
	[K=> prod( [ 1//( factorial(K[j]) * j^K[j] ) for j in 1:5]) 
	for K in keys(λ_Counts) ]
)
χ_sq = length(keys(λ_Thry))*sum(
	[ (λ_Obs[x]-λ_Thry[x])^2 / λ_Thry[x]
	for x in keys(λ_Thry)]
)
println(χ_sq)

# 2) Panels taken under equivalence relation and ordered by min-max rule
Panel_eqλs = map(x->x[2].λ_eqc,WAVES)
Y = map(eqPanel -> [count(==(i),length.(eqPanel)) for i in 1:5], Panel_eqλs )
λeq_Counts = countmap(Y)
λeq_Props = Dict([x=>λeq_Counts[x]/N for x in keys(λeq_Counts)])
λeq_Thry = Dict(
	[K=> prod( [ 1//( factorial(K[j]) *  j^K[j] ) for j in 1:5]) 
	for K in keys(λeq_Counts) ]
)
eqχ_sq = length(keys(λeq_Props))*sum(
	[ (λeq_Props[x]-λeq_Thry[x])^2 / λeq_Thry[x]
	for x in keys(λeq_Props)]
)

# 3) 
Panel_eqλs_nom = map(λ-> unique.(λ),Panel_eqλs)
Y = map(eqλ_nom -> [count(==(i),length.(eqλ_nom)) for i in 1:5],Panel_eqλs_nom)
nOrigs_pd = [ count(==(i), map(x->sum(length.(x)),Panel_eqλs_nom))/N for i in 1:5]
λeq_noM_Counts = countmap(Y)
λeq_noM_Props = Dict([x=>λeq_noM_Counts[x]/N for x in keys(λeq_noM_Counts)])
λ_Thry_Cond_nOrigs_pd = Dict(
	[ K => nOrigs_pd[sum(K)] / prod([1//(factorial(K[j]) * j^K[j]) for j in 1:sum(K)]) 
	for K in keys(λeq_noM_Props) ]
)
eq_noM_cond_χ_sq = sum(
	[ (λeq_noM_Props[x]-λ_Thry_Cond_nOrigs_pd[x])/λ_Thry_Cond_nOrigs_pd[x]
	for x in keys(λ_Thry_Cond_nOrigs_pd)]
)


