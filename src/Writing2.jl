#' Some Preliminaries:
using GRUtils # for visualization
using Statistics # for general purpose stats things
using StatsBase

# Currently whats being used?
include("SimplifyComp.jl")
include("toolkit/SymGrpAndReps.jl")

#' Here we read in the raw data. I have elected to restrict our attention to 2018 and 2019 because those are the most complete years. See Table:
#' |Year | Waves | Waves w/ 0 Judge Origins Listed | Waves w/ 5 Judge Origins Listed | Total Number of Judge Origins | Waves w/ 3 Sub Scores Listed | Waves w/ 5 Sub Scores Listed | Total number of Judge Scores |
#' | -----|-----|-----|-----|-----|-----|-----|----- |
#' | 2017 | 7328 | 5210 | 2118 | 10590 | 299 | 7029 | 36042 |
#' | 2018 | 6639 | 336 | 6303 | 31515 | 0 | 6639 | 33195 |
#' | 2019 | 7648 | 79 | 7569 | 37845 | 0 | 7648 | 38240 |
#' | -----|-----|-----|-----|-----|-----|-----|----- |
#' | All | 21615 | 5625 | 15990 | 79950 | 299 | 21316 | 107477 |

#' ## Motivation for A Different Approach


#' # The Simple Approach
#' Our goal is to determine if judges have a tendency to give higher scores to surfers that share their same nationality.

#' The straight forward approach is to test the differences in means between judges with the same nationality as the surfer and those with a different nationality. I.e. Test
#' H₀: Mean(Match Scores) - Mean(No Match Scores) = 0
#' H₁: Mean(Match Scores) - Mean(No Match Scores) != 0
HasMatch = filter(x->x.I_match==true, last.(WAVES) )
DiffInMeans = map(x->mean(x.labeledPanelBinary[:Match]) - mean(x.labeledPanelBinary[:NoMatch]), HasMatch)
ttest(X) = mean(X) / (std(X)/sqrt(length(X)))
println( ttest(DiffInMeans) )
#' The difference of mean(DiffInMeans) is significant. However, the distribution of the differences in means is less compelling. 
savefig("visuals/DistOfDiffInMeans.jl", histogram(DiffInMeans,title="Distribution of Differences in Means") )
#' !(visuals/DistOfDiffInMeans.png)

#' We may carry out this process for each country:
#' Note C = :AUS, :BRA, :FRA, :PRT, :USA, :ZAF
B = [:Match, :NoMatch]
C = sort(unique(map(x->x.athOrig, HasMatch)))
function diffInMeansDist(c::ORIG)
	I = findall(x->x.athOrig==c, HasMatch)
	t = ttest( DiffInMeans[I] )
	savefig(
		"visuals/DistOfDiffInMeansFor$(c).jl",
		histogram(DiffInMeans[I],title="Difference in Means where C=$c (n=$(length(I)) and t=$t)")
	)
	return t
end
for c in C diffInMeansDist(c) end

#' And we find that differences in means between Matching and Non-matching judges, conditional on a Nationality is signifigant for AUS, BRA, FRA, USA, ZAF.
#' Though we will use the term "matching judge(s)" throughout the paper, it is important to keep in mind that it is the Athlete's origin which determines whether a judge is a "matching judge" or a "non-matching judge".

#' ## Judging Panel Data
#' Question: What is the definition of a panel? What type of data is it?
#' For any given heat, there are 5 judges on the judging panel. Anytime a surfer rides a wave, each of them observes the way 
#' What does a panel look like? A lot of very different things that may seem very similar BUT ARE NOT. For example:
Panels = map(x->x.panel,last.(WAVES))




display(Panels[1])
display(Panels[3])
display(Panels[5])


#' # Methods
#' ## Lets explore the Data!
#' We have lots of missing Judge Origins from panels in the 2017 World Surf League season so we will omit the 2017 season ... for now (This begs an intersting question which we should return to later).
#' 

#' We have constructed a multidimensional array, aka an m-way, cross classified, contingency table. We have m classification factors:
#' - WSL Season
#' - Event
#' - Round
#' - Heat
#' - Ahtlete Origin
#' - Judge Origin
#' - Size of Partition of Panel (Max Rank)
#' - Rank of Judge
println("m = $(ndims(rnkData))")
println(size(rnkData))

function marginalBarPlot(M::Array{T,N}, d::Integer) where {T,N}
	dIndex_(h::Integer) = [ i==d ? h : Colon() for i in 1:N]
	marginal_d = [ sum( M[ dIndex_(h)... ] ) for h in 1:size(M)[d] ]
	savefig(
		"visuals/marginal_$(rnkDataVarNames[d]).png",
		barplot(String.(Symbol.(rnkDataVarRngs[d])), marginal_d )
	)
	return marginal_d
end

for d in length(rnkDataVarNames) marginalBarPlot(rnkData, d) end

#' ![Judge Origin Marginal](visuals/marginals_YR.png)
#' ![Judge Origin Marginal](visuals/marginals_EVT.png)
#' ![Judge Origin Marginal](visuals/marginals_RND.png)
#' ![Judge Origin Marginal](visuals/marginals_HEAT.png)
#' ![Judge Origin Marginal](visuals/marginals_ATH_orig.png)
#' ![Judge Origin Marginal](visuals/marginals_JUD_orig.png)
#' ![Judge Origin Marginal](visuals/marginals_MAX_RANK.png)
#' ![Judge Rank Marginal](visuals/marginals_RANK.png)


#' # Analysis
#' One thing I could change:
#' - instead giving every element of rth partition rank  give them each:
#' { Nᵣ + i for i in 1:n_i}  where Nᵣ is Σʳnᵢ
#' ... but then any answer we give loses meaning
#' 
#' I wonder: Mᵣ := max order for that wave
#' B = 1_{order(Jᵢ) = 1 }
#' T = 1_{order(Jᵢ) = Mᵣ }
#' B ∪ T only really meaninful for Mᵣ>2
#' 
#' Wondering₁: if P(T | JUD_orig==ATH_orig ) =  P( T | JUD_orig != ATH_orig )
#' Wondering₂: if Wondering₁ depends on [   ]_orig
#' Wondering₃: if P( B ∪ T | Mᵣ ) = 2/Mᵣ
#' Wondering₄: if P(T | (JUD_orig==ATH_orig,Mᵣ) ) = 1/Mᵣ

#' ## Method 1
MATCH_ORIGS = sort(unique(map(x->x.athOrig,HasMatch)))
GivenMatchByMᵣ = [
sum( [sum(rnkData[:,:,:,:,c,c,r,:]) for c in Int.(MATCH_ORIGS) ] )
for r in 1:5
]
println("Match conditional on Max Rank is:")
GivenMatchByMᵣ

TopGivenMatchByMᵣ = [
sum( [ sum(rnkData[:,:,:,:,c,c,r,1]) for c in Int.(MATCH_ORIGS)])
for r in 1:5
]
println("Judge has Max Rank given Match (by Max Rank)")
TopGivenMatchByMᵣ

D₁ = TopGivenMatchByMᵣ ./ GivenMatchByMᵣ

#' We would expect:
E₁ = [1/r for r in 1:5]

#' So we have:
χsq₁ = sum(GivenMatchByMᵣ)*sum( (D₁ .- E₁).^2 ./ E₁ )

#' Now...
GivenMatchByMᵣandOrig = [ 
sum( rnkData[:,:,:,:,c,c,r,:] )
for r in 1:5, c in Int.(MATCH_ORIGS)
]
println("Match conditional on Max Rank (rows) and Nationality (cols)")
display(GivenMatchByMᵣandOrig)

TopGivenMatchByMᵣandOrig = [
sum( rnkData[:,:,:,:,c,c,r,1] )
for r in 1:5, c in Int.(MATCH_ORIGS)
]
println("Judge has Max Rank given Match (by Max Rank (rows) by Nationality (cols))")
println(TopGivenMatchByMᵣandOrig)

println(MATCH_ORIGS)
D₂ = TopGivenMatchByMᵣandOrig ./ GivenMatchByMᵣandOrig

#' We would expect:
E₂ = [ 1/r for r in 1:5, c in Int.(MATCH_ORIGS) ]

#' So we have a total χ^2 of:
χsq₂ = sum(GivenMatchByMᵣandOrig)*sum( (D₂ .- E₂).^2 ./ E₂)

W = (D₂ .- E₂).^2 ./ E₂
χsqbyctry = [ sum(GivenMatchByMᵣandOrig[:,c])*sum(W[:,c]) for c in 1:6]

MᵣGivenMatch =[ sum(rnkData[:,:,:,:,c,c,3:5,:] ) for c in Int.(MATCH_ORIGS) ]
OrderIsMᵣ = [ sum(rnkData[:,:,:,:,c,c,3:5,1]) for c in Int.(MATCH_ORIGS) ]
println(OrderIsMᵣ ./ MᵣGivenMatch)

#' But!!!

#' ``
#' |\{ g ∈ S_d | cycles(g) = 1^{k₁},2^{k₂},…,d^{k_d} \}| = \frac{d!}{∏_{j=1}^{d} k_j!j^k_j }
#' ⟹ P(1^k_1, … , d^k_d) = \frac{\frac{d!}{∏_{j=1}^{d} k_j!j^k_j }}{d!} = \frac{1}{∏_{j=1}^{d} k_j!j^k_j }
#' ``
#' Panels is the Array of the observed Panels, each of which is ordered by score.
#' We do not have a total order for every panel
#' So a observed panel is an ordered partition.
#' Y is an array of cycle counts.
#' Below, λ, is the cycle counts.
Panel_λs = map(x->x.λ_origs,last.(WAVES))
N = length(WAVES)
Y = map(panel -> [count(==(i),length.(panel)) for i in 1:5], Panel_λs)
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

#' ... χ² is very large....
eqPanels = map(x->x.eqPanel, last.(WAVES))
eqParts = map(eqPanel -> [count(==(i),length.(last.(eqPanel))) for i in 1:5], eqPanels )
eqParts_Counts = countmap(eqParts)
eqParts_Props = Dict([x=>eqParts_Counts[x]/N for x in keys(eqParts_Counts)])
Theory_Parts_Props = Dict(
	[K=> prod( [ 1//( factorial(K[j]) *  j^K[j] ) for j in 1:5]) 
	for K in keys(eqParts_Counts) ]
)
eqχ_sq = length(keys(eqParts_Props))*sum(
	[ (eqParts_Props[x]-Theory_Parts_Props[x])^2 / Theory_Parts_Props[x]
	for x in keys(eqParts_Props)]
)

#' Now eqPanels with No mulitplicity
Y = map(eqPanels) do panel
	x = unique.(last.(panel))
	return [ count( ==(i), length.(x) ) for i in 1:5]
end
nOrigs_Counts = countmap( map(x->sum([i*x[i] for i in 1:5]), Y))
nOrigs_Props = Dict([x=>nOrigs_Counts[x]/N for x in keys(nOrigs_Counts)])
λ_Counts = countmap(Y)
λ_Props = Dict([x=>λ_Counts[x]/N for x in keys(λ_Counts)])
λ_Thry_Cond_Origs = Dict(
	[ K => nOrigs_Props[sum([i*K[i] for i in 1:5])]/N * prod( [ 1//( factorial(K[j]) *  j^K[j] ) for j in 1:5]) 
	for K in keys(λ_Props) ]
)
eqNoM_cond_χ_sq = length(keys(λ_Thry_Cond_Origs))*sum(
	[ (λ_Props[x]-λ_Thry_Cond_Origs[x])^2 / λ_Thry_Cond_Origs[x]
	for x in keys(λ_Thry_Cond_Origs)]
)

#' what marty thought of
D = Dict([wave[1] => wave[2].eqPanel for wave in WAVES])
χ_sq_byHT = []
χ_sq_byHT_nom = []
for heat in partitionBy(:heatId)
	Panels = [ x.λ_origs for x in heat[2] ]
	Panels_nom = [ unique.(λ) for λ in Panels ]
	n = length(Panels)
	λ_HT_Origs_Counts = countmap( map(x->sum(length.(x)), Panels_nom))
	λ_HT_nOrigs_Props = Dict([x=>λ_HT_Origs_Counts[x]/n for x in keys(λ_HT_Origs_Counts)])
	λ_HT = map(x -> [count(==(i),length.(x)) for i in 1:5], Panels )
	λ_HT_nom = map(x -> [count(==(i),length.(x)) for i in 1:5], Panels_nom )
	λ_HT_Counts = countmap(λ_HT)
	λ_HT_Counts_nom = countmap(λ_HT_nom)
	λ_HT_Props = Dict([x=>λ_HT_Counts[x]/n for x in keys(λ_HT_Counts)])
	λ_HT_Props_nom = Dict([x=>λ_HT_Counts_nom[x]/n for x in keys(λ_HT_Counts_nom)])
	λ_HT_Thry = Dict(
		[K=> prod( [ 1//( factorial(K[j]) *  j^K[j] ) for j in 1:5]) 
		for K in keys(λ_HT_Counts) ]
	)
	λ_HT_Thry_Cond_Origs = Dict(
		[ K => λ_HT_nOrigs_Props[sum([i*K[i] for i in 1:5])]/n * prod( [ 1//( factorial(K[j]) *  j^K[j] ) for j in 1:5]) 
		for K in keys(λ_HT_Props_nom) ]
	)
	eqχ_sq = length(keys(λ_HT_Thry))*sum(
		[ (λ_HT_Props[x]-λ_HT_Thry[x])^2 / λ_HT_Thry[x]
		for x in keys(λ_HT_Thry)]
	)
	eqχ_sq_cond_nom = length(keys(λ_HT_Thry_Cond_Origs))*sum(
	[ (λ_HT_Props_nom[x]-λ_HT_Thry_Cond_Origs[x])^2 / λ_HT_Thry_Cond_Origs[x]
	for x in keys(λ_HT_Thry_Cond_Origs)]
	)
	push!(χ_sq_byHT, eqχ_sq)
	push!(χ_sq_byHT_nom, eqχ_sq_cond_nom)
end

savefig("visuals/Partition_χ_sq_byHT_Hist.png", histogram(χ_sq_byHT))
savefig("visuals/Partition_nom_χ_sq_byHT_Hist.png", histogram(χ_sq_byHT_nom))

#' ![Judge Origin Marginal](visuals/Partition_χ_sq_byHT_Hist.png)
#' ![Judge Origin Marginal](visuals/Partition_nom_χ_sq_byHT_Hist.png)

#' The most granular object of study is the panel. It would not be justified to analyze our data a sequence of 13,872 panels because during each heat the set of judges is fixed (or at least it appears to be this way).
all( all(ht[2][i-1].panel_origs == ht[2][i].panel_origs for i in 2:length(ht[2])) for ht in partitionBy(:heatId) ) 
#' For this reason we will be mostly interested in heat-level dynamics of the panel. This is similar to the game-level approach Price takes in analyzing NBA refereeing (as opposed to simply many foul calls).
#' A panel is comprised of a set of distinct humans (Judges), J = {j₁, j₂, ..., j₅}, and their origns by the multiset ``J_c:C→ℕ``.
#' Multisets are interesting. Roman defines a multiset on a set A, as an element of A×N. Blizzard gives in broad overview of different approaches to multisets in [The Development of Multiset Theory], and takes a more normative approach in [Dedekind Multisets and Function Shells]. We will not delve into this, but they are important to understand in the construction below:
#' A Defn: A multiset on a set A, is a function f:A→ℕ.
#' Rmk: Some authors use the pair (A,f) to define a multiset. While it may bring us comfort to know the underlying set which the function is defined on, this is merely a luxury of defining a multiset. In an observational setting one does not nessisary know the underlying set, only the support of f, ie. the set of all elements mapped to a non-zero number. For example, the judges in a heat are selected from a pool of judges provided by the WSL for that event. Specifically, suppose it is march ___ 2018, Caio Ibelli rides the first wave of the WSL 2018 Mens Championship Tour, when all 5 judges provide their scores, Caio Ibelli is awarded a score of:
round(mean(sort(WAVES[1][2].judge_scores)[2:4]),digits=2)
#' Not a dramatic first wave eh? But maybe you want to know a bit more about the sub scores that resulted in a wave score of 0.43 . The panel is as follows:
WAVES[1][2].panel
#' We can certainly deduce the nationalities of 5 of the 8 Judges in the event pool.
WAVES[1][2].panel_origs
#' Where are the remaining 3 judges from? Well, we have to wait for the next heat because these judges are the panel for this heat.
#' In the next heat, Michael Rodrigues is the first to ride a wave, receiving a 0.3.
round(mean(sort(partitionBy(:heatId)[2][2][1].judge_scores)[2:4]),digits=2)
#' This feels similar to the first wave of last heat, are the judges the same?
partitionBy(:heatId)[2][2][1].panel
#' Looks like No. Okay, well are some of them the judges the same? Certainly. We could have answered that even before the Caio Ibelli took the first wave of 2018 because if there were two panels, each with five judges, and no judges were the same, then there would be at least 10 total judges for the event, a breach of WSL rulebook. So, which of the judges from heat 1 are on the panel from round 2?
#' we know there are at least:
println([AUS => 2, BRA => 2, ESP => 1, USA => 1, ZAF => 1])
#' If we assumed the panels were disjoint, we could deduce that the pool of judges for the event, at least, consisted of:
println([AUS => 3, BRA => 3, ESP => 1, USA => 2, ZAF => 1])
#` However, our information from the rulebook tells us that this certainly cannot be the case.

partitionBy(:heatId)[2][2][1].panel_origs
partitionBy(:heatId)[1][2][1].panel_origs

# are in you observe
partitionBy(:heatId)[1][2][1].panel
partitionBy(:heatId)[2][2][1].panel_origs

#' some utilities that will be of use
⊗(A::Array{T},B::Array{T}) where T<: Number = prod.(Base.product(A,B))
⊗(a::NTuple{T},b::NTuple{T}) where T  = (a...,b...)

×(A::Set,B::Set) = Set(Base.product(A,B))
×(A::Array,B::Array) = collect(Base.product(A,B))

#=
E(X::Array,i)=dropdims( sum(X,dims=setdiff(1:ndims(X),i)),dims=setdiff(1:ndims(X),i))
E(X::Array,I::NTuple)=dropdims(sum(X,dims=i),setdiff(1:ndims(X),I))
cov(X::Array,i::K,j::K) where {K<:Integer} = E(X,(i,j))-E(X,i)⊗E(X,j)
cov(X::Array,I::NTuple{K},j::K) where K<:Integer = E(X,(I...,j))-E(X,I)⊗E(X,j)
cov(X::Array,i::K,J::NTuple{K}) where K<:Integer = E(X,(i,J...))-E(X,i)⊗E(X,J)
cov(X::Array,I::NTuple{K},J::NTuple{k}) where K<:Integer =E(X,(I...,J...))-E(I)⊗E(J)
=#

# this is succinct Julia for create a multidimensional array with 
# 		Prob(arrangement of judges) = D[arrangement of judges]/sum(D)
D = zeros(Float64, (7,7,7,7,7) );

for P in Panel_λs
	c = prod( factorial.(length.(P)) )
	for a in Base.product( Sym.(P)...)
		D[ Int.(vcat(a...))... ] += 1/c
	end
end

#CHK
sum(D)
length(Panels)

#' D is currently the entire data set, This isn't partitularly accurate of a view of the data. Instead lets chop it up.

Y = Array{<:Number,5}[]
H = map(x->map(y->y.λ_origs,x),last.(partitionBy(:heatId)) );
for ht_λs in H
	Yᵢ = zeros(Float32, (7,7,7,7,7))
	for λ in ht_λs
		c = prod( factorial.(length.(λ)) )
		for a in Base.product( Sym.(λ)...)
			Yᵢ[ Int.(vcat(a...))... ] += 1/c
		end
	end
	push!(Y,Yᵢ)
end

#' First note that since there are varying numbers of waves in a heat so, the sums of Yᵢ are distributed as #waves in heat.
histogram(sum.(Y))

#'
μ = sum(Y)/ length(Y) ;
AbsErr = map(y-> abs.(y - μ), Y);
SqErr = map(y-> y.^2,AbsErr);
SqrtErr = map(y-> sqrt.(y),AbsErr);
histogram( sum.(AbsErr) )


Y = Array{<:Number,5}[]
Hts = map(x->map(y->y.λ_origs,x),last.(partitionBy(:heatId)) );
for ht_λs in H
	Yᵢ = zeros(Float32, (7,7,7,7,7))
	for λ in ht_λs
		c = prod( factorial.(length.(λ)) )
		for a in Base.product( Sym.(λ)...)
			Yᵢ[ Int.(vcat(a...))... ] += 1/c
		end
	end
	push!(Y,Yᵢ)
end

Hts = map(x->x[1]=>map(y->(round(mean(y.judge_scores),digits=2),y.λ_origs),x[2]),partitionBy(:heatId));
sort!(Hts)
# Heat x Mean judge score x panel
# since judge score is rounded to 2 decimal places, sco range needs to be 1000
H = zeros(Float32,(910,1000,7,7,7,7,7));
for (i,ht) in enumerate(Hts)
	for (s,λ) in ht[2]
		sco = Int(round(s*100))
		c = prod(factorial.(length.(λ)))
		for a in Base.product( Sym.(λ)...)
			H[i,sco,Int.(vcat(a...))... ] += 1/c
		end
	end
end
#=
scoPD = [sum(H[:,c,:,:,:,:,:]) for c in 1:100]
scoPDF = scoPD / sum(scoPD)
scoCDF = [sum(scoPDF[1:k]) for k in eachindex(scoPDF)] / sum(scoPDF)

ht_sco_PDs = sum(H,dims=[3,4,5,6,7])
sco_ht_PDs = permutedims(dropdims(ht_sco_PDs,dims=(3,4,5,6,7)),[2,1])
sco_ht_PDFs = mapslices(x->x/sum(x), sco_ht_PDs,dims=1)
sco_ht_CDFs = mapslices(x->cumsum(x/sum(x)), sco_ht_PDs,dims=1)
=#
# see this plot
plot([i for i in 0.1:0.1:10],sco_ht_CDFs)
# now I'd like to find some c that minimizes |shCFS 1 - c|
# I think this is the mean
#=
wslHts = map(x->x[1]=>map(y->(round(mean(sort(y.judge_scores)[2:4]),digits=2),y.λ_origs),x[2]),partitionBy(:heatId));
wslH = zeros(Float32,(910,1000,7,7,7));
for (i,ht) in enumerate(wslHts)
	for (s,λ) in ht[2]
		sco = Int(round(s*100))
		c = prod(factorial.(length.(λ)))
		for a in Base.product( Sym.(λ)...)
			wslH[i,sco,Int.(vcat(a...)[2:4])... ] += 1/c
		end
	end
end
wslscoPD = dropdims(sum(wslH, dims=(1,3,4,5)),dims=(1,3,4,5) )
wslscoPDF = wslscoPD/sum(wslscoPD)
wslscoCDF = cumsum(wslscoPDF)

ht_wslsco_PDs = dropdims(sum(wslH,dims=(3,4,5)),dims=(3,4,5))
wslsco_ht_PDs = permutedims(ht_wslsco_PDs, [2,1])
wslsco_ht_PDFs = mapslices(x->x/sum(x), wslsco_ht_PDs,dims=1)
wslsco_ht_CDFs = mapslices(x->cumsum(x), wslsco_ht_PDFs,dims=1)
=#

#' Now just some decompositions
PartSizeDecomp = [zeros(Float64, ntuple(i->7,k) ) for k in 1:5]
for λ in Panel_λs
	c = prod( factorial.(length.(λ)) )
	for p in λ
		k=length(p)
		for a in Sym(p)
			PartSizeDecomp[k][Int.(a)...] += 1/c
		end
	end
end
PSD = PartSizeDecomp
sum(sum.(PSD))
sum.(PSD)
sum.(PSD) ./ sum(sum.(PSD))
AltPart = AltOp.(PSD)
SymPart = SymOp.(PSD)

# ... if we are doing variance decompostion then mayeb:
#			... so my brain ended up doing whats below
#	it seems like I ought to refactor this, perhaps write some math down and see what comes.
#=
P1_2 = PSD[1]⊗PSD[1] /sum(PSD[1])^2
P1_2 *= mean(PSD[2]) / mean(P1_2)
est_c1_2 = map(c->sum(abs.(PSD[2]-c*P1_2)), -2:0.1:5)
plot(-2:0.1:5,est_c1_2, hold=true)
c1_2 = (-2:0.1:5)[findmin(est_c1_2)[2]]

P1_3 = PSD[1]⊗PSD[1]⊗PSD[1] /sum(PSD[1])^3
P1_3 *= mean(PSD[3])/ mean(P1_3)
est_c1_3 = map(c->sum(abs.(PSD[3]-c*P1_3)), -2:0.1:5)
plot(-2:0.1:5,est_c1_3)
c1_3 = (-2:0.1:5)[findmin(est_c1_3)[2]]

P1_4 = PSD[1]⊗PSD[1]⊗PSD[1]⊗PSD[1] /sum(PSD[1])^4;
P1_4 *= mean(PSD[4]) / mean(P1_4);
est_c1_4 = map(c->sum(abs.(PSD[4]-c*P1_4)), -2:0.1:5)
plot(-2:0.1:5,est_c1_4)
c1_4 = (-2:0.1:5)[findmin(est_c1_4)[2]]

P1_5 = PSD[1]⊗PSD[1]⊗PSD[1]⊗PSD[1]⊗PSD[1] /sum(PSD[1])^5;
P1_5 *= mean(PSD[5]) / mean(P1_5);
est_c1_5 = map(c->sum(abs.(PSD[5]-c*P1_5)), -2:0.1:5)
plot(-2:0.1:5,est_c1_5)
c1_5 = (-2:0.1:5)[findmin(est_c1_5)[2]]
legend("c1^2","c1^3","c1^4","c1^5")

P12_3 = (PSD[1]⊗PSD[2] + PSD[2]⊗PSD[1]) /2*(sum(PSD[1])*sum(PSD[2]))
P12_3 *= mean(PSD[3]-c1_3*P1_3)/ mean(P12_3)
est_c12_3 = map(c->sum(abs.((PSD[3]-c1_3*P1_3)-c*P12_3)), -2:0.01:5)
plot(-2:0.01:5,est_c12_3)
c12_3 = (-2:0.01:5)[findmin(est_c12_3)[2]]

P12_4 = (PSD[2]⊗PSD[2])/sum(PSD[2]) +(PSD[2]⊗PSD[1]⊗PSD[1]+PSD[1]⊗PSD[2]⊗PSD[1]+PSD[1]⊗PSD[1]⊗PSD[2])/(3*sum(PSD[1])^2 *sum(PSD[2]))
P12_4 *= mean(PSD[4]-c1_4*P1_4)/ mean(P12_4)
est_c12_4 = map(c->sum(abs.((PSD[4]-c1_4*P1_4)-c*P12_4)), -2:0.1:5)
plot(-2:0.1:5,est_c12_4)
c12_4 = (-2:0.1:5)[findmin(est_c12_4)[2]]
Figure()





diffs = map(x->length(x.panel)==1 ? 0 : maximum([x.panel[i][1] - x.panel[i-1][1] for i in 2:length(x.panel)]), last.(WAVES))
relativediffs = map(x->length(x.panel)==1 ? 0 : mean(x.judge_scores)/maximum([x.panel[i][1] - x.panel[i-1][1] for i in 2:length(x.panel)]), last.(WAVES))
histogram(relativediffs[findall(!=(Inf), relativediffs)])
relativediffs = map(x->length(x.panel)==1 ? 0 : mean(x.judge_scores)/maximum([round(x.panel[i][1] - x.panel[i-1][1],digits=2) for i in 2:length(x.panel)]), HasMatch)

relativeSpread = map(x->(maximum(x.judge_scores)-minimum(x.judge_scores))/mean(x.judge_scores), last.(WAVES))
histogram(relativeSpread)

relativeStd = map(x->std(x.judge_scores)/mean(x.judge_scores), last.(WAVES))

farfromedgesStd = map(x->std(x.judge_scores)/maximum([10-mean(x.judge_scores),mean(x.judge_scores)]), last.(WAVES))
histogram(farfromedgesStd)### SEEEE THIS

spread = map(x->(maximum(x.judge_scores)-minimum(x.judge_scores)),last.(WAVES) )
fromedge = map(x->maximum([10-mean(x.judge_scores),mean(x.judge_scores)]), last.(WAVES))
spread ./fromedge
=#