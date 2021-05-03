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


#' # Table of Contents
#' ## Abstract
#' ## Introduction
#' ## Data
#' ## The Simple Approach
#' ### Formulation
#' ### Results
#' ### Lingering Questions

#' ## Motivation for A Different Approach

#' ## Concrete
#' Does nationality influence the World Surf League Championship Tour?
#' If so, how?, to what degree?, and does this vary depending on nationality?

#' # Introduction
#' The World Surf League (WSL) is the most prominent organizer of international surf compeitions. Each year the WSL organizes a variety of "tours" which include Mens and Womens versions of Big Wave events, the "Longboard Tour", the "Qualifying Series", and the "Championship Tour" (CT). 

#' ## Format
#' Each year, the 32 highest ranked (shortboard) surfers are invited to participate in the "Championship Tour" (CT), which constists of 11 surf compeitions in 7 different countries. Each competition has 7 rounds, each consisting of 1 to 16 heats, and each heat has 2 to 3 surfers. Within a heat, a surfer may attempt to ride any number of waves, but their final heat score is the sum of their two highest scoring waves. The surfer with the highest heat score places 1st in the heat, the surfer with the next highest heat score places 2nd in the heat. In some rounds, heats consist of 3 surfers, in which case the surfer with the third highest heat score will place 3rd in the heat. Each round has a rule that determines which surfers advance, and what (round,heat) they advance to. Surfers that do not advance "exit" the event, and are given some number of points (the closer the exit is to the compeitions final round, the higher the amount of points).

#' ## How are waves scored?
#' In any given heat, there is judging panel comprised of 5 judges. Anytime a surfer attempts to ride a wave, each judge analyzes the surfer's ride with respect to the following criteria and write down a score:
#'   - Commitment and degree of difficulty
#'   - Innovative and progressive maneuvers
#'   - Combination of major maneuvers
#'   - Variety of maneuvers
#'   - Speed, power, and flow
#' Note: that Different elements of this list may be emphasized more or less depending on the location, conditions, and changes of conditions. (Chapter 13 Article 182)
#' Additionally, the General Judging Rules (Chapter 13; Article 183; Section 1,2) state that judges should be visually separated, should not discuss scores, and may not change their scores.

#' Though there are 5 judges and 5 scores, the score a surfer receives, called the "wave score", is the trimmed mean of the scores given by the 5 judges. This is very important information and particularly interesting because it distinguishes the *existence* of biased judges from the effect of biased judging (if it exists).

#' ## The Judging Process
#' The role of judges is not limited to the five judges on a panel. For each Men's CT event, there is 1 international Head Judge, 7 international judges, and 1 international priority Judge (Chapter 13, Article 179.01). And Chapter 13, article 179.17 states that "At CT Events, the number of judges from any one regional area is limited to 3". For each heat, the WSL Head Judge must assure there are at least 5 judges on the panel for each heat and that they are a subset of the 7 International Judges and 1 International Head Judge (Chapter 13, Article 179.13).

#' 5 judges are selected from a pool of 8 judges to form a panel for a heat. During that heat, when a surfer rides a wave, each judge independently observes the ride and writes down a score (based on some broadly defined criteria), which is some number in {0.0, 0.1, 0.2, ..., 9.8, 9.9, 10.0}.

#' ## Data
#' We collected some incredibly rich data on surf compeitions from the 2017, 2018, and 2019 seasons of the Mens World Championship Tour (WCT). Each year the World Surf League (WCT) holds 10 to 11 surf competitions, which are called "events". While the format of events have changed slightly between the 2017 and 2019 seasons, they are all very similar. Each event consists of 7 rounds, and within each round there are some number of heats. A heat is the level at which intra-athlete competition takes place and may consist of 2 or 3 surfers. Throughout a timed heat (usually between 22 and 35 minutes), each athlete may surf any number of waves and their "heat score" is the sum of the scores of their two highest scoring waves. 

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




#' # Definitions
#'   - **The Symmetric Group on a set, ``X``** is ``S_X := Isomorphisms(X,X)``. When |X|<∞,``S_X = \{ τ:X→X | Image(τ) = X \}  = \{τ:X→X ∣ \{τ(x) ∣ x∈X \} = X \} ``
#' 
#'   - ``Sₖ := S_{\{ 1,…,k\} }. ``
#'
#'  - When G is a Group, and 𝔽 is a Field, the **Group Algebra of G over 𝔽 **, denoted 𝔽[G], is the space of formal linear combinations of elements of G. Elements of 𝔽[G] are of the form: ``c₁g₁ + … + cₙgₙ = ∑^n_{i=1} cᵢgᵢ``, where ``cᵢ∈𝔽, gᵢ∈G``. Note that i≠j ⟹ gᵢ ≠ gⱼ because G is a set of elements, so no element occurs with multiplicity. For example, c₁g + c₂g ∉ 𝔽[G] whereas (c₁+c₂)g ∈ 𝔽[G].
#'
#'  - 𝔽[Sₖ] is comprised on formal linear combinations of elements of Sₖ. This may be understood as two different ways:
#'     - There exists a function a:Sₖ→𝔽 and ``A:=∑_{τ∈Sₖ} a(τ)τ`` ∈ 𝔽[Sₖ].
#'     - Or: `` A = [ π₁ …  π_{k!} ] \begin{bmatrix}a_{\pi_1} \\ \vdots  \\ a_{π_{k!}} \end{bmatrix} = \sum_{\tau \in Sₖ } a_τ τ ∈ 𝔽[Sₖ] ``,where τ ∈ Sₖ and a_τ ∈𝔽.
#' 
#'   - **Addition in 𝔽[Sₖ] ** is defined by ``A+B = (∑_{τ} a_τ τ) + (∑_τ b_τ τ) := ∑_τ (a_τ + b_τ) τ ``
#' 
#'   - **Scalar Multiplication in ``𝔽[Sₖ]``** is ``c(A) = c(∑_{τ} a_τ τ)  = ∑_{τ} ca_τ τ``.
#' 
#'   - Multiplication in 𝔽[Sₖ] is *defined* by:
#'     ``A*B =(∑_{τ} a_τ τ)* (∑_π b_π π)
#'     :=∑_{γ ∈ Sₖ}(∑_{τ,π | τπ=γ} a_τ b_π) γ
#'     = ∑_{γ ∈ Sₖ}(∑_{τ ∈ Sₖ} a_τ b_{τ^{-1}γ})γ``.
#' 

#' We should not over complicate ``*``. The *definition* of ``*`` may look odd, but it exactly the same as our basic understanding of multiplication:
#' ``
#' (\sum_{  \tau} a_\tau\tau) *(\sum_{\pi }   b_\pi \pi)
#' = \sum_{  \tau} (a_\tau\tau) *(\sum_{\pi }   b_\pi \pi)
#' = \sum_{  \tau} \sum_{\pi }  (a_\tau\tau) * (b_\pi \pi)
#' = \sum_{  \tau}\sum_{\pi } a_\tau b_\pi \tau \pi
#' =\sum_{\gamma \in S_d}(\sum_{  \tau,\pi | \tau\pi = \gamma  } a_\tau b_\pi) \gamma
#' =(\sum_{\tau} a_\tau \tau)* (\sum_\pi b_\pi \pi)
#' ``
#' Even though ``(∑_{τ} a_τ τ) *(∑_{π} b_π π)`` is equal to the intuitive form, ``∑_{τ}\sum_{π} a_τ b_π τπ``, the latter is not an element of the the group algebra because elements are repeated in the sum, hence our chosen definition. Also, ``*`` merely extends multiplication in the Field,``\cdot: \mathbb{F} \times \mathbb{F} \rightarrow \mathbb{F}``, and the operation in the group, ``∘``, by ``a_τ τ * b_π π = a_τ⋅b_π τ∘π = a_τ b_π τπ``, where the last equality is simply notation-reduction.
#' 
#' A measure on Sₖ, is an element of the group algebra ``ℂ[Sₖ]``.
#' A measure on Sₖ, ``F = ∑_{τ ∈ Sₖ} f_τ τ``, is a probability measure on Sₖ if and only if ``∀ τ ∈ Sₖ f_τ ≥ 0`` and ``∑_{τ ∈ Sₖ} f_τ = 1 ``.
#' 
#' A linear representation of a group G, is a group homomorphism, ``ρ : (G,∘) → (GL(V),⋅)``. A group homomorphism satisfies:
#'   - ∀ x,y ∈ G ρ(x∘y) = ρ(x)⋅ρ(y) where ⋅ is multiplication in GL(V).
#'   - ρ(e) = I, where e is the identity element in G and I is the identity element in GL(V).
#'   - ``∀ x ∈ G, ρ(x^{-1}) = ρ(x)^*``, where ``^*`` denotes involution in GL(V).

#' The representation of a measure F, is:``F̂ := ∑_{τ ∈ Sₖ} F(τ)ρ(τ)``.
#' Note: This is sometimes called the Fourier transform at a representation, I avoid that lingo. 
#' 
#' **Convolution of two functions on Sₖ ** is a binary operation ``A * B := ∑_{τ ∈ Sₖ} a(τ) b(τ^{-1}g)``.
#' 
#' Note:
#' ``
#' \widehat{A*B} = ∑_{γ} (A*B)(γ)ρ(γ) 
#' = ∑_{γ} ∑_{τ} a(γ τ^{-1})b(τ)ρ(γ )
#' = ∑_{τ} ∑_{γτ} a(γ τ τ^{-1})b(τ)ρ(γτ)
#' = ∑_{τ} ∑_{γτ} a(γ)b(τ)ρ(γ)ρ(τ)
#' = ∑_{τ} b(τ)ρ(τ) ∑_{γτ} a(γ)ρ(γ) 
#' = ∑_{τ} b(τ)ρ(τ) Â
#' = Â ∑_{τ} b(τ)ρ(τ)
#' = Â ⋅ B̂``
#' 
#' We take ρₖ the permutation representation acting on the vector space V := ℝᵏ with basis indexed by {1,…,k}. So a typical element of V is of the form: 
#' 
#' ``
#' \begin{pmatrix} a^1 \\ a^2 \\ ⋮ \\ a^k \end{pmatrix} = \begin{pmatrix} a^1 \\ 0 \\ ⋮ \\ 0 \end{pmatrix} + \begin{pmatrix} 0 \\ a^2 \\ ⋮ \\ 0 \end{pmatrix} + … + \begin{pmatrix} 0 \\  ⋮ \\ 0\\ a^k \end{pmatrix} 
#' = a^1\begin{pmatrix} 1 \\ 0 \\ ⋮ \\ 0 \end{pmatrix} + a^2\begin{pmatrix} 0 \\ 1 \\ ⋮ \\ 0 \end{pmatrix} + … + a^k\begin{pmatrix} 0 \\  ⋮ \\ 0\\ 1\end{pmatrix} 
#' = a^1 e₁ + a^2 e₂ + … a^k eₖ
#' ``
#' 
#' This could be rewritten as ``∑_{i ∈ \{1,…,n\} } a^i eᵢ``. But remember, we are interested in functions that act on the basis.
#' 
#' Many authors call ``*`` "convolution", however this terminology is superfluous.
#' 
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

#' Though the definition of a multiset on A as a function f:A→ℕ is plesant to work with, it omits some of the structure intrinsic to a multiset. Namely, consider how a multiset arrises in an applied setting. We have some knowledge or infomation about the judges, ie. some function ``k: J_{pool} → C``, which maps a judge in the judging pool to their nationality c ∈ C. How do we arrive at *the* f:C→ℕ ? By: ``f(c) = ∑_{j ∈ J} 1_{k(j)=c }``. Evidently f is the sum of simple functions of our information. So f really isn't representative of much and truly depends on some reference set or information. Additionally, if we had additional structure on J, for example a total order where j₁<j₂<…<j₅, its not entriely clear how we can extend the structure to the multiset.

#' Hence, we define: a multiset on a set A is **an** equivalence class of 𝓁(A) where:
#'   - ``𝓁(A):= ⋃_{n∈𝑵ˣ}A^{×n} `` where ×n denotes the nth cartesian power. Note that this is certainly a disjoint union.
#'   - ``M(A) := 𝓁(A)/~`` where x~y iff. ⨳x=⨳y and ``∃ τ ∈ S_{⨳x}`` s.t. τx=y, where ⨳x:=n s.t. ``x∈A^{×n}``. (NOTE: Sₙ acts on an n-tuple by permuting the **positions** of entries, not on the entries themsleves).
#' To show that M(A) is well defined, we must show that ~ is an equivalence relation:
#'   * WTS: x~x. PF: ``e ∈ S_{⨳x}`` so ex = x and ⨳x=⨳x, thus x~x.
#'   * WTS: x~y ⟹ y~x. PF: Assume x~y, ie. ⨳x=⨳y and ``∃ τ ∈ S_{⨳x}`` s.t. ``τx = y. τ∈S_{⨳x} ⟹τ^{-1} ∈ S_{⨳x} = S_{⨳y}`` so ``τ^{-1}y = τ^{-1}τx = ex = x``. And we have ⨳y=⨳x, thus y~x.
#'   * WTS [x~y ∧ y~z] ⟹ [x~z]. Assume x~y and y~z. ⨳x=⨳y and ⨳y=⨳z implies that ⨳x=⨳z and ``S_{⨳x}=S_{⨳y}=S_{⨳z}``. By assumtion ``∃ τ,γ∈S_{⨳x}`` s.t. τx=y and γy=z. So γτx = γy = z and ``γτ ∈ S_{⨳x}``. Thus, x~z.

#' Further, let us define:
#'   * ⊗:𝓁(A)×𝓁(B)→𝓁(A∪B) by ⊗(x,y) = x⊗y = (x₁,…,xₘ,y₁,…yₙ). This is well defined because ``x⊗y∈(A∪B)^{×(m+n)} ⊂ 𝓁(A∪B)``.
#'   * And define ⊗:M(A)×M(B)→M(A∪B) by ⊗(mA,mB):={⊗(x,y)|(x,y)∈ mA×mB\} = {x⊗y|(x,y)∈ mA×mB}.

#' To be clear:
#'   * M(A) is a *set of equivalence classes of 𝓁(A)* under ~, and a multiset on A, denoted mA, is an element of M(A), ie. mA is an equivalence class.
#'   * Though (mA₁,…,mAₖ) and ⊗(mA₁,…,mAₖ)=mA₁⊗…⊗mAₖ look very similar, they are not the same.
#'     * (mA₁,…,mAₖ) is a k-tuple of multisets on A which means it is an element of ``M(A)^{×k}`` and thus an element of 𝓁(M(A)) because ``M(A)^{×k} ⊂ 𝓁(M(A))``.
#'     * ⊗(mA₁,…,mAₖ)=mA₁⊗…⊗mAₖ is not in M(A), rather it is a subset of an element of M(A), specifically, ``mA₁⊗…⊗mAₖ ⊂ \{ τx ∣ x∈mA₁⊗…⊗mAₖ ∧ τ∈S_{⨳x} \} = mA ∈ M(A)``.

#' This construction has some nice properties:
#'   - Let ``V_A ``:= span{eₐ|a∈A} over ℝ (any field with char 0 should do i think). Given any function α:M(A)→ℝ, we may define ``v ∈ T(V_A)`` by ``v = ∑_{W ∈ M(A)} α(W)∑_{x∈W}vₓ``, which is completely symmetric by definition of W (``v_x = e_{x^{(1)}}⊗…⊗e_{x^{(⨳x)}}``)
#'   - Any β:mA→ℝ satisfying ``∑_{t∈mA} β(t) = 1``, defines a probability distribution. Moreover, we may represent this as ``P ∈ T(V_A)`` by ``∑_{t∈mA} β(t)vₜ``.
#'   - And more to come.

#' Defn: An *ordered partition* of a mutliset mA is a k-tuple of multisets on A,``λ=(mA^{(1)},…,mA^{(k)})``, satisfying ``⊗λ=⊗(mA^{(1)},…,mA^{(k)})=mA^{(1)}⊗…⊗mA^{(k)}⊆ mA``. By definition of ⨳, ⨳λ=k. 

#' Defn: The set of ordered partitions of mA∈M(A) is ``ℿ_{mA}:=\{λ ∈ 𝓁(M(A)) ∣ ⊗λ ⊆ mA\}``

#' Defn: An *unordered partition* of a multiset mA, is a multiset of multisets on A, where an arbitrary element satisfies ⊗v ⊂ mA. By definition of multiset,V={τv∣``τ∈S_{⨳v}``}.

#' JOJO can you show that the function defintion by blizzard is a faithful represetnation of this ?? that would be useful. Think slowley here Representaiton in a discrete set function space should be the same thone that you understand concretely in GL(V).

#' This notion of a partion is different than the typical ones. P Martin and others use a Partition Algebra first introduced by Brauer in 1913?, we are no experts in this realm but the notion seems most useful when analyzing the automorphisms of partitions of a set. Multisets cannot be accomodated (JOJO you must be more informed than this come on man).

#' A more statistical interpretation of partitions appears in a seminal paper on population genetics by Kingman in 1978 (Random Partitions in Population Genetics). More abstract interpretations come in a later paper (The Representation of Partition Structures). The partitions defined are integer partitions (which is the typical interpretation, see also Vershik's paper Statistical Mechanics of Combinitorial Partitions and Their Limit Shapes). Integer partiions are no different than **set** partitions, one can pass freeely between the two by taking the cardinality of sets, with one assumption: that the elements of the set are exchangeable. This is defined in Pitman's 1995 paper, (Exchangeable and Partially Exchangeable Random Paritions). "Partial" exchangeability is a very intuitive notion. Consider a random variable, λ, which takes values in the set of ordered set partitions of {1,…,n}, denoted ℿₙ. And consider its probabailtiy distribution, p:ℿₙ→ℝ. λ is partially exchangeable iff. ``∀ λ ∈ ℿₙ ∀ τ ∈ S_{⨳λ} p(|λ^{1}|, …, |λ^{⨳λ}|) = p(|λ^{τ^{-1}(1)}|, …, |λ^{τ^{-1}(⨳λ)}|) = τp``. In english, this means that the probabiltiy of observing some ordered parition with λ₁ blocks of size 1, λ₂ blocks of size 2, …, λₙ blocks of size n is the same regarless of the order of the blocks (to be clear ``λₖ := |\{j| |λ^{(j)}| = k \}|``). This defines *partial* exchangeablility (read part-ial exchangeability).

#' The defintion of plain old exchangeablility is distinct. A random partition ``A∈ℿₙ`` is exchangeable iff. p(A) = p(τ̂A) ∀τ ∈ Sₙ. You might be wondering what τ̂ is, and why there is a chapeau ̂ involed:

#' Let τ∈Sₖ. How does τ act on ``(A^{1},…,A^{k})`` ? By ``τ(A^{1},…,A^{k}) = (A_{τ^{-1}(1)},…,A^{τ^{-1}(1)})``, ie. moving ith argument to the τ(i)th position. So it makes sense to define ``τ⊗(A^{1},…,A^{k}) = ⊗(A^{τ^{-1}(1)},…,A^{τ^{-1}(k)})`` because we have defined it that way for functions, and ⊗ is a function.
#' In the case of the first wave of 2018, ridden by Caio Ibelli: A=({BRA, ZAF},{ESP, AUS, USA}) and A ∈ ℿ_{\{AUS,BRA,ESP,USA,ZAF}\}}. So τ:=(1 2)∈S₂, acts on A by τ({BRA, ZAF},{ESP, AUS, USA}) = ({ESP, AUS, USA},{BRA, ZAF}).

#' Let γ ∈ Sₙ. How does γ act on ``(A^{(1)},…,A^{(k)})`` ? It doesn't. However, we can *define a function* ``f_τ``, to act on ``(A_{(1)},…,A_{(k)})`` by ``f_τ(A^{(1)},…,A^{(k)}) = (τ(A^{(1)}),…,τ(A^{(k)})) = (\{τ(a) ∣a∈A^{(1)} \},…,\{τ(a)∣a∈A^{(k)} \})``, and let τ̂:=``f_τ``.

#' Pitman calls this the "natural action of permutations of {1,…,n}  on λ" and refers to earlier references by Kingman and a delightful survey of exchangability by Aldous (Exchangeability and related topics)(we have altered the quote to use {1,…,n} instead of Pitman's ``N_n``).
#' The "natural action" of permutations of a set X on a partition of X is a subjective notion. We agree Pitman's thinking here and are fans of the chapeau. 
#' In practice this looks like taking some ``γ∈S_{\{AUS,BRA,ESP,USA,ZAF\}}``, perhaps γ=``\begin{pmatrix} AUS & BRA & ESP & USA & ZAF \\ BRA & USA & ZAF & AUS & ESP \end{pmatrix}`` (this is arbitrary choice). γ̂({BRA, ZAF},{ESP, AUS, USA}) = (γ({BRA, ZAF}),γ({ESP, AUS, USA})) = ({γ(BRA),γ(ZAF)},{γ(ESP),γ(AUS),γ(USA)}) = ({USA,ESP},{ZAF,BRA,AUS}).
#' In this concrete setting, the definition of exchnagability requires that ProbabilityOf({BRA, ZAF},{ESP, AUS, USA}) = ProbabilityOf({γ(BRA),γ(ZAF)},{γ(ESP), γ(AUS),γ(USA)}) for every permutation, γ, of {AUS,BRA,ESP,USA,ZAF}. If the distribtuion over ordered panels is exchangeable, then we can assert: the probability of any arrangement of judges is equal to the the probability of that arrangement under any permutation of judge labels. Further analysis would be required to see if this holds conditional on a surfer's country. But is this the right notion of exchangability that we are after?

#' Remeber the first wave of the second heat of 2018? Michael Rodrigues started the heat of with a whopping 0.3, but taking the panel the under an equivalence relation that judge a ~ judge b iff. judge's a score = judge b's score, we ended up with an ordered partition of the panel: ([BRA],[AUS,AUS,USA],[BRA]). Our desired exchangeability is: P([BRA],[AUS,AUS,USA],[BRA]) = P(γ̂([BRA],[AUS,AUS,USA],[BRA])).

#' Taking a similar γ:=``\begin{pmatrix} AUS & BRA & USA \\ BRA & USA & AUS \end{pmatrix}`` as before, P(γ[BRA],γ[AUS,AUS,USA],γ[BRA])= P([γ(BRA)],[γ(AUS),γ(AUS),γ(USA)],[γ(BRA)]) = P([USA],[BRA,BRA,AUS],[USA]). This seems a little hard to expect, at least within our heat, because the panel for this heat constists of 2 Australian judges, 2 Brazilian judges, and 1 American judge. Or maybe we would like to be strict and require that arrangements, with their multiplicies, are equally probable under any permutation of judge origins?

#' (JOJO i need you to be able to pin down when this is a "fair" or somewhat "fair" expectation and when it is not. Does this nessesitate that countries are equally represented? Yes i think, if we are seeing ordered arrangements of our full multiset, but no if we are seeing ordered arrangements of a subset of the multiset, for example, we could have exactly 5 judges of each contry. But it is also viable if we have ≥ 5 judges from each ctry. How does the definition change when we use γ∈S_{ctry set} vs. γ∈S₅? Generalizing this seems incredibly fun. )

#' Regardless of whether this is a fair notion, you may notice that we passed the action of γ to the elements of our multiset, without justifciation. This is not justified. Nor have we defined anything that would work this way (yet).These are the difficulties with multisets that motivate our definitions and a lot of the interest in these panels.

#' Retain part-ial exchangeability as [ P(mA₁,…,mAₖ)= τP(mA₁,…,mAₖ) ∀τ∈Sₖ ].

#' Define exchangeability as [ P(⊗λ)=P(γ⊗λ) ∀γ∈Sₙ ] where n:=∑λₖ (note that γ is applies to elements of λ which are all n-tuples so action of Sₙ makes sense, and note: P(γ(λ)) = γ̂P(λ)

#' Define label exchangeability as [ P(λ) = P(β̂(λ)) ∀λ∈ℿ_C ∀β∈S_C ] where C = {AUS,BRA,ESP,FRA,PRT,USA,ZAF}, and for ω∈λ β̂ω = β̂(ω₁,…,ωₙ) = (β(ω₁),…	,β(ωₙ)), and note: ``P(β̂(λ))=f_{β̂}P(λ)``

#' I have a bag of judges from in 2018 or 2019
#' (Y_2018 bag for 2018, Y_2019 bag for 2019)
#' ((evt1 bag, …,evt11 bag), (evt1 bag,…,evt11bag))
#' ....

#' given pool = (n_AUS,n_BRA,n_ESP,n_FRA,n_PRT,n_USA,n_ZAF)
#' given panel = (k_AUS,k_BRA,k_ESP,k_FRA,k_PRT,k_USA,k_ZAF)

#=
P(Heat | Panel ) = (1/#waves∑(k_AUS!k_BRA!k_ESP!k_FRA!k_PRT!k_USA!k_ZAF / 5!)∏1/λ_i!) 
P(Heat | Panel ) = (1/(n_AUS)(n_AUS-1)…(n_AUS-k_AUS))⋅…⋅(1/(n_ZAF)(n_ZAF-1)…(n_ZAF-k_ZAF))(1/#waves∑(k_AUS!k_BRA!k_ESP!k_FRA!k_PRT!k_USA!k_ZAF / 5!)∏1/λ_i!) 
=#
#' JOJO how bout we implement a multivariate-hypergeometric model for panel composition, you can use multinomial prior, and i think you can probably get marginals that look like multinomials, but you will have some vector p corresponding to ORIG concentration to estimate, but i think we know what p is empirically for almost every event. So we should be able to assess deviations from our expectations, or just do this fully bayesian (but then we'd end up with a result for an estimate of a quantity we already know...)

#' JOJO How bout you actually just make a quick through all of the rank tests you know / can find. they really aren't too hard to code, just do it. Then we can see how well they do, and maybe improve upon them? Figure out WHY they work or don't work.

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