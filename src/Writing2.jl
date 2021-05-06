#' Some Preliminaries:
using GRUtils # for visualization
using Statistics # for general purpose stats things
using StatsBase

# Currently whats being used?
include("SimplifyComp.jl")
include("toolkit/SymGrpAndReps.jl")

#' Here we read in the raw data. I have elected to restrict our attention to 2018 and 2019 because those are the most complete years. See Table:
#' Year | Waves | Waves w/ 0 Judge Origins Listed | Waves w/ 5 Judge Origins Listed | Total Number of Judge Origins | Waves w/ 3 Sub Scores Listed | Waves w/ 5 Sub Scores Listed | Total number of Judge Scores
#'-----|-----|-----|-----|-----|-----|-----|-----
#'2017 | 7328 | 5210 | 2118 | 10590 | 299 | 7029 | 36042
#'2018 | 6639 | 336 | 6303 | 31515 | 0 | 6639 | 33195
#'2019 | 7648 | 79 | 7569 | 37845 | 0 | 7648 | 38240
#'All | 21615 | 5625 | 15990 | 79950 | 299 | 21316 | 107477


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
#'   - **The Symmetric Group on a set, ``X``** is ``S_X := Isomorphisms(X,X)``. When ``|X|<\infty``,``S_X = \{ \tau :X\rightarrow X \mid Image(\tau) = X \}  = \{\tau:X\rightarrow X \mid \{\tau(x) \mid x \in X\} = X \} ``
#' 
#'   - ``S_d := S_{\{ 1 ,\dots, d\}}. ``
#'
#'  - When G is a Group, and ``\mathbb{F}`` is a Field, the **Group Algebra of G over ``F``**, denoted ``\mathbb{F}[G]``, is the space of formal linear combinations of elements of G. Elements of ``\mathbb{F}[G]`` are of the form: ``c_1 g_1 + \dots + c_n g_n = \sum^n_{i=1} c_i g_i``, where ``c_i \in \mathbb{F}, g_i \in G``. Note that ``i\neq j \implies g_i \neq g_j`` because G is a set of elements, so no element occurs with multiplicity. For example, ``c_1 g + c_2 g \not\in \mathbb{F}[G]`` whereas ``(c_1 + c_2)g \in \mathbb{F}[G]``.
#'
#'  - ``\mathbb{F}[S_d]`` is comprised on formal linear combinations of elements of ``S_d``. This may be understood as two different ways:
#'     - There exists a function ``a:S_n \rightarrow \mathbb{F}`` and ``A :=  \sum_{\tau \in S_d} a(\tau) \tau \quad \in \mathbb{F}[S_n]``.
#'     - Or: `` A = [ \pi_1 \dots  \pi_{d!} ] \begin{bmatrix}a_{\pi_1} \\ \vdots  \\ a_{\pi_{d!}} \end{bmatrix} = \sum_{\tau \in S_d } a_\tau \tau \in \mathbb{F} [S_d] ``,where ``\tau \in S_d, a_\tau \in \mathbb{F}``.
#' 
#'   - **Addition in ``\mathbb{F}[S_n]``** is defined by ``A+B = (\sum_{\tau} a_\tau \tau) + (\sum_\tau b_\tau \tau) := \sum_\tau (a_\tau + b_\tau) \tau ``
#' 
#'   - **Scalar Multiplication in ``\mathbb{F}[S_d]``** is ``c(A) = c(\sum_{\tau} a_\tau \tau)  = \sum_{\tau} ca_\tau \tau``.
#' 
#'   - Multiplication in ``\mathbb{F}[S_d]`` is *defined* by:
#'     ``A*B =(\sum_{\tau} a_\tau \tau)* (\sum_\pi b_\pi \pi)
#'     :=\sum_{\gamma \in S_d}(\sum_{\tau,\pi | \tau\pi=\gamma} a_\tau b_\pi) \gamma
#'     = \sum_{\gamma \in S_d}(\sum_{  \tau \in S_d} a_\tau b_{\tau^{-1}\gamma} ) \gamma ``.
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
#' A measure on ``S_d``, is an element of the group algebra ``\mathbb{C}[S_n]``.
#' A measure on ``S_d``, ``F = ∑_{τ ∈ S_d} f_τ τ``, is a probability measure on ``S_d`` if and only if ``∀ τ ∈ S_d f_τ \geq 0 `` and ``∑_{τ ∈ S_d} f_τ = 1 ``.
#' 
#' A linear representation of a group G, is a group homomorphism, ``\rho : (G,\circ) \rightarrow (GL(V),\cdot)``. A group homomorphism satisfies:
#'   
#'   - `` ∀ x,y ∈ G \quad ρ(x∘y) = ρ(x) \cdot ρ(y)`` where ``\cdot`` is multiplication in ``GL(V)``.
#'   - `` ρ(e) = I ``, where e is the identity element in ``G`` and I is the identity element in ``GL(V)``.
#'   - ``∀ x ∈ G, ρ(x^{-1}) = ρ(x)^*``, where ``^*`` denotes involution in ``GL(V)``.
#' 
#' The representation of a measure F, is: `` \hat{F} := ∑_{τ ∈ S_n} P(τ) ρ(τ) ``.
#' Note: This is sometimes called the Fourier transform at a representation, I avoid that lingo. 
#' 
#' **Convolution of two functions on ``S_d``** is a binary operation `` A * B := ∑_{τ ∈ S_n} a(τ) b(τ^{-1}g)``.
#' 
#' Note:
#' ``
#' \widehat{A*B} = ∑_{γ} (A*B)(γ)ρ(γ) 
#' = ∑_{γ} ∑_{τ} a(γ τ^{-1})b(τ)\rho(γ )
#' = ∑_{τ} ∑_{γτ} a(γ τ τ^{-1})b(τ)ρ(γτ)
#' = ∑_{τ} ∑_{γτ} a(γ)b(τ)ρ(γ)ρ(τ)
#' = ∑_{τ} b(τ)ρ(τ) ∑_{γτ} a(γ)ρ(γ) 
#' = ∑_{τ} b(τ)ρ(τ) \hat{A}
#' = \hat{A} ∑_{τ} b(τ)\rho(τ)
#' = \hat{A} \cdot \hat{B}
#' ``
#' 
#' We take ``ρ_d`` the permutation representation acting on the vector space ``V := \mathbb{R}^d`` with basis indexed by ``\{1,…,d\}``. So a typical element of V is of the form: 
#' 
#' ``
#' \begin{pmatrix} a^1 \\ a^2 \\ \vdots \\ a^n \end{pmatrix} = \begin{pmatrix} a^1 \\ 0 \\ \vdots \\ 0 \end{pmatrix} + \begin{pmatrix} 0 \\ a^2 \\ \vdots \\ 0 \end{pmatrix} + \dots + \begin{pmatrix} 0 \\  \vdots \\ 0\\ a^n \end{pmatrix} 
#' = a^1\begin{pmatrix} 1 \\ 0 \\ \vdots \\ 0 \end{pmatrix} + a^2\begin{pmatrix} 0 \\ 1 \\ \vdots \\ 0 \end{pmatrix} + \dots + a^n\begin{pmatrix} 0 \\  \vdots \\ 0\\ 1\end{pmatrix} 
#' = a^1 e_1 + a^2 e_2 + \dots a^n e_n
#' ``
#' 
#' This could be rewritten as ``∑_{i ∈ \{1,…,n\} } a^i e_i``. But remember, we are interested in functions that act on the basis.
#' 
#' 
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
Panel_λs = map(x->x.λ_origs, last.(WAVES))
N = length(Panel_λs)
Ord_Parts = map(λ -> length.(λ) , Panel_λs)
Ord_Parts_Counts = countmap(Ord_Parts)
Ord_Parts_Props = Dict( [x=>Ord_Parts_Counts[x]/N  for x in keys(Ord_Parts_Counts)] )

#' ``
#' |\{ g ∈ S_d | cycles(g) = 1^{k₁},2^{k₂},…,d^{k_d} \}| = \frac{d!}{\prod_{j=1}^{d} k_j!j^k_j }
#' ⟹ P(1^k_1, … , d^k_d) = \frac{\frac{d!}{\prod_{j=1}^{d} k_j!j^k_j }}{d!} = \frac{1}{\prod_{j=1}^{d} k_j!j^k_j }
#' ``
#' Panels is the Array of the observed Panels, each of which is ordered by score.
#' We do not have a total order for every panel
#' So a observed panel is an ordered partition.
#' Y is an array of cycle counts.
#' Below, λ, is the cycle counts.
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

#' Question: Are lengths of partions asymptotically normal?
lParts = length.(eqPanels)
histogram(lParts)


#' Comprende?
D = map(eqPanels) do x
	y = Set.(last.(x))
	P = Iterators.product(Sym.(y)...)
	Perms = Array{Int8}[]
	onPan = sort([ Int(c) for c in union(y...)])
	notPan = setdiff(1:7,onPan)
	for z in P
		p = union(z...)
		perm = zeros(Int8, 7)
		perm[ notPan] .= notPan
		perm[ onPan ] .= [ Int(c) for c in p ]
		push!(Perms, perm )
	end
	return Perms
end

Reps = map(eqPanels) do x
	y = Set.(last.(x))
	P = Iterators.product(Sym.(y)...)
	Perms = Array{Int8}[]
	onPan = sort([ Int(c) for c in union(y...)])
	notPan = setdiff(1:7,onPan)
	repNotPan = [ ((i==j) & (i in notPan)) ? 1 : 0 for i in 1:7, j in 1:7]
	for z in P
		p = union(z...)
		perm = zeros(Int8, 7)
		perm[ notPan] .= notPan
		perm[ onPan ] .= [ Int(c) for c in p ]
		push!(Perms, Rep(perm))
	end
	return sum(Perms) / length(Perms) - repNotPan
end

expApp(A::Array{Float64,2}) = sum([ A^k / factorial(k) for k in 0:10] )

Reps = map(x->sum(map(y->Rep(y),x))/length(x), D)
Reps = expApp.(Reps) ./ exp(1)

#=
videofile("MomentsForMixtureElement.mp4") do
	for n in 1:length(Reps)
		draw(gcf(wireframe(prod(Reps[1:n]))))
	end
end
=#
#=
[
exp(im*tα)
for α in A
]

map(x-> map(y->Set(y[2]),x),eqPanels)
=#

#' The most granular object of study is the panel. It would not be justified to analyze our data a sequence of 13,872 panels because during each heat the set of judges is fixed (or at least it appears to be this way).
all( all(ht[2][i-1].panel_origs == ht[2][i].panel_origs for i in 2:length(ht[2])) for ht in partitionBy(:heatId) ) 
#' For this reason we will be mostly interested in heat-level dynamics of the panel. This is similar to the game-level approach Price takes in analyzing NBA refereeing (as opposed to simply many foul calls).
#' A panel is comprised of a set of distinct humans (Judges), J = {j₁, j₂, ..., j₅}, and their origns by the multiset J_c:C→𝑵.
#' Multisets are interesting. Roman defines a multiset on a set A, as an element of A×N. Blizzard gives in broad overview of different approaches to multisets in [The Development of Multiset Theory], and takes a more normative approach in [Dedekind Multisets and Function Shells]. We will not delve into this, but they are important to understand in the construction below:
#' DEFN: A multiset on a set S, is a function f:S→𝑵.
#' Rmk: Some authors use the pair (S,f) to define a multiset. While it may bring us comfort to know the underlying set which the function is defined on, this is merely a luxury of defining a multiset. In an observational setting one does not nessisary know the underlying set, only the support of f, ie. the set of all elements mapped to a non-zero number. For example, the judges in a heat are selected from a pool of judges provided by the WSL for that event. Specifically, suppose it is march ___ 2018, Caio Ibelli rides the first wave of the WSL 2018 Mens Championship Tour, when all 5 judges provide their scores, Caio Ibelli is awarded a score of:
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

#' Though the definition of a multiset on A as a function f:C→N is plesant to work with, it omits some of the structure intrinsic to a multiset. Namely, consider how a multiset arrises in an applied setting. We have some knowledge or infomation about the judges, ie. some function k: J_pool → C, which maps a judge in the judging pool to their nationality c ∈ C. How do we arrive at *the* f:C→N ? By: `` f(c) = ∑_{j ∈ J} 1_{k(j)=c} ``. Evidently f is the sum of simple functions of our information. So f really isn't representative of much and truly depends on some reference set or information. Additionally, if we had additional structure on J, for example a total order where ``j₁ < j₂ < … < j₅ ``, its not entriely clear how we can extend the structure to the multiset.

#` Hence, we define: a multiset on a set A is an equivalence class of ``𝓁(A):= ⋃_{n∈𝑵ˣ}A^{×n} `` where ×n denotes the nth cartesian power. And let ``M(A) := 𝓁(A)_{/~}`` where x~y iff. ``y ∈ S_{|x|}x`` , and |x| is the number of entries of x. (NOTE: S_n acts on an n-tuple by permuting the positions of the entries).

#' M(A) has some nice properties:
#'   - Let ``V_A := \span \{e_a|a∈A \}`` over 𝑪. Given any function ``α:M(A)→𝑪``, we may define ``v ∈ T(V_A)`` by ``v = ∑_{W ∈ M(A)} α(W)∑_{x∈W}v_x`` which is completely symmetric by definition of W.
#`   - We can define ⊗:𝓁(A)x𝓁(A)→ by

partitionBy(:heatId)[2][2][1].panel_origs
partitionBy(:heatId)[1][2][1].panel_origs

# are in you observe
partitionBy(:heatId)[1][2][1].panel
partitionBy(:heatId)[2][2][1].panel_origs


#' some utilities that will be of use
⊗(A::Array{T},B::Array{T}) where T<: Number = prod.(Base.product(A,B))

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

#' D is currently the entire data set, This isn't partitularly accurate of a view of the data. Instead lets chop shit up.

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

μ = sum(Y)/ length(Y) ;
AbsErr = map(y-> abs.(y - μ), Y);
histogram( sum.(AbsErr) ) 

SqErr = map(y -> mapreduce(prod,+,Base.product(y-μ,y-μ)), Y)


#' Asymptotically cycle lengths are asymptotically normal? So Asymptotically, partition blocks of size k will have normally distributed sub reps with params based on k? in which case by ind. normal + ind. normal = ind. normal, we have asymptotic normality of partition sizes?

#' Now just some decompositions

PartSizeDecomp = [zeros(Float64, ntuple(i->7,k) ) for k in 1:5]
for λ in Panel_λs
	c = prod( factorial.(length.(λ)) )
	for p in λ
		k=length(p)
		for a in Sym(p)
			PartSizeDecomp[k][Int.(a)...] += 1//c
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
PSD ./= sum(sum.(PSD))
D1to2 = PSD[1]*transpose(PSD[1])
D1to2 /= sum(D1to2) 
D1to2 *= sum(PSD[2])





