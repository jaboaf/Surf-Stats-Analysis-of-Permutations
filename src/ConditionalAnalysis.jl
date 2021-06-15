include("SimplifyComp.jl")
using GRUtils

vars = [:EvtOrig,:AthOrig,:I_match,:multi_ctrys,:panel_dist]

procs = []
for cond_panel_comp in partitionBy((:m_c_supp,))
	comp = cond_panel_comp[1]
	Z = zeros(Float32,(7,7,7,7,7))
	panels = map(x->x.λ_c,sort(cond_panel_comp[2]))
	A = Float32[] # || Alternating Part ||
	D = Float32[] # || Z - SymPart ||
	for (t,λ) in enumerate(panels)
		c = prod( factorial.(length.(λ)) )
		for p in Base.product(Sym.(λ)...)
			Z[Int.(vcat(p...))...] += 1/c
		end
		push!(A, sum(abs.(AltOp(Z))) /t )
		push!(D, sum(abs.(Z - SymOp(Z))) / t )
	end
	push!(procs, comp => [A,D])
end

Aprocs = map(p->p[2][1],procs)
Dprocs = map(p->p[2][2],procs)

# Plot of Antisymmetric Part, looks like it converges to 0 ish
for p in Aprocs plot(1/2 * p,ylim=(0,0.4),hold=true) end
# embedded into unit interval
for p in Aprocs plot(LinRange(0,1,length(p)),1/2 * p,ylim=(0,0.4),hold=true) end

# Plot of (Antisymmetric conv to 0 empirical process, i.e. sqrt(n)(AntiSym_n - 0)
for p in Aprocs plot(1/2 * p .* sqrt.(1:length(p)),ylim=(0,0.4),hold=true) end
# embedded into unit interval
for p in Aprocs plot(LinRange(0,1,length(p)),1/2 * p .* sqrt.(1:length(p)),ylim=(0,0.4),hold=true) end

# same as above but with D procs
for p in Dprocs plot(1/2 * p ,ylim=(0,1),hold=true) end
for p in Dprocs plot(LinRange(0,1,length(p)),1/2 * p ,ylim=(0,1),hold=true) end

for p in Dprocs plot(1/2 * p .* sqrt.(1:length(p)), ylim=(0,6),hold=true) end
for p in Dprocs plot(LinRange(0,1,length(p)),1/2 * p .* sqrt.(1:length(p)), ylim=(0,6),hold=true) end

# just interesting
for p in Dprocs plot(1/2 * p .* collect(1:length(p)), ylim=(0,6),hold=true) end


# What prop of a Dproc is Aproc?
for p in procs plot( p[2][1] ./ p[2][2],ylim=(0,1),hold=true) end
# embedded into [0,1]
for p in procs plot(LinRange(0,1,length(p[2][1])),p[2][1] ./ p[2][2],ylim=(0,1),hold=true) end


# NOTE THESE ARE FOR partion by m_c_supp

#=
Regarding a Probability algebra:
Segal (1953) defined a proabability algebra as the pair (R,E) where R is an algebra and E is a distinguished linear function that satisfy:
1) \forall a \in R, E(a^2) \geq 0, and E(a^2) = 0 only if a = 0.
2) \forall a \in R \exists \mu s.t. E(ab^2) \leq \muE(b^2) \forall b \in R
3) R has a unit, e, and E(e) = 1.

We mention this because of the signifigance of segal's works and to contrast it with our construction.
Re 1): The norm gaurentees this.
Re 2): We do not require this, however our construction satisfies this. For arbitary a, take \mu_a to be the maximal entry.
Re 3): Our construction does not adhere to this. Instead, the unit, e, is the tensor of all 1s and the expectation functional is pointwise multiplication by e.
In our setting, \mathbb{P}(e) = 1, and E(a)=1 iff a is a probability measure.
This makes a lot more sense in the commutative case because for every a, a = 1*a, so P(a) should equal P(1*a)
because P: T^d(V) \rightarrow \mathbb{R} may be defined on \{e_w\} as P(v) = ||p*v||.
so we'd hope that p*1 = p


=#





#=
Lebanon and Mao, Non-Parametric Modelling of Ranked Data
n := number of items
m := "number of samples"

1) Estimate P based on full as well as partial rankings.
2) Assign probabilities in a coherent and contradiction free manner.
3) Estimate P based on partial rankings of different types.
4) Statistical Consistency \hat{P} \rightarrow P  (conv. in prob) as m \rightarrow \infty and n \rightarro \infty
5) Stattiscal Accuracy can be slow for fully ranked data, accelerated when restricted to simpler partial rankings
6) Estimate \hat{P} and computing probabilites \hat{P}(A) should be computationally feasible

We will take the spirit of these criteria.
n := number of things
m := length of sample
K := maximal length of arrangement

- Estimate P based on Data (combiing 1 and 3)
- Assign probabilities in a coherent manner (our interpretation of 2)
	- Subadditivity
- Convergence in distribution as m \rightarrow \infty
- Computationally feasible
- Ability to extrapolate for arrangements of length d > K.

Our methods will be of a non-parametric flavor, however, they are parametric.
The statistic we would like of each data point is the ordered pair:
( integer partition of d, representative of data point)

This is not sufficient, there are multiple representativces for certain observations (those that aren't totally ordered).
A sufficent statistic is:
an ordered l-tuple of monomials (defined by maintained arbitrary order of the items) defining the number of item i in class l_j of 
equivalently, an ordered l-tuple of multisets where t_j(c) = number of item c in the jth class

... These are both larger in dimension than the non-sufficient one?





=#

#=
Begining of methodology.
We love, share, and endorse all of the enthusiasm around permutation groups in statistics (as well as other groups).
Unfortunately, the word and notion of a permutation or permutations get thrown around in settings where they aren't exactly the object at hand.
This distinction is vital to remainder of the paper.

a=(1,2,3,4,5)
a is an arrangement of {1,2,3,4,5}

b=(2,1,3,4,5)
b is an arranegment of {1,2,3,4,5}

a is not a permutation. b is not a permutation.
a is not a permutation of {1,2,3,4,5}. b is not a permutation of {1,2,3,4,5}

a is a re-arrangement of b.
b is a re-arrangement of a.

a is a re-arrangement of a.
b is a re-arrangement of b.

the map, a \mapsto b, is a permutation.
the map, b \mapsto a, is a permutation. It is the inverse of a \mapsto b.

the map, a \mapsto a, is a permutation.
the map, b \mapsto b, is a permutation.
The above 2 maps are the same.


=#


