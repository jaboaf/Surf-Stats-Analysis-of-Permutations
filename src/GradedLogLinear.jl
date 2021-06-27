include("SimplifyComp.jl")
L = [ w.λ_c for w in Tups]; # observed panel
Q = [ w.λ_eqc for w in Tups]; # observed panel

Nᵪ(X::Array{ORIG,1}) = vec([count(==(c),X) for c in JUD_ORIGS])
Nᵪ(X::NTuple{N,ORIG}) where N = map(x-> vec([count(==(c),x) for c in JUD_ORIGS]),X)
N̂ᵪ(X::NTuple{N,Array{ORIG,1}}) where N = map(Nᵪ,X)

# This is suff. statistic for ~graded~ log linear model!
# A = ⨁̂a = (⨁a¹,…,⨁a^#)
# = a⨁
A_L = ⨁(map(l->hcat(N̂ᵪ(l)...),L))
A_Q = ⨁(map(l->hcat(N̂ᵪ(l)...),Q))
# 𝔸 = A
𝔸_L = hcat([ sum([ A_L[k][:,j] for k in j:5 ]) for j in 1:5 ]...)
𝔸_Q = hcat([ sum([ A_Q[k][:,j] for k in j:5 ]) for j in 1:5 ]...)
# The number of judges observed is:
njuds = sum(𝔸) # one can check that this is length(WAVES)*5
# Judge Orig margingal is
L_JudOrig_marg = 𝔸_L*ones(5)
Q_JudOrig_marg = 𝔸_Q*ones(5)
# Block marginal is
L_Blk_marg = ones(7)' * 𝔸_L
Q_Blk_marg = ones(7)' * 𝔸_Q

# Estimated joint distribution assuming independence is
L_estJoint = L_JudOrig_marg * L_Blk_marg/njuds
Q_estJoint = Q_JudOrig_marg * Q_Blk_marg/njuds

# Chi Squared
sum((𝔸_L - L_estJoint) .^2 ./ L_estJoint)
sum((𝔸_Q - Q_estJoint) .^2 ./ Q_estJoint)



b = map(l-> ( map(c->c in vcat(l...),JUD_ORIGS),embedd(l;strong=true,minimal=true)),L);
B = [ (m,⨁(last.(b)[findall(==(m),first.(b))] ) ) for m in unique(first.(b)) ];
# N = number of judges
N = sum(𝔸)



ℙJOrig = margJudOrig/njuds
ℙPod = margPod/njuds

Ejoint = ℙJOrig * ℙPod
ℙ𝔸 = 𝔸/sum(𝔸)
