using StatsBase
include("SimplifyComp.jl")
WAVES = sort(last.(WAVES));

# COUNTRY aka JudOrigs
M = [ w.m_c for w in WAVES];

L = [ w.λ_c for w in WAVES]; # observed panel
Q = [ w.λ_eqc for w in WAVES]; # observed panel
#Abstraction from multiplicities within parts of panel
NL = [ unique.(w.λ_c) for w in WAVES ]; # i.e. map( cls -> unique(cls), w.λ_c)
NQ = [ unique.(w.λ_eqc) for w in WAVES ]; # Abstraction from multiplicities within parts of panel under equivalence

# BINARY aka Match or NoMatch
M_b = [ w.m_b for w in WAVES];

L_b = [ w.λ_b for w in WAVES];
Q_b = [ w.λ_eqb for w in WAVES];
NL_b = [ unique.(w.λ_b) for w in WAVES ];
NQ_b = [ unique.(w.λ_eqb) for w in WAVES ];

ML = [ (w.m_c_supp,w.λ_c) for w in WAVES]; 
MQ = [ (w.m_c_supp,w.λ_eqc) for w in WAVES]; 

# Interesting to note:
# W = M .⊗ L gives an error
# ERROR: MethodError: no method matching ⊗(::NTuple{7,Int64}, ::Tuple{Array{ORIG,1},Array{ORIG,1}})
# typeof(M) --> Array{NTuple{7,Int64},1}
# typeof(L) --> Array{Tuple{Array{ORIG,1},Vararg{Array{ORIG,1},N} where N},1}
# It makes sense that this error comes about. 

# note: our data is exactly Tuple{Array{Orig,1},Vararg{Array{ORIG,1},N} where N}
# the data is kind of NTuple{N,Array{Orig,1}} where N
# (bc. we always get some ordered tuple of an array of Origs, i.e. at least a 1 tuple and NTuple permits for 0-tuple)
e(c::ORIG) = vcat(zeros(Int(c)-1),1,zeros(7-Int(c)))
e(b::Bool) = [ Int(b==0); Int(b==1)] 
e(i::Integer,n::Integer) = vcat(zeros(i-1),1,zeros(n-i)) 
function embedd(t::Tuple{T,Vararg{T,N} where N} where T <:Union{Array{ORIG,1},BitArray{1}};strong=false, countinblock=false, minimal=false, normalize=false)
	if minimal==true # then embedd in VS with as many dims are needed
		orgs = unique(union(t...))
		k = length(orgs)
		idx = Dict([c => i for (i,c) in enumerate(sort(orgs)) ])
		if countinblock==true
			q = map(cls->[ count(==(c),cls)*e(idx[c],k) for c in unique(cls)],t)
		else
			q = map(cls->map(c->e(idx[c],k),cls),t)
		end
	else # embedd in VS with e_Int(label)
		if countinblock==true
			q = map(cls->[count(==(c),cls)*e(c) for c in unique(cls)],t)
		else
			q = map(cls->e.(cls),t)
		end
	end
	if normalize==true
		q = map(pod -> pod ./ sum(sum(pod)), q)
	end
	if strong==true # then do sum within pod
		q = map(pod->sum(pod),q)
		return ⊗(q...)
	else
		blocks = map(pod->SymOp(⊗(pod...)),q)
		return ⊗(blocks...)
	end
end


# idea!
### Binary 
TL_b = sum(map(x->embedd(x),L_b));
TL_bs = ⨁(map(x->embedd(x;strong=true),L_b));
TL_bm = ⨁(map(x->embedd(x;minimal=true),L_b));
TL_bsm = ⨁(map(x->embedd(x;strong=true,minimal=true),L_b));

###  CTRY
TL = sum(map(x->embedd(x),L));
TLv = ⨁(map(x->embedd(x;countinblock=true),L));
TLs = ⨁(map(x->embedd(x;strong=true),L));
TLvs = ⨁(map(x->embedd(x;strong=true,countinblock=true),L));

TQ = sum(map(x->embedd(x),L)); # only has 1 size type
TQv = ⨁(map(x->embedd(x;countinblock=true),Q));
TQs = ⨁(map(x->embedd(x;strong=true),Q));
TQvs = ⨁(map(x->embedd(x;strong=true,countinblock=true),Q));

μₐ(D::Array{Array{T,N} where N}) where T<:Number = map(x->sum(abs.(AltOp(x))),⨁(D))
μₛ(D::Array{Array{T,N} where N}) where T<:Number = map(x->sum(abs.(SymOp(x))),⨁(D))

size(TL)
size(QL)

size.(TLv)
size.(QLv)

size.(TLs)
size.(QLs)

size.(TLvs)
size.(QLvs)

barplot([sum.(TLs) sum.(TLvs) sum.(TQs) sum.(TQvs)],labels=["TLs","TLvs","TQs","TQvs"])

barplot([sum.(TL_b) sum.(TL) sum.(TQ_b) sum.(TQ)],labels=["TLe_b","TL","TQe_b","TQ"])
barplot(sum.(TL_subc),label="TLe_subc")
barplot(sum.(TQ_subc),label="TQe_subc")

barplot([μₐ(TL_b) μₐ(TL) μₐ(TQ_b) μₐ(TQ)],labels=["TLe_b","TL","TQe_b","TQ"])

barplot(μₐ(TL_subc),label="TLe_subc")
barplot(μₐ(TQ_subc),label="TQe_subc")

barplot(μₐ(TL_subc) ./ μₛ(TL_subc),label="TLe_subc")
barplot(μₐ(TQ_subc) ./ μₛ(TQ_subc),label="TQe_subc")

barplot([μₐ(TL_b) μₐ(TL) μₐ(TQ_b) μₐ(TQ)] ./ [μₛ(TL_b) μₛ(TL) μₛ(TQ_b) μₛ(TQ)],labels=["TLe_b","TL","TQe_b","TQ"])

# now you may apply vector-variate partion ideas...
# go look at vershink and kerov asymptotics
# I would like if you did the above wihtout the summing of pods
# maybe do it for:
# 1) each judge independent slot
# 2) each podXctry an independednt slot (so sum (within ctry, within pod) i.e. when 2+ judges from same ctry in same pod)


# This is suff. statistic for ~graded~ log linear model!
# A = ⨁̂a = (⨁a¹,…,⨁a^{\#})
# = a⨁
A = ⨁(map(l-> hcat( [sum(e.(pod)) for pod in l]... ),L))
# 𝔸 = A
𝔸 = hcat([ sum([ A[k][:,j] for k in j:5 ]) for j in 1:5 ]...)
b = map(l-> ( map(c->c in vcat(l...),JUD_ORIGS),embedd(l;strong=true,minimal=true)),L);
B = [ (m,⨁(last.(b)[findall(==(m),first.(b))] ) ) for m in unique(first.(b)) ];

# N = number of judges
N = sum(𝔸)

# Judge Orig margingal is
margJudOrig = 𝔸*ones(5)
# class marginal is
margPod = ones(7)' * 𝔸
# estimated joint distribution assuming independence is
estJoint = margJudOrig * margPod / njuds
# Chi Squared
sum((𝔸 - estJoint) .^2 ./ estJoint)

ℙJOrig = margJudOrig/njuds
ℙPod = margPod/njuds

Ejoint = ℙJOrig * ℙPod
ℙ𝔸 = 𝔸/sum(𝔸)

reducedL = map(L) do l
orgs = unique(union(l...))
k = length(orgs)
idx = Dict([c => i for (i,c) in enumerate(sort(orgs)) ])
q = map(cls->sum([e(idx[c],k) for c in cls])/length(cls),l)
reduce(⊗,q)
end;
GL = ⨁(reducedL);
size.(GL)

reducedQ = map(Q) do l
orgs = unique(union(l...))
k = length(orgs)
idx = Dict([c => i for (i,c) in enumerate(sort(orgs)) ])
q = map(cls->sum([e(idx[c],k) for c in cls])/length(cls),l)
reduce(⊗,q)
end;
GQ = ⨁(reducedQ);
size.(GQ)



#
DL = map(ML) do d
k = sum(d[1])
idx = cumsum(d[1])
q = [ [ (count(==(c),cls)/length(cls))*e(idx[Int(c)],k) for c in unique(cls)] for cls in d[2]]
v = mapreduce(p->reduce(⊗,p), +, Base.product( q...))
s = zeros((2,2,2,2,2,2,2))
s[1 .+ d[1]...] = 1
return s⊗v
end;

# this gives a prob measure i think.
PL = ⨁(DL);





