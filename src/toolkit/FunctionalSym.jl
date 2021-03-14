import Base: one,inv,*,^,==
include("SymGrpAndReps.jl")

struct Perm{T}
	A::Vector{T}
end

# Making perms into functors
function (p::Perm)(x) return p.A[x] end

Base.one(p::Perm) = Perm(sort(p.A))
Base.inv(p::Perm) = Perm(invperm(p.A))
Base.:*(a::Perm, b::Perm) = Perm(a.A[b.A])
function ^(p::Perm, k::Integer)
	if k > 1 return p * ^(p, k-1)
	elseif k < 0 return ^(inv(p), abs(k))
	elseif k == 1 return p
	else return one(p)
	end
end

==(p::Perm{T}, q::Perm{T}) where T  = (p.A == q.A)
sgn(p::Perm) = sgn(p.A)
function period(p::Perm)
	k=1
	while p^k != one(p)
		k += 1
	end
	return k
end

Base.:*(A::Set{Perm}, B::Set{Perm}) = Set( *(z...) for z in Base.product(A,B) )
Base.:*(X::Set{Perm}, Y::Set{Perm}, Z::Set{Perm}...) = Set( *(z...) for z in Base.product(X...) )
Base.:*(A::Array{Perm}, B::Array{Perm}) = [ *(z...) for z in Base.product(A,B) ]
Base.:*(X::Array{Perm}...) = [ *(x...) for x in Base.product(X...)]


SymGroup(S::Set) = Set( (Perm(t)) for t in Sym(S::Set) )

# Conceptually I dont like these. SymGroup() is not 
# Note: Sym calls unique!() on Arrays
SymGroup(S::Array) = Set( (Perm(t)) for t in Sym(S) )
SymGroup(d::Integer) = Set( (Perm(t)) for t in Sym(d) )


#= Sym group idea
struct SymGroup{T} where T <: Symbol
	S::Set{T}
	D Set( (Perm(t)) for t in Sym(S) )
end
=#

mutable struct GrpAlgElem{F<:Number, Perm}
	vals::Vector{F}
	perms::Vector{Perm}
end


# Group is (S::Set, op::Function) where op: S x S --> S

G = SymGroup(4)
SgnDecomp = Dict([-1 => Perm[], 1 => Perm[] ])
for p in G
	push!(SgnDecomp[sgn(p)], p)
end

PeriodDecomp = [Perm[] for i in 1:4]
for p in G
	z = period(p)
	if z in keys(PeriodDecomp) push!(PeriodDecomp[z], p)
	else PeriodDecomp[z] = [p]
	end
end

FixedDecomp = [Perm[] for i in 1:(4+1)]
for p in G
	hasfixed = false
	for i in 1:4
		if p(i) == i
			push!(FixedDecomp[i], p)
			hasfixed = true
		end
	end
	if hasfixed == false push!(FixedDecomp[end], p) end
end

⊗(A::Array{Perm{T},1}, B::Array{Perm{T},1}) where T  = [ a*b for a in A, b in B ]
⊗(A::Array{Perm{T},N}, B::Array{Perm{T},M}) where T  = [ a*b for (a,b) in Iterators.product(A,B) ]


#subsets of Gk
A = rand(G, 14) 
B = rand(G, 14)
#sgn operator
ϵ = [ sgn(a*b) == s for a in A for b in B, s in [-1,1] ] 
ε(X::Array{Perm{T}}, Z::Array{Perm{T}}...) where T = [ sgn( *(z...) ) for z in Iterators.product(X,Y...) ]

∘

SgnDecomp(A::Array{Perm{T},1}, B::Array{Perm{T},1})
sum(ϵ)

ρ(p::Perm) = Rep(p.A)

# Woudl this work? or just use GroupAlgebra elements hmmm
⊗(A::Array{Tuple{UniformScaling{T<:Number},Perm{T}},1}, B::Array{Tuple{uniformScalar,Perm{T}},1}) where T  = [ a*b for a in A for b in B ]
⊗(A::Array{Perm,1}, B::Array{Perm,1})  = [ a*b for a in A for b in B ]
⊗(A::Array{Perm{T},1}, B::Array{Perm{T},1}) where T  = [ a*b for a in A for b in B ]



