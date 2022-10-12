import Base: one,inv,*,^,==
include("SymGrpAndReps.jl")

struct Perm{T}
	A::Vector{T}
end

#struct Sym{T}
#	S::Set{T,1}
#end3

# Making perms into functors
function (p::Perm{T})(x::T) where T return p.A[x] end

Base.one(p::Perm) = Perm(sort(p.A))
Base.inv(p::Perm) = Perm(invperm(p.A))

# Definitions of product
Base.:*(a::Perm, b::Perm) = Perm(a.A[b.A])
Base.:*(A::Set{Perm{T}}, B::Set{Perm{T}}) where T = Set( a*b for a in A for b in B )
Base.:*(A::Array{Perm{T}}, B::Array{Perm{T}}) where T = [ a*b for (a,b) in Base.product(A,B) ]

⊗(A::Perm, B::Perm) = Perm([a.A...,b.A...])


function ^(p::Perm, k::Integer)
	if k > 1 return p * ^(p, k-1)
	elseif k < 0 return ^(inv(p), abs(k))
	elseif k == 1 return p
	else return one(p)
	end
end
^(a::Perm{T}, b::Perm{T}) where {T} = b * a * b^-1
^(H::Set{Perm{T}}, b::Perm{T}) where {T} = b * H * b^-1

==(p::Perm, q::Perm) = (p.A == q.A)
sgn(p::Perm) = sgn(p.A)
function period(p::Perm)
	k=1
	while p^k != one(p)
		k += 1
	end
	return k
end
	

SymGroup(S::Set) = [ (Perm(t)) for t in Sym(S::Set) ]

# Conceptually I dont like these. SymGroup() is not 
# Note: Sym calls unique!() on Arrays
SymGroup(S::Array) = [ (Perm(t)) for t in Sym(S) ]

function NormalSubgroup(H::Set{Perm})
	return all(a-> (H^a == H), Sym(first(H).A))
end


#= STILL NEED THE FOLLOWING

*(::Perm, ::Set{Perm}) & vise versa... do i want this though?
Tensor Representation, Should prolly go in SymGrpAndReps

Sym group idea
struct SymGroup{T} where T <: Symbol
	S::Set{T}
	D Set( (Perm(t)) for t in Sym(S) )
end
=#

function SgnDecomp(G::Array{Perm{T},1}) where T
	Decomp = Decomp([-1 => Perm{T}[], 1 => Perm{T}[] ])
	for p in G
		push!(Decomp[sgn(p)], p)
	end
	return Decomp
end

function FixedDecomp(G::Array{Perm{T},1}) where T
	X = one(G[1])
	Decomp = [Perm{T}[] for i in 1:(length(X.A)+1)]
	for p in G
		hasfixed = false
		for i in X.A
			if p(i) == i
				push!(Decomp[i], p)
				hasfixed = true
			end
		end
		if hasfixed == false push!(Decomp[end], p) end
	end
	return Decomp
end

function PeriodDecomp(G::Array{Perm{T}}) where T
	Decomp = Dict{Integer,Array{Perm{T},1}}()
	maxPeriod = 0
	for p in G
		z = period(p)
		if z in keys(Decomp) push!(Decomp[z], p)
		else Decomp[z] = [p]
		end
	end
	return Decomp
end