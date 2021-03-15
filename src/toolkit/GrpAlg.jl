#=
Idea:
abstract type SymGrp end
(Q: can we do: abstract type SymGrp{Set} end)
(Q: can we do: abstract type SymGrp{Set{T}} end where T)

type PermGrp <: SymGrp
type Perm <: PermGrp
=#


#=
Jojo You have to think about this more.
Q: Do you want this struct for αg, where α in F, g in G?
		OR do you want it for α = ∑α_g g ?

Ans?: the first one is simpler, but 
=#
mutable struct GrpAlgElem{F<:Number, Perm}
	vals::Array{F,0}
	perms::Array{Perm,0}
end

+(a::GrpAlgElem, b::GrpAlgElem) = 

# Group is (S::Set, op::Function) where op: S x S --> S

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

