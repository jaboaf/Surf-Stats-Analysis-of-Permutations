include("SimplifyComp.jl")

M = [ w[2].m_c for w in WAVES];
L = [ w[2].λ_c for w in WAVES];
D = [ (w[2].m_c_supp,w[2].λ_c) for w in WAVES]; 

# Interesting to note:
# W = M .⊗ L gives an error
# ERROR: MethodError: no method matching ⊗(::NTuple{7,Int64}, ::Tuple{Array{ORIG,1},Array{ORIG,1}})
# typeof(M)
# Array{NTuple{7,Int64},1}
# typeof(L)
# Array{Tuple{Array{ORIG,1},Vararg{Array{ORIG,1},N} where N},1}
# It makes sense that this error comes about.

M_b = [ w[2].m_b for w in WAVES];
L_b = [ w[2].λ_b for w in WAVES];

# W_b = M_b .⊗ L_b
# ERROR: MethodError: no method matching ⊗(::Tuple{Int64,Int64}, ::NTuple{4,BitArray{1}})

# i need an embeddeding function
VL_b = map(L_b) do l
Z = zeros((2,2,2,2,2))
c = prod( factorial.(length.(l)) )
for p in Base.product(Sym.(Array.(l))...)
	Z[ 1 .+ Int.(vcat(p...))...] += 1/c
end
return Z
end
Bs = SymOp.(VL_b)
Ba = AltOp.(VL_b)

sum(abs.(sum.(Ba)))
# returns 0...
# so maybe eq version of panel is better
#

# idea!
# let e() be the embedding function in {0,1} SS
e_b(i) = vcat(zeros(Int(i)+1 -1),1,zeros(2-(Int(i)+1))) 

Q_b = map(L_b) do l
map(cls->count(==(i),cls)*e_b(i) for i in unique(cls)],l)
end

NVL_b = map(Q_b) do q
reduce(+,map(p->reduce(⊗,p),Base.product( q...)))
end

# this has a bunch of NaN and Inf as returned values, but this is actually correct
# and its what we want actually right?
pQ_b = map(L_b) do l
Q = [ [ e_b(i)/count(==(i-1),cls) for i in 1:2] for cls in l]
end
# i don't think so actually
# instead lets make a prob measure, this is kinda an unjustified approach tbh
# interestingly: if we do e_b(i)*count(==(i),cls) we get a density which has norm diff from 1 sometimes
# if we normalize, we end up normalizing by prod( prod(counts of type in class) for class in partition)
# which just results in e_b(i_1)⊗...⊗e_b(i_k)
P_b = map(L_b) do l
q = map(cls->[ e_b(i) for i in unique(cls)],l)
reduce(⊗, map(c->SymOp(reduce(⊗,c)),q) )
end;

p_b = map(L_b) do l
q = map(cls->[ e_b(i)*count(==(i),cls)/length(cls) for i in unique(cls)],l)
reduce(+,map(p->reduce(⊗,p),Base.product( q...)))
end;


unique(size.(P_b))


# whoahhhhhhh look!!!! this is what we want~!!!!
unique(size.(NVL_b))





e_c(i) = vcat(zeros(Int(i)-1),1,zeros(7-Int(i))) 
Q = map(L) do l
q = [ [ count(==(c),cls)*e_c(c) for c in unique(cls)] for cls in l]
end;

pQ = map(L) do l
q = [ [ e_c(c)/count(==(c),cls) for c in unique(cls)] for cls in l]
end;

NVL = map(Q) do q
reduce(+,map(p->reduce(⊗,p),Base.product( q...)))
end;




# This isn't technically what were after.
# computationally we'd rather embedd using only the number ctrys there are
# can we define an order (given by enum) and attach info about the support of q in Q
unique(size.(NVL))


e(i,n) = vcat(zeros(Int(i)-1),1,zeros(n-Int(i))) 
QD = map(D) do d
k = sum(Int.(d[1]))
q = [ [ count(==(c),cls)*e(sum(Int.(d[1])[1:Int(c)]),k) for c in unique(cls)] for cls in d[2]]
v = reduce(+,map(p->reduce(⊗,p),Base.product( q...)))
s = zeros((2,2,2,2,2,2,2))
s[1 .+ Int.(d[1])...] = 1
return s⊗v
end;
# this gives a prob measure i think.
PD = reduce(+,QD)


rdQD = ⨁(QD);



# I kinda went 2 steps ahead now that i think about it
# lets take one step back
QmD = map(D) do d
k = length(unique(vcat(d[2]...)))
q = [ [ count(==(c),cls)*e(sum(Int.(d[1])[1:Int(c)]),k) for c in unique(cls)] for cls in d[2]]
v = reduce(+,map(p->reduce(⊗,p),Base.product( q...)))
return v
end;

mD = ⨁(map(q-> q/sum(q), QmD));
size.(mD)
# mmmm i like it!!!
sum(sum.(mD))

info = ()




