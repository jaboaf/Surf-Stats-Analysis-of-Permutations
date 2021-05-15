Cumulants = []

Cumulants[1] = [ E(P,i) for i in 1:5 ]
Cumulants[2] = [ E(P,(i,j))-E(P,i)⊗E(P,j) for i in 1:5,j in 1:5]
Cumulants[3] = [ E(P,(i,j,k))-(E(P,i)⊗E(P,(j,k))+E(P,j)⊗E(P,(i,k))+E(P,k)⊗E(P,(i,j)))+2*E(P,i)⊗E(P,j)⊗E(P,k) for i in 1:5,j in 1:5, k in 1:5]
Cumulants[4] = 


Ydensity = map(y->y/sum(y),Y)
devs = map(y-> y - SymOp(y),Ydensity)
plot(map(y->sum(abs.(y)),devs))

# theses are the same
plot(map(y->sum(abs.(filter(>(0),y))),devs))
oplot(map(y->sum(abs.(filter(<(0),y))),devs))
#... these are bdd. by 2
# so we have 
plot(map(y->1/2*sum(abs.(y)),devs))

#Make these densities graded (let class be dim index) and then apply sym op to each part of the graded space
# measure abs. dev see what happens



# Get Writing2.jl part

