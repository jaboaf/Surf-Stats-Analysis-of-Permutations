include("SimplifyComp.jl")

# JOJO ADD ATHLETE ORIGIN TO THIS!

Hts = map(x->x[1]=>map(y->(round(mean(y.judge_scores),digits=2),y.λ_origs),x[2]),partitionBy(:heatId));
sort!(Hts)
# Heat x Mean judge score x panel
# since judge score is rounded to 2 decimal places, sco range needs to be 1000
e_(i::Int) = vec([j==i ? 1 : 0 for j in 1:7])
Mi = [zeros(7) for r in 1:5]
Mu = [zeros(7) for r in 1:5]
M = [zeros(7) for r in 1:5]

H = Array{Float32,7}(undef,(910,1000,7,7,7,7,7));
for (i,ht) in enumerate(Hts)
	for (s,λ) in ht[2]
		sco = Int(round(s*100))
		c = prod(factorial.(length.(λ)))
		b = [1/factorial(length(λ[k])) for k in 1:length(λ) for r in 1:length(λ[k])]
		for a in Base.product( Sym.(λ)...)
			I = Int.(vcat(a...))
			M .+= e_.(I)
			Mu .+= (1/c)^(1/5) .* e_.(I)
			Mi .+=  b .* e_.(I)
			H[i,sco, I... ] += 1/c
		end
	end
end


# Is this the dual of M ?
K = [[M[r][c] for r in 1:5] for c in 1:7]
Ki = [[Mi[r][c] for r in 1:5] for c in 1:7]
Ku = [[Mu[r][c] for r in 1:5] for c in 1:7]

# C for "by Country"
b(a,b) = length(a)==length(b) ? sum([x*y for x in a for y in b]) : error("array do not have the same length")
C = sum([M[i]⊗M[j] / (sum(M[i])*sum(M[j])) for i in 1:5 for j in 1:5])
# R for "by Rank"
R = sum([K[i]⊗K[j] / (sum(K[i])*sum(K[j])) for i in 1:7 for j in 1:7])

# simply panel data prob dist.
P = E(H,3,4,5,6,7) / sum(H);
# probability vector variate
F_m = mapreduce(v->v/sum(v),⊗,M);
F_mi = mapreduce(v->v/sum(v),⊗,Mi);
F_mu = mapreduce(v->v/sum(v),⊗,Mu);



#=
#det(C)
#det(R)
println("Approximately 0...")

sum(C)
sum(R)
println("Interesting hmmmmm")

C_n= C/25
R_n= R/49

U_C = [1/7*1/7 for c in 1:7,c2 in 1:7]
U_R = [1/5*1/5 for r in 1:5,r2 in 1:5]

A = [M[r][c] for r in 1:5, c in 1:7]

ht1 = map(p-> sum([Rep(vcat(t...)) for t in Base.product(Sym.(map(cls->Int.(cls),p[2].λ_origs))...)]), filter(w->w[2].heatId=="79998",WAVES) )
ht1P = sum(map(p->p/(sum(p)/5),ht1)) # 5 is number of panel origs
# note: ht1P above, row and column sums are == 26 which is len of heat
ht1P = sum(map(p->p/(sum(p)/5),ht1))/length(ht1)

ht2 = map(p-> sum([Rep(vcat(t...)) for t in Base.product(Sym.(map(cls->Int.(cls),p[2].λ_origs))...)]), filter(w->w[2].heatId=="68816",WAVES) )
ht2P = sum(map(p->p/(sum(p)/5),ht2))
# 5 is number of panel origs, 4 is number of origs
# note: ht2P above, column sums are == 17 which is len of heat
# note: ht2P above, row sums are 17 for every row except 1st
# 1st row corresponds to AUS which two jugdes have orig of
ht2P = sum(map(p->p/(sum(p)/5),ht2))/length(ht2)

now see:
sum(ht1P,dims=(1))
sum(ht1P,dims=(2))
sum(ht2P,dims=(1))
sum(ht2P,dims=(2))

also this is what you'd expect:
sum(ht1P*ht1P',dims=(1))
sum(ht1P*ht1P',dims=(2))

now you'd expect this (because there are 2 judges from AUS):
sum(ht2P*ht2P',dims=(2))
would you expect this? (i didn't at frist but it makes sense cause weve made it symmetric):
sum(ht2P*ht2P',dims=(1))
=#