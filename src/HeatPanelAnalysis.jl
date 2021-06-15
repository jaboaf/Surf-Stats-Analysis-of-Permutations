include("PanelAnalysis.jl")

N = 910
S = SymOp(P) / sum(P);
total = sum(P)
Y = Array{<:Number,5}[]
Hts = map(x->map(y->y.λ_c,x),last.(partitionBy(:heatId)));
for ht_λs in Hts
	Yᵢ = zeros(Float32, (7,7,7,7,7))
	for λ in ht_λs
		c = prod( factorial.(length.(λ)) )
		for a in Base.product( Sym.(λ)...)
			Yᵢ[ Int.(vcat(a...))... ] += 1/c
		end
	end
	push!(Y,Yᵢ)
end

Yd = map(y->y/sum(y),Y)

W = sum.(Y)

# Total Variation
TV = map(y-> y-SymOp(y) ,Y);
ATV = map(y-> 1/2*abs.(y-SymOp(y)) ,Yd);
ATVfromS = map(y-> 1/2*abs.(y-S) ,Yd);
BV = [abs.(P-SymOp(P))];



evtPart = [findall(x->x[1].evt==E,last.(partitionBy(:heatId)) ) for E in varRng(:evt)]
evtEffects = [ 1/2*sum(abs.(Yd[h] - SymOp(sum(Yd[E]))/length(E) )) for E in evtPart for h in E]

rndPart = [findall(x->x[1].rnd==R,last.(partitionBy(:heatId)) ) for R in varRng(:rnd)]
rndEffects = [ 1/2*sum(abs.(Yd[h] - SymOp(sum(Yd[R]))/length(R) )) for R in rndPart for h in R]

sznPart = [findall(x->x[1].szn==CT,last.(partitionBy(:heatId)) ) for CT in varRng(:szn)]
sznEffects = [ 1/2*sum(abs.(Yd[h] - SymOp(sum(Yd[CT]))/length(CT) )) for CT in sznPart for h in CT]

evtOrigPart = [findall(x->x[1].evtOrig==C,last.(partitionBy(:heatId)) ) for C in varRng(:evtOrig)]
evtOrigEffects = [ 1/2*sum(abs.(Yd[h] - SymOp(sum(Yd[C]))/length(C) )) for C in evtOrigPart for h in C]

heatPart = [findall(x->x[1].heat==ht,last.(partitionBy(:heatId)) ) for ht in varRng(:heat)]
heatEffects = [ 1/2*sum(abs.(Yd[h] - SymOp(sum(Yd[H]))/length(H) )) for H in heatPart for h in H]


#=
title("Absolute Total Var for each Wave, Sorted")
plot(1/N : 1/N : 1, sort(sum.(ATV)), ylim=(0,1))
xticks(N*0.2)
xlim(0,N)
ylim(0,1)

Figure()
title("Absolute Total Variation from Conditional Symmetric Estimate")
pctile = (1/N):(1/N):1
plot(pctile,sort(sznEffects),ylim=(0,1), label="Season", hold=true)
plot(pctile,sort(evtOrigEffects), label="Event Location")
plot(pctile,sort(evtEffects), label="Event")
plot(pctile,sort(rndEffects), label="Round #")
plot(pctile,sort(heatEffects), label="Heat #")
plot(pctile,1/6 * (sort(sznEffects) .+ sort(evtEffects) .+ sort(rndEffects) .+ sort(evtOrigEffects) .+ sort(athOrigEffects) .+ sort(heatEffects) ), "x", markersize=.25,label="Avg. of ATVs")
legend()
legend(location=4)
ylabel("Abs. Total Variation")
xlabel("Percentile")

=#

KS

#=
DevTV = map(y-> abs.(y) - SymOp(abs.(y)),TV)

# Abs Total Var

DevATV = map(y-> 1/2*abs.(y - SymOp(y)),ATV);
Dev2ATV = map(y-> 1/2*abs.(y - SymOp(y)),DevATV);
Dev3ATV = map(y-> 1/2*abs.(y - SymOp(y)),Dev2ATV);
Dev4ATV = map(y-> 1/2*abs.(y - SymOp(y)),Dev3ATV);
Dev5ATV = map(y-> 1/2*abs.(y - SymOp(y)),Dev4ATV);
Dev6ATV = map(y-> 1/2*abs.(y - SymOp(y)),Dev5ATV);


DevTVanother = map(y-> y-SymOp(y),TV);
TV2 = map(y-> y - SymOp(y),DevTV)
DevTV2 = map(y-> abs.(y) - SymOp(abs.(y)),TV2)

# Absolute Deviation
plot(map(y->sum(abs.(y)),TV))

plot(map(y->sum(abs.(y)),TV))
oplot(map(y->sum(abs.(y)),DevTV))

#Pct devs
plot(map(y->sum(abs.(y)),1/2 .* TV ./ W ))
plot(map(y->sum(abs.(y)),1/2 .* DevTV ./ W ))



plot(sum.(ATV))
oplot(sum.(DevATV))
oplot(sum.(DevDevATV))
oplot(sum.(Dev3ATV))

plot(sum.(ATV) , ylim=(0,1))
oplot(sum.(DevATV) )
oplot(sum.(Dev2ATV) )
oplot(sum.(Dev3ATV) )
oplot( sum.(Dev4ATV) )
=#

#=
Note:
Taking 1/2^(k+1) * DevkATV ./W gives kth symmetric proportional error


=#

