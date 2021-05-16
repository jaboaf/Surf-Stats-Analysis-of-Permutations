include("PanelAnalysis")

N = 910
S = SymOp(P)
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

W = sum.(Y)

# Total Variation
TV = map(y-> y-SymOp(y) ,Y);
DevTV = map(y-> abs.(y) - SymOp(abs.(y)),TV)

# Abs Total Var
ATV = map(y-> abs.(y-SymOp(y)) ,Y);
DevATV = map(y-> abs.(y - SymOp(y)),ATV);
DevDevATV = map(y-> abs.(y - SymOp(y)),DevATV);
Dev3ATV = map(y-> abs.(y - SymOp(y)),DevDevATV);
Dev4ATV = map(y-> abs.(y - SymOp(y)),DevDevATV);


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

plot(1/2 .* sum.(ATV) ./ W)
oplot(1/2 .* sum.(DevATV) ./ W)
oplot(1/2 .* sum.(DevDevATV) ./ W)
oplot(1/2 .* sum.(Dev3ATV) ./ W)
oplot(1/2 .* sum.(Dev4ATV) ./ W)



Dev2TV = map(y-> y - SymOp(y),DevTV)

