include("PanelAnalysis.jl")

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

Yd = map(y->y/sum(y),Y)

W = sum.(Y)

# Total Variation
TV = map(y-> y-SymOp(y) ,Y);
ATV = map(y-> 1/2*abs.(y-SymOp(y)) ,Yd);

BV = [abs.(P-SymOp(P))];


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

