include("PanelAnalysis.jl")
using GRUtils

P = zeros(Rational, (7,7,7,7,7) );
for λ in Panel_λs
	c = prod( factorial.(length.(λ)) )
	for a in Base.product( Sym.(λ)...)
		P[ Int.(vcat(a...))... ] += 1//c
	end
end

ℙ = P/sum(P);

V = Array{Array{Rational,1}}(undef,(7,7,7,7,7));
W = Array{Array{Rational,1}}(undef,(7,7,7,7,7));
for i in eachindex(V) V[i] = Rational[] end;
V = vcat.(V,1//factorial(5) .* P);

for g in setdiff(Sym(Int8[1,2,3,4,5]),[[1,2,3,4,5]])
	V = vcat.(V, sgn(g)*1//factorial(5) .* permutedims(P,g))
end
W = map(v->abs.(v),V); 

V ./= sum(P);

plot(hcat(vec(V)...))
plot(hcat(sort.(vec(V))...))
plot(hcat(sort.(vec(V))...)*ones(length(V)))
plot(hcat(vec(V)...)*ones(length(V)))
plot(hcat(map(s-> sort(s .- mean(s)),vec(W))...))

