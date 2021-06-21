include("SimplifyComp.jl")
using GRUtils

#= Description
L = series of observed panels
Q = series of observed panels taken under equivalence

DLm = series of, sub-series of L with panel composition = β.
DQm = series of, sub-series of Q with panel composition = β.

Note: If we retain the panel composition information, β, (and have some
consistent ordering of countries) we can take each panel in Vₖˣ⁵ where
k is the number of distinct countries on the panel, i.e. |supp(β)|

Gₖ is the symmetric group on k letters
S is the symmetrization operator, F ↦ 1/(n!)∑ᵧ γF where γ ranges G_ndims(F)
A is the anti-symmetrization operator , F ↦ 1/(n!)∑ᵧ sgn(γ)γF

For each b, DLmᵦ, is a series of panels with panel composition β. I am interested in 2 statistics.
nsᵦ := (Yᵢ - S(Fᵦ) ) where Yᵢ ranges over elements of DLmᵦ
aᵦ := A(Yᵢ) where Yi ranges over elements of  DLm=b
and Fᵦ = 1/#Fᵦ ∑(DLmᵦ) is the empirical distribution
=#


W = last.(WAVES);
M = [ w.m_c for w in W];
L = [ w.λ_c for w in W]; # observed panel
Q = [ w.λ_eqc for w in W]; # observed panel under equivalence

DLm = [ m => map(x->embedd(x;minimal=true),L[findall(==(m),M)]) for m in sort(unique(M))] ;
DQm = [ m => map(x->embedd(x;minimal=true),Q[findall(==(m),M)]) for m in sort(unique(M))] ;
#D = sort(unique(M)) .=> [ map(x->embedd(x;minimal=true),L[findall(==(m),M)]) .- SymOp(sum(L[findall(==(m),M)])) for m in sort(unique(M))] ;
#=
U = Pair{NTuple{7,Int64},Array{Array{Float64,1},5}}[]
for d in DLₘ[1:2]
	Nd = length(d[2])
	S = SymOp(sum(d[2]))
	A = map(x-> x-S, d[2])
	R = vcat.( A...)
	push!(U, d[1] => R)
end
=#
Lns = Pair{NTuple{7,Int64},Array{Float64,1}}[]
La = Pair{NTuple{7,Int64},Array{Float64,1}}[]
for d in DLm
	Nd = length(d[2])
	S = SymOp(sum(d[2])/Nd)
	nonS = map(x-> 1/2*sum(abs.(x-S)), d[2])
	push!(Lns, d[1] => nonS)
	A = map(x->1/2*sum(abs.(AltOp(x))),d[2])
	push!(La, d[1] => A)
end

Qns = Pair{NTuple{7,Int64},Array{Float64,1}}[]
Qa = Pair{NTuple{7,Int64},Array{Float64,1}}[]
for d in DQm
	Nd = length(d[2])
	S = SymOp(sum(d[2])/Nd)
	nonS = map(x-> 1/2*sum(abs.(x-S)), d[2])
	push!(Qns, d[1] => nonS)
	A = map(x->1/2*sum(abs.(AltOp(x))),d[2])
	push!(Qa, d[1] => A)
end

L_b = [ w.λ_c for w in W]; # observed panel
Q_b = [ w.λ_eqc for w in W]; # observed panel under equivalence

function pltAa(Aa::Array{Array{Float64,1}};inInterval=false,tlab="",xlab="",ylab="")
	M = maximum(maximum.(Aa))
	m = minimum(minimum.(Aa))
	Figure()
	if inInterval ==true
		for a in Aa plot(LinRange(0,1,length(a)),a,hold=true,ylim=(m,M)) end
	else
		for a in Aa plot(a,hold=true,ylim=(m,M)) end
	end
	title(tlab)
	xlabel(xlab)
	ylabel(ylab)
	gcf()
end
# W[n] = 1/√n ∑ⁿ X[i]
partsum(X::Array{T,1} where T) = sqrt.(1:length(X)) .^-1 .* cumsum(X)
# W[n] = 1/√n ∑ⁿ (X[i]-μₓ)/σₓ
partsumProc(X::Array{T,1} where T) = sqrt.(1:length(X)) .^-1 .* cumsum( (X .- mean(X)) ./ std(X) )
# W[n] = 1/√n ∑ⁿ X[i]-μₓ
centsumProc(X::Array{T,1} where T) = sqrt.(1:length(X)) .^-1 .* cumsum( X .- mean(X))
# W[n] = 1/n ∑ⁿ X[i]
asymptoticmean(X::Array{T,1} where T) = collect(1:length(X)) .^-1 .* cumsum(X)

# Distributions of Lns and Qns
pltAa( last.(Lns);tlab="Non-Symmetric Part of L",xlab="observation",ylab="Lns")
pltAa( last.(Qns);tlab="Non-Symmetric Part of Q",xlab="observation",ylab="Qns" )
# Distributions of Lns and Qns embedded into [0,1]
pltAa( last.(Lns);inInterval=true, tlab="Non-Symmetric Part of L embedded into [0,1]",xlab="observation/N_β",ylab="Lns")
pltAa( last.(Qns);inInterval=true, tlab="Non-Symmetric Part of Q embedded into [0,1]",xlab="observation/N_β",ylab="Qns")
# Empirical CDF of empirical means
plot( sort(sum.(last.(Lns)) ./ length.(last.(Lns))),ylim=(0,1), title="Sorted Empricial Means of Lns")
plot( sort(sum.(last.(Qns)) ./ length.(last.(Qns))),ylim=(0,1), title="Sorted Empricial Means of Qns")

# Distributions of La and Qa
pltAa( last.(La);tlab="Anti-symmetric Part of L",xlab="observation",ylab="La")
pltAa( last.(Qa);tlab="Anti-symmetric Part of Q",xlab="observation",ylab="Qa")
# Distributions of Lns and Qns embedded into [0,1]
pltAa( last.(La);inInterval=true, tlab="Anti-symmetric Part of L embedded into [0,1]",xlab="observation/N_β",ylab="La")
pltAa( last.(Qa);inInterval=true, tlab="Anti-symmetric Part of Q embedded into [0,1]",xlab="observation/N_β",ylab="Qa")
# Empirical CDF of empirical means
plot( sort(sum.(last.(La)) ./ length.(last.(La))),ylim=(0,1), title="Sorted Empricial Means of La")
plot( sort(sum.(last.(Qa)) ./ length.(last.(Qa))),ylim=(0,1), title="Sorted Empricial Means of Qa")

# now for some visuals of cumulative sums
pltAa( cumsum.(last.(Lns));tlab ="Cumulative Sums of Lns_{β} for each β")
pltAa( cumsum.(last.(Qns));tlab ="Cumulative Sums of Wns_{β} for each β")

pltAa( cumsum.(last.(Lns));inInterval=true,tlab ="Cumulative Sums of Lns_{β} for each β, embedded in [0,1]")
pltAa( cumsum.(last.(Qns));inInterval=true,tlab ="Cumulative Sums of Wns_{β} for each β, embedded in [0,1]")

pltAa( cumsum.(last.(Qa));tlab ="Cumulative Sums of Qa_{β} for each β (same as La_{β})")
pltAa( cumsum.(last.(Qa));inInterval=true,tlab ="Cumulative Sums of Qa_{β} for each β (same as La_{β}), embedded in [0,1]")

pltAa( partsum.(last.(Lns));tlab ="Partial Sums: W_{β}[n] = (√n)^{-1} ∑^{n}_{i} Lns_{β}[i]")
pltAa( partsum.(last.(Qns));tlab ="Partial Sums: W_{β}[n] = (√n)^{-1} ∑^{n}_{i} Qns_{β}[i]")
pltAa( partsum.(last.(Lns));inInterval=true,tlab ="Partial Sums: W_{β} (t) = (√floor (tN_{β}))^{-1} ∑^{floor (tN_{β})}_{i} Lns_{β}[i]")
pltAa( partsum.(last.(Qns));inInterval=true,tlab ="Partial Sums: W_{β} (t) = (√floor (tN_{β}))^{-1} ∑^{floor (tN_{β})}_{i} Qns_{β}[i]")

pltAa( partsum.(last.(Qa));tlab ="Partial Sums: W_{β}[n] = (√n)^{-1} ∑^{n}_{i} Qa_{β}[i]")
pltAa( partsum.(last.(Qa));inInterval=true,tlab="Partial Sums: W_{β} (t) = (√floor (tN_{β}))^{-1} ∑^{floor (tN_{β})}_{i} Qa_{β}[i]")

pltAa( partsumProc.(last.(Lns));tlab="Partial Sum Process: W_{β}[n] = (√n)^{-1} ∑^{n} (Lns_{β}[i]-μ̂_{β})σ_{β}^{-1}")
pltAa( partsumProc.(last.(Qns));tlab="Partial Sum Process: W_{β}[n] = (√n)^{-1} ∑^{n} (Qns_{β}[i]-μ̂_{β})σ_{β}^{-1}")
pltAa( partsumProc.(last.(Lns));inInterval=true,tlab="Partial Sum Process: W_{β} (t) = (√floor (tN_{β}))^{-1} ∑^{floor (tN_{β})}_{i} (Lns_{β}[i]-μ̂_{β})σ_{β}^{-1}")
pltAa( partsumProc.(last.(Qns));inInterval=true,tlab="Partial Sum Process: W_{β} (t) = (√floor (tN_{β}))^{-1} ∑^{floor (tN_{β})}_{i} (Qns_{β}[i]-μ̂_{β})σ_{β}^{-1}")

pltAa( partsumProc.(filter(q->std(q)!=0,last.(Qa)));tlab="Partial Sum Process (for β where σ_{β} ≠ 0): W_{β}[n] = (√n)^{-1} ∑^{n} (Qa_{β}[i]-μ̂_{β})σ_{β}^{-1}")
pltAa( partsumProc.(filter(q->std(q)!=0,last.(Qa)));inInterval=true,tlab="Partial Sum Process (for β where σ_{β} ≠ 0): W_{β} (t) = (√floor (tN_{β}))^{-1} ∑^{floor (tN_{β})}_{i} (Qa_{β}[i]-μ̂_{β})σ_{β}^{-1}")

pltAa( centsumProc.(last.(Qa));tlab="Centered Sum Process: W_{β}[n] = (√n)^{-1} ∑^{n} (Qa_{β}[i]-μ̂_{β})")
pltAa( centsumProc.(last.(Qa));inInterval=true,tlab="Centered Sum Process: W_{β} (t) = (√floor (tN_{β}))^{-1} ∑^{floor (tN_{β})}_{i} (Qa_{β}[i]-μ̂_{β})")

histogram(std.(last.(Lns)),title="Histogram of Standard Deviations of Lns_{β} over all β",xlim=(0,0.6))
histogram(std.(last.(Qns)),title="Histogram of Standard Deviations of Qns_{β} over all β",xlim=(0,0.6))
histogram(std.(last.(Qa)),title="Histogram of Standard Deviations of Qa_{β} over all β",xlim=(0,0.6))

pltAa( asymptoticmean.(last.(Lns));tlab="Asymptotic Mean: W_{β}[n] = (n)^{-1} ∑^{n} Lns_{β}[i]")
pltAa( asymptoticmean.(last.(Qns));tlab="Asymptotic Mean: W_{β}[n] = (n)^{-1} ∑^{n} Qns_{β}[i]")
pltAa( asymptoticmean.(last.(Lns));inInterval=true,tlab="Asymptotic Mean: W_{β} (t) = (floor (tN_{β}))^{-1} ∑^{floor (tN_{β})}_{i} Lns_{β}[i]")
pltAa( asymptoticmean.(last.(Qns));inInterval=true,tlab="Asymptotic Mean: W_{β} (t) = (floor (tN_{β}))^{-1} ∑^{floor (tN_{β})}_{i} Qns_{β}[i]")

pltAa( asymptoticmean.(last.(Qa));tlab="Asymptotic Mean: W_{β}[n] = (n)^{-1} ∑^{n} Qa_{β}[i]")
pltAa( asymptoticmean.(last.(Qa));inInterval=true,tlab="Asymptotic Mean: W_{β} (t) = (floor (tN_{β}))^{-1} ∑^{floor (tN_{β})}_{i} Qa_{β}[i]")

𝔸 = [ q[1] => map(n-> 1/2*sum(abs.(AltOp( sum(q[2][1:n])/n ))),1:length(q[2])) for q in DQm];

pltAa( last.(𝔸);tlab=" W_{β}[n] = AntiSym (𝔽_{β,n}) = A (𝔽_{β,n}) ")


Figure()
for u in Uns plot(u[2],hold=true) end
gcf()

Figure()
for u in U plot(sqrt.(1:length(u[2])) .* u[2],hold=true) end
gcf()

Figure()
for u in U plot(u[2] ./ sqrt.(1:length(u[2])),hold=true) end
gcf()

for u in U plot( u[2] .^(-1),hold=true) end

for u in U plot(cumsum(u[2]) ./ collect(1:length(u[2])),hold=true) end
for u in U plot(LinRange(0,1,length(u[2])),cumsum(u[2]) ./ collect(1:length(u[2])),hold=true) end
for u in U plot( cumsum(u[2]) ./ (collect(1:length(u[2])) .^(3/2)),hold=true) end

Msupp = [ w.m_c_supp for w in W];
Dsupp = [m => map(x->embedd(x;minimal=true),L[findall(==(m),Msupp)]) for m in sort(unique(Msupp))] ;
Usupp = Pair{NTuple{7,Int64},Array{Float64,1}}[]
for d in Dsupp
	Nd = length(d[2])
	S = SymOp(sum(d[2])/Nd)
	A = map(x-> 1/2*sum(abs.(x-S)), d[2])
	push!(Usupp, d[1] => A)
end
for u in Usupp plot(u[2],hold=true) end
Figure()
for u in Usupp plot( 1 .+ u[2] .^(-1),hold=true) end
gcf()


include("PanelAnalysis.jl")


Hts = map(x->map(y->embedd(y.λ_c;minimal=true),x),last.(partitionBy((:heatId,))));
Devs = map(Hts) do H
k = length(H)
a = 1/2*sum(abs.(AltOp(sum(H)/k)))
end



