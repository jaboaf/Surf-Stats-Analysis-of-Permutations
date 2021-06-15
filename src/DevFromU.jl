include("SimplifyComp.jl")
using GRUtils

#=(:id, :szn, :evtOrig, :evt, :evtId, :rnd, :rndId,
:heat, :heatId, :athName, :athId, :athOrig,
:currentPoints, :endingPoints,
:judge_scores, :judge_origs, :panel, :panelBinary,
:labeledPanel, :labeledPanelBinary, :eqPanel,
:λ_origs, :λ_binary, :λ_c, :λ_eqc, :m_c, :m_b, :panel_scores, :panel_origs, :I_match
=#

# add m_c which is just the monomial of ctrys
# add m_binary which is monomial of match, no match
# add m_c_supp which is the support of m_c


C = Set([AUS,BRA,ESP,FRA,PRT,USA,ZAF])

panel_compositions = varRng(:m_c)
unique(sort.(panel_compositions))
unique(sort.(unique.(panel_compositions)))

# theorhetically
unqiue(sort.(collect.(Base.product(C,C,C,C,C))))
unique(sort.(unique.(collect.(Base.product(C,C,C,C,C)))))


P = zeros(Float32, (5,5,5,5,5,5,5,7,7,7,7,7));
for w in last.(WAVES)
	c = prod( factorial.(length.(w.λ_c)) )
	for t in Base.product(Sym.(w.λ_c)...)
		I = CartesianIndex((w.m_c .+1)... , Int.(vcat(t...))...)
		P[ I ] += 1/c
	end
end

plot(sort([ sum(abs.(AltOp( P[(a .+1)...,:,:,:,:,:] ))) / sum( P[(a .+1)...,:,:,:,:,:] ) for a in panel_compositions ]))


processes = []
for cond_panel in partitionBy(:m_c)
	a = cond_panel[1]
	Z = zeros(Float32,(7,7,7,7,7))
	S = map(x->x.λ_c,cond_panel[2])
	U = Float32[]
	for (j,λ) in enumerate(S)
		c = prod( factorial.(length.(λ)) )
		for t in Base.product(Sym.(λ)...)
			Z[Int.(vcat(t...))...] += 1/c
		end
		push!(U, 1/2*sum(abs.(AltOp(Z)))/j )
	end
	plot(U, hold=true, ylim=(0,0.6))
	push!(processes, a => U )
	println("done with $a, of len $(length(S))")
end

# Deviations 
for p in processes
	plot(p[2],hold=true,ylim=(0,0.6))
end

for p in processes
	plot(LinRange(0,1,length(p[2])),p[2],hold=true,ylim=(0,1))
end

for p in processes
	plot(LinRange(0,1,length(p[2])-1),diff(p[2]),hold=true,ylim=(-0.5,0.5))
end

for p in processes
plot(LinRange(0,1,length(p[2])-1),diff(sort(p[2])),hold=true,ylim=(0,0.6))
end

empProcs = map(Y-> Y[1] => sqrt.(1:length(Y[2])) .* Y[2] ,processes)



# Fake Data for comparison
include("toolkit/SymGrpAndReps.jl")
fakeprocs= []
for p in processes
	n = length(p[2])
	a = p[1]
	pan = vcat([ fill(ORIG(i),k) for (i,k) in enumerate(a) ]...)
	S = rand(Sym(5), n)
	Z = zeros(Float32,(7,7,7,7,7))
	U = Float32[]
	for (j,t) in enumerate(S)
		I = pan[t]
		Z[ Int.(I)... ] += 1
		push!(U, 1/2*sum(abs.(AltOp(Z))) )
	end
	plot(U, hold=true, ylim=(0,0.6))
	push!(fakeprocs, a => U )
	println("done with $a, of len $(length(S))")
end






