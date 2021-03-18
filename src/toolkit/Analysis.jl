include("SurfingInfo.jl"); println("Loaded SurfingInfo")
include("FunctionalSym.jl"); println("Loaded FunctionalSym")
include("OrderingUtils.jl");println("Loaded OrderingUtils")

B = Set([:Match, :NoMatch])
Sym_B = Sym(B)
B_panels = map(w->w.panelBinary, waves)

B_valids = Array{Array{Symbol,1},1}[]
for panel in B_panels
	concordant = Array{Symbol,1}[]
	for p in Sym_B
		subp = filter(z->z in keys(panel), p) 
		if isordered([ panel[c] for c in subp ])
			push!(concordant, p)
		end
	end
	push!(B_valids, concordant)
end

println(count(x-> length(x) ==0, B_valids))
println(count(x-> length(x) ==1, B_valids))
println(count(x-> length(x) ==2, B_valids))



C = union(map(w->w.panelOrigs, waves)...)
CtoI = Dict([c=>i for (i,c) in enumerate(sort(collect(C)))])
C_panels = map(w-> Dict([ Int8(CtoI[c])=>w.panel[c] for c in keys(w.panel)]), waves)
Sym_C = Sym(length(C))

Concord_Indicator = map(C_panels) do panel
	for p in Sym( Set(keys(panel)) )
		if isordered([ panel[c] for c in p ])
			return true
		end
	end
	return false
end

Concordant_panels = C_panels[ Concord_Indicator ]
global i = 0
for panel in Concordant_panels
	concordant = Array{Int8,1}[]
	sub = collect(keys(panel))
	for p in Sym_C
		subp = filter(z-> z in sub, p) 
		if isordered([ panel[c] for c in subp ])
			push!(concordant, p)
		end
	end
	if length(concordant) > 0 push!(C_maps, sum(Rep.(concordant)) )
	else push!(C_maps, zeros(Int8,7,7)) 
	end
	push!(C_prjs, sub)
	global i += 1
	if i % 100 == 0 print(i) end
end

C_sgn = Array{BitArray{2},1}[]
C_fixed = Array{BitArray{2},1}[]
C_period =  Array{BitArray{2},1}[]

C_maps = Array{<:Number,2}[]
C_prjs = Array{<:Number,1}[]

#=
for panel in C_panels
	concordant = Array{Symbol,1}[]
	for p in Sym_C
		subp = filter(z->z in keys(panel), p) 
		if isordered([ panel[c] for c in subp ])
			push!(concordant, p)
		end
	end
	push!(C_valids, concordant)
end
=#

