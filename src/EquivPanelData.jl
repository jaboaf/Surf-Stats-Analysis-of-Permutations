include("toolkit/OrderingUtils.jl")
isoDict = Dict([
    "Australia" => :AUS,
    "Basque Country" => :ESP,
    "Brazil" => :BRA,
    "Fiji" => :FJI,
    "France" => :FRA,
    "Hawaii" => :USA,
    "Indonesia" => :IDN,
    "Italy" => :ITA,
    "Japan" => :JPN,
    "New Zealand" => :NZL,
    "Portugal" => :PRT,
    "South Africa" => :ZAF,
    "Spain" => :ESP,
    "United States" => :USA
])

data = parse( open("Data/CleanAllDataCC.txt", "r"))

filter!(data) do wave
	# Simplify data to most "complete years"
	# ... we can always undo this
	# subScoOrigDefect== true when no info on judge origins
	wave[2]["evtYear"] in ["2018","2019"] && wave[2]["subScoOrigDefect"]==false
end

eqPanels = Array{Pair,1}[]
for wave in data
	judge_origs = map(x->isoDict[x], wave[2]["subScoOrig"])
	panelOrigs = Set(judge_origs)
	judge_scores = Float16.(wave[2]["subSco"])
	eqPanel = Pair{Array{Float16,1},Array{Symbol,1}}[]
	for c in panelOrigs
		I = findall(==(c), judge_origs)
		push!(eqPanel, sort(judge_scores[I])=>judge_origs[I])
	end
	sort!(eqPanel)

	if !(isordered(first.(eqPanel)))
		cond = true
		while cond
			j = findfirst(i->!(isordered(first.(eqPanel)[i:(i+1)])), 1:(length(eqPanel)-1))
			newElm = [ vcat(first.(eqPanel)[j:(j+1)]...) => vcat(last.(eqPanel)[j:(j+1)]...)]
			eqPanel = [ eqPanel[1:(j-1)]... , newElm... , eqPanel[(j+2):end]... ]
			cond = eval( !(isordered( collect(first.(eqPanel)) )) )
		end
	end
	push!(eqPanels, eqPanel)
end

