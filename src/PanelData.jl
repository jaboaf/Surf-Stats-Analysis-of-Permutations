using GRUtils
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

data = parse( open("../Data/CleanAllDataCC.txt", "r"))

filter!(data) do wave
	# Simplify data to most "complete years"
	# ... we can always undo this
	# subScoOrigDefect== true when no info on judge origins
	wave[2]["evtYear"] in ["2018","2019"] && wave[2]["subScoOrigDefect"]==false
end

panels = Array{Array{Symbol,1},1}[]
for wave in data
	judge_origs = map(x->isoDict[x], wave[2]["subScoOrig"])
	judge_scores = Float16.(wave[2]["subSco"])
	judging = Pair(judge_scores, judge_origs)
	panel = []
	for (i,s) in enumerate(sort(unique(judge_scores)))
		I = findall(x->x[1]==s, judge_scores)
		push!(panel, judge_origs[I])
	end
	push!(panels, panel)
end


