# use JSON for parsing, Use Base.parse for parsing purposes
using JSON: parse
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
# Simplify data to most "complete years"
filter!(data) do wave
	wave[2]["evtYear"] in ["2018","2019"]
end
# subScoOrigDefect== true when no info on judge origins
filter!(data) do wave
	wave[2]["subScoOrigDefect"]==false
end

WaveIds = sort(collect(keys(data)))
HeatIds = sort(unique(map( wid-> data[wid]["heatId"], WaveIds )))
RndIds = sort(unique(map( wid-> data[wid]["rndId"], WaveIds )))
EvtIds = sort(unique(map( wid-> data[wid]["evtId"], WaveIds )))

maxRnd = Dict()
for id in EvtIds
    thisEvt = filter(wid->data[wid]["evtId"]==id, WaveIds )
    maxRnd[id] = maximum(map(wid->Base.parse(Int,data[wid]["rnd"]),thisEvt))
end

for wave in data
	wave[2]["subScoOrig"] = map(x->isoDict[x], wave[2]["subScoOrig"])
	wave[2]["athOrig"] = isoDict[ wave[2]["athOrig"] ]
	wave[2]["evtOrig"] = isoDict[ wave[2]["evtOrig"] ]
	wave[2]["rnd"] = maxRnd[wave[2]["evtId"]]-Base.parse(Int,wave[2]["rnd"])+1
	wave[2]["heat"] = Base.parse(Int, wave[2]["heat"])
end

function varRng(s::String)
	rng=Set()
	for wave in data
		push!(rng, wave[2][s])
	end
	return sort(collect(rng))
end

ToInd = Dict()
# YR has dim 2
ToInd[:YR] = Dict([y=>i for (i,y) in enumerate(varRng("evtYear"))] )
# EVT has dim 10
ToInd[:EVT] = Dict([n=>i for (i,n) in enumerate(varRng("evtName"))] )
# RND has dim 7
ToInd[:RND] = Dict([r=>i for (i,r) in enumerate(varRng("rnd"))] )
# HEAT has dim 16
ToInd[:HEAT] = Dict([h=>i for (i,h) in enumerate(varRng("heat"))] )
# ATH_orig das dim 10
ToInd[:ATH_orig] = Dict([r=>i for (i,r) in enumerate(varRng("athOrig"))] )
# JUD_orig has dim 7
ToInd[:JUD_orig] = Dict([orig=>i for (i,orig) in enumerate(sort(union(varRng("subScoOrig")...)))])
# RANK has dim 5
ToInd[:RANK] = Dict([i=>i for i in 1:5])

vars = [:YR, :EVT, :RND, :HEAT, :ATH_orig, :JUD_orig, :RANK ]
info = zeros(Int8,(2,10,7,16,10,7,5));
for wave in data
	yr = ToInd[:YR][ wave[2]["evtYear"] ]
	evt = ToInd[:EVT][ wave[2]["evtName"] ]
	rnd = ToInd[:RND][ wave[2]["rnd"] ]
	heat = ToInd[:HEAT][ wave[2]["heat"] ]
	ath_orig = ToInd[:ATH_orig][ wave[2]["athOrig"] ]

	judge_origs = wave[2]["subScoOrig"]
	Scores = Float16.(wave[2]["subSco"])
	for (i,score) in enumerate(sort(unique(Scores)))
		I = findall(x->x==score, Scores)
		rank = Int(i) # to reassure julia's type inference
		for orig in judge_origs[I]
			jud_orig = Int(ToInd[:JUD_orig][ orig ]) # reassure type inference
			info[ yr, evt , rnd , ath_orig , jud_orig , rank ] += 1
		end
	end
end