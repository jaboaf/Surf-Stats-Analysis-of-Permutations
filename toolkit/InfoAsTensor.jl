using JSON: parse
include("SymGrpAndReps.jl")
include("OrderingUtils.jl")

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

WaveIds = sort(collect( filter( wid-> data[wid]["nJudOrigs"] == 5 & data[wid]["nSubScos"], keys(data)) ))
HeatIds = sort( unique( map( wid-> data[wid]["heatId"], WaveIds )))
RndIds = sort( unique( map( wid-> data[wid]["rndId"], WaveIds )))
EvtIds = sort( unique( map( wid-> data[wid]["evtId"], WaveIds )))

maxRnd = Dict()
for id in EvtIds
    thisEvt = filter(wid->data[wid]["evtId"]==id, WaveIds )
    maxRnd[id] = maximum( map( wid -> Base.parse(Int, data[wid]["rnd"] ), thisEvt))
end

EvtYears = sort( unique( map( wid-> data[wid]["evtYear"], WaveIds )))
EvtOrigs = sort( unique( map( wid-> isoDict[ data[wid]["evtOrig"] ], WaveIds)))
#EvtNames = sort( unique( map( wid-> data[wid]["evtName"], WaveIds)))
Rnds = [i for i in 1:maximum(values(maxRnd))]
#Heats = sort( unique( map( wid-> Base.parse(Int, data[wid]["heat"]), WaveIds)))
AthOrigs = sort( unique( map( wid-> isoDict[ data[wid]["athOrig"] ], WaveIds)))
Imatch = sort( unique( map( wid-> data[wid]["athOrig"] in data[wid]["subScoOrig"] ? 1 : 2, WaveIds)))

evtYrsIND = Dict([yr=>i for (i,yr) in enumerate(EvtYears)])
evtOrigsIND = Dict([eo=>i for (i,eo) in enumerate(EvtOrigs)])
athOrigsIND = Dict([ao=>i for (i,ao) in enumerate(AthOrigs)])

JustPanelOrigs = sort( map(c-> isoDict[c], collect( union( map( wid-> data[wid]["subScoOrig"], WaveIds)...))))
JustPanelToInd = Dict([orig => i for (i,orig) in enumerate(JustPanelOrigs)])

# basically constructing "concordant permutations" for each wave
function GetPermsIntIndex(wid)
	Perms = Array{Int,1}[]
	subScoOrigs = map( c-> isoDict[c], data[wid]["subScoOrig"])
	origXscos = zip(subScoOrigs, data[wid]["subSco"])
	labeledScos = Dict([c => Float16[] for c in unique(subScoOrigs)])
	for p in origXscos
		push!( labeledScos[ p[1] ], p[2])
	end
	for G in Sym(JustPanelOrigs)
		A = Array{Float16,1}[]
		for c in G if c in subScoOrigs push!(A, labeledScos[c]) end end 
		if isordered(A)
			push!(Perms, [JustPanelToInd[c] for c in G])
		end
	end
	return Perms
end

Data = zeros(Float16, (3,7,8,10,2,7,7,7,7,7,7,7) );

for htid in HeatIds
	thisHeatsIds = sort( filter( wid-> data[wid]["heatId"] == htid ,WaveIds))
	panelOrigs = sort(unique( map(x->isoDict[x], data[first(thisHeatsIds)]["subScoOrig"] )))
	M = maxRnd[ data[first(thisHeatsIds)]["evtId"] ]
	HT = Base.parse(Int, data[first(thisHeatsIds)]["heat"])
	for wid in thisHeatsIds
		evtYear = data[wid]["evtYear"]
		evtOrig = isoDict[ data[wid]["evtOrig"] ]
		evtName = data[wid]["evtName"]
		rnd = M - Base.parse(Int, data[wid]["rnd"])
		ht = HT
		athOrig = isoDict[ data[wid]["athOrig"] ]
		Imatch = (data[wid]["athOrig"] in panelOrigs ? 1 : 2)

		someINDS = [ evtYrsIND[evtYear], evtOrigsIND[evtOrig], rnd, athOrigsIND[athOrig], Imatch]
		perms = GetPermsIntIndex(wid)
		if length(perms) != 0
			for P in perms
				IND = [someINDS... P...]
				Data[ IND... ] = 1/length(perms)
			end
		else
			Data[someINDS...,:,:,:,:,:,:,:] .= 1/factorial(7)
		end
	end
end