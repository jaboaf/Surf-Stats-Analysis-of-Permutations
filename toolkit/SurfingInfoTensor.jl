using JSON: parse

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
ISOs = sort(unique(collect(values(isoDict))))

#Reading In Data
data = parse( open("../Data/CleanAllDataCC.txt", "r") )

# This is me fixing the problem that round is currently based on "depth into event" 
# as opposed to "shallowness from final", which would make data across events more comparable
EvtIds = unique!([data[wid]["evtId"] for wid in keys(data)])
maxRnd = Dict()
for id in EvtIds
    thisEvt = filter(wid->data[wid]["evtId"]==id, keys(data) )
    maxRnd[id] = maximum(map( wid -> Base.parse(Int, data[wid]["rnd"] ), collect(thisEvt)))
end

variables = [
	"evtYear", 
    "evtOrig",
    "evtName",
    "rnd",
    "heat",
    "athOrig",
    "I_match",
    "waveInHeat",
    "currentPoints",
    "endingPoints",
    "actualSco"
]
varRange = Dict([ var => [] for var in variables ] )

HeatIds = []
for wid in keys(data)
	push!( varRange[ "evtYear" ] 		, Base.parse(Int, data[wid]["evtYear"] ) )
	push!( varRange[ "evtOrig" ] 		, isoDict[ data[wid][ "evtOrig" ] ] )
	push!( varRange[ "evtName" ] 		, data[wid][ "evtName" ] )
	push!( varRange[ "rnd" ] 			, maxRnd[data[wid]["evtId"]] - Base.parse(Int, data[wid]["rnd"]) )
	push!( varRange[ "heat" ] 			, Base.parse(Int, data[wid][ "heat" ] ) )
	push!( varRange[ "athOrig" ] 		, isoDict[ data[wid]["athOrig" ] ] )
	#push!( varRange[ "currentPoints" ] 	, data[wid]["currentPoints"] )
	#push!( varRange[ "endingPoints" ] 	, data[wid]["endingPoints"] )
	push!( varRange[ "actualSco" ] 		, data[wid]["actualSco"] )
    push!(HeatIds, data[wid]["heatId"])
end

unique!(HeatIds)
for var in variables
	sort!( unique!(varRange[var]) )
end

IND = Dict([ var => Dict([val => i for (i,val) in enumerate(varRange[var])]) for var in keys(varRange) ])

# w.r.t variables (evtYear, evtOrig, evtName, rnd, heat, athOrig, I_match, waveInHt, actualSco,)
MartysTensor = Array{Real}(undef, dims=(3,7,11,1,16,10,2,1,1,1,1,7,7) )

for htid in HeatIds
    thisHeatWids = sort!( collect( filter( wid-> data[wid]["heatId"]=="80700", keys(data) )))
    for (k, id) in enumerate(thisHeatWids)
        W = data[wid]
        indexInto = [
            IND[W["evtYear"]],
            IND[W["evtOrig"]],
            IND[W["evtName"]],
            IND[W["rnd"]],
            IND[W["heat"]],
            IND[W["athOrig"]],
            k,
            IND[W["currentPoints"]],
            IND[W["endingPoints"]],
            IND[W["actualSco"]]
        ]
        origs = unique(W["subScoOrig"])
        labeledScos = Dict([varRange[origin] => Float16[] for origin in origs])
        for perm in Sym(origs)
            if isordered([ labeledScos[p] for p in perm])
                

        MartysTensor[indexInto]
    end
end

MartysTensor[]



