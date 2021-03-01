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
data = parse(open("../Data/CleanAllDataCC.txt", "r"))

# This is me fixing the problem that round is currently based on "depth into event" 
# as opposed to "shallowness from final", which would make data across events more comparable
EvtIds = unique!([data[wid]["evtId"] for wid in keys(data)])
maxRnd = Dict()
for id in EvtIds
    thisEvt = filter(wid->data[wid]["evtId"]==id, keys(data) )
    maxRnd[id] = maximum(map( wid -> parse(data[wid]["rnd"] ), collect(thisEvt)))
end

variables = [
	"evtYear", 
    "evtOrig",
    "evtName",
    "rnd",
    "heat",
    "athOrig",
    "currentPoints",
    "endingPoints",
    "actualSco"
]
varRange = Dict([ var => [] for var in variables ] )

for wid in keys(data)
	push!( varRange[ "evtYear" ] 		, parse( data[wid][ "evtYear" ] ) )
	push!( varRange[ "evtOrig" ] 		, isoDict[ data[wid][ "evtOrig" ] ] )
	push!( varRange[ "evtName" ] 		, data[wid][ "evtName" ] )
	push!( varRange[ "rnd" ] 			, maxRnd[data[wid]["evtId"]] - parse(data[wid]["rnd"]) )
	push!( varRange[ "heat" ] 			, parse( data[wid][ "heat" ] ) )
	push!( varRange[ "athOrig" ] 		, isoDict[ data[wid]["athOrig" ] ] )
	push!( varRange[ "currentPoints" ] 	, data[wid]["currentPoints"] )
	push!( varRange[ "endingPoints" ] 	, data[wid]["endingPoints"] )
	push!( varRange[ "actualSco" ] 		, data[wid]["actualSco"] )
end

for var in variables
	sort!( unique!(varRange[var]) )
end

MartysTensor = Array( dims=(3,7,11,11,16,10,1,1,1) )


