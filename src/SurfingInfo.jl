using JSON: parse
include("toolkit/Multisets.jl")
include("toolkit/SymGrpAndReps.jl")

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

#KeyVars = ["evtYears", "evtOrig", "evtName", "rnd", "heat"]
#Index = Dict([ var => Dict() for var in [ ] ])
#Ind["evtYear"]

data = parse( open("Data/CleanAllDataCC.txt", "r"))

waves = []
WaveIds = []
for wid in sort(collect(keys(data)))
    if data[wid]["nJudOrigs"] == 5 & data[wid]["nSubScos"] == 5
        origs = unique(data[wid]["subScoOrig"])
        matchIndicator = (data[wid]["athOrig"] in origs)

        labeledScos = Dict([isoDict[origin] => Float16[] for origin in origs])
        labeledScosBinary = Dict([:Match => Float16[], :NoMatch => Float16[] ])

        origScoPairs = zip(data[wid]["subScoOrig"], data[wid]["subSco"])
        for p in origScoPairs
            # Psuedo: push!( array of judge scores from country p[1], score=p[2] )
            push!(labeledScos[ isoDict[p[1]] ], p[2])
            if p[1] == data[wid]["athOrig"]
                push!(labeledScosBinary[:Match], p[2])
            else
                push!(labeledScosBinary[:NoMatch], p[2])
            end
        end

        wave = (
            id=wid,
            evtYear=data[wid]["evtYear"],
            evtOrig=isoDict[data[wid]["evtOrig"]],
            evtName=data[wid]["evtName"],
            evtId=data[wid]["evtId"],
            rnd=data[wid]["rnd"],
            rndId=data[wid]["rndId"],
            heat=data[wid]["heat"],
            heatId=data[wid]["heatId"],
            athName=data[wid]["athName"],
            athId=data[wid]["athId"],
            athOrig=isoDict[data[wid]["athOrig"]],
            currentPoints=data[wid]["currentPoints"],
            endingPoints=data[wid]["endingPoints"],
            panel=labeledScos,
            panelBinary=labeledScosBinary,
            subScos=data[wid]["subSco"],
            judgeOrigs= Multiset(map(x->isoDict[x], data[wid]["subScoOrig"])),
            panelOrigs= Set(map(x->isoDict[x], data[wid]["subScoOrig"])),
            I_match = matchIndicator
        )
        push!(waves, wave)
        push!(WaveIds, wid)
    end
end
