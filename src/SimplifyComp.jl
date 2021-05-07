using Statistics
using JSON
include("toolkit/Multisets.jl")
include("toolkit/OrderingUtils.jl")
include("toolkit/SymGrpAndReps.jl")

@enum ORIG AUS=1 BRA=2 ESP=3 FRA=4 PRT=5 USA=6 ZAF=7 IDN=8 ITA=9 JPN=10 NZL=11
@enum SZN WCT18=1 WCT19
@enum EVT BaliPro=1 BellsBeach FrancePro GoldCoast JBayOpen MargaretRiver PenichePro PipeMasters RioPro Teahupoo
const RND = [i for i in 1:7]
const HT = [i for i in 1:16]
const JUD_ORIGS = collect(instances(ORIG)[1:7])
evtDict = Dict([
	"Bali Pro" => BaliPro,
	"Bells Beach" => BellsBeach,
	"France" => FrancePro,
	"Gold Coast" => GoldCoast,
	"J-Bay Open" => JBayOpen,
	"Margaret River" => MargaretRiver,
	"Peniche Pro" => PenichePro,
	"Pipe Masters" => PipeMasters,
	"Rio Pro" => RioPro,
	"Tahiti" => Teahupoo
])
isoDict = Dict([
    "Australia" => AUS,
    "Brazil" => BRA,
    "Basque Country" => ESP,
    "Spain" => ESP,
    "France" => FRA,
    "Portugal" => PRT,
    "Hawaii" => USA,
    "United States" => USA,
    "South Africa" => ZAF,
    "Indonesia" => IDN,
    "Italy" => ITA,
    "Japan" => JPN,
    "New Zealand" => NZL
])
sznDict = Dict(["2018"=>WCT18, "2019"=>WCT19])

# varRng takes in element of field name s of a wave ∈ WAVES and returns sorted {wave.s | wave ∈ WAVES}
function varRng(s::Symbol)
	Rng=Set()
	for wave in WAVES
		push!(Rng, wave[2][s])
	end
	return sort(collect(Rng))
end

function partitionBy(s::Symbol)
	By = varRng(s)
	partition = [b => [] for b in By]
	for wave in WAVES
		i = findfirst(==(wave[2][s]), By) # note: can call field of named tuple with []
		push!(partition[i][2], wave[2] )
	end
	return partition
end

data = JSON.parse( open("data/CleanAllDataCC.txt", "r") )
filter!(data) do wave
	wave[2]["evtYear"] in ["2018","2019"] && wave[2]["subScoOrigDefect"]==false
end

WaveIds = sort(collect(keys(data)))
EvtIds = Set(map(x->data[x]["evtId"],WaveIds))
HeatIds = Set(map(x->data[x]["heatId"],WaveIds))
maxRnd = Dict()
for id in EvtIds
	thisEvt = filter(wid->data[wid]["evtId"]==id, WaveIds )
	maxRnd[id] = maximum(map(wid->Base.parse(Int,data[wid]["rnd"]), thisEvt))
end

# Modifications / transformations
for wave in data
	wave[2]["evtYear"] = sznDict[wave[2]["evtYear"]]
	wave[2]["evtOrig"] = isoDict[ wave[2]["evtOrig"] ]
	wave[2]["evtName"] = evtDict[ wave[2]["evtName"] ]
	wave[2]["rnd"] = maxRnd[wave[2]["evtId"]]-Base.parse(Int,wave[2]["rnd"])+1
	wave[2]["heat"] = Base.parse(Int, wave[2]["heat"])
	wave[2]["athOrig"] = isoDict[ wave[2]["athOrig"] ]
	wave[2]["subScoOrig"] = map(x->isoDict[x], wave[2]["subScoOrig"])
	wave[2]["subSco"] = Float16.(wave[2]["subSco"])
end

WAVES = []
for wid in WaveIds
    panelInfo =  sort(data[wid]["subSco"] .=> data[wid]["subScoOrig"])
    binaryPanelInfo = sort(data[wid]["subSco"] .=> (data[wid]["subScoOrig"] .== data[wid]["athOrig"]))

    labPan = []
    binaryLabPan = Dict([:Match => [], :NoMatch => [] ])
    for c in unique(last.(panelInfo))
    	I = findall(==(c), data[wid]["subScoOrig"])
    	push!(labPan, sort(data[wid]["subSco"][I])=>data[wid]["subScoOrig"][I])
    	if c==data[wid]["athOrig"]
    		binaryLabPan[:Match] = data[wid]["subSco"][I]
    	else
    		push!(binaryLabPan[:NoMatch], data[wid]["subSco"][I]... )
    	end
    end
    eqPan = sort(labPan)
    if !(isordered(first.(eqPan)))
		cond = true
		while cond
			j = findfirst(i->!(isordered(first.(eqPan)[i:(i+1)])), 1:(length(eqPan)-1))
			newElm = [ vcat(first.(eqPan)[j:(j+1)]...) => vcat(last.(eqPan)[j:(j+1)]...)]
			eqPan = [ eqPan[1:(j-1)]... , newElm... , eqPan[(j+2):end]... ]
			cond = eval( !(isordered( collect(first.(eqPan)) )) )
		end
	end

    partition_binary = []
    partition_origs = []
	for s in sort(unique(data[wid]["subSco"]))
		I = findall(x-> x==s, data[wid]["subSco"])
		push!(partition_binary, data[wid]["subScoOrig"][I] .== data[wid]["athOrig"] )
		push!(partition_origs, data[wid]["subScoOrig"][I])
	end
	scos = unique(first.(panelInfo))
	part = ntuple(i->last.(panelInfo[findall(x->x[1]==scos[i],panelInfo)]),length(scos))


    wave = (
        id=wid,

        szn=data[wid]["evtYear"],
        evtOrig=data[wid]["evtOrig"],
        evt=data[wid]["evtName"],
        evtId=data[wid]["evtId"],

        rnd=data[wid]["rnd"],
        rndId=data[wid]["rndId"],

        heat=data[wid]["heat"],
        heatId=data[wid]["heatId"],

        athName=data[wid]["athName"],
        athId=data[wid]["athId"],
        athOrig=data[wid]["athOrig"],

        currentPoints=data[wid]["currentPoints"],
        endingPoints=data[wid]["endingPoints"],
        
        judge_scores=data[wid]["subSco"],
        judge_origs=data[wid]["subScoOrig"],

        panel=panelInfo,
        panelBinary=binaryPanelInfo,

        labeledPanel=labPan,
        labeledPanelBinary=binaryLabPan,

        eqPanel = eqPan,

        λ_origs = partition_origs,
        λ_binary = partition_binary,
        λ_c = part,
        
        panel_scores= Multiset(data[wid]["subSco"]),
        panel_origs= Multiset(data[wid]["subScoOrig"]),
		I_match = Int(data[wid]["athOrig"] in data[wid]["subScoOrig"])
    )
    push!(WAVES, wid=>wave)
end
#=
rnkData = zeros(Rational,(2,10,7,16,11,7,5,5));
rnkDataVarRngs = [collect(instances(SZN)), collect(instances(EVT)), RND, HT, collect(instances(ORIG)), JUD_ORIGS,1:5, 1:5]
rnkDataVarNames = [:Season, :Event, :Round, :Heat, :AthOrig, :JudOrig, :MaxRank, :Rank]
for wave in last.(WAVES)
	szn = Int(wave.szn)
	evt = Int(wave.evt)
	rnd = wave.rnd
	ht = wave.heat
	ath_orig = Int(wave.athOrig)
	max_rank = length(wave.λ_origs)

	for (rank,part) in enumerate(wave.λ_origs)
		for orig in part
			rnkData[szn,evt,rnd,ht,ath_orig,Int(orig),max_rank,rank] += 1
		end
	end
end
=#

⊗(A::Array{T},B::Array{T}) where T<: Number = prod.(Base.product(A,B))
⊗(a::NTuple{T},b::NTuple{T}) where T  = (a...,b...)
×(A::Set,B::Set) = Set(Base.product(A,B))
×(A::Array,B::Array) = collect(Base.product(A,B))

#=
Hts = map(x->x[1]=>map(y->(round(mean(y.judge_scores),digits=2),y.λ_origs),x[2]),partitionBy(:heatId));
sort!(Hts)
# Heat x Mean judge score x panel
# since judge score is rounded to 2 decimal places, sco range needs to be 1000
e_(i::Int) = vec([j==i ? 1 : 0 for j in 1:7])
Mi = [zeros(7) for r in 1:5]
Mu = [zeros(7) for r in 1:5]
C = [zeros(7) for r in 1:5]
H = Array{Float32,7}(undef,(910,1000,7,7,7,7,7));
for (i,ht) in enumerate(Hts)
	for (s,λ) in ht[2]
		sco = Int(round(s*100))
		c = prod(factorial.(length.(λ)))
		b = [1/factorial(length(λ[k])) for k in 1:length(λ) for r in 1:length(λ[k])]
		C .+= map(cls->sum(e_.(cls)),λ)
		for a in Base.product( Sym.(λ)...)
			I = Int.(vcat(a...))
			Mu .+= (1/c)^(1/5) .* e_.(I)
			Mi .+=  b .* e_.(I)
			H[i,sco, I... ] += 1/c
		end
	end
end
=#
#=
# C for "by Country"
C = sum([M[i]⊗M[j] / (sum(M[i])*sum(M[j])) for i in 1:5 for j in 1:5])

# Is this the dual of M ?
K = [[M[r][c] for r in 1:5] for c in 1:7]
Ki = [[Mi[r][c] for r in 1:5] for c in 1:7]
Ku = [[Mu[r][c] for r in 1:5] for c in 1:7]
# R for "by Rank"
R = sum([K[i]⊗K[j] / (sum(K[i])*sum(K[j])) for i in 1:7 for j in 1:7])


#det(C)
#det(R)
println("Approximately 0...")

sum(C)
sum(R)
println("Interesting hmmmmm")

C_n= C/25
R_n= R/49

U_C = [1/7*1/7 for c in 1:7,c2 in 1:7]
U_R = [1/5*1/5 for r in 1:5,r2 in 1:5]

A = [M[r][c] for r in 1:5, c in 1:7]

ht1 = map(p-> sum([Rep(vcat(t...)) for t in Base.product(Sym.(map(cls->Int.(cls),p[2].λ_origs))...)]), filter(w->w[2].heatId=="79998",WAVES) )
ht1P = sum(map(p->p/(sum(p)/5),ht1)) # 5 is number of panel origs
# note: ht1P above, row and column sums are == 26 which is len of heat
ht1P = sum(map(p->p/(sum(p)/5),ht1))/length(ht1)

ht2 = map(p-> sum([Rep(vcat(t...)) for t in Base.product(Sym.(map(cls->Int.(cls),p[2].λ_origs))...)]), filter(w->w[2].heatId=="68816",WAVES) )
ht2P = sum(map(p->p/(sum(p)/5),ht2))
# 5 is number of panel origs, 4 is number of origs
# note: ht2P above, column sums are == 17 which is len of heat
# note: ht2P above, row sums are 17 for every row except 1st
# 1st row corresponds to AUS which two jugdes have orig of
ht2P = sum(map(p->p/(sum(p)/5),ht2))/length(ht2)

now see:
sum(ht1P,dims=(1))
sum(ht1P,dims=(2))
sum(ht2P,dims=(1))
sum(ht2P,dims=(2))

also this is what you'd expect:
sum(ht1P*ht1P',dims=(1))
sum(ht1P*ht1P',dims=(2))

now you'd expect this (because there are 2 judges from AUS):
sum(ht2P*ht2P',dims=(2))
would you expect this? (i didn't at frist but it makes sense cause weve made it symmetric):
sum(ht2P*ht2P',dims=(1))
=#