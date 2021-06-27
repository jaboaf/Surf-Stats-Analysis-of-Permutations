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
	Rng=Set(NamedTuple{(s,)}.(last.(WAVES)))
	return sort(collect(Rng))
end

function varRng(S::NTuple{N,Symbol}) where N
	Rng=Set(NamedTuple{S}.(last.(WAVES)))
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

function partitionBy(S::NTuple{N,Symbol}) where N
	By = varRng(S)
	partition = [b => [] for b in By]
	for wave in WAVES
		i = findfirst(==(NamedTuple{S}(wave[2])), By) # note: can call field of named tuple with []
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
	wave[2]["subSco"] = map(x->round(x,digits=1),wave[2]["subSco"])
end

WAVES = []
Tups = []
for wid in WaveIds
	c = data[wid]["subScoOrig"]
	s = data[wid]["subSco"]
    panelInfo =  sort(s.=> data[wid]["subScoOrig"])
    binaryPanelInfo = sort(s .=> (data[wid]["subScoOrig"] .== data[wid]["athOrig"]))

    invlabPan = []
    binaryLabPan = [0 => [], 1 => [] ]
    for c in unique(last.(panelInfo))
    	I = findall(==(c), data[wid]["subScoOrig"])
    	push!(invlabPan, sort(s[I])=>data[wid]["subScoOrig"][I])
    	if c==data[wid]["athOrig"]
    		push!(binaryLabPan[1][2], s[I]...)
    	else
    		push!(binaryLabPan[2][2], s[I]... )
    	end
    end
    labPan = unique.(last.(invlabPan)) .=> first.(invlabPan)

    eqPan = sort(invlabPan)
    if !(isordered(first.(eqPan)))
		cond = true
		while cond
			j = findfirst(i->!(isordered(first.(eqPan)[i:(i+1)])), 1:(length(eqPan)-1))
			newElm = [ vcat(first.(eqPan)[j:(j+1)]...) => vcat(last.(eqPan)[j:(j+1)]...)]
			eqPan = [ eqPan[1:(j-1)]... , newElm... , eqPan[(j+2):end]... ]
			cond = eval( !(isordered( collect(first.(eqPan)) )) )
		end
	end
	eqPart = ntuple(i->eqPan[i][2],length(eqPan))
	eqPartBinary = map(p-> p .== data[wid]["athOrig"],eqPart)

	scos = unique(first.(panelInfo))
	part_c = ntuple(i->last.(panelInfo[findall(x->x[1]==scos[i],panelInfo)]),length(scos))
	part_b = ntuple(i->last.(binaryPanelInfo[findall(x->x[1]==scos[i],binaryPanelInfo)]),length(scos))

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
        
        judge_scores=s,
        judge_origs=data[wid]["subScoOrig"],

        panel=panelInfo,
        panelBinary=binaryPanelInfo,

        labeledPanel=labPan,
        labeledPanelBinary=binaryLabPan,

        eqPanel = eqPan,

        λ_c = part_c,
        λ_b = part_b,
        λ_eqc = eqPart,
        λ_eqb = eqPartBinary,

        m_c = map(c->count(==(c),data[wid]["subScoOrig"]), tuple(JUD_ORIGS...) ),
        m_b = map(ind->count(==(ind), last.(binaryPanelInfo)),(0,1) ),
        m_c_supp = map(c-> c in data[wid]["subScoOrig"], tuple(JUD_ORIGS...) ),
        
        panel_scores= Multiset(s),
        panel_origs= Multiset(data[wid]["subScoOrig"]),
		I_match = Int(data[wid]["athOrig"] in data[wid]["subScoOrig"])
    )
    push!(WAVES, wid=>wave)
    push!(Tups, wave)
end

⊗(A::Array{T},B::Array{T}) where T<: Number = prod.(Base.product(A,B))
⊗(V::Vararg{Array{T}}) where T<: Number =reduce(⊗,V) 
⊗(a::NTuple{T},b::NTuple{T}) where T  = (a...,b...)
⊗(a::NTuple{N,T},b::NTuple{M,T}) where {T,N,M}  = (a...,b...)
×(A::Set,B::Set) = Set(Base.product(A,B))
×(A::Array,B::Array) = collect(Base.product(A,B))
⊕(V::Vararg{Array{T,N} where N}) where T = [ sum(filter(x->size(x)==d,V)) for d in unique(size.(V)) ] 
# maybe
# ⊕(V::Vararg{Array{T,N} where N,m}) where {T,m} = 
⨁(A::Array{Array{T,N} where N,1}) where T <: Number = [ sum(filter(x->size(x)==d,A)) for d in sort(unique(size.(A))) ]
⨁(A::Array{Array{T},1}) where T<:Number = [ sum(filter(x->size(x)==d,A)) for d in sort(unique(size.(A))) ]
⨁(A::Array{Array{T,N},1}) where {T<:Number,N} = [ sum(filter(x->size(x)==d,A)) for d in sort(unique(size.(A))) ]
⨂(A::Array{Array{T,N} where N,1}) where T <: Number = reduce(⊗,A)


E(X::Array,i::K) where {K<:Integer} =dropdims( sum(X,dims=setdiff(1:ndims(X),i)),dims=tuple(setdiff(1:ndims(X),i)...) )
E(X::Array,I::NTuple) =dropdims(sum(X,dims=setdiff(1:ndims(X),I)),dims=tuple(setdiff(1:ndims(X),I)...))
cov(X::Array,i::K,j::K) where {K<:Integer} = E(X,(i,j))-E(X,i)⊗E(X,j)
cov(X::Array,I::NTuple,j::K) where {K<:Integer} = E(X,(I...,j))-E(X,I)⊗E(X,j)
cov(X::Array,i::K,J::NTuple) where {K<:Integer} = E(X,(i,J...))-E(X,i)⊗E(X,J)
cov(X::Array,I::NTuple,J::NTuple) =E(X,(I...,J...))-E(I)⊗E(J)


#corr(X,Y) = 

