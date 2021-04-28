using JSON
doc = JSON.parse( open("data/CleanAllDataCC.txt", "r"))

filter!(doc) do wave
	# Simplify data to most "complete years"
	# ... we can always undo this
	# subScoOrigDefect== true when no info on judge origins
	wave[2]["evtYear"] in ["2018","2019"] && wave[2]["subScoOrigDefect"]==false
end

WaveIds = sort(collect(keys(doc)))
EvtIds = Set(map(x->doc[x]["evtId"],WaveIds))
maxRnd = Dict()
for id in EvtIds
    thisEvt = filter(wid->doc[wid]["evtId"]==id, WaveIds )
    maxRnd[id] = maximum(map(wid->Base.parse(Int,doc[wid]["rnd"]),thisEvt))
end

@enum Orig AUS=1 BRA=2 ESP=3 FRA=4 PRT=5 USA=6 ZAF=7 IDN=8 ITA=9 JPN=10 NZL=11
isoDict = Dict([
    "Australia" => AUS,
    "Brazil" => BRA,
    "Spain" => ESP,
    "Basque Country" => ESP,
    "France" => FRA,
    "Indonesia" => IDN,
    "Italy" => ITA,
    "Japan" => JPN,
    "New Zealand" => NZL,
    "Portugal" => PRT,
    "United States" => USA,
    "Hawaii" => USA,
    "South Africa" => ZAF
])

EvtInd(S::String) = findfirst(==(S),["Bali Pro","Bells Beach","France","Gold Coast","J-Bay Open","Margaret River","Peniche Pro","Pipe Masters","Rio Pro","Tahiti"])
OrigInd(S::Symbol) = findfirst(x->Symbol(x)==S, instances(Orig))

# YR ∈ {2018,2019}
# EVT ∈ {"Bali Pro","Bells Beach","France","Gold Coast","J-Bay Open","Margaret River","Peniche Pro","Pipe Masters","Rio Pro","Tahiti"}
# RND ∈ {1,2,3,4,5,6,7}
# HEAT ∈ {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}
# ATH_orig ∈ {:AUS,:BRA,:FRA,:IDN,:ITA,:JPN,:NZL,:PRT,:USA,:ZAF}
# JUD_orig ∈ {:AUS,:BRA,:ESP,:FRA,:PRT,:USA,:ZAF}
# JUD_score ∈ {0.1, 0.2, ..., 9.8, 9.9, 10.0}
# JUD_cls ∈ {1,2,3,4,5}
# SZ_cls ∈ {1,2,3,4,5}
# N_clss ∈ {1,2,3,4,5}

vars = [:YR,:EVT,:RND,:HEAT,:ATH_orig,:JUD_orig,:JUD_score,:cls_1,:cls_2,:cls_3,:cls_4,:cls_5]
data = zeros(Rational,(2,10,7,16,11,7,100,2,2,2,2,2))
N = 5*length(doc)
for wave in doc
	yr = (wave[2]["evtYear"]=="2018" ? 1 : 2)
	evt = EvtInd(wave[2]["evtName"])
	rnd = 7-abs(Base.parse(Int,wave[2]["rnd"])-maxRnd[wave[2]["evtId"]])
	ht = Base.parse(Int,wave[2]["heat"])
	ath_orig = Int(isoDict[wave[2]["athOrig"]])
	#=
    for i in 1:5
		jud_orig = Int(isoDict[ wave[2]["subScoOrig"][i] ])
		jud_score = Int(round(wave[2]["subSco"][i], digits=1)*10)
		data[yr,evt,rnd,ht,ath_orig,jud_orig,jud_score] += 1
	end
    =#
    jud_scores = Float16.(wave[2]["subSco"])
    jud_origs = map(c->isoDict[c], wave[2]["subScoOrig"])
    n_clss = length(unique(jud_scores))
    for (cls,score) in enumerate(sort(unique(jud_scores)))
        I = findall(x->x==score, jud_scores)
        #jud_cls = Int(cls) # to reassure julia's type inference
        #Sym(jud_origs[I])
        #sz_cls = length(I)
        for i in I
            jud_orig = Int(jud_origs[i])
            jud_score = Int(round(jud_scores[i], digits=1)*10)
            data[yr,evt,rnd,ht,ath_orig,jud_orig,jud_score,(cls==1)+1,(cls==2)+1,(cls==3)+1,(cls==4)+1,(cls==5)+1] += 1//N
        end
    end
end