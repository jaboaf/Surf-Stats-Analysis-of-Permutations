panels = Array{Pair{Float16,Array{Symbol,1}}}[]
histogram(length.(panels), title="Distribution of Number of Subspaces")


#= CURRENT VARS AVAILABLE
x "subScoOrigDefect"	true for some but not many
x "subScoDefect"		false for  2018,2019

x "nSubScos"			always 5 for 2018,2019
x "nJudOrigs"			
x "validSubScos"		always 5 gor filtered

x "matchMean" 			actual stat
x "noMatchMean"			actual stat
x "matchVar" 			actual stat
x "noMatchVar"			actual stat
x "actualScoVar"		actual stat
x "subScoMean"			actual stat
x "subScoVar"			actual stat

x "matchSubScos"	manipulation of data gen process
x "noMatchSubScos" 	manipulation of data gen process

x "atHome"			actual indicator var
x "noMatches"		actual indicator var
x "matches"			actual rand. variable

external information?
	"actualSco"			computation
	"actualScoLevel"	f: 0:0.1:10 --> Poor:Excellent
	"endingPoints"		f: Athletes --> [ numer. rng ]
	"currentPoints"		f: Athlete, Time --> [ numer. rng ]

Actual Variables:
"evtYear"
"evtOrig"		
"evtName"
"evtId"

"rnd"
"rndId"

"heat"
"heatId"

"athName"
"athId"
"athOrig"

"subSco"
"subScoOrig"
=#

# IDEA!!!
# let sym be: SymGrp{TypeA,N}(::Set{Typea})
# where typea isa TypeA
# and len of set is N

#= Simplify
"evtYear" --> 
	"2018"	--> WCT18
	"2019"	-->	WCT19
"evtName" --> ?
	"Gold Coast"
 	"Bells Beach"
 	"Bali Pro"
 	"Margaret River"
 	"Tahiti"
 	"J-Bay Open"
 	"Rio Pro"
 	"Peniche Pro"
 	"France"
 	"Pipe Masters"

"rnd" # we have 11 of these
	"Final"
	"Semifinal"
	"Quarterfinal"

"heat" # we have 16 of these

"athName"

=#

"evtYear"		

"evtName"

"rnd"

"heat"

"athOrig"

"subScoOrig"

"rank"

varRng("evtName")
varRng("rnd")
varRng("heat")
athOrig


info = BitArray()


@enum WCT begin
    WCT17=1
    WCT18
    WCT19
end
@enum EVENT begin
	"Bali Pro"			--> BaliPro=1
	"Bells Beach" 		--> BellsBeach
	"France"			--> FrancePro
	"Gold Coast" 		--> GoldCoast
	"J-Bay Open"		--> JBayOpen
	"Margaret River"	--> MargaretRiver
	"Peniche Pro"		--> PenichePro
	"Pipe Masters"		--> PipeMasters
	"Rio Pro"			--> RioPro
	"Tahiti" 			--> Teahupoo
end



struct TypArray{T,N} <: Array{T,N}

