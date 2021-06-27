# Unique size accurate
include("SimplifyComp.jl")
e(c::ORIG) = vcat(zeros(Int8,Int(c)-1),1,zeros(Int8,length(JUD_ORIGS)-Int(c)))
E_C(X::Array{ORIG,1}) = sum(map(c->Int(count(==(c),X))*e(c),JUD_ORIGS))::Array{Int64,1}
N_C(X::Array{ORIG,1}) = map(c->count(==(c),X),Tuple(JUD_ORIGS))

