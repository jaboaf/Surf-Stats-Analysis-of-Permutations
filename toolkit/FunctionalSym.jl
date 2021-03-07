include("SymGrpAndReps.jl")
F = Function[]
S(n) = [ x -> return p[x] for p in Sym(n) ] 