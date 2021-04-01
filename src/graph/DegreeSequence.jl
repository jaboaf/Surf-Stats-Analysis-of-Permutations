t = [2,1,1,2,4,1,2,1]

function graphFromDeg(T::Array)
    G = Array[]
    t = sort(copy(T), rev=true)
    dotDiagram(t)
    println("This is the degree sequence you passed (yes, we sorted it)")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    while any(t .> 0)
        rem = findall(n-> n>0, t)
        m = minimum(t[rem])
        gs = rand(permGroup(rem), m)
        append!(G, gs)
        dotDiagram(t, forceMin=1)
        println("------------------------")
        println("sampled $m perms from Sym($rem)")
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
        t .-= m
    end
    dotDiagram(t, forceMin=1)
    println("------------------------")
    println("Nothing left? Correct.")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    println("Now lets see whats under the hood: \n")
    dotDiagram(t)
    return G
end