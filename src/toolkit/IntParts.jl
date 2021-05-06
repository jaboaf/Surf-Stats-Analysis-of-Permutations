using Combinatorics

F(z) = map(integer_partitions(z)) do λ
       c = [count(==(k),λ) for k in 1:n]
       m = factorial(n)/prod([factorial(c[k])* k^c[k] for k in 1:n])
       l = length(λ)
       return θ->m*θ^l / prod([θ+k-1 for k in 1:n])
end

n = 5
Θ = [i for i in .05:0.05:2]
map(θ->mapreduce(f->f(θ),+,F(n)), Θ)
