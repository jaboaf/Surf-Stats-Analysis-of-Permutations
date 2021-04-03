function totVar(A::Array{T,2}) where T <: Number
    if size(A)[1] != size(A)[2] error("A is not square") end
    n = size(A)[1]
    v = 0
    for τ in S(n)
        v += abs( reduce(*,[ A[i,τ[i]] for i in 1:n]) )
    end
    return v
end

# Defining isordered on orderable collections of numbers
function isordered(X::Array{T}) where T <: Number
    return all([ X[i] < X[i+1] for i in 1:(length(X)-1) ])
end

function isordered(X::Tuple{T}) where T <: Number
    return all([ X[i] < X[i+1] for i in 1:(length(X)-1)])
end

function isordered(X::NTuple{N,<:Number}) where N
    return all([ X[i] < X[i+1] for i in 1:(N-1) ])
end

# Defining isordered on orderable collections of collections
# First given two arrays
function isordered(X::Array{<:Number},Y::Array{<:Number})
    return all(isordered.( Base.product(X,Y) ) )
end

# Given an array of arrays
function isordered(X::Array{Array{T,1},1}) where T<:Number
    return all( isordered.(Base.product(X...) ))
end




