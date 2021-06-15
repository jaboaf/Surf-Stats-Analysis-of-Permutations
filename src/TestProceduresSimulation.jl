include("toolkit/SymGrpAndReps.jl")
D = 5
aStat = [[] for i in 1:D]
uStat = [ [] for i in 1:D]
for d in 1:D
	Samp = rand(Sym(d),5000)
	U = zeros(Float32,ntuple(i->d,d))
	U[Samp[1]...] = 1
	U = SymOp(U);
	F = zeros(Float32, ntuple(i->d,d));
	for N in 1:5000
		F *= (N-1)/N
		F[ Samp[N]...] += 1/N
		a = sum(abs.( 1/2*(F-SymOp(F)) ) )
		u = sum(abs.( 1/2*(F-U)) )
		push!(aStat[d],a)
		push!(uStat[d],u)
	end
	println("done with dim= $d")
end

using GRUtils
for i in 1:D plot(aStat[i],ylim=(0,1),hold=true) end


d = 5
NumberOfSamples = 10
ObsPerSample = 1000
NotSYMsims = Array{Float32}(undef,(ObsPerSample,NumberOfSamples))
NotUsims = Array{Float32}(undef,(ObsPerSample,NumberOfSamples))
Altsims = Array{Float32}(undef,(ObsPerSample,NumberOfSamples))
Z = zeros(Float32, ntuple(i->d,d));
Z[1,2,3,4,5] = 1
U = SymOp(Z);
for k in 1:NumberOfSamples
	Samp = rand(Sym(d),ObsPerSample)
	F = zeros(Float32, ntuple(i->d,d));
	for n in 1:ObsPerSample
		F *= (n-1)/n
		F[ Samp[n]...] += 1/n
		NotSYMsims[n,k] = sum(abs.( 1/2*(F-SymOp(F)) ) )
		NotUsims[n,k] = sum(abs.( 1/2*(F-U) ) )
		Altsims[n,k] = sum(abs.( 1/2 * AltOp(F)))
	end
end

plot(mapslices(x-> x .* sqrt.(1:1000), NotSYMsims,dims=1))
plot(mapslices(x-> x .* sqrt.(1:1000), Altsims,dims=1))
plot(mapslices(x-> x .* sqrt.(1:1000), NotUsims,dims=1))


plot(hcat(Sims...))


