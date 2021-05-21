include("SimplifyComp.jl")
#SZN
#EVT

byEvt = partitionBy(:evt)
byEvt = map(x-> x[1] => map(y->y.λ_c, x[2]),byEvt)
P = zeros((10,7,7,7,7,7));
for by in byEvt
	evt = by[1]
	for λ in by[2]
		# note a = number of elements in ∏ Sym.(λ)...
		a = prod( factorial.(length.(λ)) )
		for t in Base.product(Sym.(λ)...)
			P[Int(evt), Int.(vcat(t...))... ] += 1/a
		end
	end
end

evtNs = [sum(P[i,:,:,:,:,:]) for i in 1:10]
for e in 1:10
	P[e,:,:,:,:,:] /= evtNs[e]
end

KS = zeros(Float32,(7,7,7,7,7));
for i in CartesianIndices(size(KS))
	t = Tuple(i)
	vals = P[:,t...]
	KS[i] = maximum(vals) - minimum(vals)
end





