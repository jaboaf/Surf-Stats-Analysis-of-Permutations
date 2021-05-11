include("SimplifyComp")

rnkData = zeros(Rational,(2,10,7,16,11,7,5,5));
rnkDataVarRngs = [collect(instances(SZN)), collect(instances(EVT)), RND, HT, collect(instances(ORIG)), JUD_ORIGS,1:5, 1:5]
rnkDataVarNames = [:Season, :Event, :Round, :Heat, :AthOrig, :JudOrig, :MaxRank, :Rank]
for wave in last.(WAVES)
	szn = Int(wave.szn)
	evt = Int(wave.evt)
	rnd = wave.rnd
	ht = wave.heat
	ath_orig = Int(wave.athOrig)
	max_rank = length(wave.λ_origs)

	for (rank,part) in enumerate(wave.λ_origs)
		for orig in part
			rnkData[szn,evt,rnd,ht,ath_orig,Int(orig),max_rank,rank] += 1
		end
	end
end