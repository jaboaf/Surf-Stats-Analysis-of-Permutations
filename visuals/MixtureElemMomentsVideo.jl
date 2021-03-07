include("../toolkit/SurfingInfo.jl")
using GRUtils

p = 12
g = sum(rand(RepSym(7), 30))
videofile("MomentsForMixtureElement.mp4") do
	for n in 0:p
		fig = Figure( (800,800) )
		subplot(2,2, (1) )
		wireframe(
			g^n ,
			title=""" Moment #$n """
		)
		subplot(2,2, (2) )
		wireframe(
			g^n - g^n * transpose(g),
			title="""Change from Moment #$(n-1) to Moment #$n """
		)
		subplot(2,2, (3) )
		wireframe(
			sum([g^k for k in 0:n]) ,
			title="""Sum of First $n Moments """
		)
		subplot(2,2, (4) )
		wireframe(
			sum([g^k - g^k *transpose(g) for k in 0:n]),
			title="""Sum of Changes of Moments up to $n """
		)
		for w in 1:10
			draw(gcf())
		end
	end
end

g *= 1/30
videofile("MomentsForMixtureElementStandardized.mp4") do
	for n in 0:p
		fig = Figure( (800,800) )
		subplot(2,2, (1) )
		wireframe(
			g^n ,
			title=""" Moment #$n """
		)
		subplot(2,2, (2) )
		wireframe(
			g^n - g^n * transpose(g),
			title="""Change from Moment #$(n-1) to Moment #$n """
		)
		subplot(2,2, (3) )
		wireframe(
			sum([g^k for k in 0:n]) ,
			title="""Sum of First $n Moments """
		)
		subplot(2,2, (4) )
		wireframe(
			sum([g^k - g^k *transpose(g) for k in 0:n]),
			title="""Sum of Changes of Moments up to $n """
		)
		for w in 1:10
			draw(gcf())
		end
	end
end