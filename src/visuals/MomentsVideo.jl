include("../toolkit/SymGrpAndReps.jl")
using GRUtils
# Note subplot(num_rows, num_columns, indices)

p = 12
g = rand(RepSym(7))
videofile("MomentsForGenericElement.mp4") do
	for n in 0:p
		fig = Figure( (800,800) )
		subplot(2,2, (1) )
		wireframe(
			g^n ,
			title=""" Moment #$n """,
			zticks=(1,1),
			zlims=(0,1)
		)
		subplot(2,2, (2) )
		wireframe(
			g^n - g^n * transpose(g),
			title="""Change from Moment #$(n-1) to Moment #$n """,
			zticks=(1,1),
			zlims=(-1,1)
		)
		subplot(2,2, (3) )
		wireframe(
			sum([g^k for k in 0:n]) ,
			title="""Sum of First $n Moments """,
			zticks=(1,p),
			zlims=(0,p)
		)
		subplot(2,2, (4) )
		wireframe(
			sum([g^k - g^k *transpose(g) for k in 0:n]),
			title="""Sum of Changes of Moments up to $n """,
			zticks=(1,1),
			zlims=(-1,1)
		)
		for w in 1:10
			draw(gcf())
		end
	end
end