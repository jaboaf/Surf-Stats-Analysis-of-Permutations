using GRUtils

# Note subplot(num_rows, num_columns, indices)

videofile("Moments.mp4") do
	for n in 0:10
		Figure()
		draw( Figure("") )
		subplot(2,2, 1)
		draw( gcf( wireframe(g^n , title="Raw Moment $n")))

		subplot(2,2, 2)
		draw( gcf( wireframe(g^n - g^n * transpose(g), title="Centralized moment")))

		subplot(2,2, 3)
		draw( gcf( wireframe( sum([g^k for k in 0:n]) , title="Cumulative Raw Moment $n")))

	    subplot(2,2, 4)
		draw( gcf( wireframe( sum([g^k - g^k *transpose(g) for k in 0:n]), title="Cumulative of first $n central moments")))

		sleep(.5)
	end
end