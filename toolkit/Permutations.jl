

struct perm{d}
	a::Array{Int8,1}
	s = prod( [ a[j] - a[i] >= 0 ? 1 : -1 for i in 1:d-1 for j in i:d ] )
end