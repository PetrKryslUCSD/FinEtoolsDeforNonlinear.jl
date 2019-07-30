module testelasticitytensor1
using FinEtools
using FinEtoolsDeforNonlinear.MatDeforNonlinearModule: totalLagrangean2current!
using FinEtoolsDeforLinear.MatDeforModule: tens4symmto6x6t!
using LinearAlgebra
using Test

function test()
	delta = (I, J) -> I == J ? 1.0 : 0.0
	lambda = 3.3
	mu = 0.156
	F = [1.02 0.03 -0.04; 0.01 0.99 -0.03; -0.01 0.02 0.95]
	C = fill(0.0, 3, 3, 3, 3)
	for I in 1:3, J in 1:3, K in 1:3, L in 1:3
		C[I, J, K, L] = lambda * delta(I, J) * delta(K, L) + 
		mu * (delta(I, K) * delta(J, L) + delta(I, L) * delta(J, K))
	end
	Cm = fill(0.0, 6, 6)
	tens4symmto6x6t!(Cm, C)
	c = fill(0.0, 6, 6)
	totalLagrangean2current!(c, Cm, F)
	c2 = fill(0.0, 3, 3, 3, 3)
	for i in 1:3, j in 1:3, k in 1:3, l in 1:3
		c2[i, j, k, l] = 0.0
		for I in 1:3, J in 1:3, K in 1:3, L in 1:3
			c2[i, j, k, l] += C[I, J, K, L] / det(F) * F[i, I] * F[j, J] * F[k, K] * F[l, L]
		end
	end
	c2m = fill(0.0, 6, 6)
	tens4symmto6x6t!(c2m, c2)
	@test norm(c - c2m) <= 1.0e-12
end
end
using .testelasticitytensor1
testelasticitytensor1.test()