@testset "is_monomial_free" begin

    # The following input gave a wrong `true` result
    # with an earlier version of `is_monomial_free`.
    using Oscar
    nr_variables = 4
    vars_strings = map( i -> "x"*string(i), 1:nr_variables )
    R, T = PolynomialRing(QQ,vars_strings)
    ideal_gens = [
        T[1]^2*T[2]*T[3]^2*T[4]^2+T[1]^2*T[2]^2*T[3],
        2*T[1]*T[2]^2*T[3]^2*T[4]+T[1]*T[3]^2*T[4]^2-T[1]^2,
        T[1]*T[2]*T[3]*T[4]-T[1]^2*T[4]^2+T[1]*T[2]-T[4],
        T[1]*T[2]*T[3]-T[2]*T[3]*T[4],
        T[1]^2*T[2]^2*T[4]^2-T[1]^2*T[3]^2*T[4]^2+T[1]*T[2]*T[3]*T[4]
    ]

    @test ! GITFans.is_monomial_free(ideal(R, ideal_gens))
    
end
