using Singular
include("/home/sebastian/Software/gc-converge-julia/gap/pkg/GAPJulia/LibGAP.jl/src/initialization.jl")
GAP.run_it("/home/sebastian/Software/gc-converge-julia/gap")

######## GLOBAL SETTNGS ################

perms = []
push!( perms, GAP.Globals.PermList(GAP.julia_to_gap([1,3,2,4,6,5,7,8,10,9])) )
push!( perms, GAP.Globals.PermList(GAP.julia_to_gap([5,7,1,6,9,2,8,4,10,3])) )
mat_list = [ 1,1,1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,0,1,1,0,0,0,-1,1,0,0,0,1,0,1,0,-1,0,0,1,0,0,0,1,1,-1,0,0,0,0,1 ]
nr_variables = 10
projected_dimension = 5

vars_strings = map( i -> "x"*string(i), 1:nr_variables )

R,T = PolynomialRing(QQ,vars_strings)

ideal_list = [
    T[5]*T[10] - T[6]*T[9] + T[7]*T[8],
    T[1]*T[9]  - T[2]*T[7] + T[4]*T[5],
    T[1]*T[8]  - T[2]*T[6] + T[3]*T[5],
    T[1]*T[10] - T[3]*T[7] + T[4]*T[6],
    T[2]*T[10] - T[3]*T[9] + T[4]*T[8],
]

function is_monomial_free(ideal_list,vars_to_zero)
    for i in vars_to_zero
        ideal_list = map(j->Singular.substitute_ith_variable(j,i,R(0)),ideal_list)
    end
    saturated = satstd(Ideal(R,ideal_list),Ideal(R,T...))
    return ngens(saturated) == 1 && saturated[1] == R(1)
end


######## ALGORITHM ##################

grp = GAP.Globals.Group(GAP.julia_to_gap(perms))

set = GAP.Globals.Combinations(GAP.julia_to_gap(convert(Array{Int64,1},1:nr_variables)))

orbs = GAP.Globals.Orbits(grp,set,GAP.Globals.OnSets)

length(orbs)

reps = map(i->orbs[i][1],1:length(orbs))

reps_julia = map( i -> GAP.gap_to_julia(Array{Int64,1}, i ), reps )

mat = transpose(reshape(mat_list,nr_variables,projected_dimension))

collector = []
for i in reps_julia
    current_mat = mat[:,i]
    if rank(current_mat) == projected_dimension && !is_monomial_free(ideal_list,i)
        push!(collector,current_mat)
    end
end





