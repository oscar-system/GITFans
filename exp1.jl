using Singular
using Polymake
import Nemo
include("/home/sebastian/Software/gc-converge-julia/gap/pkg/GAPJulia/LibGAP.jl/src/initialization.jl")
GAP.run_it("/home/sebastian/Software/gc-converge-julia/gap")

######## GLOBAL SETTNGS ################


perms_list = [ [1,3,2,4,6,5,7,8,10,9], [5,7,1,6,9,2,8,4,10,3] ]
perms = []
for i in perms_list
    push!( perms, GAP.Globals.PermList(GAP.julia_to_gap(i)) )
end
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

function is_monomial_free(poly_list,vars_to_zero)
    for i in 1:nr_variables
        if !(i in vars_to_zero)
            poly_list = map(j->Singular.substitute_variable(j,i,R(0)),poly_list)
        end
    end
    perm = [ j for j in 2:(nr_variables) ]
    push!(perm,1)
    for i in 1:nr_variables
        if !(nr_variables in vars_to_zero)
            saturated = satstd(Ideal(R,poly_list),Ideal(R,T...))
            if ngens(saturated) == 1 && saturated[1] == R(1)
                return false
            elseif reduce(R(1),saturated) == 0
                return false
            end
        end
        vars_to_zero = [ perm[j] for j in vars_to_zero ]
        poly_list = [ Singular.permute_variables(j, perm, R) for j in poly_list ]
    end
    return true
end

function unique_cone_list( cone_list )
    if length(cone_list) == 1
        return cone_list
    end
    new_list = [ cone_list[ 1 ] ]
    for i in 2:length(cone_list)
        if ! any( j -> polytope.equal_polyhedra(j,cone_list[i]),new_list)
            push!(new_list,cone_list[i])
        end
    end
    return new_list
end

function matrix_unpack(matrix_packed)
    return map( i -> matrix_packed[i,:], 1:size(matrix_packed,1) )
end

## solve x * left_side = right_side
function solve_linear_system(left_side, right_side)
    left_side_unpacked = matrix_unpack(left_side)
    right_side_unpacked = matrix_unpack(right_side)
    left_gap_matrix = GAP.julia_to_gap(left_side_unpacked,Val(true))
    right_gap_matrix = GAP.julia_to_gap(right_side_unpacked,Val(true))
    result = []
    for i in 1:size(right_side,1)
        push!(result,GAP.Globals.SolutionMat(left_gap_matrix,right_gap_matrix[i]))
    end
    result_converted = map(i->GAP.gap_to_julia(Array{Rational{BigInt},1},i),result)
    return result_converted
end

function orbit( element, generators, action, compare_func )
    orb = [ element ]
    for b in orb
        for g in generators
            c = action(b, g)
            tester = any(i->compare_func(i,c),orb)
            if !tester
                push!( orb, c )
            end
        end
    end
    return orb
end

function cone_operation(cone,matrix)
    rays = convert( Array{Rational{BigInt},2}, cone.RAYS )
    new_rays = rays * transpose(matrix)
    return perlobj("Cone", Dict("INPUT_RAYS" => new_rays ) )
end

######## ALGORITHM ##################

grp = GAP.Globals.Group(GAP.julia_to_gap(perms))

set = GAP.Globals.Combinations(GAP.julia_to_gap(convert(Array{Int64,1},1:nr_variables)))

orbs = GAP.Globals.Orbits(grp,set,GAP.Globals.OnSets)

length(orbs)

reps = map(i->(orbs[i][1],length(orbs[i])),1:length(orbs))

reps_julia = map( i -> (GAP.gap_to_julia(Array{Int64,1}, i[1] ), i[2]) , reps )

mat = reshape(mat_list,nr_variables,projected_dimension)

idmat = one(Array{Rational{BigInt},2}(undef,nr_variables,nr_variables))

collector = []
collector_non_projected = []
collector_cones = []

for (i,len) in reps_julia
    current_mat = mat[i,:]
    if rank(current_mat) == projected_dimension && is_monomial_free(ideal_list,i)
    # if is_monomial_free(ideal_list,i)
        push!(collector,current_mat)
        input_dict = Dict( "INPUT_RAYS" => current_mat )
        cone = perlobj("Cone",input_dict)
        push!(collector_cones,cone)
    end
end

collector_cones_unique = unique_cone_list(collector_cones)

perms_list_projected = []
for i in perms_list
    mat_perm = mat[i,:]
    current_solution = solve_linear_system(transpose(mat),transpose(mat_perm))
    push!( perms_list_projected, copy(hcat(current_solution...)') )
end

orbit_list = []

for current_cone in collector_cones_unique
    push!(orbit_list,orbit(current_cone,perms_list_projected,cone_operation,polytope.equal_polyhedra))
end
