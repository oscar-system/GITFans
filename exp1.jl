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

function as_permutation(element, set, action, compare )
    perm = Array{Int64,1}(undef,length(set))
    for i in 1:length(set)
        image = action(set[i],element)
        perm[i] = findfirst( j -> compare(j,image), set )
    end
    return perm
end

function get_interior_point(matrix)
    matrix_unpacked = matrix_unpack(matrix)
    result = matrix_unpacked[1]
    for i in 2:length(matrix_unpacked)
        result += matrix_unpacked[i]
    end
    return result
end

function compute_bit_list( orbits, point )
    bitset_list = map(i->BitSet(),orbits)
    pm_point = convert(pm_Vector{pm_Rational},point)
    for current_orbit_nr in 1:length(orbits)
        current_orbit = orbits[current_orbit_nr]
        for current_cone_nr in 1:length(current_orbit)
            current_cone = current_orbit[current_cone_nr]
            if call_method(current_cone,:contains,pm_point)
                push!(bitset_list[current_orbit_nr],current_cone_nr)
            end
        end
    end
    return bitset_list
end

function get_neighbor_hash(orbits, facet_point, inner_normal_vector )
    lambda = 1024
    facet_point_bl = compute_bit_list(orbits,facet_point)
    while true
        current_point = lambda*facet_point - inner_normal_vector
        ## FIXME: compute only necessary part of BL
        current_bl = compute_bit_list(orbits,current_point)
        current_symmdiff_size = 0
        if all(i->issubset(current_bl[i],facet_point_bl[i]),1:length(facet_point_bl))
            return current_bl
        end
        println(lambda)
        lambda *= 2
    end
end

## a <= b
function less_or_equal_array_bitlist(a, b)
    for i in 1:length(a)
        if a[i].bits > b[i].bits
            return false
        end
    end
    return true
end

function bitlist_oper( bitset, perm )
    return_list = Int64[]
    for i in bitset
        push!(return_list, perm[i] )
    end
    return BitSet(return_list)
end

function bitlist_oper_tuple( bitset_tuple , perm_tuple )
    return map(i->bitlist_oper(bitset_tuple[i],perm_tuple[i]),1:length(bitset_tuple))
end

function cones_from_bitlist( cone_list, bit_list_tuple )
    return_list = Any[]
    for i in 1:length(cone_list)
        for j in bit_list_tuple[i]
            push!(return_list,cone_list[i][j])
        end
    end
    return return_list
end

function find_smallest_orbit_element(elem, gens, action, comparator, leq )
    current_orbit = orbit(elem,gens,action,comparator)
    sorted_orbit = sort(current_orbit;lt=leq)
    return sorted_orbit[1]
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

new_perm_presentation = Array{Array{Array{Int64,1},1},1}(undef,0)
for current_orbit in orbit_list
    current_presentation = []
    for current_element in perms_list_projected
        push!(current_presentation,as_permutation(current_element,current_orbit,cone_operation,polytope.equal_polyhedra))
    end
    push!(new_perm_presentation,current_presentation)
end

q_cone = perlobj( "Cone", Dict( "INPUT_RAYS" => mat ) )
q_cone_rays = q_cone.RAYS
q_cone_facets = q_cone.FACETS
q_cone_facets_converted = convert(Array{Rational{BigInt},2},q_cone_facets)

q_cone_int_point = get_interior_point(convert(Array{Rational{BigInt},2},q_cone_rays))

start_hash = compute_bit_list(orbit_list,q_cone_int_point)

generators_new_perm = map( i-> Any[], 1:length(new_perm_presentation[1]))

for i in 1:length(new_perm_presentation)
    for j in 1:length(new_perm_presentation[i])
        push!(generators_new_perm[j],new_perm_presentation[i][j])
    end
end

orbit_start_hash_smallest = find_smallest_orbit_element(start_hash,generators_new_perm,bitlist_oper_tuple,==,less_or_equal_array_bitlist)

hash_list = [ orbit_start_hash_smallest ]
current_finished_index = 1

while current_finished_index <= length(hash_list)
    global current_finished_index
    current_hash = hash_list[current_finished_index]
    current_cone_list = cones_from_bitlist( orbit_list, current_hash );
    intersected_cone = call_function(:intersection, current_cone_list...)
    facets = intersected_cone.FACETS
    facets = convert(Array{Rational{BigInt},2},facets)
    facet_points = []
    for i in 1:size(facets,1)
        push!(facet_points, convert(Array{Rational{BigInt},1},call_function(:facet,intersected_cone,i-1).REL_INT_POINT))
    end
    neighbor_hashes = []
    for i in 1:length(facet_points)
        if any(i->i==0, q_cone_facets_converted*facet_points[i])
            continue
        end
        push!(neighbor_hashes,get_neighbor_hash(orbit_list,facet_points[i],facets[i,:]))
    end
    neighbor_hashes = map(i->find_smallest_orbit_element(i,generators_new_perm,bitlist_oper_tuple,==,less_or_equal_array_bitlist),neighbor_hashes)
    for i in neighbor_hashes
        if !(i in hash_list)
            push!(hash_list,i)
        end
    end
    current_finished_index += 1
end
