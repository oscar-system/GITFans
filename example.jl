############################################################################
##  The computations shown here belong to Example 5.2 in the paper at
##  https://arxiv.org/abs/1603.09241
##  (Computing GIT-Fans with Symmetry and the Mori Chamber Decomposition ...)

using GITFans
using GAP
using Polymake
using Singular


############################################################################
##  1. the input data (see the paper)

# grading matrix
Q = [
 1  1   0   0   0 ;
 1  0   1   1   0 ;
 1  0   1   0   1 ;
 1  0   0   1   1 ;
 0  1   0   0  -1 ;
 0  1   0  -1   0 ;
 0  1  -1   0   0 ;
 0  0   1   0   0 ;
 0  0   0   1   0 ;
 0  0   0   0   1 ];

# polynomial ring
R, T = Singular.PolynomialRing( Singular.QQ,
           [ "x" * string(i) for i in 1:size( Q, 1 ) ] )

# generators for the ideal
II = Singular.Ideal( R, [
    T[5]*T[10] - T[6]*T[9] + T[7]*T[8],
    T[1]*T[9]  - T[2]*T[7] + T[4]*T[5],
    T[1]*T[8]  - T[2]*T[6] + T[3]*T[5],
    T[1]*T[10] - T[3]*T[7] + T[4]*T[6],
    T[2]*T[10] - T[3]*T[9] + T[4]*T[8],
] )

# symmetries: (2,3)(5,6)(9,10), (1,5,9,10,3)(2,7,8,4,6)
perms_list = [ [1,3,2,4,6,5,7,8,10,9], [5,7,1,6,9,2,8,4,10,3] ];
perms = [ GAP.Globals.PermList(GAP.julia_to_gap(i)) for i in perms_list ];
G = GAP.Globals.Group( GAP.julia_to_gap( perms ) )


############################################################################
##  2. initialize

collector_cones = GITFans.orbit_cones( II, Q, G )
hom = GITFans.action_on_target( Q, G )
orbit_list = GITFans.orbit_cone_orbits( collector_cones, hom );

# array of homomorphisms,
# the i-th entry describes the action of G on the i-th orbit
homs = GITFans.action_on_orbit_cone_orbits( orbit_list, hom )

# the induced actions on each of the orbits
generators_new_perm = GITFans.rewrite_action_to_orbits( homs )

q_cone = Polymake.polytope.Cone( INPUT_RAYS = Q )
q_cone_rays = q_cone.RAYS
q_cone_facets = q_cone.FACETS
q_cone_facets_converted = convert(Array{Rational{BigInt},2},q_cone_facets)
q_cone_int_point = GITFans.get_interior_point(convert(Array{Rational{BigInt},2},q_cone_rays))

start_hash = GITFans.compute_bit_list(orbit_list,q_cone_int_point)
orbit_start_hash_smallest = GITFans.find_smallest_orbit_element(start_hash,generators_new_perm,GITFans.bitlist_oper_tuple,==,GITFans.less_or_equal_array_bitlist)

hash_list = [ orbit_start_hash_smallest ]


############################################################################
##  3. compute neighbors until the list is complete

current_finished_index = 1
while current_finished_index <= length(hash_list)
    global current_finished_index
    current_hash = hash_list[current_finished_index]
    current_cone_list = GITFans.cones_from_bitlist( orbit_list, current_hash );
    intersected_cone = Polymake.polytope.intersection( current_cone_list...)
    facets = intersected_cone.FACETS
    facets = convert(Array{Rational{BigInt},2},facets)
    facet_points = []
    for i in 1:size(facets,1)
        push!(facet_points, convert(Array{Rational{BigInt},1},Polymake.polytope.facet(intersected_cone,i-1).REL_INT_POINT))
    end

    neighbor_hashes = []
    for i in 1:length(facet_points)
        if any(i->i==0, q_cone_facets_converted*facet_points[i])
            continue
        end
        push!(neighbor_hashes, GITFans.get_neighbor_hash(orbit_list,facet_points[i],facets[i,:]))
    end

    neighbor_hashes = map(i->GITFans.find_smallest_orbit_element(i,generators_new_perm,GITFans.bitlist_oper_tuple,==,GITFans.less_or_equal_array_bitlist),neighbor_hashes)
    for i in neighbor_hashes
        if !(i in hash_list)
            push!(hash_list,i)
        end
    end
    current_finished_index += 1
end

