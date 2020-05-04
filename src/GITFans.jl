module GITFans

# Load the necessary Julia packages.
import Singular
import Polymake
import Nemo
import GAP


#############################################################################

# utility functions

## a <= b for two arrays `a`, `b` of bitlists
function less_or_equal_array_bitlist(a, b)
    for i in 1:length(a)
        if a[i].bits > b[i].bits
            return false
        end
    end
    return true
end;

##  Return a new bitset, consisting of the images of `bitset` under the
##  permutation `perm`.
function bitlist_oper( bitset, perm )
    return_list = Int64[]
    for i in bitset
        push!(return_list, perm[i] )
    end
    return BitSet(return_list)
end;

function bitlist_oper_tuple( bitset_tuple, perm_tuple )
    return map(i->bitlist_oper(bitset_tuple[i],perm_tuple[i]),1:length(bitset_tuple))
end;

function find_smallest_orbit_element(elem, gens, action, comparator, leq )
    current_orbit = orbit(elem,gens,action,comparator)
    sorted_orbit = sort(current_orbit;lt=leq)
    return sorted_orbit[1]
end;

function rewrite_action_to_orbits( homs )
    G = GAP.Globals.Source( homs[1] )
    Ggens = GAP.Globals.GeneratorsOfGroup( G )
    generators_new_perm = map( i-> Any[], 1:length( Ggens ) )

    for i in 1:length( homs )
      n = GAP.Globals.NrMovedPoints( GAP.Globals.Image( homs[i] ) )
      if n == 0
        n = 1
      end
      for j in 1:length( generators_new_perm )
        push!( generators_new_perm[j],
               GAP.gap_to_julia( Vector{Int},
                   GAP.Globals.ListPerm( GAP.Globals.Image( homs[i],
                       Ggens[j] ), n ) ) )
      end
    end

    return generators_new_perm
end;


#############################################################################

# user functions

"""
    is_monomial_free( I, vars_to_zero = [] )

> see Prop. 3.1 in the paper
"""
function is_monomial_free( I, vars_to_zero = [] )
    R = Singular.base_ring( I )
    nr_variables = Singular.nvars( R )
    poly_list = [ I[i] for i in 1:Singular.ngens( I ) ]
    for i in 1:nr_variables
        if (i in vars_to_zero)
            poly_list = map(j->Singular.substitute_variable(j,i,R(0)),poly_list)
        end
    end

    # perm is an n-cycle
    perm = [ j for j in 2:(nr_variables) ]
    push!(perm,1)

    for i in 1:nr_variables
        if !(nr_variables in vars_to_zero)
            saturated = Singular.satstd(
                Singular.Ideal( R, poly_list ), Singular.MaximalIdeal( R, 1 ) )
            if Singular.ngens(saturated) == 1 && saturated[1] == R(1)
                return false
            elseif Singular.reduce( R(1), saturated ) == 0
                return false
            end
        end
        vars_to_zero = [ perm[j] for j in vars_to_zero ]
        poly_list = [ Singular.permute_variables(j, perm, R) for j in poly_list ]
    end
    return true
end


"""
"""
function orbit_cones( I, Q, G = GAP.Globals.SymmetricGroup( 1 ) )
    R = Singular.base_ring( I )
    nr_variables = Singular.nvars( R )
    projected_dimension = size( Q, 2 )

    set = GAP.Globals.Combinations( GAP.julia_to_gap( 1:nr_variables ) )
    orbs = GAP.Globals.Orbits( G, set, GAP.Globals.OnSets )
    reps = [ orbs[i][1] for i in 1:length(orbs) ]
    reps_julia = [ GAP.gap_to_julia( Vector{Int64}, i ) for i in reps ]

    collector_cones = []

    for i in reps_julia
      current_mat = Q[i,:]

      if Nemo.rank( current_mat ) == projected_dimension &&
         is_monomial_free( I, setdiff( 1:nr_variables, i ) )
        cone = Polymake.polytope.Cone( INPUT_RAYS = current_mat )
        if ! any( j -> Polymake.polytope.equal_polyhedra( j, cone ),
                  collector_cones )
          push!( collector_cones, cone )
        end
      end
    end

    return collector_cones
end
#T what if some projections lie in the same orbit?
#T later we expand the orbits, do we want to check this here?


"""
    action_on_target( Q, G )

> Let `Q` be a Q-matrix, and `G` be a GAP permutation group
> that describes an action on the rows of `Q`.
> The function returns the group homomorphism from `G`
> to its induced matrix action on the column space of `Q`.
> 
> Example:
> perms_list = [ [1,3,2,4,6,5,7,8,10,9], [5,7,1,6,9,2,8,4,10,3] ];
> perms = [ GAP.Globals.PermList(GAP.julia_to_gap(i)) for i in perms_list ];
> G = GAP.Globals.Group( GAP.julia_to_gap( perms ) )
> Q = [
>  1  1   0   0   0 ;
>  1  0   1   1   0 ;
>  1  0   1   0   1 ;
>  1  0   0   1   1 ;
>  0  1   0   0  -1 ;
>  0  1   0  -1   0 ;
>  0  1  -1   0   0 ;
>  0  0   1   0   0 ;
>  0  0   0   1   0 ;
>  0  0   0   0   1 ];
> GITFans.action_on_target( Q, G )
"""
function action_on_target( Q, G )
#T rename to orbit_cone_action?
    gapQ = GAP.julia_to_gap( Q )
    gapQtr = GAP.Globals.TransposedMat( gapQ )
    gens = GAP.Globals.GeneratorsOfGroup( G )
    Qimgs = GAP.Globals.List( gens, p -> GAP.Globals.Permuted( gapQ, p ) )
    Qimgstr = GAP.Globals.List( Qimgs, GAP.Globals.TransposedMat )
    genimgs = GAP.Globals.List( Qimgstr,
                  mat -> GAP.Globals.List( mat,
                             row -> GAP.Globals.SolutionMat( gapQtr, row ) ) )

    return GAP.Globals.GroupHomomorphismByImages( G,
               GAP.Globals.Group( genimgs ), gens, genimgs )
end


"""
    orbit( point, generators, action, compare_func )

> Return the orbit of the point `point` under the action of the group that is
> generated by the elements in the array `generators`.
> The function `action` defines the action (pt,g) -> pt^g on some set.
> The function `compare_func` compares two elements of the set
> and returns `true` if the two objects are considered as equal,
> and `false` otherwise.
> The elements of the returned array are in general not sorted.
"""
function orbit( point, generators, action, compare_func )
    orb = [ point ]
    # Note that in Julia (like in GAP),
    # the 'for' loop runs also over entries that get added
    # inside the loop.
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
end;


"""
    as_permutation( element, set, action, compare )

> Return the permutation (encoded as the array of image positions)
> induced by the action of the group element `element` on the array `set`
> via the function `action`: (pt,element) -> pt^element.
> The equality of points is decided via the binary function `compare`.
"""
function as_permutation( element, set, action, compare )
    perm = Vector{Int64}(undef,length(set))
    for i in 1:length(set)
        image = action(set[i],element)
        perm[i] = findfirst( j -> compare(j,image), set )
    end
    return perm
end;

function matrix_action_on_cones( cone, matrix )
    rays = convert( Matrix{Rational{BigInt}}, cone.RAYS )
    new_rays = rays * transpose( matrix )
    return Polymake.polytope.Cone( INPUT_RAYS = new_rays )
end


function orbit_cone_orbits( cones, hom; disjoint_orbits = false )

    gap_matgrp = GAP.Globals.Image( hom )
    gens = GAP.Globals.GeneratorsOfGroup( gap_matgrp )
    gens = [ GAP.gap_to_julia( Matrix{BigInt}, gens[i] )
             for i in 1:length( gens ) ]

    act = matrix_action_on_cones
    comp = Polymake.polytope.equal_polyhedra

    result = []
    for cone in cones
      orb = orbit( cone, gens, act, comp )
      if disjoint_orbits || all( o -> all( c -> ! comp( cone, c ), o ), result )
        push!( result, orbit( cone, gens, act, comp ) )
      end
    end

    return result
end


function action_on_orbit_cone_orbits( orbits, hom )
    G = GAP.Globals.Source( hom )
    gap_matgrp = GAP.Globals.Image( hom )
    gens = GAP.Globals.GeneratorsOfGroup( gap_matgrp )
    gens = [ GAP.gap_to_julia( Matrix{BigInt}, gens[i] )
             for i in 1:length( gens ) ]

    act = matrix_action_on_cones
    comp = Polymake.polytope.equal_polyhedra

    result = [ [ as_permutation( gen, orb, act, comp ) for gen in gens ]
               for orb in orbits ]

    res = []
    for list in result
      gapperms = GAP.Globals.List( GAP.julia_to_gap( list, Val(true) ),
                                   GAP.Globals.PermList )
      img = GAP.Globals.Group( gapperms )
      push!( res, GAP.Globals.GroupHomomorphismByImages( G, img,
          GAP.Globals.GeneratorsOfGroup( G ), gapperms ) )
    end

    return res
end

function matrix_unpack( matrix_packed )
    return map( i -> matrix_packed[i,:], 1:size( matrix_packed, 1 ) )
end;

"""
    get_interior_point(cone)

> Return the sum of the rows of the matrix given by the rays of the cone.
"""
function get_interior_point( cone )
    matrix = convert(Matrix{Rational{BigInt}}, cone.RAYS)
    matrix_unpacked = matrix_unpack( matrix )
    result = matrix_unpacked[1]
    for i in 2:length( matrix_unpacked )
        result += matrix_unpacked[i]
    end
    return result
end;


"""
    compute_bit_list( orbits, point )

> Let `orbits` be an array of arrays of cones.
> Return an array of bitsets,
> containing at position i the bitset for the i-th entry of `orbits`,
> that is, `true` at the positions of those cones in `orbits[i]`
> that contain the point `point`.
"""
function compute_bit_list( orbits, point )
    bitset_list = map(i->BitSet(),orbits)
    for current_orbit_nr in 1:length(orbits)
        current_orbit = orbits[current_orbit_nr]
        for current_cone_nr in 1:length(current_orbit)
            current_cone = current_orbit[current_cone_nr]
            if Polymake.polytope.contains(current_cone,point)
                push!(bitset_list[current_orbit_nr],current_cone_nr)
            end
        end
    end
    return bitset_list
end;

"""
    cones_from_bitlist( cone_list, bit_list_tuple )

> Return the array of all cones at `true` positions in the bitsets,
> where `cone_list` is an array of arrays of cones,
> and `bit_list_tuple` is an array of bitsets.
"""
function cones_from_bitlist( cone_list, bit_list_tuple )
    return_list = Any[]
    for i in 1:length(cone_list)
        for j in bit_list_tuple[i]
            push!(return_list,cone_list[i][j])
        end
    end
    return return_list
end;

"""
    get_neighbor_hash( orbits, facet_point, inner_normal_vector )

> Return the list of bitsets describing the cones adjacent to the point
> `facet_point` in the direction of `inner_normal_vector`.
"""
function get_neighbor_hash(orbits, facet_point, inner_normal_vector )
    lambda = 1024
    facet_point_bl = compute_bit_list(orbits,facet_point)
    while true
        current_point = lambda*facet_point - inner_normal_vector
        ## FIXME: compute only necessary part of BL
        current_bl = compute_bit_list(orbits,current_point)
        if all(i->issubset(current_bl[i],facet_point_bl[i]),1:length(facet_point_bl))
            return current_bl
        end
        lambda *= 2
    end
end;

"""
    GITFan(I, Q, G)

> Return a pair containing the array of orbit cone orbits and the array of
> the corresponding GAP homomorphism objects from `G` to the induced
> permutation action on the orbits,
> where `I` is the ideal in question, `Q` is the Q-matrix,
> and `G` is a symmetry group of the problem.
"""
function GITFan(I, Q, G)
    collector_cones = orbit_cones( I, Q, G )
    hom = action_on_target( Q, G )
    orbit_list = orbit_cone_orbits( collector_cones, hom; disjoint_orbits = true );
    homs = action_on_orbit_cone_orbits( orbit_list, hom )

    return (orbit_list, homs)
end;

"""
    fan_traversal( fan_descr, Q )

> Return the array of hash values ...
"""
function fan_traversal( fan_descr, Q )

    orbit_list = fan_descr[1]
    homs = fan_descr[2]

    # the induced actions on each of the orbits
    generators_new_perm = rewrite_action_to_orbits( homs )

    q_cone = Polymake.polytope.Cone( INPUT_RAYS = Q )
    q_cone_facets_converted = convert(Matrix{Rational{BigInt}}, q_cone.FACETS)
    q_cone_int_point = get_interior_point(q_cone)

    start_hash = compute_bit_list(orbit_list, q_cone_int_point)
    orbit_start_hash_smallest = find_smallest_orbit_element(start_hash,generators_new_perm,bitlist_oper_tuple,==,less_or_equal_array_bitlist)
    hash_list = [ orbit_start_hash_smallest ]

    current_finished_index = 1
    while current_finished_index <= length(hash_list)
        current_hash = hash_list[current_finished_index]
        current_cone_list = cones_from_bitlist( orbit_list, current_hash );
        intersected_cone = Polymake.polytope.intersection( current_cone_list...)
        facets = intersected_cone.FACETS
        facets = convert(Matrix{Rational{BigInt}},facets)
        facet_points = []
        for i in 1:size(facets,1)
            push!(facet_points, convert(Vector{Rational{BigInt}},Polymake.polytope.facet(intersected_cone,i-1).REL_INT_POINT))
        end

        neighbor_hashes = []
        for i in 1:length(facet_points)
            if any(x->x==0, q_cone_facets_converted*facet_points[i])
                continue
            end
            push!(neighbor_hashes, get_neighbor_hash(orbit_list,facet_points[i],facets[i,:]))
        end

        neighbor_hashes = map(i->find_smallest_orbit_element(i,generators_new_perm,bitlist_oper_tuple,==,less_or_equal_array_bitlist),neighbor_hashes)
        for i in neighbor_hashes
            if !(i in hash_list)
                push!(hash_list,i)
            end
        end
        current_finished_index += 1
    end

    return hash_list
end;

end # module
