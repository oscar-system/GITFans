
# submodule of GITFans
module OLD

# the necessary Julia packages
import Singular
import Polymake
import Nemo
import GAP

# for using LaTex in docstrings processed by Documenter.jl
import Markdown


#############################################################################

# utility functions

## a <= b for two arrays `a`, `b` of bitsets,
## defined lexicographically
function less_or_equal_array_bitlist(a::Vector{BitSet}, b::Vector{BitSet})
    for i in 1:length(a)
        if a[i].bits < b[i].bits
            return true
        elseif a[i].bits > b[i].bits
            return false
        end
    end
    return true
end;

##  Return a new bitset, consisting of the images of `bitset` under the
##  permutation `perm`.
function bitlist_oper(bitset::BitSet, perm::Vector{Int})
    return_list = Int64[]
    for i in bitset
        push!(return_list, perm[i])
    end
    return BitSet(return_list)
end;

function bitlist_oper_tuple(bitset_tuple, perm_tuple)
    return map(i -> bitlist_oper(bitset_tuple[i], perm_tuple[i]),
               1:length(bitset_tuple))
end;

function find_smallest_orbit_element(elem, gens, action, comparator, leq)
    current_orbit = orbit(elem, gens, action, comparator)
    sorted_orbit = sort(current_orbit; lt=leq)
    return sorted_orbit[1]
end;

function rewrite_action_to_orbits(homs)
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
    is_monomial_free(I, vars_to_zero = [])

See Prop. 3.1 in [BKR](https://arxiv.org/abs/1603.09241).
"""
function is_monomial_free(I, vars_to_zero = [])
    R = Singular.base_ring( I )
    nr_variables = Singular.nvars( R )
    poly_list = [ I[i] for i in 1:Singular.ngens( I ) ]
    for i in vars_to_zero
        poly_list = map(j->Singular.substitute_variable(j, i, R(0)), poly_list)
    end

    # perm is an n-cycle
    perm = collect( 2:(nr_variables) )
    push!(perm, 1)

    for i in 1:nr_variables
        if !(nr_variables in vars_to_zero)
            saturated = Singular.satstd(
                Singular.Ideal(R, poly_list), Singular.MaximalIdeal(R, 1))
            if Singular.ngens(saturated) == 1 && saturated[1] == R(1)
                return false
            elseif Singular.reduce(R(1), saturated) == 0
                return false
            end
        end
        vars_to_zero = perm[vars_to_zero]
        poly_list = [ Singular.permute_variables(j, perm, R) for j in poly_list ]
    end
    return true
end


"""
    orbit_cones(I, Q, G = GAP.Globals.SymmetricGroup(1))

> Return orbit representatives of the group `G` on the set of those cones
> whose defining rays are given by subsets `S` of the rows of the matrix `Q`,
> such that the matrix `S` has full rank and such that the ideal `I` is
> monomial-free (see [`is_monomial_free`](@ref)) w.r.t. the variables `x_i`
> for which the `i`-th row of `Q` is not contained in `S`.
"""
function orbit_cones(I, Q, G = GAP.Globals.SymmetricGroup(1))
    nr_variables = size( Q, 1 )
    projected_dimension = size( Q, 2 )

    collector_cones = []

    # We need not consider sets of smaller size because of the rank condition.
    for k in projected_dimension:nr_variables
        set = GAP.Globals.Combinations( GAP.julia_to_gap( 1:nr_variables ), k )
        orbs = GAP.Globals.Orbits( G, set, GAP.Globals.OnSets )
        reps = [ orbs[i][1] for i in 1:length(orbs) ]
        reps_julia = [ GAP.gap_to_julia( Vector{Int64}, i ) for i in reps ]

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
    end

    return collector_cones
end
#T what if some projections lie in the same orbit?
#T later we expand the orbits, do we want to check this here?

@doc Markdown.doc"""
    action_on_target(Q, G)

Let `Q` be an $n \times m$ Q-matrix, and `G` be a GAP permutation group
on $n$ points that describes an action on the rows of `Q`.
The function returns the group homomorphism $\rho$ from `G`
to its induced matrix action on the $m$-dimensional space over the Rationals.

# Examples

```jldoctest
julia> Q = [
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

julia> perms_list = [ [1,3,2,4,6,5,7,8,10,9], [5,7,1,6,9,2,8,4,10,3] ];

julia> perms = [ GAP.Globals.PermList(GAP.julia_to_gap(i)) for i in perms_list ];

julia> G = GAP.Globals.Group( GAP.julia_to_gap( perms ) );

julia> GITFans.action_on_target( Q, G )
GAP: [ (2,3)(5,6)(9,10), (1,5,9,10,3)(2,7,8,4,6) ] -> 
[ 
  [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ], [ 0, 0, 0, 0, 1 ],
      [ 0, 0, 0, 1, 0 ] ], 
  [ [ -1, 1, -1, -1, -2 ], [ 1, 0, 1, 1, 1 ], [ 1, 0, 0, 1, 1 ], 
      [ 0, 0, 0, 0, 1 ], [ 1, 0, 1, 0, 1 ] ] ]
```

Note that $\rho(\pi)$, for $\pi \in $`G`, satisfies
`Q`$^\pi$ = `Q` * $\rho(\pi^{-1})$, where `Q`$^\pi$ is the matrix obtained
from `Q` by permuting its rows by $\pi$.

Let $Q\{J\}$ denote the matrix whose rows are those rows of $Q$ indexed by
the subset $J$ of $\{ 1, 2, \ldots, n \}$.
The cone defined by $Q{J}$ is mapped to
$Q{J} \cdot \rho(\pi^{-1}) = (Q \cdot \rho(\pi^{-1}))\{J\} = Q^\pi\{J\}$,
which is equal to $Q\{J^{\pi^{-1}}\}$.
"""
function action_on_target(Q, G)
#T rename to orbit_cone_action?
    gapQ = GAP.julia_to_gap( Q )
    gapQtr = GAP.Globals.TransposedMat( gapQ )
    gens = GAP.Globals.GeneratorsOfGroup( G )
    Qimgs = GAP.Globals.List( gens, p -> GAP.Globals.Permuted( gapQ, p ) )
    Qimgstr = GAP.Globals.List( Qimgs, GAP.Globals.TransposedMat )
    genimgs = GAP.Globals.List( Qimgstr,
                  mat -> GAP.Globals.TransposedMat( GAP.Globals.Inverse(
                             GAP.Globals.List( mat,
                             row -> GAP.Globals.SolutionMat( gapQtr, row ) ) ) ) )

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
            if ! any(i->compare_func(i, c), orb)
                push!( orb, c )
            end
        end
    end
    return orb
end;

"""
    as_permutation(element, set, action, compare)

> Return the permutation (encoded as the array of image positions)
> induced by the action of the group element `element` on the array `set`
> via the function `action`: (pt,element) -> pt^element.
> The equality of points is decided via the binary function `compare`.
"""
function as_permutation( element, set, action, compare )
    perm = Vector{Int64}(undef, length(set))
    for i in 1:length(set)
        image = action(set[i], element)
        perm[i] = findfirst(j -> compare(j, image), set)
    end
    return perm
end;

function matrix_action_on_cones( cone, matrix )
    rays = convert( Matrix{Rational{BigInt}}, cone.RAYS )
    return Polymake.polytope.Cone( INPUT_RAYS = rays * matrix )
end;

function orbit_cone_orbits(cones, hom; disjoint_orbits = false)
    gap_matgrp = GAP.Globals.Image( hom )
    gens = GAP.Globals.GeneratorsOfGroup( gap_matgrp )
    gens = [ GAP.gap_to_julia( Matrix{BigInt}, gens[i] )
             for i in 1:length( gens ) ]

    act = matrix_action_on_cones
    comp = Polymake.polytope.equal_polyhedra

    result = []
    for cone in cones
        orb = orbit(cone, gens, act, comp)
        if disjoint_orbits || all(o -> all(c -> ! comp(cone, c), o), result)
            push!(result, orbit(cone, gens, act, comp))
        end
    end

    return result
end

function action_on_orbit_cone_orbits(orbits, hom)
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
        gapperms = GAP.Globals.List(GAP.julia_to_gap(list, Val(true)),
                                    GAP.Globals.PermList)
        img = GAP.Globals.Group(gapperms)
        push!(res, GAP.Globals.GroupHomomorphismByImages(G, img,
                        GAP.Globals.GeneratorsOfGroup(G), gapperms))
    end

    return res
end


"""
    get_interior_point(cone)

Return the sum of the rows of the matrix given by the rays of the cone.
"""
get_interior_point(cone) = vec(sum(cone.RAYS, dims = 1))


"""
    compute_bit_list(orbits, point)

Let `orbits` be an array of arrays of cones.
Return an array of bitsets,
containing at position i the bitset for the i-th entry of `orbits`,
that is, `true` at the positions of those cones in `orbits[i]`
that contain the point `point`.
"""
function compute_bit_list(orbits, point)
    bitset_list = map(i->BitSet(), orbits)
    for current_orbit_nr in 1:length(orbits)
        current_orbit = orbits[current_orbit_nr]
        for current_cone_nr in 1:length(current_orbit)
            current_cone = current_orbit[current_cone_nr]
            if Polymake.polytope.contains(current_cone, point)
                push!(bitset_list[current_orbit_nr], current_cone_nr)
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
function cones_from_bitlist(cone_list, bit_list_tuple)
    return_list = Any[]
    for i in 1:length(cone_list)
        for j in bit_list_tuple[i]
            push!(return_list, cone_list[i][j])
        end
    end
    return return_list
end;


"""
    hash_to_cone(orbit_list, hash)

"""
hash_to_cone(orbit_list, hash) = Polymake.polytope.intersection(
    cones_from_bitlist(orbit_list, hash)...);


"""
    get_neighbor_hash( orbits, facet_point, inner_normal_vector )

> Return the list of bitsets describing the cones adjacent to the point
> `facet_point` in the direction of `inner_normal_vector`.
"""
function get_neighbor_hash(orbits, facet_point, inner_normal_vector)
    lambda = 1024
    facet_point_bl = compute_bit_list(orbits, facet_point)
    while true
        current_point = lambda * facet_point - inner_normal_vector
## FIXME: compute only necessary part of BL
        current_bl = compute_bit_list(orbits, current_point)
        if all(i->issubset(current_bl[i], facet_point_bl[i]), 1:length(facet_point_bl))
            return current_bl
        end
        lambda *= 2
    end
end;

"""
    orbit_cone_orbits_and_action(I, Q, G)

Return a tuple containing
the array of orbit cone orbits,
the array of the corresponding GAP homomorphism objects from `G` to the induced
permutation action on the orbits,
and the matrix `Q`.
Here `I` is the ideal in question, `Q` is the Q-matrix,
and `G` is a symmetry group of the problem.
"""
function orbit_cone_orbits_and_action(I, Q, G)
    collector_cones = orbit_cones(I, Q, G)
    hom = action_on_target(Q, G)
    orbit_list = orbit_cone_orbits(collector_cones, hom; disjoint_orbits = true);
    homs = action_on_orbit_cone_orbits(orbit_list, hom)

    return Dict( "orbit_list" => orbit_list,
                 "hom" => hom,
                 "homs" => homs,
                 "Q" => Q,
                 "G" => G )
end;

"""
    fan_traversal(oco)

Return the a pair `(hash_list, edges)` where `hash_list` is an array that
encodes orbit representatives of the maximal cones of the GIT fan described
by `oco`,
and `edges` encodes the `Set` of edges of the incidence graph of the orbits.

The input is expected to be a dictionary as computed by
[`orbit_cone_orbits_and_action`](@ref).
"""
function fan_traversal(oco)
    orbit_list = oco[ "orbit_list" ]
    homs = oco[ "homs" ]
    Q = oco[ "Q" ]

    # the induced actions on each of the orbits
    generators_new_perm = rewrite_action_to_orbits(homs)

    q_cone = Polymake.polytope.Cone(INPUT_RAYS = Q)
    q_cone_facets_converted = convert(Matrix{Rational{BigInt}}, q_cone.FACETS)
    q_cone_int_point = get_interior_point(q_cone)

    start_hash = compute_bit_list(orbit_list, q_cone_int_point)
    orbit_start_hash_smallest = find_smallest_orbit_element(
        start_hash, generators_new_perm, bitlist_oper_tuple, ==,
        less_or_equal_array_bitlist)
    hash_list = [ orbit_start_hash_smallest ]

    current_pos = 0
    edges = Set( Vector{Int}[] )

    for current_hash in hash_list
        current_pos = current_pos + 1

        # note that we run also over elements added inside the loop
        current_cone_list = cones_from_bitlist( orbit_list, current_hash )
        intersected_cone = Polymake.polytope.intersection( current_cone_list...)
        facets = intersected_cone.FACETS
        facets = convert(Matrix{Rational{BigInt}},facets)
        facet_points = Vector{Rational{BigInt}}[
                           Polymake.polytope.facet(intersected_cone, i-1).REL_INT_POINT
                           for i in 1:size(facets, 1)]

        neighbor_hashes = []
        for i in 1:length(facet_points)
            if any(x->x==0, q_cone_facets_converted*facet_points[i])
                continue
            end
            push!(neighbor_hashes, get_neighbor_hash(orbit_list, facet_points[i], facets[i, :]))
        end

        neighbor_hashes = map(i->find_smallest_orbit_element(i, generators_new_perm, bitlist_oper_tuple, ==, less_or_equal_array_bitlist), neighbor_hashes)
        for i in neighbor_hashes
            if i in hash_list
                # perhaps we have found a new incidence
                push!(edges, sort([findfirst(x->x == i, hash_list), current_pos]))
            else
                # new representative found
                push!(hash_list, i)
                # new incidence
                push!(edges, [current_pos, length(hash_list)])
            end
        end
    end

    return (hash_list, edges)
end;


orbits_of_maximal_GIT_cones(oc, hash_list) = orbit_cone_orbits(
    map(x -> hash_to_cone(oc[ "orbit_list" ], x), hash_list), oc[ "hom" ]);


"""
    hashes_to_polyhedral_fan(oc, hash_list)

...
"""
function hashes_to_polyhedral_fan(oc, hash_list)
    # translate the descriptions of the orbit repres. of maximal cones
    # to cone objects
    result_cones = map(x -> hash_to_cone(oc[ "orbit_list" ], x), hash_list)

    # expand their orbits
    expanded = orbit_cone_orbits(result_cones, oc[ "hom" ])
    maxcones = vcat( expanded... )

    # the defining rays for all maximal cones
    rays_maxcones = [ [ convert( Vector{Rational{BigInt}}, cone.RAYS[i, :] )
                        for i in 1:size( cone.RAYS, 1 ) ]
                      for cone in maxcones ]

    # the set of rays
    allrays = sort( collect( Set( vcat( rays_maxcones... ) ) ) )

    # the indices of rays that belong to each maximal cone (0-based)
    index_maxcones = [ sort( [ findfirst( x -> x == v, allrays )-1
                               for v in rays ] )
                       for rays in rays_maxcones ]

    return Polymake.fan.PolyhedralFan(INPUT_RAYS = hcat( allrays...)',
                                      INPUT_CONES = index_maxcones)
end


"""
    git_fan(a, Q, G)

Return the polymake object that represents the polyhedral fan given by
the ideal a, the grading matrix Q, and the symmetry group G.
"""
function git_fan(a, Q, G)
    oc = orbit_cone_orbits_and_action(a, Q, G)
    (hash_list, edges) = fan_traversal(oc)

    return hashes_to_polyhedral_fan(oc, hash_list)
end


"""
    edges_intersection_graph(maxcones, inter_dim)

Return the array of those pairs `[i, j]` such that `i < j` and such that
the intersection of the cones `maxcones[i]` and `maxcones[j]` has dimension
`inter_dim`.

If `maxcones` is an array of cones of dimension `inter_dim + 1`
then the returned array describes the edges of the intersection graph.
"""
function edges_intersection_graph(maxcones, inter_dim)
    edges = Vector{Int}[]
    for j in 1:length( maxcones )
        for i in 1:(j-1)
            if Polymake.polytope.dim( Polymake.polytope.intersection(
                   maxcones[i], maxcones[j] ) ) == inter_dim
                push!( edges, [ i, j ] )
            end
        end
    end

    return edges
end;

end
